"""Align module library."""

from __future__ import annotations

import csv
import gzip
import logging
import tempfile
import traceback
from functools import partial
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Any, BinaryIO, Iterable, Iterator, Literal, NamedTuple, Sequence, cast, overload

try:
    # Try the modern location first
    from typing import Self
except ImportError:
    # Fallback to the extension library
    from typing_extensions import Self

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta
from pyhmmer import hmmalign, hmmsearch
from pyhmmer.easel import AA, DNA, Alphabet, DigitalSequence, DigitalSequenceBlock, SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

from ..exception import EmptyWarning, SeqtypeError
from ..external import Muscle
from . import SeqTypes, _abc
from ._utils import guess_seqtype, is_gzip_file, load_msa

__all__ = [
    "SampleSeqs",
    "SampleList",
    "SearchHit",
    "OrthologSeqs",
    "OrthologList",
    "run_hmmsearch",
    "run_hmmalign",
    "run_muscle",
    "bp_mrtrans",
]
logger = logging.getLogger(__name__)

_abc.FileWrapperABC.register(HMM)


class HMMMarkerSet(_abc.DataListABC[HMM]):  # type: ignore[type-var]
    """Manage a collection of pyhmmer HMM profiles with optional cutoffs for model-specific bitscore thresholds.

    Provides methods to initialize, access, and manipulate HMM data.
    """

    _bound_class = HMM

    @overload
    def __init__(self, data: str | Path) -> None: ...
    @overload
    def __init__(self, data: str | Path, cutoff_file: str | Path) -> None: ...
    def __init__(self, data: str | Path | None = None, cutoff_file: str | Path | None = None) -> None:
        """Initialize the HMMMarkerSet with data and optional cutoffs.

        Args:
            data (str | Path | None): HMM profiles or paths to load.
            cutoff_file (str | Path | None): File with model-specific bitscore cutoffs.

        Raises:
            TypeError: If an item cannot be converted to an HMM.
            RuntimeError: If the cutoff_file is not a valid file.
        """
        data_tuple = ()
        if isinstance(data, str):
            data = Path(data)
        if isinstance(data, Path):
            if data.is_dir():
                data_tuple = tuple(data.iterdir())
            else:
                data_tuple = (data,)

        self._data: list[HMM] = []
        if not data_tuple and cutoff_file:
            raise RuntimeError("Cannot specify cutoff_file without data.")

        for d in data_tuple:
            with HMMFile(d) as hmm_profile:  # type: ignore[type-var]
                profile_name = d.stem
                hmm = hmm_profile.read()
                if hmm:
                    hmm.name = profile_name
                    self.append(hmm)

        if cutoff_file:
            self.set_cutoffs(cutoff_file)

    @overload
    def __getitem__(self, key: int) -> HMM: ...
    @overload
    def __getitem__(self, key: str) -> HMM: ...
    @overload
    def __getitem__(self, key: slice) -> Self: ...
    def __getitem__(self, key: int | slice | str) -> HMM | Self:
        """Retrieve HMM profile(s) by index or name.

        Args:
            key (int | slice | str): Index, slice, or name of the profile.

        Returns:
            HMM | HMMMarkerSet: The HMM profile or a subset of profiles.

        Raises:
            KeyError: If the specified name is not found.
        """
        if isinstance(key, str):
            if key not in self.names:
                raise KeyError(f"{key}: Sample not found.")
            return self._data[self.names.index(key)]
        return super().__getitem__(key)

    def have_cutoffs(self, verbose: bool = True) -> bool:
        """Check if all HMM profiles have model-specific bitscore cutoffs.

        Args:
            verbose (bool): Log warnings if cutoffs are missing.

        Returns:
            bool: True if all profiles have cutoffs, otherwise False.
        """
        no_cutoffs = [hmm.name for hmm in self if not hmm.cutoffs.trusted]
        if no_cutoffs:
            if verbose:
                logger.warning(
                    "The following markers do not contain model-specific bitscore cutoff: %s",
                    ", ".join(no_cutoffs),
                )
            return False
        if verbose:
            logger.info("All markers contain model-specific bitscore cutoff.")
        return True

    @overload
    def set_cutoffs(self, cutoffs: str | Path) -> None: ...
    @overload
    def set_cutoffs(self, cutoffs: dict[str, float]) -> None: ...
    def set_cutoffs(self, cutoffs: str | Path | dict[str, float]) -> None:
        """Set model-specific bitscore cutoffs for HMM profiles.

        Args:
            cutoffs_dict (dict[str, float] | None): Dictionary of cutoffs.
        """
        cutoffs_dict: dict[str, float] = {}
        if isinstance(cutoffs, (str, Path)):
            cutoff_file = Path(cutoffs)
            if not cutoff_file.exists():
                raise FileNotFoundError(f"{cutoff_file}")
            if not cutoff_file.is_file():
                raise FileNotFoundError(f"{cutoff_file} is not a file.")
            with open(cutoff_file) as f:
                for line in csv.reader(f, delimiter="\t"):
                    if line[0].startswith("#"):
                        continue
                    cutoffs_dict[line[0]] = float(line[1])
        else:
            try:
                cutoffs_dict = dict(cutoffs)
            except TypeError:
                raise TypeError("Invalid cutoffs format.")

        for hmm, cutoff in cutoffs_dict.items():
            self[hmm].cutoffs.trusted = (float(cutoff),) * 2


class SampleSeqs(_abc.SeqFileWrapperABC):
    """Converts the plain or bgzipped peptide/CDS fasta file to pyhmmer compatible format.

    This class mainly used for storing sequences for hmmsearch.
    """

    @overload
    def __init__(self, file: str | Path, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...
    @overload
    def __init__(self, file: str | Path, name: str, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...
    def __init__(self, file: str | Path, name: str | None = None, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None:
        """Initialize a SampleSeqs object.

        Initialize with a path of a fasta file and a representative name (optional). The basename of the file will be used as the
        name if no name specified. The fasta won't be loaded immediately unless calling the load method.

        Args:
            file (str | Path): The path to the fasta file (plain or bgzipped).
            name (str, optional): A representative name for the sequences. Defaults to None.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Examples:
            Load a peptide fasta:
            >>> SampleSeqs("data/Actinomucor_elegans_CBS_100.09.aa.fa")

            Load a peptide fasta with a given name:
            >>> SampleSeqs("data/Actinomucor_elegans_CBS_100.09.aa.fa", "Actinomucor_elegans")

            Load a bgzf compressed cds fasta:
            >>> SampleSeqs("data/Actinomucor_elegans_CBS_100.09.cds.fa.gz", "Actinomucor_elegans")
        """
        super().__init__(file, name, seqtype=seqtype)

    @_abc.check_loaded
    def __len__(self) -> int:
        """Return the number of sequences contain in this object.

        Returns:
            int: The number of sequences in this object.
        """
        return len(self._data)

    @_abc.check_loaded
    def __iter__(self) -> Iterator[DigitalSequence]:
        """Returns an iterator over the sequences.

        Returns:
            Iterator[DigitalSequence]: An iterator for the sequences.
        """
        return iter(self._data)

    def load(self) -> None:
        """Load the sequences for hmmsearch.

        Raises:
            SeqtypeError: If RNA sequences are detected.
        """
        self._data = DigitalSequenceBlock(Alphabet.amino(), [])
        seqtype_conversion = {
            SeqTypes.PEP: Alphabet.amino(),
            SeqTypes.DNA: Alphabet.dna(),
            SeqTypes.RNA: Alphabet.rna(),
        }
        f = gzip.open(self.file) if is_gzip_file(self.file) else open(self.file, "rb")
        with SequenceFile(cast(BinaryIO, f), digital=True, alphabet=seqtype_conversion[self.seqtype]) as sf:
            seqblock: DigitalSequenceBlock = sf.read_block()
        f.close()

        if self.seqtype == SeqTypes.PEP:
            logger.debug("%s contains %s sequences.", self.file, self.seqtype)
            self._process_pep_seqs(cast("DigitalSequenceBlock[AA]", seqblock))

        elif self.seqtype == SeqTypes.DNA:
            logger.debug("%s contains %s sequences.", self.file, self.seqtype)
            self._process_cds_seqs(cast("DigitalSequenceBlock[DNA]", seqblock))

        else:
            raise SeqtypeError(f"{self.file} contains rna sequences, which is not supported. Please convert them to dna first.")

    @_abc.load_data
    def search(self, hmms: HMMMarkerSet, *, evalue: float = 1e-10, threads: int = 1) -> list[SearchHit]:
        """Run the hmmsearch process using the pyhmmer library.

        This supports multithreaded running. Empirical testing shows that efficiency drops off after more than 4 CPU threads are
        used so considering limit the threads up to 4.

        Args:
            hmms (HMMMarkerSet): A set of HMM profiles.
            evalue (float): E-value threshold for filtering. Defaults to 1e-10.
            threads (int): Number of threads to use. Defaults to 1.

        Returns:
            list[SearchHit]: A list of search hits.
        """
        return run_hmmsearch(self, hmms, evalue=evalue, threads=threads)

    def _guess_seqtype(self) -> Literal["dna", "pep", "rna", "NaN"]:
        """Guess and return the sequence type."""
        f = gzip.open(self.file, "rt") if is_gzip_file(self.file) else open(self.file)

        for r in SeqIO.FastaIO.SimpleFastaParser(f):
            seqtype = guess_seqtype(r[1], ignore_failed=True)
            if seqtype:
                break
        f.close()
        return seqtype

    def _process_pep_seqs(self, seqblock: DigitalSequenceBlock[AA]) -> None:
        """Process the peptide sequences."""
        while seqblock:
            seq = seqblock.pop()
            seq.description = self.name
            self._data.append(seq)

    def _process_cds_seqs(self, seqblock: DigitalSequenceBlock[DNA]) -> None:
        """Process the CDS sequences.

        Logs information about problematic sequences with invalid lengths.
        """
        problematic_seqs_name = []
        original_size = len(seqblock)
        while seqblock:
            seq = seqblock.pop(0)
            seq.description = self.name
            try:
                self._data.append(seq.translate())
            except ValueError:
                problematic_seqs_name.append(seq.name)

        if problematic_seqs_name:
            logger.warning(
                "%s: %s out of %s seqs have invalid length.",
                self.file,
                len(problematic_seqs_name),
                original_size,
            )
            logger.info("Problematic seqs: %s", ", ".join(problematic_seqs_name))


class SampleList(_abc.SeqDataListABC[SampleSeqs]):
    """A wrapper that stores all the SampleSeqs for an analysis."""

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, data: Sequence[str | Path | SampleSeqs]) -> None: ...
    @overload
    def __init__(self, data: Sequence[str | Path | SampleSeqs], names: Sequence[str]) -> None: ...
    @overload
    def __init__(self, data: Sequence[str | Path | SampleSeqs], *, seqtype: Literal["dna", "pep", "AUTO"]) -> None: ...
    @overload
    def __init__(
        self, data: Sequence[str | Path | SampleSeqs], names: Sequence[str], *, seqtype: Literal["dna", "pep", "AUTO"]
    ) -> None: ...
    @overload
    def __init__(
        self, data: Sequence[str | Path | SampleSeqs], names: Sequence[str], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO"
    ) -> None: ...
    def __init__(
        self,
        data: Sequence[str | Path | SampleSeqs] = (),
        names: Sequence[str] = (),
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ) -> None:
        """Initializes the object and stores data into a list.

        Args:
            data (Sequence[str | Path | SampleSeqs] | None, optional): A sequence of data items.
            names (Sequence[str] | None, optional): A sequence of names corresponding to the data items.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to a SampleSeqs.
            KeyError: If the item already exists.
        """
        if not hasattr(self, "_bound_class"):
            self._bound_class = SampleSeqs
        super().__init__(data, names, seqtype=seqtype)

    @overload
    def __getitem__(self, key: int) -> SampleSeqs: ...
    @overload
    def __getitem__(self, key: str) -> SampleSeqs: ...
    @overload
    def __getitem__(self, key: slice) -> Self: ...
    def __getitem__(self, key: int | slice | str) -> SampleSeqs | Self:
        """Retrieves an item or subset of items by name, index, or slice.

        Args:
            key (str | int | slice): The key to retrieve.

        Returns:
            SampleSeqs | SampleList: The corresponding item or subset of items.
        """
        return super().__getitem__(key)

    @_abc.load_data
    def search(self, hmms: HMMMarkerSet, *, evalue: float = 1e-10, jobs: int = 1, threads: int = 1) -> list[SearchHit]:
        """Search for HMM markers in each sample using multiple processes or threads.

        Args:
            hmms (HMMMarkerSet): The set of HMM markers to search for.
            evalue (float, optional): E-value threshold for HMM search. Defaults to 1e-10.
            jobs (int, optional): Number of parallel processes. Defaults to 1.
            threads (int, optional): Number of threads per process. Defaults to 1.

        Returns:
            SearchHitsManager: The manager containing search results.

        Raises:
            RuntimeError: If the number of jobs is outside the valid range.

        Examples:
            >>> sample_list = SampleList(["sample1.fasta", "sample2.fasta"])
            >>> hmmset = HMMMarkerSet(["marker1.hmm", "marker2.hmm"])
            >>> result = sample_list.search(hmmset, evalue=1e-5, jobs=2)
        """
        if not 0 < jobs <= len(self):
            raise RuntimeError(f"jobs = {jobs}: jobs should be between 1 and {len(self)}")

        func = partial(_search_helper, hmms=hmms, evalue=evalue, threads=threads)

        step_size = min(max(1, len(self) // 200 * 10), 50)
        total = len(self)
        search_res = []
        try:
            # if logger.parent.getEffectiveLevel() > 10:
            #     logger.setLevel("WARNING")
            if jobs <= 1:
                logger.debug("Sequential mode with %s threads.", threads)
                results = map(func, self)
                for i, r in enumerate(results, 1):
                    search_res.append(r)
                    if i % step_size == 0 or i == total:
                        logger.info("Progress: %d / %d", i, total)
            else:
                logger.debug("Multiprocesses mode with %s jobs and %s threads for each.", jobs, threads)
                with ThreadPool(processes=jobs) as pool:
                    results = pool.imap(func, self)
                    for i, r in enumerate(results, 1):
                        search_res.append(r)
                        if i % step_size == 0 or i == total:
                            logger.info("Progress: %d / %d", i, total)
        except ValueError:
            logger.error("Input sequence reading error. Did you assign the wrong seqtype?\n%s", traceback.format_exc())
            raise
        except Exception:
            logger.error("%s", traceback.format_exc())
            raise

        search_res = [hit for res in search_res for hit in res]
        return search_res


class SearchHit(NamedTuple):
    """A named tuple representing a search hit.

    Attributes:
        hmm (str): The HMM identifier associated with the search hit.
        sample (str): The name or identifier of the sample.
        sequence (str): The sequence associated with the search hit.
    """

    hmm: str
    sample: SampleSeqs
    sequence: str


class SearchHitsManager:
    """A list-like object that saves the hmmsearch hits and facilitates sequence retrieval."""

    __slots__ = ("_data", "_orthologs", "_samples", "_mfa_dir")

    def __init__(self, hits: list[SearchHit] | None = None) -> None:
        """Initialize the manager with optional hits and sample list.

        Args:
            hits (list[SearchHit] | None): Optional list of search hits.
        """
        hits = hits or []
        self._data: list[SearchHit] = []
        self._orthologs: dict[str, list[int]] = {}  # Based on hits
        self._samples: dict[SampleSeqs, list[int]] = {}  # Based on hits
        self._mfa_dir: tempfile.TemporaryDirectory | None = None

        for hit in hits:
            self.add(hit)

    def __del__(self) -> None:
        """Clean up when the object is destroyed."""
        self.unload()

    def __repr__(self) -> str:
        """Return a string representation of the manager.

        Returns:
            str: The string representation of the manager.
        """
        return type(self).__qualname__ + f"(hits={len(self)}; samples={len(self._samples)}; orthologs={len(self._orthologs)}"

    @_abc.check_class
    def __eq__(self, other: object) -> bool:
        """Check if this manager is equal to another manager.

        Args:
            other (SearchHitsManager): Another SearchHitsManager to compare.

        Returns:
            bool: True if the managers are equal, False otherwise.
        """
        if not isinstance(other, type(self)):
            return NotImplemented
        return (self._orthologs == other._orthologs) and (self._samples == other._samples)

    @overload
    def __getitem__(self, key: int) -> SearchHit: ...
    @overload
    def __getitem__(self, key: slice) -> Self: ...
    @overload
    def __getitem__(self, key: list) -> Self: ...
    def __getitem__(self, key: int | slice | list) -> SearchHit | Self:
        """Retrieve search hits using an index, slice, or list of indices.

        Args:
            key (int | slice | list): Index, slice, or list of indices to retrieve.

        Returns:
            SearchHit | SearchHitsManager: The corresponding hits.
        """
        if isinstance(key, slice):
            return self.__class__(self._data[key])
        elif isinstance(key, list):
            return self.__class__([self._data[x] for x in key])
        else:
            return self._data[key]

    def __len__(self) -> int:
        """Get the number of hits in the manager.

        Returns:
            int: The number of hits.
        """
        return len(self._data)

    def __iter__(self) -> Iterator[SearchHit]:
        """Iterate over the search hits.

        Returns:
            Iterator[SearchHit]: An iterator over the hits.
        """
        return iter(self._data)

    def __getstate__(self) -> dict:
        """Prepare the object for pickling.

        Returns:
            dict: The state of the object.
        """
        self.samplelist.unload()
        self.unload()
        return {slot: getattr(self, slot) for slot in self.__slots__}

    def __setstate__(self, state: dict[str, Any]):
        """Restore the object from a pickled state.

        Args:
            state (dict[str, Any]): The state to restore.
        """
        for slot, value in state.items():
            setattr(self, slot, value)

    def add(self, hit: SearchHit) -> None:
        """Add a new search hit to the manager.

        Args:
            hit (SearchHit): The search hit to add.
        """
        idx = len(self._data)
        self._data.append(hit)
        self._orthologs.setdefault(hit.hmm, []).append(idx)
        self._samples.setdefault(hit.sample, []).append(idx)

    def update(self, hits: Iterable[SearchHit]) -> None:
        """Update new search hits to the manager.

        Args:
            hit (Iterator[SearchHit]): The search hit to add.
        """
        for hit in hits:
            self.add(hit)

    @property
    def samplelist(self) -> SampleList:
        """Get the associated SampleList.

        Returns:
            SampleList[SampleSeqs]: The associated SampleList.
        """
        return SampleList(tuple(self._samples.keys()))

    @property
    def orthologs(self) -> dict[str, Path]:
        """Get the orthologs.

        Returns:
            dict[str, tuple | Path]: The orthologs.
        """
        if self._mfa_dir:
            mfa_dir_path = Path(self._mfa_dir.name)
            return {hmm: mfa_dir_path / f"{hmm}.fa" for hmm in self._orthologs.keys()}
        else:
            raise RuntimeError("Please run the load method first.")

    @overload
    def filter(self, *, min_taxa: int) -> Self: ...
    @overload
    def filter(self, *, drop_samples: Sequence) -> Self: ...
    def filter(self, *, min_taxa: int = 0, drop_samples: Sequence[str] = ()) -> Self:
        """Filter hits by minimum taxa or dropping specific samples.

        Args:
            min_taxa (int, optional): Minimum number of taxa. Defaults to 0.
            drop_samples (Sequence[str] | None, optional): Samples to drop. Defaults to None.

        Returns:
            SearchHitsManager: A new manager with filtered hits.

        Raises:
            EmptyWarning: If no hits are left after filtering.
        """
        drop_samples = tuple(set(drop_samples))
        selected_idx = [x for x in range(len(self))]
        if min_taxa and drop_samples:
            return self.filter(drop_samples=drop_samples).filter(min_taxa=min_taxa)
        elif drop_samples:
            selected_idx = [idx for idx, data in enumerate(self) if data.sample.name not in drop_samples]
        elif min_taxa:
            selected_idx = sorted([idx for indices in self._orthologs.values() for idx in indices if len(indices) >= min_taxa])
        if not selected_idx:
            raise EmptyWarning("No hit left after filtering.")

        return self[selected_idx]

    def load(self) -> None:
        """Retrieve the sequences for each alignment."""
        if not self._mfa_dir:
            self._mfa_dir = tempfile.TemporaryDirectory(prefix="phyling_align_mfa_")
        try:
            with tempfile.TemporaryDirectory() as faidx_dir:
                faidx_dir = Path(faidx_dir)
                msa_dir_path = Path(self._mfa_dir.name)
                loaded_seqs: list[SeqRecord] = [None] * len(self)  # type: ignore
                for sample, indices in self._samples.items():
                    self._to_seqrecord(sample, indices, loaded_seqs, faidx_dir)
                for hmm, indices in self._orthologs.items():
                    self._write_to_file(loaded_seqs, indices, msa_dir_path / f"{hmm}.fa")
        except Exception:
            self.unload()
            raise

    def unload(self) -> None:
        """Release the memory allocated by loaded sequences."""
        if self._mfa_dir:
            self._mfa_dir.cleanup()
            self._mfa_dir = None

    def _to_seqrecord(self, sample: SampleSeqs, indices: list[int], loaded_seqs: list, path: Path) -> None:
        """Convert hits to sequence records.

        Args:
            sample (str): Sample name.
            indices (list[int]): Indices of the hits.
            loaded_seqs (list): A fixed-size list that used to store the SeqRecord objects.
            path (Path): Path for temporary storage.
        """
        try:
            seqs_handler = Fasta(
                sample.file,
                indexname=Path(path) / f"{sample.name}.fai",
                gzi_indexname=Path(path) / f"{sample.name}.gzi",
            )
            for idx in indices:
                data = self._data[idx]
                seqrec = SeqRecord(
                    seq=Seq(str(seqs_handler[data.sequence])),
                    id=sample.name,
                    name=sample.name,
                    description=data.sequence,
                )
                loaded_seqs[idx] = seqrec
        except Exception as e:
            raise type(e)(f"{sample.file}: {e}").with_traceback(e.__traceback__)

    def _write_to_file(self, seqrecords: list[SeqRecord], indices: list[int], path: Path) -> None:
        """Write sequence records to a file.

        Args:
            seqrecords (list[SeqRecord]): The sequence records.
            indices (list[int]): Indices of the records to write.
            path (Path): Path to the output file.
        """
        ortholog_seqs = [seqrecords[idx] for idx in indices]
        SeqIO.write(ortholog_seqs, path, "fasta")


class OrthologSeqs(SampleSeqs):
    """Converts the plain or bgzipped peptide/CDS fasta file to pyhmmer compatible format.

    This class mainly used for storing sequences that are considered orthologs and provide alignment features.
    """

    __slots__ = ("_data_cds",)

    @overload
    def __init__(self, file: str | Path, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...
    @overload
    def __init__(self, file: str | Path, name: str, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...
    def __init__(self, file: str | Path, name: str | None = None, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None:
        """Initialize a OrthologSeqs object.

        Initialize with a path of a fasta file and a representative name (optional). The basename of the file will be used as the
        name if no name specified. The fasta won't be loaded immediately unless calling the load method.

        Args:
            file (str | Path): The path to the fasta file (plain or bgzipped).
            name (str, optional): A representative name for the sequences. Defaults to None.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Examples:
            Load a peptide fasta:
            >>> OrthologSeqs("data/101133at4751.fa")

            >>> Load a peptide fasta with a given name:
            OrthologSeqs("data/101133at4751.fa", "hmm_101133at4751")

            Load a bgzf compressed cds fasta:
            >>> OrthologSeqs("data/101133at4751.fa.bgzf", "hmm_101133at4751")
        """
        if name:
            super().__init__(file, name, seqtype=seqtype)
        else:
            super().__init__(file, seqtype=seqtype)

    def load(self) -> None:
        """Load the sequences for hmmalign.

        Raises:
            SeqtypeError: If RNA sequences are detected.
        """
        self._data_cds: list[SeqRecord] = []
        super().load()

    def _process_pep_seqs(self, seqblock: DigitalSequenceBlock[AA]) -> None:
        """Process the peptide sequences."""
        while seqblock:
            seq = seqblock.pop()
            self._data.append(seq)

    def _process_cds_seqs(self, seqblock: DigitalSequenceBlock[DNA]) -> None:
        """Process the CDS sequences.

        The CDS sequences will being saved to internal attribute _data_cds. Logs information about problematic sequences with
        invalid lengths.
        """
        problematic_seqs_name = []
        original_size = len(seqblock)
        while seqblock:
            seq = seqblock.pop(0)
            try:
                self._data.append(seq.translate())
                self._data_cds.append(
                    SeqRecord(
                        Seq(seq.textize().sequence),
                        id=seq.name,
                        name=seq.name,
                        description=seq.description,
                    )
                )
            except ValueError:
                problematic_seqs_name.append(seq.name)

        if problematic_seqs_name:
            logger.warning(
                "%s: %s out of %s seqs have invalid length.",
                self.file,
                len(problematic_seqs_name),
                original_size,
            )
            logger.debug("Problematic seqs: %s", ", ".join(problematic_seqs_name))

    def search(self, hmms: HMMMarkerSet, *, evalue: float = 1e-10, threads: int = 1) -> list[SearchHit]:
        """Placeholder for the search method in the superclass.

        Raises:
            NotImplementedError: When execute the method.
        """
        raise NotImplementedError(f"Method is not implemented in {type(self).__qualname__} class.")

    @overload
    def align(self, method: Literal["hmmalign"], *, hmm: HMM) -> MultipleSeqAlignment: ...
    @overload
    def align(self, method: Literal["muscle"], *, threads: int = 1) -> MultipleSeqAlignment: ...
    @_abc.load_data
    def align(
        self,
        method: Literal["hmmalign", "muscle"],
        *,
        hmm: HMM | None = None,
        threads: int = 1,
    ) -> MultipleSeqAlignment:
        """Align sequences from each sample using hmmalign or Muscle.

        Args:
            method (Literal["hmmalign", "muscle"]): The alignment method to use.
                - "hmmalign": Uses HMM alignment with provided marker set.
                - "muscle": Uses Muscle alignment.
            hmm (HMM): The HMM model for "hmmalign". Required for this method.
            threads (int, optional): Number of threads for "muscle". Defaults to 1.

        Returns:
            MultipleSeqAlignment: The resulting multiple sequence alignment.

        Example:
            >>> ortholog = OrthologSeqs("ortholog.fasta")
            >>> hmm = HMM("model.hmm")

            Align sequences using hmmalign:

            >>> alignment = ortholog.align(method="hmmalign", hmm)

            Align sequences using muscle:
            >>> alignment = ortholog.align(method="muscle", threads=4)
        """
        if method == "hmmalign":
            if hmm:
                pep_msa = run_hmmalign(self, hmm)
            else:
                raise ValueError('hmmalign required "hmm" argument.')
        elif method == "muscle":
            pep_msa = run_muscle(self, threads=threads)
        else:
            raise ValueError('Argument method only accepts "hmmalign" or "muscle".')
        if self.seqtype == SeqTypes.PEP:
            return pep_msa
        else:
            return bp_mrtrans(pep_msa, self._data_cds)


class OrthologList(_abc.SeqDataListABC[OrthologSeqs]):
    """A wrapper that stores all the OrthologSeqs for an analysis."""

    @overload
    def __init__(self, data: Sequence[str | Path | OrthologSeqs], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...
    @overload
    def __init__(
        self, data: Sequence[str | Path | OrthologSeqs], names: Sequence[str], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO"
    ) -> None: ...
    def __init__(
        self,
        data: Sequence[str | Path | OrthologSeqs] = (),
        names: Sequence[str] = (),
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ) -> None:
        """Initializes the object and stores data into a list.

        Args:
            data (Sequence[str | Path | OrthologSeqs] | None, optional): A sequence of data items.
            names (Sequence[str], optional): A sequence of names corresponding to the data items.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to a OrthologSeqs.
            KeyError: If the item already exists.
        """
        if not hasattr(self, "_bound_class"):
            self._bound_class = OrthologSeqs
        super().__init__(data, names, seqtype=seqtype)

    @overload
    def __getitem__(self, key: int) -> OrthologSeqs: ...
    @overload
    def __getitem__(self, key: str) -> OrthologSeqs: ...
    @overload
    def __getitem__(self, key: slice) -> Self: ...
    def __getitem__(self, key: int | slice | str) -> OrthologSeqs | Self:
        """Retrieves an item or subset of items by name, index, or slice.

        Args:
            key (str | int | slice): The key to retrieve.

        Returns:
            OrthologSeqs | OrthologList: The corresponding item or subset of items.
        """
        return super().__getitem__(key)

    @overload
    def align(self, method: Literal["hmmalign"], *, hmms: HMMMarkerSet) -> list[MultipleSeqAlignment]: ...
    @overload
    def align(self, method: Literal["hmmalign"], *, hmms: HMMMarkerSet, jobs: int) -> list[MultipleSeqAlignment]: ...
    @overload
    def align(self, method: Literal["hmmalign"], *, hmms: HMMMarkerSet, jobs: int, threads) -> list[MultipleSeqAlignment]: ...
    @overload
    def align(self, method: Literal["muscle"]) -> list[MultipleSeqAlignment]: ...
    @overload
    def align(self, method: Literal["muscle"], *, jobs: int) -> list[MultipleSeqAlignment]: ...
    @overload
    def align(self, method: Literal["muscle"], *, jobs: int, threads: int) -> list[MultipleSeqAlignment]: ...
    def align(
        self,
        method: Literal["hmmalign", "muscle"],
        *,
        hmms: HMMMarkerSet | None = None,
        jobs: int = 1,
        threads: int = 1,
    ) -> list[MultipleSeqAlignment]:
        """Align the orthologous sequences using multiple processes or threads.

        Args:
            method (Literal["hmmalign", "muscle"]): The alignment method to use.
                - "hmmalign": Uses HMM alignment with provided marker set.
                - "muscle": Uses Muscle alignment.
            hmms (HMMMarkerSet, optional): The set of HMM markers for "hmmalign". Required for this method.
            jobs (int, optional): Number of parallel processes. Defaults to 1.
            threads (int, optional): Number of threads per process. Defaults to 1.

        Returns:
            list[MultipleSeqAlignment]: A list of `MultipleSeqAlignment` objects containing alignment results.

        Raises:
            RuntimeError: If the number of jobs is outside the valid range.

        Examples:
            Align sequences using hmmalign with a set of HMM markers:

            >>> ortho_list = OrthologList(["marker1.fasta", "marker2.fasta"])
            >>> hmmset = HMMMarkerSet(["marker1.hmm", "marker2.hmm"])
            >>> msa = ortho_list.align(method="hmmalign", hmms=hmmset, jobs=2)

            Align sequences using muscle:

            >>> msa = ortho_list.align(method="muscle", jobs=2, threads=4)
        """
        if not 0 < jobs <= len(self):
            raise RuntimeError(f"jobs = {jobs}: jobs should be between 1 and {len(self)}")

        func = partial(_align_helper, hmms=hmms, method=method, threads=threads)

        step_size = min(max(10, len(self) // 200 * 50), 500)
        total = len(self)
        msa = []

        try:
            if jobs <= 1:
                logger.debug("Sequential mode with %s threads.", threads)
                results = map(func, self)
                for i, r in enumerate(results, 1):
                    msa.append(r)
                    if i % step_size == 0 or i == total:
                        logger.info("Progress: %d / %d", i, total)
            else:
                logger.debug("Multiprocesses mode with %s jobs and %s threads for each.", jobs, threads)
                with ThreadPool(processes=jobs) as pool:
                    results = pool.imap(func, self)
                    for i, r in enumerate(results, 1):
                        msa.append(r)
                        if i % step_size == 0 or i == total:
                            logger.info("Progress: %d / %d", i, total)
        except Exception:
            logger.error("%s", traceback.format_exc())
            raise
        return msa


@_abc.check_loaded
def run_hmmsearch(sample: SampleSeqs, hmms: HMMMarkerSet, *, evalue: float = 1e-10, threads: int = 1) -> list[SearchHit]:
    """Run hmmsearch on a sample sequence against a set of HMM markers.

    Args:
        sample (SampleSeqs): The sample sequences to search within.
        hmms (HMMMarkerSet): The set of HMM markers to search for.
        evalue (float, optional): E-value threshold for HMM search. Defaults to 1e-10.
        jobs (int, optional): Number of parallel processes. Defaults to 1.
        threads (int, optional): Number of threads per process. Defaults to 1.

    Returns:
        list[SearchHit]: A list of search hits, each representing the best match for a marker.

    Raises:
        RuntimeError: If the checkpoint fails to load properly.

    Example:
        >>> sample = SampleSeqs("sample.fasta")
        >>> hmmset = HMMMarkerSet(["marker.hmm"])
        >>> hits = run_hmmsearch(sample, hmmset, evalue=1e-5, threads=4)
    """
    if hmms.have_cutoffs(verbose=False):
        cutoffs_args: dict[str, Any] = {"bit_cutoffs": "trusted"}
    else:
        cutoffs_args: dict[str, Any] = {"E": evalue}
    r = []
    for hits in hmmsearch(hmms, sample, cpus=threads, **cutoffs_args):
        hmm = hits.query.name
        reported: list[SearchHit] = []
        idx = 0
        while idx < min(2, len(hits)):
            hit = hits[idx]
            if hit.reported:
                reported.append(SearchHit(hmm, sample, hit.name))
            idx += 1
        if len(reported) > 1:
            continue
        r.extend(reported)

    logger.debug("Hmmsearch on %s is done.", sample.name)
    return r


@_abc.check_loaded
def run_hmmalign(ortholog: OrthologSeqs, hmm: HMM) -> MultipleSeqAlignment:
    """Align a set of protein sequences against an HMM using hmmalign.

    Args:
        ortholog (OrthologSeqs): The ortholog sequences to align.
        hmm (HMM): The HMM model for alignment.

    Returns:
        MultipleSeqAlignment: The resulting multiple sequence alignment.

    Example:
        >>> ortholog = OrthologSeqs("ortholog.fasta")
        >>> hmm = HMM("model.hmm")
        >>> alignment = run_hmmalign(ortholog, hmm)
    """
    f = tempfile.NamedTemporaryFile()
    try:
        hmmalign(hmm, ortholog).write(cast(BinaryIO, f), "afa")
        f.seek(0)
        alignment = load_msa(Path(f.name))
        alignment.annotations = {"seqtype": SeqTypes.PEP}
        logger.debug("Hmmalign on %s is done.", ortholog.name)
        return alignment
    finally:
        f.close()


@_abc.check_loaded
def run_muscle(ortholog: OrthologSeqs, *, threads: int = 1) -> MultipleSeqAlignment:
    """Run Muscle to perform multiple sequence alignment.

    Args:
        ortholog (OrthologSeqs): The ortholog sequences to align.
        threads (int, optional): Number of threads for parallel processing. Defaults to 1.

    Returns:
        MultipleSeqAlignment: The resulting multiple sequence alignment.

    Raises:
        BinaryNotFoundError: If the Muscle binary is not found.
        RuntimeError: If Muscle fails to complete successfully.

    Example:
        >>> ortholog = OrthologSeqs("ortholog.fasta")
        >>> alignment = run_muscle(ortholog, threads=4)
    """
    with tempfile.NamedTemporaryFile() as f_aln:
        with tempfile.NamedTemporaryFile() as f_pep:
            SeqIO.write(
                (
                    SeqRecord(
                        Seq(seq.textize().sequence),
                        id=seq.name,
                        name=seq.name,
                        description=ortholog.name,
                    )
                    for seq in ortholog
                ),
                f_pep.name,
                "fasta",
            )
            f_pep.seek(0)
            runner = Muscle(f_pep.name, f_aln.name, threads=threads)
            try:
                runner.run()
            except Exception as e:
                raise RuntimeError(f"Muscle failed on peptide fasta translated from {ortholog.file}:\n{e}")
        f_aln.seek(0)
        alignment = load_msa(Path(f_aln.name))
    alignment.annotations = {"seqtype": SeqTypes.PEP}
    logger.debug("Muscle on %s is done.", ortholog.name)
    return alignment


def bp_mrtrans(pep_msa: MultipleSeqAlignment, cds_recs: list[SeqRecord]) -> MultipleSeqAlignment:
    """Transform protein MSA to mRNA/cDNA coordinates.

    Args:
        pep_msa (MultipleSeqAlignment): Protein sequence alignment.
        cds_rec (list[SeqRecord]): List of corresponding nucleotide sequences.

    Returns:
        MultipleSeqAlignment: Alignment in cDNA coordinates.

    Example:
        >>> cds_msa = bp_mrtrans(pep_msa, cds_records)
    """
    stop_codons = {"TAA", "TAG", "TGA"}
    codon_size = 3
    gap = "-"
    codon_gap = codon_size * gap

    cds_msa = MultipleSeqAlignment([], annotations={"seqtype": SeqTypes.DNA})
    for pep_rec, cds_rec in zip(pep_msa, cds_recs):
        pep_rec: SeqRecord
        dna_idx = 0
        dna_align: list[str] = []
        info = {
            "id": cds_rec.id if hasattr(cds_rec, "id") else pep_rec.id,
            "name": cds_rec.name if hasattr(cds_rec, "name") else pep_rec.name,
            "description": cds_rec.description if hasattr(cds_rec, "description") else pep_rec.description,
        }
        cds_seq: str = str(cds_rec.seq).upper()
        for align_idx in range(len(pep_rec)):
            codon = cds_seq[dna_idx : dna_idx + codon_size]

            if pep_rec[align_idx] == gap or dna_idx >= len(cds_seq):
                dna_align.append(codon_gap)
                if codon not in stop_codons:
                    continue
            else:
                dna_align.append(codon)
            dna_idx += codon_size
        cds_msa.append(SeqRecord(Seq(("".join(dna_align))), **info))
    return cds_msa


def _search_helper(
    instance: SampleSeqs,
    hmms: HMMMarkerSet,
    evalue: float = 1e-10,
    threads: int = 1,
) -> list[SearchHit]:
    """Helper function to perform hmmsearch.

    Args:
        instance (SampleSeqs): The sample sequences to search within.
        hmms (HMMMarkerSet): The set of HMM markers to search for.
        evalue (float, optional): E-value threshold for HMM search. Defaults to 1e-10.
        threads (int, optional): Number of threads per process. Defaults to 1.

    Returns:
        list[SearchHit]: A list of search hits, each representing the best match for a marker.
    """
    r = instance.search(hmms=hmms, evalue=evalue, threads=threads)
    return r


def _align_helper(
    instance: OrthologSeqs,
    hmms: HMMMarkerSet,
    method: Literal["hmmalign", "muscle"],
    threads: int = 1,
) -> MultipleSeqAlignment:
    """Helper function to perform alignment using a specified method.

    Args:
        instance (OrthologSeqs): The sequences to align.
        method (str): Alignment method, either 'hmmalign' or 'muscle'.
        threads (int, optional): Number of threads for 'muscle'. Defaults to 1.

    Returns:
        MultipleSeqAlignment: Resulting alignment.
    """
    if method == "hmmalign":
        r = instance.align(method=method, hmm=hmms[instance.name])
    else:
        r = instance.align(method=method, threads=threads)
    return r
