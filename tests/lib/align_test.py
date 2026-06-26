"""Tests for the align module library."""

from __future__ import annotations

import logging
import pickle
from multiprocessing.pool import ThreadPool
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyhmmer.plan7 import HMM

from phyling.exception import EmptyWarning, SeqtypeError
from phyling.lib.align import (
    HMMMarkerSet,
    OrthologList,
    OrthologSeqs,
    SampleList,
    SampleSeqs,
    SearchHit,
    SearchHitsManager,
    bp_mrtrans,
    run_hmmalign,
    run_hmmsearch,
    run_muscle,
)


@pytest.fixture(scope="module")
def mock_hmm() -> MagicMock:
    """Generates a base HMM mock object."""
    hmm = MagicMock()
    hmm.name = "HMM1"
    return hmm


@pytest.fixture(scope="module")
def mock_hmms_no_cutoff(mock_hmm) -> MagicMock:
    """HMMMarkerSet container mock holding mock_hmm with have_cutoffs = False."""
    hmms = MagicMock(spec=HMMMarkerSet)
    hmms.have_cutoffs.return_value = False

    container_data = [mock_hmm]
    hmms.__iter__.return_value = iter(container_data)
    hmms.__getitem__.side_effect = lambda index: container_data[index]
    hmms.__len__.return_value = len(container_data)

    return hmms


@pytest.fixture(scope="module")
def real_hmms(path_hmm_dir) -> HMMMarkerSet:
    return HMMMarkerSet(path_hmm_dir)


@pytest.fixture(scope="module")
def mock_hmms_with_cutoff(mock_hmm) -> MagicMock:
    """HMMMarkerSet container mock holding mock_hmm with have_cutoffs = True."""
    hmms = MagicMock(spec=HMMMarkerSet)
    hmms.have_cutoffs.return_value = True

    container_data = [mock_hmm]
    hmms.__iter__.return_value = iter(container_data)
    hmms.__getitem__.side_effect = lambda index: container_data[index]
    hmms.__len__.return_value = len(container_data)

    return hmms


@pytest.fixture
def mock_sampleseq_factory():
    """Returns a function that generates custom SampleSeqs mock instances."""

    def _make(name: str) -> MagicMock:
        mock = MagicMock(spec=SampleSeqs)
        mock.name = name
        mock.file = f"{name}.fasta"
        mock._data = ("data",)
        return mock

    return _make


@pytest.fixture
def mock_sampleseqs_a(mock_sampleseq_factory):
    return mock_sampleseq_factory("Sample_A")


@pytest.fixture
def mock_sampleseqs_b(mock_sampleseq_factory):
    return mock_sampleseq_factory("Sample_B")


@pytest.fixture
def mock_sampleseqs_c(mock_sampleseq_factory):
    return mock_sampleseq_factory("Sample_C")


@pytest.fixture
def mock_orthologseqs():
    name = "ortho_A"
    mock = MagicMock(spec=OrthologSeqs)
    mock.name = name
    mock.file = f"{name}.fasta"
    mock._data = ("data",)
    return mock


@pytest.fixture
def mock_alignment():
    alignment = MagicMock(spec=MultipleSeqAlignment)
    alignment.annotations = {}
    return alignment


@pytest.fixture
def mock_tempfile_factory(monkeypatch: pytest.MonkeyPatch, tmp_path):
    """A factory to dynamically create mock tempfiles supporting single or sequential context manager tracking."""
    mock_ctor = MagicMock()
    mock_files = []
    monkeypatch.setattr("phyling.lib.align.tempfile.NamedTemporaryFile", mock_ctor)

    def _make(file_name: str):
        mock_file = MagicMock()
        mock_file.name = str(tmp_path / file_name)

        # Setup __enter__ to return itself so it works inside 'with' blocks natively
        mock_file.__enter__.return_value = mock_file

        mock_files.append(mock_file)
        mock_ctor.side_effect = mock_files

        return {
            "constructor": mock_ctor,
            "instance": mock_file,
        }

    return _make


# ---------------------------------------------------------------------------
# TestHMMMarkerSet
# ---------------------------------------------------------------------------


@pytest.fixture
def path_hmm_file(path_hmm_dir) -> Path:
    return path_hmm_dir / "10at10240.hmm"


class TestHMMMarkerSet:
    def test_init_with_real_single_file(self, path_hmm_file: Path):
        """Verify loading a single real HMM file."""
        assert path_hmm_file.exists(), f"Missing test file: {path_hmm_file}"

        marker_set = HMMMarkerSet(data=path_hmm_file)

        # Assertions based on real file data
        assert len(marker_set) == 1
        assert marker_set[0].name == "10at10240"  # stem of the filename
        # Check that it's a real pyhmmer HMM object
        assert hasattr(marker_set[0], "consensus")

    def test_init_with_real_directory(self, path_hmm_dir: Path):
        """Verify loading all HMMs from the real directory."""
        assert path_hmm_dir.is_dir(), f"Missing test directory: {path_hmm_dir}"

        marker_set = HMMMarkerSet(data=path_hmm_dir)

        # Check that we loaded multiple files
        assert len(marker_set) > 0
        # Ensure they are all named correctly
        for hmm in marker_set:
            assert hmm.name is not None

    def test_init_with_cutoff_file(self, path_hmm_dir: Path, path_cutoff_file: Path):
        """Verify that bitscore cutoffs are applied correctly from the real scores file."""
        assert path_cutoff_file.exists(), f"Missing cutoff file: {path_cutoff_file}"

        # Initialize with both data and cutoffs
        marker_set = HMMMarkerSet(data=path_hmm_dir, cutoff_file=path_cutoff_file)

        # Check if the cutoff was actually applied to the HMM object
        # have_cutoffs returns True only if .cutoffs.trusted is set
        assert marker_set.have_cutoffs(verbose=True) is True

        # Specific bitscore check (if you know the value for 10at10240)
        # The code sets trusted as a tuple: (float, float)
        cutoff_val = marker_set["10at10240"].cutoffs.trusted
        assert isinstance(cutoff_val, tuple)
        assert len(cutoff_val) == 2

    def test_init_runtime_error_no_data(self, path_cutoff_file: Path):
        """Verify check for cutoff_file without data."""
        with pytest.raises(RuntimeError, match="Cannot specify cutoff_file without data"):
            HMMMarkerSet(data=None, cutoff_file=path_cutoff_file)  # type: ignore

    def test_getitem_by_name_real(self, path_hmm_file: Path):
        """Verify string-based retrieval using real data."""
        marker_set = HMMMarkerSet(data=path_hmm_file)

        # This should return the HMM object directly
        hmm = marker_set["10at10240"]
        assert hmm.name == "10at10240"

        with pytest.raises(KeyError):
            _ = marker_set["non_existent_marker"]

    def test_set_cutoffs_with_dict(self, path_hmm_file: Path):
        """Test the dictionary overload of set_cutoffs."""
        marker_set = HMMMarkerSet(data=path_hmm_file)
        manual_cutoffs = {"10at10240": 450.0}

        marker_set.set_cutoffs(manual_cutoffs)

        assert marker_set["10at10240"].cutoffs.trusted == (450.0, 450.0)

    def test_set_cutoffs_file_not_found(self, path_hmm_file: Path):
        """Verify FileNotFoundError logic."""
        marker_set = HMMMarkerSet(data=path_hmm_file)
        with pytest.raises(FileNotFoundError):
            marker_set.set_cutoffs("non_existent_file.tsv")

    def test_set_cutoffs_invalid_type(self, path_hmm_file: Path):
        """Verify TypeError for invalid input."""
        marker_set = HMMMarkerSet(data=path_hmm_file)
        with pytest.raises(TypeError, match="Invalid cutoffs format"):
            # Passing an integer instead of path or dict
            marker_set.set_cutoffs([1, 2, 3])  # type: ignore


# ---------------------------------------------------------------------------
# TestSampleSeqs
# ---------------------------------------------------------------------------


@pytest.fixture
def path_single_pep_fasta(path_pep_fasta) -> Path:
    return path_pep_fasta[0]


@pytest.fixture
def path_single_cds_fasta(path_cds_fasta) -> Path:
    return path_cds_fasta[0]


class TestSampleSeqs:
    def test_init_pep(self, path_single_pep_fasta: Path):
        """Test initialization and lazy loading of peptide sequences."""
        sample = SampleSeqs(path_single_pep_fasta, name="MPXV_PEP", seqtype="pep")
        assert sample.name == "MPXV_PEP"
        assert sample.file == path_single_pep_fasta.absolute()

        # Test __len__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            len(sample)

        # Test __iter__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            [x for x in sample]

        sample.load()
        assert len(sample) == len([x for x in sample])
        # Verify digital sequences are amino acids (Alphabet.amino codes are usually 0-20+)
        assert sample._data.alphabet.is_amino()  # type: ignore
        assert sample.seqtype == "pep"

    def test_init_cds(self, path_single_cds_fasta: Path):
        """Test that CDS sequences are correctly translated to protein on load."""
        sample = SampleSeqs(path_single_cds_fasta, name="MPXV_CDS", seqtype="dna")
        assert sample.name == "MPXV_CDS"
        assert sample.file == path_single_cds_fasta.absolute()

        # Test __len__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            len(sample)

        # Test __iter__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            [x for x in sample]

        sample.load()
        assert len(sample) == len([x for x in sample])
        # Even though input is DNA, _data should store translated Amino Acid sequences
        assert sample._data.alphabet.is_amino()  # type: ignore
        assert sample.seqtype == "dna"

    def test_init_auto(self, path_single_pep_fasta: Path, path_single_cds_fasta: Path):
        """Test guess seqtype from sequences."""
        sample = SampleSeqs(path_single_pep_fasta, name="MPXV_PEP")
        assert sample.seqtype == "pep"
        sample = SampleSeqs(path_single_cds_fasta, name="MPXV_CDS")
        assert sample.seqtype == "dna"

    def test_eq_lt(self, path_single_pep_fasta: Path, path_data_dir: Path):
        a = SampleSeqs(path_single_pep_fasta)
        b = SampleSeqs(path_single_pep_fasta)
        c = SampleSeqs(path_data_dir / "pep" / "Monkeypox_virus.faa.gz")
        assert a == b < c

    def test_eq_lt_seqtypeerror(self, path_single_pep_fasta: Path, path_single_cds_fasta: Path, path_data_dir: Path):
        a = SampleSeqs(path_single_pep_fasta)
        b = SampleSeqs(path_single_cds_fasta)
        c = SampleSeqs(path_data_dir / "pep" / "Monkeypox_virus.faa.gz")
        with pytest.raises(SeqtypeError, match="Items represent different seqtypes"):
            a == b
        with pytest.raises(SeqtypeError, match="Items represent different seqtypes"):
            assert b > c

    def test_search_mock_data(
        self, path_single_pep_fasta: Path, mock_hmms_no_cutoff: HMMMarkerSet, monkeypatch: pytest.MonkeyPatch
    ):
        captured_args = {}

        def mock_run_hmmsearch(sample, hmms, evalue, threads):
            # Record what was passed to verify it later
            captured_args["sample"] = sample
            captured_args["hmms"] = hmms
            captured_args["evalue"] = evalue
            captured_args["threads"] = threads
            # Return a dummy list to simulate SearchHits
            return ["hit1", "hit2"]

        monkeypatch.setattr("phyling.lib.align.run_hmmsearch", mock_run_hmmsearch)

        # 4. Initialize objects
        sample = SampleSeqs(path_single_pep_fasta, name="test_sample", seqtype="pep")

        # 5. Execute search (with custom parameters)
        results = sample.search(mock_hmms_no_cutoff, evalue=0.001, threads=4)

        # 6. Assertions
        assert results == ["hit1", "hit2"]
        assert captured_args["evalue"] == 0.001
        assert captured_args["threads"] == 4
        assert captured_args["sample"] == sample
        assert captured_args["hmms"] == mock_hmms_no_cutoff

    def test_search_real_data(self, path_single_pep_fasta: Path, real_hmms: HMMMarkerSet):
        """Run search using standard E-value filtering (no bitscore cutoffs)."""
        sample = SampleSeqs(path_single_pep_fasta, seqtype="pep")

        # run_hmmsearch should return a list of SearchHit
        hits = sample.search(real_hmms, evalue=1e-5)

        assert isinstance(hits, list)
        assert isinstance(hits[0], SearchHit)

    def test_cds_problematic_reporting(self, path_data_dir: Path, caplog: pytest.LogCaptureFixture):
        """Test that the class logs warnings for invalid CDS lengths."""
        # If you have a known 'bad' CDS file, use it here.
        # Otherwise, this just ensures the process finishes.
        sample = SampleSeqs(path_data_dir / "cds" / "Monkeypox_virus_with_bad_seq.fna", seqtype="dna")
        sample.load()

        # check if logger captured the warning.
        assert "seqs have invalid length" in caplog.text


# ---------------------------------------------------------------------------
# TestSampleList
# ---------------------------------------------------------------------------


class TestSampleList:
    def mock_helper(self, sample, hmms, evalue, threads):
        # Return a dummy list of "hits"
        return [f"hit_from_{sample.name}"]

    def test_sample_list_init_with_paths(self, path_pep_fasta: list[Path]):
        """Verify initialization with a list of file paths and test __getitem__."""
        sl = SampleList(path_pep_fasta)
        assert len(sl) == len(path_pep_fasta)

        subset = sl[1]
        assert isinstance(subset, SampleSeqs)
        subset = sl[0:2]
        assert isinstance(subset, SampleList)
        assert len(subset) == 2

        for i in range(len(sl)):
            assert isinstance(sl[i], SampleSeqs)
            assert sl[i].file == path_pep_fasta[i].absolute()
            assert sl[i].name == path_pep_fasta[i].name
            assert sl[i] == sl[path_pep_fasta[i].name]

        names = tuple(f"name_{i}" for i in range(len(path_pep_fasta)))
        sl = SampleList(path_pep_fasta, names=names)
        assert len(sl) == len(path_pep_fasta)
        for i in range(len(sl)):
            assert sl[i].file == path_pep_fasta[i].absolute()
            assert sl[i].name == names[i]
            assert sl[i] == sl[names[i]]

    def test_init_typeerror(self):
        with pytest.raises(TypeError, match="int cannot be converted to SampleSeqs"):
            SampleList([1, 2, 3])  # type: ignore

    def test_init_names_and_data_mismatched(self, path_pep_fasta: list[Path]):
        """Ensure error if names provided without data or with different length."""
        names = tuple(f"name_{i}" for i in range(len(path_pep_fasta)))
        with pytest.raises(RuntimeError, match="Received no data with names specified."):
            SampleList(data=[], names=names)

        with pytest.raises(RuntimeError, match="Data and names have different length."):
            SampleList(data=path_pep_fasta, names=names + ("extra_name",))

    def test_sample_list_init_mixed_seqtypes(self, path_pep_fasta: list[Path]) -> None:
        """Verify Error with files with mixed seqtypes."""
        mixed_seqtypes_files = path_pep_fasta + [Path("tests/data/cds/Monkeypox_virus.fna.gz")]
        with pytest.raises(SeqtypeError, match="Items represent different seqtypes"):
            SampleList(mixed_seqtypes_files)

    def test_search_sequential_logic(
        self, path_pep_fasta: list[Path], mock_hmms_with_cutoff: HMMMarkerSet, monkeypatch, caplog: pytest.LogCaptureFixture
    ):
        """Test sequential search mode logic."""
        sl = SampleList(path_pep_fasta)

        monkeypatch.setattr("phyling.lib.align._search_helper", self.mock_helper)

        with caplog.at_level("DEBUG"):
            results = sl.search(mock_hmms_with_cutoff, jobs=1, threads=2)

        assert "Sequential mode" in caplog.text
        assert f"Progress: {len(path_pep_fasta)} / {len(path_pep_fasta)}" in caplog.text
        assert f"hit_from_{path_pep_fasta[0].name}" in results

    def test_search_parallel_logic(
        self, path_pep_fasta: list[Path], mock_hmms_with_cutoff: HMMMarkerSet, monkeypatch, caplog: pytest.LogCaptureFixture
    ):
        """Test search orchestration using monkeypatch to avoid actual pyhmmer calls."""

        # 1. Setup SampleList
        sl = SampleList(path_pep_fasta)

        # Patch the helper used by 'partial' in the search method
        monkeypatch.setattr("phyling.lib.align._search_helper", self.mock_helper)

        # 3. Test Parallel Mode (jobs > 1)
        with caplog.at_level("DEBUG"):
            results = sl.search(mock_hmms_with_cutoff, jobs=3, threads=1)

        assert len(results) == len(path_pep_fasta)
        assert "Multiprocesses mode" in caplog.text
        assert f"Progress: {len(path_pep_fasta)} / {len(path_pep_fasta)}" in caplog.text
        assert f"hit_from_{path_pep_fasta[0].name}" in results

    def test_search_invalid_jobs(self, path_pep_fasta: list[Path], mock_hmms_with_cutoff: HMMMarkerSet):
        """Ensure RuntimeError when jobs are out of bounds."""
        sl = SampleList(path_pep_fasta)
        with pytest.raises(RuntimeError, match="jobs should be between 1 and"):
            sl.search(mock_hmms_with_cutoff, jobs=len(path_pep_fasta) + 1)


# ---------------------------------------------------------------------------
# TestSearchHitsManager
# ---------------------------------------------------------------------------


@pytest.fixture
def standard_hits(mock_sampleseqs_a, mock_sampleseqs_b) -> list[SearchHit]:
    return [
        SearchHit(hmm="HMM_1", sample=mock_sampleseqs_a, seqid="seq_01"),
        SearchHit(hmm="HMM_1", sample=mock_sampleseqs_b, seqid="seq_02"),
        SearchHit(hmm="HMM_2", sample=mock_sampleseqs_a, seqid="seq_03"),
    ]


@pytest.fixture(scope="class")
def real_search_hits(path_pep_fasta: list[Path], real_hmms) -> list[SearchHit]:
    sl = SampleList(path_pep_fasta[:3])
    results = sl.search(real_hmms, evalue=1e-10, jobs=1)
    return results


class TestSearchHitsManager:
    def test_init_and_update_encoding(self, standard_hits: list[SearchHit]):
        """Verify that initialization and batch updating correctly map string/object data to integers."""
        manager = SearchHitsManager(standard_hits)

        assert len(manager) == 3
        # Check that identical HMMs/Samples resolve to the same underlying integer ID
        assert manager._hmm_arr[0] == manager._hmm_arr[1]  # Both HMM_1
        assert manager._hmm_arr[0] != manager._hmm_arr[2]  # HMM_1 vs HMM_2

        assert manager._sample_arr[0] == manager._sample_arr[2]  # Both Sample_A
        assert manager._sample_arr[0] != manager._sample_arr[1]  # Sample_A vs Sample_B

        for i, hit in enumerate(standard_hits):
            assert manager._seqid_arr[i] == hit.seqid

    def test_getitem_single_element(
        self, standard_hits: list[SearchHit], mock_sampleseqs_a: SampleSeqs, mock_sampleseqs_b: SampleSeqs
    ):
        """Ensure element access correctly reverses integers back to namedtuple objects on-the-fly."""
        manager = SearchHitsManager(standard_hits)

        hit_0 = manager[0]
        assert isinstance(hit_0, SearchHit)
        assert hit_0.hmm == "HMM_1"
        assert hit_0.sample is mock_sampleseqs_a
        assert hit_0.seqid == "seq_01"

        hit_1 = manager[1]
        assert hit_1.hmm == "HMM_1"
        assert hit_1.sample is mock_sampleseqs_b
        assert hit_1.seqid == "seq_02"

    def test_getitem_lazy_slicing(self, standard_hits: list[SearchHit]):
        """Verify slicing yields a new manager sharing the master dictionary maps (lazy evaluation)."""
        manager = SearchHitsManager(standard_hits)

        sliced = manager[0:2]
        assert isinstance(sliced, SearchHitsManager)
        assert len(sliced) == 2

        # Lazy behavior check: Slicing should retain the parent dictionary mappings
        assert sliced._orthologs == manager._orthologs
        assert sliced._samples == manager._samples

    def test_pickle_serialization_loop(self, real_search_hits: list[SearchHit]):
        """Verify that __getstate__ and __setstate__ preserve all state and numpy arrays during serialization."""
        # 1. Setup original manager and populate data
        original_manager = SearchHitsManager(real_search_hits)

        # 2. Serialize (Pickle) the object to bytes (triggers __getstate__)
        pickled_data = pickle.dumps(original_manager)

        # 3. Deserialize (Unpickle) back to a new object (triggers __setstate__)
        unpickled_manager = pickle.loads(pickled_data)

        # 4. Assert structural integrity and object equivalence
        assert isinstance(unpickled_manager, SearchHitsManager)
        assert len(unpickled_manager) == len(original_manager)

        # Verify primary maps match exactly
        assert unpickled_manager._orthologs == original_manager._orthologs
        # Object references inside dictionaries should evaluate identically by key value
        assert list(unpickled_manager._samples.keys())[0].name == list(original_manager._samples.keys())[0].name

        # Verify NumPy array values and dtypes are preserved exactly
        assert np.array_equal(unpickled_manager._hmm_arr, original_manager._hmm_arr)
        assert np.array_equal(unpickled_manager._sample_arr, original_manager._sample_arr)
        assert np.array_equal(unpickled_manager._seqid_arr, original_manager._seqid_arr)

        assert unpickled_manager._hmm_arr.dtype == original_manager._hmm_arr.dtype
        assert unpickled_manager._seqid_arr.dtype == original_manager._seqid_arr.dtype

    def test_pickle_empty_manager(self):
        """Ensure an unpopulated or fresh manager can be serialized without throwing errors."""
        manager = SearchHitsManager()

        pickled_data = pickle.dumps(manager)
        unpickled_manager = pickle.loads(pickled_data)

        assert len(unpickled_manager._hmm_arr) == 0
        assert unpickled_manager._orthologs == {}

    def test_filter_drop_samples(
        self, standard_hits: list[SearchHit], mock_sampleseqs_a: SampleSeqs, mock_sampleseqs_b: SampleSeqs
    ):
        """Verify vectorized filtering correctly handles sample exclusions and drops metadata via compact."""
        manager = SearchHitsManager(standard_hits)

        filtered = manager.filter(drop_samples=["Sample_B"])

        assert len(filtered) == 2
        assert mock_sampleseqs_b not in filtered._samples
        rev_samp = {v: k for k, v in manager._samples.items()}
        assert all(rev_samp[s_idx] is mock_sampleseqs_a for s_idx in filtered._sample_arr)

    def test_filter_min_taxa_logic(
        self, mock_sampleseqs_a: SampleSeqs, mock_sampleseqs_b: SampleSeqs, mock_sampleseqs_c: SampleSeqs
    ):
        """Confirm the ortholog-aware taxa counter targets and purges markers below threshold limits."""
        hits = [
            # HMM_1 satisfies min_taxa=2 (found across 2 different samples)
            SearchHit(hmm="HMM_1", sample=mock_sampleseqs_a, seqid="s1"),
            SearchHit(hmm="HMM_1", sample=mock_sampleseqs_b, seqid="s2"),
            # HMM_2 fails min_taxa=2 (2 hits, but both trapped inside a single sample)
            SearchHit(hmm="HMM_2", sample=mock_sampleseqs_c, seqid="s3"),
            SearchHit(hmm="HMM_2", sample=mock_sampleseqs_c, seqid="s4"),
        ]
        manager = SearchHitsManager(hits)

        filtered = manager.filter(min_taxa=2)

        assert len(filtered) == 2
        assert "HMM_1" in filtered._orthologs
        assert "HMM_2" not in filtered._orthologs

    def test_filter_empty_raises_warning(self, standard_hits: list[SearchHit]):
        """Ensure an EmptyWarning gets triggered when filtering strips every row out."""
        manager = SearchHitsManager(standard_hits)

        with pytest.raises(EmptyWarning):
            manager.filter(drop_samples=["Sample_A", "Sample_B"])

    # @pytest.mark.dependency(depends=["search_integration"])
    def test_load_allocates_temporary_directory(self, real_search_hits: list[SearchHit]):
        """Verify that calling load provisions a physical temporary directory on disk."""
        manager = SearchHitsManager(real_search_hits)

        # Verify that prior to loading, no live directory tracking handles exist
        assert manager._mfa_dir is None

        # 1. Trigger the load sequence
        manager.load()

        # 2. Check that the temporary directory object exists
        assert manager._mfa_dir is not None

        # 3. Check that the physical path actually exists on the filesystem
        mfa_path = Path(manager._mfa_dir.name)
        assert mfa_path.exists()
        assert mfa_path.is_dir()

    # @pytest.mark.dependency(depends=["search_integration"])
    def test_unload_and_destruction_cleans_filesystem(self, real_search_hits: list[SearchHit]):
        """Verify that unload and object garbage collection completely purge disk files."""
        manager = SearchHitsManager(real_search_hits)

        manager.load()
        assert manager._mfa_dir is not None
        mfa_path = Path(manager._mfa_dir.name)
        assert mfa_path.exists()

        # 1. Trigger explicit unload
        manager.unload()

        # 2. Verify track state reset and physical path deletion
        assert manager._mfa_dir is None
        assert not mfa_path.exists()

    # @pytest.mark.dependency(depends=["search_integration"])
    def test_destructor_cleanup_fallback(self, real_search_hits: list[SearchHit]):
        """Ensure the __del__ magic method acts as a safety valve to prevent storage leaks."""
        manager = SearchHitsManager(real_search_hits)

        manager.load()
        assert manager._mfa_dir is not None
        mfa_path = Path(manager._mfa_dir.name)
        assert mfa_path.exists()

        # 1. Force object dereference to trigger __del__ garbage collection loop
        del manager

        # 2. Verify that the physical temporary directory was automatically scrubbed
        assert not mfa_path.exists()

    def test_compact_removes_ghosts_and_reindexes(
        self, standard_hits: list[SearchHit], mock_sampleseqs_a: SampleSeqs, mock_sampleseqs_b: SampleSeqs
    ):
        """Test that compact() completely drops unused categories and shifts remaining keys to 0..N."""
        manager = SearchHitsManager(standard_hits)

        # Slice down to just the hit containing HMM_2 and Sample_A
        sliced = manager[2:3]
        sliced.compact()

        # Verification of purged metadata
        assert "HMM_1" not in sliced._orthologs
        assert "HMM_2" in sliced._orthologs
        assert mock_sampleseqs_a in sliced._samples
        assert mock_sampleseqs_b not in sliced._samples
        assert len(sliced._samples) == 1

        # Verification of re-indexed arrays (must map to contiguous 0..N sequence)
        assert sliced._hmm_arr[0] == 0
        assert sliced._sample_arr[0] == 0

        # Test on-the-fly reconstruction after compaction
        reconstructed = sliced[0]
        assert reconstructed.hmm == "HMM_2"
        assert reconstructed.sample is mock_sampleseqs_a

    def test_compact_empty_manager(self):
        """Verify that calling compact on an empty manager behaves gracefully without errors."""
        manager = SearchHitsManager()
        manager.compact()

        assert manager._orthologs == {}
        assert manager._samples == {}
        assert len(manager._hmm_arr) == 0


# ---------------------------------------------------------------------------
# TestOrthologSeqs
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_alignment_runners(monkeypatch, mock_alignment):
    """Stub out the actual alignment processing functions."""
    monkeypatch.setattr("phyling.lib.align.run_hmmalign", lambda self, hmm: mock_alignment)
    monkeypatch.setattr("phyling.lib.align.run_muscle", lambda self, threads: mock_alignment)
    monkeypatch.setattr("phyling.lib.align.bp_mrtrans", lambda pep_msa, cds_data: mock_alignment)


@pytest.fixture
def path_single_pep_mfa(path_pep_mfa: list[Path]) -> Path:
    return path_pep_mfa[0]


@pytest.fixture
def path_single_cds_mfa(path_cds_mfa: list[Path]) -> Path:
    return path_cds_mfa[0]


class TestOrthologSeqs:
    def test_init_with_explicit_name(self, path_single_pep_mfa: Path):
        """Verify that initialization preserves an explicitly supplied name mapping."""
        ortho = OrthologSeqs(file=path_single_pep_mfa, name="custom_name", seqtype="pep")
        assert ortho.file == path_single_pep_mfa.absolute()
        assert ortho.name == "custom_name"
        assert ortho.seqtype == "pep"

    def test_init_with_auto_name(self, path_single_cds_mfa: Path):
        """Ensure that omitting a name falls back securely to the file name."""
        ortho = OrthologSeqs(file=path_single_cds_mfa, seqtype="dna")
        assert ortho.file == path_single_cds_mfa.absolute()
        assert ortho.name == path_single_cds_mfa.name
        assert ortho.seqtype == "dna"

    def test_load_initializes_empty_cds_attribute(self, path_single_cds_mfa: Path):
        """Confirm the load lifecycle properly sets up data slots before reading structures."""
        ortho = OrthologSeqs(path_single_cds_mfa)
        # Before load, slot shouldn't exist or isn't initialized
        with pytest.raises(AttributeError):
            _ = ortho._data_cds

        ortho.load()
        assert isinstance(ortho._data_cds, list)
        assert len(ortho._data_cds) == 5

    def test_process_pep_seqs(self, path_single_pep_mfa: Path):
        """Ensure raw peptide items pass sequentially straight into the primary array data stack."""
        ortho = OrthologSeqs(path_single_pep_mfa, seqtype="pep")
        ortho.load()

        assert len(ortho._data) == 5
        assert "Goatpox_virus.faa.gz" in [x.name for x in ortho._data]

    def test_process_cds_seqs_successful(self, path_single_cds_mfa: Path):
        """Verify clean CDS text reading converts data to translated peptides and stores back-mapped CDS templates."""
        ortho = OrthologSeqs(path_single_cds_mfa, seqtype="dna")
        ortho.load()

        assert len(ortho._data) == 5
        assert "Goatpox_virus.fna.gz" in [x.name for x in ortho._data]
        assert len(ortho._data_cds) == 5
        assert isinstance(ortho._data_cds[0], SeqRecord)
        assert "Goatpox_virus.fna.gz" in [x.id for x in ortho._data_cds]

    def test_process_cds_seqs_with_invalid_lengths(self, caplog: pytest.LogCaptureFixture):
        """Test that translation anomalies trigger validation exceptions and get gracefully logged as warnings."""
        ortho = OrthologSeqs("tests/data/cds/Monkeypox_virus_with_bad_seq.fna", seqtype="dna")
        with caplog.at_level("WARNING"):
            ortho.load()

        assert "have invalid length" in caplog.text

    def test_search_raises_not_implemented_error(self, path_single_pep_mfa: Path, mock_hmms_with_cutoff: HMMMarkerSet):
        """Enforce the constraint that search operations must explicitly fail on this subclass."""
        ortho = OrthologSeqs(path_single_pep_mfa)
        with pytest.raises(NotImplementedError) as exc_info:
            ortho.search(mock_hmms_with_cutoff)
        assert "Method is not implemented" in str(exc_info.value)

    def test_align_hmmalign_missing_argument(self, path_single_pep_mfa: Path):
        """Ensure executing an HMM-driven alignment architecture without a reference profile fails safely."""
        ortho = OrthologSeqs(path_single_pep_mfa)
        with pytest.raises(ValueError) as exc_info:
            ortho.align(method="hmmalign", hmm=None)  # type: ignore
        assert 'required "hmm" argument' in str(exc_info.value)

    def test_align_invalid_method_argument(self, path_single_pep_mfa: Path):
        """Enforce string-literal verification constraints over unsupported foreign alignment backends."""
        ortho = OrthologSeqs(path_single_pep_mfa)
        with pytest.raises(ValueError) as exc_info:
            ortho.align(method="unsupported_engine")  # type: ignore
        assert "Argument method only accepts" in str(exc_info.value)

    @pytest.mark.parametrize("method", ["hmmalign", "muscle"])
    def test_align_execution_flows_pep(
        self, path_single_pep_mfa: Path, method: str, mock_hmms_with_cutoff: HMMMarkerSet, mock_alignment_runners
    ):
        """Validate cross-combinations of text-types and backend align engines compute successfully."""
        ortho = OrthologSeqs(path_single_pep_mfa, seqtype="pep")
        ortho._data_cds = []

        res = ortho.align(method=method, hmm=mock_hmms_with_cutoff, threads=2)  # type: ignore

        assert isinstance(res, MultipleSeqAlignment)

    @pytest.mark.parametrize("method", ["hmmalign", "muscle"])
    def test_align_execution_flows_cds(
        self, path_single_cds_mfa: Path, method: str, mock_hmms_with_cutoff: HMMMarkerSet, mock_alignment_runners
    ):
        """Validate cross-combinations of text-types and backend align engines compute successfully."""
        ortho = OrthologSeqs(path_single_cds_mfa, seqtype="dna")
        ortho._data_cds = []

        res = ortho.align(method=method, hmm=mock_hmms_with_cutoff, threads=2)  # type: ignore

        assert isinstance(res, MultipleSeqAlignment)


# ---------------------------------------------------------------------------
# TestOrthologList
# ---------------------------------------------------------------------------


@pytest.fixture
def spy_align_env(monkeypatch: pytest.MonkeyPatch, mock_alignment):
    monkeypatch.setattr("phyling.lib.align._align_helper", lambda ortholog_seq, hmms, method, threads: mock_alignment)


class TestOrthologList:
    def test_init_bound_class_assignment(self, path_pep_mfa: list[Path]):
        """Verify that initialization sets the correct internal target class type."""
        # Partially patch init to let hasattr check fire normally
        ortho_list = OrthologList(data=path_pep_mfa, seqtype="pep")
        assert ortho_list._bound_class is OrthologSeqs

    def test_getitem_passthrough(self, path_pep_mfa: list[Path]):
        """Ensure that indexing or fetching keys safely proxies up to the base ABC container."""
        ortho_list = OrthologList(path_pep_mfa)

        # Accessing via index or text label should yield an OrthologSeqs interface instance
        result = ortho_list[0]
        assert isinstance(result, OrthologSeqs)

    @pytest.mark.parametrize(
        "jobs",
        [0, 10],
    )
    def test_align_jobs_out_of_bounds_raises_error(self, path_pep_mfa: list[Path], jobs):
        """Enforce range validation bounds on process task distribution counters."""
        ortho_list = OrthologList(path_pep_mfa)

        # Length is mocked to 2. Setting jobs to 0 or 3 should trigger a validation crash
        with pytest.raises(RuntimeError) as exc_info:
            ortho_list.align(method="muscle", jobs=jobs)
        assert "jobs should be between 1 and" in str(exc_info.value)

    def test_align_sequential_mode_flow(
        self, path_pep_mfa: list[Path], mock_hmms_with_cutoff: HMMMarkerSet, spy_align_env, caplog: pytest.LogCaptureFixture
    ):
        """Verify that processing runs sequentially when jobs=1 and updates logs."""
        ortho_list = OrthologList(path_pep_mfa)

        with caplog.at_level(logging.DEBUG):
            alignments = ortho_list.align(method="hmmalign", hmms=mock_hmms_with_cutoff, jobs=1, threads=2)

        assert len(alignments) == 4
        assert isinstance(alignments[0], MultipleSeqAlignment)
        assert "Sequential mode with 2 threads." in caplog.text
        assert "Progress: 4 / 4" in caplog.text

    def test_align_parallel_threadpool_mode_flow(
        self, path_pep_mfa: list[Path], spy_align_env, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ):
        """Confirm multiprocess threading pool contexts instantiate and process cleanly."""
        ortho_list = OrthologList(path_pep_mfa)

        monkeypatch.setattr(ThreadPool, "imap", lambda self, func, iterable: map(func, iterable))

        with caplog.at_level(logging.DEBUG):
            alignments = ortho_list.align(method="muscle", jobs=2, threads=4)

        assert len(alignments) == 4
        assert "Multiprocesses mode with 2 jobs and 4 threads" in caplog.text
        assert "Progress: 4 / 4" in caplog.text

    def test_align_exception_handling_logs_and_re_raises(
        self, path_pep_mfa: list[Path], monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
    ):
        """Verify that pipeline execution failures register logs correctly before throwing exceptions."""
        ortho_list = OrthologList(path_pep_mfa)

        def crashing_helper(*args, **kwargs):
            raise ValueError("MUSCLE binary execution failure")

        monkeypatch.setattr("phyling.lib.align._align_helper", crashing_helper)

        with pytest.raises(ValueError, match="MUSCLE binary execution failure"):
            with caplog.at_level(logging.ERROR):
                ortho_list.align(method="muscle", jobs=1)

        # Confirm traceback structural details captured inside error logger
        assert "Traceback" in caplog.text


# ---------------------------------------------------------------------------
# TestRunHmmsearch
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_tophits_factory():
    """A factory fixture to generate list-capable TopHits mocks dynamically."""

    def _create_tophits(query_name: str, hits_list: list) -> MagicMock:
        mock_tophits = MagicMock(wraps=hits_list)

        # 2. Add the sub-mock query property
        mock_tophits.query = MagicMock()
        mock_tophits.query.name = query_name

        # 3. Fix Python's handling of the list length under a 'wraps' mock
        mock_tophits.__len__.return_value = len(hits_list)
        mock_tophits.__getitem__.side_effect = lambda idx: hits_list[idx]

        return mock_tophits

    return _create_tophits


@pytest.fixture
def spy_hmmsearch(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Fixtures that patches hmmsearch with a mock and returns it for controller tracking."""
    mock_search_func = MagicMock()
    monkeypatch.setattr("phyling.lib.align.hmmsearch", mock_search_func)
    return mock_search_func


class TestRunHmmsearch:
    def test_run_hmmsearch_uses_evalue_when_no_cutoffs(
        self, mock_sampleseqs_a: SampleSeqs, mock_hmms_no_cutoff: HMMMarkerSet, mock_tophits_factory, spy_hmmsearch
    ):
        """Verify that E-value thresholds are applied when trusted cutoffs are missing."""
        mock_hit_1 = MagicMock(name="seq_1", reported=True)
        mock_hit_2 = MagicMock(name="seq_2", reported=False)
        fake_tophits = mock_tophits_factory(query_name="HMM1", hits_list=[mock_hit_1, mock_hit_2])
        spy_hmmsearch.return_value = [fake_tophits]

        results = run_hmmsearch(mock_sampleseqs_a, mock_hmms_no_cutoff, evalue=1e-5, threads=2)

        spy_hmmsearch.assert_called_once_with(mock_hmms_no_cutoff, mock_sampleseqs_a, cpus=2, E=1e-5)
        assert len(results) == 1

    def test_run_hmmsearch_uses_trusted_cutoffs(
        self, mock_sampleseqs_a: SampleSeqs, mock_hmms_with_cutoff: HMMMarkerSet, mock_tophits_factory, spy_hmmsearch
    ):
        """Ensure bit_cutoffs flag shifts to 'trusted' when reference cutoffs are present."""
        mock_hit_1 = MagicMock(name="seq_1", reported=False)
        mock_hit_2 = MagicMock(name="seq_2", reported=True)
        fake_tophits = mock_tophits_factory(query_name="HMM1", hits_list=[mock_hit_1, mock_hit_2])
        spy_hmmsearch.return_value = [fake_tophits]

        results = run_hmmsearch(mock_sampleseqs_a, mock_hmms_with_cutoff, evalue=1e-5, threads=2)

        spy_hmmsearch.assert_called_once_with(mock_hmms_with_cutoff, mock_sampleseqs_a, cpus=2, bit_cutoffs="trusted")
        assert len(results) == 1

    def test_run_hmmsearch_skips_unreported_hits(
        self, mock_sampleseqs_a: SampleSeqs, mock_hmms_no_cutoff: HMMMarkerSet, mock_tophits_factory, spy_hmmsearch
    ):
        """Ensure that hits that fail the inclusion threshold ('reported=False') are dropped."""
        mock_hit_1 = MagicMock(name="seq_1", reported=False)
        mock_hit_2 = MagicMock(name="seq_2", reported=False)
        fake_tophits = mock_tophits_factory(query_name="HMM1", hits_list=[mock_hit_1, mock_hit_2])
        spy_hmmsearch.return_value = [fake_tophits]

        results = run_hmmsearch(mock_sampleseqs_a, mock_hmms_no_cutoff)
        assert len(results) == 0

    def test_run_hmmsearch_skips_multi_reported_hits(
        self, mock_sampleseqs_a: SampleSeqs, mock_hmms_no_cutoff: HMMMarkerSet, mock_tophits_factory, spy_hmmsearch
    ):
        """Confirm that multiple reported hits for a single profile discard the entire marker sequence block."""
        # Mock a marker matching two separate loci on the genome (paralogs/non-specific hits)
        mock_hit_1 = MagicMock(name="seq_1", reported=True)
        mock_hit_2 = MagicMock(name="seq_2", reported=True)
        fake_tophits = mock_tophits_factory(query_name="HMM1", hits_list=[mock_hit_1, mock_hit_2])
        spy_hmmsearch.return_value = [fake_tophits]

        results = run_hmmsearch(mock_sampleseqs_a, mock_hmms_no_cutoff)

        # Due to 'if len(reported) > 1: continue', this entire marker group should be dropped
        assert len(results) == 0

    def test_run_hmmsearch_logs_completion(
        self,
        mock_sampleseqs_a: SampleSeqs,
        mock_hmms_no_cutoff: HMMMarkerSet,
        mock_tophits_factory,
        spy_hmmsearch,
        caplog: pytest.LogCaptureFixture,
    ):
        """Verify that execution completions trigger expected diagnostic logging messages."""
        mock_hit_1 = MagicMock(name="seq_1", reported=True)
        mock_hit_2 = MagicMock(name="seq_2", reported=False)
        fake_tophits = mock_tophits_factory(query_name="HMM1", hits_list=[mock_hit_1, mock_hit_2])
        spy_hmmsearch.return_value = [fake_tophits]

        with caplog.at_level(logging.DEBUG):
            _ = run_hmmsearch(mock_sampleseqs_a, mock_hmms_no_cutoff)

        assert f"Hmmsearch on {mock_sampleseqs_a.name} is done." in caplog.text


# ---------------------------------------------------------------------------
# TestRunHmmalign
# ---------------------------------------------------------------------------


@pytest.fixture
def spy_hmmalign_env(monkeypatch: pytest.MonkeyPatch, mock_alignment, mock_tempfile_factory):
    """A lightweight patch for hmmalign and load_msa."""
    mock_hmmalign = MagicMock()

    # 2. Apply patches to the target paths
    monkeypatch.setattr("phyling.lib.align.hmmalign", mock_hmmalign)
    mock_load_msa = MagicMock()
    mock_load_msa.return_value = mock_alignment
    monkeypatch.setattr("phyling.lib.align.load_msa", mock_load_msa)
    mock_tempfile = mock_tempfile_factory("hmmalign")

    return {
        "runner": mock_hmmalign,
        "alignment_result": mock_alignment,
        "tempfile_constructor": mock_tempfile["constructor"],
        "tempfile_instance": mock_tempfile["instance"],
    }


class TestRunHmmalign:
    def test_run_hmmalign_execution_flow(
        self, spy_hmmalign_env, mock_orthologseqs: OrthologSeqs, mock_hmm: HMM, caplog: pytest.LogCaptureFixture
    ):
        """Verify the full execution matrix: invocation, temporary file consumption, and tagging."""
        with caplog.at_level(logging.DEBUG):
            result = run_hmmalign(mock_orthologseqs, mock_hmm)

        spy_hmmalign_env["runner"].assert_called_once()
        spy_hmmalign_env["tempfile_constructor"].assert_called_once()
        spy_hmmalign_env["tempfile_instance"].seek.assert_called_once_with(0)
        spy_hmmalign_env["tempfile_instance"].close.assert_called_once()

        assert isinstance(result, MultipleSeqAlignment)
        assert result.annotations == {"seqtype": "pep"}
        assert f"Hmmalign on {mock_orthologseqs.name} is done." in caplog.text

    def test_run_hmmalign_cleans_up_file_on_exception(
        self, spy_hmmalign_env, monkeypatch: pytest.MonkeyPatch, mock_orthologseqs: OrthologSeqs, mock_hmm: HMM
    ):
        """Guarantee that NamedTemporaryFile descriptors close and drop even if the pipeline crashes."""
        spy_hmmalign_env["runner"].side_effect = RuntimeError("Execution failed")
        with pytest.raises(RuntimeError):
            run_hmmalign(mock_orthologseqs, mock_hmm)

        spy_hmmalign_env["tempfile_constructor"].assert_called_once()
        spy_hmmalign_env["tempfile_instance"].seek.assert_not_called()
        spy_hmmalign_env["tempfile_instance"].close.assert_called_once()


# ---------------------------------------------------------------------------
# TestRunMuscle
# ---------------------------------------------------------------------------


@pytest.fixture
def spy_muscle_env(monkeypatch: pytest.MonkeyPatch, mock_alignment, mock_tempfile_factory):
    # 3. Mock Muscle runner and load_msa
    mock_muscle_runner = MagicMock()
    mock_muscle_class = MagicMock(return_value=mock_muscle_runner)

    monkeypatch.setattr("phyling.lib.align.Muscle", mock_muscle_class)
    mock_load_msa = MagicMock()
    mock_load_msa.return_value = mock_alignment
    monkeypatch.setattr("phyling.lib.align.load_msa", mock_load_msa)
    mock_f_aln = mock_tempfile_factory("muscle")
    mock_f_pep = mock_tempfile_factory("pep_fasta")

    return {
        "runner": mock_muscle_runner,
        "alignment_result": mock_alignment,
        "tempfile_constructor": mock_f_pep["constructor"],
        "f_aln_instance": mock_f_aln["instance"],
        "f_pep_instance": mock_f_pep["instance"],
    }


class TestRunMuscle:
    def test_run_muscle_execution_flow(
        self, spy_muscle_env, mock_orthologseqs: OrthologSeqs, mock_hmm: HMM, caplog: pytest.LogCaptureFixture
    ):
        """Verify standard successful execution flow, parameter mapping, and typing assertions."""
        with caplog.at_level(logging.DEBUG):
            result = run_muscle(mock_orthologseqs, threads=4)

        spy_muscle_env["runner"].run.assert_called_once
        assert spy_muscle_env["tempfile_constructor"].call_count == 2
        spy_muscle_env["f_aln_instance"].seek.assert_called_once_with(0)
        spy_muscle_env["f_pep_instance"].seek.assert_called_once_with(0)
        spy_muscle_env["f_aln_instance"].__exit__.assert_called_once()
        spy_muscle_env["f_pep_instance"].__exit__.assert_called_once()

        assert isinstance(result, MultipleSeqAlignment)
        assert result.annotations == {"seqtype": "pep"}
        assert f"Muscle on {mock_orthologseqs.name} is done." in caplog.text

    def test_run_muscle_exception_handling_and_re_raise(
        self,
        spy_muscle_env,
        monkeypatch: pytest.MonkeyPatch,
        mock_orthologseqs: OrthologSeqs,
        mock_hmm: HMM,
        caplog: pytest.LogCaptureFixture,
    ):
        """Confirm errors encountered at runtime wrap the crash trace in a localized error explanation."""
        spy_muscle_env["runner"].run.side_effect = RuntimeError("Execution failed")

        # Enforce that the runner bubbles up into an informative RuntimeError container wrap
        with pytest.raises(RuntimeError) as exc_info:
            run_muscle(mock_orthologseqs, threads=2)

        spy_muscle_env["runner"].run.assert_called_once()
        assert spy_muscle_env["tempfile_constructor"].call_count == 2
        spy_muscle_env["f_aln_instance"].seek.assert_not_called()
        spy_muscle_env["f_pep_instance"].seek.assert_called_once_with(0)
        spy_muscle_env["f_aln_instance"].__exit__.assert_called_once()
        spy_muscle_env["f_pep_instance"].__exit__.assert_called_once()

        assert f"Muscle failed on peptide fasta translated from {mock_orthologseqs.file}" in str(exc_info.value)


class TestBpMrtrans:
    def test_bp_mrtrans_successful_translation(self):
        """Verify back-translation maps valid codons 1:1 against non-gapped alignment templates."""
        # Protein: M   A   Y   A
        pep_rec = SeqRecord(Seq("MAYA"), id="gene_01", name="n1", description="desc1")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # DNA matching exactly 4 codons: ATG (M), GCT (A), TAC (Y), GCA (A)
        cds_rec = SeqRecord(Seq("ATGGCCTACGCA"), id="gene_01", name="n1", description="desc1")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)

        assert isinstance(result, MultipleSeqAlignment)
        assert result.annotations == {"seqtype": "dna"}
        assert len(result) == 1
        assert str(result[0].seq) == "ATGGCCTACGCA"
        assert result[0].id == "gene_01"

    def test_bp_mrtrans_with_gaps_in_protein(self):
        """Ensure protein gaps ('-') project 3 nucleotide gaps into the cDNA sequence."""
        # Protein alignment with a deletion insertion gap: M  -  Y
        pep_rec = SeqRecord(Seq("M-Y"), id="gene_02")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # The underlying coding sequence contains only the real structural coordinates: ATG (M), TAC (Y)
        cds_rec = SeqRecord(Seq("ATGTAC"), id="gene_02")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)

        # Expected out: ATG (from M), --- (from gap), TAC (from Y)
        assert str(result[0].seq) == "ATG---TAC"

    def test_bp_mrtrans_handles_trailing_stop_codons(self):
        """Verify that stop codons encountered concurrently with alignment boundaries parse correctly."""
        # Protein sequence layout: M  A  Y
        pep_rec = SeqRecord(Seq("MAY"), id="gene_03")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # DNA sequence ending in a valid stop codon (TAA)
        cds_rec = SeqRecord(Seq("ATGGCCTACTAA"), id="gene_03")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)

        # The function handles a trailing stop codon matching a boundary by expanding it out if present
        assert str(result[0].seq) == "ATGGCCTAC"

    def test_bp_mrtrans_metadata_fallback_logic(self):
        """Confirm that if SeqRecord properties are missing from the CDS, they fall back to the peptide's."""
        pep_rec = SeqRecord(Seq("M"), id="pep_id", name="pep_name", description="pep_desc")
        pep_msa = MultipleSeqAlignment([pep_rec])

        cds_rec = SeqRecord(Seq("ATG"))
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)  # type: ignore

        # Verify fallback properties mapped successfully
        assert result[0].id == "pep_id"
        assert result[0].name == "pep_name"
        assert result[0].description == "pep_desc"

    def test_bp_mrtrans_case_insensitivity_normalization(self):
        """Ensure lower-case coding letters normalize seamlessly into upper-case alignments."""
        pep_rec = SeqRecord(Seq("M"), id="gene_04")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # Lower-case string sequence representation
        cds_rec = SeqRecord(Seq("atg"), id="gene_04")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)
        assert str(result[0].seq) == "ATG"
