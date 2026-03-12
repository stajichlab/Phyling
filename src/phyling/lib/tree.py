"""Filter and tree modules library."""

from __future__ import annotations

import gzip
import logging
import tempfile
import traceback
from functools import wraps
from itertools import product
from multiprocessing.pool import ThreadPool
from pathlib import Path
from typing import Callable, Literal, Sequence, TypeVar, overload

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord

from .. import AVAIL_CPUS
from . import SeqTypes, TreeMethods, TreeOutputFiles, _abc
from ._utils import CheckAttrs, Timer, guess_seqtype, is_gzip_file

__all__ = ["MFA2Tree", "MFA2TreeList"]
logger = logging.getLogger(__name__)
_R = TypeVar("Return")
_C = TypeVar("Callable", bound=Callable[..., _R])


def _check_attributes(*attrs: str) -> Callable[[_C], _C]:
    """Decorator to ensure specific attributes are initialized before executing the function.

    Args:
        *attrs: Attribute names to check in the instance.

    Raises:
        AttributeError: If any specified attribute is `None` in the instance.
    """
    var_mapping = {
        "_tree": "build",
        "_toverr": "compute_toverr",
        "_saturation": "compute_saturation",
    }
    invalid_attrs = [attr for attr in attrs if attr not in var_mapping]
    if invalid_attrs:
        raise AttributeError(f"Invalid attribute names: {invalid_attrs}")

    def decorator(func: _C) -> _C:
        @wraps(func)
        def wrapper(instance, *args, **kwargs):
            """Validate variable inequality and execute the wrapped function."""
            none_attrs = CheckAttrs.is_none(instance, *attrs)
            for var in sorted(none_attrs, key=lambda x: list(var_mapping.keys()).index(x)):
                raise AttributeError(f"Please run the {var_mapping[var]} method first.")
            return func(instance, *args, **kwargs)

        return wrapper

    return decorator


def _check_single_file(func: _C) -> _C:
    """Decorator to ensure that the method is only called on a list with multiple MFA files.

    This decorator checks if the `MFA2TreeList` instance contains more than one MFA file. If there is only one file, it raises a
    `RuntimeError`. Otherwise, it calls the wrapped function.

    Args:
        func (Callable): The function to be wrapped.

    Returns:
        Callable: The decorated function.

    Raises:
        RuntimeError: If the `MFA2TreeList` instance contains fewer than two MFA files.
    """

    @wraps(func)
    def wrapper(self: MFA2TreeList, *args, **kwargs):
        if len(self) < 2:
            raise RuntimeError(f"{func.__name__} cannot be run since only found a single mfa.")
        return func(self, *args, **kwargs)

    return wrapper


class MFA2Tree(_abc.SeqFileWrapperABC):
    """Convert peptide/DNA multiple sequence alignment to Biopython Tree object.

    Optionally take a partition file to build tree with partition mode with RAxML-NG and IQ-TREE.
    """

    __slots__ = ("_logger", "_method", "_tree", "_toverr", "_saturation")

    @overload
    def __init__(self, file: str | Path, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...

    @overload
    def __init__(self, file: str | Path, name: str, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...

    def __init__(self, file: str | Path, name: str | None = None, *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None:
        """Initialize a MFA2Tree object.

        Initialize with a path of a multiple sequence alignment file (in fasta format) and a representative name (optional). The
        basename of the file will be used as the name if no name specified. The fasta won't be loaded immediately unless calling
        the load method.

        Args:
            file (str | Path): The path to the MSA file (plain or bgzipped).
            name (str, optional): A representative name for the sequences. Defaults to None.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Examples:
            Load a peptide fasta:
            >>> MFA2Tree("data/101133at4751.aa.mfa")

            Load a peptide fasta with a given name:
            >>> MFA2Tree("data/101133at4751.aa.mfa", "msa_101133at4751")

            Load a bgzf compressed cds fasta:
            >>> MFA2Tree("data/101133at4751.cds.mfa", "msa_101133at4751")
        """
        super().__init__(file, name, seqtype=seqtype)
        self._logger: logging.Logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._method: str = None
        self._tree: Tree = None
        self._toverr: float = None
        self._saturation: float = None

    @_abc.check_loaded
    def __len__(self) -> int:
        """Return the length of the MSA.

        Returns:
            int: The length of the MSA.
        """
        return self._data.get_alignment_length()

    @_abc.check_class
    def __gt__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is larger than another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is larger, otherwise False.
        """
        return self.toverr > other.toverr

    @_abc.check_class
    def __ge__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is larger than or equal to another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is larger or equal, otherwise False.
        """
        return self.toverr >= other.toverr

    @_abc.check_class
    def __lt__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is smaller than another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is smaller, otherwise False.
        """
        return self.toverr < other.toverr

    @_abc.check_class
    def __le__(self, other: MFA2Tree) -> bool:
        """Compare if the current object is smaller than or equal to another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the current object's toverr is smaller or equal, otherwise False.
        """
        return self.toverr <= other.toverr

    @_abc.check_class
    def __eq__(self, other: MFA2Tree) -> bool:
        """Check if the current object is equal to another object.

        Args:
            other (MFA2Tree): Another instance of the same class.

        Returns:
            bool: True if the objects have the same toverr, otherwise False.
        """
        return self.toverr == other.toverr

    def load(self) -> None:
        """Load the MultipleSeqAlignment.

        Assign the sequence type and its name after the MSA being loaded.
        """
        self._data = AlignIO.read(self.file, format="fasta")
        self._data.annotations["seqtype"] = self.seqtype
        self._data.annotations["seqname"] = self.name

    @property
    @_check_attributes("_tree")
    def tree(self) -> Tree:
        """Return the Biopython tree object.

        Raises:
            AttributeError: If tree haven't built.
        """
        return self._tree

    @property
    @_check_attributes("_tree")
    def method(self) -> str:
        """Return the method used for tree building.

        Raises:
            AttributeError: If tree haven't built.
        """
        return self._method

    @property
    @_check_attributes("_tree", "_toverr")
    def toverr(self) -> float:
        """Return the treeness/RCVs.

        Raises:
            AttributeError: If tree haven't built.
            AttributeError: If the toverr haven't calculated.
        """
        return self._toverr

    @property
    @_check_attributes("_tree", "_saturation")
    def saturation(self) -> float:
        """Return the saturation.

        Raises:
            AttributeError: If tree haven't built.
            AttributeError: If the saturation haven't calculated.
        """
        return self._saturation

    @overload
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        output: None,
        model: str | Path = "AUTO",
        *,
        noml: bool = False,
        bs: int = 0,
        scfl: int = 0,
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
    ) -> Tree: ...

    @overload
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        output: str | Path,
        model: str | Path = "AUTO",
        *,
        noml: bool = False,
        bs: int = 0,
        scfl: int = 0,
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
    ) -> Path: ...

    @_abc.load_data
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        output: str | Path | None = None,
        model: str | Path = "AUTO",
        *,
        noml: bool = False,
        bs: int = 0,
        scfl: int = 0,
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
        **kwargs,
    ) -> Tree | Path:
        """Build a phylogenetic tree using the specified method.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
            output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
            model (str | Path | None, optional): String of model or path to a partition file for model partitioning. Defaults to
                "AUTO" that use ModelFinder to identify the best model.
            noml (bool): Disable Maximum-likelihood estimation to speed up tree building. Only available for FastTree. Defaults to
                False.
            bs (int): UFBoot value. Must be at least 1000 if set. Defaults to 0 (disabled).
            scfl (int): Number of quartets to randomly sample around each internal branch when computing the site concordance
                factor (sCF). A minimum of 100 is recommended for stable results. Defaults to 0 (disabled).
            seed (int): Seed number for reproducibility. Defaults to -1 (random seed).
            threads (int): Number of threads to use. Only applicable to RAxML-NG and IQ-TREE. Defaults to -1
                (autodetect).
            max_threads (int): Maximum number of threads used during autodetect. Only applicable when `threads` is set to -1.
                Defaults to maximum available CPUs.

        Returns:
            Tree | Path:
                - If `output` is None (default), returns a `Bio.Phylo.BaseTree.Tree` object.
                - If `output` is str or Path, returns a `pathlib.Path` to the output tree file.
        """
        from ..external import (
            FastTree,
            Iqtree,
            ModelFinder,
            Raxml,
        )

        method = method.upper()
        if method not in [m.name for m in TreeMethods]:
            raise KeyError(f"Invalid method: {method}. Expected one of: {', '.join([m.name.lower() for m in TreeMethods])}")

        tempdir = None
        try:
            if not output:
                tempdir = tempfile.TemporaryDirectory()
                output = tempdir.name

            output = Path(output)
            output.mkdir(exist_ok=True)

            if model == "AUTO":
                modelfinder = ModelFinder(
                    self.file,
                    output / "modelfinder",
                    seqtype=self.seqtype,
                    method=method.lower(),
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                )
                self._logger.info("Find the best-fit model...")
                modelfinder.run()
                model = modelfinder.result
                self._logger.info(f"Best-fit model: {model}")
            elif Path(model).is_file():
                # Partitioning analysis
                if method == TreeMethods.FT.name:
                    raise ValueError(f"{TreeMethods.FT.method} does not support partitioning analysis.")

            if method == TreeMethods.FT.name:
                runner = FastTree(
                    self.file,
                    output / f"{self.file.name}.nw",
                    seqtype=self.seqtype,
                    model=model,
                    noml=noml,
                    **kwargs,
                )
            elif method == TreeMethods.RAXML.name:
                runner = Raxml(
                    self.file,
                    output / self.file.name,
                    seqtype=self.seqtype,
                    model=model,
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                    **kwargs,
                )
            elif method == TreeMethods.IQTREE.name:
                runner = Iqtree(
                    self.file,
                    output / self.file.name,
                    seqtype=self.seqtype,
                    model=model,
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                    **kwargs,
                )
            else:
                raise NotImplementedError(f"Method not implemented yet: {TreeMethods[method].name}.")

            self._logger.info(f"Build tree by {runner._prog}...")
            runner.run()
            self._method = method
            self._tree = Phylo.read(runner.result, "newick")

            if bs > 0:
                self._logger.info("Bootstrapping...")
                tree_file = bootstrap(
                    self,
                    runner.result,
                    output / "ufboot",
                    runner.model,
                    bs=bs,
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                )
                self._tree = Phylo.read(tree_file, "newick")

            if scfl > 0:
                self._logger.info("Calculating site concordance factor...")
                tree_file = branch_concordance(
                    self,
                    tree_file,
                    output / "concord",
                    runner.model,
                    scfl=scfl,
                    seed=seed,
                    threads=threads,
                    threads_max=threads_max,
                )
                self._tree = Phylo.read(tree_file, "newick")

        finally:
            if tempdir:
                tempdir.cleanup()
        self._logger.info("Tree building on %s is done.", self.name)

        if tempdir:
            return self.tree
        else:
            return runner.result

    @_check_attributes("_tree")
    def compute_toverr(self) -> None:
        """Calculate the treeness/RCV of the tree by PhyKIT implementation."""
        from ..external._libphykit import compute_toverr

        self._toverr = compute_toverr(self.file, self.tree)

    @_check_attributes("_tree")
    def compute_saturation(self) -> None:
        """Calculate the saturation of the tree by PhyKIT implementation."""
        from ..external._libphykit import Saturation

        self._saturation = Saturation().compute_saturation(self._data, self.tree)

    def _guess_seqtype(self) -> str:
        """Guess and return the sequence type."""
        f = gzip.open(self.file, "rt") if is_gzip_file(self.file) else open(self.file)

        for r in FastaIO.SimpleFastaParser(f):
            seqtype = guess_seqtype(r[1], ignore_failed=True)
            if seqtype:
                f.close()
                return seqtype


class MFA2TreeList(_abc.SeqDataListABC[MFA2Tree]):
    """A wrapper for managing MFA2Tree objects, providing tree-building, concatenation and consensus tree calculations.

    This class allows working with multiple MFA2Tree instances, providing methods to build trees, compute metrics, and handle MSA
    concatenation efficiently with parallelization support.
    """

    _bound_class = MFA2Tree

    @overload
    def __init__(self, data: Sequence[str | Path | MFA2Tree], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO") -> None: ...

    @overload
    def __init__(
        self, data: Sequence[str | Path | MFA2Tree], names: Sequence[str], *, seqtype: Literal["dna", "pep", "AUTO"] = "AUTO"
    ) -> None: ...

    def __init__(
        self,
        data: Sequence[str | Path | MFA2Tree] | None = None,
        names: Sequence[str] | None = None,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ) -> None:
        """Initialize the MFA2TreeList wrapper object.

        Args:
            data (Sequence[str | Path | MFA2Tree] | None, optional): A sequence of data items.
            names (Sequence[str] | None, optional): A sequence of names corresponding to the data items.
            seqtype (Literal["dna", "pep", "AUTO"]): The sequence type of the file. Defaults to AUTO.

        Raises:
            RuntimeError: If names are provided but data is not.
            TypeError: If a data item cannot be converted to a MFA2Tree.
            KeyError: If the item already exists.
        """
        super().__init__(data, names, seqtype=seqtype)

    @overload
    def __getitem__(self, key: str) -> MFA2Tree: ...

    @overload
    def __getitem__(self, key: int) -> MFA2Tree: ...

    @overload
    def __getitem__(self, key: slice) -> MFA2TreeList: ...

    def __getitem__(self, key: str | int | slice) -> MFA2Tree | MFA2TreeList:
        """Retrieves an item or subset of items by name, index, or slice.

        Args:
            key (str | int | slice): The key to retrieve.

        Returns:
            MFA2Tree | MFA2TreeList: The corresponding item or subset of items.
        """
        return super().__getitem__(key)

    @property
    def trees(self) -> list[Tree]:
        """List of phylogenetic trees built from multiple sequence alignments.

        Returns:
            list[Tree]: A list of Biopython Tree objects.
        """
        return [mfa.tree for mfa in self]

    @property
    def method(self) -> str:
        """Return the method used for tree building.

        Raises:
            RuntimeError: If multiple methods are detected among the trees.

        Returns:
            str: The tree building method.
        """
        method = {mfa.method for mfa in self}
        if len(method) > 1:
            raise RuntimeError("Not all the trees built with the same method.")
        return next(iter(method))

    @property
    def toverrs(self) -> list[float]:
        """List of treeness/RCV values for each MFA2Tree.

        Returns:
            list[float]: Treeness/RCV values.
        """
        return [mfa.toverr for mfa in self]

    @property
    def saturations(self) -> list[float]:
        """List of saturation values for each MFA2Tree.

        Returns:
            list[float]: Saturation values.
        """
        return [mfa.saturation for mfa in self]

    def sort(self) -> None:
        """Sort the MFA2Tree objects by treeness/RCV in descending order."""
        super().sort(reverse=True)

    @Timer.timer
    def build(
        self,
        method: Literal["ft", "raxml", "iqtree"],
        model: str | Path = "AUTO",
        *,
        noml: bool = False,
        bs: int = 0,
        scfl: int = 0,
        seed: int = -1,
        jobs: int = AVAIL_CPUS,
        threads: int = 1,
    ) -> None:
        """Build trees for each MFA2Tree object.

        Args:
            method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
            noml (bool): Disable Maximum-likelihood estimation to speed up tree building. Only available for FastTree. Defaults to
                False.
            bs (int): UFBoot value. Must be at least 1000 if set. Defaults to 0 (disabled).
            scfl (int): Number of quartets to randomly sample around each internal branch when computing the site concordance
                factor (sCF). A minimum of 100 is recommended for stable results. Defaults to 0 (disabled).
            seed (int): Seed number for reproducibility. Defaults to -1 (random seed).
            jobs (int, optional): Number of parallel processes. Defaults to maximum available CPUs.
            threads (int, optional): Number of threads per process. Defaults to 1.
        """
        step_size = min(max(10, len(self) // 200 * 50), 500)
        total = len(self)
        params_per_task = [(mfa2tree, method, None, model, noml, bs, scfl, seed, threads) for mfa2tree in self]
        try:
            if len(self) == 1 or jobs <= 1:
                logger.debug("Sequential mode with %s threads.", threads)
                for i, params in enumerate(params_per_task, 1):
                    _build_helper(*params)
                    if i % step_size == 0 or i == total:
                        logger.info("Progress: %d / %d", i, total)
            else:
                logger.debug("Multiprocesses mode with %s jobs and %s threads for each.", jobs, threads)
                with ThreadPool(jobs) as pool:
                    async_tasks = [pool.apply_async(_build_helper, p) for p in params_per_task]
                    for i, r in enumerate(async_tasks, 1):
                        r.get()
                        if i % step_size == 0 or i == total:
                            logger.info("Progress: %d / %d", i, total)
        except Exception:
            logger.error("%s", traceback.format_exc())
            raise

    @Timer.timer
    def compute_toverr(self, *, threads: int = 1) -> None:
        """Compute treeness/RCV for each tree.

        Args:
            threads (int, optional): Number of parallel threads. Defaults to 1.
        """
        logger.info("Calculating treeness/RCV...")
        if len(self) == 1 or threads == 1:
            logger.debug("Sequential mode with maximum %s threads.", threads)
            for mfa2tree in self:
                _compute_toverr_helper(mfa2tree)
        else:
            logger.debug("Multiprocesses mode: %s jobs are run concurrently.", threads)
            with ThreadPool(threads) as pool:
                pool.map(_compute_toverr_helper, self)
        logger.info("Done.")

    @Timer.timer
    def compute_saturation(self, *, threads: int = 1) -> None:
        """Compute saturation values for each tree.

        Args:
            threads (int, optional): Number of parallel threads. Defaults to 1.
        """
        logger.info("Calculating saturation...")
        if len(self) == 1 or threads == 1:
            logger.debug("Sequential mode with maximum %s threads.", threads)
            for mfa2tree in self:
                _compute_saturation_helper(mfa2tree)
        else:
            logger.debug("Multiprocesses mode: %s jobs are run concurrently.", threads)
            with ThreadPool(threads) as pool:
                pool.map(_compute_saturation_helper, self)
        logger.info("Done.")

    @overload
    def get_consensus_tree(self, output: None, *, seed: int = -1, threads: int = 1) -> Tree: ...

    @overload
    def get_consensus_tree(self, output: str | Path, *, seed: int = -1, threads: int = 1) -> Path: ...

    @Timer.timer
    @_check_single_file
    def get_consensus_tree(self, output: str | Path | None = None, *, seed: int = -1, threads: int = 1) -> Tree | Path:
        """Compute a consensus tree using ASTRAL.

        Args:
            output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
            seed (int): Seed number for reproducibility. Defaults to -1 (random seed).
            threads (int): Number of threads to use. Applicable for RAxML-NG and IQ-TREE. Defaults to 1.

        Returns:
            Tree | Path:
                - If `output` is None (default), returns a `Bio.Phylo.BaseTree.Tree` object.
                - If `output` is str or Path, returns a `pathlib.Path` to the output tree file.

        Raises:
            BinaryNotFoundError: If ASTRAL is not installed.
            RuntimeError: If ASTRAL fails.
        """
        from ..external import Astral

        tempdir = None
        try:
            if not output:
                tempdir = tempfile.TemporaryDirectory()
                output = tempdir.name

            output = Path(output)
            output.mkdir(parents=True, exist_ok=True)

            with open(output / "trees.nw", "w") as f:
                Phylo.write(self.trees, f, "newick")
            runner = Astral(output / "trees.nw", output / "consensus.nw", seed=seed, threads=threads)
            runner.run()
        finally:
            if tempdir:
                tempdir.cleanup()
        logger.info("Consensus tree is done.")
        if tempdir:
            return Phylo.read(runner.result, "newick")
        else:
            return runner.result

    @Timer.timer
    @_check_single_file
    def concat(self, output: str | Path, *, threads: int = 1) -> tuple[Path, Path]:
        """Concatenate selected MSAs into a single file and generate an MFA2Tree.

        Args:
            output (str | Path): Path to save the concatenated MSA file.
            partition (bool, optional): Whether to create a partition file. Defaults to False.
            threads (int, optional): Number of threads for parallel processing. Defaults to 1.

        Returns:
            MFA2Tree: An MFA2Tree object representing the concatenated result.

        Raises:
            TypeError: If the output is not a string or Path.
        """
        from ..external import PartitionRecord, Partitions, RaxmlHandler

        self.load()
        try:
            # Sort by MSA length with descending order (long -> short)
            if isinstance(output, (str, Path)):
                output = Path(output)
            else:
                raise TypeError("Argument output only accepts str or Path.")
            mfa2treelist = sorted(self, key=lambda msa: len(msa), reverse=True)

            # Get sample names and create an empty MSA holder
            concat_alignments = MultipleSeqAlignment([])
            samples = set()
            msa_list = []
            for msa in mfa2treelist:
                msa_list.append(msa._data)
                for rec in msa._data:
                    samples.add(rec.id)
            samples = tuple(samples)
            for sample in samples:
                concat_alignments.append(SeqRecord(Seq(""), id=sample, description=""))
            concat_alignments.sort()
        finally:
            self.unload()

        # Fill missing taxon on some of the MSAs
        logger.info("Filling missing taxon...")
        if threads > 1:
            with ThreadPool(threads) as pool:
                msa_list = pool.starmap(_fill_missing_taxon, product([samples], msa_list))
        else:
            msa_list = [_fill_missing_taxon(samples, msa) for msa in msa_list]
        logger.info("Done.")

        # Concatenate
        logger.info("Concatenating selected MSAs...")

        start_idx = 0
        partition_info = Partitions("raxml")
        model = "GTR" if self.seqtype == SeqTypes.DNA else "LG"
        for msa in msa_list:
            concat_alignments += msa
            end_idx = start_idx + msa.get_alignment_length()
            partition_info.add(
                PartitionRecord(
                    name=msa.annotations["seqname"],
                    start=start_idx + 1,
                    end=end_idx,
                    model=model,
                )
            )
            start_idx = end_idx

        for seq in concat_alignments:
            seq.description = ""
        logger.info("Done.")

        output.mkdir(exist_ok=True)
        concat_file = output / TreeOutputFiles.CONCAT
        with open(concat_file, "w") as f:
            SeqIO.write(concat_alignments, f, format="fasta")
            logger.info("Concatenated fasta is output to %s", concat_file)

        partition_file = output / TreeOutputFiles.PARTITION
        with RaxmlHandler(partition_file, "w") as f:
            f.write(partition_info)
        logger.info("Partition file is output to %s", partition_file)
        return concat_file, partition_file


@overload
def bootstrap(
    mfa2tree: MFA2Tree,
    tree: str | Path,
    output: None,
    model: str | Path = "AUTO",
    *,
    bs: int = 0,
    seed: int = -1,
    threads: int = -1,
    threads_max: int = AVAIL_CPUS,
) -> Tree: ...


@overload
def bootstrap(
    mfa2tree: MFA2Tree,
    tree: str | Path,
    output: str | Path,
    model: str | Path = "AUTO",
    *,
    bs: int = 0,
    seed: int = -1,
    threads: int = -1,
    threads_max: int = AVAIL_CPUS,
) -> Path: ...


def bootstrap(
    mfa2tree: MFA2Tree,
    tree: str | Path,
    output: str | Path | None = None,
    model: str | Path = "AUTO",
    *,
    bs: int = 0,
    seed: int = -1,
    threads: int = -1,
    threads_max: int = AVAIL_CPUS,
) -> Tree | Path:
    """
    Performs bootstrap analysis on a phylogenetic tree using UFBoot.

    Args:
        mfa2tree (MFA2Tree): The object contains sequence alignment.
        tree: Path to the tree file.
        output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
        model (str | Path): String of model or path to a partition file for model partitioning. Defaults to "AUTO" that use
            ModelFinder to identify the best model.
        bs (int): UFBoot value. Must be at least 1000 if set. Defaults to 0 (disabled).
        seed (int): Seed number for reproducibility. Defaults to -1 (random seed).
        threads (int): Number of threads to use. Defaults to -1 (autodetect).
        threads_max (int): Maximum number of threads used during autodetect. Defaults to maximum available CPUs.

    Returns:
        Tree | Path:
            - If `output` is None (default), returns a `Bio.Phylo.BaseTree.Tree` object.
            - If `output` is str or Path, returns a `pathlib.Path` to the output tree file.
    """
    from ..external import UFBoot

    tempdir = None
    try:
        if not output:
            tempdir = tempfile.TemporaryDirectory()
            output = tempdir.name

        output = Path(output)
        output.mkdir(exist_ok=True)

        runner = UFBoot(mfa2tree.file, tree, output, model=model, bs=bs, seed=seed, threads=threads, threads_max=threads_max)
        runner.run()

    finally:
        if tempdir:
            tempdir.cleanup()
    logger.info("Bootstrapping on %s is done.", mfa2tree.name)
    if tempdir:
        return Phylo.read(runner.result, "newick")
    else:
        return runner.result


@overload
def branch_concordance(
    mfa2tree: MFA2Tree,
    tree: str | Path,
    output: None,
    model: str | Path = "AUTO",
    *,
    scfl: int = 0,
    seed: int = -1,
    threads: int = -1,
    threads_max: int = AVAIL_CPUS,
) -> Tree: ...


@overload
def branch_concordance(
    mfa2tree: MFA2Tree,
    tree: str | Path,
    output: str | Path,
    model: str | Path = "AUTO",
    *,
    scfl: int = 0,
    seed: int = -1,
    threads: int = -1,
    threads_max: int = AVAIL_CPUS,
) -> Path: ...


def branch_concordance(
    mfa2tree: MFA2Tree,
    tree: str | Path,
    output: str | Path | None = None,
    model: str | Path = "AUTO",
    *,
    scfl: int = 0,
    seed: int = -1,
    threads: int = -1,
    threads_max: int = AVAIL_CPUS,
) -> Tree | Path:
    """
    Performs site concordance factor calculation on a phylogenetic tree using UFBoot.

    Args:
        mfa2tree (MFA2Tree): The object contains sequence alignment.
        tree: Path to the tree file.
        output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
        model (str | Path): String of model or path to a partition file for model partitioning. Defaults to "AUTO" that use
            ModelFinder to identify the best model.
        scfl (int): Number of quartets to randomly sample around each internal branch when computing the site concordance factor
            (sCF). A minimum of 100 is recommended for stable results. Defaults to 0 (disabled).
        seed (int): Seed number for reproducibility. Defaults to -1 (random seed).
        threads (int): Number of threads to use. Defaults to -1 (autodetect).
        threads_max (int): Maximum number of threads used during autodetect. Defaults to maximum available CPUs.

    Returns:
        Tree | Path:
            - If `output` is None (default), returns a `Bio.Phylo.BaseTree.Tree` object.
            - If `output` is str or Path, returns a `pathlib.Path` to the output tree file.
    """
    from ..external import Concordance

    tempdir = None
    try:
        if not output:
            tempdir = tempfile.TemporaryDirectory()
            output = tempdir.name

        output = Path(output)
        output.mkdir(exist_ok=True)

        runner = Concordance(
            mfa2tree.file, tree, output, model=model, scfl=scfl, seed=seed, threads=threads, threads_max=threads_max
        )
        runner.run()

    finally:
        if tempdir:
            tempdir.cleanup()
    logger.info("Site concordance factor calculation on %s is done.", mfa2tree.name)
    if tempdir:
        return Phylo.read(runner.result, "newick")
    else:
        return runner.result


def _compute_toverr_helper(instance: MFA2Tree) -> None:
    """Helper function to run the `compute_toverr` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `compute_toverr` method.
    """
    instance.compute_toverr()


def _compute_saturation_helper(instance: MFA2Tree) -> None:
    """Helper function to run the `compute_saturation` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `compute_saturation` method.
    """
    instance.compute_saturation()


def _build_helper(
    instance: MFA2Tree,
    method: Literal["ft", "raxml", "iqtree"],
    output: str | Path | None = None,
    model: str | Path = "AUTO",
    noml: bool = False,
    bs: int = 0,
    scfl: int = 0,
    seed: int = -1,
    threads: int = 1,
) -> Tree:
    """Helper function to run the `build` method on an MFA2Tree instance.

    Args:
        instance (MFA2Tree): The MFA2Tree instance on which to run the `build` method.
        method (Literal["ft", "raxml", "iqtree"]): Tree-building method.
                - "ft": Use FastTree.
                - "raxml": Use RAxML-NG.
                - "iqtree": Use IQ-TREE.
        output (str | Path | None, optional): Path to save the resulting tree file. If None, a temporary path is used.
        model (str | Path | None, optional): String of model or path to a partition file for model partitioning. Defaults to
            "AUTO" that use ModelFinder to identify the best model.
        noml (bool): Disable Maximum-likelihood estimation to speed up tree building. Only available for FastTree. Defaults to
            False.
        bs (int): UFBoot value. Defaults to 0 (disabled).
        scfl (int): Site concordance factor threshold for IQ-TREE. Must be at least 1000 if set. Defaults to 0 (disabled).
        seed (int): Seed number for reproducibility. Defaults to -1 (random seed).
        threads (int): Maximum number of threads to use. Only applicable to RAxML-NG and IQ-TREE. Defaults to 1.

    Returns:
        Tree: The resulting phylogenetic tree.
    """
    mfa2tree_logger = logging.getLogger("phyling.lib.tree.MFA2Tree")
    original_level = mfa2tree_logger.level
    mfa2tree_logger.setLevel(logging.WARNING)
    try:
        instance.build(method, output, model, noml=noml, bs=bs, scfl=scfl, seed=seed, threads_max=threads)
    finally:
        mfa2tree_logger.setLevel(original_level)


def _fill_missing_taxon(samples: Sequence[str], msa: MultipleSeqAlignment) -> MultipleSeqAlignment:
    """Add gap-only sequences to the alignment for missing taxa.

    Args:
        samples (Sequence[str]): List of all expected taxa.
        msa (MultipleSeqAlignment): The existing multiple sequence alignment.

    Returns:
        MultipleSeqAlignment: The updated alignment with added gap-only sequences for missing taxa.
    """
    missing = set(samples) - {seq.id for seq in msa}
    for sample in missing:
        msa.append(
            SeqRecord(
                Seq("-" * msa.get_alignment_length()),
                id=sample,
                description="",
            )
        )
        msa.sort()
    return msa
