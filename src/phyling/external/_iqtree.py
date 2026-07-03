"""IQ-Tree utilities"""

from __future__ import annotations

import warnings

from phyling.exception import SeqtypeError

warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import gzip
import re
from pathlib import Path
from typing import Literal, overload

from Bio import SeqIO

from .. import AVAIL_CPUS
from ..lib import SeqTypes, TreeMethods
from ..lib._utils import check_binary, guess_seqtype, is_gzip_file
from ._abc import BinaryWrapper, TreeToolWrapper
from ._models import (
    ALL_MODELS,
    DNA_MODELS,
    INVARIANT_CODES,
    PEP_MODELS,
    STATIONARY_CODES,
    NexusHandler,
    RaxmlHandler,
)


def _get_bin():
    return check_binary(
        TreeMethods.IQTREE.method, TreeMethods.IQTREE.bins, "bioconda::iqtree", "https://github.com/iqtree/iqtree3"
    )


class ModelFinder(BinaryWrapper[Path]):
    _prog: str = "ModelFinder"

    @overload
    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        method: Literal["ft", "raxml", "iqtree"] = "iqtree",
        seed: int = -1,
        threads: int = 1,
        threads_max: int = AVAIL_CPUS,
    ) -> None: ...
    @overload
    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        partition_file: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        method: Literal["ft", "raxml", "iqtree"] = "iqtree",
        seed: int = -1,
        threads: int = 1,
        threads_max: int = AVAIL_CPUS,
    ) -> None: ...
    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        partition_file: str | Path | None = None,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        method: Literal["ft", "raxml", "iqtree"] = "iqtree",
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
    ) -> None:
        super().__init__(file, output)

        if partition_file is not None:
            if not isinstance(partition_file, (str, Path)):
                raise TypeError(f"Argument partition only accepts str or Path. Got {type(partition_file)}.")
            if method == "ft":
                raise ValueError(f"Partitioning analysis is not allowed when using {TreeMethods.FT.method}.")
            partition_file = Path(partition_file)
        method_idx = list(TreeMethods).index(TreeMethods[method.upper()])
        if seqtype == "AUTO" or seqtype != "DNA" or seqtype != "AA":
            f = gzip.open(file, "rt") if is_gzip_file(file) else open(file)

            for r in SeqIO.FastaIO.SimpleFastaParser(f):
                seqtype_ = guess_seqtype(r[1])
                if seqtype_ in (SeqTypes.DNA, SeqTypes.PEP):
                    break
            if seqtype_ == SeqTypes.RNA:
                raise SeqtypeError(f"Invalid seqtype: {seqtype}.")
            f.close()
            seqtype = seqtype_
        if seqtype == SeqTypes.DNA:
            seqtype_ = "DNA"
            # Find the support models for each tool and map to the name in IQTree format
            mset = ",".join(DNA_MODELS[DNA_MODELS[:, method_idx] != "", 2])
        else:
            seqtype_ = "AA"
            mset = ",".join(PEP_MODELS[PEP_MODELS[:, method_idx] != "", 2])
        self._method = method

        self._construct_cmd(
            partition_file=partition_file, seqtype=seqtype_, mset=mset, seed=seed, threads=threads, threads_max=threads_max
        )

    def _post_run(self) -> None:
        if not self._output:
            raise RuntimeError("No output file was generated.")
        if self._output.suffix == ".nex":  # partitioning analysis
            if self._method == "raxml":
                with NexusHandler(self._output) as fi, RaxmlHandler(self._output.with_suffix(""), mode="w") as fo:
                    part_info = fi.read()["sets"]
                    part_info = part_info.convert_to("raxml")
                    fo.write(part_info)
                self._result = self._output.with_suffix("")
            else:
                self._result = self._output
        else:
            model, params = [], []
            with gzip.open(self._output, "rt") as f:
                for line in f.read().strip("\n").split("\n"):
                    best_model_prefix = "best_model_BIC: "
                    if line.startswith(best_model_prefix):
                        best_model = line.lstrip(best_model_prefix)
                        model, *params = best_model.split("+")
                        method_idx = list(TreeMethods).index(TreeMethods[self._method.upper()])
                        # Convert model in IQTree name to the name of each tool
                        model = ALL_MODELS[ALL_MODELS[:, 2] == model, method_idx].tolist()
                        for param in params:
                            if self._method == "ft":
                                if param.startswith(("F", "I")):
                                    params
                            else:
                                if param.startswith("F"):
                                    param = STATIONARY_CODES[STATIONARY_CODES[:, 2] == param, method_idx][0].item()
                                elif param.startswith("F"):
                                    param = INVARIANT_CODES[INVARIANT_CODES[:, 2] == param, method_idx][0].item()
            self._result = "+".join(model + params)

    def _construct_cmd(
        self,
        *,
        partition_file: Path | None,
        seqtype: Literal["AA", "DNA"],
        mset: str,
        seed: int,
        threads: int,
        threads_max: int,
    ):
        self._cmd = [
            _get_bin(),
            "-s",
            str(self._file.absolute()),
            "--prefix",
            str(self._output.absolute()),
            "-T",
            str(threads) if threads >= 1 else "AUTO",
            "--threads-max",
            str(threads_max),
            "-m",
            "TESTONLY",
        ]
        self._cmd.extend(["--seqtype", seqtype])
        self._cmd.extend(["--mset", mset])
        if partition_file:
            self._cmd.extend(["-p", str(partition_file)])
            self._output = self._output.with_suffix(".best_scheme.nex")
        else:
            self._output = self._output.with_suffix(".model.gz")
        if seed >= 0:
            self._cmd.extend(["--seed", str(seed)])


class Iqtree(TreeToolWrapper[Literal["DNA", "AA", "AUTO"]]):
    _prog: str = TreeMethods.IQTREE.method
    _ALLOWED_SEQTYPES: tuple[str, ...] = ("DNA", "AA", "AUTO")

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
        model: str = "AUTO",
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
    ) -> None:
        super().__init__(file, output, seqtype=seqtype, model=model)

        self._construct_cmd(seed=seed, threads=threads, threads_max=threads_max)

    def _post_run(self) -> None:
        if not self._output:
            raise RuntimeError("No output file was generated.")
        model_file = self._output.with_suffix(".best_model.nex")

        if model_file.is_file():
            self._model = str(model_file)
        else:
            with open(self._output.with_suffix(".iqtree")) as f:
                if match := re.search(r"alisim simulated_MSA .* (\-m) \"(.*)\" ", f.read()):
                    self._model = match[2]

    def _construct_cmd(self, *, seed: int, threads: int, threads_max: int):
        self._cmd = [
            _get_bin(),
            "-s",
            str(self._file.absolute()),
            "--prefix",
            str(self._output.absolute()),
            "-T",
            str(threads) if threads >= 1 else "AUTO",
            "--threads-max",
            str(threads_max),
        ]
        self._cmd.extend(["--seqtype", self._seqtype])
        if seed >= 0:
            self._cmd.extend(["--seed", str(seed)])
        if Path(self._model).is_file():
            self._cmd.extend(["-p", self._model])
        else:
            self._cmd.extend(["-m", self._model])

        self._output = self._output.with_suffix(".treefile")


class UFBoot(TreeToolWrapper[Literal["DNA", "AA", "AUTO"]]):
    _prog: str = "UFBoot"
    _ALLOWED_SEQTYPES: tuple[str, ...] = ("DNA", "AA", "AUTO")

    def __init__(
        self,
        file: str | Path,
        tree: str | Path,
        output: str | Path,
        *,
        model: str = "AUTO",
        bs: int = 1000,
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
    ):
        super().__init__(file, output, seqtype="AUTO", model=model)

        tree = Path(tree)
        if not tree.exists():
            raise FileNotFoundError(f"{tree}")
        if not tree.is_file():
            raise RuntimeError(f"{tree} is not a file.")

        self._construct_cmd(tree=tree, bs=bs, seed=seed, threads=threads, threads_max=threads_max)

    def _construct_cmd(self, *, tree: Path, bs: int, seed: int, threads: int, threads_max: int):
        self._cmd = [
            _get_bin(),
            "-s",
            str(self._file.absolute()),
            "--prefix",
            str(self._output.absolute()),
            "-t",
            str(tree.absolute()),
            "-B",
            str(bs),
            "-bnni",
            "-T",
            str(threads) if threads >= 1 else "AUTO",
            "--threads-max",
            str(threads_max),
        ]
        if seed >= 0:
            self._cmd.extend(["--seed", str(seed)])
        if Path(self._model).is_file():
            self._cmd.extend(["-p", self._model])
        else:
            self._cmd.extend(["-m", self._model])

        self._output = self._output.with_suffix(".treefile")


class Concordance(TreeToolWrapper[Literal["DNA", "AA", "AUTO"]]):
    _prog: str = "Branch concordance calculation"
    _ALLOWED_SEQTYPES: tuple[str, ...] = ("DNA", "AA", "AUTO")

    def __init__(
        self,
        file: str | Path,
        tree: str | Path,
        output: str | Path,
        *,
        model: str = "AUTO",
        scfl: int = 100,
        seed: int = -1,
        threads: int = -1,
        threads_max: int = AVAIL_CPUS,
    ):
        super().__init__(file, output, seqtype="AUTO", model=model)

        tree = Path(tree)
        if not tree.exists():
            raise FileNotFoundError(f"{tree}")
        if not tree.is_file():
            raise RuntimeError(f"{tree} is not a file.")

        self._construct_cmd(tree=tree, scfl=scfl, seed=seed, threads=threads, threads_max=threads_max)

    def _construct_cmd(self, *, tree: Path, scfl: int, seed: int, threads: int, threads_max: int):
        self._cmd = [
            _get_bin(),
            "-s",
            str(self._file.absolute()),
            "--prefix",
            str(self._output.absolute()),
            "-te",
            str(tree.absolute()),
            "--scfl",
            str(scfl),
            "-T",
            str(threads) if threads >= 1 else "AUTO",
            "--threads-max",
            str(threads_max),
        ]
        if seed >= 0:
            self._cmd.extend(["--seed", str(seed)])
        if Path(self._model).is_file():
            self._cmd.extend(["-p", self._model])
        else:
            self._cmd.extend(["-m", self._model])

        self._output = self._output.with_suffix(".cf.tree")
