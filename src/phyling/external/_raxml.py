"""RAxML utilities"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from .. import AVAIL_CPUS
from ..lib import TreeMethods
from ..lib._utils import check_binary
from ._abc import TreeToolWrapper
from ._models import NexusHandler, RaxmlHandler

RAXML_BIN = check_binary(
    TreeMethods.RAXML.method, TreeMethods.RAXML.bins, "bioconda::raxml-ng", "https://github.com/amkozlov/raxml-ng"
)


class Raxml(TreeToolWrapper):
    """Runs RAxML-NG to build a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree (MFA2Tree): The MFA2Tree object containing alignment data.
        output (str | Path | None, optional): Path to save the resulting tree file. If not provided, a temporary path is used.
        partition_file (str | Path | None, optional): Path to a partition file for model partitioning. Optional.
        bs (int): Bootstrap value. Defaults to 50.
        threads (int): Number of threads to use for RAxML-NG computation.
        capture_cmd (bool): If True, returns the RAxML-NG command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple of the Tree object and the command string.
    """

    _prog: str = TreeMethods.RAXML.method

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep", "AUTO"],
        model: str = "AUTO",
        seed: int = -1,
        threads: int = 1,
        threads_max: int = AVAIL_CPUS,
    ):
        """Instantiate RAxML-NG runner.

        Args:
            file (str | Path): Path of the MSA file.
            output (str | Path): Directory to save the results.
            seqtype (str | Path | None, optional): The sequence type of the file.
            model (int): Model to use for phylogeny inference.
            threads (int): Number of threads to use.

        Returns:
            If capture_cmd is False (default), returns a Tree object.
            If capture_cmd is True, returns a tuple of the Tree object and the command string.
        """
        super().__init__(file, output, seqtype=seqtype, model=model, seed=seed, threads=threads, threads_max=threads_max)

    def _post_run(self):
        if not self._output:
            raise RuntimeError("No output file was generated.")
        model_file = self._output.with_suffix(".bestModel")

        with RaxmlHandler(model_file) as fi:
            part_info = fi.read()
            part_info = part_info.convert_to("iqtree")

        if len(part_info) > 1:
            with NexusHandler(f"{model_file}.nex", mode="w") as fo:
                fo.write({"sets": part_info})
            self._model = f"{model_file}.nex"
        else:
            info = part_info[0]
            self._model = "+".join([info.model, info.stationary, info.invariant, info.gamma])

    def _construct_cmd(
        self,
        file: Path,
        output: Path,
        *,
        seqtype: Literal["DNA", "AA"] | None,
        model: str,
        seed: int,
        threads: int,
        threads_max: int,
    ):
        self._cmd = [
            RAXML_BIN,
            "--msa",
            str(file),
            "--prefix",
            str(output),
            "--model",
            model,
            "--threads",
            str(threads) if threads >= 1 else f"auto{{{threads_max}}}",
        ]
        if seqtype:
            self._cmd.extend(["--data-type", seqtype])
        if seed >= 0:
            self._cmd.extend(["--seed", str(seed)])

        self._output = Path(f"{output}.raxml.bestTree")
