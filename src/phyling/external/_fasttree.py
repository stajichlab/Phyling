"""FastTree utilities"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Literal

from ..lib import TreeMethods
from ..lib._utils import check_binary
from ._abc import TreeToolWrapper
from ._models import DNA_MODELS, PEP_MODELS


def _get_bin():
    return check_binary(
        TreeMethods.FT.method, TreeMethods.FT.bins, "bioconda::fasttree", "https://github.com/morgannprice/fasttree"
    )


class FastTree(TreeToolWrapper[Literal["DNA", "AA"]]):
    """Runs FastTree to build a phylogenetic tree from the given MFA2Tree object.

    Args:
        mfa2tree (MFA2Tree): The MFA2Tree object containing multiple sequence alignment data.
        capture_cmd (bool): If True, returns the FastTree command along with the resulting tree.

    Returns:
        If capture_cmd is False (default), returns a Tree object.
        If capture_cmd is True, returns a tuple containing the Tree object and the command string.
    """

    _prog: str = TreeMethods.FT.method
    _cmd_log = "stderr"
    _ALLOWED_SEQTYPES: tuple[str, ...] = ("DNA", "AA")

    def __init__(
        self,
        file: str | Path,
        output: str | Path,
        *,
        seqtype: Literal["dna", "pep"],
        model: str = "AUTO",
        seed: int = -1,
        noml: bool = False,
    ):
        super().__init__(file, output, seqtype=seqtype, model=model)

        self._construct_cmd(seed=seed, noml=noml)

    def _construct_cmd(self, *, noml: bool, seed: int):
        self._cmd = [
            _get_bin(),
            "-nosupport",
            "-out",
            str(self._output),
            str(self._file),
        ]
        model, *params = self._model.split("+")
        if self._seqtype == "DNA":
            self._cmd.insert(1, "-nt")
            if model.upper() not in DNA_MODELS[DNA_MODELS[:, 0] != "", 0]:
                raise ValueError(f"Model {model} is not supported in {TreeMethods.FT.method} with {self._seqtype} alignments.")
            if model.upper() == "GTR":
                self._cmd.insert(3, f"-{model.lower()}")
        else:
            if model.upper() not in PEP_MODELS[PEP_MODELS[:, 0] != "", 0]:
                raise ValueError(f"Model {model} is not supported in {TreeMethods.FT.method} with {self._seqtype} alignments.")
            if model.upper() in ("LG", "WAG"):
                self._cmd.insert(2, f"-{model.lower()}")

        if re.match("G", "+".join(params)):
            self._cmd.insert(-1, "-gamma")

        if noml:
            self._cmd.insert(-1, "-noml")

        if seed >= 0:
            self._cmd[-1:-1] = ["-seed", str(seed)]
