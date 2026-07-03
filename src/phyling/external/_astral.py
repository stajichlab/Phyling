"""Astral utilities"""

from __future__ import annotations

from pathlib import Path

from ..lib._utils import check_binary
from ._abc import BinaryWrapper


def _get_bin():
    return check_binary("ASTRAL", ("astral",), "bioconda::aster", "https://github.com/chaoszhang/ASTER")


class Astral(BinaryWrapper[Path]):
    """Compute a consensus tree using ASTRAL.

    Returns:
        Tree: The consensus tree.

    Raises:
        BinaryNotFoundError: If ASTRAL is not installed.
        RuntimeError: If ASTRAL fails.
    """

    _prog: str = "ASTRAL"
    _cmd_log = "stderr"

    def __init__(self, file: str | Path, output: str | Path, *, seed: int = -1, threads: int = 1):
        super().__init__(file, output)

        self._construct_cmd(seed=seed, threads=threads)

    def _construct_cmd(self, *, seed: int, threads: int):
        self._cmd = [_get_bin(), "--output", str(self._output), "--thread", str(threads), str(self._file)]

        if seed >= 0:
            self._cmd[-1:-1] = ["--seed", str(seed)]
