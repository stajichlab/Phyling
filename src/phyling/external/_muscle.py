"""Muscle utilities"""

from __future__ import annotations

from pathlib import Path

from ..lib._utils import check_binary
from ._abc import BinaryWrapper


def _get_bin():
    return check_binary("Muscle", ("muscle",), "bioconda::muscle", "https://github.com/rcedgar/muscle")


class Muscle(BinaryWrapper[Path]):
    _prog: str = "Muscle"

    def __init__(self, file: str | Path, output: str | Path, *, threads: int = 1):
        super().__init__(file, output)
        self._construct_cmd(threads=threads)

    def _construct_cmd(self, *, threads: int = 1):
        self._cmd = [
            _get_bin(),
            "-align",
            str(self._file),
            "-output",
            str(self._output),
            "-threads",
            str(threads),
        ]
