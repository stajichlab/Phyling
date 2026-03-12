"""Muscle utilities"""

from __future__ import annotations

from pathlib import Path

from ..lib._utils import check_binary
from ._abc import BinaryWrapper

MUSCLE_BIN = check_binary("Muscle", ("muscle",), "bioconda::muscle", "https://github.com/rcedgar/muscle")


class Muscle(BinaryWrapper):
    _prog: str = "Muscle"

    def __init__(self, file: str | Path, output: str | Path, *, threads: int = 1):
        super().__init__(file, output, threads=threads)
        print(self._output)

    def _construct_cmd(self, file: Path, output: Path, *, threads: int = 1):
        self._cmd = [
            MUSCLE_BIN,
            "-align",
            str(file),
            "-output",
            str(output),
            "-threads",
            str(threads),
        ]
