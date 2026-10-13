"""Tests for the partition models and file handlers."""

from __future__ import annotations

from pathlib import Path

import pytest

from phyling.external._models import NexusHandler


def _write_nexus(path: Path, charsets: str) -> Path:
    path.write_text(f"#nexus\nbegin sets;\n{charsets}  charpartition mymodel =\n    LG+G4: p1;\nend;\n")
    return path


class TestNexusHandler:
    @pytest.mark.parametrize("charset", ["p1 = 1-100", "p1 = concat.fa: 1-100"])
    def test_read_single_range(self, tmp_path: Path, charset: str):
        nexus = _write_nexus(tmp_path / "part.nex", f"  charset {charset};\n")
        with NexusHandler(nexus) as fi:
            (record,) = list(fi.read()["sets"])
        assert (record.name, record.start, record.end, record.model) == ("p1", 1, 100, "LG")

    @pytest.mark.parametrize("charset", ["p1 = 1-100 201-300", "p1 = 1-100\\3"])
    def test_read_unsupported_charset_raises(self, tmp_path: Path, charset: str):
        nexus = _write_nexus(tmp_path / "part.nex", f"  charset {charset};\n")
        with NexusHandler(nexus) as fi:
            with pytest.raises(RuntimeError, match="Only a single continuous range"):
                fi.read()
