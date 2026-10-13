"""Tests for the external tree tool wrappers."""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from phyling.external import _fasttree, _iqtree, _raxml
from phyling.external._models import PartitionRecord


class TestModelFinder:
    @staticmethod
    def _post_run_with_result(tmp_path: Path, method: str, best_model: str) -> _iqtree.ModelFinder:
        output = tmp_path / "out.model.gz"
        with gzip.open(output, "wt") as f:
            f.write(f"some_line: 1\nbest_model_BIC: {best_model}\n")
        runner = object.__new__(_iqtree.ModelFinder)
        runner._output = output
        runner._method = method
        runner._post_run()
        return runner

    @pytest.mark.parametrize(
        "method, best_model, expected",
        [
            ("iqtree", "BLOSUM62+G4", "BLOSUM62+G4"),  # lstrip used to turn this into "LOSUM62"
            ("raxml", "LG+FQ+I+G4", "LG+FE+I+G4"),  # stationary code converted to the RAxML-NG name
            ("ft", "LG+F+I+G4", "LG+G4"),  # FastTree has no F/I codes
        ],
    )
    def test_converts_best_model(self, tmp_path, method, best_model, expected):
        assert self._post_run_with_result(tmp_path, method, best_model)._result == expected

    def test_keeps_given_seqtype(self, tmp_path, monkeypatch):
        """A given seqtype is used as is instead of being guessed from the file."""
        msa = tmp_path / "aln.fa"
        msa.write_text(">a\nACGT\n")
        captured = {}
        monkeypatch.setattr(_iqtree.ModelFinder, "_construct_cmd", lambda self, **kw: captured.update(kw))

        _iqtree.ModelFinder(msa, tmp_path / "out", seqtype="pep")

        assert captured["seqtype"] == "AA"


class TestFastTree:
    def test_adds_gamma_when_g_is_not_first(self, tmp_path, monkeypatch):
        monkeypatch.setattr(_fasttree, "_get_bin", lambda: "FastTree")
        (tmp_path / "aln.fa").write_text(">a\nACGT\n")

        runner = _fasttree.FastTree(tmp_path / "aln.fa", tmp_path / "out.nw", seqtype="dna", model="GTR+I+G4")

        assert "-gamma" in runner._cmd


class TestRaxml:
    def test_model_skips_empty_parts(self, tmp_path, monkeypatch):
        part_info = MagicMock()
        part_info.convert_to.return_value = [PartitionRecord(name="p1", start=1, end=10, model="LG", gamma="G4")]
        handler = MagicMock()
        handler.return_value.__enter__.return_value.read.return_value = part_info
        monkeypatch.setattr(_raxml, "RaxmlHandler", handler)
        runner = object.__new__(_raxml.Raxml)
        runner._output = tmp_path / "out.raxml"

        runner._post_run()

        assert runner._model == "LG+G4"
