"""Tests for the output precheck module."""

from __future__ import annotations

from unittest.mock import MagicMock

from phyling.lib.align import SampleList
from phyling.lib.tree import MFA2Tree, MFA2TreeList
from phyling.pipeline._outputprecheck import AlignPrecheck, FilterPrecheck


class TestAlignPrecheck:
    def test_cutoff_change_removes_old_outputs(self, tmp_path):
        (tmp_path / AlignPrecheck.ckp).write_bytes(b"")
        stale = tmp_path / "old_marker.mfa"
        stale.write_text(">a\nACGT\n")
        data = MagicMock(spec=SampleList)
        data.seqtype = "pep"
        data.checksums = {}
        prev_data = MagicMock()
        prev_data.samplelist.seqtype = "pep"
        prev_data.samplelist.checksums = {}
        precheck = AlignPrecheck(tmp_path, data, markerset="m", markerset_cutoff="new", method="hmmalign")

        precheck._determine_rerun({"markerset": "m", "markerset_cutoff": "old", "method": "hmmalign"}, prev_data)

        assert not stale.exists()
        assert (tmp_path / AlignPrecheck.ckp).exists()

    def test_params_are_not_shared_between_instances(self, tmp_path):
        data = MagicMock(spec=SampleList)
        AlignPrecheck(tmp_path, data, markerset="first", markerset_cutoff=1, method="hmmalign")
        second = AlignPrecheck(tmp_path, data, method="muscle")

        assert second.params["markerset"] is None
        assert AlignPrecheck.params["markerset"] is None


class TestFilterPrecheck:
    @staticmethod
    def _msa(name: str, checksum: int, toverr: float | None = 0.5) -> MagicMock:
        msa = MagicMock(spec=MFA2Tree)
        msa.name, msa.checksum, msa.toverr = name, checksum, toverr
        return msa

    def test_reruns_when_inputs_change(self, tmp_path):
        (tmp_path / FilterPrecheck.ckp).write_bytes(b"")
        (tmp_path / "other_file").write_text("")
        data = MagicMock(spec=MFA2TreeList)
        data.seqtype = "pep"
        data.__iter__.return_value = [self._msa("m1", 2)]
        prev_data = MagicMock()
        prev_data.seqtype = "pep"
        prev_data.__iter__.return_value = [self._msa("m1", 1)]
        prev_data.__contains__.return_value = True
        prev_data.__getitem__.return_value = self._msa("m1", 1)
        precheck = FilterPrecheck(tmp_path, data, top_n_toverr=10)

        remained, completed = precheck._determine_rerun({"top_n_toverr": 10}, prev_data)

        assert [msa.name for msa in remained] == ["m1"]
