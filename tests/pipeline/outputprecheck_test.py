"""Tests for the output precheck module."""

from __future__ import annotations

from unittest.mock import MagicMock

from phyling.lib.align import SampleList
from phyling.pipeline._outputprecheck import AlignPrecheck


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
