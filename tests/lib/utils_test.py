"""Tests for the _utils module library."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from phyling.exception import SeqtypeError
from phyling.lib._utils import guess_file_seqtype


class TestGuessFileSeqtype:
    def test_skips_unrecognized_records(self, tmp_path: Path):
        fasta = tmp_path / "seqs.fa"
        fasta.write_text(">a\n12345\n>b\nMKLV\n")
        assert guess_file_seqtype(fasta) == "pep"

    def test_gzipped_file(self, tmp_path: Path):
        fasta = tmp_path / "seqs.fa.gz"
        with gzip.open(fasta, "wt") as f:
            f.write(">a\nACGT\n")
        assert guess_file_seqtype(fasta) == "dna"

    @pytest.mark.parametrize(
        "content, match",
        [
            ("", "No sequences found"),
            (">a\n12345\n", "Cannot determine seqtype"),
            (">a\nACGU\n", "Invalid seqtype"),
        ],
    )
    def test_raises(self, tmp_path: Path, content: str, match: str):
        fasta = tmp_path / "seqs.fa"
        fasta.write_text(content)
        with pytest.raises(SeqtypeError, match=match):
            guess_file_seqtype(fasta)
