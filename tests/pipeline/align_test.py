"""Tests for the align pipeline module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from Bio.AlignIO import MultipleSeqAlignment

import phyling.pipeline.align as align_mod
from phyling.lib.align import OrthologList, SampleList, SearchHitsManager
from phyling.pipeline.align import _args_check, _search_threads_check, align

# ---------------------------------------------------------------------------
# _args_check
# ---------------------------------------------------------------------------


class TestArgsCheck:
    def test_returns_tuple_of_paths(self, path_pep_fasta: list[Path], path_hmm_dir: Path):
        inputs_out, ms_out, ev_out, method_out = _args_check(path_pep_fasta, path_hmm_dir, 1e-10, "hmmalign")
        assert isinstance(inputs_out, tuple)
        assert all(isinstance(p, Path) for p in inputs_out)
        assert ms_out.absolute() == path_hmm_dir.absolute()
        assert ev_out == 1e-10
        assert method_out == "hmmalign"

    def test_with_hmm_name(self, path_pep_fasta: list[Path], path_hmm_dir: Path):
        markerset = "poxviridae_odb10"
        inputs_out, ms_out, ev_out, method_out = _args_check(path_pep_fasta, markerset, 1e-10, "muscle")
        assert isinstance(inputs_out, tuple)
        assert all(isinstance(p, Path) for p in inputs_out)
        assert ms_out.absolute() == path_hmm_dir.absolute()
        assert ev_out == 1e-10
        assert method_out == "muscle"

    def test_single_directory_input(self, path_pep_fasta_dir: Path, path_hmm_dir: Path):
        inputs_out, *_ = _args_check(path_pep_fasta_dir, path_hmm_dir, 1e-10, "hmmalign")
        assert len(inputs_out) == len(tuple(path_pep_fasta_dir.iterdir()))

    def test_empty_directory_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Empty input"):
            _args_check(tmp_path, "poxviridae_odb10", 1e-10, "hmmalign")

    def test_single_file_input(self, path_pep_fasta: list[Path]):
        """Pass a path of a single file rather then a single dir"""
        with pytest.raises(ValueError, match="Requires at least 4 samples"):
            _args_check(path_pep_fasta[0], "poxviridae_odb10", 1e-10, "hmmalign")

    def test_raises_on_fewer_than_four_samples(self, path_pep_fasta: list[Path], path_hmm_dir: Path):
        files = path_pep_fasta[:3]
        with pytest.raises(ValueError, match="Requires at least 4 samples"):
            _args_check(files, path_hmm_dir, 1e-10, "hmmalign")

    def test_raises_on_missing_markerset(self, path_pep_fasta: list[Path]):
        markerset = "nonexistent_markerset"
        with pytest.raises(FileNotFoundError, match=f"Markerset folder does not exist: {markerset}"):
            _args_check(path_pep_fasta, markerset, 1e-10, "hmmalign")

    def test_raises_on_evalue_gte_1(self, path_pep_fasta: list[Path]):
        with pytest.raises(ValueError, match="Invalid evalue"):
            _args_check(path_pep_fasta, "poxviridae_odb10", 1.0, "hmmalign")

    def test_raises_on_invalid_method(self, path_pep_fasta: list[Path]):
        with pytest.raises(ValueError, match="Invalid method"):
            _args_check(path_pep_fasta, "poxviridae_odb10", 1e-10, "invalid")


# ---------------------------------------------------------------------------
# _search_threads_check
# ---------------------------------------------------------------------------


class TestSearchThreadsCheck:
    @pytest.mark.parametrize(
        "n_samples, threads, expected_jobs, expected_tpj",
        [
            (4, 1, 1, 1),
            (4, 4, 1, 4),  # threads < 8: jobs=1, tpj=threads
            (1, 16, 1, 16),  # n_samples==1: jobs=1, tpj=threads
            (10, 16, 4, 4),  # threads >= 8, n_samples > 1: tpj=4, jobs=threads/4
            (2, 16, 2, 4),  # jobs would be 4 but capped to n_samples=2
        ],
    )
    def test_thread_distribution(self, n_samples, threads, expected_jobs, expected_tpj):
        jobs, tpj = _search_threads_check(n_samples, threads)
        assert tpj == expected_tpj
        assert jobs == expected_jobs

    def test_jobs_capped_by_n_samples(self):
        jobs, _ = _search_threads_check(n_samples=2, threads=32)
        assert jobs <= 2


# ---------------------------------------------------------------------------
# align
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_samplelist(monkeypatch: pytest.MonkeyPatch):
    samplelist = MagicMock(spec=SampleList)
    samplelist.seqtype = "pep"
    container_data = [MagicMock()]
    samplelist.__iter__.return_value = iter(container_data)
    samplelist.__getitem__.side_effect = lambda index: container_data[index]
    samplelist.__len__.return_value = len(container_data)

    monkeypatch.setattr(align_mod, "SampleList", MagicMock(return_value=samplelist))
    return samplelist


@pytest.fixture
def mock_align_precheck(mock_samplelist, monkeypatch: pytest.MonkeyPatch, tmp_path):
    manager = MagicMock(spec=SearchHitsManager)
    manager.orthologs = {"HMM1": "/mock/path/HMM1.fa"}
    manager.filter.return_value = manager

    precheck = MagicMock()
    precheck.output = tmp_path
    precheck.precheck.return_value = (mock_samplelist, manager)

    monkeypatch.setattr(align_mod, "AlignPrecheck", MagicMock(return_value=precheck))

    return {
        "AlignPrecheck": precheck,
        "remaining_samples": mock_samplelist,
        "manager": manager,
    }


@pytest.fixture
def mock_orthologs(monkeypatch: pytest.MonkeyPatch):
    msa_mock = MagicMock(spec=MultipleSeqAlignment)

    orthologs = MagicMock(spec=OrthologList)
    orthologs.names = ["HMM1"]
    orthologs.align.return_value = [msa_mock]

    monkeypatch.setattr(align_mod, "OrthologList", MagicMock(return_value=orthologs))
    return {
        "orthologs": orthologs,
        "msa": msa_mock,
    }


@pytest.fixture
def mock_align_env(mock_orthologs, mock_hmms_no_cutoff, monkeypatch: pytest.MonkeyPatch):
    mock_trim = MagicMock(return_value=mock_orthologs["msa"])
    mock_seqio_write = MagicMock()

    monkeypatch.setattr(align_mod, "HMMMarkerSet", MagicMock(return_value=mock_hmms_no_cutoff))
    monkeypatch.setattr(align_mod, "trim_gaps", mock_trim)
    monkeypatch.setattr("Bio.SeqIO.write", mock_seqio_write)

    return {
        "orthologs": mock_orthologs["orthologs"],
        "msa": mock_orthologs["msa"],
        "trim": mock_trim,
    }


class TestAlignDataFlow:
    """Test the align() function by monkeypatching all heavy external calls."""

    def test_align_full_data_flow(self, mock_align_precheck, mock_align_env, path_pep_fasta: list[Path], tmp_path):
        align(
            path_pep_fasta,
            tmp_path,
            markerset="poxviridae_odb10",
            seqtype="pep",
            method="hmmalign",
            non_trim=False,
            threads=1,
        )

        # Verify key data-flow calls
        mock_align_precheck["AlignPrecheck"].precheck.assert_called_once()
        mock_align_precheck["remaining_samples"].search.assert_called_once()
        mock_align_precheck["manager"].update.assert_called_once()
        mock_align_precheck["AlignPrecheck"].save_checkpoint.assert_called_once_with(mock_align_precheck["manager"])
        mock_align_precheck["manager"].filter.assert_called_once_with(min_taxa=4)
        mock_align_env["orthologs"].align.assert_called_once()
        mock_align_env["trim"].assert_called_once()
        mock_align_env["msa"].sort.assert_called_once()

    def test_align_no_trim(self, mock_align_precheck, mock_align_env, path_pep_fasta: list[Path], tmp_path):
        """When non_trim=True, trim_gaps should not be called."""
        align(
            path_pep_fasta,
            tmp_path,
            markerset="poxviridae_odb10",
            seqtype="pep",
            method="hmmalign",
            non_trim=True,
            threads=1,
        )
        mock_align_env["trim"].assert_not_called()

    def test_align_raises_when_all_orthologs_filtered(
        self, mock_align_precheck, mock_align_env, path_pep_fasta: list[Path], tmp_path
    ):
        """RuntimeError is raised when filter() removes all orthologs."""
        from phyling.exception import EmptyWarning

        manager_obj = mock_align_precheck["manager"]
        manager_obj.filter.side_effect = EmptyWarning("empty")

        with pytest.raises(RuntimeError, match="All orthologs were gone"):
            align(
                path_pep_fasta,
                tmp_path,
                markerset="poxviridae_odb10",
                seqtype="pep",
                method="hmmalign",
                non_trim=False,
                threads=1,
            )

    def test_align_skips_search_when_no_remaining(
        self, mock_align_precheck, mock_align_env, path_pep_fasta: list[Path], tmp_path
    ):
        """When precheck returns no remaining_samples, search() is never called."""

        precheck_obj = mock_align_precheck["AlignPrecheck"]
        manager_obj = mock_align_precheck["manager"]
        precheck_obj.precheck.return_value = (None, manager_obj)

        align(
            path_pep_fasta,
            tmp_path,
            markerset="poxviridae_odb10",
            seqtype="pep",
            method="hmmalign",
            non_trim=False,
            threads=1,
        )

        # search() should never have been called on anything
        manager_obj.update.assert_not_called()
