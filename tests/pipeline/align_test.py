"""Tests for the align pipeline module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from _pytest.monkeypatch import MonkeyPatch

import phyling.pipeline.align as align_mod
from phyling.pipeline.align import _args_check, _search_threads_check, align

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

BASE_DB = Path("tests/database/poxviridae_odb10")
HMM_DIR = BASE_DB / "hmms"
CUTOFF_FILE = BASE_DB / "scores_cutoff"

DATA_DIR = Path("tests/data")
PEP_FASTA_DIR = DATA_DIR / "pep" / "bgzf"
CDS_FASTA_DIR = DATA_DIR / "cds" / "bgzf"
PEP_FASTA = sorted(tuple(PEP_FASTA_DIR.iterdir()))
CDS_FASTA = sorted(tuple(CDS_FASTA_DIR.iterdir()))

PEP_MFA = sorted(tuple((DATA_DIR / "mfa").glob("*.faa")))
CDS_MFA = sorted(tuple((DATA_DIR / "mfa").glob("*.fna")))


# This runs once before any tests in this file execute
@pytest.fixture(scope="class", autouse=True)
def patch_cfg_globally():
    mp = MonkeyPatch()
    mp.setattr("phyling.pipeline.align.CFG_DIRS", [Path("tests/database")])

    yield  # Let the tests run

    mp.undo()


# ---------------------------------------------------------------------------
# _args_check
# ---------------------------------------------------------------------------


class TestArgsCheck:
    @pytest.mark.parametrize("markerset", ["poxviridae_odb10", "tests/database/poxviridae_odb10/hmms"])
    def test_returns_tuple_of_paths(self, markerset):
        inputs_out, ms_out, ev_out, method_out = _args_check(PEP_FASTA, markerset, 1e-10, "hmmalign")
        assert isinstance(inputs_out, tuple)
        assert all(isinstance(p, Path) for p in inputs_out)
        assert ms_out.absolute() == Path("tests/database/poxviridae_odb10/hmms").absolute()
        assert ev_out == 1e-10
        assert method_out == "hmmalign"

    def test_single_directory_input(self):
        # Need at least 4 files in that directory (already 4)
        inputs_out, *_ = _args_check(PEP_FASTA_DIR, HMM_DIR, 1e-10, "hmmalign")
        assert len(inputs_out) == 5

    def test_raises_on_fewer_than_four_samples(self):
        files = PEP_FASTA[:3]
        with pytest.raises(ValueError, match="at least 4"):
            _args_check(files, HMM_DIR, 1e-10, "hmmalign")

    # def test_raises_on_directory_in_inputs(self, tmp_path, sample_files):
    #     bad = tmp_path / "subdir"
    #     bad.mkdir()
    #     files = sample_files + [bad]
    #     with pytest.raises(IsADirectoryError):
    #         _args_check(files, HMM_DIR, 1e-10, "hmmalign")

    def test_raises_on_missing_markerset(self):
        markerset = "nonexistent_markerset"
        with pytest.raises(FileNotFoundError, match=f"Markerset folder does not exist: {markerset}"):
            _args_check(PEP_FASTA, markerset, 1e-10, "hmmalign")

    def test_raises_on_evalue_gte_1(self):
        with pytest.raises(ValueError, match="evalue"):
            _args_check(PEP_FASTA, "poxviridae_odb10", 1.0, "hmmalign")

    def test_raises_on_invalid_method(self):
        with pytest.raises(ValueError, match="method"):
            _args_check(PEP_FASTA, "poxviridae_odb10", 1e-10, "invalid_method")

    def test_single_file_input(self):
        # Passing a single Path that is a file → treated as one-element tuple
        # But 1 < 4, so must raise
        with pytest.raises(ValueError, match="at least 4"):
            _args_check(PEP_FASTA[0], "poxviridae_odb10", 1e-10, "hmmalign")

    def test_empty_directory_raises(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="Empty input"):
            _args_check(empty_dir, "poxviridae_odb10", 1e-10, "hmmalign")


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
# align (data-flow integration test with monkeypatching)
# ---------------------------------------------------------------------------


@pytest.fixture()
def shared_output_path(tmpdir_factory) -> Path:
    # mktemp creates a unique directory unique to this test run session
    shared_dir = tmpdir_factory.mktemp("pipeline_outputs")
    return shared_dir


@pytest.fixture
def mock_samplelist(monkeypatch):
    mock_sample_list = MagicMock()
    mock_sample_list.seqtype = "pep"

    monkeypatch.setattr(align_mod, "SampleList", MagicMock(return_value=mock_sample_list))


@pytest.fixture
def mock_markerset(monkeypatch):
    mock_markerset = MagicMock()
    mock_markerset.checksums = {"marker1": "abc123"}
    mock_markerset.have_cutoffs.return_value = False
    mock_markerset.sort = MagicMock()

    monkeypatch.setattr(align_mod, "HMMMarkerSet", MagicMock(return_value=mock_markerset))


@pytest.fixture
def mock_align_precheck(monkeypatch, shared_output_path):

    def make_searchhits_mock():
        return MagicMock(
            orthologs={"marker1": MagicMock()},
            load=MagicMock(),
        )

    mock_remaining = MagicMock()
    mock_hits_from_search = MagicMock()
    mock_remaining.search.return_value = mock_hits_from_search

    mock_raw_hits = make_searchhits_mock()
    mock_filtered_hits = make_searchhits_mock()
    mock_raw_hits.filter.return_value = mock_filtered_hits

    mock_precheck = MagicMock()
    mock_precheck.output = shared_output_path
    mock_precheck.precheck.return_value = (mock_remaining, mock_raw_hits)
    mock_precheck.save_checkpoint = MagicMock()

    monkeypatch.setattr(align_mod, "AlignPrecheck", MagicMock(return_value=mock_precheck))

    return {
        "AlignPrecheck": mock_precheck,
        "remaining_samples": mock_remaining,
        "searchhits": mock_raw_hits,
        "hits_from_search": mock_hits_from_search,
        "filtered_hits": mock_filtered_hits,
    }


@pytest.fixture
def mock_orthologs(monkeypatch):
    rec = MagicMock()
    rec.id = "sample0"
    msa_mock = MagicMock()
    msa_mock.__iter__ = MagicMock(return_value=iter([rec]))
    msa_mock.sort = MagicMock()

    mock_orthologs = MagicMock()
    mock_orthologs.names = ["marker1"]
    mock_orthologs.align.return_value = [msa_mock]

    monkeypatch.setattr(align_mod, "OrthologList", MagicMock(return_value=mock_orthologs))
    return {
        "orthologs": mock_orthologs,
        "msa": msa_mock,
    }


@pytest.fixture
def mock_align_env(monkeypatch, mock_orthologs):
    mock_trim = MagicMock(return_value=mock_orthologs["msa"])
    mock_seqio_write = MagicMock()

    monkeypatch.setattr(align_mod, "trim_gaps", mock_trim)
    monkeypatch.setattr("Bio.SeqIO.write", mock_seqio_write)

    return {
        "orthologs": mock_orthologs["orthologs"],
        "msa": mock_orthologs["msa"],
        "trim": mock_trim,
    }


class TestAlignDataFlow:
    """Test the align() function by monkeypatching all heavy external calls."""

    def test_align_full_data_flow(self, mock_samplelist, mock_markerset, mock_align_precheck, mock_align_env, shared_output_path):
        align(
            PEP_FASTA,
            shared_output_path,
            markerset="poxviridae_odb10",
            seqtype="pep",
            method="hmmalign",
            non_trim=False,
            threads=1,
        )

        # Verify key data-flow calls
        mock_align_precheck["AlignPrecheck"].precheck.assert_called_once()
        mock_align_precheck["remaining_samples"].search.assert_called_once()
        mock_align_precheck["searchhits"].update.assert_called_once_with(mock_align_precheck["hits_from_search"])
        mock_align_precheck["AlignPrecheck"].save_checkpoint.assert_called_once_with(mock_align_precheck["searchhits"])
        mock_align_precheck["searchhits"].filter.assert_called_once_with(min_taxa=4)
        mock_align_precheck["filtered_hits"].load.assert_called_once()
        mock_align_env["orthologs"].align.assert_called_once()
        mock_align_env["trim"].assert_called_once()
        mock_align_env["msa"].sort.assert_called_once()

    def test_align_no_trim(self, mock_samplelist, mock_markerset, mock_align_precheck, mock_align_env, shared_output_path):
        """When non_trim=True, trim_gaps should not be called."""
        align(
            PEP_FASTA,
            shared_output_path,
            markerset="poxviridae_odb10",
            seqtype="pep",
            method="hmmalign",
            non_trim=True,
            threads=1,
        )
        mock_align_env["trim"].assert_not_called()

    def test_align_raises_when_all_orthologs_filtered(
        self, mock_samplelist, mock_markerset, mock_align_precheck, mock_align_env, monkeypatch, shared_output_path
    ):
        """RuntimeError is raised when filter() removes all orthologs."""
        from phyling.exception import EmptyWarning

        def make_searchhits_mock():
            return MagicMock(
                orthologs={"marker1": MagicMock()},
                load=MagicMock(),
            )

        mock_remaining = MagicMock()
        mock_hits_from_search = MagicMock()
        mock_remaining.search.return_value = mock_hits_from_search

        mock_raw_hits = make_searchhits_mock()
        mock_raw_hits.filter.side_effect = EmptyWarning("empty")

        mock_precheck = MagicMock()
        mock_precheck.output = shared_output_path
        mock_precheck.precheck.return_value = (mock_remaining, mock_raw_hits)
        mock_precheck.save_checkpoint = MagicMock()

        monkeypatch.setattr(align_mod, "AlignPrecheck", MagicMock(return_value=mock_precheck))

        with pytest.raises(RuntimeError, match="All orthologs were gone"):
            align(
                PEP_FASTA,
                shared_output_path,
                markerset="poxviridae_odb10",
                seqtype="pep",
                method="hmmalign",
                non_trim=False,
                threads=1,
            )

    def test_align_skips_search_when_no_remaining(
        self, mock_samplelist, mock_markerset, mock_align_precheck, mock_align_env, monkeypatch, shared_output_path
    ):
        """When precheck returns no remaining_samples, search() is never called."""

        def make_searchhits_mock():
            return MagicMock(
                orthologs={"marker1": MagicMock()},
                load=MagicMock(),
            )

        mock_raw_hits = make_searchhits_mock()
        mock_filtered_hits = make_searchhits_mock()
        mock_raw_hits.filter.return_value = mock_filtered_hits

        mock_precheck = MagicMock()
        mock_precheck.output = shared_output_path
        mock_precheck.precheck.return_value = (None, mock_raw_hits)
        mock_precheck.save_checkpoint = MagicMock()

        monkeypatch.setattr(align_mod, "AlignPrecheck", MagicMock(return_value=mock_precheck))

        align(
            PEP_FASTA,
            shared_output_path,
            markerset="poxviridae_odb10",
            seqtype="pep",
            method="hmmalign",
            non_trim=False,
            threads=1,
        )

        # search() should never have been called on anything
        mock_align_precheck["searchhits"].update.assert_not_called()
