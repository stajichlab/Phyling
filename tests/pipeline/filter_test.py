"""Tests for the filter pipeline module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import phyling.pipeline.filter as filter_mod
from phyling.lib.tree import MFA2Tree, TreeOutputFiles
from phyling.pipeline.filter import filter

# ---------------------------------------------------------------------------
# filter
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_treelist(path_pep_msa: list[Path], monkeypatch: pytest.MonkeyPatch):
    def make_dynamic_mock_tree(file_input):
        path_obj = Path(file_input)

        mock_tree = MagicMock(spec=MFA2Tree)
        mock_tree.name = path_obj.name
        mock_tree.file = path_obj
        mock_tree.toverr = 0.5

        mock_terminal = MagicMock()
        mock_terminal.name = path_obj.stem
        mock_tree.tree.get_terminals.return_value = [mock_terminal]
        return mock_tree

    mock_trees = [make_dynamic_mock_tree(msa) for msa in path_pep_msa]

    mock_list = MagicMock()
    mock_list.seqtype = "pep"

    mock_list.__len__.return_value = len(mock_trees)
    mock_list.__iter__.return_value = iter(mock_trees)
    mock_list.__getitem__.side_effect = lambda idx: mock_trees[idx] if isinstance(idx, (int, slice)) else mock_trees

    monkeypatch.setattr(filter_mod, "MFA2TreeList", MagicMock(return_value=mock_list))

    return mock_list


@pytest.fixture
def mock_filter_precheck(mock_treelist, monkeypatch: pytest.MonkeyPatch, tmp_path):
    mock_precheck = MagicMock()
    mock_precheck.output = tmp_path
    mock_precheck.precheck.return_value = (mock_treelist, mock_treelist)

    monkeypatch.setattr(filter_mod, "FilterPrecheck", MagicMock(return_value=mock_precheck))

    return {
        "FilterPrecheck": mock_precheck,
        "mfa2treelist": mock_treelist,
    }


class TestFilter:
    """Test the filter() function by monkeypatching all heavy external calls."""

    def test_filter_raises_invalid_inputs_number(self, path_pep_msa: list[Path], tmp_path) -> None:
        with pytest.raises(ValueError, match="Fewer than 3 inputs"):
            filter(path_pep_msa[0], tmp_path, top_n_toverr=2)
        with pytest.raises(ValueError, match="Fewer than 3 inputs"):
            filter(path_pep_msa[0:2], tmp_path, top_n_toverr=2)

    def test_filter_raises_when_top_n_equals_input_count(self, path_pep_msa: list[Path], tmp_path):
        """SystemExit raised when top_n_toverr equals number of inputs."""
        with pytest.raises(SystemExit, match="Argument top_n_toverr is equal to the number of inputs. Do not need filtering"):
            filter(path_pep_msa, tmp_path, top_n_toverr=len(path_pep_msa))

    def test_filter_raises_invalid_top_n(self, path_pep_msa: list[Path], tmp_path):
        with pytest.raises(ValueError, match="should between 2"):
            filter(path_pep_msa, tmp_path, top_n_toverr=1)
        with pytest.raises(ValueError, match="should between 2"):
            filter(path_pep_msa, tmp_path, top_n_toverr=10)
        with pytest.raises(ValueError, match="can only be 2 since there are only 3 inputs"):
            filter(path_pep_msa[:3], tmp_path, top_n_toverr=4)

    def test_filter_full_data_flow(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """Test the full data flow of filter()."""
        filter(path_pep_msa, tmp_path, top_n_toverr=3)

        mock_filter_precheck["FilterPrecheck"].precheck.assert_called_once()
        mock_filter_precheck["mfa2treelist"].build.assert_called_once()
        mock_filter_precheck["mfa2treelist"].compute_toverr.assert_called_once()
        mock_filter_precheck["mfa2treelist"].extend.assert_called_once()
        mock_filter_precheck["mfa2treelist"].sort.assert_called_once()
        mock_filter_precheck["FilterPrecheck"].save_checkpoint.assert_called_once_with(mock_filter_precheck["mfa2treelist"])

    def test_filter_skips_build_when_no_remained(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """When precheck returns no remained, build() and compute_toverr() are never called."""
        precheck_obj = mock_filter_precheck["FilterPrecheck"]
        mock_treelist = mock_filter_precheck["mfa2treelist"]
        precheck_obj.precheck.return_value = (None, mock_treelist)
        filter(path_pep_msa, tmp_path, top_n_toverr=3)

        mock_treelist.build.assert_not_called()
        mock_treelist.compute_toverr.assert_not_called()
        mock_treelist.extend.assert_not_called()

    def test_filter_treeness_file_written(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """Verify treeness output file is created."""
        filter(path_pep_msa, tmp_path, top_n_toverr=3)

        treeness_file = tmp_path / TreeOutputFiles.TREENESS
        assert treeness_file.exists()

    def test_filter_treeness_file_contains_selected_and_filtered(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """Verify treeness file content includes selected and filtered sections."""
        filter(path_pep_msa, tmp_path, top_n_toverr=3)

        treeness_file = tmp_path / TreeOutputFiles.TREENESS
        content = treeness_file.read_text()
        assert "top 3" in content

    def test_filter_symlinks_created(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """Verify symlinks are created for selected MSA files."""
        filter(path_pep_msa, tmp_path, top_n_toverr=3)

        for i in range(3):
            expected_symlink = tmp_path / path_pep_msa[i].name
            assert expected_symlink.exists() or expected_symlink.is_symlink()
        assert (tmp_path / path_pep_msa[-1].name).exists() is False

    def test_filter_build_called_with_noml_true_when_ml_false(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """When ml=False (default), build should be called with noml=True."""
        filter(path_pep_msa, tmp_path, top_n_toverr=3, ml=False, threads=2)

        call_kwargs = mock_filter_precheck["mfa2treelist"].build.call_args.kwargs
        assert call_kwargs.get("noml") is True

    def test_filter_with_ml_flag(self, mock_filter_precheck, path_pep_msa: list[Path], tmp_path):
        """Test that ml=True is passed correctly."""
        filter(path_pep_msa, tmp_path, top_n_toverr=3, ml=True, threads=2)

        call_kwargs = mock_filter_precheck["mfa2treelist"].build.call_args.kwargs
        assert call_kwargs.get("noml") is False
