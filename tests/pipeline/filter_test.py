"""Tests for the filter pipeline module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import phyling.pipeline.filter as filter_mod
from phyling.lib.tree import TreeOutputFiles
from phyling.pipeline.filter import filter

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DATA_DIR = Path("tests/data")
MFA_DIR = DATA_DIR / "mfa"
PEP_MFA = sorted(tuple((MFA_DIR).glob("*.faa")))


# ---------------------------------------------------------------------------
# filter (data-flow integration test with monkeypatching)
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_treelist_gen(monkeypatch):
    """A factory fixture that preserves real file names and paths on mock elements."""

    def _generate_env(inputs, seqtype):
        # A helper function to build a mock tree with real name/file attributes
        def make_dynamic_mock_tree(file_input):
            path_obj = Path(file_input)

            mock_tree = MagicMock()
            mock_tree.name = path_obj.name
            mock_tree.file = path_obj
            mock_tree.toverr = 0.5

            mock_terminal = MagicMock()
            mock_terminal.name = path_obj.stem
            mock_tree.tree.get_terminals.return_value = [mock_terminal]
            return mock_tree

        mock_trees = [make_dynamic_mock_tree(inp) for inp in inputs]

        # 3. Create a top-level MagicMock to represent the MFA2TreeList container
        mock_list = MagicMock()
        mock_list.seqtype = seqtype

        # Configure the mock container to behave like a standard Python list
        mock_list.__len__ = MagicMock(return_value=len(mock_trees))
        mock_list.__iter__ = MagicMock(return_value=iter(mock_trees))

        # Handle index lookup and slicing dynamically
        mock_list.__getitem__ = MagicMock(
            side_effect=lambda idx: mock_trees[idx] if isinstance(idx, (int, slice)) else mock_trees
        )

        # 4. Patch the class constructor inside the production module
        monkeypatch.setattr(filter_mod, "MFA2TreeList", MagicMock(return_value=mock_list))

        return mock_list

    return _generate_env


@pytest.fixture
def mock_filter_precheck_gen(monkeypatch):
    def _generate_env(output, mfa2treelist, remains=1):

        remains = len(mfa2treelist) if remains > len(mfa2treelist) else remains

        mock_remained = MagicMock()
        mock_remained.__bool__.return_value = bool(remains)
        mock_remained.__iter__.side_effect = lambda: iter(mfa2treelist[:remains])
        mock_remained.__len__ = MagicMock(return_value=len(mfa2treelist[:remains]))

        mock_completed = MagicMock()
        mock_completed.__iter__.side_effect = lambda: iter(mfa2treelist)
        mock_completed.__len__ = MagicMock(return_value=len(mfa2treelist))
        mock_completed.__getitem__.side_effect = lambda idx: mfa2treelist[idx]

        mock_precheck = MagicMock()
        mock_precheck.output = output
        mock_precheck.precheck.return_value = (mock_remained, mock_completed)
        mock_precheck.save_checkpoint = MagicMock()

        monkeypatch.setattr(filter_mod, "FilterPrecheck", MagicMock(return_value=mock_precheck))

        return {
            "FilterPrecheck": mock_precheck,
            "remained": mock_remained,
            "completed": mock_completed,
        }

    return _generate_env


class TestFilter:
    """Test the filter() function by monkeypatching all heavy external calls."""

    @pytest.mark.parametrize("inputs, match", [[PEP_MFA[0], "Fewer than 3 inputs"], [PEP_MFA[0:2], "Fewer than 3 inputs"]])
    def test_filter_raises_invalid_inputs_number(self, inputs, match, tmp_path) -> None:
        with pytest.raises(ValueError, match=match):
            filter(
                inputs,
                tmp_path,
                top_n_toverr=2,
            )

    def test_filter_raises_when_top_n_equals_input_count(self, tmp_path):
        """SystemExit raised when top_n_toverr equals number of inputs."""
        with pytest.raises(SystemExit, match="Argument top_n_toverr is equal to the number of inputs. Do not need filtering"):
            filter(
                PEP_MFA,
                tmp_path,
                top_n_toverr=len(PEP_MFA),
            )

    @pytest.mark.parametrize(
        "inputs, top_n_toverr, match",
        [
            [PEP_MFA, 1, "should between 2"],
            [PEP_MFA[:3], 4, "can only be 2 since there are only 3 inputs"],
            [PEP_MFA, 10, "should between 2"],
        ],
    )
    def test_filter_raises_invalid_top_n(self, inputs, top_n_toverr, match, tmp_path):
        with pytest.raises(ValueError, match=match):
            filter(
                inputs,
                tmp_path,
                top_n_toverr=top_n_toverr,
            )

    def test_filter_full_data_flow(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """Test the full data flow of filter()."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck = mock_filter_precheck_gen(tmp_path, mock_treelist, 1)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
        )

        mock_filter_precheck["FilterPrecheck"].precheck.assert_called_once()
        mock_filter_precheck["remained"].build.assert_called_once()
        mock_filter_precheck["remained"].compute_toverr.assert_called_once()
        mock_filter_precheck["completed"].extend.assert_called_once_with(mock_filter_precheck["remained"])
        mock_filter_precheck["completed"].sort.assert_called_once()
        mock_filter_precheck["FilterPrecheck"].save_checkpoint.assert_called_once_with(mock_filter_precheck["completed"])

    def test_filter_skips_build_when_no_remained(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """When precheck returns no remained, build() and compute_toverr() are never called."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck = mock_filter_precheck_gen(tmp_path, mock_treelist, 0)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
        )

        mock_filter_precheck["remained"].build.assert_not_called()
        mock_filter_precheck["remained"].compute_toverr.assert_not_called()
        mock_filter_precheck["completed"].extend.assert_not_called()

    def test_filter_treeness_file_written(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """Verify treeness output file is created."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck_gen(tmp_path, mock_treelist)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
        )

        treeness_file = tmp_path / TreeOutputFiles.TREENESS
        assert treeness_file.exists()

    def test_filter_treeness_file_contains_selected_and_filtered(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """Verify treeness file content includes selected and filtered sections."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck_gen(tmp_path, mock_treelist)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
        )

        treeness_file = tmp_path / TreeOutputFiles.TREENESS
        content = treeness_file.read_text()
        assert "top 3" in content

    def test_filter_symlinks_created(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """Verify symlinks are created for selected MSA files."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck_gen(tmp_path, mock_treelist)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
        )

        for i in range(3):
            expected_symlink = tmp_path / PEP_MFA[i].name
            assert expected_symlink.exists() or expected_symlink.is_symlink()
        assert (tmp_path / PEP_MFA[-1].name).exists() is False

    def test_filter_build_called_with_noml_true_when_ml_false(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """When ml=False (default), build should be called with noml=True."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck = mock_filter_precheck_gen(tmp_path, mock_treelist)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
            ml=False,
            threads=2,
        )

        call_kwargs = mock_filter_precheck["remained"].build.call_args.kwargs
        assert call_kwargs.get("noml") is True

    def test_filter_with_ml_flag(self, mock_treelist_gen, mock_filter_precheck_gen, tmp_path):
        """Test that ml=True is passed correctly."""
        mock_treelist = mock_treelist_gen(PEP_MFA, "pep")
        mock_filter_precheck = mock_filter_precheck_gen(tmp_path, mock_treelist)

        filter(
            PEP_MFA,
            tmp_path,
            top_n_toverr=3,
            ml=True,
            threads=2,
        )

        call_kwargs = mock_filter_precheck["remained"].build.call_args.kwargs
        assert call_kwargs.get("noml") is False
