"""Tests for the tree pipeline module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import phyling.pipeline.tree as tree_mod
from phyling.lib import TreeOutputFiles
from phyling.lib.tree import MFA2Tree, MFA2TreeList
from phyling.pipeline.tree import _input_check, _validate_partition, tree


@pytest.fixture
def mock_mfa2treelist(monkeypatch: pytest.MonkeyPatch, tmp_path):
    mock_list = MagicMock(spec=MFA2TreeList)
    mock_list.seqtype = "pep"
    container_data = [MagicMock()]
    mock_list.__iter__.return_value = iter(container_data)
    mock_list.__getitem__.side_effect = lambda index: container_data[index]
    mock_list.__len__.return_value = len(container_data)

    monkeypatch.setattr(tree_mod, "MFA2TreeList", MagicMock(return_value=mock_list))

    fake_consensus_file = tmp_path / "consensus.nwk"
    fake_consensus_file.write_text("consensus_tree", encoding="utf-8")

    mock_list.get_consensus_tree.return_value = fake_consensus_file

    fake_concat_file = tmp_path / "concat.faa"
    fake_concat_file.write_text("concat_faa", encoding="utf-8")
    fake_partition_file = tmp_path / "partition.txt"
    fake_partition_file.write_text("partition", encoding="utf-8")

    mock_list.concat.return_value = (fake_concat_file, fake_partition_file)

    return mock_list


@pytest.fixture
def mock_concat_mfa2tree(monkeypatch: pytest.MonkeyPatch, tmp_path):
    fake_tree_file = tmp_path / "concat.nwk"
    fake_tree_file.write_text("concat_tree", encoding="utf-8")
    mfa2tree = MagicMock(spec=MFA2Tree)
    mfa2tree.file = Path(fake_tree_file)
    mfa2tree.build.return_value = Path(fake_tree_file)

    monkeypatch.setattr(tree_mod, "MFA2Tree", MagicMock(return_value=mfa2tree))

    return mfa2tree


# ---------------------------------------------------------------------------
# tree
# ---------------------------------------------------------------------------


class TestTree:
    """Test the tree() function by monkeypatching all heavy external calls."""

    def test_tree_raises_with_single_input(self, path_pep_msa: list[Path], tmp_path):
        """Raises ValueError when only 1 MSA file is given."""
        with pytest.raises(ValueError, match="Found only 1 MSA fasta"):
            tree(path_pep_msa[0], tmp_path)

    def test_tree_raises_with_empty_dir(self, tmp_path):
        """Raises FileNotFoundError when input directory is empty."""
        empty_dir = tmp_path / "empty"
        output_dir = tmp_path / "output"
        empty_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="Empty input directory"):
            tree(empty_dir, output_dir)

    def test_tree_raises_partition_without_concat(self, path_pep_msa: list[Path], tmp_path):
        """Raises ValueError when partition=True but concat=False."""
        with pytest.raises(ValueError, match="Partition is not allowed in consensus mode"):
            tree(path_pep_msa, tmp_path, method="raxml", partition=True, concat=False)

    def test_tree_raises_partition_with_ft(self, path_pep_msa: list[Path], tmp_path):
        """Raises ValueError when partition=True with FastTree method."""
        with pytest.raises(ValueError, match="Partition is not allowed with"):
            tree(path_pep_msa, tmp_path, method="ft", partition=True, concat=True)

    def test_tree_raises_inputs_from_different_folders(self, path_pep_msa: list[Path], tmp_path):
        """Raises RuntimeError when inputs come from different directories."""
        other_dir = tmp_path / "other"
        other_dir.mkdir()
        other_file = other_dir / "other.faa"
        other_file.touch()
        mixed_inputs = [path_pep_msa[0], other_file]
        with pytest.raises(RuntimeError, match="The inputs aren't in the same folder"):
            tree(mixed_inputs, tmp_path / "output")

    def test_tree_consensus_mode(self, mock_mfa2treelist, path_pep_msa: list[Path], tmp_path):
        """Test consensus mode (concat=False) data flow."""
        tree(path_pep_msa, tmp_path, method="ft", concat=False, bs=0, scfl=0)

        mock_mfa2treelist.build.assert_called_once()
        mock_mfa2treelist.get_consensus_tree.assert_called_once()

        output_tree = tmp_path / TreeOutputFiles.TREE_NW
        assert output_tree.exists()
        assert output_tree.read_text(encoding="utf-8") == "consensus_tree"

    def test_tree_concat_mode(self, mock_mfa2treelist, path_pep_msa: list[Path], mock_concat_mfa2tree, tmp_path):
        """Test concatenation mode (concat=True) data flow."""
        tree(path_pep_msa, tmp_path, method="ft", concat=True, partition=False, bs=0, scfl=0)

        mock_mfa2treelist.concat.assert_called_once()
        mock_concat_mfa2tree.build.assert_called_once()

        output_tree = tmp_path / TreeOutputFiles.TREE_NW
        assert output_tree.exists()
        assert output_tree.read_text(encoding="utf-8") == "concat_tree"

    def test_tree_concat_with_partition(
        self, mock_mfa2treelist, path_pep_msa: list[Path], mock_concat_mfa2tree, monkeypatch, tmp_path
    ):
        """Test concatenation mode with partition=True."""
        # Create a modelfinder result file
        modelfinder_result = tmp_path / "raxml" / "modelfinder" / "partition_result.txt"
        modelfinder_result.parent.mkdir(parents=True, exist_ok=True)
        modelfinder_result.write_text("model_result", encoding="utf-8")

        mock_modelfinder = MagicMock()
        mock_modelfinder.result = str(modelfinder_result)
        monkeypatch.setattr(tree_mod, "ModelFinder", MagicMock(return_value=mock_modelfinder))

        tree(path_pep_msa, tmp_path, method="raxml", concat=True, partition=True, bs=0, scfl=0)

        mock_modelfinder.run.assert_called_once()
        mock_concat_mfa2tree.build.assert_called_once()

        output_tree = tmp_path / TreeOutputFiles.TREE_NW
        assert output_tree.exists()
        assert output_tree.read_text(encoding="utf-8") == "concat_tree"

    def test_tree_figure_output(self, path_tree_file: Path, mock_mfa2treelist, path_pep_msa: list[Path], monkeypatch, tmp_path):
        """Test that figure output is created when figure=True."""
        mock_phylo = MagicMock()
        mock_phylo.read.side_effect = lambda path, *args: path.read_text(encoding="utf-8")

        monkeypatch.setattr(tree_mod, "Phylo", mock_phylo)

        tree(path_pep_msa, tmp_path, method="ft", concat=False, bs=0, scfl=0, figure=True)

        mock_phylo.draw.assert_called_once()
        output_fig = tmp_path / TreeOutputFiles.TREE_IMG
        assert output_fig.exists()

    def test_tree_consensus_bs_warning(self, mock_mfa2treelist, path_pep_msa: list[Path], caplog, tmp_path):
        """Test that a warning is issued when bs > 0 in consensus mode."""
        import logging

        with caplog.at_level(logging.WARNING):
            tree(path_pep_msa, tmp_path, method="ft", concat=False, bs=100, scfl=0)

        assert any(
            "Bootstrap and concordance factor calculation are disabled with consensus mode." in record.message
            for record in caplog.records
        )

    def test_validate_partition_invalid_type(self):
        """_validate_partition raises ValueError when partition is not bool."""
        with pytest.raises(ValueError, match="boolean value"):
            _validate_partition("yes", "RAXML", True)   # pyrefly: ignore [bad-argument-type]

    def test_validate_partition_ft_raises(self):
        """_validate_partition raises ValueError when method is FT with partition."""
        with pytest.raises(ValueError, match="Partition is not allowed with"):
            _validate_partition(True, "FT", True)

    @pytest.mark.parametrize("method", ["RAXML", "IQTREE"])
    def test_validate_partition_no_concat_raises(self, method):
        """_validate_partition raises ValueError when concat=False with partition."""
        with pytest.raises(ValueError, match="Partition is not allowed in consensus mode"):
            _validate_partition(True, method, False)

    @pytest.mark.parametrize("method", ["RAXML", "IQTREE"])
    def test_validate_partition_valid(self, method):
        """_validate_partition returns True when valid."""
        result = _validate_partition(True, method, True)
        assert result is True

    def test_input_check_single_file_raises(self, path_pep_msa: list[Path]):
        """_input_check raises ValueError when only 1 file is given."""
        with pytest.raises(ValueError, match="Found only 1 MSA fasta"):
            _input_check(path_pep_msa[0])

    def test_input_check_multiple_files(self, path_pep_msa: list[Path]):
        """_input_check returns tuple when multiple files are given."""
        result = _input_check(path_pep_msa)
        assert isinstance(result, tuple)
        assert len(result) == len(path_pep_msa)
