"""Tests for the tree pipeline module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

import phyling.pipeline.tree as tree_mod
from phyling.lib import TreeOutputFiles
from phyling.pipeline.tree import tree

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DATA_DIR = Path("tests/data")
MFA_DIR = DATA_DIR / "mfa"
PEP_MFA = sorted(tuple((MFA_DIR).glob("*.faa")))


# ---------------------------------------------------------------------------
# tree (data-flow integration test with monkeypatching)
# ---------------------------------------------------------------------------


@pytest.fixture()
def shared_output_path(tmpdir_factory) -> Path:
    shared_dir = tmpdir_factory.mktemp("tree_outputs")
    return shared_dir


@pytest.fixture
def mock_mfa2treelist_gen(monkeypatch):
    """A factory fixture that returns a mock MFA2TreeList."""

    def _generate_env(inputs, seqtype):
        mock_list = MagicMock()
        mock_list.seqtype = seqtype
        mock_list.__len__ = MagicMock(return_value=len(inputs))
        mock_list.__iter__ = MagicMock(return_value=iter(inputs))

        monkeypatch.setattr(tree_mod, "MFA2TreeList", MagicMock(return_value=mock_list))

        return mock_list

    return _generate_env


@pytest.fixture
def mock_mfa2tree_gen(monkeypatch):
    """A factory fixture that returns a mock MFA2Tree."""

    def _generate_env(tree_file):
        mock_tree = MagicMock()
        mock_tree.file = Path(tree_file)
        mock_tree._data = ["seq"] * 10
        mock_tree.__len__ = MagicMock(return_value=5000)
        mock_tree.build.return_value = Path(tree_file)

        monkeypatch.setattr(tree_mod, "MFA2Tree", MagicMock(return_value=mock_tree))

        return mock_tree

    return _generate_env


class TestTree:
    """Test the tree() function by monkeypatching all heavy external calls."""

    def test_tree_raises_with_single_input(self, tmp_path):
        """Raises ValueError when only 1 MSA file is given."""
        with pytest.raises(ValueError, match="Found only 1 MSA fasta"):
            tree(
                PEP_MFA[0],
                tmp_path,
            )

    def test_tree_raises_with_empty_dir(self, tmp_path):
        """Raises FileNotFoundError when input directory is empty."""
        empty_dir = tmp_path / "empty"
        output_dir = tmp_path / "output"
        empty_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="Empty input directory"):
            tree(
                empty_dir,
                output_dir,
            )

    def test_tree_raises_partition_without_concat(self, tmp_path):
        """Raises ValueError when partition=True but concat=False."""
        with pytest.raises(ValueError, match="Partition is not allowed in consensus mode"):
            tree(
                PEP_MFA,
                tmp_path,
                method="raxml",
                partition=True,
                concat=False,
            )

    def test_tree_raises_partition_with_ft(self, tmp_path):
        """Raises ValueError when partition=True with FastTree method."""
        with pytest.raises(ValueError, match="Partition is not allowed with"):
            tree(
                PEP_MFA,
                tmp_path,
                method="ft",
                partition=True,
                concat=True,
            )

    def test_tree_raises_inputs_from_different_folders(self, tmp_path):
        """Raises RuntimeError when inputs come from different directories."""
        other_dir = tmp_path / "other"
        other_dir.mkdir()
        other_file = other_dir / "other.faa"
        other_file.touch()
        mixed_inputs = [PEP_MFA[0], other_file]
        with pytest.raises(RuntimeError, match="The inputs aren't in the same folder"):
            tree(
                mixed_inputs,
                tmp_path / "output",
            )

    def test_tree_consensus_mode(self, mock_mfa2treelist_gen, tmp_path):
        """Test consensus mode (concat=False) data flow."""
        fake_tree_file = tmp_path / "consensus.nwk"
        fake_tree_file.write_text("consensus_tree", encoding="utf-8")

        mock_list = mock_mfa2treelist_gen(PEP_MFA, "pep")
        mock_list.get_consensus_tree.return_value = fake_tree_file

        tree(
            PEP_MFA,
            tmp_path,
            method="ft",
            concat=False,
            bs=0,
            scfl=0,
        )

        mock_list.build.assert_called_once()
        mock_list.get_consensus_tree.assert_called_once()

        output_tree = tmp_path / TreeOutputFiles.TREE_NW
        assert output_tree.exists()
        assert output_tree.read_text(encoding="utf-8") == "consensus_tree"

    def test_tree_concat_mode(self, mock_mfa2treelist_gen, mock_mfa2tree_gen, tmp_path):
        """Test concatenation mode (concat=True) data flow."""
        fake_concat_file = tmp_path / "concat.faa"
        fake_concat_file.write_text("concat_faa", encoding="utf-8")
        fake_partition_file = tmp_path / "partition.txt"
        fake_partition_file.write_text("partition", encoding="utf-8")
        fake_tree_file = tmp_path / "concat.nwk"
        fake_tree_file.write_text("concat_tree", encoding="utf-8")

        mock_list = mock_mfa2treelist_gen(PEP_MFA, "pep")
        mock_list.concat.return_value = (fake_concat_file, fake_partition_file)

        mock_ct = mock_mfa2tree_gen(fake_tree_file)
        mock_ct.build.return_value = fake_tree_file

        tree(
            PEP_MFA,
            tmp_path,
            method="ft",
            concat=True,
            partition=False,
            bs=0,
            scfl=0,
        )

        mock_list.concat.assert_called_once()
        mock_ct.build.assert_called_once()

        output_tree = tmp_path / TreeOutputFiles.TREE_NW
        assert output_tree.exists()
        assert output_tree.read_text(encoding="utf-8") == "concat_tree"

    def test_tree_concat_with_partition(self, mock_mfa2treelist_gen, mock_mfa2tree_gen, monkeypatch, tmp_path):
        """Test concatenation mode with partition=True."""
        fake_concat_file = tmp_path / "concat.faa"
        fake_concat_file.write_text("concat_faa", encoding="utf-8")
        fake_partition_file = tmp_path / "partition.txt"
        fake_partition_file.write_text("partition", encoding="utf-8")
        fake_tree_file = tmp_path / "concat.nwk"
        fake_tree_file.write_text("concat_tree", encoding="utf-8")

        mock_list = mock_mfa2treelist_gen(PEP_MFA, "pep")
        mock_list.concat.return_value = (fake_concat_file, fake_partition_file)

        mock_ct = mock_mfa2tree_gen(fake_tree_file)
        mock_ct.build.return_value = fake_tree_file

        # Create a modelfinder result file
        modelfinder_result = tmp_path / "raxml" / "modelfinder" / "partition_result.txt"
        modelfinder_result.parent.mkdir(parents=True, exist_ok=True)
        modelfinder_result.write_text("model_result", encoding="utf-8")

        mock_modelfinder = MagicMock()
        mock_modelfinder.result = str(modelfinder_result)
        monkeypatch.setattr(tree_mod, "ModelFinder", MagicMock(return_value=mock_modelfinder))

        tree(
            PEP_MFA,
            tmp_path,
            method="raxml",
            concat=True,
            partition=True,
            bs=0,
            scfl=0,
        )

        mock_modelfinder.run.assert_called_once()
        mock_ct.build.assert_called_once()

        output_tree = tmp_path / TreeOutputFiles.TREE_NW
        assert output_tree.exists()
        assert output_tree.read_text(encoding="utf-8") == "concat_tree"

    def test_tree_figure_output(self, mock_mfa2treelist_gen, monkeypatch, tmp_path):
        """Test that figure output is created when figure=True."""
        fake_tree_file = tmp_path / "tree.nwk"
        fake_tree_file.write_text("tree", encoding="utf-8")

        mock_list = mock_mfa2treelist_gen(PEP_MFA, "pep")
        mock_list.get_consensus_tree.return_value = fake_tree_file

        mock_phylo = MagicMock()
        mock_phylo.read.side_effect = lambda path, *args: path.read_text(encoding="utf-8")

        monkeypatch.setattr("phyling.pipeline.tree.Phylo", mock_phylo)

        tree(
            PEP_MFA,
            tmp_path,
            method="ft",
            concat=False,
            bs=0,
            scfl=0,
            figure=True,
        )

        output_fig = tmp_path / TreeOutputFiles.TREE_IMG
        assert output_fig.exists()

    def test_tree_consensus_bs_warning(self, mock_mfa2treelist_gen, caplog, tmp_path):
        """Test that a warning is issued when bs > 0 in consensus mode."""
        import logging

        fake_tree_file = tmp_path / "consensus.nwk"
        fake_tree_file.write_text("(A,B);")

        mock_list = mock_mfa2treelist_gen(PEP_MFA, "pep")
        mock_list.get_consensus_tree.return_value = fake_tree_file

        with caplog.at_level(logging.WARNING):
            tree(
                PEP_MFA,
                tmp_path,
                method="ft",
                concat=False,
                bs=100,
                scfl=0,
            )

        assert any(
            "Bootstrap and concordance factor calculation are disabled with consensus mode." in record.message
            for record in caplog.records
        )

    def test_validate_partition_invalid_type(self, tmp_path):
        """_validate_partition raises ValueError when partition is not bool."""
        from phyling.pipeline.tree import _validate_partition

        with pytest.raises(ValueError, match="boolean value"):
            _validate_partition("yes", "RAXML", True)

    def test_validate_partition_ft_raises(self, tmp_path):
        """_validate_partition raises ValueError when method is FT with partition."""
        from phyling.pipeline.tree import _validate_partition

        with pytest.raises(ValueError, match="Partition is not allowed with"):
            _validate_partition(True, "FT", True)

    def test_validate_partition_no_concat_raises(self, tmp_path):
        """_validate_partition raises ValueError when concat=False with partition."""
        from phyling.pipeline.tree import _validate_partition

        with pytest.raises(ValueError, match="Partition is not allowed in consensus mode"):
            _validate_partition(True, "RAXML", False)

    def test_validate_partition_valid(self, tmp_path):
        """_validate_partition returns True when valid."""
        from phyling.pipeline.tree import _validate_partition

        result = _validate_partition(True, "RAXML", True)
        assert result is True

    def test_input_check_single_file_raises(self, tmp_path):
        """_input_check raises ValueError when only 1 file is given."""
        from phyling.pipeline.tree import _input_check

        with pytest.raises(ValueError, match="Found only 1 MSA fasta"):
            _input_check(PEP_MFA[0])

    def test_input_check_multiple_files(self, tmp_path):
        """_input_check returns tuple when multiple files are given."""
        from phyling.pipeline.tree import _input_check

        result = _input_check(PEP_MFA)
        assert isinstance(result, tuple)
        assert len(result) == len(PEP_MFA)
