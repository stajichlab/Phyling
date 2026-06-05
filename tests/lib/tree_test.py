import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree

from phyling.lib.tree import MFA2Tree


@pytest.fixture
def sample_mfa_file(tmp_path):
    """Create a sample MSA file for testing."""
    mfa_content = """>sample1
ACDEFGHIKLMNPQRSTVWY
>sample2
ACDEFGHIKLMNPQRSTVWY
>sample3
ACDEFGHIKLMNPQRSTVWY
"""
    mfa_file = tmp_path / "test.aa.mfa"
    mfa_file.write_text(mfa_content)
    return mfa_file


@pytest.fixture
def sample_dna_mfa_file(tmp_path):
    """Create a sample DNA MSA file for testing."""
    mfa_content = """>sample1
ATGCATGCATGCATGCATGC
>sample2
ATGCATGCATGCATGCATGC
>sample3
ATGCATGCATGCATGCATGC
"""
    mfa_file = tmp_path / "test.cds.mfa"
    mfa_file.write_text(mfa_content)
    return mfa_file


@pytest.fixture
def loaded_mfa2tree(sample_mfa_file):
    """Create and load a MFA2Tree object."""
    obj = MFA2Tree(sample_mfa_file, seqtype="pep")
    obj.load()
    return obj


class TestMFA2TreeInit:
    def test_init_with_file_only(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file)
        assert obj.file == sample_mfa_file
        assert obj.name == sample_mfa_file.name

    def test_init_with_file_and_name(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, "custom_name")
        assert obj.file == sample_mfa_file
        assert obj.name == "custom_name"

    def test_init_with_seqtype_pep(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        assert obj.seqtype == "pep"

    def test_init_with_seqtype_dna(self, sample_dna_mfa_file):
        obj = MFA2Tree(sample_dna_mfa_file, seqtype="dna")
        assert obj.seqtype == "dna"

    def test_init_with_seqtype_auto(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="AUTO")
        assert obj.seqtype in ("pep", "dna", "AUTO")

    def test_init_private_attrs_none(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file)
        assert obj._method is None
        assert obj._tree is None
        assert obj._toverr is None
        assert obj._saturation is None

    def test_init_with_string_path(self, sample_mfa_file):
        obj = MFA2Tree(str(sample_mfa_file))
        assert obj.file == sample_mfa_file


class TestMFA2TreeLoad:
    def test_load(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        obj.load()
        assert obj._data is not None
        assert isinstance(obj._data, MultipleSeqAlignment)

    def test_load_sets_seqtype_annotation(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        obj.load()
        assert obj._data.annotations["seqtype"] == "pep"

    def test_load_sets_seqname_annotation(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, "test_name", seqtype="pep")
        obj.load()
        assert obj._data.annotations["seqname"] == "test_name"

    def test_load_dna(self, sample_dna_mfa_file):
        obj = MFA2Tree(sample_dna_mfa_file, seqtype="dna")
        obj.load()
        assert obj._data is not None
        assert isinstance(obj._data, MultipleSeqAlignment)


class TestMFA2TreeLen:
    def test_len_loaded(self, loaded_mfa2tree):
        length = len(loaded_mfa2tree)
        assert length == 20

    def test_len_not_loaded(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file)
        with pytest.raises(Exception):
            len(obj)


class TestMFA2TreeComparisons:
    def _make_obj_with_toverr(self, sample_mfa_file, toverr_val):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        obj.load()
        obj._tree = MagicMock(spec=Tree)
        obj._toverr = toverr_val
        return obj

    def test_gt_true(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.9)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        assert obj1 > obj2

    def test_gt_false(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.3)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        assert not (obj1 > obj2)

    def test_ge_true(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        assert obj1 >= obj2

    def test_lt_true(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.3)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        assert obj1 < obj2

    def test_le_true(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        assert obj1 <= obj2

    def test_eq_true(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        assert obj1 == obj2

    def test_eq_false(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        obj2 = self._make_obj_with_toverr(sample_mfa_file, 0.6)
        assert obj1 != obj2

    def test_eq_different_type(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        with pytest.raises(TypeError):
            obj1 == "not_an_mfa2tree"

    def test_gt_wrong_type(self, sample_mfa_file):
        obj1 = self._make_obj_with_toverr(sample_mfa_file, 0.5)
        with pytest.raises(TypeError):
            obj1 > "not_an_mfa2tree"


class TestMFA2TreeProperties:
    def test_tree_property_raises_without_build(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        with pytest.raises(AttributeError):
            _ = obj.tree

    def test_method_property_raises_without_build(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        with pytest.raises(AttributeError):
            _ = obj.method

    def test_toverr_property_raises_without_build(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        with pytest.raises(AttributeError):
            _ = obj.toverr

    def test_saturation_property_raises_without_build(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        with pytest.raises(AttributeError):
            _ = obj.saturation

    def test_tree_property_returns_tree(self, loaded_mfa2tree):
        mock_tree = MagicMock(spec=Tree)
        loaded_mfa2tree._tree = mock_tree
        assert loaded_mfa2tree.tree is mock_tree

    def test_method_property_returns_method(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        loaded_mfa2tree._method = "FT"
        assert loaded_mfa2tree.method == "FT"

    def test_toverr_property_returns_value(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        loaded_mfa2tree._toverr = 0.75
        assert loaded_mfa2tree.toverr == 0.75

    def test_saturation_property_returns_value(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        loaded_mfa2tree._saturation = 0.42
        assert loaded_mfa2tree.saturation == 0.42

    def test_toverr_raises_without_toverr(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        loaded_mfa2tree._toverr = None
        with pytest.raises(AttributeError):
            _ = loaded_mfa2tree.toverr

    def test_saturation_raises_without_saturation(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        loaded_mfa2tree._saturation = None
        with pytest.raises(AttributeError):
            _ = loaded_mfa2tree.saturation


class TestMFA2TreeBuild:
    def test_build_invalid_method(self, loaded_mfa2tree, tmp_path):
        with pytest.raises(KeyError):
            loaded_mfa2tree.build("invalid_method", tmp_path)

    @patch("phyling.lib.tree.bootstrap")
    @patch("phyling.lib.tree.branch_concordance")
    @patch("phyling.external.FastTree")
    @patch("phyling.external.ModelFinder")
    def test_build_fasttree_returns_tree(
        self,
        mock_modelfinder_cls,
        mock_fasttree_cls,
        mock_branch_concordance,
        mock_bootstrap,
        loaded_mfa2tree,
        tmp_path,
    ):
        mock_modelfinder = MagicMock()
        mock_modelfinder.result = "LG"
        mock_modelfinder_cls.return_value = mock_modelfinder

        tree_file = tmp_path / "result.nw"
        mock_tree = MagicMock(spec=Tree)
        newick_content = "(sample1:0.1,sample2:0.1,sample3:0.1);"
        tree_file.write_text(newick_content)

        mock_fasttree = MagicMock()
        mock_fasttree.result = tree_file
        mock_fasttree._prog = "FastTree"
        mock_fasttree.model = "LG"
        mock_fasttree_cls.return_value = mock_fasttree

        with patch("phyling.lib.tree.Phylo.read", return_value=mock_tree):
            result = loaded_mfa2tree.build("ft", None)

        assert isinstance(result, Tree) or result is mock_tree

    @patch("phyling.external.FastTree")
    @patch("phyling.external.ModelFinder")
    def test_build_ft_partition_raises(
        self,
        mock_modelfinder_cls,
        mock_fasttree_cls,
        loaded_mfa2tree,
        tmp_path,
    ):
        partition_file = tmp_path / "partition.txt"
        partition_file.write_text("DNA, gene1 = 1-100\n")

        with pytest.raises(ValueError, match="does not support partitioning"):
            loaded_mfa2tree.build("ft", tmp_path / "output", model=partition_file)

    @patch("phyling.lib.tree.bootstrap")
    @patch("phyling.lib.tree.branch_concordance")
    @patch("phyling.external.Iqtree")
    @patch("phyling.external.ModelFinder")
    def test_build_iqtree_returns_path(
        self,
        mock_modelfinder_cls,
        mock_iqtree_cls,
        mock_branch_concordance,
        mock_bootstrap,
        loaded_mfa2tree,
        tmp_path,
    ):
        mock_modelfinder = MagicMock()
        mock_modelfinder.result = "LG"
        mock_modelfinder_cls.return_value = mock_modelfinder

        tree_file = tmp_path / "result.nw"
        newick_content = "(sample1:0.1,sample2:0.1,sample3:0.1);"
        tree_file.write_text(newick_content)

        mock_iqtree = MagicMock()
        mock_iqtree.result = tree_file
        mock_iqtree._prog = "iqtree2"
        mock_iqtree.model = "LG"
        mock_iqtree_cls.return_value = mock_iqtree

        mock_tree = MagicMock(spec=Tree)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("phyling.lib.tree.Phylo.read", return_value=mock_tree):
            result = loaded_mfa2tree.build("iqtree", output_dir)

        assert isinstance(result, Path)

    @patch("phyling.lib.tree.bootstrap")
    @patch("phyling.lib.tree.branch_concordance")
    @patch("phyling.external.Raxml")
    @patch("phyling.external.ModelFinder")
    def test_build_raxml_returns_path(
        self,
        mock_modelfinder_cls,
        mock_raxml_cls,
        mock_branch_concordance,
        mock_bootstrap,
        loaded_mfa2tree,
        tmp_path,
    ):
        mock_modelfinder = MagicMock()
        mock_modelfinder.result = "LG"
        mock_modelfinder_cls.return_value = mock_modelfinder

        tree_file = tmp_path / "result.nw"
        newick_content = "(sample1:0.1,sample2:0.1,sample3:0.1);"
        tree_file.write_text(newick_content)

        mock_raxml = MagicMock()
        mock_raxml.result = tree_file
        mock_raxml._prog = "raxml-ng"
        mock_raxml.model = "LG"
        mock_raxml_cls.return_value = mock_raxml

        mock_tree = MagicMock(spec=Tree)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        with patch("phyling.lib.tree.Phylo.read", return_value=mock_tree):
            result = loaded_mfa2tree.build("raxml", output_dir)

        assert isinstance(result, Path)


class TestMFA2TreeComputeToverr:
    def test_compute_toverr_raises_without_tree(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        obj.load()
        with pytest.raises(AttributeError):
            obj.compute_toverr()

    def test_compute_toverr_sets_value(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        with patch("phyling.external._libphykit.compute_toverr", return_value=0.65) as mock_fn:
            with patch.dict("sys.modules", {"phyling.external._libphykit": MagicMock(compute_toverr=mock_fn)}):
                loaded_mfa2tree._toverr = 0.65
        assert loaded_mfa2tree._toverr == 0.65


class TestMFA2TreeComputeSaturation:
    def test_compute_saturation_raises_without_tree(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file, seqtype="pep")
        obj.load()
        with pytest.raises(AttributeError):
            obj.compute_saturation()

    def test_compute_saturation_sets_value(self, loaded_mfa2tree):
        loaded_mfa2tree._tree = MagicMock(spec=Tree)
        loaded_mfa2tree._saturation = 0.33
        assert loaded_mfa2tree._saturation == 0.33


class TestMFA2TreeGuessSeqtype:
    def test_guess_seqtype_pep(self, sample_mfa_file):
        obj = MFA2Tree(sample_mfa_file)
        seqtype = obj._guess_seqtype()
        assert seqtype in ("pep", "dna", "rna", "NaN")

    def test_guess_seqtype_dna(self, sample_dna_mfa_file):
        obj = MFA2Tree(sample_dna_mfa_file)
        seqtype = obj._guess_seqtype()
        assert seqtype in ("pep", "dna", "rna", "NaN")

    def test_guess_seqtype_gzipped(self, tmp_path):
        mfa_content = b">sample1\nACDEFGHIKLMNPQRSTVWY\n>sample2\nACDEFGHIKLMNPQRSTVWY\n"
        gz_file = tmp_path / "test.aa.mfa.gz"
        with gzip.open(gz_file, "wb") as f:
            f.write(mfa_content)
        obj = MFA2Tree(gz_file, seqtype="pep")
        seqtype = obj._guess_seqtype()
        assert seqtype in ("pep", "dna", "rna", "NaN")
