"""Tests for the tree module library."""

from __future__ import annotations

import random
import re
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from _pytest.monkeypatch import MonkeyPatch
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyling.lib.tree import MFA2Tree, MFA2TreeList, bootstrap, branch_concordance, fill_missing_taxon

DATA_DIR = Path("tests/data")
PEP_MSA = tuple((DATA_DIR / "msa").glob("*.faa"))
CDS_MSA = tuple((DATA_DIR / "msa").glob("*.fna"))


@pytest.fixture
def mfa2tree_obj():
    """Create and load a MFA2Tree object."""
    obj = MFA2Tree(PEP_MSA[0], seqtype="pep")
    return obj


def make_tree_with_toverr(msa, toverr_val):
    obj = MFA2Tree(msa)
    obj.load()
    obj._tree = MagicMock(spec=Tree)
    obj._toverr = toverr_val
    return obj


class TestMFA2Tree:
    pep_msa = PEP_MSA[0]
    cds_msa = CDS_MSA[0]

    def test_init_with_file_only(self):
        obj = MFA2Tree(self.pep_msa)
        assert obj.file == self.pep_msa.absolute()
        assert obj.name == self.pep_msa.name

    def test_init_with_file_and_name(self):
        obj = MFA2Tree(self.pep_msa, "custom_name")
        assert obj.file == self.pep_msa.absolute()
        assert obj.name == "custom_name"

    def test_init_with_seqtype_pep(self):
        obj = MFA2Tree(self.pep_msa, seqtype="pep")
        assert obj.seqtype == "pep"

    def test_init_with_seqtype_dna(self):
        obj = MFA2Tree(self.cds_msa, seqtype="dna")
        assert obj.seqtype == "dna"

    @pytest.mark.parametrize("input, expect", [["pep_msa", "pep"], ["cds_msa", "dna"]])
    def test_init_with_seqtype_auto(self, input, expect):
        obj = MFA2Tree(getattr(self, input), seqtype="AUTO")
        assert obj.seqtype == expect

    def test_init_private_attrs_none(self):
        obj = MFA2Tree(self.pep_msa)
        assert obj._method is None
        assert obj._tree is None
        assert obj._toverr is None
        assert obj._saturation is None

    def test_init_with_string_path(self):
        obj = MFA2Tree(str(self.pep_msa))
        assert obj.file == self.pep_msa.absolute()

    def test_load(self, mfa2tree_obj):
        mfa2tree_obj.load()
        assert mfa2tree_obj._data is not None
        assert isinstance(mfa2tree_obj._data, MultipleSeqAlignment)

    def test_load_dna(self):
        obj = MFA2Tree(self.cds_msa, seqtype="dna")
        obj.load()
        assert obj._data is not None
        assert isinstance(obj._data, MultipleSeqAlignment)

    def test_load_sets_seqtype_annotation(self, mfa2tree_obj):
        mfa2tree_obj.load()
        assert mfa2tree_obj._data.annotations["seqtype"] == "pep"

    def test_load_sets_seqname_annotation(self, mfa2tree_obj):
        mfa2tree_obj.load()
        assert mfa2tree_obj._data.annotations["seqname"] == PEP_MSA[0].name

    def test_len_not_loaded(self, mfa2tree_obj):
        with pytest.raises(Exception):
            len(mfa2tree_obj)
        mfa2tree_obj.load()
        length = len(mfa2tree_obj)
        assert length == 1345

    def test_gt_true(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.9)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.5)
        assert obj1 > obj2

    def test_gt_false(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.3)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.5)
        assert not (obj1 > obj2)

    def test_ge_true(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.5)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.5)
        assert obj1 >= obj2

    def test_lt_true(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.3)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.5)
        assert obj1 < obj2

    def test_le_true(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.5)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.5)
        assert obj1 <= obj2

    def test_eq_true(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.5)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.5)
        assert obj1 == obj2

    def test_eq_false(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.5)
        obj2 = make_tree_with_toverr(self.pep_msa, 0.6)
        assert obj1 != obj2

    def test_eq_different_type(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.5)
        with pytest.raises(TypeError):
            obj1 == "not_an_mfa2tree"

    def test_gt_wrong_type(self):
        obj1 = make_tree_with_toverr(self.pep_msa, 0.5)
        with pytest.raises(TypeError):
            obj1 > "not_an_mfa2tree"  # type: ignore

    def test_tree_property_raises_without_build(self, mfa2tree_obj):
        with pytest.raises(AttributeError):
            _ = mfa2tree_obj.tree

    def test_method_property_raises_without_build(self, mfa2tree_obj):
        with pytest.raises(AttributeError):
            _ = mfa2tree_obj.method

    def test_toverr_property_raises_without_build(self, mfa2tree_obj):
        with pytest.raises(AttributeError):
            _ = mfa2tree_obj.toverr

    def test_saturation_property_raises_without_build(self, mfa2tree_obj):
        with pytest.raises(AttributeError):
            _ = mfa2tree_obj.saturation

    def test_tree_property_returns_tree(self, mfa2tree_obj):
        mock_tree = MagicMock(spec=Tree)
        mfa2tree_obj._tree = mock_tree
        assert mfa2tree_obj.tree is mock_tree

    def test_method_property_returns_method(self, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mfa2tree_obj._method = "FT"
        assert mfa2tree_obj.method == "FT"

    def test_toverr_property_returns_value(self, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mfa2tree_obj._toverr = 0.75
        assert mfa2tree_obj.toverr == 0.75

    def test_saturation_property_returns_value(self, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mfa2tree_obj._saturation = 0.42
        assert mfa2tree_obj.saturation == 0.42

    def test_toverr_raises_without_toverr(self, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mfa2tree_obj._toverr = None
        with pytest.raises(AttributeError):
            _ = mfa2tree_obj.toverr

    def test_saturation_raises_without_saturation(self, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mfa2tree_obj._saturation = None
        with pytest.raises(AttributeError):
            _ = mfa2tree_obj.saturation


@pytest.fixture(scope="class")
def class_monkeypatch():
    """A class-scoped monkeypatch fixture."""
    mp = MonkeyPatch()
    yield mp
    mp.undo()


@pytest.fixture(scope="class")
def mock_build_env(class_monkeypatch, tmpdir_factory):
    """Mocks all external runners, utilities, and file readers inside tree_builder."""
    shared_dir = tmpdir_factory.mktemp("tree_data")

    # 2. Create the MagicMocks
    mock_modelfinder = MagicMock()
    mock_modelfinder.return_value.result = "model"

    mock_fasttree = MagicMock()
    ft_tree = shared_dir / "ft_tree.nw"
    ft_tree.write_text("ft_tree", encoding="utf-8")
    mock_fasttree.return_value.result = ft_tree

    mock_raxml = MagicMock()
    raxml_tree = shared_dir / "raxml_tree.nw"
    raxml_tree.write_text("raxml_tree", encoding="utf-8")
    mock_raxml.return_value.result = raxml_tree

    mock_iqtree = MagicMock()
    iqtree_tree = shared_dir / "iqtree_tree.nw"
    iqtree_tree.write_text("iqtree_tree", encoding="utf-8")
    mock_iqtree.return_value.result = iqtree_tree

    mock_bootstrap = MagicMock()
    bootstrap_file = shared_dir / "ufboot.nw"
    bootstrap_file.write_text("ufboot", encoding="utf-8")
    mock_bootstrap.return_value = bootstrap_file

    mock_concordance = MagicMock()
    mock_phylo = MagicMock()
    mock_phylo.read.side_effect = lambda path, format: MagicMock(
        spec=Tree,
        _mock_file_content=path.read_text(encoding="utf-8"),  # Optional: store the text inside the mock if needed
    )
    concordance_file = shared_dir / "concord.nw"
    concordance_file.write_text("concord", encoding="utf-8")
    mock_concordance.return_value = concordance_file

    class_monkeypatch.setattr("phyling.external.ModelFinder", mock_modelfinder)
    class_monkeypatch.setattr("phyling.external.FastTree", mock_fasttree)
    class_monkeypatch.setattr("phyling.external.Raxml", mock_raxml)
    class_monkeypatch.setattr("phyling.external.Iqtree", mock_iqtree)
    class_monkeypatch.setattr("phyling.lib.tree.bootstrap", mock_bootstrap)
    class_monkeypatch.setattr("phyling.lib.tree.branch_concordance", mock_concordance)
    class_monkeypatch.setattr("phyling.lib.tree.Phylo", mock_phylo)


class TestMFA2TreeBuild:
    pep_msa = PEP_MSA[0]
    cds_msa = CDS_MSA[0]

    def test_build_invalid_method(self, mock_build_env, mfa2tree_obj, tmp_path):
        with pytest.raises(KeyError):
            mfa2tree_obj.build("invalid_method", tmp_path)

    @pytest.mark.parametrize("method", ["ft", "raxml", "iqtree"])
    def test_build_returns_tree(self, mock_build_env, mfa2tree_obj, method):
        result = mfa2tree_obj.build(method, None)
        print(type(result))
        assert isinstance(result, Tree)
        assert result._mock_file_content == f"{method}_tree"  # type: ignore

    @pytest.mark.parametrize("method", ["ft", "raxml", "iqtree"])
    def test_build_returns_path(self, mock_build_env, mfa2tree_obj, method, tmp_path):
        result = mfa2tree_obj.build(method, tmp_path)
        assert isinstance(result, Path)
        assert result.name == f"{method}_tree.nw"

    def test_build_ft_partition_raises(self, mock_build_env, mfa2tree_obj, tmp_path):
        partition_file = tmp_path / "partition.txt"
        partition_file.touch()

        with pytest.raises(ValueError, match="does not support partitioning"):
            mfa2tree_obj.build("ft", tmp_path / "output", model=partition_file)

    def test_build_with_bs(self, mock_build_env, mfa2tree_obj, tmp_path):
        result = mfa2tree_obj.build("ft", tmp_path, bs=1000)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "ufboot"

    def test_build_with_scfl(self, mock_build_env, mfa2tree_obj, tmp_path):
        result = mfa2tree_obj.build("ft", tmp_path, scfl=100)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "concord"


class TestMFA2TreeComputeToverr:
    def test_compute_toverr_raises_without_tree(self, mfa2tree_obj):
        mfa2tree_obj.load()
        with pytest.raises(AttributeError):
            mfa2tree_obj.compute_toverr()

    def test_compute_toverr_sets_value(self, monkeypatch, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mock_compute_toverr = MagicMock()
        val = random.random()
        mock_compute_toverr.return_value = val
        monkeypatch.setattr("phyling.external._libphykit.compute_toverr", mock_compute_toverr)
        mfa2tree_obj.compute_toverr()
        assert mfa2tree_obj.toverr == val


class TestMFA2TreeComputeSaturation:
    def test_compute_saturation_raises_without_tree(self, mfa2tree_obj):
        mfa2tree_obj.load()
        with pytest.raises(AttributeError):
            mfa2tree_obj.compute_saturation()

    def test_compute_saturation_sets_value(self, monkeypatch, mfa2tree_obj):
        mfa2tree_obj._tree = MagicMock(spec=Tree)
        mock_saturation_class = MagicMock()
        val = random.random()
        mock_saturation_instance = mock_saturation_class.return_value
        mock_saturation_instance.compute_saturation.return_value = val
        monkeypatch.setattr("phyling.external._libphykit.Saturation", mock_saturation_class)
        mfa2tree_obj.compute_saturation()
        assert mfa2tree_obj.saturation == val


@pytest.fixture
def mfa2treelist_obj():
    return MFA2TreeList(PEP_MSA, seqtype="pep")


class TestMFA2TreeList:
    def test_init_empty(self):
        obj = MFA2TreeList()
        assert len(obj) == 0

    @pytest.mark.parametrize("msa", [PEP_MSA, CDS_MSA])
    def test_init_with_files(self, msa):
        obj = MFA2TreeList(msa)
        assert len(obj) == 2

    def test_init_with_names(self):
        names = ["name1", "name2"]
        obj = MFA2TreeList(PEP_MSA, names, seqtype="pep")
        assert len(obj) == 2
        for mfa2tree in obj:
            assert mfa2tree.name in names

    def test_getitem_by_index(self, mfa2treelist_obj):
        item = mfa2treelist_obj[0]
        assert isinstance(item, MFA2Tree)

    def test_getitem_by_slice(self, mfa2treelist_obj):
        subset = mfa2treelist_obj[0:1]
        assert isinstance(subset, MFA2TreeList)
        assert len(subset) == 1

    def test_getitem_by_name(self):
        names = ["name1", "name2"]
        obj = MFA2TreeList(PEP_MSA, names, seqtype="pep")
        item = obj["name1"]
        assert isinstance(item, MFA2Tree)
        assert item.file == PEP_MSA[0].absolute()

    def test_trees_property_raises_without_build(self, mfa2treelist_obj):
        with pytest.raises(AttributeError):
            _ = mfa2treelist_obj.trees

    def test_trees_property_returns_list(self, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
        trees = mfa2treelist_obj.trees
        assert isinstance(trees, list)
        assert len(trees) == 2
        assert isinstance(trees[0], Tree)

    def test_method_property_returns_method(self, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
            mfa._method = "FT"
        assert mfa2treelist_obj.method == "FT"

    def test_method_property_raises_with_multiple_methods(self, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        methods = ["FT", "iqtree"]
        for mfa, method in zip(mfa2treelist_obj, methods):
            mfa._tree = mock_tree
            mfa._method = method
        with pytest.raises(RuntimeError):
            _ = mfa2treelist_obj.method

    def test_toverrs_property_returns_list(self, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        expect = []
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
            val = random.random()
            expect.append(val)
            mfa._toverr = val
        toverrs = mfa2treelist_obj.toverrs
        assert isinstance(toverrs, list)
        assert len(toverrs) == 2
        assert toverrs == expect

    def test_saturations_property_returns_list(self, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        expect = []
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
            val = random.random()
            expect.append(val)
            mfa._saturation = val
        saturations = mfa2treelist_obj.saturations
        assert isinstance(saturations, list)
        assert len(saturations) == 2
        assert saturations == expect

    def test_sort(self, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        expect = []
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
            val = random.random()
            expect.append(val)
            mfa._toverr = val
        mfa2treelist_obj.sort()
        toverrs = mfa2treelist_obj.toverrs
        assert toverrs == sorted(expect, reverse=True)

    def test_build_ft(
        self,
        mock_build_env,
        mfa2treelist_obj,
    ):
        mfa2treelist_obj.load()

        mfa2treelist_obj.build("ft")
        for mfa2tree in mfa2treelist_obj:
            assert isinstance(mfa2tree.tree, Tree)

    @pytest.mark.parametrize("jobs, msg", [[1, "Sequential mode"], [2, "Multiprocesses mode"]])
    def test_build_flow_success(self, monkeypatch, mfa2treelist_obj, jobs, msg, caplog: pytest.LogCaptureFixture):
        def mock_helper(obj, method, output, model, noml, bs, scfl, seed, threads):
            return

        monkeypatch.setattr("phyling.lib.tree._build_helper", mock_helper)

        mfa2treelist_obj.load()
        with caplog.at_level("DEBUG"):
            mfa2treelist_obj.build("ft", jobs=jobs, threads=2)

        assert f"{msg}" in caplog.text
        assert f"Progress: {len(PEP_MSA)} / {len(PEP_MSA)}" in caplog.text

    def test_build_flow_failed(self, monkeypatch, mfa2treelist_obj):
        def mock_helper(obj, method, output, model, noml, bs, scfl, seed, threads):
            raise RuntimeError("Fake runtime error")

        monkeypatch.setattr("phyling.lib.tree._build_helper", mock_helper)

        mfa2treelist_obj.load()
        with pytest.raises(RuntimeError, match="Fake runtime error"):
            mfa2treelist_obj.build("ft", jobs=1, threads=2)

    def test_compute_toverr_raises_without_tree(self, mfa2treelist_obj):
        mfa2treelist_obj.load()
        with pytest.raises(AttributeError):
            mfa2treelist_obj.compute_toverr()

    def test_compute_toverr(self, monkeypatch, mfa2treelist_obj):
        mock_compute_toverr = MagicMock()
        val = random.random()
        mock_compute_toverr.return_value = val
        monkeypatch.setattr("phyling.external._libphykit.compute_toverr", mock_compute_toverr)

        mock_tree = MagicMock(spec=Tree)
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
        mfa2treelist_obj.compute_toverr()
        for mfa in mfa2treelist_obj:
            assert mfa.toverr == val

    @pytest.mark.parametrize("threads, msg", [[1, "Sequential mode"], [2, "Multiprocesses mode"]])
    def test_compute_toverr_flow(self, monkeypatch, mfa2treelist_obj, threads, msg, caplog: pytest.LogCaptureFixture):
        def mock_helper(obj):
            return

        monkeypatch.setattr("phyling.lib.tree._compute_toverr_helper", mock_helper)

        mfa2treelist_obj.load()
        with caplog.at_level("DEBUG"):
            mfa2treelist_obj.compute_toverr(threads=threads)

        assert f"{msg}" in caplog.text

    def test_compute_saturation_raises_without_tree(self, mfa2treelist_obj):
        mfa2treelist_obj.load()
        with pytest.raises(AttributeError):
            mfa2treelist_obj.compute_saturation()

    def test_compute_saturation(self, monkeypatch, mfa2treelist_obj):
        mock_saturation_class = MagicMock()
        val = random.random()
        mock_saturation_instance = mock_saturation_class.return_value
        mock_saturation_instance.compute_saturation.return_value = val
        monkeypatch.setattr("phyling.external._libphykit.Saturation", mock_saturation_class)

        mock_tree = MagicMock(spec=Tree)
        for mfa in mfa2treelist_obj:
            mfa._tree = mock_tree
        mfa2treelist_obj.compute_saturation()
        for mfa in mfa2treelist_obj:
            assert mfa.saturation == val

    @pytest.mark.parametrize("threads, msg", [[1, "Sequential mode"], [2, "Multiprocesses mode"]])
    def test_compute_saturation_flow(self, monkeypatch, mfa2treelist_obj, threads, msg, caplog: pytest.LogCaptureFixture):
        def mock_helper(obj):
            return

        monkeypatch.setattr("phyling.lib.tree._compute_saturation_helper", mock_helper)

        mfa2treelist_obj.load()
        with caplog.at_level("DEBUG"):
            mfa2treelist_obj.compute_saturation(threads=threads)

        assert f"{msg}" in caplog.text


@pytest.fixture
def mock_get_consensus_tree_env(monkeypatch, tmp_path):
    mock_astral = MagicMock()
    result_file = tmp_path / "consensus.nw"
    result_file.write_text("consensus", encoding="utf-8")
    mock_astral.return_value.result = result_file

    mock_phylo = MagicMock()
    mock_phylo.read.side_effect = lambda path, format: MagicMock(
        spec=Tree,
        _mock_file_content=path.read_text(encoding="utf-8"),  # Optional: store the text inside the mock if needed
    )

    monkeypatch.setattr("phyling.external.Astral", mock_astral)
    monkeypatch.setattr("phyling.lib.tree.Phylo", mock_phylo)


class TestMFA2TreeListConsensusTree:
    def test_returns_tree(self, mock_get_consensus_tree_env, mfa2treelist_obj):
        mock_tree = MagicMock(spec=Tree)
        for mfa2tree in mfa2treelist_obj:
            mfa2tree._tree = mock_tree

        result = mfa2treelist_obj.get_consensus_tree(None)
        assert isinstance(result, Tree)

    def test_get_consensus_tree_returns_path(self, mock_get_consensus_tree_env, mfa2treelist_obj, tmp_path):
        mock_tree = MagicMock(spec=Tree)
        for mfa2tree in mfa2treelist_obj:
            mfa2tree._tree = mock_tree

        result = mfa2treelist_obj.get_consensus_tree(tmp_path)
        assert isinstance(result, Path)


class TestMFA2TreeListConcat:
    def test_concat_returns_paths(self, mfa2treelist_obj, tmp_path):
        result = mfa2treelist_obj.concat(tmp_path)
        assert len(result) == 2
        assert result[0].name == "concat_alignments.mfa"
        assert result[1].name == "concat_alignments.partition"
        assert len(re.findall(r"(^|\n)>", Path(result[0]).read_text())) == 5


@pytest.fixture(scope="class")
def mock_bootstrap_env(class_monkeypatch, tmpdir_factory):
    shared_tmp = tmpdir_factory.mktemp("tree_data")

    mock_ufboot = MagicMock()
    ufboot_tree = shared_tmp / "ufboot_tree.nw"
    ufboot_tree.write_text("ufboot_tree", encoding="utf-8")
    mock_ufboot.return_value.result = ufboot_tree

    mock_phylo = MagicMock()
    mock_phylo.read.side_effect = lambda path, format: MagicMock(
        spec=Tree,
        _mock_file_content=path.read_text(encoding="utf-8"),
    )

    class_monkeypatch.setattr("phyling.external.UFBoot", mock_ufboot)
    class_monkeypatch.setattr("phyling.lib.tree.Phylo", mock_phylo)


class TestBootstrap:
    def test_bootstrap_returns_tree(self, mock_bootstrap_env, mfa2tree_obj):
        mfa2tree_obj.load()
        result = bootstrap(mfa2tree_obj, "tests/data/tree.nw")  # type: ignore
        assert isinstance(result, Tree)
        assert result._mock_file_content == "ufboot_tree"  # type: ignore

    def test_bootstrap_returns_path(self, mock_bootstrap_env, mfa2tree_obj, tmp_path):
        mfa2tree_obj.load()
        result = bootstrap(mfa2tree_obj, "tests/data/tree.nw", output=tmp_path)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "ufboot_tree"


@pytest.fixture(scope="class")
def mock_branch_concordance_env(class_monkeypatch, tmpdir_factory):
    shared_tmp = tmpdir_factory.mktemp("tree_data")
    concordance_tree = shared_tmp / "concordance_tree.nw"

    mock_concordance_instance = MagicMock()
    mock_concordance_instance.result = concordance_tree

    def mock_run():
        concordance_tree.write_text("concordance_tree", encoding="utf-8")

    mock_concordance_instance.run.side_effect = mock_run

    mock_concordance_class = MagicMock()
    mock_concordance_class.return_value = mock_concordance_instance

    mock_phylo = MagicMock()
    mock_phylo.read.side_effect = lambda path, format: MagicMock(
        spec=Tree,
        _mock_file_content=path.read_text(encoding="utf-8"),
    )

    class_monkeypatch.setattr("phyling.external.Concordance", mock_concordance_class)
    class_monkeypatch.setattr("phyling.lib.tree.Phylo", mock_phylo)


class TestBranchConcordance:
    def test_branch_concordance_returns_tree(self, mock_branch_concordance_env, mfa2tree_obj, tmp_path):
        mfa2tree_obj.load()
        result = branch_concordance(mfa2tree_obj, "tests/data/tree.nw")  # type: ignore
        assert isinstance(result, Tree)
        assert result._mock_file_content == "concordance_tree"  # type: ignore

    def test_branch_concordance_returns_path(self, mock_branch_concordance_env, mfa2tree_obj, tmp_path):
        mfa2tree_obj.load()
        result = branch_concordance(mfa2tree_obj, "tests/data/tree.nw", output=tmp_path)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "concordance_tree"


class TestFillMissingTaxon:
    def test_fill_missing_taxon_adds_gap_sequences(self):
        msa = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="taxon1", description=""),
                SeqRecord(Seq("ACGT"), id="taxon2", description=""),
            ]
        )
        samples = ["taxon1", "taxon2", "taxon3"]
        result = fill_missing_taxon(samples, msa)
        ids = [seq.id for seq in result]
        assert "taxon3" in ids
        taxon3_seq = next(seq for seq in result if seq.id == "taxon3")
        assert str(taxon3_seq.seq) == "----"

    def test_fill_missing_taxon_no_missing(self):
        msa = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="taxon1", description=""),
                SeqRecord(Seq("ACGT"), id="taxon2", description=""),
            ]
        )
        samples = ["taxon1", "taxon2"]
        result = fill_missing_taxon(samples, msa)
        assert len(result) == 2

    def test_fill_missing_taxon_returns_sorted(self):
        msa = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="taxon3", description=""),
                SeqRecord(Seq("ACGT"), id="taxon1", description=""),
            ]
        )
        samples = ["taxon1", "taxon2", "taxon3"]
        result = fill_missing_taxon(samples, msa)
        ids = [seq.id for seq in result]
        assert ids == sorted(ids)

    def test_fill_missing_taxon_gap_length_matches_alignment(self):
        msa = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGTACGT"), id="taxon1", description=""),
                SeqRecord(Seq("ACGTACGT"), id="taxon2", description=""),
            ]
        )
        samples = ["taxon1", "taxon2", "taxon3"]
        result = fill_missing_taxon(samples, msa)
        taxon3_seq = next(seq for seq in result if seq.id == "taxon3")
        assert len(str(taxon3_seq.seq)) == 8
        assert str(taxon3_seq.seq) == "--------"

    def test_fill_missing_taxon_multiple_missing(self):
        msa = MultipleSeqAlignment(
            [
                SeqRecord(Seq("ACGT"), id="taxon1", description=""),
            ]
        )
        samples = ["taxon1", "taxon2", "taxon3", "taxon4"]
        result = fill_missing_taxon(samples, msa)
        assert len(result) == 4
        ids = [seq.id for seq in result]
        assert "taxon2" in ids
        assert "taxon3" in ids
        assert "taxon4" in ids
