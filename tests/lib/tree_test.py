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


@pytest.fixture
def real_mfa2tree(path_pep_msa) -> MFA2Tree:
    """Create and load a MFA2Tree object."""
    return MFA2Tree(path_pep_msa[0], seqtype="pep")


@pytest.fixture
def mfa2treelist_obj() -> MFA2TreeList:
    def mock_mfa2tree_factory(file, seqtype="pep"):
        mfa2tree = MagicMock(spec=MFA2Tree)
        mfa2tree.name = Path(file).name
        mfa2tree.file = Path(file).absolute()
        mfa2tree.seqtype = seqtype
        return mfa2tree

    obj_1 = mock_mfa2tree_factory("/mock/path/MSA1.faa", seqtype="pep")
    obj_2 = mock_mfa2tree_factory("/mock/path/MSA2.faa", seqtype="pep")
    obj_1.tree = MagicMock(spec=Tree)
    obj_2.tree = MagicMock(spec=Tree)
    return MFA2TreeList([obj_1, obj_2], seqtype="pep")


@pytest.fixture(scope="class")
def class_monkeypatch():
    """A class-scoped monkeypatch fixture."""
    mp = MonkeyPatch()
    yield mp
    mp.undo()


@pytest.fixture(scope="class")
def class_temp_dir(tmpdir_factory):
    shared_tmp = tmpdir_factory.mktemp("tree_data")
    return shared_tmp


@pytest.fixture(scope="class")
def spy_phylo_read(class_monkeypatch: pytest.MonkeyPatch):
    def fake_read(path, format):
        mock_tree = MagicMock(spec=Tree)
        content = path.read_text(encoding="utf-8")
        mock_tree.__str__.return_value = content

        return mock_tree

    mock_phylo = MagicMock()
    mock_phylo.read.side_effect = fake_read
    class_monkeypatch.setattr("phyling.lib.tree.Phylo", mock_phylo)


@pytest.fixture(scope="class")
def mock_bootstrap(class_temp_dir, spy_phylo_read, class_monkeypatch: pytest.MonkeyPatch):
    ufboot = MagicMock()
    ufboot_tree = class_temp_dir / "ufboot.nw"
    ufboot_tree.write_text("ufboot", encoding="utf-8")
    ufboot.return_value.result = ufboot_tree

    class_monkeypatch.setattr("phyling.external.UFBoot", ufboot)


@pytest.fixture(scope="class")
def mock_branch_concordance(class_temp_dir, spy_phylo_read, class_monkeypatch: pytest.MonkeyPatch):
    concordance_tree = class_temp_dir / "concord.nw"

    mock_concordance_instance = MagicMock()
    mock_concordance_instance.result = concordance_tree

    def mock_run():
        concordance_tree.write_text("concord", encoding="utf-8")

    mock_concordance_instance.run.side_effect = mock_run

    mock_concordance_class = MagicMock()
    mock_concordance_class.return_value = mock_concordance_instance

    class_monkeypatch.setattr("phyling.external.Concordance", mock_concordance_class)


@pytest.fixture(scope="class")
def spy_build_env(class_monkeypatch, mock_bootstrap, mock_branch_concordance, spy_phylo_read, class_temp_dir):
    """Mocks all external runners, utilities, and file readers inside tree_builder."""
    # 2. Create the MagicMocks
    mock_modelfinder = MagicMock()
    mock_modelfinder.return_value.result = "model"

    mock_fasttree = MagicMock()
    ft_tree = class_temp_dir / "ft_tree.nw"
    ft_tree.write_text("ft_tree", encoding="utf-8")
    mock_fasttree.return_value.result = ft_tree

    mock_raxml = MagicMock()
    raxml_tree = class_temp_dir / "raxml_tree.nw"
    raxml_tree.write_text("raxml_tree", encoding="utf-8")
    mock_raxml.return_value.result = raxml_tree

    mock_iqtree = MagicMock()
    iqtree_tree = class_temp_dir / "iqtree_tree.nw"
    iqtree_tree.write_text("iqtree_tree", encoding="utf-8")
    mock_iqtree.return_value.result = iqtree_tree

    class_monkeypatch.setattr("phyling.external.ModelFinder", mock_modelfinder)
    class_monkeypatch.setattr("phyling.external.FastTree", mock_fasttree)
    class_monkeypatch.setattr("phyling.external.Raxml", mock_raxml)
    class_monkeypatch.setattr("phyling.external.Iqtree", mock_iqtree)


# ---------------------------------------------------------------------------
# TestMFA2Tree
# ---------------------------------------------------------------------------


def make_tree_with_toverr(msa, toverr_val) -> MFA2Tree:
    obj = MFA2Tree(msa)
    obj.load()
    obj._tree = MagicMock(spec=Tree)
    obj._toverr = toverr_val
    return obj


class TestMFA2Tree:
    def test_init_with_file_only(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(path_pep_msa[0])
        assert obj.file == path_pep_msa[0].absolute()
        assert obj.name == path_pep_msa[0].name

    def test_init_with_file_and_name(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(path_pep_msa[0], "custom_name")
        assert obj.file == path_pep_msa[0].absolute()
        assert obj.name == "custom_name"

    def test_init_with_seqtype_pep(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(path_pep_msa[0], seqtype="pep")
        assert obj.seqtype == "pep"

    def test_init_with_seqtype_dna(self, path_cds_msa: list[Path]):
        obj = MFA2Tree(path_cds_msa[0], seqtype="dna")
        assert obj.seqtype == "dna"

    def test_init_with_seqtype_auto_pep(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(path_pep_msa[0], seqtype="AUTO")
        assert obj.seqtype == "pep"

    def test_init_with_seqtype_auto_cds(self, path_cds_msa: list[Path]):
        obj = MFA2Tree(path_cds_msa[0], seqtype="AUTO")
        assert obj.seqtype == "dna"

    def test_init_private_attrs_none(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(path_pep_msa[0])
        assert obj._method is None
        assert obj._tree is None
        assert obj._toverr is None
        assert obj._saturation is None

    def test_init_with_string_path(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(str(path_pep_msa[0]))
        assert obj.file == path_pep_msa[0].absolute()

    def test_load_pep(self, path_pep_msa: list[Path]):
        obj = MFA2Tree(path_pep_msa[0], seqtype="pep")
        obj.load()
        assert obj._data is not None
        assert isinstance(obj._data, MultipleSeqAlignment)
        assert obj._data.annotations["seqtype"] == "pep"
        assert obj._data.annotations["seqname"] == path_pep_msa[0].name

    def test_load_cds(self, path_cds_msa: list[Path]):
        obj = MFA2Tree(path_cds_msa[0], seqtype="dna")
        obj.load()
        assert obj._data is not None
        assert isinstance(obj._data, MultipleSeqAlignment)
        assert obj._data.annotations["seqtype"] == "dna"
        assert obj._data.annotations["seqname"] == path_cds_msa[0].name

    def test_len(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(Exception):
            len(real_mfa2tree)
        real_mfa2tree.load()
        length = len(real_mfa2tree)
        assert length > 0

    def test_gt_true(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.9)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        assert obj1 > obj2

    def test_gt_false(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.3)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        assert not (obj1 > obj2)

    def test_ge_true(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        assert obj1 >= obj2

    def test_lt_true(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.3)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        assert obj1 < obj2

    def test_le_true(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        assert obj1 <= obj2

    def test_eq_true(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        assert obj1 == obj2

    def test_eq_false(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        obj2 = make_tree_with_toverr(path_pep_msa[0], 0.6)
        assert obj1 != obj2

    def test_eq_different_type(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        with pytest.raises(TypeError):
            obj1 == "not_an_mfa2tree"

    def test_gt_wrong_type(self, path_pep_msa: list[Path]):
        obj1 = make_tree_with_toverr(path_pep_msa[0], 0.5)
        with pytest.raises(TypeError):
            obj1 > "not_an_mfa2tree"  # type: ignore

    def test_tree_property_raises_without_build(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(AttributeError):
            real_mfa2tree.tree

    def test_method_property_raises_without_build(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(AttributeError):
            real_mfa2tree.method

    def test_toverr_property_raises_without_build(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(AttributeError):
            real_mfa2tree.toverr

    def test_saturation_property_raises_without_build(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(AttributeError):
            real_mfa2tree.saturation

    def test_tree_property_returns_tree(self, real_mfa2tree: MFA2Tree):
        mock_tree = MagicMock(spec=Tree)
        real_mfa2tree._tree = mock_tree
        assert real_mfa2tree.tree is mock_tree

    def test_method_property_returns_method(self, real_mfa2tree: MFA2Tree):
        real_mfa2tree._tree = MagicMock(spec=Tree)
        real_mfa2tree._method = "FT"
        assert real_mfa2tree.method == "FT"

    def test_toverr_property_returns_value(self, real_mfa2tree: MFA2Tree):
        val = random.random()
        real_mfa2tree._tree = MagicMock(spec=Tree)
        real_mfa2tree._toverr = val
        assert real_mfa2tree.toverr == val

    def test_saturation_property_returns_value(self, real_mfa2tree: MFA2Tree):
        val = random.random()
        real_mfa2tree._tree = MagicMock(spec=Tree)
        real_mfa2tree._saturation = val
        assert real_mfa2tree.saturation == val

    def test_toverr_raises_without_toverr(self, real_mfa2tree: MFA2Tree):
        real_mfa2tree._tree = MagicMock(spec=Tree)
        real_mfa2tree._toverr = None
        with pytest.raises(AttributeError):
            _ = real_mfa2tree.toverr

    def test_saturation_raises_without_saturation(self, real_mfa2tree: MFA2Tree):
        real_mfa2tree._tree = MagicMock(spec=Tree)
        real_mfa2tree._saturation = None
        with pytest.raises(AttributeError):
            _ = real_mfa2tree.saturation


class TestMFA2TreeBuild:
    def test_build_invalid_method(self, spy_build_env, real_mfa2tree: MFA2Tree, tmp_path):
        with pytest.raises(KeyError):
            real_mfa2tree.build("invalid_method", tmp_path)  # pyrefly: ignore [bad-argument-type]

    @pytest.mark.parametrize("method", ["ft", "raxml", "iqtree"])
    def test_build_returns_tree(self, spy_build_env, real_mfa2tree: MFA2Tree, method):
        result = real_mfa2tree.build(method)
        assert isinstance(result, Tree)
        assert f"{method}_tree" in str(result)

    @pytest.mark.parametrize("method", ["ft", "raxml", "iqtree"])
    def test_build_returns_path(self, spy_build_env, real_mfa2tree: MFA2Tree, method, tmp_path):
        result = real_mfa2tree.build(method, tmp_path)
        assert isinstance(result, Path)
        assert result.name == f"{method}_tree.nw"

    def test_build_ft_partition_raises(self, spy_build_env, real_mfa2tree: MFA2Tree, tmp_path):
        partition_file = tmp_path / "partition.txt"
        partition_file.touch()

        with pytest.raises(ValueError, match="does not support partitioning"):
            real_mfa2tree.build("ft", tmp_path / "output", model=partition_file)

    def test_build_with_bs(self, spy_build_env, real_mfa2tree: MFA2Tree, tmp_path):
        result = real_mfa2tree.build("ft", tmp_path, bs=1000)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "ufboot"

    def test_build_with_scfl(self, spy_build_env, real_mfa2tree: MFA2Tree, tmp_path):
        result = real_mfa2tree.build("ft", tmp_path, scfl=100)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "concord"


class TestMFA2TreeComputeToverr:
    def test_compute_toverr_raises_without_tree(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(AttributeError):
            real_mfa2tree.compute_toverr()

    def test_compute_toverr_sets_value(self, monkeypatch: pytest.MonkeyPatch, real_mfa2tree: MFA2Tree):
        real_mfa2tree._tree = MagicMock(spec=Tree)
        mock_compute_toverr = MagicMock()
        val = random.random()
        mock_compute_toverr.return_value = val
        monkeypatch.setattr("phyling.external._libphykit.compute_toverr", mock_compute_toverr)
        real_mfa2tree.compute_toverr()
        assert real_mfa2tree.toverr == val


class TestMFA2TreeComputeSaturation:
    def test_compute_saturation_raises_without_tree(self, real_mfa2tree: MFA2Tree):
        with pytest.raises(AttributeError):
            real_mfa2tree.compute_saturation()

    def test_compute_saturation_sets_value(self, monkeypatch: pytest.MonkeyPatch, real_mfa2tree: MFA2Tree):
        real_mfa2tree._tree = MagicMock(spec=Tree)
        mock_saturation_class = MagicMock()
        val = random.random()
        mock_saturation_instance = mock_saturation_class.return_value
        mock_saturation_instance.compute_saturation.return_value = val
        monkeypatch.setattr("phyling.external._libphykit.Saturation", mock_saturation_class)
        real_mfa2tree.compute_saturation()
        assert real_mfa2tree.saturation == val


# ---------------------------------------------------------------------------
# TestMFA2TreeList
# ---------------------------------------------------------------------------


class TestMFA2TreeList:
    def test_init_empty(self):
        obj = MFA2TreeList()
        assert len(obj) == 0

    def test_init_pep(self, path_pep_msa: list[Path]):
        obj = MFA2TreeList(path_pep_msa)
        assert len(obj) == 2
        assert obj.seqtype == "pep"

    def test_init_cds(self, path_cds_msa: list[Path]):
        obj = MFA2TreeList(path_cds_msa)
        assert len(obj) == 2
        assert obj.seqtype == "dna"

    def test_init_with_names(self, path_pep_msa: list[Path]):
        names = ["name1", "name2"]
        obj = MFA2TreeList(path_pep_msa, names, seqtype="pep")
        assert len(obj) == 2
        for mfa2tree in obj:
            assert mfa2tree.name in names

    def test_getitem_by_index(self, mfa2treelist_obj: MFA2TreeList):
        item = mfa2treelist_obj[0]
        assert isinstance(item, MFA2Tree)

    def test_getitem_by_slice(self, mfa2treelist_obj: MFA2TreeList):
        subset = mfa2treelist_obj[0:1]
        assert isinstance(subset, MFA2TreeList)
        assert len(subset) == 1

    def test_getitem_by_name(self, path_pep_msa: list[Path]):
        names = ["name1", "name2"]
        obj = MFA2TreeList(path_pep_msa, names, seqtype="pep")
        item = obj["name1"]
        assert isinstance(item, MFA2Tree)
        assert item.file == path_pep_msa[0].absolute()

    def test_trees_property_raises_without_build(self, path_pep_msa: list[Path]):
        obj = MFA2TreeList(path_pep_msa, seqtype="pep")
        with pytest.raises(AttributeError):
            obj.trees

    def test_trees_property_returns_list(self, mfa2treelist_obj: MFA2TreeList):
        trees = mfa2treelist_obj.trees
        assert isinstance(trees, list)
        assert len(trees) == 2
        assert isinstance(trees[0], Tree)

    def test_sort(self, path_pep_msa: list[Path]):
        obj = MFA2TreeList(path_pep_msa, seqtype="pep")
        mock_tree = MagicMock(spec=Tree)
        expect = []
        for mfa in obj:
            mfa._tree = mock_tree
            val = random.random()
            expect.append(val)
            mfa._toverr = val
        obj.sort()
        toverrs = obj.toverrs
        assert toverrs == sorted(expect, reverse=True)

    def test_method_property_returns_method(self, mfa2treelist_obj: MFA2TreeList):
        for mfa in mfa2treelist_obj:
            mfa.method = "FT"  # pyrefly: ignore [read-only]
        assert mfa2treelist_obj.method == "FT"

    def test_method_property_raises_with_multiple_methods(self, mfa2treelist_obj: MFA2TreeList):
        methods = ["FT", "iqtree"]
        for mfa, method in zip(mfa2treelist_obj, methods):
            mfa.method = method  # pyrefly: ignore [read-only]
        with pytest.raises(RuntimeError):
            mfa2treelist_obj.method

    def test_toverrs_property_returns_list(self, mfa2treelist_obj: MFA2TreeList):
        expect = []
        for mfa in mfa2treelist_obj:
            val = random.random()
            expect.append(val)
            mfa.toverr = val  # pyrefly: ignore [read-only]
        toverrs = mfa2treelist_obj.toverrs
        assert isinstance(toverrs, list)
        assert len(toverrs) == 2
        assert toverrs == expect

    def test_saturations_property_returns_list(self, mfa2treelist_obj: MFA2TreeList):
        expect = []
        for mfa in mfa2treelist_obj:
            val = random.random()
            expect.append(val)
            mfa.saturation = val  # pyrefly: ignore [read-only]
        saturations = mfa2treelist_obj.saturations
        assert isinstance(saturations, list)
        assert len(saturations) == 2
        assert saturations == expect

    @pytest.mark.parametrize("jobs, msg", [[1, "Sequential mode"], [2, "Multiprocesses mode"]])
    def test_build_flow_success(self, mfa2treelist_obj: MFA2TreeList, jobs: int, msg: str, caplog: pytest.LogCaptureFixture):
        with caplog.at_level("DEBUG"):
            mfa2treelist_obj.build("ft", jobs=jobs, threads=2)

        assert f"{msg}" in caplog.text
        assert f"Progress: {len(mfa2treelist_obj)} / {len(mfa2treelist_obj)}" in caplog.text
        for mfa in mfa2treelist_obj:
            mfa.build.assert_called_once()  # pyrefly: ignore [missing-attribute]

    def test_build_flow_failed(self, monkeypatch: pytest.MonkeyPatch, mfa2treelist_obj: MFA2TreeList):
        def crashing_helper(*args, **kwargs):
            raise RuntimeError("Fake runtime error")

        monkeypatch.setattr("phyling.lib.tree._build_helper", crashing_helper)

        with pytest.raises(RuntimeError, match="Fake runtime error"):
            mfa2treelist_obj.build("ft", jobs=1, threads=2)

    @pytest.mark.parametrize("threads, msg", [[1, "Sequential mode"], [2, "Multiprocesses mode"]])
    def test_compute_toverr_flow(self, mfa2treelist_obj: MFA2TreeList, threads: int, msg: str, caplog: pytest.LogCaptureFixture):
        with caplog.at_level("DEBUG"):
            mfa2treelist_obj.compute_toverr(threads=threads)

        assert f"{msg}" in caplog.text
        for mfa in mfa2treelist_obj:
            mfa.compute_toverr.assert_called_once()  # pyrefly: ignore [missing-attribute]

    @pytest.mark.parametrize("threads, msg", [[1, "Sequential mode"], [2, "Multiprocesses mode"]])
    def test_compute_saturation_flow(
        self, mfa2treelist_obj: MFA2TreeList, threads: int, msg: str, caplog: pytest.LogCaptureFixture
    ):
        with caplog.at_level("DEBUG"):
            mfa2treelist_obj.compute_saturation(threads=threads)

        assert f"{msg}" in caplog.text
        for mfa in mfa2treelist_obj:
            mfa.compute_saturation.assert_called_once()  # pyrefly: ignore [missing-attribute]


@pytest.fixture
def mock_get_consensus_tree_env(spy_phylo_read, tmp_path, monkeypatch: pytest.MonkeyPatch):
    mock_astral = MagicMock()
    result_file = tmp_path / "consensus.nw"
    result_file.write_text("consensus", encoding="utf-8")
    mock_astral.return_value.result = result_file

    monkeypatch.setattr("phyling.external.Astral", mock_astral)


class TestMFA2TreeListConsensus:
    def test_returns_tree(self, mock_get_consensus_tree_env, mfa2treelist_obj: MFA2TreeList):
        mock_tree = MagicMock(spec=Tree)
        for mfa2tree in mfa2treelist_obj:
            mfa2tree._tree = mock_tree

        result = mfa2treelist_obj.get_consensus_tree(None)
        assert isinstance(result, Tree)

    def test_get_consensus_tree_returns_path(self, mock_get_consensus_tree_env, mfa2treelist_obj: MFA2TreeList, tmp_path):
        mock_tree = MagicMock(spec=Tree)
        for mfa2tree in mfa2treelist_obj:
            mfa2tree._tree = mock_tree

        result = mfa2treelist_obj.get_consensus_tree(tmp_path)
        assert isinstance(result, Path)


class TestMFA2TreeListConcat:
    def test_concat_returns_paths(self, path_pep_msa: list[Path], tmp_path):
        obj = MFA2TreeList(path_pep_msa, seqtype="pep")
        result = obj.concat(tmp_path)
        assert len(result) == 2
        assert result[0].name == "concat_alignments.mfa"
        assert result[1].name == "concat_alignments.partition"
        assert len(re.findall(r"(^|\n)>", Path(result[0]).read_text())) == 5


class TestBootstrap:
    def test_bootstrap_returns_tree(self, mock_bootstrap, real_mfa2tree: MFA2Tree, path_tree_file: Path):
        result = bootstrap(real_mfa2tree, path_tree_file)
        assert isinstance(result, Tree)
        assert str(result) == "ufboot"  # type: ignore

    def test_bootstrap_returns_path(self, mock_bootstrap, real_mfa2tree: MFA2Tree, path_tree_file: Path, tmp_path):
        result = bootstrap(real_mfa2tree, path_tree_file, tmp_path)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "ufboot"


class TestBranchConcordance:
    def test_branch_concordance_returns_tree(
        self, mock_branch_concordance, real_mfa2tree: MFA2Tree, path_tree_file: Path, tmp_path
    ):
        result = branch_concordance(real_mfa2tree, path_tree_file)  # type: ignore
        assert isinstance(result, Tree)
        assert str(result) == "concord"  # type: ignore

    def test_branch_concordance_returns_path(
        self, mock_branch_concordance, real_mfa2tree: MFA2Tree, path_tree_file: Path, tmp_path
    ):
        result = branch_concordance(real_mfa2tree, path_tree_file, output=tmp_path)
        assert isinstance(result, Path)
        assert result.read_text(encoding="utf-8") == "concord"


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
