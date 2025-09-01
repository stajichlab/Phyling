from __future__ import annotations

import logging
import os
from argparse import ArgumentParser
from itertools import product
from pathlib import Path
from shutil import copytree, rmtree

import pytest
from Bio import AlignIO

from phyling.libphyling import ALIGN_METHODS, FileExts, TreeMethods, TreeOutputFiles
from phyling.pipeline import align, download, filter, tree


@pytest.fixture
def parser() -> ArgumentParser:
    """Fixture to create a new parser for each test."""
    return ArgumentParser(add_help=False)


@pytest.fixture
def make_empty_dir():
    created_dirs = []

    def _make_dir(base_folder: Path, name: str = "test_dir"):
        test_dir = Path(base_folder) / name
        test_dir.mkdir(parents=True, exist_ok=True)
        created_dirs.append(test_dir)
        return test_dir

    yield _make_dir

    # cleanup after test
    for d in created_dirs:
        if d.exists():
            rmtree(d)


def is_subset(d1: dict, d2: dict):
    return all(k in d1 and d1[k] == v for k, v in d2.items())


@pytest.fixture
def copied_pep_align_outputs(tmp_path: Path, persistent_tmp_path: Path) -> Path:
    copytree(persistent_tmp_path / "pep_align", tmp_path, dirs_exist_ok=True)
    return tmp_path


@pytest.fixture
def copied_pep_align_filtered(tmp_path: Path, persistent_tmp_path: Path) -> Path:
    copytree(persistent_tmp_path / "pep_align_filtered", tmp_path, dirs_exist_ok=True)
    return tmp_path


@pytest.mark.order(1)
class TestDownload:
    def test_menu(self, parser: ArgumentParser):
        download.menu(parser)
        true_param_pairs = vars(parser.parse_args(["poxviridae_odb10"]))
        expected_param_pairs = {"markerset": "poxviridae_odb10"}
        assert is_subset(true_param_pairs, expected_param_pairs)

    def test_download_list(self):
        download.download("list")

    def test_download_with_valid_name(self, capsys: pytest.CaptureFixture):
        with capsys.disabled():
            download.download("tevenvirinae_odb10")
        download.download("list")
        captured: str = capsys.readouterr().out
        captured = captured.split("Datasets available on local:")
        assert "tevenvirinae_odb10" in captured[1]

    def test_download_with_invalid_name(self):
        with pytest.raises(RuntimeError) as excinfo:
            download.download("InvalidName")
        assert "Markerset not available" in str(excinfo.value)

    @pytest.mark.dependency(name="database_update")
    def test_download_update(self, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.INFO):
            download.download("bclasvirinae_odb10")
        assert "Local markerset outdated." in caplog.text

    @pytest.mark.dependency("database_update")
    def test_download_already_up_to_date(self, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.INFO):
            download.download("bclasvirinae_odb10")
        assert "Markerset already exists and up to date." in caplog.text


@pytest.mark.dependency("database_update")
@pytest.mark.order(2)
class TestAlign:
    inputs_pep = (
        "tests/data/pep/bgzf",
        [
            "tests/data/pep/bgzf/Anomala_cuprea_entomopoxvirus.faa.gz",
            "tests/data/pep/bgzf/Canarypox_virus.faa.gz",
            "tests/data/pep/bgzf/Cowpox_virus.faa.gz",
            "tests/data/pep/bgzf/Goatpox_virus.faa.gz",
            "tests/data/pep/bgzf/Variola_virus.faa.gz",
        ],
    )
    inputs_cds = (
        "tests/data/cds/bgzf",
        [
            "tests/data/cds/bgzf/Anomala_cuprea_entomopoxvirus.fna.gz",
            "tests/data/cds/bgzf/Canarypox_virus.fna.gz",
            "tests/data/cds/bgzf/Cowpox_virus.fna.gz",
            "tests/data/cds/bgzf/Goatpox_virus.fna.gz",
            "tests/data/cds/bgzf/Variola_virus.fna.gz",
        ],
    )
    additional_inputs_pep = ["tests/data/pep/Monkeypox_virus.faa.gz"]
    prob_inputs_cds = ["tests/data/cds/Monkeypox_virus_with_bad_seq.fna"]
    nontrim = (True, False)

    def test_menu(self, parser: ArgumentParser):
        align.menu(parser)
        true_param_pairs = vars(
            parser.parse_args(
                [
                    "--input_dir",
                    "I",
                    "--output",
                    "O",
                    "--markerset",
                    "M",
                    "--evalue",
                    "1e-10",
                    "--method",
                    "hmmalign",
                    "--non_trim",
                    "--threads",
                    "3",
                ]
            )
        )
        expected_param_pairs = {
            "inputs": Path("I"),
            "output": Path("O"),
            "markerset": Path("M"),
            "evalue": 1e-10,
            "method": "hmmalign",
            "non_trim": True,
            "threads": 3,
        }
        assert is_subset(true_param_pairs, expected_param_pairs)

    @pytest.mark.dependency(name="align")
    @pytest.mark.parametrize("inputs, method", list(product(inputs_pep + inputs_cds, ALIGN_METHODS)))
    def test_align(self, inputs: str | list[str], tmp_path: Path, method: str, persistent_tmp_path: Path):
        align.align(
            inputs,
            tmp_path,
            markerset="tests/database/poxviridae_odb10/hmms",
            method=method,
            non_trim=False,
            threads=1,
        )
        samples = set()
        for file in tmp_path.glob(f"*.{FileExts.ALN}"):
            for rec in AlignIO.read(file, "fasta"):
                samples.add(rec.id)
        expected_samples = set(
            [
                "Anomala_cuprea_entomopoxvirus",
                "Canarypox_virus",
                "Cowpox_virus",
                "Goatpox_virus",
                "Variola_virus",
            ]
        )
        assert samples == expected_samples

        if method == "hmmalign":
            if inputs == self.inputs_pep[0]:
                copytree(tmp_path, persistent_tmp_path / "pep_align")
            if inputs == self.inputs_cds[0]:
                copytree(tmp_path, persistent_tmp_path / "cds_align")

    @pytest.mark.parametrize("inputs, non_trim", list(product((inputs_pep[0], inputs_cds[0]), nontrim)))
    def test_align_non_trim(self, inputs: str, tmp_path: Path, non_trim: bool):
        align.align(
            inputs,
            tmp_path,
            markerset="tests/database/poxviridae_odb10/hmms",
            method="hmmalign",
            non_trim=non_trim,
            threads=1,
        )
        samples = set()
        for file in tmp_path.glob(f"*.{FileExts.ALN}"):
            for rec in AlignIO.read(file, "fasta"):
                samples.add(rec.id)
        expected_samples = set(
            [
                "Anomala_cuprea_entomopoxvirus",
                "Canarypox_virus",
                "Cowpox_virus",
                "Goatpox_virus",
                "Variola_virus",
            ]
        )
        assert samples == expected_samples

    def test_align_dir_as_inputs(self, make_empty_dir, tmp_path: Path):
        make_empty_dir(self.inputs_pep[0], "empty_dir")
        with pytest.raises(IsADirectoryError) as excinfo:
            align.align(self.inputs_pep[0], tmp_path, markerset="poxviridae_odb10")
        assert "Input files contain directory" in str(excinfo.value)

    def test_align_too_few_inputs(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            align.align(self.inputs_pep[1][1:4], tmp_path, markerset="poxviridae_odb10")
        assert "Requires at least 4 samples" in str(excinfo.value)

    def test_align_inputs_cds_invalid_length(self, tmp_path: Path, caplog: pytest.LogCaptureFixture):
        with caplog.at_level(logging.DEBUG):
            align.align(self.inputs_cds[1][1:4] + self.prob_inputs_cds, tmp_path, markerset="poxviridae_odb10")
        assert "seqs have invalid length." in caplog.text
        assert "lcl|NC_003310.1_cds_NP_536428.1_1" in caplog.text

    def test_align_invalid_markerset(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError) as excinfo:
            align.align(self.inputs_pep[0], tmp_path, markerset="InvalidName")
        assert "Markerset folder does not exist" in str(excinfo.value)

    def test_align_invalid_evalue(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            align.align(self.inputs_pep[0], tmp_path, markerset="poxviridae_odb10", evalue=1)
        assert "Invalid evalue" in str(excinfo.value)

    def test_align_invalid_method(self, tmp_path: Path):
        with pytest.raises(ValueError) as excinfo:
            align.align(self.inputs_pep[0], tmp_path, markerset="poxviridae_odb10", method="InvalidName")
        assert "Invalid method" in str(excinfo.value)

    @pytest.mark.dependency(depends=["align"])
    def test_rerun(self, copied_pep_align_outputs: Path):
        prev_output = copied_pep_align_outputs
        inputs = self.inputs_pep[1][1:] + self.additional_inputs_pep  # Replace Anomala_cuprea with Monkeypox

        align.align(
            inputs,
            prev_output,
            markerset="poxviridae_odb10",
            method="hmmalign",
            non_trim="False",
            threads=1,
        )
        samples = set()
        for file in prev_output.glob(f"*.{FileExts.ALN}"):
            for rec in AlignIO.read(file, "fasta"):
                samples.add(rec.id)
        expected_samples = set(
            [
                "Monkeypox_virus",
                "Canarypox_virus",
                "Cowpox_virus",
                "Goatpox_virus",
                "Variola_virus",
            ]
        )
        assert samples == expected_samples


@pytest.mark.parametrize(
    "thread_param, expected",
    [
        ((1, 1), (1, 1)),
        ((4, 1), (1, 1)),
        ((4, 7), (1, 7)),
        ((4, 8), (2, 4)),
        ((4, 16), (4, 4)),
        ((8, 96), (8, 4)),
    ],
)
def test_search_threads_check(thread_param, expected):
    assert align._search_threads_check(*thread_param) == expected


@pytest.mark.dependency(depends=["align"])
@pytest.mark.order(3)
class TestFilter:
    inputs_pep = (
        "pep_align",
        [
            "pep_align/14at10240.aa.mfa",
            "pep_align/10at10240.aa.mfa",
            "pep_align/22at10240.aa.mfa",
            "pep_align/43at10240.aa.mfa",
            "pep_align/48at10240.aa.mfa",
            "pep_align/53at10240.aa.mfa",
        ],
    )
    inputs_cds = (
        "cds_align",
        [
            "cds_align/10at10240.cds.mfa",
            "cds_align/14at10240.cds.mfa",
            "cds_align/22at10240.cds.mfa",
            "cds_align/43at10240.cds.mfa",
            "cds_align/48at10240.cds.mfa",
            "cds_align/53at10240.cds.mfa",
        ],
    )
    additional_inputs_pep = ["pep_align/74at10240.aa.mfa", "pep_align/81at10240.aa.mfa", "pep_align/118at10240.aa.mfa"]
    invalid_top_n = (1, 7)

    def test_menu(self, parser: ArgumentParser):
        filter.menu(parser)
        true_param_pairs = vars(
            parser.parse_args(["--input_dir", "I", "--output", "O", "--top_n_toverr", "10", "--threads", "3"])
        )
        expected_param_pairs = {
            "inputs": Path("I"),
            "output": Path("O"),
            "top_n_toverr": 10,
            "threads": 3,
        }
        assert is_subset(true_param_pairs, expected_param_pairs)

    @pytest.mark.parametrize("inputs", inputs_pep + inputs_cds)
    def test_filter(self, inputs: str | list[str], tmp_path: Path, persistent_tmp_path: Path):
        if isinstance(inputs, str):
            inputs = persistent_tmp_path / inputs
        else:
            inputs = [persistent_tmp_path / f for f in inputs]
        filter.filter(inputs, tmp_path, top_n_toverr=5, threads=os.cpu_count())
        treeness = tmp_path / TreeOutputFiles.TREENESS
        assert treeness.is_file()
        assert len(open(treeness).read().split("#")[1].strip().split("\n")) == 6
        assert sum(1 for _ in tmp_path.glob(f"*.{FileExts.ALN}")) == 5
        if inputs == [persistent_tmp_path / f for f in self.inputs_pep[1]]:
            copytree(tmp_path, persistent_tmp_path / "pep_align_filtered")

    @pytest.mark.parametrize("inputs", (inputs_pep[1][0:2], inputs_cds[1][0:2]))
    def test_filter_too_few_inputs(self, inputs: list[str], tmp_path: Path, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in inputs]
        with pytest.raises(ValueError) as excinfo:
            filter.filter(inputs, tmp_path, top_n_toverr=5)
        assert "Fewer than 3 inputs" in str(excinfo.value)

    @pytest.mark.parametrize("inputs, top_n", list(product((inputs_pep[1], inputs_pep[1][0:3]), invalid_top_n)))
    def test_filter_invalid_top_n(self, inputs: list[str], tmp_path: Path, top_n: int, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in inputs]
        with pytest.raises(ValueError) as excinfo:
            filter.filter(inputs, tmp_path, top_n_toverr=top_n)
        assert "Argument top_n_toverr out of range" in str(excinfo.value)

    def test_filter_top_n_equal_to_inputs(self, tmp_path: Path, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in self.inputs_pep[1]]
        with pytest.raises(SystemExit) as excinfo:
            filter.filter(inputs, tmp_path, top_n_toverr=6)
        assert "Argument top_n_toverr is equal to the number of inputs" in str(excinfo.value)

    def test_rerun(self, copied_pep_align_filtered: Path, persistent_tmp_path: Path):
        prev_output = copied_pep_align_filtered
        inputs = self.inputs_pep[1][1:] + self.additional_inputs_pep
        inputs = [persistent_tmp_path / f for f in inputs]
        filter.filter(inputs, prev_output, top_n_toverr=6)
        treeness = prev_output / TreeOutputFiles.TREENESS
        assert treeness.is_file()
        assert len(open(treeness).read().split("#")[1].strip().split("\n")) == 7
        assert sum(1 for _ in prev_output.glob(f"*.{FileExts.ALN}")) == 6


@pytest.mark.dependency(depends=["align"])
@pytest.mark.order(4)
class TestTree:
    inputs_pep = (
        "pep_align",
        [
            "pep_align/10at10240.aa.mfa",
            "pep_align/14at10240.aa.mfa",
            "pep_align/43at10240.aa.mfa",
        ],
    )
    inputs_cds = (
        "cds_align",
        [
            "cds_align/10at10240.cds.mfa",
            "cds_align/14at10240.cds.mfa",
            "cds_align/22at10240.cds.mfa",
        ],
    )
    method = tuple(m.name.lower() for m in TreeMethods)
    concat = (False, True)
    partition = (False, True)
    figure = (False, True)

    def test_menu(self, parser: ArgumentParser):
        tree.menu(parser)
        true_param_pairs = vars(
            parser.parse_args(
                [
                    "--input_dir",
                    "I",
                    "--output",
                    "O",
                    "--method",
                    "iqtree",
                    "--concat",
                    "--partition",
                    "--threads",
                    "3",
                ]
            )
        )
        expected_param_pairs = {
            "inputs": Path("I"),
            "output": Path("O"),
            "method": "iqtree",
            "concat": True,
            "partition": True,
            "threads": 3,
        }
        assert is_subset(true_param_pairs, expected_param_pairs)

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method, concat", list(product(inputs_pep + inputs_cds, method, concat)))
    def test_tree(self, inputs: str | list[str], tmp_path: Path, method: str, concat: bool, persistent_tmp_path: Path):
        if isinstance(inputs, str):
            inputs = persistent_tmp_path / inputs
        else:
            inputs = [persistent_tmp_path / f for f in inputs]
        tree.tree(inputs, tmp_path, method=method, concat=concat, threads=1)
        concat_file = tmp_path / TreeOutputFiles.CONCAT
        if concat:
            assert concat_file.is_file()
        else:
            assert not concat_file.exists()
        treefile = tmp_path / TreeOutputFiles.TREE_NW
        assert treefile.is_file()

    @pytest.mark.slow
    @pytest.mark.parametrize("inputs, method", list(product((inputs_pep[1], inputs_cds[1]), method[1:])))
    def test_tree_partition(self, inputs: list[str], tmp_path: Path, method: str, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in inputs]
        tree.tree(inputs, tmp_path, method=method, concat=True, partition=True, threads=1)
        concat_file = tmp_path / TreeOutputFiles.CONCAT
        partition_file = tmp_path / TreeOutputFiles.PARTITION
        treefile = tmp_path / TreeOutputFiles.TREE_NW
        assert concat_file.is_file()
        assert partition_file.is_file()
        assert treefile.is_file()

    @pytest.mark.parametrize("method", method[1:])
    def test_partition_no_concat(self, tmp_path: Path, method: str, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in self.inputs_pep[1]]
        with pytest.raises(ValueError) as excinfo:
            tree.tree(inputs, tmp_path, method=method, concat=False, partition=True)
        assert "Partition is not allowed in consensus mode" in str(excinfo.value)

    def test_partition_invalid_method(self, tmp_path: Path, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in self.inputs_pep[1]]
        with pytest.raises(ValueError) as excinfo:
            tree.tree(inputs, tmp_path, method="ft", concat=True, partition=True)
        assert f"Partition is not allowed with {TreeMethods.FT.method}" in str(excinfo.value)

    @pytest.mark.parametrize("inputs, method", list(product((inputs_pep[1][0], inputs_cds[1][0]), method)))
    def test_tree_single_input(self, inputs: str, tmp_path: Path, method: str, persistent_tmp_path: Path):
        inputs = persistent_tmp_path / inputs
        with pytest.raises(ValueError) as excinfo:
            tree.tree(inputs, tmp_path, method=method)
        assert "Found only 1 MSA fasta" in str(excinfo.value)

    def test_tree_figure(self, tmp_path: Path, persistent_tmp_path: Path):
        inputs = [persistent_tmp_path / f for f in self.inputs_pep[1]]
        tree.tree(inputs, tmp_path, figure=True)
        img_file = tmp_path / TreeOutputFiles.TREE_IMG
        assert img_file.is_file()
