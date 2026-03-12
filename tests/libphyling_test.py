# from __future__ import annotations

# import logging
# import shutil
# from copy import deepcopy
# from pathlib import Path

# import pytest
# from Bio.Align import MultipleSeqAlignment
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from pyhmmer.easel import DigitalSequence, DigitalSequenceBlock
# from pyhmmer.plan7 import HMM

# from phyling.align._markerset import HMMMarkerSet
# from phyling.align.align import SearchHitsManager
# from phyling.align.align import SampleList
# from phyling.align.align import SampleSeqs
# from phyling.pipeline.align import main, search, trim
# import phyling.exception as exception
# from phyling.pipeline._outputprecheck import AlignPrecheck


# @pytest.fixture(scope="class")
# def copy_lib_ckp(shared_tmpdir_class):
#     shutil.copy("tests/data/pep_align/.align.ckp", shared_tmpdir_class / ".align.ckp")
#     AlignPrecheck.setup(folder=shared_tmpdir_class)


# def need_search(samplelist: SampleList) -> int:
#     need_search = 0
#     for sample in samplelist:
#         if sample.is_scanned is False:
#             need_search += 1
#     return need_search


# class TestSampleSeqs:
#     def test_init_pep(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         assert a.name == "Monkeypox_virus"
#         assert a.seqtype == "pep"
#         assert a.is_scanned is False

#     def test_init_cds(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.fna.gz")
#         assert a.name == "Monkeypox_virus"
#         assert a.seqtype == "cds"
#         assert a.is_scanned is False

#     def test_init_cds_with_bad_seq(self, caplog: pytest.LogCaptureFixture):
#         with caplog.at_level(logging.WARNING):
#             SampleSeqs("tests/data/Monkeypox_virus_with_bad_seq.fna")
#         assert "seqs has invalid length" in caplog.text

#     def test_len(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         assert type(len(a)) is int

#     def test_eq_lt(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         b = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         c = SampleSeqs("tests/data/pep/Anomala_cuprea_entomopoxvirus.faa.gz")
#         assert a == b > c

#     def test_eq_lt_seqtypeerror(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         b = SampleSeqs("tests/data/Monkeypox_virus.fna.gz")
#         c = SampleSeqs("tests/data/pep/Anomala_cuprea_entomopoxvirus.faa.gz")
#         with pytest.raises(exception.SeqtypeError, match="Items represent different seqtypes"):
#             a == b
#         with pytest.raises(exception.SeqtypeError, match="Items represent different seqtypes"):
#             assert b > c

#     def test_getitem(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         assert type(a["pep"]) == DigitalSequenceBlock
#         assert type(a["pep", b"NP_536428.1"]) == DigitalSequence
#         assert a["pep", b"NP_536428.1"].alphabet.is_amino() is True

#         b = SampleSeqs("tests/data/Monkeypox_virus.fna.gz")
#         assert type(b["pep"]) == DigitalSequenceBlock
#         assert type(b["pep", b"lcl|NC_003310.1_cds_NP_536428.1_1"]) == DigitalSequence
#         assert b["pep", b"lcl|NC_003310.1_cds_NP_536428.1_1"].alphabet.is_amino() is True
#         assert type(b["cds"]) == DigitalSequenceBlock
#         assert type(b["cds", b"lcl|NC_003310.1_cds_NP_536428.1_1"]) == DigitalSequence
#         assert b["cds", b"lcl|NC_003310.1_cds_NP_536428.1_1"].alphabet.is_dna() is True

#     def test_getitem_errors(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         with pytest.raises(IndexError, match="Too many keys are given"):
#             a["pep", b"NP_536428.1", "Invalid"]
#         with pytest.raises(AttributeError, match="Cannot obtain the CDS sequences from a peptide inputs"):
#             a["cds"]
#         with pytest.raises(KeyError, match='First key only accepts "pep" or "cds"'):
#             a["rna"]

#     def test_scanned(self):
#         a = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         a.scanned()
#         assert a.is_scanned is True


# class TestSampleList:
#     files = tuple(file for file in Path("tests/data/pep").iterdir())

#     def test_init(self):
#         a = SampleList(self.files)
#         assert a.seqtype == "pep"
#         assert len(a) == len(self.files)
#         b = SampleList(Path("tests/data/cds").iterdir())
#         assert b.seqtype == "cds"
#         assert len(b) == len(self.files)

#     def test_init_gunzip(self):
#         a = SampleList(Path("tests/data/pep_gunzip").iterdir())
#         assert a.seqtype == "pep"

#     def test_init_typeerror(self):
#         with pytest.raises(TypeError, match="only accepts list of str/Path/SampleSeqs"):
#             SampleList([1, 2, 3])

#     def test_init_seqtypeerror(self):
#         with pytest.raises(exception.SeqtypeError, match="Inputs contains more than one seqtypes"):
#             SampleList(self.files + tuple(file for file in Path("tests/data/cds").iterdir()))

#     def test_init_identicalkeyerror(self):
#         with pytest.raises(exception.IdenticalKeyError, match="The following files share the same sample name"):
#             SampleList(self.files + tuple(file for file in Path("tests/data/pep_gunzip").iterdir()))

#     def test_getitem(self):
#         a = SampleList(self.files)
#         assert type(a["Cowpox_virus"]) == SampleSeqs
#         assert type(a[0]) == SampleSeqs
#         r = a[1:4]
#         assert type(r) == SampleList
#         assert len(r) == 3

#     def test_getitem_keyerror(self):
#         a = SampleList(self.files)
#         with pytest.raises(KeyError, match="Sample not found"):
#             a["Invalid_name"]

#     def test_append(self):
#         a = SampleList(self.files)
#         b = SampleSeqs("tests/data/Monkeypox_virus.faa.gz")
#         a.append(b)
#         assert len(a) == 6

#     def test_append_typeerror(self):
#         a = SampleList(self.files)
#         b = 1
#         with pytest.raises(TypeError, match="Can only append SampleSeqs object"):
#             a.append(b)

#     def test_append_seqtypeerror(self):
#         a = SampleList(self.files)
#         b = SampleSeqs("tests/data/Monkeypox_virus.fna.gz")
#         with pytest.raises(exception.SeqtypeError, match="Item represents different seqtype"):
#             a.append(b)

#     def test_append_identicalkeyerror(self):
#         a = SampleList(self.files)
#         b = SampleSeqs("tests/data/pep_gunzip/Cowpox_virus.faa")
#         with pytest.raises(exception.IdenticalKeyError, match="already exists"):
#             a.append(b)

#     def test_pop(self):
#         a = SampleList(self.files)
#         a.pop(3)
#         assert len(a) == 4

#     def test_update(self):
#         prev = SampleList(Path("tests/data/pep_gunzip").iterdir())
#         [file.scanned() for file in prev]
#         cur = SampleList(self.files)
#         assert cur.update(prev) is None
#         assert len(cur) == 5
#         assert need_search(cur) == 0

#     def test_update_addsample(self):
#         prev = SampleList(Path("tests/data/pep_gunzip").iterdir())
#         [file.scanned() for file in prev]
#         cur = SampleList(self.files + ("tests/data/Monkeypox_virus.faa.gz",))
#         assert cur.update(prev) is None
#         assert len(cur) == 6
#         assert need_search(cur) == 1

#     def test_update_dropsample(self):
#         prev = SampleList(Path("tests/data/pep_gunzip").iterdir())
#         prev.pop(-1)
#         prev.append(SampleSeqs("tests/data/Monkeypox_virus.faa.gz"))
#         [file.scanned() for file in prev]
#         cur = SampleList(self.files)
#         assert cur.update(prev) == ["Monkeypox_virus"]
#         assert len(cur) == 5
#         assert need_search(cur) == 1

#     def test_update_typeerror(self):
#         prev = [1, 2, 3]
#         cur = SampleList(self.files)
#         with pytest.raises(TypeError, match="Can only update with SampleList object"):
#             cur.update(prev)

#     def test_update_seqtypeerror(self):
#         prev = SampleList(Path("tests/data/cds").iterdir())
#         cur = SampleList(self.files)
#         with pytest.raises(exception.SeqtypeError, match="Item represents different seqtype"):
#             cur.update(prev)


# class TestHMMMarkerSet:
#     markerset_path = Path("tests/database/poxviridae_odb10/hmms")
#     cutoff_path = markerset_path.parent / "scores_cutoff"

#     def test_init(self):
#         markerset = HMMMarkerSet(self.markerset_path)
#         assert markerset.have_cutoff is False
#         markerset = HMMMarkerSet(self.markerset_path, self.cutoff_path)
#         assert markerset.have_cutoff is True

#     def test_init_typerror(self):
#         with pytest.raises(TypeError, match='Argument "folder" only accepts of str or Path'):
#             HMMMarkerSet([1, 2, 3])

#     def test_init_cutoff_not_found(self, caplog: pytest.LogCaptureFixture):
#         with pytest.raises(FileNotFoundError, match="HMM cutoff file not found"):
#             HMMMarkerSet(self.markerset_path, "Invalid_path")
#         with caplog.at_level(logging.WARNING):
#             HMMMarkerSet(self.markerset_path, "Invalid_path", raise_err=False)
#         assert "HMM cutoff file not found. Will use evalue instead" in caplog.text

#     def test_init_cutoff_not_match(self, caplog: pytest.LogCaptureFixture):
#         with pytest.raises(KeyError, match="HMM cutoff file doesn't match the markerset"):
#             HMMMarkerSet(self.markerset_path, "tests/database/alphaherpesvirinae_odb10/scores_cutoff")
#         with caplog.at_level(logging.WARNING):
#             HMMMarkerSet(self.markerset_path, "tests/database/alphaherpesvirinae_odb10/scores_cutoff", raise_err=False)
#         assert "HMM cutoff file doesn't match the markerset. Will use evalue instead" in caplog.text

#     def test_getitem(self):
#         markerset = HMMMarkerSet(self.markerset_path)
#         assert type(markerset[b"14at10240"]) == HMM
#         assert type(markerset[0]) == HMM
#         r = markerset[1:4]
#         assert len(r) == 3
#         assert type(r[0]) == HMM

#     def test_getitem_keyerror(self):
#         markerset = HMMMarkerSet(self.markerset_path)
#         with pytest.raises(KeyError, match="Sample not found"):
#             markerset[b"8at10293"]


# class TestOrthologs:
#     samplelist_pep = SampleList([Path("tests/data/pep/Cowpox_virus.faa.gz"), Path("tests/data/pep/Goatpox_virus.faa.gz")])
#     samplelist_cds = SampleList([Path("tests/data/cds/Cowpox_virus.fna.gz"), Path("tests/data/cds/Goatpox_virus.fna.gz")])

#     def test_init(self):
#         d = {b"hmm1": {("sample_1", b"seq_1"), ("sample_2", b"seq_2")}}
#         SearchHitsManager(d)

#     def test_init_typeerror(self):
#         d = {b"hmm1": {("sample_1", b"seq_1"), ("sample_2", 3)}}
#         with pytest.raises(TypeError):
#             SearchHitsManager(d)

#     def test_eq(self):
#         d1 = {b"hmm1": {("sample_1", b"seq_1"), ("sample_2", b"seq_2")}}
#         d2 = {b"hmm2": {("sample_1", b"seq_3"), ("sample_2", b"seq_3")}}
#         a = SearchHitsManager(d1)
#         b = SearchHitsManager(d1)
#         c = SearchHitsManager(d2)
#         assert a == b != c

#     def test_eq_typeerror(self):
#         d = {b"hmm1": {("sample_1", b"seq_1"), ("sample_2", b"seq_2")}}
#         with pytest.raises(TypeError):
#             SearchHitsManager(d) == d

#     def test_setitem(self):
#         a = SearchHitsManager()
#         a[b"hmm1"] = ("sample_1", b"seq_2")

#     def test_setitem_typeerror(self):
#         a = SearchHitsManager()
#         with pytest.raises(TypeError):
#             a[b"hmm1"] = {"sample_1", "seq_2"}
#         with pytest.raises(TypeError):
#             a[b"hmm1"] = ("sample_1", "seq_2")
#         with pytest.raises(TypeError):
#             a[b"hmm1"] = ("sample_1", 1)

#     def test_setitem_indexerror(self):
#         a = SearchHitsManager()
#         with pytest.raises(IndexError):
#             a[b"hmm1"] = ("sample_1", b"seq_2", b"seq_3")

#     def test_update(self):
#         d1 = {b"hmm1": {("sample_1", b"seq_1"), ("sample_2", b"seq_2")}}
#         d2 = {
#             b"hmm1": {("sample_1", b"seq_2"), ("sample_2", b"seq_3")},
#             b"hmm2": {("sample_1", b"seq_1"), ("sample_3", b"seq_3")},
#         }
#         r = {
#             b"hmm1": {("sample_1", b"seq_1"), ("sample_2", b"seq_2"), ("sample_1", b"seq_2"), ("sample_2", b"seq_3")},
#             b"hmm2": {("sample_1", b"seq_1"), ("sample_3", b"seq_3")},
#         }
#         a = SearchHitsManager(d1)
#         a.update(d2)
#         assert a.data == r

#     def test_query_pep(self):
#         ortho = SearchHitsManager({b"455at10240": {("Cowpox_virus", b"NP_619891.1"), ("Goatpox_virus", b"YP_001293259.1")}})
#         ortho.map(self.samplelist_pep)
#         r = ortho.query(b"455at10240")
#         assert type(r) is list
#         assert len(r) == 2
#         assert isinstance(r[0], DigitalSequence) and isinstance(r[1], DigitalSequence)
#         assert r[0].alphabet.is_amino() is True

#     def test_query_cds(self):
#         ortho = SearchHitsManager(
#             {
#                 b"455at10240": {
#                     ("Cowpox_virus", b"lcl|NC_003663.2_cds_NP_619891.1_108"),
#                     ("Goatpox_virus", b"lcl|NC_004003.1_cds_YP_001293259.1_65"),
#                 }
#             }
#         )
#         ortho.map(self.samplelist_cds)
#         r = ortho.query(b"455at10240")
#         assert type(r) is list
#         assert len(r) == 2
#         assert isinstance(r[0], DigitalSequence) and isinstance(r[1], DigitalSequence)
#         assert r[0].alphabet.is_amino() is True
#         r = ortho.query(b"455at10240", seqtype="cds")
#         assert type(r) is list
#         assert len(r) == 2
#         assert isinstance(r[0], DigitalSequence) and isinstance(r[1], DigitalSequence)
#         assert r[0].alphabet.is_dna() is True

#     def test_query_keyerror(self):
#         ortho = SearchHitsManager({b"455at10240": {("Cowpox_virus", b"NP_619891.1"), ("Goatpox_virus", b"YP_001293259.1")}})
#         ortho.map(self.samplelist_pep)
#         with pytest.raises(KeyError):
#             ortho.query(b"Invalid_name")

#     def test_query_seqtype_error(self):
#         ortho = SearchHitsManager({b"455at10240": {("Cowpox_virus", b"NP_619891.1"), ("Goatpox_virus", b"YP_001293259.1")}})
#         ortho.map(self.samplelist_pep)
#         with pytest.raises(AttributeError):
#             ortho.query(b"455at10240", seqtype="cds")

#     def test_map(self):
#         ortho = SearchHitsManager({b"455at10240": {("Cowpox_virus", b"NP_619891.1"), ("Goatpox_virus", b"YP_001293259.1")}})
#         ortho.map(self.samplelist_pep)
#         assert ortho.is_mapped is True

#     def test_filter_min_taxa(self):
#         ortho = SearchHitsManager(
#             {
#                 b"455at10240": {("Cowpox_virus", b"NP_619891.1"), ("Goatpox_virus", b"YP_001293259.1")},
#                 b"43at10240": {("Cowpox_virus", b"NP_619906.1")},
#                 b"14at10240": {("Goatpox_virus", b"YP_001293307.1")},
#             }
#         )
#         ortho.map(self.samplelist_pep)
#         r = ortho.filter(min_taxa=2)
#         assert type(r) == SearchHitsManager
#         assert len(r) == 1
#         assert b"455at10240" in r
#         assert r.is_mapped is False

#     def test_filter_droplist(self):
#         ortho = SearchHitsManager(
#             {
#                 b"455at10240": {("Cowpox_virus", b"NP_619891.1"), ("Goatpox_virus", b"YP_001293259.1")},
#                 b"43at10240": {("Cowpox_virus", b"NP_619906.1")},
#                 b"14at10240": {("Goatpox_virus", b"YP_001293307.1")},
#             }
#         )
#         ortho.map(self.samplelist_pep)
#         r = ortho.filter(droplist=["Goatpox_virus"])
#         assert type(r) == SearchHitsManager
#         assert len(r) == 2
#         assert len(r[b"455at10240"]) == 1
#         assert r.is_mapped is False


# class TestOutputPrecheck:
#     inputs = SampleList(Path("tests/data/pep").iterdir())
#     markerset_path = Path("tests/database/poxviridae_odb10/hmms")
#     cutoff_path = markerset_path.parent / "scores_cutoff"
#     markerset = HMMMarkerSet(markerset_path, cutoff_path)
#     params = {
#         "inputs": inputs.checksum,
#         "markerset": markerset.checksum,
#         "markerset_cutoff": "markerset cutoff" if markerset.have_cutoff else 1e-10,
#         "method": "hmmalign",
#         "non_trim": False,
#     }

#     def test_precheck_new(self, tmp_path: Path):
#         assert (tmp_path / "output").exists() is False
#         AlignPrecheck.setup(folder=tmp_path / "output")
#         _, r = AlignPrecheck.precheck(self.params, self.inputs)
#         assert r is None
#         assert (tmp_path / "output").is_dir() is True

#     def test_precheck_folder_exists(self, tmp_path: Path):
#         AlignPrecheck.setup(folder=tmp_path)
#         _, r = AlignPrecheck.precheck(self.params, self.inputs)
#         assert r is None
#         assert (tmp_path).is_dir() is True

#     def test_precheck_notadirectoryerror(self, tmp_path: Path):
#         (tmp_path / "output").touch()
#         AlignPrecheck.setup(folder=tmp_path / "output")
#         with pytest.raises(NotADirectoryError, match="already existed but not a folder"):
#             AlignPrecheck.precheck(self.params, self.inputs)

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_dropsample(self, caplog: pytest.LogCaptureFixture):
#         inputs = deepcopy(self.inputs)
#         popped_sample = inputs.pop(-1)
#         params = self.params.copy()
#         params.update(inputs=inputs.checksum)
#         with caplog.at_level(logging.INFO):
#             updated_inputs, _ = AlignPrecheck.precheck(params, inputs)
#         assert f"Remove hits corresponding to {popped_sample.name} from orthologs" in caplog.text
#         assert need_search(updated_inputs) == 0

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_addsample(self):
#         inputs = deepcopy(self.inputs)
#         inputs.append(SampleSeqs(Path("tests/data/Monkeypox_virus.faa.gz")))
#         params = self.params.copy()
#         params.update(inputs=inputs.checksum)
#         updated_inputs, _ = AlignPrecheck.precheck(params, inputs)
#         assert need_search(updated_inputs) == 1

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_change_method(self):
#         inputs = deepcopy(self.inputs)
#         params = self.params.copy()
#         params.update(method="muscle")
#         updated_inputs, _ = AlignPrecheck.precheck(params, inputs)
#         assert need_search(updated_inputs) == 0

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_change_non_trim(self):
#         inputs = deepcopy(self.inputs)
#         params = self.params.copy()
#         params.update(non_trim=True)
#         updated_inputs, _ = AlignPrecheck.precheck(params, inputs)
#         assert need_search(updated_inputs) == 0

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_change_hmm_cutoff(self):
#         inputs = deepcopy(self.inputs)
#         markerset = HMMMarkerSet(self.markerset_path)
#         params = self.params.copy()
#         params.update(markerset=markerset.checksum, markerset_cutoff=1e-10)
#         updated_inputs, _ = AlignPrecheck.precheck(params, inputs)
#         assert need_search(updated_inputs) == 5

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_params_error(self):
#         params = {"Invalid_key": "Invalid_value"}
#         with pytest.raises(KeyError, match="Params should contain keys"):
#             AlignPrecheck.precheck(params, self.inputs)

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_seqtype_error(self):
#         inputs = SampleList(Path("tests/data/cds").iterdir())
#         params = self.params.copy()
#         params.update(inputs=inputs.checksum)
#         with pytest.raises(SystemExit, match="Seqtype is changed. Aborted."):
#             AlignPrecheck.precheck(params, inputs)

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_file_not_changed_error(self):
#         with pytest.raises(SystemExit, match="Files not changed and parameters are identical to the previous run. Aborted."):
#             AlignPrecheck.precheck(self.params, self.inputs)

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_precheck_markerset_error(self):
#         markerset = HMMMarkerSet(Path("tests/database/alphaherpesvirinae_odb10/hmms"))
#         params = self.params.copy()
#         params.update(markerset=markerset.checksum)
#         with pytest.raises(SystemExit, match="Markerset is changed. Aborted."):
#             AlignPrecheck.precheck(params, self.inputs)

#     @pytest.mark.usefixtures("copy_lib_ckp")
#     def test_load_checkpoint(self):
#         params, samplelist, orthologs = AlignPrecheck.load_checkpoint()
#         assert isinstance(params, dict)
#         assert isinstance(samplelist, SampleList)
#         assert isinstance(orthologs, SearchHitsManager)

#     def test_save_checkpoint(self, tmp_path):
#         AlignPrecheck.setup(folder=tmp_path)

#         AlignPrecheck.save_checkpoint(self.params, self.inputs, SearchHitsManager())
#         a, b, c = AlignPrecheck.load_checkpoint()
#         assert isinstance(a, dict)
#         assert b == self.inputs
#         assert isinstance(c, SearchHitsManager)


# @pytest.mark.slow
# @pytest.mark.usefixtures("copy_lib_ckp")
# class TestSearch:
#     inputs = SampleList(Path("tests/data/pep").iterdir())
#     markerset_path = Path("tests/database/poxviridae_odb10/hmms")
#     cutoff_path = markerset_path.parent / "scores_cutoff"
#     markerset = HMMMarkerSet(markerset_path, cutoff_path)

#     def test_search_with_cutoff(self):
#         inputs = SampleList((Path("tests/data/pep")).iterdir())
#         orthologs = search(inputs, self.markerset, evalue=1e-50)
#         assert len(orthologs) == 18

#         hit_count = 0
#         for hits in orthologs.data.values():
#             hit_count += len(hits)
#         assert hit_count == 90

#     def test_search_with_evalue(self):
#         inputs = SampleList((Path("tests/data/pep")).iterdir())
#         markerset = HMMMarkerSet(self.markerset_path)
#         orthologs = search(inputs, markerset, evalue=1e-50)
#         assert len(orthologs) == 16

#         hit_count = 0
#         for hits in orthologs.data.values():
#             hit_count += len(hits)
#         assert hit_count == 79

#     def test_search_checkpoint_dropsample(self):
#         inputs = SampleList(Path("tests/data/pep").glob("*pox*"))
#         params = {
#             "inputs": inputs.checksum,
#             "markerset": self.markerset.checksum,
#             "markerset_cutoff": "markerset cutoff" if self.markerset.have_cutoff else 1e-10,
#             "method": "hmmalign",
#             "non_trim": False,
#         }
#         inputs, orthologs = AlignPrecheck.precheck(params, inputs)
#         orthologs = search(inputs, self.markerset, orthologs=orthologs)
#         collection = set()
#         for hits in orthologs.data.values():
#             for sample, _ in hits:
#                 collection.add(sample)
#         assert len(collection) == 4

#     def test_search_checkpoint_addsample(self, caplog: pytest.LogCaptureFixture):
#         inputs = [file for file in Path("tests/data/pep").iterdir()]
#         inputs.append(Path("tests/data/Monkeypox_virus.faa.gz"))
#         inputs = SampleList(inputs)
#         params = {
#             "inputs": inputs.checksum,
#             "markerset": self.markerset.checksum,
#             "markerset_cutoff": "markerset cutoff" if self.markerset.have_cutoff else 1e-10,
#             "method": "hmmalign",
#             "non_trim": False,
#         }
#         inputs, orthologs = AlignPrecheck.precheck(params, inputs)
#         with caplog.at_level(logging.INFO):
#             orthologs = search(inputs, self.markerset, orthologs=orthologs)
#         assert "hmmsearch on Monkeypox_virus.faa.gz is done" in caplog.text

#         collection = set()
#         for hits in orthologs.data.values():
#             for sample, _ in hits:
#                 collection.add(sample)
#         assert len(collection) == 6


# @pytest.mark.usefixtures("copy_lib_ckp")
# class TestAlign:
#     markerset_path = Path("tests/database/poxviridae_odb10/hmms")
#     cutoff_path = markerset_path.parent / "scores_cutoff"
#     markerset = HMMMarkerSet(markerset_path, cutoff_path)
#     inputs_pep = SampleList(["tests/data/pep/Cowpox_virus.faa.gz", "tests/data/pep/Goatpox_virus.faa.gz"])
#     ortho_pep = search(inputs_pep, markerset)

#     inputs_cds = SampleList(["tests/data/cds/Cowpox_virus.fna.gz", "tests/data/cds/Goatpox_virus.fna.gz"])
#     ortho_cds = search(inputs_cds, markerset)

#     @pytest.mark.parametrize("markerset, method", list(zip((markerset, None), ("hmmalign", "muscle"))))
#     def test_align_pep(self, markerset, method):
#         ortho = deepcopy(self.ortho_pep)
#         ortho.map(self.inputs_pep)
#         r = main(ortho, markerset, method=method)
#         assert len(r) == 1
#         assert type(r[0][0]) == MultipleSeqAlignment

#     @pytest.mark.parametrize("markerset, method", list(zip((markerset, None), ("hmmalign", "muscle"))))
#     def test_align_cds(self, markerset, method):
#         ortho = deepcopy(self.ortho_cds)
#         ortho.map(self.inputs_cds)
#         r = main(ortho, markerset, method=method)
#         assert len(r) == 2
#         assert type(r[0][0]) == MultipleSeqAlignment
#         assert type(r[1][0]) == MultipleSeqAlignment


# class Testtrim:
#     pep_msa = MultipleSeqAlignment(
#         [
#             SeqRecord(Seq(rec), id=id)
#             for rec, id in zip(
#                 (
#                     "MI-T-C",
#                     "M-AT-C",
#                     "MI-TG-",
#                     "M-GT-C",
#                     "M-AT--",
#                     "MI-T-C",
#                     "M-AT-C",
#                     "MI-T--",
#                     "M-G--C",
#                     "M-AT--",
#                 ),
#                 (
#                     "Species_A",
#                     "Species_B",
#                     "Species_C",
#                     "Species_D",
#                     "Species_E",
#                     "Species_F",
#                     "Species_G",
#                     "Species_H",
#                     "Species_I",
#                     "Species_J",
#                 ),
#             )
#         ]
#     )
#     cds_msa = MultipleSeqAlignment(
#         [
#             SeqRecord(Seq(rec), id=id)
#             for rec, id in zip(
#                 (
#                     "ATGATC---ACG---TGT",
#                     "ATG---GCTACT---TGT",
#                     "ATGATA---ACAGGA---",
#                     "ATG---GGTACT---TGC",
#                     "ATG---GCTACG------",
#                     "ATGATC---ACG---TGT",
#                     "ATG---GCTACT---TGT",
#                     "ATGATA---ACA------",
#                     "ATG---GGT------TGC",
#                     "ATG---GCTACG------",
#                 ),
#                 (
#                     "Species_A",
#                     "Species_B",
#                     "Species_C",
#                     "Species_D",
#                     "Species_E",
#                     "Species_F",
#                     "Species_G",
#                     "Species_H",
#                     "Species_I",
#                     "Species_J",
#                 ),
#             )
#         ]
#     )

#     def test_trim_pep(self):
#         (result_msa,) = trim([self.pep_msa])
#         assert result_msa.get_alignment_length() == 5
#         assert str(result_msa[0].seq) == "MI-TC"

#     def test_trim_cds(self):
#         (result_msa,) = trim([self.pep_msa], [self.cds_msa])
#         assert result_msa.get_alignment_length() == 15
#         assert str(result_msa[0].seq) == "ATGATC---ACGTGT"
