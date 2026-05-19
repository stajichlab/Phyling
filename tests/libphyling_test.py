from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyling.exception import EmptyWarning, SeqtypeError
from phyling.lib import SeqTypes
from phyling.lib.align import HMMMarkerSet, OrthologSeqs, SampleList, SampleSeqs, SearchHit, SearchHitsManager

BASE_DB = Path("tests/database/poxviridae_odb10")
HMM_FILE = BASE_DB / "hmms" / "10at10240.hmm"
HMM_DIR = BASE_DB / "hmms"
CUTOFF_FILE = BASE_DB / "scores_cutoff"


class TestHMMMarkerSet:
    hmm_file = BASE_DB / "hmms" / "10at10240.hmm"

    def test_init_with_real_single_file(self):
        """Verify loading a single real HMM file."""
        assert self.hmm_file.exists(), f"Missing test file: {self.hmm_file}"

        marker_set = HMMMarkerSet(data=self.hmm_file)

        # Assertions based on real file data
        assert len(marker_set) == 1
        assert marker_set[0].name == "10at10240"  # stem of the filename
        # Check that it's a real pyhmmer HMM object
        assert hasattr(marker_set[0], "consensus")

    def test_init_with_real_directory(self):
        """Verify loading all HMMs from the real directory."""
        assert HMM_DIR.is_dir(), f"Missing test directory: {HMM_DIR}"

        marker_set = HMMMarkerSet(data=HMM_DIR)

        # Check that we loaded multiple files
        assert len(marker_set) > 0
        # Ensure they are all named correctly
        for hmm in marker_set:
            assert hmm.name is not None

    def test_init_with_cutoff_file(self):
        """Verify that bitscore cutoffs are applied correctly from the real scores file."""
        assert CUTOFF_FILE.exists(), f"Missing cutoff file: {CUTOFF_FILE}"

        # Initialize with both data and cutoffs
        marker_set = HMMMarkerSet(data=HMM_DIR, cutoff_file=CUTOFF_FILE)

        # Check if the cutoff was actually applied to the HMM object
        # have_cutoffs returns True only if .cutoffs.trusted is set
        assert marker_set.have_cutoffs(verbose=True) is True

        # Specific bitscore check (if you know the value for 10at10240)
        # The code sets trusted as a tuple: (float, float)
        cutoff_val = marker_set["10at10240"].cutoffs.trusted
        assert isinstance(cutoff_val, tuple)
        assert len(cutoff_val) == 2

    def test_init_runtime_error_no_data(self):
        """Verify check for cutoff_file without data."""
        with pytest.raises(RuntimeError, match="Cannot specify cutoff_file without data"):
            HMMMarkerSet(data=None, cutoff_file=CUTOFF_FILE)

    def test_getitem_by_name_real(self):
        """Verify string-based retrieval using real data."""
        marker_set = HMMMarkerSet(data=self.hmm_file)

        # This should return the HMM object directly
        hmm = marker_set["10at10240"]
        assert hmm.name == "10at10240"

        with pytest.raises(KeyError):
            _ = marker_set["non_existent_marker"]

    def test_set_cutoffs_with_dict(self):
        """Test the dictionary overload of set_cutoffs."""
        marker_set = HMMMarkerSet(data=self.hmm_file)
        manual_cutoffs = {"10at10240": 450.0}

        marker_set.set_cutoffs(manual_cutoffs)

        assert marker_set["10at10240"].cutoffs.trusted == (450.0, 450.0)

    def test_set_cutoffs_file_not_found(self):
        """Verify FileNotFoundError logic."""
        marker_set = HMMMarkerSet(data=self.hmm_file)
        with pytest.raises(FileNotFoundError):
            marker_set.set_cutoffs("non_existent_file.tsv")

    def test_set_cutoffs_invalid_type(self):
        """Verify TypeError for invalid input."""
        marker_set = HMMMarkerSet(data=self.hmm_file)
        with pytest.raises(TypeError, match="Invalid cutoffs format"):
            # Passing an integer instead of path or dict
            marker_set.set_cutoffs([1, 2, 3])


# --- Fixtures ---


@pytest.fixture(scope="module")
def hmms_no_cutoff() -> HMMMarkerSet:
    """Load HMMMarkerSet without bitscore cutoffs."""
    return HMMMarkerSet(data=HMM_DIR)


@pytest.fixture(scope="module")
def hmms_with_cutoff() -> HMMMarkerSet:
    """Load HMMMarkerSet with bitscore cutoffs."""
    return HMMMarkerSet(data=HMM_DIR, cutoff_file=CUTOFF_FILE)


class TestSampleSeqs:
    # --- Paths ---
    data_dir = Path("tests/data")
    pep_fasta = data_dir / "pep" / "Monkeypox_virus.faa.gz"
    cds_fasta = data_dir / "cds" / "Monkeypox_virus.fna.gz"

    def test_init_pep(self):
        """Test initialization and lazy loading of peptide sequences."""
        sample = SampleSeqs(self.pep_fasta, name="MPXV_PEP", seqtype="pep")
        assert sample.name == "MPXV_PEP"
        assert sample.file == self.pep_fasta.absolute()

        # Test __len__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            len(sample)

        # Test __iter__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            [x for x in sample]

        sample.load()
        assert len(sample) == len([x for x in sample])
        # Verify digital sequences are amino acids (Alphabet.amino codes are usually 0-20+)
        assert sample._data.alphabet.is_amino()
        assert sample.seqtype == "pep"

    def test_init_cds(self):
        """Test that CDS sequences are correctly translated to protein on load."""
        sample = SampleSeqs(self.cds_fasta, name="MPXV_CDS", seqtype="dna")
        assert sample.name == "MPXV_CDS"
        assert sample.file == self.cds_fasta.absolute()

        # Test __len__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            len(sample)

        # Test __iter__ before loading
        with pytest.raises(RuntimeError, match="Data is not loaded yet. Please run the load method to load it first."):
            [x for x in sample]

        sample.load()
        assert len(sample) == len([x for x in sample])
        # Even though input is DNA, _data should store translated Amino Acid sequences
        assert sample._data.alphabet.is_amino()
        assert sample.seqtype == "dna"

    def test_init_auto(self):
        """Test guess seqtype from sequences."""
        sample = SampleSeqs(self.pep_fasta, name="MPXV_PEP")
        assert sample.seqtype == "pep"
        sample = SampleSeqs(self.cds_fasta, name="MPXV_CDS")
        assert sample.seqtype == "dna"

    def test_eq_lt(self):
        a = SampleSeqs(self.pep_fasta)
        b = SampleSeqs(self.pep_fasta)
        c = SampleSeqs(self.data_dir / "pep" / "bgzf" / "Anomala_cuprea_entomopoxvirus.faa.gz")
        assert a == b > c

    def test_eq_lt_seqtypeerror(self):
        a = SampleSeqs(self.pep_fasta)
        b = SampleSeqs(self.cds_fasta)
        c = SampleSeqs(self.data_dir / "pep" / "bgzf" / "Anomala_cuprea_entomopoxvirus.faa.gz")
        with pytest.raises(SeqtypeError, match="Items represent different seqtypes"):
            a == b
        with pytest.raises(SeqtypeError, match="Items represent different seqtypes"):
            assert b > c

    def test_search_mock_data(self, hmms_no_cutoff: HMMMarkerSet, monkeypatch):
        captured_args = {}

        def mock_run_hmmsearch(sample, hmms, evalue, threads):
            # Record what was passed to verify it later
            captured_args["sample"] = sample
            captured_args["hmms"] = hmms
            captured_args["evalue"] = evalue
            captured_args["threads"] = threads
            # Return a dummy list to simulate SearchHits
            return ["hit1", "hit2"]

        monkeypatch.setattr("phyling.lib.align.run_hmmsearch", mock_run_hmmsearch)

        # 4. Initialize objects
        sample = SampleSeqs(self.pep_fasta, name="test_sample", seqtype="pep")

        # 5. Execute search (with custom parameters)
        results = sample.search(hmms_no_cutoff, evalue=0.001, threads=4)

        # 6. Assertions
        assert results == ["hit1", "hit2"]
        assert captured_args["evalue"] == 0.001
        assert captured_args["threads"] == 4
        assert captured_args["sample"] == sample
        assert captured_args["hmms"] == hmms_no_cutoff

    def test_search_real_data(self, hmms_no_cutoff: HMMMarkerSet):
        """Run search using standard E-value filtering (no bitscore cutoffs)."""
        sample = SampleSeqs(self.pep_fasta, seqtype="pep")

        # run_hmmsearch should return a list of SearchHit
        hits = sample.search(hmms_no_cutoff, evalue=1e-5)

        assert isinstance(hits, list)
        assert isinstance(hits[0], SearchHit)

    def test_cds_problematic_reporting(self, caplog: pytest.LogCaptureFixture):
        """Test that the class logs warnings for invalid CDS lengths."""
        # If you have a known 'bad' CDS file, use it here.
        # Otherwise, this just ensures the process finishes.
        sample = SampleSeqs(self.data_dir / "cds" / "Monkeypox_virus_with_bad_seq.fna", seqtype="dna")
        sample.load()

        # check if logger captured the warning.
        assert "seqs have invalid length" in caplog.text


class TestSampleList:
    files = tuple(file for file in Path("tests/data/pep/bgzf").iterdir())
    names = tuple(f"name_{i}" for i in range(len(files)))

    def mock_helper(self, sample, hmms, evalue, threads):
        # Return a dummy list of "hits"
        return [f"hit_from_{sample.name}"]

    def test_sample_list_init_with_paths(self):
        """Verify initialization with a list of file paths and test __getitem__."""
        sl = SampleList(self.files)
        assert len(sl) == len(self.files)

        subset = sl[1]
        assert isinstance(subset, SampleSeqs)
        subset = sl[0:2]
        assert isinstance(subset, SampleList)
        assert len(subset) == 2

        for i in range(len(sl)):
            assert isinstance(sl[i], SampleSeqs)
            assert sl[i].file == self.files[i].absolute()
            assert sl[i].name == self.files[i].name
            assert sl[i] == sl[self.files[i].name]

        sl = SampleList(self.files, names=self.names)
        assert len(sl) == len(self.files)
        for i in range(len(sl)):
            assert sl[i].file == self.files[i].absolute()
            assert sl[i].name == self.names[i]
            assert sl[i] == sl[self.names[i]]

    def test_init_typeerror(self):
        with pytest.raises(TypeError, match="int cannot be converted to SampleSeqs"):
            SampleList([1, 2, 3])

    def test_init_names_and_data_mismatched(self):
        """Ensure error if names provided without data or with different length."""
        with pytest.raises(RuntimeError, match="Received no data with names specified."):
            SampleList(data=[], names=self.names)

        with pytest.raises(RuntimeError, match="Data and names have different length."):
            SampleList(data=self.files, names=self.names + ("extra_name",))

    def test_sample_list_init_mixed_seqtypes(self) -> None:
        """Verify Error with files with mixed seqtypes."""
        mixed_seqtypes_files = self.files + (Path("tests/data/cds/Monkeypox_virus.fna.gz"),)
        with pytest.raises(SeqtypeError, match="Items represent different seqtypes"):
            SampleList(mixed_seqtypes_files)

    # --- Search Orchestration Tests (Mocked) ---

    def test_search_sequential_logic(self, hmms_with_cutoff: HMMMarkerSet, monkeypatch, caplog: pytest.LogCaptureFixture):
        """Test sequential search mode logic."""
        sl = SampleList(self.files)

        monkeypatch.setattr("phyling.lib.align._search_helper", self.mock_helper)

        with caplog.at_level("DEBUG"):
            results = sl.search(hmms_with_cutoff, jobs=1, threads=2)

        assert "Sequential mode" in caplog.text
        assert f"Progress: {len(self.files)} / {len(self.files)}" in caplog.text
        assert f"hit_from_{self.files[0].name}" in results

    def test_search_parallel_logic(self, hmms_with_cutoff: HMMMarkerSet, monkeypatch, caplog: pytest.LogCaptureFixture):
        """Test search orchestration using monkeypatch to avoid actual pyhmmer calls."""

        # 1. Setup SampleList
        sl = SampleList(self.files)

        # Patch the helper used by 'partial' in the search method
        monkeypatch.setattr("phyling.lib.align._search_helper", self.mock_helper)

        # 3. Test Parallel Mode (jobs > 1)
        with caplog.at_level("DEBUG"):
            results = sl.search(hmms_with_cutoff, jobs=3, threads=1)

        assert len(results) == len(self.files)
        assert "Multiprocesses mode" in caplog.text
        assert f"Progress: {len(self.files)} / {len(self.files)}" in caplog.text
        assert f"hit_from_{self.files[0].name}" in results

    def test_search_invalid_jobs(self, hmms_with_cutoff: HMMMarkerSet):
        """Ensure RuntimeError when jobs are out of bounds."""
        sl = SampleList(self.files)
        with pytest.raises(RuntimeError, match="jobs should be between 1 and"):
            sl.search(hmms_with_cutoff, jobs=len(self.files) + 1)

    # --- Integration Test (Real Search) ---

    @pytest.mark.slow
    def test_search_integration_real_data(self, hmms_with_cutoff: HMMMarkerSet):
        """A real search against the poxviridae database (No mocks)."""
        # Load real samples
        sl = SampleList(self.files)

        # Perform search
        results = sl.search(hmms_with_cutoff, evalue=1e-10, jobs=1)

        assert isinstance(results, list)
        assert isinstance(results[0], SearchHit)


class MockSampleSeqs:
    def __init__(self, name: str, file: str = "mock.fasta"):
        self.name = name
        self.file = file
        self._data: tuple = ("data",)

    def unload(self):
        self._data = ()


@pytest.fixture(scope="class")
def sample_a():
    return MockSampleSeqs("Sample_A")


@pytest.fixture(scope="class")
def sample_b():
    return MockSampleSeqs("Sample_B")


@pytest.fixture(scope="class")
def sample_c():
    return MockSampleSeqs("Sample_C")


@pytest.fixture(scope="class")
def standard_hits(sample_a, sample_b):
    return [
        SearchHit(hmm="HMM_1", sample=sample_a, seqid="seq_01"),
        SearchHit(hmm="HMM_1", sample=sample_b, seqid="seq_02"),
        SearchHit(hmm="HMM_2", sample=sample_a, seqid="seq_03"),
    ]


class TestSearchHitsManager:
    def test_init_and_update_encoding(self, standard_hits):
        """Verify that initialization and batch updating correctly map string/object data to integers."""
        manager = SearchHitsManager(standard_hits)

        assert len(manager) == 3
        # Check that identical HMMs/Samples resolve to the same underlying integer ID
        assert manager._hmm_arr[0] == manager._hmm_arr[1]  # Both HMM_1
        assert manager._hmm_arr[0] != manager._hmm_arr[2]  # HMM_1 vs HMM_2

        assert manager._sample_arr[0] == manager._sample_arr[2]  # Both Sample_A
        assert manager._sample_arr[0] != manager._sample_arr[1]  # Sample_A vs Sample_B

        for i, hit in enumerate(standard_hits):
            assert manager._seqid_arr[i] == hit.seqid

    def test_getitem_single_element(self, standard_hits, sample_a, sample_b):
        """Ensure element access correctly reverses integers back to namedtuple objects on-the-fly."""
        manager = SearchHitsManager(standard_hits)

        hit_0 = manager[0]
        assert isinstance(hit_0, SearchHit)
        assert hit_0.hmm == "HMM_1"
        assert hit_0.sample is sample_a
        assert hit_0.seqid == "seq_01"

        hit_1 = manager[1]
        assert hit_1.hmm == "HMM_1"
        assert hit_1.sample is sample_b
        assert hit_1.seqid == "seq_02"

    def test_getitem_lazy_slicing(self, standard_hits):
        """Verify slicing yields a new manager sharing the master dictionary maps (lazy evaluation)."""
        manager = SearchHitsManager(standard_hits)

        sliced = manager[0:2]
        assert isinstance(sliced, SearchHitsManager)
        assert len(sliced) == 2

        # Lazy behavior check: Slicing should retain the parent dictionary mappings
        assert sliced._orthologs == manager._orthologs
        assert sliced._samples == manager._samples

    def test_pickle_serialization_loop(self, standard_hits, sample_a):
        """Verify that __getstate__ and __setstate__ preserve all state and numpy arrays during serialization."""
        # 1. Setup original manager and populate data
        original_manager = SearchHitsManager(standard_hits)

        # 2. Serialize (Pickle) the object to bytes (triggers __getstate__)
        pickled_data = pickle.dumps(original_manager)

        # 3. Deserialize (Unpickle) back to a new object (triggers __setstate__)
        unpickled_manager = pickle.loads(pickled_data)

        # 4. Assert structural integrity and object equivalence
        assert isinstance(unpickled_manager, SearchHitsManager)
        assert len(unpickled_manager) == len(original_manager)

        # Verify primary maps match exactly
        assert unpickled_manager._orthologs == original_manager._orthologs
        # Object references inside dictionaries should evaluate identically by key value
        assert list(unpickled_manager._samples.keys())[0].name == sample_a.name

        # Verify NumPy array values and dtypes are preserved exactly
        assert np.array_equal(unpickled_manager._hmm_arr, original_manager._hmm_arr)
        assert np.array_equal(unpickled_manager._sample_arr, original_manager._sample_arr)
        assert np.array_equal(unpickled_manager._seqid_arr, original_manager._seqid_arr)

        assert unpickled_manager._hmm_arr.dtype == original_manager._hmm_arr.dtype
        assert unpickled_manager._seqid_arr.dtype == original_manager._seqid_arr.dtype

    def test_pickle_empty_manager(self):
        """Ensure an unpopulated or fresh manager can be serialized without throwing errors."""
        manager = SearchHitsManager()

        pickled_data = pickle.dumps(manager)
        unpickled_manager = pickle.loads(pickled_data)

        assert len(unpickled_manager._hmm_arr) == 0
        assert unpickled_manager._orthologs == {}

    def test_compact_removes_ghosts_and_reindexes(self, standard_hits, sample_a, sample_b):
        """Test that compact() completely drops unused categories and shifts remaining keys to 0..N."""
        manager = SearchHitsManager(standard_hits)

        # Slice down to just the hit containing HMM_2 and Sample_A
        sliced = manager[2:3]
        sliced.compact()

        # Verification of purged metadata
        assert "HMM_1" not in sliced._orthologs
        assert "HMM_2" in sliced._orthologs
        assert sample_a in sliced._samples
        assert sample_b not in sliced._samples
        assert len(sliced._samples) == 1

        # Verification of re-indexed arrays (must map to contiguous 0..N sequence)
        assert sliced._hmm_arr[0] == 0
        assert sliced._sample_arr[0] == 0

        # Test on-the-fly reconstruction after compaction
        reconstructed = sliced[0]
        assert reconstructed.hmm == "HMM_2"
        assert reconstructed.sample is sample_a

    def test_filter_drop_samples(self, standard_hits, sample_a, sample_b):
        """Verify vectorized filtering correctly handles sample exclusions and drops metadata via compact."""
        manager = SearchHitsManager(standard_hits)

        filtered = manager.filter(drop_samples=["Sample_B"])

        assert len(filtered) == 2
        assert sample_b not in filtered._samples
        rev_samp = {v: k for k, v in manager._samples.items()}
        assert all(rev_samp[s_idx] is sample_a for s_idx in filtered._sample_arr)

    def test_filter_min_taxa_logic(self, sample_a, sample_b, sample_c):
        """Confirm the ortholog-aware taxa counter targets and purges markers below threshold limits."""
        hits = [
            # HMM_1 satisfies min_taxa=2 (found across 2 different samples)
            SearchHit(hmm="HMM_1", sample=sample_a, seqid="s1"),
            SearchHit(hmm="HMM_1", sample=sample_b, seqid="s2"),
            # HMM_2 fails min_taxa=2 (2 hits, but both trapped inside a single sample)
            SearchHit(hmm="HMM_2", sample=sample_c, seqid="s3"),
            SearchHit(hmm="HMM_2", sample=sample_c, seqid="s4"),
        ]
        manager = SearchHitsManager(hits)

        filtered = manager.filter(min_taxa=2)

        assert len(filtered) == 2
        assert "HMM_1" in filtered._orthologs
        assert "HMM_2" not in filtered._orthologs

    def test_filter_empty_raises_warning(self, standard_hits):
        """Ensure an EmptyWarning gets triggered when filtering strips every row out."""
        manager = SearchHitsManager(standard_hits)

        with pytest.raises(EmptyWarning):
            manager.filter(drop_samples=["Sample_A", "Sample_B"])

    def test_compact_empty_manager(self):
        """Verify that calling compact on an empty manager behaves gracefully without errors."""
        manager = SearchHitsManager()
        manager.compact()

        assert manager._orthologs == {}
        assert manager._samples == {}
        assert len(manager._hmm_arr) == 0


class MockSequence:
    def __init__(self, name: str, sequence: str, description: str = ""):
        self.name = name
        self.sequence = sequence
        self.description = description

    def translate(self):
        # Returns a mock translated object
        return MockSequence(self.name, "MOCK_AMINO_ACID_SEQ")

    def textize(self):
        return self


class MockDNASeq:
    """Mock for pyhmmer's DigitalSequence for DNA."""

    def __init__(self, name: str, sequence: str, description: str = ""):
        self.name = name
        self.sequence = sequence
        self.description = description

    def translate(self):
        # Simplistic translation mock returning a dummy object
        return f"TRANSLATED_{self.sequence}"

    def textize(self):
        class Textized:
            def __init__(self, seq):
                self.sequence = seq

        return Textized(self.sequence)


class MockDigitalSequenceBlock(list):
    """Mock for pyhmmer's DigitalSequenceBlock container."""

    pass


# -----------------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------------


@pytest.fixture(scope="class")
def dummy_faa(tmp_path_factory):
    """Creates a real FAA file in a class-scoped temporary directory."""
    # Create a unique directory under the class session
    tmp_dir = tmp_path_factory.mktemp("fasta_class_data")
    fasta_file = tmp_dir / "marker_001.faa"

    content = ">Species_A\nMITC\n>Species_B\nMATC\n>Species_C\nMITG\n"
    fasta_file.write_text(content, encoding="utf-8")
    return fasta_file


@pytest.fixture(scope="class")
def dummy_fna(tmp_path_factory):
    """Creates a real FNA file in a class-scoped temporary directory."""
    tmp_dir = tmp_path_factory.mktemp("fasta_class_data")
    fasta_file = tmp_dir / "marker_001.fna"

    content = ">Species_A\nATGATCACGTGT\n>Species_B\nATGGCTACTTGT\n>Species_C\nATGATAACAGGA\n"
    fasta_file.write_text(content, encoding="utf-8")
    return fasta_file


@pytest.fixture
def mock_load(monkeypatch):
    """Stub the base SampleSeqs.load() so it doesn't trigger true file reads."""
    monkeypatch.setattr(SampleSeqs, "load", lambda self: None)


@pytest.fixture
def mock_alignment_runners(monkeypatch):
    """Stub out the actual alignment processing functions."""
    mock_msa = MultipleSeqAlignment([])
    monkeypatch.setattr("phyling.lib.align.run_hmmalign", lambda self, hmm: mock_msa)
    monkeypatch.setattr("phyling.lib.align.run_muscle", lambda self, threads: mock_msa)
    monkeypatch.setattr("phyling.lib.align.bp_mrtrans", lambda pep_msa, cds_data: mock_msa)
    return mock_msa


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
class TestOrthologSeqs:
    def test_init_with_explicit_name(self, dummy_faa):
        """Verify that initialization preserves an explicitly supplied name mapping."""
        ortho = OrthologSeqs(file=dummy_faa, name="custom_name", seqtype="pep")
        assert ortho.file == dummy_faa
        assert ortho.name == "custom_name"
        assert ortho.seqtype == SeqTypes.PEP

    def test_init_with_auto_name(self, tmp_path):
        """Ensure that omitting a name falls back securely to the file's base stem name."""
        real_file = tmp_path / "marker_002.fa"
        real_file.touch()

        ortho = OrthologSeqs(file=real_file, seqtype="dna")
        assert ortho.file == real_file
        assert ortho.name == "marker_002.fa"
        assert ortho.seqtype == SeqTypes.DNA

    def test_load_initializes_empty_cds_attribute(self, dummy_fna, mock_load):
        """Confirm the load lifecycle properly sets up data slots before reading structures."""
        ortho = OrthologSeqs(dummy_fna)
        # Before load, slot shouldn't exist or isn't initialized
        with pytest.raises(AttributeError):
            _ = ortho._data_cds

        ortho.load()
        assert isinstance(ortho._data_cds, list)
        assert len(ortho._data_cds) == 0

    def test_process_pep_seqs(self, dummy_faa):
        """Ensure raw peptide items pass sequentially straight into the primary array data stack."""
        ortho = OrthologSeqs(dummy_faa, seqtype="pep")
        ortho._data = []

        seq_block = MockDigitalSequenceBlock(["pep_seq1", "pep_seq2"])
        ortho._process_pep_seqs(seq_block)

        assert "pep_seq1" in ortho._data
        assert "pep_seq2" in ortho._data

    def test_process_cds_seqs_successful(self, dummy_fna):
        """Verify clean CDS text reading converts data to translated peptides and stores back-mapped CDS templates."""
        ortho = OrthologSeqs(dummy_fna, seqtype="dna")
        ortho._data = []
        ortho._data_cds = []

        seq_block = MockDigitalSequenceBlock([MockDNASeq(name="gene_1", sequence="ATG", description="desc_1")])

        ortho._process_cds_seqs(seq_block)

        assert ortho._data[0] == "TRANSLATED_ATG"
        assert len(ortho._data_cds) == 1
        assert isinstance(ortho._data_cds[0], SeqRecord)
        assert ortho._data_cds[0].id == "gene_1"
        assert ortho._data_cds[0].seq == Seq("ATG")

    def test_process_cds_seqs_with_invalid_lengths(self, dummy_fna, caplog):
        """Test that translation anomalies trigger validation exceptions and get gracefully logged as warnings."""
        ortho = OrthologSeqs(dummy_fna, seqtype="dna")
        ortho._data = []
        ortho._data_cds = []

        broken_seq = MockDNASeq(name="broken_gene", sequence="AT")

        def raise_value_error():
            raise ValueError("Invalid codon length")

        broken_seq.translate = raise_value_error

        seq_block = MockDigitalSequenceBlock([broken_seq])

        ortho._process_cds_seqs(seq_block)

        assert len(ortho._data) == 0
        assert len(ortho._data_cds) == 0
        assert "have invalid length" in caplog.text

    def test_search_raises_not_implemented_error(self, dummy_faa):
        """Enforce the constraint that search operations must explicitly fail on this subclass."""
        ortho = OrthologSeqs(dummy_faa)
        with pytest.raises(NotImplementedError) as exc_info:
            ortho.search("mock_hmms")
        assert "Method is not implemented" in str(exc_info.value)

    def test_align_hmmalign_missing_argument(self, dummy_faa):
        """Ensure executing an HMM-driven alignment architecture without a reference profile fails safely."""
        ortho = OrthologSeqs(dummy_faa)
        with pytest.raises(ValueError) as exc_info:
            ortho.align(method="hmmalign", hmm=None)
        assert 'required "hmm" argument' in str(exc_info.value)

    def test_align_invalid_method_argument(self, dummy_faa):
        """Enforce string-literal verification constraints over unsupported foreign alignment backends."""
        ortho = OrthologSeqs(dummy_faa)
        with pytest.raises(ValueError) as exc_info:
            ortho.align(method="unsupported_engine")  # type: ignore
        assert "Argument method only accepts" in str(exc_info.value)

    @pytest.mark.parametrize(
        "fasta_fixture, seqtype, method",
        [
            ("dummy_faa", SeqTypes.PEP, "hmmalign"),
            ("dummy_faa", SeqTypes.PEP, "muscle"),
            ("dummy_fna", SeqTypes.DNA, "hmmalign"),
            ("dummy_fna", SeqTypes.DNA, "muscle"),
        ],
    )
    def test_align_execution_flows(self, fasta_fixture, seqtype, method, request, mock_alignment_runners):
        """Validate cross-combinations of text-types and backend align engines compute successfully."""
        dummy_fasta = request.getfixturevalue(fasta_fixture)
        ortho = OrthologSeqs(dummy_fasta)
        ortho._seqtype = seqtype
        ortho._data_cds = []

        mock_hmm = "mock_hmm_profile_object"
        res = ortho.align(method=method, hmm=mock_hmm, threads=2)

        assert isinstance(res, MultipleSeqAlignment)


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
