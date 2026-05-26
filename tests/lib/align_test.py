from __future__ import annotations

import logging
import pickle
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyling.exception import BinaryNotFoundError, EmptyWarning, SeqtypeError
from phyling.lib import SeqTypes, _abc
from phyling.lib.align import (
    HMMMarkerSet,
    OrthologList,
    OrthologSeqs,
    SampleList,
    SampleSeqs,
    SearchHit,
    SearchHitsManager,
    bp_mrtrans,
    run_hmmalign,
    run_hmmsearch,
    run_muscle,
)

BASE_DB = Path("tests/database/poxviridae_odb10")
HMM_DIR = BASE_DB / "hmms"
CUTOFF_FILE = BASE_DB / "scores_cutoff"

DATA_DIR = Path("tests/data")
PEP_FASTA_DIR = DATA_DIR / "pep" / "bgzf"
CDS_FASTA_DIR = DATA_DIR / "cds" / "bgzf"
PEP_FASTA = tuple(PEP_FASTA_DIR.iterdir())
CDS_FASTA = tuple(CDS_FASTA_DIR.iterdir())

PEP_MFA = tuple((DATA_DIR / "mfa").glob("*.faa"))
CDS_MFA = tuple((DATA_DIR / "mfa").glob("*.fna"))


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


@pytest.fixture(scope="module")
def hmms_no_cutoff() -> HMMMarkerSet:
    """Load HMMMarkerSet without bitscore cutoffs."""
    return HMMMarkerSet(data=HMM_DIR)


@pytest.fixture(scope="module")
def hmms_with_cutoff() -> HMMMarkerSet:
    """Load HMMMarkerSet with bitscore cutoffs."""
    return HMMMarkerSet(data=HMM_DIR, cutoff_file=CUTOFF_FILE)


class TestSampleSeqs:
    pep_fasta = PEP_FASTA[0]
    cds_fasta = CDS_FASTA[0]

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
        c = SampleSeqs(DATA_DIR / "pep" / "Monkeypox_virus.faa.gz")
        assert a == b < c

    def test_eq_lt_seqtypeerror(self):
        a = SampleSeqs(self.pep_fasta)
        b = SampleSeqs(self.cds_fasta)
        c = SampleSeqs(DATA_DIR / "pep" / "Monkeypox_virus.faa.gz")
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
        sample = SampleSeqs(DATA_DIR / "cds" / "Monkeypox_virus_with_bad_seq.fna", seqtype="dna")
        sample.load()

        # check if logger captured the warning.
        assert "seqs have invalid length" in caplog.text


class TestSampleList:
    files = PEP_FASTA
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

    # @pytest.mark.slow
    # @pytest.mark.dependency(name="search_integration")
    # def test_search_integration_real_data(self, hmms_with_cutoff: HMMMarkerSet, real_search_hits):
    #     """A real search against the poxviridae database (No mocks)."""
    #     # Load real samples
    #     sl = SampleList(self.files)

    #     # Perform search
    #     results = sl.search(hmms_with_cutoff, evalue=1e-10, jobs=1)

    #     assert isinstance(results, list)
    #     assert isinstance(results[0], SearchHit)


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


@pytest.fixture(scope="class")
def real_search_hits(hmms_with_cutoff):
    sl = SampleList(PEP_FASTA[:3])
    results = sl.search(hmms_with_cutoff, evalue=1e-10, jobs=1)
    return results


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

    # @pytest.mark.dependency(depends=["search_integration"])
    def test_load_allocates_temporary_directory(self, real_search_hits):
        """Verify that calling load provisions a physical temporary directory on disk."""
        if not real_search_hits:
            pytest.fail("Fixture real_search_hits were not properly captured.")

        manager = SearchHitsManager(real_search_hits)

        # Verify that prior to loading, no live directory tracking handles exist
        assert manager._mfa_dir is None

        # 1. Trigger the load sequence
        manager.load()

        # 2. Check that the temporary directory object exists
        assert manager._mfa_dir is not None

        # 3. Check that the physical path actually exists on the filesystem
        mfa_path = Path(manager._mfa_dir.name)
        assert mfa_path.exists()
        assert mfa_path.is_dir()

    # @pytest.mark.dependency(depends=["search_integration"])
    def test_unload_and_destruction_cleans_filesystem(self, real_search_hits):
        """Verify that unload and object garbage collection completely purge disk files."""
        manager = SearchHitsManager(real_search_hits)

        manager.load()
        assert manager._mfa_dir is not None
        mfa_path = Path(manager._mfa_dir.name)
        assert mfa_path.exists()

        # 1. Trigger explicit unload
        manager.unload()

        # 2. Verify track state reset and physical path deletion
        assert manager._mfa_dir is None
        assert not mfa_path.exists()

    # @pytest.mark.dependency(depends=["search_integration"])
    def test_destructor_cleanup_fallback(self, real_search_hits):
        """Ensure the __del__ magic method acts as a safety valve to prevent storage leaks."""
        manager = SearchHitsManager(real_search_hits)

        manager.load()
        assert manager._mfa_dir is not None
        mfa_path = Path(manager._mfa_dir.name)
        assert mfa_path.exists()

        # 1. Force object dereference to trigger __del__ garbage collection loop
        del manager

        # 2. Verify that the physical temporary directory was automatically scrubbed
        assert not mfa_path.exists()

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


@pytest.fixture
def mock_alignment_runners(monkeypatch):
    """Stub out the actual alignment processing functions."""
    mock_msa = MultipleSeqAlignment([])
    monkeypatch.setattr("phyling.lib.align.run_hmmalign", lambda self, hmm: mock_msa)
    monkeypatch.setattr("phyling.lib.align.run_muscle", lambda self, threads: mock_msa)
    monkeypatch.setattr("phyling.lib.align.bp_mrtrans", lambda pep_msa, cds_data: mock_msa)
    return mock_msa


class TestOrthologSeqs:
    pep_mfa = PEP_MFA[0]
    cds_mfa = CDS_MFA[0]

    def test_init_with_explicit_name(self):
        """Verify that initialization preserves an explicitly supplied name mapping."""
        ortho = OrthologSeqs(file=self.pep_mfa, name="custom_name", seqtype="pep")
        assert ortho.file == self.pep_mfa.absolute()
        assert ortho.name == "custom_name"
        assert ortho.seqtype == SeqTypes.PEP

    def test_init_with_auto_name(self):
        """Ensure that omitting a name falls back securely to the file name."""
        ortho = OrthologSeqs(file=self.cds_mfa, seqtype="dna")
        assert ortho.file == self.cds_mfa.absolute()
        assert ortho.name == "10at10240.fna"
        assert ortho.seqtype == SeqTypes.DNA

    def test_load_initializes_empty_cds_attribute(self):
        """Confirm the load lifecycle properly sets up data slots before reading structures."""
        ortho = OrthologSeqs(self.cds_mfa)
        # Before load, slot shouldn't exist or isn't initialized
        with pytest.raises(AttributeError):
            _ = ortho._data_cds

        ortho.load()
        assert isinstance(ortho._data_cds, list)
        assert len(ortho._data_cds) == 5

    def test_process_pep_seqs(self):
        """Ensure raw peptide items pass sequentially straight into the primary array data stack."""
        ortho = OrthologSeqs(self.pep_mfa, seqtype="pep")
        ortho.load()

        assert len(ortho._data) == 5
        assert "Goatpox_virus.faa.gz" in [x.name for x in ortho._data]

    def test_process_cds_seqs_successful(self):
        """Verify clean CDS text reading converts data to translated peptides and stores back-mapped CDS templates."""
        ortho = OrthologSeqs(self.cds_mfa, seqtype="dna")
        ortho.load()

        assert len(ortho._data) == 5
        assert "Goatpox_virus.fna.gz" in [x.name for x in ortho._data]
        assert len(ortho._data_cds) == 5
        assert isinstance(ortho._data_cds[0], SeqRecord)
        assert "Goatpox_virus.fna.gz" in [x.id for x in ortho._data_cds]

    def test_process_cds_seqs_with_invalid_lengths(self, caplog):
        """Test that translation anomalies trigger validation exceptions and get gracefully logged as warnings."""
        ortho = OrthologSeqs(self.cds_mfa, seqtype="dna")
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

    def test_search_raises_not_implemented_error(self, hmms_with_cutoff):
        """Enforce the constraint that search operations must explicitly fail on this subclass."""
        ortho = OrthologSeqs(self.pep_mfa)
        with pytest.raises(NotImplementedError) as exc_info:
            ortho.search(hmms_with_cutoff)
        assert "Method is not implemented" in str(exc_info.value)

    def test_align_hmmalign_missing_argument(self):
        """Ensure executing an HMM-driven alignment architecture without a reference profile fails safely."""
        ortho = OrthologSeqs(self.pep_mfa)
        with pytest.raises(ValueError) as exc_info:
            ortho.align(method="hmmalign", hmm=None)
        assert 'required "hmm" argument' in str(exc_info.value)

    def test_align_invalid_method_argument(self):
        """Enforce string-literal verification constraints over unsupported foreign alignment backends."""
        ortho = OrthologSeqs(self.pep_mfa)
        with pytest.raises(ValueError) as exc_info:
            ortho.align(method="unsupported_engine")  # type: ignore
        assert "Argument method only accepts" in str(exc_info.value)

    @pytest.mark.parametrize(
        "fasta_attr_name, seqtype, method",
        [
            ("pep_mfa", SeqTypes.PEP, "hmmalign"),
            ("pep_mfa", SeqTypes.PEP, "muscle"),
            ("cds_mfa", SeqTypes.DNA, "hmmalign"),
            ("cds_mfa", SeqTypes.DNA, "muscle"),
        ],
    )
    def test_align_execution_flows(self, fasta_attr_name, seqtype, method, hmms_with_cutoff, mock_alignment_runners):
        """Validate cross-combinations of text-types and backend align engines compute successfully."""
        mfa = getattr(self, fasta_attr_name)
        ortho = OrthologSeqs(mfa)
        ortho._seqtype = seqtype
        ortho._data_cds = []

        res = ortho.align(method=method, hmm=hmms_with_cutoff, threads=2)

        assert isinstance(res, MultipleSeqAlignment)


# We mock the global alignment helper function that `partial` targets
def mock_align_helper(ortholog_seq, hmms, method, threads):
    """Stub function that mimics individual sequence alignments."""
    return MultipleSeqAlignment([])


@pytest.fixture
def mock_align_dependencies(monkeypatch):
    """Patches base class methods and the alignment worker engine."""
    # Prevent base class from trying to read or unpack actual files during init/getitem
    monkeypatch.setattr(_abc.SeqDataListABC, "__init__", lambda self, data, names, seqtype: None)

    # Define a clean mock list behaviors so len() and iteration work on the ABC wrapper
    monkeypatch.setattr(_abc.SeqDataListABC, "__len__", lambda self: 2)

    # Mock individual slice retrieval return templates
    mock_seq_instance = OrthologSeqs.__new__(OrthologSeqs)
    monkeypatch.setattr(_abc.SeqDataListABC, "__getitem__", lambda self, key: mock_seq_instance)

    # Patch the private target alignment worker function used by partial
    monkeypatch.setattr("phyling.lib.align._align_helper", mock_align_helper)

    return mock_seq_instance


class TestOrthologList:
    pep_mfa = PEP_MFA
    cds_mfa = CDS_MFA

    def test_init_bound_class_assignment(self):
        """Verify that initialization sets the correct internal target class type."""
        # Partially patch init to let hasattr check fire normally
        ortho_list = OrthologList(data=self.pep_mfa, seqtype="pep")
        assert ortho_list._bound_class is OrthologSeqs

    def test_getitem_passthrough(self):
        """Ensure that indexing or fetching keys safely proxies up to the base ABC container."""
        ortho_list = OrthologList(self.pep_mfa)

        # Accessing via index or text label should yield an OrthologSeqs interface instance
        result = ortho_list[0]
        assert isinstance(result, OrthologSeqs)

    @pytest.mark.parametrize(
        "jobs",
        [0, 3],
    )
    def test_align_jobs_out_of_bounds_raises_error(self, jobs):
        """Enforce range validation bounds on process task distribution counters."""
        ortho_list = OrthologList(self.pep_mfa)

        # Length is mocked to 2. Setting jobs to 0 or 3 should trigger a validation crash
        with pytest.raises(RuntimeError) as exc_info:
            ortho_list.align(method="muscle", jobs=jobs)
        assert "jobs should be between 1 and 2" in str(exc_info.value)

    def test_align_sequential_mode_flow(self, hmms_with_cutoff, monkeypatch, caplog):
        """Verify that processing runs sequentially when jobs=1 and updates logs."""
        ortho_list = OrthologList(self.pep_mfa)

        monkeypatch.setattr("phyling.lib.align._align_helper", mock_align_helper)

        with caplog.at_level(logging.DEBUG):
            alignments = ortho_list.align(method="hmmalign", hmms=hmms_with_cutoff, jobs=1, threads=2)

        assert len(alignments) == 2
        assert isinstance(alignments[0], MultipleSeqAlignment)
        assert "Sequential mode with 2 threads." in caplog.text
        assert "Progress: 2 / 2" in caplog.text

    def test_align_parallel_threadpool_mode_flow(self, monkeypatch, caplog):
        """Confirm multiprocess threading pool contexts instantiate and process cleanly."""
        ortho_list = OrthologList(self.pep_mfa)

        monkeypatch.setattr("phyling.lib.align._align_helper", mock_align_helper)

        # Stub the multi-worker imap to execute sequentially under test to isolate pool calls
        def mock_imap(self, func, iterable):
            return map(func, iterable)

        monkeypatch.setattr(ThreadPool, "imap", mock_imap)

        with caplog.at_level(logging.DEBUG):
            alignments = ortho_list.align(method="muscle", jobs=2, threads=4)

        assert len(alignments) == 2
        assert "Multiprocesses mode with 2 jobs and 4 threads" in caplog.text
        assert "Progress: 2 / 2" in caplog.text

    def test_align_exception_handling_logs_and_re_raises(self, monkeypatch, caplog):
        """Verify that pipeline execution failures register logs correctly before throwing exceptions."""
        ortho_list = OrthologList(self.pep_mfa)

        # Force the worker call to throw an unexpected engine failure
        def crashing_helper(*args, **kwargs):
            raise ValueError("MUSCLE binary execution failure")

        monkeypatch.setattr("phyling.lib.align._align_helper", crashing_helper)

        with pytest.raises(ValueError, match="MUSCLE binary execution failure"):
            with caplog.at_level(logging.ERROR):
                ortho_list.align(method="muscle", jobs=1)

        # Confirm traceback structural details captured inside error logger
        assert "Traceback" in caplog.text


class MockHmmsearchReturn:
    """A helper object to capture call context and control mock pyhmmer hits."""

    def __init__(self):
        self.captured_kwargs = {}
        self.return_hits = []

    def __call__(self, hmms, sample, cpus, **kwargs):
        self.captured_kwargs = kwargs
        return self.return_hits


@pytest.fixture
def mock_hmmsearch(monkeypatch):
    """Fixture that patches hmmsearch and returns a controller to verify calls/return data."""
    f = MockHmmsearchReturn()
    monkeypatch.setattr("phyling.lib.align.hmmsearch", f)
    return f


class MockQuery:
    def __init__(self, name: str):
        self.name = name


class MockHit:
    def __init__(self, name: str, reported: bool = True):
        self.name = name
        self.reported = reported


class MockTopHits(list):
    """Mocks pyhmmer.plan7.TopHits."""

    def __init__(self, query_name: str, hits_list: list[MockHit]):
        super().__init__(hits_list)
        self.query = MockQuery(query_name)


class TestRunHmmsearch:
    def test_run_hmmsearch_uses_evalue_when_no_cutoffs(self, sample_a, hmms_no_cutoff, mock_hmmsearch):
        """Verify that E-value thresholds are applied when trusted cutoffs are missing."""
        mock_hmmsearch.return_hits = []

        results = run_hmmsearch(sample_a, hmms_no_cutoff, evalue=1e-5, threads=2)

        assert results == []
        assert mock_hmmsearch.captured_kwargs["E"] == 1e-5
        assert "bit_cutoffs" not in mock_hmmsearch.captured_kwargs

    def test_run_hmmsearch_uses_trusted_cutoffs(self, sample_a, hmms_with_cutoff, mock_hmmsearch):
        """Ensure bit_cutoffs flag shifts to 'trusted' when reference cutoffs are present."""
        mock_hmmsearch.return_hits = []

        _ = run_hmmsearch(sample_a, hmms_with_cutoff, threads=4)

        assert mock_hmmsearch.captured_kwargs["bit_cutoffs"] == "trusted"
        assert "E" not in mock_hmmsearch.captured_kwargs

    def test_run_hmmsearch_skips_unreported_hits(self, sample_a, hmms_no_cutoff, monkeypatch):
        """Ensure that hits that fail the inclusion threshold ('reported=False') are dropped."""
        mock_tophits = MockTopHits(query_name="HMM1", hits_list=[MockHit(name="seq_1", reported=False)])
        monkeypatch.setattr("phyling.lib.align.hmmsearch", lambda *args, **kwargs: [mock_tophits])

        results = run_hmmsearch(sample_a, hmms_no_cutoff)
        assert len(results) == 0

    def test_run_hmmsearch_skips_multi_reported_hits(self, sample_a, hmms_no_cutoff, monkeypatch):
        """Confirm that multiple reported hits for a single profile discard the entire marker sequence block."""
        # Mock a marker matching two separate loci on the genome (paralogs/non-specific hits)
        mock_tophits = MockTopHits(
            query_name="HMM1",
            hits_list=[
                MockHit(name="seq_1", reported=True),
                MockHit(name="seq_2", reported=True),
            ],
        )
        monkeypatch.setattr("phyling.lib.align.hmmsearch", lambda *args, **kwargs: [mock_tophits])

        results = run_hmmsearch(sample_a, hmms_no_cutoff)

        # Due to 'if len(reported) > 1: continue', this entire marker group should be dropped
        assert len(results) == 0

    def test_run_hmmsearch_logs_completion(self, sample_a, hmms_no_cutoff, monkeypatch, caplog):
        """Verify that execution completions trigger expected diagnostic logging messages."""
        monkeypatch.setattr("phyling.lib.align.hmmsearch", lambda *args, **kwargs: [])

        with caplog.at_level(logging.DEBUG):
            _ = run_hmmsearch(sample_a, hmms_no_cutoff)

        assert f"Hmmsearch on {sample_a.name} is done." in caplog.text


@pytest.fixture(scope="module")
def get_orthologseqs():
    """Intercepts load_msa to return a genuine BioPython MSA without reading disk."""
    pep_mfa = PEP_MFA[0]
    ortho = OrthologSeqs(pep_mfa)
    ortho.load()
    return ortho


class MockHmmalignReturn:
    """Spies on hmmalign calls and intercepts data flow to simulate file writing."""

    def __init__(self):
        self.captured_hmm = None
        self.captured_ortholog = None
        self.alignment_to_write = ">SeqA\nMAYA\n>SeqB\nM-YA\n"

    def __call__(self, hmm, ortholog):
        self.captured_hmm = hmm
        self.captured_ortholog = ortholog
        return self  # Return self so it can chain with the .write() call

    def write(self, target_file, file_format):
        """Simulates writing the pyhmmer alignment to the active temp file handle."""
        assert file_format == "afa"
        # Since NamedTemporaryFile is open in binary mode ('w+b'), write encoded bytes
        target_file.write(self.alignment_to_write.encode("utf-8"))


@pytest.fixture
def mock_hmmalign(monkeypatch):
    """Fixture that intercepts and stubs the global hmmalign execution block."""
    f = MockHmmalignReturn()
    monkeypatch.setattr("phyling.lib.align.hmmalign", f)
    return f


class TestRunHmmalign:
    def test_run_hmmalign_execution_flow(self, get_orthologseqs, hmms_no_cutoff, mock_hmmalign):
        """Verify the full execution matrix: invocation, temporary file consumption, and tagging."""
        # 1. Instantiate light stub versions of your pipeline inputs
        mock_hmm = hmms_no_cutoff[0]

        # 2. Execute the function
        result = run_hmmalign(get_orthologseqs, mock_hmm)

        # 3. Assertions checking argument pass-through to pyhmmer
        assert mock_hmmalign.captured_hmm is mock_hmm
        assert mock_hmmalign.captured_ortholog is get_orthologseqs

        # 4. Assertions checking BioPython annotation metadata mutations
        assert isinstance(result, MultipleSeqAlignment)
        assert result.annotations == {"seqtype": SeqTypes.PEP}

    def test_run_hmmalign_logs_on_completion(self, get_orthologseqs, hmms_no_cutoff, mock_hmmalign, caplog):
        """Ensure successful task runs trigger correct diagnostic logger events."""
        mock_hmm = hmms_no_cutoff[0]

        with caplog.at_level(logging.DEBUG):
            _ = run_hmmalign(get_orthologseqs, mock_hmm)

        assert f"Hmmalign on {get_orthologseqs.name} is done." in caplog.text

    def test_run_hmmalign_cleans_up_file_on_exception(self, get_orthologseqs, hmms_no_cutoff, monkeypatch):
        """Guarantee that NamedTemporaryFile descriptors close and drop even if the pipeline crashes."""
        mock_hmm = hmms_no_cutoff[0]

        # Force a terminal crash directly inside the core call block
        def broken_hmmalign(*args, **kwargs):
            raise RuntimeError("pyhmmer segfault or memory failure")

        monkeypatch.setattr("phyling.lib.align.hmmalign", broken_hmmalign)

        # Intercept NamedTemporaryFile to spy on its life cycle
        import tempfile

        original_temp_file = tempfile.NamedTemporaryFile
        captured_file_handle = None

        def mock_tempfile(*args, **kwargs):
            nonlocal captured_file_handle
            captured_file_handle = original_temp_file(*args, **kwargs)
            return captured_file_handle

        monkeypatch.setattr(tempfile, "NamedTemporaryFile", mock_tempfile)

        # Verify the crash bubbles up properly to your process context
        with pytest.raises(RuntimeError, match="pyhmmer segfault"):
            run_hmmalign(get_orthologseqs, mock_hmm)

        # Assert that the finally block executed, executing an unconditional cleanup path
        assert captured_file_handle is not None
        # In Python's tempfile module, a closed file handle updates its internal close state flag
        assert captured_file_handle.closed


class MockMuscleRunner:
    """Spy/Stub container for mocking the underlying binary wrapper process."""

    def __init__(self, infile, outfile, threads):
        self.infile = infile
        self.outfile = outfile
        self.threads = threads
        self.run_called = False

    def run(self):
        self.run_called = True


@pytest.fixture
def mock_muscle_runner(monkeypatch):
    """Intercepts and records instantiations/executions of the Muscle wrapper."""
    captured_instances = []

    from phyling.lib.align import Muscle

    def mock_init(self, infile, outfile, threads):
        runner = MockMuscleRunner(infile, outfile, threads)
        captured_instances.append(runner)
        self._mock_impl = runner

    def mock_run(self):
        self._mock_impl.run()

    monkeypatch.setattr(Muscle, "__init__", mock_init)
    # Using a secondary patch block so that `runner.run()` references our logic hook safely
    monkeypatch.setattr(Muscle, "run", mock_run)

    return captured_instances


@pytest.fixture
def mock_load_msa(monkeypatch):
    """Stub out load_msa to bypass physical file reading and provide a baseline MSA."""
    dummy_msa = MultipleSeqAlignment([])
    monkeypatch.setattr("phyling.lib.align.load_msa", lambda path: dummy_msa)
    return dummy_msa


class TestRunMuscle:
    def test_run_muscle_successful_flow(self, get_orthologseqs, mock_muscle_runner, mock_load_msa):
        """Verify standard successful execution flow, parameter mapping, and typing assertions."""
        # 1. Fire execution path with custom thread requirements
        result = run_muscle(get_orthologseqs, threads=4)

        # 2. Check that the binary wrapper constructor was safely hit
        assert len(mock_muscle_runner) == 1
        runner_spy = mock_muscle_runner[0]

        assert runner_spy.threads == 4
        assert runner_spy.run_called is True

        # 3. Verify BioPython wrapper metadata annotations were pinned cleanly
        assert isinstance(result, MultipleSeqAlignment)
        assert result.annotations == {"seqtype": SeqTypes.PEP}

    def test_run_muscle_logs_on_completion(self, get_orthologseqs, mock_muscle_runner, caplog, mock_load_msa):
        """Ensure completion triggers accurate downstream debugging trace alerts."""
        with caplog.at_level(logging.DEBUG):
            _ = run_muscle(get_orthologseqs, threads=1)

        assert f"Muscle on {get_orthologseqs.name} is done." in caplog.text

    def test_run_muscle_exception_handling_and_re_raise(self, get_orthologseqs, monkeypatch):
        """Confirm errors encountered at runtime wrap the crash trace in a localized error explanation."""

        # Force muscle execution pass to throw a process terminal fault or syntax break
        def crashing_run(self):
            raise OSError("Access violation or binary crash")

        monkeypatch.setattr("phyling.lib.align.Muscle.run", crashing_run)

        # Enforce that the runner bubbles up into an informative RuntimeError container wrap
        with pytest.raises(RuntimeError) as exc_info:
            run_muscle(get_orthologseqs, threads=2)

        assert f"Muscle failed on peptide fasta translated from {get_orthologseqs.file}" in str(exc_info.value)
        assert "Access violation or binary crash" in str(exc_info.value)

    def test_run_muscle_handles_missing_binary_exception(self, get_orthologseqs, monkeypatch):
        """Enforce that missing installation binaries fail elegantly before creating runtime corruptions."""

        def missing_binary_init(self, *args, **kwargs):
            raise BinaryNotFoundError("The 'muscle' executable could not be resolved on PATH environmental vectors.")

        monkeypatch.setattr("phyling.lib.align.Muscle", missing_binary_init)

        with pytest.raises(BinaryNotFoundError, match="could not be resolved on PATH"):
            run_muscle(get_orthologseqs, threads=1)


class TestBpMrtrans:
    def test_bp_mrtrans_successful_translation(self):
        """Verify back-translation maps valid codons 1:1 against non-gapped alignment templates."""
        # Protein: M   A   Y   A
        pep_rec = SeqRecord(Seq("MAYA"), id="gene_01", name="n1", description="desc1")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # DNA matching exactly 4 codons: ATG (M), GCT (A), TAC (Y), GCA (A)
        cds_rec = SeqRecord(Seq("ATGGCCTACGCA"), id="gene_01", name="n1", description="desc1")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)

        assert isinstance(result, MultipleSeqAlignment)
        assert result.annotations == {"seqtype": SeqTypes.DNA}
        assert len(result) == 1
        assert str(result[0].seq) == "ATGGCCTACGCA"
        assert result[0].id == "gene_01"

    def test_bp_mrtrans_with_gaps_in_protein(self):
        """Ensure protein gaps ('-') project 3 nucleotide gaps into the cDNA sequence."""
        # Protein alignment with a deletion insertion gap: M  -  Y
        pep_rec = SeqRecord(Seq("M-Y"), id="gene_02")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # The underlying coding sequence contains only the real structural coordinates: ATG (M), TAC (Y)
        cds_rec = SeqRecord(Seq("ATGTAC"), id="gene_02")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)

        # Expected out: ATG (from M), --- (from gap), TAC (from Y)
        assert str(result[0].seq) == "ATG---TAC"

    def test_bp_mrtrans_handles_trailing_stop_codons(self):
        """Verify that stop codons encountered concurrently with alignment boundaries parse correctly."""
        # Protein sequence layout: M  A  Y
        pep_rec = SeqRecord(Seq("MAY"), id="gene_03")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # DNA sequence ending in a valid stop codon (TAA)
        cds_rec = SeqRecord(Seq("ATGGCCTACTAA"), id="gene_03")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)

        # The function handles a trailing stop codon matching a boundary by expanding it out if present
        assert str(result[0].seq) == "ATGGCCTAC"

    def test_bp_mrtrans_metadata_fallback_logic(self):
        """Confirm that if SeqRecord properties are missing from the CDS, they fall back to the peptide's."""
        pep_rec = SeqRecord(Seq("M"), id="pep_id", name="pep_name", description="pep_desc")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # Construct a bare record missing standard ID attributes
        class MinimalSeqRecord:
            def __init__(self, seq):
                self.seq = seq

        cds_rec = MinimalSeqRecord(Seq("ATG"))
        cds_recs = [cds_rec]  # type: ignore

        result = bp_mrtrans(pep_msa, cds_recs)

        # Verify fallback properties mapped successfully
        assert result[0].id == "pep_id"
        assert result[0].name == "pep_name"
        assert result[0].description == "pep_desc"

    def test_bp_mrtrans_case_insensitivity_normalization(self):
        """Ensure lower-case coding letters normalize seamlessly into upper-case alignments."""
        pep_rec = SeqRecord(Seq("M"), id="gene_04")
        pep_msa = MultipleSeqAlignment([pep_rec])

        # Lower-case string sequence representation
        cds_rec = SeqRecord(Seq("atg"), id="gene_04")
        cds_recs = [cds_rec]

        result = bp_mrtrans(pep_msa, cds_recs)
        assert str(result[0].seq) == "ATG"
