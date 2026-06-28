from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from pyhmmer.plan7 import HMM

TEST_DATA_DIR = Path("tests/data")
TEST_DB_DIR = Path("tests/database")
MARKERSET_DIR = TEST_DB_DIR / "poxviridae_odb10"


project_root = Path(__file__).parent.parent
database_path = project_root / TEST_DB_DIR
# Inject it directly into the process environment
# THIS MUST BE DONE BEFORE ANY PHYLING IMPORT!!!
os.environ["PHYLING_DB"] = str(database_path.resolve())


@pytest.fixture(scope="session")
def fake_online_metadata() -> Path:
    return TEST_DB_DIR / "fake_online_metadata"


@pytest.fixture(scope="session")
def fake_local_metadata() -> Path:
    return TEST_DB_DIR / ".metadata"


@pytest.fixture(scope="session")
def path_markerset() -> Path:
    return MARKERSET_DIR


@pytest.fixture(scope="session")
def path_hmm_dir() -> Path:
    return MARKERSET_DIR / "hmms"


@pytest.fixture(scope="session")
def path_cutoff_file() -> Path:
    return MARKERSET_DIR / "scores_cutoff"


@pytest.fixture(scope="session")
def path_data_dir() -> Path:
    return TEST_DATA_DIR


@pytest.fixture(scope="session")
def path_pep_fasta_dir() -> Path:
    return TEST_DATA_DIR / "pep" / "bgzf"


@pytest.fixture(scope="session")
def path_cds_fasta_dir() -> Path:
    return TEST_DATA_DIR / "cds" / "bgzf"


@pytest.fixture(scope="session")
def path_pep_fasta(path_pep_fasta_dir) -> list[Path]:
    return sorted(tuple(path_pep_fasta_dir.iterdir()))


@pytest.fixture(scope="session")
def path_cds_fasta(path_cds_fasta_dir) -> list[Path]:
    return sorted(tuple(path_cds_fasta_dir.iterdir()))


@pytest.fixture(scope="session")
def path_pep_mfa() -> list[Path]:
    return sorted(tuple((TEST_DATA_DIR / "mfa").glob("*.faa")))


@pytest.fixture(scope="session")
def path_cds_mfa() -> list[Path]:
    return sorted(tuple((TEST_DATA_DIR / "mfa").glob("*.fna")))


@pytest.fixture(scope="session")
def path_pep_msa() -> list[Path]:
    return sorted(tuple((TEST_DATA_DIR / "msa").glob("*.aa.mfa")))


@pytest.fixture(scope="session")
def path_cds_msa() -> list[Path]:
    return sorted(tuple((TEST_DATA_DIR / "msa").glob("*.cds.mfa")))


@pytest.fixture(scope="session")
def path_tree_file() -> Path:
    return TEST_DATA_DIR / "tree.nw"


# ---------------------------------------------------------------------------
# Mock fixtures
# ---------------------------------------------------------------------------
from phyling.lib.align import HMMMarkerSet


@pytest.fixture(scope="module")
def mock_hmm() -> MagicMock:
    """Generates a base HMM mock object."""
    hmm = MagicMock(spec=HMM)
    hmm.name = "HMM1"
    return hmm


@pytest.fixture(scope="module")
def mock_hmms_no_cutoff(mock_hmm) -> MagicMock:
    """HMMMarkerSet container mock holding mock_hmm with have_cutoffs = False."""
    hmms = MagicMock(spec=HMMMarkerSet)
    hmms.have_cutoffs.return_value = False
    hmms.checksums = {"HMM1": "abc123"}

    container_data = [mock_hmm]
    hmms.__iter__.return_value = iter(container_data)
    hmms.__getitem__.side_effect = lambda index: container_data[index]
    hmms.__len__.return_value = len(container_data)

    return hmms


@pytest.fixture(scope="module")
def mock_hmms_with_cutoff(mock_hmm) -> MagicMock:
    """HMMMarkerSet container mock holding mock_hmm with have_cutoffs = True."""
    hmms = MagicMock(spec=HMMMarkerSet)
    hmms.have_cutoffs.return_value = True
    hmms.checksums = {"HMM1": "abc123"}

    container_data = [mock_hmm]
    hmms.__iter__.return_value = iter(container_data)
    hmms.__getitem__.side_effect = lambda index: container_data[index]
    hmms.__len__.return_value = len(container_data)

    return hmms
