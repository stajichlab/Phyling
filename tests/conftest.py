from __future__ import annotations

from pathlib import Path

import pytest

TEST_DATA_DIR = Path("tests/data")
TEST_DB_DIR = Path("tests/database")
MARKERSET_DIR = TEST_DB_DIR / "poxviridae_odb10"


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
    return sorted(tuple((TEST_DATA_DIR / "msa").glob("*.faa")))


@pytest.fixture(scope="session")
def path_cds_msa() -> list[Path]:
    return sorted(tuple((TEST_DATA_DIR / "msa").glob("*.fna")))
