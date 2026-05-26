from __future__ import annotations

import shutil
import tarfile
from pathlib import Path

import pytest

import phyling


class PackMetadata:
    def __init__(self, persistent_tmp_path: Path) -> None:
        self.cfg_dir = phyling.CFG_DIRS[0]
        self.temp_tar = persistent_tmp_path / "temp.tgz"

    def pack_metadata(self):
        with tarfile.open(self.temp_tar, "w") as tar:
            if self.cfg_dir.exists() and self.cfg_dir.is_dir():
                tar.add(self.cfg_dir, arcname=self.cfg_dir.name)

    def unpack_metadata(self):
        shutil.rmtree(self.cfg_dir)
        with tarfile.open(self.temp_tar, "r") as tar:
            tar.extractall(self.cfg_dir.parent, filter="fully_trusted")
        self.temp_tar.unlink()


@pytest.fixture(scope="session")
def persistent_tmp_path(tmp_path_factory: pytest.TempPathFactory):
    return tmp_path_factory.mktemp("persistent_tmpdir")


@pytest.fixture(scope="session", autouse=True)
def hide_metadata(persistent_tmp_path: Path):
    obj = PackMetadata(persistent_tmp_path)
    obj.pack_metadata()
    yield
    obj.unpack_metadata()


@pytest.fixture(scope="class")
def class_shared_tmpdir(tmp_path_factory: pytest.TempPathFactory):
    return tmp_path_factory.mktemp("shared_tmpdir")


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="run slow tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
