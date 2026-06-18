"""Tests for the download module library."""

from __future__ import annotations

import csv
import hashlib
import tarfile
from io import BytesIO
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from phyling.lib.download import BuscoParser, _fetch_url

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

FAKE_TSV = (
    "bacteria_odb10\t2020-03-06\tabc123\tlineages\n"
    "fungi_odb10\t2021-01-01\tdef456\tlineages\n"
    "archaea_odb10\t2019-05-10\t111aaa\tother\n"  # should be ignored (not lineages)
)

FAKE_LOCAL_METADATA = [
    ["bacteria_odb10", "abc123"],
    ["fungi_odb10", "oldmd5"],  # outdated
]


def _make_tar_gz(name: str = "markerset_dir") -> bytes:
    """Create an in-memory .tar.gz archive containing a single directory."""
    buf = BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        info = tarfile.TarInfo(name=name)
        info.type = tarfile.DIRTYPE
        info.mode = 0o755
        tf.addfile(info)

        # Add a dummy file inside
        content = b"dummy hmm content"
        finfo = tarfile.TarInfo(name=f"{name}/dummy.hmm")
        finfo.size = len(content)
        tf.addfile(finfo, BytesIO(content))

    return buf.getvalue()


def _write_metadata(directory: Path, rows: list[list[str]]) -> None:
    """Write a metadata TSV file to *directory*."""
    meta = directory / ".metadata"
    with open(meta, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)


@pytest.fixture()
def cfg_dir(tmp_path: Path) -> Path:
    """Return a temporary directory that acts as the primary config dir."""
    return tmp_path / "cfg"


@pytest.fixture()
def cfg_dir_with_metadata(cfg_dir: Path) -> Path:
    cfg_dir.mkdir(parents=True)
    _write_metadata(cfg_dir, FAKE_LOCAL_METADATA)
    # Create the markerset directories so the parser sees them as valid
    (cfg_dir / "bacteria_odb10").mkdir()
    (cfg_dir / "fungi_odb10").mkdir()
    return cfg_dir


# Patch _fetch_url globally for most tests so we never hit the network.
@pytest.fixture()
def mock_fetch(monkeypatch):
    """Patch _fetch_url to return fake TSV data for the metadata endpoint."""

    def _fake_fetch(url: str, timeout: int = 30) -> bytes:
        if "file_versions.tsv" in url:
            return FAKE_TSV.encode()
        # For markerset downloads return a valid tar.gz
        markerset = url.split("/")[-1].split(".")[0]
        return _make_tar_gz(markerset)

    monkeypatch.setattr("phyling.lib.download._fetch_url", _fake_fetch)
    return _fake_fetch


# ---------------------------------------------------------------------------
# _fetch_url
# ---------------------------------------------------------------------------


class TestFetchUrl:
    def test_returns_bytes_on_success(self, monkeypatch, tmp_path):
        content = b"hello world"
        mock_response = MagicMock()
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = MagicMock(return_value=False)
        mock_response.read.return_value = content

        mock_urlopen = MagicMock()
        mock_urlopen.return_value = mock_response

        monkeypatch.setattr("phyling.lib.download.urlopen", mock_urlopen)
        result = _fetch_url("http://example.com/file.tsv")

        assert result == content

    def test_raises_http_error(self, monkeypatch):
        from urllib.error import HTTPError

        mock_urlopen = MagicMock(
            side_effect=HTTPError(url="http://x.com", code=404, msg="Not Found", hdrs={}, fp=None)  # type: ignore
        )

        monkeypatch.setattr("phyling.lib.download.urlopen", mock_urlopen)

        with pytest.raises(HTTPError):
            _fetch_url("http://x.com/missing")

    def test_raises_url_error(self, monkeypatch):
        from urllib.error import URLError

        mock_urlopen = MagicMock(side_effect=URLError("Name or service not known"))
        monkeypatch.setattr("phyling.lib.download.urlopen", mock_urlopen)

        with pytest.raises(URLError):
            _fetch_url("http://unreachable.invalid/")


# ---------------------------------------------------------------------------
# BuscoParser – initialisation
# ---------------------------------------------------------------------------


class TestBuscoParserInit:
    def test_online_metadata_parsed(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        with BuscoParser(cfg_dir) as bp:
            assert "bacteria_odb10" in bp.online
            assert "fungi_odb10" in bp.online
            # 'other' type should NOT be present
            assert "archaea_odb10" not in bp.online

    def test_local_metadata_empty_when_no_file(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        with BuscoParser(cfg_dir) as bp:
            assert bp.local == []

    def test_local_metadata_loaded(self, cfg_dir_with_metadata, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            local_names = [name.replace(" [Outdated]", "") for name in bp.local]
            assert "bacteria_odb10" in local_names
            assert "fungi_odb10" in local_names

    def test_outdated_markerset_flagged(self, cfg_dir_with_metadata, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            assert any("[Outdated]" in name for name in bp.local)

    def test_up_to_date_markerset_not_flagged(self, cfg_dir_with_metadata, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            assert "bacteria_odb10" in bp.local  # md5 matches → no flag

    def test_no_network_graceful(self, cfg_dir, monkeypatch):
        """When the network is unavailable online metadata should be empty."""
        from urllib.error import URLError

        monkeypatch.setattr(
            "phyling.lib.download._fetch_url",
            lambda *a, **kw: (_ for _ in ()).throw(URLError("down")),
        )
        cfg_dir.mkdir(parents=True)
        bp = BuscoParser(cfg_dir)
        assert bp.online == []
        bp.close()


# ---------------------------------------------------------------------------
# BuscoParser.online / local properties
# ---------------------------------------------------------------------------


class TestBuscoParserProperties:
    def test_online_returns_sorted_list(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        with BuscoParser(cfg_dir) as bp:
            assert bp.online == sorted(bp.online)

    def test_local_returns_sorted_list(self, cfg_dir_with_metadata, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            stripped = [n.replace(" [Outdated]", "") for n in bp.local]
            assert stripped == sorted(stripped)


# ---------------------------------------------------------------------------
# BuscoParser.download
# ---------------------------------------------------------------------------


class TestBuscoParserDownload:
    def test_download_new_markerset(self, monkeypatch, cfg_dir):
        cfg_dir.mkdir(parents=True)
        markerset = "bacteria_odb10"
        tar_data = _make_tar_gz(markerset)
        with BuscoParser(cfg_dir) as bp:
            monkeypatch.setattr("phyling.lib.download._fetch_url", lambda *a, **kw: tar_data)
            # Use a real md5 so the check passes
            bp._online_metadata[markerset]["md5"] = hashlib.md5(tar_data).hexdigest()
            bp.download(markerset)

            assert markerset in bp._local_metadata
            assert (cfg_dir / markerset).is_dir()

    def test_download_skips_up_to_date(self, cfg_dir_with_metadata, mock_fetch, caplog):
        import logging

        with BuscoParser(cfg_dir_with_metadata) as bp:
            markerset = "bacteria_odb10"
            # Ensure md5 matches so it is considered up-to-date
            bp._online_metadata[markerset]["md5"] = "abc123"

            with caplog.at_level(logging.INFO, logger="phyling"):
                bp.download(markerset)

            assert any("up to date" in msg.lower() for msg in caplog.messages)

    def test_download_updates_outdated(self, monkeypatch, cfg_dir_with_metadata):
        cfg_dir = cfg_dir_with_metadata
        markerset = "fungi_odb10"
        tar_data = _make_tar_gz(markerset)
        with BuscoParser(cfg_dir) as bp:
            monkeypatch.setattr("phyling.lib.download._fetch_url", lambda *a, **kw: tar_data)
            new_md5 = hashlib.md5(tar_data).hexdigest()
            bp._online_metadata[markerset]["md5"] = new_md5
            bp.download(markerset)

            assert bp._local_metadata[markerset]["md5"] == new_md5

    def test_download_raises_on_unknown_markerset(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        with BuscoParser(cfg_dir) as bp:
            with pytest.raises(KeyError, match="not found in BUSCO database"):
                bp.download("nonexistent_odb10")

    def test_download_raises_on_md5_mismatch(self, monkeypatch, cfg_dir):
        cfg_dir.mkdir(parents=True)
        markerset = "bacteria_odb10"
        tar_data = _make_tar_gz(markerset)
        with BuscoParser(cfg_dir) as bp:
            monkeypatch.setattr("phyling.lib.download._fetch_url", lambda *a, **kw: tar_data)
            bp._online_metadata[markerset]["md5"] = "intentionally_wrong_md5"
            with pytest.raises(ValueError, match="MD5 mismatch"):
                bp.download(markerset)


# ---------------------------------------------------------------------------
# BuscoParser – metadata persistence
# ---------------------------------------------------------------------------


class TestBuscoParserMetadataPersistence:
    def test_metadata_written_on_change(self, monkeypatch, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        markerset = "bacteria_odb10"
        tar_data = _make_tar_gz(markerset)
        correct_md5 = hashlib.md5(tar_data).hexdigest()

        with BuscoParser(cfg_dir) as bp:
            monkeypatch.setattr("phyling.lib.download._fetch_url", lambda *a, **kw: tar_data)
            bp._online_metadata[markerset]["md5"] = correct_md5
            bp.download(markerset)

        meta_file = cfg_dir / ".metadata"
        assert meta_file.is_file()
        content = meta_file.read_text()
        assert markerset in content

    def test_metadata_not_written_when_unchanged(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        with BuscoParser(cfg_dir) as _:
            pass  # no downloads

        meta_file = cfg_dir / ".metadata"
        assert not meta_file.exists()

    def test_metadata_written_only_for_primary_dir(self, monkeypatch, tmp_path, mock_fetch):
        """Entries from a secondary (global) dir must not appear in primary metadata."""
        global_dir = tmp_path / "global"
        global_dir.mkdir()
        local_dir = tmp_path / "local"
        local_dir.mkdir()

        # Seed global dir with a markerset
        (global_dir / "bacteria_odb10").mkdir()
        _write_metadata(global_dir, [["bacteria_odb10", "abc123"]])

        markerset = "fungi_odb10"
        tar_data = _make_tar_gz(markerset)
        correct_md5 = hashlib.md5(tar_data).hexdigest()

        with BuscoParser(local_dir, global_dir) as bp:
            monkeypatch.setattr("phyling.lib.download._fetch_url", lambda *a, **kw: tar_data)
            bp._online_metadata[markerset]["md5"] = correct_md5
            bp.download(markerset)

        meta_file = local_dir / ".metadata"
        assert meta_file.is_file()
        content = meta_file.read_text()
        # Only the newly downloaded markerset should be in local metadata
        assert "fungi_odb10" in content
        assert "bacteria_odb10" not in content


# ---------------------------------------------------------------------------
# BuscoParser – context manager protocol
# ---------------------------------------------------------------------------


class TestBuscoParserContextManager:
    def test_enter_returns_self(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        bp = BuscoParser(cfg_dir)
        result = bp.__enter__()
        assert result is bp
        bp.close()

    def test_exit_returns_false(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        bp = BuscoParser(cfg_dir)
        bp.__enter__()
        assert bp.__exit__(None, None, None) is False

    def test_exceptions_propagate(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        with pytest.raises(RuntimeError, match="test error"):
            with BuscoParser(cfg_dir):
                raise RuntimeError("test error")

    def test_close_equivalent_to_exit(self, cfg_dir, mock_fetch):
        cfg_dir.mkdir(parents=True)
        bp = BuscoParser(cfg_dir)
        # Should not raise
        bp.close()


# ---------------------------------------------------------------------------
# BuscoParser – multiple config dirs
# ---------------------------------------------------------------------------


class TestBuscoParserMultipleDirs:
    def test_local_dir_overrides_global(self, tmp_path, mock_fetch):
        """Local directory metadata should take precedence over global."""
        global_dir = tmp_path / "global"
        global_dir.mkdir()
        local_dir = tmp_path / "local"
        local_dir.mkdir()

        # Both dirs have the same markerset but different md5s
        (global_dir / "bacteria_odb10").mkdir()
        _write_metadata(global_dir, [["bacteria_odb10", "global_md5"]])

        (local_dir / "bacteria_odb10").mkdir()
        _write_metadata(local_dir, [["bacteria_odb10", "local_md5"]])

        with BuscoParser(local_dir, global_dir) as bp:
            assert bp._local_metadata["bacteria_odb10"]["md5"] == "local_md5"
            assert bp._local_metadata["bacteria_odb10"]["path"] == local_dir

    def test_global_markerset_visible_when_not_in_local(self, tmp_path, mock_fetch):
        global_dir = tmp_path / "global"
        global_dir.mkdir()
        local_dir = tmp_path / "local"
        local_dir.mkdir()

        (global_dir / "bacteria_odb10").mkdir()
        _write_metadata(global_dir, [["bacteria_odb10", "abc123"]])

        with BuscoParser(local_dir, global_dir) as bp:
            local_names = [n.replace(" [Outdated]", "") for n in bp.local]
            assert "bacteria_odb10" in local_names
