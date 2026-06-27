"""Tests for the download module library."""

from __future__ import annotations

import shutil
import tarfile
from io import BytesIO
from pathlib import Path
from unittest.mock import MagicMock

import pytest

import phyling.lib.download as download_mod
from phyling.lib.download import BuscoParser, _fetch_url


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


@pytest.fixture
def cfg_dir(tmp_path: Path) -> Path:
    """Return a temporary directory that acts as the primary config dir."""
    p = tmp_path / "global_config"
    p.mkdir()
    return p


@pytest.fixture
def cfg_dir_with_metadata(fake_local_metadata: Path, cfg_dir: Path) -> Path:
    # Create the markerset directories so the parser sees them as valid
    shutil.copy(fake_local_metadata, cfg_dir / ".metadata")
    for line in (cfg_dir / ".metadata").read_text(encoding="utf-8").splitlines():
        dataset_name = line.split("\t")[0]
        (cfg_dir / dataset_name).mkdir()
    return cfg_dir


@pytest.fixture
def mock_fetch(fake_online_metadata: Path, monkeypatch: pytest.MonkeyPatch):
    """Patch _fetch_url to return fake TSV data for the metadata endpoint."""

    def _fake_fetch(url: str, timeout: int = 30) -> bytes:
        if "file_versions.tsv" in url:
            return fake_online_metadata.read_bytes()
        # For markerset downloads return a valid tar.gz
        markerset = url.split("/")[-1].split(".")[0]
        return _make_tar_gz(markerset)

    monkeypatch.setattr("phyling.lib.download._fetch_url", _fake_fetch)
    return _fake_fetch


@pytest.fixture
def mock_hashlib_gen():

    def gen(fake_hash):
        mock_hash = MagicMock()
        mock_hash.md5.return_value.hexdigest.return_value = fake_hash
        return mock_hash

    return gen


# ---------------------------------------------------------------------------
# _fetch_url
# ---------------------------------------------------------------------------


class TestFetchUrl:
    def test_returns_bytes_on_success(self, monkeypatch: pytest.MonkeyPatch):
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

    def test_raises_http_error(self, monkeypatch: pytest.MonkeyPatch):
        from urllib.error import HTTPError

        mock_urlopen = MagicMock(
            side_effect=HTTPError(url="http://x.com", code=404, msg="Not Found", hdrs={}, fp=None)  # type: ignore
        )

        monkeypatch.setattr("phyling.lib.download.urlopen", mock_urlopen)

        with pytest.raises(HTTPError):
            _fetch_url("http://x.com/missing")

    def test_raises_url_error(self, monkeypatch: pytest.MonkeyPatch):
        from urllib.error import URLError

        mock_urlopen = MagicMock(side_effect=URLError("Name or service not known"))
        monkeypatch.setattr("phyling.lib.download.urlopen", mock_urlopen)

        with pytest.raises(URLError):
            _fetch_url("http://unreachable.invalid/")


# ---------------------------------------------------------------------------
# BuscoParser
# ---------------------------------------------------------------------------


class TestBuscoParserInit:
    def test_online_metadata_parsed(self, cfg_dir: Path, mock_fetch):
        with BuscoParser(cfg_dir) as bp:
            assert "bacteria_odb10" in bp.online
            assert "bclasvirinae_odb10" in bp.online
            # 'other' type should NOT be present
            assert "archaea_odb10" not in bp.online

    def test_local_metadata_empty_when_no_file(self, cfg_dir: Path, mock_fetch):
        with BuscoParser(cfg_dir) as bp:
            assert bp.local == []

    def test_local_metadata_loaded(self, cfg_dir_with_metadata: Path, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            local_names = [name.replace(" [Outdated]", "") for name in bp.local]
            assert "poxviridae_odb10" in local_names
            assert "bclasvirinae_odb10" in local_names

    def test_up_to_date_markerset(self, cfg_dir_with_metadata: Path, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            assert "poxviridae_odb10" in bp.local

    def test_outdated_markerset(self, cfg_dir_with_metadata: Path, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            assert "bclasvirinae_odb10 [Outdated]" in bp.local

    def test_no_network_graceful(self, cfg_dir: Path, monkeypatch: pytest.MonkeyPatch):
        """When the network is unavailable online metadata should be empty."""
        from urllib.error import URLError

        monkeypatch.setattr(
            "phyling.lib.download._fetch_url",
            lambda *a, **kw: (_ for _ in ()).throw(URLError("down")),
        )
        bp = BuscoParser(cfg_dir)
        assert bp.online == []
        bp.close()


class TestBuscoParserProperties:
    def test_online_returns_sorted_list(self, cfg_dir: Path, mock_fetch):
        with BuscoParser(cfg_dir) as bp:
            assert bp.online == sorted(bp.online)

    def test_local_returns_sorted_list(self, cfg_dir_with_metadata: Path, mock_fetch):
        with BuscoParser(cfg_dir_with_metadata) as bp:
            stripped = [n.replace(" [Outdated]", "") for n in bp.local]
            assert stripped == sorted(stripped)


class TestBuscoParserDownload:
    def test_download_new_markerset(self, monkeypatch: pytest.MonkeyPatch, mock_hashlib_gen, mock_fetch, cfg_dir: Path):
        markerset = "bacteria_odb10"
        with BuscoParser(cfg_dir) as bp:
            mock_hash = mock_hashlib_gen(bp._online_metadata[markerset]["md5"])
            monkeypatch.setattr(download_mod, "hashlib", mock_hash)
            bp.download(markerset)

            assert markerset in bp._local_metadata
            assert (cfg_dir / markerset).is_dir()

    def test_download_skips_up_to_date(self, cfg_dir_with_metadata: Path, mock_fetch, caplog: pytest.LogCaptureFixture):
        import logging

        with BuscoParser(cfg_dir_with_metadata) as bp:
            markerset = "poxviridae_odb10"

            with caplog.at_level(logging.INFO, logger="phyling"):
                bp.download(markerset)

            assert any("up to date" in msg.lower() for msg in caplog.messages)

    def test_download_updates_outdated(
        self, monkeypatch: pytest.MonkeyPatch, mock_hashlib_gen, mock_fetch, cfg_dir_with_metadata: Path
    ):
        markerset = "bclasvirinae_odb10"
        with BuscoParser(cfg_dir_with_metadata) as bp:
            mock_hash = mock_hashlib_gen(bp._online_metadata[markerset]["md5"])
            monkeypatch.setattr(download_mod, "hashlib", mock_hash)
            bp.download(markerset)

            assert bp._local_metadata[markerset]["md5"] == bp._online_metadata[markerset]["md5"]

    def test_download_raises_on_unknown_markerset(self, cfg_dir, mock_fetch):
        with BuscoParser(cfg_dir) as bp:
            with pytest.raises(KeyError, match="not found in BUSCO database"):
                bp.download("nonexistent_odb10")

    def test_download_raises_on_md5_mismatch(self, monkeypatch: pytest.MonkeyPatch, mock_hashlib_gen, mock_fetch, cfg_dir: Path):
        markerset = "bacteria_odb10"
        with BuscoParser(cfg_dir) as bp:
            mock_hash = mock_hashlib_gen("intentionally_wrong_md5")
            monkeypatch.setattr(download_mod, "hashlib", mock_hash)
            with pytest.raises(ValueError, match="MD5 mismatch"):
                bp.download(markerset)


class TestBuscoParserMetadataPersistence:
    def test_metadata_written_on_change(self, monkeypatch: pytest.MonkeyPatch, mock_hashlib_gen, cfg_dir: Path, mock_fetch):
        markerset = "bacteria_odb10"
        hash = "bff9523df89980699c01e8f31e490496"
        mock_hash = mock_hashlib_gen(hash)
        monkeypatch.setattr(download_mod, "hashlib", mock_hash)

        with BuscoParser(cfg_dir) as bp:
            bp.download(markerset)

        meta_file = cfg_dir / ".metadata"
        assert meta_file.is_file()
        content = meta_file.read_text()
        assert markerset in content

    def test_metadata_not_written_when_unchanged(self, cfg_dir: Path, mock_fetch):
        with BuscoParser(cfg_dir) as _:
            pass  # no downloads

        meta_file = cfg_dir / ".metadata"
        assert not meta_file.exists()


class TestBuscoParserContextManager:
    def test_enter_returns_self(self, cfg_dir: Path, mock_fetch):
        bp = BuscoParser(cfg_dir)
        result = bp.__enter__()
        assert result is bp
        bp.close()

    def test_exit_returns_false(self, cfg_dir: Path, mock_fetch):
        bp = BuscoParser(cfg_dir)
        bp.__enter__()
        assert bp.__exit__(None, None, None) is False

    def test_exceptions_propagate(self, cfg_dir: Path, mock_fetch):
        with pytest.raises(RuntimeError, match="test error"):
            with BuscoParser(cfg_dir):
                raise RuntimeError("test error")

    def test_close_equivalent_to_exit(self, cfg_dir: Path, mock_fetch):
        bp = BuscoParser(cfg_dir)
        # Should not raise
        bp.close()


class TestBuscoParserMultipleDirs:
    def test_local_dir_overrides_global(self, tmp_path, mock_fetch, cfg_dir_with_metadata: Path):
        """Local directory metadata should take precedence over global."""
        upstream_cfg_dir = cfg_dir_with_metadata
        user_cfg_dir = tmp_path / "local_config"
        shutil.copytree(upstream_cfg_dir, user_cfg_dir)
        upstream_cfg_meta = upstream_cfg_dir / ".metadata"
        user_cfg_meta = user_cfg_dir / ".metadata"
        assert upstream_cfg_meta.is_file()
        assert user_cfg_meta.is_file()

        # Both dirs have the same markerset but different md5s
        markerset = "bacteria_odb10"
        (upstream_cfg_dir / markerset).mkdir(exist_ok=True)
        with open(upstream_cfg_meta, "a", encoding="utf-8") as f:
            f.write(f"{markerset}\tglobal_md5\n")

        (user_cfg_dir / markerset).mkdir(exist_ok=True)
        with open(user_cfg_meta, "a", encoding="utf-8") as f:
            f.write(f"{markerset}\tlocal_md5\n")

        with BuscoParser(user_cfg_dir, upstream_cfg_dir) as bp:
            assert bp._local_metadata["bacteria_odb10"]["md5"] == "local_md5"
            assert bp._local_metadata["bacteria_odb10"]["path"] == user_cfg_dir

    def test_global_markerset_visible_when_not_in_local(self, tmp_path, mock_fetch, cfg_dir_with_metadata: Path):
        upstream_cfg_dir = cfg_dir_with_metadata
        user_cfg_dir = tmp_path / "local_config"
        shutil.copytree(upstream_cfg_dir, user_cfg_dir)
        upstream_cfg_meta = upstream_cfg_dir / ".metadata"
        user_cfg_meta = user_cfg_dir / ".metadata"
        assert upstream_cfg_meta.is_file()
        assert user_cfg_meta.is_file()

        # Add a new markerset to upstream_cfg_dir
        markerset = "bacteria_odb10"
        hash = "bff9523df89980699c01e8f31e490496"
        (upstream_cfg_dir / markerset).mkdir(exist_ok=True)
        with open(upstream_cfg_meta, "a", encoding="utf-8") as f:
            f.write(f"{markerset}\t{hash}\n")

        with BuscoParser(user_cfg_dir, upstream_cfg_dir) as bp:
            local_names = [n.replace(" [Outdated]", "") for n in bp.local]
            assert "bacteria_odb10" in local_names

    def test_metadata_written_only_for_local(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path, mock_fetch, mock_hashlib_gen, cfg_dir_with_metadata: Path
    ):
        """Entries from a secondary (global) dir must not appear in primary metadata."""
        upstream_cfg_dir = cfg_dir_with_metadata
        user_cfg_dir = tmp_path / "local_config"
        shutil.copytree(upstream_cfg_dir, user_cfg_dir)
        upstream_cfg_meta = upstream_cfg_dir / ".metadata"
        user_cfg_meta = user_cfg_dir / ".metadata"
        assert upstream_cfg_meta.is_file()
        assert user_cfg_meta.is_file()

        markerset = "bclasvirinae_odb10"
        hash = "51a9dce485c577e575d411d10c566820"
        mock_hash = mock_hashlib_gen(hash)
        monkeypatch.setattr(download_mod, "hashlib", mock_hash)

        with BuscoParser(user_cfg_dir, upstream_cfg_dir) as bp:
            bp.download(markerset)

        assert hash in user_cfg_meta.read_text(encoding="utf-8")
        assert hash not in upstream_cfg_meta.read_text(encoding="utf-8")
