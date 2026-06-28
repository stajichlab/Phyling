"""Tests for the download pipeline module."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

import phyling.pipeline.download as download_mod
from phyling.pipeline.download import _wrapper, download


@pytest.fixture
def mock_busco_parser(monkeypatch: pytest.MonkeyPatch):
    mock_metadata = MagicMock()
    mock_metadata.online = ["fungi_odb10", "bacteria_odb10", "poxviridae_odb10"]
    mock_metadata.local = []
    mock_metadata.__enter__ = MagicMock(return_value=mock_metadata)
    mock_metadata.__exit__ = MagicMock(return_value=False)

    monkeypatch.setattr(download_mod, "BuscoParser", MagicMock(return_value=mock_metadata))
    return mock_metadata


# ---------------------------------------------------------------------------
# download
# ---------------------------------------------------------------------------


class TestDownload:
    def test_list_option_empty_local(self, mock_busco_parser, monkeypatch: pytest.MonkeyPatch):
        mock_wrap = MagicMock()
        monkeypatch.setattr(download_mod, "_wrapper", mock_wrap)

        download("list")

        assert mock_wrap.call_count == 1
        mock_busco_parser.download.assert_not_called()

    def test_list_option_has_local(self, mock_busco_parser, monkeypatch: pytest.MonkeyPatch):
        mock_busco_parser.local = ["fungi_odb10"]

        mock_wrap = MagicMock()
        monkeypatch.setattr(download_mod, "_wrapper", mock_wrap)

        download("list")

        assert mock_wrap.call_count == 2
        mock_busco_parser.download.assert_not_called()

    def test_download_markerset_not_in_local(self, mock_busco_parser):
        download("fungi_odb10")
        mock_busco_parser.download.assert_called_once_with("fungi_odb10")

    def test_download_markerset_already_in_local(self, mock_busco_parser):
        mock_busco_parser.local = ["fungi_odb10"]
        download("fungi_odb10")
        mock_busco_parser.download.assert_called_once_with("fungi_odb10")

    def test_download_invalid_markerset_raises(self, mock_busco_parser):
        with pytest.raises(RuntimeError):
            download("nonexistent_odb10")

    def test_url_error_raises(self, mock_busco_parser):
        from urllib.error import URLError

        mock_busco_parser.__enter__ = MagicMock(side_effect=URLError("connection refused"))
        mock_busco_parser.__exit__ = MagicMock(return_value=False)

        with pytest.raises(URLError):
            download("list")


# ---------------------------------------------------------------------------
# _wrapper
# ---------------------------------------------------------------------------


class TestWrapper:
    def test_prints_message(self, capsys: pytest.CaptureFixture):
        _wrapper(["item1", "item2", "item3"], col=2, col_width=15, msg="Test header:")
        captured = capsys.readouterr()
        assert "Test header:" in captured.out

    def test_no_message(self, capsys: pytest.CaptureFixture):
        _wrapper(["item1", "item2"], col=2, col_width=15, msg=None)
        captured = capsys.readouterr()
        assert "item1" in captured.out

    def test_single_item(self, capsys: pytest.CaptureFixture):
        _wrapper(["only_item"], col=3, col_width=20, msg="Single:")
        captured = capsys.readouterr()
        assert "only_item" in captured.out

    def test_empty_list(self, capsys: pytest.CaptureFixture):
        _wrapper([], col=3, col_width=20, msg="Empty:")
        captured = capsys.readouterr()
        assert "Empty:" in captured.out

    def test_items_are_printed(self, capsys: pytest.CaptureFixture):
        items = ["fungi_odb10", "bacteria_odb10", "archaea_odb10"]
        _wrapper(items, col=2, col_width=20)
        captured = capsys.readouterr()
        for item in items:
            assert item in captured.out
