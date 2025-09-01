from __future__ import annotations
import argparse
import importlib
import logging
from pathlib import Path

import pytest

import phyling
import phyling.__main__ as main


def test_cfg_dirs_with_env():
    assert phyling.CFG_DIRS[0].resolve().absolute() == Path("tests/database").resolve().absolute()


def test_cfg_dirs_wo_env(monkeypatch):
    monkeypatch.setenv("PHYLING_DB", "")
    importlib.reload(phyling)
    assert phyling.CFG_DIRS[0] == Path.home() / ".phyling"


def test_main_help():
    with pytest.raises(SystemExit) as excinfo:
        assert main.main(["--help"]) == main.main(["-h"])
        assert "usage: phyling" in str(excinfo.value)


def test_main_version(capsys: pytest.CaptureFixture):
    with pytest.raises(SystemExit) as excinfo:
        assert main.main(["--version"]) == main.main(["-V"])
        # captured: str = capsys.readouterr().out
        # captured = captured.split("\n")
        # assert captured[0] == captured[1] == phyling.VERSION
        assert phyling.VERSION in str(excinfo)


def test_main_verbose(caplog: pytest.LogCaptureFixture):
    with caplog.at_level(logging.DEBUG):
        assert main.main(["download", "list", "-v"]) == main.main(["download", "list", "--verbose"])

    assert "Debug mode enabled" in caplog.text


def test_main_FileHandler(monkeypatch, tmp_path: Path):
    def dummy_func(**kwargs):
        return 0

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verbose", action="store_true")
    parser.set_defaults(func=dummy_func)
    monkeypatch.setattr(main, "menu", lambda: parser)

    assert main.main(["--output", str(tmp_path)]) == 0

    assert (tmp_path / "log.txt").exists()

def test_main_KeyboardInterrupt(monkeypatch, caplog: pytest.LogCaptureFixture):
    def dummy_func(**kwargs):
        raise KeyboardInterrupt

    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true")
    parser.set_defaults(func=dummy_func)
    monkeypatch.setattr(main, "menu", lambda: parser)

    with caplog.at_level(logging.WARNING):
        assert main.main(["--verbose"]) == 1

    assert "Terminated by user." in caplog.text
