from __future__ import annotations

import importlib
import logging
import sys
from argparse import Namespace
from pathlib import Path

import pytest

import phyling
from phyling.main import main, run_cli


class TestModuleInit:
    def test_cfg_dirs_with_env(self):
        assert phyling.CFG_DIRS[0].resolve().absolute() == Path("tests/database").resolve().absolute()

    def test_cfg_dirs_wo_env(self, monkeypatch):
        monkeypatch.setenv("PHYLING_DB", "")
        importlib.reload(phyling)
        assert phyling.CFG_DIRS[0] == Path.home() / ".phyling"


class TestRunCli:
    @pytest.mark.parametrize("exit_code", [0, 1])
    def test_exit_code(self, monkeypatch: pytest.MonkeyPatch, exit_code: int):
        monkeypatch.setattr(sys, "argv", ["phyling", "download", "list"])
        monkeypatch.setattr("phyling.main.main", lambda args: exit_code)

        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        assert excinfo.value.code == exit_code

    def test_version(self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
        # Force main to return exit code 0
        monkeypatch.setattr("phyling.main.main", lambda args: 0)

        monkeypatch.setattr(sys, "argv", ["phyling", "-V"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out1 = capsys.readouterr().out
        assert excinfo.value.code == 0

        monkeypatch.setattr(sys, "argv", ["phyling", "--version"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out2 = capsys.readouterr().out

        assert excinfo.value.code == 0
        assert out1.strip() == out2.strip() == phyling.VERSION

    def test_help(self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
        # Force main to return exit code 0
        monkeypatch.setattr("phyling.main.main", lambda args: 0)

        monkeypatch.setattr(sys, "argv", ["phyling", "-h"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out1 = capsys.readouterr().out
        assert excinfo.value.code == 0

        monkeypatch.setattr(sys, "argv", ["phyling", "--help"])
        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        out2 = capsys.readouterr().out

        assert excinfo.value.code == 0
        assert out1 == out2
        assert "usage: phyling" in out1

    def test_missing_args(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["phyling"])
        # Force main to return exit code 0
        monkeypatch.setattr("phyling.main.main", lambda args: 0)

        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        assert excinfo.value.code == 2

    def test_invalid_args(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["phyling", "invalid"])
        # Force main to return exit code 0
        monkeypatch.setattr("phyling.main.main", lambda args: 0)

        with pytest.raises(SystemExit) as excinfo:
            run_cli()

        assert excinfo.value.code == 2


@pytest.fixture
def args_factory() -> Namespace:
    """Returns a function that creates a mock args object with defaults."""

    def create_args(**kwargs):
        # Define your standard defaults here
        defaults = {
            "module": "fake module",
            "output": "test",
            "verbose": True,
            "func": lambda **kwargs: 0,
        }
        defaults.update(kwargs)
        return Namespace(**defaults)

    return create_args


@pytest.fixture
def mock_main_dependencies(monkeypatch: pytest.MonkeyPatch):
    """Fixture to mock all heavy dependencies in the refine module."""
    # Mock Path.mkdir and Path.unlink to prevent actual disk changes
    monkeypatch.setattr(Path, "mkdir", lambda *args, **kwargs: None)
    monkeypatch.setattr(Path, "rename", lambda self, target: None)


class TestMain:
    def test_main_success(self, args_factory, mock_main_dependencies):
        """Test that main returns 0 on a successful run."""
        args = args_factory()
        result = main(args)
        assert result == 0

    def test_main_FileHandler(self, args_factory, tmp_path: Path, mock_main_dependencies):

        args = args_factory()
        args.output = tmp_path

        main(args)

        assert (tmp_path / "log.txt").exists()

    def test_main_keyboard_interrupt(self, args_factory, mock_main_dependencies):
        """Test that main catches Ctrl+C and returns 130."""

        def mock_interrupt(**kwargs):
            raise KeyboardInterrupt

        args = args_factory()
        args.func = mock_interrupt

        result = main(args)
        assert result == 130

    def test_main_filenotfound(self, args_factory, mock_main_dependencies):
        """Test that main catches FileNotFound and returns 2."""

        def mock_interrupt(**kwargs):
            raise FileNotFoundError

        args = args_factory()
        args.func = mock_interrupt

        result = main(args)
        assert result == 2

    def test_main_exception(self, args_factory, mock_main_dependencies):
        """Test that main catches other errors and returns 1."""

        def mock_interrupt(**kwargs):
            raise RuntimeError

        args = args_factory()
        args.func = mock_interrupt

        result = main(args)
        assert result == 1
