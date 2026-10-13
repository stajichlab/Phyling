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

    def test_cfg_dirs_wo_env(self, monkeypatch: pytest.MonkeyPatch):
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
        assert out1.strip() == out2.strip() == phyling.__version__

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

    return create_args  # type: ignore


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

    def test_log_captures_pipeline_records(self, args_factory, tmp_path: Path):
        """Records from phyling.pipeline.* and phyling.lib.* are written to the log file."""
        output = tmp_path / "out"

        def func(**kwargs):
            output.mkdir()
            logging.getLogger("phyling.pipeline.fake").warning("from pipeline")

        main(args_factory(output=output, func=func))

        assert "from pipeline" in (output / "log.txt").read_text()

    def test_log_written_outside_existing_output_during_run(self, args_factory, tmp_path: Path):
        """The log file does not appear in the output folder while the module runs."""
        output = tmp_path / "out"
        output.mkdir()
        seen = {}

        def func(**kwargs):
            seen["files"] = list(output.iterdir())

        main(args_factory(output=output, func=func))

        assert seen["files"] == []
        assert (output / "log.txt").exists()

    def test_failed_run_keeps_log_outside_output(self, args_factory, tmp_path: Path):
        output = tmp_path / "out"

        def func(**kwargs):
            output.mkdir()
            raise RuntimeError("boom")

        assert main(args_factory(output=output, func=func)) == 1
        assert list(output.iterdir()) == []
        assert "boom" in (tmp_path / "out.log").read_text()


class TestMainMenu:
    def test_align_threads_default_is_positive(self, monkeypatch: pytest.MonkeyPatch):
        from phyling import _menu

        monkeypatch.setattr(_menu, "AVAIL_CPUS", 2)

        args = _menu.main_menu(["align", "-I", "x", "-m", "y"])

        assert args.threads == 2

    def test_module_help_text(self, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
        from phyling import _menu

        with pytest.raises(SystemExit):
            _menu.main_menu(["-h"])
        out = " ".join(capsys.readouterr().out.split())

        assert "download Download HMM markers" in out
        assert "align Run multiple sequence alignments" in out
