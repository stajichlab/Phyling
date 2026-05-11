#!/usr/bin/env python3
"""The Phyling CLI menu and sub-process execution."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from ._menu import Args, main_menu

LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


def run_cli() -> None:
    """Runs the orthoSynAssign CLI entry point."""
    parsed_args = main_menu(sys.argv[1:])
    print(parsed_args)
    sys.exit(main(parsed_args))


def main(args: Args) -> int:
    """Main function to parse command-line arguments and execute the corresponding functionality."""

    _setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    logger.debug("Debug mode enabled.")
    logger.debug(vars(args))
    file_handler = None

    try:
        output_path = getattr(args, "output", None)
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            log_file = output_path.with_suffix(".log") if not output_path.is_dir() else output_path / "log.txt"

            file_handler = logging.FileHandler(log_file, mode="w", delay=True)
            file_handler.setFormatter(logging.Formatter(LOG_FORMAT))
            logger.addHandler(file_handler)
            logger.debug("Log file initialized.")

        try:
            args.func(**vars(args))

        except KeyboardInterrupt:
            logger.warning("Terminated by user.")
            return 130

        except FileNotFoundError as e:
            logger.error("An error occurred: %s", e)
            logger.debug("Traceback details:", exc_info=True)
            return 2

        except Exception as e:
            logger.error("An error occurred: %s", e)
            logger.debug("Traceback details:", exc_info=True)
            return 1

    finally:
        if file_handler:
            logger.removeHandler(file_handler)
            file_handler.close()
            if Path(args.output).is_dir():
                Path(file_handler.baseFilename).rename(args.output / "log.txt")

    return 0


def _setup_logging(verbose: bool = False) -> None:
    """Configure logging for the application.

    Args:
        verbose (bool): Flag to enable verbose logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format=LOG_FORMAT,
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


if __name__ == "__main__":
    run_cli()
