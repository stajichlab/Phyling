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
    sys.exit(main(parsed_args))


def main(args: Args) -> int:
    """Main function to parse command-line arguments and execute the corresponding functionality."""

    _setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    logger.debug("Debug mode enabled.")
    logger.debug(vars(args))
    file_handler = None
    succeeded = False

    try:
        output_path = getattr(args, "output", None)
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            # Write the log next to the output folder during the run so it neither makes the folder look non-empty to the
            # output precheck nor gets removed by the rerun cleanup. It is moved into the folder once the run succeeds.
            log_file = output_path.parent / f"{output_path.name}.log"

            file_handler = logging.FileHandler(log_file, mode="w", delay=True)
            file_handler.setFormatter(logging.Formatter(LOG_FORMAT))
            # Attach to the package logger so records from phyling.pipeline.* and phyling.lib.* are captured too
            logging.getLogger("phyling").addHandler(file_handler)
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

        succeeded = True

    finally:
        if file_handler:
            logging.getLogger("phyling").removeHandler(file_handler)
            file_handler.close()
            source_log = Path(file_handler.baseFilename)
            # Keep a failed run's log outside the output folder so it does not block the next attempt
            if succeeded and source_log.exists() and Path(args.output).is_dir():
                source_log.replace(Path(args.output) / "log.txt")

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
