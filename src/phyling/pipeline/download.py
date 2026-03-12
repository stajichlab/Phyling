"""Help to download/update BUSCO v5 markerset to a local folder.

First it checks whether the metadata file is exist under the config folder ~/.phyling. A missing or outdated file will trigger the
module to download/update the metadata.

Passing "list" to markerset argument will list all the available/already downloaded markersets. Passing a valid name to the
markerset argument will download the markerset to the config folder ~/.phyling/HMM.
"""

from __future__ import annotations

import argparse
import shutil
from urllib.error import URLError

from .. import CFG_DIRS
from ..lib._utils import Timer
from ..lib.download import BuscoParser

__all__ = ["download"]


def menu(parser: argparse.ArgumentParser) -> None:
    """Menu for download module."""
    parser.add_argument("markerset", metavar='HMM markerset or "list"', help="Name of the HMM markerset")
    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    parser.set_defaults(func=download)


@Timer.timer
def download(markerset: str, **kwargs) -> None:
    """A pipeline that list the available BUSCO HMM markerset and download it when specifying."""
    try:
        with BuscoParser(*CFG_DIRS) as metadata:
            markerset_list = metadata.online + metadata.local
            if markerset == "list":
                width, _ = shutil.get_terminal_size((80, 24))
                col = width // 40
                markerset_list = [markerset_list[x : x + col] for x in range(0, len(markerset_list), col)]
                col_width = max(len(word) for row in markerset_list for word in row) + 3  # padding

                if metadata.online:
                    msg = "Datasets available online:"
                    _wrapper(metadata.online, col=col, col_width=col_width, msg=msg)

                if metadata.local:
                    msg = "Datasets available on local:"
                    _wrapper(metadata.local, col=col, col_width=col_width, msg=msg)
            elif markerset in metadata.online:
                metadata.download(markerset)
            else:
                raise RuntimeError(
                    'Markerset not available: %s Please check it again with "list" option.',
                    markerset,
                )
    except URLError as e:
        raise URLError("Connection lost or URL currently not available: %s", e)


def _wrapper(item_list: list[str], col: int, col_width: int, msg: str | None = None) -> None:
    """Adjust databases display according to the terminal size."""
    item_list = [item_list[x : x + col] for x in range(0, len(item_list), col)]
    if msg:
        print(msg)
        print()
    for row in item_list:
        # Print the database list
        print(" ".join(word.ljust(col_width) for word in row))
    print()
