"""Help to download/update BUSCO v5 markerset to a local folder."""

from __future__ import annotations

import shutil
from urllib.error import URLError

from .. import CFG_DIRS
from ..lib._utils import Timer
from ..lib.download import BuscoParser

__all__ = ["download"]


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
        raise URLError(e.reason)
        raise URLError("Connection lost or URL currently not available") from e


def _wrapper(item_list: list[str], col: int, col_width: int, msg: str | None = None) -> None:
    """Adjust databases display according to the terminal size."""
    items = [item_list[x : x + col] for x in range(0, len(item_list), col)]
    if msg:
        print(msg)
        print()
    for row in items:
        # Print the database list
        print(" ".join(word.ljust(col_width) for word in row))
    print()
