"""Help to download/update BUSCO v5 markerset to a local folder."""

from __future__ import annotations

import logging
import shutil
from urllib.error import URLError

from .. import CFG_DIRS
from ..lib._utils import Timer
from ..lib.download import BuscoParser

__all__ = ["download"]

logger = logging.getLogger(__name__)


@Timer.timer
def download(markerset: str, **kwargs) -> None:
    """A pipeline that list the available BUSCO HMM markerset and download it when specifying."""
    CFG_DIRS[0].mkdir(parents=True, exist_ok=True)
    try:
        with BuscoParser(*CFG_DIRS) as metadata:
            markerset_list = metadata.online + metadata.local
            if markerset == "list":
                width, _ = shutil.get_terminal_size((80, 24))
                if not markerset_list:
                    logger.warning("No markersets available online or on local.")
                    return
                col = max(width // 40, 1)
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
            elif metadata.online_error:
                raise RuntimeError(f"Could not retrieve online metadata: {metadata.online_error}") from metadata.online_error
            else:
                raise RuntimeError(f'Markerset not available: {markerset}. Please check it again with "list" option.')
    except URLError as e:
        raise URLError(f"Connection lost or URL currently not available: {e.reason}") from e


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
