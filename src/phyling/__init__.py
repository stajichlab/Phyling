"""Phyling - Phylogenomic reconstruction from genomes."""

import logging
import os
from importlib.metadata import PackageNotFoundError, metadata
from pathlib import Path

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


try:
    _meta = metadata("phyling")
    __version__ = _meta["Version"]
    __author__ = _meta["Author-email"]
except PackageNotFoundError:  # running from an uninstalled source tree
    __version__ = "0+unknown"
    __author__ = ""

# Config folders from $PHYLING_DB; fall back to $HOME/.phyling when the first one is not writable.
# The folder is created lazily by the download pipeline.
CFG_DIRS: list[Path] = [Path(p) for p in os.getenv("PHYLING_DB", "").split(":") if p]
if not (CFG_DIRS and os.access(CFG_DIRS[0], os.W_OK)):
    CFG_DIRS.insert(0, Path.home() / ".phyling")

# Available CPUs
AVAIL_CPUS = int(os.environ.get("SLURM_CPUS_ON_NODE") or os.cpu_count() or 1)
