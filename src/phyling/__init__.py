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

# Config folders from $PHYLING_DB; fall back to $XDG_CACHE_HOME/phyling when the first one is not writable.
# The folder is created lazily by the download pipeline.
# Per the XDG spec, an unset, empty or relative $XDG_CACHE_HOME means $HOME/.cache
_xdg_cache_home = Path(os.getenv("XDG_CACHE_HOME") or Path.home() / ".cache")
DEFAULT_CFG_DIR = (_xdg_cache_home if _xdg_cache_home.is_absolute() else Path.home() / ".cache") / "phyling"
# Deprecated location, still searched for existing markersets until support is removed
LEGACY_CFG_DIR = Path.home() / ".phyling"

CFG_DIRS: list[Path] = [Path(p) for p in os.getenv("PHYLING_DB", "").split(":") if p]
if not (CFG_DIRS and os.access(CFG_DIRS[0], os.W_OK)):
    CFG_DIRS.insert(0, DEFAULT_CFG_DIR)
    if LEGACY_CFG_DIR.is_dir() and LEGACY_CFG_DIR not in CFG_DIRS:
        CFG_DIRS.append(LEGACY_CFG_DIR)


def warn_legacy_cfg_dir(*markersets: str) -> None:
    """Warn that the given markersets are still read from the deprecated $HOME/.phyling folder."""
    logger.warning(
        "Markerset(s) %s in %s are deprecated and support will be removed in a future release. Please move them to %s.",
        ", ".join(markersets),
        LEGACY_CFG_DIR,
        DEFAULT_CFG_DIR,
    )


# Available CPUs
AVAIL_CPUS = int(os.environ.get("SLURM_CPUS_ON_NODE") or os.cpu_count() or 1)
