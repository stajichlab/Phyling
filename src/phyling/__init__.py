"""Phyling - Phylogenomic reconstruction from genomes."""

import logging
import os
from importlib.metadata import metadata
from pathlib import Path

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


VERSION = metadata("phyling")["Version"]
AUTHOR = metadata("phyling")["Author-email"]

# Create config folder in $HOME/.phyling
CFG_DIRS: list[Path] = []
phyling_db_loc = os.getenv("PHYLING_DB")
if phyling_db_loc:
    CFG_DIRS.extend([Path(path) for path in phyling_db_loc.split(":")])
if CFG_DIRS and os.access(CFG_DIRS[0], os.W_OK):
    pass
else:
    cfg_dir = Path.home() / ".phyling"
    cfg_dir.mkdir(exist_ok=True)
    CFG_DIRS.insert(0, cfg_dir)

# Available CPUs
AVAIL_CPUS = int(os.environ.get("SLURM_CPUS_ON_NODE") or os.cpu_count() or 1)
