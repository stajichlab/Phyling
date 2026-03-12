"""Construct a phylogenetic tree by the selected multiple sequence alignment (MSA) results.

By default the consensus tree method will be employed which use a 50% cutoff to represent the majority of all the trees. You can
use the -c/--concat option to concatenate the MSA and build a single tree instead. Note that enable the -p/--partition option will
also output a partition file that compatible to RAxML-NG and IQ-TREE.

For the tree building step, the FastTree will be used as default algorithm. Users can switch to the RAxML-NG or IQ-TREE by
specifying the -m/--method raxml/iqtree.

Once the tree is built, an ASCII figure representing the tree will be displayed, and a treefile in Newick format will be generated
as output. Additionally, users can choose to obtain a matplotlib-style figure using the -f/--figure option.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import time
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
from Bio import Phylo

from .. import AVAIL_CPUS
from ..external import ModelFinder
from ..lib import FileExts, SeqTypes, TreeMethods, TreeOutputFiles
from ..lib._utils import Timer, check_threads
from ..lib.tree import MFA2Tree, MFA2TreeList

logger = logging.getLogger(__name__)


def menu(parser: argparse.ArgumentParser) -> None:
    """Menu for tree module."""
    req_args = parser.add_argument_group("Required arguments")
    input_type = req_args.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        dest="inputs",
        metavar=("file", "files"),
        nargs="+",
        type=Path,
        help="Multiple sequence alignment fasta of the markers",
    )
    input_type.add_argument(
        "-I",
        "--input_dir",
        dest="inputs",
        metavar="directory",
        type=Path,
        help="Directory containing multiple sequence alignment fasta of the markers",
    )
    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument(
        "-o",
        "--output",
        metavar="directory",
        type=Path,
        default="phyling-tree-%s" % time.strftime("%Y%m%d-%H%M%S", time.gmtime()),
        help="Output directory of the newick treefile (default: phyling-tree-[YYYYMMDD-HHMMSS] (UTC timestamp))",
    )
    opt_args.add_argument(
        "--seqtype",
        choices=["pep", "dna", "AUTO"],
        default="AUTO",
        help="Input data sequence type",
    )
    opt_args.add_argument(
        "-M",
        "--method",
        choices=[m.name.lower() for m in TreeMethods],
        default=TreeMethods.FT.name.lower(),
        help="Algorithm used for tree building. (default: %(default)s)\n"
        + "Available options:\n"
        + "\n".join(f"{m.name.lower()}: {m.method}" for m in TreeMethods),
    )
    opt_args.add_argument(
        "-c",
        "--concat",
        action="store_true",
        help="Concatenated alignment results",
    )
    opt_args.add_argument(
        "-p",
        "--partition",
        action="store_true",
        help="Partitioned analysis by sequence. Only works when --concat enabled.",
    )
    opt_args.add_argument("-f", "--figure", action="store_true", help="Generate a matplotlib tree figure")
    opt_args.add_argument(
        "--seed",
        type=int,
        default=-1,
        help="Seed number for stochastic elements during inferences. (default: %(default)s to generate randomly)",
    )
    opt_args.add_argument(
        "-t",
        "--threads",
        type=int,
        default=AVAIL_CPUS,
        help="Threads for tree construction",
    )
    opt_args.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    parser.set_defaults(func=tree)


@Timer.timer
@check_threads
def tree(
    inputs: str | Path | list[str | Path],
    output: str | Path,
    *,
    method: Literal["ft", "raxml", "iqtree"] = "ft",
    seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    bs: int = 1000,
    scfl: int = 100,
    concat: bool = False,
    partition: bool = False,
    figure: bool = False,
    seed: int = -1,
    threads: int = 1,
    **kwargs,
) -> None:
    """A pipeline that build phylogenetic tree through either FastTree, RAxML-NG or IQ-TREE."""

    inputs = _input_check(inputs)
    logger.info("Found %s MSA fasta.", len(inputs))

    mfa2treelist = MFA2TreeList(inputs, seqtype=seqtype)
    partition = _validate_partition(partition, method, concat)

    output = Path(output)
    output.mkdir(exist_ok=True)

    if concat:
        logger.info("Inference with concatenation mode...")
        concat_file, partition_file = mfa2treelist.concat(output=output, threads=threads)
        concat_tree = MFA2Tree(concat_file, seqtype=seqtype)
        concat_tree.load()
        threads_max = threads
        if len(concat_tree._data) > 100 and len(concat_tree) > 10000:
            threads = min(min(threads, len(concat_tree._data)), 32)  # Use at most 32 cpus
        else:
            threads = -1  # Auto detect optimal threads
        concat_tree.unload()

        if partition:
            modelfinder_runner = ModelFinder(
                concat_tree.file,
                output / method / "modelfinder",
                partition_file,
                seqtype=concat_tree.seqtype,
                method=method,
                seed=seed,
                threads=threads,
                threads_max=threads_max,
            )
            modelfinder_runner.run()
            new_partition_file = modelfinder_runner.result
            partition_file.unlink()
            partition_file.symlink_to(new_partition_file.parent.relative_to(partition_file.parent) / new_partition_file.name)
            tree = concat_tree.build(
                method, output / method, partition_file, bs=bs, scfl=scfl, seed=seed, threads=threads, threads_max=threads_max
            )
        else:
            tree = concat_tree.build(
                method, output / method, "AUTO", bs=bs, scfl=scfl, seed=seed, threads=threads, threads_max=threads_max
            )

    else:
        logger.info("Inference with consensus mode...")
        if bs > 0 or scfl > 0:
            logger.warning("Bootstrap and concordance factor calculation are disabled with consensus mode.")
        mfa2treelist.build(method, "LG" if mfa2treelist.seqtype == SeqTypes.PEP else "GTR", bs=0, scfl=0, jobs=threads, threads=1)
        logger.info("Tree building done.")
        tree = mfa2treelist.get_consensus_tree(output / "astral", seed=seed)

    """Output the tree in newick format and figure."""
    output_tree = output / TreeOutputFiles.TREE_NW
    logger.info("Output tree to %s", output_tree)
    shutil.copy(tree, output_tree)

    if figure:
        fig, ax = plt.subplots(figsize=(20, 12))
        output_fig = output / TreeOutputFiles.TREE_IMG
        logger.info("Output figure to %s", output_fig)
        Phylo.draw(Phylo.read(tree, "newick"), do_show=False, axes=ax)
        fig.savefig(output_fig)

    logger.info(f"{__name__.split('.')[-1].capitalize()} module done.")


def _input_check(inputs: str | Path | list) -> tuple[Path]:
    """Check and adjust the arguments passed in."""
    if isinstance(inputs, list):
        inputs = tuple(Path(file) for file in inputs)
        input_dir = {file.parent for file in inputs}
        if len(input_dir) > 1:
            raise RuntimeError("The inputs aren't in the same folder, which indicates it might come from different analysis.")
    else:
        inputs = Path(inputs)
        if inputs.is_file():
            inputs = (inputs,)
        else:
            inputs = tuple(file for file in inputs.glob(f"*.{FileExts.ALN}"))
            if not inputs:
                raise FileNotFoundError("Empty input directory.")

    if len(inputs) == 1:
        raise ValueError("Found only 1 MSA fasta. Please directly build tree with your desired tree building software.")
    return inputs


def _validate_partition(partition: bool, method: str, concat: bool) -> bool:
    """Validate and preprocess partition settings."""
    method = method.upper()
    if not isinstance(partition, bool):
        raise ValueError("Argument partition only accepts boolean value.")
    if partition:
        if method not in (TreeMethods.RAXML.name, TreeMethods.IQTREE.name):
            raise ValueError(f"Partition is not allowed with {TreeMethods[method].method}.")
        if not concat:
            raise ValueError("Partition is not allowed in consensus mode.")

    return partition
