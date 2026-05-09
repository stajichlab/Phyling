"""Argparse menu for Phyling."""

from __future__ import annotations

import argparse
import re
import textwrap
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Callable, Protocol, cast

from . import AUTHOR, AVAIL_CPUS, VERSION
from . import __name__ as _module_name
from .lib import ALIGN_METHODS, TreeMethods
from .pipeline.align import align
from .pipeline.download import download
from .pipeline.filter import filter
from .pipeline.tree import tree


class Args(Protocol):
    """Arguments shared by ALL entry points.

    Attributes:
        verbose (bool): Flag to enable verbose logging.
    """

    module: str
    output: Path
    verbose: bool
    func: Callable


class BaseModule(ABC):
    """The blueprint for every PHYLing module."""

    help: str
    description: str
    func: Callable

    @classmethod
    @abstractmethod
    def add_args(cls, parser: argparse.ArgumentParser) -> None:
        """Add module-specific arguments to the parser."""
        parser.set_defaults(func=cls.func)


class CustomHelpFormatter(argparse.HelpFormatter):
    """Custom help formatter for argparse to enhance text wrapping and default value display.

    This formatter adjusts the text wrapping for better readability and ensures that the default value of an argument is displayed
    in the help message if not already present.
    """

    def _fill_text(self, text, width, indent):
        """Format the class/function docstring with appropriate wrapping."""
        text = [self._whitespace_matcher.sub(" ", paragraph.strip()) for paragraph in text.split("\n\n") if paragraph.strip()]
        return "\n\n".join([textwrap.fill(line, width) for line in text])

    def _split_lines(self, text, width):
        """Enable multi-line display in argument help message."""
        text = [self._whitespace_matcher.sub(" ", line.strip()) for line in text.split("\n") if line.strip()]
        return [wrapped_line for line in text for wrapped_line in textwrap.wrap(line, width)]

    def _get_help_string(self, action):
        """Allow additional message after default parameter displayed."""
        help = action.help
        pattern = r"\(default: .+\)"
        if re.search(pattern, action.help) is None:
            if action.default not in [argparse.SUPPRESS, None, False]:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += " (default: %(default)s)"
        return help


def main_menu(argv=None) -> Args:
    parser = argparse.ArgumentParser(
        prog=_module_name,
        formatter_class=CustomHelpFormatter,
        description="""Phyling - Phylogenomic reconstruction from genomes.

        Phyling comprises 4 modules - download, align, filter and tree. The download module can be used to download HMM markerset
        from BUSCO. The align module is the core element of this package which generate multiple sequence alignment among the
        orthologs found across samples. The filter module calculates treeness/RCV scores to filter out the uninformative alignment
        results. The tree module help to build a phylogenetic tree by different algorithms.
        """,
        epilog=f"Written by {AUTHOR}",
        add_help=False,
    )

    # Configuration for each module
    modules = {"download": Download, "align": Align, "filter": Filter, "tree": Tree}

    subparsers = parser.add_subparsers(title="Modules", dest="module")

    for name, module in modules.items():
        subparser = subparsers.add_parser(
            name,
            formatter_class=CustomHelpFormatter,
            help=module.help,
            description=module.description,
            add_help=False,
        )
        # Add module-specific arguments
        module.add_args(subparser)

        # Add shared help/options to each subparser if needed
        # (Alternatively, keep them only on the parent)

    # Global options
    opt_args = parser.add_argument_group("Options")
    opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
    opt_args.add_argument("-V", "--version", action="version", version=VERSION)

    parsed_args = parser.parse_args(argv)
    if not hasattr(parsed_args, "func"):
        parser.error("No module specified. Please choose a module from the list above.")

    return cast(Args, parsed_args)


class Download(BaseModule):
    help = "Run multiple sequence alignments"
    description = """Help to download/update BUSCO v5 markerset to a local folder.

    First it checks whether the metadata file is exist under the config folder ~/.phyling. A missing or outdated file will trigger
    the module to download/update the metadata.

    Passing "list" to markerset argument will list all the available/already downloaded markersets. Passing a valid name to the
    markerset argument will download the markerset to the config folder ~/.phyling/HMM.
    """
    func = download

    @classmethod
    def add_args(cls, parser: argparse.ArgumentParser) -> None:
        """Menu for download module."""
        parser.add_argument("markerset", metavar='HMM markerset or "list"', help="Name of the HMM markerset")
        opt_args = parser.add_argument_group("Options")
        opt_args.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")
        opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
        super().add_args(parser)


class Align(BaseModule):
    help = "Download HMM markers"
    description = """Perform multiple sequence alignment (MSA) on orthologous sequences that match the hmm markers across samples.

    Initially, hmmsearch is used to match the samples against a given markerset and report the top hit of each sample for each hmm
    marker, representing "orthologs" across all samples. In order to build a tree, minimum of 4 samples should be used. If the
    bitscore cutoff file is present in the hmms folder, it will be used as the cutoff. Otherwise, an evalue of 1e-10 will be used
    as the default cutoff.

    Sequences corresponding to orthologs found in more than 4 samples are extracted from each input. These sequences then undergo
    MSA with hmmalign or muscle. The resulting alignments are further trimmed using clipkit by default. You can use the --non_trim
    option to skip the trimming step. Finally, The alignment results are output separately for each hmm marker.
    """
    func = align

    @classmethod
    def add_args(cls, parser: argparse.ArgumentParser) -> None:
        """Menu for align module."""
        req_args = parser.add_argument_group("Required arguments")
        input_type = req_args.add_mutually_exclusive_group(required=True)
        input_type.add_argument(
            "-i",
            "--inputs",
            dest="inputs",
            metavar=("file", "files"),
            nargs="+",
            type=Path,
            help="Query pepetide/cds fasta or gzipped fasta",
        )
        input_type.add_argument(
            "-I",
            "--input_dir",
            dest="inputs",
            metavar="directory",
            type=Path,
            help="Directory containing query pepetide/cds fasta or gzipped fasta",
        )
        req_args.add_argument(
            "-m",
            "--markerset",
            metavar="directory",
            type=Path,
            required=True,
            help="Directory of the HMM markerset",
        )
        opt_args = parser.add_argument_group("Options")
        opt_args.add_argument(
            "-o",
            "--output",
            metavar="directory",
            type=Path,
            default="phyling-align-%s" % time.strftime("%Y%m%d-%H%M%S", time.gmtime()),
            help="Output directory of the alignment results (default: phyling-align-[YYYYMMDD-HHMMSS] (UTC timestamp))",
        )
        opt_args.add_argument(
            "--seqtype",
            choices=["dna", "pep", "AUTO"],
            default="AUTO",
            help="Input data sequence type",
        )
        opt_args.add_argument(
            "-E",
            "--evalue",
            metavar="float",
            type=float,
            default=1e-10,
            help="Hmmsearch reporting threshold (default: %(default)s, "
            + "only being used when bitscore cutoff file is not available)",
        )
        opt_args.add_argument(
            "-M",
            "--method",
            choices=ALIGN_METHODS,
            default="hmmalign",
            help="Program used for multiple sequence alignment",
        )
        opt_args.add_argument(
            "--non_trim",
            action="store_true",
            help="Report non-trimmed alignment results",
        )
        opt_args.add_argument(
            "-t",
            "--threads",
            type=int,
            default=AVAIL_CPUS // 4 * 4,
            help="Threads for hmmsearch and the number of parallelized jobs in MSA step. "
            + "Better be multiple of 4 if using more than 8 threads",
        )
        opt_args.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")
        opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")
        super().add_args(parser)


class Filter(BaseModule):
    help = "Filter MSA results"
    description = """Filter the multiple sequence alignment (MSA) results for tree module.

    The align step usually reports a lot of markers but many of them are uninformative or susceptible to composition bias. The
    Treeness/RCV value computed by PhyKIT is used to estimate how informative the markers are. By default the -n/--top_n_toverr is
    set to 50 to select only the top 50 markers.
    """
    func = filter

    @classmethod
    def add_args(cls, parser: argparse.ArgumentParser) -> None:
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
        req_args.add_argument(
            "-n",
            "--top_n_toverr",
            type=int,
            required=True,
            help="Select the top n markers based on their treeness/RCV for final tree building",
        )
        opt_args = parser.add_argument_group("Options")
        opt_args.add_argument(
            "-o",
            "--output",
            metavar="directory",
            type=Path,
            default="phyling-filter-%s" % time.strftime("%Y%m%d-%H%M%S", time.gmtime()),
            help="Output directory of the treeness.tsv and selected MSAs "
            + "(default: phyling-tree-[YYYYMMDD-HHMMSS] (UTC timestamp))",
        )
        opt_args.add_argument(
            "--seqtype",
            choices=["pep", "dna", "AUTO"],
            default="AUTO",
            help="Input data sequence type",
        )
        opt_args.add_argument(
            "--ml",
            action="store_true",
            help="Use maximum-likelihood estimation during tree building",
        )
        opt_args.add_argument(
            "-t",
            "--threads",
            type=int,
            default=AVAIL_CPUS,
            help="Threads for filtering",
        )
        opt_args.add_argument("-v", "--verbose", action="store_true", help="Verbose mode for debug")
        opt_args.add_argument("-h", "--help", action="help", help="show this help message and exit")


class Tree(BaseModule):
    help = "Build a phylogenetic tree"
    description = """Construct a phylogenetic tree by the selected multiple sequence alignment (MSA) results.

    By default the consensus tree method will be employed which use a 50% cutoff to represent the majority of all the trees. You
    can use the -c/--concat option to concatenate the MSA and build a single tree instead. Note that enable the -p/--partition
    option will also output a partition file that compatible to RAxML-NG and IQ-TREE.

    For the tree building step, the FastTree will be used as default algorithm. Users can switch to the RAxML-NG or IQ-TREE by
    specifying the -m/--method raxml/iqtree.

    Once the tree is built, an ASCII figure representing the tree will be displayed, and a treefile in Newick format will be
    generated as output. Additionally, users can choose to obtain a matplotlib-style figure using the -f/--figure option.
    """
    func = tree

    @classmethod
    def add_args(cls, parser: argparse.ArgumentParser) -> None:
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
