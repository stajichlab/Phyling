"""Perform multiple sequence alignment (MSA) on orthologous sequences that match the hmm markers across samples.

Initially, hmmsearch is used to match the samples against a given markerset and report the top hit of each sample for each hmm
marker, representing "orthologs" across all samples. In order to build a tree, minimum of 4 samples should be used. If the
bitscore cutoff file is present in the hmms folder, it will be used as the cutoff. Otherwise, an evalue of 1e-10 will be used as
the default cutoff.

Sequences corresponding to orthologs found in more than 4 samples are extracted from each input. These sequences then undergo MSA
with hmmalign or muscle. The resulting alignments are further trimmed using clipkit by default. You can use the --non_trim option
to skip the trimming step. Finally, The alignment results are output separately for each hmm marker.
"""

from __future__ import annotations

import argparse
import re
import time
from multiprocessing import Manager, Pool
from pathlib import Path
from typing import Literal

from Bio import SeqIO

from .. import AVAIL_CPUS, CFG_DIRS
from ..exception import EmptyWarning
from ..external._libclipkit import trim_gaps
from ..libphyling import ALIGN_METHODS, FileExts, SeqTypes
from ..libphyling._utils import Timer, check_threads
from ..libphyling.align import HMMMarkerSet, OrthologList, SampleList
from . import logger
from ._outputprecheck import AlignPrecheck


def menu(parser: argparse.ArgumentParser) -> None:
    """Menu for align module."""
    req_args = parser.add_argument_group("Required arguments")
    input_type = req_args.add_mutually_exclusive_group(required=True)
    input_type.add_argument(
        "-i",
        "--inputs",
        dest="inputs",
        metavar=("file", "files"),
        nargs="+",
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
        help="Hmmsearch reporting threshold (default: %(default)s, only being used when bitscore cutoff file is not available)",
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
    parser.set_defaults(func=align)


@Timer.timer
@check_threads
def align(
    inputs: str | Path | list[str | Path],
    output: str | Path,
    *,
    markerset: str | Path,
    seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    evalue: float = 1e-10,
    method: Literal["hmmalign", "muscle"] = "hmmalign",
    non_trim: bool = False,
    threads: int = 1,
    **kwargs,
) -> None:
    """A pipeline that do hmmsearch to identify orthologs and align them through hmmalign or MUSCLE."""

    inputs, markerset, evalue, method = _args_check(inputs, markerset, evalue, method)

    logger.info("Found %s samples.", len(inputs))
    names = [
        re.sub(
            r"(\.(aa|pep|cds|fna|faa))?\.(fasta|fas|faa|fna|seq|fa)(\.gz)?",
            "",
            sample.name,
        )
        for sample in inputs
    ]
    samples = SampleList(inputs, names, seqtype=seqtype)

    logger.info("Loading markerset from %s...", markerset)
    hmmmarkerset = HMMMarkerSet(markerset, markerset.parent / "scores_cutoff")
    hmmmarkerset.sort(key=lambda x: x.name)
    logger.debug("Load markerset done.")

    # Params for precheck
    params = {
        "markerset": tuple(hmmmarkerset.checksums.keys()),
        "markerset_cutoff": "markerset cutoff" if hmmmarkerset.have_cutoffs() else evalue,
        "method": method,
    }
    # Precheck and load checkpoint if it exist
    output_precheck = AlignPrecheck(output, samples, **params)
    remaining_samples, searchhits = output_precheck.precheck()
    if remaining_samples:
        logger.info("Search start...")
        jobs, threads_per_job = _search_threads_check(len(remaining_samples), threads)
        hits = remaining_samples.search(hmmmarkerset, evalue=evalue, jobs=jobs, threads=threads_per_job)
        searchhits.update(hits)
        logger.info("Search done.")

    output_precheck.save_checkpoint(searchhits)
    logger.debug("Save checkpoint done.")

    try:
        searchhits = searchhits.filter(min_taxa=4)
    except EmptyWarning:
        raise RuntimeError(
            "All orthologs were gone after filtering. Please confirm whether the inputs have sufficient "
            "number of sample or if the correct HMM markers were being used."
        )
    logger.debug("Filter hits done.")

    searchhits.load()
    orthologs = OrthologList(searchhits.orthologs.values(), searchhits.orthologs.keys(), seqtype=samples.seqtype)

    logger.info("Alignment start...")
    msa_list = orthologs.align(method=method, hmms=hmmmarkerset if method == "hmmalign" else None, jobs=threads)
    logger.info("Alignment done.")

    if not non_trim:
        logger.info("Trimming start...")
        if threads == 1:
            msa_list = [trim_gaps(msa) for msa in msa_list]
        else:
            manager = Manager()
            with Pool(threads) as pool:
                msa_list = pool.map(trim_gaps, manager.list(msa_list))
        logger.info("Trimming done.")

    logger.info("Output individual fasta to folder %s...", output)
    ext = FileExts.CDS_ALN if samples.seqtype == SeqTypes.DNA else FileExts.PEP_ALN

    for msa, hmm in zip(msa_list, orthologs.names):
        msa.sort()
        with open(output_precheck.output / f"{hmm}.{ext}", "w") as f:
            SeqIO.write(msa, f, format="fasta")

    logger.info(f"{__name__.split('.')[-1].capitalize()} module done.")


def _args_check(
    inputs: str | Path | list[str | Path],
    hmmmarkerset: str | Path,
    evalue: float,
    method: str,
) -> tuple[tuple[Path], Path, float, str]:
    """Check and adjust the arguments passed in."""
    if isinstance(inputs, list):
        inputs = tuple(Path(sample) for sample in inputs)
    else:
        inputs = Path(inputs)
        if inputs.is_file():
            inputs = (inputs,)
        else:
            inputs = tuple(Path(inputs).iterdir())
            if not inputs:
                raise FileNotFoundError("Empty input directory.")
    for sample in inputs:
        if sample.is_dir():
            raise IsADirectoryError("Input files contain directory.")
    if len(inputs) < 4:
        raise ValueError("Requires at least 4 samples.")

    hmmmarkerset = Path(hmmmarkerset)
    if not hmmmarkerset.exists():
        for cfg_dir in CFG_DIRS:
            if Path(cfg_dir, hmmmarkerset).exists():
                hmmmarkerset = Path(cfg_dir, hmmmarkerset)
                if hmmmarkerset.name != "hmms" and [f for f in hmmmarkerset.glob("hmms") if f.is_dir()]:
                    hmmmarkerset = hmmmarkerset / "hmms"
    if not hmmmarkerset.exists():
        raise FileNotFoundError(f"Markerset folder does not exist: {hmmmarkerset} - did you download BUSCO?")

    if evalue >= 1:
        raise ValueError(f"Invalid evalue: {evalue}")
    if method not in ALIGN_METHODS:
        raise ValueError(f"Invalid method: {method}")

    return inputs, hmmmarkerset, evalue, method


def _search_threads_check(n_samples: int, threads: int):
    threads_per_job = threads if (threads < 8 or n_samples == 1) else 4
    jobs = threads // threads_per_job
    if jobs > n_samples:
        jobs = n_samples
    return jobs, threads_per_job
