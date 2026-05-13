"""Perform multiple sequence alignment (MSA) on orthologous sequences that match the hmm markers across samples."""

from __future__ import annotations

import logging
import re
from multiprocessing import Manager, Pool
from pathlib import Path
from typing import Literal

from Bio import SeqIO

from .. import CFG_DIRS
from ..exception import EmptyWarning
from ..external._libclipkit import trim_gaps
from ..lib import ALIGN_METHODS, FileExts, SeqTypes
from ..lib._utils import Timer, check_threads
from ..lib.align import HMMMarkerSet, OrthologList, SampleList
from ._outputprecheck import AlignPrecheck

logger = logging.getLogger(__name__)


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

    inputs_, markerset_, evalue_, method_ = _args_check(inputs, markerset, evalue, method)
    output = Path(output)

    logger.info("Found %s samples.", len(inputs_))
    names = [
        re.sub(
            r"(\.(aa|pep|cds|fna|faa))?\.(fasta|fas|faa|fna|seq|fa)(\.gz)?",
            "",
            sample.name,
        )
        for sample in inputs_
    ]
    samples = SampleList(inputs_, names, seqtype=seqtype)

    logger.info("Loading markerset from %s...", markerset_)
    hmmmarkerset = HMMMarkerSet(markerset_, markerset_.parent / "scores_cutoff")
    hmmmarkerset.sort(key=lambda x: x.name)
    logger.debug("Load markerset done.")

    # Params for precheck
    params = {
        "markerset": tuple(hmmmarkerset.checksums.keys()),
        "markerset_cutoff": "markerset cutoff" if hmmmarkerset.have_cutoffs() else evalue_,
        "method": method_,
    }
    # Precheck and load checkpoint if it exist
    output_precheck = AlignPrecheck(output, samples, **params)
    remaining_samples, searchhits = output_precheck.precheck()
    if remaining_samples:
        logger.info("Search start...")
        jobs, threads_per_job = _search_threads_check(len(remaining_samples), threads)
        hits = remaining_samples.search(hmmmarkerset, evalue=evalue_, jobs=jobs, threads=threads_per_job)
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
    orthologs = OrthologList(list(searchhits.orthologs.values()), list(searchhits.orthologs.keys()), seqtype=seqtype)

    logger.info("Alignment start...")
    align_kwargs = {}
    align_kwargs["jobs"] = threads
    msa_list = orthologs.align(method=method_, hmms=hmmmarkerset if method_ == "hmmalign" else None, jobs=threads)  # type: ignore
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
) -> tuple[tuple[Path, ...], Path, float, str]:
    """Check and adjust the arguments passed in."""
    if isinstance(inputs, list):
        inputs_tuple = tuple(Path(sample) for sample in inputs)
    else:
        inputs = Path(inputs)
        if inputs.is_file():
            inputs_tuple = (inputs,)
        else:
            inputs_tuple = tuple(Path(inputs).iterdir())
            if not inputs_tuple:
                raise FileNotFoundError("Empty input directory.")
    for sample in inputs_tuple:
        if sample.is_dir():
            raise IsADirectoryError("Input files contain directory.")
    if len(inputs_tuple) < 4:
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

    return inputs_tuple, hmmmarkerset, evalue, method


def _search_threads_check(n_samples: int, threads: int):
    threads_per_job = threads if (threads < 8 or n_samples == 1) else 4
    jobs = threads // threads_per_job
    if jobs > n_samples:
        jobs = n_samples
    return jobs, threads_per_job
