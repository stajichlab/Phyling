"""Filter the multiple sequence alignment (MSA) results for tree module."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

from ..lib import FileExts, SeqTypes, TreeMethods
from ..lib._utils import Timer, check_threads
from ..lib.tree import MFA2TreeList, TreeOutputFiles
from ._outputprecheck import FilterPrecheck

logger = logging.getLogger(__name__)


@Timer.timer
@check_threads
def filter(
    inputs: str | Path | list[str | Path],
    output: str | Path,
    top_n_toverr: int,
    *,
    seqtype: Literal["dna", "pep", "AUTO"] = "AUTO",
    ml: bool = False,
    threads: int = 1,
    **kwargs,
) -> None:
    """A pipeline that filter the multiple sequence alignment results through their treeness/RCVs."""

    inputs_ = _input_check(inputs)
    output = Path(output)
    if not 1 < top_n_toverr < len(inputs_):
        if top_n_toverr == len(inputs_):
            raise SystemExit("Argument top_n_toverr is equal to the number of inputs. Do not need filtering.")
        elif len(inputs_) == 3:
            detail_msg = "can only be 2 since there are only 3 inputs"
        else:
            detail_msg = f"should between 2 to {len(inputs_) - 1}"
        raise ValueError(f"Argument top_n_toverr out of range. ({detail_msg})")

    logger.info("Found %s MSA fasta.", len(inputs_))

    mfa2treelist = MFA2TreeList(data=inputs_, seqtype=seqtype)

    # Params for precheck
    params = {"top_n_toverr": top_n_toverr}

    # Precheck and load checkpoint if it exist
    output_precheck = FilterPrecheck(output, mfa2treelist, **params)
    remained_mfa2treelist, completed_mfa2treelist = output_precheck.precheck()

    if remained_mfa2treelist:
        logger.info(
            "Use %s to generate trees and filter by the rank of their toverr.",
            TreeMethods.FT.method,
        )
        remained_mfa2treelist.build("ft", "LG" if mfa2treelist.seqtype == SeqTypes.PEP else "JC", noml=not (ml), threads=threads)
        remained_mfa2treelist.compute_toverr(threads=threads)
        completed_mfa2treelist.extend(remained_mfa2treelist)
    completed_mfa2treelist.sort()
    all_samples = set()
    for mfa2tree in completed_mfa2treelist:
        all_samples.update([tip.name for tip in mfa2tree.tree.get_terminals()])

    # Generate treeness tsv
    logger.info("Output selected fasta to folder %s...", output)
    treeness_file = output / TreeOutputFiles.TREENESS
    with open(treeness_file, "w") as f:
        f.write(f"# The MSA fasta which the toverr within top {top_n_toverr} are selected:\n")
        picked_samples = set()
        for mfa2tree in completed_mfa2treelist[:top_n_toverr]:
            f.write("\t".join([mfa2tree.name, str(mfa2tree.toverr)]) + "\n")
            picked_samples.update([tip.name for tip in mfa2tree.tree.get_terminals()])
        missing_samples = all_samples - picked_samples
        if missing_samples:
            logger.warning(
                "The following samples are totally missed in the selected markers: %s", ", ".join(sorted(missing_samples))
            )
            logger.warning(
                "They will not show up in the final tree. "
                + "If the presence of all samples is important please consider adjusting the argument --top_n_toverr/-n."
            )
        if completed_mfa2treelist[top_n_toverr:]:
            f.write("# The MSA fasta below are filtered out:\n")
            for mfa2tree in completed_mfa2treelist[top_n_toverr:]:
                f.write("\t".join([mfa2tree.name, str(mfa2tree.toverr)]) + "\n")

    # Symlink to seletced MSAs
    output = Path(output)
    files = [mfa2tree.file for mfa2tree in completed_mfa2treelist[:top_n_toverr]]
    for file in files:
        (output / file.name).symlink_to(file.absolute())
    logger.debug("Create symlinks done.")

    output_precheck.save_checkpoint(completed_mfa2treelist)
    logger.debug("Save checkpoint done.")

    logger.info(f"{__name__.split('.')[-1].capitalize()} module done.")


def _input_check(inputs: str | Path | list[str | Path]) -> tuple[Path, ...]:
    """Check and adjust the arguments passed in."""
    if isinstance(inputs, list):
        inputs_tuple = tuple(Path(file) for file in inputs)
        input_dir = {file.parent for file in inputs_tuple}
        if len(input_dir) > 1:
            raise RuntimeError("The inputs aren't in the same folder, which indicates it might come from different analysis.")
    else:
        inputs = Path(inputs)
        if inputs.is_file():
            inputs_tuple = (inputs,)
        else:
            inputs_tuple = tuple(file for file in inputs.glob(f"*.{FileExts.ALN}"))
            if not inputs_tuple:
                raise FileNotFoundError("Empty input directory.")

    if len(inputs_tuple) < 3:
        raise ValueError("Fewer than 3 inputs. Please directly build tree with your desired tree building software.")
    return inputs_tuple
