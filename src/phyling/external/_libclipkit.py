"""Clipkit utilities."""

from __future__ import annotations

import numpy as np
from Bio.Align import MultipleSeqAlignment
from clipkit.modes import TrimmingMode
from clipkit.msa import MSA

from phyling.lib import SeqTypes


def trim_gaps(msa: MultipleSeqAlignment, gaps: float = 0.9) -> MultipleSeqAlignment:
    """Trim columns from an MSA based on gappyness ratio.

    Args:
        msa (MultipleSeqAlignment): Multiple sequence alignment to trim.
        gaps (float, optional): Gappyness ratio for trimming. Defaults to 0.9.

    Returns:
        MultipleSeqAlignment: Trimmed alignment.

    Raises:
        ValueError: If gaps is not between 0 and 1.

    Example:
        >>> trimmed_msa = trim_gaps(msa, gaps=0.8)
    """
    if gaps > 1 or gaps < 0:
        raise ValueError('The argument "gaps" should be a float between 0 to 1.')
    infoList = [{"id": rec.id, "name": rec.name, "description": rec.description} for rec in msa]
    cds_msa = None
    if msa.annotations["seqtype"] == SeqTypes.DNA:
        pep_msa = MultipleSeqAlignment(
            [rec.translate(gap="-") for rec in msa],
            annotations={"seqtype": SeqTypes.PEP},
        )
        cds_msa = msa
    else:
        pep_msa = msa
    clipkit_pep_msa = MSA.from_bio_msa(pep_msa, gap_chars="-")
    clipkit_pep_msa.trim(mode=TrimmingMode.gappy, gap_threshold=gaps)
    pep_msa = clipkit_pep_msa.to_bio_msa()

    if cds_msa:
        if clipkit_pep_msa._site_positions_to_trim.size > 0:
            clipkit_cds_msa = MSA.from_bio_msa(cds_msa)
            pep_trimList_expanded = np.expand_dims(clipkit_pep_msa._site_positions_to_trim, axis=1)
            cds_site_positions_to_trim = (pep_trimList_expanded * np.array([3]) + np.array([0, 1, 2])).flatten()
            clipkit_cds_msa.trim(site_positions_to_trim=cds_site_positions_to_trim)
            cds_msa = clipkit_cds_msa.to_bio_msa()
        msa = cds_msa
    else:
        msa = pep_msa

    for new_rec, rec in zip(msa, infoList):
        new_rec.id, new_rec.name, new_rec.description = (
            rec["id"],
            rec["name"],
            rec["description"],
        )

    return msa
