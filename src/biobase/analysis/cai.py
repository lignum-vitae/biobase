"""
Codon Adaptation Index (CAI)

This module provides a function to compute the Codon Adaptation Index (CAI),
a measure of how well a coding sequence’s codon usage matches a reference set
of highly expressed genes in a given organism.

Public module variables:
(None)

Key components:
1. CAI computation based on per–amino-acid synonymous codon weights
2. Robust handling of DNA/RNA input (DNA 'T' converted to RNA 'U')
3. Exclusion of stop codons from CAI
4. Reference normalization (RNA uppercase) and per-family maxima

Notes:
- CAI is the geometric mean of per-codon weights, where each weight is the
  codon’s reference frequency divided by the maximum reference frequency among
  synonymous codons for the same amino acid.
- Stop codons are excluded from both the product and the length L.
- Reference counts are organism-specific and must be supplied by the caller.
"""

from math import log, exp
from typing import Mapping
from biobase.constants.amino_acid import CODON_TABLE


def _build_family_max(
    ref_counts: Mapping[str, int | float]
) -> tuple[dict[str, float], dict[str, float]]:
    """
    Normalize reference counts to RNA codons and compute per–amino-acid maxima.

    Returns
    -------
    ref_rna : dict[str, float]
        Reference codon counts with RNA keys (T->U, uppercase).
    family_max : dict[str, float]
        Amino acid -> maximum reference count among its synonymous codons.
        STOP codons are ignored.
    """
    # Normalize reference to RNA uppercase
    ref_rna: dict[str, float] = {
        codon.upper().replace("T", "U"): float(v) for codon, v in ref_counts.items()
    }

    family_max: dict[str, float] = {}
    for codon, aa in CODON_TABLE.items():
        if aa == "STOP":
            continue
        f = ref_rna.get(codon, 0.0)
        if f > family_max.get(aa, 0.0):
            family_max[aa] = f

    return ref_rna, family_max


def cai(
    seq: str,
    ref_counts: Mapping[str, int | float],
) -> float:
    """
    Calculate the Codon Adaptation Index (CAI) for a coding sequence.

    This function computes the geometric mean of per-codon weights, where each
    weight is the codon’s reference frequency divided by the maximum reference
    frequency among synonymous codons for the same amino acid. DNA input is
    converted to RNA by replacing 'T' with 'U'. Stop codons are ignored.

    Parameters
    ----------
    seq : str
        DNA or RNA coding sequence. Whitespace/newlines are ignored.
        DNA bases 'T' are converted to RNA 'U' for codon lookups.
    ref_counts : Mapping[str, int | float]
        Reference codon usage counts or frequencies (organism-specific).
        Keys may be DNA or RNA codons; they are normalized internally to
        RNA uppercase (T->U).

    Returns
    -------
    float
        CAI value in [0.0, 1.0]. Returns 0.0 if no valid codons contribute.

    Raises
    ------
    ValueError
        If the sequence contains invalid codons (after normalization).

    Examples
    --------
    >>> ref = {"AAA": 80, "AAG": 20, "GCC": 50, "GCU": 10, "GCA": 20, "GCG": 20}
    >>> seq = "AAA AAG GCC GCU"
    >>> round(cai(seq, ref), 4)
    0.4729
    """
    # 1) Normalize sequence: uppercase, strip spaces/newlines, DNA->RNA
    s = seq.upper().replace(" ", "").replace("\n", "")
    s = s.replace("T", "U")

    # 2) Split into full codons only
    codons = [s[i:i+3] for i in range(0, len(s) - (len(s) % 3), 3)]

    # 3) Validate codons against the genetic code
    invalid = [c for c in codons if c not in CODON_TABLE]
    if invalid:
        raise ValueError(f"Invalid codons found in sequence: {invalid}")

    # 4) Build reference lookups
    ref_rna, family_max = _build_family_max(ref_counts)

    # 5) Geometric mean of per-codon weights (w = f_i / max_f_in_family)
    total_log = 0.0
    L = 0

    for codon in codons:
        aa = CODON_TABLE[codon]
        if aa == "STOP":
            continue  # exclude stops from CAI

        f_i = ref_rna.get(codon)
        max_f = family_max.get(aa)

        # Skip if this codon/family has no coverage in the reference
        if not f_i or not max_f:
            continue

        w_i = f_i / max_f  # in (0, 1]
        if w_i <= 0.0:
            continue

        total_log += log(w_i)
        L += 1

    return exp(total_log / L) if L > 0 else 0.0
