import pytest
from biobase.analysis.cai import cai


def test_basic_cai_example():
    """Verify CAI on a small, hand-computable example."""
    # Reference counts (synthetic)
    ref_counts = {
        "AAA": 80,  # Lysine (K)
        "AAG": 20,
        "GCC": 50,  # Alanine (A)
        "GCU": 10,
        "GCA": 20,
        "GCG": 20,
    }
    # Sequence: AAA AAG GCC GCU (K, K, A, A)
    seq = "AAA AAG GCC GCU"
    result = cai(seq, ref_counts)
    # Expected CAI = (1.0 * 0.25 * 1.0 * 0.2) ** (1/4) ≈ 0.4729
    expected = 0.4729
    assert abs(result - expected) < 1e-4, f"Expected ~{expected}, got {result}"


def test_handles_dna_input_t_to_u_conversion():
    """DNA input should be converted to RNA (T->U) and processed normally."""
    ref_counts = {"AAA": 10, "AAG": 5, "UUU": 8}  # include UUU for middle TTT
    seq = "AAATTTAAA"  # becomes "AAA UUU AAA"
    result = cai(seq, ref_counts)
    assert isinstance(result, float)
    assert 0.0 <= result <= 1.0


def test_empty_and_stop_only_sequences_return_zero():
    """Empty sequence or sequences that only contribute STOP codons should yield 0.0."""
    ref_counts = {"AAA": 10}
    assert cai("", ref_counts) == 0.0
    # "UAA" is STOP; ignored → no contributing codons → 0.0
    assert cai("UAA", ref_counts) == 0.0
    # STOP mixed with partial trailing codon that doesn't form a full triplet
    assert cai("UAAU", ref_counts) == 0.0


def test_trailing_partial_codon_is_ignored():
    """A trailing incomplete triplet must be dropped."""
    ref_counts = {"AAA": 10}
    # "AAAAA" -> "AAA" + "AA" (drop "AA")
    assert cai("AAAAA", ref_counts) > 0.0  # contributed by the first "AAA"


def test_invalid_codon_raises_value_error():
    """Invalid codons after normalization must raise ValueError."""
    ref_counts = {"AAA": 10}
    with pytest.raises(ValueError):
        cai("XYZ", ref_counts)  # not in codon table
    # Also ensure invalid length triplets aren't the issue here: use a full invalid triplet
    with pytest.raises(ValueError):
        cai("ABU", ref_counts)


def test_mixed_case_and_whitespace_are_handled():
    """Upper/lower case and whitespace/newlines should be normalized away."""
    ref_counts = {"AAA": 10, "AAG": 5}
    seq1 = "aAa \n aAg"
    seq2 = "AAA AAG"
    assert cai(seq1, ref_counts) == cai(seq2, ref_counts)


def test_unknown_codon_in_reference_is_skipped_not_error():
    """
    If a codon appears in the sequence but is missing in ref_counts,
    it is skipped (not treated as zero) for CAI computation.
    """
    # Only AAA has reference coverage; AAG missing on purpose.
    ref_counts = {"AAA": 100}
    seq = "AAA AAG AAA"  # AAG contributes no weight; both AAAs do
    v = cai(seq, ref_counts)
    assert 0.0 < v <= 1.0


def test_stop_codons_are_ignored_when_mixed():
    """STOP codons mixed into a valid sequence should be ignored."""
    ref_counts = {"AAA": 50, "AAG": 50}
    # UAA (STOP) between two valid codons
    seq = "AAA UAA AAG"
    v = cai(seq, ref_counts)
    assert 0.0 < v <= 1.0
