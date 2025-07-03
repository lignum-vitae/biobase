import pytest
from biobase.analysis import Dna
from biobase.constants.nucleic_acid import (
    MOLECULAR_WEIGHT,
    RNA_COMPLEMENTS,
    DNA_COMPLEMENTS,
    IUPAC_NUCLEOTIDES,
)


class TestMolecularWeights:
    """Tests for nucleotide molecular weights"""

    def test_molecular_weight_values(self):
        """Test that all standard nucleotides have correct molecular weights"""
        expected_weights = {
            "A": 135.13,
            "C": 111.10,
            "G": 151.13,
            "T": 126.11,
            "U": 112.09,
        }
        for nuc, weight in expected_weights.items():
            assert MOLECULAR_WEIGHT[nuc] == pytest.approx(weight, rel=1e-2)

    def test_molecular_weight_completeness(self):
        """Test that molecular weights exist for all standard nucleotides"""
        standard_nucs = set("ACGTU")
        assert all(nuc in MOLECULAR_WEIGHT for nuc in standard_nucs)

    def test_molecular_weight_type(self):
        """Test that all molecular weights are floats"""
        assert all(isinstance(weight, float) for weight in MOLECULAR_WEIGHT.values())

    def test_molecular_weight_positive(self):
        """Test that all molecular weights are positive"""
        assert all(weight > 0 for weight in MOLECULAR_WEIGHT.values())


class TestComplements:
    """Tests for nucleotide complement mappings"""

    def test_dna_complements(self):
        """Test DNA complement pairs"""
        expected_pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
        for nuc, complement in expected_pairs.items():
            assert DNA_COMPLEMENTS[nuc] == complement
            # Test reverse complement exists
            assert DNA_COMPLEMENTS[complement] == nuc

    def test_rna_complements(self):
        """Test RNA complement pairs"""
        expected_pairs = {"A": "U", "U": "A", "C": "G", "G": "C"}
        for nuc, complement in expected_pairs.items():
            assert RNA_COMPLEMENTS[nuc] == complement
            # Test reverse complement exists
            assert RNA_COMPLEMENTS[complement] == nuc


class TestIUPACNucleotides:
    """Tests for IUPAC nucleotide codes"""

    def test_iupac_standard_nucleotides(self):
        """Test that standard nucleotides are in IUPAC set"""
        standard_nucs = set("ACGTU")
        assert all(nuc in IUPAC_NUCLEOTIDES for nuc in standard_nucs)

    def test_iupac_degenerate_nucleotides(self):
        """Test that degenerate nucleotides are in IUPAC set"""
        degenerate_nucs = set("RYKMSWBDHVN")
        assert all(nuc in IUPAC_NUCLEOTIDES for nuc in degenerate_nucs)


class TestDataTypeConsistency:
    """Tests for data type consistency across constants"""

    def test_molecular_weight_type(self):
        """Test that MOLECULAR_WEIGHT is a dict with str keys and float values"""
        assert isinstance(MOLECULAR_WEIGHT, dict)
        assert all(isinstance(k, str) for k in MOLECULAR_WEIGHT.keys())
        assert all(isinstance(v, float) for v in MOLECULAR_WEIGHT.values())

    def test_complement_types(self):
        """Test that complement dictionaries have str keys and values"""
        for complement_dict in [DNA_COMPLEMENTS, RNA_COMPLEMENTS, IUPAC_NUCLEOTIDES]:
            assert isinstance(complement_dict, dict)
            assert all(isinstance(k, str) for k in complement_dict.keys())
            assert all(isinstance(v, str) for v in complement_dict.values())

    def test_iupac_nucleotides_type(self):
        """Test that IUPAC_NUCLEOTIDES is a dictionary of strings"""
        assert isinstance(IUPAC_NUCLEOTIDES, dict)
        assert all(isinstance(nuc, str) for nuc in IUPAC_NUCLEOTIDES)


class TestEdgeCases:
    """Tests for edge cases and potential issues"""

    def test_case_sensitivity(self):
        """Test that all constants use uppercase letters"""
        assert all(k.isupper() for k in MOLECULAR_WEIGHT.keys())
        assert all(k.isupper() for k in DNA_COMPLEMENTS.keys())
        assert all(k.isupper() for k in RNA_COMPLEMENTS.keys())
        assert all(nuc.isupper() for nuc in IUPAC_NUCLEOTIDES)

    def test_invalid_nucleotides(self):
        """Test that invalid nucleotides are not in any constants"""
        invalid_chars = set("EFILPQZ0123456789")
        assert not any(c in MOLECULAR_WEIGHT for c in invalid_chars)
        assert not any(c in DNA_COMPLEMENTS for c in invalid_chars)
        assert not any(c in RNA_COMPLEMENTS for c in invalid_chars)
        assert not any(c in IUPAC_NUCLEOTIDES for c in invalid_chars)


class TestBiologicalConsistency:
    """Tests for biological consistency of the data"""

    def test_dna_rna_consistency(self):
        """Test consistency between DNA and RNA data"""
        # T and U should have different weights
        assert MOLECULAR_WEIGHT["T"] != MOLECULAR_WEIGHT["U"]

        # T and U should complement to A
        assert DNA_COMPLEMENTS["T"] == "A"
        assert RNA_COMPLEMENTS["U"] == "A"

    def test_complement_symmetry(self):
        """Test that complement relationships are symmetric"""
        for complement_dict in [DNA_COMPLEMENTS, RNA_COMPLEMENTS]:
            for nuc, complement in complement_dict.items():
                if nuc != "N":  # N is self-complementary
                    assert complement_dict[complement] == nuc

    def test_molecular_weight_ratios(self):
        """Test that molecular weight ratios are biologically sensible"""
        # Purines (A, G) should be heavier than pyrimidines (C, T, U)
        assert MOLECULAR_WEIGHT["A"] > MOLECULAR_WEIGHT["C"]
        assert MOLECULAR_WEIGHT["G"] > MOLECULAR_WEIGHT["C"]
        assert MOLECULAR_WEIGHT["A"] > MOLECULAR_WEIGHT["T"]
        assert MOLECULAR_WEIGHT["G"] > MOLECULAR_WEIGHT["T"]
        assert MOLECULAR_WEIGHT["A"] > MOLECULAR_WEIGHT["U"]
        assert MOLECULAR_WEIGHT["G"] > MOLECULAR_WEIGHT["U"]


@pytest.mark.parametrize("seq, expected", [("aTcG", "ATCG")])
def test_validate_dna_sequence_valid(seq, expected):
    assert Dna._validate_dna_sequence(seq) == expected


@pytest.mark.parametrize("seq", ["ATCU", "XYZ", "123", "ACGTU"])
def test_validate_dna_sequence_invalid_chars(seq):
    with pytest.raises(ValueError, match="Invalid DNA nucleotides"):
        Dna._validate_dna_sequence(seq)


@pytest.mark.parametrize(
    "seq", ["", None]
)  # None should raise TypeError inside isinstance check
def test_validate_dna_sequence_empty_or_none(seq):
    with pytest.raises(ValueError):
        Dna._validate_dna_sequence(seq)


def test_validate_dna_sequence_non_str():
    with pytest.raises(ValueError, match="Expected string input"):
        Dna._validate_dna_sequence(123)  # type: ignore[arg-type]


def test_transcribe_basic():
    assert Dna.transcribe("aTcG") == "AUCG"


@pytest.mark.parametrize(
    "seq, expected_rev, expected_no_rev",
    [
        ("ATCG", "CGAT", "TAGC"),
        ("GGCCAA", "TTGGCC", "CCGGTT"),
        ("atgc", "GCAT", "TACG"),
    ],
)
def test_complement_dna(seq, expected_rev, expected_no_rev):
    assert Dna.complement_dna(seq) == expected_rev
    assert Dna.complement_dna(seq, reverse=False) == expected_no_rev


def test_transcribe_invalid_input_raises():
    with pytest.raises(ValueError):
        Dna.transcribe("ATBX")


def test_complement_dna_empty_raises():
    with pytest.raises(ValueError):
        Dna.complement_dna("")


def test_gc_content_invalid_raises():
    with pytest.raises(ValueError):
        Dna.calculate_gc_content("NNNN")
