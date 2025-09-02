import pytest
from biobase.analysis import find_motifs


# Test data fixtures
@pytest.fixture
def single_sequence():
    return "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture
def fasta_dict():
    return {
        ">SP001": "ACDEFCDEFCDEFGHIKLMN",  # has matches for "CDE" that span indexes [(1, 4), (5, 8), (9, 12)]
        ">SP002": "MNPQRSTVWYACDEFGHIKL",  # has match for "CDE" that span indexes [(11, 14)]
        ">SP003": "AAAAAAAAAAAAAAAAAA12",  # invalid: contains "1", "2"
        ">SP004": "GGGGGGGGGGGGGGGGGGGG",  # no match
        ">SP005": "HHHHHHHHHHHHHHHHH@#$",  # invalid: contains "@", "#", "$"
        ">SP006": "DDDDDDDDDDDDDDDDDDDD",  # no match
        ">SP007": "CDEFGHCDEFKLCDEFPQRS",  # has matches for "CDE" that span indexes [(0, 3), (6, 9), (12, 15)]
        ">SP008": "LLLLLLLLLLLLLLLLLLLL",  # no match
        ">SP009": "KKKKKKKKKKKK123KKKKK",  # invalid: contains "1", "2", "3"
        ">SP010": "CDEACDEDCDEFAAAAAAAA",  # has matches for "CDE" that span indexes [(0, 3), (4, 7), (8, 11)]
    }


class TestSingleSequence:
    """Tests for single sequence input"""

    def test_valid_sequence_with_match(self, single_sequence):
        result = find_motifs(single_sequence, "DEF")
        assert result == [(2, 5)]

    def test_valid_sequence_no_match(self):
        result = find_motifs("GGGGGGGGGGGGGGGGGGGG", "CDE")
        assert result == []

    def test_empty_sequence(self):
        with pytest.raises(ValueError, match="empty"):
            find_motifs("", "CDE")

    def test_invalid_sequence_chars(self):
        with pytest.raises(ValueError, match="Invalid"):
            find_motifs("ACDEF123GHIKLMNPQRSTVWY", "CDE")

    def test_adjacent_matches(self):
        result = find_motifs("CDEFDEFGHI", "DEF")
        assert result == [(1, 4), (4, 7)]

    def test_overlapping_matches(self):
        result = find_motifs("CEDEDEFGHI", "EDE")
        assert result == [(1, 4), (3, 6)]

    def test_empty_pattern(self):
        with pytest.raises(ValueError, match="empty"):
            find_motifs("ACDEFGHIKLMNPQRSTVWY", "")

    def test_pattern_longer_than_sequence(self):
        result = find_motifs("CDE", "CDEFG")
        assert result == []

    @pytest.mark.parametrize(
        "sequence,pattern,expected",
        [
            ("ACDEFGHIKLMNPQRSTVWY", "CDE", [(1, 4)]),
            ("CDEFGHIKLMNPQRSTVWY", "CDE", [(0, 3)]),
            ("ACDEFCDEFCDEF", "CDE", [(1, 4), (5, 8), (9, 12)]),
            (
                "AAAAAAAAAA",
                "AAA",
                [(0, 3), (1, 4), (2, 5), (3, 6), (4, 7), (5, 8), (6, 9), (7, 10)],
            ),
        ],
    )
    def test_various_patterns(self, sequence, pattern, expected):
        result = find_motifs(sequence, pattern)
        assert result == expected


class TestFASTADictionary:
    """Tests for FASTA dictionary input"""

    def test_valid_fasta_dict(self, fasta_dict):
        result_dict, invalid_dict, no_matches = find_motifs(fasta_dict, "CDE")

        # Check sequences with matches
        assert result_dict[">SP001"] == [(1, 4), (5, 8), (9, 12)]
        assert result_dict[">SP002"] == [(11, 14)]
        assert result_dict[">SP007"] == [(0, 3), (6, 9), (12, 15)]
        assert result_dict[">SP010"] == [(0, 3), (4, 7), (8, 11)]

        # Check invalid sequences
        assert ">SP003" in invalid_dict
        assert ">SP005" in invalid_dict
        assert ">SP009" in invalid_dict

        # Check sequences with no matches
        assert sorted(no_matches) == sorted([">SP004", ">SP006", ">SP008"])

    def test_empty_fasta_dict(self):
        with pytest.raises(ValueError, match="empty"):
            find_motifs({}, "CDE")

    def test_single_entry_fasta(self):
        fasta = {">SP001": "ACDEFGHIKLMNPQRSTVWY"}
        result_dict, invalid_dict, no_matches = find_motifs(fasta, "CDE")
        assert result_dict[">SP001"] == [(1, 4)]
        assert not invalid_dict
        assert not no_matches

    def test_all_invalid_sequences(self):
        fasta = {
            ">SP001": "123456789",
            ">SP002": "@#$%^&*",
        }
        result_dict, invalid_dict, no_matches = find_motifs(fasta, "CDE")
        assert not result_dict
        assert len(invalid_dict) == 2
        assert not no_matches

    def test_all_no_matches(self):
        fasta = {
            ">SP001": "GGGGGGGGGG",
            ">SP002": "AAAAAAAAAA",
        }
        result_dict, invalid_dict, no_matches = find_motifs(fasta, "CDE")
        assert not result_dict
        assert not invalid_dict
        assert len(no_matches) == 2

    def test_empty_sequence_in_fasta(self):
        fasta = {">SP001": ""}
        with pytest.raises(ValueError):
            find_motifs(fasta, "CDE")


class TestExtendedAminoAcids:
    """Tests for extended amino acid codes"""

    def test_extended_valid_sequence(self):
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        result = find_motifs(sequence, "DEF", ext=True)
        assert result == [(2, 5)]

    def test_extended_invalid_sequence(self):
        with pytest.raises(ValueError):
            find_motifs("ACDEFGHIKLMNPQRSTVWYBJZX123", "CDE", ext=True)

    def test_extended_fasta(self):
        fasta = {
            ">SP001": "ACDEFGHIKLMNPQRSTVWYUUOU",
            ">SP002": "DUOUACDEFGHIKLMNPQRSTVWY",
        }
        result_dict, invalid_dict, no_matches = find_motifs(fasta, "CDE", ext=True)
        assert result_dict[">SP001"] == [(1, 4)]
        assert result_dict[">SP002"] == [(5, 8)]
        assert not invalid_dict
        assert not no_matches


class TestEdgeCases:
    """Tests for edge cases and error handling"""

    def test_invalid_input_type(self):
        with pytest.raises(Exception):  # Could be TypeError or ValueError
            find_motifs(123, "CDE")

    def test_invalid_pattern_type(self):
        with pytest.raises(ValueError):
            find_motifs("ACDEF", 123)

    def test_none_input(self):
        with pytest.raises(ValueError):
            find_motifs("", "CDE")

    def test_none_pattern(self):
        with pytest.raises(ValueError):
            find_motifs("ACDEF", None)

    @pytest.mark.parametrize("invalid_char", ["1", "2", "@", "#", "$", " ", "\n", "\t"])
    def test_specific_invalid_chars(self, invalid_char):
        with pytest.raises(ValueError, match="Invalid"):
            find_motifs(f"ACDEF{invalid_char}GHIKL", "CDE")
