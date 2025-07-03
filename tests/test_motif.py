import pytest
from biobase.analysis import find_motifs


# Test data fixtures
@pytest.fixture
def single_sequence():
    return "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture
def fasta_dict():
    return {
        ">SP001": "ACDEFCDEFCDEFGHIKLMN",  # has matches for "CDE" at positions 2, 6, 10
        ">SP002": "MNPQRSTVWYACDEFGHIKL",  # has match for "CDE" at position 12
        ">SP003": "AAAAAAAAAAAAAAAAAA12",  # invalid: contains "1", "2"
        ">SP004": "GGGGGGGGGGGGGGGGGGGG",  # no match
        ">SP005": "HHHHHHHHHHHHHHHHH@#$",  # invalid: contains "@", "#", "$"
        ">SP006": "DDDDDDDDDDDDDDDDDDDD",  # no match
        ">SP007": "CDEFGHCDEFKLCDEFPQRS",  # has matches for "CDE" at positions 1, 7, 13
        ">SP008": "LLLLLLLLLLLLLLLLLLLL",  # no match
        ">SP009": "KKKKKKKKKKKK123KKKKK",  # invalid: contains "1", "2", "3"
        ">SP010": "CDEACDEBCDEFAAAAAAAA",  # has matches for "CDE" at positions 1, 5, 9
    }


class TestSingleSequence:
    """Tests for single sequence input"""

    def test_valid_sequence_with_match(self, single_sequence):
        result = find_motifs(single_sequence, "DEF")
        assert result == [3]

    def test_valid_sequence_no_match(self):
        result = find_motifs("GGGGGGGGGGGGGGGGGGGG", "CDE")
        assert result == []

    def test_empty_sequence(self):
        with pytest.raises(ValueError, match="empty"):
            find_motifs("", "CDE")

    def test_invalid_sequence_chars(self):
        with pytest.raises(ValueError, match="Invalid"):
            find_motifs("ACDEF123GHIKLMNPQRSTVWY", "CDE")

    def test_overlapping_matches(self):
        result = find_motifs("CDEFDEFGHI", "DEF")
        assert result == [2, 5]

    def test_empty_pattern(self):
        with pytest.raises(ValueError, match="empty"):
            find_motifs("ACDEFGHIKLMNPQRSTVWY", "")

    def test_pattern_longer_than_sequence(self):
        result = find_motifs("CDE", "CDEFG")
        assert result == []

    @pytest.mark.parametrize(
        "sequence,pattern,expected",
        [
            ("ACDEFGHIKLMNPQRSTVWY", "CDE", [2]),
            ("CDEFGHIKLMNPQRSTVWY", "CDE", [1]),
            ("ACDEFCDEFCDEF", "CDE", [2, 6, 10]),
            (
                "AAAAAAAAAAAAAAAAAAAA",
                "AAA",
                [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
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
        assert result_dict[">SP001"] == [2, 6, 10]
        assert result_dict[">SP002"] == [12]
        assert result_dict[">SP007"] == [1, 7, 13]
        assert result_dict[">SP010"] == [1, 5, 9]

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
        assert result_dict[">SP001"] == [2]
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
        assert result == [3]

    def test_extended_invalid_sequence(self):
        with pytest.raises(ValueError):
            find_motifs("ACDEFGHIKLMNPQRSTVWYBJZX123", "CDE", ext=True)

    def test_extended_fasta(self):
        fasta = {
            ">SP001": "ACDEFGHIKLMNPQRSTVWYBJZX",
            ">SP002": "BJZXACDEFGHIKLMNPQRSTVWY",
        }
        result_dict, invalid_dict, no_matches = find_motifs(fasta, "CDE", ext=True)
        assert result_dict[">SP001"] == [2]
        assert result_dict[">SP002"] == [6]


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


class TestPerformance:
    """Basic performance tests"""

    def test_long_sequence(self):
        # Test with a 10,000 character sequence
        sequence = "ACDEFGHIKL" * 1000
        result = find_motifs(sequence, "CDE")
        assert len(result) == 1000
        assert result[0] == 2

    def test_long_fasta(self):
        # Test with 1000 sequences
        fasta = {f">SP{i:03d}": "ACDEFGHIKL" * 10 for i in range(1000)}
        result_dict, invalid_dict, no_matches = find_motifs(fasta, "CDE")
        assert len(result_dict) == 1000
        assert not invalid_dict
        assert not no_matches
