import pytest
from pathlib import Path
from biobase.matrix import _Matrix, Blosum, Pam, Identity, Match


# Fixtures for common test setups
@pytest.fixture
def matrix_dir():
    project_root = Path(__file__).parent.parent.resolve()
    return (project_root / "src" / "biobase" / "matrix" / "matrices").resolve()


@pytest.fixture
def blosum45():
    return Blosum(45)


@pytest.fixture
def blosum62():
    return Blosum(62)


@pytest.fixture
def pam250():
    return Pam(250)


@pytest.fixture
def identity():
    return Identity(0)


@pytest.fixture
def match():
    return Match()


# Base Matrix Class Tests
class TestMatrixBase:
    def test_init_default_folder(self):
        matrix = _Matrix()
        assert matrix.folder.name == "matrices"
        assert matrix.folder.parent.name == "matrix"
        assert matrix.folder.parent.parent.name == "biobase"
        assert matrix.folder.parent.parent.parent.name == "src"
        assert matrix.matrix_data is None
        assert matrix.matrix_name is None
        assert matrix.version is None

    def test_init_custom_folder(self, tmp_path):
        matrix = _Matrix(matrix_folder=tmp_path)
        assert matrix.folder == tmp_path

    def test_available_matrices(self):
        matrices = _Matrix.available_matrices()
        assert "BLOSUM62" in matrices
        assert "PAM250" in matrices
        assert "IDENTITY0" in matrices
        assert "MATCH" in matrices

    def test_select_matrix_valid(self):
        matrix = _Matrix()
        matrix.select_matrix("Blosum", 62)
        assert matrix.matrix_name == "BLOSUM"
        assert matrix.version == 62

    def test_select_matrix_invalid(self):
        matrix = _Matrix()
        with pytest.raises(ValueError):
            matrix.select_matrix("INVALID", 999)

    def test_str_representation(self):
        matrix = _Matrix()
        assert str(matrix) == "No matrix selected"
        matrix.select_matrix("Blosum", 62)
        assert str(matrix) == "BLOSUM62 Matrix"


# Blosum Matrix Tests
class TestBLOSUM:
    def test_blosum_init(self, blosum62):
        assert blosum62.matrix_name == "BLOSUM"
        assert blosum62.version == 62
        assert blosum62.matrix_data is not None

    def test_blosum_invalid_version(self):
        with pytest.raises(ValueError):
            Blosum(999)

    def test_blosum_lookups(self, blosum62):
        assert blosum62["A"]["A"] > 0  # Match should be positive
        assert isinstance(blosum62["A"]["W"], int)  # Should return integer score

    def test_blosum_invalid_lookup(self, blosum62):
        with pytest.raises(KeyError):
            blosum62["J"]

    @pytest.mark.parametrize("version", [45, 50, 62, 80, 90])
    def test_blosum_versions(self, version):
        matrix = Blosum(version)
        assert matrix.version == version
        assert matrix["A"]["A"] is not None


# Pam Matrix Tests
class TestPAM:
    def test_pam_init(self, pam250):
        assert pam250.matrix_name == "PAM"
        assert pam250.version == 250
        assert pam250.matrix_data is not None

    def test_pam_invalid_version(self):
        with pytest.raises(ValueError):
            Pam(999)

    def test_pam_lookups(self, pam250):
        assert isinstance(pam250["A"]["A"], int)
        assert isinstance(pam250["W"]["C"], int)

    @pytest.mark.parametrize("version", [30, 70, 250])
    def test_pam_versions(self, version):
        matrix = Pam(version)
        assert matrix.version == version
        assert matrix["A"]["A"] is not None


# Identity Matrix Tests
class TestIDENTITY:
    def test_identity_init(self, identity):
        assert identity.matrix_name == "IDENTITY"
        assert identity.version == 0
        assert identity.matrix_data is not None

    def test_identity_matches(self, identity):
        assert identity["A"]["A"] == 1  # Match
        assert identity["A"]["C"] == 0  # Mismatch

    @pytest.mark.parametrize(
        "aa1,aa2,expected", [("A", "A", 1), ("C", "C", 1), ("A", "C", 0), ("W", "Y", 0)]
    )
    def test_identity_scores(self, identity, aa1, aa2, expected):
        assert identity[aa1][aa2] == expected


# Match Matrix Tests
class TestMATCH:
    def test_match_init(self, match):
        assert match.matrix_name == "MATCH"
        assert match.matrix_data is not None

    def test_match_scores(self, match):
        assert match["A"]["A"] == 1  # Match
        assert match["A"]["C"] == -1  # Mismatch

    @pytest.mark.parametrize(
        "aa1,aa2,expected",
        [("A", "A", 1), ("C", "C", 1), ("A", "C", -1), ("W", "Y", -1)],
    )
    def test_match_all_combinations(self, match, aa1, aa2, expected):
        assert match[aa1][aa2] == expected


# File Loading Tests
class TestFileLoading:
    def test_missing_file(self, tmp_path):
        matrix = _Matrix(matrix_folder=tmp_path)
        matrix.select_matrix("Blosum", 62)
        with pytest.raises(RuntimeError):
            matrix.load_json_matrix()

    def test_invalid_json(self, tmp_path):
        # Create invalid JSON file
        file_path = tmp_path / "BLOSUM62.json"
        file_path.write_text("invalid json")

        matrix = _Matrix(matrix_folder=tmp_path)
        matrix.select_matrix("Blosum", 62)
        with pytest.raises(Exception):  # Could be JSON decode error
            matrix.load_json_matrix()


# Integration Tests
class TestIntegration:
    def test_matrix_chain_access(self, blosum62):
        # Test chained access works correctly
        score = blosum62["A"]["W"]
        assert isinstance(score, int)

    def test_matrix_symmetry(self, blosum62):
        # Test matrix is symmetric
        assert blosum62["A"]["W"] == blosum62["W"]["A"]

    @pytest.mark.parametrize("matrix_class", [Blosum, Pam, Identity, Match])
    def test_matrix_initialization(self, matrix_class):
        print(matrix_class.__name__)
        if matrix_class in [Blosum, Pam]:
            matrix = matrix_class(matrix_class.matrices[matrix_class.__name__][0])
        elif matrix_class == Identity:
            matrix = matrix_class(0)
        else:
            matrix = matrix_class()
        assert matrix.matrix_data is not None


class TestAPIConsistency:
    """Ensure consistent behavior across all matrix types"""

    def test_common_interface(self, blosum62, pam250, identity, match):
        """All matrices should support the same basic operations"""
        for matrix_obj in [blosum62, pam250, identity, match]:
            # Basic attributes
            assert hasattr(matrix_obj, "matrix_name")
            assert hasattr(matrix_obj, "matrix_data")

            # Common methods
            assert hasattr(matrix_obj, "available_matrices")
            assert callable(getattr(matrix_obj, "available_matrices"))

            # Dictionary-like access
            assert matrix_obj["A"]["A"] is not None

            # String representation
            assert str(matrix_obj).endswith("Matrix")

    def test_error_handling_consistency(self, blosum62, pam250, identity, match):
        """All matrices should handle errors consistently"""
        for matrix_obj in [blosum62, pam250, identity, match]:
            # Invalid amino acid
            with pytest.raises(KeyError):
                matrix_obj["J"]

            # Invalid types
            with pytest.raises(Exception):  # Could be TypeError or ValueError
                matrix_obj[1]

            with pytest.raises(Exception):
                matrix_obj[None]

    def test_chaining_consistency(self, blosum62, pam250, identity, match):
        """All matrices should handle chained access consistently"""
        matrices = [blosum62, pam250, identity, match]
        for matrix in matrices:
            # First access should return matrix-like object
            assert hasattr(matrix["A"], "__getitem__")
            # Second access should return score
            assert isinstance(matrix["A"]["A"], int)


class TestBioinformaticsIntegration:
    """Test integration with common bioinformatics workflows"""

    def test_pairwise_alignment_scoring(self, blosum62):
        """Test using matrix for pairwise alignment scoring"""
        seq1 = "ARND"
        seq2 = "ARCD"

        # Calculate alignment score
        score = sum(blosum62[a][b] for a, b in zip(seq1, seq2))
        assert isinstance(score, int)

        # Verify specific alignment scores
        assert score > blosum62["N"]["C"]  # Mismatch score
        assert score < sum(blosum62[a][a] for a in seq1)  # Perfect match score

    def test_multiple_sequence_scoring(self, blosum62):
        """Test scoring multiple sequence alignments"""
        sequences = ["ARND", "ARCD", "ARAD"]

        # Score all sequences against first sequence
        scores = []
        for seq in sequences[1:]:
            score = sum(blosum62[a][b] for a, b in zip(sequences[0], seq))
            scores.append(score)

        assert all(isinstance(score, int) for score in scores)
        assert len(scores) == len(sequences) - 1

    def test_conservation_analysis(self, blosum62, pam250, identity, match):
        """Test using matrices for sequence conservation analysis"""
        for matrix_obj in [blosum62, pam250, identity, match]:
            alignment = [
                "ARNDCEQ",
                "ARNDCEQ",  # Identical
                "ARNDCEK",  # One difference
                "ARNDCWA",  # Multiple differences
            ]

            # Calculate conservation scores
            scores = []
            for pos in range(len(alignment[0])):
                column = [seq[pos] for seq in alignment]
                # Score each position against consensus (first sequence)
                col_scores = [matrix_obj[column[0]][aa] for aa in column[1:]]
                scores.append(sum(col_scores))

            assert len(scores) == len(alignment[0])
            assert scores[0] > scores[-1]  # First position more conserved than last

    def test_substitution_analysis(self, blosum45, pam250):
        """Test analyzing amino acid substitutions"""
        # Compare scoring between matrices
        aa_pairs = [("A", "G"), ("W", "F"), ("D", "E")]

        for aa1, aa2 in aa_pairs:
            blosum_score = blosum45[aa1][aa2]
            pam_score = pam250[aa1][aa2]

            # Scores should correlate (similar amino acids score higher in both matrices)
            print(blosum_score, pam_score)
            assert (blosum_score >= 0) == (pam_score >= 0)

    def test_profile_generation(self, blosum62):
        """Test generating position-specific scoring profiles"""
        sequence = "ARNDCEQ"
        profile = {}

        # Generate simple position-specific scoring profile
        for pos, aa in enumerate(sequence):
            profile[pos] = {
                target_aa: blosum62[aa][target_aa]
                for target_aa in "ACDEFGHIKLMNPQRSTVWY"
            }

        assert len(profile) == len(sequence)
        assert all(len(pos_scores) == 20 for pos_scores in profile.values())
        assert all(
            isinstance(score, int)
            for pos_scores in profile.values()
            for score in pos_scores.values()
        )
