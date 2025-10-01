import json
import importlib.resources
import pytest
from biobase.matrix import _Matrix


MATRIX_FILES = [
    "BLOSUM62.json",
    "IDENTITY0.json",
    "MATCH.json",
    "PAM50.json",
]


@pytest.mark.parametrize("filename", MATRIX_FILES)
def test_load_json_from_package(filename):
    package = "biobase.matrix.matrices"
    try:
        with importlib.resources.files(package).joinpath(filename).open("r") as file:
            data = json.load(file)
    except FileNotFoundError:
        pytest.fail(f"File not found in package resources: {filename}")
    except Exception as e:
        pytest.fail(f"Failed to load JSON {filename}: {e}")

    # Basic sanity check: data should be a dict
    assert isinstance(data, dict), f"JSON content is not a dict for {filename}"

    # Optional: check some expected keys, e.g. amino acid letters like 'A'
    assert "A" in data, f"'A' key missing in JSON matrix {filename}"


@pytest.mark.parametrize("filename", MATRIX_FILES)
def test_matrix_json_file_exists_and_loads(filename):
    package = "biobase.matrix.matrices"
    file = importlib.resources.files(package).joinpath(filename)
    assert file.is_file(), f"{filename} not found in package"


@pytest.mark.parametrize("filename", MATRIX_FILES)
def test_matrix_structure_and_symmetry(filename):
    package = "biobase.matrix.matrices"
    file = importlib.resources.files(package).joinpath(filename)

    with file.open("r") as f:
        matrix = json.load(f)

    # Check that all keys are single-character amino acids
    for key in matrix:
        assert isinstance(key, str) and len(key) == 1, f"Invalid key: {key}"
        for inner_key, score in matrix[key].items():
            assert isinstance(inner_key, str)
            assert isinstance(score, (int, float))

    # Check symmetry: score[A][B] == score[B][A]
    for a in matrix:
        for b in matrix[a]:
            assert b in matrix
            assert a in matrix[b]
            assert matrix[a][b] == matrix[b][a], f"Asymmetry at {a}, {b}"


def test_matrix_class_loads_matrix():
    matrix = _Matrix()
    matrix.select_matrix("BLOSUM", 62)
    matrix.load_json_matrix()

    assert matrix.matrix_data is not None
    assert isinstance(matrix.matrix_data, dict)
    assert "A" in matrix.matrix_data


def test_matrix_raises_on_missing_file():
    matrix = _Matrix()
    with pytest.raises(ValueError):
        matrix.select_matrix("FAKE", 999)
