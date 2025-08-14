import json
import importlib.resources
from pathlib import Path


def main():
    blosum = Blosum(62)
    pam = Pam(250)
    identity0 = Identity(0)
    match = Match()

    print(blosum["A"]["A"])
    print(blosum["A"])
    print(blosum)
    print(pam["A"]["A"])
    print(pam)
    print(identity0["A"]["A"])
    print(identity0["C"]["A"])
    print(identity0)
    print(match["A"]["A"])
    print(match["C"]["A"])
    print(match)
    print(Match.available_matrices())
    print(blosum.available_matrices())


class _Matrix:
    matrices = {
        "Blosum": [30, 35, 40, 45, 50, 55, 60, 62, 65, 70, 75, 80, 85, 90],
        "Pam": [
            10,
            30,
            50,
            60,
            70,
            80,
            90,
            100,
            110,
            120,
            160,
            170,
            200,
            250,
            300,
            400,
            450,
            500,
        ],
        "Identity": [-10000, 0],
        "Match": [""],
    }
    # Get the project root directory (src/biobase)
    PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()

    # Common subdirectories
    # Ensure load_json function matches this folder layout
    MATRICES_DIR = PROJECT_ROOT / "biobase" / "matrix" / "matrices"

    def __init__(self, matrix_folder: str | Path | None = None) -> None:
        """
        Initialize a Matrix object with a specified matrix folder.

        Parameters:
        - matrix_folder (str | Path): Path to the folder containing matrix files.
                                     Defaults to 'PROJECT_ROOT/src/biobase/matrix/matrices'

        Returns:
        - None

        Example:
        >>> matrix = _Matrix()  # Uses default folder
        >>> matrix = _Matrix("path/to/matrices")  # Uses custom folder
        """
        self.folder = (
            self.MATRICES_DIR if matrix_folder is None else Path(matrix_folder)
        )
        self.matrix_data = None
        self.matrix_name = None
        self.version = None

    @classmethod
    def available_matrices(cls) -> list[str]:
        """
        Get a list of all available scoring matrices.

        Returns:
        - list[str]: List of available matrices in format "NAME{version}"

        Example:
        >>> _Matrix.available_matrices()
        ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']
        """
        return [
            f"{name.upper()}{version}"
            for name, version_list in cls.matrices.items()
            for version in version_list
        ]

    def select_matrix(self, matrix_name: str, version: int | str) -> None:
        """
        Select a specific scoring matrix by name and version.

        Parameters:
        - matrix_name (str): Name of the matrix (e.g., "Blosum", "Pam")
        - version (int): Version number of the matrix

        Raises:
        - ValueError: If the requested matrix/version combination is not available

        Example:
        >>> matrix = _Matrix()
        >>> matrix.select_matrix("Blosum", 62)  # Selects BLOSUM62
        >>> matrix.select_matrix("Pam", 999)  # raises ValueError
        ValueError: Only BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, and PAM250 are currently supported matrices
        """
        matrix_name = matrix_name.upper()
        available = self.available_matrices()
        if f"{matrix_name}{version}" in available:
            self.matrix_name = matrix_name
            self.version = version
            return
        raise ValueError(
            f"Only {', '.join(available[:-1])}, and {available[-1]} are currently supported matrices"
        )

    def load_json_matrix(self) -> None:
        """
        Load the selected matrix data from its JSON file.

        The matrix must be selected using select_matrix()

        Raises:
        - RuntimeError: If the matrix file is not found
        - ValueError: If no matrix has been selected

        Example:
        >>> matrix = _Matrix()
        >>> matrix.select_matrix("Blosum", 62)
        >>> matrix.load_json_matrix()  # Loads BLOSUM62.json
        """
        filename = f"{self.matrix_name}{self.version}.json"
        if self.folder is not None:
            # Load from local filesystem path (dev)
            json_file_path = self.folder / filename

            if not json_file_path.exists():
                raise RuntimeError(f"File not found: {json_file_path}")

            with open(json_file_path) as file:
                self.matrix_data = json.load(file)
        else:
            # Load from package resources (pip installation)
            package = "biobase.matrix.matrices"
            with (
                importlib.resources.files(package).joinpath(filename).open("r") as file
            ):
                self.matrix_data = json.load(file)

    def __getitem__(self, key: str):
        """
        Access matrix values using dictionary-style lookup.

        Allows for chained indexing to access nested values in the matrix.

        Parameters:
        - key (str): Matrix key (typically an amino acid letter)

        Returns:
        - Union[dict, int]: Either a sub-matrix (dict) or score value (int)

        Raises:
        - ValueError: If no matrix data is loaded
        - KeyError: If the key is not found in the matrix

        Example:
        >>> matrix = _Matrix()
        >>> matrix.select_matrix("Blosum", 62)
        >>> matrix.load_json_matrix()
        >>> matrix["A"]["A"]  # Get score for A-A match
        4
        """
        if not self.matrix_data:
            raise ValueError("No matrix data loaded")

        if key not in self.matrix_data:
            raise KeyError(f"Key '{key}' not found in matrix")

        return self.matrix_data[key]

    def __str__(self):
        """
        Get a string representation of the matrix.

        Returns:
        - str: String in format "NAME{version} Matrix" or "No matrix selected"

        Example:
        >>> matrix = _Matrix()
        >>> str(matrix)
        'No matrix selected'
        >>> matrix.select_matrix("Blosum", 62)
        >>> str(matrix)
        'BLOSUM62 Matrix'
        """
        return (
            f"{self.matrix_name}{self.version} Matrix"
            if self.matrix_name
            else "No matrix selected"
        )


class _BaseMatrixClass(_Matrix):
    def __init__(self, matrix_name: str, version: int | str, matrix_folder) -> None:
        super().__init__(matrix_folder)
        self.select_matrix(matrix_name, version)
        self.load_json_matrix()


class Blosum(_BaseMatrixClass):
    def __init__(self, version: int, matrix_folder: str | None = None) -> None:
        """
        Initialize a Blosum (BLOcks SUbstitution Matrix) scoring matrix.

        Blosum matrices are amino acid substitution matrices based on observed alignments.
        Higher numbers (e.g., BLOSUM80) are designed for comparing closely related sequences,
        while lower numbers (e.g., BLOSUM45) are for more divergent sequences.

        Parameters:
        - version (int): Blosum version number (ex. 45, 50, 62, 80, or 90)
        - matrix_folder (str | Path): Path to matrix files. Defaults to Matrix.default_matrix_folder

        Raises:
        - ValueError: If version is not one of the available version numbers
        - RuntimeError: If matrix file is not found

        Example:
        >>> blosum62 = Blosum(62)
        >>> blosum62["A"]["A"]  # Score for matching Alanine-Alanine
        4
        >>> blosum62["W"]["C"]  # Score for substituting Tryptophan with Cysteine
        -2
        """
        super().__init__("Blosum", version, matrix_folder)


class Pam(_BaseMatrixClass):
    def __init__(self, version: int, matrix_folder: str | None = None) -> None:
        """
        Initialize a Pam (Point Accepted Mutation) scoring matrix.

        Pam matrices are amino acid substitution matrices based on evolutionary distance.
        Lower numbers (e.g., PAM30) are for closely related sequences,
        while higher numbers (e.g., PAM250) are for more divergent sequences.

        Parameters:
        - version (int): Pam version number ex. (30, 70, or 250)
        - matrix_folder (str | Path): Path to matrix files. Defaults to Matrix.default_matrix_folder

        Raises:
        - ValueError: If version is not one of the available version numbers
        - RuntimeError: If matrix file is not found

        Example:
        >>> pam250 = Pam(250)
        >>> pam250["A"]["A"]  # Score for matching Alanine-Alanine
        2
        >>> pam250["W"]["C"]  # Score for substituting Tryptophan with Cysteine
        -8
        """
        super().__init__("Pam", version, matrix_folder)


class Identity(_BaseMatrixClass):
    def __init__(self, version: int, matrix_folder: str | None = None) -> None:
        """
        Initialize an Identity scoring matrix.

        The Identity matrix is a simple scoring matrix that gives:
        - A positive score of typically 1 for matching amino acids
        - A negative score of typically 0 for mismatching amino acids
        This matrix treats all mismatches equally, unlike Blosum or Pam matrices.

        Parameters:
        - version (int): Identity version number (-10000, or 0)
        - matrix_folder (str | Path): Path to matrix files. Defaults to Matrix.default_matrix_folder

        Raises:
        - RuntimeError: If matrix file is not found

        Example:
        >>> identity = Identity(0)
        >>> identity["A"]["A"]  # Score for matching Alanine-Alanine
        1
        >>> identity["W"]["C"]  # Score for any mismatch
        0
        """
        super().__init__("Identity", version, matrix_folder)


class Match(_BaseMatrixClass):
    def __init__(self, matrix_folder: str | None = None) -> None:
        """
        Initialize a Match scoring matrix.

        The Match matrix is a binary scoring matrix that gives:
        - A score of +1 for matching amino acids
        - A score of -1 for mismatching amino acids
        This matrix is useful for calculating sequence identity and in cases
        where you only want to count exact matches.

        Parameters:
        - matrix_folder (str | Path): Path to matrix files. Defaults to Matrix.default_matrix_folder

        Raises:
        - RuntimeError: If matrix file is not found

        Example:
        >>> match = Match()
        >>> match["A"]["A"]  # Score for matching Alanine-Alanine
        1
        >>> match["W"]["C"]  # Score for any mismatch
        -1
        """
        super().__init__("Match", "", matrix_folder)


def text_matrix_to_json(
    input_matrix_path: str | Path, output_matrix_path: str | Path, matrix_name: str
) -> None:
    r"""
    Convert a text matrix file to JSON format.

    Parameters:
    - input_matrix_path (str | Path): Path to the input text matrix file
    - output_matrix_path (str | Path): Path to the output JSON matrix file
    - matrix_name (str): Name of matrix

    Prints:
    - FileNotFoundError: If the matrix file is not found
    - ValueError: If the file path is empty
    - File Successfully created: If file is successfully created

    Example:
    >>> # Get the project root directory (src/biobase)
    >>> PROJECT_ROOT = Path(__file__).parent.parent.resolve()

    >>> # Common subdirectories
    >>> MATRICES_DIR = PROJECT_ROOT / "biobase" / "matrices"
    >>> TEXT_MATRICES_DIR = MATRICES_DIR / "text_matrices"
    >>> chosen_matrix = "PAM70"
    >>> matrix_input = TEXT_MATRICES_DIR / chosen_matrix
    >>> matrix_output = MATRICES_DIR / chosen_matrix
    >>> text_matrix_to_json(matrix_input, matrix_output, chosen_matrix)
    File successfully created: JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\PAM70.json

    >>> matrices_to_download = ["PAM300", "PAM400", "PAM450", "BLOSUM30", "BLOSUM85"]
    >>> for matrix in matrices_to_download:
    >>>     chosen_matrix = matrix
    >>>     matrix_input = TEXT_MATRICES_DIR / chosen_matrix
    >>>     matrix_output = MATRICES_DIR / chosen_matrix
    >>>     text_matrix_to_json(matrix_input, matrix_output, chosen_matrix)
    File Successfully created: JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\PAM300.json
    File Successfully created: JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\PAM400.json
    File Successfully created: JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\PAM450.json
    File Successfully created: JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\BLOSUM30.json
    File successfully created: JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\BLOSUM85.json
    """

    if not input_matrix_path:
        print("ValueError: Empty file path provided.")

    input_path = Path(f"{input_matrix_path}.txt")
    if not input_path.exists():
        print(f"FileNotFoundError: Matrix file not found: {input_path}")

    with open(input_path) as input_file:
        raw_lines = input_file.readlines()

    # Filter out comments and split lines into tokens
    matrix_lines = [line.split() for line in raw_lines if not line.startswith("#")]
    # First line contains amino acid labels
    amino_acid_labels = matrix_lines[0]

    scoring_matrix = {}
    for _, row_data in enumerate(matrix_lines[1:]):
        # First element is the row label, rest are scores
        row_label = row_data[0]
        row_scores = [int(score) for score in row_data[1:]]
        # Create dictionary mapping amino acids to their scores
        scoring_matrix[row_label] = dict(zip(amino_acid_labels, row_scores))

    output_path = Path(f"{output_matrix_path}.json")
    with open(output_path, "w") as output_file:
        json.dump(scoring_matrix, output_file, indent=4)

    print(f"File Successfully created: JSON file created at: {output_path}")


if __name__ == "__main__":
    main()
