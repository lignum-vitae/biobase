import json
from pathlib import Path

class _Matrix:
    matrices = {
           "BLOSUM": [45, 50, 62, 80, 90],
           "PAM": [30, 70, 250]
           }
    default_matrix_folder = Path("src/biobase/matrices").resolve()
    def __init__(self, matrix_folder: str = None)->None:
        """
        Initialize a Matrix object with a specified matrix folder.

        Parameters:
        - matrix_folder (str | Path): Path to the folder containing matrix files.
                                     Defaults to './matrices'

        Returns:
        - None

        Example:
        >>> matrix = _Matrix()  # Uses default folder
        >>> matrix = _Matrix("path/to/matrices")  # Uses custom folder
        """
        self.folder = self.default_matrix_folder if matrix_folder == None else matrix_folder
        self.matrix_data = None
        self.matrix_name = None
        self.version = None

    @classmethod
    def available_matrices(cls)->list[str]:
        """
        Get a list of all available scoring matrices.

        Returns:
        - list[str]: List of available matrices in format "NAME{version}"

        Example:
        >>> _Matrix.available_matrices()
        ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250']
        """
        return [f"{name}{version}" for name, version_list in cls.matrices.items() for version in version_list]

    def select_matrix(self, matrix_name: str, version: int)->None:
        """
        Select a specific scoring matrix by name and version.

        Parameters:
        - matrix_name (str): Name of the matrix (e.g., "BLOSUM", "PAM")
        - version (int): Version number of the matrix

        Raises:
        - ValueError: If the requested matrix/version combination is not available

        Example:
        >>> matrix = _Matrix()
        >>> matrix.select_matrix("BLOSUM", 62)  # Selects BLOSUM62
        >>> matrix.select_matrix("PAM", 999)  # raises ValueError
        ValueError: Only BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, and PAM250 are currently supported matrices
        """
        matrix_name = matrix_name.upper()
        available = self.available_matrices()
        if f"{matrix_name}{version}" in available:
            self.matrix_name = matrix_name
            self.version = version
            return
        raise ValueError(f"Only {', '.join(available[:-1])}, and {available[-1]} are currently supported matrices")

    def load_json_matrix(self)->None:
        """
        Load the selected matrix data from its JSON file.

        The matrix must be selected using select_matrix() before calling this method.

        Raises:
        - RuntimeError: If the matrix file is not found
        - ValueError: If no matrix has been selected

        Example:
        >>> matrix = _Matrix()
        >>> matrix.select_matrix("BLOSUM", 62)
        >>> matrix.load_json_matrix()  # Loads BLOSUM62.json
        """
        filename = f"{self.matrix_name}{self.version}.json"
        json_file_path  = self.folder / filename

        if not json_file_path.exists():
            raise RuntimeError("File not found")

        with open(json_file_path) as file:
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
        >>> matrix.select_matrix("BLOSUM", 62)
        >>> matrix.load_json_matrix()
        >>> matrix["A"]["A"]  # Get score for A-A match
        4
        """
        if not self.matrix_data:
            raise ValueError("No matrix data loaded")


        if key not in self.matrix_data:
            raise KeyError(f"Key '{key}' not found in matrix")

        # Handle chained indexing: if the value is a dictionary, return another Matrix-like object
        sub_matrix = self.matrix_data[key]
        if isinstance(sub_matrix, dict):
            # Return the current matrix with the sub-matrix data, simulating the chained access
            new_matrix = self.__class__(self.version)
            new_matrix.matrix_data = sub_matrix
            new_matrix.matrix_name = self.matrix_name
            new_matrix.version = self.version
            return new_matrix
        return sub_matrix

    def __str__(self):
        """
        Get a string representation of the matrix.

        Returns:
        - str: String in format "NAME{version} Matrix" or "No matrix selected"

        Example:
        >>> matrix = _Matrix()
        >>> str(matrix)
        'No matrix selected'
        >>> matrix.select_matrix("BLOSUM", 62)
        >>> str(matrix)
        'BLOSUM62 Matrix'
        """
        return f"{self.matrix_name}{self.version} Matrix" if self.matrix_name else "No matrix selected"

class _BaseMatrixClass(_Matrix):
    def __init__(self, matrix_name: str, version: int, matrix_folder) -> None:
        super().__init__(matrix_folder)
        self.select_matrix(matrix_name, version)
        self.load_json_matrix()

class BLOSUM(_BaseMatrixClass):
    def __init__(self, version: int, matrix_folder: str = None) -> None:
        """
        Initialize a BLOSUM (BLOcks SUbstitution Matrix) scoring matrix.

        BLOSUM matrices are amino acid substitution matrices based on observed alignments.
        Higher numbers (e.g., BLOSUM80) are designed for comparing closely related sequences,
        while lower numbers (e.g., BLOSUM45) are for more divergent sequences.

        Parameters:
        - version (int): BLOSUM version number (45, 50, 62, 80, or 90)
        - matrix_folder (str | Path): Path to matrix files. Defaults to Matrix.default_matrix_folder

        Raises:
        - ValueError: If version is not one of: 45, 50, 62, 80, 90
        - RuntimeError: If matrix file is not found

        Example:
        >>> blosum62 = BLOSUM(62)
        >>> blosum62["A"]["A"]  # Score for matching Alanine-Alanine
        4
        >>> blosum62["W"]["C"]  # Score for substituting Tryptophan with Cysteine
        -2
        """
        super().__init__("BLOSUM", version, matrix_folder)

class PAM(_BaseMatrixClass):
    def __init__(self, version: int, matrix_folder: str = None) -> None:
        """
        Initialize a PAM (Point Accepted Mutation) scoring matrix.

        PAM matrices are amino acid substitution matrices based on evolutionary distance.
        Lower numbers (e.g., PAM30) are for closely related sequences,
        while higher numbers (e.g., PAM250) are for more divergent sequences.

        Parameters:
        - version (int): PAM version number (30, 70, or 250)
        - matrix_folder (str | Path): Path to matrix files. Defaults to Matrix.default_matrix_folder

        Raises:
        - ValueError: If version is not one of: 30, 70, 250
        - RuntimeError: If matrix file is not found

        Example:
        >>> pam250 = PAM(250)
        >>> pam250["A"]["A"]  # Score for matching Alanine-Alanine
        2
        >>> pam250["W"]["C"]  # Score for substituting Tryptophan with Cysteine
        -8
        """
        super().__init__("PAM", version, matrix_folder)

def text_matrix_to_json(input_matrix_path: str | Path, output_matrix_path: str | Path, matrix_name: str) -> None:
    """
    Convert a text matrix file to JSON format.

    Parameters:
    - input_matrix_path (str | Path): Path to the input text matrix file
    - output_matrix_path (str | Path): Path to the output JSON matrix file
    - matrix_name (str): Name of matrix

    Raises:
    - FileNotFoundError: If the matrix file is not found
    - ValueError: If the file path is empty

    Example:
    >>> chosen_matrix = "PAM70"
    >>> matrix_input = f"src/biobase/matrices/text_matrices/{chosen_matrix}"
    >>> matrix_output = f"src/biobase/matrices/{chosen_matrix}"
    >>> text_matrix_to_json(matrix_input, matrix_output, chosen_matrix)
    JSON file created at: C:\REST\OF\ABSOLUTE\PATH\src\biobase\matrices\PAM70.json
    """
    if not input_matrix_path:
        raise ValueError("Empty file path provided.")

    input_path = Path(f"{input_matrix_path}.txt").resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Matrix file not found: {input_path}")

    with open(input_path) as input_file:
        raw_lines = input_file.readlines()

    # Filter out comments and split lines into tokens
    matrix_lines = [line.split() for line in raw_lines if not line.startswith("#")]
    # First line contains amino acid labels
    amino_acid_labels = matrix_lines[0]

    scoring_matrix = {}
    for row_index, row_data in enumerate(matrix_lines[1:]):
        # First element is the row label, rest are scores
        row_label = row_data[0]
        row_scores = [int(score) for score in row_data[1:]]
        # Create dictionary mapping amino acids to their scores
        scoring_matrix[row_label] = dict(zip(amino_acid_labels, row_scores))

    output_path = Path(f"{output_matrix_path}.json").resolve()
    with open(output_path, 'w') as output_file:
        json.dump(scoring_matrix, output_file, indent=4)

    print(f"JSON file created at: {output_path}")

if __name__ == "__main__":
    blosum = BLOSUM(62)
    pam = PAM(250)

    print(blosum["A"]["A"])
    print(blosum)
    print(pam["A"]["A"])
    print(pam)
    print(_Matrix.available_matrices())
