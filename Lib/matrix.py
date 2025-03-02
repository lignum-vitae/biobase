import json
from pathlib import Path

class Matrix:
    matrices = {
           "BLOSUM": [45, 50, 62, 80, 90], 
           "PAM": [30, 70, 250]
           }
    def __init__(self)->None:
        self.folder = Path("./matrices")
        self.matrix_data = None
        self.matrix_name = None
        self.version = None

    @classmethod
    def available_matrices(cls)->list[str]:
        return [f"{name}{version}" for name, version_list in cls.matrices.items() for version in version_list]

    def select_matrix(self, matrix_name: str, version: int)->None:
        matrix_name = matrix_name.upper()
        available = self.available_matrices()
        if f"{matrix_name}{version}" in available:
            self.matrix_name = matrix_name
            self.version = version
            return
        raise ValueError(f"Only {', '.join(available[:-1])}, and {available[-1]} are currently supported matrices")

    def load_json_matrix(self)->None:
        """
        Change this to appropriate directory later
        """
        filename = f"{self.matrix_name}{self.version}.json"
        json_file_path  = self.folder / filename

        if not json_file_path.exists():
            raise RuntimeError("File not found")

        with open(json_file_path) as file:
            self.matrix_data = json.load(file)

    def __getitem__(self, key: str):
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
        return f"{self.matrix_name}{self.version} Matrix" if self.matrix_name else "No matrix selected"

class BaseMatrixClass(Matrix):
    def __init__(self, matrix_name: str, version: int) -> None:
        super().__init__()
        self.select_matrix(matrix_name, version)
        self.load_json_matrix()

class BLOSUM(BaseMatrixClass):
    def __init__(self, version: int) -> None:
        super().__init__("BLOSUM", version)

class PAM(BaseMatrixClass):
    def __init__(self, version: int) -> None:
        super().__init__("PAM", version)

if __name__ == "__main__":
    blosum = BLOSUM(62)
    pam = PAM(250)

    print(blosum["A"]["A"])
    print(blosum)
    print(pam["A"]["A"])
    print(pam)
    print(Matrix.available_matrices())
