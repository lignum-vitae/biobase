"""A collection of DNA and RNA constants.

Public module variables:

DNA -- a string containing DNA nucleotides
RNA -- a string containing RNA nucleotides
COMPLEMENTS -- a dictionary containing complements for each DNA nucleotide
NUCLEOTIDE_NAMES -- a list containing the name of each nucleotide
MOLECULAR_WEIGHT -- a dictionary containing the molecular weights for each nucleotide
"""

DNA = 'ATCG'
RNA = 'AUCG'
NUCLEOTIDES = 'ATCGU'
COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C"}
NUCLEOTIDE_NAMES = ["Adenine", "Thymine", "Cytosine", "Guanine", "Uracil"]
MOLECULAR_WEIGHT = {"A":135.13, "T":126.12, "C":111.10, "G":151.13, "U":112.09}

class Nucleotides():
    nuc_molecular_weight = MOLECULAR_WEIGHT

    @staticmethod
    def _validate_nucleotide(nucs: str, is_single_nucleotide: bool) -> str:
        if is_single_nucleotide and len(nucs) != 1:
            raise ValueError("Only a single nucleotide character can be analyzed.")
        if not isinstance(nucs, str):
            raise ValueError("Only strings can be analyzed.")
        nucs = nucs.upper()
        if any(nuc not in 'ATCGU' for nuc in nucs):
            raise ValueError("Only valid nucleotides (A, T, C, G, U) can be analyzed.")
        return nucs

    @classmethod
    def molecular_weight(cls, nuc: str) -> float | None:
        """
        Returns the molecular weight of a single nucleotide.

        This method retrieves the molecular weight of a nucleotide from the `MOLECULAR_WEIGHT` dictionary.
        The nucleotide is validated before being used.

        Parameters:
        - nuc (str): A single nucleotide character (A, T, C, G, or U).

        Returns:
        - float: The molecular weight of the nucleotide.
        - none: If there is no input sequence.

        Raises:
        - ValueError: If the nucleotide is not valid or if input sequence is not a string.
        - ValueError: If input string is longer than one character.

        Example:
        >>> Nucleotides.molecular_weight("A")
        135.13
        >>> Nucleotides.molecular_weight("G")
        151.13
        >>> Nucleotides.molecular_weight("X")
        Traceback (most recent call last):
        ValueError: Only valid nucleotides (A, T, C, G, U) can be analyzed.
        """
        if not nuc:
            return None
        nuc = cls._validate_nucleotide(nuc, is_single_nucleotide=True)
        return cls.nuc_molecular_weight[nuc]

    @classmethod
    def  cumulative_molecular_weight(cls, nucs: str) -> float | None:
        """
        Returns the cumulative molecular weight of a nucleotide sequence.

        This method computes the sum of the molecular weights for each nucleotide in the provided sequence.
        The sequence is validated before being processed.

        Parameters:
        - nucs (str): A string representing a sequence of nucleotides (e.g., "ACGT").

        Returns:
        - float: The cumulative molecular weight of the nucleotide sequence.
        - none: If there is no input sequence.

        Raises:
        - ValueError: If the nucleotide is not valid or if input sequence is not a string.

        Example:
        >>> Nucleotides.cumulative_molecular_weight("ACGT")
        523.48
        >>> Nucleotides.cumulative_molecular_weight("AATTGG")
        824.76
        >>> Nucleotides.cumulative_molecular_weight("AXT")
        Traceback (most recent call last):
        ValueError: Only valid nucleotides (A, T, C, G, U) can be analyzed.
        """
        if not nucs:
            return None
        nucs = cls._validate_nucleotide(nucs, is_single_nucleotide=False)
        return sum(cls.nuc_molecular_weight[nuc] for nuc in nucs)

class Dna:
    complements = COMPLEMENTS

    @staticmethod
    def _validate_dna_sequence(dna_seq: str) -> str:
        """Helper method to validate if a DNA sequence is valid."""
        if not isinstance(dna_seq, str):
            raise ValueError("Only strings can be analyzed.")
        dna_seq = dna_seq.upper()
        if not all(letter in 'ATCG' for letter in dna_seq):
            raise ValueError("Only valid DNA nucleotides (A, T, C, G) can be used.")
        return dna_seq

    @classmethod
    def transcribe(cls, dna_seq: str) -> str | None:
        """
        Transcribes a DNA sequence into RNA.

        This method converts a given DNA sequence into an RNA sequence by replacing 'T' with 'U'.
        The DNA sequence is validated before the transcription is performed.

        Parameters:
        - dna_seq (str): A string representing the DNA sequence to be transcribed.

        Returns:
        - str: The transcribed RNA sequence, where 'T' is replaced by 'U'.
        - none: If there is no input sequence.

        Raises:
        - ValueError: If the DNA sequence is invalid or if the DNA sequence is not a string.

        Example:
        >>> Dna.transcribe("ATCG")
        'AUCG'
        >>> Dna.transcribe("ATCGATCG")
        'AUCGAUCG'
        >>> Dna.transcribe("ATCGX")
        Traceback (most recent call last):
        ValueError: Only valid DNA nucleotides (A, T, C, G) can be used.
        """
        if not dna_seq:
            return None
        dna_seq = cls._validate_dna_sequence(dna_seq)
        return dna_seq.replace('T', 'U')

    @classmethod
    def complement_dna(cls, dna_seq: str, reverse: bool = True) -> str | None:
        """
        Returns the complement of a DNA sequence.

        This method returns the complement of the given DNA sequence. The default behavior is to reverse the 
        sequence before complementing, but this can be changed by setting `reverse` to `False`.

        Parameters:
        - dna_seq (str): A string representing the DNA sequence to be complemented.
        - reverse (bool): Flag indicating whether to reverse the DNA sequence before returning the complement 
                          (default is `True`).

        Returns:
        - str: The complement of the DNA sequence, optionally reversed.
        - none: If there is no input sequence.

        Raises:
        - ValueError: If the DNA sequence is invalid or if the DNA sequence is not a string.

        Example:
        >>> Dna.complement_dna("ATCG")
        'CGAT'
        >>> Dna.complement_dna("ATCG", reverse=False)
        'TAGC'
        >>> Dna.complement_dna("ATCGX")
        Traceback (most recent call last):
        ValueError: Only valid DNA nucleotides (A, T, C, G) can be used.
        """
        if not dna_seq:
            return None
        dna_seq = cls._validate_dna_sequence(dna_seq)
        if not reverse:
            return "".join([cls.complements[x] for x in dna_seq])
        return "".join([cls.complements[x] for x in dna_seq[::-1]])

if __name__ == "__main__":
    print(Nucleotides.cumulative_molecular_weight("acgtagtgat"))
    print(Nucleotides.cumulative_molecular_weight("AATTGG"))
    print(Nucleotides.cumulative_molecular_weight("ACTG"))
    print(Nucleotides.molecular_weight("u"))

    print(Dna.transcribe("acccggtccatcatcattca"))
    print(Dna.complement_dna("acccggtccatcatcattca"))
    print(Dna.complement_dna("acccggtccatcatcattca", reverse=False))
