# standard library
import math
import re
from typing import Iterator

# internal dependencies
from biobase.constants.nucleic_acid import DNA_COMPLEMENTS, MOLECULAR_WEIGHT


def main():
    print(Nucleotides.cumulative_molecular_weight("acgtagtgat"))
    print(Nucleotides.cumulative_molecular_weight("AATTGG"))
    print(Nucleotides.cumulative_molecular_weight("ACTG"))
    print(Nucleotides.cumulative_molecular_weight("AU"))
    print(Nucleotides.molecular_weight("u"))

    print(Dna.transcribe("acccggtccatcatcattca"))
    print(Dna.complement("acccggtccatcatcattca"))
    print(Dna.complement("acccggtccatcatcattca", reverse=False))
    # Homopolymer, exterme low entropy case, exp value = 0.0
    print(Dna.entropy("AAAAAAA"))
    print(Dna.entropy("ACGTACGT"))  # Equal distribution, exp value = 2.0
    print(Dna.entropy("AAACCCGG"))  # Mixed, exp value = 1.561278124459133

    for s, e in Dna.find_orfs("ccatgccctaaatggggtag") :
        print(s, e)


class Nucleotides:
    VALID_NUCLEOTIDES = frozenset(
        "ATCGU"
    )  # frozenset more efficient for repeated lookups
    nuc_molecular_weight = MOLECULAR_WEIGHT

    @staticmethod
    def _validate_nucleotide(nucs: str, is_single_nucleotide: bool = False) -> str:
        """
        Validate a nucleotide sequence or single nucleotide.

        Parameters:
        - nucs (str): Nucleotide sequence to validate
        - is_single_nucleotide (bool): If True, ensures input is exactly one nucleotide

        Returns:
        - str: Uppercase validated sequence

        Raises:
        - ValueError: If sequence is invalid or empty

        Example:
        >>> Nucleotides._validate_nucleotide("atcgu")
        'ATCGU'
        >>> Nucleotides._validate_nucleotide("X")  # raises ValueError
        ValueError: Invalid nucleotides found: ['X']
        """
        if not nucs:
            raise ValueError("Empty nucleotide sequence provided.")
        if not isinstance(nucs, str):
            raise ValueError(f"Expected string input, got {type(nucs).__name__}")
        if is_single_nucleotide and len(nucs) != 1:
            raise ValueError(
                "Expected single nucleotide, got sequence of length {len(nucs)}"
            )

        nucs = nucs.upper()
        invalids = set(nucs) - Nucleotides.VALID_NUCLEOTIDES
        if invalids:
            raise ValueError(f"Invalid nucleotides found: {invalids}")

        return nucs

    @classmethod
    def molecular_weight(cls, nuc: str) -> float:
        """
        Calculate molecular weight of a single nucleotide.

        Parameters:
        - nuc (str): Single nucleotide character

        Returns:
        - float: Molecular weight in g/mol

        Raises:
        - ValueError: If input is invalid or empty

        Example:
        >>> Nucleotides.molecular_weight("A")
        135.13
        >>> Nucleotides.molecular_weight("u")
        112.09
        """
        nuc = cls._validate_nucleotide(nuc, is_single_nucleotide=True)
        return cls.nuc_molecular_weight[nuc]

    @classmethod
    def cumulative_molecular_weight(cls, nucs: str) -> float:
        """
        Calculate cumulative molecular weight of a nucleotide sequence.

        Parameters:
        - nucs (str): Nucleotide sequence

        Returns:
        - float: Total molecular weight in g/mol

        Raises:
        - ValueError: If input is invalid or empty

        Example:
        >>> Nucleotides.cumulative_molecular_weight("ATCG")
        523.48
        >>> Nucleotides.cumulative_molecular_weight("AU")
        247.22
        """
        nucs = cls._validate_nucleotide(nucs)
        return sum(cls.nuc_molecular_weight[nuc] for nuc in nucs)


class Dna:
    VALID_DNA = frozenset("ATCG")
    complements = DNA_COMPLEMENTS

    @staticmethod
    def _validate_dna_sequence(dna_seq: str) -> str:
        """
        Validate a DNA sequence.

        Parameters:
        - dna_seq (str): DNA sequence to validate

        Returns:
        - str: Uppercase validated sequence

        Raises:
        - ValueError: If sequence is invalid or empty

        Example:
        >>> Dna._validate_dna_sequence("atcg")
        'ATCG'
        >>> Dna._validate_dna_sequence("ATCU")  # raises ValueError
        ValueError: Invalid DNA nucleotides found: ['U']
        """
        if not dna_seq:
            raise ValueError("Empty sequence provided.")
        if not isinstance(dna_seq, str):
            raise ValueError(f"Expected string input, got {type(dna_seq).__name__}")

        dna_seq = dna_seq.upper()
        invalids = set(dna_seq) - Dna.VALID_DNA
        if invalids:
            raise ValueError(f"Invalid DNA nucleotides found: {sorted(invalids)}")

        return dna_seq

    @classmethod
    def transcribe(cls, dna_seq: str) -> str:
        """
        Transcribe DNA to RNA.

        Parameters:
        - dna_seq (str): DNA sequence

        Returns:
        - str: RNA sequence

        Raises:
        - ValueError: If sequence is invalid or empty

        Example:
        >>> Dna.transcribe("ATCG")
        'AUCG'
        >>> Dna.transcribe("atcg")
        'AUCG'
        """
        dna_seq = cls._validate_dna_sequence(dna_seq)
        return dna_seq.replace("T", "U")

    @classmethod
    def complement(cls, dna_seq: str, reverse: bool = False) -> str:
        """
        Generate DNA complement sequence.

        Parameters:
        - dna_seq (str): DNA sequence
        - reverse (bool): If True, return reverse complement

        Returns:
        - str: Complement sequence

        Raises:
        - ValueError: If sequence is invalid or empty

        Example:
        >>> Dna.complement("ATCG")
        'CGAT'  # reverse complement
        >>> Dna.complement("ATCG", reverse=False)
        'TAGC'  # complement only
        """
        dna_seq = cls._validate_dna_sequence(dna_seq)
        complement = "".join(cls.complements[x] for x in dna_seq)
        return complement[::-1] if reverse else complement

    @classmethod
    def calculate_gc_content(cls, dna_seq: str) -> float:
        """
        Calculate the GC content percentage of a DNA sequence.

        This function calculates the percentage of G and C nucleotides in a DNA sequence.
        The sequence is converted to uppercase before calculation.

        Parameters:
        - sequence (str): A DNA sequence string

        Returns:
        - float: The percentage of GC content (0-100)

        Raises:
        - ValueError: If sequence is invalid or empty

        Example:
        >>> calculate_gc_content("ATGC")
        50.0
        >>> calculate_gc_content("GCGC")
        100.0
        """
        sequence = cls._validate_dna_sequence(dna_seq)
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100 if sequence else 0.0

    @classmethod
    def calculate_at_content(cls, dna_seq: str) -> float:
        """
        Calculate the AT content percentage of a DNA sequence.

        This function calculates the percentage of A and T nucleotides in a DNA sequence.
        The sequence is converted to uppercase before calculation.

        Parameters:
        - sequence (str): A DNA sequence string

        Returns:
        - float: The percentage of AT content (0-100)

        Raises:
        - ValueError: If sequence is invalid or empty

        Example:
        >>> calculate_at_content("ATGC")
        50.0
        >>> calculate_at_content("ATAT")
        100.0
        """
        sequence = cls._validate_dna_sequence(dna_seq)
        at_count = sequence.count("A") + sequence.count("T")
        return (at_count / len(sequence)) * 100 if sequence else 0.0

    @classmethod
    def entropy(cls, dna_sequence: str) -> float:
        """
        Calculate the Shannon entropy of a DNA sequence.

        This function computes the Shannon entropy, which measures the
        uncertainty or randomness in the nucleotide composition of a DNA sequence.
        Shannon entropy is always non-negative for discrete distributions,
        with 0 indicating a completely uniform sequence
        and higher values indicating more diversity.

        Parameters:
        - dna_sequence (str): A DNA sequence string

        Returns:
        - float: The Shannon entropy of the sequence

        Raises:
        - ValueError: If the sequence is invalid or empty

        Example:
        >>> entropy("AAAAAAA")
        0.0
        >>> entropy("ACGTACGT")
        2.0
        >>> entropy("AAACCCGG")
        1.561278124459133
        """
        dna_sequence = cls._validate_dna_sequence(dna_sequence)
        # Calculate the proportion of bases in the sequence
        counts = [dna_sequence.count(nuc) for nuc in cls.VALID_DNA]
        seq_length = len(dna_sequence)
        entropy = 0
        for count in counts:
            if count > 0:  # Avoid log(0) because it is undefined
                p = count / seq_length
                entropy -= p * math.log2(p)
        return entropy

    # Compile once and use it many times 
    _ORF_PATTERN: re.Pattern = re.compile(
        r'atg(?:[atgc]{3})*?(?:taa|tag|tga)',
        re.IGNORECASE
    )
    @classmethod
    # Using iterator here is more memory efficient, instead of collecting we stream the results
    def find_orfs(cls, dna_sequence : str) -> Iterator[tuple[int, int]] :
        """
        Yeilds a stream of ORFs found in DNA sequence.
        ORF or open reading frame is defined as a spans of DNA sequences
        between the start (ATG) and stop codons (TAA | TAG | TGA).
        Parameters:
        - dna_sequence (str) 
        Yeilds:
        A tuple(start, end) for each ORF found. Indices are 0 based and 
        follows Python slice convention (start is inclusive, end is exclusive)
        Example:
            >>> for s, e in Dna.find_orfs("ccatgccctaaatggggtag") :
                    print(s,e)
            2 11
            11 18
            optional, print the sequence as well
            >>> seq = "ccatgccctaaatggggtag"
            >>> for s, e in Dna.find_orfs(seq) :
                    print(s,e)
            2 11 atgccctaa
            11 18 atggggtag
        No exceptions are raised by this function itself. However, validation
        errors might be raised by cls._validate_dna_sequence if the input is invalid
        """
        dna_sequence = cls._validate_dna_sequence(dna_sequence)
        for m in cls._ORF_PATTERN.finditer(dna_sequence) :
            # return one value at a time
            yield m.span()
        

if __name__ == "__main__":
    main()
