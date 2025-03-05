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
    def _validate_nucleotide(nucs: str, is_single_nucleotide: bool) -> None:
        if is_single_nucleotide and len(nucs) != 1:
            raise ValueError("Only a single nucleotide character can be analyzed.")
        if not isinstance(nucs, str):
            raise ValueError("Only strings can be analyzed.")
        nucs = nucs.upper()
        if any(nuc not in 'ATCGU' for nuc in nucs):
            raise ValueError("Only valid nucleotides (A, T, C, G, U) can be analyzed.")
        return nucs

    @classmethod
    def molecular_weight(cls, nuc: str) -> float:
        if not nuc:
            return None
        nuc = cls._validate_nucleotide(nuc, is_single_nucleotide=True)
        return cls.nuc_molecular_weight[nuc]

    @classmethod
    def  cumulative_molecular_weight(cls, nucs: str) -> float:
        if not nucs:
            return None
        nucs = cls._validate_nucleotide(nucs, is_single_nucleotide=False)
        return sum(cls.nuc_molecular_weight[nuc] for nuc in nucs)

class Dna:
    complements = COMPLEMENTS

    @staticmethod
    def _validate_dna_sequence(dna_seq: str) -> None:
        """Helper method to validate if a DNA sequence is valid."""
        if not isinstance(dna_seq, str):
            raise ValueError("Only strings can be analyzed.")
        dna_seq = dna_seq.upper()
        if not all(letter in 'ATCG' for letter in dna_seq):
            raise ValueError("Only valid DNA nucleotides (A, T, C, G) can be used.")
        return dna_seq


    @classmethod
    def transcribe(cls, dna_seq: str) -> str:
        if not dna_seq:
            return None
        dna_seq = cls._validate_dna_sequence(dna_seq)
        return dna_seq.replace('T', 'U')

    @classmethod
    def complement_dna(cls, dna_seq: str, reverse: bool = True) -> str:
        if not dna_seq:
            return None
        dna_seq = cls._validate_dna_sequence(dna_seq)
        if not reverse:
            return "".join([cls.complements[x] for x in dna_seq])
        return "".join([cls.complements[x] for x in dna_seq[::-1]])

if __name__ == "__main__":
    print(Nucleotides.cumulative_molecular_weight("acgtagtgat"))
    print(Nucleotides.molecular_weight("u"))

    print(Dna.transcribe("acccggtccatcatcattca"))
    print(Dna.complement_dna("acccggtccatcatcattca"))
    print(Dna.complement_dna("acccggtccatcatcattca", reverse=False))
