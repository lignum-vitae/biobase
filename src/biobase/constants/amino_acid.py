"""
Amino Acid Constants and Methods

This module provides fundamental constants and methods for protein analysis.

Public module variables:
ONE_LETTER_CODES -- string containing all amino acid one letter codes (ACDEFGHIKLMNPQRSTVWY)
ONE_LETTER_CODES_EXT -- same as above plus selenocysteine (U) and pyrrolysine (O)
THREE_LETTER_CODES -- set containing all amino acid three letter codes
THREE_LETTER_CODES_EXT -- same as above plus selenocysteine (Sec) and pyrrolysine (Pyl)
AMINO_ACID_NAMES -- set containing all amino acid full names
AMINO_ACID_NAMES_EXT -- same as above plus selenocysteine and pyrrolysine
MONO_MASS -- dictionary containing monoisotopic masses for each amino acid
MONO_MASS_EXT -- same as above plus selenocysteine and pyrrolysine
CODONS_PER_AA -- dictionary containing the number of codons per amino acid
CODON_TABLE -- dictionary containing RNA codon to amino acid conversions
CODONS -- list containing all possible RNA codons

Key components:
1. Amino acid identifiers (one letter, three letter, full names)
2. Extended amino acids (selenocysteine, pyrrolysine)
3. Physical properties (monoisotopic masses)
4. Genetic code (codons, codon usage)

Constants are organized by:
- Standard amino acid codes (ONE_LETTER_CODES, THREE_LETTER_CODES, AMINO_ACID_NAMES)
- Extended amino acids (*_EXT versions include pyrrolysine and selenocysteine)
- Physical properties (MONO_MASS, MONO_MASS_EXT)
- Genetic code (CODONS_PER_AA, CODON_TABLE, CODONS)

Notes:
- Masses are monoisotopic (most abundant isotope)
- Extended amino acids (U, O) are rare but naturally occurring
- Codon usage varies between organisms
- Stop codons are represented as 'Stop' in the codon table
"""

from itertools import product

ONE_LETTER_CODES = "ACDEFGHIKLMNPQRSTVWY"
ONE_LETTER_CODES_EXT = ONE_LETTER_CODES + "OU"

THREE_LETTER_CODES = {
    "Ala",
    "Arg",
    "Asn",
    "Asp",
    "Cys",
    "Glu",
    "Gln",
    "Gly",
    "His",
    "Ile",
    "Leu",
    "Lys",
    "Met",
    "Phe",
    "Pro",
    "Ser",
    "Thr",
    "Trp",
    "Tyr",
    "Val",
}
THREE_LETTER_CODES_EXT = {"Pyl", "Sec"}
THREE_LETTER_CODES_EXT.update(THREE_LETTER_CODES)

AMINO_ACID_NAMES = {
    "Alanine",
    "Arginine",
    "Asparagine",
    "Aspartic acid",
    "Cysteine",
    "Glutamic acid",
    "Glutamine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Leucine",
    "Lysine",
    "Methionine",
    "Phenylalanine",
    "Proline",
    "Serine",
    "Threonine",
    "Tryptophan",
    "Tyrosine",
    "Valine",
}

AMINO_ACID_NAMES_EXT = {"Pyrrolysine", "Selenocysteine"}
AMINO_ACID_NAMES_EXT.update(AMINO_ACID_NAMES)

# [1]
MONO_MASS = {
    "A": 71.037113805,
    "C": 103.009184505,
    "D": 115.026943065,
    "E": 129.042593135,
    "F": 147.068413945,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084064015,
    "K": 128.094963050,
    "L": 113.084064015,
    "M": 131.040484645,
    "N": 114.042927470,
    "P": 97.052763875,
    "Q": 128.058577540,
    "R": 156.101111050,
    "S": 87.032028435,
    "T": 101.047678505,
    "V": 99.068413945,
    "W": 186.079312980,
    "Y": 163.063328575,
}

MONO_MASS_EXT = {"O": 237.147726925, "U": 150.953633405}
MONO_MASS_EXT.update(MONO_MASS)

CODONS = ["".join(x) for x in product("AUCG", repeat=3)]
CODONS_PER_AA = {
    "A": 4,
    "W": 1,
    "D": 2,
    "F": 2,
    "S": 6,
    "Stop": 3,
    "P": 4,
    "E": 2,
    "V": 4,
    "R": 6,
    "M": 1,
    "Q": 2,
    "L": 6,
    "K": 2,
    "I": 3,
    "C": 2,
    "H": 2,
    "T": 4,
    "G": 4,
    "N": 2,
    "Y": 2,
}

# [2]
CODON_TABLE = {
    "UUU": "F",
    "UCU": "S",
    "UAU": "Y",
    "UGU": "C",
    "UUC": "F",
    "UCC": "S",
    "UAC": "Y",
    "UGC": "C",
    "UUA": "L",
    "UCA": "S",
    "UAA": "STOP",
    "UGA": "STOP",
    "UUG": "L",
    "UCG": "S",
    "UAG": "STOP",
    "UGG": "W",
    "CUU": "L",
    "CCU": "P",
    "CAU": "H",
    "CGU": "R",
    "CUC": "L",
    "CCC": "P",
    "CAC": "H",
    "CGC": "R",
    "CUA": "L",
    "CCA": "P",
    "CAA": "Q",
    "CGA": "R",
    "CUG": "L",
    "CCG": "P",
    "CAG": "Q",
    "CGG": "R",
    "AUU": "I",
    "ACU": "T",
    "AAU": "N",
    "AGU": "S",
    "AUC": "I",
    "ACC": "T",
    "AAC": "N",
    "AGC": "S",
    "AUA": "I",
    "ACA": "T",
    "AAA": "K",
    "AGA": "R",
    "AUG": "M",
    "ACG": "T",
    "AAG": "K",
    "AGG": "R",
    "GUU": "V",
    "GCU": "A",
    "GAU": "D",
    "GGU": "G",
    "GUC": "V",
    "GCC": "A",
    "GAC": "D",
    "GGC": "G",
    "GUA": "V",
    "GCA": "A",
    "GAA": "E",
    "GGA": "G",
    "GUG": "V",
    "GCG": "A",
    "GAG": "E",
    "GGG": "G",
}
