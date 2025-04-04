"""A collection of protein constants.

Public module variables:

ONE_LETTER_CODES -- a string containing all amino acid one letter codes
THREE_LETTER_CODES -- a list containing all amino acid three letter codes
AMINO_ACID_NAMES -- a list containing all amino acid names
CODONS_PER_AA -- a dictionary containing the number of CODONS that can result in specified amino acid
CODON_TABLE -- a dictionary containing all RNA codon to amino acid conversions
CODONS -- a list containing all CODONS
MONO_MASS --  a dictionary containing the mass of all amino acids
*_ext -- same as the above with the addition of selenocysteine and pyrrolysine
"""
from itertools import product

ONE_LETTER_CODES = 'ACDEFGHIKLMNPQRSTVWY'
ONE_LETTER_CODES_EXT = ONE_LETTER_CODES + 'OU'

THREE_LETTER_CODES = {"Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
                      "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
                      "Tyr", "Val"}
THREE_LETTER_CODES_EXT = {"Pyl", "Sec"}
THREE_LETTER_CODES_EXT.update(THREE_LETTER_CODES)

AMINO_ACID_NAMES = {"Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", 
                    "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                    "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"}

AMINO_ACID_NAMES_EXT = {"Pyrrolysine", "Selenocysteine"}
AMINO_ACID_NAMES_EXT.update(AMINO_ACID_NAMES)

# [1]
MONO_MASS ={"A": 71.037113805, "C":103.009184505, "D":115.026943065, "E":129.042593135,
            "F":147.068413945, "G": 57.021463735, "H":137.058911875, "I":113.084064015,
            "K":128.094963050, "L":113.084064015, "M":131.040484645, "N":114.042927470,
            "P": 97.052763875, "Q":128.058577540, "R":156.101111050, "S": 87.032028435,
            "T":101.047678505, "V": 99.068413945, "W":186.079312980, "Y":163.063328575}

CODONS = ["".join(x) for x in product('AUCG', repeat=3)]
CODONS_PER_AA ={'A': 4, 'W': 1, 'D': 2, 'F': 2, 'S': 6, 'Stop': 3, 'P': 4,
                'E': 2, 'V': 4, 'R': 6, 'M': 1, 'Q': 2, 'L': 6, 'K': 2, 
                'I': 3, 'C': 2, 'H': 2, 'T': 4, 'G': 4, 'N': 2, 'Y': 2}

# [2]
CODON_TABLE =  {"UUU":"F",    "UCU":"S",    "UAU":"Y",    "UGU":"C",
                "UUC":"F",    "UCC":"S",    "UAC":"Y",    "UGC":"C",
                "UUA":"L",    "UCA":"S",    "UAA":"STOP", "UGA":"STOP",
                "UUG":"L",    "UCG":"S",    "UAG":"STOP", "UGG":"W",
                "CUU":"L",    "CCU":"P",    "CAU":"H",    "CGU":"R",
                "CUC":"L",    "CCC":"P",    "CAC":"H",    "CGC":"R",
                "CUA":"L",    "CCA":"P",    "CAA":"Q",    "CGA":"R",
                "CUG":"L",    "CCG":"P",    "CAG":"Q",    "CGG":"R",
                "AUU":"I",    "ACU":"T",    "AAU":"N",    "AGU":"S",
                "AUC":"I",    "ACC":"T",    "AAC":"N",    "AGC":"S",
                "AUA":"I",    "ACA":"T",    "AAA":"K",    "AGA":"R",
                "AUG":"M",    "ACG":"T",    "AAG":"K",    "AGG":"R",
                "GUU":"V",    "GCU":"A",    "GAU":"D",    "GGU":"G",
                "GUC":"V",    "GCC":"A",    "GAC":"D",    "GGC":"G",
                "GUA":"V",    "GCA":"A",    "GAA":"E",    "GGA":"G",
                "GUG":"V",    "GCG":"A",    "GAG":"E",    "GGG":"G"}

MONO_MASS_EXT ={"O":237.147726925, "U":150.953633405}
MONO_MASS_EXT.update(MONO_MASS)
