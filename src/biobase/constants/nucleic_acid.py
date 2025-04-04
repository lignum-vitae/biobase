"""A collection of DNA and RNA constants.

Public module variables:

DNA -- a string containing DNA nucleotides
RNA -- a string containing RNA nucleotides
NUCLEOTIDES -- a string containing all DNA and RNA nucleotides
NUCLEOTIDE_NAMES -- a list containing the name of each nucleotide
MOLECULAR_WEIGHT -- a dictionary containing the molecular weights for each nucleotide
DNA_COMPLEMENTS -- a dictionary containing complements for each DNA nucleotide
RNA_COMPLEMENTS -- a dictionary containing complements for each RNA nucleotide
IUPAC_NUCLEOTIDES -- a dictionary containing the nucleotide equivalent for each IUPAC code
"""

DNA = 'ATCG'
RNA = 'AUCG'
NUCLEOTIDES = 'ATCGU'
NUCLEOTIDE_NAMES = ["Adenine", "Thymine", "Cytosine", "Guanine", "Uracil"]
MOLECULAR_WEIGHT = {"A":135.13, "T":126.12, "C":111.10, "G":151.13, "U":112.09}

DNA_COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C"}
RNA_COMPLEMENTS = {"A":"U", "U":"A", "C":"G", "G":"C"}

# [4]
IUPAC_NUCLEOTIDES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    "R": "[AG]",  # Purine
    "Y": "[CT]",  # Pyrimidine
    "M": "[AC]",  # Amino
    "K": "[GT]",  # Keto
    "S": "[CG]",  # Strong
    "W": "[AT]",  # Weak
    "H": "[ACT]", # not G
    "B": "[CGT]", # not A
    "V": "[ACG]", # not T
    "D": "[AGT]", # not C
    "N": "[ACGTU]" # any
}

# Melting temperature constants
GC_CONTENT_FACTOR = 41  # For GC content calculation
SALT_CORRECTION_FACTOR = 16.6  # For salt correction

# Common sequence motifs in IUPAC notation
# Lower case letters are variable. Upper case letters are conserved.
KOZAK_SEQUENCE = "gccgccRccAUGG"  # Consensus Kozak sequence. First gcc is of unknown significance [5]
TATA_BOX = "TATAWAW"       # TATA box consensus [6]
POLY_A_SIGNAL = "AAUAAA"  # Polyadenylation signal [7]

# DNA/RNA structure
PURINE_BASES = "AG"
PYRIMIDINE_BASES = "CTU"
WEAK_PAIRS = "AT"    # 2 hydrogen bonds
STRONG_PAIRS = "GC"  # 3 hydrogen bonds

# Codon tables
START_CODONS = {"ATG"}  # Standard start codon
STOP_CODONS = {"TAA", "TAG", "TGA"}  # Standard stop codons

# Common restriction sites
RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC"
}
