"""
Nucleic Acid Constants and Methods

This module provides fundamental constants and methods for DNA and RNA analysis.

Public module variables:
DNA -- string containing DNA nucleotides (ATCG)
RNA -- string containing RNA nucleotides (AUCG)
NUCLEOTIDES -- string containing all DNA and RNA nucleotides (ATCGU)
NUCLEOTIDE_NAMES -- list containing the name of each nucleotide
MOLECULAR_WEIGHT -- dictionary containing molecular weights for each nucleotide
DNA_COMPLEMENTS -- dictionary containing complements for each DNA nucleotide
RNA_COMPLEMENTS -- dictionary containing complements for each RNA nucleotide
IUPAC_NUCLEOTIDES -- dictionary containing nucleotide patterns for IUPAC codes
SEQUENCE_MOTIFS -- dictionary containing common sequence motifs and their metadata
PURINE_BASES -- string containing purine bases
PYRIMIDINE_BASES -- string containing pyrimidine bases
WEAK_PAIRS -- base pairs with 2 hydrogen bonds
STRONG_PAIRS -- base pairs with 3 hydrogen bonds
START_CODONS -- set of standard start codons
STOP_CODONS -- set of standard stop codons
RESTRICTION_SITES -- dictionary of common restriction enzyme recognition sequences

Key components:
1. Nucleotide sequences (DNA, RNA)
2. Molecular properties (molecular weights, complements)
3. IUPAC codes and their meanings
4. Common sequence motifs and regulatory elements
5. Structural properties (base pairing, GC content)
6. Restriction sites and enzymes

Constants are organized by:
- Basic nucleotides (DNA, RNA, NUCLEOTIDES)
- Molecular properties (MOLECULAR_WEIGHT)
- Complementarity (DNA_COMPLEMENTS, RNA_COMPLEMENTS)
- IUPAC nomenclature (IUPAC_NUCLEOTIDES)
- Sequence motifs (KOZAK_SEQUENCE, TATA_BOX, POLY_A_SIGNAL)
- Structural features (PURINE_BASES, PYRIMIDINE_BASES)
- Genetic elements (START_CODONS, STOP_CODONS)
- Restriction sites (RESTRICTION_SITES)

Notes:
- IUPAC codes allow for ambiguous base representation
- Some sequence motifs may vary between organisms
- Restriction sites are typically palindromic
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

# Common sequence motifs in IUPAC notation
SEQUENCE_MOTIFS = {
    'kozak': { # [5]
        'sequence': 'gccgccRccAUGG', # First gcc is of unknown significance
        'description': 'Consensus Kozak sequence for translation initiation',
    },
    'tata_box': { # [6]
        'sequence': 'TATAWAW',
        'description': 'Core promoter element for transcription',
    },
    'polya_signal': { # [7]
        'sequence': 'AAUAAA',
        'description': 'Polyadenylation signal sequence',
    },
    'splice_donor': { # [9]
        'sequence': 'MAG|GTRAGT',
        'description': '5\' splice site consensus',
    },
    'splice_acceptor': { # [9]
        'sequence': 'CAG|G',
        'description': '3\' splice site consensus',
    }
}

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
