"""A collection of protein constants.


Public module variables:


one_letter_codes -- a string containing all amino acid one letter codes
three_letter_codes -- a list containing all amino acid three letter codes
amino_acid_names -- a list containing all amino acid names
codons_per_aa -- a dictionary containing the number of codons that can result in specified amino acid
codon_table -- a dictionary containing all RNA codon to amino acid conversions
codons -- a list containing all codons
mono_mass --  a dictionary containing the mass of all amino acids
*_ext -- same as the above with the addition of selenocysteine and pyrrolysine
"""

from itertools import product
import re

one_letter_codes = 'ACDEFGHIKLMNPQRSTVWY'
three_letter_codes = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
                      "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
                      "Tyr", "Val"]
amino_acid_names = ["Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", 
                    "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                    "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"]

one_letter_codes_ext = 'ACDEFGHIKLMNOPQRSTUVWY'
three_letter_codes_ext = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", 
                          "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Pyl", 
                          "Sec", "Ser", "Thr", "Trp", "Tyr", "Val"]
amino_acid_names_ext = ["Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", 
                        "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                        "Phenylalanine", "Proline", "Pyrrolysine", "Selenocysteine", "Serine", "Threonine", 
                        "Tryptophan", "Tyrosine", "Valine"]

codons = ["".join(x) for x in product('AUCG', repeat=3)]
codons_per_aa ={'A': 4, 'W': 1, 'D': 2, 'F': 2, 'S': 6, 'Stop': 3, 'P': 4,
                'E': 2, 'V': 4, 'R': 6, 'M': 1, 'Q': 2, 'L': 6, 'K': 2, 
                'I': 3, 'C': 2, 'H': 2, 'T': 4, 'G': 4, 'N': 2, 'Y': 2}

codon_table =  {"UUU":"F",    "UCU":"S",    "UAU":"Y",    "UGU":"C",
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

mono_mass ={"A": 71.037113805, "C":103.009184505, "D":115.026943065, "E":129.042593135,
            "F":147.068413945, "G": 57.021463735, "H":137.058911875, "I":113.084064015,
            "K":128.094963050, "L":113.084064015, "M":131.040484645, "N":114.042927470,
            "P": 97.052763875, "Q":128.058577540, "R":156.101111050, "S": 87.032028435,
            "T":101.047678505, "V": 99.068413945, "W":186.079312980, "Y":163.063328575}

mono_mass_ext ={"A": 71.037113805, "C":103.009184505, "D":115.026943065, "E":129.042593135,
                "F":147.068413945, "G": 57.021463735, "H":137.058911875, "I":113.084064015,
                "K":128.094963050, "L":113.084064015, "M":131.040484645, "N":114.042927470,
                "O":237.147726925, "P": 97.052763875, "Q":128.058577540, "R":156.101111050,
                "S": 87.032028435, "T":101.047678505, "U":150.953633405, "V": 99.068413945,
                "W":186.079312980, "Y":163.063328575}

def motif_finder(protein: str, pattern: str) -> list[int]:
    """
    :param protein: a string of amino acids
    :param pattern: python-flavoured regex string for desired pattern or plain text for exact matches
    :return: List of all start positions of motif (including overlapping positions)
    """
    aa_codes = 'ACDEFGHIKLMNPQRSTVWY'
    if not all(aa in aa_codes for aa in protein):
        raise ValueError("Invalid protein sequence used in motif finder")
    return [m.start(0)+1 for m in re.finditer(f"(?={pattern})", protein)] #positive lookahead match is zero-width

def fasta_motif_finder(fasta_dict: dict[str, str], pattern: str) -> dict[str, list[str]]:
    """
    :param fasta_dict: dictionary of fasta values from get_uniprot_fasta
    :param pattern: python-flavoured regex string for desired pattern or plain text for exact matches
    :return: dictionary containing uniprot_id and list of all positions matching motif (including overlapping positions)
    """
    result_dict = {}
    for uniprot_id, value in fasta_dict.items():
        matches = [str(m.start(0)+1) for m in re.finditer(f"(?={pattern})", value)] #positive lookahead match is zero-width
        if matches: #only adds sequences with desired motif to output dictionary
            result_dict[uniprot_id] = matches
    return result_dict
