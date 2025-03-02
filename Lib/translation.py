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
one_letter_codes_ext = one_letter_codes + 'OU'

three_letter_codes = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
                      "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
                      "Tyr", "Val"]
three_letter_codes_ext = ["Pyl", "Sec"]
three_letter_codes_ext.extend(three_letter_codes)

amino_acid_names = ["Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", 
                    "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                    "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"]

amino_acid_names_ext = ["Pyrrolysine", "Selenocysteine"]
amino_acid_names_ext.extend(amino_acid_names)

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

mono_mass_ext ={"O":237.147726925, "U":150.953633405}
mono_mass_ext.update(mono_mass)

def motif_finder(protein: str, pattern: str, ext: bool = False) -> list[int]:
    """
    :param protein: a string of amino acids
    :param pattern: python-flavoured regex string for desired pattern or plain text for exact matches
    :return: List of all start positions of motif (including overlapping positions)
    """
    """
    Find all occurrences of a specified motif (pattern) in a given protein sequence.

    This function uses regular expressions to locate all positions of a motif (including overlapping matches) 
    within a protein sequence. It raises a ValueError if the protein sequence contains invalid characters 
    (i.e., characters not belonging to the standard amino acid one-letter codes).

    Parameters:
    - protein (str): A string representing the protein sequence, where each character is an amino acid's 
                     one-letter code (e.g., "ACDEFGHIKLMNPQRSTVWY").
    - pattern (str): A string representing the motif (pattern) to search for in the protein sequence.
                     It can be a plain text for exact matches or a Python-flavoured regular expression 
                     string to specify complex patterns.

    Returns:
    - list[int]: A list of 1-based start positions (indices) where the motif is found in the protein sequence.
                 If no motif matches, an empty list is returned.

    Raises:
    - ValueError: If the protein sequence contains characters not present in the standard amino acid one-letter codes.
    - ValueError: If any parameter is left blank
    - ValueError: If pattern is not a stirng

    Example:
    >>> motif_finder("ACDEFGHIKLMNPQRSTVWY", "CDE")
    [2]
    >>> motif_finder("ACDEFGHIKLMNPQRSTVWY", "A.*F")
    [1]
    """
    aa_codes = one_letter_codes if not ext else one_letter_codes_ext
    aa_codes = set(aa_codes)
    if not protein:
        raise ValueError("The input fasta dictionary is empty. Please provide valid data.")
    if not isinstance(pattern, str) or not pattern:
        raise ValueError("The pattern must be a non-empty string.")
    if not all(aa in aa_codes for aa in protein):
        raise ValueError("Invalid protein sequence used in motif finder")
    return [m.start(0)+1 for m in re.finditer(f"(?={pattern})", protein)] #positive lookahead match is zero-width

def fasta_motif_finder(fasta_dict: dict[str, str], pattern: str, ext: bool = False) -> dict[str, list[str]]:
    """
    :param fasta_dict: dictionary of fasta values from get_uniprot_fasta
    :param pattern: python-flavoured regex string for desired pattern or plain text for exact matches
    :return: dictionary containing fasta_id and list of all positions matching motif (including overlapping positions)
    """
    """
    Find all occurrences of a specified motif (pattern) in protein sequences from a FASTA dictionary.

    This function searches for a motif in each protein sequence in a FASTA-format dictionary, where keys 
    represent protein identifiers (e.g., Uniprot IDs) and values are protein sequences. It returns a 
    dictionary of protein identifiers as keys and a list of start positions (1-based) of the motif as values. 
    If no motif is found in a particular sequence, that sequence is omitted from the result.

    Parameters:
    - fasta_dict (dict[str, str]): A dictionary where keys are protein identifiers (e.g., Uniprot IDs) 
                                   and values are protein sequences (strings) containing one-letter amino acid codes.
    - pattern (str): A string representing the motif (pattern) to search for in the protein sequences. 
                     It can be a Python-flavoured regular expression string or plain text for exact matches.

    Returns:
    - dict[str, list[str]]: A dictionary with protein identifiers as keys and a list of start positions 
                             (1-based) where the motif is found in the protein sequence as values. Only sequences 
                             with at least one match are included in the output.

    Raises:
    - ValueError: If the protein sequence contains characters not present in the standard amino acid one-letter codes.
    - ValueError: If any parameter is left blank
    - ValueError: If pattern is not a stirng

    Example:
    >>> fasta_dict = {"P12345": "ACDEFGHIKLMNPQRSTVWY", "Q67890": "NNNNNNNNNNNN"}
    >>> fasta_motif_finder(fasta_dict, "CDE")
    {"P12345": [2]}
    >>> fasta_motif_finder(fasta_dict, "A.*F")
    {"P12345": [1]}
    """
    aa_codes = one_letter_codes if not ext else one_letter_codes_ext
    aa_codes = set(aa_codes)
    if not fasta_dict:
        raise ValueError("The input fasta dictionary is empty. Please provide valid data.")
    if not isinstance(pattern, str) or not pattern:
        raise ValueError("The pattern must be a non-empty string.")
    for fasta_id, value in fasta_dict.items():
        # Ensure that the sequence in fasta_dict only contains valid amino acids
        invalid_aa = [aa for aa in value if aa not in aa_codes]
        if invalid_aa:
            raise ValueError(f"Invalid amino acids found in the sequence for {fasta_id}: {', '.join(invalid_aa)}.")

    result_dict = {}
    for fasta_id, value in fasta_dict.items():
        matches = [str(m.start(0)+1) for m in re.finditer(f"(?={pattern})", value)] #positive lookahead match is zero-width
        if matches: #only adds sequences with desired motif to output dictionary
            result_dict[fasta_id] = matches
    return result_dict
