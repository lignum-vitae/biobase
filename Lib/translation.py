"""A collection of protein constants.

Public module variables:

one_letter_codes -- a string containing all amino acid one letter codes
three_letter_codes -- a list containing all amino acid three letter codes
amino_acids -- a list containing all all amino acid names
codons_per_aa -- a dictionary containing the number of codons that can result in specified amino acid
codon_table -- a dictionary containing all RNA to amino acid conversions
"""
one_letter_codes = 'ACDEFGHIKLMNPQRSTVWY'
three_letter_codes = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]
amino_acid_names = ["Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Glutamic acid", 
                    "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine",
                    "Phenylalanine", "Proline", "Serine", "Threonine", "Tryptophan", "Tyrosine", "Valine"]

codons_per_aa = {'A': 4, 'W': 1, 'D': 2, 'F': 2, 'S': 6, 'Stop': 3, 'P': 4,
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
             


