"""A collection of DNA and RNA constants.


Public module variables:


dna -- a string containing DNA nucleotides
rna -- a string containing RNA nucleotides
complements -- a dictionary containing complements for each DNA nucleotide
"""

dna = 'ATCG'
rna = 'AUCG'
complements = {"A":"T", "T":"A", "C":"G", "G":"C"}

def transcribe(dna_seq: str) -> str:
    if not all(letter in 'ATCG' for letter in dna_seq):
        raise ValueError("Invalid DNA sequence used as argument in transcribe")
    return dna_seq.replace('T', 'U')


def complement_dna(dna_seq: str, reverse: bool = True) -> str:
    if not all(letter in 'ATCG' for letter in dna_seq):
        raise ValueError("Invalid DNA sequence used as argument in transcribe")
    replacements = {"A":"T", "T":"A", "C":"G", "G":"C"}
    if not reverse:
        return "".join([replacements[x] for x in dna_seq])
    return "".join([replacements[x] for x in dna_seq[::-1]])
