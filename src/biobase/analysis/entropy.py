import math

# For the validation of the sequence
from biobase.analysis.nucleic_analysis import Dna


def main() -> None :
    print(entropy("AAAAAAA")) # Homopolymer, exterme low entropy case, exp value = 0.0
    print(entropy("ACGTACGT")) # Equal distribution, exp value = 2.0
    print(entropy("AAACCCGG")) # Mixed, exp value = 1.561278124459133


def entropy(dna_sequence: str) -> float:
    dna_sequence = Dna._validate_dna_sequence(dna_sequence)
    # Calculate the proportion of bases in the sequence
    counts = [dna_sequence.count(nuc) for nuc in Dna.VALID_DNA]
    seq_length = len(dna_sequence)
    entropy = 0
    for count in counts : 
        if count > 0 : # Avoid log(0) because it is undefined
            p = count / seq_length
            entropy -= p * math.log2(p)
    return entropy

if __name__ == "__main__" :
    main()
