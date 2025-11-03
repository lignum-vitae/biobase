from biobase.analysis import nucleic_analysis, motif
from biobase.analysis import cai as cai_module  # avoids naming conflict with function

# Functions
find_motifs = motif.find_motifs
cai = cai_module.cai
ref_counts_from_sequences = cai_module.ref_counts_from_sequences
ref_freqs_from_sequences = cai_module.ref_freqs_from_sequences

# Nucleic classes
Dna = nucleic_analysis.Dna
Nucleotides = nucleic_analysis.Nucleotides
