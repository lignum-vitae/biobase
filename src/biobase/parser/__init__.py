from biobase.parser import fasta, fastq, genbank

# Fasta parser functions
FastaParser = fasta.FastaParser
FastaFileParser = fasta.FastaFileParser
fasta_parser = fasta.fasta_parser
fasta_file_parser = fasta.fasta_file_parser

# Fastq parser functions
FastqParser = fastq.FastqParser
FastqFileParser = fastq.FastqFileParser
fastq_parser = fastq.fastq_parser
fastq_file_parser = fastq.fastq_file_parser

# Record Class
FastaRecord = fasta.FastaRecord
FastqRecord = fastq.FastqRecord

# GenBank Parser Classes
GenBankRecord = genbank.GenBankRecord
Feature = genbank.Feature
Locus = genbank.Locus
Definition = genbank.Definition
Accession = genbank.Accession
Version = genbank.Version
Origin = genbank.Origin
Features = genbank.Features # The container class for all features
