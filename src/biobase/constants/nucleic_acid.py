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

# ===========================
# RESTRICTION_SITES
# ===========================
# This dictionary maps enzymes to their corresponding restriction recognition site sequences.
# Enzymes recognize specific short sequences in DNA to make cuts.

RESTRICTION_SITES = {
    "AatII": "GACGTC",
    "Acc65I": "GGTACC",
    "AccI": "GTMKAC",
    "AciI": "CCGC",
    "AclI": "AACGTT",
    "AcuI": "CTGAAGNNNNNNNNNNNNNNNN",
    "AfeI": "AGCGCT",
    "AflII": "CTTAAG",
    "AflIII": "ACRYGT",
    "AgeI": "ACCGGT",
    "AhdI": "GACNNNNNGTC",
    "AleI-v2": "CACNNNNGTG",
    "AluI": "AGCT",
    "AlwI": "GGATCNNNNN",
    "AlwNI": "CAGNNNCTG",
    "ApaI": "GGGCCC",
    "ApaLI": "GTGCAC",
    "ApoI": "RAATTY",
    "AscI": "GGCGCGCC",
    "AseI": "ATTAAT",
    "AsiSI": "GCGATCGC",
    "AvaI": "CYCGRG",
    "BsoBI": "CYCGRG",
    "AvaII": "GGWCC",
    "AvrII": "CCTAGG",
    "BaeGI": "GKGCMC",
    "BaeI": "NNNNNNNNNNNNNNNACNNNNGTAYCNNNNNNNNNNNN",
    "BamHI": "GGATCC",
    "BanI": "GGYRCC",
    "BanII": "GRGCYC",
    "BbsI": "GAAGACNNNNNN",
    "BbvCI": "CCTCAGC",
    "BbvI": "GCAGCNNNNNNNNNNNN",
    "BccI": "CCATCNNNNN",
    "BceAI": "ACGGCNNNNNNNNNNNNNN",
    "BcgI": "NNNNNNNNNNNNCGANNNNNNTGCNNNNNNNNNNNN",
    "BciVI": "GTATCCNNNNNN",
    "BclI": "TGATCA",
    "BfaI": "CTAG",
    "BfuAI": "ACCTGCNNNNNNNN",
    "BspMI": "ACCTGCNNNNNNNN",
    "BglI": "GCCNNNNNGGC",
    "BglII": "AGATCT",
    "BlpI": "GCTNAGC",
    "BmgBI": "CACGTC",
    "BmrI": "ACTGGGNNNNN",
    "BmtI": "GCTAGC",
    "BpmI": "CTGGAGNNNNNNNNNNNNNNNN",
    "Bpu10I": "CCTNAGC",
    "BpuEI": "CTTGAGNNNNNNNNNNNNNNNN",
    "BsaAI": "YACGTR",
    "BsaBI": "GATNNNNATC",
    "BsaHI": "GRCGYC",
    "BsaIv2": "GGTCTCNNNNN",
    "BsaJI": "CCNNGG",
    "BsaWI": "WCCGGW",
    "BsaXI": "NNNNNNNNNNNNACNNNNNCTCCNNNNNNNNNN",
    "BseRI": "GAGGAGNNNNNNNNNN",
    "BseYI": "CCCAGC",
    "BsgI": "GTGCAGNNNNNNNNNNNNNNNN",
    "BsiEI": "CGRYCG",
    "BsiHKAI": "GWGCWC",
    "BsiWI": "CGTACG",
    "BslI": "CCNNNNNNNGG",
    "BsmAI": "GTCTCNNNNN",
    "BcoDI": "GTCTCNNNNN",
    "BsmBI-v2": "CGTCTCNNNNN",
    "BsmFI": "GGGACNNNNNNNNNNNNNN",
    "BsmI": "GAATGCN",
    "Bsp1286I": "GDGCHC",
    "BspCNI": "CTCAGNNNNNNNNN",
    "BspEI": "TCCGGA",
    "BspHI": "TCATGA",
    "BsrBI": "CCGCTC",
    "BsrDI": "GCAATGNN",
    "BsrFI-v2": "RCCGGY",
    "BsrGI": "TGTACA",
    "BsrI": "ACTGGN",
    "BssHII": "GCGCGC",
    "BssSI-v2": "CACGAG",
    "BstAPI": "GCANNNNNTGC",
    "BstBI": "TTCGAA",
    "BstEII": "GGTNACC",
    "BstNI": "CCWGG",
    "BstUI": "CGCG",
    "BstXI": "CCANNNNNNTGG",
    "BstYI": "RGATCY",
    "BstZ17I": "GTATAC",
    "Bsu36I": "CCTNAGG",
    "BtgI": "CCRYGG",
    "BtgZI": "GCGATGNNNNNNNNNNNNNN",
    "BtsCI": "GGATGNN",
    "BtsI-v2": "GCAGTGNN",
    "BtsIMutI": "CAGTGNN",
    "Cac8I": "GCNNGC",
    "ClaI": "ATCGAT",
    "BspDI": "ATCGAT",
    "CspCI": "NNNNNNNNNNNNNCAANNNNNGTGGNNNNNNNNNNNN",
    "CviKI-1": "RGCY",
    "CviQI": "GTAC",
    "DdeI": "CTNAG",
    "DpnI": "GATC",
    "DraI": "TTTAAA",
    "DraIII": "CACNNNGTG",
    "DrdI": "GACNNNNNNGTC",
    "EaeI": "YGGCCR",
    "EagI": "CGGCCG",
    "EarI": "CTCTTCNNNN",
    "EciI": "GGCGGANNNNNNNNNNN",
    "Eco53kI": "GAGCTC",
    "EcoNI": "CCTNNNNNAGG",
    "EcoO109I": "RGGNCCY",
    "EcoP15I": "CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNN",
    "EcoRI": "GAATTC",
    "EcoRV": "GATATC",
    "Esp3I": "CGTCTCNNNNN",
    "FatI": "CATG",
    "FauI": "CCCGCNNNNNN",
    "Fnu4HI": "GCNGC",
    "FokI": "GGATGNNNNNNNNNNNNN",
    "FseI": "GGCCGGCC",
    "FspI": "TGCGCA",
    "HaeII": "RGCGCY",
    "HaeIII": "GGCC",
    "HgaI": "GACGCNNNNNNNNNN",
    "HhaI": "GCGC",
    "HincII": "GTYRAC",
    "HindIII": "AAGCTT",
    "HinfI": "GANTC",
    "HinP1I": "GCGC",
    "HpaI": "GTTAAC",
    "HpaII": "CCGG",
    "MspI": "CCGG",
    "HphI": "GGTGANNNNNNNN",
    "Hpy166II": "GTNNAC",
    "Hpy188I": "TCNGA",
    "Hpy188III": "TCNNGA",
    "Hpy99I": "CGWCG",
    "HpyAV": "CCTTCNNNNNN",
    "HpyCH4III": "ACNGT",
    "HpyCH4IV": "ACGT",
    "HpyCH4V": "TGCA",
    "I-CeuI": "TAACTATAACGGTCCTAAGGTAGCGAA",
    "I-SceI": "TAGGGATAACAGGGTAAT",
    "KasI": "GGCGCC",
    "KpnI": "GGTACC",
    "MboII": "GAAGANNNNNNNN",
    "MfeI": "CAATTG",
    "MluCI": "AATT",
    "MluI": "ACGCGT",
    "MlyI": "GAGTCNNNNN",
    "MmeI": "TCCRACNNNNNNNNNNNNNNNNNNNN",
    "MnlI": "CCTCNNNNNNN",
    "MscI": "TGGCCA",
    "MseI": "TTAA",
    "MslI": "CAYNNNNRTG",
    "MspA1I": "CMGCKG",
    "MspJI": "CNNRNNNNNNNNNNNNN",
    "MwoI": "GCNNNNNNNGC",
    "NaeI": "GCCGGC",
    "NarI": "GGCGCC",
    "Nb.BbvCI": "CCTCAGC",
    "Nb.BsmI": "GAATGC",
    "Nb.BsrDI": "GCAATG",
    "Nb.BssSI": "CACGAG",
    "Nb.BtsI": "GCAGTG",
    "NciI": "CCSGG",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "NgoMIV": "GCCGGC",
    "NheI": "GCTAGC",
    "NlaIII": "CATG",
    "NlaIV": "GGNNCC",
    "NmeAIII": "GCCGAGNNNNNNNNNNNNNNNNNNNNN",
    "NotI": "GCGGCCGC",
    "NruI": "TCGCGA",
    "NsiI": "ATGCAT",
    "NspI": "RCATGY",
    "Nt.AlwI": "GGATCNNNN",
    "Nt.BbvCI": "CCTCAGC",
    "Nt.BsmAI": "GTCTCN",
    "Nt.BspQI": "GCTCTTCN",
    "Nt.BstNBI": "GAGTCNNNN",
    "Nt.CviPII": "CCD",
    "PacI": "TTAATTAA",
    "PaqCI": "CACCTGCNNNNNNNN",
    "PciI": "ACATGT",
    "PflFI": "GACNNNGTC",
    "Tth111I": "GACNNNGTC",
    "PflMI": "CCANNNNNTGG",
    "PI-PspI": "TGGCAAACAGCTATTATGGGTATTATGGGT",
    "PI-SceI": "ATCTATGTCGGGTGCGGAGAAAGAGGTAAT",
    "PleI": "GAGTCNNNNN",
    "PluTI": "GGCGCC",
    "PmeI": "GTTTAAAC",
    "PmlI": "CACGTG",
    "PpuMI": "RGGWCCY",
    "PshAI": "GACNNNNGTC",
    "PsiI-v2": "TTATAA",
    "PspGI": "CCWGG",
    "PspOMI": "GGGCCC",
    "PspXI": "VCTCGAGB",
    "PstI": "CTGCAG",
    "PvuI": "CGATCG",
    "PvuII": "CAGCTG",
    "RsaI": "GTAC",
    "RsrII": "CGGWCCG",
    "SacI": "GAGCTC",
    "SacII": "CCGCGG",
    "SalI": "GTCGAC",
    "SapI": "GCTCTTCNNNN",
    "BspQI": "GCTCTTCNNNN",
    "Sau3AI": "GATC",
    "MboI": "GATC",
    "DpnII": "GATC",
    "Sau96I": "GGNCC",
    "SbfI": "CCTGCAGG",
    "ScaI": "AGTACT",
    "ScrFI": "CCNGG",
    "SexAI": "ACCWGGT",
    "SfaNI": "GCATCNNNNNNNNN",
    "SfcI": "CTRYAG",
    "SfiI": "GGCCNNNNNGGCC",
    "SfoI": "GGCGCC",
    "SgrAI": "CRCCGGYG",
    "SmaI": "CCCGGG",
    "SmlI": "CTYRAG",
    "SnaBI": "TACGTA",
    "SpeI": "ACTAGT",
    "SphI": "GCATGC",
    "SrfI": "GCCCGGGC",
    "SspI": "AATATT",
    "StuI": "AGGCCT",
    "StyD4I": "CCNGG",
    "StyI": "CCWWGG",
    "SwaI": "ATTTAAAT",
    "TaqI-v2": "TCGA",
    "TfiI": "GAWTC",
    "TseI": "GCWGC",
    "ApeKI": "GCWGC",
    "Tsp45I": "GTSAC",
    "TspRI": "NNCASTGNN",
    "XbaI": "TCTAGA",
    "XcmI": "CCANNNNNNNNNTGG",
    "XhoI": "CTCGAG",
    "PaeR7I": "CTCGAG",
    "XmaI": "CCCGGG",
    "TspMI": "CCCGGG",
    "XmnI": "GAANNNNTTC",
    "ZraI": "GACGTC"
}

# ===========================
# RESTRICTION_SITES_REV
# ===========================
# This dictionary contains the reverse complement of each restriction site in RESTRICTION_SITES.
# Reverse complements are important because many enzymes recognize palindromic sequences, 
# where the sequence reads the same forward and backward on the opposite strand.

RESTRICTION_SITES_REV = {
    "AatII": "GACGTC",
    "Acc65I": "GGTACC",
    "AccI": "GTMKAC",
    "AciI": "GCGG",
    "AclI": "AACGTT",
    "AcuI": "NNNNNNNNNNNNNNNNCTTCAG",
    "AfeI": "AGCGCT",
    "AflII": "CTTAAG",
    "AflIII": "ACRYGT",
    "AgeI": "ACCGGT",
    "AhdI": "GACNNNNNGTC",
    "AleI-v2": "CACNNNNGTG",
    "AluI": "AGCT",
    "AlwI": "NNNNNGATCC",
    "AlwNI": "CAGNNNCTG",
    "ApaI": "GGGCCC",
    "ApaLI": "GTGCAC",
    "ApoI": "RAATTY",
    "AscI": "GGCGCGCC",
    "AseI": "ATTAAT",
    "AsiSI": "GCGATCGC",
    "AvaI": "CYCGRG",
    "BsoBI": "CYCGRG",
    "AvaII": "GGWCC",
    "AvrII": "CCTAGG",
    "BaeGI": "GKGCMC",
    "BaeI": "NNNNNNNNNNNNGRTACNNNNGTNNNNNNNNNNNNNNN",
    "BamHI": "GGATCC",
    "BanI": "GGYRCC",
    "BanII": "GRGCYC",
    "BbsI": "NNNNNNGTCTTC",
    "BbvCI": "GCTGAGG",
    "BbvI": "NNNNNNNNNNNNGCTGC",
    "BccI": "NNNNNGATGG",
    "BceAI": "NNNNNNNNNNNNNNGCCGT",
    "BcgI": "NNNNNNNNNNNNGCANNNNNNTCGNNNNNNNNNNNN",
    "BciVI": "NNNNNNGGATAC",
    "BclI": "TGATCA",
    "BfaI": "CTAG",
    "BfuAI": "NNNNNNNNGCAGGT",
    "BspMI": "NNNNNNNNGCAGGT",
    "BglI": "GCCNNNNNGGC",
    "BglII": "AGATCT",
    "BlpI": "GCTNAGC",
    "BmgBI": "GACGTG",
    "BmrI": "NNNNNCCCAGT",
    "BmtI": "GCTAGC",
    "BpmI": "NNNNNNNNNNNNNNNNCTCCAG",
    "Bpu10I": "GCTNAGG",
    "BpuEI": "NNNNNNNNNNNNNNNNCTCAAG",
    "BsaAI": "YACGTR",
    "BsaBI": "GATNNNNATC",
    "BsaHI": "GRCGYC",
    "BsaIv2": "NNNNNGAGACC",
    "BsaJI": "CCNNGG",
    "BsaWI": "WCCGGW",
    "BsaXI": "NNNNNNNNNNGGAGNNNNNGTNNNNNNNNNNNN",
    "BseRI": "NNNNNNNNNNCTCCTC",
    "BseYI": "GCTGGG",
    "BsgI": "NNNNNNNNNNNNNNNNCTGCAC",
    "BsiEI": "CGRYCG",
    "BsiHKAI": "GWGCWC",
    "BsiWI": "CGTACG",
    "BslI": "CCNNNNNNNGG",
    "BsmAI": "NNNNNGAGAC",
    "BcoDI": "NNNNNGAGAC",
    "BsmBI-v2": "NNNNNGAGACG",
    "BsmFI": "NNNNNNNNNNNNNNGTCCC",
    "BsmI": "NGCATTC",
    "Bsp1286I": "GDGCHC",
    "BspCNI": "NNNNNNNNNCTGAG",
    "BspEI": "TCCGGA",
    "BspHI": "TCATGA",
    "BsrBI": "GAGCGG",
    "BsrDI": "NNCATTGC",
    "BsrFI-v2": "RCCGGY",
    "BsrGI": "TGTACA",
    "BsrI": "NCCAGT",
    "BssHII": "GCGCGC",
    "BssSI-v2": "CTCGTG",
    "BstAPI": "GCANNNNNTGC",
    "BstBI": "TTCGAA",
    "BstEII": "GGTNACC",
    "BstNI": "CCWGG",
    "BstUI": "CGCG",
    "BstXI": "CCANNNNNNTGG",
    "BstYI": "RGATCY",
    "BstZ17I": "GTATAC",
    "Bsu36I": "CCTNAGG",
    "BtgI": "CCRYGG",
    "BtgZI": "NNNNNNNNNNNNNNCATCGC",
    "BtsCI": "NNCATCC",
    "BtsI-v2": "NNCACTGC",
    "BtsIMutI": "NNCACTG",
    "Cac8I": "GCNNGC",
    "ClaI": "ATCGAT",
    "BspDI": "ATCGAT",
    "CspCI": "NNNNNNNNNNNNCCACNNNNNTTGNNNNNNNNNNNNN",
    "CviKI-1": "RGCY",
    "CviQI": "GTAC",
    "DdeI": "CTNAG",
    "DpnI": "GATC",
    "DraI": "TTTAAA",
    "DraIII": "CACNNNGTG",
    "DrdI": "GACNNNNNNGTC",
    "EaeI": "YGGCCR",
    "EagI": "CGGCCG",
    "EarI": "NNNNGAAGAG",
    "EciI": "NNNNNNNNNNNTCCGCC",
    "Eco53kI": "GAGCTC",
    "EcoNI": "CCTNNNNNAGG",
    "EcoO109I": "RGGNCCY",
    "EcoP15I": "NNNNNNNNNNNNNNNNNNNNNNNNNNNCTGCTG",
    "EcoRI": "GAATTC",
    "EcoRV": "GATATC",
    "Esp3I": "NNNNNGAGACG",
    "FatI": "CATG",
    "FauI": "NNNNNNGCGGG",
    "Fnu4HI": "GCNGC",
    "FokI": "NNNNNNNNNNNNNCATCC",
    "FseI": "GGCCGGCC",
    "FspI": "TGCGCA",
    "HaeII": "RGCGCY",
    "HaeIII": "GGCC",
    "HgaI": "NNNNNNNNNNGCGTC",
    "HhaI": "GCGC",
    "HincII": "GTYRAC",
    "HindIII": "AAGCTT",
    "HinfI": "GANTC",
    "HinP1I": "GCGC",
    "HpaI": "GTTAAC",
    "HpaII": "CCGG",
    "MspI": "CCGG",
    "HphI": "NNNNNNNNTCACC",
    "Hpy166II": "GTNNAC",
    "Hpy188I": "TCNGA",
    "Hpy188III": "TCNNGA",
    "Hpy99I": "CGWCG",
    "HpyAV": "NNNNNNGAAGG",
    "HpyCH4III": "ACNGT",
    "HpyCH4IV": "ACGT",
    "HpyCH4V": "TGCA",
    "I-CeuI": "TTCGCTACCTTAGGACCGTTATAGTTA",
    "I-SceI": "ATTACCCTGTTATCCCTA",
    "KasI": "GGCGCC",
    "KpnI": "GGTACC",
    "MboII": "NNNNNNNNTCTTC",
    "MfeI": "CAATTG",
    "MluCI": "AATT",
    "MluI": "ACGCGT",
    "MlyI": "NNNNNGACTC",
    "MmeI": "NNNNNNNNNNNNNNNNNNNNGTYGGA",
    "MnlI": "NNNNNNNGAGG",
    "MscI": "TGGCCA",
    "MseI": "TTAA",
    "MslI": "CAYNNNNRTG",
    "MspA1I": "CMGCKG",
    "MspJI": "NNNNNNNNNNNNNYNNG",
    "MwoI": "GCNNNNNNNGC",
    "NaeI": "GCCGGC",
    "NarI": "GGCGCC",
    "Nb.BbvCI": "GCTGAGG",
    "Nb.BsmI": "GCATTC",
    "Nb.BsrDI": "CATTGC",
    "Nb.BssSI": "CTCGTG",
    "Nb.BtsI": "CACTGC",
    "NciI": "CCSGG",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "NgoMIV": "GCCGGC",
    "NheI": "GCTAGC",
    "NlaIII": "CATG",
    "NlaIV": "GGNNCC",
    "NmeAIII": "NNNNNNNNNNNNNNNNNNNNNCTCGGC",
    "NotI": "GCGGCCGC",
    "NruI": "TCGCGA",
    "NsiI": "ATGCAT",
    "NspI": "RCATGY",
    "Nt.AlwI": "NNNNGATCC",
    "Nt.BbvCI": "GCTGAGG",
    "Nt.BsmAI": "NGAGAC",
    "Nt.BspQI": "NGAAGAGC",
    "Nt.BstNBI": "NNNNGACTC",
    "Nt.CviPII": "HGG",
    "PacI": "TTAATTAA",
    "PaqCI": "NNNNNNNNGCAGGTG",
    "PciI": "ACATGT",
    "PflFI": "GACNNNGTC",
    "Tth111I": "GACNNNGTC",
    "PflMI": "CCANNNNNTGG",
    "PI-PspI": "ACCCATAATACCCATAATAGCTGTTTGCCA",
    "PI-SceI": "ATTACCTCTTTCTCCGCACCCGACATAGAT",
    "PleI": "NNNNNGACTC",
    "PluTI": "GGCGCC",
    "PmeI": "GTTTAAAC",
    "PmlI": "CACGTG",
    "PpuMI": "RGGWCCY",
    "PshAI": "GACNNNNGTC",
    "PsiI-v2": "TTATAA",
    "PspGI": "CCWGG",
    "PspOMI": "GGGCCC",
    "PspXI": "VCTCGAGB",
    "PstI": "CTGCAG",
    "PvuI": "CGATCG",
    "PvuII": "CAGCTG",
    "RsaI": "GTAC",
    "RsrII": "CGGWCCG",
    "SacI": "GAGCTC",
    "SacII": "CCGCGG",
    "SalI": "GTCGAC",
    "SapI": "NNNNGAAGAGC",
    "BspQI": "NNNNGAAGAGC",
    "Sau3AI": "GATC",
    "MboI": "GATC",
    "DpnII": "GATC",
    "Sau96I": "GGNCC",
    "SbfI": "CCTGCAGG",
    "ScaI": "AGTACT",
    "ScrFI": "CCNGG",
    "SexAI": "ACCWGGT",
    "SfaNI": "NNNNNNNNNGATGC",
    "SfcI": "CTRYAG",
    "SfiI": "GGCCNNNNNGGCC",
    "SfoI": "GGCGCC",
    "SgrAI": "CRCCGGYG",
    "SmaI": "CCCGGG",
    "SmlI": "CTYRAG",
    "SnaBI": "TACGTA",
    "SpeI": "ACTAGT",
    "SphI": "GCATGC",
    "SrfI": "GCCCGGGC",
    "SspI": "AATATT",
    "StuI": "AGGCCT",
    "StyD4I": "CCNGG",
    "StyI": "CCWWGG",
    "SwaI": "ATTTAAAT",
    "TaqI-v2": "TCGA",
    "TfiI": "GAWTC",
    "TseI": "GCWGC",
    "ApeKI": "GCWGC",
    "Tsp45I": "GTSAC",
    "TspRI": "NNCASTGNN",
    "XbaI": "TCTAGA",
    "XcmI": "CCANNNNNNNNNTGG",
    "XhoI": "CTCGAG",
    "PaeR7I": "CTCGAG",
    "XmaI": "CCCGGG",
    "TspMI": "CCCGGG",
    "XmnI": "GAANNNNTTC",
    "ZraI": "GACGTC"
}

# ===========================
# RESTRICTION_SITES_CUTTING
# ===========================
# This dictionary maps restriction enzyme names to their cutting patterns.
# Cutting information is represented by specific symbols:
# - "|" cut to both strands (blunt cut)
# - "/" for cuts on the given strand (sticky ends)
# - "^" for cuts on the complementary strand (sticky ends)
# - Some enzymes may produce a single "nick" (either "/" or "^").
# - Rare enzymes may cut at two locations on each strand, resulting in 4 total cutting symbols.

RESTRICTION_SITES_CUTTING = {
    "AatII": "G^ACGT/C",
    "Acc65I": "G/GTAC^C",
    "AccI": "GT/MK^AC",
    "AciI": "C/CG^C",
    "AclI": "AA/CG^TT",
    "AcuI": "CTGAAGNNNNNNNNNNNNNN^NN/",
    "AfeI": "AGC|GCT",
    "AflII": "C/TTAA^G",
    "AflIII": "A/CRYG^T",
    "AgeI": "A/CCGG^T",
    "AhdI": "GACNN^N/NNGTC",
    "AleI-v2": "CACNN|NNGTG",
    "AluI": "AG|CT",
    "AlwI": "GGATCNNNN/N^",
    "AlwNI": "CAG^NNN/CTG",
    "ApaI": "G^GGCC/C",
    "ApaLI": "G/TGCA^C",
    "ApoI": "R/AATT^Y",
    "AscI": "GG/CGCG^CC",
    "AseI": "AT/TA^AT",
    "AsiSI": "GCG^AT/CGC",
    "AvaI": "C/YCGR^G",
    "BsoBI": "C/YCGR^G",
    "AvaII": "G/GWC^C",
    "AvrII": "C/CTAG^G",
    "BaeGI": "G^KGCM/C",
    "BaeI": "^NNNNN/NNNNNNNNNNACNNNNGTAYCNNNNNNN^NNNNN/",
    "BamHI": "G/GATC^C",
    "BanI": "G/GYRC^C",
    "BanII": "G^RGCY/C",
    "BbsI": "GAAGACNN/NNNN^",
    "BbvCI": "CC/TCA^GC",
    "BbvI": "GCAGCNNNNNNNN/NNNN^",
    "BccI": "CCATCNNNN/N^",
    "BceAI": "ACGGCNNNNNNNNNNNN/NN^",
    "BcgI": "^NN/NNNNNNNNNNCGANNNNNNTGCNNNNNNNNNN^NN/",
    "BciVI": "GTATCCNNNNN^N/",
    "BclI": "T/GATC^A",
    "BfaI": "C/TA^G",
    "BfuAI": "ACCTGCNNNN/NNNN^",
    "BspMI": "ACCTGCNNNN/NNNN^",
    "BglI": "GCCN^NNN/NGGC",
    "BglII": "A/GATC^T",
    "BlpI": "GC/TNA^GC",
    "BmgBI": "CAC|GTC",
    "BmrI": "ACTGGGNNNN^N/",
    "BmtI": "G^CTAG/C",
    "BpmI": "CTGGAGNNNNNNNNNNNNNN^NN/",
    "Bpu10I": "CC/TNA^GC",
    "BpuEI": "CTTGAGNNNNNNNNNNNNNN^NN/",
    "BsaAI": "YAC|GTR",
    "BsaBI": "GATNN|NNATC",
    "BsaHI": "GR/CG^YC",
    "BsaIv2": "GGTCTCN/NNNN^",
    "BsaJI": "C/CNNG^G",
    "BsaWI": "W/CCGG^W",
    "BsaXI": "^NNN/NNNNNNNNNACNNNNNCTCCNNNNNNN^NNN/",
    "BseRI": "GAGGAGNNNNNNNN^NN/",
    "BseYI": "C/CCAG^C",
    "BsgI": "GTGCAGNNNNNNNNNNNNNN^NN/",
    "BsiEI": "CG^RY/CG",
    "BsiHKAI": "G^WGCW/C",
    "BsiWI": "C/GTAC^G",
    "BslI": "CCNN^NNN/NNGG",
    "BsmAI": "GTCTCN/NNNN^",
    "BcoDI": "GTCTCN/NNNN^",
    "BsmBI-v2": "CGTCTCN/NNNN^",
    "BsmFI": "GGGACNNNNNNNNNN/NNNN^",
    "BsmI": "GAATG^CN/",
    "Bsp1286I": "G^DGCH/C",
    "BspCNI": "CTCAGNNNNNNN^NN/",
    "BspEI": "T/CCGG^A",
    "BspHI": "T/CATG^A",
    "BsrBI": "CCG|CTC",
    "BsrDI": "GCAATG^NN/",
    "BsrFI-v2": "R/CCGG^Y",
    "BsrGI": "T/GTAC^A",
    "BsrI": "ACTG^GN/",
    "BssHII": "G/CGCG^C",
    "BssSI-v2": "C/ACGA^G",
    "BstAPI": "GCAN^NNN/NTGC",
    "BstBI": "TT/CG^AA",
    "BstEII": "G/GTNAC^C",
    "BstNI": "CC/W^GG",
    "BstUI": "CG|CG",
    "BstXI": "CCAN^NNNN/NTGG",
    "BstYI": "R/GATC^Y",
    "BstZ17I": "GTA|TAC",
    "Bsu36I": "CC/TNA^GG",
    "BtgI": "C/CRYG^G",
    "BtgZI": "GCGATGNNNNNNNNNN/NNNN^",
    "BtsCI": "GGATG^NN/",
    "BtsI-v2": "GCAGTG^NN/",
    "BtsIMutI": "CAGTG^NN/",
    "Cac8I": "GCN|NGC",
    "ClaI": "AT/CG^AT",
    "BspDI": "AT/CG^AT",
    "CspCI": "^NN/NNNNNNNNNNNCAANNNNNGTGGNNNNNNNNNN^NN/",
    "CviKI-1": "RG|CY",
    "CviQI": "G/TA^C",
    "DdeI": "C/TNA^G",
    "DpnI": "GA|TC",
    "DraI": "TTT|AAA",
    "DraIII": "CAC^NNN/GTG",
    "DrdI": "GACNN^NN/NNGTC",
    "EaeI": "Y/GGCC^R",
    "EagI": "C/GGCC^G",
    "EarI": "CTCTTCN/NNN^",
    "EciI": "GGCGGANNNNNNNNN^NN/",
    "Eco53kI": "GAG|CTC",
    "EcoNI": "CCTNN/N^NNAGG",
    "EcoO109I": "RG/GNC^CY",
    "EcoP15I": "CAGCAGNNNNNNNNNNNNNNNNNNNNNNNNN/NN^",
    "EcoRI": "G/AATT^C",
    "EcoRV": "GAT|ATC",
    "Esp3I": "CGTCTCN/NNNN^",
    "FatI": "/CATG^",
    "FauI": "CCCGCNNNN/NN^",
    "Fnu4HI": "GC/N^GC",
    "FokI": "GGATGNNNNNNNNN/NNNN^",
    "FseI": "GG^CCGG/CC",
    "FspI": "TGC|GCA",
    "HaeII": "R^GCGC/Y",
    "HaeIII": "GG|CC",
    "HgaI": "GACGCNNNNN/NNNNN^",
    "HhaI": "G^CG/C",
    "HincII": "GTY|RAC",
    "HindIII": "A/AGCT^T",
    "HinfI": "G/ANT^C",
    "HinP1I": "G/CG^C",
    "HpaI": "GTT|AAC",
    "HpaII": "C/CG^G",
    "MspI": "C/CG^G",
    "HphI": "GGTGANNNNNNN^N/",
    "Hpy166II": "GTN|NAC",
    "Hpy188I": "TC^N/GA",
    "Hpy188III": "TC/NN^GA",
    "Hpy99I": "^CGWCG/",
    "HpyAV": "CCTTCNNNNN^N/",
    "HpyCH4III": "AC^N/GT",
    "HpyCH4IV": "A/CG^T",
    "HpyCH4V": "TG|CA",
    "I-CeuI": "TAACTATAACGGTC^CTAA/GGTAGCGAA",
    "I-SceI": "TAGGG^ATAA/CAGGGTAAT",
    "KasI": "G/GCGC^C",
    "KpnI": "G^GTAC/C",
    "MboII": "GAAGANNNNNNN^N/",
    "MfeI": "C/AATT^G",
    "MluCI": "/AATT^",
    "MluI": "A/CGCG^T",
    "MlyI": "GAGTCNNNNN|",
    "MmeI": "TCCRACNNNNNNNNNNNNNNNNNN^NN/",
    "MnlI": "CCTCNNNNNN^N/",
    "MscI": "TGG|CCA",
    "MseI": "T/TA^A",
    "MslI": "CAYNN|NNRTG",
    "MspA1I": "CMG|CKG",
    "MspJI": "CNNRNNNNNNNNN/NNNN^",
    "MwoI": "GCNN^NNN/NNGC",
    "NaeI": "GCC|GGC",
    "NarI": "GG/CG^CC",
    "Nb.BbvCI": "CCTCA^GC",
    "Nb.BsmI": "GAATG^C",
    "Nb.BsrDI": "GCAATG^",
    "Nb.BssSI": "CACGA^G",
    "Nb.BtsI": "GCAGTG^",
    "NciI": "CC/S^GG",
    "NcoI": "C/CATG^G",
    "NdeI": "CA/TA^TG",
    "NgoMIV": "G/CCGG^C",
    "NheI": "G/CTAG^C",
    "NlaIII": "^CATG/",
    "NlaIV": "GGN|NCC",
    "NmeAIII": "GCCGAGNNNNNNNNNNNNNNNNNNN^NN/",
    "NotI": "GC/GGCC^GC",
    "NruI": "TCG|CGA",
    "NsiI": "A^TGCA/T",
    "NspI": "R^CATG/Y",
    "Nt.AlwI": "GGATCNNNN/",
    "Nt.BbvCI": "CC/TCAGC",
    "Nt.BsmAI": "GTCTCN/",
    "Nt.BspQI": "GCTCTTCN/",
    "Nt.BstNBI": "GAGTCNNNN/",
    "Nt.CviPII": "/CCD",
    "PacI": "TTA^AT/TAA",
    "PaqCI": "CACCTGCNNNN/NNNN^",
    "PciI": "A/CATG^T",
    "PflFI": "GACN/N^NGTC",
    "Tth111I": "GACN/N^NGTC",
    "PflMI": "CCAN^NNN/NTGG",
    "PI-PspI": "TGGCAAACAGCTA^TTAT/GGGTATTATGGGT",
    "PI-SceI": "ATCTATGTCGG^GTGC/GGAGAAAGAGGTAAT",
    "PleI": "GAGTCNNNN/N^",
    "PluTI": "G^GCGC/C",
    "PmeI": "GTTT|AAAC",
    "PmlI": "CAC|GTG",
    "PpuMI": "RG/GWC^CY",
    "PshAI": "GACNN|NNGTC",
    "PsiI-v2": "TTA|TAA",
    "PspGI": "/CCWGG^",
    "PspOMI": "G/GGCC^C",
    "PspXI": "VC/TCGA^GB",
    "PstI": "C^TGCA/G",
    "PvuI": "CG^AT/CG",
    "PvuII": "CAG|CTG",
    "RsaI": "GT|AC",
    "RsrII": "CG/GWC^CG",
    "SacI": "G^AGCT/C",
    "SacII": "CC^GC/GG",
    "SalI": "G/TCGA^C",
    "SapI": "GCTCTTCN/NNN^",
    "BspQI": "GCTCTTCN/NNN^",
    "Sau3AI": "/GATC^",
    "MboI": "/GATC^",
    "DpnII": "/GATC^",
    "Sau96I": "G/GNC^C",
    "SbfI": "CC^TGCA/GG",
    "ScaI": "AGT|ACT",
    "ScrFI": "CC/N^GG",
    "SexAI": "A/CCWGG^T",
    "SfaNI": "GCATCNNNNN/NNNN^",
    "SfcI": "C/TRYA^G",
    "SfiI": "GGCCN^NNN/NGGCC",
    "SfoI": "GGC|GCC",
    "SgrAI": "CR/CCGG^YG",
    "SmaI": "CCC|GGG",
    "SmlI": "C/TYRA^G",
    "SnaBI": "TAC|GTA",
    "SpeI": "A/CTAG^T",
    "SphI": "G^CATG/C",
    "SrfI": "GCCC|GGGC",
    "SspI": "AAT|ATT",
    "StuI": "AGG|CCT",
    "StyD4I": "/CCNGG^",
    "StyI": "C/CWWG^G",
    "SwaI": "ATTT|AAAT",
    "TaqI-v2": "T/CG^A",
    "TfiI": "G/AWT^C",
    "TseI": "G/CWG^C",
    "ApeKI": "G/CWG^C",
    "Tsp45I": "/GTSAC^",
    "TspRI": "^NNCASTGNN/",
    "XbaI": "T/CTAG^A",
    "XcmI": "CCANNNN^N/NNNNTGG",
    "XhoI": "C/TCGA^G",
    "PaeR7I": "C/TCGA^G",
    "XmaI": "C/CCGG^G",
    "TspMI": "C/CCGG^G",
    "XmnI": "GAANN|NNTTC",
    "ZraI": "GAC|GTC"
}

# ===========================
# RESTRICTION_SITES_CUTTING_REV
# ===========================
# This dictionary contains reverse complement cutting information.
# The cutting symbols ("/" and "^") are adjusted to reflect the reverse complement strand.

RESTRICTION_SITES_CUTTING_REV = {
    "AatII": "G^ACGT/C",
    "Acc65I": "G/GTAC^C",
    "AccI": "GT/MK^AC",
    "AciI": "G/CG^G",
    "AclI": "AA/CG^TT",
    "AcuI": "^NN/NNNNNNNNNNNNNNCTTCAG",
    "AfeI": "AGC|GCT",
    "AflII": "C/TTAA^G",
    "AflIII": "A/CRYG^T",
    "AgeI": "A/CCGG^T",
    "AhdI": "GACNN^N/NNGTC",
    "AleI-v2": "CACNN|NNGTG",
    "AluI": "AG|CT",
    "AlwI": "/N^NNNNGATCC",
    "AlwNI": "CAG^NNN/CTG",
    "ApaI": "G^GGCC/C",
    "ApaLI": "G/TGCA^C",
    "ApoI": "R/AATT^Y",
    "AscI": "GG/CGCG^CC",
    "AseI": "AT/TA^AT",
    "AsiSI": "GCG^AT/CGC",
    "AvaI": "C/YCGR^G",
    "BsoBI": "C/YCGR^G",
    "AvaII": "G/GWC^C",
    "AvrII": "C/CTAG^G",
    "BaeGI": "G^KGCM/C",
    "BaeI": "^NNNNN/NNNNNNNGRTACNNNNGTNNNNNNNNNN^NNNNN/",
    "BamHI": "G/GATC^C",
    "BanI": "G/GYRC^C",
    "BanII": "G^RGCY/C",
    "BbsI": "/NNNN^NNGTCTTC",
    "BbvCI": "GC/TGA^GG",
    "BbvI": "/NNNN^NNNNNNNNGCTGC",
    "BccI": "/N^NNNNGATGG",
    "BceAI": "/NN^NNNNNNNNNNNNGCCGT",
    "BcgI": "^NN/NNNNNNNNNNGCANNNNNNTCGNNNNNNNNNN^NN/",
    "BciVI": "^N/NNNNNGGATAC",
    "BclI": "T/GATC^A",
    "BfaI": "C/TA^G",
    "BfuAI": "/NNNN^NNNNGCAGGT",
    "BspMI": "/NNNN^NNNNGCAGGT",
    "BglI": "GCCN^NNN/NGGC",
    "BglII": "A/GATC^T",
    "BlpI": "GC/TNA^GC",
    "BmgBI": "GAC|GTG",
    "BmrI": "^N/NNNNCCCAGT",
    "BmtI": "G^CTAG/C",
    "BpmI": "^NN/NNNNNNNNNNNNNNCTCCAG",
    "Bpu10I": "GC/TNA^GG",
    "BpuEI": "^NN/NNNNNNNNNNNNNNCTCAAG",
    "BsaAI": "YAC|GTR",
    "BsaBI": "GATNN|NNATC",
    "BsaHI": "GR/CG^YC",
    "BsaIv2": "/NNNN^NGAGACC",
    "BsaJI": "C/CNNG^G",
    "BsaWI": "W/CCGG^W",
    "BsaXI": "^NNN/NNNNNNNGGAGNNNNNGTNNNNNNNNN^NNN/",
    "BseRI": "^NN/NNNNNNNNCTCCTC",
    "BseYI": "G/CTGG^G",
    "BsgI": "^NN/NNNNNNNNNNNNNNCTGCAC",
    "BsiEI": "CG^RY/CG",
    "BsiHKAI": "G^WGCW/C",
    "BsiWI": "C/GTAC^G",
    "BslI": "CCNN^NNN/NNGG",
    "BsmAI": "/NNNN^NGAGAC",
    "BcoDI": "/NNNN^NGAGAC",
    "BsmBI-v2": "/NNNN^NGAGACG",
    "BsmFI": "/NNNN^NNNNNNNNNNGTCCC",
    "BsmI": "^NG/CATTC",
    "Bsp1286I": "G^DGCH/C",
    "BspCNI": "^NN/NNNNNNNCTGAG",
    "BspEI": "T/CCGG^A",
    "BspHI": "T/CATG^A",
    "BsrBI": "GAG|CGG",
    "BsrDI": "^NN/CATTGC",
    "BsrFI-v2": "R/CCGG^Y",
    "BsrGI": "T/GTAC^A",
    "BsrI": "^NC/CAGT",
    "BssHII": "G/CGCG^C",
    "BssSI-v2": "C/TCGT^G",
    "BstAPI": "GCAN^NNN/NTGC",
    "BstBI": "TT/CG^AA",
    "BstEII": "G/GTNAC^C",
    "BstNI": "CC/W^GG",
    "BstUI": "CG|CG",
    "BstXI": "CCAN^NNNN/NTGG",
    "BstYI": "R/GATC^Y",
    "BstZ17I": "GTA|TAC",
    "Bsu36I": "CC/TNA^GG",
    "BtgI": "C/CRYG^G",
    "BtgZI": "/NNNN^NNNNNNNNNNCATCGC",
    "BtsCI": "^NN/CATCC",
    "BtsI-v2": "^NN/CACTGC",
    "BtsIMutI": "^NN/CACTG",
    "Cac8I": "GCN|NGC",
    "ClaI": "AT/CG^AT",
    "BspDI": "AT/CG^AT",
    "CspCI": "^NN/NNNNNNNNNNCCACNNNNNTTGNNNNNNNNNNN^NN/",
    "CviKI-1": "RG|CY",
    "CviQI": "G/TA^C",
    "DdeI": "C/TNA^G",
    "DpnI": "GA|TC",
    "DraI": "TTT|AAA",
    "DraIII": "CAC^NNN/GTG",
    "DrdI": "GACNN^NN/NNGTC",
    "EaeI": "Y/GGCC^R",
    "EagI": "C/GGCC^G",
    "EarI": "/NNN^NGAAGAG",
    "EciI": "^NN/NNNNNNNNNTCCGCC",
    "Eco53kI": "GAG|CTC",
    "EcoNI": "CCTNN/N^NNAGG",
    "EcoO109I": "RG/GNC^CY",
    "EcoP15I": "/NN^NNNNNNNNNNNNNNNNNNNNNNNNNCTGCTG",
    "EcoRI": "G/AATT^C",
    "EcoRV": "GAT|ATC",
    "Esp3I": "/NNNN^NGAGACG",
    "FatI": "/CATG^",
    "FauI": "/NN^NNNNGCGGG",
    "Fnu4HI": "GC/N^GC",
    "FokI": "/NNNN^NNNNNNNNNCATCC",
    "FseI": "GG^CCGG/CC",
    "FspI": "TGC|GCA",
    "HaeII": "R^GCGC/Y",
    "HaeIII": "GG|CC",
    "HgaI": "/NNNNN^NNNNNGCGTC",
    "HhaI": "G^CG/C",
    "HincII": "GTY|RAC",
    "HindIII": "A/AGCT^T",
    "HinfI": "G/ANT^C",
    "HinP1I": "G/CG^C",
    "HpaI": "GTT|AAC",
    "HpaII": "C/CG^G",
    "MspI": "C/CG^G",
    "HphI": "^N/NNNNNNNTCACC",
    "Hpy166II": "GTN|NAC",
    "Hpy188I": "TC^N/GA",
    "Hpy188III": "TC/NN^GA",
    "Hpy99I": "^CGWCG/",
    "HpyAV": "^N/NNNNNGAAGG",
    "HpyCH4III": "AC^N/GT",
    "HpyCH4IV": "A/CG^T",
    "HpyCH4V": "TG|CA",
    "I-CeuI": "TTCGCTACC^TTAG/GACCGTTATAGTTA",
    "I-SceI": "ATTACCCTG^TTAT/CCCTA",
    "KasI": "G/GCGC^C",
    "KpnI": "G^GTAC/C",
    "MboII": "^N/NNNNNNNTCTTC",
    "MfeI": "C/AATT^G",
    "MluCI": "/AATT^",
    "MluI": "A/CGCG^T",
    "MlyI": "^/NNNNNGACTC",
    "MmeI": "^NN/NNNNNNNNNNNNNNNNNNGTYGGA",
    "MnlI": "^N/NNNNNNGAGG",
    "MscI": "TGG|CCA",
    "MseI": "T/TA^A",
    "MslI": "CAYNN|NNRTG",
    "MspA1I": "CMG|CKG",
    "MspJI": "/NNNN^NNNNNNNNNYNNG",
    "MwoI": "GCNN^NNN/NNGC",
    "NaeI": "GCC|GGC",
    "NarI": "GG/CG^CC",
    "Nb.BbvCI": "GC/TGAGG",
    "Nb.BsmI": "G/CATTC",
    "Nb.BsrDI": "/CATTGC",
    "Nb.BssSI": "C/TCGTG",
    "Nb.BtsI": "/CACTGC",
    "NciI": "CC/S^GG",
    "NcoI": "C/CATG^G",
    "NdeI": "CA/TA^TG",
    "NgoMIV": "G/CCGG^C",
    "NheI": "G/CTAG^C",
    "NlaIII": "^CATG/",
    "NlaIV": "GGN|NCC",
    "NmeAIII": "^NN/NNNNNNNNNNNNNNNNNNNCTCGGC",
    "NotI": "GC/GGCC^GC",
    "NruI": "TCG|CGA",
    "NsiI": "A^TGCA/T",
    "NspI": "R^CATG/Y",
    "Nt.AlwI": "^NNNNGATCC",
    "Nt.BbvCI": "GCTGA^GG",
    "Nt.BsmAI": "^NGAGAC",
    "Nt.BspQI": "^NGAAGAGC",
    "Nt.BstNBI": "^NNNNGACTC",
    "Nt.CviPII": "HGG^",
    "PacI": "TTA^AT/TAA",
    "PaqCI": "/NNNN^NNNNGCAGGTG",
    "PciI": "A/CATG^T",
    "PflFI": "GACN/N^NGTC",
    "Tth111I": "GACN/N^NGTC",
    "PflMI": "CCAN^NNN/NTGG",
    "PI-PspI": "ACCCATAATACCC^ATAA/TAGCTGTTTGCCA",
    "PI-SceI": "ATTACCTCTTTCTCC^GCAC/CCGACATAGAT",
    "PleI": "/N^NNNNGACTC",
    "PluTI": "G^GCGC/C",
    "PmeI": "GTTT|AAAC",
    "PmlI": "CAC|GTG",
    "PpuMI": "RG/GWC^CY",
    "PshAI": "GACNN|NNGTC",
    "PsiI-v2": "TTA|TAA",
    "PspGI": "/CCWGG^",
    "PspOMI": "G/GGCC^C",
    "PspXI": "VC/TCGA^GB",
    "PstI": "C^TGCA/G",
    "PvuI": "CG^AT/CG",
    "PvuII": "CAG|CTG",
    "RsaI": "GT|AC",
    "RsrII": "CG/GWC^CG",
    "SacI": "G^AGCT/C",
    "SacII": "CC^GC/GG",
    "SalI": "G/TCGA^C",
    "SapI": "/NNN^NGAAGAGC",
    "BspQI": "/NNN^NGAAGAGC",
    "Sau3AI": "/GATC^",
    "MboI": "/GATC^",
    "DpnII": "/GATC^",
    "Sau96I": "G/GNC^C",
    "SbfI": "CC^TGCA/GG",
    "ScaI": "AGT|ACT",
    "ScrFI": "CC/N^GG",
    "SexAI": "A/CCWGG^T",
    "SfaNI": "/NNNN^NNNNNGATGC",
    "SfcI": "C/TRYA^G",
    "SfiI": "GGCCN^NNN/NGGCC",
    "SfoI": "GGC|GCC",
    "SgrAI": "CR/CCGG^YG",
    "SmaI": "CCC|GGG",
    "SmlI": "C/TYRA^G",
    "SnaBI": "TAC|GTA",
    "SpeI": "A/CTAG^T",
    "SphI": "G^CATG/C",
    "SrfI": "GCCC|GGGC",
    "SspI": "AAT|ATT",
    "StuI": "AGG|CCT",
    "StyD4I": "/CCNGG^",
    "StyI": "C/CWWG^G",
    "SwaI": "ATTT|AAAT",
    "TaqI-v2": "T/CG^A",
    "TfiI": "G/AWT^C",
    "TseI": "G/CWG^C",
    "ApeKI": "G/CWG^C",
    "Tsp45I": "/GTSAC^",
    "TspRI": "^NNCASTGNN/",
    "XbaI": "T/CTAG^A",
    "XcmI": "CCANNNN^N/NNNNTGG",
    "XhoI": "C/TCGA^G",
    "PaeR7I": "C/TCGA^G",
    "XmaI": "C/CCGG^G",
    "TspMI": "C/CCGG^G",
    "XmnI": "GAANN|NNTTC",
    "ZraI": "GAC|GTC"
}
