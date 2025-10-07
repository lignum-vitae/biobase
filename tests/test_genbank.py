from pathlib import Path

import pytest

# Import the updated classes from your parser module
from biobase.parser.genbank import (
    Accession,
    Definition,
    Features,
    GenBankParser,
    GenBankRecord,
    Locus,
    Origin,
    SingleFeature,
    Version,
)

SAMPLE_GENBANK = """LOCUS       NC_012532               1079 bp    RNA     linear   VRL 28-JUL-2016
DEFINITION  Zika virus, complete genome.
ACCESSION   NC_012532
VERSION     NC_012532.1  GI:254753235
KEYWORDS    .
SOURCE      Zika virus
  ORGANISM  Zika virus
            Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Flasuviricetes;
            Amarillovirales; Flaviviridae; Flavivirus.
REFERENCE   1  (bases 1 to 1079)
  AUTHORS   Kuno,G. and Chang,G.J.
  TITLE     Full-length sequencing and genomic characterization of Bagaza,
            Kedougou, and Zika viruses
  JOURNAL   Arch. Virol. 152 (4), 687-696 (2007)
   PUBMED   17195954
REFERENCE   2  (bases 1 to 1079)
  AUTHORS   Lanciotti,R.S.
  TITLE     Direct Submission
  JOURNAL   Submitted (22-AUG-2009) Division of Vector-Borne Diseases,
            Centers for Disease Control and Prevention, 3156 Rampart Rd, Fort
            Collins, CO 80521, USA
FEATURES             Location/Qualifiers
     source          1..1079
                     /organism="Zika virus"
                     /mol_type="genomic RNA"
                     /strain="MR 766"
                     /db_xref="taxon:64320"
     gene            1..102
                     /gene="ANK"
     CDS             1..102
                     /gene="ANK"
                     /codon_start=1
                     /product="anchored capsid protein C"
                     /protein_id="YP_009228357.1"
                     /translation="MSNNQQKGGRLLQPSQ"
ORIGIN
        1 agttgttgat ctgtgtgaat cagactgcga cagtcatggt aacagcagca ggaagaggca
       61 ggacgcttgc agccaagtcag cagctacagc cctcgcaacg c
//
"""


@pytest.fixture
def sample_gbk_file(tmp_path: Path) -> str:
    """A pytest fixture to write the sample GenBank content to a temporary file."""
    gbk_file = tmp_path / "test.gbk"
    gbk_file.write_text(SAMPLE_GENBANK)
    return str(gbk_file)


# --- Tests for the GenBankParser class ---


def test_genbank_parser_read(sample_gbk_file: str):
    """Test that the GenBankParser can read a file correctly."""
    parser = GenBankParser(sample_gbk_file)
    content = parser.read()
    assert content.startswith("LOCUS")
    assert "Zika virus" in content


def test_genbank_parser_split_into_blocks(sample_gbk_file: str):
    """Test the internal block splitting logic of the parser."""
    parser = GenBankParser(sample_gbk_file)
    content = parser.read()
    blocks = list(parser._split_into_blocks(content))

    assert len(blocks) == 10

    keys = [key for key, block in blocks]
    assert keys.count("REFERENCE") == 2
    assert keys.count("LOCUS") == 1

    # This is the corrected list, with "REFERENCE" appearing twice.
    expected_keys = [
        "LOCUS",
        "DEFINITION",
        "ACCESSION",
        "VERSION",
        "KEYWORDS",
        "SOURCE",
        "REFERENCE",  # First occurrence
        "REFERENCE",  # Second occurrence
        "FEATURES",
        "ORIGIN",
    ]

    assert keys == expected_keys


# --- Tests for the GenBankRecord and component classes ---


def test_genbank_record_instantiation(sample_gbk_file: str):
    """Test that the main GenBankRecord class can be created and has key attributes."""
    record = GenBankRecord(sample_gbk_file)
    assert record is not None
    assert isinstance(record, GenBankRecord)
    assert "LOCUS" in record.entries
    assert "ORIGIN" in record.entries

    # Test top-level attributes
    assert record.id == "NC_012532"
    assert record.name == "Zika virus, complete genome."
    assert record.seq.startswith("agttgttgat")

    assert repr(record).startswith("<GenBankRecord for 'NC_012532'")


def test_locus_parsing(sample_gbk_file: str):
    """Verify that the LOCUS block is parsed correctly with the updated logic."""
    record = GenBankRecord(sample_gbk_file)
    locus = record.entries["LOCUS"]
    assert isinstance(locus, Locus)
    assert locus.name == "NC_012532"
    assert locus.length == 1079
    assert locus.molecule_type == "RNA"
    assert locus.topology == "linear"
    assert locus.date == "28-JUL-2016"


def test_definition_parsing(sample_gbk_file: str):
    """Verify that the DEFINITION block is parsed correctly."""
    record = GenBankRecord(sample_gbk_file)
    definition = record.entries["DEFINITION"]
    assert isinstance(definition, Definition)
    assert definition.info == "Zika virus, complete genome."


def test_accession_parsing(sample_gbk_file: str):
    """Verify that the ACCESSION block is parsed correctly."""
    record = GenBankRecord(sample_gbk_file)
    accession = record.entries["ACCESSION"]
    assert isinstance(accession, Accession)
    assert accession.info == "NC_012532"


def test_version_parsing(sample_gbk_file: str):
    """Verify that the VERSION block is parsed correctly."""
    record = GenBankRecord(sample_gbk_file)
    version = record.entries["VERSION"]
    assert isinstance(version, Version)
    assert version.version == "NC_012532.1"
    assert version.gi == "254753235"


def test_origin_parsing(sample_gbk_file: str):
    """Verify that the ORIGIN sequence is extracted and cleaned correctly."""
    record = GenBankRecord(sample_gbk_file)
    origin = record.entries["ORIGIN"]
    assert isinstance(origin, Origin)

    expected_seq = "agttgttgatctgtgtgaatcagactgcgacagtcatggtaacagcagcaggaagaggcaggacgcttgcagccaagtcagcagctacagccctcgcaacgc"
    sequence = origin.sequence

    assert len(sequence) == len(expected_seq)
    assert sequence == expected_seq
    assert "1" not in sequence
    assert " " not in sequence
    assert sequence.islower()


def test_features_parsing(sample_gbk_file: str):
    """Verify the FEATURES block, using the new SingleFeature class."""
    record = GenBankRecord(sample_gbk_file)
    features = record.entries["FEATURES"]
    assert isinstance(features, Features)

    assert len(features.entries) == 3

    # Test Feature 1: source
    f1 = features.entries[0]
    assert isinstance(f1, SingleFeature)  # Check for the new class name
    assert f1.key == "source"
    assert f1.location == "1..1079"
    assert len(f1.qualifiers) == 4
    assert f1.qualifiers["organism"] == "Zika virus"
    assert f1.qualifiers["db_xref"] == "taxon:64320"

    # Test Feature 2: gene
    f2 = features.entries[1]
    assert f2.key == "gene"
    assert f2.location == "1..102"
    assert f2.qualifiers["gene"] == "ANK"

    # Test Feature 3: CDS
    f3 = features.entries[2]
    assert f3.key == "CDS"
    assert f3.location == "1..102"
    assert len(f3.qualifiers) == 5
    assert f3.qualifiers["product"] == "anchored capsid protein C"
    assert f3.qualifiers["protein_id"] == "YP_009228357.1"
    assert f3.qualifiers["translation"] == "MSNNQQKGGRLLQPSQ"


def test_unhandled_and_duplicate_entries(sample_gbk_file: str):
    """
    Test that an entry without a dedicated parser class (e.g., KEYWORDS)
    is stored as a raw string, and that duplicate entries (REFERENCE) are
    concatenated.
    """
    record = GenBankRecord(sample_gbk_file)

    # KEYWORDS has no special class, so it should be a string
    keywords = record.entries["KEYWORDS"]
    assert isinstance(keywords, str)
    assert keywords.strip() == "KEYWORDS    ."

    # The two REFERENCE blocks should be concatenated into a single string
    reference = record.entries["REFERENCE"]
    assert isinstance(reference, str)
    assert "AUTHORS   Kuno,G. and Chang,G.J." in reference
    assert "AUTHORS   Lanciotti,R.S." in reference
    assert reference.count("REFERENCE") == 2
