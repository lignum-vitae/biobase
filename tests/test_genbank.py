from pathlib import Path

import pytest

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

SAMPLE_GENBANK_1 = """LOCUS       NC_012532           1079 bp    RNA     linear   VRL 28-JUL-2016
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
  PUBMED    17195954
REFERENCE   2  (bases 1 to 1079)
  AUTHORS   Lanciotti,R.S.
  TITLE     Direct Submission
  JOURNAL   Submitted (22-AUG-2009) Division of Vector-Borne Diseases,
            Centers for Disease Control and Prevention, 3156 Rampart Rd, Fort
            Collins, CO 80521, USA
FEATURES            Location/Qualifiers
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

SAMPLE_GENBANK_2 = """LOCUS       ADF90000            50 bp    DNA     circular INV 01-JAN-2023
DEFINITION  A second test record.
ACCESSION   ADF90000
VERSION     ADF90000.1  GI:100000000
KEYWORDS    second; test.
ORIGIN
        1 cgatcggatc gattcggact ggatcgatcg atcggatcga tcggatcgga
//
"""

# --- Pytest Fixtures ---


@pytest.fixture
def sample_gbk_file_single(tmp_path: Path) -> str:
    """A pytest fixture for a single-record GenBank file."""
    gbk_file = tmp_path / "test_single.gbk"
    gbk_file.write_text(SAMPLE_GENBANK_1)
    return str(gbk_file)


@pytest.fixture
def sample_gbk_file_multi(tmp_path: Path) -> str:
    """A pytest fixture for a multi-record GenBank file."""
    gbk_file = tmp_path / "test_multi.gbk"
    gbk_file.write_text(SAMPLE_GENBANK_1 + "\n" + SAMPLE_GENBANK_2)
    return str(gbk_file)


@pytest.fixture
def first_record(sample_gbk_file_single: str) -> GenBankRecord:
    """Fixture to get the GenBankRecord object for the first sample."""
    parser = GenBankParser(sample_gbk_file_single)
    # The parser is now an iterator, so we can use next() or list access
    return next(iter(parser))


# --- Tests for the GenBankParser class ---


def test_genbank_parser_read_all(sample_gbk_file_single: str):
    """Test that the GenBankParser can read a file correctly (read_all renamed)."""
    parser = GenBankParser(sample_gbk_file_single)
    content = parser.read_all()  # Updated function name
    assert content.startswith("LOCUS")
    assert "Zika virus" in content


def test_genbank_parser_split_into_blocks(sample_gbk_file_single: str):
    """Test the internal block splitting logic of the parser on a single record."""
    parser = GenBankParser(sample_gbk_file_single)
    content = parser.read_all()
    # Need to pass only the record part to _split_into_blocks
    blocks = list(parser._split_into_blocks(content.split("//")[0].strip()))

    assert len(blocks) == 10

    keys = [key for key, block in blocks]
    assert keys.count("REFERENCE") == 2
    assert keys.count("LOCUS") == 1

    expected_keys = [
        "LOCUS",
        "DEFINITION",
        "ACCESSION",
        "VERSION",
        "KEYWORDS",
        "SOURCE",
        "REFERENCE",
        "REFERENCE",
        "FEATURES",
        "ORIGIN",
    ]
    assert keys == expected_keys


def test_genbank_parser_split_into_records_single(sample_gbk_file_single: str):
    """Test splitting a file with a single record."""
    parser = GenBankParser(sample_gbk_file_single)
    content = parser.read_all()
    records = list(parser._split_into_records(content))
    assert len(records) == 1
    assert records[0].strip().endswith("//")


def test_genbank_parser_split_into_records_multi(sample_gbk_file_multi: str):
    """Test splitting a file with multiple records."""
    parser = GenBankParser(sample_gbk_file_multi)
    content = parser.read_all()
    records = list(parser._split_into_records(content))
    assert len(records) == 2
    assert records[0].strip().startswith("LOCUS       NC_012532")
    assert records[1].strip().startswith("LOCUS       ADF90000")


def test_genbank_parser_iter_multi_record(sample_gbk_file_multi: str):
    """Test the __iter__ method to ensure all records are yielded."""
    parser = GenBankParser(sample_gbk_file_multi)
    records: list[GenBankRecord] = list(parser)

    assert len(records) == 2

    # Check first record
    r1 = records[0]
    assert isinstance(r1, GenBankRecord)
    assert r1.id == "NC_012532"
    assert "LOCUS" in r1.entries and isinstance(r1.entries["LOCUS"], Locus)
    assert r1.entries["LOCUS"].length == 1079

    # Check second record
    r2 = records[1]
    assert isinstance(r2, GenBankRecord)
    assert r2.id == "ADF90000"
    assert r2.name == "A second test record."
    assert r2.entries["LOCUS"].length == 50
    assert len(r2.seq) == 50
    assert r2.entries["LOCUS"].topology == "circular"


# --- Tests for the GenBankRecord and component classes (Updated to use fixture) ---


def test_genbank_record_instantiation(first_record: GenBankRecord):
    """Test that the main GenBankRecord class can be created and has key attributes."""
    record = first_record
    assert record is not None
    assert isinstance(record, GenBankRecord)
    assert "LOCUS" in record.entries
    assert "ORIGIN" in record.entries

    # Test top-level attributes
    assert record.id == "NC_012532"
    assert record.name == "Zika virus, complete genome."
    assert record.seq.startswith("agttgttgat")

    assert repr(record).startswith("<GenBankRecord for 'NC_012532'")


def test_locus_parsing(first_record: GenBankRecord):
    """Verify that the LOCUS block is parsed correctly with the updated logic."""
    locus = first_record.entries["LOCUS"]
    assert isinstance(locus, Locus)
    assert locus.name == "NC_012532"
    assert locus.length == 1079
    assert locus.molecule_type == "RNA"
    assert locus.topology == "linear"
    assert locus.date == "28-JUL-2016"


def test_definition_parsing(first_record: GenBankRecord):
    """Verify that the DEFINITION block is parsed correctly."""
    definition = first_record.entries["DEFINITION"]
    assert isinstance(definition, Definition)
    assert definition.info == "Zika virus, complete genome."


def test_accession_parsing(first_record: GenBankRecord):
    """Verify that the ACCESSION block is parsed correctly."""
    accession = first_record.entries["ACCESSION"]
    assert isinstance(accession, Accession)
    assert accession.info == "NC_012532"


def test_version_parsing(first_record: GenBankRecord):
    """Verify that the VERSION block is parsed correctly."""
    version = first_record.entries["VERSION"]
    assert isinstance(version, Version)
    assert version.version == "NC_012532.1"
    assert version.gi == "254753235"


def test_origin_parsing(first_record: GenBankRecord):
    """Verify that the ORIGIN sequence is extracted and cleaned correctly."""
    origin = first_record.entries["ORIGIN"]
    assert isinstance(origin, Origin)

    expected_seq = "agttgttgatctgtgtgaatcagactgcgacagtcatggtaacagcagcaggaagaggcaggacgcttgcagccaagtcagcagctacagccctcgcaacgc"
    sequence = origin.sequence

    assert len(sequence) == len(expected_seq)
    assert sequence == expected_seq
    assert sequence.islower()


def test_features_parsing(first_record: GenBankRecord):
    """Verify the FEATURES block, using the new SingleFeature class."""
    features = first_record.entries["FEATURES"]
    assert isinstance(features, Features)

    assert len(features.entries) == 3

    # Test Feature 1: source
    f1 = features.entries[0]
    assert isinstance(f1, SingleFeature)
    assert f1.key == "source"
    assert f1.location == "1..1079"
    assert len(f1.qualifiers) == 4
    assert f1.qualifiers["organism"] == "Zika virus"
    assert f1.qualifiers["db_xref"] == "taxon:64320"

    # Test Feature 3: CDS (including /translation)
    f3 = features.entries[2]
    assert f3.key == "CDS"
    assert f3.location == "1..102"
    assert len(f3.qualifiers) == 5
    assert f3.qualifiers["product"] == "anchored capsid protein C"
    assert f3.qualifiers["translation"] == "MSNNQQKGGRLLQPSQ"


def test_unhandled_and_duplicate_entries(first_record: GenBankRecord):
    """
    Test that an entry without a dedicated parser class (e.g., KEYWORDS)
    is stored as a raw string, and that duplicate entries (REFERENCE) are
    concatenated.
    """
    record = first_record

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
