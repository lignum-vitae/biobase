import pytest
from pathlib import Path
from io import StringIO

from biobase.read import fasta_parser, fasta_file_parser

SAMPLE_FASTA = """>CAA39742.1 cytochrome b (mitochondrion) [Sus scrofa]
MTNIRKSHPLMKIINNAFIDLPAPSNISSWWNFGSLLGICLILQILTGLFLAMHYTSDTTTAFSSVTHIC

>BAA85863.1 cytochrome b, partial (mitochondrion) [Rattus rattus]
MTNIRKSHPLIKIINHSFIDLPAPSNISSWWNFGSLLGVCLMVQIITGLFLAMHYTSDTLTAFSSVTHIC
"""


def test_fasta_parser_multi_record():
    records = fasta_parser(SAMPLE_FASTA)
    assert len(records) == 2

    # Test first record
    r1 = records[0]
    assert r1.id == "CAA39742.1"
    assert r1.name == "cytochrome b (mitochondrion) [Sus scrofa]"
    assert r1.seq.startswith("MTNIRKSHPLMKIINNAF")

    # Test second record
    r2 = records[1]
    assert r2.id == "BAA85863.1"
    assert r2.name == "cytochrome b, partial (mitochondrion) [Rattus rattus]"
    assert r2.seq.startswith("MTNIRKSHPLIKIINHSF")


def test_fasta_record_str_and_repr():
    record = fasta_parser(SAMPLE_FASTA)[0]
    s = str(record)
    r = repr(record)
    assert record.id in s
    assert record.name in s
    assert "FastaRecord" in r
    assert str(len(record.seq)) in r  # seq_len appears in repr


def test_fasta_parser_single_record():
    single_fasta = """>NP_000257 TP53 protein [Homo sapiens]
MEEPQSDPSV"""
    records = fasta_parser(single_fasta)
    assert len(records) == 1
    r = records[0]
    assert r.id == "NP_000257"
    assert r.name == "TP53 protein [Homo sapiens]"
    assert r.seq == "MEEPQSDPSV"


def test_fasta_file_parser(tmp_path):
    # Create a temporary FASTA file
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(SAMPLE_FASTA)

    records = fasta_file_parser(str(fasta_file))
    assert len(records) == 2
    assert records[0].id == "CAA39742.1"
    assert records[1].id == "BAA85863.1"


def test_invalid_fasta_raises_value_error(tmp_path):
    bad_file = tmp_path / "bad.fasta"
    bad_file.write_text("This is not a FASTA file")

    with pytest.raises(ValueError):
        fasta_file_parser(str(bad_file))

    with pytest.raises(ValueError):
        fasta_parser("No headers here, just text")
