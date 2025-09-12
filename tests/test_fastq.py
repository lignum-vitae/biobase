import numpy as np
import pytest

from biobase.parser import (
    FastqFileParser,
    FastqParser,
    FastqRecord,
    fastq_parser,
    fastq_file_parser,
)

SAMPLE_FASTQ = """@2fa9ee19-5c51-4281-abdd-eac8663f9b49 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=109831 ch=33 start_time=2019-09-12T12:05:03Z
CGGTAGCCAGCTGCGTTCAGTATGGAAGATTTGATTTGTTTAGCGATCGCCATACTACCGTGACAAGAAAGTTGTCAGTCTTTGTGACTTGCCTGTCGCTCTATCTTCCAGACTCCTTGGTCCGTGTTCAATCCCGGTAGTAGCGACGGGCGGTGTATGTATTATCAGCGCAACAGAAACAAAGACACC
+
+&&-&%$%%$$$#)33&0$&%$''*''%$#%$%#+-5/---*&&%$%&())(&$#&,'))5769*+..*&(+28./#&1228956:7674';:;80.8>;91;>?B=%.**==?(/'($$$$*'&'**%&/));807;3A=;88>=?9498++0%"%%%%'#&5/($0.$2%&0'))*'%**&)(.%&&
@1f9ca490-2f25-484a-8972-d60922b41f6f runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=106343 ch=28 start_time=2019-09-12T12:05:07Z
GATGCATACTTCGTTCGATTTCGTTTCAACTGGACAACCTACCGTGACAAAGAAAGTTGTCGATGCTTTGTGACTTGCTGTCCTCTATCTTCAGACTCCTTGGTCCATTTCAAGACCAAACAATCAGTAGTAGCGACGGGCGGTGTGGCAATATCGCTTTCAACGAAACACAAAGAAT
+
&%&%''&'+,005<./'%*-)%(#$'$#$%&&'('$$..74483=.0412$*/,)9/194/+('%%(+1+'')+,-&,>;%%.*@@D>>?)3%%296070717%%16;<'<236800%(,734(0$7769879@;?8)09:+/4'1+**7<<4.4,%%(.)##%&'(&&%*++'&#%$
"""


def test_fastq_parser_class_multi_record():
    records = list(FastqParser(SAMPLE_FASTQ))
    assert len(records) == 2
    # Record 1
    r1 = records[0]
    assert r1.id.startswith("2fa9ee19-5c51")
    assert "runid=" in r1.id
    assert r1.seq.startswith("CGGTAGCCAGCTGCG")
    assert len(r1.seq) == len(r1.quality)  # always true for FASTQ
    assert r1.length() == len(r1.seq)

    # Check phred conversion doesn’t crash
    scores1 = r1.phred_scores()
    assert isinstance(scores1, np.ndarray)
    assert scores1.min() >= 0  # Phred scores can’t be negative
    assert scores1.max() < 100  # sanity check

    # Record 2
    r2 = records[1]
    assert r2.id.startswith("1f9ca490-2f25")
    assert r2.seq.startswith("GATGCATACTTCGTT")
    assert len(r2.seq) == len(r2.quality)
    assert r2.length() == len(r2.seq)

    scores2 = r2.phred_scores()
    assert isinstance(scores2, np.ndarray)
    assert scores2.min() >= 0
    assert scores2.max() < 100


def test_fastq_parser_multi_record():
    records = fastq_parser(SAMPLE_FASTQ)
    assert len(records) == 2
    # Record 1
    r1 = records[0]
    assert r1.id.startswith("2fa9ee19-5c51")
    assert "runid=" in r1.id
    assert r1.seq.startswith("CGGTAGCCAGCTGCG")
    assert len(r1.seq) == len(r1.quality)  # always true for FASTQ
    assert r1.length() == len(r1.seq)


def test_fastq_record_and_repr():
    # repr and print have the same behavior for the FastqParser - only repr is defined in the class
    # They include the id of the read (first line without @) and the read seq length
    record = list(FastqParser(SAMPLE_FASTQ))[0]
    s = str(record)
    r = repr(record)
    assert record.id in s, r
    assert str(record.length()) in r


def test_fastq_parser_single_record():
    single_fastq = """@2fa9ee19-5c51-4281-abdd-eac86
CGGTAGCCAGCTGCGTTCAGTATG
+
%%%+++'''@@@???<<<??????"""
    records = list(FastqParser(single_fastq))
    assert len(records) == 1
    r: FastqRecord = records[0]
    assert r.id == "2fa9ee19-5c51-4281-abdd-eac86"
    assert r.seq == "CGGTAGCCAGCTGCGTTCAGTATG"


def test_fastq_file_parser_class(tmp_path):
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text(SAMPLE_FASTQ)
    parser = FastqFileParser(str(fastq_file))  # Produce a read object and can be listed
    records = list(parser)
    assert len(records) == 2
    # Record 1
    r1 = records[0]
    assert r1.id.startswith("2fa9ee19-5c51")
    assert "runid=" in r1.id
    # Record 2
    r2 = records[1]
    assert r2.id.startswith("1f9ca490-2f25")
    assert r2.seq.startswith("GATGCATACTTCGTT")
    assert len(r2.seq) == len(r2.quality)
    assert r2.length() == len(r2.seq)


def test_fastq_file_parser(tmp_path):
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text(SAMPLE_FASTQ)
    records = fastq_file_parser(
        str(fastq_file)
    )  # Produce a read object and can be listed
    assert len(records) == 2
    # Record 1
    r1 = records[0]
    assert r1.id.startswith("2fa9ee19-5c51")
    assert "runid=" in r1.id
    # Record 2
    r2 = records[1]
    assert r2.id.startswith("1f9ca490-2f25")
    assert r2.seq.startswith("GATGCATACTTCGTT")
    assert len(r2.seq) == len(r2.quality)
    assert r2.length() == len(r2.seq)


def test_to_fasta():
    parser = FastqParser(SAMPLE_FASTQ)
    fasta_record_arr = parser.to_fasta()

    assert "2fa9ee19-5c51" in fasta_record_arr[0].id
    assert "1f9ca490-2f25" in fasta_record_arr[1].id

def test_to_fasta_file(tmp_path):
    out_path = tmp_path / "out.fasta"
    # Get the first record
    parser = FastqParser(SAMPLE_FASTQ)
    parser.to_fasta_file(out_path)
    fasta_content = out_path.read_text()
    # make sure the headers of two reads is there
    assert ">2fa9ee19-5c51" in fasta_content
    assert ">1f9ca490-2f25" in fasta_content

    assert "CGGTA" in fasta_content
    assert "GATGC" in fasta_content


def test_invalid_fastq_header_raises_value_eror():
    bad_fasta = """This is a bad fastq"""
    parser = FastqParser(bad_fasta)
    with pytest.raises(ValueError):
        list(parser)
