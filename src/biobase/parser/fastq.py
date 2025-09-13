# Stndard library
from typing import Iterator

# Third-party dependencies
import numpy as np

# Internal dependencies
from biobase.analysis import Dna
from .fasta import FastaRecord


def main() -> None:
    fastq_seq = """@2fa9ee19-5c51-4281-abdd-eac8663f9b49 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=109831 ch=33 start_time=2019-09-12T12:05:03Z
    CGGTAGCCAGCTGCGTTCAGTATGGAAGATTTGATTTGTTTAGCGATCGCCATACTACCGTGACAAGAAAGTTGTCAGTCTTTGTGACTTGCCTGTCGCTCTATCTTCCAGACTCCTTGGTCCGTGTTCAATCCCGGTAGTAGCGACGGGCGGTGTATGTATTATCAGCGCAACAGAAACAAAGACACC
    +
    +&&-&%$%%$$$#)33&0$&%$''*''%$#%$%#+-5/---*&&%$%&())(&$#&,'))5769*+..*&(+28./#&1228956:7674';:;80.8>;91;>?B=%.**==?(/'($$$$*'&'**%&/));807;3A=;88>=?9498++0%"%%%%'#&5/($0.$2%&0'))*'%**&)(.%&&
    @1f9ca490-2f25-484a-8972-d60922b41f6f runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=106343 ch=28 start_time=2019-09-12T12:05:07Z
    GATGCATACTTCGTTCGATTTCGTTTCAACTGGACAACCTACCGTGACAAAGAAAGTTGTCGATGCTTTGTGACTTGCTGTCCTCTATCTTCAGACTCCTTGGTCCATTTCAAGACCAAACAATCAGTAGTAGCGACGGGCGGTGTGGCAATATCGCTTTCAACGAAACACAAAGAAT
    +
    &%&%''&'+,005<./'%*-)%(#$'$#$%&&'('$$..74483=.0412$*/,)9/194/+('%%(+1+'')+,-&,>;%%.*@@D>>?)3%%296070717%%16;<'<236800%(,734(0$7769879@;?8)09:+/4'1+**7<<4.4,%%(.)##%&'(&&%*++'&#%$
    @06936a64-6c08-40e9-8a10-0fbc74812c89 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=83531 ch=23 start_time=2019-09-12T12:03:50Z
    GTTTTGTCGCTGCGTTCAGTTTATGGGTGCGGGTGTTATGATGCTTCGCTTTACGTGACAAGAAAGTTAGTAGATTGTCTTTATGTTTCTGTGGTGCTGATATTGCCACACCGCCCGATAGCTCTACCGATTGAAACACGGACCAAGGAATCGGAAATGTAGGCGAGCAGGCCGTCCTGAACACCCATTAACTTTCTTGTC
    +
    $&'((&%$$$.$2/=-*#'.2'&&##$$#$#&&(&+-%'(%&#"###""$$%#)%,+)+&'(,&%*((%%&%$%'+),,+,,&%$')1+*$.&+6*+(*%(&'*(''&%*+,*)('%#$$$%,$&&'&)))12)*&/*,364$%$%))$'')#%%&%$#$%$$#('$(%$%$%%$$*$&$%)''%%$$&'&$)+2++,)&%
    @d6a555a1-d8dd-4e55-936f-ade7c78d9d38 runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=112978 ch=97 start_time=2019-09-12T12:03:49Z
    CGTATGCTTTGAGATTCATTCAGGAGGCGGGTATTTGCTCGATCATACCATACGTGGCAAGAAAGTTGTCAGTGTCTTTGTGTTTCTCTGTGGTGCGCGATATTGCCACGCCCGTCGCTACACCGATTGAAACACGGACCGAAGTCTGAAGATAGAGCGACGAGCGAAGTCACAAAGGAACTAGAGCAACTTTTTATC
    +
    #$%%%%''(($$%$*-&%$%)%*'%(+($)(%$,.)##$&$$#$$&('(%&%%%%#$$%(&*('('+18/(6?65510+))'--*&&$$$,*+;/+%%&&''13&%&%(133<;9=/.2*$*657,0*&(237'85;A1/$$%'7:;;:<2:..%$)0,*.)(1)1&&1+-$$,-&(-&&####%%98:AHFEB4(%,
    @91ca9c6c-12fe-4255-83cc-96ba4d39ac4b runid=f53ee40429765e7817081d4bcdee6c1199c2f91d sampleid=18S_amplicon read=110811 ch=113 start_time=2019-09-12T12:04:28Z
    CGGTGTACTTCGTTCCAGCTAGATTTGGGTGCATGACCATACCGTGACAAGAAAGTTGTCGGTATCTTTGTGTTTCTGTTGGTGCTGATATTGCCGCACCGCCCGTCGCTACACCGATTGTTCTGTTGGTCTTGAAACACGGACCAGGGTCTAGAGCAG
    +
    %$&'$&'&&'0,42%*$&&%$%$#$)$*+,'($&))(*$%$%'-8644(()-&'%&*'')%*('579:?.*,9:+)1-9.'(7491:7,(52.11'7;:<===E@;>448,,(%*.''*,%&$-.;<:;66138/**,2?8<:**'%&)%&#$&&,,'&"""
    fastq: FastqParser = FastqParser(fastq_seq)
    print(f"Total reads: {fastq.count_reads()}")
    for read in fastq:
        print(read)
        print(" AvgQ:", read.average_quality())
        print(" Fasta:\n", read.convert_to_fasta())


class FastqRecord:
    def __init__(self, id: str, seq: str, separator: str, quality: str) -> None:
        # Validation is done at the file level
        self.id = id
        self.seq = seq
        self.separator = separator
        self.quality = quality

    # Behavior in case of print or repr
    def __repr__(self) -> str:
        return f"FastqRecord(id={self.id!r}, seq_len={len(self.seq)})"

    def __str__(self) -> str:
        preview = self.seq[:20] + ("..." if len(self.seq) > 20 else "")
        return f"{self.id} | {preview}"

    # FastqRecord level utility
    def length(self) -> int:
        return len(self.seq)

    def convert_to_fasta(self) -> str:
        return f">{self.id}\n{self.seq}"

    def phred_scores(self) -> np.ndarray:
        return np.fromiter((ord(ch) - 33 for ch in self.quality), dtype=np.int16)

    def average_quality(self) -> float:
        scores: np.ndarray = self.phred_scores()
        return sum(scores) / scores.size if scores.size > 0 else 0.0


class FastqParserBase:
    # Operations on all reads
    def count_reads(self) -> int:
        return sum(1 for _ in self)

    def filter_reads(self, min_avg_quality: float) -> Iterator[FastqRecord]:
        for read in self:
            if read.average_quality() >= min_avg_quality:
                yield read

    def to_fasta(self) -> list[FastaRecord]:
        return list(self.to_fasta_iter())

    def to_fasta_iter(self) -> Iterator[FastaRecord]:
        for read in self:
            yield FastaRecord(read.id, read.seq)

    def to_fasta_file(self, out_path: str) -> None:
        with open(out_path, "w") as file:
            for read in self:
                file.write(read.convert_to_fasta() + "\n")

    def read_lengths(self) -> np.ndarray:
        lengths: list[int] = []
        for read in self:
            lengths.append(read.length())
        # Using numpy arrays as the return because it will be more efficient for the further analysis
        lengths_array: np.ndarray = np.array(lengths)
        return lengths_array


class FastqFileParser(FastqParserBase):
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath

    def __iter__(self) -> Iterator[FastqRecord]:
        with open(self.filepath, "r") as f:
            for line_number, line in enumerate(f, start=1):
                line = line.strip()
                # Keep track on where we are on the file
                pos = line_number % 4

                # Just in case there are no empty lines before the first line
                if pos == 1:
                    if not line.startswith("@"):
                        raise ValueError(
                            f"Invalid Fastq format Expecetd '@' at line {line_number} got {line!r}"
                        )
                    read_identifier = line.lstrip("@")
                elif pos == 2:
                    nucleotide_sequence = line.strip()
                    nucleotide_sequence = Dna._validate_dna_sequence(line)
                # Again to make sure that there are no empty lines in between
                elif pos == 3:
                    if not line.startswith("+"):
                        raise ValueError(
                            f"Invalid Fastq format Expecetd '+' at line {line_number} got {line!r}"
                        )
                    separator = line.strip()
                elif pos == 0:
                    read_quality = line.strip()
                    yield FastqRecord(
                        read_identifier, nucleotide_sequence, separator, read_quality
                    )

    def count_reads(self) -> int:
        return super().count_reads()

    def filter_reads(self, min_avg_quality: float) -> Iterator[FastqRecord]:
        yield from super().filter_reads(min_avg_quality)

    def to_fasta(self) -> list[FastaRecord]:
        return super().to_fasta()

    def to_fasta_iter(self) -> Iterator[FastaRecord]:
        yield from super().to_fasta_iter()

    def to_fasta_file(self, out_path: str) -> None:
        return super().to_fasta_file(out_path)

    def read_lengths(self) -> np.ndarray:
        return super().read_lengths()


class FastqParser(FastqParserBase):
    def __init__(self, reads: str) -> None:
        self.reads = reads

    def __iter__(self) -> Iterator[FastqRecord]:
        for line_number, line in enumerate(self.reads.splitlines(), start=1):
            line = line.strip()
            # Keep track on where we are on the file
            pos = line_number % 4

            # Just in case there are no empty lines before the first line
            if pos == 1:
                if not line.startswith("@"):
                    raise ValueError(
                        f"Invalid Fastq format Expecetd '@' at line {line_number} got {line!r}"
                    )
                read_identifier = line.lstrip("@")
            elif pos == 2:
                nucleotide_sequence = line.strip()
                nucleotide_sequence = Dna._validate_dna_sequence(line)
            # Again to make sure that there are no empty lines in between
            elif pos == 3:
                if not line.startswith("+"):
                    raise ValueError(
                        f"Invalid Fastq format Expecetd '+' at line {line_number} got {line!r}"
                    )
                separator = line.strip()
            elif pos == 0:
                read_quality = line.strip()
                yield FastqRecord(
                    read_identifier, nucleotide_sequence, separator, read_quality
                )

    def count_reads(self) -> int:
        return super().count_reads()

    def filter_reads(self, min_avg_quality: float) -> Iterator[FastqRecord]:
        yield from super().filter_reads(min_avg_quality)

    def to_fasta(self) -> list[FastaRecord]:
        return super().to_fasta()

    def to_fasta_iter(self) -> Iterator[FastaRecord]:
        yield from super().to_fasta_iter()

    def to_fasta_file(self, out_path: str) -> None:
        return super().to_fasta_file(out_path)

    def read_lengths(self) -> np.ndarray:
        return super().read_lengths()


def fastq_parser(fastq: str) -> list[FastqRecord]:
    return list(FastqParser(fastq))


def fastq_file_parser(file_path: str) -> list[FastqRecord]:
    return list(FastqFileParser(file_path))


if __name__ == "__main__":
    main()
