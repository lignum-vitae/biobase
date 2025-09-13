# Stndard library
from typing import Iterator


def main():
    fasta_seqs = """
>CAA39742.1 cytochrome b (mitochondrion) [Sus scrofa]
MTNIRKSHPLMKIINNAFIDLPAPSNISSWWNFGSLLGICLILQILTGLFLAMHYTSDTTTAFSSVTHIC
RDVNYGWVIRYLHANGASMFFICLFIHVGRGLYYGSYMFLETWNIGVVLLFTVMATAFMGYVLPWGQMSF
WGATVITNLLSAIPYIGTDLVEWIWGGFSVDKATLTRFFAFHFILPFIITALAAVHLMFLHETGSNNPTG
ISSDMDKIPFHPYYTIKDILGALFMMLILLILVLFSPDLLGDPDNYTPANPLNTPPHIKPEWYFLFAYAI
LRSIPNKLGGVLALVASILILILMPMLHTSKQRGMMFRPLSQCLFWMLVADLITLTWIGGQPVEHPFIII
GQLASILYFLIILVLMPITSIIENNLLKW

>BAA85863.1 cytochrome b, partial (mitochondrion) [Rattus rattus]
MTNIRKSHPLIKIINHSFIDLPAPSNISSWWNFGSLLGVCLMVQIITGLFLAMHYTSDTLTAFSSVTHIC
RDVNYGWLIRYLHANGASMFFICLFLHVGRGMYYGSYTFLETWNIGIILLFAVMATAFMGYVLPWGQMSF
WGATVITNLLSAIPYIGTTLVEWIWGGFSVDKATLTRFFAFHFILPFIIAALAIVHLLFLHETGSNNPTG
LNSDADKIPFHPYYTIKDLLGVFMLLLFLMTLVLFFPDLLGDPDNYTPANPLNTPPHIKPEWYFLFAYAI
LRSIPNKLGGVVALVLSILILAFLPFLHTSKQRSLTFRPITQILYWILVANLFILTWIGGQPVEHPFIII
GQLASISYFSIILILMPISGIIEDKMLKWN
"""
    records = fasta_parser(fasta_seqs)
    for r in records:
        print(r.id)
        print(r.seq)


class FastaRecord:
    def __init__(self, header, sequence) -> None:
        parts = header.lstrip(">").split(maxsplit=1)
        self.id = parts[0]
        self.name = parts[1] if len(parts) > 1 else ""
        self.seq = sequence

    def __repr__(self) -> str:
        return (
            f"FastaRecord(id={self.id!r}, name={self.name!r}, seq_len={len(self.seq)})"
        )

    def __str__(self) -> str:
        preview = self.seq[:20] + ("..." if len(self.seq) > 20 else "")
        return f"{self.id} | {self.name} | {preview}"

    def length(self) -> int:
        return len(self.seq)


class FastaFileParser:
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath

    def __iter__(self) -> Iterator[FastaRecord]:
        with open(self.filepath, "r") as file:
            temp = [x.strip() for x in file.readlines()]
        fast_index = [i for i, j in enumerate(temp) if ">" in j]
        if len(fast_index) == 0:
            raise ValueError("Failed to parse file due to improper fasta format")
        elif len(fast_index) == 1:
            header = temp[fast_index[0]]
            seq = "".join(temp[fast_index[0] + 1 :])
            yield FastaRecord(header, seq)
        elif len(fast_index) >= 2:
            for i, temp_index in enumerate(fast_index[:-1]):
                header = temp[temp_index]
                seq = "".join(temp[temp_index + 1 : fast_index[i + 1]])
                yield FastaRecord(header, seq)

            header = temp[fast_index[-1]]
            seq = "".join(temp[fast_index[-1] + 1 :])
            yield FastaRecord(header, seq)


class FastaParser:
    def __init__(self, reads: str) -> None:
        self.reads = reads

    def __iter__(self) -> Iterator[FastaRecord]:
        temp = [x.strip() for x in self.reads.split("\n")]
        fast_index = [i for i, j in enumerate(temp) if ">" in j]
        if len(fast_index) == 0:
            raise ValueError("Failed to parse file due to improper fasta format")
        elif len(fast_index) == 1:
            header = temp[fast_index[0]]
            seq = "".join(temp[fast_index[0] + 1 :])
            yield FastaRecord(header, seq)
        elif len(fast_index) >= 2:
            for i, temp_index in enumerate(fast_index[:-1]):
                header = temp[temp_index]
                seq = "".join(temp[temp_index + 1 : fast_index[i + 1]])
                yield FastaRecord(header, seq)
            header = temp[fast_index[-1]]
            seq = "".join(temp[fast_index[-1] + 1 :])
            yield FastaRecord(header, seq)


def fasta_parser(fasta: str) -> list[FastaRecord]:
    return list(FastaParser(fasta))


def fasta_file_parser(file_path: str) -> list[FastaRecord]:
    return list(FastaFileParser(file_path))


if __name__ == "__main__":
    main()
