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
    def __init__(self, header, sequence):
        parts = header.lstrip(">").split(maxsplit=1)
        self.id = parts[0]
        self.name = parts[1] if len(parts) > 1 else ""
        self.seq = sequence

    def __repr__(self):
        return (
            f"FastaRecord(id={self.id!r}, name={self.name!r}, seq_len={len(self.seq)})"
        )

    def __str__(self):
        preview = self.seq[:20] + ("..." if len(self.seq) > 20 else "")
        return f"{self.id} | {self.name} | {preview}"


def fasta_file_parser(file_path: str) -> list[FastaRecord]:
    with open(file_path, "r") as file:
        temp = [x.strip() for x in file.readlines()]
    fast_index = [i for i, j in enumerate(temp) if ">" in j]
    if len(fast_index) == 1:
        header = temp[fast_index[0]]
        seq = "".join(temp[fast_index[0] + 1 :])
        return [FastaRecord(header, seq)]
    if len(fast_index) >= 2:
        records = []
        for i, temp_index in enumerate(fast_index[:-1]):
            header = temp[temp_index]
            seq = "".join(temp[temp_index + 1 : fast_index[i + 1]])
            records.append(FastaRecord(header, seq))

        header = temp[fast_index[-1]]
        seq = "".join(temp[fast_index[-1] + 1 :])
        records.append(FastaRecord(header, seq))
        return records
    raise ValueError("Failed to parse file due to improper fasta format")


def fasta_parser(fasta: str) -> list[FastaRecord]:
    temp = [x.strip() for x in fasta.split("\n")]
    fast_index = [i for i, j in enumerate(temp) if ">" in j]
    if len(fast_index) == 1:
        header = temp[fast_index[0]]
        seq = "".join(temp[fast_index[0] + 1 :])
        return [FastaRecord(header, seq)]
    if len(fast_index) >= 2:
        records = []
        for i, temp_index in enumerate(fast_index[:-1]):
            header = temp[temp_index]
            seq = "".join(temp[temp_index + 1 : fast_index[i + 1]])
            records.append(FastaRecord(header, seq))

        header = temp[fast_index[-1]]
        seq = "".join(temp[fast_index[-1] + 1 :])
        records.append(FastaRecord(header, seq))
        return records
    raise ValueError("Failed to parse file due to improper fasta format")


if __name__ == "__main__":
    main()
