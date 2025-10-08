import re
from collections.abc import Generator
from pathlib import Path
from typing import Any, Iterator


def main() -> None:
    """
    Main function to demonstrate the GenBank parser.
    It downloads a sample file, parses it using the new iterable parser,
    prints the contents of all records, and cleans up the file afterward.
    """
    # import requests
    # import os

    # # URL for a sample GenBank file with multiple records (e.g., several accessions)
    # multi_record_ids = "NG_017013.2,DQ366810.1,X52541.1"
    # url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={multi_record_ids}&rettype=gb&retmode=text"
    # file_path = "temp_multi_record.gb"

    # try:
    #     # 1. Download the file
    #     print(f"Downloading sample file with multiple records from NCBI: {multi_record_ids}...")
    #     response = requests.get(url, timeout=10)
    #     response.raise_for_status()  # Raise an exception for bad status codes

    #     with open(file_path, "w") as f:
    #         f.write(response.text)
    #     print(f"File saved to '{file_path}'")

    #     # 2. Parse the file and iterate over records using the GenBankParser
    #     print("\n--- Parsing GenBank Records ---")
    #     parser = GenBankParser(file_path)

    #     total_records = 0

    #     # Iterate over the parser, which yields GenBankRecord objects
    #     for i, record in enumerate(parser, 1):
    #         total_records += 1
    #         print(f"\n--- Record {i}: {record!r} ---")

    #         # Print core attributes
    #         print(f"  ID: {record.id}")
    #         print(f"  Name: {record.name}")

    #         # Demonstrate access to parsed entries
    #         if "LOCUS" in record.entries:
    #             locus = record.entries["LOCUS"]
    #             print(f"  Locus Name: {locus.name}")
    #             print(f"  Sequence Length (from LOCUS): {locus.length} bp")

    #         if "ORIGIN" in record.entries:
    #             seq = record.entries["ORIGIN"].sequence
    #             print(f"  Sequence Length (Actual): {len(seq)} bp")
    #             print(f"  Sequence Preview: {seq[:30]}...")

    #         if "FEATURES" in record.entries:
    #             features = record.entries["FEATURES"]
    #             print(f"  Number of features: {len(features.entries)}")
    #             if features.entries:
    #                 print(f"  First feature: {features.entries[0]}")

    #     print(f"\nSuccessfully parsed {total_records} records.")

    # except requests.exceptions.RequestException as e:
    #     print(f"Error downloading file: {e}")
    # except Exception as e:
    #     print(f"An unexpected error occurred during parsing: {e}")
    # finally:
    #     # 3. Clean up the downloaded file
    #     if os.path.exists(file_path):
    #         os.remove(file_path)
    #         print(f"\nCleaned up temporary file '{file_path}'.")


class Locus:
    _MOLECULE_TYPE_LIST: list[str] = ["DNA", "RNA", "PROTEIN"]

    def __init__(self, line: str) -> None:
        self._raw_line: str = line.strip()
        self._parts: list[str] = self._raw_line.split()
        self.name: str = ""
        self.length: int = 0
        self.molecule_type: str = ""
        self.topology: str = ""
        self.date: str = ""
        self._set_info()

    def _set_info(self):
        if not self._parts:
            return

        # LOCUS name is always second
        self.name: str = self._parts[1] if len(self._parts) > 1 else ""

        # Try to find the number as the length
        for i, token in enumerate(self._parts):
            if token.isdigit():
                self.length: int = int(token)
                break
        # Identify molcule type
        for token in self._parts:
            if token.upper() in self._MOLECULE_TYPE_LIST:
                self.molecule_type = token.upper()
                break
        # Detect topology
        for token in self._parts:
            if token.lower() in ("linear", "circular"):
                self.topology = token.lower()

        self.date: str = self._parts[-1] if len(self._parts) > 2 else ""

    def __repr__(self) -> str:
        return (
            f"Locus(name='{self.name}', length={self.length}, "
            f"molecule_type='{self.molecule_type}', topology='{self.topology}', "
            f"date='{self.date}')"
        )


class Definition:
    def __init__(self, info: str) -> None:
        self.info: str = info.replace("DEFINITION", "").strip()

    def __repr__(self) -> str:
        preview: str = f"{self.info[:70]}..." if len(self.info) > 70 else self.info
        return f"<Definition info='{preview}'>"


class Accession:
    def __init__(self, info: str) -> None:
        self.info: str = " ".join(info.split()[1:])

    def __repr__(self) -> str:
        return f"<Accession: {self.info}>"


class Version:
    def __init__(self, info: str) -> None:
        parts: list[str] = info.split()
        self.version: str = parts[1]
        self.gi: str | None = None

        for part in parts:
            if part.startswith("GI:"):
                self.gi = part.split(":")[1]
                break

    def __repr__(self) -> str:
        return f"<Version: {self.version} GI:{self.gi}>"


class Origin:
    """Representation of the ORIGIN entry in a GenBank format"""

    def __init__(self, raw_text: str) -> None:
        self._raw_text: str = raw_text

    @property
    def sequence(self) -> str:
        lines: list[str] = self._raw_text.splitlines()
        seq_lines: list[str] = []

        for line in lines:
            if line.startswith("ORIGIN") or line.startswith("//"):
                continue
            clean: str = "".join(ch for ch in line if ch.isalpha())
            seq_lines.append(clean)
        return "".join(seq_lines).lower()

    def __repr__(self) -> str:
        seq: str = self.sequence
        preview: str = seq[:20] + "..." if len(seq) > 20 else seq
        return f"<Origin sequence='{preview}' length={len(seq)}>"


class SingleFeature:
    """Representation of one of FEATURES entries, like a 'gene' or a 'CDS' block."""

    def __init__(self, key: str, location: str, qualifiers: dict[str, str]) -> None:
        self.key = key
        self.location = location
        self.qualifiers = qualifiers

    def __repr__(self) -> str:
        return f"Feature(key={self.key}, location={self.location}, qualifiers={len(self.qualifiers)})"


class Features:
    """Parses and stores the FEATURES block of a GenBank file"""

    def __init__(self, info: str) -> None:
        self.info: str = info
        self.entries: list[SingleFeature] = []
        self._parse_features()

    def __repr__(self) -> str:
        return f"<Features containing {len(self.entries)} entries>"

    def _parse_features(self) -> None:
        """Extract each feature line block (e.g., 'gene', 'CDS') and its qualifiers.

        Expected structure in GenBank:
            FEATURE_KEY    LOCATION
                            /qualifier1="..."
                            /qualifier2="..."
        """
        # self.info have all the block of FEATURES
        current_key: str = ""
        current_loc: str = ""
        current_quals: dict[str, str] = {}

        for line in self.info.splitlines()[1:]:
            # A new feature appears when a non-space appears in the first 21 columns
            # GenBank Uses fixed-width indentation
            if line[:21].strip() and not line.strip().startswith("/"):
                # Store the previous feature
                if current_key:
                    self.entries.append(
                        SingleFeature(current_key, current_loc, current_quals)
                    )
                parts = line.strip().split(None, 1)  # split feature and location
                current_key = parts[0]
                current_loc = parts[1] if len(parts) > 1 else ""
                current_quals = {}  # reset qualifiers for the new feature

            # This line is a qualifier for the new feature
            elif line.strip().startswith("/"):
                qualifier = line.strip()[1:]  # to remove leading slash
                # usually qualifier and values are seperated by "="
                if "=" in qualifier:
                    key, value = qualifier.split("=", 1)
                    value = value.strip('""')  # also remove quotes
                    current_quals[key] = value
                else:
                    # if qualifier doesn't have "="
                    current_quals[qualifier] = ""
            else:
                # This means that this line isn't a new feature nor a new qualifier
                # i.e it is a continuation of a previous qualifier
                # so add the current line to the last qualifier
                if current_quals:
                    last_key = list(current_quals.keys())[-1]
                    current_quals[last_key] += " " + line.strip()
        # append the last feature
        if current_key:
            self.entries.append(SingleFeature(current_key, current_loc, current_quals))


class GenBankRecord:
    """Represents a parsed GenBank record with entries"""

    _entry_classes: dict[str, type] = {
        "LOCUS": Locus,
        "DEFINITION": Definition,
        "ACCESSION": Accession,
        "FEATURES": Features,
        "ORIGIN": Origin,
        "VERSION": Version,
    }

    def __init__(
        self, entries: dict[str, Any], source_filepath: Path | None = None
    ) -> None:
        # The parser now provides the 'entries' dict directly.
        self.entries: dict[str, Any] = entries
        self._source_filepath = source_filepath

        # Access common fields, handling potential missing keys gracefully
        self.name: str = (
            self.entries.get("DEFINITION").info
            if "DEFINITION" in self.entries
            else "N/A"
        )
        self.seq: str = (
            self.entries.get("ORIGIN").sequence if "ORIGIN" in self.entries else ""
        )
        self.id: str = (
            self.entries.get("ACCESSION").info if "ACCESSION" in self.entries else "N/A"
        )

    def __repr__(self) -> str:
        locus_name = self.entries["LOCUS"].name if "LOCUS" in self.entries else "N/A"
        return f"<GenBankRecord for '{locus_name}'>"


class GenBankParser:
    """A parser that reads a GenBANK file and splits it into entry blocks"""

    _ENTRY_PATTERN = re.compile(r"^[A-Z]+(\s|$)")  # Line starts with all caps
    _RECORD_SEPARATOR = "//"

    def __init__(self, filepath: str | Path) -> None:
        self.filepath = Path(filepath)

    def read_all(self) -> str:
        with open(self.filepath) as f:
            return f.read()

    def _split_into_blocks(
        self, file_contents: str
    ) -> Generator[tuple[str, str], None, None]:
        """Yeild (key, block) pairs from a GenBank record lazily while preserving duplicates"""

        current_key: str | None = None
        buffer: list[str] = []

        for line in file_contents.splitlines():
            # Detect header line
            if self._ENTRY_PATTERN.match(line):
                # restore the block if there is already a present key
                if current_key:
                    yield current_key, "\n".join(buffer).strip()
                    buffer.clear()
                current_key = line.split()[0]
                buffer.append(line)
            else:
                # continuation line — belongs to current block
                buffer.append(line)
        # Save the last block
        if current_key:
            yield current_key, "\n".join(buffer).strip()

    def _split_into_records(self, file_contents: str) -> Generator[str, None, None]:
        """Splits the file content into multiple GenBank record blocks"""
        for record_block in file_contents.split(f"\n{self._RECORD_SEPARATOR}\n"):
            clean_block = record_block.strip()
            if clean_block:
                # Re-add the seprator if it was present
                if not clean_block.endswith(self._RECORD_SEPARATOR):
                    clean_block += f"\n{self._RECORD_SEPARATOR}"
                yield clean_block

    def _parse_record(self, record_text: str) -> GenBankRecord:
        """Parse a single record text block and returns a GenBank record object."""
        entries: dict[str, Any] = {}
        entry_classes = GenBankRecord._entry_classes

        for key, block in self._split_into_blocks(record_text):
            # if the entry maps to an entry in the _entry_classes else create new one
            entry_cls = entry_classes.get(key)
            if key in entries:
                if isinstance(entries[key], str) and isinstance(block, str):
                    # append repeated entries
                    entries[key] += "\n" + block
                else:
                    entries[key] = block
            else:
                # instantiate using the contructors in the _entry_classes
                # entry_cls here is a class (I have to remind myself that it is)
                # keep raw text if there is no mapping entry in the _entry_classes
                entries[key] = entry_cls(block) if entry_cls else block
        return GenBankRecord(entries, source_filepath=self.filepath)

    def __iter__(self) -> Iterator[GenBankRecord]:
        "Allows iterating over the records in the file"
        file_contents = self.read_all()
        for record_text in self._split_into_records(file_contents):
            # Parse each record block and yield the GenBankRecord object
            if (
                record_text.strip() != self._RECORD_SEPARATOR
            ):  # Skip empty blocks from splitting
                yield self._parse_record(record_text)


if __name__ == "__main__":
    main()
