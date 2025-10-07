import os
import re
from pathlib import Path

# import requests


def main() -> None:
    """
    Main function to demonstrate the GenBank parser.
    It downloads a sample file, parse it, prints the contents, and clean it afterwards
    """

    # URL for a sample GenBank file (human p53 gene) from NCBI
#    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NG_017013.2&rettype=gb&retmode=text"
#    file_path = "temp_record.gb"
#    try:       # Download the file
#        print("Downloading sample file from NCBI...")
#        response = requests.get(url)
#        response.raise_for_status()  # Raise an exception for bad status codes
#        with open(file_path, "w") as f:
#            f.write(response.text)
#        print(f"File saved to '{file_path}'")
#
#        # Parse the file and print results
#        print("\n--- Parsing GenBank Record ---")
#        g1 = GenBankRecord(file_path)
#        print(f"Record Object: {g1}\n")
#
#        for key, info in g1.entries.items():
#            print(f"{key}\n{info}\n===============================================")
#
#        print("\n--- Example of accessing parsed data ---")
#        if "LOCUS" in g1.entries:
#            print(f"Locus Name: {g1.entries['LOCUS'].name}")
#        if "ORIGIN" in g1.entries:
#            print(f"Sequence Length: {len(g1.entries['ORIGIN'].sequence)}")
#        if "FEATURES" in g1.entries:
#            print(f"Number of features: {len(g1.entries['FEATURES'].entries)}")
#            if g1.entries["FEATURES"].entries:
#                print(f"First feature: {g1.entries['FEATURES'].entries[0]}")
#
#    except requests.exceptions.RequestException as e:
#        print(f"Error downloading file: {e}")
#    except Exception as e:
#        print(f"An error occurred: {e}")
#    finally:
#        # Clean up the downloaded file
#        if os.path.exists(file_path):
#            os.remove(file_path)
#            print(f"\nCleaned up '{file_path}'.")


class Locus:
    _MOLECULE_TYPE_LIST : list[str] = ["DNA", "RNA", "PROTEIN"]
    def __init__(self, line: str) -> None:
        self._raw_line : str = line.strip() 
        self._parts: list[str] = self._raw_line.split()
        self.name : str = ""
        self.length : int = 0
        self.molecule_type : str = ""
        self.topology : str = ""
        self.date : str = "" 
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
        # Identify molcule type - date is not treated as one
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


class Feature:
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
        self.entries: list[Feature] = []
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
                        Feature(current_key, current_loc, current_quals)
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
            self.entries.append(Feature(current_key, current_loc, current_quals))


class GenBankRecord:
    _ENTRY_PATTERN = re.compile(r"^[A-Z]+(\s|$)")  # Line starts with all caps

    _entry_classes = {
        "LOCUS": Locus,
        "DEFINITION": Definition,
        "ACCESSION": Accession,
        "FEATURES": Features,
        "ORIGIN": Origin,
        "VERSION": Version,
    }

    def __init__(self, filepath: str | Path) -> None:
        self._filepath = filepath
        self.entries = {}  # empty because will be filled with the function
        self._parse_entries()
        # TODO: faetures goes to the Features class

    def __repr__(self) -> str:
        locus_name = self.entries["LOCUS"].name if "LOCUS" in self.entries else "N/A"
        return f"<GenBankRecord for '{locus_name}' from '{self._filepath}'>"

    def _split_into_blocks(self, file_contents: str) -> list[tuple[str, str]]:
        """Split a GenBank reord by into (key, block) pairs while preserving duplicates"""
        blocks: list[tuple[str, str]] = []
        current_key: str | None = None
        buffer: list[str] = []

        for line in file_contents.splitlines():
            # Detect header line
            if self._ENTRY_PATTERN.match(line):
                # restore the block if there is already a present key
                if current_key:
                    blocks.append((current_key, "\n".join(buffer).strip()))
                    buffer.clear()
                current_key = line.split()[0]
                buffer.append(line)
            else:
                # continuation line — belongs to current block
                buffer.append(line)
        # Save the last block
        if current_key:
            blocks.append((current_key, "\n".join(buffer).strip()))

        return blocks

    def _parse_entries(self):  # returns a dict
        """Parse all GenBank entries and instantiate known classes."""
        with open(self._filepath) as f:
            text: str = f.read()

        for key, block in self._split_into_blocks(text):
            # if the entry maps to an entry in the _entry_classes else create new one
            cls = self._entry_classes.get(key)
            if key in self.entries:
                if isinstance(self.entries[key], str) and isinstance(block, str):
                    # append repeated entries
                    self.entries[key] += "\n" + block
                else:
                    self.entries[key] = block
            else:
                # instantiate using the contructors in the _entry_classes
                # cls here is a class (I have to remind myself that it is)
                # keep raw text if there is no mapping entry in the _entry_classes
                self.entries[key] = cls(block) if cls else block


if __name__ == "__main__":
    main()
