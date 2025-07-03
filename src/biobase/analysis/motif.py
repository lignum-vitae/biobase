# external dependency
import re

# internal dependency
from biobase.constants.amino_acid import ONE_LETTER_CODES, ONE_LETTER_CODES_EXT


def main():
    print(find_motifs("ACDEFGHIKLMNPQRSTVWY", "CDE"))
    print(find_motifs("GGGGGGGGGGGGGGGGGGGG", "CDE"))
    test_dict = {
        ">SP001": "ACDEFCDEFCDEFGHIKLMN",  # has matches for "CDE" at positions 2, 6, 10
        ">SP002": "MNPQRSTVWYACDEFGHIKL",  # has match for "CDE" at position 12
        ">SP003": "AAAAAAAAAAAAAAAAAA12",  # invalid: contains "1", "2"
        ">SP004": "GGGGGGGGGGGGGGGGGGGG",  # no match
        ">SP005": "HHHHHHHHHHHHHHHHH@#$",  # invalid: contains "@", "#", "$"
        ">SP006": "DDDDDDDDDDDDDDDDDDDD",  # no match
        ">SP007": "CDEFGHCDEFKLCDEFPQRS",  # has matches for "CDE" at positions 1, 7, 13
        ">SP008": "LLLLLLLLLLLLLLLLLLLL",  # no match
        ">SP009": "KKKKKKKKKKKK123KKKKK",  # invalid: contains "1", "2", "3"
        ">SP010": "CDEACDEBCDEFAAAAAAAA",  # has matches for "CDE" at positions 1, 5, 9
    }
    print(find_motifs(test_dict, "CDE"))


def find_motifs(
    sequence: str | dict[str, str], pattern: str, ext: bool = False
) -> list[int] | tuple[dict[str, list[int]], dict[str, set[str]], list[str]] | None:
    """
    Find all occurrences of a specified motif (pattern) in protein sequence(s).

    This function uses regular expressions to locate all positions of a motif (including overlapping matches)
    within protein sequence(s). It can handle either a single protein sequence or a dictionary of sequences
    in FASTA format.

    Parameters:
    - sequence: Either:
        - str: A protein sequence where each character is an amino acid's one-letter code
        - dict[str, str]: A FASTA dictionary where keys are protein identifiers and values are sequences
    - pattern (str): A string representing the motif to search for. Can be plain text for exact matches
                    or a Python-flavoured regular expression string for complex patterns.
    - ext (bool): If True, uses extended amino acid codes. Defaults to False.

    Returns:
    - For single sequence (str input):
        list[int]: A list of 1-based start positions where the motif is found
    - For FASTA dictionary (dict input):
        tuple containing:
        - dict[str, list[int]]: Dictionary mapping sequence IDs to lists of motif positions
        - dict[str, set[str]]: Dictionary mapping sequence IDs to sets of invalid characters found
        - list[str]: List of sequence IDs that had no matches (but were valid sequences)

    Raises:
    - ValueError: If sequence(s) contain invalid characters
    - ValueError: If any input is empty
    - ValueError: If pattern is not a string

    Examples:
    >>> find_motifs("ACDEFGHIKLMNPQRSTVWY", "CDE")
    [2]
    >>> find_motifs({"P12345": "ACDEFGHIKLMNPQRSTVWY"}, "CDE")
    {"P12345": [2]}
    """
    if not isinstance(pattern, str) or not pattern:
        raise ValueError("The pattern must be a non-empty string.")

    aa_codes = set(ONE_LETTER_CODES_EXT if ext else ONE_LETTER_CODES)

    # Handle single sequence string
    if isinstance(sequence, str):
        if not sequence:
            raise ValueError("The input protein sequence is empty.")
        invalid_chars = set(sequence) - aa_codes
        if invalid_chars:
            raise ValueError(
                f"Invalid protein sequence used in motif finder. Invalid characters: {sorted(invalid_chars)}"
            )
        return [m.start(0) + 1 for m in re.finditer(f"(?={pattern})", sequence)]

    # Handle FASTA dictionary
    if isinstance(sequence, dict):
        if (
            not sequence
            or any(x == "" for x in sequence.values())
            or any(not isinstance(x, str) for x in sequence.values())
        ):
            raise ValueError("The input FASTA dictionary is empty.")

        result_dict, invalid_ids = {}, {}
        non_matches = []
        for seq_id, seq in sequence.items():
            invalid_chars = set(seq) - aa_codes
            matches = [m.start(0) + 1 for m in re.finditer(f"(?={pattern})", seq)]
            if matches:  # only include sequences with matches
                result_dict[seq_id] = matches
            elif invalid_chars:
                invalid_ids[seq_id] = invalid_chars
            else:
                non_matches.append(seq_id)
        return result_dict, invalid_ids, non_matches
    raise ValueError("The input must be a non-empty string or FASTA dictionary.")


if __name__ == "__main__":
    main()
