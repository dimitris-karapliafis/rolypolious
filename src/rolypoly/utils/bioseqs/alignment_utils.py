"""Alignment and mapping utility functions."""

import re


def calculate_percent_identity(cigar_string: str, num_mismatches: int) -> float:
    """Calculate sequence identity percentage from CIGAR string and edit distance.

    Computes the percentage identity between aligned sequences using the CIGAR
    string from an alignment and the number of mismatches (NM tag).

    Args:
        cigar_string (str): CIGAR string from sequence alignment
        num_mismatches (int): Number of mismatches (edit distance)

    Returns:
        float: Percentage identity between sequences (0-100)

    Note:
        The calculation considers matches (M), insertions (I), deletions (D),
        and exact matches (=) from the CIGAR string.

    Example:
         print(calculate_percent_identity("100M", 0))
         100.0
         print(calculate_percent_identity("100M", 2))
         98.0
    """

    cigar_tuples = re.findall(r"(\d+)([MIDNSHPX=])", cigar_string)
    matches = sum(int(length) for length, op in cigar_tuples if op in {"M", "=", "X"})
    total_length = sum(
        int(length) for length, op in cigar_tuples if op in {"M", "I", "D", "=", "X"}
    )
    return (matches - num_mismatches) / total_length * 100 