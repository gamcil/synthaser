#!/usr/bin/env python3


def count(fasta):
    """Counts sequences in an open FASTA file handle.

    Iterates file and counts header lines. Then, seeks to start
    of the file and returns the count.

    Parameters:
        fasta (file pointer): An open file handle corresponding to a FASTA file.
    Returns:
        count (int): Total number of sequences in the file.
    """
    count = 0
    for line in fasta:
        if line.startswith(">"):
            count += 1
    fasta.seek(0)
    return count


def parse(fasta):
    """Parses an open FASTA file for sequences.

    For example, given a FASTA file `fasta.faa` containing:

    ::

        >sequence
        ACGTACGTACGT

    This file can be parsed:

    >>> with open('fasta.faa') as handle:
    ...     parse_fasta(handle)
    {"sequence": "ACGTACGTACGT"}

    Parameters:
        fasta (file pointer): Open FASTA file handle
    Returns:
        sequences (dict): Sequences in the FASTA file, keyed on sequence headers.
    """
    sequences = {}
    for line in fasta:
        try:
            line = line.decode().strip()
        except AttributeError:
            line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            sequences[header] = ""
        else:
            sequences[header] += line
    return sequences


def wrap(sequence, limit=80):
    """Wraps sequences to `limit` characters per line.

    Parameters:
        sequence (str): Sequence to be wrapped.
        limit (int): Total characters per line.
    Returns:
        (str): Sequence wrapped to maximum `limit` characters per line.
    """
    return "\n".join(sequence[i : i + limit] for i in range(0, len(sequence), limit))


def create(header, sequence, limit=80):
    """Creates a FASTA format string from a header and sequence.

    For example:

    >>> fasta = create_fasta('header', 'AAAAABBBBBCCCCC', wrap=5)
    >>> print(fasta)
    >header
    AAAAA
    BBBBB
    CCCCC

    Parameters:
        header (str): Name to use in FASTA definition line (i.e. >header).
        sequence (str): The sequence corresponding to the `header`.
        limit (int): Total characters per line before sequence is wrapped.
    Returns:
        (str): FASTA format string.
    """
    return ">{}\n{}".format(header, wrap(sequence, limit=limit))
