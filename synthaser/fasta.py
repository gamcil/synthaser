#!/usr/bin/env python3

from Bio import SeqIO


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


def wrap(sequence, limit=80):
    """Wraps sequences to `limit` characters per line.

    Parameters:
        sequence (str): Sequence to be wrapped.
        limit (int): Total characters per line.
    Returns:
        (str): Sequence wrapped to maximum `limit` characters per line.
    """
    return "\n".join(sequence[i: i + limit] for i in range(0, len(sequence), limit))


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


def parse(handle):
    return {
        record.name: str(record.seq)
        for record in SeqIO.parse(handle, 'fasta')
    }
