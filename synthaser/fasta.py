#!/usr/bin/env python3

"""
This module contains helper functions for handling FASTA format files/strings.
"""


def count(fasta):
    """Count sequences in an open FASTA file handle.

    Iterates each line, and counts header lines (start with `>`). Then, seeks to start
    of the file and returns the count.

    Parameters
    ----------
    fasta : open file handle
        An open file handle corresponding to a FASTA file.

    Returns
    -------
    count : int
        Total number of sequences in the file.
    """
    count = 0
    for line in fasta:
        if line.startswith(">"):
            count += 1
    fasta.seek(0)
    return count


def parse(fasta):
    """Parse an open FASTA file for sequences.

    For example, given a FASTA file `fasta.faa` containing:

    ::
        >sequence
        ACGTACGTACGT

    This file can be parsed:

    >>> with open('fasta.faa') as handle:
    ...     parse_fasta(handle)
    {"sequence": "ACGTACGTACGT"}

    Parameters
    ----------
    fasta : str
        Either an open file handle of a FASTA file, or newline split string (e.g. read
        in via readlines()) that can be iterated over.

    Returns
    -------
    sequences : dict
        Sequences in the FASTA file, keyed on sequence headers.
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
    """Wrap sequences to `limit` characters per line.

    Parameters
    ----------
    sequence : str
        Sequence to be wrapped.
    limit : int
        Total characters per line.

    Returns
    -------
    str
        Sequence wrapped to maximum `limit` characters per line.
    """
    return "\n".join(sequence[i : i + limit] for i in range(0, len(sequence), limit))


def create(header, sequence, limit=80):
    """Create a FASTA format string from a header and sequence.

    For example:

    >>> fasta = create_fasta('header', 'AAAAABBBBBCCCCC', wrap=5)
    >>> print(fasta)
    >header
    AAAAA
    BBBBB
    CCCCC

    Parameters
    ----------
    header : str
        Name to use in FASTA definition line (i.e. >header).
    sequence : str
        The sequence corresponding to the `header`.
    limit : int
        The number of characters per line for wrapping the given `sequence`.
        This function will call `wrap_fasta`.

    Returns
    -------
    str
        FASTA format string.
    """
    return ">{}\n{}".format(header, wrap(sequence, limit=limit))
