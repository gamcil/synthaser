#!/usr/bin/env python3


import logging
import re

from synthaser.models import Domain, Synthase, hits_overlap
from synthaser.classify import classify_synthase


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class ResultParser:
    """The ResultParser class handles the parsing and filtering of CDSearch results.

    First, instantiate a `ResultParser` object:

    >>> rp = ResultParser()

    Then, pass an open results file handle to `ResultParser.parse_results()`:

    >>> with open('results.tsv') as handle:
    ...     figure = rp.parse_results(results)

    A Figure object is returned which contains fully instantiated Synthase and Domain
    models:

    >>> figure
    HR-PKS
    ------
    KS-AT-DH-ER-KR-ACP-TE
    ...

    These can then be used to directly instantiate Domain and Synthase objects using
    their `from_dict()` methods.

    Parameters
    ----------
    domains : dict
        Dictionary mapping actual domain names from the conserved domain database (CDD)
        to their broader biosynthetic type. This can be altered through the method
        `set_domains` to add more categories or domains.
    """

    default_domains = {
        "KS": ["PKS_KS", "PKS", "CLF", "KAS_I_II", "CHS_like", "KAS_III"],
        "AT": ["PKS_AT", "Acyl_transf_1"],
        "ER": ["PKS_ER", "enoyl_red"],
        "KR": ["KR", "PKS_KR"],
        "TE": ["Thioesterase", "Aes"],
        "TR": ["Thioester-redct", "SDR_e1"],
        "MT": [
            "Methyltransf_11",
            "Methyltransf_12",
            "Methyltransf_23",
            "Methyltransf_25",
            "Methyltransf_31",
            "AdoMet_MTases",
        ],
        "DH": ["PKS_DH", "PS-DH"],
        "PT": ["PT_fungal_PKS"],
        "ACP": ["PKS_PP", "PP-binding", "AcpP"],
        "SAT": ["SAT"],
        "C": ["Condensation"],
        "A": ["A_NRPS", "AMP-binding"],
        # "E": ["NRPS-para261"],
    }

    def __init__(self, domains=None):
        self.domains = self.default_domains.copy()
        if domains:
            self.set_domains(domains)

    def set_domains(self, domains):
        """Update the internal domain type dictionary with new domains.

        Parameters
        ----------
        domains : dict
            Dictionary of conserved domain names keyed on domain type. Valid keys are
            those in `ResultParser.domains`; if a provided key is not in this set, a
            ValueError is thrown. Values should be tuples or lists of str corresponding
            to names of conserved domains in the CDD.

        Raises
        ------
        TypeError
            If `domains` is not of type `dict`.
        ValueError
            If a domain type is specified that is not in `ResultParser.domains`.
        TypeError
            If dictionary values are not lists or tuples.
        TypeError
            If a given value list element is not a string.
        """
        if not isinstance(domains, dict):
            raise TypeError("Expected dict")
        for key, value in domains.items():
            if key not in self.domains:
                raise ValueError(f"Invalid domain '{key}'")
            if not isinstance(value, (tuple, list)):
                raise TypeError("Expected domains in list or tuple")
            if any(not isinstance(v, str) for v in value):
                raise TypeError("Expected list/tuple of str")
            self.domains[key].extend(d for d in value if d not in self.domains[key])

    def parse_row(self, row):
        """Parse a domain hit from a row in a CD-search results file.

        For example, a typical row might looks like:

        >>> print(row)
        Q#1 - >AN6791.2\tspecific\t225858\t9\t1134\t0\t696.51\tCOG3321\tPksD\t-\tcl09938

        Using this function will generate:

        >>> rp.parse_row(row)
        PksD [KS] 9-1134

        Parameters
        ----------
        row : str
            Tab-separated row from a CDSearch results file.

        Returns
        -------
        Domain

        Raises
        ------
        ValueError
            If the domain in this row is not in the DOMAINS dictionary.
        """
        _, _, _, start, end, _, _, _, domain, *_ = row.split("\t")
        for domain_type, domains in self.domains.items():
            if domain not in domains:
                continue
            return Domain(
                type=domain_type, domain=domain, start=int(start), end=int(end)
            )
        raise ValueError(f"'{domain}' not a synthaser key domain")

    def parse(self, results_handle, query_handle=None, sequences=None):
        """Parse CD-Search results.

        Parameters
        ----------
        results_handle: open file handle
            An open CD-Search results file handle. If you used the website to analyse your
            sequences, the file you should download is Domain hits, Data mode: Full, ASN
            text. When using a `CDSearch` object, this format is automatically selected.

        query_handle: open file handle, optional
            An open file handle for the sequences used in the CD-search run. Sequences
            can be added later via `Figure.add_query_sequences()`.

        Returns
        -------
        list:
            A list of Synthase objects parsed from the results file.
        """
        query_regex = re.compile(r"Q#\d+? - [>]?(.+?)\t")

        _query = ""
        _synthase = None
        synthases = []

        if not sequences:
            if query_handle:
                sequences = parse_fasta(query_handle) if query_handle else {}
            else:
                sequences = {}

        for row in results_handle:
            try:
                row = row.decode()
            except AttributeError:
                pass  # in case rows are unicode

            if not row.startswith("Q#") or row.isspace():
                continue

            query = query_regex.search(row).group(1)

            if query != _query:
                if _synthase:
                    synthases.append(_synthase)
                _query = query
                _synthase = Synthase(
                    header=query,
                    sequence=sequences[query] if query in sequences else "",
                )
            try:
                domain = self.parse_row(row)
            except ValueError:
                continue

            _synthase.domains.append(domain)

        if _synthase:
            synthases.append(_synthase)

        for synthase in synthases:
            synthase.filter_overlapping_domains()
            try:
                classify_synthase(synthase)
            except ValueError:
                log.warning("Failed to classify %s", synthase.header)
            synthase.rename_nrps_domains()

        return synthases


def parse_fasta(fasta):
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


def parse_results(results_file):
    """Parse a results file at a specified path.

    Instantiates a ResultParser object and returns the output of its `parse()` method.
    """
    rp = ResultParser()
    with open(results_file) as results:
        return rp.parse(results)
