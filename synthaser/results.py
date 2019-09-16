#!/usr/bin/env python3


import logging
import re
from collections import defaultdict
from operator import attrgetter

from synthaser import fasta, classify
from synthaser.models import Domain, Synthase


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

    special_rules : dict
        Special classification rules that will be tested in `filter_domains`. These
        should be lambda functions that take two variables: the container Domain (best
        scoring of an overlap group) and other Domains in the same group. If True, the
        container Domain type is set to the type the rule is keyed on.
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
        "E": ["NRPS-para261"],
    }

    special_rules = {
        "E": lambda a, b: len(a) > len(b) and a.type == "C" and b.type == "E"
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
        _, _, _, start, end, evalue, _, _, domain, *_ = row.split("\t")
        for domain_type, domains in self.domains.items():
            if domain not in domains:
                continue
            return Domain(
                type=domain_type,
                domain=domain,
                start=int(start),
                end=int(end),
                evalue=float(evalue),
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

        sequences : dict
            Sequences for query sequences, keyed on query name.

        Returns
        -------
        list:
            A list of Synthase objects parsed from the results file.
        """
        results = self.parse_table(results_handle)
        return self.build_synthases(
            results, query_handle=query_handle, sequences=sequences
        )

    def parse_table(self, results_handle):
        """Parse a CD-Search results table and instantiate Domain objects for each hit.

        Parameters
        ----------
        results_handle : open file handle
            Open file handle corresponding to a CD-Search results file.

        Returns
        -------
        results : dict
            Lists of Domain objects keyed on the query they were found in.
        """
        query_regex = re.compile(r"Q#\d+? - [>]?(.+?)\t")
        results = defaultdict(list)
        for row in results_handle:
            try:
                row = row.decode()
            except AttributeError:
                pass  # in case rows are unicode
            if not row.startswith("Q#") or row.isspace():
                continue
            query = query_regex.search(row).group(1)
            try:
                domain = self.parse_row(row)
            except ValueError:
                continue
            results[query].append(domain)
        return dict(results)

    def build_synthases(self, results, query_handle=None, sequences=None):
        """Build Synthase objects from a parsed results dictionary.

        Parameters
        ----------
        results : dict
            Grouped Domains; output from `parse_table`.
        query_handle: open file handle, optional
            Open file handle corresponding to the FASTA file used as a query. Parsed for
            sequences.
        sequences: dict, optional
            Sequences for each query sequence keyed on query name.

        Returns
        -------
        list
            Synthase objects containing all Domain objects found in the CD-Search.
        """
        if not sequences:
            if query_handle:
                sequences = fasta.parse(query_handle) if query_handle else {}
            else:
                sequences = {}
        return [
            self.new_synthase(
                name, domains, sequences[name] if name in sequences else ""
            )
            for name, domains in results.items()
        ]

    def new_synthase(self, name, domains, sequence=""):
        """Instantiate and classify a new Synthase object."""
        synthase = Synthase(
            header=name, sequence=sequence, domains=self.filter_domains(domains)
        )
        classify.classify(synthase)
        return synthase

    def filter_domains(self, domains):
        """Filter overlapping Domain objects and test special rules."""
        return [
            self.apply_special_rules(group) for group in group_overlapping_hits(domains)
        ]

    def apply_special_rules(self, group):
        """Test an overlapping Domain group for special rules.

        This function uses the rules stored in `special_rules`, which are lambdas that
        take two variables. It sorts the group by e-value, then tests each rule using
        the container (first, best scoring group) against all other Domains in the
        group.

        If any test is True, the container type is set to the rule key and returned.
        Otherwise, this function will return the container Domain with no modification.

        Parameters
        ----------
        group : list
            Overlapping `Domain` objects, the yielded product from `group_overlapping_hits`.

        Returns
        -------
        Domain
            Highest scoring `Domain` in the group. If any special rules have been satisfied,
            the type of this `Domain` will be set to that rule (e.g. Condensation ->
            Epimerization).
        """
        container, *_group = sorted(group, key=attrgetter("evalue"))
        for domain in _group:
            for new_type, rule in self.special_rules.items():
                if rule(container, domain):
                    container.type = new_type
                    return container
        return container


def parse_results(results_file):
    """Parse a results file at a specified path.

    Instantiates a ResultParser object and returns the output of its `parse()` method.
    """
    rp = ResultParser()
    with open(results_file) as results:
        return rp.parse(results)


def hits_overlap(a, b, threshold=0.9):
    """Return True if Domain overlap is greater than threshold * domain size.

    Parameters
    ----------
    a : Domain
        First Domain object.
    b : Domain
        Second Domain object.
    threshold : int
        Minimum percentage to classify two Domains as overlapping. By default,
        `threshold` is set to 0.9, i.e. two Domains are considered as overlapping if
        the total amount of overlap is greater than 90% of either Domain hit.

    Returns
    -------
    bool
        True if Domain overlap exceeds threshold, False if not.
    """
    start, end = max(a.start, b.start), min(a.end, b.end)
    overlap = max(0, end - start)
    a_threshold = threshold * (a.end - a.start)
    b_threshold = threshold * (b.end - b.start)
    return overlap >= a_threshold or overlap >= b_threshold


def group_overlapping_hits(domains, threshold=0.9):
    """Iterator that groups Domain objects based on overlapping locations.

    Parameters
    ----------
    domains : list, tuple
        Domain objects to be grouped.
    threshold : float
        See hits_overlap().

    Yields
    ------
    group : list
        A group of overlapping Domain objects, as computed by hits_overlap().
    """
    domains.sort(key=attrgetter("start"))
    i, total = 0, len(domains)
    while i < total:
        current = domains[i]  # grab current hit
        group = [current]  # start group
        if i == total - 1:  # if current hit is the last, yield
            yield group
            break
        for j in range(i + 1, total):  # iterate rest
            future = domains[j]  # grab next hit
            if hits_overlap(current, future, threshold):
                group.append(future)  # add if contained
            else:
                yield group  # else yield to iterator
                break
            if j == total - 1:  # if reached the end, yield
                yield group
        i += len(group)  # move index ahead of last group
