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
        "ACP": ["PKS_PP", "PP-binding", "AcpP", "AcpS", "acpS", "ACPS"],
        "SAT": ["SAT"],
        "C": ["Condensation"],
        "A": ["A_NRPS", "AMP-binding"],
        "E": ["NRPS-para261"],
        "cAT": ["Carn_acyltransf"],
        "Ox": ["mcbC-like_oxidoreductase"],  # NRPS oxydation domain
        "KAT": ["RimI"],  # N-epsilon-Lysine acetyltransferase
    }

    pssm_lengths = {
        "PKS_KS": 298,
        "PKS": 421,
        "CLF": 399,
        "KAS_I_II": 406,
        "CHS_like": 361,
        "KAS_III": 320,
        "PKS_AT": 298,
        "Acyl_transf_1": 319,
        "PKS_ER": 287,
        "enoyl_red": 293,
        "KR": 180,
        "PKS_KR": 287,
        "Thioesterase": 224,
        "Aes": 312,
        "Thioester-redct": 367,
        "SDR_e1": 290,
        "Methyltransf_11": 95,
        "Methyltransf_12": 98,
        "Methyltransf_23": 162,
        "Methyltransf_25": 97,
        "Methyltransf_31": 150,
        "AdoMet_MTases": 107,
        "PKS_DH": 167,
        "PS-DH": 289,
        "PT_fungal_PKS": 324,
        "PKS_PP": 86,
        "PP-binding": 67,
        "AcpP": 80,
        "AcpS": 127,
        "acpS": 125,
        "ACPS": 104,
        "SAT": 239,
        "Condensation": 455,
        "A_NRPS": 444,
        "AMP-binding": 361,
        "NRPS-para261": 153,
        "Carn_acyltransf": 575,
        "mcbC-like_oxidoreductase": 180,
        "RimI": 177,
    }

    bitscore_cutoffs = {
        "PKS_KS": 241.079,
        "PKS": 167.35,
        "CLF": 503.431,
        "KAS_I_II": 285.201,
        "CHS_like": 242.515,
        "KAS_III": 212.4,
        "PKS_AT": 201.477,
        "Acyl_transf_1": 342.53,
        "PKS_ER": 250.768,
        "enoyl_red": 129.997,
        "KR": 149.632,
        "PKS_KR": 83.3005,
        "Thioesterase": 157.139,
        "Aes": 59.5636,
        "Thioester-redct": 305.493,
        "SDR_e1": 229.845,
        "Methyltransf_11": 53.8265,
        "Methyltransf_12": 40.8165,
        "Methyltransf_23": 61.2786,
        "Methyltransf_25": 34.8626,
        "Methyltransf_31": 67.449,
        "AdoMet_MTases": 29.3203,
        "PKS_DH": 76.8814,
        "PS-DH": 128.631,
        "PT_fungal_PKS": 202.466,
        "PKS_PP": 33.3777,
        "PP-binding": 28.2759,
        "AcpP": 26.8736,
        "AcpS": 88.0777,
        "acpS": 87.108,
        "ACPS": 31.8159,
        "SAT": 129.59,
        "Condensation": 345.862,
        "A_NRPS": 356.838,
        "AMP-binding": 185.114,
        "NRPS-para261": 159.361,
        "Carn_acyltransf": 255.925,
        "mcbC-like_oxidoreductase": 62.4272,
        "RimI": 43.8319,
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
        coverage_cutoff : int
            Minimum percent coverage to save a conserved domain hit. This is multiplied
            by the PSSM length defined by NCBI and then checked against the actual
            length of the given hit.
        bitscore_cutoff : int
            Percent multiplier for determining bitscore cutoff. The PSSM threshold
            bitscore, as calculated by NCBI, is multiplied by this value and is then
            checked against the bitscore of the hit.

        Returns
        -------
        Domain

        Raises
        ------
        ValueError
            If the domain in this row is not in the DOMAINS dictionary.
        """
        _, _, _, start, end, evalue, bitscore, _, domain, *_ = row.split("\t")
        for domain_type, domains in self.domains.items():
            if domain not in domains:
                continue
            return Domain(
                type=domain_type,
                domain=domain,
                start=int(start),
                end=int(end),
                evalue=float(evalue),
                bitscore=float(bitscore),
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
        domains = self.filter_domains(domains)

        if not domains:
            log.error("No domains remain after filtering for %s", name)
            return

        synthase = Synthase(
            header=name, sequence=sequence, domains=self.filter_domains(domains)
        )
        try:
            classify.classify(synthase)
        except ValueError:
            print(f"Failed to classify {synthase}")
        return synthase

    def detect_fragmented_domain(self, one, two, coverage_pct=0.5, tolerance_pct=0.1):
        """Detect if two adjacent domains are likely a single domain.

        This is useful in cases where a domain is detected with multiple small hits. For
        example, an NRPS may have two adjacent condensation (C) domain hits that are
        both individually too small and low-scoring, but should likely just be merged.

        If two hits are close enough together, such that the distance between the start
        of the first and end of the second is within some tolerance (default +-10%) of the
        total length of a domains PSSM, this function will return True.
        """
        pssm_length = self.pssm_lengths[one.domain]
        coverage = pssm_length * coverage_pct
        tolerance = pssm_length * tolerance_pct
        one_length, two_length = len(one), len(two)

        return (
            one_length < coverage
            and two_length < coverage
            and pssm_length - tolerance
            <= two.end - one.start
            <= pssm_length + tolerance
            and one_length + two_length > coverage
        )

    def filter_domains(self, domains):
        """Filter overlapping Domain objects and test special rules.

        Special rules are tested again here, in case they are missed within overlap
        groups. For example, the NRPS-para261 domain is not always entirely contained by
        a condensation domain, so should be caught by this pass.
        """
        filtered = [
            self.filter_domain_group(group) for group in group_overlapping_hits(domains)
        ]

        i, total = 1, len(filtered)
        while i < total:
            if i + 1 == total:
                break
            previous, current = filtered[i - 1 : i + 1]

            # When domains are likely together, e.g. two small C domain hits right next
            # to each other
            if previous.domain == current.domain and self.detect_fragmented_domain(
                previous, current
            ):
                previous.end = current.end
                del filtered[i]
                continue

            for new_type, rule in self.special_rules.items():
                if rule(previous, current):
                    previous.type = new_type
                    del filtered[i]
                    break
            i += 1

        # Final filter, get rid of any obviously wrong, small hits; do here, not when
        # parsing table, so we don't discard potentially fragmented single domains
        return [
            domain
            for domain in filtered
            if len(domain) > 0.2 * self.pssm_lengths[domain.domain]
        ]

    def filter_domain_group(self, group):
        """Select the best domain from a collection of overlapping domains.

        This function tests rules stored in `special_rules`, which are lambdas that
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
        container, *_group = sorted(group, key=lambda d: d.evalue)
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


def hits_overlap(a, b, threshold=0.2):
    """Return True if Domain overlap is greater than threshold * domain size.

    The specified domains typically should share no overlap, tending to be discrete with
    inter-domain gaps. Subsequently, the default threshold is set low. However, it is
    not 0, as to accomodate SOME level of overlap.

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


def group_overlapping_hits(domains, threshold=0.2):
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
