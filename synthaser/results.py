#!/usr/bin/env python3

"""
This module stores functions for parsing CD-Search output.

It provides a single public function, `parse`, which takes an open file handle
corresponding to a CD-Search hit table and returns a list of fully characterized
Synthase objects, i.e.:

>>> from synthaser import results
>>> with open('results.txt') as handle:
...     synthases = results.parse(handle)
>>> synthases
[AN6791.2 KS-AT-DH-MT-ER-KR-ACP, ... ]

Domain associations are specified in the `DOMAINS` dictionary, which maps entries in
NCBI's Conserved Domain Database (CDD) to broader megasynthase domain categories.
For example:

>>> results.DOMAINS['KS']
['PKS_KS', 'PKS', 'FabB', 'CLF', ... ]

The functions in this module are agnostic to the contents of this dictionary. Thus,
existing domain types can be extended with new CDD entries, or new domain types can be
created, just by adding to this dictionary.

When doing this, you must add corresponding entries to `CD_LENGTHS` and
`BITSCORE_THRESHOLDS`, which store the lengths of each CD, and the bitscore threshold
used by CD-Search to register the CD as a specific hit, respectively.

For example, to add the carnitine O-acyltransferase (cAT) domain, recently characterized
in an HR-PKS system (DOI: 10.1002/anie.201705237), we perform the following three steps:

>>> results.DOMAINS['cAT'] = ['Carn_acyltransf']
>>> results.CD_LENGTHS['Carn_acyltransf'] = 575
>>> results.BITSCORE_THRESHOLDS['Carn_acyltransf'] = 255.925

Future searches in this environment will then report any instances of this domain.
Currently, a fairly extensive list of fungal PKS domains is hard-coded into this
dictionary. Suggestions for more are welcomed via GitHub issues or pull requests.

Finally, rules can be created for filtering of adjacent domains in the `ADJACENCY_RULES`
dictionary. These rules take the form of lambda functions that accept two Domain objects
as input and return True/False. This is useful since the CDD does not always have exact
CDs that match certain synthase domain types. For example, NRPS epimerization (E)
domains are always detected as condensation (C) domains. However, every NRPS containing
an E domain (in my observations) also report hits for the NRPS-para261 domain family
either within (near the end), or immediately after, the C domain hit. Thus, I created an
adjacency rule to specifically check for this:

>>> results.ADJACENCY_RULES['E'] = lambda a, b: len(a) > len(b) and a.type == "C" and b.type == "E"

Here I define a lambda function which will take two Domain objects and return True if
1) they are C and E domains, and 2) the C is bigger than the E.
This is checked within both overlapping Domain groups (i.e. look for containment) and
against the entire list of Domain objects in a Synthase.
"""


import logging
import json
import re

from pathlib import Path

from collections import defaultdict
from operator import attrgetter

from synthaser.models import Domain


LOG = logging.getLogger(__name__)


DOMAINS = {}


def load_domain_json(json_file):
    """Load domain information from JSON.

    Expects a file with format:

    .. code-block::

        {
            "PKS_KS": {
                "length": 298,
                "bitscore": 241.079,
                "type": "KS"
            },
            ...
        }
    """
    with open(json_file) as fp:
        rules = json.load(fp)
        DOMAINS.clear()
        DOMAINS.update(rules)


# Load defaults, stored in synthaser/domains.json
load_domain_json(Path(__file__).parent / "domains.json")

ADJACENCY_RULES = {
    "E": lambda a, b: len(a) > len(b) and a.type == "C" and b.type == "E"
}


def _domain_from_row(row):
    """Parse a domain hit from a row in a CD-search results file.

    For example, a typical row might looks like:

    >>> print(row)
    Q#1 - >AN6791.2\tspecific\t225858\t9\t1134\t0\t696.51\tCOG3321\tPksD\t-\tcl09938

    Using this function will generate:

    >>> domain_from_row(row)
    PksD [KS] 9-1134

    Parameters
    ----------
    row : str
        Tab-separated row from a CDSearch results file

    Returns
    -------
    models.Domain
        Instance of the Domain class containing information about this hit

    Raises
    ------
    ValueError
        If the domain in this row is not in the DOMAINS dictionary.
    """
    _, _, _, start, end, evalue, bitscore, _, domain, *_ = row.split("\t")
    if domain not in DOMAINS:
        raise ValueError(f"'{domain}' not a synthaser key domain")
    return Domain(
        type=DOMAINS[domain]["type"],
        domain=domain,
        start=int(start),
        end=int(end),
        evalue=float(evalue),
        bitscore=float(bitscore),
    )


def _parse_cdsearch_table(handle):
    """Parse a CD-Search results table and instantiate Domain objects for each hit.

    Parameters
    ----------
    handle : open file handle
        Open file handle corresponding to a CD-Search results file.

    Returns
    -------
    results : dict
        Lists of Domain objects keyed on the query they were found in.
    """
    query_regex = re.compile(r"Q#\d+? - [>]?(.+?)\t")
    results = defaultdict(list)
    for row in handle:
        try:
            row = row.decode()
        except AttributeError:
            pass  # in case rows are unicode
        if not row.startswith("Q#") or row.isspace():
            continue
        query = query_regex.search(row).group(1)
        try:
            domain = _domain_from_row(row)
        except ValueError:
            continue
        results[query].append(domain)
    return dict(results)


def _filter_results(results, **kwargs):
    """Build Synthase objects from a parsed results dictionary.

    Any additional kwargs are passed to `_filter_domains`.

    Parameters
    ----------
    results : dict
        Grouped Domains; output from `_parse_cdsearch_table`.

    Returns
    -------
    synthases : list
        Synthase objects containing all Domain objects found in the CD-Search.
    """
    _results = {}
    for name, domains in results.items():
        _domains = _filter_domains(domains, **kwargs)
        if not _domains:
            LOG.error("No domains remain after filtering for %s", name)
        _results[name] = _domains
    return _results


def _is_fragmented_domain(one, two, coverage_pct=0.5, tolerance_pct=0.1):
    """Detect if two adjacent domains are likely a single domain.

    This is useful in cases where a domain is detected with multiple small hits. For
    example, an NRPS may have two adjacent condensation (C) domain hits that are
    both individually too small and low-scoring, but should likely just be merged.

    If two hits are close enough together, such that the distance between the start
    of the first and end of the second is within some tolerance (default +-10%) of the
    total length of a domains PSSM, this function will return True.

    Parameters
    ----------
    one : models.Domain
        Domain instance
    two : models.Domain
        Domain instance
    coverage_pct : float
        Conserved domain hit percentage coverage threshold. A hit is considered
        truncated if its total length is less than `coverage_pct` * CD length.
    tolerance_pct : float
        Percentage of CD length to use when calculating acceptable lower/upper bounds
        for combined domains.

    Returns
    -------
    True
        If the Domain instances are likely fragmented and should be combined.
    False
        If the Domain instances should be separate.
    """
    if one.type != two.type:
        raise ValueError("Expected Domain instances of same type")

    pssm_length = DOMAINS[one.domain]["length"]
    coverage = pssm_length * coverage_pct
    tolerance = pssm_length * tolerance_pct
    one_length, two_length = len(one), len(two)

    return (
        one_length < coverage
        and two_length < coverage
        and pssm_length - tolerance <= two.end - one.start <= pssm_length + tolerance
        and one_length + two_length > coverage
    )


def _filter_domains(domains, by="evalue", coverage_pct=0.5, tolerance_pct=0.1):
    """Filter overlapping Domain objects and test adjcency rules.

    Adjacency rules are tested again here, in case they are missed within overlap
    groups. For example, the NRPS-para261 domain is not always entirely contained by
    a condensation domain, so should be caught by this pass.

    Parameters
    ----------
    domains : list, tuple, set
        Collection of models.Domain instances to be filtered.
    by : str
        Which measure to use in `_filter_domain_group` (def. 'evalue').
    coverage_pct : float
        Conserved domain coverage percentage threshold to use in `_is_fragmented_domain`.
    tolerance_pct : float
        CD length tolerance percentage threshold to use in `_is_fragmented_domain`.

    Returns
    -------
    list
        models.Domain instances remaining after filtering.
    """

    domains = [
        _filter_domain_group(group, by) for group in _group_overlapping_hits(domains)
    ]

    i, total = 1, len(domains)
    while i < total:
        if i + 1 == total:
            break
        previous, current = domains[i - 1 : i + 1]

        # When domains are likely together, e.g. two small C domain hits right next
        # to each other, or multiple Methyltransf_X domains
        if previous.type == current.type and _is_fragmented_domain(
            previous, current, coverage_pct, tolerance_pct
        ):
            previous.end = current.end
            del domains[i]
            continue

        for new_type, rule in ADJACENCY_RULES.items():
            if rule(previous, current):
                previous.type = new_type
                del domains[i]
                break
        i += 1

    for domain in domains:
        if len(domain) < 0.6 * DOMAINS[domain.domain]["length"]:
            domain.truncated = True

    return domains

    # return [
    #     # Final filter, get rid of any obviously wrong, small hits; do here, not when
    #     # parsing table, so we don't discard potentially fragmented single domains
    #     domain
    #     for domain in domains
    #     if len(domain) > 0.2 * CD_LENGTHS[domain.domain]
    # ]


def _filter_domain_group(group, by="evalue"):
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
    by: str
        Measure to use when determining the best domain of the group. Choices:
        'bitscore': return domain with highest bitscore (relative to its threshold)
        'evalue': return domain with lowest E-value
        'length': return longest domain hit

    Returns
    -------
    Domain
        Highest scoring `Domain` in the group. If any special rules have been satisfied,
        the type of this `Domain` will be set to that rule (e.g. Condensation ->
        Epimerization).
    """
    key_functions = {
        "bitscore": (lambda d: d.bitscore / DOMAINS[d.domain]["bitscore"], True),
        "evalue": (lambda d: d.evalue, False),
        "length": (lambda d: d.end - d.start, True),
    }

    if by not in key_functions:
        raise ValueError("Expected 'bitscore', 'evalue' or 'length'")

    key, reverse = key_functions[by]

    container, *_group = sorted(group, key=key, reverse=reverse)

    for domain in _group:
        for new_type, rule in ADJACENCY_RULES.items():
            if rule(container, domain):
                container.type = new_type
                return container

    return container


def _hits_overlap(a, b, threshold=0.2):
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


def _group_overlapping_hits(domains, threshold=0.2):
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
            if _hits_overlap(current, future, threshold):
                group.append(future)  # add if contained
            else:
                yield group  # else yield to iterator
                break
            if j == total - 1:  # if reached the end, yield
                yield group
        i += len(group)  # move index ahead of last group


def parse(handle, **kwargs):
    """Parse CD-Search results.

    Any additional kwargs are passed to `synthases_from_results`.

    Parameters
    ----------
    handle: open file handle
        An open CD-Search results file handle. If you used the website to analyse your
        sequences, the file you should download is Domain hits, Data mode: Full, ASN
        text. When using a `CDSearch` object, this format is automatically selected.

    Returns
    -------
    list:
        A list of Synthase objects parsed from the results file.
    """
    return _filter_results(_parse_cdsearch_table(handle), **kwargs)
