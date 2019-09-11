#!/usr/bin/env python3


from operator import attrgetter
import json


class Synthase:
    """The Synthase class stores a query protein sequence, its hit domains, and the
    methods for filtering and classifying.

    Parameters
    ----------
    header : str
        Name of this Synthase. This must be equal to what is used in NCBI CD-search.
    sequence : str
        Amino acid sequence of this Synthase.
    domains : list
        Conserved domain hits in this Synthase.
    type : str
        Type of synthase; 'PKS', 'NRPS' or 'Hybrid'
    subtype : str
        Subtype of synthase, e.g. HR-PKS.
    """

    __slots__ = ("header", "sequence", "domains", "type", "subtype")

    def __init__(
        self, header=None, sequence=None, domains=None, type=None, subtype=None
    ):
        self.header = header
        self.sequence = sequence
        self.domains = domains if domains else []
        self.type = type if type else ""
        self.subtype = subtype if subtype else ""

    def __repr__(self):
        return f"{self.header}\t{self.architecture}"

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.header == other.header
        raise NotImplementedError

    def filter_overlapping_domains(self):
        """Filter overlapping Domains on this Synthase, saving best of each group.

        Parameters
        ----------
        synthase : dict
            Dictionary representation of a query synthase.
        """
        self.domains = [
            max(group, key=lambda x: x.end - x.start)
            for group in group_overlapping_hits(self.domains)
        ]

    def classify(self):
        """Assign a type and subtype to this Synthase.

        Refer to ``assign_type`` and ``assign_subtype`` for description of
        classification rules.
        """
        domains = self.domain_types
        self.type = assign_type(domains)
        self.subtype = assign_subtype(self.type, domains)

    def rename_nrps_domains(self):
        """Replace domain types in Hybrid and NRPS Synthases.

        The acyl carrier protein (ACP) domain in PKSs is homologous to the thioester
        domain of the peptide carrier protein (PCP) domain in NRPSs, and as such, both
        PKS and NRPS will report the same conserved domain hit. In NRPS, it is
        convention to name these T, i.e.

        ::

            A-ACP-C --> A-T-C

        In hybrid PKS-NRPS, this replacement is made in the NRPS module of the synthase.
        Thus, this function looks for a condensation (C) domain that typically signals
        the beginning of such a module, and replaces any ACP with T after that domain.

        An example PKS-NRPS domain architecture may resemble:

        ::

            KS-AT-DH-ER-KR-ACP-C-A-T-R

        Lastly, thioester reductase (TR) domains are generally written as R in NRPS,
        thus the replacement here.
        """
        if not self.type or self.type == "PKS":
            return

        start, replace = 0, {"ACP": "T", "TR": "R"}

        if self.type == "Hybrid":
            for start, domain in enumerate(self.domains):
                if domain.type == "C":
                    break

        for domain in self.domains[start:]:
            if domain.type in replace:
                domain.type = replace[domain.type]

    def to_dict(self):
        return {
            "header": self.header,
            "sequence": self.sequence,
            "domains": [domain.to_dict() for domain in self.domains],
            "type": self.type,
            "subtype": self.subtype,
        }

    @classmethod
    def from_dict(cls, dic):
        synthase = cls()
        for key, value in dic.items():
            if key == "domains":
                synthase.domains = [Domain(**domain) for domain in value]
            else:
                setattr(synthase, key, value)
        return synthase

    def to_json(self):
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_file):
        return Synthase.from_dict(json.load(json_file))

    @property
    def sequence_length(self):
        return len(self.sequence)

    @property
    def architecture(self):
        """Return the domain architecture of this synthase as a hyphen separated string."""
        return "-".join(domain.type for domain in self.domains)

    @property
    def domain_types(self):
        return set(domain.type for domain in self.domains)


class Domain:
    """Store a conserved domain hit."""

    __slots__ = ("type", "domain", "start", "end")

    def __init__(self, type=None, domain=None, start=None, end=None):
        self.type = type
        self.domain = domain
        self.start = start
        self.end = end

    def __repr__(self):
        return f"{self.domain} [{self.type}] {self.start}-{self.end}"

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return (
                self.type == other.type
                and self.domain == other.domain
                and self.start == other.start
                and self.end == other.end
            )
        raise NotImplementedError

    def slice(self, sequence):
        """Slice segment of sequence using the position of this Domain."""
        return sequence[self.start - 1 : self.end]

    def to_dict(self):
        """Serialise this object to dict of its attributes.

        For example, if we define a Domain:

        >>> domain = Domain(type='KS', domain='PksD', start=9, end=1143)

        We can serialise it to a Python dictionary:

        >>> domain.to_dict()
        {"type": "KS", "domain": "PksD", "start": 9, "end": 1143}
        """
        return {
            "type": self.type,
            "domain": self.domain,
            "start": self.start,
            "end": self.end,
        }

    def to_json(self):
        """Serialise this object to JSON.

        This function calls json.dumps() on Domain.to_dict().
        """
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_file):
        return cls(**json.load(json_file))


def assign_type(domains):
    """Determine the broad biosynthetic type of a Synthase.

    Classification rules
    --------------------
    Hybrid (PKS-NRPS):
        Both a beta-ketoacyl synthase (KS) and adenylation (A) domains
    Polyketide synthase (PKS):
        KS domain
    Nonribosomal peptide synthase (NRPS):
        A domain

    Parameters
    ----------
    synthase : Synthase
        A Synthase object with domain hits.

    Returns
    -------
    str
        The biosynthetic type of the given Synthase (hybrid, pks, nrps).

    Raises
    ------
    ValueError
        If no identifying domain (KS, A) is found.
    """
    if {"KS", "A"}.issubset(domains):
        return "Hybrid"
    if "KS" in domains:
        return "PKS"
    if "A" in domains:
        return "NRPS"
    raise ValueError("Could not find an identifying domain")


def assign_subtype(type, domains):
    """Determine the biosynthetic subtype of a Synthase.

    Subtypes are determined by the following rules:

    Polyketide synthase (PKS):

    1) Highly-reducing (HR-PKS): enoyl-reductase (ER), keto-reductase (KR) and
       dehydratase (DH)
    2) Partially-reducing (PR-PKS): any, but not all, reducing domains used to
       classify a HR-PKS
    3) Non-reducing (NR-PKS): no reducing domains, but beta-ketoacyl synthase (KS)
       and acyltransferase (AT) present
    4) PKS-like: at least a KS domain

    Nonribosomal peptide synthetase (NRPS):

    1) NRPS: full NRPS module, consisting of adenylation (A), peptidyl-carrier (PCP,
       aka T) and condensation (C) domains
    2) NRPS-like: at least an A domain

    If a non-PKS/NRPS Synthase is supplied, then this function will return its `type`
    attribute.

    Parameters
    ----------
    synthase : Synthase
        A Synthase object with domain hits.

    Returns
    -------
    str
        The biosynthetic subtype of the given Synthase.

    Raises
    ------
    ValueError
        Synthase has type other than "pks" or "nrps" or no subtype could be assigned.

    """
    if type == "PKS":
        subtypes = [
            ("HR-PKS", all, {"ER", "KR", "DH"}),
            ("PR-PKS", any, {"ER", "KR", "DH"}),
            ("NR-PKS", all, {"KS", "AT"}),
            ("PKS-like", any, {"KS"}),
        ]
    elif type == "NRPS":
        subtypes = [("NRPS", all, {"A", "T", "C"}), ("NRPS-like", any, {"A"})]
    else:
        return type
    for subtype, function, required in subtypes:
        if function(domain in domains for domain in required):
            return subtype
    return "Other"


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
