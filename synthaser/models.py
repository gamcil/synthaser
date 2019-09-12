#!/usr/bin/env python3


from collections import defaultdict
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

        Refer to `assign_type` and `assign_subtype` for description of
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
        convention to name these T, i.e.::

            A-ACP-C --> A-T-C

        In hybrid PKS-NRPS, this replacement is made in the NRPS module of the synthase.
        Thus, this function looks for a condensation (C) domain that typically signals
        the beginning of such a module, and replaces any ACP with T after that domain.

        An example PKS-NRPS domain architecture may resemble::

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

    def extract_domains(self):
        """Extract all domains in this synthase.

        For example, given a Synthase:

        >>> synthase = Synthase(
        ...     header='synthase',
        ...     sequence='ACGT...',  # length 100
        ...     domains=[
        ...         Domain(type='KS', domain='PKS_KS', start=1, end=20),
        ...         Domain(type='AT', domain='PKS_AT', start=50, end=70)
        ...     ]
        ... )

        Then, we can call this function to extract the domain sequences:

        >>> synthase.extract_domains()
        {'KS':['ACGT...'], 'AT':['ACGT...']}

        Returns
        -------
        dict
            Sliced sequences for each domain in this synthase, keyed on domain type.

        Raises
        ------
        ValueError
            If the `Synthase` has no `Domain` objects.
        ValueError
            If the `sequence` attribute is empty.
        """
        if not self.domains:
            raise ValueError("Synthase has no domains")
        if not self.sequence:
            raise ValueError("Synthase has no sequence")
        domains = defaultdict(list)
        for domain in self.domains:
            domains[domain.type].append(domain.slice(self.sequence))
        return dict(domains)

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
        """Slice segment of sequence using the position of this Domain.

        Given a Domain:

        >>> domain = Domain(type='KS', subtype='PKS_KS', start=10, end=20)

        And its corresponding Synthase sequence:

        >>> synthase.sequence
        'ACGTACGTACACGTACGTACACGTACGTAC'

        We can extract the Domain:

        >>> domain.slice(synthase.sequence)
        'CGTACGTACA'
        """
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

    Synthases are assigned types based on the following rules:

    Hybrid (PKS-NRPS):
        Both beta-ketoacyl synthase (KS) and adenylation (A) domains
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


def wrap_fasta(sequence, limit=80):
    """Wrap FASTA record to 80 characters per line.

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


def create_fasta(header, sequence, wrap=80):
    """Create a FASTA format string from a header and sequence.

    For example:

    >>> create_fasta('header', 'AAAAABBBBBCCCCC', wrap=5)
    '>header\\nAAAAA\\nBBBBB\\nCCCCC'

    Parameters
    ----------
    header : str
        Name to use in FASTA definition line (i.e. >header).
    sequence : str
        The sequence corresponding to the `header`.
    wrap : int
        The number of characters per line for wrapping the given `sequence`.
        This function will call `wrap_fasta`.

    Returns
    -------
    str
        FASTA format string.
    """
    return ">{}\n{}".format(header, wrap_fasta(sequence, limit=wrap))


def extract_all_domains(synthases):
    """Extract all domain sequences in a list of `Synthase` objects.

    For example, given a list of `Synthase` objects:

    >>> synthases = [Synthase(header='one', ...), Synthase(header='two', ...)]

    Then, the output of this function may resemble:

    >>> domains = extract_all_domains(synthases)
    >>> domains
    {'KS': [('one_KS_1', 'IAIA...'), ('two_KS_1', 'IAIE...')], 'AT': [...]}

    We can easily write these to file in FASTA format. For example, to write all KS
    domain sequences to file, we open a file handle for writing and build a multiFASTA
    using `create_fasta`:

    >>> with open('output.faa', 'w') as out:
    ...     multifasta = '\\n'.join(
    ...         create_fasta(header, sequence)
    ...         for header, sequence in domains['KS'].items()
    ...     )
    ...     out.write(multifasta)

    Parameters
    ----------
    synthases : list
        A list of `Synthase` objects with non-empty `sequence` and `domains`
        attributes.

    Returns
    -------
    combined : dict
        Dictionary of extracted domain sequences keyed on domain type. Each domain is
        represented by a tuple consisting of a unique header, in the format
        `Synthase.header_Domain.type_index` where `index` is the index of that
        Domain in the Synthase, and the extracted sequence.
    """
    combined = defaultdict(list)
    for synthase in synthases:
        for type, sequences in synthase.extract_domains().items():
            combined[type].extend(
                (f"{synthase.header}_{type}_{i}", sequence)
                for i, sequence in enumerate(sequences)
            )
    return dict(combined)
