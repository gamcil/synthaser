#!/usr/bin/env python3
"""
Parse CD-search results to extract domain architectures of synthases.

Can also extract domain sequences for each synthase (--extract)
and generate a SVG figure of all synthases and their domain architecture
(--visual).

Author: Cameron Gilchrist
Date: 2018-06-12
"""

import json
import logging
import re

from itertools import groupby
from operator import attrgetter
from tempfile import NamedTemporaryFile as NTF

from synthaser.cdsearch import CDSearch

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


COLOURS = {
    "ACP": "#084BC6",
    "KS": "#08B208",
    "SAT": "#808080",
    "KR": "#089E4B",
    "MT": "#00ff00",
    "ER": "#089E85",
    "AT": "#DC0404",
    "DH": "#B45F04",
    "PT": "#999900",
    "TE": "#750072",
    "TR": "#9933ff",
    "T": "#084BC6",
    "R": "#9933ff",
    "C": "#393989",
    "A": "#56157F",
}

DOMAINS = {
    "KS": ["PKS_KS", "PKS"],
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
}


class Synthase:
    """The Synthase class stores a query protein sequence, its hit domains, and the
    methods for generating its SVG representation.

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

    def filter_overlapping_domains(self):
        """Filter overlapping Domains on this Synthase, saving best of each group."""
        self.domains = [
            max(group, key=lambda x: x.end - x.start)
            for group in group_overlapping_hits(self.domains)
        ]

    @property
    def sequence_length(self):
        return len(self.sequence)

    @property
    def architecture(self):
        """Return the domain architecture of this synthase as a hyphen separated string.

        Also performs any final type-specific cleanup (e.g. NRPS ACP domain -> T).
        """
        if self.type != "PKS":
            start, replace = 0, {"ACP": "T", "TR": "R"}

            if self.type == "Hybrid":
                for start, domain in enumerate(self.domains):
                    if domain.type == "C":
                        break
                    start += 1

            for domain in self.domains[start:]:
                if domain.type in replace:
                    domain.type = replace[domain.type]

        return "-".join(domain.type for domain in self.domains)

    @property
    def gradient(self):
        """Create a linearGradient SVG element based on this objects domain architecture.

        For example, if we define a Synthase object as follows:

        >>> synthase = Synthase(
        ...     header='synthase',
        ...     sequence='ACGACG...',  # length 100
        ...     domains=[Domain(type='KS', start=50, end=100)],
        ... )

        The generated gradient will be:

        >>> print(synthase.gradient)
        <linearGradient id="synthase_doms" x1="0%" y1="0%" x2="100%" y2="0%">
        <stop offset="50%" stop-color="white"/>
        <stop offset="50%" stop-color="#08B208"/>
        <stop offset="100%" stop-color="#08B208"/>
        <stop offset="100%" stop-color="white"/>
        </linearGradient>

        Note that a generated linearGradient will have 4 stops for every Domain object
        in the ``domains`` attribute; this ensures that each Domain colour has a hard
        edge instead of blending into the next Domain.
        """
        if not self.sequence:
            raise ValueError("Synthase has no sequence")

        stops = []
        sequence_length = len(self.sequence)

        for domain in self.domains:
            start_pct = int(domain.start / sequence_length * 100)
            end_pct = int(domain.end / sequence_length * 100)
            colour = COLOURS[domain.type]
            stops.append(
                f'<stop offset="{start_pct}%" stop-color="white"/>\n'
                f'<stop offset="{start_pct}%" stop-color="{colour}"/>\n'
                f'<stop offset="{end_pct}%" stop-color="{colour}"/>\n'
                f'<stop offset="{end_pct}%" stop-color="white"/>'
            )

        return (
            '<linearGradient id="{}_doms" x1="0%" y1="0%" x2="100%" y2="0%">\n{}\n'
            "</linearGradient>"
            "".format(self.header, "\n".join(stops))
        )

    def polygon(self, scale_factor=1, info_fsize=12, arrow_height=14):
        """Build SVG representation of one synthase.

        Length is determined by the supplied scale factor. Then, pairs of X and Y coordinates
        are calculated to represent each point in the synthase arrow. Finally, an SVG
        polygon feature is built with extra information above in a text feature, e.g.

        ::

            synthase, 100aa, KS-AT
            A----------B
            |           \\
            |            C
            |           /
            E----------D


        For example, if we define a Synthase object as follows:

        >>> synthase = Synthase(
        ...     header='synthase',
        ...     sequence='ACGACG...',  # length 100
        ...     domains=[Domain(type='KS', start=50, end=100)],
        ... )

        The generated polygon will be:

        >>> polygon = synthase.polygon()
        >>> print(polygon)
        <text dominant-baseline="hanging" font-size="12">synthase, 100aa, KS</text>'
        <polygon
            id="synthase" points="0,10.8,90,10.8,100,17.8,90,24.8,0,24.8"
            fill="url(#synthase_doms)" stroke="black" stroke-width="1.5"
        />

        Note that the ``fill`` attribute takes the form ``header``_doms; this is the id
        of the gradient for this Synthase.

        Parameters
        ----------
        scale_factor : int
            Scaling factor to multiply the Synthase sequence length by.

        Returns
        -------
        str
            <text> and <polygon> SVG features representing this Synthase. The fill for
            the polygon corresponds to the URL that is generated by the gradient property.
        """
        sequence_length = len(self.sequence)
        scaled_length = scale_factor * sequence_length
        info_fsize_scaled = info_fsize * 0.9
        bottom_y = info_fsize_scaled + arrow_height
        middle_y = info_fsize_scaled + arrow_height / 2

        ax, ay = 0, info_fsize_scaled
        bx, by = scaled_length - 10, info_fsize_scaled
        cx, cy = scaled_length, middle_y
        dx, dy = scaled_length - 10, bottom_y
        ex, ey = 0, bottom_y

        points = f"{ax},{ay},{bx},{by},{cx},{cy},{dx},{dy},{ex},{ey}"
        information = f"{self.header}, {sequence_length}aa, {self.architecture}"

        return (
            f'<text dominant-baseline="hanging" font-size="{info_fsize}">{information}</text>'
            f'<polygon id="{self.header}" points="{points}"'
            f' fill="url(#{self.header}_doms)" stroke="black"'
            ' stroke-width="1.5"/>'
        )


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

    @classmethod
    def from_cdsearch_row(cls, row):
        """Instantiate a new Domain from a row in a CD-search results file.

        For example, a typical row might looks like:

        >>> print(row)
        Q#1 - >AN6791.2\tspecific\t225858\t9\t1134\t0\t696.51\tCOG3321\tPksD\t - \tcl09938

        We can then instantiate a Domain object from this row, like so:

        >>> domain = Domain.from_cdsearch_row(row)
        >>> domain
        'PksD [KS] 9-1134'

        Parameters
        ----------
        row : str
            Tab-separated row from a CDSearch results file.

        Raises
        ------
        ValueError
            If the domain in this row is not in the DOMAINS dictionary.
        """
        _, _, _, start, end, _, _, _, domain, *_ = row.split("\t")
        for domain_type, domains in DOMAINS.items():
            if domain in domains:
                return cls(
                    type=domain_type, domain=domain, start=int(start), end=int(end)
                )
        raise ValueError(f"'{domain}' not a synthaser key domain")

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
        """Serialise this object to JSON."""
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_file):
        return cls(**json.load(json_file))


class Figure:
    """The Figure class acts as a repository for all Synthase objects, as well as
    holding all methods required for generating the final SVG figure.

    Attributes
    ----------
    synthases : list
        Synthase objects to be drawn.
    """

    def __init__(self, synthases=None):
        self.synthases = synthases if synthases else []

    def __repr__(self):
        return "\n\n".join(
            "{}\n{}\n{}".format(
                subtype,
                "-" * len(subtype),
                "\n".join(str(synthase) for synthase in group),
            )
            for subtype, group in self.iterate_synthase_types()
        )

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.synthases == other.synthases
        raise NotImplementedError

    def __getitem__(self, key):
        for synthase in self.synthases:
            if synthase.header == key:
                return synthase
        raise KeyError(f"No synthase with header '{key}'")

    def scale_factor(self, width):
        """Calculate the scale factor for drawing synthases.

        Parameters
        ----------
        width : int
            Total width, in pixels, of the SVG figure.

        Returns
        -------
        float :
            Scaling factor that will be used when calculating the width of each Synthase polygon.

        Raises
        ------
        ValueError
            If `width` is a negative number.
        """
        if width < 0:
            raise ValueError("Width must be greater than 0")
        largest = max(self.synthases, key=attrgetter("sequence_length"))
        return (width - 2) / largest.sequence_length

    def sort_synthases_by_length(self):
        """Sort Synthase objects by length of their sequences or domain architecture.

        Synthases are only sorted by domain architecture in the absence of sequences.
        """
        if any(not synthase.sequence for synthase in self.synthases):
            self.synthases.sort(key=lambda s: len(s.architecture), reverse=True)
        else:
            self.synthases.sort(key=attrgetter("sequence_length"), reverse=True)

    def iterate_synthase_types(self):
        """Group synthases by their types and yield.

        Synthases are first reverse sorted in-place by their sequence or architecture
        length. They are then sorted by their `type` and `subtype` attributes, grouped
        by `subtype` and yielded.

        Yields
        ------
        (str, list)
            A subtype and a list of Synthase objects with that subtype.
        """
        self.sort_synthases_by_length()
        self.synthases.sort(key=attrgetter("type", "subtype"))
        for subtype, group in groupby(self.synthases, key=attrgetter("subtype")):
            yield subtype, list(group)

    @staticmethod
    def build_polygon_block(
        subtype,
        synthases,
        scale_factor,
        arrow_spacing,
        arrow_height,
        info_fsize,
        header_fsize,
    ):
        """Generate the SVG for a block of Synthase objects.

        Parameters
        ----------
        subtype : str
            The subtype of the Synthase objects supplied to this method.
        synthases : list
            Synthase objects of a certain subtype to be visualised.

        See Figure.visualise() for description of other parameters.

        Returns
        -------
        block : str
            SVG of this Synthase subtype block.
        offset : int
            Cumulative total offset in this block. This is returned so the following
            block can be positioned below the one generated here.
        """

        block = (
            '<text dominant-baseline="hanging"'
            f' font-size="{header_fsize}"'
            f' font-weight="bold">{subtype}</text>'
        )
        offset = header_fsize
        for synthase in synthases:
            polygon = synthase.polygon(
                scale_factor, info_fsize=info_fsize, arrow_height=arrow_height
            )
            block += f'\n<g transform="translate(1,{offset})">\n{polygon}\n</g>'
            offset += info_fsize + arrow_height + 4 + arrow_spacing
        return block, offset

    def visualise(
        self,
        arrow_height=12,
        arrow_spacing=4,
        block_spacing=16,
        header_fsize=15,
        info_fsize=12,
        width=600,
    ):
        """Build the SVG figure.

        Parameters
        ----------
        arrow_height : int
            The height, in pixels, of the generated polygon for each Synthase.
        arrow_spacing : int
            Vertical spacing, in pixels, to insert between each Synthase polygon.
        block_spacing : int
            Vertical spacing, in pixels, to insert between each subtype block of
            Synthases.
        header_fsize : int
            Font size of Synthase type headers.
        info_fsize : int
            Font size of Synthase information subheaders.
        width : int
            Width, in pixels, of the final generated SVG.

        Returns
        -------
        str
            Final SVG figure.
        """
        scale_factor = self.scale_factor(width)

        blocks = ""
        offset = 3

        for subtype, synthases in self.iterate_synthase_types():
            log.info("Subtype=%s, %i synthases", subtype, len(synthases))
            block, height = self.build_polygon_block(
                subtype,
                synthases,
                scale_factor,
                arrow_spacing,
                arrow_height,
                info_fsize,
                header_fsize,
            )
            blocks += f'<g transform="translate(0,{offset})">\n{block}\n</g>'
            offset += height + block_spacing

        return '<svg width="{}" height="{}">\n{}\n{}\n</svg>'.format(
            width,
            offset - block_spacing - arrow_spacing - 4,
            "\n".join(synthase.gradient for synthase in self.synthases),
            blocks,
        )

    def to_json(self):
        """Serialise Figure to JSON."""
        return json.dumps([synthase.to_dict() for synthase in self.synthases])

    @classmethod
    def from_json(cls, json_file):
        """Load Figure from an open JSON file handle."""
        return cls([Synthase.from_dict(record) for record in json.load(json_file)])

    @classmethod
    def from_cdsearch_results(cls, results_handle, query_handle=None):
        """Instantiate a new Figure from CD-search results file.

        results_handle should be an open file handle.

        Parameters
        ----------
        results_handle: open file handle
            An open CD-Search results file handle. If you used the website to analyse your
            sequences, the file you should download is Domain hits, Data mode: Full, ASN
            text. When using a CDSearch object, this format is automatically selected.

        query_handle: open file handle, optional
            An open file handle for the sequences used in the CD-search run. Sequences
            can be added later via Figure.add_query_sequences().

        Returns
        -------
        Figure:
            A Figure object built from the CDSearch results.
        """
        figure = cls()
        _query = ""
        synthase = None
        pattern = re.compile(r"Q#\d+? - [>]?(.+?)\t")

        for row in results_handle:
            try:
                row = row.decode()
            except AttributeError:
                pass  # in case rows are unicode

            if not row.startswith("Q#") or row.isspace():
                continue

            query = pattern.search(row).group(1)
            if query != _query:
                if synthase:
                    figure.synthases.append(synthase)
                synthase = Synthase(query)
                _query = query

            try:
                domain = Domain.from_cdsearch_row(row)
            except ValueError:
                continue

            synthase.domains.append(domain)

        if synthase:  # add the final synthase from the loop
            figure.synthases.append(synthase)

        for synthase in figure.synthases:
            synthase.filter_overlapping_domains()

            try:
                synthase.type = assign_type(synthase)
            except ValueError:
                log.warning("%s has no identifying domains", synthase.header)

            try:
                synthase.subtype = assign_subtype(synthase)
            except ValueError:
                log.warning(
                    "%s (type: %s) was not assigned a subtype",
                    synthase.header,
                    synthase.type,
                )

        if query_handle:
            figure.add_query_sequences(query_handle)

        figure.sort_synthases_by_length()
        return figure

    @classmethod
    def from_cdsearch(cls, query_file, **kwargs):
        """Convenience function to directly instantiate a Figure from a new CDSearch job.

        Parameters
        ----------
        query_file : str
            Path to a FASTA file containing query sequences to be analysed.

        All additional keyword arguments are passed to CDSearch.run().

        Returns
        -------
        Figure
            Figure built from the CDSearch query.
        """
        cd = CDSearch()
        with NTF() as results:
            cd.run(query_file=query_file, output_file_handle=results, **kwargs)
            results.seek(0)
            figure = cls.from_cdsearch_results(results)

        with open(query_file) as query_handle:
            figure.add_query_sequences(query_handle=query_handle)

        return figure

    def add_query_sequences(self, query_handle=None, sequences=None):
        """Add sequences from query FASTA file to the Figure.

        Parameters
        ----------
        query_handle : open file handle
            Open file handle of a FASTA file containing sequences corresponding to the
            Synthases in this objects `synthases` attribute.

        sequences : dict
            A pre-populated dictionary containing sequences corresponding to the
            Synthases in this objects `synthases` attribute.
        """
        if query_handle and not sequences:
            sequences = parse_fasta(query_handle)
        for header, sequence in sequences.items():
            try:
                self[header].sequence = sequence
            except KeyError as exc:
                raise KeyError(
                    f"Could not match '{header}' to synthase in results"
                ) from exc


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


def parse_fasta(fasta):
    """Parse an open FASTA file for sequences.

    Parameters
    ----------
    fasta : str
        Either an open file handle of a FASTA file, or newline split string (e.g. read
        in via readlines()) that can be iterated over.

    Returns
    -------
    sequences : dict
        All sequences in the FASTA file keyed on their header lines. For example,
            >sequence
            ACGTACGTACGT

        Will be stored as:
            {"sequence": "ACGTACGTACGT"}
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


def assign_type(synthase):
    """Determine the broad biosynthetic type of a Synthase.

    Classification rules
    --------------------
    Hybrid (PKS-NRPS): both a beta-ketoacyl synthase (KS) and adenylation (A) domains
    Polyketide synthase (PKS): a KS domain
    Nonribosomal peptide synthase (NRPS): an A domain

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
    types = set(domain.type for domain in synthase.domains)
    if {"KS", "A"}.issubset(types):
        return "Hybrid"
    if "KS" in types:
        return "PKS"
    if "A" in types:
        return "NRPS"
    raise ValueError("Could not find an identifying domain")


def assign_subtype(synthase):
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
    types = set(domain.type for domain in synthase.domains)
    if synthase.type == "PKS":
        subtypes = [
            ("HR-PKS", all, {"ER", "KR", "DH"}),
            ("PR-PKS", any, {"ER", "KR", "DH"}),
            ("NR-PKS", all, {"KS", "AT"}),
            ("PKS-like", any, {"KS"}),
        ]
    elif synthase.type == "NRPS":
        subtypes = [("NRPS", all, {"A", "T", "C"}), ("NRPS-like", any, {"A"})]
    else:
        return synthase.type

    for subtype, function, required in subtypes:
        if function(domain in types for domain in required):
            return subtype
    return "Other"
