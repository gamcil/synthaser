"""
This module stores the classes used throughout `synthaser`.

The `Domain` class represents a conserved domain hit. It stores the broader domain type,
the specific conserved domain profile name (from CDD), as well as its position in its
parent synthase sequence and score from the search. It also provides methods for slicing
the corresponding sequence and serialisation. We can instantiate a `Domain` object like
so:

>>> from synthaser.models import Domain
>>> domain = Domain(
...     type='KS',
...     domain='PKS_KS',
...     start=756,
...     end=1178,
...     evalue=0.0,
...     bitscore=300
... )

and get its sequence given the parent `Synthase` object sequence:

>>> domain.slice(synthase.sequence)
'MPIAVGM..'

Likewise, the `Synthase` class stores information about a synthase, including its name,
amino acid sequence, `Domain` instances and its biosynthetic type and subtype. It also
contains methods for generating the domain architecture, extraction of domain sequences
and more. For example, we can instantiate a new `Synthase` object like so:

>>> from synthaser.models import Synthase
>>> synthase = Synthase(
...     header='SEQ001.1',
...     sequence='MASGTC...',
...     domains=[
...         Domain(type='KS'),
...         Domain(type='AT'),
...         Domain(type='DH'),
...         Domain(type='ER'),
...         Domain(type='KR'),
...         Domain(type='ACP'),
...     ],
...     type='Type I PKS',
...     subtype='HR-PKS'
... )

Then, we can generate the domain architecture:

>>> synthase.architecture
'KS-AT-DH-ER-KR-ACP'

Or extract all of the domain sequences:

>>> synthase.extract_domains()
{
    "KS_0": "MPIAVGM...",
    "AT_0": "VFTGQGA...",
    "DH_0": "DLLGVPV...",
    "ER_0": "DVEIQVS...",
    "KR_0": "IAENMCS...",
    "ACP_0": "ASTTVAQ..."
}

The object can also be serialised to JSON (note the `Domain` object works the same way):

>>> js = synthase.to_json()
>>> with open('synthase.json', 'w') as handle:
...     handle.write(js)

and subsequently loaded from JSON:

>>> with open('synthase.json') as handle:
...     synthase = Synthase.from_json(handle)

This will internally convert the `Synthase` object, as well as any `Domain` objects it
contains, to dictionaries, before converting to JSON using the builtin json library
and writing to file. When loading up from JSON, this process is reversed, and the
entries in the file are converted back to Python objects.
"""

from collections import defaultdict, UserList
from itertools import groupby
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

    def __init__(
        self, header=None, sequence=None, domains=None, type=None, subtype=None
    ):
        self.header = header
        self.sequence = sequence
        self.domains = domains if domains else []
        self.type = type if type else ""
        self.subtype = subtype if subtype else ""

    def __str__(self):
        return f"{self.header}\t{self.architecture}"

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.header == other.header
        raise NotImplementedError

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

    def to_fasta(self):
        return f">{self.header}\n{self.sequence}"

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
        return [domain.type for domain in self.domains]


class Domain:
    """Store a conserved domain hit."""

    def __init__(
        self, type=None, domain=None, start=None, end=None, evalue=None, bitscore=None
    ):
        self.type = type
        self.domain = domain
        self.start = start
        self.end = end
        self.evalue = evalue
        self.bitscore = bitscore

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

    def __len__(self):
        return self.end - self.start

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


class SynthaseContainer(UserList):
    """Simple container class for Synthase objects.

    The purpose of this class is to facilitate batch actions on Synthase objects, i.e.
    serialisation, extraction of domain sequences, iteration over type/subtype, and
    printing summaries.
    """

    def __init__(self, synthases):
        UserList.__init__(self)
        self.extend(synthases)

    def __add__(self, container):
        if not isinstance(container, SynthaseContainer):
            raise TypeError("Expected SynthaseContainer object")
        copy = self.copy()
        copy.extend(container)
        return copy

    def __str__(self):
        return "\n\n".join(
            "{}\n{}\n{}".format(
                subtype,
                "-" * len(subtype),
                "\n".join(str(synthase) for synthase in group),
            )
            for subtype, group in self.subtypes()
        )

    def append(self, synthase):
        if not isinstance(synthase, Synthase):
            raise TypeError("Expected Synthase object")
        self.data.append(synthase)

    def extend(self, synthases):
        for synthase in synthases:
            self.append(synthase)

    def get(self, header):
        for synthase in self:
            if synthase.header == header:
                return synthase
        raise KeyError(f"No Synthase object with header: '{header}'")

    def to_json(self, handle, **kwargs):
        """Serialise this container to JSON."""
        dicts = [synthase.to_dict() for synthase in self]
        json.dump(dicts, handle, **kwargs)

    @classmethod
    def from_json(cls, handle):
        """Load Synthase objects from JSON."""
        return cls(Synthase.from_dict(entry) for entry in json.load(handle))

    def _attr_iter(self, attr):
        """Iterate over Synthase objects inside the SynthaseContainer, grouped by some
        attribute.
        """
        if attr not in ("type", "subtype"):
            raise ValueError("Expected 'type' or 'subtype'")
        synthases = sorted(self.data, key=lambda s: getattr(s, attr))
        for key, group in groupby(synthases, key=lambda s: getattr(s, attr)):
            yield key, type(self)(group)

    def subtypes(self):
        """Iterate over Synthase subtypes.

        Yields
        ------
        subtype : str
            Subtype of Synthase objects in the group
        synthases : SynthaseContainer
            SynthaseContainer object containing Synthase objects
        """
        yield from self._attr_iter("subtype")

    def types(self):
        """Iterate over Synthase types.

        Yields
        ------
        type : str
            Type of Synthase objects in the group
        synthases : SynthaseContainer
            SynthaseContainer object containing Synthase objects
        """
        yield from self._attr_iter("type")

    def extract_domains(self):
        """Extract domain sequence from Synthase objects in this container.

        For example, given a `SynthaseContainer` containing `Synthase` objects:

        >>> synthases = [Synthase(header='one', ...), Synthase(header='two', ...)]
        >>> container = SynthaseContainer(synthases)

        Then, the output of this function may resemble:

        >>> container.extract_domains()
        {'KS': [('one_KS_1', 'IAIA...'), ('two_KS_1', 'IAIE...')], 'AT': [...]}
        """
        combined = defaultdict(list)
        for synthase in self:
            for type, sequences in synthase.extract_domains().items():
                combined[type].extend(
                    (f"{synthase.header}_{type}_{i}", sequence)
                    for i, sequence in enumerate(sequences)
                )
        return dict(combined)

    def add_sequences(self, sequences):
        """Add amino acid sequence to Synthase objects in this container."""
        for header, sequence in sequences.items():
            self.get(header).sequence = sequence

    def to_fasta(self):
        return "\n".join(synthase.to_fasta() for synthase in self)

    @classmethod
    def from_sequences(cls, sequences):
        """Build a SynthaseContainer from a dictionary of query sequences."""
        return cls(
            Synthase(header=key, sequence=value) for key, value in sequences.items()
        )
