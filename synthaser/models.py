import json
import logging

from collections import defaultdict, UserList


LOG = logging.getLogger(__name__)


class Serialiser:
    def to_dict(self):
        raise NotImplementedError

    @classmethod
    def from_dict(cls, d):
        raise NotImplementedError

    def to_json(self, fp=None, **kwargs):
        if fp:
            json.dump(self.to_dict(), fp, **kwargs)
        else:
            return json.dumps(self.to_dict(), **kwargs)

    @classmethod
    def from_json(cls, js, **kwargs):
        if isinstance(js, str):
            d = json.loads(js, **kwargs)
        else:
            d = json.load(js, **kwargs)
        return cls.from_dict(d)


class Synthase(Serialiser):
    """The Synthase class stores a query protein sequence, its hit domains, and the
    methods for filtering and classifying.

    Attributes:
        header (str): Synthase name.
        sequence (str): Amino acid sequence of this Synthase.
        domains (list): Conserved domain hits in this Synthase.
        classification (list): All classification rules satisfied.
    """

    def __init__(
        self,
        header=None,
        sequence=None,
        domains=None,
        classification=None
    ):
        self.header = header
        self.sequence = sequence
        self.domains = domains if domains else []
        self.classification = classification if classification else []

    def __str__(self):
        return f"{self.header}\t{self.architecture}"

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.header == other.header
        raise NotImplementedError

    def contains(self, classes=None, types=None, families=None):
        """Checks if Synthase contains given classifications, domain
        families or types."""
        return (
            classes and not set(classes).isdisjoint(self.classification)
            or types and not set(types).isdisjoint(d.type for d in self.domains)
            or families and (
                not set(families).isdisjoint(d.accession for d in self.domains)
                or not set(families).isdisjoint(d.domain for d in self.domains)
            )
        )

    def extract_domains(self, types=None, families=None):
        """Extract specific domain type/family sequences from this Synthase.
        """
        if not self.domains:
            raise ValueError("Synthase has no domains")
        if not self.sequence:
            raise ValueError("Synthase has no sequence")

        # If nothing specified, extract all
        if not (types or families):
            return self.extract_all_domains()

        output = defaultdict(lambda: defaultdict(list))

        for domain in self.domains:
            # Slice from parent sequence
            sequence = domain.slice(self.sequence)

            # Domain types
            if types and domain.type in types:
                output["types"][domain.type].append(sequence)

            # Specific domain families
            if families and (
                domain.domain in families
                or domain.accession in families
            ):
                output["families"][domain.domain].append(sequence)

        return output

    def extract_all_domains(self):
        """Extracts all domain sequences from this synthase.

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

        >>> synthase.extract_all_domains()
        {'KS':['ACGT...'], 'AT':['ACGT...']}

        Returns:
            dict: Sequences for each domain in this synthase keyed on domain type.
        Raises:
            ValueError: If the Synthase has no Domain objects.
            ValueError: If the sequence attribute is empty.
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
            "classification": self.classification,
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

    def to_long(self, delimiter=","):
        fields = [
            self.header,
            str(self.sequence_length),
            self.architecture,
            ", ".join(self.classification)
        ]
        return delimiter.join(fields)

    @property
    def sequence_length(self):
        return len(self.sequence)

    @property
    def architecture(self):
        return "-".join(str(domain) for domain in self.domains)

    @property
    def domain_types(self):
        return [domain.type for domain in self.domains]


class Domain(Serialiser):
    """A conserved domain hit.

    Attributes:
        type (str): Broader domain type (e.g. KS)
        domain (str): Specific CDD family (e.g. PKS_KS)
        start (int): Start of domain hit in parent sequence
        end (int): End of domain hit in parent sequence
        evalue (float): Domain hit E-value
        bitscore (float): Domain hit bitscore
        partial (str): Domain hit partiality ('C', 'N' or 'NC')
        accession (str): CDD accession of domain family
        superfamily (str): CDD accession of domain superfamily
    """

    def __init__(
        self,
        pssm=None,
        type=None,
        domain=None,
        start=None,
        end=None,
        evalue=None,
        bitscore=None,
        partial=False,
        accession=None,
        superfamily=None,
    ):
        self.pssm = pssm
        self.type = type
        self.domain = domain
        self.start = start
        self.end = end
        self.evalue = evalue
        self.bitscore = bitscore
        self.partial = partial
        self.accession = accession
        self.superfamily = superfamily

    def __str__(self):
        return self.type

    # TODO: rework, or move this to e.g. domain1.equals(domain2)
    #       since this messes with object comparison e.g. domain1 in [...]
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
        """Slices segment of sequence using the position of this Domain.

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
        return {
            "pssm": self.pssm,
            "type": self.type,
            "domain": self.domain,
            "start": self.start,
            "end": self.end,
            "evalue": self.evalue,
            "bitscore": self.bitscore,
            "partial": self.partial,
            "accession": self.accession,
            "superfamily": self.superfamily,
        }


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
        previous = []
        lines = []
        for synthase in sorted(
            self, key=lambda s: (s.classification, -s.sequence_length)
        ):
            if previous != synthase.classification:
                previous = synthase.classification
                if lines:
                    lines[-1] += "\n"
                line = " --> ".join(previous)
                under = "-" * len(line)
                lines.extend([line, under])
            lines.append(str(synthase))
        return "\n".join(lines)

    def to_long(self, delimiter=",", headers=True):
        """Generate summary of the container in long data format.

        For example:

        ::

            Synthase  Length (aa)  Architecture        Classification
            SEQ001.1  1000         KS-AT-DH-ER-KR-ACP  PKS, Type I, Highly-reducing

        NOTE: actual output is character delimited, not human readable.
        """
        synthases = [synthase.to_long(delimiter) for synthase in self]
        if headers:
            _headers = ("Synthase", "Length (aa)", "Architecture", "Classification")
            synthases.insert(0, delimiter.join(_headers))
        return "\n".join(synthases)

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

    def filter(self, classes=None, types=None, families=None):
        filtered = [
            synthase
            for synthase in self
            if synthase.contains(
                classes=classes,
                types=types,
                families=families,
            )
        ]
        return type(self)(synthases=filtered)

    def extract_synthases(self, classes=None, types=None, families=None):
        """Bin entire synthase sequences."""
        result = {
            "types": defaultdict(list),
            "classes": defaultdict(list),
            "families": defaultdict(list),
        }
        for synthase in self:
            for key, values in zip(
                ["classes", "types", "families"],
                [classes, types, families],
            ):
                if not values:
                    continue
                for value in values:
                    if synthase.contains(**{key: [value]}):
                        result[key][value].append(synthase)
        return result

    def extract_domains(self, classes=None, types=None, families=None, by="sequence"):
        """Extract domain sequences from Synthase objects in this container.

        For example, given a `SynthaseContainer` containing `Synthase` objects:

        >>> synthases = [Synthase(header='one', ...), Synthase(header='two', ...)]
        >>> container = SynthaseContainer(synthases)

        Then, the output of this function may resemble:

        >>> container.extract_domains()
        {'KS': [('one_KS_1', 'IAIA...'), ('two_KS_1', 'IAIE...')], 'AT': [...]}
        """
        if by == "sequence":
            result = {}
        elif by == "query":
            result = {
                "types": defaultdict(dict),
                "families": defaultdict(dict),
            }
        else:
            raise ValueError("Expected by='sequence' or by='query'")

        for synthase in self:
            if classes and not synthase.contains(classes=classes):
                LOG.warning("Sequence not one of: %s, skipping", classes)
                continue
            if not synthase.domains:
                LOG.warning("%s has no domains, skipping", synthase.header)
                continue
            sequences = synthase.extract_domains(
                types=types,
                families=families
            )
            if not sequences:
                LOG.warning("No domains extracted (%s), skipping", synthase.header)
                continue
            if by == "sequence":
                result[synthase.header] = sequences
            elif by == "query":
                for key in result:
                    if key not in sequences:
                        continue
                    for label, extracts in sequences[key].items():
                        if label not in result[key]:
                            result[key][label] = {}
                        result[key][label][synthase.header] = extracts
        return result

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
