#!/usr/bin/env python3

# TODO:
#   1. Add scale bar (1kb)

"""
This module handles the generation of the SVG figure.

All functionality is available through use of the `Figure` class, which can be
instantiated simply with a collection of `models.Synthase` instances:

>>> from synthaser.models import Synthase
>>> from synthase.figure import Figure
>>> synthases = [
...     models.Synthase(...),
...     models.Synthase(...),
...     models.Synthase(...),
...     ...
...     ]
>>> figure = Figure(synthases)

In order to generate the figure, the object needs to know the lengths of each synthase
sequence. To do this, add the sequences using `add_query_sequences`:

>>> with open('synthases.fa') as handle:
>>>     figure.add_query_sequences(handle)

Now, we can create the SVG:

>>> svg = figure.visualise()
>>> with open('synthases.svg', 'w') as handle:
...     handle.write(svg)

Internally, the `Figure` class calls a number of methods to assemble an SVG-format
string, which is then written to the file.

Synthases are drawn as `<polygon>` elements, which are scaled to the correct width given
the document width and length of the largest sequence in the figure.

The domain colouring is done via `<linearGradient>` elements, which define colour fill
blocks based on the start and end of each domain object. This means that the shapes can
be resized after generating the figure (e.g. in InkScape), and the domain colours will always
correctly correspond to the synthase length.

The colour scheme can be manipulated by changing the `colours` attribute, which is a
dictionary mapping `Domain` type to hexadecimal colour codes.
"""


import json
import logging
import re

from itertools import groupby
from operator import attrgetter

from synthaser import fasta, results
from synthaser.models import Synthase
from synthaser.ncbi import CDSearch, efetch_sequences

LOG = logging.getLogger(__name__)


class Figure:
    """The Figure class acts as a repository for all Synthase objects, as well as
    holding all methods required for generating the final SVG figure.

    Attributes
    ----------
    synthases : list
        Synthase objects to be drawn.
    colours: dict
        Colourscheme to use when visualising domain architecture of Synthases.
    """

    default_colours = {
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
        "E": "#FFA500",
        "cAT": "#FF6600",
    }

    default_config = {
        "arrow_height": 12,
        "arrow_spacing": 4,
        "block_spacing": 16,
        "header_fsize": 15,
        "info_fsize": 12,
        "width": 600,
    }

    def __init__(self, synthases=None, colours=None, config=None):
        self.synthases = synthases if synthases else []

        self.colours = self.default_colours.copy()
        if colours:
            self.set_colours(colours)

        self.config = self.default_config.copy()
        if config:
            self.set_config(config)

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

    def get_by(self, attribute, key):
        if attribute not in ("type", "subtype"):
            raise ValueError("Expected 'type' or 'subtype'")
        return [
            synthase
            for synthase in self.synthases
            if getattr(synthase, attribute) == key
        ]

    @property
    def headers(self):
        return [synthase.header for synthase in self.synthases]

    def set_colours(self, colours):
        """Change colour hex codes for domain types.

        Valid domain types are:

        `ACP, KS, SAT, KR, MT, ER, AT, DH, PT, TE, TR, T, R, C, A`

        Thus, a valid dictionary might look like:

        >>> colours = {
        ...     "ACP": "#000000",
        ...     "KS": "#FFFFFF"
        ... }

        Parameters
        ----------
        colours : dict
            A dictionary of colour hex codes keyed on domain type. These are what the
            Figure class will use when creating the gradient fill of each Synthase
            polygon.

        Raises
        ------
        TypeError
            If `colours` is not a dictionary.
        KeyError
            If a key in `colours` is not a valid domain type.
        ValueError
            If a value in `colours` is not a valid hex code.
        """
        if not isinstance(colours, dict):
            raise TypeError("Expected dict")

        for key, value in colours.items():
            if key not in self.colours:
                raise KeyError(f"Invalid domain '{key}'")
            if not validate_colour(value):
                raise ValueError(f"'{key}' is not a valid hex code")
            self.colours[key] = value

    def set_config(
        self,
        arrow_height=None,
        arrow_spacing=None,
        block_spacing=None,
        header_fsize=None,
        info_fsize=None,
        width=None,
    ):
        """Set visualisation parameters on this Figure."""
        params = {
            "arrow_height": arrow_height,
            "arrow_spacing": arrow_spacing,
            "block_spacing": block_spacing,
            "header_fsize": header_fsize,
            "info_fsize": info_fsize,
            "width": width,
        }
        for key, value in params.items():
            if value:
                self.config[key] = value

    def calculate_scale_factor(self):
        """Calculate the scale factor for drawing synthases.

        The scale factor is calculated such that the largest Synthase in the Figure will
        match the value of `width`.

        For example, given two Synthases:

        >>> a = Synthase(sequence='ACGT...')  # sequence_length == 2000
        >>> b = Synthase(sequence='ACGT...')  # sequence_length == 1000

        We can instantiate a Figure and compute the scaling factor for a document of
        `width` 1000:

        >>> figure = Figure(synthases=[a, b])
        >>> figure.scale_factor(1000)
        0.499

        i.e. The larger Synthase with `sequence_length` 2000 is multipled by 0.499 to
        scale it to the width of the Figure.

        Note that the scaling factor is calculated with a slight offset (2) to account
        for the borders of each polygon being drawn outside of their strict width and
        height.

        Parameters
        ----------
        width : int
            Total width, in pixels, of the SVG figure.

        Returns
        -------
        float
            Scaling factor that will be used when calculating the width of each Synthase polygon.

        Raises
        ------
        ValueError
            If the Synthase objects in this object have empty `sequence` attributes.
        ValueError
            If `width` is a negative number.
        """
        if any(not synthase.sequence for synthase in self.synthases):
            raise ValueError("Synthases in this Figure have no sequences")
        if self.config["width"] < 0:
            raise ValueError("Width must be greater than 0")
        largest = max(self.synthases, key=attrgetter("sequence_length"))
        return (self.config["width"] - 2) / largest.sequence_length

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

    def generate_synthase_gradient(self, synthase):
        """Create a linearGradient SVG element representing the domain architecture of a
        Synthase.

        For example, given a Synthase:

        >>> synthase = Synthase(
        ...     header='synthase',
        ...     sequence='ACGACG...',  # length 100
        ...     domains=[Domain(type='KS', start=50, end=100)],
        ... )

        The generated gradient will be:

        ::

        <linearGradient id="synthase_doms" x1="0%" y1="0%" x2="100%" y2="0%">
        <stop offset="50%" stop-color="white"/>
        <stop offset="50%" stop-color="#08B208"/>
        <stop offset="100%" stop-color="#08B208"/>
        <stop offset="100%" stop-color="white"/>
        </linearGradient>

        Note that a generated linearGradient will have 4 stops for every Domain object
        in the `domains` attribute; this ensures that each Domain colour has a hard
        edge instead of blending into the next Domain.

        Also note that the `id` parameter of the generated linearGradient takes the
        form `header_doms`. This allows the Synthase polygon `fill` attribute to
        reference this linearGradient.

        Parameters
        ----------
        synthase : Synthase
            A Synthase oject. Must have a non-empty `sequence` attribute, as this is
            used to calculate the relative positioning of each domain.

        Returns
        -------
        str
            A string containing the linearGradient SVG element to use as the fill for
            the Synthase polygon.

        Raises
        ------
        ValueError
            If the supplied `synthase` has an empty `sequence` attribute.
        """
        if not synthase.sequence:
            raise ValueError("Synthase has no sequence")

        stops = []
        sequence_length = len(synthase.sequence)

        for domain in synthase.domains:
            start_pct = int(domain.start / sequence_length * 100)
            end_pct = int(domain.end / sequence_length * 100)
            colour = self.colours[domain.type]
            stops.append(
                f'<stop offset="{start_pct}%" stop-color="white"/>\n'
                f'<stop offset="{start_pct}%" stop-color="{colour}"/>\n'
                f'<stop offset="{end_pct}%" stop-color="{colour}"/>\n'
                f'<stop offset="{end_pct}%" stop-color="white"/>'
            )

        return (
            '<linearGradient id="{}_doms" x1="0%" y1="0%" x2="100%" y2="0%">\n{}\n'
            "</linearGradient>"
            "".format(synthase.header, "\n".join(stops))
        )

    def generate_synthase_polygon(self, synthase, scale_factor=1):
        """Build SVG representation of one synthase.

        Length is determined by the supplied scale factor. Then, (x, y) coordinates
        are calculated to represent each point in the synthase arrow. Finally, an SVG
        polygon feature is built with extra information above in a text feature, e.g.

        ::

            synthase, 100aa, KS-AT
            A----------B
            |           \\
            |            C
            |           /
            E----------D


        For example, given a Synthase:

        >>> synthase = Synthase(
        ...     header='synthase',
        ...     sequence='ACGACG...',  # length 100
        ...     domains=[Domain(type='KS', start=50, end=100)],
        ... )

        The generated polygon will be (equivalent to):

        >>> figure.generate_synthase_polygon(synthase)
        <text dominant-baseline="hanging" font-size="12">synthase, 100aa, KS</text>'
        <polygon
            id="synthase" points="0,10.8,90,10.8,100,17.8,90,24.8,0,24.8"
            fill="url(#synthase_doms)" stroke="black" stroke-width="1.5"
        />

        Note that the `fill` attribute takes the form `header_doms`; this is the id
        of the linearGradient for this Synthase.

        Parameters
        ----------
        synthase : Synthase
            A Synthase oject. Must have a non-empty `sequence` attribute, as this is
            used to calculate the coordinates in the polygon.

        scale_factor : int
            Scaling factor to multiply the Synthase sequence length by.

        Returns
        -------
        str
            <text> and <polygon> SVG features representing this Synthase. The fill for
            the polygon corresponds to the linearGradient element generated using
            Figure.generate_synthase_gradient().

        Raises
        ------
        ValueError
            If the supplied `synthase` has an empty `sequence` attribute.
        """
        if not synthase.sequence:
            raise ValueError("Synthase has no sequence")

        sequence_length = len(synthase.sequence)
        scaled_length = scale_factor * sequence_length
        info_fsize_scaled = self.config["info_fsize"] * 0.9
        bottom_y = info_fsize_scaled + self.config["arrow_height"]
        middle_y = info_fsize_scaled + self.config["arrow_height"] / 2

        ax, ay = 0, info_fsize_scaled
        bx, by = scaled_length - 10, info_fsize_scaled
        cx, cy = scaled_length, middle_y
        dx, dy = scaled_length - 10, bottom_y
        ex, ey = 0, bottom_y

        points = f"{ax},{ay},{bx},{by},{cx},{cy},{dx},{dy},{ex},{ey}"
        information = f"{synthase.header}, {sequence_length}aa, {synthase.architecture}"

        return (
            f'<text dominant-baseline="hanging" font-size="{self.config["info_fsize"]}"'
            f">{information}</text>"
            f'<polygon id="{synthase.header}" points="{points}"'
            f' fill="url(#{synthase.header}_doms)" stroke="black"'
            ' stroke-width="1.5"/>'
        )

    def build_polygon_block(self, subtype, synthases, scale_factor):
        """Generate the SVG for a block of Synthase objects of a specified subtype.

        See Figure.visualise() for description of other parameters.

        Parameters
        ----------
        subtype : str
            The subtype of the Synthase objects supplied to this method.

        synthases : list
            Synthase objects of a certain subtype to be visualised.

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
            f' font-size="{self.config["header_fsize"]}"'
            f' font-weight="bold">{subtype}</text>'
        )
        offset = self.config["header_fsize"]
        for synthase in synthases:
            polygon = self.generate_synthase_polygon(
                synthase, scale_factor=scale_factor
            )
            block += f'\n<g transform="translate(1,{offset})">\n{polygon}\n</g>'
            offset += (
                self.config["info_fsize"]
                + self.config["arrow_height"]
                + 4
                + self.config["arrow_spacing"]
            )
        return block, offset

    def visualise(self):
        """Construct the SVG figure.

        This function wraps all the necessary methods in the Figure class to generate
        the final SVG.

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
        scale_factor = self.calculate_scale_factor()

        blocks = ""
        offset = 3

        for subtype, synthases in self.iterate_synthase_types():
            block, height = self.build_polygon_block(subtype, synthases, scale_factor)
            blocks += f'<g transform="translate(0,{offset})">\n{block}\n</g>'
            offset += height + self.config["block_spacing"]

        gradients = [
            self.generate_synthase_gradient(synthase) for synthase in self.synthases
        ]

        return '<svg width="{}" height="{}">\n{}\n{}\n</svg>'.format(
            self.config["width"],
            offset - self.config["block_spacing"] - self.config["arrow_spacing"] - 4,
            "\n".join(gradients),
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
    def _local_cdsearch(cls, query_file, results_file=None, **kwargs):
        """Launch a new CD-Search run from a query file.
        """
        if not results_file:
            response = CDSearch(query_file=query_file, **kwargs)
            figure = cls(results.parse(response.text.split("\n")))
        else:
            with open(results_file) as _results:
                figure = cls(results.parse(_results))

        with open(query_file) as handle:
            figure.add_query_sequences(handle)

        return figure

    @classmethod
    def _remote_cdsearch(cls, query_ids, **kwargs):
        """Launch a new CD-Search run from a collection of query IDs."""
        sequences = efetch_sequences(query_ids)
        response = CDSearch(query_ids=query_ids, **kwargs)
        figure = cls(results.parse(response.text.split("\n")))
        figure.add_query_sequences(sequences=sequences)
        return figure

    @classmethod
    def from_cdsearch(
        cls, query_file=None, results_file=None, query_ids=None, **kwargs
    ):
        """Convenience function to directly instantiate a Figure from a new CDSearch job.

        Internally calls `_local_cdsearch` or `_remote_cdsearch` classmethods if
        query_file or query_ids is specified, respectively.

        For example:

        >>> figure = Figure.from_cdsearch(query_ids=['Q5BEJ6.1'])
        >>> figure
        HR-PKS
        ------
        Q5BEJ6.1    SAT-KS-AT-PT-ACP-MT-TR
        >>> figure['Q5BEJ6.1'].sequence
        'MTRASASGSGHEASTVFLFGPHVGTFTKASMDKLVRPLSQSPQRD...'

        Parameters
        ----------
        query_file : str
            Path to a FASTA file containing query sequences to be analysed.

        Returns
        -------
        Figure
            Figure built from the CDSearch query.
        """
        if query_file:
            return cls._local_cdsearch(query_file, results_file=results_file, **kwargs)
        elif query_ids:
            return cls._remote_cdsearch(query_ids, **kwargs)
        else:
            raise ValueError("Expected query_file or query_ids")

    def add_query_sequences(self, query_handle=None, sequences=None, ncbi=False):
        """Add sequences from query FASTA file to the Figure.

        Parameters
        ----------
        query_handle : open file handle
            Open file handle of a FASTA file containing sequences corresponding to the
            Synthases in this objects `synthases` attribute.

        sequences : dict
            A pre-populated dictionary containing sequences corresponding to the
            Synthases in this objects `synthases` attribute.

        Raises
        ------
        ValueError
            If no sequences could be obtained (i.e. `sequences` is empty).
        """
        if ncbi:
            sequences = efetch_sequences(self.headers)
        elif query_handle:
            sequences = fasta.parse(query_handle)

        if not sequences:
            raise ValueError("No sequences were obtained")

        for header, sequence in sequences.items():
            try:
                self[header].sequence = sequence
            except KeyError:
                LOG.error("Could not find match for %s, skipping", header)


def validate_colour(colour):
    """Check that a supplied colour is a valid hex code."""
    hex_regex = re.compile(r"^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$")
    if hex_regex.search(colour):
        return True
    return False
