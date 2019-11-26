"""
Plot domain architectures using matplotlib.
"""


import math

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.patches import Patch

plt.style.use("default")


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
    "E": "#FFA500",
    "cAT": "#FF3300",
}


def get_domain_information(synthase):
    """Generate arrays for plotting synthase domains.
    Return
    ------
    lengths : list
        Total domain sequence length (end - start).
    starts : list
        Domain start positions
    labels : list
        Domain types
    colour : list
        Domain colours, as defined in `COLOURS`
    """
    lengths, starts, labels, colour = [], [], [], []
    for domain in synthase.domains:
        lengths.append(domain.end - domain.start)
        starts.append(domain.start)
        labels.append(domain.type)
        colour.append(COLOURS[domain.type])
    return lengths, starts, labels, colour


def plot_gene_arrow(ax, synthase):
    """Add a gene arrow and domains to the axis."""

    # Plot the gene arrow background
    gene = ax.arrow(
        0,
        synthase._index,
        synthase.sequence_length,
        0,
        label=synthase.header,
        width=0.6,
        head_width=0.6,
        head_length=50,
        length_includes_head=True,
        zorder=1,
        facecolor=(1, 1, 1, 1),
        edgecolor=(0, 0, 0, 0),
    )

    # Get information for plotting domains as barh
    lengths, starts, labels, colours = get_domain_information(synthase)

    # Plot the domains in this synthase above the background
    doms = ax.barh(
        synthase._index,
        lengths,
        left=starts,
        label=labels,
        color=colours,
        height=1,
        zorder=2,
        edgecolor=(0, 0, 0, 0),
        linewidth=0,
    )

    # Domains are set larger than the width of the gene arrow, so clip them
    for dom in doms:
        dom.set_clip_path(gene)

    # Plot arrow again for borders
    ax.arrow(
        0,
        synthase._index,
        synthase.sequence_length,
        0,
        label=synthase.header,
        width=0.6,
        head_width=0.6,
        head_length=50,
        length_includes_head=True,
        zorder=3,
        facecolor=(1, 1, 1, 0),
        edgecolor=(0, 0, 0, 1),
        linestyle="-",
        linewidth=0.6,
    )


def subtype_iter(synthases):
    """Iterate subtype groups in a collection of Synthase objects."""
    _synthases = sorted(synthases, key=lambda s: s.subtype)
    _group, _subtype = "", []
    for synthase in _synthases:
        if synthase.subtype != _subtype:
            if _group:
                yield _group
            _group, _subtype = [synthase], synthase.subtype
        else:
            _group.append(synthase)
    yield _group


def assign_synthase_indices(groups):
    """Assign hidden y indices to synthases.

    The indices are staggered by synthase subtype group for visual separation in the
    plot.
    """
    index = 1
    for group in groups:
        for synthase in group:
            setattr(synthase, "_index", index)
            index += 1
        index += 1


def generate_yticklabels(groups):
    """Generate yticklabels from synthase groups.

    The 'groupings' in the plot are made by separating each block with an empty row (y
    value). Therefore, we assign the labels for each empty row as the subtype for the
    subsequent synthase group to match the indices of plotted genes.
    """
    labels = []
    for group in groups:
        labels.append(group[0].subtype)
        labels.extend([synthase.header for synthase in group])
    return labels


def format_yticklabels(ax, types):
    """Set different formatting for synthase types and headers in yticklabels."""
    for label in ax.get_yticklabels():
        if label.get_text() in types:
            label.set(color="black", fontsize="medium", weight="bold")
        else:
            label.set(color="grey", fontsize="small")


def generate_legend_elements(synthases):
    """Generate a list of unique Patch legend elements.

    First groups elements with common colours (e.g. R and TR, T and ACP), then generates
    `matplotlib.pyplot.Patch` instances with a comma-separated string of domain types as
    its label.
    """
    elements = {}
    for synthase in synthases:
        for domain in synthase.domains:
            colour = COLOURS[domain.type]
            if colour not in elements:
                elements[colour] = {domain.type}
            else:
                elements[colour].add(domain.type)

    return [
        Patch(label=", ".join(domains), facecolor=colour)
        for colour, domains in elements.items()
    ]


def generate_figure(ax, synthases):
    """Generate a synthaser plot on a given matplotlib axes."""

    # Set the title
    ax.set_title(
        f"Domain architectures of {len(synthases)} synth(et)ases", fontsize="large"
    )

    # Sort synthases, get subtype groups and domain types
    synthases = sorted(synthases, key=lambda s: s.sequence_length, reverse=True)
    groups = sorted(
        subtype_iter(synthases), key=lambda g: g[0].sequence_length, reverse=True
    )
    types = set(group[0].subtype for group in groups)

    # Initialise the figure and adjust axes
    # Sets initial height of the figure as 0.25 * total synthases + groups (inches).
    total = len(synthases) + len(groups) - 1
    label = generate_yticklabels(groups)
    ytick = np.arange(total + 1)

    ax.set_yticks(ytick)
    ax.set_yticklabels(label)
    ax.set_ylim([total + 1, -1])

    # Round x limit up to nearest 500
    ax.set_xlim([0, math.ceil(synthases[0].sequence_length / 500) * 500])

    # Adjust ticks
    ax.tick_params(axis="x", which="major", width=1, length=5, labelbottom=True)
    ax.tick_params(axis="y", which="both", width=0, length=0)

    # Hide left, top and right borders
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Format yticklabels
    format_yticklabels(ax, types)

    # Assign unique indices to each synthase
    assign_synthase_indices(groups)

    # Plot arrows
    for synthase in synthases:
        plot_gene_arrow(ax, synthase)


def generate_legend(ax, synthases):
    # Create legend elements for each domain actually in the synthases
    legend_elements = generate_legend_elements(synthases)

    # Make the legend
    legend = ax.legend(
        handles=legend_elements,
        fancybox=True,
        columnspacing=1,
        loc=8,
        bbox_to_anchor=(0.5, 0),
        bbox_transform=ax.figure.transFigure,
        ncol=8,
        fontsize="smaller",
    )

    return legend


def initial_figure_size(synthases, width=6):
    """Calculate initial figure height and subplot spacing."""

    # Number of legend elements is the # of unique domain colours
    legend_elements = len(
        set(
            COLOURS[domain.type]
            for synthase in synthases
            for domain in synthase.domains
        )
    )
    legend_nrows = math.ceil(legend_elements / 8)

    # Approximate size (inches) for plot elements
    synths = len(synthases)
    groups = len(set(s.subtype for s in synthases))
    arrows = 0.15 * (synths + groups)
    legend = 1.4 + 0.4 * (legend_nrows - 1)  # add extra for >1 row legends

    # Initial figure size based on elements
    return width, arrows + legend


def plot(synthases, file=None, dpi=300):
    """Plot a collection of synthases, show or save to file.

    Note that the format of the file saved is determined by the extension specified in
    the str given to `file`.
    For example, to save as png:

    >>> plot(synthases, file='figure.png', dpi=600)

    or as svg:

    >>> plot(synthases, file='figure.svg')

    The `dpi` argument is ignored if not saving to file, and effectively ignored by
    matplotlib if given with a file format that doesn't use it (e.g. svg).
    """

    width, height = initial_figure_size(synthases)

    figure, axes = plt.subplots(figsize=(width, height), tight_layout=True)

    generate_figure(axes, synthases)
    generate_legend(axes, synthases)

    plt.tight_layout()

    if not file:
        plt.show()
    else:
        plt.savefig(file, dpi=dpi, bbox_inches="tight")
