"""
Plot domain architectures using matplotlib.

Very much an experiment, and not actually used by synthaser.
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
            label.set_color("black")
            label.set_fontsize("medium")
            label.set_weight("bold")
        else:
            label.set_color("grey")
            label.set_fontsize("small")


def plot(synthases):
    """Generate a synthaser plot."""

    synthases = sorted(synthases, key=lambda s: s.sequence_length, reverse=True)
    groups = sorted(
        subtype_iter(synthases), key=lambda g: g[0].sequence_length, reverse=True
    )
    types = set(group[0].subtype for group in groups)

    # Initialise the figure
    fig, ax = plt.subplots()

    # Adjust axes
    _total = len(synthases) + len(groups) - 1
    ax.set_yticks(np.arange(_total + 1))
    ax.set_yticklabels(generate_yticklabels(groups))
    ax.set_ylim([_total + 1, -1])
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

    # Add titles
    plt.suptitle("synthaser", fontsize=18, weight="bold")
    plt.title(f"Domain architectures of {len(synthases)} synth(et)ases", fontsize=15)

    # Create legend elements for each domain actually in the synthases
    legend_elements = [
        Patch(label=domain, facecolor=COLOURS[domain])
        for domain in set(
            domain.type for synthase in synthases for domain in synthase.domains
        )
    ]

    # Make the legend
    ax.legend(
        handles=legend_elements,
        fancybox=True,
        columnspacing=1,
        loc="best",
        ncol=2,
        fontsize="smaller",
    )

    plt.show()
