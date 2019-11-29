"""
Plot domain architectures using matplotlib.
"""


import math

from collections import namedtuple

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path

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
    index = 0
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
    for i, group in enumerate(groups):
        if i > 0:
            labels.append("")
        labels.extend([synthase.header for synthase in group])
    return labels


def annotate_groups(ax, groups):
    """Add synthase type annotation to axes."""
    for group in groups:
        # Grab the first and last Synthase in the group
        first, last = group[0], group[-1]

        # Calculate coordinates for forming the bar
        top_x, top_y = first.sequence_length + 60, first._index
        out_x = top_x + 50
        bot_y = last._index

        # MOVETO start coordinates, then LINETO for each subsequent point
        points = [(top_x, top_y), (out_x, top_y), (out_x, bot_y), (top_x, bot_y)]
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO]

        # Form the Patch and add it to the axes
        patch = PathPatch(Path(points, codes), facecolor="none", lw=1)
        ax.add_patch(patch)

        # Now, add the text next to bar
        mid_y = bot_y + (top_y - bot_y) / 2
        ax.text(out_x + 50, mid_y, first.subtype)


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


def draw_figure(ax, synthases):
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

    # Initialise the figure and adjust axes
    total = len(synthases) + len(groups) - 2
    label = generate_yticklabels(groups)
    ytick = np.arange(total + 1)

    # Adjust the yticks; all synthases are plotted as a unique index, with an empty
    # index used as inter-group spaces
    ax.set_yticks(ytick)
    ax.set_yticklabels(label)
    ax.set_ylim([total + 1, -1])

    # Round x limit up to nearest 500 for prettier scale
    ax.set_xlim([0, math.ceil(synthases[0].sequence_length / 500) * 500])

    # Adjust x tick formatting, hide y ticks
    ax.tick_params(axis="x", which="major", width=1, length=5, labelbottom=True)
    ax.tick_params(
        axis="y", which="both", width=0, length=0, labelsize="small", labelcolor="grey"
    )

    # Hide left, top and right borders
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Assign unique indices to each synthase
    assign_synthase_indices(groups)

    # Plot arrows
    for synthase in synthases:
        plot_gene_arrow(ax, synthase)

    # Annotate groups
    annotate_groups(ax, groups)


def draw_legend(ax, synthases):
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


def adjust_subplots(event):
    """Callback function to adjust subplot padding upon resize.

    First calls tight_layout to adjust side margins, then calculates top/bottom margins
    dynamically using title and legend heights, respectively.

    This is bound to `resize_event`, such that resizing the matplotlib figure viewer
    automatically scales the figure.
    """
    # Initial tight layout
    plt.tight_layout()

    # Get current figure/axes, make sure its drawn to enable window_extent
    figure, axes = plt.gcf(), plt.gca()
    figure.canvas.draw()

    # Adjust top/bottom spacing based on legend, title size
    legend = axes.get_legend().get_window_extent()
    figure.subplots_adjust(
        top=(event.height - 26) / event.height,
        bottom=(legend.y1 - legend.y0 + 40) / event.height,
    )


def calculate_initial_figsize(synthases):
    """Calculate a rough initial figure size based on the synthases.

    Default width = 7.
    For height, sum:
        (total synthases + total groups - 1) * 0.2
        total rows in legend * 0.6
        0.3

    Try to account for total y indices in the plot (synthases + inter-group spaces)
    as well as variable row number legends, and a bit of extra padding for the title.
    """
    synths = len(synthases)
    groups = len(set(s.subtype for s in synthases))
    dtypes = len(set(COLOURS[d.type] for s in synthases for d in s.domains))
    lgrows = math.ceil(dtypes / 8)
    return 7, (synths + groups - 1) * 0.2 + lgrows * 0.6 + 0.3


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

    rcParams["savefig.dpi"] = dpi

    # Rough guess at initial figure size
    width, height = calculate_initial_figsize(synthases)

    # Set up the figure and plot
    figure, axes = plt.subplots(figsize=(width, height))
    draw_figure(axes, synthases)
    draw_legend(axes, synthases)

    # Adjust subplots with mock event, bind function to resize event for dynamic adjustments
    Event = namedtuple("Event", ["width", "height"])
    event = Event(width * figure.dpi, height * figure.dpi)
    adjust_subplots(event)
    figure.canvas.mpl_connect("resize_event", adjust_subplots)

    if not file:
        plt.show()
    else:
        plt.savefig(file)
