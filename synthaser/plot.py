"""
Plot domain architectures using matplotlib.

Very much an experiment, and not actually used by synthaser.
"""


import math

from matplotlib import pyplot as plt
from matplotlib.patches import Patch


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
    "cAT": "#000000",
}


def gene_arrow(ax, index, synthase):
    """Add a gene arrow and domains to the axis."""
    ax.plot(
        [0, len(synthase.sequence)],
        [index, index],
        label=synthase.header,
        color="#d3d3d3",
        linewidth=10,
    )
    for domain in synthase.domains:
        ax.plot(
            [domain.start, domain.end],
            [index, index],
            label=domain.type,
            color=COLOURS[domain.type],
            linewidth=10,
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


def setup_axes(ax, synthases):
    ax.set_title(synthases[0].subtype, fontsize=14)

    # Remove top, left and right borders
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_edgecolor(None)

    # Remove y-axis ticks
    ax.get_yaxis().set_ticks(list(range(len(synthases))))
    ax.get_yaxis().set_ticklabels([synthase.header for synthase in synthases])
    ax.tick_params(axis="both", which="both", length=0)


def plot(synthases):
    """Generate a synthaser plot."""

    synthases.sort(key=lambda s: len(s.sequence))
    domains = set(domain.type for synthase in synthases for domain in synthase.domains)
    subtypes = set(synthase.subtype for synthase in synthases)

    # Set up figure and axes
    fig, ax = plt.subplots(
        nrows=len(subtypes),
        sharex=True,
        gridspec_kw={
            "height_ratios": [len(group) for group in subtype_iter(synthases)]
        },
    )

    # Add titles
    plt.title("synthaser", fontsize=14)
    plt.suptitle(f"Domain architectures of {len(synthases)} synth(et)ases", fontsize=10)

    # Plot synthases
    for index, group in enumerate(subtype_iter(synthases)):
        setup_axes(ax[index], group)
        for y_pos, synthase in enumerate(group):
            gene_arrow(ax[index], y_pos, synthase)

    # Scale x-axis ticks to largest sequence
    lengths = [len(synthase.sequence) for synthase in synthases]

    # Adjust x-axis formatting
    ax[-1].spines["bottom"].set_visible(True)
    ax[-1].xaxis.set_ticks_position("bottom")
    ax[-1].tick_params(axis="x", which="major", width=1.00, length=5)
    ax[-1].set_xlim([0, math.ceil(max(lengths) / 500) * 500])

    # Create domain legend
    legend_elements = [
        Patch(label=domain, facecolor=COLOURS[domain]) for domain in domains
    ]
    fig.legend(
        handles=legend_elements, fancybox=True, columnspacing=1, loc="best", ncol=2
    )

    plt.show()
