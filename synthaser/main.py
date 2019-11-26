"""
CLI, main routine
"""

import argparse
import logging
import sys

from pathlib import Path

from synthaser import ncbi, __version__, plot, models

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

LOG = logging.getLogger("synthaser")
LOG.setLevel(logging.INFO)


def synthaser(
    query_file=None,
    query_ids=None,
    figure=None,
    svg_dpi=300,
    json_file=None,
    cdsid=None,
    db=None,
    smode=None,
    useid1=None,
    compbasedadj=None,
    filter=None,
    evalue=None,
    maxhit=None,
    dmode=None,
):
    """Run synthaser."""

    LOG.info("Starting synthaser")

    # Set flag to prevent re-writing the JSON we load from
    _json_loaded = False

    if json_file and Path(json_file).exists():
        LOG.info("Specified JSON file exists, attempting to read...")
        with open(json_file) as fp:
            synthases = models.SynthaseContainer.from_json(fp)
        _json_loaded = True
    else:
        ncbi.set_search_params(
            db=db,
            smode=smode,
            useid1=useid1,
            compbasedadj=compbasedadj,
            filter=filter,
            evalue=evalue,
            maxhit=maxhit,
            dmode=dmode,
        )

        synthases = ncbi.CDSearch(
            query_file=query_file, query_ids=query_ids, cdsid=cdsid
        )

    print(synthases, flush=True)

    if json_file and not _json_loaded:
        LOG.info("Serialising Figure to JSON: %s", json_file)
        with open(json_file, "w") as fp:
            synthases.to_json(fp)

    if figure:
        if figure == "show":
            # Opens matplotlib viewer
            LOG.info("Generating synthaser plot...")
            plot.plot(synthases)
        else:
            LOG.info("Writing SVG to: %s", figure)
            plot.plot(synthases, file=figure, dpi=svg_dpi)

    LOG.info("Finished synthaser")


def get_arguments(args):
    parser = argparse.ArgumentParser("synthaser")
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    inputs = parser.add_argument_group("Input")
    inputs.add_argument(
        "-qf",
        "--query_file",
        help="Path to FASTA file containing query synthase sequences",
    )
    inputs.add_argument(
        "-qi",
        "--query_ids",
        nargs="+",
        help="Collection of NCBI sequence identifiers corresponding to"
        " query synthases",
    )

    outputs = parser.add_argument_group("Output")
    outputs.add_argument(
        "-fig",
        "--figure",
        help="Generate a synthaser plot. If 'show' is specified,"
        " will open the matplotlib viewer to allow adjustments. Otherwise, will use"
        " the value provided as a file name for directly saving to SVG.",
    )
    outputs.add_argument(
        "-dpi", "--svg_dpi", help="DPI to use when saving figure (def. 300)"
    )
    outputs.add_argument(
        "-json",
        "--json_file",
        help="Serialise Synthase objects to JSON. If this is specified, the synthases"
        " can be loaded from this file using the synthaser Python API.",
    )

    search = parser.add_argument_group("CD-Search parameters")
    search.add_argument(
        "--db",
        choices=["cdd", "pfam", "smart", "tigrfam", "cog", "kog"],
        help="Name of the database to search (def. cdd)",
    )
    search.add_argument(
        "--cdsid",
        help="CD-Search run ID, e.g. QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0."
        " If provided, will attempt to retrieve results instead of starting a new"
        " search.",
    )
    search.add_argument(
        "--smode",
        choices=["auto", "prec", "live"],
        help="CD-Search search mode (def. auto)",
    )
    search.add_argument(
        "--useid1",
        choices=["true", "false"],
        help="Let NCBI find protein ID's in their archival database if they are not"
        " recognized as current Protein Entrez entries (def. true)",
    )
    search.add_argument(
        "--compbasedadj",
        choices=["0", "1"],
        help="Use composition-corrected scoring (def. 1)",
    )
    search.add_argument(
        "--filter",
        choices=["true", "false"],
        help="Filter out compositionally biased regions from query sequences (def. true)",
    )
    search.add_argument(
        "--evalue",
        type=float,
        help="Maximum E-value threshold (def. 3.0). Note that by default this is very"
        " generous, to catch domains that typically have low scores (e.g. 'SAT' domains)",
    )
    search.add_argument(
        "--maxhit",
        type=float,
        help="Maximum number of conserved domain hits to return from search (def. 500)",
    )
    search.add_argument(
        "--dmode",
        default="full",
        choices=["full", "rep", "std"],
        help="What level of CD-Search hits to report (def. full)",
    )

    args = parser.parse_args(args)

    if not any([args.query_ids, args.query_file, args.json_file]):
        raise ValueError("Expected query_ids, query_file or json_file")

    return args


def main():
    args = get_arguments(sys.argv[1:])

    synthaser(
        query_file=args.query_file,
        query_ids=args.query_ids,
        figure=args.figure,
        svg_dpi=args.svg_dpi,
        json_file=args.json_file,
        cdsid=args.cdsid,
        db=args.db,
        smode=args.smode,
        useid1=args.useid1,
        compbasedadj=args.compbasedadj,
        filter=args.filter,
        evalue=args.evalue,
        maxhit=args.maxhit,
        dmode=args.dmode,
    )


if __name__ == "__main__":
    main()
