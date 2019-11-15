"""
CLI, main routine
"""

import argparse
import logging
import sys

from synthaser import figure, ncbi

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

LOG = logging.getLogger("synthaser")
LOG.setLevel(logging.INFO)


def synthaser(
    query_file=None,
    query_ids=None,
    svg_file=None,
    json_file=None,
    db=None,
    smode=None,
    useid1=None,
    compbasedadj=None,
    filter=None,
    evalue=None,
    maxhit=None,
    dmode=None,
    arrow_height=None,
    arrow_spacing=None,
    block_spacing=None,
    header_fsize=None,
    info_fsize=None,
    width=None,
):
    """Run synthaser."""

    LOG.info("Starting synthaser")

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

    fig = figure.Figure.from_cdsearch(query_file=query_file, query_ids=query_ids)

    fig.set_config(
        arrow_height=arrow_height,
        arrow_spacing=arrow_spacing,
        block_spacing=block_spacing,
        header_fsize=header_fsize,
        info_fsize=info_fsize,
        width=width,
    )

    if json_file:
        LOG.info("Serialising Figure to JSON: %s", json_file)
        json_file.write(fig.to_json())

    if svg_file:
        LOG.info("Writing SVG to: %s", svg_file)
        svg_file.write(fig.visualise())

    print(fig)

    LOG.info("Finished synthaser")


def get_arguments():
    parser = argparse.ArgumentParser("synthaser")

    _inputs = parser.add_argument_group("Input")
    inputs = _inputs.add_mutually_exclusive_group(required=True)
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
        "-svg",
        "--svg_file",
        type=argparse.FileType("w"),
        help="Write SVG figure to this file name",
    )
    outputs.add_argument(
        "-json",
        "--json_file",
        type=argparse.FileType("w"),
        help="Serialise Figure object to JSON",
    )

    search = parser.add_argument_group("CD-Search parameters")
    search.add_argument(
        "--db",
        choices=["cdd", "pfam", "smart", "tigrfam", "cog", "kog"],
        help="Name of the database to search (def. cdd)",
    )
    search.add_argument(
        "--cdsid",
        help="ID of a CD-Search run, e.g. QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0",
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

    visual = parser.add_argument_group("Visualisation options")
    visual.add_argument(
        "--arrow_height", type=int, help="Height (px) of Synthase polygons"
    )
    visual.add_argument(
        "--arrow_spacing",
        type=int,
        help="Vertical spacing (px) to insert between each Synthase polygon",
    )
    visual.add_argument(
        "--block_spacing",
        type=int,
        help="Vertical spacing (px) to insert between each block of Synthase polygons",
    )
    visual.add_argument(
        "--header_fsize",
        type=int,
        help="Font size of Synthase type headers (e.g. Type I PKS)",
    )
    visual.add_argument(
        "--info_fsize",
        type=int,
        help="Font size of Synthase information headers (e.g. SEQ001, 2000aa,"
        " KS-AT-DH-ER-KR-ACP)",
    )
    visual.add_argument(
        "--width",
        type=int,
        help="Width (px) of generated SVG. The longest sequence in the Figure will"
        " be set to this width, and all other Synthases will be scaled accordingly",
    )

    return parser.parse_args()


def main():
    args = get_arguments(sys.argv[1:])

    synthaser(
        query_file=args.query_file,
        query_ids=args.query_ids,
        svg_file=args.svg_file,
        json_file=args.json_file,
        db=args.db,
        smode=args.smode,
        useid1=args.useid1,
        compbasedadj=args.compbasedadj,
        filter=args.filter,
        evalue=args.evalue,
        maxhit=args.maxhit,
        dmode=args.dmode,
        arrow_height=args.arrow_height,
        arrow_spacing=args.arrow_spacing,
        block_spacing=args.block_spacing,
        header_fsize=args.header_fsize,
        info_fsize=args.info_fsize,
        width=args.width,
    )


if __name__ == "__main__":
    main()
