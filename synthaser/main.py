"""
CLI, main routine
"""

import argparse
import logging
import sys

from pathlib import Path

from synthaser import search, rpsblast, __version__, plot, models

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
    long_form=False,
    output=None,
    cdsid=None,
    mode="remote",
    database=None,
    smode=None,
    useid1=None,
    compbasedadj=None,
    filter=None,
    evalue=None,
    maxhit=None,
    dmode=None,
    domain_file=None,
    results_file=None,
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
        synthases = search.search(
            mode=mode,
            query_file=query_file,
            query_ids=query_ids,
            cdsid=cdsid,
            domain_file=domain_file,
            results_file=results_file,
            database=database,
            smode=smode,
            useid1=useid1,
            compbasedadj=compbasedadj,
            filter=filter,
            evalue=evalue,
            maxhit=maxhit,
            dmode=dmode,
            domain_file=domain_file,
            results_file=results_file,
        )

    if long_form:
        print(synthases.to_long(), flush=True, file=output)
    else:
        print(synthases, flush=True, file=output)

    if json_file and not _json_loaded:
        LOG.info("Serialising synthases to JSON: %s", json_file)
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
    parser = argparse.ArgumentParser(
        "synthaser",
        epilog="Cameron L.M. Gilchrist 2019",
        description="synthaser: a Python toolkit for analysing domain architecture of"
        " secondary metabolite megasynth (et) ases with NCBI CD-Search.",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    subparsers = parser.add_subparsers(dest="subcommand")

    getdb = subparsers.add_parser(
        "getdb",
        help="Download a CDD database for local searches",
        description="Download a pre-formatted rpsblast database."
        " For full description of the available databases, see:"
        " https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CDSource",
        epilog="Cameron L.M. Gilchrist 2019",
    )
    getdb.add_argument(
        "database",
        choices=["Cdd", "Cdd_NCBI", "Cog", "Kog", "Pfam", "Prk", "Smart", "Tigr"],
        help="Database to be downloaded",
    )
    getdb.add_argument(
        "folder",
        help="Folder where database is to be saved. Will save a .tar.gz file, and"
        " extract its contents to a folder of the same name.",
    )

    getseq = subparsers.add_parser(
        "getseq",
        help="Download sequences from NCBI",
        description="Download sequences from NCBI in FASTA format. This utility will"
        " accept either a file containing newline separated sequence identifiers,"
        " or directly on the command line separated by spaces.",
        epilog="Cameron L.M. Gilchrist 2019",
    )
    getseq.add_argument(
        "sequence_ids",
        nargs="+",
        help="Collection of NCBI sequence identifiers to retrieve",
    )
    getseq.add_argument(
        "-o",
        "--output",
        nargs="?",
        default=sys.stdout,
        type=argparse.FileType("w"),
        help="Where to print output (def. stdout)",
    )

    sub_search = subparsers.add_parser(
        "search",
        help="Run a synthaser search",
        description="Run a synthaser search.",
        epilog="Cameron L.M. Gilchrist 2019",
    )

    inputs = sub_search.add_argument_group("Input")
    inputs.add_argument(
        "-qf",
        "--query_file",
        type=argparse.FileType("r"),
        help="Path to FASTA file containing query synthase sequences",
    )
    inputs.add_argument(
        "-qi",
        "--query_ids",
        nargs="+",
        help="Collection of NCBI sequence identifiers corresponding to"
        " query synthases",
    )

    outputs = sub_search.add_argument_group("Output")
    outputs.add_argument(
        "-fig",
        "--figure",
        help="Generate a synthaser plot. If 'show' is specified,"
        " will open the matplotlib viewer to allow adjustments. Otherwise, will use"
        " the value provided as a file name for directly saving to SVG.",
    )
    outputs.add_argument(
        "-dpi",
        "--svg_dpi",
        help="DPI to use when saving figure (def. 300)",
        default=300,
    )
    outputs.add_argument(
        "-json",
        "--json_file",
        help="Serialise Synthase objects to JSON. If this is specified, the synthases"
        " can be loaded from this file using the synthaser Python API.",
    )
    outputs.add_argument(
        "-o",
        "--output",
        nargs="?",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Save domain architecture summary to file",
    )
    outputs.add_argument(
        "-lf",
        "--long_form",
        action="store_true",
        help="Return output in long data format",
    )

    grp_search = sub_search.add_argument_group("Search options")
    grp_search.add_argument(
        "--mode",
        choices=["local", "remote"],
        default="remote",
        help="Specifies synthaser search mode (def. 'remote'). If 'local', will run"
        " rpsblast against a local database.",
    )
    grp_search.add_argument(
        "-db",
        "--database",
        help="Name of the database to search (def. cdd). If --mode is local, this should"
        " be the name of a valid rpsblast database",
    )

    grp_cdsearch = sub_search.add_argument_group("CD-Search parameters")
    grp_cdsearch.add_argument(
        "-rf",
        "--results_file",
        help="Path to results file from a previous search, or path to save results of"
        " a new search",
    )
    grp_cdsearch.add_argument(
        "--cdsid",
        help="CD-Search run ID, e.g. QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0."
        " If provided, will attempt to retrieve results instead of starting a new"
        " search.",
    )
    grp_cdsearch.add_argument(
        "--smode",
        choices=["auto", "prec", "live"],
        help="CD-Search search mode (def. auto)",
    )
    grp_cdsearch.add_argument(
        "--useid1",
        choices=["true", "false"],
        help="Let NCBI find protein ID's in their archival database if they are not"
        " recognized as current Protein Entrez entries (def. true)",
    )
    grp_cdsearch.add_argument(
        "--compbasedadj",
        choices=["0", "1"],
        help="Use composition-corrected scoring (def. 1)",
    )
    grp_cdsearch.add_argument(
        "--filter",
        choices=["true", "false"],
        help="Filter out compositionally biased regions from query sequences (def. true)",
    )
    grp_cdsearch.add_argument(
        "--evalue",
        type=float,
        help="Maximum E-value threshold (def. 3.0). Note that by default this is very"
        " generous, to catch domains that typically have low scores (e.g. 'SAT' domains)",
    )
    grp_cdsearch.add_argument(
        "--maxhit",
        type=float,
        help="Maximum number of conserved domain hits to return from search (def. 500)",
    )
    grp_cdsearch.add_argument(
        "--dmode",
        default="full",
        choices=["full", "rep", "std"],
        help="What level of CD-Search hits to report (def. full)",
    )

    grp_other = sub_search.add_argument_group("Other arguments")
    grp_other.add_argument(
        "-df",
        "--domain_file",
        type=argparse.FileType(),
        help="JSON file containing alternate domain targets.",
    )

    args = parser.parse_args(args)

    if not args.subcommand:
        parser.print_help()
        sys.exit()

    if (
        args.subcommand == "search"
        and args.mode == "remote"
        and args.database not in ("cdd", "pfam", "smart", "tigrfam", "cog", "kog")
    ):
        raise ValueError("Expected 'cdd', 'pfam', 'smart', 'tigrfam', 'cog' or 'kog'")

    if args.subcommand == "search" and not any(
        [args.query_ids, args.query_file, args.json_file]
    ):
        raise ValueError("Expected query_ids, query_file or json_file")

    return args


def main():
    args = get_arguments(sys.argv[1:])

    if args.subcommand == "getdb":
        rpsblast.getdb(args.database, args.folder)

    elif args.subcommand == "getseq":
        container = search.prepare_input(query_ids=args.sequence_ids)
        print(container.to_fasta(), file=args.output)

    elif args.subcommand == "search":
        synthaser(
            query_file=args.query_file,
            query_ids=args.query_ids,
            figure=args.figure,
            svg_dpi=args.svg_dpi,
            json_file=args.json_file,
            output=args.output,
            long_form=args.long_form,
            cdsid=args.cdsid,
            mode=args.mode,
            database=args.database,
            smode=args.smode,
            useid1=args.useid1,
            compbasedadj=args.compbasedadj,
            filter=args.filter,
            evalue=args.evalue,
            maxhit=args.maxhit,
            dmode=args.dmode,
            domain_file=args.domain_file,
            results_file=args.results_file,
        )


if __name__ == "__main__":
    main()
