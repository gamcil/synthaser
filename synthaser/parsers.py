"""Argument parsers."""


import argparse
import sys

from synthaser import __version__


def add_rules_group(parser):
    group = parser.add_argument_group("Other arguments")
    group.add_argument(
        "-df",
        "--domain_file",
        type=argparse.FileType(),
        help="JSON file containing alternate domain targets",
    )
    group.add_argument(
        "-cf",
        "--classify_file",
        type=argparse.FileType(),
        help="JSON file containing custom classification rules"
    )


def add_cdsearch_group(parser):
    group = parser.add_argument_group("CD-Search parameters")
    group.add_argument(
        "-rf",
        "--results_file",
        help="Path to results file from a previous search, or path to save results of"
        " a new search",
    )
    group.add_argument(
        "--cdsid",
        help="CD-Search run ID, e.g. QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0."
        " If provided, will attempt to retrieve results instead of starting a new"
        " search.",
    )
    group.add_argument(
        "--smode",
        choices=["auto", "prec", "live"],
        help="CD-Search search mode (def. auto)",
    )
    group.add_argument(
        "--useid1",
        choices=["true", "false"],
        help="Let NCBI find protein ID's in their archival database if they are not"
        " recognized as current Protein Entrez entries (def. true)",
    )
    group.add_argument(
        "--compbasedadj",
        choices=["0", "1"],
        help="Use composition-corrected scoring (def. 1)",
    )
    group.add_argument(
        "--filter",
        choices=["true", "false"],
        help="Filter out compositionally biased regions from query sequences (def. true)",
    )
    group.add_argument(
        "--evalue",
        type=float,
        help="Maximum E-value threshold (def. 3.0). Note that by default this is very"
        " generous, to catch domains that typically have low scores (e.g. 'SAT' domains)",
    )
    group.add_argument(
        "--maxhit",
        type=float,
        help="Maximum number of conserved domain hits to return from search (def. 500)",
    )
    group.add_argument(
        "--dmode",
        default="full",
        choices=["full", "rep", "std"],
        help="What level of CD-Search hits to report (def. full)",
    )


def add_searchopts_group(parser):
    group = parser.add_argument_group("Search options")
    group.add_argument(
        "-m",
        "--mode",
        choices=["local", "remote"],
        default="remote",
        help="Specifies synthaser search mode (def. 'remote'). If 'local', will run"
        " rpsblast against a local database.",
    )
    group.add_argument(
        "-db",
        "--database",
        help="Name of the database to search (def. cdd). If --mode is local, this should"
        " be the name of a valid rpsblast database",
    )


def add_output_group(parser):
    group = parser.add_argument_group("Output")
    group.add_argument(
        "-p",
        "--plot",
        nargs="?",
        const=True,
        default=False,
        help="Generate a synthaser plot"
    )
    group.add_argument(
        "-json",
        "--json_file",
        help="Serialise Synthase objects to JSON. If this is specified, the synthases"
        " can be loaded from this file using the synthaser Python API.",
    )
    group.add_argument(
        "-o",
        "--output",
        nargs="?",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Save domain architecture summary to file",
    )
    group.add_argument(
        "-lf",
        "--long_form",
        action="store_true",
        help="Return output in long data format",
    )


def add_input_group(parser):
    group = parser.add_argument_group("Input")
    group.add_argument(
        "-qf",
        "--query_file",
        type=argparse.FileType("r"),
        help="Path to FASTA file containing query synthase sequences",
    )
    group.add_argument(
        "-qi",
        "--query_ids",
        nargs="+",
        help="Collection of NCBI sequence identifiers corresponding to"
        " query synthases",
    )


def add_search_subparser(subparsers):
    parser = subparsers.add_parser(
        "search",
        help="Run a synthaser search",
        description="Run a synthaser search.",
        epilog="Cameron L.M. Gilchrist 2019",
    )
    add_input_group(parser)
    add_output_group(parser)
    add_searchopts_group(parser)
    add_cdsearch_group(parser)
    add_rules_group(parser)


def add_getseq_subparser(subparsers):
    p = subparsers.add_parser(
        "getseq",
        help="Download sequences from NCBI",
        description="Download sequences from NCBI in FASTA format. This utility will"
        " accept either a file containing newline separated sequence identifiers,"
        " or directly on the command line separated by spaces.",
        epilog="Cameron L.M. Gilchrist 2019",
    )
    p.add_argument(
        "sequence_ids",
        nargs="+",
        help="Collection of NCBI sequence identifiers to retrieve",
    )
    p.add_argument(
        "-o",
        "--output",
        nargs="?",
        default=sys.stdout,
        type=argparse.FileType("w"),
        help="Where to print output (def. stdout)",
    )


def add_getdb_subparser(subparsers):
    p = subparsers.add_parser(
        "getdb",
        help="Download a CDD database for local searches",
        description="Download a pre-formatted rpsblast database."
        " For full description of the available databases, see:"
        " https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CDSource",
        epilog="Cameron L.M. Gilchrist 2019",
    )
    p.add_argument(
        "database",
        choices=["Cdd", "Cdd_NCBI", "Cog", "Kog", "Pfam", "Prk", "Smart", "Tigr"],
        help="Database to be downloaded",
    )
    p.add_argument(
        "folder",
        help="Folder where database is to be saved. Will save a .tar.gz file, and"
        " extract its contents to a folder of the same name.",
    )


def get_parser():
    parser = argparse.ArgumentParser(
        "synthaser",
        epilog="Cameron L.M. Gilchrist 2020",
        description="synthaser: a Python toolkit for analysing domain architecture of"
        " secondary metabolite megasynth (et) ases with NCBI CD-Search.",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    subparsers = parser.add_subparsers(dest="subcommand")
    add_getdb_subparser(subparsers)
    add_getseq_subparser(subparsers)
    add_search_subparser(subparsers)
    return parser


def parse_args(args):
    parser = get_parser()
    args = parser.parse_args(args)

    if not args.subcommand:
        parser.print_help()
        parser.exit(1)

    if (
        args.subcommand == "search"
        and args.mode == "remote"
        and args.database not in (None, "cdd", "pfam", "smart", "tigrfam", "cog", "kog")
    ):
        raise ValueError("Expected 'cdd', 'pfam', 'smart', 'tigrfam', 'cog' or 'kog'")

    if args.subcommand == "search" and not any(
        [args.query_ids, args.query_file, args.json_file]
    ):
        raise ValueError("Expected query_ids, query_file or json_file")

    return args
