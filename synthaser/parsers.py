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
        epilog="Usage examples\n--------------\n"
        "Analyse sequences in a FASTA file and generate a plot:\n"
        "  $ synthaser search -qf sequences.fa -p\n\n"
        "Analyse sequences from the NCBI Protein database and save the search:\n"
        "  $ synthaser search -qi Q0CJ59.1 CAA39295.1 -json session.json\n\n"
        "Use custom domain and classification rule files:\n"
        "  $ synthaser search -qf sequences.fa \\\n"
        "      -df domain_rules.json \\\n"
        "      -cf classification_rules.json\n\n"
        "Download a CDD database and do a local search:\n"
        "  $ synthaser getdb Smart mydatabases\n"
        "  $ synthaser search -qf sequences.fa \\\n"
        "      -m local -db mydatabases/Smart_LE/Smart\n\n"
        "Cameron L.M. Gilchrist, 2020.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
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


def add_domain_file_argument(parser):
    parser.add_argument(
        "domain_file",
        type=argparse.FileType("r"),
        help="CDD family file generated by rulegen parse",
    )


def add_domain_rule_argument(parser):
    parser.add_argument(
        "rule_file",
        help="synthaser Domain rule file"
    )


def add_convert_subparser(subparsers):
    text = "Convert CDD family names to accessions and vice versa."
    parser = subparsers.add_parser(
        "convert",
        help=text,
        description=text,
        epilog="Cameron L.M. Gilchrist, 2020."
    )
    add_domain_file_argument(parser)
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument(
        "-n",
        "--names",
        nargs="+",
        help="Domain family names"
    )
    inputs.add_argument(
        "-a",
        "--accessions",
        nargs="+",
        help="Domain family accessions"
    )


def add_update_domain_subparser(subparsers):
    text = "Create domain rule files, add new rules or update existing rules."
    parser = subparsers.add_parser("update", help=text, description=text)
    add_domain_file_argument(parser)
    add_domain_rule_argument(parser)

    single = parser.add_argument_group("Add single rule")
    single.add_argument("--type", help="Domain type")
    single.add_argument("--name", help="Domain name")
    single.add_argument("--families", help="CDD family names", nargs="+")

    multi = parser.add_argument_group("Add multiple rules")
    multi.add_argument("--file", type=argparse.FileType("r"), help="Skeleton file")


def add_remove_rule_subparser(subparsers):
    parser = subparsers.add_parser(
        "remove",
        help="Remove domain rules, or specific CDD families from rules.",
        description="Remove domain rules, or specific CDD families from rules."
        " Either remove multiple entire rules, or specify a single rule and multiple"
        " CDD family accessions using --families."
    )
    add_domain_rule_argument(parser)
    parser.add_argument("--rules", default=[], nargs="+", help="Rule classes")
    parser.add_argument("--families", default=[], nargs="+", help="CDD family accessions")


def add_download_cdd_subparser(subparsers):
    text = "Download CDD family information from the NCBI."
    parser = subparsers.add_parser(
        "download",
        help=text,
        description=text,
    )
    parser.add_argument("domain_file", help="Output CDD family file")
    parser.add_argument("--folder", help="Output folder")
    parser.add_argument("--indent", default=2, help="JSON indent level (def. 2)")


def add_summary_subparser(subparsers):
    text = "Print out domain rules."
    parser = subparsers.add_parser(
        "summary",
        help=text,
        description=text + " If no rules or families are specified, this module"
        " will print a summary of the current rules.",
    )
    add_domain_rule_argument(parser)
    parser.add_argument("--rules", default=[], nargs="+", help="Rule classes")
    parser.add_argument("--families", default=[], nargs="+", help="CDD family accessions")


def add_domain_rules_subparser(subparsers):
    text = "Create or manipulate domain rule files"
    parser = subparsers.add_parser(
        "domains",
        help=text,
        description=text,
        epilog="Usage examples\n--------------\n"
        "Download and build CDD family database:\n"
        "  $ synthaser domains download cdd.json\n\n"
        "Create a new domain rule file from a skeleton:\n"
        "  $ synthaser domains update cdd.json myrules.json \\\n"
        "      --file skeleton.txt\n\n"
        "Add the PKS_KS domain family as a KS rule:\n"
        "  $ synthaser domains update cdd.json myrules.json \\\n"
        "      --type KS --name PKS_KS\n\n"
        "...and remove it:\n"
        "  $ synthaser domains remove myrules.json --families PKS_KS\n\n"
        "Get the CDD accession for PKS_KS (smart00825):\n"
        "  $ synthaser domains convert cdd.json --names PKS_KS\n\n"
        "Print a summary of a domain rule file:\n"
        "  $ synthaser domains summary myrules.json\n\n"
        "Cameron L.M. Gilchrist, 2020.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="subcommand")
    add_download_cdd_subparser(subparsers)
    add_update_domain_subparser(subparsers)
    add_remove_rule_subparser(subparsers)
    add_convert_subparser(subparsers)
    add_summary_subparser(subparsers)


def add_extract_subparser(subparsers):
    parser = subparsers.add_parser(
        "extract",
        help="Extract domain/synthase sequences from search results",
        description="Extract domain/synthase sequences from search results.",
        epilog="Usage examples\n--------------\n"
        "Extract KS, A and TE domain sequences:\n"
        "  $ synthaser extract session.json out_ --types KS A TE\n"
        "  Output: out_KS.faa out_A.faa out_TE.faa\n\n"
        "Extract NRPS and non-reducing PKS sequences:\n"
        "  $ synthaser extract session.json out_ \\\n"
        "      --mode synthase \\\n"
        "      --classes Non-reducing NRPS\n"
        "  Output: out_Non-reducing.faa out_NRPS.faa\n\n"
        "Extract PKS_KS domains (CDD) only from highly-reducing PKSs:\n"
        "  $ synthaser extract session.json out_ \\\n"
        "      --families PKS_KS \\\n"
        "      --classes Highly-reducing\n"
        "  Output: out_PKS_KS.faa\n\n"
        "Cameron L.M. Gilchrist, 2020.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("session", help="Synthaser session file")
    parser.add_argument("prefix", help="Output file prefix")
    parser.add_argument(
        "-m",
        "--mode",
        default="domain",
        choices=["domain", "synthase"],
        help="Extract domain sequences or whole synthases from a session file"
    )
    parser.add_argument("--types", nargs="+", help="Domain types")
    parser.add_argument("--classes", nargs="+", help="Sequence classifications")
    parser.add_argument("--families", nargs="+", help="CDD families")


def add_genbank_subparser(subparsers):
    parser = subparsers.add_parser(
        "genbank",
        help="Extract protein sequences from GenBank files for analysis",
        description="Extract protein sequences from GenBank files."
        " To extract PKS or NRPS sequences from antiSMASH GenBank files,"
        " use the --antismash option."
    )
    parser.add_argument("genbank", help="GenBank file")
    parser.add_argument("output", help="Output file name")
    parser.add_argument(
        "--antismash",
        action="store_true",
        help="Extract PKS/NRPS sequences from an antiSMASH file"
    )


def get_parser():
    parser = argparse.ArgumentParser(
        "synthaser",
        description="synthaser: a Python toolkit for analysing domain architecture of"
        " secondary metabolite megasynth (et) ases with NCBI CD-Search.",
        epilog="Cameron L.M. Gilchrist 2020",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    subparsers = parser.add_subparsers(dest="command")
    add_search_subparser(subparsers)
    add_domain_rules_subparser(subparsers)
    add_getdb_subparser(subparsers)
    add_getseq_subparser(subparsers)
    add_extract_subparser(subparsers)
    add_genbank_subparser(subparsers)
    return parser


def parse_args(args):
    parser = get_parser()
    args = parser.parse_args(args)

    if not args.command:
        parser.print_help()
        parser.exit(1)

    if (
        args.command == "search"
        and args.mode == "remote"
        and args.database not in (None, "cdd", "pfam", "smart", "tigrfam", "cog", "kog")
    ):
        raise ValueError("Expected 'cdd', 'pfam', 'smart', 'tigrfam', 'cog' or 'kog'")

    if args.command == "search" and not any(
        [args.query_ids, args.query_file, args.json_file]
    ):
        raise ValueError("Expected query_ids, query_file or json_file")

    return args
