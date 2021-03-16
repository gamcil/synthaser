"""
CLI, main routine
"""

import logging
import sys

from pathlib import Path

from synthaser import (
    search,
    models,
    parsers,
    download,
)
from synthaser.plot import plot_synthases

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

LOG = logging.getLogger("synthaser")
LOG.setLevel(logging.INFO)


def synthaser(
    query_file=None,
    query_ids=None,
    plot=None,
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
    classify_file=None,
    results_file=None,
):
    """Run synthaser."""
    # Set flag to prevent re-writing the JSON we load from
    _json_loaded = False

    if json_file and Path(json_file).exists():
        LOG.info("Specified JSON file exists, attempting to read...")
        with open(json_file) as fp:
            synthases = models.SynthaseContainer.from_json(fp)
        _json_loaded = True
    else:
        try:
            synthases = search.search(
                mode=mode,
                query_file=query_file,
                query_ids=query_ids,
                cdsid=cdsid,
                domain_file=domain_file,
                classify_file=classify_file,
                results_file=results_file,
                database=database,
                smode=smode,
                useid1=useid1,
                compbasedadj=compbasedadj,
                filter=filter,
                evalue=evalue,
                maxhit=maxhit,
                dmode=dmode,
            )
        except ValueError:
            LOG.exception("Search failed! Exiting...")
            return
    if long_form:
        print(synthases.to_long(), flush=True, file=output)
    else:
        print(synthases, flush=True, file=output)

    if json_file and not _json_loaded:
        LOG.info("Serialising synthases to JSON: %s", json_file)
        with open(json_file, "w") as fp:
            synthases.to_json(fp)

    if plot:
        LOG.info("Generating synthaser plot...")
        plot = None if plot is True else plot
        plot_synthases(synthases, plot)


def main():
    LOG.info("Starting synthaser")

    args = parsers.parse_args(sys.argv[1:])

    if args.command == "getdb":
        download.getdb(args.database, args.folder)

    elif args.command == "getseq":
        container = search.prepare_input(query_ids=args.sequence_ids)
        print(container.to_fasta(), file=args.output)

    elif args.command == "genbank":
        from synthaser import genbank
        genbank.convert(
            args.genbank,
            antismash=args.antismash
        )

    elif args.command == "search":
        synthaser(
            query_file=args.query_file,
            query_ids=args.query_ids,
            plot=args.plot,
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
            classify_file=args.classify_file,
            results_file=args.results_file,
        )

    elif args.command == "extract":
        from synthaser import extract
        with open(args.session) as fp:
            synthases = models.SynthaseContainer.from_json(fp)
        extract.extract(
            synthases,
            args.prefix,
            types=args.types,
            classes=args.classes,
            families=args.families,
            mode=args.mode,
        )

    LOG.info("Done!")


if __name__ == "__main__":
    main()
