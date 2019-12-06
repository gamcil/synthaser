"""
This module contains routines for performing local/remote searches.
"""

import logging

from pathlib import Path

from synthaser import rpsblast, ncbi, fasta, results
from synthaser.models import SynthaseContainer
from synthaser.classify import classify


LOG = logging.getLogger(__name__)

SEARCH_HISTORY = []


def history():
    """Print out summary of previously saved CD-Search runs.

    The parameters and results of prior searches are stored in `SEARCH_HISTORY`, which
    is simply a list of dictionaries:

    >>> ncbi.SEARCH_HISTORY
    [{'cdsid': 'QM3-qcdsearch-C7558163144A288-147D9ECB39A09F8B',
    'parameters': {'db': 'cdd', 'smode': 'auto', 'useid1': 'true', 'compbasedadj': '1',
    'filter': 'true', 'evalue': '3.0', 'maxhit': '500', 'dmode': 'full', 'tdata':
    'hits'}, 'results': <Response [200]>}]

    This function simply formats the entries in `SEARCH_HISTORY` nicely:

    >>> ncbi.history()
    1.      Run ID: QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0
        Parameters:
                    db: cdd
                 smode: auto
                useid1: true
          compbasedadj: 1
                filter: true
                evalue: 3.0
                maxhit: 500
                 dmode: full
                 tdata: hits

    Raises
    ------
    ValueError
        If `SEARCH_HISTORY` is empty (i.e. no searches have been run).
    """
    if not SEARCH_HISTORY:
        raise ValueError("No searches have been run")

    for index, run in enumerate(SEARCH_HISTORY, 1):
        mode = run["mode"]
        params = "\n".join(
            f"{key}: {value}"
            for key, value in run.items()
            if key not in {"results", "mode"}
        )
        print(f"{index}. {mode}\n{params}")


def _container_from_query_file(handle):
    """Build SynthaseContainer from FASTA file handle."""
    return SynthaseContainer.from_sequences(fasta.parse(handle))


def _container_from_query_ids(ids):
    """Build SynthaseContainer from query ID file or collection.

    First checks if `ids` is an iterable; if so, fetch sequences from NCBI and return a
    new SynthaseContainer. Otherwise, expects a file containing a collection of IDs
    each on a new line.
    """
    if not hasattr(ids, "__iter__"):
        raise ValueError("Expected iterable")

    if len(ids) == 1 and Path(ids[0]).exists():
        # This is a file
        with Path(ids[0]).open() as fp:
            _ids = [line.strip() for line in fp]
        return SynthaseContainer.from_sequences(ncbi.efetch_sequences(_ids))

    # Otherwise, expect nargs with IDs
    return SynthaseContainer.from_sequences(ncbi.efetch_sequences(ids))


def prepare_input(query_ids=None, query_file=None):
    """Generate a SynthaseContainer from either query IDs or a query file.

    Returns
    -------
    models.SynthaseContainer
        Collection of `models.Synthase` objects representing query sequences.

    Raises
    ------
    ValueError
        Neither `query_ids` nor `query_file` provided
    ValueError
        Too many sequences were provided (NCBI limits searches at 4000 sequences)
    """
    if query_ids:
        container = _container_from_query_ids(query_ids)
    elif query_file:
        container = _container_from_query_file(query_file)
    else:
        raise ValueError("Expected 'query_ids' or 'query_file'")

    if len(container) > 4000:
        raise ValueError("Too many sequences (NCBI limit = 4000)")

    return container


def search(
    mode="remote",
    query_ids=None,
    query_file=None,
    domain_file=None,
    results_file=None,
    cdsid=None,
    delay=20,
    max_retries=-1,
    database=None,
    cpu=2,
):
    """..."""

    query = prepare_input(query_ids, query_file)

    if domain_file:
        LOG.info("Reading domain rules from: %s", domain_file.name)
        results.load_domain_json(domain_file)

    try:
        # If results_file is specified, first assume it's an actual results file
        with open(results_file) as rf:
            LOG.info("Reading results from: %s", results_file)
            for header, domains in results.parse(rf).items():
                query.get(header).domains = domains

    except (TypeError, FileNotFoundError):
        # Otherwise, user wants to start a search and save results under that name
        # OR just hasn't specified a results_file -> TypeError
        if mode == "remote":
            handle = _remote(
                query,
                output=results_file,
                cdsid=cdsid,
                delay=delay,
                max_retries=max_retries,
            )
        elif mode == "local":
            handle = _local(query, database, cpu=cpu, output=results_file)
        else:
            raise ValueError("Expected 'remote' or 'local'")

        LOG.info("Parsing results for domains...")
        for header, domains in results.parse(handle).items():
            query.get(header).domains = domains

    LOG.info("Classifying synthases...")
    for synthase in query:
        classify(synthase)

    return query


def _remote(query, cdsid=None, delay=20, max_retries=-1, output=None):
    """Launch a new CD-Search job.

    This function takes a query FASTA file or collection of NCBI sequence
    identifiers, sends them to the CD-Search server and then polls for response until
    one is returned.

    >>> results = CDSearch(query_file='path/to/query.fasta')
    >>> results = CDSearch(query_ids=['XP_000000001.1', ... ])
    >>> results
    <Response [200]>

    Search parameters are specified in `SEARCH_PARAMS`. For example, to adjust the
    e-value cutoff used in the search, simply change it there:

    >>> SEARCH_PARAMS['evalue'] = 2.0

    Then, future searches will use the updated value.
    """
    if not cdsid:
        LOG.info("Launching new CD-Search run")
        cdsid = ncbi.launch(query)

    LOG.info("Run ID: %s", cdsid)
    LOG.info("Polling NCBI for results...")
    response = ncbi.retrieve(cdsid, delay=delay, max_retries=max_retries)

    SEARCH_HISTORY.append(
        {"cdsid": cdsid, "query": query, "results": results, **ncbi.SEARCH_PARAMS}
    )

    if output:
        LOG.info("Writing CD-Search results table to %s", output)
        with open(output, "w") as out:
            out.write(response.text)

    return response.text.split("\n")


def _local(query, database, cpu=2, output=None, domain_file=None):
    """Run RPSBLAST."""
    LOG.info("Starting RPSBLAST")
    process = rpsblast.search(query.to_fasta().encode(), database, cpu)

    entry = {
        "mode": "rpsblast",
        "query": query,
        "database": database,
        "results": process.stdout,
    }

    SEARCH_HISTORY.append(entry)

    if output:
        LOG.info("Writing CD-Search results table to %s", output)
        with open(output, "wb") as handle:
            handle.write(process.stdout)

    return process.stdout.splitlines()
