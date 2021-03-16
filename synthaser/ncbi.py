#!/usr/bin/env python3

import time
import logging
import re

import requests

from Bio import Entrez, SeqIO

from synthaser import fasta

LOG = logging.getLogger(__name__)


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
CDSEARCH_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"

SEARCH_PARAMS = {
    "db": "cdd",
    "smode": "auto",
    "useid1": "true",
    "compbasedadj": "1",
    "filter": "true",
    "evalue": "3.0",
    "maxhit": "500",
    "dmode": "full",
    "tdata": "hits",
}


def launch(query):
    """Launches a new CDSearch run.

    Arguments:
        query (Synthase, SynthaseContainer):
            Synthase objects to be searched. This could either be a single Synthase
            object or a SynthaseContainer; other objects could be used as long as they
            implement a to_fasta method.
    Returns:
        cdsid (str):
            CDSearch ID (CDSID) corresponding to the new run. This takes the form:
            QM3-qcdsearch-XXXXXXXXXXXXXXXX-YYYYYYYYYYYYYYY.
    Raises:
        AttributeError: query has no to_fasta method
        AttributeError: No CDSID was returned from NCBI
    """
    try:
        files = {"queries": query.to_fasta()}
    except AttributeError:
        LOG.exception("Expected Synthase or SynthaseContainer")
        raise
    response = requests.post(CDSEARCH_URL, params=SEARCH_PARAMS, files=files)
    try:
        cdsid = re.search(r"#cdsid\t(.+?)\n", response.text).group(1)
    except AttributeError:
        LOG.exception("Search failed; no search ID returned")
        raise
    return cdsid


def check(cdsid):
    """Checks the status of a running CD-search job.

    CD-Search runs are assigned a unique search ID, which typically take the form::

        QM3-qcdsearch-xxxxxxxxxxx-yyyyyyyyyyy

    This function queries NCBI for the status of a running CD-Search job
    corresponding to the search ID specified by ``cdsid``.

    >>> response = check('QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0')

    If the job has finished, this function will return the ``requests.Response`` object
    which contains the run results. If the job is still running, this function will
    return None. If an error is encountered, a ValueError will be thrown with the
    corresponding error code and message.

    Arguments:
        cdsid (str): CD-search identifier (CDSID).
    Returns:
        True: If the job has completed and is ready for download
        False: If the job is still running
    Raises:
        ValueError:
            If the returned results file has a successful status code but is actually
            empty (i.e. contains no results), perhaps due to an invalid query.
        ValueError: When a status code of 1, 2, 4 or 5 is returned from the request.
    """
    response = requests.post(
        CDSEARCH_URL,
        params={"cdsid": cdsid, "tdata": "hits"}
    )
    match = re.search(r"#status\s+([\d])", response.text)
    if not match:
        raise ValueError("Failed to read search status")
    code = match.group(1)
    if code == "0":
        if response.text.endswith("Superfamily\n"):
            raise ValueError("Empty results file; perhaps invalid query?")
        return True
    if code == "3":
        return False
    errors = {
        "1": "Invalid search ID",
        "2": "No effective input (usually no query proteins or search ID specified)",
        "4": "Queue manager (qman) service error",
        "5": "Data is corrupted or no longer available (cache cleaned, etc)",
    }
    raise ValueError(f"Request failed; NCBI returned code {code} ({errors[code]})")


def get_results(cdsid):
    """Downloads results corresponding to a CDSID.

    Arguments:
        cdsid (str): CD-Search identifier
    Returns:
        requests.Response: Response object containing search results
    Raises:
        ValueError: If response has bad status code
    """
    params = {
        "tdata": "hits",
        "cddefl": "false",
        "qdefl": "false",
        "dmode": "full",
        "clonly": "false",
        "cdsid": cdsid,
    }
    response = requests.post(CDSEARCH_URL, params)
    if not response.ok:
        raise ValueError("Failed to retrieve results!")
    return response


def retrieve(cdsid, max_retries=-1, delay=20):
    """Poll CDSearch for results.

    This method queries the NCBI for results from a CDSearch job corresponding to
    the supplied cdsid. If max_retries is -1, this function will check for
    results every delay interval until something is returned.

    If you wish to save the results of a CD-Search run to file, you can supply an
    open file handle via the output parameter:

    >>> with open('results.tsv', 'w') as results:
    ...     retrieve(
    ...         'QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0',
    ...         output=results
    ...     )

    This function returns the Response object returned by check():

    >>> response = retrieve('QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0')
    >>> print(response.text)
    #Batch CD-search tool	NIH/NLM/NCBI
    #cdsid	QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0
    #datatype	hitsFull Results
    #status	0
    ...

    Arguments:
        cdsid (str):
            CD-search job ID. Looks like ``QM3-qcdsearch-xxxxxxxxxxx-yyyyyyyyyyy``.
        output (file pointer):
            Save results to a given open file handle instead of a local file. This
            facilitates usage of e.g. tempfile objects.
        max_retries (int):
            Maximum number of retries for checking job completion. If -1 is given, this
            function will keep paging for results until something is returned.
        delay (int):
            Number of seconds to wait between each request to the NCBI. The wait time is
            re-calculated to this value each time, based on the time taken by the
            previous request. By default, this is set to 20; giving a value less than 10
            will result in a ValueError being thrown.
    Returns:
        (requests.models.Response): Response returned by the check()
    Raises:
        ValueError: If delay is less than 10.
        ValueError: If no Response is returned by check()
    """
    if delay < 10:
        raise ValueError("Delay must be at least 10s")
    retries, previous = 0, 0
    finished = False
    while True:
        current = time.time()
        wait = previous + delay - current
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current
        LOG.info("Checking search status...")
        finished = check(cdsid)
        if finished:
            LOG.info("Search successfully completed!")
            break
        if max_retries > 0 and retries == max_retries:
            LOG.error("Maximum retry limit (%i) exceeded, breaking", max_retries)
            break
        retries += 1
    if not finished:
        raise ValueError("No results were returned")
    return get_results(cdsid)


def efetch_sequences(headers):
    """Retrieve protein sequences from NCBI for supplied accessions.

    This function uses EFetch from the NCBI E-utilities to retrieve the sequences for
    all synthases specified in headers. It then calls fasta.parse to parse the
    returned response; note that extra processing has to occur because the returned
    FASTA will contain a full sequence description in the header line after the
    accession.

    Arguments:
        headers (list): A collection of NCBI sequence identifiers (accession, GI, etc)
    Returns:
        sequences (dict): Sequences downloaded from NCBI
    """
    try:
        handle = Entrez.efetch(
            db="protein",
            id=headers,
            rettype="fasta",
            retmode="text",
        )
    except IOError:
        LOG.exception("Failed to fetch sequences")
        raise
    return fasta.parse(handle)

    
def set_search_params(
    database=None,
    smode=None,
    useid1=None,
    compbasedadj=None,
    filter=None,
    evalue=None,
    maxhit=None,
    dmode=None,
):
    """Set CD-Search search parameters.

    All search parameters are stored in SEARCH_PARAMS; this can either be edited
    directly, or through this function, prior to a search.

    Arguments:
        database (str):
            Name of search database. Available options are 'cdd' (default), 'pfam',
            'smart', 'tigrfam', 'cog' and 'kog'. Only applies when smode is live.
        smode (str):
            Search mode; 'auto' (automatic), 'prec' (precalculated only) or
            'live' (live searches).
        useid1 (str): Search archived sequences ('true' or 'false')
        compbasedadj (str): Composition-corrected scoring ('0' or '1')
        filter (str): Filter out compositionally biased regions ('true' or 'false')
        evalue (float): E-value cutoff
        maxhit (int): Maximum number of hits per query
        dmode (str): Data mode of output ('rep', 'std', or 'full')

    For a full description of parameters, refer to the NCBI's documentation_.

    .. _documentation: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBSearchParameters
    """
    params = {
        "db": database,
        "smode": smode,
        "useid1": useid1,
        "compbasedadj": compbasedadj,
        "filter": filter,
        "evalue": evalue,
        "maxhit": maxhit,
        "dmode": dmode,
    }
    for key, value in params.items():
        if not value:
            continue
        SEARCH_PARAMS[key] = value
