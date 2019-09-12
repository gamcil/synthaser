#!/usr/bin/env python3


import time
import logging
import re

import requests


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class CDSearch:
    """The CDSearch class handles web API calls to NCBI CD-Search.

    Basic usage
    -----------
    First, instantiate a new CDSearch object:

    >>> cd = CDSearch()

    Then, call ``CDSearch.search()`` on either a query FASTA file or a list of NCBI
    identifiers. The results are stored in a ``requests.Response`` object:

    >>> response = cd.search(query_file='/path/to/queries.fasta')
    >>> response = cd.search(query_ids=['XP_664395.1', 'XP_663604.1'])
    >>> print(response.text)
    #Batch CD-search tool	NIH/NLM/NCBI
    #cdsid	QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0
    #datatype	hitsFull Results
    #status	0
    ...

    An open file handle can be passed via ``output`` to save results to a file:

    >>> with open('results.tsv', 'w') as out:
    ...     response = cd.run(query_file='/path/to/queries.fasta', output=out)

    The results will be saved to ``results.tsv``, and the Response object returned by the
    NCBI API call is returned.

    Parameters
    ----------
    config_file : str, optional
        Path to a custom configuration JSON file. If this is not specified, the default
        config distributed with synthaser (`cdsearch.json`) will be loaded.
    """

    ncbi = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"

    parameters = {
        "db": "cdd",
        "smode": "auto",
        "useid1": "true",
        "compbasedadj": "0",
        "filter": "false",
        "evalue": "0.01",
        "maxhit": "500",
        "dmode": "full",
        "tdata": "hits",
    }

    def __init__(self, parameters=None):
        if parameters:
            self.set_search_parameters(parameters)

    def set_search_parameters(self, config):
        """Validate and set parameters in the ``config`` attribute.

        The supplied config dictionary should contain valid CD-search parameter values
        keyed on the parameter name. Valid keywords are::

            db, smode, useid1, compbasedadj, filter, evalue, maxhit

        For a full description of CD-Search parameters, see here_.

        To ensure that the returned results file is of the correct format for use in
        ``synthaser``, the parameters 'dmode', 'tdata', 'ainfmt', 'qdefl' and 'cddefl'
        are reserved. Including these in ``config`` will result in a ValueError being
        thrown.

        .. _here: https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBWebAPI_submit

        Parameters
        ----------
        config : dict
            Dictionary of CD-search parameters.

        Raises
        ------
        ValueError
            If a reserved or otherwise invalid keyword is encountered in the supplied
            `config` dictionary.
        """
        valid = {"db", "smode", "useid1", "compbasedadj", "filter", "evalue", "maxhit"}
        reserved = {"dmode", "tdata", "ainfmt", "qdefl", "cddefl"}
        for key, value in config.items():
            if key not in valid:
                if key in reserved:
                    raise ValueError(f"{key} is a reserved keyword")
                raise ValueError(f"{key} is not a valid keyword")
            self.parameters[key] = value

    def new_search(self, query_file=None, query_ids=None):
        """Launch a new CD-Search run using a FASTA file or list of NCBI identifiers.

        To analyse sequences in a local FASTA file, supply the path to the file via
        `query_file`:

        >>> cdsid = cd.new_search(query_file='/path/to/queries.fasta')

        Otherwise, supply a list of valid NCBI identifiers (accessions, GI numbers) via
        `query_ids`:

        >>> cdsid = cd.new_search(query_ids=['XP_664395.1', 'XP_663604.1'])

        Parameters
        ----------
        query_file : str, optional
            Protein amino acid sequences in FASTA format.

        query_ids : list, optional
            Valid NCBI sequence identifiers (accession, GI, etc) of query proteins.

        Returns
        -------
        cdsid : str
            Search ID assigned by NCBI. Used for checking status and retrieving results.

        Raises
        ------
        FileNotFoundError
            If no file exists at the path specified by ``query_file``. This is
            propagated from calling ``open(query_file)``.
        TypeError
            If ``query_ids`` is not a list or tuple.
        TypeError
            If values in ``query_ids`` are not of type ``int`` or ``str``.
        ValueError
            If neither ``query_ids`` or ``query_file`` are specified.
        AttributeError
            If Response returned by the POST request contains no search ID.
        """
        if query_file and not query_ids:
            with open(query_file) as query_handle:
                response = requests.post(
                    self.ncbi, files={"queries": query_handle}, data=self.parameters
                )

        elif query_ids and not query_file:
            if not isinstance(query_ids, (list, tuple)):
                raise TypeError("Query IDs must be given in a list or tuple")

            if not all(isinstance(x, (int, str)) for x in query_ids):
                raise TypeError("Expected a list/tuple of int/str")

            data = self.parameters.copy()
            data["queries"] = "\n".join(query_ids)  # line feed -> %0A
            response = requests.post(self.ncbi, params=data)
        else:
            raise ValueError("A query must be specified with query_file OR query_ids")

        try:
            return re.search(r"#cdsid\t(.+?)\n", response.text).group(1)
        except AttributeError as exc:
            raise AttributeError("Search failed; no search ID returned") from exc

    def check_status(self, cdsid):
        """Check the status of a running CD-search job.

        CD-Search runs are assigned a unique search ID, which typically take the form::

            QM3-qcdsearch-xxxxxxxxxxx-yyyyyyyyyyy

        This function queries NCBI for the status of a running CD-Search job
        corresponding to the search ID specified by ``cdsid``.

        >>> response = cd.check_status('QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0')

        If the job has finished, this function will return the ``requests.Response`` object
        which contains the run results. If the job is still running, this function will
        return None. If an error is encountered, a ValueError will be thrown with the
        corresponding error code and message.

        Parameters
        ----------
        cdsid : str
            CD-search job ID.

        Returns
        -------
        requests.models.Response
            If the job has successfully completed (Response contains status code 0).
            This is the object returned by ``requests.get()``; if the search has
            completed, this object will contain the search results in its ``text`` and
            ``contents`` attributes.
        None
            If the job is still running (Response contains status code 3).

        Raises
        ------
        ValueError
            If the returned results file has a successful status code but is actually
            empty (i.e. contains no results), perhaps due to an invalid query.
        ValueError
            When a status code of 1, 2, 4 or 5 is returned from the request.
        """
        response = requests.get(self.ncbi, params={"cdsid": cdsid})
        code = re.search(r"#status\t([012345])[\n\t]", response.text).group(1)
        if code == "0":
            if response.text.endswith("Superfamily\n"):
                raise ValueError("Empty results file; perhaps invalid query?")
            return response
        if code == "3":
            return None
        errors = {
            "1": "Invalid search ID",
            "2": "No effective input (usually no query proteins or search ID specified)",
            "4": "Queue manager (qman) service error",
            "5": "Data is corrupted or no longer available (cache cleaned, etc)",
        }
        raise ValueError(f"Request failed; NCBI returned code {code} ({errors[code]})")

    def retrieve_results(self, cdsid, output=None, max_retries=-1, delay=20):
        """Retrieve the results of a CD-search job.

        This method queries the NCBI for results from a CDSearch job corresponding to
        the supplied ``cdsid``. If ``max_retries`` is -1, this function will check for
        results every ``delay`` interval until something is returned.

        If you wish to save the results of a CD-Search run to file, you can supply an
        open file handle via the ``output`` parameter:

        >>> cd = CDSearch()
        >>> with open('results.tsv', 'w') as results:
        ...     cd.retrieve_results(
        ...         'QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0',
        ...         output=results
        ...     )

        This function returns the ``Response`` object returned by
        ``CDSearch.check_status()``:

        >>> cd = CDSearch()
        >>> response = cd.retrieve_results('QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0')
        >>> print(response.text)
        #Batch CD-search tool	NIH/NLM/NCBI
        #cdsid	QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0
        #datatype	hitsFull Results
        #status	0
        ...

        Parameters
        ----------
        cdsid : str
            CD-search job ID. Looks like ``QM3-qcdsearch-xxxxxxxxxxx-yyyyyyyyyyy``.

        output : open file handle
            Save results to a given open file handle instead of a local file. This
            facilitates usage of e.g. tempfile objects.

        max_retries : int
            Maximum number of retries for checking job completion. If -1 is given, this
            function will keep paging the NCBI for results until something is returned.

        delay : int
            Number of seconds to wait between each request to the NCBI. The wait time is
            re-calculated to this value each time, based on the time taken by the
            previous request. By default, this is set to 20; giving a value less than 10
            will result in a ValueError being thrown.

        Returns
        -------
        response : requests.Response
            Response returned by the ``requests.get()`` in ``CDSearch.check_status()``.

        Raises
        -------
        ValueError
            If ``delay`` is < 10.
        ValueError
            If no Response model is returned by self.check_status()
        """
        if delay < 10:
            raise ValueError("Delay must be at least 10s")
        retries, previous = 0, 0
        while True:
            current = time.time()
            wait = previous + delay - current
            if wait > 0:
                time.sleep(wait)
                previous = current + wait
            else:
                previous = current
            response = self.check_status(cdsid)
            if response or (max_retries > 0 and retries == max_retries):
                break
            retries += 1
        if not response:
            raise ValueError("No results were returned")
        if output:
            output.write(response.text)
        return response

    def search(
        self, query_file=None, query_ids=None, output=None, delay=20, max_retries=-1
    ):
        """Convenience function to start a new batch CD-search job and retrieve results.

        This function first launches a new CD-Search run via ``CDSearch.new_search()``,
        then attempts to retrieve the results from the resulting ``cdsid`` via
        ``CDSearch.retrieve_results()``.

        Refer to ``CDSearch.new_search()`` and ``CDSearch.retrieve_results()`` for full
        description of parameters.

        >>> with open('results.tsv', 'w') as results:
        ...     response = cd.search(query_ids=['AN6791.2', 'AN6000.2'], output=results)
        """
        cdsid = self.new_search(query_file=query_file, query_ids=query_ids)

        return self.retrieve_results(
            cdsid, delay=delay, max_retries=max_retries, output=output
        )
