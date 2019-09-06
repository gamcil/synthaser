#!/usr/bin/env python3


"""
Module for handling web API calls to NCBI Batch CD-search.

1. Send fasta file to web API to create new query
2. Track job ID
3. Retrieve download for parsing

Cameron Gilchrist
"""


import json
import time
import logging
import re

from pathlib import Path

import requests


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


class CDSearch:
    """Run, check status and retrieve results of a batch CD-search run."""

    base_url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?"
    base_cfg = Path(__file__).resolve().parent / "cdsearch.json"

    def __init__(self, config_file=None):
        self.load_config(config_file=config_file)

    def load_config(self, config_file=None):
        """Load a configuration file.

        If `config_file` is not specified, this will load the default settings stored in
        'cdsearch.json'.
        """
        config_file = Path(config_file) if config_file else self.base_cfg

        if not config_file.exists():
            raise FileNotFoundError(f"No file exists at {config_file}")

        with config_file.open() as cfg:
            self.configure(json.load(cfg))

    def configure(self, config):
        """Validate and set parameters in configuration file.

        Parameters
        ----------
        config : dict
            Dictionary of synthaser and CD-search parameters.
            Schema should match (default cdsearch.json):
                {
                    "db": "cdd",
                    "smode": "auto",
                    "useid1": "true",
                    "compbasedadj": "0",
                    "filter": "false",
                    "evalue": "0.01",
                    "maxhit": "500"
                }
        """
        valid = {"db", "smode", "useid1", "compbasedadj", "filter", "evalue", "maxhit"}
        reserved = {"dmode", "tdata", "ainfmt", "qdefl", "cddefl"}
        for key, value in config.items():
            if key not in valid:
                if key in reserved:
                    raise ValueError(f"{key} is a reserved keyword")
                raise ValueError(f"{key} is not a valid keyword")

        self.config = config
        self.config["dmode"] = "full"
        self.config["tdata"] = "hits"

    def new_search(self, query_file=None, query_ids=None):
        """Run a new search against a supplied query file.

        Accepted input is either a list of NCBI sequence identifiers (accession, GI), or
        a FASTA formatted file containing amino acid sequences of query proteins.

        Parameters
        ----------
        query_file : str, optional
            Protein amino acid sequences in FASTA format.

        query_ids : list, optional
            Valid NCBI sequence identifiers (accession, GI, etc) of query proteins.

        Returns
        -------
        cdsid : str
            Job ID of the search. Used for checking status and retrieving results later.
        """
        if query_file and not query_ids:
            query_path = Path(query_file)

            if not query_path.exists():
                raise FileNotFoundError(f"No file exists at {query_file}")

            with query_path.open() as query_handle:
                response = requests.post(
                    self.base_url, files={"queries": query_handle}, data=self.config
                )
        elif query_ids and not query_file:
            if not isinstance(query_ids, (list, tuple)):
                raise ValueError("Query IDs must be given in a list or tuple")

            if not all(isinstance(x, (int, str)) for x in query_ids):
                raise ValueError("Expected a list/tuple of str")

            data = self.config.copy()
            data["queries"] = "\n".join(query_ids)  # line feed -> %0A
            response = requests.post(self.base_url, params=data)
        else:
            raise ValueError("A query must be specified with query_file OR query_ids")

        try:
            return re.search(r"#cdsid\t(.+?)\n", response.text).group(1)
        except AttributeError as exc:
            raise AttributeError("Search failed; no job ID returned") from exc

    def check_status(self, cdsid):
        """Check the status of a running CD-search job.

        Parameters
        ----------
        cdsid : int
            CD-search job ID.

        Returns
        -------
        requests.models.Response, None
            If the job has successfully completed (as indicated by status code 0), the
            Response object used to check will be returned, as this Response will also
            contain the job results.

            If the job is still running (status code 3), None will be returned.

        Raises
        ------
        ValueError when a status code of 1, 2, 4 or 5 is returned from the request.
        """
        response = requests.get(self.base_url, params={"cdsid": cdsid})
        code = re.search(r"#status\t([012345])[\n\t]", response.text).group(1)
        if code == "0":
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

    def retrieve_results(self, cdsid, output=None, check_interval=10, max_retries=20):
        """Retrieve the results of a CD-search job.

        Parameters
        ----------
        cdsid : str
            CD-search job ID. This looks like:
                QM3-qcdsearch-xxxxxxxxxxx-yyyyyyyyyyy

        output : open file handle
            Save results to a given open file handle instead of a local file. This
            facilitates usage of e.g. tempfile objects.

        check_interval : int
            Total seconds to wait before checking for job completion.

        max_retries : int
            Maximum number of retries for checking job completion.

        Raises
        -------
        ValueError
            If no Response model is returned by self.check_status()
        """
        log.info("Checking %s has finished", cdsid)
        retries = 0
        while retries <= max_retries:
            if retries > 0:
                log.info("Retry %i/%i", retries, max_retries)
            response = self.check_status(cdsid)
            if response:
                break
            retries += 1
            time.sleep(check_interval)

        if not response:
            raise ValueError("No results were returned")

        if output:
            log.info("Writing results to supplied file handle")
            output.write(response.text)

        return response

    def run(
        self,
        query_file=None,
        query_ids=None,
        output=None,
        force=False,
        check_interval=5,
        max_retries=10,
    ):
        """Convenience function to start a new batch CD-search job and retrieve results.

        Refer to self.new_search() and self.retrieve_results() for full description of
        parameters.
        """
        log.info("Starting new CD-search")
        cdsid = self.new_search(query_file=query_file, query_ids=query_ids)

        log.info("Retrieving results")
        log.debug(cdsid)
        self.retrieve_results(
            cdsid,
            check_interval=check_interval,
            max_retries=max_retries,
            output=output,
            force=force,
        )
