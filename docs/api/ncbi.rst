.. _ncbi_module:

:mod:`synthaser.ncbi`
-------------------------

This module handles all interaction with NCBI.

Given a collection of `Synthase` objects, a workflow might look like:

1. Launch new CD-Search run

>>> cdsid = ncbi.launch(synthases)
>>> cdsid
QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0

The query sequences are sent to the batch CD-Search API, where a new run is started and
assigned a unique CD-Search identifier (CDSID) which can be used to check on search
progress.

Note that search parameters for this search are specified by the values in `SEARCH_PARAMS`:

>>> ncbi.SEARCH_PARAMS
{
  'db': 'cdd',
  'smode': 'auto',
  'useid1': 'true',
  'compbasedadj': '1',
  'filter': 'true',
  'evalue': '3.0',
  'maxhit': '500',
  'dmode': 'full',
  'tdata': 'hits'
}

which can be freely edited, either directly or by using the `set_search_params`
function.

2. Poll CD-Search API for results using the CDSID

>>> response = ncbi.retrieve(cdsid)

This function repeatedly polls the API at regular intervals until either results or an
error has occurred. Internally, this function calls `check()`, which takes a CDSID and
sends a single request to the API. It returns a `Response` object (from the `requests`
library), which will have any search content saved in its `text` or `content`
properties.

>>> print(response.text)
#Batch CD-search tool	NIH/NLM/NCBI
#cdsid	QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0
#datatype	hitsFull Results
#status	0
#Start time	2019-09-03T04:21:23	Run time	0:00:04:23
#status	success

3. Parse results and create `Synthase` objects

>>> from synthaser import results
>>> handle = results.text.split("\n")
>>> synthases = results.parse(handle)


Additionally, this module provides `efetch_sequences`, a function for fetching sequences
from NCBI from a collection of accessions. For example:

>>> ncbi.efetch_sequences(['CBF71467.1', 'XP_681681.1'])
{'CBF71467.1': 'MQSAGMHRATA...', 'XP_681681.1': 'MQDLIAIVGSA...'}

The accessions are sent to the NCBI's Entrez API, which returns the sequences in FASTA
format. They are parsed using `fasta.parse`, and the resulting dictionary is returned.

.. automodule:: synthaser.ncbi
        :members:
