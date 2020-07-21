.. _results_module:

:mod:`synthaser.results`
-------------------------

This module stores functions for parsing CD-Search output.

All functionality is provided by `parse`, which takes an open file handle
corresponding to a CD-Search hit table and returns a list of fully characterized
Synthase objects, i.e.:

>>> from synthaser import results
>>> with open('results.txt') as handle:
...     synthases = results.parse(handle)
>>> synthases
[AN6791.2 KS-AT-DH-MT-ER-KR-ACP, ... ]

`synthaser` uses the `results.DOMAINS` dictionary to control which domain families, and
the quality thresholds (length, bitscore) they must meet, to save in any given search.
This can be edited directly using `update_domains`, or loaded from a JSON file using
`load_domain_json`. An entry in this dictionary may look like:

.. code-block: python

        "PKS_KS": {
                "type": "KS"  # broader domain type
                "length": 298,  # CDD domain family profile length
                "bitscore": 241.079,  # threshold bitscore for 'specific hit'
        }

For further details on how to obtain these values and use a custom domain file,
please refer to the user guide.

.. automodule:: synthaser.results
        :members:
