.. _search_module:

:mod:`synthaser.search`
-------------------------

This module serves as the starting point for `synthaser`, preparing input and
dispatching it to either local or remote searches.

In any given search, input can either be a FASTA file or a collection of NCBI sequence
identifiers. The `prepare_input` function is used to generate a `SynthaseContainer`
object from either source which can then be used as a query. For example:

>>> sc1 = search.prepare_input(query_ids=["SEQ001.1", "SEQ002.1"])
>>> sc2 = search.prepare_input(query_file="my_sequences.fasta")

If `query_ids` are used, the sequences are first retrieved using NCBI Entrez using
`ncbi.efetch_sequences()`.

A full `synthaser` search can be performed using the `search` function. This prepares
the input (ids or FASTA) as above, then launches local and remote searches using the
`ncbi` and `rpsblast` modules, respectively. Results are then parsed using the `results`
module, and classified using the `classify` module. Lastly, the `SynthaseContainer`
object which was created inside this function is returned.

>>> sc = search.search(query_file="my_sequences.fasta")

To use custom domain and classification rules, simply provide the paths to each file to
the `search` function:

>>> sc = search.search(
...     query_file="my_sequences.fasta",
...     domain_file="my_domains.json",
...     classify_file="my_rules.json",
... )

Previous searches are stored in the `SEARCH_HISTORY` variable, and can be summarised
using the `history` function:

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

.. automodule:: synthaser.search
        :members:
