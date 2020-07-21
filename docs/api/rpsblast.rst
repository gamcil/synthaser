.. _rpsblast_module:

:mod:`synthaser.rpsblast`
-------------------------

This module provides functionality for setting up and performing local `synthaser` searches using `RPS-BLAST` and `rpsbproc`. `RPS-BLAST` (Reverse PSI-BLAST) searches
query sequences against databases of domain family profiles, and `rpsbproc` is used to
post-process the raw results into something resembling results from an online CD-Search
run. If `synthaser` cannot find either program on the system `$PATH`, it will raise an
exception. For details on installing `RPS-BLAST` and `rpsbproc`, please refer to the user
guide.

A basic search can be performed using the `search` function:

>>> rpsblast.search("sequences.fasta", "Cdd_LE", cpu=4)

This will automatically search the sequences in ``sequences.fasta`` against the
``Cdd_LE`` using `RPS-BLAST` and process the raw results using `rpsbproc`, resulting in
`Response` object which be readily parsed like in a remote CD-Search.

A profile database can be downloaded using the `download_database` function, e.g.:

>>> path = rpsblast.download_database("my_folder", flavour="Cdd")

This will connect to the NCBI's FTP and download the "Cdd" database (the complete
database). The downloaded file will be a .tar archive, which can be extracted using
`untar`:

>>> untarred_path = rpsblast.untar(path)

Alternatively, just use `getdb` to do both steps at once:

>>> rpsblast.getdb("Cdd", "myfolder")

.. automodule:: synthaser.rpsblast
        :members:
