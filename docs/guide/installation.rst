.. _installation:

Installation
============

This section of the documentation covers the installation of `synthaser`.


Python version
--------------

``synthaser`` is written using Python 3, and should work with any version above 3.3.

Dependencies
------------

These packages are automatically installed when installing cblaster:

- requests_
- biopython_

Other dependencies
------------------

- RPS-BLAST is the search tool used in local cblaster searches
- rpsbproc is used to post-process RPS-BLAST results to remove redundant hits and
  fill in information about domain families like in the web CD-Search tool

Installation
------------

1. (Optional) Create a new virtual environment

.. code-block:: python

        python3 -m virtualenv venv
        source venv/bin/activate

This will create (and activate) a sandboxed environment where you can install
Python packages separately to those available on your system. This isn't necessarily
required, but is recommended.

2. Install ``synthaser``

The easiest way to obtain ``synthaser`` is to install it directly from PyPI using ``pip``:

.. code-block:: sh

        pip install synthaser

This will install ``synthaser``, as well as all of its required dependencies.
Alternatively, you could clone the cblaster repository from GitHub and
install it like so:

.. code-block:: sh

        git clone https://www.github.com/gamcil/synthaser
        cd synthaser
        pip install .

This will download the latest version of cblaster and install it from the downloaded
folder, rather than from PyPI.

``synthaser`` should now be available directly on your terminal:

::

        $ synthaser -h
        usage: synthaser [-h] [--version] {getdb,getseq,search} ...     
        synthaser: a Python toolkit for analysing domain architecture of secondary metabolite megasynth (et) ases with NCBI CD-Search.

        positional arguments:
          {getdb,getseq,search}
            getdb               Download a CDD database for local searches
            getseq              Download sequences from NCBI
            search              Run a synthaser search

        optional arguments:
          -h, --help            show this help message and exit
          --version             show program's version number and exit

        Cameron L.M. Gilchrist 2020


Installing RPS-BLAST and rpsbproc
---------------------------------

``RPS-BLAST`` is a distributed in the NCBI's BLAST+ toolkit. This can be acquired
either directly from NCBI's FTP_ or from your distributions repositories, for example
in Ubuntu: ``sudo apt install ncbi-blast+``.

To install ``rpsbproc``, follow these steps:

1. Acquire the relevant archive for your system from the `CDD FTP`__.

__ ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/

2. Extract the contents

3. Acquire the data files required by rpsbproc either by running
``utils/getcdddata.sh``, or directly from the FTP as detailed by the README
(see: domain-annotation files). The program does NOT require you to
download all of the domain databases. So, if doing the former, you can
chttps://www.circles.life/au/plan/ancel the run after the necessary files are in ``data/``, then delete ``db/``
and the database ``.tar.gz`` files.

4. Make sure the rpsbproc binary file is on your system $PATH.
This is a requirement of ``synthaser``, as it will throw an error if it
cannot find ``rpsbproc`` directly on the $PATH (i.e. accessable in terminal
just by typing 'rpsbproc').


.. _requests: https://requests.readthedocs.io/en/master/
.. _biopython: https://biopython.org/
.. _numpy: https://numpy.org/
.. _scipy: https://scipy.org/
.. _PySimpleGUI: https://pysimplegui.readthedocs.io/en/latest/
.. _genome2json: https://github.com/gamcil/genome2json
.. _diamond: https://github.com/bbuchfink/diamond
.. _FTP: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
