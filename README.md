# synthaser
[![Build Status](https://travis-ci.org/gamcil/synthaser.svg?branch=master)](https://travis-ci.org/gamcil/synthaser)
[![Coverage Status](https://coveralls.io/repos/github/gamcil/synthaser/badge.svg?branch=master)](https://coveralls.io/github/gamcil/synthaser?branch=master&service=github)
[![Documentation Status](https://readthedocs.org/projects/synthaser/badge/?version=latest)](https://synthaser.readthedocs.io/en/latest/?badge=latest)

## Process
`synthaser` parses the results of a batch NCBI conserved domain search and determines
the domain architecture of secondary metabolite synthases.

## Installation
Install from PyPI via pip:
```sh
$ pip install synthaser
```

or clone the repo and install locally:
```sh
$ git clone https://www.github.com/gamcil/synthaser
$ cd synthaser
$ pip install -e .
```

## Dependencies
`synthaser` is written for Python 3.6+ and has been tested on Linux (Ubuntu 18.04) and
Windows (10). The only external Python dependency is `requests`, which is used for
querying the NCBI's APIs.

## Usage
A search can be launched as simply as:
```sh
$ synthaser -qi <accessions> OR synthaser -qf query.fasta
```

For example, performing a `synthaser` run on the cichorine PKS:
```sh
$ synthaser -qi CBF69451.1
[11:13:42] INFO - Starting synthaser
[11:13:44] INFO - Launching new CDSearch run on IDs: ['CBF69451.1']
[11:13:45] INFO - Run ID: QM3-qcdsearch-14C5BC063AA03DDE-15B11AB00918AED0
[11:13:45] INFO - Polling NCBI for results...
[11:13:45] INFO - Checking search status...
[11:14:05] INFO - Checking search status...
[11:14:06] INFO - Search successfully completed!
NR-PKS
------
CBF69451.1      SAT-KS-AT-PT-ACP-ACP-MT-TE
[11:14:06] INFO - Finished synthaser
```

`synthaser` can also produce an SVG representation of the query synthases. For example,
we could take the CDSID (CD-Search ID) of the previous run, and provide the `--svg` flag:

```sh
$ synthaser -qi CBF69451.1 \
    --cdsid QM3-qcdsearch-14C5BC063AA03DDE-15B11AB00918AED0 \
    --svg figure.svg
```

The generated figure is then saved in `figure.svg`, and looks like:

<img
  src="https://raw.githubusercontent.com/gamcil/synthaser/master/img/cichorine_svg.png"
  width="600"
>

`synthaser` can also start batch searches, either by providing more than one sequence in
a query FASTA file (`-qf`), or more than one NCBI accession (`-qi`).

For example, searching PKS sequences from *A. nidulans*:

```sh
$ synthaser -qf sequences.fasta --json nidulans.svg
```

Produces:

<img
  src="https://raw.githubusercontent.com/gamcil/synthaser/master/img/anid_pks.png"
  width="600"
>

Refer to `synthaser --help` for all tweakable parameters for generating the SVG.

## Citations
If you found `synthaser` helpful, please cite:

```sh
1. <pending>
```
