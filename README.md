# synthaser
[![Coverage Status](https://coveralls.io/repos/github/gamcil/synthaser/badge.svg?branch=master)](https://coveralls.io/github/gamcil/synthaser?branch=master&service=github)
![Tests passing](https://github.com/gamcil/synthaser/actions/workflows/python-app.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/synthaser/badge/?version=latest)](https://synthaser.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/synthaser.svg)](https://badge.fury.io/py/synthaser)

## Process
`synthaser` parses the results of a batch NCBI conserved domain search and determines
the domain architecture of secondary metabolite synthases.

## Installation
Install from PyPI using pip:

```sh
$ pip install --user synthaser
```

or clone the repo and install locally:

```sh
$ git clone https://www.github.com/gamcil/synthaser
$ cd synthaser
$ pip install .
```

Finally, configure synthaser with your e-mail address or NCBI API key (used when making requests to NCBI servers), for example:

```sh
$ synthaser config --email your@email.com
```

## Dependencies
`synthaser` is written in pure Python (3.6+), and requires only the following dependencies for
remote searches:
- `requests`, for interaction with the NCBI's CD-Search API
- `biopython`, for retrieving sequences from NCBI Entrez

If you want to do local searches, you'll need:
- `RPS-BLAST`, for performing local domain searches
- `rpsbproc`, for formatting RPS-BLAST results like CD-Search

These can be obtained from the [NCBI FTP](ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/).

## Usage
A full `synthaser` search can be performed as simply as:

```sh
$ synthaser search -qf sequences.fasta
```

Where `sequences.fasta` is a FASTA format file containing the protein sequences
that you would like to search.

For a full listing of available arguments, enter:

```sh
$ synthaser -h
```

### Visualising your results
`synthaser` is capable of generating fully-interactive, annotated visualisations
so you can easily explore your results. All that is required is one
extra argument:

```sh
$ synthaser search -qf sequences.fasta -p
```

This will generate a figure like so:

<img src="./img/anid_pks.png"
	width="400"
	alt="Example synthaser output">

[Click here](https://synthaser.readthedocs.io/en/latest/_static/anid.html) to play around with the full version of this example.

### Saving your search session
`synthaser` allows you to save your search results such that they can be easily
reloaded for further visualisation or exploration without having to fully re-do
the search.

To do this, use the `--json_file` command:

```sh
$ synthaser search -qf sequences.fasta --json_file sequences.json
```

This will save all of your results, in JSON format, to the file
`sequences.json`. Then, loading this session back into `synthaser`, is as easy
as:

```sh
$ synthaser search --json_file sequences.json ...
```

### Using your own rules
Though `synthaser` was originally designed to analyse secondary metabolite synthases,
it can easily be repurposed to analyse the domain architectures of any type of protein sequence.

Under the hood, `synthaser` uses a central rule file which contains:
1. Domain types, containing specific families to save in CD-Search results, corresponding to domain 'islands';
2. Rules for classifying the sequences based on domain architecture predictions; and
3. A hierarchy which determines the order of evaluation for the rules.

We distribute our fungal megasynthase rule file as the default, but providing your own rule file
is as simple as:

```sh
$ synthaser search -qf sequences.fasta --rule_file my_rules.json
```

We also provide a web application for assembling your own rule files, which can be
[found here](https://gamcil.github.io/synthaser/).

For a detailed explanation of how the rule file works, as well as API documentation,
please refer to the [documentation](https://synthaser.readthedocs.io/en/latest/).

## Citations
If you found `synthaser` helpful, please cite:

```text
Gilchrist, C. L., & Chooi, Y. H. (2021).
Synthaser: a CD-Search enabled Python toolkit for analysing domain architecture of fungal secondary metabolite megasynth (et) ases.
Fungal Biology and Biotechnology, 8(1), 1-19.
```
