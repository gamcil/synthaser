#!/usr/bin/env python3

"""
Test suite for synthaser.py

Cameron Gilchrist
"""

from pathlib import Path
from operator import attrgetter

import pytest

from synthaser.models import Synthase, Domain
from synthaser.figure import Figure
from synthaser.results import ResultParser, parse_fasta, parse_results


TEST_DIR = Path(__file__).resolve().parent


def test_parse_fasta(tmp_path):
    """Test synthaser.parse_fasta()."""
    fasta_file = tmp_path / "fasta.faa"
    fasta_str = (
        ">TEST1\nACGTACGTACGTACGT\nACGTGCACGTGCACGTGCA\n"
        ">TEST2\nACGTGCAGTAGCTGAC\nAGACGGATTGCCAGTGACA"
    )
    fasta_file.write_text(fasta_str)

    with fasta_file.open() as fasta_handle:
        fasta = parse_fasta(fasta_handle)

    assert fasta == {
        "TEST1": "ACGTACGTACGTACGTACGTGCACGTGCACGTGCA",
        "TEST2": "ACGTGCAGTAGCTGACAGACGGATTGCCAGTGACA",
    }


@pytest.fixture
def resultparser():
    return ResultParser()


def test_ResultParser_set_domains_input(resultparser):
    with pytest.raises(TypeError):
        resultparser.set_domains(["list"])
    with pytest.raises(ValueError):
        resultparser.set_domains({"invalid_key": 1})
    with pytest.raises(TypeError):
        resultparser.set_domains({"KS": "not list or tuple"})
    with pytest.raises(TypeError):
        resultparser.set_domains({"KS": ["asd", 2, 3]})


def test_ResultParser_set_domains(resultparser):
    # "KS": ["PKS_KS", "PKS"],
    # "AT": ["PKS_AT", "Acyl_transf_1"],
    domains = {
        "KS": ["PKS_KS", "PKS", "test_KS"],  # duplication
        "AT": ["test_AT"],  # purely new domain
    }
    resultparser.set_domains(domains)
    assert resultparser.domains["KS"] == ["PKS_KS", "PKS", "test_KS"]
    assert resultparser.domains["AT"] == ["PKS_AT", "Acyl_transf_1", "test_AT"]


def test_ResultParser_parse_row(resultparser):
    row = "\t\t\t0\t100\t\t\t\tPKS_KS\t\t"
    domain = resultparser.parse_row(row)
    assert domain.start == 0
    assert domain.end == 100
    assert domain.type == "KS"
    assert domain.domain == "PKS_KS"


def test_ResultParser_parse_row_not_key_domain(resultparser):
    row = "\t\t\t0\t100\t\t\t\tTESTING\t\t"
    with pytest.raises(ValueError):
        resultparser.parse_row(row)


@pytest.fixture
def domains():
    return [
        Domain(start=0, end=90, type="KS", domain="PKS_KS"),
        Domain(start=10, end=80, type="KS", domain="PKS"),
        Domain(start=100, end=200, type="AT", domain="PKS_AT"),
        Domain(start=130, end=190, type="AT", domain="Acyl_transf_1"),
    ]


@pytest.fixture
def synthase(domains):
    synth = Synthase(
        header="test", sequence="A" * 200, domains=domains, type="PKS", subtype="NR-PKS"
    )
    synth.filter_overlapping_domains()
    return synth


@pytest.fixture
def anid():
    """Returns a Figure from query FASTA and CDSearch results table in tests directory."""
    with open(TEST_DIR / "anid.json") as js:
        return Figure.from_json(js)


def test_ResultParser_parse(resultparser, anid):
    anid_tsv = TEST_DIR / "anid.tsv"
    with anid_tsv.open() as tsv:
        anid_synthases = resultparser.parse(tsv)
    anid_synthases.sort(key=attrgetter("header"))
    anid.synthases.sort(key=attrgetter("header"))
    assert anid_synthases == anid.synthases


def test_parse_results(anid):
    anid_tsv = str(TEST_DIR / "anid.tsv")
    anid_synthases = parse_results(anid_tsv)
    anid_synthases.sort(key=attrgetter("header"))
    anid.synthases.sort(key=attrgetter("header"))
    assert anid_synthases == anid.synthases
