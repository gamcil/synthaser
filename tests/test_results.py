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
from synthaser.results import (
    ResultParser,
    parse_results,
    hits_overlap,
    group_overlapping_hits,
)


TEST_DIR = Path(__file__).resolve().parent


@pytest.mark.parametrize(
    "a,b,threshold,result",
    [((0, 100), (10, 110), 0.9, True), ((0, 100), (80, 180), 0.9, False)],
)
def test_hits_overlap(a, b, threshold, result):
    dom_a = Domain(start=a[0], end=a[1])
    dom_b = Domain(start=b[0], end=b[1])
    assert hits_overlap(dom_a, dom_b, threshold=threshold) == result


def test_group_overlapping_hits():
    domains = [
        Domain(start=0, end=100),
        Domain(start=10, end=110),
        Domain(start=90, end=200),
        Domain(start=180, end=290),
        Domain(start=180, end=300),
    ]

    groups = [group for group in group_overlapping_hits(domains)]

    assert groups == [domains[0:2], [domains[2]], domains[3:]]


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
    domains = {
        "KS": ["PKS_KS", "PKS", "test_KS"],  # duplication
        "AT": ["test_AT"],  # purely new domain
    }
    resultparser.set_domains(domains)
    assert resultparser.domains["KS"] == [
        "PKS_KS",
        "PKS",
        "CLF",
        "KAS_I_II",
        "CHS_like",
        "KAS_III",
        "test_KS",
    ]
    assert resultparser.domains["AT"] == ["PKS_AT", "Acyl_transf_1", "test_AT"]


def test_ResultParser_parse_row(resultparser):
    row = "\t\t\t0\t100\t0\t\t\tPKS_KS\t\t"
    domain = resultparser.parse_row(row)
    assert domain.start == 0
    assert domain.end == 100
    assert domain.evalue == 0.0
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


@pytest.fixture
def results():
    return {
        "one": [
            Domain(type="KS", domain="PKS_KS", start=0, end=100, evalue=0.0),
            Domain(type="KS", domain="CLF", start=0, end=100, evalue=5.11327e-141),
        ],
        "two": [Domain(type="A", domain="A_NRPS", start=0, end=100, evalue=0.0)],
    }


def test_ResultParser_parse_table(resultparser, tmp_path, results):
    results_file = tmp_path / "results.tsv"
    results_file.write_text(
        "Q#1 - >one\t\t\t0\t100\t0\t\t\tPKS_KS\t\t\n"
        "Q#1 - >one\t\t\t0\t100\t5.11327e-141\t\t\tCLF\t\t\n"
        "Q#1 - >one\t\t\t0\t100\t5.11327e-141\t\t\tNot_a_key_domain\t\t\n"
        "Q#2 - >two\t\t\t0\t100\t0\t\t\tA_NRPS\t\t\n"
    )
    with results_file.open() as r:
        assert resultparser.parse_table(r) == results


def test_ResultParser_build_synthases(resultparser, results):
    sequences = {"one": "ACGT"}
    assert resultparser.build_synthases(results, sequences=sequences) == [
        Synthase(header="one", type="Type I PKS", subtype="PKS-like", sequence="ACGT"),
        Synthase(header="two", type="NRPS", subtype="NRPS-like"),
    ]


@pytest.mark.parametrize(
    "header,sequence,domains,result",
    [
        (
            "one",
            "ACGT",
            [
                Domain(type="KS", domain="PKS_KS", start=0, end=90, evalue=0.0),
                Domain(type="KS", domain="PKS", start=0, end=100, evalue=1.0),
                Domain(type="AT", domain="PKS_AT", start=100, end=200, evalue=0.0),
            ],
            Synthase(
                header="one",
                sequence="ACGT",
                type="Type I PKS",
                subtype="PKS-like",
                domains=[
                    Domain(type="KS", domain="PKS_KS", start=0, end=90, evalue=0.0),
                    Domain(type="AT", domain="PKS_AT", start=100, end=200, evalue=0.0),
                ],
            ),
        )
    ],
)
def test_ResultParser_new_synthase(resultparser, header, sequence, domains, result):
    assert resultparser.new_synthase(header, domains, sequence) == result


def test_ResultParser_filter_domains(resultparser):
    domains = [
        Domain(type="C", start=0, end=100, evalue=1.0),
        Domain(type="C", start=10, end=110, evalue=0.0),
        Domain(type="E", start=80, end=100, evalue=0.1),
        Domain(type="KS", start=90, end=200, evalue=0.0),
        Domain(type="KS", start=100, end=200, evalue=1.0),
    ]
    assert resultparser.filter_domains(domains) == [
        Domain(type="E", start=10, end=110, evalue=0.0),
        Domain(type="KS", start=90, end=200, evalue=0.0),
    ]


def test_ResultParser_apply_special_rules(resultparser):
    domains = [
        Domain(type="C", start=0, end=100, evalue=1.0),
        Domain(type="C", start=10, end=110, evalue=0.0),
        Domain(type="E", start=80, end=100, evalue=0.1),
    ]
    assert resultparser.apply_special_rules(domains[:2]) == Domain(
        type="C", start=10, end=110, evalue=0.0
    )
    assert resultparser.apply_special_rules(domains) == Domain(
        type="E", start=10, end=110, evalue=0.0
    )


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
