#!/usr/bin/env python3

"""
Test suite for synthaser.py

Cameron Gilchrist
"""

from pathlib import Path

import pytest

from synthaser.models import Synthase, Domain
from synthaser import results


TEST_DIR = Path(__file__).resolve().parent


@pytest.mark.parametrize(
    "a,b,threshold,result",
    [((0, 100), (10, 110), 0.9, True), ((0, 100), (80, 180), 0.9, False)],
)
def test_hits_overlap(a, b, threshold, result):
    dom_a = Domain(start=a[0], end=a[1])
    dom_b = Domain(start=b[0], end=b[1])
    assert results._hits_overlap(dom_a, dom_b, threshold=threshold) == result


def test_group_overlapping_hits():
    domains = [
        Domain(start=0, end=100),
        Domain(start=10, end=110),
        Domain(start=90, end=200),
        Domain(start=180, end=290),
        Domain(start=180, end=300),
    ]

    groups = [group for group in results._group_overlapping_hits(domains)]

    assert groups == [domains[0:2], [domains[2]], domains[3:]]


def test_parse_row():
    row = "\t\t\t0\t100\t0\t0\t\tPKS_KS\t\t"
    domain = results._domain_from_row(row)
    assert domain.start == 0
    assert domain.end == 100
    assert domain.evalue == 0.0
    assert domain.type == "KS"
    assert domain.domain == "PKS_KS"


def test_parse_row_not_key_domain():
    row = "\t\t\t0\t100\t0\t\t\tTESTING\t\t"
    with pytest.raises(ValueError):
        results._domain_from_row(row)


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
def test_results():
    return {
        "one": [
            Domain(type="KS", domain="PKS_KS", start=0, end=100, evalue=0.0),
            Domain(type="KS", domain="CLF", start=0, end=100, evalue=5.11327e-141),
        ],
        "two": [Domain(type="A", domain="A_NRPS", start=0, end=100, evalue=0.0)],
    }


def test_parse_cdsearch_table(tmp_path, test_results):
    results_file = tmp_path / "results.tsv"
    results_file.write_text(
        "Q#1 - >one\t\t\t0\t100\t0\t0\t\tPKS_KS\t\t\n"
        "Q#1 - >one\t\t\t0\t100\t5.11327e-141\t0\t\tCLF\t\t\n"
        "Q#1 - >one\t\t\t0\t100\t5.11327e-141\t0\t\tNot_a_key_domain\t\t\n"
        "Q#2 - >two\t\t\t0\t100\t0\t0\t\tA_NRPS\t\t\n"
    )
    with results_file.open() as r:
        assert results._parse_cdsearch_table(r) == test_results


def test_filter_domains():
    domains = [
        Domain(type="C", domain="Condensation", start=0, end=100, evalue=1.0),
        Domain(type="C", domain="Condensation", start=10, end=110, evalue=0.0),
        Domain(type="E", domain="NRPS-para261", start=80, end=100, evalue=0.1),
        Domain(type="KS", domain="PKS", start=90, end=200, evalue=0.0),
        Domain(type="KS", domain="PKS_KS", start=100, end=200, evalue=1.0),
    ]
    assert results._filter_domains(domains) == [
        Domain(type="E", domain="Condensation", start=10, end=110, evalue=0.0),
        Domain(type="KS", domain="PKS", start=90, end=200, evalue=0.0),
    ]


def test_filter_results(test_results):
    assert results._filter_results(test_results) == {
        "one": [Domain(type="KS", domain="PKS_KS", start=0, end=100, evalue=0.0)],
        "two": [Domain(type="A", domain="A_NRPS", start=0, end=100, evalue=0.0)],
    }
