#!/usr/bin/env python3

"""
Test suite for synthaser.py

Cameron Gilchrist
"""

from pathlib import Path

import pytest

from synthaser.models import Domain
from synthaser import results


TEST_DIR = Path(__file__).resolve().parent


def test_group_overlapping_hits():
    domains = [
        Domain(start=0, end=100),
        Domain(start=10, end=110),
        Domain(start=90, end=200),
        Domain(start=220, end=290),
        Domain(start=210, end=300),
    ]

    groups = [group for group in results.group_overlapping_hits(domains)]

    print(groups)
    assert groups == [domains[0:3], domains[3:][::-1]]


def test_parse_row():
    row = "\t\t\t0\t100\t0\t0\tsmart00825\tPKS_KS\t\t"
    domain = results.domain_from_row(row)
    assert domain.start == 0
    assert domain.end == 100
    assert domain.evalue == 0.0
    assert domain.type == "KS"
    assert domain.domain == "PKS_KS"


def test_parse_row_not_key_domain():
    row = "\t\t\t0\t100\t0\t\t\tTESTING\t\t"
    with pytest.raises(ValueError):
        results.domain_from_row(row)


@pytest.fixture
def domains():
    return [
        Domain(start=0, end=90, type="KS", domain="PKS_KS"),
        Domain(start=10, end=80, type="KS", domain="PKS"),
        Domain(start=100, end=200, type="AT", domain="PKS_AT"),
        Domain(start=130, end=190, type="AT", domain="Acyl_transf_1"),
    ]


@pytest.fixture
def test_results():
    return {
        "one": [
            Domain(
                type="KS",
                domain="PKS_KS",
                start=0,
                end=100,
                evalue=0.0,
                bitscore=500,
                accession="smart00825",
                pssm=214836,
                superfamily="cl09938",
            ),
            Domain(
                pssm=238428,
                accession="cd00832",
                superfamily="cl09938",
                type="KS",
                domain="CLF",
                start=0,
                end=100,
                evalue=5.11327e-141,
                bitscore=500,
            ),
        ],
        "two": [
            Domain(
                type="A",
                domain="A_NRPS",
                start=0,
                end=100,
                evalue=0.0,
                bitscore=500,
                pssm=341253,
                accession="cd05930",
                superfamily="cl17068",
            )
        ],
    }


def test_parse_cdsearch_table(tmp_path, test_results):
    results_file = tmp_path / "results.tsv"
    results_file.write_text(
        "Q#1 - >one\t\t\t0\t100\t0\t0\tsmart00825\tPKS_KS\t\t\n"
        "Q#1 - >one\t\t\t0\t100\t5.11327e-141\t0\tcd00832\tCLF\t\t\n"
        "Q#1 - >one\t\t\t0\t100\t5.11327e-141\t0\t\tNot_a_key_domain\t\t\n"
        "Q#2 - >two\t\t\t0\t100\t0\t0\tcd05930\tA_NRPS\t\t\n"
    )
    with results_file.open() as r:
        assert results.parse_cdsearch(r) == test_results


def test_filter_domains():
    domains = [
        Domain(
            type="C", domain="Condensation", start=0, end=100, evalue=0.0, bitscore=999
        ),
        Domain(
            type="C", domain="Condensation", start=10, end=100, evalue=1.0, bitscore=0
        ),
        Domain(
            type="E", domain="NRPS-para261", start=80, end=100, evalue=0.1, bitscore=999
        ),
        Domain(type="KS", domain="PKS", start=110, end=200, evalue=0.0, bitscore=999),
        Domain(
            type="KS", domain="PKS_KS", start=105, end=200, evalue=0.0, bitscore=999
        ),
    ]
    result = results.filter_domains(domains)
    expected = [
        Domain(
            type="E", domain="Condensation", start=0, end=100, evalue=0.0, bitscore=999
        ),
        Domain(
            type="KS", domain="PKS_KS", start=105, end=200, evalue=0.0, bitscore=999
        ),
    ]
    for a, b in zip(result, expected):
        assert a.domain == b.domain
        assert a.start == b.start
        assert a.end == b.end


def test_filter_results(test_results):
    assert results.filter_results(test_results) == {
        "one": [
            Domain(
                type="KS", domain="PKS_KS", start=0, end=100, evalue=0.0, bitscore=500
            )
        ],
        "two": [
            Domain(
                type="A", domain="A_NRPS", start=0, end=100, evalue=0.0, bitscore=500
            )
        ],
    }
