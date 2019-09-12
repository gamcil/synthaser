#!/usr/bin/env python3

"""
Test suite for models.py
"""

from pathlib import Path

import pytest

from synthaser.models import (
    Synthase,
    Domain,
    assign_type,
    assign_subtype,
    hits_overlap,
    group_overlapping_hits,
)


TEST_DIR = Path(__file__).resolve().parent


@pytest.mark.parametrize(
    "types,result", [(["KS", "A"], "Hybrid"), (["KS"], "PKS"), (["A"], "NRPS")]
)
def test_assign_type(types, result):
    with pytest.raises(ValueError):
        assign_type([])
    assert assign_type(types) == result


@pytest.mark.parametrize(
    "stype,dtypes,result",
    [
        ("PKS", ["ER", "KR", "DH"], "HR-PKS"),
        ("PKS", ["KR", "DH"], "PR-PKS"),
        ("PKS", ["DH"], "PR-PKS"),
        ("PKS", ["KR"], "PR-PKS"),
        ("PKS", ["KS", "AT"], "NR-PKS"),
        ("PKS", ["KS"], "PKS-like"),
        ("PKS", ["AT"], "Other"),
        ("NRPS", ["A", "T", "C"], "NRPS"),
        ("NRPS", ["A", "T"], "NRPS-like"),
        ("NRPS", ["A"], "NRPS-like"),
        ("Test", [], "Test"),
    ],
)
def test_assign_subtype(stype, dtypes, result):
    assert assign_subtype(stype, dtypes) == result


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


def test_Domain_slice():
    sequence = "ACGCAGCAGTGCAGTGAGACGATGA"
    domain = Domain(start=10, end=20)
    assert domain.slice(sequence) == "TGCAGTGAGAC"


def test_Domain_serialisation(tmp_path):
    domain = Domain(start=0, end=100, type="KS", domain="PKS_KS")
    result_dict = {"start": 0, "end": 100, "type": "KS", "domain": "PKS_KS"}
    assert domain.to_dict() == result_dict

    json_file = tmp_path / "json"
    json_file.write_text(domain.to_json())

    with json_file.open() as js:
        from_json = Domain.from_json(js)

    assert from_json.start == 0
    assert from_json.end == 100
    assert from_json.type == "KS"
    assert from_json.domain == "PKS_KS"


@pytest.fixture
def domains():
    return [
        Domain(start=1, end=90, type="KS", domain="PKS_KS"),
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


def test_domain_repr(domains):
    assert repr(domains[0]) == "PKS_KS [KS] 1-90"


def test_domain_eq(domains):
    assert domains[0] == domains[0]
    assert domains[0] != domains[1]
    with pytest.raises(NotImplementedError):
        domains[0] == 1


def test_synthase_repr(synthase):
    assert repr(synthase) == "test\tKS-AT"


def test_synthase_eq(synthase):
    synthase2 = Synthase(header="test")
    assert synthase == synthase2

    synthase2.header = "test2"
    assert synthase != synthase2

    with pytest.raises(NotImplementedError):
        synthase == 1


def test_Synthase_serialisation(synthase, domains, tmp_path):
    """Test serialisation of a Synthase to dict and from/to JSON."""
    result_dict = {
        "header": "test",
        "sequence": "A" * 200,
        "domains": [
            {"start": 1, "end": 90, "type": "KS", "domain": "PKS_KS"},
            {"start": 100, "end": 200, "type": "AT", "domain": "PKS_AT"},
        ],
        "type": "PKS",
        "subtype": "NR-PKS",
    }

    assert synthase.to_dict() == result_dict

    json_file = tmp_path / "json"
    json_file.write_text(synthase.to_json())

    with json_file.open() as js:
        from_json = Synthase.from_json(js)

    assert from_json.header == "test"
    assert from_json.sequence == "A" * 200
    assert from_json.domains == [domains[0], domains[2]]
    assert from_json.type == "PKS"
    assert from_json.subtype == "NR-PKS"


def test_Synthase_extract_domains(synthase, domains):
    with pytest.raises(ValueError):
        synthase.domains = []
        synthase.extract_domains()

    synthase.sequence = "A" * 10 + "B" * 70 + "C" * 10 + "D" * 110
    synthase.domains = domains[:3]
    assert synthase.extract_domains() == {
        "KS": ["A" * 10 + "B" * 70 + "C" * 10, "A" + "B" * 70],
        "AT": ["D" * 101],
    }

    with pytest.raises(ValueError):
        synthase.sequence = ""
        synthase.extract_domains()


def test_Synthase_filter_overlapping_domains(domains, synthase):
    """Test domain filtering; function called in fixture function."""
    assert synthase.domains == [domains[0], domains[2]]


def test_Synthase_rename_nrps_domains(synthase):
    """Test computation of domain architecture based on Synthase type.
    For example, non-PKS synthases should replace ACP->T and TR->R; 'Hybrid' type
    should only replace after the condensation domain, 'NRPS' should replace all.
    """
    synthase.type = "Hybrid"
    nrps_module = [
        Domain(start=200, end=250, type="ACP", domain="PKS_PP"),
        Domain(start=250, end=350, type="C", domain="Condensation"),
        Domain(start=350, end=450, type="A", domain="A_NRPS"),
        Domain(start=450, end=550, type="ACP", domain="PKS_PP"),
        Domain(start=550, end=650, type="TR", domain="Thioester-redct"),
    ]
    synthase.domains.extend(nrps_module)
    assert synthase.architecture == "KS-AT-ACP-C-A-ACP-TR"
    synthase.rename_nrps_domains()
    assert synthase.architecture == "KS-AT-ACP-C-A-T-R"

    synthase.type = "NRPS"
    synthase.domains = nrps_module
    synthase.rename_nrps_domains()
    assert synthase.architecture == "T-C-A-T-R"


def test_Synthase_architecture(synthase):
    assert synthase.architecture == "KS-AT"
