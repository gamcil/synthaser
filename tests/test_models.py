#!/usr/bin/env python3

"""
Test suite for models.py
"""

import pytest

from synthaser.models import Synthase, Domain, extract_all_domains


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
    return Synthase(
        header="test",
        sequence="A" * 200,
        domains=[domains[0], domains[2]],
        type="PKS",
        subtype="NR-PKS",
    )


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


def test_Synthase_architecture(synthase):
    assert synthase.architecture == "KS-AT"


def test_extract_all_domains():
    one = Synthase(
        header="one", sequence="AAAAABBBBB", domains=[Domain(type="KS", start=1, end=5)]
    )
    two = Synthase(
        header="two",
        sequence="BBBBBAAAAA",
        domains=[Domain(type="KS", start=6, end=10)],
    )
    assert extract_all_domains([one, two]) == {
        "KS": [("one_KS_0", "AAAAA"), ("two_KS_0", "AAAAA")]
    }
