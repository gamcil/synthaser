#!/usr/bin/env python3


"""
Tests for classify.py
"""


import pytest


from synthaser import classify
from synthaser.models import Domain, Synthase


@pytest.mark.parametrize(
    "domains,result",
    [
        ([Domain(type="KS"), Domain(type="KS")], "multi-modular PKS"),
        ([Domain(type="KS", domain="PKS")], "Type I PKS"),
        ([Domain(type="KS", domain="CLF")], "trans-AT PKS"),
        ([Domain(type="KS", domain="CHS_like")], "Type III PKS"),
    ],
)
def test_assign_PKS_type(domains, result):
    assert classify.assign_PKS_type(domains) == result


@pytest.mark.parametrize(
    "types,result",
    [
        (["ER", "KR", "DH"], "HR-PKS"),
        (["KR", "DH"], "PR-PKS"),
        (["KS", "AT", "DH"], "PR-PKS"),
        (["KS", "AT"], "NR-PKS"),
        (["KS"], "PKS-like"),
    ],
)
def test_assign_T1PKS_subtype(types, result):
    assert classify.assign_T1PKS_subtype(types) == result


@pytest.mark.parametrize(
    "types,result",
    [
        (["A", "T", "C"], "NRPS"),
        (["A", "ACP", "C"], "NRPS"),
        (["A", "T"], "NRPS-like"),
        (["T", "C"], "NRPS-like"),
        (["A"], "NRPS-like"),
    ],
)
def test_assign_NRPS_type(types, result):
    assert classify.assign_NRPS_type(types) == result


@pytest.mark.parametrize(
    "domains,result",
    [
        ([Domain(type="KS"), Domain(type="KS")], "multi-modular PKS"),
        ([Domain(type="KS"), Domain(type="A")], "Hybrid"),
        ([Domain(type="KS", domain="PKS")], "Type I PKS"),
        ([Domain(type="A"), Domain(type="T"), Domain(type="C")], "NRPS"),
        ([Domain(type="A")], "NRPS-like"),
    ],
)
def test_assign_broad_type(domains, result):
    assert classify.assign_broad_type(domains) == result


def test_assign_broad_type_no_domains():
    with pytest.raises(ValueError):
        classify.assign_broad_type([Domain(type="C")])
        classify.assign_broad_type([])


def test_classify_synthase():
    synthase = Synthase(domains=[Domain(type="KS", domain="PKS")])

    classify.classify(synthase)
    assert synthase.type == "Type I PKS"
    assert synthase.subtype == "PKS-like"

    synthase.domains.append(Domain(type="A"))
    classify.classify(synthase)
    assert synthase.type == "Hybrid"
    assert synthase.subtype == "Hybrid"


@pytest.mark.parametrize(
    "type,domains,result",
    [
        ("Hybrid", [Domain(t) for t in ("KS", "AT")], "KS-AT"),
        (
            "Hybrid",
            [Domain(t) for t in ("KS", "AT", "ACP", "C", "A", "ACP", "TR")],
            "KS-AT-ACP-C-A-T-R",
        ),
        ("NRPS", [Domain(t) for t in ("ACP", "C", "A", "ACP", "TR")], "T-C-A-T-R"),
    ],
)
def test_rename_NRPS_domains(type, domains, result):
    synthase = Synthase(type=type, domains=domains)
    classify.rename_NRPS_domains(synthase)
    assert synthase.architecture == result


def test_rename_NRPS_domains_bad_type():
    synthase = Synthase(type="Type I PKS")
    with pytest.raises(ValueError):
        classify.rename_NRPS_domains(synthase)
