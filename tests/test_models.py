#!/usr/bin/env python3

"""
Test suite for models.py
"""

import pytest

from synthaser.models import SynthaseContainer, Synthase, Domain


def test_Domain_slice():
    sequence = "ACGCAGCAGTGCAGTGAGACGATGA"
    domain = Domain(start=10, end=20)
    assert domain.slice(sequence) == "TGCAGTGAGAC"


def test_Domain_serialisation(tmp_path):
    domain = Domain(
        start=0,
        end=100,
        type="KS",
        domain="PKS_KS",
        evalue=0.01,
        bitscore=100,
        accession="smart00825",
        pssm=214836,
        superfamily="cl09938",
    )
    result_dict = {
        "start": 0,
        "end": 100,
        "type": "KS",
        "domain": "PKS_KS",
        "evalue": 0.01,
        "bitscore": 100,
        "accession": "smart00825",
        "pssm": 214836,
        "superfamily": "cl09938",
    }
    assert domain.to_dict() == result_dict

    json_file = tmp_path / "json"

    with json_file.open("w") as fp:
        domain.to_json(fp)

    with json_file.open() as js:
        from_json = Domain.from_json(js)

    for key, value in result_dict.items():
        assert getattr(from_json, key) == value


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
        classification=None
    )


def test_domain_str(domains):
    assert str(domains[0]) == "KS"


def test_domain_eq(domains):
    assert domains[0] == domains[0]
    assert domains[0] != domains[1]
    with pytest.raises(TypeError):
        domains[0] == 1


def test_synthase_str(synthase):
    assert str(synthase) == "test\tKS-AT"


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
            {
                "start": 1,
                "end": 90,
                "type": "KS",
                "domain": "PKS_KS",
                "evalue": None,
                "bitscore": None,
                "pssm": None,
                "superfamily": None,
                "accession": None,
            },
            {
                "start": 100,
                "end": 200,
                "type": "AT",
                "domain": "PKS_AT",
                "evalue": None,
                "bitscore": None,
                "pssm": None,
                "superfamily": None,
                "accession": None,
            },
        ],
        "classification": []
    }

    assert synthase.to_dict() == result_dict

    json_file = tmp_path / "json"
    json_file.write_text(synthase.to_json())

    with json_file.open() as js:
        from_json = Synthase.from_json(js)

    assert from_json.header == "test"
    assert from_json.sequence == "A" * 200
    assert from_json.domains == [domains[0], domains[2]]


def test_Synthase_extract_domains(synthase, domains):
    with pytest.raises(ValueError):
        synthase.domains = []
        synthase.extract_domains()

    synthase.sequence = "A" * 10 + "B" * 70 + "C" * 10 + "D" * 110
    synthase.domains = domains[:3]

    assert synthase.extract_domains() == {
        "KS": ["A" * 10 + "B" * 70 + "C" * 10, "A" + "B" * 70],
        "AT": ["D" * 101]
    }
    with pytest.raises(ValueError):
        synthase.sequence = ""
        synthase.extract_domains()


def test_Synthase_architecture(synthase):
    assert synthase.architecture == "KS-AT"


@pytest.fixture
def sc():
    return SynthaseContainer(
        [
            Synthase(
                header="one",
                sequence="AAAAABBBBB",
                domains=[Domain(type="KS", start=1, end=5)],
                classification=None,
            ),
            Synthase(
                header="two",
                sequence="BBBBBAAAAA",
                domains=[Domain(type="KS", start=6, end=10)],
                classification=None,
            ),
        ]
    )


def test_SynthaseContainer_extract_domains(sc):
    assert sc.extract_domains() == {'one': {'KS': ['AAAAA']}, 'two': {'KS': ['AAAAA']}}


def test_SynthaseContainer_add_typeerror(sc):
    with pytest.raises(TypeError):
        sc += 1


def test_SynthaseContainer_add(sc):
    one, two = sc[:1], sc[1:]
    assert one + two == sc


def test_SynthaseContainer_append_typerror(sc):
    with pytest.raises(TypeError):
        sc.append(1)


def test_SynthaseContainer_append(sc):
    new = Synthase(header="test", sequence="abc")
    sc.append(new)
    assert len(sc) == 3


def test_SynthaseContainer_extend_typeerror(sc):
    with pytest.raises(TypeError):
        sc.extend([1])


def test_SynthaseContainer_extend(sc):
    new = Synthase(header="test", sequence="abc")
    sc.extend([new])
    assert len(sc) == 3


def test_SynthaseContainer_get_keyerror(sc):
    with pytest.raises(KeyError):
        sc.get("fake")


def test_SynthaseContainer_get(sc):
    assert sc.get("two") == sc[1]


def test_SynthaseContainer_to_fasta(sc):
    assert sc.to_fasta() == ">one\nAAAAABBBBB\n>two\nBBBBBAAAAA"


def test_SynthaseContainer_from_sequences(sc):
    sequences = {"one": "AAAAABBBBB", "two": "BBBBBAAAAA"}
    sc2 = SynthaseContainer.from_sequences(sequences)
    assert (sc[0].header, sc[0].sequence) == (sc2[0].header, sc2[0].sequence)
    assert (sc[1].header, sc[1].sequence) == (sc2[1].header, sc2[1].sequence)


def test_SynthaseContainer_str(sc):
    print(sc)
    assert str(sc) == "one\tKS\ntwo\tKS"


def test_SynthaseContainer_json_serialisation(sc, tmp_path):
    d = tmp_path / "test.json"

    with d.open("w") as fp:
        sc.to_json(fp)

    with d.open() as fp:
        sc2 = SynthaseContainer.from_json(fp)

    assert sc == sc2
