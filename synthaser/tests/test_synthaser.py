#!/usr/bin/env python3

"""
Test suite for synthaser.py

Cameron Gilchrist
"""

from pathlib import Path

import pytest

import synthaser


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
        fasta = synthaser.parse_fasta(fasta_handle)

    assert fasta == {
        "TEST1": "ACGTACGTACGTACGTACGTGCACGTGCACGTGCA",
        "TEST2": "ACGTGCAGTAGCTGACAGACGGATTGCCAGTGACA",
    }


def test_wrap_fasta():
    sequence = "GAGAACGTCGACGTCGATCGATCTAGCTGACAGCTAGCTA"
    wrapped = synthaser.wrap_fasta(sequence, limit=10)
    assert wrapped == "GAGAACGTCG\nACGTCGATCG\nATCTAGCTGA\nCAGCTAGCTA"


@pytest.mark.parametrize(
    "types,result", [(["KS", "A"], "Hybrid"), (["KS"], "PKS"), (["A"], "NRPS")]
)
def test_assign_type(types, result):
    synthase = synthaser.Synthase()

    with pytest.raises(ValueError):
        synthaser.assign_type(synthase)

    for typ in types:
        synthase.domains.append(synthaser.Domain(type=typ))

    assert synthaser.assign_type(synthase) == result


@pytest.mark.parametrize(
    "stype,dtypes,result",
    [
        ("PKS", ["ER", "KR", "DH"], "HR-PKS"),
        ("PKS", ["KR", "DH"], "PR-PKS"),
        ("PKS", ["DH"], "PR-PKS"),
        ("PKS", ["KR"], "PR-PKS"),
        ("PKS", ["KS", "AT"], "NR-PKS"),
        ("PKS", ["KS"], "Other"),
        ("PKS", ["AT"], "Other"),
        ("NRPS", ["A", "T", "C"], "NRPS"),
        ("NRPS", ["A", "T"], "NRPS-like"),
        ("NRPS", ["A"], "NRPS-like"),
        ("Test", [], "Test"),
    ],
)
def test_assign_subtype(stype, dtypes, result):
    synthase = synthaser.Synthase(type=stype)

    for t in dtypes:
        synthase.domains.append(synthaser.Domain(type=t))

    assert synthaser.assign_subtype(synthase) == result


@pytest.mark.parametrize(
    "a,b,threshold,result",
    [((0, 100), (10, 110), 0.9, True), ((0, 100), (80, 180), 0.9, False)],
)
def test_hits_overlap(a, b, threshold, result):
    dom_a = synthaser.Domain(start=a[0], end=a[1])
    dom_b = synthaser.Domain(start=b[0], end=b[1])
    assert synthaser.hits_overlap(dom_a, dom_b, threshold=threshold) == result


def test_group_overlapping_hits():
    domains = [
        synthaser.Domain(start=0, end=100),
        synthaser.Domain(start=10, end=110),
        synthaser.Domain(start=90, end=200),
        synthaser.Domain(start=180, end=290),
        synthaser.Domain(start=180, end=300),
    ]

    groups = [group for group in synthaser.group_overlapping_hits(domains)]

    assert groups == [domains[0:2], [domains[2]], domains[3:]]


def test_Domain_from_cdsearch_row():
    row = "\t\t\t0\t100\t\t\t\tPKS_KS\t\t"
    domain = synthaser.Domain.from_cdsearch_row(row)
    assert domain.start == 0
    assert domain.end == 100
    assert domain.type == "KS"
    assert domain.domain == "PKS_KS"


def test_Domain_slice():
    sequence = "ACGCAGCAGTGCAGTGAGACGATGA"
    domain = synthaser.Domain(start=10, end=20)
    assert domain.slice(sequence) == "TGCAGTGAGAC"


def test_Domain_serialisation(tmp_path):
    domain = synthaser.Domain(start=0, end=100, type="KS", domain="PKS_KS")
    result_dict = {"start": 0, "end": 100, "type": "KS", "domain": "PKS_KS"}
    assert domain.to_dict() == result_dict

    json_file = tmp_path / "json"
    json_file.write_text(domain.to_json())

    with json_file.open() as js:
        from_json = synthaser.Domain.from_json(js)

    assert from_json.start == 0
    assert from_json.end == 100
    assert from_json.type == "KS"
    assert from_json.domain == "PKS_KS"


@pytest.fixture
def domains():
    return [
        synthaser.Domain(start=0, end=90, type="KS", domain="PKS_KS"),
        synthaser.Domain(start=10, end=80, type="KS", domain="PKS"),
        synthaser.Domain(start=100, end=200, type="AT", domain="PKS_AT"),
        synthaser.Domain(start=130, end=190, type="AT", domain="Acyl_transf_1"),
    ]


@pytest.fixture
def synthase(domains):
    synth = synthaser.Synthase(
        header="test", sequence="A" * 200, domains=domains, type="PKS", subtype="NR-PKS"
    )
    synth.filter_overlapping_domains()
    return synth


def test_Synthase_serialisation(synthase, domains, tmp_path):
    """Test serialisation of a Synthase to dict and from/to JSON."""
    result_dict = {
        "header": "test",
        "sequence": "A" * 200,
        "domains": [
            {"start": 0, "end": 90, "type": "KS", "domain": "PKS_KS"},
            {"start": 100, "end": 200, "type": "AT", "domain": "PKS_AT"},
        ],
        "type": "PKS",
        "subtype": "NR-PKS",
    }

    assert synthase.to_dict() == result_dict

    json_file = tmp_path / "json"
    json_file.write_text(synthase.to_json())

    with json_file.open() as js:
        from_json = synthaser.Synthase.from_json(js)

    assert from_json.header == "test"
    assert from_json.sequence == "A" * 200
    assert from_json.domains == [domains[0], domains[2]]
    assert from_json.type == "PKS"
    assert from_json.subtype == "NR-PKS"


def test_Synthase_filter_overlapping_domains(domains, synthase):
    """Test domain filtering; function called in fixture function."""
    assert synthase.domains == [domains[0], domains[2]]


def test_Synthase_architecture(synthase):
    """Test computation of domain architecture based on Synthase type.
    For example, non-PKS synthases should replace ACP->T and TR->R; 'Hybrid' type
    should only replace after the condensation domain, 'NRPS' should replace all.
    """
    assert synthase.architecture == "KS-AT"

    synthase.type = "Hybrid"
    synthase.domains.extend(
        [
            synthaser.Domain(start=200, end=250, type="ACP", domain="PKS_PP"),
            synthaser.Domain(start=250, end=350, type="C", domain="Condensation"),
            synthaser.Domain(start=350, end=450, type="A", domain="A_NRPS"),
            synthaser.Domain(start=450, end=550, type="ACP", domain="PKS_PP"),
            synthaser.Domain(start=550, end=650, type="TR", domain="Thioester-redct"),
        ]
    )
    assert synthase.architecture == "KS-AT-ACP-C-A-T-R"

    synthase.type = "NRPS"
    synthase.domains = synthase.domains[-5:]
    assert synthase.architecture == "T-C-A-T-R"


def test_Synthase_gradient(synthase):
    """Test generation of Synthase gradient based on its Domains.

    Each domain should have 4 stops, the first and last having stop-color="white" to
    create strict edges, the middle two having the colour of the domain as assigned in
    synthaser.COLOURS. The offset attribute should be an integer indicating at what
    percentage of the Synthase length a Domain starts and stops.
    """
    gradient = (
        '<linearGradient id="test_doms" x1="0%" y1="0%" x2="100%" y2="0%">\n'
        f'<stop offset="0%" stop-color="white"/>\n'
        f'<stop offset="0%" stop-color="#08B208"/>\n'
        f'<stop offset="45%" stop-color="#08B208"/>\n'
        f'<stop offset="45%" stop-color="white"/>\n'
        f'<stop offset="50%" stop-color="white"/>\n'
        f'<stop offset="50%" stop-color="#DC0404"/>\n'
        f'<stop offset="100%" stop-color="#DC0404"/>\n'
        f'<stop offset="100%" stop-color="white"/>\n'
        "</linearGradient>"
    )
    assert synthase.gradient == gradient


def test_Synthase_polygon(synthase):
    """Test generation of Synthase polygon SVG element.

    This comprises of two SVG elements:
        <text>
            Holds Synthase information. Should have font-size corresponding to argument
            given to Synthase.polygon() and information string generated by
            Synthase.architecture property.

        <polygon>
            Represents the gene arrow. Defines 10 points (coordinate pairs for each of
            the 5 anchors of the polygon) calculated based on scale_factor, sequence
            length and arrow_height.
    """
    # length=200, scale_factor=1, arrow_height=14, info_fsize=12
    assert synthase.polygon(scale_factor=1) == (
        f'<text dominant-baseline="hanging" font-size="12">test, 200aa, KS-AT</text>'
        f'<polygon id="test" points="0,10.8,190,10.8,200,17.8,190,24.8,0,24.8"'
        f' fill="url(#test_doms)" stroke="black"'
        ' stroke-width="1.5"/>'
    )


@pytest.fixture
def figure():
    synthases = [
        synthaser.Synthase(
            header="hrpks1", sequence="A" * 1000, type="PKS", subtype="HR-PKS"
        ),
        synthaser.Synthase(
            header="hrpks2", sequence="A" * 1100, type="PKS", subtype="HR-PKS"
        ),
        synthaser.Synthase(
            header="nrpks", sequence="A" * 800, type="PKS", subtype="NR-PKS"
        ),
        synthaser.Synthase(
            header="hybrid", sequence="A" * 2000, type="Hybrid", subtype="Hybrid"
        ),
        synthaser.Synthase(
            header="nrps", sequence="A" * 1500, type="NRPS", subtype="NRPS"
        ),
    ]

    figure = synthaser.Figure(synthases=synthases)
    return figure


def test_Figure_serialisation(figure, tmp_path):
    """Test serialisation of Figure to/from JSON."""
    json_file = tmp_path / "json"
    json_file.write_text(figure.to_json())
    with json_file.open() as js:
        figure2 = synthaser.Figure.from_json(js)
    assert figure == figure2


def test_Figure_get(figure):
    with pytest.raises(KeyError):
        figure["not_in_figure"]
    figure["hrpks1"]


@pytest.mark.parametrize("width,result", [(1000, 0.499), (500, 0.249)])
def test_Figure_scale_factor(figure, width, result):
    assert figure.scale_factor(width) == result


def test_Figure_iterate_synthase_types(figure):
    """Test grouping of Synthases on this Figure by their subtypes.
    Have to place function call on RHS since Synthases are sorted in-place.
    """
    assert [
        ("Hybrid", [figure.synthases[3]]),
        ("NRPS", [figure.synthases[4]]),
        ("HR-PKS", [figure.synthases[1], figure.synthases[0]]),
        ("NR-PKS", [figure.synthases[2]]),
    ] == list(figure.iterate_synthase_types())


@pytest.fixture
def anid():
    """Returns a Figure from query FASTA and CDSearch results table in tests directory.
    """
    results = TEST_DIR / "anid.tsv"
    fasta = TEST_DIR / "anid.faa"
    with results.open() as res, fasta.open() as faa:
        figure = synthaser.Figure.from_cdsearch_results(res, query_handle=faa)
    return figure


def test_Figure_from_cdsearch_results(anid):
    """Test instantiation of a Figure from CDSearch results.
    Loads a Figure from JSON generated from the test FASTA/CDSearch results and compares
    to the `anid` Figure fixture.
    """
    test_figure_json = TEST_DIR / "anid.json"
    with test_figure_json.open() as js:
        test_figure = synthaser.Figure.from_json(js)
    assert str(anid) == str(test_figure)
    assert anid == test_figure


def test_Figure_visualise(anid):
    """Test generation of SVG from a Figure object.

    Reads in the precomputed SVG in the tests directory and compares to the output of
    Figure.visualise() using the anid fixture.
    """
    test_svg = TEST_DIR / "anid.svg"
    with test_svg.open() as handle:
        svg = handle.read()
    assert anid.visualise() == svg


def test_Figure_add_query_sequences(figure):
    """Test adding amino acid sequences to the Figure."""
    with pytest.raises(KeyError):
        figure.add_query_sequences(sequences={"not_in_figure": "ACGT"})
    figure.add_query_sequences(sequences={"hrpks1": "ACGT", "nrps": "ACGT"})
    assert figure["hrpks1"].sequence == "ACGT"
