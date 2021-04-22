#!/usr/bin/env python3

"""
Tests for fasta.py
"""

import pytest

from synthaser import fasta


def test_wrap():
    sequence = "GAGAACGTCGACGTCGATCGATCTAGCTGACAGCTAGCTA"
    wrapped = fasta.wrap(sequence, limit=10)
    assert wrapped == "GAGAACGTCG\nACGTCGATCG\nATCTAGCTGA\nCAGCTAGCTA"


@pytest.mark.parametrize(
    "header,sequence,wrap,result",
    [
        ("header", "AAAAABBBBBCCCCC", 5, ">header\nAAAAA\nBBBBB\nCCCCC"),
        (12345, "AAAAABBBBBCCCCC", 10, ">12345\nAAAAABBBBB\nCCCCC"),
    ],
)
def test_create(header, sequence, wrap, result):
    assert fasta.create(header, sequence, limit=wrap) == result


def test_parse(tmp_path):
    """Test synthaser.parse_fasta()."""
    fasta_file = tmp_path / "fasta.faa"
    fasta_file.write_text(
        ">TEST1\nACGTACGTACGTACGT\nACGTGCACGTGCACGTGCA\n"
        ">TEST2\nACGTGCAGTAGCTGAC\nAGACGGATTGCCAGTGACA\n\n"
        ">TEST3 test space\n\nACGTGCACGTGCA"
    )

    with fasta_file.open() as fasta_handle:
        parsed = fasta.parse(fasta_handle)

    assert parsed == {
        "TEST1": "ACGTACGTACGTACGTACGTGCACGTGCACGTGCA",
        "TEST2": "ACGTGCAGTAGCTGACAGACGGATTGCCAGTGACA",
        "TEST3": "ACGTGCACGTGCA",
    }


@pytest.mark.parametrize("limit", [0, 50, 39])
def test_count(tmp_path, limit):
    fasta_file = tmp_path / "fasta.faa"
    fasta_file.write_text(">\n" * limit)
    with fasta_file.open() as fasta_handle:
        assert fasta.count(fasta_handle) == limit
