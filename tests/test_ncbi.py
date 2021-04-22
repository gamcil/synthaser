#!/usr/bin/env python3

"""
Unit tests for cdsearch.py
"""

from pathlib import Path

import pytest
import requests
import requests_mock

from synthaser import ncbi
from synthaser.models import SynthaseContainer, Synthase

TEST_DIR = Path(__file__).resolve().parent


@pytest.fixture
def parameters():
    return {
        "db": "pfam",
        "smode": "prec",
        "useid1": "false",
        "compbasedadj": "1",
        "filter": "true",
        "evalue": "0.05",
        "maxhit": "500",
    }


def test_launch():
    query_file = TEST_DIR / "anid.faa"
    synthase = Synthase(header="one", sequence="abcdef")
    container = SynthaseContainer([synthase])
    with requests_mock.Mocker() as m:
        text = (
            "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
            "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
            "#datatype\thitsFull Results\n#status\t3\tmsg\tJob is still running\n'"
        )
        m.post(ncbi.CDSEARCH_URL, text=text)
        cdsid = "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8"
        assert ncbi.launch(synthase) == cdsid
        assert ncbi.launch(container) == cdsid


def test_CDSearch_check_status_empty_file():
    cdsid = "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8"
    text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
        "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
        "#datatype\thitsFull Results\n#status\t0\n"
        "#Start time\t2019-09-03T04:21:23\tRun time\t0:00:04:23\n"
        "#status\tsuccess\n\n"
        "Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore"
        "\tAccession\tShort\tname\tIncomplete\tSuperfamily\n"
    )
    with requests_mock.Mocker() as mock, pytest.raises(ValueError):
        mock.post(ncbi.CDSEARCH_URL, text=text)
        ncbi.check(cdsid)


def test_CDSearch_check_status():
    # mock responses
    # test 0 -> returns response, 3 -> returns None, ValueError otherwise
    cdsid = "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8"
    text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
        "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
        "#datatype\thitsFull Results\n#status\t{}\tmsg\tJob is still running\n'"
    )
    with requests_mock.Mocker() as m:
        # Returns True if job has completed
        m.post(ncbi.CDSEARCH_URL, text=text.format("0"))
        response = ncbi.check(cdsid)
        assert response is True

        # Returns False if job is still running
        m.post(ncbi.CDSEARCH_URL, text=text.format("3"))
        response = ncbi.check(cdsid)
        assert response is False

        # Raises ValueError if something went wrong with the search
        with pytest.raises(ValueError):
            m.post(ncbi.CDSEARCH_URL, text=text.format("1"))
            response = ncbi.check(cdsid)


def test_CDSearch_retrieve_results():
    cdsid = "QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0"
    text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
        "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
        "#datatype\thitsFull Results\n#status\t{}\tmsg\tJob is still running\n'"
    )
    with requests_mock.Mocker() as m:
        # Case when bad status code received, raises ValueError
        m.post(ncbi.CDSEARCH_URL, text=text.format("1"))
        with pytest.raises(ValueError):
            response = ncbi.retrieve(cdsid)

        # Case when status is 3 (still running), loop ends but no response
        m.post(ncbi.CDSEARCH_URL, text=text.format("3"))
        with pytest.raises(ValueError):
            response = ncbi.retrieve(cdsid, delay=0, max_retries=1)

    anid_tsv = TEST_DIR / "anid.tsv"
    with anid_tsv.open() as anid:
        text = anid.read()

    with requests_mock.Mocker() as m:
        m.post(ncbi.CDSEARCH_URL, text=text)

        # Successful run, test save to file handle
        response = ncbi.retrieve(cdsid)

        assert response.text == text


def test_CDSearch_retrieve_results_no_response(monkeypatch):
    def no_response(cdsid):
        return None

    monkeypatch.setattr(ncbi, "check", no_response)

    with pytest.raises(ValueError):
        ncbi.retrieve("test", delay=0, max_retries=1)


def test_CDSearch_retrieve_results_input():
    with pytest.raises(ValueError):
        ncbi.retrieve("test", delay=5)


def test_efetch_sequences_bad_response(monkeypatch):
    from Bio import Entrez

    def mocked_efetch(db, id, rettype, retmode):
        raise IOError

    monkeypatch.setattr(Entrez, "efetch", mocked_efetch)

    with pytest.raises(IOError):
        ncbi.efetch_sequences(["header"])


def test_efetch_sequences(monkeypatch):
    from Bio import Entrez
    from io import StringIO

    def mocked_efetch(db, id, rettype, retmode):
        return StringIO(">sequence with long definition\nACGT\n>sequence2 test\nACGT")

    monkeypatch.setattr(Entrez, "efetch", mocked_efetch)

    assert ncbi.efetch_sequences(["sequence", "sequence2"]) == {
        "sequence": "ACGT",
        "sequence2": "ACGT",
    }
