#!/usr/bin/env python3

"""
Unit tests for cdsearch.py
"""

from pathlib import Path

import pytest
import requests
import requests_mock

from synthaser import ncbi

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


@pytest.mark.parametrize(
    "params,exception",
    [
        ({}, ValueError),
        ({"query_file": "doesnt_exist"}, FileNotFoundError),
        ({"query_ids": "not a list or tuple"}, TypeError),
        ({"query_ids": [["not a list/tuple of int/str"]]}, TypeError),
    ],
)
def test_CDSearch_new_search_input_exceptions(params, exception):
    with pytest.raises(exception):
        ncbi.CDSearch(**params)


def test_CDSearch_new_search_bad_response():
    query_ids = ["AN6791.2", "AN6000.2", "AN6431.2"]

    with requests_mock.Mocker() as m, pytest.raises(AttributeError):
        # In case of empty response, or response where regex pattern is not matched
        m.post(ncbi.CDSEARCH_URL, text="")
        ncbi.CDSearch(query_ids=query_ids)


def test_search_query_ids():
    query_ids = ["AN6791.2", "AN6000.2", "AN6431.2"]
    with requests_mock.Mocker() as m:
        text = (
            "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
            "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
            "#datatype\thitsFull Results\n#status\t3\tmsg\tJob is still running\n'"
        )
        m.post(ncbi.CDSEARCH_URL, text=text)
        cdsid = "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8"
        assert ncbi._search(query_ids=query_ids) == cdsid


def test_search_query_file():
    query_file = TEST_DIR / "anid.faa"
    with requests_mock.Mocker() as m:
        text = (
            "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
            "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
            "#datatype\thitsFull Results\n#status\t3\tmsg\tJob is still running\n'"
        )
        m.post(ncbi.CDSEARCH_URL, text=text)
        cdsid = "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8"
        assert ncbi._search(query_file=query_file) == cdsid


def test_CDSearch_search_query_file_over_limit(tmp_path):
    fasta_file = tmp_path / "fasta.faa"
    fasta_file.write_text(">\n" * 4001)
    with pytest.raises(ValueError):
        with fasta_file.open() as handle:
            ncbi._search_query_file(handle)


def test_CDSearch_search_query_ids_over_limit():
    query_ids = ["id"] * 4001
    with pytest.raises(ValueError):
        ncbi._search_query_ids(query_ids)


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
    with requests_mock.Mocker() as m, pytest.raises(ValueError):
        m.get(ncbi.CDSEARCH_URL, text=text)
        ncbi._check_status(cdsid)


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
        # Returns Response if job has completed
        m.get(ncbi.CDSEARCH_URL, text=text.format("0"))
        response = ncbi._check_status(cdsid)
        assert response is not None

        # Returns None if job is still running
        m.get(ncbi.CDSEARCH_URL, text=text.format("3"))
        response = ncbi._check_status(cdsid)
        assert response is None

        # Raises ValueError if something went wrong with the search
        with pytest.raises(ValueError):
            m.get(ncbi.CDSEARCH_URL, text=text.format("1"))
            response = ncbi._check_status(cdsid)


def test_CDSearch_retrieve_results():
    cdsid = "QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0"
    text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
        "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
        "#datatype\thitsFull Results\n#status\t{}\tmsg\tJob is still running\n'"
    )
    with requests_mock.Mocker() as m:
        # Case when bad status code received, raises ValueError
        m.get(ncbi.CDSEARCH_URL, text=text.format("1"))
        with pytest.raises(ValueError):
            response = ncbi._retrieve_results(cdsid)

        # Case when status is 3 (still running), loop ends but no response
        m.get(ncbi.CDSEARCH_URL, text=text.format("3"))
        with pytest.raises(ValueError):
            response = ncbi._retrieve_results(cdsid, delay=0, max_retries=1)

    anid_tsv = TEST_DIR / "anid.tsv"
    with anid_tsv.open() as anid:
        text = anid.read()

    with requests_mock.Mocker() as m:
        m.get(ncbi.CDSEARCH_URL, text=text)

        # Successful run, test save to file handle
        response = ncbi._retrieve_results(cdsid)

        assert response.text == text


def test_CDSearch_retrieve_results_no_response(monkeypatch):
    def no_response(cdsid):
        return None

    monkeypatch.setattr(ncbi, "_check_status", no_response)

    with pytest.raises(ValueError):
        ncbi._retrieve_results("test", delay=0, max_retries=1)


def test_CDSearch_retrieve_results_input():
    with pytest.raises(ValueError):
        ncbi._retrieve_results("test", delay=5)


def test_efetch_sequences_bad_response():
    with requests_mock.Mocker() as m, pytest.raises(requests.HTTPError):
        m.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
            status_code="400",
        )
        ncbi.efetch_sequences(["header"])


def test_efetch_sequences():
    headers = ["sequence", "sequence2"]
    with requests_mock.Mocker() as m:
        m.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
            text=">sequence with long definition\nACGT\n>sequence2 test\nACGT",
        )
        assert ncbi.efetch_sequences(headers) == {
            "sequence": "ACGT",
            "sequence2": "ACGT",
        }
