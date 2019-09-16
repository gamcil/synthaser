#!/usr/bin/env python3

"""
Unit tests for cdsearch.py
"""

from pathlib import Path

import pytest
import requests
import requests_mock

from synthaser.ncbi import CDSearch, efetch_sequences

TEST_DIR = Path(__file__).resolve().parent


@pytest.fixture
def cdsearch():
    """Returns a CDSearch instance, loading from default config file."""
    return CDSearch()


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


def test_CDSearch_set_search_parameters_invalid(cdsearch):
    with pytest.raises(ValueError):
        cdsearch.set_search_parameters({"id": 1})


def test_CDSearch_set_search_parameters_reserved(cdsearch):
    with pytest.raises(ValueError):
        cdsearch.set_search_parameters({"dmode": 1})


def test_CDSearch_set_search_parameters(cdsearch, parameters):
    assert cdsearch.parameters != parameters
    cdsearch.set_search_parameters(parameters)
    parameters.update({"dmode": "full", "tdata": "hits"})
    assert cdsearch.parameters == parameters


@pytest.mark.parametrize(
    "params,exception",
    [
        ({}, ValueError),
        ({"query_file": "doesnt_exist"}, FileNotFoundError),
        ({"query_ids": "not a list or tuple"}, TypeError),
        ({"query_ids": [["not a list/tuple of int/str"]]}, TypeError),
    ],
)
def test_CDSearch_new_search_input_exceptions(cdsearch, params, exception):
    with pytest.raises(exception):
        cdsearch.new_search(**params)


def test_CDSearch_new_search(cdsearch):
    query_ids = ["AN6791.2", "AN6000.2", "AN6431.2"]
    query_file = TEST_DIR / "anid.faa"

    with requests_mock.Mocker() as m, pytest.raises(AttributeError):
        # In case of empty response, or response where regex pattern is not matched
        m.post("https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?", text="")
        cdsearch.new_search(query_ids=query_ids)

    with requests_mock.Mocker() as m:
        text = (
            "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
            "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
            "#datatype\thitsFull Results\n#status\t3\tmsg\tJob is still running\n'"
        )
        m.post("https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?", text=text)
        cdsid = "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8"
        assert cdsearch.new_search(query_ids=query_ids) == cdsid
        assert cdsearch.new_search(query_file=query_file) == cdsid


def test_CDSearch_search_query_file_over_limit(cdsearch, tmp_path):
    fasta_file = tmp_path / "fasta.faa"
    fasta_file.write_text(">\n" * 4001)
    with pytest.raises(ValueError):
        cdsearch._search_query_file(fasta_file)


def test_CDSearch_search_query_ids_over_limit(cdsearch):
    query_ids = ["id"] * 4001
    with pytest.raises(ValueError):
        cdsearch._search_query_ids(query_ids)


def test_CDSearch_check_status_empty_file(cdsearch):
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
        m.get(cdsearch.ncbi, text=text)
        cdsearch.check_status(cdsid)


def test_CDSearch_check_status(cdsearch):
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
        m.get(cdsearch.ncbi, text=text.format("0"))
        response = cdsearch.check_status(cdsid)
        assert response is not None

        # Returns None if job is still running
        m.get(cdsearch.ncbi, text=text.format("3"))
        response = cdsearch.check_status(cdsid)
        assert response is None

        # Raises ValueError if something went wrong with the search
        with pytest.raises(ValueError):
            m.get(cdsearch.ncbi, text=text.format("1"))
            response = cdsearch.check_status(cdsid)


def test_CDSearch_retrieve_results(cdsearch, tmp_path):
    cdsid = "QM3-qcdsearch-B4BAD4B59BC5B80-3E7CFCD3F93E21D0"
    text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n#cdsid\t"
        "QM3-qcdsearch-5C11C54FC403F91-E792C9BEE9818D8\n"
        "#datatype\thitsFull Results\n#status\t{}\tmsg\tJob is still running\n'"
    )
    with requests_mock.Mocker() as m:
        # Case when bad status code received, raises ValueError
        m.get(cdsearch.ncbi, text=text.format("1"))
        with pytest.raises(ValueError):
            response = cdsearch.retrieve_results(cdsid)

        # Case when status is 3 (still running), loop ends but no response
        m.get(cdsearch.ncbi, text=text.format("3"))
        with pytest.raises(ValueError):
            response = cdsearch.retrieve_results(cdsid, delay=0, max_retries=1)

    anid_tsv = TEST_DIR / "anid.tsv"
    with anid_tsv.open() as anid:
        text = anid.read()

    with requests_mock.Mocker() as m:
        m.get(cdsearch.ncbi, text=text)
        results = tmp_path / "results"

        # Successful run, test save to file handle
        with results.open("w") as out:
            response = cdsearch.retrieve_results(cdsid, output=out)

        assert response.text == text

        with results.open() as res:
            assert res.read() == text


def test_CDSearch_retrieve_results_no_response(cdsearch, monkeypatch):
    def no_response(cdsid):
        return None

    monkeypatch.setattr(cdsearch, "check_status", no_response)

    with pytest.raises(ValueError):
        cdsearch.retrieve_results("test", delay=0, max_retries=1)


def test_CDSearch_retrieve_results_input(cdsearch):
    with pytest.raises(ValueError):
        cdsearch.retrieve_results("test", delay=5)


def test_efetch_sequences_bad_response():
    with requests_mock.Mocker() as m, pytest.raises(requests.HTTPError):
        m.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
            status_code="400",
        )
        efetch_sequences(["header"])


def test_efetch_sequences():
    headers = ["sequence", "sequence2"]
    with requests_mock.Mocker() as m:
        m.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
            text=">sequence with long definition\nACGT\n>sequence2 test\nACGT",
        )
        assert efetch_sequences(headers) == {"sequence": "ACGT", "sequence2": "ACGT"}
