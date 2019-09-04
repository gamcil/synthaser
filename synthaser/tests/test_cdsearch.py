#!/usr/bin/env python3

"""
Test suite for cdsearch.py

Cameron Gilchrist
"""

from pathlib import Path

import json
import pytest

import requests
import requests_mock

from synthaser.cdsearch import CDSearch

TEST_DIR = Path(__file__).resolve().parent


@pytest.fixture
def cdsearch():
    """Returns a CDSearch instance, loading from default config file."""
    return CDSearch()


@pytest.fixture
def config():
    return {
        "settings": {
            "db": "pfam",
            "smode": "prec",
            "useid1": "false",
            "compbasedadj": "1",
            "filter": "true",
            "evalue": "0.05",
            "maxhit": "500",
        },
        "output_folder": "..",
    }


def make_tmp_config(tmp_path, config_dict):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "config.json"
    p.write_text(json.dumps(config_dict))
    return p


def test_default_CDSearch_init(cdsearch):
    assert cdsearch.config == {
        "db": "cdd",
        "smode": "auto",
        "useid1": "true",
        "compbasedadj": "0",
        "filter": "false",
        "evalue": "0.01",
        "maxhit": "500",
        "dmode": "full",
        "tdata": "hits",
    }
    assert cdsearch.output_folder == Path("")


def test_CDSearch_load_config(tmp_path, config, cdsearch):
    default_config = cdsearch.config.copy()

    cdsearch.load_config()
    assert cdsearch.config == default_config

    with pytest.raises(IOError):
        cdsearch.load_config("fake_path")

    test_config = make_tmp_config(tmp_path, config)
    cdsearch.load_config(test_config)

    config["settings"]["dmode"] = "full"
    config["settings"]["tdata"] = "hits"
    assert cdsearch.config == config["settings"]


def test_CDSearch_configure(cdsearch):
    invalid = {"settings": {"id": 1}}
    with pytest.raises(ValueError):
        cdsearch.configure(invalid)

    reserved = {"settings": {"dmode": 1}}
    with pytest.raises(ValueError):
        cdsearch.configure(reserved)

    output = {"settings": {}, "output_folder": "test"}
    cdsearch.configure(output)
    assert cdsearch.output_folder == Path("test")


def test_CDSearch_new_search(cdsearch):
    with pytest.raises(ValueError):
        cdsearch.new_search()

    with pytest.raises(FileNotFoundError):
        cdsearch.new_search(query_file="test_query")

    with pytest.raises(ValueError):
        cdsearch.new_search(query_ids="not a list or tuple")

    with pytest.raises(ValueError):
        cdsearch.new_search(query_ids=[["not a list/tuple of int/str"]])

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
        m.get(cdsearch.base_url, text=text.format("0"))
        response = cdsearch.check_status(cdsid)
        assert response is not None

        # Returns None if job is still running
        m.get(cdsearch.base_url, text=text.format("3"))
        response = cdsearch.check_status(cdsid)
        assert response is None

        # Raises ValueError if something went wrong with the search
        with pytest.raises(ValueError):
            m.get(cdsearch.base_url, text=text.format("1"))
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
        m.get(cdsearch.base_url, text=text.format("1"))
        with pytest.raises(ValueError):
            response = cdsearch.retrieve_results(cdsid)

        # Case when status is 3 (still running), returns None, loop ends
        m.get(cdsearch.base_url, text=text.format("3"))
        with pytest.raises(ValueError):
            response = cdsearch.retrieve_results(cdsid, check_interval=0, max_retries=1)

    anid_tsv = TEST_DIR / "anid.tsv"
    with anid_tsv.open() as anid:
        text = anid.read()

    with requests_mock.Mocker() as m:
        m.get(cdsearch.base_url, text=text)
        results = tmp_path / "results"

        # Successful run, test save to file handle
        with results.open("w") as out:
            response = cdsearch.retrieve_results(cdsid, output=out)

        assert response.text == text

        with results.open() as res:
            assert res.read() == text
