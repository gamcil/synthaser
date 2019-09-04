#!/usr/bin/env python3

"""
Test suite for cdsearch.py

Cameron Gilchrist
"""

from pathlib import Path

import json
import pytest

import requests

from synthaser.cdsearch import CDSearch


@pytest.fixture
def cdsearch():
    """Returns a CDSearch instance, loading from default config file."""
    return CDSearch()


ALT_CONFIG = {
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


def make_tmp_config(tmp_path, config):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "config.json"
    p.write_text(json.dumps(config))
    return p


def test_default_CDSearch_init(cdsearch):
    assert cdsearch.base_cfg == Path("cdsearch.json")
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


def test_CDSearch_load_config(tmp_path, cdsearch):
    default_config = cdsearch.config.copy()

    cdsearch.load_config()
    assert cdsearch.config == default_config

    with pytest.raises(IOError):
        cdsearch.load_config("fake_path")

    alt_config = ALT_CONFIG.copy()
    test_config = make_tmp_config(tmp_path, config=alt_config)
    cdsearch.load_config(test_config)

    alt_config["settings"]["dmode"] = "full"
    alt_config["settings"]["tdata"] = "hits"
    assert cdsearch.config == alt_config["settings"]


def test_CDSearch_configure(cdsearch):
    """."""
    # test invalid keys
    # test reserved keys
    # test output folder
    return


def test_CDSearch_new_search(cdsearch):
    """."""
    # test query IDs
    # test query file, IOError
    # no params valueerror
    # mock returned status in response.text
    return


def test_CDSearch_check_status(cdsearch):
    # mock responses
    # test 0 -> returns response, 3 -> returns None, ValueError otherwise
    return


def test_CDSearch_retrieve_results(cdsearch):
    # mock loop
    # check file output to handle or custom path
    return
