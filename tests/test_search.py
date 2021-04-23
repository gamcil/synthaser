"""
Test suite for search.py
"""


import pytest

import requests_mock

from synthaser import ncbi, search, models


def test_container_from_query_file(tmp_path):
    fasta = tmp_path / "test.fa"

    fasta.write_text(">one\nabcdef")

    with fasta.open() as fp:
        sc = search._container_from_query_file(fp)

    assert sc == models.SynthaseContainer(
        [models.Synthase(header="one", sequence="abcdef")]
    )


def test_container_from_query_ids_valueerror():
    with pytest.raises(ValueError):
        search._container_from_query_ids(1)


def test_container_from_query_ids_iterable(monkeypatch):
    def mocked_efetch(headers):
        return {"one": "ACGT", "two": "ACGT"}


    monkeypatch.setattr(ncbi, "efetch_sequences", mocked_efetch)

    ids = ["one", "two"]
    sc = search._container_from_query_ids(ids)

    assert sc == models.SynthaseContainer(
        [
            models.Synthase(header="one", sequence="acgt"),
            models.Synthase(header="two", sequence="acgt")
        ]
    )


def test_container_from_query_ids_file(tmp_path, monkeypatch):
    ids = tmp_path / "ids.txt"
    ids.write_text("one")

    def response(_):
        return {"one": "abcdef"}

    monkeypatch.setattr(ncbi, "efetch_sequences", response)

    sc = search._container_from_query_ids([ids.name])

    assert sc == models.SynthaseContainer(
        [models.Synthase(header="one", sequence="abcdef")]
    )


def test_prepare_input_valueerror():
    with pytest.raises(ValueError):
        search.prepare_input()
