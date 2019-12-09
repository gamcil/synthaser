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


def test_container_from_query_ids_iterable():
    ids = ["one"]

    with requests_mock.Mocker() as m:
        m.post(ncbi.EFETCH_URL, text=">one\nabcdef")
        sc = search._container_from_query_ids(ids)

    assert sc == models.SynthaseContainer(
        [models.Synthase(header="one", sequence="abcdef")]
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


def test_prepare_input_too_many_sequences(monkeypatch):
    def response(ids):
        return models.SynthaseContainer([models.Synthase()] * 4001)

    monkeypatch.setattr(search, "_container_from_query_ids", response)

    with pytest.raises(ValueError):
        search.prepare_input(query_ids=["test"])


def test_search_bad_response(monkeypatch):
    query_ids = ["test"]

    def response(ids):
        return models.SynthaseContainer([models.Synthase(header="test")])

    monkeypatch.setattr(search, "_container_from_query_ids", response)

    with requests_mock.Mocker() as m, pytest.raises(AttributeError):
        # RuntimeError In case of empty response, or response where regex pattern is not matched
        m.post(ncbi.CDSEARCH_URL, text="")
        search.search(mode="remote", query_ids=query_ids)
