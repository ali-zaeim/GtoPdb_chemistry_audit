"""Tests for the notebook's *Retrieval* functions: downloads and API retrieval, cached on disk.

No request leaves the machine: the `session` fixture stands in for `http_session()`, answering each request with the
next response queued by the test and recording what was asked. The replies are shaped like those of the real APIs,
cut down to the fields the code reads.
"""

import datetime as dt
import gzip
import hashlib
import json

import pandas as pd
import pytest
import requests

import gtp_chemistry_problem as retrieval
from tests.inchis import ALANINE, L_ALANINE

NOW = "2026-11-12T09:30:00+00:00"
ALANINE_KEY = "QNAYBMKLOCPYGJ-UHFFFAOYSA-N"
L_ALANINE_KEY = "QNAYBMKLOCPYGJ-REOHCLBHSA-N"


class FakeResponse:
    """The parts of `requests.Response` that the retrieval code uses."""

    def __init__(self, json=None, status_code=200, content=b"", headers=None):
        self.status_code = status_code
        self.content = content
        self.headers = headers or {}
        self._json = json

    def json(self):
        return self._json

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} error", response=self)


class FakeSession:
    """Stands in for the `requests.Session` of `http_session()`."""

    def __init__(self):
        self.responses = []  # queued by the test, in the order the requests are expected
        self.requests = []
        self.sleeps = []

    def reply(self, *responses):
        self.responses.extend(responses)

    def get(self, url, **kwargs):
        return self._answer("GET", url, kwargs)

    def post(self, url, **kwargs):
        return self._answer("POST", url, kwargs)

    def _answer(self, method, url, kwargs):
        self.requests.append({"method": method, "url": url, **kwargs})
        return self.responses.pop(0)


@pytest.fixture
def session(monkeypatch):
    """The fake session, with the clock fixed at NOW and the pauses between PubChem requests recorded, not waited."""
    fake = FakeSession()
    monkeypatch.setattr(retrieval, "http_session", lambda: fake)
    monkeypatch.setattr(retrieval, "utc_now", lambda: NOW)
    monkeypatch.setattr(retrieval.time, "sleep", fake.sleeps.append)
    yield fake
    assert fake.responses == [], "responses queued but never requested"


def assert_rows(table, rows):
    """`table` holds `rows`, with None for a missing value, as pandas reads them into a table."""
    pd.testing.assert_frame_equal(table, pd.DataFrame(rows))


# --- Session and downloads -------------------------------------------------------------------------------------------


def test_http_session():
    session = retrieval.http_session()
    assert session.headers["User-Agent"] == "gtopdb-xref-audit (marimo notebook)"
    retry = session.get_adapter("https://www.ebi.ac.uk/chembl/api/data/molecule.json").max_retries
    assert retry.total == 3
    assert retry.backoff_factor == 0.5
    assert set(retry.status_forcelist) == {429, 500, 502, 503, 504}
    # POST requests (PubChem, UniChem) are retried too
    assert retry.allowed_methods is None


def test_utc_now():
    now = dt.datetime.fromisoformat(retrieval.utc_now())
    assert now.utcoffset() == dt.timedelta(0)
    assert now.microsecond == 0
    assert abs(dt.datetime.now(dt.UTC) - now) < dt.timedelta(minutes=1)


def test_download_to_cache(session, tmp_path):
    path = tmp_path / "raw" / "ligands.csv"
    content = b"a,b\n1,2\n"
    session.reply(FakeResponse(content=content, headers={"Last-Modified": "Mon, 12 Oct 2026 10:00:00 GMT"}))
    log = retrieval.download_to_cache("https://example.org/ligands.csv", path)
    assert session.requests == [{"method": "GET", "url": "https://example.org/ligands.csv", "timeout": 300}]
    assert path.read_bytes() == content
    assert log == {
        "downloaded_from": "https://example.org/ligands.csv",
        "downloaded_at_utc": NOW,
        "last_modified": "Mon, 12 Oct 2026 10:00:00 GMT",
        "file_size_bytes": 8,
        "sha256": hashlib.sha256(content).hexdigest(),
    }
    assert json.loads((tmp_path / "raw" / "ligands.csv.log.json").read_text()) == log


def test_download_to_cache_reuses_the_file(session, tmp_path):
    path = tmp_path / "ligands.csv"
    session.reply(FakeResponse(content=b"first"))
    first = retrieval.download_to_cache("https://example.org/ligands.csv", path)
    # cached: nothing requested
    assert retrieval.download_to_cache("https://example.org/ligands.csv", path) == first
    assert len(session.requests) == 1
    session.reply(FakeResponse(content=b"second"))
    refreshed = retrieval.download_to_cache("https://example.org/ligands.csv", path, refresh=True)
    assert len(session.requests) == 2
    assert path.read_bytes() == b"second"
    assert refreshed["sha256"] == hashlib.sha256(b"second").hexdigest()


def test_download_to_cache_without_its_log(session, tmp_path):
    # a file whose download log is missing is downloaded again, as the log is what the function returns
    path = tmp_path / "ligands.csv"
    path.write_bytes(b"old")
    session.reply(FakeResponse(content=b"new"))
    log = retrieval.download_to_cache("https://example.org/ligands.csv", path)
    assert path.read_bytes() == b"new"
    assert log["file_size_bytes"] == 3


def test_seed_cache(session, tmp_path):
    cache = tmp_path / "raw"
    cache.mkdir()
    (cache / "present.parquet").write_bytes(b"local")
    session.reply(
        FakeResponse(content=b"from github"),
        FakeResponse(status_code=404),  # not in the repository: left to the API calls
    )
    copied = retrieval.seed_cache(
        cache, "https://example.org/cache/", ["present.parquet", "ligands.csv", "missing.parquet"]
    )
    assert copied == ["ligands.csv"]
    # only the missing files are requested, and the local one is kept
    assert [request["url"] for request in session.requests] == [
        "https://example.org/cache/ligands.csv",
        "https://example.org/cache/missing.parquet",
    ]
    assert (cache / "present.parquet").read_bytes() == b"local"
    assert (cache / "ligands.csv").read_bytes() == b"from github"
    assert sorted(path.name for path in cache.iterdir()) == ["ligands.csv", "present.parquet"]


def test_seed_cache_with_nothing_missing(session, tmp_path):
    (tmp_path / "ligands.csv").write_bytes(b"local")
    assert retrieval.seed_cache(tmp_path, "https://example.org/cache/", ["ligands.csv"]) == []
    assert session.requests == []


def test_seed_cache_creates_the_directory(session, tmp_path):
    session.reply(FakeResponse(content=b"x"))
    assert retrieval.seed_cache(tmp_path / "a" / "b", "https://example.org/", ["f.json"]) == ["f.json"]
    assert (tmp_path / "a" / "b" / "f.json").read_bytes() == b"x"


def test_download_to_cache_http_error(session, tmp_path):
    path = tmp_path / "ligands.csv"
    session.reply(FakeResponse(status_code=404))
    with pytest.raises(requests.HTTPError):
        retrieval.download_to_cache("https://example.org/ligands.csv", path)
    assert list(tmp_path.iterdir()) == []


def test_download_gtp_file(session, tmp_path):
    csv = (
        '"# GtoPdb Version: 2026.3.1 - published: 2026-10-12"\n'
        '"Ligand ID","Name","PubChem CID","INN"\n'
        '"6069","10,10-difluoro TXA<sub>2</sub>","73755171",""\n'
    )
    session.reply(FakeResponse(content=csv.encode()))
    table, version, log = retrieval.download_gtp_file("ligands.csv", tmp_path)
    assert session.requests[0]["url"] == "https://www.guidetopharmacology.org/DATA/ligands.csv"
    assert version == "GtoPdb Version: 2026.3.1 - published: 2026-10-12"
    # every field as text, and an empty one as ""
    assert table.to_dict("records") == [
        {"Ligand ID": "6069", "Name": "10,10-difluoro TXA<sub>2</sub>", "PubChem CID": "73755171", "INN": ""}
    ]
    assert log["downloaded_from"] == "https://www.guidetopharmacology.org/DATA/ligands.csv"


def test_download_unichem_mapping(session, tmp_path):
    session.reply(FakeResponse(content=gzip.compress(b"From src:'1'\tTo src:'4'\nCHEMBL290916\t11994\nCHEMBL178090\t8754\n")))
    # the file names the lower source ID first, whichever order the sources are given in
    pairs, log = retrieval.download_unichem_mapping("GtoPdb", "ChEMBL", tmp_path)
    url = "https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src4.txt.gz"
    assert session.requests[0]["url"] == url
    assert log["downloaded_from"] == url
    assert (tmp_path / "unichem_src1src4.txt.gz").exists()
    assert pairs.to_dict("list") == {"ChEMBL": ["CHEMBL290916", "CHEMBL178090"], "GtoPdb": ["11994", "8754"]}


def test_unichem_compound(session, tmp_path):
    session.reply(
        FakeResponse(
            json={
                "compounds": [
                    {
                        "uci": 304,
                        "inchi": {"inchi": L_ALANINE, "keyLayer": "N"},
                        "standardInchiKey": L_ALANINE_KEY,
                        "sources": [],
                    }
                ],
                "response": "Success",
            }
        )
    )
    compound = retrieval.unichem_compound("GtoPdb", "720", tmp_path)
    assert compound == {"uci": 304, "inchi": L_ALANINE, "inchikey": L_ALANINE_KEY, "retrieved_at": NOW}
    assert session.requests == [
        {
            "method": "POST",
            "url": "https://www.ebi.ac.uk/unichem/api/v1/compounds",
            "json": {"type": "sourceID", "sourceID": 4, "compound": "720"},
            "timeout": 60,
        }
    ]
    assert json.loads((tmp_path / "unichem_compound_GtoPdb_720.json").read_text()) == compound
    # cached
    assert retrieval.unichem_compound("GtoPdb", "720", tmp_path) == compound
    assert len(session.requests) == 1


def test_fetch_unichem_gtopdb_batch(session):
    session.reply(
        FakeResponse(
            json={
                "compounds": [
                    {
                        "uci": 304,
                        "inchi": {"inchi": ALANINE, "keyLayer": "N"},
                        "standardInchiKey": ALANINE_KEY,
                        "sources": [
                            {"id": 4, "compoundId": "720"},
                            {"id": 4, "compoundId": "9"},
                            {"id": 4, "compoundId": "31"},
                            {"id": 1, "compoundId": "CHEMBL279597"},
                        ],
                    }
                ],
                "response": "Success",
            }
        ),
        # an ID UniChem does not hold: 200, not 404
        FakeResponse(json={"compounds": [], "response": "Not found"}),
    )
    assert_rows(
        retrieval.fetch_unichem_gtopdb_batch(["720", "999999"]),
        [
            {
                "query_id": "720",
                "found": True,
                "uci": 304,
                "inchi": ALANINE,
                "inchikey": ALANINE_KEY,
                # the other GtoPdb IDs under the UCI, not the other sources' IDs
                "other_gtopdb_ids_same_uci": "31|9",
                "retrieved_at": NOW,
            },
            {
                "query_id": "999999",
                "found": False,
                "uci": None,
                "inchi": None,
                "inchikey": None,
                "other_gtopdb_ids_same_uci": "",
                "retrieved_at": NOW,
            },
        ],
    )
    assert [request["json"] for request in session.requests] == [
        {"type": "sourceID", "sourceID": 4, "compound": "720"},
        {"type": "sourceID", "sourceID": 4, "compound": "999999"},
    ]


# --- Cached look-ups -------------------------------------------------------------------------------------------------


def test_ids_to_fetch():
    gtp_xrefs = pd.DataFrame(
        {
            "xref_db": ["ChEMBL", "ChEMBL", "PubChem", "ChEMBL", "ChEMBL"],
            "xref_id_normalised": ["CHEMBL2", None, "2244", "CHEMBL1", "CHEMBL2"],
        }
    )
    # GtoPdb's own usable IDs, once each
    assert retrieval.ids_to_fetch("ChEMBL", gtp_xrefs) == ["CHEMBL1", "CHEMBL2"]
    assert retrieval.ids_to_fetch("PubChem", gtp_xrefs) == ["2244"]
    assert retrieval.ids_to_fetch("ChEBI", gtp_xrefs) == []


class LengthAPI:
    """A batch API for `fetch_cached` that returns the length of each ID, and fails for the IDs in `failing`."""

    def __init__(self, failing=()):
        self.batches = []
        self.failing = set(failing)

    def __call__(self, ids):
        self.batches.append(ids)
        if self.failing & set(ids):
            raise requests.HTTPError("503 error")
        return pd.DataFrame({"query_id": ids, "found": True, "length": [len(i) for i in ids], "retrieved_at": NOW})


def cached_ids(path):
    return sorted(pd.read_parquet(path).query_id)


def test_fetch_cached(tmp_path):
    path = tmp_path / "raw" / "lengths.parquet"
    api = LengthAPI()
    ids = ["eeeee", "a", "bb", "dddd", "ccc"]
    records = retrieval.fetch_cached(path, ids, api, batch_size=2, workers=3, title="Lengths")
    # in sorted batches
    assert sorted(api.batches) == [["a", "bb"], ["ccc", "dddd"], ["eeeee"]]
    assert records.query_id.tolist() == ["a", "bb", "ccc", "dddd", "eeeee"]
    assert records.length.tolist() == [1, 2, 3, 4, 5]
    assert cached_ids(path) == ["a", "bb", "ccc", "dddd", "eeeee"]


def test_fetch_cached_fetches_only_new_ids(tmp_path):
    path = tmp_path / "lengths.parquet"
    retrieval.fetch_cached(path, ["a", "bb", "ccc"], LengthAPI(), batch_size=10)
    api = LengthAPI()
    records = retrieval.fetch_cached(path, ["ffffff", "a"], api, batch_size=10)
    assert api.batches == [["ffffff"]]
    # only the rows asked for, cached ones first
    assert records.query_id.tolist() == ["a", "ffffff"]
    assert cached_ids(path) == ["a", "bb", "ccc", "ffffff"]


def test_fetch_cached_refresh(tmp_path):
    path = tmp_path / "lengths.parquet"
    retrieval.fetch_cached(path, ["a", "bb"], LengthAPI(), batch_size=10)
    api = LengthAPI()
    records = retrieval.fetch_cached(path, ["bb"], api, batch_size=10, refresh=True)
    assert api.batches == [["bb"]]
    assert records.query_id.tolist() == ["bb"]
    # the cache starts again from what was fetched
    assert cached_ids(path) == ["bb"]


def test_fetch_cached_nothing_to_fetch(tmp_path):
    path = tmp_path / "lengths.parquet"
    records = retrieval.fetch_cached(path, [], LengthAPI(), batch_size=10)
    assert list(records.columns) == ["query_id", "found", "retrieved_at"]
    assert records.empty
    assert not path.exists()


def test_fetch_cached_resumes_after_a_failure(tmp_path):
    path = tmp_path / "lengths.parquet"
    # one ID per batch, saved every two batches: the third fails
    with pytest.raises(requests.HTTPError):
        retrieval.fetch_cached(path, ["a", "b", "c", "d"], LengthAPI(failing={"c"}), batch_size=1)
    assert cached_ids(path) == ["a", "b"]
    api = LengthAPI()
    records = retrieval.fetch_cached(path, ["a", "b", "c", "d"], api, batch_size=1)
    assert api.batches == [["c"], ["d"]]
    assert records.query_id.tolist() == ["a", "b", "c", "d"]


def test_retrieved_on():
    records = pd.DataFrame({"retrieved_at": ["2026-10-29T12:44:34+00:00", "2026-11-02T08:00:00+00:00"]})
    assert retrieval.retrieved_on(records) == "2026-11-02"
    assert retrieval.retrieved_on(records.iloc[:0]) == ""


# --- ChEMBL ----------------------------------------------------------------------------------------------------------


def test_fetch_chembl_batch(session):
    session.reply(
        FakeResponse(
            json={
                "molecules": [
                    {
                        "molecule_chembl_id": "CHEMBL279597",
                        "pref_name": "ALANINE",
                        "molecule_structures": {"standard_inchi": L_ALANINE, "standard_inchi_key": L_ALANINE_KEY},
                        "molecule_hierarchy": {"parent_chembl_id": "CHEMBL279597"},
                    },
                    # a biologic: no structure, and no hierarchy
                    {
                        "molecule_chembl_id": "CHEMBL1201580",
                        "pref_name": None,
                        "molecule_structures": None,
                        "molecule_hierarchy": None,
                    },
                ],
                "page_meta": {},
            }
        )
    )
    records = retrieval.fetch_chembl_batch(["CHEMBL279597", "CHEMBL1201580", "CHEMBL0"])
    assert session.requests == [
        {
            "method": "GET",
            "url": "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
            "params": {
                "molecule_chembl_id__in": "CHEMBL279597,CHEMBL1201580,CHEMBL0",
                "only": "molecule_chembl_id,pref_name,molecule_structures,molecule_hierarchy",
                "limit": 1000,
            },
            "timeout": 30,
        }
    ]
    missing = {"xref_name": None, "inchi": None, "inchikey": None, "chembl_parent_id": None, "retrieved_at": NOW}
    assert_rows(
        records,
        [
            {
                "query_id": "CHEMBL279597",
                "found": True,
                "xref_name": "ALANINE",
                "inchi": L_ALANINE,
                "inchikey": L_ALANINE_KEY,
                "chembl_parent_id": "CHEMBL279597",
                "retrieved_at": NOW,
            },
            {"query_id": "CHEMBL1201580", "found": True, **missing},
            # not returned by the API
            {"query_id": "CHEMBL0", "found": False, **missing},
        ],
    )


def test_fetch_chembl_id_lookup_batch(session):
    session.reply(
        FakeResponse(
            json={
                "chembl_id_lookups": [
                    {"chembl_id": "CHEMBL10041", "entity_type": "COMPOUND", "status": "INACTIVE", "last_active": 28},
                    {"chembl_id": "CHEMBL1121", "entity_type": "DOCUMENT", "status": "ACTIVE", "last_active": None},
                ]
            }
        )
    )
    records = retrieval.fetch_chembl_id_lookup_batch(["CHEMBL10041", "CHEMBL1121", "CHEMBL99999999"])
    assert session.requests == [
        {
            "method": "GET",
            "url": "https://www.ebi.ac.uk/chembl/api/data/chembl_id_lookup.json",
            "params": {"chembl_id__in": "CHEMBL10041,CHEMBL1121,CHEMBL99999999", "limit": 1000},
            "timeout": 60,
        }
    ]
    assert_rows(
        records,
        [
            {
                "query_id": "CHEMBL10041",
                "found": True,
                "entity_type": "COMPOUND",
                "status": "INACTIVE",
                "last_active": 28,
                "retrieved_at": NOW,
            },
            {
                "query_id": "CHEMBL1121",
                "found": True,
                "entity_type": "DOCUMENT",
                "status": "ACTIVE",
                "last_active": None,
                "retrieved_at": NOW,
            },
            {
                "query_id": "CHEMBL99999999",
                "found": False,
                "entity_type": None,
                "status": None,
                "last_active": None,
                "retrieved_at": NOW,
            },
        ],
    )


def test_fetch_chembl_by_inchikey(session):
    session.reply(
        FakeResponse(
            json={
                "molecules": [
                    {
                        "molecule_chembl_id": "CHEMBL279597",
                        "molecule_structures": {"standard_inchi": L_ALANINE, "standard_inchi_key": L_ALANINE_KEY},
                    },
                    {
                        "molecule_chembl_id": "CHEMBL1200617",
                        "molecule_structures": {"standard_inchi": L_ALANINE, "standard_inchi_key": L_ALANINE_KEY},
                    },
                    # without a structure, so without an InChIKey to match
                    {"molecule_chembl_id": "CHEMBL1201580", "molecule_structures": None},
                ]
            }
        )
    )
    records = retrieval.fetch_chembl_by_inchikey([L_ALANINE_KEY, ALANINE_KEY])
    assert session.requests == [
        {
            "method": "GET",
            "url": "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
            "params": {
                "molecule_structures__standard_inchi_key__in": f"{L_ALANINE_KEY},{ALANINE_KEY}",
                "only": "molecule_chembl_id,molecule_structures",
                "limit": 1000,
            },
            "timeout": 60,
        }
    ]
    assert_rows(
        records,
        [
            {"query_id": L_ALANINE_KEY, "found": True, "xref_id": "CHEMBL279597", "inchi": L_ALANINE, "retrieved_at": NOW},
            {"query_id": L_ALANINE_KEY, "found": True, "xref_id": "CHEMBL1200617", "inchi": L_ALANINE, "retrieved_at": NOW},
            {"query_id": ALANINE_KEY, "found": False, "xref_id": None, "inchi": None, "retrieved_at": NOW},
        ],
    )


# --- ChEBI -----------------------------------------------------------------------------------------------------------


def test_fetch_chebi_batch(session):
    session.reply(
        FakeResponse(
            json={
                "CHEBI:16977": {
                    "exists": True,
                    "id_type": "PRIMARY_ID",
                    "primary_chebi_id": "CHEBI:16977",
                    "data": {
                        "ascii_name": "L-alanine",
                        "name": "L-alanine",
                        "stars": 3,
                        "default_structure": {"standard_inchi": L_ALANINE, "standard_inchi_key": L_ALANINE_KEY},
                    },
                },
                # a secondary ID; the name only where there is no ASCII name
                "CHEBI:13069": {
                    "exists": True,
                    "id_type": "SECONDARY_ID",
                    "primary_chebi_id": "CHEBI:16977",
                    "data": {"ascii_name": None, "name": "L-alanine", "stars": 3, "default_structure": None},
                },
                # a class, whose structure fields are empty
                "CHEBI:33709": {
                    "exists": True,
                    "id_type": "PRIMARY_ID",
                    "primary_chebi_id": "CHEBI:33709",
                    "data": {
                        "ascii_name": "amino acid",
                        "stars": 3,
                        "default_structure": {"standard_inchi": "", "standard_inchi_key": ""},
                    },
                },
                "CHEBI:999999999": {"exists": False, "id_type": None, "primary_chebi_id": None, "data": None},
            }
        )
    )
    ids = ["CHEBI:16977", "CHEBI:13069", "CHEBI:33709", "CHEBI:999999999", "CHEBI:0"]
    records = retrieval.fetch_chebi_batch(ids)
    assert session.requests == [
        {
            "method": "GET",
            "url": "https://www.ebi.ac.uk/chebi/backend/api/public/compounds/",
            "params": {"chebi_ids": ",".join(ids)},
            "timeout": 60,
        }
    ]
    no_structure = {"inchi": None, "inchikey": None, "retrieved_at": NOW}
    missing = {"xref_name": None, "chebi_id_type": None, "chebi_primary_id": None, "chebi_stars": None, **no_structure}
    assert_rows(
        records,
        [
            {
                "query_id": "CHEBI:16977",
                "found": True,
                "xref_name": "L-alanine",
                "chebi_id_type": "PRIMARY_ID",
                "chebi_primary_id": "CHEBI:16977",
                "chebi_stars": 3,
                "inchi": L_ALANINE,
                "inchikey": L_ALANINE_KEY,
                "retrieved_at": NOW,
            },
            {
                "query_id": "CHEBI:13069",
                "found": True,
                "xref_name": "L-alanine",
                "chebi_id_type": "SECONDARY_ID",
                "chebi_primary_id": "CHEBI:16977",
                "chebi_stars": 3,
                **no_structure,
            },
            {
                "query_id": "CHEBI:33709",
                "found": True,
                "xref_name": "amino acid",
                "chebi_id_type": "PRIMARY_ID",
                "chebi_primary_id": "CHEBI:33709",
                "chebi_stars": 3,
                **no_structure,
            },
            {"query_id": "CHEBI:999999999", "found": False, **missing},
            # left out of the reply
            {"query_id": "CHEBI:0", "found": False, **missing},
        ],
    )


def test_fetch_chebi_by_inchikey(session):
    session.reply(
        FakeResponse(
            json={
                "results": [
                    {"_source": {"chebi_accession": "CHEBI:16977", "inchi": L_ALANINE, "inchikey": L_ALANINE_KEY, "ascii_name": "L-alanine"}},
                    # the search is full-text: an entry that only mentions the key is left out
                    {"_source": {"chebi_accession": "CHEBI:1", "inchi": ALANINE, "inchikey": ALANINE_KEY, "ascii_name": "x"}},
                    {"_source": {"chebi_accession": "CHEBI:57972", "inchikey": L_ALANINE_KEY, "name": "L-alaninate"}},
                ]
            }
        ),
        FakeResponse(status_code=404),
        FakeResponse(json={}),
    )
    records = retrieval.fetch_chebi_by_inchikey([L_ALANINE_KEY, ALANINE_KEY, "XXXXXXXXXXXXXX-UHFFFAOYSA-N"])
    url = "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/"
    assert session.requests == [
        {"method": "GET", "url": url, "params": {"term": key}, "timeout": 60}
        for key in [L_ALANINE_KEY, ALANINE_KEY, "XXXXXXXXXXXXXX-UHFFFAOYSA-N"]
    ]
    assert_rows(
        records,
        [
            {
                "query_id": L_ALANINE_KEY,
                "found": True,
                "xref_id": "CHEBI:16977",
                "inchi": L_ALANINE,
                "xref_name": "L-alanine",
                "retrieved_at": NOW,
            },
            {
                "query_id": L_ALANINE_KEY,
                "found": True,
                "xref_id": "CHEBI:57972",
                "inchi": None,
                "xref_name": "L-alaninate",
                "retrieved_at": NOW,
            },
            {"query_id": ALANINE_KEY, "found": False, "xref_id": None, "inchi": None, "retrieved_at": NOW},
            {"query_id": "XXXXXXXXXXXXXX-UHFFFAOYSA-N", "found": False, "xref_id": None, "inchi": None, "retrieved_at": NOW},
        ],
    )


def test_fetch_chebi_by_inchikey_http_error(session):
    session.reply(FakeResponse(status_code=500))
    with pytest.raises(requests.HTTPError):
        retrieval.fetch_chebi_by_inchikey([L_ALANINE_KEY])


# --- PubChem ---------------------------------------------------------------------------------------------------------


def test_fetch_pubchem_batch(session):
    session.reply(
        FakeResponse(
            json={
                "PropertyTable": {
                    "Properties": [
                        {"CID": 5950, "InChI": L_ALANINE, "InChIKey": L_ALANINE_KEY, "Title": "L-Alanine"},
                        # a structure too large for an InChI
                        {"CID": 16129778, "Title": "Large peptide"},
                    ]
                }
            }
        )
    )
    records = retrieval.fetch_pubchem_batch(["5950", "16129778", "999999999"])
    assert session.requests == [
        {
            "method": "POST",
            "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/InChI,InChIKey,Title/JSON",
            "data": {"cid": "5950,16129778,999999999"},
            "timeout": 60,
        }
    ]
    # PubChem asks for no more than 5 requests per second
    assert session.sleeps == [0.25]
    assert_rows(
        records,
        [
            {
                "query_id": "5950",
                "found": True,
                "xref_name": "L-Alanine",
                "inchi": L_ALANINE,
                "inchikey": L_ALANINE_KEY,
                "retrieved_at": NOW,
            },
            {
                "query_id": "16129778",
                "found": True,
                "xref_name": "Large peptide",
                "inchi": None,
                "inchikey": None,
                "retrieved_at": NOW,
            },
            # left out of the reply
            {"query_id": "999999999", "found": False, "xref_name": None, "inchi": None, "inchikey": None, "retrieved_at": NOW},
        ],
    )


def test_fetch_pubchem_sid_cids(session):
    session.reply(
        FakeResponse(
            json={
                "InformationList": {
                    "Information": [
                        {"SID": 178102690, "CID": [73755171]},
                        {"SID": 135650000, "CID": [5950, 602]},
                        # a deposit not merged into any compound
                        {"SID": 381120000},
                    ]
                }
            }
        )
    )
    records = retrieval.fetch_pubchem_sid_cids(["178102690", "135650000", "381120000", "1"])
    assert session.requests == [
        {
            "method": "POST",
            "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/cids/JSON",
            "data": {"sid": "178102690,135650000,381120000,1"},
            "timeout": 120,
        }
    ]
    assert session.sleeps == [0.25]
    assert_rows(
        records,
        [
            {"query_id": "178102690", "found": True, "cids": "73755171", "retrieved_at": NOW},
            {"query_id": "135650000", "found": True, "cids": "5950|602", "retrieved_at": NOW},
            {"query_id": "381120000", "found": False, "cids": "", "retrieved_at": NOW},
            {"query_id": "1", "found": False, "cids": "", "retrieved_at": NOW},
        ],
    )


def test_fetch_pubchem_by_inchikey(session):
    session.reply(
        FakeResponse(
            json={
                "PropertyTable": {
                    "Properties": [
                        {"CID": 5950, "InChI": L_ALANINE, "Title": "L-Alanine"},
                        {"CID": 7311724, "InChI": L_ALANINE},
                    ]
                }
            }
        ),
        # PUGREST.NotFound: no compound has this InChIKey
        FakeResponse(status_code=404),
    )
    records = retrieval.fetch_pubchem_by_inchikey([L_ALANINE_KEY, ALANINE_KEY])
    assert session.requests == [
        {
            "method": "GET",
            "url": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{key}/property/InChI,Title/JSON",
            "timeout": 60,
        }
        for key in [L_ALANINE_KEY, ALANINE_KEY]
    ]
    assert session.sleeps == [0.25, 0.25]
    assert_rows(
        records,
        [
            {
                "query_id": L_ALANINE_KEY,
                "found": True,
                "xref_id": "5950",
                "inchi": L_ALANINE,
                "xref_name": "L-Alanine",
                "retrieved_at": NOW,
            },
            {
                "query_id": L_ALANINE_KEY,
                "found": True,
                "xref_id": "7311724",
                "inchi": L_ALANINE,
                "xref_name": None,
                "retrieved_at": NOW,
            },
            {"query_id": ALANINE_KEY, "found": False, "xref_id": None, "inchi": None, "retrieved_at": NOW},
        ],
    )


def test_fetch_pubchem_by_inchikey_http_error(session):
    session.reply(FakeResponse(status_code=503))
    with pytest.raises(requests.HTTPError):
        retrieval.fetch_pubchem_by_inchikey([L_ALANINE_KEY])


def pug_compound(*props):
    """A PubChem standardisation reply with the given (label, value) properties."""
    return {"PC_Compounds": [{"props": [{"urn": {"label": label}, "value": {"sval": value}} for label, value in props]}]}


def test_fetch_pubchem_standardised(session):
    session.reply(
        FakeResponse(json=pug_compound(("IUPAC Name", "(2S)-2-aminopropanoic acid"), ("InChI", L_ALANINE))),
        # PUGREST.BadRequest: PubChem cannot standardise this structure
        FakeResponse(status_code=400),
        FakeResponse(json=pug_compound(("IUPAC Name", "unnamed"))),
    )
    records = retrieval.fetch_pubchem_standardised([L_ALANINE, "InChI=1S/CH5/h1H5", ALANINE])
    assert session.requests == [
        {
            "method": "POST",
            "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/standardize/inchi/JSON",
            "data": {"inchi": inchi},
            "timeout": 60,
        }
        for inchi in [L_ALANINE, "InChI=1S/CH5/h1H5", ALANINE]
    ]
    assert session.sleeps == [0.25, 0.25, 0.25]
    assert_rows(
        records,
        [
            {"query_id": L_ALANINE, "found": True, "standardised_inchi": L_ALANINE, "retrieved_at": NOW},
            {"query_id": "InChI=1S/CH5/h1H5", "found": False, "standardised_inchi": None, "retrieved_at": NOW},
            # no InChI among the properties
            {"query_id": ALANINE, "found": False, "standardised_inchi": None, "retrieved_at": NOW},
        ],
    )


def test_fetch_pubchem_standardised_http_error(session):
    session.reply(FakeResponse(status_code=500))
    with pytest.raises(requests.HTTPError):
        retrieval.fetch_pubchem_standardised([L_ALANINE])
