"""Tests for the notebook's *Identifiers* functions: normalising ChEMBL, ChEBI and PubChem IDs, what is wrong with a
cross-reference's identifier, and GtoPdb's HTML names."""

import pytest

from gtp_chemistry_problem import IdentifierIssue as Issue
from gtp_chemistry_problem import identifier_issues, normalise_xref, plain_text

NAN = float("nan")


@pytest.mark.parametrize(
    ("name_html", "text"),
    [
        ("&alpha;-bungarotoxin", "α-bungarotoxin"),
        ("[<sup>3</sup>H]NMS", "[3H]NMS"),
    ],
)
def test_plain_text(name_html, text):
    assert plain_text(name_html) == text


# --- normalise_xref --------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("db", "raw_id", "normalised"),
    [
        pytest.param("ChEMBL", " CHEMBL25\t", "CHEMBL25", id="whitespace"),
        pytest.param("PubChem", "\ufeff22\u200d44\ufeff", "2244", id="invisible-characters"),
        pytest.param("ChEBI", "ChEBI:15365", "CHEBI:15365", id="prefix-case"),
        # PubChem reads CIDs as signed 32-bit integers
        pytest.param("PubChem", "2147483647", "2147483647", id="largest-cid"),
        pytest.param("PubChem", "2147483648", None, id="cid-beyond-32-bits"),
        # the whole field must be one ID
        pytest.param("ChEMBL", "CHEMBL25, CHEMBL26", None, id="two-ids"),
        pytest.param("ChEBI", "15365", None, id="chebi-without-prefix"),
        pytest.param("PubChem", "02244", None, id="leading-zero"),
        pytest.param("PubChem", "", None, id="empty"),
    ],
)
def test_normalise_xref(db, raw_id, normalised):
    assert normalise_xref(db, raw_id) == normalised


# --- identifier_issues -----------------------------------------------------------------------------------------------


def xref_row(**changes):
    """The columns of a row of the notebook's `xref_table` that `identifier_issues` reads: by default a well-written
    ChEMBL ID that the molecule API returns, of a parent compound."""
    row = {
        "xref_db": "ChEMBL",
        "xref_id_as_written": "CHEMBL1",
        "ids_in_field": 1,
        "xref_found": True,
        # from ChEMBL's ID lookup, only for IDs the molecule API does not return
        "chembl_id_issued": NAN,
        "chembl_entity_type": NAN,
        "chembl_status": NAN,
        "chembl_last_active": NAN,
        "chembl_parent_id": NAN,
        "chebi_id_type": NAN,
        "chebi_primary_id": NAN,
    }
    return row | changes


def assert_issues(row, first_difference, expected):
    """The issues of `row` are `expected`, in the same order."""
    assert list(identifier_issues(row, first_difference).items()) == list(expected.items())


@pytest.mark.parametrize(
    ("db", "raw_id", "expected"),
    [
        pytest.param("ChEMBL", "CHEMBL25", {}, id="well-written"),
        pytest.param(
            "ChEMBL", " CHEMBL25\t", {Issue.WHITESPACE: "leading/trailing whitespace"}, id="whitespace"
        ),
        # each invisible character named once, in order
        pytest.param(
            "PubChem",
            "\ufeff22\u200d44\ufeff",
            {Issue.INVISIBLE_CHARACTER: "invisible character U+200D, U+FEFF"},
            id="several-invisible-characters",
        ),
        pytest.param(
            "ChEBI", "ChEBI:15365", {Issue.PREFIX_CASE: "non-standard prefix case 'ChEBI'"}, id="prefix-case"
        ),
        pytest.param(
            "ChEMBL",
            "\u200bchembl25 ",
            {
                Issue.INVISIBLE_CHARACTER: "invisible character U+200B",
                Issue.WHITESPACE: "leading/trailing whitespace",
                Issue.PREFIX_CASE: "non-standard prefix case 'chembl'",
            },
            id="every-problem",
        ),
        pytest.param(
            "ChEMBL", "CHEMBL25, CHEMBL26", {Issue.NOT_VALID: "not a valid ChEMBL identifier"}, id="two-ids"
        ),
        pytest.param(
            "ChEBI",
            "chebi 15365",
            {
                Issue.PREFIX_CASE: "non-standard prefix case 'chebi'",
                Issue.NOT_VALID: "not a valid ChEBI identifier",
            },
            id="prefix-case-and-not-valid",
        ),
    ],
)
def test_identifier_issues_how_written(db, raw_id, expected):
    # an ID that cannot be read is not looked up, so it is not reported as missing from the linked database either
    found = normalise_xref(db, raw_id) is not None
    assert_issues(xref_row(xref_db=db, xref_id_as_written=raw_id, xref_found=found), "identical", expected)


def test_identifier_issues_several_ids_in_field():
    # GtoPdb separates the IDs of one field with "|": each is a cross-reference of its own
    assert_issues(
        xref_row(ids_in_field=3), "identical", {Issue.SEVERAL_IDS_IN_FIELD: "3 IDs in one ChEMBL field"}
    )


NOT_RETURNED = {"xref_found": False, "chembl_id_issued": True, "chembl_entity_type": "COMPOUND"}


@pytest.mark.parametrize(
    ("changes", "expected"),
    [
        pytest.param(
            {"xref_db": "ChEBI", "xref_id_as_written": "CHEBI:1", "xref_found": False},
            {Issue.CHEBI_NOT_FOUND: "does not exist in ChEBI (deleted, or never issued)"},
            id="chebi",
        ),
        pytest.param(
            {"xref_db": "PubChem", "xref_id_as_written": "1", "xref_found": False},
            {Issue.PUBCHEM_NO_STRUCTURE: "PubChem returns no structure for this CID"},
            id="pubchem",
        ),
        # ChEMBL IDs that the molecule API does not return are explained by ChEMBL's ID lookup
        pytest.param(
            {"xref_found": False, "chembl_id_issued": False},
            {Issue.CHEMBL_NEVER_ISSUED: "never issued by ChEMBL (likely a typo)"},
            id="chembl-never-issued",
        ),
        pytest.param(
            {**NOT_RETURNED, "chembl_entity_type": "DOCUMENT", "chembl_status": "ACTIVE"},
            {Issue.CHEMBL_NOT_A_COMPOUND: "ChEMBL document ID, not a compound"},
            id="chembl-document",
        ),
        pytest.param(
            {**NOT_RETURNED, "chembl_status": "INACTIVE", "chembl_last_active": 28.0},
            {Issue.CHEMBL_RETIRED: "retired ChEMBL compound (last in ChEMBL 28)"},
            id="chembl-retired",
        ),
        pytest.param(
            {**NOT_RETURNED, "chembl_status": "INACTIVE"},
            {Issue.CHEMBL_RETIRED: "retired ChEMBL compound"},
            id="chembl-retired-release-unknown",
        ),
        pytest.param(
            {**NOT_RETURNED, "chembl_status": "ACTIVE", "chembl_last_active": 36},
            {Issue.CHEMBL_NOT_RETURNED: "active ChEMBL compound that the molecule API does not return"},
            id="chembl-active",
        ),
    ],
)
def test_identifier_issues_not_found(changes, expected):
    assert_issues(xref_row(**changes), "not compared", expected)


def test_identifier_issues_secondary_chebi_id():
    row = xref_row(
        xref_db="ChEBI",
        xref_id_as_written="ChEBI:1",
        ids_in_field=2,
        chebi_id_type="SECONDARY_ID",
        chebi_primary_id="CHEBI:2",
    )
    assert_issues(
        row,
        "same compound",
        {
            Issue.PREFIX_CASE: "non-standard prefix case 'ChEBI'",
            Issue.SEVERAL_IDS_IN_FIELD: "2 IDs in one ChEBI field",
            Issue.CHEBI_SECONDARY_ID: "secondary ChEBI ID (primary CHEBI:2)",
        },
    )


@pytest.mark.parametrize(
    ("parent", "first_difference", "expected"),
    [
        # a salt or child form with the GtoPdb structure is the right link, so its parent is only named when they differ
        pytest.param("CHEMBL2", "identical", {}, id="identical"),
        pytest.param("CHEMBL2", "not compared", {}, id="not-compared"),
        pytest.param(
            "CHEMBL2",
            "formula",
            {Issue.CHEMBL_CHILD_FORM: "ChEMBL salt/child form (parent CHEMBL2)"},
            id="structures-differ",
        ),
        # a parent compound is its own parent
        pytest.param("CHEMBL1", "stereo", {}, id="own-parent"),
    ],
)
def test_identifier_issues_chembl_child_form(parent, first_difference, expected):
    assert_issues(xref_row(chembl_parent_id=parent), first_difference, expected)
