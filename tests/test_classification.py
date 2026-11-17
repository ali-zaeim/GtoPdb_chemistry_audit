"""Tests for the notebook's *Structure comparison and suggested actions* functions: how a linked record's structure
compares with the GtoPdb ligand's, layer by layer (`classify_inchi_pair`), when that comparison is made and the
identifier issues reported with it (`assess_xref`), the labels and suggested actions of the curation file, and the
check of a GtoPdb SMILES against its InChI.

The InChIs are named in `tests/inchis.py`, which gives their sources.
"""

import pytest
from rdkit import Chem

from gtp_chemistry_problem import (
    MORE_STEREO_IN,
    STEREO_DEFINED_DIFFERENTLY,
    Action,
    assess_xref,
    classify_inchi_pair,
    difference_label,
    inorganic_structure_action,
    morgan_similarity,
    smiles_inchi_problem,
    suggested_action,
)
from gtp_chemistry_problem import IdentifierIssue as Issue
from tests.inchis import (
    ACETONE,
    ALANINE,
    ALANINE_HCL,
    AZ2158_CHEMBL,
    AZ2158_GTOPDB,
    BUT_1_ENE,
    BUT_2_ENE,
    BUTANE,
    BUTANEDIOL_ONE_CENTRE,
    BUTANEDIOL_TWO_CENTRES,
    CAFFEINE,
    CALCIUM_ACETATE,
    CATION_16C,
    CATION_16F,
    CISPLATIN,
    D_ALANINATE,
    D_ALANINE,
    D_ALANINE_HBR,
    D_ALANINE_HCL,
    DACH,
    DIMETHYL_ETHER,
    ELAMIPRETIDE,
    ETHANOL,
    ETHANOL_1_D,
    GLYCINATE,
    GLYCINE,
    GLYCINE_HBR,
    GLYCINE_HCL,
    HEXENEDIOL_CENTRES_ONLY,
    HEXENEDIOL_CENTRES_ONLY_INVERTED,
    HEXENEDIOL_E_AND_CENTRES,
    HEXENEDIOL_E_ONLY,
    IBUPROFEN,
    IBUPROFEN_METHYL_ESTER,
    L_ALANINE,
    L_ALANINE2_SULFATE,
    L_ALANINE_1_13C,
    L_ALANINE_D3,
    L_ALANINE_HCL,
    LEVALBUTEROL_D_TARTRATE,
    LEVALBUTEROL_L_TARTRATE,
    MESO_TARTARIC,
    MESO_TARTARIC_13C,
    MESO_TARTARIC_13C_OTHER_END,
    NEUTRAL_16G,
    O_CHLOROTOLUENE,
    OXALIPLATIN,
    P_CHLOROTOLUENE,
    PHENYLETHYLAMINE_R_TARTRATE,
    PHENYLETHYLAMINE_TARTRATE,
    PT_DACH_DICHLORIDE,
    R_ETHANOL_1_D,
    R_PENTENEDIOL,
    RICASETRON_GTOPDB,
    RICASETRON_PUBCHEM,
    S_ETHANOL_1_D,
    S_HYDROXYPENTANONE,
    S_PENTENEDIOL,
    SODIUM,
    SODIUM_ACETATE,
    SODIUM_ION,
    TARTARIC_D,
    TARTARIC_L,
    TARTARIC_RELATIVE,
    THEOPHYLLINE,
    THREONINE_C2_ONLY,
    THREONINE_C3_ONLY,
    THREONINE_C3_ONLY_INVERTED,
    THREONINIUM_C3_ONLY,
    TRANS_BUT_2_ENE,
    TRIENE_A,
    TRIENE_B,
    TRIHYDROXYGLUTARIC_2R_4R,
    TRIHYDROXYGLUTARIC_2R_4S,
    UREA,
    UREA_D2,
    YTTRIUM,
)

# missing values in the notebook's tables, which reach these functions through `DataFrame.to_dict("records")`
NAN = float("nan")

# --- classify_inchi_pair ---------------------------------------------------------------------------------------------

RESULT_KEYS = {
    "first_difference",
    "difference",
    "difference_detail",
    "stereo_relation",
    "more_stereo_in",
    "missing_stereo",
    "tanimoto_similarity",
}


def case(test_id, gtp_inchi, xref_inchi, difference, detail="", stereo_relation=""):
    """The expected comparison; the first difference is the layer that `difference` names last."""
    expected = {
        "first_difference": difference.lower().rpartition("differs in ")[2],
        "difference": difference,
        "difference_detail": detail,
        "stereo_relation": stereo_relation,
    }
    return pytest.param(gtp_inchi, xref_inchi, expected, id=test_id)


FORMULA = "Differs in formula"
CONNECTIONS = "Shares formula; differs in connections"
HYDROGENS = "Shares formula and connections; differs in hydrogens"
CHARGE = "Shares formula; differs in charge"
PROTONS = "Shares formula, connections and hydrogens; differs in protons"
STEREO = "Shares formula, connections and hydrogens; differs in stereo"
ISOTOPES = "Shares formula, connections and hydrogens; differs in isotopes"
ISOTOPES_AFTER_STEREO = "Shares formula, connections, hydrogens and stereo; differs in isotopes"

CASES = [
    # Formula: both formulas and their components, whether the atoms and connections are the same, and whether the
    # records have the same parent after standardisation
    case(
        "formula-salt-repeated-components",
        L_ALANINE2_SULFATE,
        L_ALANINE,
        FORMULA,
        "2C3H7NO2.H2O4S vs C3H7NO2; GtoPdb 3 component(s), linked record 1; same parent after standardisation",
    ),
    # the parents are compared without stereo, then with it, so a salt of the enantiomer is not the same parent
    case(
        "formula-salt-and-enantiomer",
        L_ALANINE,
        D_ALANINE_HCL,
        FORMULA,
        "C3H7NO2 vs C3H7NO2.ClH; GtoPdb 1 component(s), linked record 2; same parent after standardisation, apart "
        "from stereo",
    ),
    case(
        "formula-salt-and-xref-has-no-stereo",
        L_ALANINE,
        ALANINE_HCL,
        FORMULA,
        "C3H7NO2 vs C3H7NO2.ClH; GtoPdb 1 component(s), linked record 2; same parent after standardisation, apart "
        "from stereo",
    ),
    case(
        "formula-other-counterion",
        GLYCINE_HCL,
        GLYCINE_HBR,
        FORMULA,
        "C2H5NO2.ClH vs C2H5NO2.BrH; same parent after standardisation",
    ),
    case(
        "formula-other-counterion-and-enantiomer",
        L_ALANINE_HCL,
        D_ALANINE_HBR,
        FORMULA,
        "C3H7NO2.ClH vs C3H7NO2.BrH; same parent after standardisation, apart from stereo",
    ),
    case(
        "formula-salt-counterion-metals-are-not-complexes",
        SODIUM_ACETATE,
        CALCIUM_ACETATE,
        FORMULA,
        "C2H4O2.Na vs 2C2H4O2.Ca; GtoPdb 2 component(s), linked record 3; same parent after standardisation",
    ),
    # hydrogens missing on one side: a double bond left out, or a radical
    case("formula-same-atoms-and-connections", BUTANE, BUT_2_ENE, FORMULA, "C4H10 vs C4H8; same atoms and connections"),
    case("formula-similar-compound", IBUPROFEN, IBUPROFEN_METHYL_ESTER, FORMULA, "C13H18O2 vs C14H20O2"),
    case("formula-close-analogue", CAFFEINE, THEOPHYLLINE, FORMULA, "C8H10N4O2 vs C7H8N4O2"),
    # with the bonds to the metal removed, neither the parents nor the fingerprints describe metal complexes
    case(
        "formula-metal-complexes-sharing-a-ligand",
        OXALIPLATIN,
        PT_DACH_DICHLORIDE,
        FORMULA,
        "C6H14N2.C2H2O4.Pt vs C6H14N2.2ClH.Pt; GtoPdb 3 component(s), linked record 4; metal complexes, compared by "
        "formula only",
    ),
    # a metal complex has no parent, so none is compared
    case(
        "formula-metal-complex-vs-its-free-ligand",
        OXALIPLATIN,
        DACH,
        FORMULA,
        "C6H14N2.C2H2O4.Pt vs C6H14N2; GtoPdb 3 component(s), linked record 1",
    ),
    case("formula-organic-compound-vs-metal", ELAMIPRETIDE, YTTRIUM, FORMULA, "C32H49N9O5 vs Y"),
    # Connections: the part of each structure outside their largest common substructure
    case(
        "connections-positional-isomer",
        O_CHLOROTOLUENE,
        P_CHLOROTOLUENE,
        CONNECTIONS,
        "GtoPdb *c1ccccc1Cl; linked record *c1ccc(Cl)cc1",
    ),
    case(
        "connections-not-isolated",
        ETHANOL,
        DIMETHYL_ETHER,
        CONNECTIONS,
        "the atoms connected differently could not be isolated automatically",
    ),
    # Hydrogens: how the parents after standardisation compare, also as RDKit canonical tautomers
    case("hydrogens-tautomer", S_HYDROXYPENTANONE, S_PENTENEDIOL, HYDROGENS, "same parent as RDKit canonical tautomers"),
    case(
        "hydrogens-tautomer-and-enantiomer",
        S_HYDROXYPENTANONE,
        R_PENTENEDIOL,
        HYDROGENS,
        "same parent as RDKit canonical tautomers, apart from stereo",
    ),
    case(
        "hydrogens-moved-double-bond",
        BUT_1_ENE,
        BUT_2_ENE,
        HYDROGENS,
        "different parents, also as RDKit canonical tautomers",
    ),
    case("hydrogens-same-parent", CATION_16C, CATION_16F, HYDROGENS, "same parent after standardisation"),
    # the protons differ too, but after the hydrogens, where the comparison stops
    case("hydrogens-before-protons", CATION_16C, NEUTRAL_16G, HYDROGENS, "same parent after standardisation"),
    # Charge and protons: the net charge of each record, and how many protons GtoPdb has more or fewer
    case("charge-ion-vs-element", SODIUM_ION, SODIUM, CHARGE, "net charge GtoPdb +1, linked record 0"),
    case(
        "protons-glycine-vs-glycinate",
        GLYCINE,
        GLYCINATE,
        PROTONS,
        "net charge GtoPdb 0, linked record -1; GtoPdb has 1 proton more",
    ),
    # the stereo differs too (an enantiomer), but after the protons
    case(
        "protons-before-stereo",
        L_ALANINE,
        D_ALANINATE,
        PROTONS,
        "net charge GtoPdb 0, linked record -1; GtoPdb has 1 proton more",
    ),
    case(
        "protons-before-stereo-defined-differently",
        THREONINE_C2_ONLY,
        THREONINIUM_C3_ONLY,
        PROTONS,
        "net charge GtoPdb 0, linked record +1; GtoPdb has 1 proton fewer",
    ),
    # Stereo: compared centre by centre only where both records define the same centres and double bonds
    case("stereo-enantiomer-by-m-flag-tartaric", TARTARIC_L, TARTARIC_D, STEREO, "enantiomer", "enantiomer"),
    case(
        "stereo-diastereomer-L-vs-meso-tartaric",
        TARTARIC_L,
        MESO_TARTARIC,
        STEREO,
        "epimer or diastereomer",
        "epimer or diastereomer",
    ),
    case("stereo-EZ-triene", TRIENE_A, TRIENE_B, STEREO, "E/Z isomer", "E/Z isomer"),
    # both leave C2 undefined and define C3 inverted: an enantiomer only if C2 is inverted too
    case(
        "stereo-same-centre-undefined-on-both-the-rest-inverted",
        THREONINE_C3_ONLY,
        THREONINE_C3_ONLY_INVERTED,
        STEREO,
        "enantiomer or diastereomer",
        "enantiomer or diastereomer",
    ),
    case(
        "stereo-component-after-repeated-ones",
        LEVALBUTEROL_L_TARTRATE,
        LEVALBUTEROL_D_TARTRATE,
        STEREO,
        "epimer or diastereomer",
        "epimer or diastereomer",
    ),
    # Records that define different centres are not compared centre by centre: the detail counts what each defines.
    # One record defines more stereo when it defines every centre the other does, written the same way
    case(
        "stereo-undefined-centre-gtp-defines-more",
        BUTANEDIOL_TWO_CENTRES,
        BUTANEDIOL_ONE_CENTRE,
        STEREO,
        "GtoPdb defines more stereo (stereocentres defined: GtoPdb 2 of 2, linked record 1 of 2)",
        "GtoPdb defines more stereo",
    ),
    case(
        "stereo-undefined-centre-xref-defines-more",
        BUTANEDIOL_ONE_CENTRE,
        BUTANEDIOL_TWO_CENTRES,
        STEREO,
        "linked record defines more stereo (stereocentres defined: GtoPdb 1 of 2, linked record 2 of 2)",
        "linked record defines more stereo",
    ),
    # the linked record lists nothing for a component it defines nothing in, so the count is out of both lists
    case(
        "stereo-salt-with-one-component-undefined",
        PHENYLETHYLAMINE_R_TARTRATE,
        PHENYLETHYLAMINE_TARTRATE,
        STEREO,
        "GtoPdb defines more stereo (stereocentres defined: GtoPdb 3 of 3, linked record 2 of 3)",
        "GtoPdb defines more stereo",
    ),
    case(
        "stereo-xref-has-none",
        L_ALANINE,
        ALANINE,
        STEREO,
        "linked record has no stereo",
        "linked record has no stereo",
    ),
    case(
        "stereo-gtp-has-none",
        ALANINE,
        D_ALANINE,
        STEREO,
        "GtoPdb structure has no stereo",
        "GtoPdb structure has no stereo",
    ),
    case(
        "stereo-double-bond-only-xref-has-none",
        TRANS_BUT_2_ENE,
        BUT_2_ENE,
        STEREO,
        "linked record has no stereo",
        "linked record has no stereo",
    ),
    # neither record defines everything the other does
    case(
        "stereo-each-defines-a-different-centre",
        THREONINE_C2_ONLY,
        THREONINE_C3_ONLY,
        STEREO,
        "stereo defined differently on each side (stereocentres defined: GtoPdb 1 of 2, linked record 1 of 2)",
        STEREO_DEFINED_DIFFERENTLY,
    ),
    # a record that defines no centre or no double bond has no /t or /b layer
    case(
        "stereo-two-centres-vs-one-double-bond",
        HEXENEDIOL_CENTRES_ONLY,
        HEXENEDIOL_E_ONLY,
        STEREO,
        "stereo defined differently on each side (stereocentres defined: GtoPdb 2 of 2, linked record 0 of 2; "
        "double bonds defined: GtoPdb 0 of 1, linked record 1 of 1)",
        STEREO_DEFINED_DIFFERENTLY,
    ),
    # the same centres defined, but only one record lists the undefined C3, which next to it may be numbered
    # differently, so the parities are not compared
    case(
        "stereo-same-centres-defined-one-listed-undefined-on-one-side",
        TRIHYDROXYGLUTARIC_2R_4R,
        TRIHYDROXYGLUTARIC_2R_4S,
        STEREO,
        "stereo defined differently on each side (stereocentres defined: GtoPdb 2 of 3, linked record 2 of 3; "
        "1 of the 2 defined on both written inverted)",
        STEREO_DEFINED_DIFFERENTLY,
    ),
    # a centre written inverted next to one defined on one side only may be InChI's numbering, so no conflict is named,
    # whether every shared centre is inverted or only some
    case(
        "stereo-centres-inverted-double-bond-on-one-side",
        HEXENEDIOL_E_AND_CENTRES,
        HEXENEDIOL_CENTRES_ONLY_INVERTED,
        STEREO,
        "stereo defined differently on each side (stereocentres defined: GtoPdb 2 of 2, linked record 2 of 2; "
        "double bonds defined: GtoPdb 1 of 1, linked record 0 of 1; 2 of the 2 defined on both written inverted)",
        STEREO_DEFINED_DIFFERENTLY,
    ),
    # ricasetron and its PubChem record: both bridgeheads look inverted, but defining the record's endo/exo centre 13
    # gives exactly the GtoPdb InChI, since InChI numbers the symmetric bridgeheads by that centre
    case(
        "stereo-ricasetron-bridgeheads-numbered-by-an-undefined-centre",
        RICASETRON_GTOPDB,
        RICASETRON_PUBCHEM,
        STEREO,
        "stereo defined differently on each side (stereocentres defined: GtoPdb 3 of 3, linked record 2 of 3; "
        "2 of the 2 defined on both written inverted)",
        STEREO_DEFINED_DIFFERENTLY,
    ),
    # AZ2158 and its ChEMBL record: GtoPdb defines one end of a cyclohexane cis/trans pair (20) and not the other (19),
    # which on its own describes no configuration
    case(
        "stereo-AZ2158-half-of-a-ring-cis-trans-pair",
        AZ2158_GTOPDB,
        AZ2158_CHEMBL,
        STEREO,
        "stereo defined differently on each side (stereocentres defined: GtoPdb 2 of 3, linked record 3 of 3; "
        "1 of the 2 defined on both written inverted)",
        STEREO_DEFINED_DIFFERENTLY,
    ),
    # with the same labels, the stereo that the labels create is compared
    case("stereo-created-by-isotopes-enantiomers", S_ETHANOL_1_D, R_ETHANOL_1_D, STEREO, "enantiomer", "enantiomer"),
    case(
        "stereo-created-by-isotopes-xref-has-none",
        S_ETHANOL_1_D,
        ETHANOL_1_D,
        STEREO,
        "linked record has no stereo",
        "linked record has no stereo",
    ),
    case(
        "stereo-label-position-sets-enantiomer",
        MESO_TARTARIC_13C,
        MESO_TARTARIC_13C_OTHER_END,
        STEREO,
        "enantiomer",
        "enantiomer",
    ),
    # with different labels, the main stereo layers (the unlabelled structure) are compared, before the isotopes
    case("stereo-before-isotopes", L_ALANINE_D3, D_ALANINE, STEREO, "enantiomer", "enantiomer"),
    # Isotopes: which record carries labels
    case("isotopes-exchangeable-D-only", UREA, UREA_D2, ISOTOPES, "linked record labelled, GtoPdb unlabelled"),
    case("isotopes-different-labels", L_ALANINE_D3, L_ALANINE_1_13C, ISOTOPES_AFTER_STEREO, "different labelling"),
    case(
        "isotopes-label-makes-meso-chiral-vs-unlabelled",
        MESO_TARTARIC_13C,
        MESO_TARTARIC,
        ISOTOPES_AFTER_STEREO,
        "GtoPdb labelled, linked record unlabelled",
    ),
]


def test_identical():
    assert classify_inchi_pair(GLYCINE, GLYCINE) == {
        "first_difference": "identical",
        "difference": "Identical",
        "difference_detail": "",
        "stereo_relation": "",
        "more_stereo_in": None,
        "missing_stereo": None,
        "tanimoto_similarity": 1.0,
    }


@pytest.mark.parametrize(("gtp_inchi", "xref_inchi", "expected"), CASES)
def test_classification(gtp_inchi, xref_inchi, expected):
    result = classify_inchi_pair(gtp_inchi, xref_inchi)
    assert result.keys() == RESULT_KEYS
    assert {key: result[key] for key in expected} == expected


@pytest.mark.parametrize(
    ("gtp_inchi", "xref_inchi", "more_stereo_in", "missing"),
    [
        pytest.param(L_ALANINE, ALANINE, "GtoPdb", "R/S centre", id="gtp-defines-more"),
        pytest.param(ALANINE, L_ALANINE, "linked record", "R/S centre", id="linked-record-defines-more"),
        pytest.param(HEXENEDIOL_E_AND_CENTRES, HEXENEDIOL_CENTRES_ONLY, "GtoPdb", "E/Z", id="double-bond"),
        # the isotopes differ too, but after the stereo
        pytest.param(L_ALANINE_D3, ALANINE, "GtoPdb", "R/S centre", id="before-an-isotope-difference"),
        # the stereo is not compared when an earlier layer differs
        pytest.param(D_ALANINATE, ALANINE, None, None, id="after-a-protons-difference"),
        pytest.param(L_ALANINE_HCL, ALANINE, None, None, id="after-a-formula-difference"),
        # stereo defined on both sides
        pytest.param(L_ALANINE, D_ALANINE, None, None, id="enantiomer"),
        pytest.param(THREONINE_C2_ONLY, THREONINE_C3_ONLY, None, None, id="defined-differently"),
    ],
)
def test_stereo_defined_on_one_side_only(gtp_inchi, xref_inchi, more_stereo_in, missing):
    result = classify_inchi_pair(gtp_inchi, xref_inchi)
    assert (result["more_stereo_in"], result["missing_stereo"]) == (more_stereo_in, missing)


def test_similarity_is_the_morgan_similarity():
    # of the two structures (below 1, so not of one structure with itself)
    result = classify_inchi_pair(O_CHLOROTOLUENE, P_CHLOROTOLUENE)
    assert result["tanimoto_similarity"] == morgan_similarity(O_CHLOROTOLUENE, P_CHLOROTOLUENE) < 1


def test_non_standard_inchi_is_refused_even_when_both_are_the_same():
    with pytest.raises(ValueError, match="not a standard InChI"):
        classify_inchi_pair(TARTARIC_RELATIVE, TARTARIC_RELATIVE)


def test_stereo_written_two_ways_is_refused():
    # L-tartaric acid with the mirror-image /t layer and the opposite /m flag: the same two centres, but InChI only
    # writes the /t layer that sorts first (Fig. 25), so the layers differ only if one was not written by InChI
    rewritten = "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1+,2+/m0/s1"
    with pytest.raises(ValueError, match="no defined centre or double bond"):
        classify_inchi_pair(TARTARIC_L, rewritten)


# --- assess_xref -----------------------------------------------------------------------------------------------------


def xref_row(**changes):
    """A row of the notebook's `xref_table`: by default a well-written ChEMBL ID of a parent compound, with a structure
    on both sides. `tests/test_identifiers.py` tests the identifier issues themselves."""
    row = {
        "xref_db": "ChEMBL",
        "xref_id_as_written": "CHEMBL1",
        "ids_in_field": 1,
        "xref_id_normalised": "CHEMBL1",
        "xref_found": True,
        "chembl_id_issued": NAN,
        "chembl_entity_type": NAN,
        "chembl_status": NAN,
        "chembl_last_active": NAN,
        "chebi_id_type": NAN,
        "chebi_primary_id": NAN,
        "chembl_parent_id": "CHEMBL1",
        "gtp_inchi": L_ALANINE,
        "xref_inchi": L_ALANINE,
    }
    return row | changes


@pytest.mark.parametrize(
    ("changes", "difference", "identifier_issue", "identifier_issue_kinds"),
    [
        # an ID that cannot be read is not looked up, so only how it is written is reported
        pytest.param(
            {
                "xref_id_as_written": "CHEMBL 1",
                "xref_id_normalised": NAN,
                "ids_in_field": 2,
                "xref_found": False,
                "xref_inchi": NAN,
            },
            "Identifier not usable",
            "not a valid ChEMBL identifier; 2 IDs in one ChEMBL field",
            "not_valid|several_ids_in_field",
            id="identifier-not-usable",
        ),
        pytest.param(
            {
                "xref_found": False,
                "xref_inchi": NAN,
                "chembl_parent_id": NAN,
                "chembl_id_issued": True,
                "chembl_entity_type": "COMPOUND",
                "chembl_status": "INACTIVE",
                "chembl_last_active": 28.0,
            },
            "Identifier not resolvable",
            "retired ChEMBL compound (last in ChEMBL 28)",
            "chembl_retired",
            id="identifier-not-resolvable",
        ),
        pytest.param({"gtp_inchi": NAN}, "GtoPdb ligand has no structure", "", "", id="gtp-has-no-structure"),
        pytest.param({"xref_inchi": NAN}, "Linked record has no structure", "", "", id="xref-has-no-structure"),
        # a salt or child form is only reported when its structure differs, which is unknown here
        pytest.param(
            {"chembl_parent_id": "CHEMBL2", "xref_inchi": NAN},
            "Linked record has no structure",
            "",
            "",
            id="child-form-not-compared",
        ),
    ],
)
def test_assess_xref_not_compared(changes, difference, identifier_issue, identifier_issue_kinds):
    assert assess_xref(xref_row(**changes)) == {
        "first_difference": "not compared",
        "difference": difference,
        "difference_detail": "",
        "stereo_relation": "",
        "more_stereo_in": None,
        "missing_stereo": None,
        "tanimoto_similarity": None,
        "metal_differs": False,
        "identifier_issue": identifier_issue,
        "identifier_issue_kinds": identifier_issue_kinds,
    }


@pytest.mark.parametrize(
    ("changes", "identifier_issue", "identifier_issue_kinds"),
    [
        # the default row links a ChEMBL record that is its own parent: the comparison reports the different structure,
        # but there is no other parent to suggest, so the identifier itself has no issue
        pytest.param({"xref_inchi": D_ALANINE}, "", "", id="own-parent-differs"),
        # a salt or child form with the GtoPdb structure is the right link, so its parent is only named when they differ
        pytest.param({"chembl_parent_id": "CHEMBL2"}, "", "", id="child-form-identical"),
        pytest.param(
            {"chembl_parent_id": "CHEMBL2", "xref_inchi": D_ALANINE},
            "ChEMBL salt/child form (parent CHEMBL2)",
            "chembl_child_form",
            id="child-form-differs",
        ),
        pytest.param(
            {
                "xref_db": "ChEBI",
                "xref_id_as_written": "CHEBI:1",
                "xref_id_normalised": "CHEBI:1",
                "chembl_parent_id": NAN,
                "chebi_id_type": "SECONDARY_ID",
                "chebi_primary_id": "CHEBI:2",
            },
            "secondary ChEBI ID (primary CHEBI:2)",
            "chebi_secondary_id",
            id="secondary-chebi-id",
        ),
        pytest.param(
            {"xref_id_as_written": " CHEMBL1", "ids_in_field": 2},
            "leading/trailing whitespace; 2 IDs in one ChEMBL field",
            "whitespace|several_ids_in_field",
            id="identifier-written-badly",
        ),
    ],
)
def test_assess_xref_compared(changes, identifier_issue, identifier_issue_kinds):
    row = xref_row(**changes)
    assert assess_xref(row) == {
        **classify_inchi_pair(row["gtp_inchi"], row["xref_inchi"]),
        "metal_differs": False,
        "identifier_issue": identifier_issue,
        "identifier_issue_kinds": identifier_issue_kinds,
    }


def test_assess_xref_compares_whatever_the_identifier_issues():
    # a link to another structure is compared even when its ID is also badly written: reporting only the identifier
    # issue would hide that the structures differ
    assessed = assess_xref(xref_row(xref_id_as_written="chembl1", ids_in_field=2, xref_inchi=ACETONE))
    assert assessed["first_difference"] == "formula"
    assert assessed["identifier_issue"] == "non-standard prefix case 'chembl'; 2 IDs in one ChEMBL field"


@pytest.mark.parametrize(
    ("changes", "metal_differs"),
    [
        # standard InChI removes the bonds to the metal, so a metal-containing structure is flagged whenever the
        # linked record does not have exactly the GtoPdb InChI, whatever the layer
        pytest.param({"gtp_inchi": OXALIPLATIN, "xref_inchi": PT_DACH_DICHLORIDE}, True, id="metal-complexes-differ"),
        pytest.param({"gtp_inchi": OXALIPLATIN, "xref_inchi": DACH}, True, id="metal-complex-vs-its-free-ligand"),
        pytest.param({"gtp_inchi": ELAMIPRETIDE, "xref_inchi": YTTRIUM}, True, id="organic-compound-vs-metal"),
        pytest.param({"gtp_inchi": OXALIPLATIN, "xref_inchi": OXALIPLATIN}, False, id="metal-complexes-identical"),
        pytest.param({"gtp_inchi": OXALIPLATIN, "xref_inchi": NAN}, True, id="linked-record-has-no-structure"),
        # a GtoPdb ligand without a structure (section 6.5) and an unresolvable ID (section 4) are reported elsewhere
        pytest.param({"gtp_inchi": NAN, "xref_inchi": CISPLATIN}, False, id="gtp-has-no-structure"),
        pytest.param(
            {
                "gtp_inchi": OXALIPLATIN,
                "xref_inchi": NAN,
                "xref_found": False,
                "chembl_parent_id": NAN,
                "chembl_id_issued": True,
                "chembl_entity_type": "COMPOUND",
                "chembl_status": "INACTIVE",
                "chembl_last_active": 28.0,
            },
            False,
            id="identifier-not-resolvable",
        ),
        # sodium and calcium are salt counterions, not metals that standard InChI disconnects
        pytest.param({"gtp_inchi": SODIUM_ACETATE, "xref_inchi": CALCIUM_ACETATE}, False, id="salt-counterions"),
        pytest.param({"xref_inchi": D_ALANINE}, False, id="no-metal"),
    ],
)
def test_assess_xref_metal_differs(changes, metal_differs):
    assert assess_xref(xref_row(**changes))["metal_differs"] is metal_differs


# --- difference_label ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("first_difference", "shared", "label"),
    [
        ("identical", ["formula", "connections"], "Identical"),
        ("not compared", [], "Not compared"),
        ("formula", [], "Differs in formula"),
        ("connections", ["formula"], "Shares formula; differs in connections"),
        ("hydrogens", ["formula", "connections"], "Shares formula and connections; differs in hydrogens"),
        (
            "isotopes",
            ["formula", "connections", "hydrogens", "stereo"],
            "Shares formula, connections, hydrogens and stereo; differs in isotopes",
        ),
    ],
)
def test_difference_label(first_difference, shared, label):
    assert difference_label(first_difference, shared) == label


# --- suggested_action ------------------------------------------------------------------------------------------------


def differs(layer, detail="", stereo_relation=""):
    """The comparison fields of a row of the notebook's `curation` table whose structures first differ in `layer`."""
    return {
        "first_difference": layer,
        "difference_detail": detail,
        "stereo_relation": stereo_relation,
        "more_stereo_in": MORE_STEREO_IN.get(stereo_relation, NAN),
        "missing_stereo": NAN,
    }


def curation_row(**changes):
    """A row of the notebook's `curation` table: by default a ChEMBL record without the stereocentre that GtoPdb
    defines, with no identifier issue, no other record found with the GtoPdb structure and no other evidence."""
    row = {
        "xref_db": "ChEMBL",
        "xref_id_normalised": "CHEMBL1",
        "xref_name": "ALANINE",
        "identifier_issue_kinds": "",
        **differs("stereo", "linked record has no stereo", "linked record has no stereo"),
        "missing_stereo": "R/S centre",
        "gtp_inchi": L_ALANINE,
        "xref_inchi": ALANINE,
        "radioactive": False,
        "external_consensus_differs": False,
        "metal_differs": False,
        "only_gtp_deposit_agrees": False,
        "exact_record_already_linked": False,
        "other_gtp_ligands_with_xref_structure": "",
        "xref_id_with_gtp_structure": NAN,
        "xref_name_with_gtp_structure": NAN,
        "chebi_primary_id": NAN,
        "chembl_parent_id": NAN,
        "chembl_parent_name": NAN,
        "pubchem_standardisation_removes": NAN,
    }
    return row | changes


# another ChEMBL record, with exactly the GtoPdb structure
WITH_GTP_STRUCTURE = {"xref_id_with_gtp_structure": "CHEMBL2", "xref_name_with_gtp_structure": "L-ALANINE"}
IDENTICAL = differs("identical")
NOT_COMPARED = differs("not compared")
# a difference in each layer, as the notebook's sections 5.1-5.6 list them
ANALOGUE = differs("formula", "C13H18O2 vs C14H20O2")
POSITIONAL_ISOMER = differs("connections", "GtoPdb *c1ccccc1Cl; linked record *c1ccc(Cl)cc1")
TAUTOMER = differs("hydrogens", "same parent as RDKit canonical tautomers")
MOVED_DOUBLE_BOND = differs("hydrogens", "different parents, also as RDKit canonical tautomers")
CHARGES_ELSEWHERE = differs("charge", "same net charge, different component charge layers")
CONJUGATE_BASE = differs("protons", "net charge GtoPdb 0, linked record -1; GtoPdb has 1 proton more")
ENANTIOMER = differs("stereo", "enantiomer", "enantiomer")
DEFINED_DIFFERENTLY = differs("stereo", "stereo defined differently on each side (...)", STEREO_DEFINED_DIFFERENTLY)
LINKED_RECORD_DEFINES_MORE = {
    **differs("stereo", "linked record defines more stereo (...)", "linked record defines more stereo"),
    "missing_stereo": "R/S centre",
}
# a labelled GtoPdb ligand linked to the unlabelled compound, and the reverse
UNLABELLED_LINK = {
    **differs("isotopes", "GtoPdb labelled, linked record unlabelled"),
    "radioactive": True,
    "gtp_inchi": L_ALANINE_1_13C,
    "xref_inchi": L_ALANINE,
}
LABELLED_LINK = {**differs("isotopes", "linked record labelled, GtoPdb unlabelled"), "xref_inchi": L_ALANINE_1_13C}

NO_EXACT_RECORD = "ChEMBL has no record with exactly the GtoPdb structure"
REPLACE_OR_CHANGE = (
    "Replace with (L-ALANINE, CHEMBL2), which has exactly the GtoPdb structure, or change the GtoPdb structure if the "
    "current link is the intended compound"
)
CONSENSUS = "Check the GtoPdb structure: all linked databases agree with each other but not with GtoPdb"
DEPOSIT_NOTE = (
    "Note: the only record that agrees with the GtoPdb structure is the PubChem compound that GtoPdb's own deposit is "
    "merged into"
)
NOT_SEARCHED_NOTE = (
    "Note: GtoPdb draws this radioactive ligand without its isotope labels, so the linked database was not searched "
    "for a record with the labelled structure"
)


def issue_kinds(*kinds):
    """`identifier_issue_kinds` as `assess_xref` writes it."""
    return "|".join(kinds)


def expect(text, *actions, xref_id=None):
    """`suggested_action`'s result: the `Action` of each step, the ID suggested for linking, and the sentences."""
    return {"action": "|".join(actions), "suggested_xref_id": xref_id, "suggested_action": text}


def action(test_id, changes, expected):
    return pytest.param(changes, expected, id=test_id)


ACTIONS = [
    # --- Identifier steps, which depend only on the identifier issues
    action(
        "identifier-text-to-fix",
        {**IDENTICAL, "identifier_issue_kinds": Issue.INVISIBLE_CHARACTER},
        expect("Fix the identifier text: (ALANINE, CHEMBL1)", Action.FIX_ID, xref_id="CHEMBL1"),
    ),
    action(
        "identifier-text-not-usable",
        {
            **NOT_COMPARED,
            **WITH_GTP_STRUCTURE,
            "xref_id_normalised": NAN,
            "identifier_issue_kinds": issue_kinds(Issue.PREFIX_CASE, Issue.NOT_VALID),
        },
        expect(
            "Fix the identifier text; (L-ALANINE, CHEMBL2) has the GtoPdb structure", Action.FIX_ID, xref_id="CHEMBL2"
        ),
    ),
    action(
        "chembl-id-never-issued",
        {**NOT_COMPARED, "identifier_issue_kinds": Issue.CHEMBL_NEVER_ISSUED},
        expect("Fix the identifier: it was never issued (likely a typo)", Action.FIX_ID),
    ),
    action(
        "chembl-id-not-a-compound",
        {**NOT_COMPARED, **WITH_GTP_STRUCTURE, "identifier_issue_kinds": Issue.CHEMBL_NOT_A_COMPOUND},
        expect(
            "Replace: this ChEMBL ID belongs to another kind of ChEMBL record, not a compound; (L-ALANINE, CHEMBL2) "
            "has the GtoPdb structure",
            Action.REPLACE_LINK,
            xref_id="CHEMBL2",
        ),
    ),
    action(
        "retired-chembl-id-record-with-gtp-structure",
        {**NOT_COMPARED, **WITH_GTP_STRUCTURE, "identifier_issue_kinds": Issue.CHEMBL_RETIRED},
        expect(
            "Replace the retired ChEMBL compound ID with (L-ALANINE, CHEMBL2), which has the GtoPdb structure",
            Action.REPLACE_LINK,
            xref_id="CHEMBL2",
        ),
    ),
    action(
        "retired-chembl-id-no-record-with-gtp-structure",
        {**NOT_COMPARED, "identifier_issue_kinds": Issue.CHEMBL_RETIRED},
        expect(
            "Replace or remove the retired ChEMBL compound ID: no current ChEMBL compound has exactly the GtoPdb "
            "structure",
            Action.REMOVE_LINK,
        ),
    ),
    action(
        "retired-chembl-id-gtp-has-no-structure",
        {**NOT_COMPARED, "identifier_issue_kinds": Issue.CHEMBL_RETIRED, "gtp_inchi": NAN},
        expect("Replace the retired ChEMBL compound ID with the current one", Action.REPLACE_LINK),
    ),
    action(
        "chebi-id-deleted",
        {**NOT_COMPARED, "xref_db": "ChEBI", "identifier_issue_kinds": Issue.CHEBI_NOT_FOUND},
        expect("Remove or replace the obsolete identifier", Action.REMOVE_LINK),
    ),
    action(
        "pubchem-cid-without-structure",
        {**NOT_COMPARED, "xref_db": "PubChem", "identifier_issue_kinds": Issue.PUBCHEM_NO_STRUCTURE},
        expect("Remove or replace the obsolete identifier", Action.REMOVE_LINK),
    ),
    # every identifier issue gets its step
    action(
        "secondary-chebi-id-and-two-ids-in-one-field",
        {
            **IDENTICAL,
            "xref_db": "ChEBI",
            "identifier_issue_kinds": issue_kinds(Issue.CHEBI_SECONDARY_ID, Issue.SEVERAL_IDS_IN_FIELD),
            "chebi_primary_id": "CHEBI:2",
        },
        expect(
            "Replace with the primary ChEBI ID CHEBI:2. Check whether both IDs are intended",
            Action.REPLACE_LINK,
            Action.FIX_ID,
            xref_id="CHEBI:2",
        ),
    ),
    action(
        "two-ids-in-one-field",
        {**IDENTICAL, "identifier_issue_kinds": Issue.SEVERAL_IDS_IN_FIELD},
        expect("Check whether both IDs are intended", Action.FIX_ID),
    ),
    action(
        "three-ids-in-one-field",
        {**IDENTICAL, "identifier_issue_kinds": Issue.SEVERAL_IDS_IN_FIELD, "ids_in_field": 3},
        expect("Check whether all 3 IDs are intended", Action.FIX_ID),
    ),
    # --- Identifier and structure steps together: neither hides the other
    action(
        "identifier-text-and-hydrogens",
        {**TAUTOMER, "identifier_issue_kinds": Issue.PREFIX_CASE},
        expect(
            "Fix the identifier text: (ALANINE, CHEMBL1). Keep the link (differs in hydrogens: same parent as RDKit "
            f"canonical tautomers); {NO_EXACT_RECORD}",
            Action.FIX_ID,
            Action.KEEP_LINK,
            xref_id="CHEMBL1",
        ),
    ),
    action(
        "secondary-chebi-id-and-linked-record-defines-more-stereo",
        {
            **LINKED_RECORD_DEFINES_MORE,
            "xref_db": "ChEBI",
            "identifier_issue_kinds": Issue.CHEBI_SECONDARY_ID,
            "chebi_primary_id": "CHEBI:2",
        },
        expect(
            "Replace with the primary ChEBI ID CHEBI:2. Check the GtoPdb structure: the linked record defines stereo "
            "that GtoPdb leaves undefined (R/S centre)",
            Action.REPLACE_LINK,
            Action.CHECK_GTP_STRUCTURE,
            xref_id="CHEBI:2",
        ),
    ),
    action(
        "chembl-child-form-and-formula",
        {
            **ANALOGUE,
            "identifier_issue_kinds": Issue.CHEMBL_CHILD_FORM,
            "chembl_parent_id": "CHEMBL3",
            "chembl_parent_name": "ALANINE",
        },
        expect(
            "Consider linking the ChEMBL parent compound (ALANINE, CHEMBL3). Check which structure is intended "
            f"(differs in formula: C13H18O2 vs C14H20O2); {NO_EXACT_RECORD}",
            Action.REPLACE_LINK,
            Action.CHECK_WHICH,
            xref_id="CHEMBL3",
        ),
    ),
    action(
        "chembl-child-form-parent-not-retrieved",
        {**ANALOGUE, "identifier_issue_kinds": Issue.CHEMBL_CHILD_FORM},
        expect(
            "Consider linking the ChEMBL parent compound. Check which structure is intended (differs in formula: "
            f"C13H18O2 vs C14H20O2); {NO_EXACT_RECORD}",
            Action.REPLACE_LINK,
            Action.CHECK_WHICH,
        ),
    ),
    # a record with exactly the GtoPdb structure is the one to link, so other records suggested are left out...
    action(
        "chembl-child-form-whose-parent-has-the-gtp-structure",
        {
            **ANALOGUE,
            **WITH_GTP_STRUCTURE,
            "identifier_issue_kinds": Issue.CHEMBL_CHILD_FORM,
            "chembl_parent_id": "CHEMBL2",
            "chembl_parent_name": "L-ALANINE",
        },
        expect(
            f"{REPLACE_OR_CHANGE} (the current link differs in formula: C13H18O2 vs C14H20O2)",
            Action.REPLACE_OR_CHANGE,
            xref_id="CHEMBL2",
        ),
    ),
    # ...but the identifier text is still fixed, in case the link stays
    action(
        "identifier-text-and-record-with-gtp-structure",
        {**TAUTOMER, **WITH_GTP_STRUCTURE, "identifier_issue_kinds": Issue.PREFIX_CASE},
        expect(
            "Fix the identifier text: (ALANINE, CHEMBL1). "
            f"{REPLACE_OR_CHANGE} (the current link differs in hydrogens: same parent as RDKit canonical tautomers)",
            Action.FIX_ID,
            Action.REPLACE_OR_CHANGE,
            xref_id="CHEMBL2",
        ),
    ),
    # and removing the link leaves nothing to fix
    action(
        "two-ids-in-one-field-and-the-other-has-the-gtp-structure",
        {
            **ANALOGUE,
            **WITH_GTP_STRUCTURE,
            "exact_record_already_linked": True,
            "identifier_issue_kinds": Issue.SEVERAL_IDS_IN_FIELD,
        },
        expect(
            "Remove the cross-reference: GtoPdb also links (L-ALANINE, CHEMBL2), which has the GtoPdb structure",
            Action.REMOVE_LINK,
        ),
    ),
    # --- Any layer: a record with exactly the GtoPdb structure is linked instead, or the GtoPdb structure changed
    action(
        "record-with-gtp-structure-already-linked",
        {**ENANTIOMER, **WITH_GTP_STRUCTURE, "exact_record_already_linked": True},
        expect(
            "Remove the cross-reference: GtoPdb also links (L-ALANINE, CHEMBL2), which has the GtoPdb structure",
            Action.REMOVE_LINK,
        ),
    ),
    *[
        action(
            f"{test_id}-record-with-gtp-structure",
            {**fields, **WITH_GTP_STRUCTURE},
            expect(f"{REPLACE_OR_CHANGE} (the current link {how})", Action.REPLACE_OR_CHANGE, xref_id="CHEMBL2"),
        )
        for test_id, fields, how in [
            ("formula", ANALOGUE, "differs in formula: C13H18O2 vs C14H20O2"),
            ("connections", POSITIONAL_ISOMER, "differs in connections: GtoPdb *c1ccccc1Cl; linked record *c1ccc(Cl)cc1"),
            ("hydrogens-same-parent", TAUTOMER, "differs in hydrogens: same parent as RDKit canonical tautomers"),
            (
                "hydrogens-different-parents",
                MOVED_DOUBLE_BOND,
                "differs in hydrogens: different parents, also as RDKit canonical tautomers",
            ),
            ("charge", CHARGES_ELSEWHERE, "differs in charge: same net charge, different component charge layers"),
            ("protons", CONJUGATE_BASE, "differs in protons: net charge GtoPdb 0, linked record -1; GtoPdb has 1 proton more"),
            ("opposite-configuration", ENANTIOMER, "differs in stereo: enantiomer"),
            ("stereo-defined-differently", DEFINED_DIFFERENTLY, "differs in stereo: stereo defined differently on each side (...)"),
            ("gtp-defines-more-stereo", {}, "differs in stereo: linked record has no stereo"),
            ("isotopes", UNLABELLED_LINK, "differs in isotopes: GtoPdb labelled, linked record unlabelled"),
        ]
    ],
    # PubChem titles a compound without a name "CID n": that is left out, as is a name missing at the end
    action(
        "three-pubchem-records-with-gtp-structure",
        {
            **POSITIONAL_ISOMER,
            "xref_db": "PubChem",
            "xref_id_with_gtp_structure": "1|2|3",
            "xref_name_with_gtp_structure": "CID 1|Two",
        },
        expect(
            "Replace with PubChem CID 1 or (Two, PubChem CID 2) or PubChem CID 3, which all have exactly the GtoPdb "
            "structure, or change the GtoPdb structure if the current link is the intended compound (the current link "
            "differs in connections: GtoPdb *c1ccccc1Cl; linked record *c1ccc(Cl)cc1)",
            Action.REPLACE_OR_CHANGE,
            xref_id="1|2|3",
        ),
    ),
    # the evidence becomes notes
    action(
        "record-with-gtp-structure-another-ligand-has-the-linked-structure",
        {**ENANTIOMER, **WITH_GTP_STRUCTURE, "other_gtp_ligands_with_xref_structure": "7"},
        expect(
            f"{REPLACE_OR_CHANGE} (the current link differs in stereo: enantiomer). Note: the current link has the "
            "structure of GtoPdb ligand 7",
            Action.REPLACE_OR_CHANGE,
            xref_id="CHEMBL2",
        ),
    ),
    action(
        "record-with-gtp-structure-external-databases-agree-against-gtp",
        {**ENANTIOMER, **WITH_GTP_STRUCTURE, "external_consensus_differs": True},
        expect(
            f"{REPLACE_OR_CHANGE} (the current link differs in stereo: enantiomer). Note: all linked databases agree "
            "with each other but not with GtoPdb",
            Action.REPLACE_OR_CHANGE,
            xref_id="CHEMBL2",
        ),
    ),
    action(
        "record-with-gtp-structure-only-gtp-deposit-agrees",
        {
            **ENANTIOMER,
            "xref_db": "PubChem",
            "only_gtp_deposit_agrees": True,
            "xref_id_with_gtp_structure": "5950",
            "xref_name_with_gtp_structure": "L-alanine",
        },
        expect(
            "Replace with (L-alanine, PubChem CID 5950), which has exactly the GtoPdb structure, or change the GtoPdb "
            "structure if the current link is the intended compound (the current link differs in stereo: enantiomer). "
            f"{DEPOSIT_NOTE}",
            Action.REPLACE_OR_CHANGE,
            xref_id="5950",
        ),
    ),
    # --- Any layer without such a record: the linked record is another GtoPdb ligand's structure...
    action(
        "another-ligand-has-the-linked-structure",
        {**ANALOGUE, "other_gtp_ligands_with_xref_structure": "7, 8", "external_consensus_differs": True},
        expect(
            "Remove the cross-reference: the linked record has the structure of GtoPdb ligands 7, 8 instead",
            Action.REMOVE_LINK,
        ),
    ),
    action(
        "labelled-ligand-linked-to-another-ligands-other-compound",
        {**ANALOGUE, "gtp_inchi": L_ALANINE_1_13C, "other_gtp_ligands_with_xref_structure": "7"},
        expect(
            "Remove the cross-reference: the linked record has the structure of GtoPdb ligand 7 instead. Note: the "
            "linked record is unlabelled",
            Action.REMOVE_LINK,
        ),
    ),
    action(
        "unlabelled-ligand-linked-to-another-ligands-labelled-record",
        {**LABELLED_LINK, "other_gtp_ligands_with_xref_structure": "7"},
        expect(
            "Remove the cross-reference: the linked record has the structure of GtoPdb ligand 7 instead",
            Action.REMOVE_LINK,
        ),
    ),
    # except a labelled ligand's link to the unlabelled compound, which is often a GtoPdb ligand too
    action(
        "labelled-ligand-linked-to-another-ligands-unlabelled-record",
        {**UNLABELLED_LINK, "other_gtp_ligands_with_xref_structure": "7"},
        expect(
            "Keep the link (the linked record is the unlabelled compound, GtoPdb ligand 7); ChEMBL has no record with "
            "the GtoPdb isotope labels",
            Action.KEEP_LINK,
        ),
    ),
    # ...or the linked databases agree with each other but not with GtoPdb
    action(
        "external-databases-agree-against-gtp",
        {**TAUTOMER, "external_consensus_differs": True},
        expect(CONSENSUS, Action.CHECK_GTP_STRUCTURE),
    ),
    # a cross-reference with the GtoPdb structure, or whose structures were not compared, has no structure step
    action("identical-whatever-the-other-databases", {**IDENTICAL, "external_consensus_differs": True}, expect("")),
    action("identical", IDENTICAL, expect("")),
    # --- 5.1 Formula, 5.2 connections and the charge of 5.4: which structure is intended
    action(
        "formula",
        ANALOGUE,
        expect(
            f"Check which structure is intended (differs in formula: C13H18O2 vs C14H20O2); {NO_EXACT_RECORD}",
            Action.CHECK_WHICH,
        ),
    ),
    action(
        "connections",
        POSITIONAL_ISOMER,
        expect(
            "Check which structure is intended (differs in connections: GtoPdb *c1ccccc1Cl; linked record "
            f"*c1ccc(Cl)cc1); {NO_EXACT_RECORD}",
            Action.CHECK_WHICH,
        ),
    ),
    action(
        "charge",
        CHARGES_ELSEWHERE,
        expect(
            "Check which structure is intended (differs in charge: same net charge, different component charge layers); "
            f"{NO_EXACT_RECORD}",
            Action.CHECK_WHICH,
        ),
    ),
    # --- 5.3 Hydrogens: kept with the same parent after standardisation, otherwise which structure is intended
    action(
        "hydrogens-same-parent",
        TAUTOMER,
        expect(
            f"Keep the link (differs in hydrogens: same parent as RDKit canonical tautomers); {NO_EXACT_RECORD}",
            Action.KEEP_LINK,
        ),
    ),
    action(
        "hydrogens-same-parent-apart-from-stereo",
        differs("hydrogens", "same parent as RDKit canonical tautomers, apart from stereo"),
        expect(
            "Keep the link (differs in hydrogens: same parent as RDKit canonical tautomers, apart from stereo); "
            f"{NO_EXACT_RECORD}",
            Action.KEEP_LINK,
        ),
    ),
    action(
        "hydrogens-different-parents",
        MOVED_DOUBLE_BOND,
        expect(
            "Check which structure is intended (differs in hydrogens: different parents, also as RDKit canonical "
            f"tautomers); {NO_EXACT_RECORD}",
            Action.CHECK_WHICH,
        ),
    ),
    action(
        "hydrogens-without-detail",
        differs("hydrogens"),
        expect(f"Check which structure is intended (differs in hydrogens); {NO_EXACT_RECORD}", Action.CHECK_WHICH),
    ),
    # --- 5.4 Protons: kept
    action(
        "protons",
        CONJUGATE_BASE,
        expect(
            "Keep the link (differs in protons: net charge GtoPdb 0, linked record -1; GtoPdb has 1 proton more); "
            f"{NO_EXACT_RECORD}",
            Action.KEEP_LINK,
        ),
    ),
    # --- 5.5 Stereo
    action(
        "opposite-configuration",
        ENANTIOMER,
        expect(
            "Check stereochemistry: ChEMBL has no record with the GtoPdb stereochemistry, so review the GtoPdb "
            "structure",
            Action.CHECK_GTP_STRUCTURE,
        ),
    ),
    action(
        "stereo-defined-differently",
        DEFINED_DIFFERENTLY,
        expect(
            "Review: the records define different stereocentres, so the InChIs cannot tell whether they are the same "
            "stereoisomer",
            Action.REVIEW,
        ),
    ),
    action(
        "gtp-defines-more-stereo",
        {},
        expect(
            "Keep the link (the linked record leaves stereo undefined: R/S centre); ChEMBL has no record with the "
            "GtoPdb stereochemistry",
            Action.KEEP_LINK,
        ),
    ),
    action(
        "gtp-defines-more-stereo-pubchem-standardisation-removes-it",
        {"xref_db": "PubChem", "missing_stereo": "ring cis/trans", "pubchem_standardisation_removes": "ring cis/trans"},
        expect(
            "Keep the link (the linked record leaves stereo undefined: ring cis/trans); PubChem's standardisation "
            "removes this stereo (ring cis/trans) from the GtoPdb structure and reproduces the linked InChI",
            Action.KEEP_LINK,
        ),
    ),
    action(
        "linked-record-defines-more-stereo",
        LINKED_RECORD_DEFINES_MORE,
        expect(
            "Check the GtoPdb structure: the linked record defines stereo that GtoPdb leaves undefined (R/S centre)",
            Action.CHECK_GTP_STRUCTURE,
        ),
    ),
    # --- 5.6 Isotopes
    action(
        "labelled-ligand-linked-to-the-unlabelled-compound",
        UNLABELLED_LINK,
        expect(
            "Keep the link (the linked record is the unlabelled compound); ChEMBL has no record with the GtoPdb "
            "isotope labels",
            Action.KEEP_LINK,
        ),
    ),
    action(
        "radioactive-ligand-drawn-without-its-isotope",
        {**LABELLED_LINK, "radioactive": True},
        expect(
            "Check the GtoPdb structure: this radioactive ligand is drawn without the isotope labels that the linked "
            "record carries",
            Action.CHECK_GTP_STRUCTURE,
        ),
    ),
    action(
        "unlabelled-ligand-linked-to-a-labelled-record",
        LABELLED_LINK,
        expect(
            "Check which structure is intended (differs in isotopes: linked record labelled, GtoPdb unlabelled); "
            f"{NO_EXACT_RECORD}",
            Action.CHECK_WHICH,
        ),
    ),
    # --- Notes, which never change the action
    action(
        "only-gtp-deposit-agrees",
        {**ANALOGUE, "only_gtp_deposit_agrees": True},
        expect(
            f"Check which structure is intended (differs in formula: C13H18O2 vs C14H20O2); {NO_EXACT_RECORD}. "
            f"{DEPOSIT_NOTE}",
            Action.CHECK_WHICH,
        ),
    ),
    action(
        "identifier-text-hydrogens-and-only-gtp-deposit-agrees",
        {**TAUTOMER, "identifier_issue_kinds": Issue.PREFIX_CASE, "only_gtp_deposit_agrees": True},
        expect(
            "Fix the identifier text: (ALANINE, CHEMBL1). Keep the link (differs in hydrogens: same parent as RDKit "
            f"canonical tautomers); {NO_EXACT_RECORD}. {DEPOSIT_NOTE}",
            Action.FIX_ID,
            Action.KEEP_LINK,
            xref_id="CHEMBL1",
        ),
    ),
    # a radioactive ligand drawn without its isotope has the unlabelled compound's InChIKey, so it was not searched for
    action(
        "radioactive-ligand-drawn-without-its-isotope-not-searched",
        {**ANALOGUE, "radioactive": True},
        expect(
            f"Check which structure is intended (differs in formula: C13H18O2 vs C14H20O2). {NOT_SEARCHED_NOTE}",
            Action.CHECK_WHICH,
        ),
    ),
    # an opposite configuration does not claim that the linked database has no record of the GtoPdb stereoisomer
    action(
        "radioactive-ligand-drawn-without-its-isotope-opposite-configuration",
        {**ENANTIOMER, "radioactive": True},
        expect(
            "Check stereochemistry: the linked record has the opposite configuration, so review the GtoPdb structure. "
            f"{NOT_SEARCHED_NOTE}",
            Action.CHECK_GTP_STRUCTURE,
        ),
    ),
    # a kept link still says why
    action(
        "radioactive-ligand-drawn-without-its-isotope-kept",
        {**TAUTOMER, "radioactive": True},
        expect(
            f"Keep the link (differs in hydrogens: same parent as RDKit canonical tautomers). {NOT_SEARCHED_NOTE}",
            Action.KEEP_LINK,
        ),
    ),
    action(
        "labelled-ligand-linked-to-an-unlabelled-analogue",
        {**ANALOGUE, "gtp_inchi": L_ALANINE_1_13C},
        expect(
            f"Check which structure is intended (differs in formula: C13H18O2 vs C14H20O2); {NO_EXACT_RECORD}. Note: "
            "the linked record is unlabelled",
            Action.CHECK_WHICH,
        ),
    ),
]


# --- 5.7 Metal-containing structures: standard InChI removes the bonds to the metal, so whatever the layer the curator
# checks which structure is intended, unless an exact record decides it
METAL_STEP = (
    "Check which structure is intended (differs in charge: same net charge, different component charge layers): standard "
    "InChI removes the bonds to metals, so it cannot tell a different compound from the same one drawn differently"
)
ACTIONS += [
    action(
        "metal-containing",
        {**CHARGES_ELSEWHERE, "metal_differs": True},
        expect(f"{METAL_STEP}; {NO_EXACT_RECORD}", Action.CHECK_WHICH),
    ),
    action(
        "metal-containing-record-with-gtp-structure",
        {**CHARGES_ELSEWHERE, **WITH_GTP_STRUCTURE, "metal_differs": True},
        expect(
            f"{REPLACE_OR_CHANGE} (the current link differs in charge: same net charge, different component charge layers)",
            Action.REPLACE_OR_CHANGE,
            xref_id="CHEMBL2",
        ),
    ),
    action(
        "metal-containing-linked-record-without-structure",
        {**NOT_COMPARED, "metal_differs": True, "xref_inchi": NAN},
        expect("Compare the structures by hand: the linked record has no structure to compare with", Action.REVIEW),
    ),
]


@pytest.mark.parametrize(("changes", "expected"), ACTIONS)
def test_suggested_action(changes, expected):
    assert suggested_action(curation_row(**changes)) == expected


# --- inorganic_structure_action --------------------------------------------------------------------------------------

ADD = "Add this SMILES, InChI and InChIKey"


@pytest.mark.parametrize(
    ("linked_cids", "cids_with_structure", "cid_names", "expected"),
    [
        pytest.param([], [], {}, f"{ADD} (PubChem has no compound with this structure)", id="no-pubchem-compound"),
        pytest.param(["5460341"], ["5460341", "271"], {}, ADD, id="linked-compound-has-the-structure"),
        pytest.param(
            ["104727"],
            ["5460341"],
            {"5460341": "Calcium"},
            f"{ADD}, and replace the PubChem link with (Calcium, PubChem CID 5460341), which has this structure",
            id="linked-compound-has-another-structure",
        ),
        pytest.param(
            ["104727"],
            ["5460341", "271"],
            {"5460341": "Calcium", "271": NAN},
            f"{ADD}, and replace the PubChem link with (Calcium, PubChem CID 5460341) or PubChem CID 271, which have "
            "this structure",
            id="two-compounds-have-the-structure",
        ),
        pytest.param(
            [], ["5460341"], {"5460341": "Calcium"}, f"{ADD}, and link (Calcium, PubChem CID 5460341)", id="none-linked"
        ),
    ],
)
def test_inorganic_structure_action(linked_cids, cids_with_structure, cid_names, expected):
    assert inorganic_structure_action(linked_cids, cids_with_structure, cid_names) == expected


# --- smiles_inchi_problem --------------------------------------------------------------------------------------------

SAME_STRUCTURE = "Make the SMILES and InChI describe the same structure"


@pytest.mark.parametrize(
    ("smiles", "inchi", "problem", "action"),
    [
        pytest.param("CC(N)C(=O)O", ALANINE, "", "", id="smiles-reproduces-the-inchi"),
        pytest.param(
            "CC(N)C(=O)O",
            L_ALANINE,
            "the InChI defines stereo that the SMILES lacks (R/S centre)",
            "Add the missing stereo to the SMILES",
            id="inchi-defines-more-stereo",
        ),
        pytest.param(
            "C/C=C/C",
            BUT_2_ENE,
            "the SMILES defines stereo that the InChI lacks (E/Z)",
            "Regenerate the InChI from the SMILES",
            id="smiles-defines-more-stereo",
        ),
        pytest.param(
            "C[C@H](N)C(=O)O",
            L_ALANINE_D3,
            "isotope layers differ: InChI /i1D3, SMILES none",
            "Make the SMILES and InChI describe the same isotopes",
            id="isotopes-differ",
        ),
        # any other difference is described as in the curation file, with the InChI as GtoPdb and the SMILES as the
        # linked record
        pytest.param(
            "C[C@@H](N)C(=O)O", L_ALANINE, f"{STEREO}: enantiomer", SAME_STRUCTURE, id="enantiomer"
        ),
        pytest.param(
            "C[C@H](N)C(=O)[O-]",
            L_ALANINE,
            f"{PROTONS}: net charge InChI 0, SMILES -1; InChI has 1 proton more",
            SAME_STRUCTURE,
            id="protons",
        ),
        pytest.param(
            "Cc1ccc(Cl)cc1",
            O_CHLOROTOLUENE,
            f"{CONNECTIONS}: InChI *c1ccccc1Cl; SMILES *c1ccc(Cl)cc1",
            SAME_STRUCTURE,
            id="connections",
        ),
        pytest.param(
            "C=C(C)O", ACETONE, f"{HYDROGENS}: same parent as RDKit canonical tautomers", SAME_STRUCTURE, id="tautomer"
        ),
    ],
)
def test_smiles_inchi_problem(smiles, inchi, problem, action):
    # the InChI generated from the SMILES comes with the problem
    smiles_inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
    assert smiles_inchi_problem(smiles, inchi) == (smiles_inchi, problem, action)


@pytest.mark.parametrize(
    ("smiles", "problem", "action"),
    [
        pytest.param("CC(N", "SMILES cannot be parsed", "Fix the SMILES", id="unparsable"),
        # RDKit names the atom whose bonds exceed its valence
        pytest.param(
            "C(C)(C)(C)(C)C",
            "invalid SMILES: Explicit valence for atom # 0 C, 5, is greater than permitted",
            "Fix the SMILES: it is not a valid structure",
            id="pentavalent-carbon",
        ),
        # InChI has no notation for an attachment point
        pytest.param("*C", "no InChI can be generated from the SMILES", "Fix the SMILES", id="attachment-point"),
    ],
)
def test_smiles_inchi_problem_without_inchi(smiles, problem, action):
    assert smiles_inchi_problem(smiles, ALANINE) == (None, problem, action)
