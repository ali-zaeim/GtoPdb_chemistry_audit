"""Tests for the notebook's *InChI layers* functions: reading the layers and stereochemistry of standard InChIs.

`stereo_difference`, which compares the stereo of two records, is tested through `classify_inchi_pair` in
`test_classification.py`. The InChIs are named in `tests/inchis.py`, which gives their sources.
"""

import pytest

from gtp_chemistry_problem import (
    STEREO_LAYERS,
    charge_difference,
    comparable_stereo,
    compare_layers,
    defined_centres,
    defined_double_bonds,
    first_difference_keys,
    fixed_hydrogens,
    formula_components,
    inchi_parts,
    isotope_difference,
    isotope_layer,
    labelled_stereo,
    layer_components,
    net_charge,
    parse_standard_inchi,
    share_skeleton,
    stereo_definitions,
    stereo_difference_keys,
    stereo_layers_of,
    stereo_listing,
    strip_layers,
)
from tests.inchis import (
    ACETONE,
    ALANINE,
    BENZENE,
    BENZENE_LABELLED,
    BUT_2_ENE,
    BUTANE,
    BUTANEDIOL_ONE_CENTRE,
    CATION_16C,
    CISPLATIN,
    D_ALANINE2_SULFATE,
    DIMETHYL_ETHER,
    ETHANOL,
    GLYCINATE,
    GLYCINE,
    GLYCINIUM,
    HEXAAMMINE_RUTHENIUM_CHLORIDE,
    HEXADIENE_ONE_DOUBLE_BOND,
    HEXENEDIOL_CENTRES_ONLY_INVERTED,
    HEXENEDIOL_E_AND_CENTRES,
    L_ALANINE,
    L_ALANINE2_SULFATE,
    L_ALANINE_1_13C,
    L_ALANINE_D3,
    MESO_TARTARIC,
    MESO_TARTARIC_13C,
    MESO_TARTARIC_13C_OTHER_END,
    PHENYLETHYLAMINE_R_TARTRATE,
    PHENYLETHYLAMINE_TARTRATE,
    PROPEN_2_OL,
    RU360_GTOPDB,
    RU360_PUBCHEM,
    SODIUM,
    SODIUM_ION,
    TARTARIC_D,
    TARTARIC_L,
    TARTARIC_RELATIVE,
    TEA_D_LACTATE,
    TRIENE_A,
    TRIENE_B,
    UREA,
    UREA_D2,
)

# --- Reading and rebuilding the layers -------------------------------------------------------------------------------


def test_parse_standard_inchi():
    assert parse_standard_inchi(TARTARIC_L) == {
        "formula": "C4H6O6",
        "c": "5-1(3(7)8)2(6)4(9)10",
        "h": "1-2,5-6H,(H,7,8)(H,9,10)",
        "t": "1-,2-",
        "m": "1",
        "s": "1",
    }


def test_parse_standard_inchi_keeps_the_isotopic_sublayers_apart():
    # the label makes meso-tartaric acid chiral: /m and /s exist only in the isotopic layer, next to the main /t
    assert parse_standard_inchi(MESO_TARTARIC_13C) == {
        "formula": "C4H6O6",
        "c": "5-1(3(7)8)2(6)4(9)10",
        "h": "1-2,5-6H,(H,7,8)(H,9,10)",
        "t": "1-,2+",
        "i": "3+1",
        "im": "1",
        "is": "1",
    }
    # exchangeable D only: an empty /i layer with an isotopic /h sublayer
    assert parse_standard_inchi(UREA_D2) == {"formula": "CH4N2O", "c": "2-1(3)4", "h": "(H4,2,3,4)", "i": "", "ih": "D2"}


def test_parse_standard_inchi_refuses_a_non_standard_inchi():
    with pytest.raises(ValueError, match="not a standard InChI"):
        parse_standard_inchi(TARTARIC_RELATIVE)


@pytest.mark.parametrize(
    ("inchi", "drop", "expected"),
    [
        pytest.param(TARTARIC_L, STEREO_LAYERS, "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)", id="stereo"),
        # stereo dropped from the isotopic sublayers too
        pytest.param(
            MESO_TARTARIC_13C,
            STEREO_LAYERS,
            "InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/i3+1",
            id="stereo-with-isotopic-sublayers",
        ),
        # the isotopic layer dropped with its sublayers, leaving the unlabelled meso compound
        pytest.param(MESO_TARTARIC_13C, ("i",), MESO_TARTARIC, id="isotopes-with-their-stereo"),
        pytest.param(UREA_D2, ("t",), UREA_D2, id="isotopic-hydrogen-kept"),
    ],
)
def test_strip_layers(inchi, drop, expected):
    assert strip_layers(inchi, drop) == expected


@pytest.mark.parametrize(
    ("inchi", "written"),
    [
        (HEXENEDIOL_E_AND_CENTRES, "/b4-3+ /t5-,6+ /m0 /s1"),
        (GLYCINE, "none"),
    ],
)
def test_stereo_layers_of(inchi, written):
    assert stereo_layers_of(inchi) == written


@pytest.mark.parametrize(
    ("layer", "components"),
    [
        (";1-,2-", ["", "1-,2-"]),
        # n* repeats a component
        ("2*12-;1-,2-", ["12-", "12-", "1-,2-"]),
    ],
)
def test_layer_components(layer, components):
    assert layer_components(layer) == components


@pytest.mark.parametrize(("formula", "components"), [("C4H6O6", 1), ("3ClH.6H2N.Ru", 10)])
def test_formula_components(formula, components):
    assert formula_components(formula) == components


@pytest.mark.parametrize(
    ("inchi_a", "inchi_b", "shared"),
    [
        pytest.param(BUTANE, BUT_2_ENE, True, id="hydrogens-aside"),
        pytest.param(ETHANOL, DIMETHYL_ETHER, False, id="other-connectivity"),
    ],
)
def test_share_skeleton(inchi_a, inchi_b, shared):
    assert share_skeleton(parse_standard_inchi(inchi_a), parse_standard_inchi(inchi_b)) is shared


# --- Stereo ----------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("inchi", "centres"),
    [
        # Fig. 25: the /m flag inverts the parities written in /t
        pytest.param(TARTARIC_D, {(0, "1"): "-", (0, "2"): "-"}, id="m0-as-written"),
        pytest.param(TARTARIC_L, {(0, "1"): "+", (0, "2"): "+"}, id="m1-inverted"),
        pytest.param(MESO_TARTARIC, {(0, "1"): "-", (0, "2"): "+"}, id="no-m-flag"),
        pytest.param(BUTANEDIOL_ONE_CENTRE, {(0, "3"): "-"}, id="undefined-centre-left-out"),
        # one flag per component, "." for one without stereo
        pytest.param(TEA_D_LACTATE, {(1, "2"): "+"}, id="achiral-first-component"),
        pytest.param(D_ALANINE2_SULFATE, {(0, "2"): "+", (1, "2"): "+"}, id="repeated-components"),
    ],
)
def test_defined_centres(inchi, centres):
    assert defined_centres(parse_standard_inchi(inchi)) == centres


def test_defined_double_bonds():
    assert defined_double_bonds(parse_standard_inchi(TRIENE_B)) == {(0, "7-4"): "-", (0, "8-5"): "+", (0, "9-6"): "-"}
    assert defined_double_bonds(parse_standard_inchi(HEXADIENE_ONE_DOUBLE_BOND)) == {(0, "5-3"): "+"}


@pytest.mark.parametrize(
    ("inchi", "listing"),
    [
        pytest.param(BUTANEDIOL_ONE_CENTRE, {("t", 0, "3"): True, ("t", 0, "4"): False}, id="undefined-centre-listed"),
        pytest.param(
            HEXADIENE_ONE_DOUBLE_BOND,
            {("b", 0, "5-3"): True, ("b", 0, "6-4"): False},
            id="undefined-double-bond-listed",
        ),
        # a component without a defined centre lists none of its undefined ones ("/t;1-,2-")
        pytest.param(
            PHENYLETHYLAMINE_TARTRATE, {("t", 1, "1"): True, ("t", 1, "2"): True}, id="component-listing-nothing"
        ),
    ],
)
def test_stereo_listing(inchi, listing):
    assert stereo_listing(parse_standard_inchi(inchi)) == listing


@pytest.mark.parametrize(
    ("gtp_inchi", "xref_inchi", "counts"),
    [
        # counted out of the centres either record lists
        (PHENYLETHYLAMINE_R_TARTRATE, PHENYLETHYLAMINE_TARTRATE, "stereocentres defined: GtoPdb 3 of 3, linked record 2 of 3"),
        (
            HEXENEDIOL_E_AND_CENTRES,
            HEXENEDIOL_CENTRES_ONLY_INVERTED,
            (
                "stereocentres defined: GtoPdb 2 of 2, linked record 2 of 2; double bonds defined: GtoPdb 1 of 1, "
                "linked record 0 of 1"
            ),
        ),
    ],
)
def test_stereo_definitions(gtp_inchi, xref_inchi, counts):
    assert stereo_definitions(parse_standard_inchi(gtp_inchi), parse_standard_inchi(xref_inchi)) == counts


@pytest.mark.parametrize(
    ("inchi_a", "inchi_b", "centres", "double_bonds"),
    [
        pytest.param(L_ALANINE, ALANINE, {(0, "2")}, set(), id="defined-on-one-side"),
        pytest.param(
            HEXENEDIOL_E_AND_CENTRES, HEXENEDIOL_CENTRES_ONLY_INVERTED, {(0, "5"), (0, "6")}, {(0, "4-3")}, id="both"
        ),
        pytest.param(TARTARIC_L, TARTARIC_L, set(), set(), id="same"),
    ],
)
def test_stereo_difference_keys(inchi_a, inchi_b, centres, double_bonds):
    assert stereo_difference_keys(parse_standard_inchi(inchi_a), parse_standard_inchi(inchi_b)) == (centres, double_bonds)


# --- Isotopes --------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("inchi", "written"), [(BENZENE_LABELLED, "/i1+1,4+1D"), (UREA_D2, "/i/hD2"), (BENZENE, "none")]
)
def test_isotope_layer(inchi, written):
    assert isotope_layer(parse_standard_inchi(inchi)) == written


@pytest.mark.parametrize(
    ("gtp_inchi", "xref_inchi", "difference"),
    [
        (BENZENE_LABELLED, BENZENE, "GtoPdb labelled, linked record unlabelled"),
        (BENZENE, BENZENE_LABELLED, "linked record labelled, GtoPdb unlabelled"),
        (L_ALANINE_D3, L_ALANINE_1_13C, "different labelling"),
    ],
)
def test_isotope_difference(gtp_inchi, xref_inchi, difference):
    assert isotope_difference(parse_standard_inchi(gtp_inchi), parse_standard_inchi(xref_inchi)) == difference


def test_labelled_stereo():
    # the isotopic sublayers add /m and /s to the main /t layer
    assert labelled_stereo(parse_standard_inchi(MESO_TARTARIC_13C)) == {"t": "1-,2+", "m": "1", "s": "1"}


def test_comparable_stereo_with_the_same_labels_includes_their_stereo():
    # labelling one carboxyl or the other gives the two enantiomers
    assert comparable_stereo(parse_standard_inchi(MESO_TARTARIC_13C), parse_standard_inchi(MESO_TARTARIC_13C_OTHER_END)) == (
        {"t": "1-,2+", "m": "1", "s": "1"},
        {"t": "1-,2+", "m": "0", "s": "1"},
    )


def test_comparable_stereo_with_different_labels_is_the_unlabelled_stereo():
    # read against the unlabelled meso /t, the isotopic /m1 would invert both centres
    assert comparable_stereo(parse_standard_inchi(MESO_TARTARIC_13C), parse_standard_inchi(MESO_TARTARIC)) == (
        {"t": "1-,2+"},
        {"t": "1-,2+"},
    )


# --- Comparing two records layer by layer ----------------------------------------------------------------------------

UP_TO_STEREO = ["formula", "connections", "hydrogens"]


@pytest.mark.parametrize(
    ("inchi_a", "inchi_b", "shared", "first"),
    [
        pytest.param(L_ALANINE, L_ALANINE, [*UP_TO_STEREO, "stereo"], None, id="identical"),
        pytest.param(BUTANE, BUT_2_ENE, [], "formula", id="formula"),
        pytest.param(ETHANOL, DIMETHYL_ETHER, ["formula"], "connections", id="connections"),
        # a layer neither record has is not listed as shared: glycine and glycinate have no /q layer
        pytest.param(GLYCINE, GLYCINATE, UP_TO_STEREO, "protons", id="protons"),
        pytest.param(L_ALANINE, ALANINE, UP_TO_STEREO, "stereo", id="stereo"),
        # the comparison stops at the first difference, so the isotopes that also differ are not reached
        pytest.param(L_ALANINE_D3, ALANINE, UP_TO_STEREO, "stereo", id="stops-at-the-first-difference"),
        pytest.param(L_ALANINE_D3, L_ALANINE_1_13C, [*UP_TO_STEREO, "stereo"], "isotopes", id="isotopes"),
        # with the same labels, the stereo they create is compared as stereo
        pytest.param(MESO_TARTARIC_13C, MESO_TARTARIC_13C_OTHER_END, UP_TO_STEREO, "stereo", id="labelled-stereo"),
        # with different labels, the unlabelled stereo is compared, and is the same
        pytest.param(MESO_TARTARIC_13C, MESO_TARTARIC, [*UP_TO_STEREO, "stereo"], "isotopes", id="unlabelled-stereo"),
    ],
)
def test_compare_layers(inchi_a, inchi_b, shared, first):
    assert compare_layers(parse_standard_inchi(inchi_a), parse_standard_inchi(inchi_b)) == (shared, first)


def test_inchi_parts_are_the_layers_as_written():
    parts = inchi_parts(MESO_TARTARIC_13C)
    assert [key for key, _ in parts] == ["formula", "c", "h", "t", "i", "im", "is"]
    assert "InChI=1S/" + "/".join(text for _, text in parts) == MESO_TARTARIC_13C


@pytest.mark.parametrize(
    ("inchi_a", "inchi_b", "keys"),
    [
        pytest.param(L_ALANINE, L_ALANINE, set(), id="identical"),
        # every later layer differs too, but only the first is given
        pytest.param(BUTANE, BUT_2_ENE, {"formula"}, id="formula"),
        pytest.param(ETHANOL, DIMETHYL_ETHER, {"c"}, id="connections"),
        pytest.param(ACETONE, PROPEN_2_OL, {"h"}, id="hydrogens"),
        pytest.param(SODIUM, SODIUM_ION, {"q"}, id="charge"),
        pytest.param(GLYCINE, GLYCINATE, {"p"}, id="protons"),
        # only the stereo layers that differ
        pytest.param(TARTARIC_L, TARTARIC_D, {"m"}, id="inversion-flag"),
        pytest.param(TRIENE_A, TRIENE_B, {"b"}, id="double-bond"),
        pytest.param(L_ALANINE, ALANINE, {"t", "m", "s"}, id="stereo-on-one-side"),
        # with the same labels, the stereo they create is compared as stereo
        pytest.param(MESO_TARTARIC_13C, MESO_TARTARIC_13C_OTHER_END, {"im"}, id="labelled-stereo"),
        pytest.param(L_ALANINE_D3, L_ALANINE_1_13C, {"i"}, id="isotopes"),
        pytest.param(UREA_D2, UREA, {"i", "ih"}, id="isotopic-hydrogens"),
    ],
)
def test_first_difference_keys(inchi_a, inchi_b, keys):
    assert first_difference_keys(parse_standard_inchi(inchi_a), parse_standard_inchi(inchi_b)) == keys


# --- Charge and hydrogens --------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("inchi", "charge"),
    [
        (GLYCINATE, -1),
        # /q;;2*-1;+4 counts the repeated NH2- twice, and /p-2 removes the protons of 2 HCl
        (CISPLATIN, 0),
        # charges only, no /p layer
        (RU360_GTOPDB, 5),
    ],
)
def test_net_charge(inchi, charge):
    assert net_charge(parse_standard_inchi(inchi)) == charge


@pytest.mark.parametrize(
    ("gtp_inchi", "xref_inchi", "difference"),
    [
        (GLYCINE, GLYCINATE, "net charge GtoPdb 0, linked record -1; GtoPdb has 1 proton more"),
        (GLYCINATE, GLYCINIUM, "net charge GtoPdb -1, linked record +1; GtoPdb has 2 protons fewer"),
        # protons are only counted between equal formulas
        (CISPLATIN, HEXAAMMINE_RUTHENIUM_CHLORIDE, "net charge GtoPdb 0, linked record -3"),
        (RU360_GTOPDB, RU360_PUBCHEM, "same net charge, different component charge layers"),
        (GLYCINE, GLYCINE, ""),
    ],
)
def test_charge_difference(gtp_inchi, xref_inchi, difference):
    assert charge_difference(parse_standard_inchi(gtp_inchi), parse_standard_inchi(xref_inchi)) == difference


@pytest.mark.parametrize(
    ("inchi", "hydrogens"),
    [
        # "/h3-4H,1-2H3,(H,10,11,13,14)": ranges of atoms, and a mobile H shared by four atoms, left out
        pytest.param(CATION_16C, {(0, 3): 1, (0, 4): 1, (0, 1): 3, (0, 2): 3}, id="ranges"),
        # "/h2*2H,4H2,1H3,(H,5,6);(H2,1,2,3,4)": two alanines, and a sulfuric acid with only mobile H
        pytest.param(
            L_ALANINE2_SULFATE,
            {(0, 2): 1, (0, 4): 2, (0, 1): 3, (1, 2): 1, (1, 4): 2, (1, 1): 3},
            id="repeated-components",
        ),
    ],
)
def test_fixed_hydrogens(inchi, hydrogens):
    assert fixed_hydrogens(parse_standard_inchi(inchi)) == hydrogens
