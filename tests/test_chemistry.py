"""Tests for the notebook's *RDKit helpers* functions: the RDKit structure helpers.

The InChIs are named in `tests/inchis.py`, which gives their sources.
"""

from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import inchi as rdinchi

import gtp_chemistry_problem as chemistry
from gtp_chemistry_problem import (
    common_substructure,
    connectivity_difference,
    fragment_smiles,
    has_metal,
    heavy_and_hydrogen_count,
    inchi_matching_key,
    inorganic_smiles,
    is_radioisotope,
    isotopes_written,
    missing_stereo,
    mol_from_inchi,
    morgan_similarity,
    parent_inchi,
    radical_count,
    ring_systems,
)
from tests.inchis import (
    ACETIC_ACID,
    ACETONE,
    ALANINE,
    ANISOLE,
    AZ2158_CHEMBL,
    AZ2158_GTOPDB,
    BENZYL_ALCOHOL,
    BUTANEDIOL_ONE_CENTRE,
    BUTANEDIOL_TWO_CENTRES,
    CAFFEINE,
    CALCIUM_ACETATE,
    CISPLATIN,
    DIMETHYL_ETHER,
    ELAMIPRETIDE,
    ETHANE,
    ETHANOL,
    ETHANOL_1_D,
    ETHYL_RADICAL,
    GLYCINE,
    HEXENEDIOL_CENTRES_ONLY,
    HEXENEDIOL_E_AND_CENTRES,
    HEXENEDIOL_E_ONLY,
    L_ALANINE,
    L_ALANINE_D3,
    L_ALANINE_HCL,
    O_CHLOROTOLUENE,
    OXALIPLATIN,
    P_CHLOROTOLUENE,
    PROPEN_2_OL,
    S_ETHANOL_1_D,
    SODIUM_ACETATE,
    THEOPHYLLINE,
    YTTRIUM,
)

NAN = float("nan")

# --- Parsing ---------------------------------------------------------------------------------------------------------


def test_mol_from_inchi():
    assert mol_from_inchi(GLYCINE).GetNumAtoms() == 5


@pytest.mark.parametrize("inchi", ["", NAN, "InChI=1S/garbage"])
def test_mol_from_inchi_missing_or_unparsable(inchi):
    assert mol_from_inchi(inchi) is None


@pytest.mark.parametrize(("smiles", "count"), [("CCO", 9), ("C(C", None), (NAN, None)])
def test_heavy_and_hydrogen_count(smiles, count):
    assert heavy_and_hydrogen_count(smiles) == count


@pytest.mark.parametrize(
    ("smiles", "inchi"),
    [
        pytest.param("CCO", ETHANOL, id="matches"),
        pytest.param("COC", None, id="other-structure"),
        pytest.param("C(C", None, id="unparsable"),
    ],
)
def test_inchi_matching_key(smiles, inchi):
    assert inchi_matching_key(smiles, rdinchi.InchiToInchiKey(ETHANOL)) == inchi


@pytest.mark.parametrize(
    ("name", "smiles"),
    [
        ("Ca2+", "[Ca+2]"),
        ("Cl-", "[Cl-]"),
        ("CO", "[C-]#[O+]"),
        # the charge comes after the count in a GtoPdb name
        ("Ca+2", None),
    ],
)
def test_inorganic_smiles(name, smiles):
    assert inorganic_smiles(name) == smiles


@pytest.mark.parametrize(
    ("inchi", "metal"),
    [
        pytest.param(CISPLATIN, True, id="platinum"),
        # alkali and alkaline-earth metals are salt counterions
        pytest.param(SODIUM_ACETATE, False, id="sodium"),
        pytest.param(CALCIUM_ACETATE, False, id="calcium"),
    ],
)
def test_has_metal(inchi, metal):
    assert has_metal(inchi) is metal


# --- Standardised parents --------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("inchi", "options", "parent"),
    [
        pytest.param(L_ALANINE_HCL, {}, ALANINE, id="salt-without-stereo"),
        pytest.param(L_ALANINE_HCL, {"keep_stereo": True}, L_ALANINE, id="salt-with-stereo"),
        pytest.param(L_ALANINE_D3, {"keep_stereo": True}, L_ALANINE, id="without-isotopes"),
        pytest.param(SODIUM_ACETATE, {}, ACETIC_ACID, id="sodium-salt"),
        pytest.param(PROPEN_2_OL, {"canonical_tautomer": True}, ACETONE, id="canonical-tautomer"),
        pytest.param(PROPEN_2_OL, {}, PROPEN_2_OL, id="tautomer-kept"),
        # with the metal disconnected, the largest fragment would be a ligand
        pytest.param(OXALIPLATIN, {}, None, id="metal-complex"),
        pytest.param(None, {}, None, id="missing"),
    ],
)
def test_parent_inchi(inchi, options, parent):
    assert parent_inchi(inchi, **options) == parent


# --- Stereo one record lacks -----------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("more_inchi", "less_inchi", "kinds"),
    [
        pytest.param(BUTANEDIOL_TWO_CENTRES, BUTANEDIOL_ONE_CENTRE, "R/S centre", id="centre"),
        pytest.param(HEXENEDIOL_E_AND_CENTRES, HEXENEDIOL_CENTRES_ONLY, "E/Z", id="double-bond-next-to-centres"),
        pytest.param(HEXENEDIOL_E_AND_CENTRES, HEXENEDIOL_E_ONLY, "R/S centre", id="centres-next-to-a-double-bond"),
        # AZ2158: C19 of the 1,4-disubstituted cyclohexane has two symmetry-equivalent ring neighbours
        pytest.param(AZ2158_CHEMBL, AZ2158_GTOPDB, "ring cis/trans", id="ring-cis-trans"),
        # the centre that the deuterium creates is in the isotopic sublayers, compared between the same labels
        pytest.param(S_ETHANOL_1_D, ETHANOL_1_D, "R/S centre", id="centre-created-by-a-label"),
        pytest.param(L_ALANINE, L_ALANINE, "", id="none"),
    ],
)
def test_missing_stereo(more_inchi, less_inchi, kinds):
    assert missing_stereo(more_inchi, less_inchi) == kinds


# --- Similarity, radicals and isotopes -------------------------------------------------------------------------------


def test_morgan_similarity():
    assert morgan_similarity(GLYCINE, GLYCINE) == 1.0
    # rounded to 3 decimals
    assert morgan_similarity(CAFFEINE, THEOPHYLLINE) == 0.378


def test_morgan_similarity_unparsable():
    assert morgan_similarity(GLYCINE, "InChI=1S/garbage") is None
    assert morgan_similarity("InChI=1S/garbage", GLYCINE) is None


@pytest.mark.parametrize(
    ("structure", "radicals"),
    [
        (ETHYL_RADICAL, 1),
        ("[CH2][CH2]", 2),
        # unparsable
        ("C((", 0),
    ],
)
def test_radical_count(structure, radicals):
    assert radical_count(structure) == radicals


@pytest.mark.parametrize(
    ("smiles", "name_html", "isotopes"),
    [
        pytest.param("[2H]C([2H])([2H])[C@H](N)C(=O)O", "alanine-d3", ["2H"], id="in-the-smiles"),
        pytest.param("CC(N)C(=O)O", "[<sup>13</sup>C,<sup>15</sup>N]alanine", ["13C", "15N"], id="in-the-name"),
        pytest.param(None, "Phe(3<sup>125</sup>I)", ["125I"], id="no-smiles"),
        pytest.param("CC", "<sup>11</sup>C-raclopride", ["11C"], id="start-of-the-name"),
        # superscripts after a residue or before a residue name are sequence positions
        pytest.param("CC", "[Thr<sup>4</sup>]oxytocin", [], id="residue-position"),
        pytest.param("CC", "[<sup>2</sup>Ala]GLP-1", [], id="position-before-a-residue"),
    ],
)
def test_isotopes_written(smiles, name_html, isotopes):
    assert isotopes_written(smiles, name_html) == isotopes


@pytest.mark.parametrize(
    ("isotope", "radioactive"),
    [
        ("3H", True),
        ("2H", False),
        # not an element symbol
        ("11Xx", None),
        ("999C", None),
        ("40K", True),
        ("39K", False),
    ],
)
def test_is_radioisotope(isotope, radioactive):
    assert is_radioisotope(isotope) is radioactive


# --- Connectivity ----------------------------------------------------------------------------------------------------


def test_ring_systems():
    # naphthalene (fused) and a cyclopropane (separate)
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1CC1CC1")
    assert sorted(map(sorted, ring_systems(mol))) == [list(range(10)), [11, 12, 13]]
    assert ring_systems(Chem.MolFromSmiles("CCO")) == []


@pytest.mark.parametrize(
    ("atoms", "smiles"),
    [
        pytest.param({0, 7}, "*C.*Cl", id="two-pieces"),
        pytest.param(set(range(8)), "Cc1ccccc1Cl", id="whole-molecule"),
    ],
)
def test_fragment_smiles(atoms, smiles):
    assert fragment_smiles(Chem.MolFromSmiles("Cc1ccccc1Cl"), atoms) == smiles


def test_common_substructure():
    mols, pattern = common_substructure(O_CHLOROTOLUENE, P_CHLOROTOLUENE)
    assert [mol.GetNumAtoms() for mol in mols] == [8, 8]
    # toluene: the chlorine is bonded to a different ring atom in each
    assert pattern.GetNumAtoms() == 7
    assert all(mol.HasSubstructMatch(pattern) for mol in mols)


@pytest.mark.parametrize(
    ("inchi_a", "inchi_b", "options"),
    [
        pytest.param(GLYCINE, None, {}, id="unparsable"),
        pytest.param(YTTRIUM, ETHANE, {}, id="nothing-in-common"),
        pytest.param(ELAMIPRETIDE, ELAMIPRETIDE, {"max_atoms": 45}, id="too-large"),
    ],
)
def test_common_substructure_not_found(inchi_a, inchi_b, options):
    assert common_substructure(inchi_a, inchi_b, **options) is None


def test_common_substructure_timed_out(monkeypatch):
    # a search cancelled at its timeout returns the best match so far, which may be far from the largest
    cancelled = SimpleNamespace(numAtoms=3, canceled=True, smartsString="[#6]-[#6]-[#6]")
    monkeypatch.setattr(chemistry.rdFMCS, "FindMCS", lambda *args, **kwargs: cancelled)
    assert common_substructure(O_CHLOROTOLUENE, P_CHLOROTOLUENE) is None


@pytest.mark.parametrize(
    ("inchi_a", "inchi_b", "difference"),
    [
        # the chlorine alone looks the same on both sides, so it is shown with the ring it is attached to
        pytest.param(
            O_CHLOROTOLUENE, P_CHLOROTOLUENE, "GtoPdb *c1ccccc1Cl; linked record *c1ccc(Cl)cc1", id="positional-isomer"
        ),
        # the parts outside the common phenyl already differ: hydroxymethyl against methoxy
        pytest.param(BENZYL_ALCOHOL, ANISOLE, "GtoPdb *CO; linked record *OC", id="other-substituent"),
        pytest.param(ETHANOL, DIMETHYL_ETHER, None, id="not-isolated"),
        pytest.param(YTTRIUM, ETHANE, None, id="nothing-in-common"),
    ],
)
def test_connectivity_difference(inchi_a, inchi_b, difference):
    assert connectivity_difference(inchi_a, inchi_b) == difference
