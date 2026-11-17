"""Tests for the notebook's *UniChem checks* functions: how UniChem holds GtoPdb's structures.

The InChIs are named in `tests/inchis.py`, which gives their sources.
"""

import pytest

from gtp_chemistry_problem import unichem_status, unichem_structure_source
from tests.inchis import (
    ALANINE,
    ALANINE_HCL,
    D_ALANINE,
    L_ALANINE,
    L_ALANINE_HCL,
)


@pytest.mark.parametrize(
    ("gtp_inchi", "unichem_inchi", "status"),
    [
        (L_ALANINE, float("nan"), "not in UniChem"),
        (L_ALANINE, L_ALANINE, "identical"),
        (L_ALANINE, ALANINE, "stereo layers missing in UniChem"),
        (L_ALANINE, D_ALANINE, "other difference"),
    ],
)
def test_unichem_status(gtp_inchi, unichem_inchi, status):
    assert unichem_status(gtp_inchi, unichem_inchi) == status


@pytest.mark.parametrize(
    ("gtp_inchi", "smiles_inchi", "unichem_inchi", "source"),
    [
        # the SMILES lacks the salt the InChI has; UniChem holds either one without its stereo
        pytest.param(L_ALANINE_HCL, L_ALANINE, ALANINE, "SMILES", id="smiles"),
        pytest.param(L_ALANINE_HCL, L_ALANINE, ALANINE_HCL, "InChI", id="inchi"),
        # stereo is compared out, so UniChem's stereo does not matter
        pytest.param(L_ALANINE_HCL, L_ALANINE, L_ALANINE, "SMILES", id="smiles-with-stereo"),
        # the SMILES lacks only stereo: both give the same structure without it
        pytest.param(L_ALANINE, ALANINE, ALANINE, "either (the two differ only in stereo)", id="either"),
        pytest.param(L_ALANINE_HCL, L_ALANINE, "InChI=1S/CH4/h1H4", "neither", id="neither"),
        pytest.param(L_ALANINE, None, ALANINE, "InChI (the SMILES is invalid)", id="invalid-smiles"),
        pytest.param(L_ALANINE, None, ALANINE_HCL, "neither (the SMILES is invalid)", id="invalid-smiles-neither"),
        pytest.param(L_ALANINE, ALANINE, None, "not in UniChem", id="not-held"),
    ],
)
def test_unichem_structure_source(gtp_inchi, smiles_inchi, unichem_inchi, source):
    assert unichem_structure_source(gtp_inchi, smiles_inchi, unichem_inchi) == source
