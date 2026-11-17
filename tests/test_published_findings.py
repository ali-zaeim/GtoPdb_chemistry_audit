"""Reconcile the reported findings with the rows shipped for the same release."""

import json
from pathlib import Path

import pandas as pd
import pytest

RESULTS = Path(__file__).resolve().parents[1] / "data/processed_data/gtp_xref_audit"


@pytest.mark.parametrize("summary_file", sorted(RESULTS.glob("*__summary.json")), ids=lambda path: path.stem)
def test_published_findings(summary_file):
    summary = json.loads(summary_file.read_text())
    prefix = summary_file.name.removesuffix("summary.json")
    xrefs = pd.read_csv(RESULTS / f"{prefix}xref_structure_comparison.csv", low_memory=False)
    comparable = xrefs.gtp_inchi.notna() & xrefs.xref_inchi.notna()
    exact = comparable & xrefs.gtp_inchi.eq(xrefs.xref_inchi)
    assert len(xrefs) == summary["xrefs_total"]
    assert comparable.sum() == summary["xrefs_compared"]
    assert exact.sum() == summary["xrefs_identical"]
    different = comparable & ~exact
    assert (different & xrefs.metal_differs).sum() == summary["metal_compared"]
    assert different.sum() == sum(summary["xrefs_by_first_difference"].values()) + summary["metal_compared"]

    # Count mappings directly from the exported IDs, including stereochemistry in isotopic sublayers.
    stereo_links = xrefs[exact & xrefs.gtp_inchi.str.contains(r"/[btms]", na=False)]
    unichem = pd.read_csv(RESULTS / f"{prefix}unichem_gtopdb_structures.csv", dtype=str).set_index("ligand_id")
    mapped = 0
    for link in stereo_links.itertuples():
        ids = unichem.at[str(link.ligand_id), f"linked_{link.xref_db.lower()}"]
        if isinstance(ids, str) and link.xref_id_normalised in ids.split("|"):
            mapped += 1
    assert len(stereo_links) == summary["stereo_exact_links"]
    assert mapped == summary["stereo_exact_links_in_unichem"]
