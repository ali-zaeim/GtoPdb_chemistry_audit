# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo>=0.25.2",
#     "pandas>=3.1.0",
#     "plotly>=7.1.0",
#     "pyarrow>=26.0.0",
#     "rdkit>=2026.9.1",
#     "requests>=2.34.2",
#     "urllib3>=2.8.0",
# ]
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(
    width="full",
    app_title="GtoPdb cross-reference structure audit",
    css_file="gtp_chemistry_problem.css",
)

with app.setup:
    import datetime as dt
    import hashlib
    import html
    import json
    import re
    import time
    import unicodedata
    from concurrent.futures import ThreadPoolExecutor
    from enum import StrEnum
    from pathlib import Path
    from typing import NamedTuple

    import marimo as mo
    import pandas as pd
    import plotly.graph_objects as go
    import requests
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import inchi as rdinchi
    from rdkit.Chem import rdDepictor, rdFingerprintGenerator, rdFMCS
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    RDLogger.DisableLog("rdApp.*")

    # constants used by the notebook and by the functions in the appendix
    # Ground states marked stbl in NUBASE2020, Kondev et al., Chin. Phys. C 45, 030001 (2021).
    STABLE_NUCLIDES = frozenset(
        "1H 2H 3He 4He 6Li 7Li 9Be 10B 11B 12C 13C 14N 15N 16O 17O 18O 19F 20Ne 21Ne 22Ne 23Na 24Mg 25Mg 26Mg "  # noqa: SIM905
        "27Al 28Si 29Si 30Si 31P 32S 33S 34S 35Cl 36S 36Ar 37Cl 38Ar 39K 40Ar 40Ca 41K 42Ca 43Ca 44Ca 45Sc "
        "46Ca 46Ti 47Ti 48Ti 49Ti 50Ti 50Cr 51V 52Cr 53Cr 54Cr 54Fe 55Mn 56Fe 57Fe 58Fe 58Ni 59Co 60Ni 61Ni "
        "62Ni 63Cu 64Ni 64Zn 65Cu 66Zn 67Zn 68Zn 69Ga 70Zn 70Ge 71Ga 72Ge 73Ge 74Ge 74Se 75As 76Se 77Se 78Se "
        "78Kr 79Br 80Se 80Kr 81Br 82Kr 83Kr 84Kr 84Sr 85Rb 86Kr 86Sr 87Sr 88Sr 89Y 90Zr 91Zr 92Zr 92Mo 93Nb "
        "94Zr 94Mo 95Mo 96Mo 96Ru 97Mo 98Mo 98Ru 99Ru 100Ru 101Ru 102Ru 102Pd 103Rh 104Ru 104Pd 105Pd 106Pd "
        "106Cd 107Ag 108Pd 108Cd 109Ag 110Pd 110Cd 111Cd 112Cd 112Sn 113In 114Cd 114Sn 115Sn 116Sn 117Sn "
        "118Sn 119Sn 120Sn 120Te 121Sb 122Sn 122Te 123Sb 123Te 124Sn 124Te 124Xe 125Te 126Te 126Xe 127I 128Xe "
        "129Xe 130Xe 130Ba 131Xe 132Xe 132Ba 133Cs 134Xe 134Ba 135Ba 136Ba 136Ce 137Ba 138Ba 138Ce 139La "
        "140Ce 141Pr 142Ce 142Nd 143Nd 144Sm 145Nd 146Nd 148Nd 149Sm 150Sm 152Sm 153Eu 154Sm 154Gd 155Gd "
        "156Gd 156Dy 157Gd 158Gd 158Dy 159Tb 160Gd 160Dy 161Dy 162Dy 162Er 163Dy 164Dy 164Er 165Ho 166Er "
        "167Er 168Er 168Yb 169Tm 170Er 170Yb 171Yb 172Yb 173Yb 174Yb 175Lu 176Yb 176Hf 177Hf 178Hf 179Hf "
        "180Hf 181Ta 182W 183W 184W 185Re 186W 187Os 188Os 189Os 190Os 191Ir 192Os 192Pt 193Ir 194Pt 195Pt "
        "196Pt 196Hg 197Au 198Pt 198Hg 199Hg 200Hg 201Hg 202Hg 203Tl 204Hg 204Pb 205Tl 206Pb 207Pb 208Pb "
        .split()
    )
    STEREO_LAYERS = ("b", "t", "m", "s")
    STEREO_LABELS = {True: "GtoPdb structure has stereo", False: "no stereo"}
    # `stereo_relation` values that make a conflict. They are only given when both records define the same
    # stereocentres and double bonds
    STEREO_CONFLICT_RELATIONS = frozenset(
        {
            "enantiomer",
            "epimer or diastereomer",
            "E/Z isomer",
            # every defined centre inverted, but the same stereocentres or double bonds undefined on both sides
            "enantiomer or diastereomer",
        }
    )
    # `stereo_relation` values for stereo defined on one side only -> the record that defines more
    MORE_STEREO_IN = {
        "GtoPdb defines more stereo": "GtoPdb",
        "linked record has no stereo": "GtoPdb",
        "linked record defines more stereo": "linked record",
        "GtoPdb structure has no stereo": "linked record",
    }
    # `stereo_relation` when each side defines stereo that the other leaves undefined, or one defines more and a centre
    # both define is written inverted
    STEREO_DEFINED_DIFFERENTLY = "stereo defined differently on each side"
    # the layers of a standard InChI, in the order it writes them; two records are compared layer by layer up to the first
    # that differs
    LAYERS = ("formula", "connections", "hydrogens", "charge", "protons", "stereo", "isotopes")
    # non-metals, metalloids, and the alkali and alkaline-earth metals (salt counterions); any other makes a metal complex
    SALT_AND_NONMETAL_ELEMENTS = frozenset(
        "H He B C N O F Ne Si P S Cl Ar Ge As Se Br Kr Sb Te I Xe At Rn "  # noqa: SIM905 (easier to read as one string)
        "Li Na K Rb Cs Fr Be Mg Ca Sr Ba Ra".split()
    )
    XREF_COLUMNS = {"ChEMBL": "ChEMBl ID", "ChEBI": "Chebi ID", "PubChem": "PubChem CID"}
    XREF_PATTERNS = {"ChEMBL": r"CHEMBL\d+", "ChEBI": r"CHEBI:\d+", "PubChem": r"[1-9]\d{0,9}"}
    XREF_URLS = {
        "ChEMBL": "https://www.ebi.ac.uk/chembl/explore/compound/{}",
        "ChEBI": "https://www.ebi.ac.uk/chebi/{}",
        "PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/{}",
    }
    GTP_LIGAND_URL = "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId={}"
    # CSS for the text cells of `linked_table`, in camelCase as `mo.ui.table(style_cell=...)` takes it
    TABLE_CELL_STYLE = {"padding": "6px 16px 6px 8px", "lineHeight": "1.35"}
    # `linked_table` column widths in pixels, the same in every table; other text columns are sized by their content
    COLUMN_WIDTHS = {
        # who and what is linked
        "ligand_id": 60,
        "name": 110,
        "type": 110,
        "xref_db": 60,
        "xref_id_as_written": 110,
        "xref_name": 160,
        "chebi_stars": 60,
        "radioactive": 80,
        # how the structures differ
        "identifier_issue": 160,
        "difference": 180,
        "difference_detail": 200,
        "stereo_relation": 150,
        "missing_stereo": 120,
        "pubchem_standardisation_removes": 120,
        "tanimoto_similarity": 90,
        "external databases agree": 110,
        # what to do
        "action": 110,
        "suggested_action": 280,
        "xref_id_with_gtp_structure": 150,
        "chembl_parent_id": 120,
        "other_gtp_ligands_with_xref_structure": 110,
        # structures as text
        "inchis": 300,
        "smiles": 220,
        "inchi": 220,
        "gtp_inchi": 220,
        "xref_inchi": 220,
        "inchikey": 150,
    }
    UNICHEM_SOURCE_IDS = {"GtoPdb": 4, "ChEMBL": 1, "ChEBI": 7, "PubChem": 22}
    # the cache committed to the GitHub repository, copied where the local cache is missing (section 1), for example
    # when the notebook runs on molab without the rest of the repository
    GITHUB_CACHE_URL = (
        "https://raw.githubusercontent.com/ali-zaeim/GtoPdb_chemistry_audit/main/data/raw_data/gtp_xref_audit/"
    )
    CACHE_FILES = (
        "ligands.csv",
        "ligands.csv.log.json",
        "ligand_id_mapping.csv",
        "ligand_id_mapping.csv.log.json",
        "unichem_src1src4.txt.gz",
        "unichem_src1src4.txt.gz.log.json",
        "unichem_src4src7.txt.gz",
        "unichem_src4src7.txt.gz.log.json",
        "unichem_src4src22.txt.gz",
        "unichem_src4src22.txt.gz.log.json",
        "unichem_gtopdb_compounds.parquet",
        "unichem_compound_ChEMBL_CHEMBL1069.json",
        "chembl_molecules.parquet",
        "chembl_id_lookup.parquet",
        "chebi_compounds.parquet",
        "pubchem_compounds.parquet",
        "pubchem_sid_cids.parquet",
        "chembl_by_inchikey.parquet",
        "chebi_by_inchikey.parquet",
        "pubchem_by_inchikey.parquet",
        "pubchem_standardised.parquet",
        "chembl_suggested_records.parquet",
        "pubchem_inorganic_by_inchikey.parquet",
    )
    # where a table draws its structures, the drawings' captions give these
    SHOWN_IN_DRAWINGS = ("ligand_id", "name", "xref_db", "xref_id_as_written", "xref_name", "gtp_inchi", "xref_inchi")
    # the records that `suggested_action` names, with links
    NAMED_IN_ACTION = ("xref_id_with_gtp_structure", "chembl_parent_id", "other_gtp_ligands_with_xref_structure")
    # the IDs a suggested action names: ChEMBL and ChEBI IDs, 'PubChem CID n', and 'GtoPdb ligand n' or 'ligands n, m'
    ACTION_IDS = re.compile(r"(CHEMBL\d+)|(CHEBI:\d+)|(PubChem CID )(\d+)|(GtoPdb ligands? )(\d+(?:, \d+)*)")


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Guide to Pharmacology ligand cross-references: structure discrepancies

    This notebook checks the ChEMBL, ChEBI and PubChem identifiers that the IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb) assigns to its ligands. For each cross-reference, it compares the standard InChI of the GtoPdb ligand with the one held by the linked database, and suggests an action for curators where they differ. It also explains why GtoPdb ligands do not match their ChEMBL, ChEBI and PubChem records in UniChem.
    """)
    return


@app.cell(hide_code=True)
def _(gtp_version, summary):
    layer_counts = summary["xrefs_by_first_difference"]
    mo.md(f"""
    ## Summary

    *{gtp_version}; latest source retrieval {summary["retrieved_on"]}.*

    **Cross-references** (sections 3 to 5). An *exact record* is a record in the linked database with exactly the GtoPdb structure: the one to link if the GtoPdb structure is right.

    - **{summary["xrefs_total"]:,}** ChEMBL, ChEBI and PubChem cross-references on **{summary["ligands_with_xrefs"]:,}** ligands were checked. **{summary["xrefs_compared"]:,}** could be compared, and **{summary["xrefs_identical"]:,}** of these match the GtoPdb structure exactly.
    - **{summary["identifier_issues"]:,}** identifiers have a problem of their own: malformed, retired or secondary, or a ChEMBL salt or child form where the parent may be intended (section 4). ChEMBL has an exact record for {summary["chembl_unresolved_with_suggestion"]:,} of the {summary["chembl_unresolved"]:,} ChEMBL IDs that do not resolve to a compound.
    - **{summary["formula_or_connections"]:,}** cross-references on {summary["formula_or_connections_ligands"]:,} ligands differ in formula or connections ({summary["formula_or_connections_by_db"]}; sections 5.1 and 5.2). An exact record exists for {summary["formula_or_connections_with_gtp_structure_id"]:,} of them.
    - **{layer_counts["stereo"]:,}** first differ in stereo (section 5.5): **{summary["stereo_conflicts"]:,}** have an opposite configuration, **{summary["stereo_defined_differently"]:,}** define stereo differently on each side, and **{summary["stereo_one_side"]:,}** define it on one side only. An exact record exists for {summary["stereo_conflicts_with_gtp_stereo_id"]:,} of the opposite configurations, and for {summary["stereo_better_record"]:,} of the one-sided ones.
    - **{layer_counts["hydrogens"]:,}** first differ in hydrogens, **{layer_counts["charge"] + layer_counts["protons"]:,}** in charge or protons, and **{layer_counts["isotopes"]:,}** in isotopes (sections 5.3, 5.4 and 5.6).
    - **{summary["metal_differs"]:,}** cross-references on {summary["metal_differs_ligands"]:,} ligands involve a metal-containing structure and require separate review. In **{summary["metal_compared"]}** the InChIs differ. The other **{summary["metal_differs"] - summary["metal_compared"]}** lack a linked InChI (section 5.7).
    - On **{summary["consensus_ligands"]:,}** ligands, all comparable records from at least two external databases share an InChI that differs from GtoPdb’s. These ligands merit review against the original evidence.
    - **{summary["pubchem_agree_from_gtp_deposit"]:,}** of the {summary["pubchem_agree"]:,} PubChem cross-references that agree with GtoPdb point to the compound that GtoPdb's own deposit is merged into, so the CID match alone does not independently confirm the GtoPdb structure (section 3). On {summary["ligands_only_confirmed_by_gtp_deposit"]:,} ligands, this is the only agreeing record while another database disagrees.

    **GtoPdb's own structures** (section 6)

    - **{summary["radioactive_without_isotope"]:,}** ligands flagged as radioactive have an InChI with no isotope layer. Across all ligands, **{summary["duplicate_keys"]:,}** InChIKeys occur more than once. **{summary["duplicate_keys_with_radioactive_ligand"]:,}** of those groups include a ligand flagged as radioactive. The isotope check also identifies {summary["unflagged_radioactive"]:,} unflagged ligands with a radioisotope and {summary["radioactive_with_stable_isotopes_only"]:,} flagged ligands whose detected labels are all stable (section 6.1).
    - **{summary["smiles_inchi_issues"]:,}** ligands have a SMILES that is invalid ({summary["invalid_smiles"]:,}) or does not reproduce their own InChI. The website draws structures from the SMILES (section 6.4).
    - **{summary["inchi_rebuilt"]:,}** large peptides and oligonucleotides have an InChIKey but no InChI, so the InChI was rebuilt from the SMILES (section 6.3).
    - **{summary["inorganic_without_structure"]:,}** inorganic ligands have no SMILES, InChI or InChIKey. For **{summary["inorganic_structures"]:,}** of them, such as Na<sup>+</sup>, Ca<sup>2+</sup>, Cl<sup>-</sup> and NO, the name is the formula, so section 6.5 gives the structure. {summary["inorganic_pubchem_link_changes"]:,} of these have a suggested PubChem link added or replaced.

    **UniChem** (section 7)

    - UniChem's own records (from its API) hold **{summary["unichem_with_stereo"]:,}** GtoPdb ligands that have stereochemistry, and {"**none** of them keeps" if summary["unichem_with_stereo_kept"] == 0 else f"only **{summary['unichem_with_stereo_kept']:,}** of them keep"} any of it: **{summary["unichem_stereo_stripped"]:,}** are GtoPdb's structure with every stereo layer removed. Of **{summary["stereo_exact_links"]:,}** cross-references with identical stereo-defined InChIs, **{summary["stereo_exact_links_in_unichem"]:,}** appear in the UniChem mapping files. UniChem holds **{summary["unichem_ucis_merging_stereoisomers"]:,}** groups of GtoPdb ligands differing only in stereo specification under shared UCIs (sections 7 and 7.2).
    - Where GtoPdb's SMILES and InChI differ in more than stereo, UniChem holds the SMILES structure for {summary["unichem_structure_follows_smiles"]} of {summary["unichem_structure_decisive"]} ligands (section 7.1). Whether the stereo is missing from the file GtoPdb supplies or lost when UniChem loads it cannot be seen from public data; GtoPdb can check the file it sends.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### How to reproduce

    Everything comes from public sources. GtoPdb data come from its download files. UniChem data come from its flat-file mappings and from its API, which gives the structure UniChem holds for each GtoPdb ligand (section 7). ChEMBL, ChEBI and PubChem data come from their public APIs.

    Downloads and API responses are cached next to the notebook under `data/raw_data/gtp_xref_audit/`, so later runs only fetch identifiers that are not cached yet. Where that cache is missing, for example when the notebook runs on molab without the rest of its repository, it first copies the cache committed to the [GitHub repository](https://github.com/ali-zaeim/GtoPdb_chemistry_audit) (about 17 MB), so it reproduces the published results without querying the APIs. Tick **Re-download all source data** only if you need to start from scratch: that downloads about 10 MB from GtoPdb and 130 KB of UniChem mapping files, then makes a few hundred ChEMBL, ChEBI and PubChem requests (a few minutes) and one UniChem request per GtoPdb ligand (about 14,000; around 10 minutes). Result tables are written to `data/processed_data/gtp_xref_audit/` (section 8).
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## 1. Data retrieval
    """)
    return


@app.cell
def _(refresh):
    notebook_dir = mo.notebook_dir()
    CACHE_DIR = notebook_dir / "data" / "raw_data" / "gtp_xref_audit"
    OUTPUT_DIR = notebook_dir / "data" / "processed_data" / "gtp_xref_audit"
    # a missing cache (the notebook on its own, as on molab) starts from the one in the GitHub repository
    copied_from_github = [] if refresh.value else seed_cache(CACHE_DIR, GITHUB_CACHE_URL, CACHE_FILES)
    mo.callout(
        mo.md(
            f"Started from the cached data in the [GitHub repository](https://github.com/ali-zaeim/GtoPdb_chemistry_audit) "
            f"({len(copied_from_github)} files), so the results are those published there. Tick **Re-download all "
            "source data** to fetch current data instead."
        ),
        kind="info",
    ) if copied_from_github else None
    return CACHE_DIR, OUTPUT_DIR


@app.cell
def _():
    refresh = mo.ui.checkbox(
        label="Re-download all source data (otherwise cached files and API responses are reused)"
    )
    refresh
    return (refresh,)


@app.cell
def _(CACHE_DIR, refresh):
    gtp_ligands_raw, gtp_version, gtp_ligands_log = download_gtp_file(
        "ligands.csv", CACHE_DIR, refresh.value
    )
    gtp_mapping_raw, gtp_mapping_version, _ = download_gtp_file(
        "ligand_id_mapping.csv", CACHE_DIR, refresh.value
    )
    assert gtp_version == gtp_mapping_version, (gtp_version, gtp_mapping_version)
    gtp_version_number = re.search(r"Version: ([\d.]+)", gtp_version).group(1)
    return (
        gtp_ligands_log,
        gtp_ligands_raw,
        gtp_mapping_raw,
        gtp_version,
        gtp_version_number,
    )


@app.cell
def _(gtp_ligands_raw):
    gtp_ligands = (
        gtp_ligands_raw.rename(
            columns={
                "Ligand ID": "ligand_id",
                "Name": "name_html",
                "Type": "type",
                "SMILES": "smiles",
                "InChI": "inchi",
                "InChIKey": "inchikey",
            }
        )
        .assign(
            ligand_id=lambda d: d.ligand_id.astype(int),
            name=lambda d: d.name_html.map(plain_text),
            radioactive=lambda d: d.Radioactive.eq("yes"),
            smiles=lambda d: d.smiles.mask(d.smiles.eq("")),
            inchi=lambda d: d.inchi.mask(d.inchi.eq("")),
            inchikey=lambda d: d.inchikey.mask(d.inchikey.eq("")),
        )[["ligand_id", "name", "type", "radioactive", "smiles", "inchi", "inchikey"]]
        .assign(
            # GtoPdb gives an InChIKey but no InChI for its largest ligands (section 6.3), so InChIs are rebuilt from the SMILES and kept only if they match the InChIKey
            needs_inchi=lambda d: d.inchi.isna() & d.inchikey.notna(),
            inchi=lambda d: [
                inchi_matching_key(smiles, key) if rebuild else inchi
                for inchi, smiles, key, rebuild in zip(
                    d.inchi, d.smiles, d.inchikey, d.needs_inchi
                )
            ],
        )
        .assign(inchi_rebuilt=lambda d: d.needs_inchi & d.inchi.notna())
        .drop(columns="needs_inchi")
    )
    return (gtp_ligands,)


@app.cell
def _(gtp_mapping_raw):
    xref_rows = []
    for record in gtp_mapping_raw.to_dict("records"):
        for db, column in XREF_COLUMNS.items():
            values = [value for value in record[column].split("|") if value.strip()]
            if db == "PubChem":
                # GtoPdb writes 0 as the PubChem CID of ligands with a PubChem substance (SID) but no compound (CID)
                values = [value for value in values if value != "0"]
            for raw_id in values:
                # what is wrong with how the ID is written is reported with its other issues, by `assess_xref`
                xref_rows.append(
                    {
                        "ligand_id": int(record["Ligand id"]),
                        "xref_db": db,
                        "xref_id_as_written": raw_id,
                        "xref_id_normalised": normalise_xref(db, raw_id),
                        "ids_in_field": len(values),
                    }
                )
    gtp_xrefs = pd.DataFrame(xref_rows)
    return (gtp_xrefs,)


@app.cell
def _(CACHE_DIR, refresh):
    unichem_pairs, unichem_mapping_logs = {}, {}
    for source_name in ("ChEMBL", "ChEBI", "PubChem"):
        mapping, unichem_mapping_logs[source_name] = download_unichem_mapping(
            "GtoPdb", source_name, CACHE_DIR, refresh.value
        )
        unichem_pairs[source_name] = mapping.rename(
            columns={"GtoPdb": "gtopdb_id", source_name: "xref_id"}
        )
    return unichem_mapping_logs, unichem_pairs


@app.cell
def _(CACHE_DIR, gtp_ligands, refresh):
    # the UCI and standard InChI UniChem holds for each GtoPdb ligand (section 7); its API takes one ID per request
    unichem_gtp_records = fetch_cached(
        CACHE_DIR / "unichem_gtopdb_compounds.parquet",
        sorted(gtp_ligands.ligand_id.astype(str)),
        fetch_unichem_gtopdb_batch,
        batch_size=25,
        workers=6,
        refresh=refresh.value,
        title="UniChem API",
    )
    return (unichem_gtp_records,)


@app.cell
def _(CACHE_DIR, gtp_xrefs, refresh):
    chembl_records = fetch_cached(
        CACHE_DIR / "chembl_molecules.parquet",
        ids_to_fetch("ChEMBL", gtp_xrefs),
        fetch_chembl_batch,
        batch_size=200,
        workers=4,
        refresh=refresh.value,
        title="ChEMBL API",
    )
    return (chembl_records,)


@app.cell
def _(CACHE_DIR, chembl_records, gtp_xrefs, refresh):
    # IDs the ChEMBL API does not return are checked against ChEMBL's ID lookup, which also knows retired IDs
    # and IDs of other ChEMBL entities (documents, targets, assays)
    unfound_chembl_ids = sorted(
        set(gtp_xrefs.loc[gtp_xrefs.xref_db == "ChEMBL", "xref_id_normalised"].dropna())
        & set(chembl_records.loc[~chembl_records.found, "query_id"])
    )
    chembl_id_lookup = fetch_cached(
        CACHE_DIR / "chembl_id_lookup.parquet",
        unfound_chembl_ids,
        fetch_chembl_id_lookup_batch,
        batch_size=200,
        workers=2,
        refresh=refresh.value,
        title="ChEMBL ID lookup",
    )
    return (chembl_id_lookup,)


@app.cell
def _(CACHE_DIR, gtp_xrefs, refresh):
    chebi_records = fetch_cached(
        CACHE_DIR / "chebi_compounds.parquet",
        ids_to_fetch("ChEBI", gtp_xrefs),
        fetch_chebi_batch,
        batch_size=200,
        workers=6,
        refresh=refresh.value,
        title="ChEBI API",
    )
    return (chebi_records,)


@app.cell
def _(CACHE_DIR, gtp_xrefs, refresh):
    pubchem_records = fetch_cached(
        CACHE_DIR / "pubchem_compounds.parquet",
        ids_to_fetch("PubChem", gtp_xrefs),
        fetch_pubchem_batch,
        batch_size=200,
        workers=2,
        refresh=refresh.value,
        title="PubChem PUG-REST",
    )
    return (pubchem_records,)


@app.cell
def _(CACHE_DIR, gtp_mapping_raw, refresh):
    # the substance GtoPdb deposited in PubChem for each ligand (PubChem SID), and the compound it is merged into
    gtp_pubchem_sids = (
        gtp_mapping_raw.assign(
            ligand_id=lambda d: d["Ligand id"].astype(int),
            sid=lambda d: d["PubChem SID"].str.split("|"),
        )
        .explode("sid")
        .assign(sid=lambda d: d.sid.str.strip())[["ligand_id", "sid"]]
    )
    gtp_pubchem_sids = gtp_pubchem_sids[gtp_pubchem_sids.sid.str.fullmatch(r"\d+")]
    gtp_sid_cids = fetch_cached(
        CACHE_DIR / "pubchem_sid_cids.parquet",
        sorted(set(gtp_pubchem_sids.sid)),
        fetch_pubchem_sid_cids,
        batch_size=500,
        workers=1,
        refresh=refresh.value,
        title="PubChem SID look-up",
    )
    return gtp_pubchem_sids, gtp_sid_cids


@app.cell
def _(
    chebi_records,
    chembl_id_lookup,
    chembl_records,
    chembl_suggested_records,
    gtp_ligands,
    gtp_ligands_log,
    gtp_sid_cids,
    inorganic_pubchem_lookup,
    pubchem_records,
    pubchem_standardisation,
    records_by_inchikey,
    unichem_gtp_records,
    unichem_mapping_logs,
    unichem_pairs,
):
    retrieval_summary = pd.DataFrame(
        [
            {
                "source": "GtoPdb ligands.csv + ligand_id_mapping.csv (download)",
                "rows": len(gtp_ligands),
                "IDs found": len(gtp_ligands),
                "with structure (InChI)": int(gtp_ligands.inchi.notna().sum()),
                "retrieved": gtp_ligands_log["downloaded_at_utc"][:10],
            },
            *[
                {
                    "source": f"UniChem GtoPdb–{name} mapping ({Path(unichem_mapping_logs[name]['downloaded_from']).name}, published {unichem_mapping_logs[name]['last_modified'][5:16]})",
                    "rows": len(pairs),
                    "IDs found": int(pairs.gtopdb_id.nunique()),
                    "with structure (InChI)": None,
                    "retrieved": unichem_mapping_logs[name]["downloaded_at_utc"][:10],
                }
                for name, pairs in unichem_pairs.items()
            ],
            {
                "source": "UniChem API (the structure UniChem holds for each GtoPdb ligand, section 7)",
                "rows": len(unichem_gtp_records),
                "IDs found": int(unichem_gtp_records.found.sum()),
                "with structure (InChI)": int(unichem_gtp_records.inchi.notna().sum()),
                "retrieved": retrieved_on(unichem_gtp_records),
            },
            *[
                {
                    "source": name,
                    "rows": len(records),
                    "IDs found": int(records.found.sum()),
                    "with structure (InChI)": int(records.inchi.notna().sum()),
                    "retrieved": retrieved_on(records),
                }
                for name, records in [
                    ("ChEMBL API (GtoPdb ChEMBL IDs)", chembl_records),
                    ("ChEBI API (GtoPdb ChEBI IDs)", chebi_records),
                    ("PubChem PUG-REST (GtoPdb CIDs)", pubchem_records),
                ]
            ],
            {
                "source": "ChEMBL ID lookup (GtoPdb ChEMBL IDs the molecule API does not return)",
                "rows": len(chembl_id_lookup),
                "IDs found": int(chembl_id_lookup.found.sum()),
                "with structure (InChI)": None,
                "retrieved": retrieved_on(chembl_id_lookup),
            },
            *[
                {
                    "source": f"{db} InChIKey look-up (records with the GtoPdb structure, sections 4 and 5)",
                    "rows": int(records.query_id.nunique()),
                    "IDs found": int(records.xref_id.notna().sum()),
                    "with structure (InChI)": int(records.inchi.notna().sum()),
                    "retrieved": retrieved_on(records),
                }
                for db, records in records_by_inchikey.groupby("xref_db")
            ],
            {
                "source": "PubChem SID look-up (the compound GtoPdb's own deposit is merged into, section 3)",
                "rows": len(gtp_sid_cids),
                "IDs found": int(gtp_sid_cids.found.sum()),
                "with structure (InChI)": None,
                "retrieved": retrieved_on(gtp_sid_cids),
            },
            {
                "source": "PubChem standardisation service (GtoPdb structures of PubChem links less specific than GtoPdb, section 5.5.3)",
                "rows": len(pubchem_standardisation),
                "IDs found": int(pubchem_standardisation.found.sum()),
                "with structure (InChI)": int(
                    pubchem_standardisation.standardised_inchi.notna().sum()
                ),
                "retrieved": retrieved_on(pubchem_standardisation),
            },
            {
                "source": "PubChem InChIKey look-up (structures of inorganic ligands, section 6.5)",
                "rows": int(inorganic_pubchem_lookup.query_id.nunique()),
                "IDs found": int(inorganic_pubchem_lookup.xref_id.notna().sum()),
                "with structure (InChI)": int(inorganic_pubchem_lookup.inchi.notna().sum()),
                "retrieved": retrieved_on(inorganic_pubchem_lookup),
            },
            {
                "source": "ChEMBL API (ChEMBL records named in suggested actions: IDs with the GtoPdb structure and parent compounds)",
                "rows": len(chembl_suggested_records),
                "IDs found": int(chembl_suggested_records.found.sum()),
                "with structure (InChI)": int(
                    chembl_suggested_records.get("inchi", pd.Series(dtype=object)).notna().sum()
                ),
                "retrieved": retrieved_on(chembl_suggested_records),
            },
        ]
    ).astype({"with structure (InChI)": "Int64"})
    retrieval_summary
    return (retrieval_summary,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    *Rows* are the records in each download, or the IDs sent to each API. The UniChem mapping files have one row per pair of IDs that share a UniChem compound, so a GtoPdb ligand appears once for every ChEBI ID or PubChem CID it shares a structure with. For these files, *IDs found* counts distinct GtoPdb ligands. A ligand has several partners when ChEBI keeps as separate entries tautomers that standard InChI merges (it merges some tautomers, not all: [Dhaked et al.](https://doi.org/10.1021/acs.jcim.9b01080)), or when PubChem holds the same compound under several CIDs (vendor products, grades, E/Z names). **A few PubChem groups mix stereoisomers, because UniChem holds both the GtoPdb ligand and some of its CIDs with their stereo removed (section 7). ChEMBL registers each structure once, so its mapping has one row per ligand.**
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(rf"""
    ## 2. How structures were compared

    A standard InChI writes a structure as a series of layers, always in the same order, each adding detail to the ones before: formula (`/C9H8O4`), connections (`/c`), hydrogens (`/h`), charge (`/q`), protons (`/p`), stereo (`/b`, `/t`, `/m`, `/s`) and isotopes (`/i`). Each cross-reference is compared by reading the two InChIs from the left and stopping at the first layer that differs, since the later layers describe a structure the two records no longer share. Section 5 is organised by that layer. Standard InChI removes the bonds to metals, so a cross-reference where either structure contains a metal is a category of its own whenever the InChIs are not identical (section 5.7).

    {mo.image(
        "https://www.inchi-trust.org/wp/wp-content/uploads/2022/12/InChI-KEY-2048x728.png",
        alt="A standard InChI with each layer labelled, and its InChIKey with each block labelled",
        width="100%",
        style={"max-width": "760px", "align-self": "center"},
        caption='The layers of a standard InChI and the blocks of its InChIKey (image: <a href="https://www.inchi-trust.org/">InChI Trust</a>)',
    )}

    The GtoPdb InChI is the `InChI` column of `ligands.csv`, or the InChI rebuilt from the SMILES where that column is empty (section 6.3). The linked record's InChI is the standard InChI that its own database provides.

    Three terms recur below:

    - **Linked record**: the ChEMBL, ChEBI or PubChem entry that a GtoPdb cross-reference points to.
    - **Exact record**: a record in the linked database whose standard InChI is exactly GtoPdb's, found by searching the GtoPdb InChIKey (`xref_id_with_gtp_structure`). It is the record to link if the GtoPdb structure is right.
    - **Same parent after standardisation**: the two records have the same largest organic fragment once it is neutralised and its isotope labels removed, either directly or as RDKit canonical tautomers. *Apart from stereo* means that the two parents differ only in stereo.

    Each cross-reference gets a suggested action (`suggested_action`). It has one step for each identifier issue (section 4), then one for the structure, from the first of these rules that applies:

    1. GtoPdb already links an exact record as well: *remove link*.
    2. An exact record exists: *replace link or change GtoPdb structure*. Either the exact record is the right link, or the current link is right and the GtoPdb structure should change.
    3. The linked record has the structure of another GtoPdb ligand (`other_gtp_ligands_with_xref_structure`): *remove link*, since the identifier probably belongs to that ligand.
    4. Two or more external databases agree with each other but not with GtoPdb: *check GtoPdb structure*.
    5. Otherwise, the first layer that differs sets the step, or a metal in either structure does:

    | Section | First layer that differs | Step |
    |---|---|---|
    | 5.1, 5.2 | Formula or connections | *check which structure is intended* |
    | 5.3.1 | Hydrogens, same parent after standardisation | *keep link* |
    | 5.3.2 | Hydrogens, different parents | *check which structure is intended* |
    | 5.4 | Charge | *check which structure is intended* |
    | 5.4 | Protons | *keep link* |
    | 5.5.1 | Stereo, opposite configuration | *check GtoPdb structure* |
    | 5.5.2 | Stereo defined differently on each side | *review* |
    | 5.5.3 | Stereo, GtoPdb defines more | *keep link* |
    | 5.5.3 | Stereo, the linked record defines more | *check GtoPdb structure* |
    | 5.6 | Isotopes, only GtoPdb labelled | *keep link* |
    | 5.6 | Isotopes, only the linked record labelled | *check GtoPdb structure* for a radioactive ligand, otherwise *check which structure is intended* |
    | 5.7 | Any layer, with a metal in either structure | *check which structure is intended* |
    | 5.7 | None: the linked record of a metal-containing ligand has no structure | *review* |

    A kept link always says why: how it differs, and whether the lookup found an exact record. Notes add context without changing the action. They say when:

    - an exact record is suggested, but the current link has another GtoPdb ligand's structure, or the external databases agree with each other but not with GtoPdb;
    - the only record that agrees with GtoPdb is the PubChem compound of GtoPdb's own deposit (section 3);
    - a radioactive ligand drawn without its isotope could not be searched for an exact record (section 5.6);
    - the linked record of a labelled ligand is unlabelled.

    Section 8 describes the columns of the output tables.
    """)
    return


@app.cell
def _(chebi_records, chembl_records, pubchem_records):
    external_records = pd.concat(
        [
            chembl_records.assign(xref_db="ChEMBL"),
            chebi_records.assign(xref_db="ChEBI"),
            pubchem_records.assign(xref_db="PubChem"),
        ],
        ignore_index=True,
    ).rename(
        columns={
            "query_id": "xref_id",
            "found": "xref_found",
            "inchi": "xref_inchi",
            "inchikey": "xref_inchikey",
        }
    )
    # ChEBI returns an empty InChI for entries without a defined structure (classes, R groups): no structure
    external_records = external_records.assign(
        xref_inchi=lambda d: d.xref_inchi.mask(d.xref_inchi.eq("")),
        xref_inchikey=lambda d: d.xref_inchikey.mask(d.xref_inchikey.eq("")),
    )
    return (external_records,)


@app.cell
def _(chembl_id_lookup, external_records, gtp_ligands, gtp_xrefs):
    ligands_by_inchikey = (
        gtp_ligands.dropna(subset=["inchikey"]).groupby("inchikey").ligand_id.agg(list)
    )

    xref_table = gtp_xrefs.merge(
        gtp_ligands.rename(
            columns={
                "inchi": "gtp_inchi",
                "inchikey": "gtp_inchikey",
                "inchi_rebuilt": "gtp_inchi_rebuilt",
            }
        ),
        on="ligand_id",
        how="left",
    ).merge(
        external_records.drop(columns="retrieved_at").rename(
            columns={"xref_id": "xref_id_normalised"}
        ),
        on=["xref_db", "xref_id_normalised"],
        how="left",
    ).merge(
        # why ChEMBL's molecule API does not return an ID: never issued, another kind of record, or retired
        chembl_id_lookup.drop(columns="retrieved_at")
        .rename(
            columns={
                "query_id": "xref_id_normalised",
                "found": "chembl_id_issued",
                "entity_type": "chembl_entity_type",
                "status": "chembl_status",
                "last_active": "chembl_last_active",
            }
        )
        .assign(xref_db="ChEMBL"),
        on=["xref_db", "xref_id_normalised"],
        how="left",
    )
    xref_table["xref_found"] = xref_table.xref_found.fillna(False).astype(bool)
    assessed = pd.DataFrame(
        [assess_xref(row) for row in xref_table.to_dict("records")], index=xref_table.index
    )
    comparison = pd.concat([xref_table, assessed], axis=1).assign(
        other_gtp_ligands_with_xref_structure=lambda d: [
            ", ".join(str(i) for i in ligands_by_inchikey.get(key, []) if i != own)
            if isinstance(key, str)
            else ""
            for key, own in zip(d.xref_inchikey, d.ligand_id)
        ],
    )
    return (comparison,)


@app.cell
def _(comparison):
    # cross-references for which the linked database is searched for a record with exactly the GtoPdb structure:
    # every structure that differs, except a linked record that defines more stereo than GtoPdb, and ChEMBL IDs that
    # do not resolve to a compound (section 4)
    lookup_rows = comparison[
        (
            (comparison.first_difference.isin(LAYERS) & ~comparison.more_stereo_in.eq("linked record"))
            | ((comparison.xref_db == "ChEMBL") & ~comparison.xref_found)
        )
        # radioactive ligands stored without their isotope are left out: their InChIKey only finds the unlabelled compound
        & ~(comparison.radioactive & ~comparison.gtp_inchi.fillna("").str.contains("/i"))
    ].dropna(subset=["gtp_inchikey"])
    return (lookup_rows,)


@app.cell
def _(CACHE_DIR, lookup_rows, refresh):
    records_by_inchikey = pd.concat(
        [
            fetch_cached(
                CACHE_DIR / f"{db.lower()}_by_inchikey.parquet",
                sorted(set(lookup_rows.loc[lookup_rows.xref_db == db, "gtp_inchikey"])),
                fetch,
                batch_size=10,
                workers=1,
                refresh=refresh.value,
                title=f"{db} InChIKey look-up",
            ).assign(xref_db=db)
            for db, fetch in (
                ("ChEMBL", fetch_chembl_by_inchikey),
                ("ChEBI", fetch_chebi_by_inchikey),
                ("PubChem", fetch_pubchem_by_inchikey),
            )
        ],
        ignore_index=True,
    ).reindex(
        columns=["xref_db", "query_id", "found", "xref_id", "inchi", "xref_name", "retrieved_at"]
    )
    return (records_by_inchikey,)


@app.cell
def _(lookup_rows, records_by_inchikey):
    # IDs in the linked database whose standard InChI is exactly GtoPdb's, stereo included
    hits_by_inchikey = {
        key: list(zip(group.xref_id, group.inchi))
        for key, group in records_by_inchikey.dropna(subset=["xref_id"]).groupby(
            ["xref_db", "query_id"]
        )
    }
    xref_ids_with_gtp_structure = pd.Series(
        {
            index: "|".join(matches)
            for index, db, inchikey, gtp_inchi, own_id in zip(
                lookup_rows.index,
                lookup_rows.xref_db,
                lookup_rows.gtp_inchikey,
                lookup_rows.gtp_inchi,
                lookup_rows.xref_id_normalised,
            )
            if (
                matches := [
                    xref_id
                    for xref_id, inchi in hits_by_inchikey.get((db, inchikey), [])
                    if inchi == gtp_inchi and xref_id != own_id
                ]
            )
        },
        dtype=object,
    )
    return (xref_ids_with_gtp_structure,)


@app.cell
def _(CACHE_DIR, comparison, refresh, xref_ids_with_gtp_structure):
    # PubChem links less specific than GtoPdb where PubChem has no CID with the GtoPdb structure: PubChem's
    # standardisation is run on the GtoPdb structure, to measure which of its stereo PubChem removes
    standardise_rows = comparison[comparison.more_stereo_in.eq("GtoPdb")]
    standardise_rows = standardise_rows[
        (standardise_rows.xref_db == "PubChem")
        & ~standardise_rows.index.isin(xref_ids_with_gtp_structure.index)
    ]
    pubchem_standardisation = fetch_cached(
        CACHE_DIR / "pubchem_standardised.parquet",
        sorted(set(standardise_rows.gtp_inchi)),
        fetch_pubchem_standardised,
        batch_size=10,
        workers=1,
        refresh=refresh.value,
        title="PubChem standardisation of GtoPdb structures",
    ).reindex(columns=["query_id", "found", "standardised_inchi", "retrieved_at"])
    pubchem_standardised = pubchem_standardisation.set_index("query_id").standardised_inchi
    standardised_inchi = standardise_rows.gtp_inchi.map(pubchem_standardised)
    # Attribute the linked record's missing stereo to standardisation only when the service
    # reproduces that record's full InChI, including the centres it still defines.
    reproduces_link = standardised_inchi.eq(standardise_rows.xref_inchi)
    pubchem_standardisation_removes = standardise_rows.loc[reproduces_link, "missing_stereo"]
    return pubchem_standardisation, pubchem_standardisation_removes


@app.cell
def _(
    CACHE_DIR,
    comparison,
    records_by_inchikey,
    refresh,
    xref_ids_with_gtp_structure,
):
    # names of the records that suggested actions point to: IDs with the GtoPdb structure and ChEMBL parent compounds.
    # ChEBI and PubChem names come with the InChIKey look-up; ChEMBL records are fetched here, with the parents
    suggested_ids = (
        pd.DataFrame(
            {
                "xref_db": comparison.loc[xref_ids_with_gtp_structure.index, "xref_db"],
                "xref_id": xref_ids_with_gtp_structure.str.split("|"),
            }
        )
        .explode("xref_id")
        .groupby("xref_db")
        .xref_id.agg(set)
        .to_dict()
    )
    chembl_parents = comparison.chembl_parent_id.where(
        comparison.chembl_parent_id.ne(comparison.xref_id_normalised)
    ).dropna()
    chembl_suggested_records = fetch_cached(
        CACHE_DIR / "chembl_suggested_records.parquet",
        sorted(suggested_ids.get("ChEMBL", set()) | set(chembl_parents)),
        fetch_chembl_batch,
        batch_size=100,
        workers=1,
        refresh=refresh.value,
        title="ChEMBL suggested records",
    ).assign(xref_db="ChEMBL")
    record_names = {
        (db, xref_id): name
        for db, xref_id, name in [
            *zip(
                records_by_inchikey.xref_db,
                records_by_inchikey.xref_id,
                records_by_inchikey.xref_name,
            ),
            *zip(
                chembl_suggested_records.xref_db,
                chembl_suggested_records.query_id,
                chembl_suggested_records.get("xref_name", []),
            ),
        ]
        if isinstance(name, str)
    }
    xref_names_with_gtp_structure = pd.Series(
        {
            index: "|".join(
                record_names.get((comparison.at[index, "xref_db"], xref_id), "")
                for xref_id in ids.split("|")
            )
            for index, ids in xref_ids_with_gtp_structure.items()
        },
        dtype=object,
    )
    chembl_parent_names = chembl_parents.map(lambda parent: record_names.get(("ChEMBL", parent)))
    return (
        chembl_parent_names,
        chembl_suggested_records,
        xref_names_with_gtp_structure,
    )


@app.cell
def _(comparison, gtp_pubchem_sids, gtp_sid_cids):
    # Identify CIDs associated with the ligand’s own SID independently of the InChI comparison.
    gtp_deposit_cids = (
        gtp_pubchem_sids.merge(gtp_sid_cids, left_on="sid", right_on="query_id")
        .assign(cid=lambda d: d.cids.str.split("|"))
        .explode("cid")
        .groupby("ligand_id")
        .cid.agg(set)
    )
    pubchem_cid_from_gtp_deposit = pd.Series(
        [
            db == "PubChem" and xref_id in gtp_deposit_cids.get(ligand_id, set())
            for ligand_id, db, xref_id in zip(
                comparison.ligand_id, comparison.xref_db, comparison.xref_id_normalised
            )
        ],
        index=comparison.index,
    )
    agreeing = comparison.first_difference == "identical"
    # ligands whose structure only that PubChem compound agrees with
    only_gtp_deposit_agrees = comparison.ligand_id.isin(
        set(comparison.loc[agreeing & pubchem_cid_from_gtp_deposit, "ligand_id"])
        - set(comparison.loc[agreeing & ~pubchem_cid_from_gtp_deposit, "ligand_id"])
    )
    return only_gtp_deposit_agrees, pubchem_cid_from_gtp_deposit


@app.cell
def _(
    chembl_parent_names,
    comparison,
    only_gtp_deposit_agrees,
    pubchem_cid_from_gtp_deposit,
    pubchem_standardisation_removes,
    xref_ids_with_gtp_structure,
    xref_names_with_gtp_structure,
):
    structure_compared = comparison[comparison.first_difference != "not compared"]
    consensus = structure_compared.groupby("ligand_id").agg(
        n_databases=("xref_db", "nunique"),
        external_structures=("xref_inchi", lambda inchis: sorted(set(inchis))),
        gtp_inchi=("gtp_inchi", "first"),
    )
    consensus_against_gtp = consensus[
        (consensus.n_databases >= 2)
        & (consensus.external_structures.map(len) == 1)
        & (consensus.external_structures.str[0] != consensus.gtp_inchi)
    ].index
    # IDs GtoPdb links to each ligand, to tell a wrong extra link from a link to replace
    ids_linked_by_gtp = (
        comparison.dropna(subset=["xref_id_normalised"])
        .groupby(["ligand_id", "xref_db"])
        .xref_id_normalised.agg(set)
    )
    curation = (
        comparison.assign(
            xref_id_with_gtp_structure=xref_ids_with_gtp_structure,
            xref_name_with_gtp_structure=xref_names_with_gtp_structure,
            chembl_parent_name=chembl_parent_names,
            pubchem_standardisation_removes=pubchem_standardisation_removes,
            pubchem_cid_from_gtp_deposit=pubchem_cid_from_gtp_deposit,
            only_gtp_deposit_agrees=only_gtp_deposit_agrees,
            external_consensus_differs=comparison.ligand_id.isin(consensus_against_gtp),
        )
        .assign(
            exact_record_already_linked=lambda d: [
                isinstance(ids, str)
                and bool(set(ids.split("|")) & ids_linked_by_gtp.get((ligand, db), set()))
                for ids, ligand, db in zip(d.xref_id_with_gtp_structure, d.ligand_id, d.xref_db)
            ],
        )
    )
    # `action`, `suggested_xref_id` and `suggested_action`
    curation = curation.join(
        pd.DataFrame([suggested_action(row) for row in curation.to_dict("records")], index=curation.index)
    )
    return consensus_against_gtp, curation


@app.cell(hide_code=True)
def _(curation):
    pubchem_agreeing = curation[
        (curation.xref_db == "PubChem") & (curation.first_difference == "identical")
    ]
    mo.md(f"""
    ## 3. How every cross-reference compares

    The first InChI layer that differs from the GtoPdb structure, for every cross-reference, by linked database. *Metal-containing* counts comparable cross-references in section 5.7 whose InChIs differ. *Not compared* means that the GtoPdb ligand or the linked record has no structure (section 6.5 gives one for most inorganic ligands), or that the identifier could not be resolved (section 4).

    **Agreement with PubChem is weaker evidence than it looks.** GtoPdb deposits its own ligands in PubChem, and PubChem groups substance records (SIDs) whose structures standardise successfully under the corresponding compound (CID). For {int(pubchem_agreeing.pubchem_cid_from_gtp_deposit.sum()):,} of the {len(pubchem_agreeing):,} PubChem cross-references that agree with GtoPdb, the CID is the one that GtoPdb's own deposit (`PubChem SID` in GtoPdb's download) is merged into. A matching CID alone cannot establish independent corroboration. Other depositors may also contribute to the CID. Where such a CID is the only record that agrees with GtoPdb and another database disagrees, the suggested action says so.
    """)
    return


@app.cell
def _(curation):
    # metal-containing structures that differ are a category of their own (section 5.7), whatever the layer
    first_difference_table = (
        pd.crosstab(
            curation.first_difference.mask(curation.metal_differs & curation.first_difference.isin(LAYERS), "metal-containing"),
            curation.xref_db,
            margins=True,
            margins_name="Total",
        )
        .reindex(["identical", *LAYERS, "metal-containing", "not compared", "Total"], fill_value=0)
        .rename_axis(index="first difference", columns=None)
    )
    first_difference_table
    return


@app.cell
def _(curation):
    _differ = curation[curation.first_difference.isin(LAYERS)]
    issue_counts = (
        # metal-containing structures as one category (section 5.7), whatever the layer
        _differ.assign(
            first_difference=_differ.first_difference.mask(_differ.metal_differs, "metal-containing"),
            stereo_relation=_differ.stereo_relation.mask(_differ.metal_differs, ""),
        )
        .groupby(["first_difference", "stereo_relation", "xref_db"])
        .size()
        .rename("count")
        .reset_index()
    )
    return (issue_counts,)


@app.cell
def _(issue_counts):
    # in InChI order from the top, so reversed: Plotly draws the first horizontal category at the bottom
    layer_order = [
        layer for layer in reversed([*LAYERS, "metal-containing"]) if layer in set(issue_counts.first_difference)
    ]
    surface, ink, muted, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9"
    discrepancy_chart = go.Figure()
    for linked_db, colour in (("ChEMBL", "#2a78d6"), ("ChEBI", "#eb6834"), ("PubChem", "#1baf7a")):
        db_counts = (
            issue_counts[issue_counts.xref_db == linked_db].groupby("first_difference")["count"].sum()
        )
        discrepancy_chart.add_bar(
            name=linked_db,
            x=db_counts.values,
            y=db_counts.index,
            orientation="h",
            marker={"color": colour, "line": {"color": surface, "width": 2}},
            hovertemplate=f"<b>%{{y}}</b><br>{linked_db}: %{{x}} cross-references<extra></extra>",
        )
    discrepancy_chart.update_layout(
        barmode="stack",
        bargap=0.45,
        title={
            "text": "Where cross-referenced structures first differ from GtoPdb",
            "subtitle": {
                "text": "Cross-references whose InChI differs from the GtoPdb ligand, by the first InChI layer that differs, "
                "with metal-containing structures apart",
                "font": {"color": muted, "size": 12},
            },
            "x": 0,
            "xanchor": "left",
            "font": {"color": ink, "size": 15},
        },
        legend={
            "orientation": "h",
            "x": 0,
            "y": 1.0,
            "yanchor": "bottom",
            "title": {"text": "Linked database"},
            "traceorder": "normal",
        },
        font={
            "family": "system-ui, -apple-system, Segoe UI, sans-serif",
            "color": muted,
            "size": 12,
        },
        paper_bgcolor=surface,
        plot_bgcolor=surface,
        height=140 + 30 * len(layer_order),
        margin={"l": 10, "r": 20, "t": 110, "b": 50},
        xaxis={
            "title": "Cross-references",
            "gridcolor": grid,
            "zeroline": False,
            "showline": False,
        },
        yaxis={
            "title": "First layer that differs",
            "categoryorder": "array",
            "categoryarray": layer_order,
            "tickfont": {"color": ink},
            "ticks": "",
            "ticksuffix": "  ",
            "showgrid": False,
        },
        hoverlabel={"bgcolor": "white", "font": {"color": ink}},
    )
    discrepancy_chart
    return


@app.cell
def _(issue_counts):
    # each stereo relation on its own row, in InChI order with the total last
    issue_table = issue_counts.pivot_table(
        index=["first_difference", "stereo_relation"],
        columns="xref_db",
        values="count",
        fill_value=0,
        aggfunc="sum",
        margins=True,
        margins_name="Total",
    )
    _order = [*LAYERS, "metal-containing", "Total"]
    issue_table = issue_table.reindex(
        sorted(issue_table.index, key=lambda key: (_order.index(key[0]), key[1]))
    ).rename_axis(index=["first difference", "stereo relation"], columns=None)
    issue_table
    return


@app.cell(hide_code=True)
def _(curation, curation_issues, inorganic_links):
    # The curation issues retain their original row indices. Section 6.5 lists
    # every cross-reference of each ligand in inorganic_links.
    _remaining = curation[
        curation.first_difference.eq("not compared")
        & ~curation.index.isin(curation_issues.index)
        & ~curation.ligand_id.isin(inorganic_links.ligand_id)
    ].copy()
    _remaining["missing_inchi"] = "Linked record only"
    _remaining.loc[_remaining.gtp_inchi.isna(), "missing_inchi"] = "GtoPdb only"
    _remaining.loc[_remaining.gtp_inchi.isna() & _remaining.xref_inchi.isna(), "missing_inchi"] = "Both records"
    _counts = _remaining.missing_inchi.value_counts()
    mo.vstack([
        mo.md(f"""
        ### 3.1 Cross-references missing a standard InChI

        These **{len(_remaining):,} cross-references** on **{_remaining.ligand_id.nunique():,} ligands** are the remaining links
        for which a structure comparison was impossible. The standard InChI is missing from GtoPdb alone for
        **{_counts.get("GtoPdb only", 0):,}** links, from the linked record alone for **{_counts.get("Linked record only", 0):,}**,
        and from both records for **{_counts.get("Both records", 0):,}**.
        """),
        linked_table(
            _remaining.sort_values(["missing_inchi", "xref_db", "ligand_id"])[
                ["ligand_id", "name", "type", "xref_db", "xref_id_as_written", "xref_name", "missing_inchi"]
            ].reset_index(drop=True)
        ),
    ])
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## 4. Identifier issues

    This section lists identifiers that:

    - cannot be used as written (hidden characters, non-standard case, wrong database), or share a field with another ID;
    - no longer exist, or are secondary ChEBI IDs;
    - point to a ChEMBL salt or child form whose structure differs from GtoPdb's, where the parent compound may be the intended link.

    A ChEMBL child form with exactly the GtoPdb structure, such as a salt or isotopologue that GtoPdb lists as its own ligand, is the right link and is not listed. `difference` shows how the linked structure compares with GtoPdb's; a cross-reference whose structure also differs is listed in section 5 as well.

    A ChEMBL ID that the molecule API does not return is checked in ChEMBL's ID lookup. This tells a retired compound (with the last ChEMBL release that contained it) from another kind of ChEMBL record, such as a document, and from an ID that was never issued. ChEBI's public data cannot tell a deleted ID from one that was never issued. For ChEMBL IDs that do not resolve to a compound, the suggested action names the exact record where there is one.

    In the tables of sections 4 and 5, `structures` draws the GtoPdb structure above the linked record's, each captioned with its ID, linked to the record, and its name. `inchis` gives both InChIs, GtoPdb's first, with the layer that differs first highlighted. The IDs in `suggested_action` are links too. The columns with the names and IDs, the other records the suggested action names, and each InChI on its own are hidden: the table's search still finds them, and its column panel shows them again.
    """)
    return


@app.cell
def _(curation):
    identifier_problems = (
        curation[curation.identifier_issue != ""]
        .sort_values(["xref_db", "identifier_issue", "ligand_id"])
        .assign(
            xref_id_as_written=lambda d: d.xref_id_as_written.map(
                lambda value: value if value.isascii() else ascii(value)
            ),
            # the parent compound, only where linking it is suggested
            chembl_parent_id=lambda d: d.chembl_parent_id.where(
                d.identifier_issue_kinds.str.split("|").map(
                    lambda kinds: IdentifierIssue.CHEMBL_CHILD_FORM in kinds
                )
            ),
        )[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "identifier_issue",
                "difference",
                "action",
                "suggested_action",
                "type",
                "chebi_stars",
                "chembl_parent_id",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(identifier_problems)
    return


@app.cell(hide_code=True)
def _(curation):
    chebi_link_rows = curation[curation.xref_db == "ChEBI"]
    chebi_other_compound_rows = chebi_link_rows[
        chebi_link_rows.first_difference.isin(["formula", "connections"]) & ~chebi_link_rows.metal_differs
    ]
    mo.md(f"""
    ## 5. Structure differences, layer by layer

    Every cross-reference whose structure differs from GtoPdb's, by the first InChI layer that differs (section 2), except those with a metal in either structure, which are listed together in 5.7. A cross-reference that also has an identifier issue is listed in section 4 too. Where there is an exact record, the suggested action names it.

    `chebi_stars` is ChEBI's curation level: 3 means annotated by the ChEBI curation team, 2 means provided by a third party and not checked by ChEBI. {int(chebi_link_rows.chebi_stars.eq(2).sum())} of GtoPdb's {len(chebi_link_rows):,} ChEBI links ({chebi_link_rows.chebi_stars.eq(2).mean():.0%}) point to 2-star entries, but {int(chebi_other_compound_rows.chebi_stars.eq(2).sum())} of the {len(chebi_other_compound_rows)} that differ in formula or connections ({chebi_other_compound_rows.chebi_stars.eq(2).mean():.0%}) do.

    ### 5.1 Differs in formula

    `difference_detail` gives both formulas, then what the records still share: the number of components where it differs (a salt or mixture), *same atoms and connections* where only the hydrogen count differs, and *same parent after standardisation*. A hydrogen count differs when a double bond is left out or a hydrogen is missing (drawn as a radical): VU0010010 (ligand 7633), for example, is drawn in GtoPdb without the C2=C3 double bond of its thieno[2,3-b]pyridine. `tanimoto_similarity` is the fingerprint similarity of the two structures.

    Where the linked record has the structure of another GtoPdb ligand, the suggested action names that ligand: the identifier probably belongs to it. Without an exact record, the suggested action removes the link.
    """)
    return


@app.cell
def _(curation):
    # the cross-references listed by layer in sections 5.1 to 5.6; metal-containing structures are listed in 5.7
    layer_rows = curation[~curation.metal_differs]
    return (layer_rows,)


@app.cell
def _(layer_rows):
    formula_differences = (
        layer_rows[layer_rows.first_difference.eq("formula")]
        .sort_values(["xref_db", "tanimoto_similarity", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "action",
                "suggested_action",
                "tanimoto_similarity",
                "xref_id_with_gtp_structure",
                "chebi_stars",
                "other_gtp_ligands_with_xref_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(formula_differences)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 5.2 Differs in connections

    The records have the same formula, but their atoms are connected differently, for example a substituent at another position. `difference_detail` gives the part of each structure outside their largest common substructure, as SMILES with `*` where it bonds to the shared part. In `structures`, those atoms are highlighted in orange and the bonds joining them to the shared part in red. Structures above 60 heavy atoms, such as peptides, are too large for this search and are marked as not isolated automatically.
    """)
    return


@app.cell
def _(layer_rows):
    connection_differences = (
        layer_rows[layer_rows.first_difference.eq("connections")]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "action",
                "suggested_action",
                "tanimoto_similarity",
                "xref_id_with_gtp_structure",
                "chebi_stars",
                "other_gtp_ligands_with_xref_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(connection_differences)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 5.3 Differs in hydrogens

    The records have the same atoms connected the same way, but place hydrogens, and with them double bonds, differently. The atoms whose hydrogen count differs are highlighted. `difference_detail` says whether the records have the same parent after standardisation (section 2).

    #### 5.3.1 Same parent after standardisation

    Where the parents are the same only as RDKit canonical tautomers, the records are tautomers that standard InChI does not merge. Without an exact record, the link is kept.
    """)
    return


@app.cell
def _(layer_rows):
    hydrogen_tautomers = (
        layer_rows[layer_rows.first_difference.eq("hydrogens") & layer_rows.difference_detail.map(same_parent)]
        .sort_values(["difference_detail", "xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "action",
                "suggested_action",
                "xref_id_with_gtp_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(hydrogen_tautomers)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    #### 5.3.2 Different parents

    The parents differ even as RDKit canonical tautomers. The records may be tautomers that neither standard InChI nor RDKit recognises, such as an N-oxide drawn as an N-hydroxy imine, or different isomers, such as one with a double bond moved. The structures alone cannot tell which, so the curator checks which structure is intended.
    """)
    return


@app.cell
def _(layer_rows):
    hydrogen_other_isomers = (
        layer_rows[layer_rows.first_difference.eq("hydrogens") & ~layer_rows.difference_detail.map(same_parent)]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "action",
                "suggested_action",
                "xref_id_with_gtp_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(hydrogen_other_isomers)
    return


@app.cell(hide_code=True)
def _(charge_differences, gtp_ligands):
    _keyed = gtp_ligands.dropna(subset=["inchikey"])
    _flags = _keyed.assign(core=_keyed.inchikey.str[:25], flag=_keyed.inchikey.str[-1])
    _two_states = _flags.groupby("core").flag.nunique().gt(1).sum()
    _protons_changed = int(_flags.flag.ne("N").sum())
    gtp_ion_rows = charge_differences[
        charge_differences.gtp_inchi.map(lambda inchi: net_charge(parse_standard_inchi(inchi))).ne(0)
        & charge_differences.xref_inchi.map(lambda inchi: net_charge(parse_standard_inchi(inchi))).eq(0)
    ]
    _ions_by_db = ", ".join(f"{db} {n}" for db, n in gtp_ion_rows.xref_db.value_counts().items())
    _anions = int(
        gtp_ion_rows.gtp_inchi.map(lambda inchi: int(parse_standard_inchi(inchi).get("p") or 0) < 0).sum()
    )
    mo.md(f"""
    ### 5.4 Differs in charge or protons

    The records carry different charges (`/q`) or have protons added or removed (`/p`): an acid against its conjugate base, or an ion against its neutral form. `difference_detail` gives the net charge of each record and how many more or fewer protons GtoPdb has. Without an exact record, a link that differs in protons is kept, and for one that differs in charge the curator checks which structure is intended.

    Before linking an exact record instead, check that the GtoPdb structure is the intended form. ChEBI keeps an acid and its conjugate base as separate entries on purpose, and GtoPdb's acetoacetic acid (ligand 4522), for example, is drawn as the acetoacetate anion.

    GtoPdb mostly stores one form per compound: {_protons_changed:,} of its {len(_flags):,} ligands with an InChIKey are drawn with protons added or removed (InChIKey protonation flag other than N), and only {_two_states} compounds are stored in two protonation states. In {len(gtp_ion_rows)} of the {len(charge_differences)} cross-references below, GtoPdb draws an ion and the linked record is neutral ({_ions_by_db}). {_anions} of these GtoPdb ions are drawn with protons removed and no counterion, like an acid drawn as its anion.
    """)
    return


@app.cell
def _(layer_rows):
    charge_differences = (
        layer_rows[layer_rows.first_difference.isin(["charge", "protons"])]
        .sort_values(["first_difference", "xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference",
                "difference_detail",
                "action",
                "suggested_action",
                "xref_id_with_gtp_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(charge_differences)
    return (charge_differences,)


@app.cell(hide_code=True)
def _():
    mo.vstack(
        [
            mo.md(r"""
    ### 5.5 Differs in stereo

    The records disagree about stereochemistry in one of three ways (`stereo_relation`):

    - **Opposite configuration (5.5.1):** both define the same stereocentres and double bonds, and at least one is inverted.
    - **Defined differently on each side (5.5.2):** each defines a centre that the other leaves undefined, or one defines more and a centre both define is written inverted.
    - **Defined on one side only (5.5.3):** one record defines every centre that the other does, and more, all written the same way.

    Configurations are compared only in the first case, because next to an undefined centre InChI can write the same configuration inverted. Stereo created by isotope labels, such as a CHT group made chiral by tritium, is compared here when both records carry the same labels.
    """),
            mo.accordion(
                {
                    "Why InChI can write the same configuration inverted": mo.md(r"""
    InChI writes a configuration relative to its numbering of the neighbouring atoms. For symmetry-equivalent atoms, such as the bridgeheads of a tropane or the two ends of a ring cis/trans pair in a 1,4-disubstituted cyclohexane, that numbering depends on the other centres. Next to an undefined centre, the same configuration can therefore be written inverted. Ricasetron (ligand 2302) and its PubChem record look like opposite configurations at both bridgeheads, yet defining the record's endo/exo centre gives exactly the GtoPdb InChI.

    Where both records define the same centres, the descriptors are compared centre by centre after applying the InChI `/m` inversion flag.
    """)
                }
            ),
            mo.md(r"""
    #### 5.5.1 Opposite configuration

    The records define the same stereocentres and double bonds, and at least one has the opposite configuration: *enantiomer*, *epimer or diastereomer*, *E/Z isomer*, or *enantiomer or diastereomer* where both leave the same centres undefined (an enantiomer if those are inverted too). Either the cross-reference points at the wrong stereoisomer, or the GtoPdb stereochemistry is wrong. Agreement among the comparable external records flags the GtoPdb structure and its links for review. Establishing the intended stereoisomer requires primary evidence.
    """),
        ]
    )
    return


@app.cell
def _(layer_rows):
    stereo_conflicts = (
        layer_rows[layer_rows.stereo_relation.isin(STEREO_CONFLICT_RELATIONS)]
        .sort_values(["stereo_relation", "ligand_id", "xref_db"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "xref_id_with_gtp_structure",
                "stereo_relation",
                "action",
                "suggested_action",
                "chebi_stars",
                "external_consensus_differs",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .rename(columns={"external_consensus_differs": "external databases agree"})
        .reset_index(drop=True)
    )
    linked_table(stereo_conflicts)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    #### 5.5.2 Stereo defined differently on each side

    Each record defines a stereocentre or double bond that the other leaves undefined, or one defines more and a centre both define is written inverted. Next to an undefined centre, InChI can write the same configuration inverted (5.5), so the configurations are not compared: the two may be the same stereoisomer or different ones. The suggested action is to review them.

    `difference_detail` gives how many centres each record defines and how many are written inverted. It counts the centres and double bonds that either InChI lists; InChI does not list the undefined ones of a molecule or salt component that defines none.
    """)
    return


@app.cell
def _(layer_rows):
    stereo_defined_differently = (
        layer_rows[layer_rows.stereo_relation.eq(STEREO_DEFINED_DIFFERENTLY)]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "action",
                "suggested_action",
                "xref_id_with_gtp_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(stereo_defined_differently)
    return


@app.cell(hide_code=True)
def _(layer_rows):
    _one_sided = layer_rows[layer_rows.more_stereo_in.notna()]
    mo.md(f"""
    #### 5.5.3 Stereo defined on one side only

    In {len(_one_sided)} cross-references, one record defines stereo that the other leaves undefined, and every centre both define is written the same way. GtoPdb defines more in {int(_one_sided.more_stereo_in.eq("GtoPdb").sum())} and the linked record in {int(_one_sided.more_stereo_in.eq("linked record").sum())}; the table below shows that the direction depends on the database. These differences explain why exact InChIKey matching between GtoPdb and the linked databases fails.

    The missing stereo is a stereocentre (R/S), a double bond (E/Z), or ring cis/trans: the relative configuration of ring atoms that are not stereocentres, such as C1 and C4 of a 1,4-disubstituted cyclohexane. `difference_detail` says which record defines more, and `missing_stereo` what the other lacks. The centres defined on one side only are highlighted.
    """)
    return


@app.cell
def _(layer_rows):
    one_sided_stereo_table = pd.crosstab(
        [layer_rows.more_stereo_in, layer_rows.missing_stereo],
        layer_rows.xref_db,
        margins=True,
        margins_name="Total",
    ).rename_axis(index=["more specific record", "stereo the other record lacks"], columns=None)
    one_sided_stereo_table
    return


@app.cell(hide_code=True)
def _(layer_rows):
    _with_record = layer_rows.more_stereo_in.eq("GtoPdb") & layer_rows.xref_id_with_gtp_structure.notna()
    mo.md(f"""
    **GtoPdb defines more, and an exact record exists ({int(_with_record.sum())}).** ChEMBL, ChEBI and PubChem all keep a flat structure and its stereo-defined forms as separate records, so the suggested action is to link the exact record instead.
    """)
    return


@app.cell
def _(layer_rows):
    better_stereo_records = (
        layer_rows[layer_rows.xref_id_with_gtp_structure.notna() & layer_rows.more_stereo_in.eq("GtoPdb")]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "missing_stereo",
                "xref_id_with_gtp_structure",
                "action",
                "suggested_action",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(better_stereo_records)
    return


@app.cell(hide_code=True)
def _(layer_rows):
    _without_record = layer_rows[layer_rows.more_stereo_in.eq("GtoPdb") & layer_rows.xref_id_with_gtp_structure.isna()]
    _pubchem = _without_record[_without_record.xref_db == "PubChem"]
    _removed = _pubchem.pubchem_standardisation_removes.dropna()
    # their InChIKey is the unlabelled compound's, so it was not looked up
    _not_searched = _without_record.radioactive & ~_without_record.gtp_inchi.str.contains("/i")
    _removed_by_kind = ", ".join(f"{kind} {n}" for kind, n in _removed.value_counts().items())
    mo.vstack(
        [
            mo.md(f"""
    **GtoPdb defines more, and no exact record was identified ({len(_without_record)}).** This includes {int(_not_searched.sum())} links where the radioactive ligand was stored without its isotope and the lookup was skipped (section 5.6). For {len(_removed)} of the {len(_pubchem)} PubChem links, PubChem's own standardisation removes the missing stereo from the GtoPdb structure ({_removed_by_kind}; `pubchem_standardisation_removes`), reproducing the linked InChI.
    """),
            mo.accordion(
                {
                    "Which stereo PubChem's standardisation removes": mo.md(r"""
    PubChem's standardisation service was run on each GtoPdb structure. PubChem defines stereocentres by OEChem symmetry classes, which ignore ring cis/trans configuration. It also leaves E/Z undefined where the nitrogen of a double bond carries a carbon substituent (azo compounds, N-substituted imines), and stereo undefined at a phosphorus bonded to more than one O or S that may be charged or double-bonded, as in phosphorothioates ([Hähnke et al. 2018](https://doi.org/10.1186/s13321-018-0293-8)).
    """)
                }
            ),
        ]
    )
    return


@app.cell
def _(layer_rows):
    stereo_without_gtp_record = (
        layer_rows[layer_rows.more_stereo_in.eq("GtoPdb") & layer_rows.xref_id_with_gtp_structure.isna()]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "missing_stereo",
                "pubchem_standardisation_removes",
                "action",
                "suggested_action",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(stereo_without_gtp_record)
    return


@app.cell(hide_code=True)
def _(layer_rows):
    _linked_more = layer_rows[layer_rows.more_stereo_in.eq("linked record")]
    mo.md(f"""
    **The linked record defines more ({len(_linked_more)}).** ChEMBL, ChEBI or PubChem define stereo that GtoPdb leaves undefined. In {int(_linked_more.missing_stereo.str.contains("R/S").sum())} of these, a stereocentre is undefined in GtoPdb, which a deliberately drawn racemate could explain; the rest is ring cis/trans or E/Z geometry, which has no racemate. The suggested action is to check the GtoPdb structure.
    """)
    return


@app.cell
def _(layer_rows):
    stereo_more_in_linked_record = (
        layer_rows[layer_rows.more_stereo_in.eq("linked record")]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "missing_stereo",
                "action",
                "suggested_action",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(stereo_more_in_linked_record)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 5.6 Differs in isotopes

    Only one record carries isotope labels, or the two carry different ones; `difference_detail` says which. A radioactive ligand that GtoPdb stores without its isotope (section 6.1) has the InChIKey of the unlabelled compound, so the linked database could not be searched for an exact record. Without an exact record, a labelled ligand's link to the unlabelled compound is kept, even where that compound is also a GtoPdb ligand, and the suggested action says that it is unlabelled.
    """)
    return


@app.cell
def _(layer_rows):
    isotope_differences = (
        layer_rows[layer_rows.first_difference.eq("isotopes")]
        .sort_values(["xref_db", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "radioactive",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference_detail",
                "action",
                "suggested_action",
                "xref_id_with_gtp_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(isotope_differences)
    return


@app.cell(hide_code=True)
def _(curation):
    _metal = curation[curation.metal_differs]
    _no_structure = int(_metal.first_difference.eq("not compared").sum())
    mo.md(f"""
    ### 5.7 Metal-containing structures

    Standard InChI removes the bonds to metals, writing a complex as its ligands and separate metal ions, and where the bonds to the metal were drawn changes the charges and hydrogens it writes. Auranofin (ligand 6306), for example, is written with Au<sup>+</sup> and a thiolate in GtoPdb, but with neutral gold and a thiol in PubChem, so the two first differ in charge. The layers therefore cannot tell a different compound from the same one drawn differently.

    Every cross-reference where either structure contains a metal and the linked record does not have exactly the GtoPdb InChI is listed here instead of under its layer: {len(_metal)} cross-references on {_metal.ligand_id.nunique()} ligands. The alkali and alkaline-earth metals of salts, such as Na<sup>+</sup> and Ca<sup>2+</sup>, do not count. `difference` still gives the first layer that differs. For {_no_structure} of them, the linked record has no structure at all, so structural agreement is unknown. Without an exact record, the suggested action is to check which structure is intended.

    Cross-references with an identical InChI are not listed, although for a metal complex an identical standard InChI does not confirm the geometry: cisplatin and transplatin share one.
    """)
    return


@app.cell
def _(curation):
    metal_differences = (
        curation[curation.metal_differs]
        .sort_values(
            ["first_difference", "xref_db", "ligand_id"],
            key=lambda column: (
                column.map([*LAYERS, "not compared"].index) if column.name == "first_difference" else column
            ),
        )[
            [
                "ligand_id",
                "name",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "difference",
                "difference_detail",
                "action",
                "suggested_action",
                "xref_id_with_gtp_structure",
                "gtp_inchi",
                "xref_inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(metal_differences)
    return


@app.cell(hide_code=True)
def _(gtp_ligands):
    radioactive_count = int(gtp_ligands.radioactive.sum())
    mo.md(f"""
    ## 6. Structure issues within GtoPdb

    ### 6.1 Radioactive ligands stored without isotopes

    GtoPdb flags {radioactive_count:,} ligands as radioactive (its Labelled flag also covers fluorescent probes, which carry no isotope). Many have an InChI without an `/i` isotope layer. Omitting a label can give a radioligand the same InChIKey as its unlabelled parent. This table lists every InChIKey used by more than one GtoPdb ligand.
    """)
    return


@app.cell
def _(gtp_ligands):
    radioactive_without_isotope = gtp_ligands[
        gtp_ligands.radioactive
        & gtp_ligands.inchi.notna()
        & ~gtp_ligands.inchi.fillna("").str.contains("/i")
    ]
    duplicate_structures = (
        gtp_ligands.dropna(subset=["inchikey"])
        .groupby("inchikey")
        .filter(lambda group: len(group) > 1)
        .groupby("inchikey")
        .agg(
            ligands=("ligand_id", "size"),
            ligand_ids=("ligand_id", lambda ids: ", ".join(map(str, sorted(ids)))),
            names=("name", lambda names: " | ".join(names)),
            includes_radioactive_ligand=("radioactive", "any"),
        )
        .sort_values(["includes_radioactive_ligand", "ligands"], ascending=False)
        .reset_index()
    )
    linked_table(duplicate_structures)
    return duplicate_structures, radioactive_without_isotope


@app.cell(hide_code=True)
def _(radioactive_flag_issues):
    mo.md(f"""
    The flag can also disagree with isotope labels detected in SMILES or superscripted in names. Ground-state stability follows Kondev et al., [NUBASE2020](https://doi.org/10.1088/1674-1137/abddae). {int((~radioactive_flag_issues.radioactive).sum())} ligands have a radioisotope but are not flagged as radioactive, and {int(radioactive_flag_issues.radioactive.sum())} flagged as radioactive have only stable labels among those detected, such as <sup>13</sup>C, <sup>15</sup>N or <sup>19</sup>F, which are not radioactive.
    """)
    return


@app.cell
def _(gtp_ligands, gtp_ligands_raw):
    # names are read too: Ac-[13C,15N]PGP has no structure, so only its name shows its labels
    names_html = dict(zip(gtp_ligands_raw["Ligand ID"].astype(int), gtp_ligands_raw.Name))
    flag_checked = gtp_ligands.assign(
        isotopes=lambda d: [
            isotopes_written(smiles, names_html[ligand_id])
            for smiles, ligand_id in zip(d.smiles, d.ligand_id)
        ],
        radioisotope=lambda d: d.isotopes.map(lambda found: any(map(is_radioisotope, found))),
        stable_only=lambda d: d.isotopes.map(
            lambda found: bool(found) and all(is_radioisotope(label) is False for label in found)
        ),
    )
    radioactive_flag_issues = (
        flag_checked[
            (~flag_checked.radioactive & flag_checked.radioisotope)
            | (
                flag_checked.radioactive & flag_checked.stable_only
            )
        ]
        .assign(isotopes=lambda d: d.isotopes.map(", ".join))
        .sort_values(["radioactive", "ligand_id"])[
            ["ligand_id", "name", "type", "radioactive", "isotopes", "smiles"]
        ]
        .reset_index(drop=True)
    )
    linked_table(radioactive_flag_issues)
    return (radioactive_flag_issues,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 6.2 GtoPdb structures containing radicals

    These atoms have unpaired electrons. This is usually a drawing error, such as a missing hydrogen or charge. The atoms that carry them are highlighted in orange in the drawings, which are made from GtoPdb's SMILES.
    """)
    return


@app.cell
def _(gtp_ligands):
    radical_ligands = (
        gtp_ligands[gtp_ligands.smiles.notna()]
        .assign(radical_electrons=lambda d: d.smiles.map(radical_count))
        .query("radical_electrons > 0")[
            ["ligand_id", "name", "type", "radical_electrons", "smiles"]
        ]
        .reset_index(drop=True)
    )
    drawings = [
        radical_figure(smiles, f"GtoPdb {ligand_id}")
        for ligand_id, smiles in zip(radical_ligands.ligand_id, radical_ligands.smiles)
    ]
    linked_table(
        radical_ligands.assign(structure=[mo.Html(figure) for figure, _ in drawings])[
            ["ligand_id", "name", "type", "radical_electrons", "structure", "smiles"]
        ],
        column_widths={"structure": max(width for _, width in drawings) + 24},
    )
    return


@app.cell(hide_code=True)
def _(rebuilt_inchi):
    mo.md(f"""
    ### 6.3 Ligands with an InChIKey but no InChI

    GtoPdb gives {len(rebuilt_inchi)} ligands a SMILES and an InChIKey but leaves their InChI empty. They are large peptides and oligonucleotides with {rebuilt_inchi.atoms.min():,}–{rebuilt_inchi.atoms.max():,} atoms (counting hydrogens), among the largest structures in GtoPdb. For {int(rebuilt_inchi.inchikey_reproduced.sum())} of them, the standard InChI rebuilt from the SMILES has exactly GtoPdb's InChIKey, so the audit uses that InChI for them. Without it, their cross-references could not be compared.
    """)
    return


@app.cell
def _(gtp_ligands):
    rebuilt_inchi = (
        gtp_ligands[gtp_ligands.inchi_rebuilt | (gtp_ligands.inchi.isna() & gtp_ligands.inchikey.notna())]
        .assign(
            atoms=lambda d: pd.array(
                [heavy_and_hydrogen_count(smiles) for smiles in d.smiles], dtype="Int64"
            ),
            inchikey_reproduced=lambda d: d.inchi.notna(),
        )[["ligand_id", "name", "type", "atoms", "inchikey", "inchikey_reproduced"]]
        .sort_values("atoms")
        .reset_index(drop=True)
    )
    linked_table(rebuilt_inchi)
    return (rebuilt_inchi,)


@app.cell(hide_code=True)
def _(smiles_checked, smiles_inchi_issues):
    # one row per ligand (the table repeats a ligand for each of its PubChem CIDs)
    issue_ligands = smiles_inchi_issues.drop_duplicates("ligand_id")
    # PubChem links of the ligands whose SMILES gives an InChI (an invalid SMILES cannot match anything)
    smiles_follow = smiles_inchi_issues.loc[
        ~smiles_inchi_issues.problem.str.startswith("invalid SMILES")
        & smiles_inchi_issues.pubchem_cid_follows.ne(""),
        "pubchem_cid_follows",
    ]
    mo.md(f"""
    ### 6.4 SMILES that disagree with the InChI

    The GtoPdb website draws each structure from its SMILES (with the CDK); the InChI and InChIKey are separate fields. Regenerating the standard InChI from every SMILES with RDKit reproduces GtoPdb's own InChI for {len(smiles_checked) - len(issue_ligands):,} of {len(smiles_checked):,} ligands. The others are listed below, once per PubChem CID. {int(issue_ligands.problem.str.startswith("invalid SMILES").sum())} SMILES are not valid structures (RDKit's message names the atom whose bonds exceed its valence). The rest describe a different structure from the InChI, mostly double-bond or ring stereo that the InChI defines and the SMILES lacks. `pubchem_cid_follows` shows whether GtoPdb's PubChem compound has the structure of the SMILES or of the InChI: for {int(smiles_follow.eq("SMILES").sum())} of the {len(smiles_follow)} PubChem links of ligands whose SMILES gives an InChI, it is the SMILES, consistent with GtoPdb's PubChem deposit being made from the SMILES.

    This check cannot find a SMILES that is wrong in the same way as its InChI. Diadenosine diphosphate (ligand 1065) is an example: its SMILES, InChI and IUPAC name all describe the 5',3'-linked compound, while its synonyms and its ChEBI entry are the 5',5'-linked one (section 5.2).
    """)
    return


@app.cell
def _(curation, gtp_ligands):
    # GtoPdb's SMILES against its own InChI (the rebuilt InChIs of section 6.3 come from the SMILES, so they are skipped)
    smiles_checked = gtp_ligands[
        gtp_ligands.inchi.notna() & gtp_ligands.smiles.notna() & ~gtp_ligands.inchi_rebuilt
    ]
    smiles_results = pd.DataFrame(
        [
            smiles_inchi_problem(smiles, inchi)
            for smiles, inchi in zip(smiles_checked.smiles, smiles_checked.inchi)
        ],
        columns=["smiles_inchi", "problem", "suggested_action"],
        index=smiles_checked.index,
    )
    pubchem_link_inchis = curation.loc[
        (curation.xref_db == "PubChem") & curation.xref_inchi.notna(),
        ["ligand_id", "xref_id_normalised", "xref_inchi"],
    ]
    smiles_inchi_issues = (
        smiles_checked.join(smiles_results)
        .query("problem != ''")
        .merge(pubchem_link_inchis, on="ligand_id", how="left")
        .assign(
            pubchem_cid_follows=lambda d: [
                ""
                if not isinstance(xref_inchi, str)
                else "SMILES"
                if xref_inchi == smiles_inchi
                else "InChI"
                if xref_inchi == inchi
                else "neither"
                for xref_inchi, smiles_inchi, inchi in zip(d.xref_inchi, d.smiles_inchi, d.inchi)
            ]
        )
        .rename(columns={"xref_id_normalised": "pubchem_cid"})
        .sort_values(["suggested_action", "ligand_id"])[
            [
                "ligand_id",
                "name",
                "type",
                "problem",
                "suggested_action",
                "pubchem_cid",
                "pubchem_cid_follows",
                "smiles",
                "inchi",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(smiles_inchi_issues)
    return smiles_checked, smiles_inchi_issues


@app.cell(hide_code=True)
def _(
    gtp_ligands,
    inorganic_links,
    inorganic_structures,
    inorganic_without_structure,
):
    inorganic_complexes = inorganic_without_structure[
        ~inorganic_without_structure.ligand_id.isin(inorganic_structures.ligand_id)
    ]
    # per ligand with cross-references: whether every linked record has the structure shown
    inorganic_links_agree = inorganic_links.groupby("ligand_id").has_structure.all()
    mo.md(f"""
    ### 6.5 Inorganic ligands without a structure

    GtoPdb has {int(gtp_ligands.type.eq("Inorganic").sum())} inorganic ligands. {len(inorganic_without_structure)} of them have no SMILES, InChI or InChIKey, so none of their cross-references can be compared. {len(inorganic_structures)} of these are monatomic ions, small molecules or simple salts whose GtoPdb name is their formula, so the structure follows from the name. Monatomic ions such as Ca<sup>2+</sup> and Cl<sup>-</sup> are read from the name; the other formulas are written out as SMILES in `inorganic_smiles` (appendix). The standard InChI and InChIKey were generated from these SMILES with RDKit.

    For {int(inorganic_links_agree.sum())} of these ligands, every ChEMBL, ChEBI and PubChem record that GtoPdb links already has exactly this InChI. The others have no cross-reference ({len(inorganic_structures) - len(inorganic_links_agree)}) or link a record with a different structure ({int((~inorganic_links_agree).sum())}; see `linked_records`). `pubchem_cid_with_structure` lists the PubChem compounds with exactly the InChI shown, found by InChIKey.

    The other {len(inorganic_complexes)} ({", ".join(sorted(inorganic_complexes.name, key=str.lower))}) are metal complexes. Standard InChI disconnects bonds to metals, so it cannot distinguish, for example, cisplatin from transplatin. They are not given a structure here.
    """)
    return


@app.cell
def _(CACHE_DIR, curation, gtp_ligands, refresh):
    inorganic_without_structure = gtp_ligands[
        gtp_ligands.type.eq("Inorganic")
        & gtp_ligands.smiles.isna()
        & gtp_ligands.inchi.isna()
        & gtp_ligands.inchikey.isna()
    ]
    inorganic_proposed = (
        inorganic_without_structure.assign(smiles=lambda d: d.name.map(inorganic_smiles))
        .dropna(subset=["smiles"])
        .assign(
            inchi=lambda d: [Chem.MolToInchi(Chem.MolFromSmiles(smiles)) for smiles in d.smiles],
            inchikey=lambda d: d.inchi.map(rdinchi.InchiToInchiKey),
        )
    )
    # PubChem compounds with exactly the proposed InChI (looked up by InChIKey), and their names
    inorganic_pubchem_lookup = fetch_cached(
        CACHE_DIR / "pubchem_inorganic_by_inchikey.parquet",
        sorted(inorganic_proposed.inchikey),
        fetch_pubchem_by_inchikey,
        batch_size=10,
        refresh=refresh.value,
        title="PubChem InChIKey look-up (inorganic ligands)",
    ).reindex(columns=["query_id", "found", "xref_id", "inchi", "xref_name", "retrieved_at"])
    pubchem_hits = inorganic_pubchem_lookup.dropna(subset=["xref_id"])
    proposed_inchi = dict(zip(inorganic_proposed.inchikey, inorganic_proposed.inchi))
    cids_with_structure = (
        pubchem_hits[pubchem_hits.inchi.eq(pubchem_hits.query_id.map(proposed_inchi))]
        .groupby("query_id")
        .xref_id.agg(lambda ids: sorted(ids, key=int))
    )
    cid_names = dict(zip(pubchem_hits.xref_id, pubchem_hits.xref_name))
    # GtoPdb's ChEMBL, ChEBI and PubChem links of these ligands, against the structure from the name
    inorganic_links = (
        curation.loc[
            curation.ligand_id.isin(inorganic_proposed.ligand_id),
            ["ligand_id", "xref_db", "xref_id_as_written", "xref_id_normalised", "xref_inchi"],
        ]
        .merge(inorganic_proposed[["ligand_id", "inchi"]], on="ligand_id")
        .assign(has_structure=lambda d: d.xref_inchi.eq(d.inchi))
    )
    linked_records = (
        inorganic_links.assign(
            record=lambda d: [
                f"{db} {xref_id}: "
                + (
                    "this structure"
                    if same
                    else inchi
                    if isinstance(inchi, str)
                    else "no structure"
                )
                for db, xref_id, inchi, same in zip(
                    d.xref_db, d.xref_id_as_written, d.xref_inchi, d.has_structure
                )
            ]
        )
        .groupby("ligand_id")
        .record.agg(" | ".join)
    )
    linked_cids = (
        inorganic_links[inorganic_links.xref_db.eq("PubChem")]
        .groupby("ligand_id")
        .xref_id_normalised.agg(list)
    )
    inorganic_structures = (
        inorganic_proposed.assign(
            pubchem_cid_with_structure=lambda d: [
                "|".join(cids_with_structure.get(key, [])) for key in d.inchikey
            ],
            linked_records=lambda d: d.ligand_id.map(linked_records).fillna("none"),
            suggested_action=lambda d: [
                inorganic_structure_action(
                    linked_cids.get(ligand_id, []), cids_with_structure.get(key, []), cid_names
                )
                for ligand_id, key in zip(d.ligand_id, d.inchikey)
            ],
        )
        # PubChem links to add or replace first, then by name
        .sort_values("name", key=lambda names: names.str.lower())
        .sort_values(
            "suggested_action", key=lambda actions: ~actions.str.contains("link"), kind="stable"
        )[
            [
                "ligand_id",
                "name",
                "smiles",
                "inchi",
                "inchikey",
                "pubchem_cid_with_structure",
                "linked_records",
                "suggested_action",
            ]
        ]
        .reset_index(drop=True)
    )
    linked_table(inorganic_structures)
    return (
        inorganic_links,
        inorganic_pubchem_lookup,
        inorganic_structures,
        inorganic_without_structure,
    )


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## 7. Why GtoPdb ligands do not match in UniChem

    UniChem identifies a compound (UCI) by its full standard InChI and InChIKey. If UniChem held GtoPdb's structures as GtoPdb publishes them, a GtoPdb ligand and its correct ChEMBL, ChEBI and PubChem records would share a UCI.

    UniChem's API gives, for any compound ID of a source, the UCI and the standard InChI that UniChem holds for it. Every GtoPdb ligand was looked up (section 1), and the InChI UniChem holds was compared with the InChI in GtoPdb's `ligands.csv`:
    """)
    return


@app.cell
def _(gtp_ligands, unichem_gtp_records, unichem_pairs):
    # the IDs UniChem's mapping files pair with each GtoPdb ligand (section 7.2)
    unichem_partner_ids = pd.DataFrame(
        {
            f"linked_{name.lower()}": pairs.groupby("gtopdb_id").xref_id.agg(
                lambda ids: "|".join(sorted(ids))
            )
            for name, pairs in unichem_pairs.items()
        }
    )
    unichem_gtp = (
        gtp_ligands[gtp_ligands.inchi.notna()]
        .assign(gtopdb_id=lambda d: d.ligand_id.astype(str))
        .merge(
            unichem_gtp_records.drop(columns="retrieved_at").rename(
                columns={
                    "query_id": "gtopdb_id",
                    "found": "in_unichem",
                    "inchi": "unichem_inchi",
                    "inchikey": "unichem_inchikey",
                }
            ),
            on="gtopdb_id",
            how="left",
        )
        .merge(unichem_partner_ids, left_on="gtopdb_id", right_index=True, how="left")
        .assign(
            in_unichem=lambda d: d.in_unichem.eq(True),
            uci=lambda d: d.uci.astype("Int64"),
            gtp_has_stereo=lambda d: [strip_layers(i, STEREO_LAYERS) != i for i in d.inchi],
            unichem_vs_gtp=lambda d: [
                unichem_status(gtp_inchi, unichem_inchi)
                for gtp_inchi, unichem_inchi in zip(d.inchi, d.unichem_inchi)
            ],
        )
    )
    return (unichem_gtp,)


@app.cell
def _(unichem_gtp):
    unichem_status_table = pd.crosstab(
        unichem_gtp.unichem_vs_gtp,
        unichem_gtp.gtp_has_stereo.map(STEREO_LABELS),
        margins=True,
        margins_name="Total",
    ).rename_axis(index="UniChem structure vs GtoPdb", columns=None)
    unichem_status_table
    return


@app.cell
def _(unichem_gtp):
    # the GtoPdb ligands with stereochemistry that UniChem holds: which of their layers UniChem keeps
    held_stereo_ligands = unichem_gtp[unichem_gtp.gtp_has_stereo & unichem_gtp.in_unichem]
    gtp_layers = held_stereo_ligands.inchi.map(parse_standard_inchi)
    unichem_layers = held_stereo_ligands.unichem_inchi.map(parse_standard_inchi)
    unichem_layer_table = pd.DataFrame(
        [
            {
                "InChI layer": label,
                "GtoPdb ligands with the layer": sum(key in layers for layers in gtp_layers),
                "of these, with the layer in UniChem": sum(
                    key in gtp and key in held for gtp, held in zip(gtp_layers, unichem_layers)
                ),
            }
            for key, label in [
                ("t", "/t tetrahedral stereo"),
                ("b", "/b double-bond (E/Z) stereo"),
                ("i", "/i isotopes"),
            ]
        ]
    )
    unichem_layer_table
    return


@app.cell(hide_code=True)
def _(unichem_gtp):
    held_ligands = unichem_gtp[unichem_gtp.in_unichem]
    held_with_stereo = held_ligands[held_ligands.gtp_has_stereo]
    stereo_kept = sum(strip_layers(i, STEREO_LAYERS) != i for i in held_with_stereo.unichem_inchi)
    flat_key_share = (
        held_with_stereo.unichem_inchikey.str.split("-").str[1].eq("UHFFFAOYSA").mean()
    )
    not_held = unichem_gtp[~unichem_gtp.in_unichem]
    not_held_names = ", ".join(
        f"{name} (ligand {ligand_id})" for ligand_id, name in zip(not_held.ligand_id, not_held.name)
    )
    other_differences = int((unichem_gtp.unichem_vs_gtp == "other difference").sum())
    mo.md(f"""
    UniChem holds {len(held_ligands):,} of the {len(unichem_gtp):,} GtoPdb ligands that have an InChI{f"; it does not hold {not_held_names}" if len(not_held) else ""}. Among the {int((~held_ligands.gtp_has_stereo).sum()):,} held ligands without stereo, {int((~held_ligands.gtp_has_stereo & held_ligands.unichem_vs_gtp.eq("identical")).sum()):,} have identical InChIs. Of the **{len(held_with_stereo):,}** with stereochemistry, {"**none** keeps" if stereo_kept == 0 else f"only **{stereo_kept:,}** keep"} any stereo layer in UniChem: tetrahedral centres and E/Z double bonds are both gone, while isotope labels are kept (second table). The {other_differences} *other differences* are explained in section 7.1.

    **UniChem does use full InChIKeys.** The second block of a standard InChIKey hashes the stereo and isotope layers. For {flat_key_share:.1%} of the GtoPdb ligands with stereochemistry, UniChem's key has `UHFFFAOYSA` there, the hash of "no stereo or isotope layers", because UniChem holds no stereo for them and they carry no isotope label.
    """)
    return


@app.cell
def _(CACHE_DIR, chembl_records, refresh, unichem_gtp):
    gtp_valsartan = unichem_gtp.set_index("ligand_id").loc[3937]
    chembl_valsartan = chembl_records.set_index("query_id").loc["CHEMBL1069"]
    unichem_chembl_valsartan = unichem_compound("ChEMBL", "CHEMBL1069", CACHE_DIR, refresh.value)
    valsartan_example = pd.DataFrame(
        [
            {
                "record": "GtoPdb ligand 3937 (ligands.csv)",
                "uci": None,
                "standard InChIKey": gtp_valsartan.inchikey,
                "stereo layers": stereo_layers_of(gtp_valsartan.inchi),
            },
            {
                "record": "GtoPdb ligand 3937 in UniChem",
                "uci": gtp_valsartan.uci,
                "standard InChIKey": gtp_valsartan.unichem_inchikey,
                "stereo layers": stereo_layers_of(gtp_valsartan.unichem_inchi),
            },
            {
                "record": "CHEMBL1069 in ChEMBL",
                "uci": None,
                "standard InChIKey": chembl_valsartan.inchikey,
                "stereo layers": stereo_layers_of(chembl_valsartan.inchi),
            },
            {
                "record": "CHEMBL1069 in UniChem",
                "uci": unichem_chembl_valsartan["uci"],
                "standard InChIKey": unichem_chembl_valsartan["inchikey"],
                "stereo layers": stereo_layers_of(unichem_chembl_valsartan["inchi"]),
            },
        ]
    ).astype({"uci": "Int64"})
    valsartan_example
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Valsartan shows what this does to UniChem's links. ChEMBL's CHEMBL1069 is held in UniChem with its stereocentre (`/t22-/m0/s1`), but GtoPdb ligand 3937, the same compound, is held without it. The two therefore get different InChIKeys and different UCIs, and UniChem cannot link them.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 7.1 What UniChem's GtoPdb structures are built from

    UniChem publishes the structures it holds, not the files its sources supply (it asks them for an SDF file or a SMILES table: [UniChem provider examples](https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/examples/README)). GtoPdb publishes each structure twice, as a SMILES and as an InChI, and for a few ligands the two disagree (section 6.4). These ligands show which of the two structures UniChem holds. The comparison leaves out stereo, since UniChem holds none for GtoPdb.
    """)
    return


@app.cell
def _(smiles_inchi_issues, unichem_gtp):
    unichem_inchi_of = unichem_gtp.set_index("ligand_id").unichem_inchi
    unichem_structure_check = (
        smiles_inchi_issues.drop_duplicates("ligand_id")
        .assign(
            smiles_inchi=lambda d: [
                smiles_inchi_problem(smiles, inchi)[0] for smiles, inchi in zip(d.smiles, d.inchi)
            ],
            unichem_holds=lambda d: [
                unichem_structure_source(inchi, smiles_inchi, unichem_inchi_of.get(ligand_id))
                for ligand_id, inchi, smiles_inchi in zip(d.ligand_id, d.inchi, d.smiles_inchi)
            ],
        )
        .sort_values(["unichem_holds", "ligand_id"])[
            ["ligand_id", "name", "problem", "unichem_holds"]
        ]
        .reset_index(drop=True)
    )
    linked_table(unichem_structure_check)
    return (unichem_structure_check,)


@app.cell(hide_code=True)
def _(unichem_structure_check):
    # only where the SMILES and InChI differ in more than stereo can just one of them match UniChem
    decisive = unichem_structure_check[
        unichem_structure_check.unichem_holds.isin(["SMILES", "InChI", "neither"])
    ]
    follow_smiles = decisive[decisive.unichem_holds == "SMILES"]
    decisive_names = ", ".join(
        f"{name} (ligand {ligand_id})" for ligand_id, name in zip(decisive.ligand_id, decisive.name)
    )
    mo.md(f"""
    For {len(decisive)} of these ligands, the SMILES and InChI differ in more than stereo{f": {decisive_names}" if len(decisive) else ""}. Only one of the two can then match UniChem's structure, and it is the SMILES for **{len(follow_smiles)} of {len(decisive)}**. The other ligands cannot tell the two apart: their SMILES and InChI differ only in stereo, or the SMILES is invalid and UniChem holds the InChI's structure.

    This comparison identifies which published representation matches in the informative cases. The ChEMBL record for valsartan retains its stereo (CHEMBL1069 above). Whether the stereo is already missing from the file GtoPdb supplies or is lost when UniChem loads it cannot be seen from the public data.

    **Suggested fix:** GtoPdb should check that the file it supplies to UniChem carries the stereo: wedge bonds with 2D coordinates, or atom parities, in an SDF file, or isomeric SMILES in a SMILES table. GtoPdb should then ask UniChem to reload source 4.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### 7.2 Effect on UniChem's links

    For cross-references whose structure is identical to GtoPdb's (section 3, *identical*), this table shows whether UniChem's mapping files pair the GtoPdb ligand with the same identifier. Where GtoPdb has stereochemistry, UniChem typically links GtoPdb to the unspecified-stereo entry (for valsartan, PubChem CID 5650 rather than 60846) or to nothing.
    """)
    return


@app.cell
def _(curation, unichem_gtp):
    unichem_links = curation[curation.first_difference == "identical"].merge(
        unichem_gtp[
            ["ligand_id", "gtp_has_stereo", "linked_chembl", "linked_chebi", "linked_pubchem"]
        ],
        on="ligand_id",
        how="inner",
    )
    unichem_links["linked_by_unichem"] = [
        xref_id in str(row[f"linked_{db.lower()}"]).split("|")
        for xref_id, db, row in zip(
            unichem_links.xref_id_normalised,
            unichem_links.xref_db,
            unichem_links.to_dict("records"),
        )
    ]
    unichem_link_table = (
        unichem_links.groupby(
            [
                "xref_db",
                unichem_links.gtp_has_stereo.map(STEREO_LABELS),
            ]
        )
        .linked_by_unichem.agg(curated_links="size", also_linked_by_unichem="sum")
        .assign(share=lambda d: (d.also_linked_by_unichem / d.curated_links).round(3))
        .rename_axis(index=["linked database", "GtoPdb ligand"])
    )
    unichem_link_table
    return (unichem_links,)


@app.cell
def _(unichem_gtp):
    # UCIs that UniChem files more than one GtoPdb ligand under
    shared_uci_ligands = unichem_gtp[
        unichem_gtp.in_unichem & unichem_gtp.uci.duplicated(keep=False)
    ].assign(gtp_inchi_without_stereo=lambda d: d.inchi.map(lambda i: strip_layers(i, STEREO_LAYERS)))
    shared_ucis = shared_uci_ligands.groupby("uci").agg(
        gtp_inchis=("inchi", "nunique"),
        gtp_inchis_without_stereo=("gtp_inchi_without_stereo", "nunique"),
        with_radioactive=("radioactive", "any"),
    )
    # GtoPdb keeps these ligands apart only by their stereo
    stereo_merged_ucis = shared_ucis.index[
        (shared_ucis.gtp_inchis > 1) & (shared_ucis.gtp_inchis_without_stereo == 1)
    ]
    merged_stereoisomers = (
        shared_uci_ligands[shared_uci_ligands.uci.isin(stereo_merged_ucis)]
        .assign(gtp_stereo=lambda d: d.inchi.map(stereo_layers_of))
        .sort_values(["uci", "ligand_id"])[["uci", "ligand_id", "name", "inchikey", "gtp_stereo"]]
        .reset_index(drop=True)
    )
    return merged_stereoisomers, shared_uci_ligands, shared_ucis


@app.cell(hide_code=True)
def _(merged_stereoisomers, shared_uci_ligands, shared_ucis):
    same_gtp_inchi = shared_ucis[shared_ucis.gtp_inchis == 1]
    other_merges = int((shared_ucis.gtp_inchis_without_stereo > 1).sum())
    with_undefined = int(
        merged_stereoisomers.groupby("uci").gtp_stereo.agg(lambda s: s.eq("none").any()).sum()
    )
    n_stereo_merged = merged_stereoisomers.uci.nunique()
    mo.md(f"""
    **GtoPdb ligands that UniChem holds as one compound.** UniChem files {len(shared_ucis):,} UCIs under more than one GtoPdb ligand ({len(shared_uci_ligands):,} ligands). In {len(same_gtp_inchi)} of them GtoPdb's own InChIs are identical too; {int(same_gtp_inchi.with_radioactive.sum())} of these include a ligand flagged as radioactive, the duplicates of section 6.1. In the other **{n_stereo_merged}**, GtoPdb's InChIs differ only in stereo, so UniChem merges stereoisomers that GtoPdb keeps apart{f" ({other_merges} UCIs hold GtoPdb ligands that differ in more than stereo)" if other_merges else ""}. In {with_undefined} of them one ligand leaves the stereo undefined and another defines it, such as (±)-YM796 (ligand 285) and (−)-YM796 (ligand 287); in the other {n_stereo_merged - with_undefined}, every ligand carries stereo information and the group contains at least two distinct InChIs. They are listed below, one row per ligand.
    """)
    return


@app.cell
def _(merged_stereoisomers):
    linked_table(merged_stereoisomers)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## 8. Output files

    Every run writes the tables below. `curation_issues.csv` is the file for curators: one row per cross-reference whose identifier has a problem or whose structure differs from GtoPdb's, with links to both records. Rows follow the order of sections 4 and 5: identifier issues on structures that could not be compared or that agree with GtoPdb first, then the layers in InChI order, then the metal-containing structures. Its main columns:

    | Column | Meaning |
    |---|---|
    | `first_difference` | The first InChI layer that differs (section 2) |
    | `difference` | The layers both records share up to that point and the one that differs, such as *Shares formula, connections and hydrogens; differs in stereo*. Layers neither record has, such as the charge of two neutral compounds, are not named |
    | `difference_detail` | How that layer differs; each part of section 5 says what it gives |
    | `stereo_relation`, `more_stereo_in`, `missing_stereo` | How the stereo differs, which record defines more, and what the other lacks (section 5.5) |
    | `metal_differs` | Whether either structure contains a metal and the InChIs are not identical (section 5.7) |
    | `xref_id_with_gtp_structure` | The exact record: a record in the linked database with exactly the GtoPdb structure |
    | `action` | The kind of each step of the suggested action, such as *keep link* or *check GtoPdb structure* |
    | `suggested_xref_id` | The ID to link instead, if any |
    | `suggested_action` | The steps and notes, as sentences |
    | `tanimoto_similarity` | Morgan (radius 2) fingerprint similarity of the two structures |
    | `other_gtp_ligands_with_xref_structure` | Other GtoPdb ligands whose structure matches the linked record |
    | `external_consensus_differs` | Whether two or more external databases agree with each other but not with GtoPdb |
    """)
    return


@app.cell
def _(curation):
    curation_issues = (
        curation[
            curation.first_difference.isin(LAYERS)
            | curation.metal_differs
            | curation.identifier_issue.ne("")
            | curation.xref_id_with_gtp_structure.notna()
        ]
        .assign(
            gtp_url=lambda d: d.ligand_id.map(GTP_LIGAND_URL.format),
            xref_url=lambda d: [
                XREF_URLS[db].format(i) if isinstance(i, str) else ""
                for db, i in zip(d.xref_db, d.xref_id_normalised)
            ],
        )
        # in the order of sections 4 and 5, with the metal-containing structures of 5.7 last
        .sort_values(
            ["metal_differs", "first_difference", "xref_db", "ligand_id"],
            key=lambda column: (
                column.map(["not compared", "identical", *LAYERS].index)
                if column.name == "first_difference"
                else column
            ),
        )[
            [
                "ligand_id",
                "name",
                "type",
                "xref_db",
                "xref_id_as_written",
                "xref_name",
                "chebi_stars",
                "first_difference",
                "difference",
                "difference_detail",
                "stereo_relation",
                "more_stereo_in",
                "missing_stereo",
                "pubchem_standardisation_removes",
                "metal_differs",
                "xref_id_with_gtp_structure",
                "xref_name_with_gtp_structure",
                "identifier_issue",
                "action",
                "suggested_xref_id",
                "suggested_action",
                "tanimoto_similarity",
                "external_consensus_differs",
                "other_gtp_ligands_with_xref_structure",
                "gtp_inchikey",
                "xref_inchikey",
                "gtp_inchi",
                "xref_inchi",
                "gtp_url",
                "xref_url",
            ]
        ]
    )
    return (curation_issues,)


@app.cell
def _(
    comparison,
    consensus_against_gtp,
    curation,
    duplicate_structures,
    gtp_ligands,
    gtp_version,
    inorganic_structures,
    inorganic_without_structure,
    merged_stereoisomers,
    radioactive_flag_issues,
    radioactive_without_isotope,
    retrieval_summary,
    smiles_checked,
    smiles_inchi_issues,
    unichem_gtp,
    unichem_links,
    unichem_structure_check,
):
    compared = curation[curation.first_difference != "not compared"]
    # the layer counts leave out the metal-containing structures, which have their own (section 5.7)
    layered = curation[~curation.metal_differs]
    formula_or_connections = layered.first_difference.isin(["formula", "connections"])
    with_stereo = unichem_gtp[unichem_gtp.gtp_has_stereo & unichem_gtp.in_unichem]
    summary = {
        "gtopdb_version": gtp_version,
        "retrieved_on": retrieval_summary.retrieved.str[:10].max(),
        "xrefs_total": len(comparison),
        "ligands_with_xrefs": int(comparison.ligand_id.nunique()),
        "xrefs_compared": len(compared),
        "xrefs_identical": int((compared.first_difference == "identical").sum()),
        "xrefs_by_first_difference": {
            layer: int((layered.first_difference == layer).sum()) for layer in LAYERS
        },
        "metal_differs": int(curation.metal_differs.sum()),
        "metal_compared": int(compared.metal_differs.sum()),
        "metal_differs_ligands": int(curation.loc[curation.metal_differs, "ligand_id"].nunique()),
        "pubchem_agree": int(
            ((curation.xref_db == "PubChem") & (curation.first_difference == "identical")).sum()
        ),
        "pubchem_agree_from_gtp_deposit": int(
            (
                curation.pubchem_cid_from_gtp_deposit & (curation.first_difference == "identical")
            ).sum()
        ),
        "ligands_only_confirmed_by_gtp_deposit": int(
            curation.loc[
                curation.only_gtp_deposit_agrees & curation.first_difference.isin(LAYERS),
                "ligand_id",
            ].nunique()
        ),
        "consensus_ligands": len(consensus_against_gtp),
        "identifier_issues": int((curation.identifier_issue != "").sum()),
        "chembl_unresolved": int(((curation.xref_db == "ChEMBL") & ~curation.xref_found).sum()),
        "chembl_unresolved_with_suggestion": int(
            (
                (curation.xref_db == "ChEMBL")
                & ~curation.xref_found
                & curation.xref_id_with_gtp_structure.notna()
            ).sum()
        ),
        "stereo_defined_differently": int(layered.stereo_relation.eq(STEREO_DEFINED_DIFFERENTLY).sum()),
        "stereo_one_side": int(layered.more_stereo_in.notna().sum()),
        "stereo_better_record": int(
            (
                layered.more_stereo_in.eq("GtoPdb") & layered.xref_id_with_gtp_structure.notna()
            ).sum()
        ),
        "formula_or_connections": int(formula_or_connections.sum()),
        "formula_or_connections_ligands": int(layered.loc[formula_or_connections, "ligand_id"].nunique()),
        "formula_or_connections_by_db": ", ".join(
            f"{db} {n}"
            for db, n in layered.loc[formula_or_connections, "xref_db"].value_counts().sort_index().items()
        ),
        "formula_or_connections_with_gtp_structure_id": int(
            (formula_or_connections & layered.xref_id_with_gtp_structure.notna()).sum()
        ),
        "stereo_conflicts": int(layered.stereo_relation.isin(STEREO_CONFLICT_RELATIONS).sum()),
        "stereo_conflicts_with_gtp_stereo_id": int(
            (
                layered.stereo_relation.isin(STEREO_CONFLICT_RELATIONS)
                & layered.xref_id_with_gtp_structure.notna()
            ).sum()
        ),
        "radioactive_without_isotope": len(radioactive_without_isotope),
        "duplicate_keys": len(duplicate_structures),
        "duplicate_keys_with_radioactive_ligand": int(duplicate_structures.includes_radioactive_ligand.sum()),
        "unflagged_radioactive": int((~radioactive_flag_issues.radioactive).sum()),
        "radioactive_with_stable_isotopes_only": int(radioactive_flag_issues.radioactive.sum()),
        "inchi_rebuilt": int(gtp_ligands.inchi_rebuilt.sum()),
        "smiles_checked": len(smiles_checked),
        "smiles_inchi_issues": int(smiles_inchi_issues.ligand_id.nunique()),
        "invalid_smiles": int(
            smiles_inchi_issues.drop_duplicates("ligand_id")
            .problem.str.startswith("invalid SMILES")
            .sum()
        ),
        "inorganic_without_structure": len(inorganic_without_structure),
        "inorganic_structures": len(inorganic_structures),
        "inorganic_pubchem_link_changes": int(
            inorganic_structures.suggested_action.str.contains("link").sum()
        ),
        "unichem_with_stereo": len(with_stereo),
        "stereo_exact_links": int(unichem_links.gtp_has_stereo.sum()),
        "stereo_exact_links_in_unichem": int((unichem_links.gtp_has_stereo & unichem_links.linked_by_unichem).sum()),
        "unichem_stereo_stripped": int(
            (with_stereo.unichem_vs_gtp == "stereo layers missing in UniChem").sum()
        ),
        "unichem_with_stereo_kept": sum(
            strip_layers(i, STEREO_LAYERS) != i for i in with_stereo.unichem_inchi
        ),
        "unichem_held": int(unichem_gtp.in_unichem.sum()),
        "unichem_ucis_merging_stereoisomers": int(merged_stereoisomers.uci.nunique()),
        "unichem_structure_follows_smiles": int(
            (unichem_structure_check.unichem_holds == "SMILES").sum()
        ),
        "unichem_structure_decisive": int(
            unichem_structure_check.unichem_holds.isin(["SMILES", "InChI", "neither"]).sum()
        ),
    }
    return (summary,)


@app.cell
def _(
    OUTPUT_DIR,
    curation,
    curation_issues,
    duplicate_structures,
    gtp_version_number,
    inorganic_structures,
    merged_stereoisomers,
    smiles_inchi_issues,
    summary,
    unichem_gtp,
):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    prefix = f"gtopdb_{gtp_version_number}"
    output_tables = {
        f"{prefix}__curation_issues.csv": curation_issues,
        f"{prefix}__xref_structure_comparison.csv": curation.drop(columns=["smiles"]),
        f"{prefix}__unichem_gtopdb_structures.csv": unichem_gtp.drop(columns=["smiles"]),
        f"{prefix}__duplicate_structures.csv": duplicate_structures,
        f"{prefix}__smiles_inchi_issues.csv": smiles_inchi_issues,
        f"{prefix}__inorganic_structures.csv": inorganic_structures,
        f"{prefix}__unichem_merged_stereoisomers.csv": merged_stereoisomers,
    }
    for filename, table in output_tables.items():
        table.to_csv(OUTPUT_DIR / filename, index=False)
    summary_path = OUTPUT_DIR / f"{prefix}__summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    written_rows = {OUTPUT_DIR / name: len(table) for name, table in output_tables.items()}
    written_files = pd.DataFrame(
        [
            {"file": str(path.relative_to(mo.notebook_dir())), "rows": rows}
            for path, rows in (written_rows | {summary_path: None}).items()
        ]
    )
    written_files
    return (prefix,)


@app.cell
def _(curation_issues, prefix):
    mo.download(
        data=curation_issues.to_csv(index=False).encode(),
        filename=f"{prefix}__curation_issues.csv",
        mimetype="text/csv",
        label="Download curation issues (CSV)",
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Appendix: the functions used above

    The notebook defines every function it uses here, so it runs on its own (in molab, or as a single file). They are marimo top-level functions, so `from gtp_chemistry_problem import classify_inchi_pair` imports them, as the tests in `tests/` do. Shared constants, such as the InChI layers and the identifier patterns, are in the setup cell at the top.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### Retrieval

    Downloads and API calls (GtoPdb, UniChem, ChEMBL, ChEBI, PubChem), cached under `data/raw_data/`. Where that cache is missing, `seed_cache` first copies the one committed to the GitHub repository.
    """)
    return


@app.function
def http_session():
    """requests session with retry/back-off for the public APIs used here."""
    session = requests.Session()
    session.headers["User-Agent"] = "gtopdb-xref-audit (marimo notebook)"
    retry = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=None,
    )
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


@app.function
def utc_now():
    """Current UTC time as an ISO-8601 string."""
    return dt.datetime.now(dt.UTC).isoformat(timespec="seconds")


@app.function
def seed_cache(cache_dir, base_url, filenames):
    """Copy the `filenames` missing from `cache_dir` from `base_url` (the cache committed to the GitHub repository), so
    that the notebook on its own reproduces the published results without querying every API again. A file that
    cannot be fetched is left to the downloads and API calls. Returns the names copied."""
    missing = [name for name in filenames if not (cache_dir / name).exists()]
    if not missing:
        return []
    cache_dir.mkdir(parents=True, exist_ok=True)
    session = http_session()
    copied = []
    with mo.status.progress_bar(
        total=len(missing), title="Copying the cached data from GitHub", remove_on_exit=True
    ) as bar:
        for name in missing:
            try:
                response = session.get(base_url + name, timeout=120)
                response.raise_for_status()
            except requests.RequestException:
                bar.update()
                continue
            # written under another name first, so that an interrupted copy leaves no partial file
            partial = cache_dir / f"{name}.part"
            partial.write_bytes(response.content)
            partial.replace(cache_dir / name)
            copied.append(name)
            bar.update()
    return copied


@app.function
def download_to_cache(url, path, refresh=False):
    """Download `url` to `path` unless cached (with its download log), and return its download log."""
    log_path = path.with_name(path.name + ".log.json")
    if refresh or not path.exists() or not log_path.exists():
        with mo.status.spinner(title=f"Downloading {url}"):
            response = http_session().get(url, timeout=300)
        response.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(response.content)
        log_path.write_text(
            json.dumps(
                {
                    "downloaded_from": url,
                    "downloaded_at_utc": utc_now(),
                    "last_modified": response.headers.get("Last-Modified"),
                    "file_size_bytes": len(response.content),
                    "sha256": hashlib.sha256(response.content).hexdigest(),
                },
                indent=2,
            )
        )
    return json.loads(log_path.read_text())


@app.function
def download_gtp_file(filename, cache_dir, refresh=False):
    """A GtoPdb data file as (table, version line, download log)."""
    log = download_to_cache(
        f"https://www.guidetopharmacology.org/DATA/{filename}", cache_dir / filename, refresh
    )
    with (cache_dir / filename).open(encoding="utf-8") as handle:
        version = handle.readline().strip().strip('"').lstrip("# ")
    table = pd.read_csv(cache_dir / filename, skiprows=1, dtype=str, keep_default_na=False)
    return table, version, log


@app.function
def download_unichem_mapping(source_a, source_b, cache_dir, refresh=False):
    """(ID pairs of two sources that share a UniChem UCI, columns named after the sources; download log)."""
    names = {UNICHEM_SOURCE_IDS[source_a]: source_a, UNICHEM_SOURCE_IDS[source_b]: source_b}
    low, high = sorted(names)
    filename = f"src{low}src{high}.txt.gz"
    path = cache_dir / f"unichem_{filename}"
    log = download_to_cache(
        f"https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id{low}/{filename}",
        path,
        refresh,
    )
    pairs = pd.read_csv(
        path, sep="\t", skiprows=1, header=None, names=[names[low], names[high]], dtype=str
    )
    return pairs, log


@app.function
def unichem_compound(source_name, compound_id, cache_dir, refresh=False):
    """UniChem API look-up of one compound's UCI, standard InChI and InChIKey, cached as JSON."""
    path = cache_dir / f"unichem_compound_{source_name}_{compound_id}.json"
    if refresh or not path.exists():
        response = http_session().post(
            "https://www.ebi.ac.uk/unichem/api/v1/compounds",
            json={
                "type": "sourceID",
                "sourceID": UNICHEM_SOURCE_IDS[source_name],
                "compound": compound_id,
            },
            timeout=60,
        )
        response.raise_for_status()
        compound = response.json()["compounds"][0]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {
                    "uci": compound["uci"],
                    "inchi": compound["inchi"]["inchi"],
                    "inchikey": compound["standardInchiKey"],
                    "retrieved_at": utc_now(),
                },
                indent=2,
            )
        )
    return json.loads(path.read_text())


@app.function
def fetch_unichem_gtopdb_batch(ligand_ids):
    """UniChem API: the UCI, standard InChI and InChIKey UniChem holds for each GtoPdb ligand ID, and the other GtoPdb
    ligands UniChem files under the same UCI. The API looks up one ID per request."""
    session = http_session()
    rows = []
    for ligand_id in ligand_ids:
        response = session.post(
            "https://www.ebi.ac.uk/unichem/api/v1/compounds",
            json={"type": "sourceID", "sourceID": UNICHEM_SOURCE_IDS["GtoPdb"], "compound": ligand_id},
            timeout=60,
        )
        response.raise_for_status()
        # an ID UniChem does not hold returns 200 with {"compounds": [], "response": "Not found"}
        compounds = response.json().get("compounds") or [{}]
        compound = compounds[0]
        rows.append(
            {
                "query_id": ligand_id,
                "found": bool(compound),
                "uci": compound.get("uci"),
                "inchi": (compound.get("inchi") or {}).get("inchi"),
                "inchikey": compound.get("standardInchiKey"),
                "other_gtopdb_ids_same_uci": "|".join(
                    sorted(
                        source["compoundId"]
                        for source in compound.get("sources", [])
                        if source["id"] == UNICHEM_SOURCE_IDS["GtoPdb"]
                        and source["compoundId"] != ligand_id
                    )
                ),
                "retrieved_at": utc_now(),
            }
        )
    return pd.DataFrame(rows)


@app.function
def ids_to_fetch(xref_db, gtp_xrefs):
    """IDs of one database to look up: GtoPdb's own usable cross-references, once each."""
    return sorted(set(gtp_xrefs.loc[gtp_xrefs.xref_db == xref_db, "xref_id_normalised"].dropna()))


@app.function
def fetch_cached(path, ids, fetch_batch, *, batch_size, workers=1, refresh=False, title=""):
    """Rows from `fetch_batch` for `ids`, fetching only those not cached at `path`. Progress is saved as it goes,
    so a run that fails resumes when rerun."""
    frames = [] if refresh or not path.exists() else [pd.read_parquet(path)]
    cached_ids = set(frames[0]["query_id"]) if frames else set()
    todo = sorted(set(ids) - cached_ids)
    batches = [todo[start : start + batch_size] for start in range(0, len(todo), batch_size)]
    if batches:
        path.parent.mkdir(parents=True, exist_ok=True)
        group = workers * 2  # batches between saves
        pool = ThreadPoolExecutor(max_workers=workers)
        try:
            with mo.status.progress_bar(
                total=len(batches), title=title, remove_on_exit=True
            ) as bar:
                for start in range(0, len(batches), group):
                    for frame in pool.map(fetch_batch, batches[start : start + group]):
                        frames.append(frame)
                        bar.update()
                    pd.concat(frames, ignore_index=True).to_parquet(path, index=False)
        finally:
            # on a failure, cancel the queued requests
            pool.shutdown(wait=False, cancel_futures=True)
    if not frames:
        return pd.DataFrame(
            {
                "query_id": pd.Series(dtype=str),
                "found": pd.Series(dtype=bool),
                "retrieved_at": pd.Series(dtype=str),
            }
        )
    table = pd.concat(frames, ignore_index=True)
    return table[table.query_id.isin(set(ids))].reset_index(drop=True)


@app.function
def retrieved_on(records):
    """Date of the latest record fetched into `records`, or "" when nothing was fetched."""
    return str(records.retrieved_at.max())[:10] if len(records) else ""


@app.function
def fetch_chembl_batch(chembl_ids):
    """ChEMBL REST API: standard InChI, preferred name and parent molecule for a batch of ChEMBL IDs."""
    response = http_session().get(
        "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
        params={
            "molecule_chembl_id__in": ",".join(chembl_ids),
            "only": "molecule_chembl_id,pref_name,molecule_structures,molecule_hierarchy",
            "limit": 1000,
        },
        timeout=30,
    )
    response.raise_for_status()
    found = {
        molecule["molecule_chembl_id"]: molecule for molecule in response.json()["molecules"]
    }
    retrieved = utc_now()
    rows = []
    for chembl_id in chembl_ids:
        molecule = found.get(chembl_id) or {}
        structure = molecule.get("molecule_structures") or {}
        hierarchy = molecule.get("molecule_hierarchy") or {}
        rows.append(
            {
                "query_id": chembl_id,
                "found": bool(molecule),
                "xref_name": molecule.get("pref_name"),
                "inchi": structure.get("standard_inchi"),
                "inchikey": structure.get("standard_inchi_key"),
                "chembl_parent_id": hierarchy.get("parent_chembl_id"),
                "retrieved_at": retrieved,
            }
        )
    return pd.DataFrame(rows)


@app.function
def fetch_chembl_id_lookup_batch(chembl_ids):
    """ChEMBL ID lookup: entity type and status of every ID ChEMBL has issued, including retired ones."""
    response = http_session().get(
        "https://www.ebi.ac.uk/chembl/api/data/chembl_id_lookup.json",
        params={"chembl_id__in": ",".join(chembl_ids), "limit": 1000},
        timeout=60,
    )
    response.raise_for_status()
    found = {item["chembl_id"]: item for item in response.json()["chembl_id_lookups"]}
    retrieved = utc_now()
    return pd.DataFrame(
        [
            {
                "query_id": chembl_id,
                "found": chembl_id in found,
                "entity_type": (found.get(chembl_id) or {}).get("entity_type"),
                "status": (found.get(chembl_id) or {}).get("status"),
                "last_active": (found.get(chembl_id) or {}).get("last_active"),
                "retrieved_at": retrieved,
            }
            for chembl_id in chembl_ids
        ]
    )


@app.function
def fetch_chebi_batch(chebi_ids):
    """ChEBI API: existence, primary/secondary status, name, stars and InChI for a batch of ChEBI IDs."""
    response = http_session().get(
        "https://www.ebi.ac.uk/chebi/backend/api/public/compounds/",
        params={"chebi_ids": ",".join(chebi_ids)},
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()
    retrieved = utc_now()
    rows = []
    for chebi_id in chebi_ids:
        entry = payload.get(chebi_id) or {}
        data = entry.get("data") or {}
        structure = data.get("default_structure") or {}
        rows.append(
            {
                "query_id": chebi_id,
                "found": bool(entry.get("exists")),
                "xref_name": data.get("ascii_name") or data.get("name"),
                "chebi_id_type": entry.get("id_type"),
                "chebi_primary_id": entry.get("primary_chebi_id"),
                "chebi_stars": data.get("stars"),
                # empty for entries without a defined structure (classes, R groups)
                "inchi": structure.get("standard_inchi") or None,
                "inchikey": structure.get("standard_inchi_key") or None,
                "retrieved_at": retrieved,
            }
        )
    return pd.DataFrame(rows)


@app.function
def fetch_pubchem_batch(cids):
    """PubChem PUG-REST: InChI, InChIKey and title for a batch of CIDs."""
    response = http_session().post(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/InChI,InChIKey,Title/JSON",
        data={"cid": ",".join(cids)},
        timeout=60,
    )
    response.raise_for_status()
    found = {
        str(record["CID"]): record for record in response.json()["PropertyTable"]["Properties"]
    }
    time.sleep(0.25)  # PubChem asks for no more than 5 requests per second
    retrieved = utc_now()
    rows = []
    for cid in cids:
        record = found.get(cid) or {}
        rows.append(
            {
                "query_id": cid,
                "found": bool(record),  # even without an InChI (very large structures)
                "xref_name": record.get("Title"),
                "inchi": record.get("InChI"),
                "inchikey": record.get("InChIKey"),
                "retrieved_at": retrieved,
            }
        )
    return pd.DataFrame(rows)


@app.function
def fetch_pubchem_sid_cids(sids):
    """PubChem PUG-REST: the compound(s) (CIDs) each deposited substance (SID) is merged into."""
    response = http_session().post(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/cids/JSON",
        data={"sid": ",".join(sids)},
        timeout=120,
    )
    response.raise_for_status()
    found = {
        str(entry["SID"]): [str(cid) for cid in entry.get("CID", [])]
        for entry in response.json()["InformationList"]["Information"]
    }
    time.sleep(0.25)  # PubChem asks for no more than 5 requests per second
    retrieved = utc_now()
    return pd.DataFrame(
        [
            {
                "query_id": sid,
                "found": bool(found.get(sid)),
                "cids": "|".join(found.get(sid, [])),
                "retrieved_at": retrieved,
            }
            for sid in sids
        ]
    )


@app.function
def fetch_pubchem_by_inchikey(inchikeys):
    """PubChem PUG-REST: CID, InChI and title of every compound with each InChIKey (one request per key)."""
    session = http_session()
    rows = []
    for inchikey in inchikeys:
        response = session.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/InChI,Title/JSON",
            timeout=60,
        )
        hits = []
        if response.status_code != 404:  # no compound has this InChIKey
            response.raise_for_status()
            hits = [
                {
                    "xref_id": str(record["CID"]),
                    "inchi": record.get("InChI"),
                    "xref_name": record.get("Title"),
                }
                for record in response.json()["PropertyTable"]["Properties"]
            ]
        rows += inchikey_hits(inchikey, hits, utc_now())
        time.sleep(0.25)  # PubChem asks for no more than 5 requests per second
    return pd.DataFrame(rows)


@app.function
def fetch_pubchem_standardised(inchis):
    """PubChem's standardisation of each InChI, as a PubChem compound would hold it (one request per InChI)."""
    session = http_session()
    rows = []
    for inchi in inchis:
        response = session.post(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/standardize/inchi/JSON",
            data={"inchi": inchi},
            timeout=60,
        )
        standardised = None
        if response.status_code != 400:  # 400: PubChem cannot standardise this structure
            response.raise_for_status()
            props = response.json()["PC_Compounds"][0]["props"]
            standardised = next(
                (
                    prop["value"]["sval"]
                    for prop in props
                    if prop["urn"].get("label") == "InChI"
                ),
                None,
            )
        rows.append(
            {
                "query_id": inchi,
                "found": standardised is not None,
                "standardised_inchi": standardised,
                "retrieved_at": utc_now(),
            }
        )
        time.sleep(0.25)  # PubChem asks for no more than 5 requests per second
    return pd.DataFrame(rows)


@app.function
def fetch_chembl_by_inchikey(inchikeys):
    """ChEMBL REST API: ChEMBL ID and standard InChI of every molecule with each standard InChIKey."""
    response = http_session().get(
        "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
        params={
            "molecule_structures__standard_inchi_key__in": ",".join(inchikeys),
            "only": "molecule_chembl_id,molecule_structures",
            "limit": 1000,
        },
        timeout=60,
    )
    response.raise_for_status()
    hits = {}
    for molecule in response.json()["molecules"]:
        structure = molecule.get("molecule_structures") or {}
        hits.setdefault(structure.get("standard_inchi_key"), []).append(
            {
                "xref_id": molecule["molecule_chembl_id"],
                "inchi": structure.get("standard_inchi"),
            }
        )
    retrieved = utc_now()
    return pd.DataFrame(
        [
            row
            for inchikey in inchikeys
            for row in inchikey_hits(inchikey, hits.get(inchikey, []), retrieved)
        ]
    )


@app.function
def fetch_chebi_by_inchikey(inchikeys):
    """ChEBI API search: accession, InChI and name of every entry with each InChIKey (one request per key)."""
    session = http_session()
    rows = []
    for inchikey in inchikeys:
        response = session.get(
            "https://www.ebi.ac.uk/chebi/backend/api/public/es_search/",
            params={"term": inchikey},
            timeout=60,
        )
        hits = []
        if response.status_code != 404:  # no entry matches the term
            response.raise_for_status()
            hits = [
                {
                    "xref_id": source["chebi_accession"],
                    "inchi": source.get("inchi"),
                    "xref_name": source.get("ascii_name") or source.get("name"),
                }
                for source in (entry["_source"] for entry in response.json().get("results", []))
                if source.get("inchikey") == inchikey  # the search is full-text
            ]
        rows += inchikey_hits(inchikey, hits, utc_now())
    return pd.DataFrame(rows)


@app.function
def inchikey_hits(inchikey, hits, retrieved):
    """Rows of an InChIKey look-up: one per hit, or one empty row if there are none."""
    return [
        {"query_id": inchikey, "found": True, **hit, "retrieved_at": retrieved} for hit in hits
    ] or [
        {
            "query_id": inchikey,
            "found": False,
            "xref_id": None,
            "inchi": None,
            "retrieved_at": retrieved,
        }
    ]


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### Identifiers

    ChEMBL, ChEBI and PubChem identifiers: normalising them and finding what is wrong with them. Also GtoPdb's HTML names.
    """)
    return


@app.class_definition
class IdentifierIssue(StrEnum):
    """Kinds of identifier issue, as found by `identifier_issues`."""

    # how the ID is written in GtoPdb
    INVISIBLE_CHARACTER = "invisible_character"
    WHITESPACE = "whitespace"
    PREFIX_CASE = "prefix_case"
    NOT_VALID = "not_valid"
    SEVERAL_IDS_IN_FIELD = "several_ids_in_field"
    # what the linked database holds under the ID
    CHEMBL_NEVER_ISSUED = "chembl_never_issued"
    CHEMBL_NOT_A_COMPOUND = "chembl_not_a_compound"
    CHEMBL_RETIRED = "chembl_retired"
    CHEMBL_NOT_RETURNED = "chembl_not_returned"
    CHEBI_NOT_FOUND = "chebi_not_found"
    CHEBI_SECONDARY_ID = "chebi_secondary_id"
    PUBCHEM_NO_STRUCTURE = "pubchem_no_structure"
    # found by comparing the structures
    CHEMBL_CHILD_FORM = "chembl_child_form"


@app.function
def plain_text(name_html):
    """A GtoPdb ligand name without its HTML markup (<sup>, &alpha;)."""
    return html.unescape(re.sub(r"<[^>]+>", "", name_html))


@app.function
def without_invisible_characters(text):
    """`text` without format characters (Unicode category Cf), such as zero-width spaces and byte-order marks."""
    return "".join(char for char in text if unicodedata.category(char) != "Cf")


@app.function
def normalise_xref(db, raw_id):
    """The ID as `db` writes it (cleaned and upper-cased), or None if it is not a valid `db` identifier."""
    cleaned = without_invisible_characters(raw_id).strip().upper()
    # PubChem reads CIDs as signed 32-bit integers: a larger one fails the whole batch
    valid = re.fullmatch(XREF_PATTERNS[db], cleaned) is not None and (
        db != "PubChem" or int(cleaned) < 2**31
    )
    return cleaned if valid else None


@app.function
def identifier_issues(row, first_difference):
    """{`IdentifierIssue`: description} for one cross-reference (a row of `xref_table`), in the order reported.

    Checks how the ID is written, then what the linked database holds under it, then whether a ChEMBL ID is a salt
    or child form of a structure that differs (`first_difference`, from `classify_inchi_pair`)."""
    db, raw_id = row["xref_db"], row["xref_id_as_written"]
    issues = {}
    hidden = sorted(
        {f"U+{ord(char):04X}" for char in raw_id if unicodedata.category(char) == "Cf"}
    )
    if hidden:
        issues[IdentifierIssue.INVISIBLE_CHARACTER] = f"invisible character {', '.join(hidden)}"
    if raw_id != raw_id.strip():
        issues[IdentifierIssue.WHITESPACE] = "leading/trailing whitespace"
    prefix = re.match(r"[A-Za-z]+", without_invisible_characters(raw_id).strip())
    if prefix and prefix.group() != prefix.group().upper():
        issues[IdentifierIssue.PREFIX_CASE] = f"non-standard prefix case '{prefix.group()}'"
    xref_id = normalise_xref(db, raw_id)
    if xref_id is None:
        issues[IdentifierIssue.NOT_VALID] = f"not a valid {db} identifier"
    if row["ids_in_field"] > 1:
        issues[IdentifierIssue.SEVERAL_IDS_IN_FIELD] = (
            f"{row['ids_in_field']} IDs in one {db} field"
        )
    if xref_id is None:  # unreadable IDs are not looked up
        return issues

    if not row["xref_found"]:
        if db == "ChEBI":
            issues[IdentifierIssue.CHEBI_NOT_FOUND] = (
                "does not exist in ChEBI (deleted, or never issued)"
            )
        elif db == "PubChem":
            issues[IdentifierIssue.PUBCHEM_NO_STRUCTURE] = (
                "PubChem returns no structure for this CID"
            )
        # the ID lookup, unlike the molecule API, knows retired IDs and those of targets, assays, documents
        elif not row["chembl_id_issued"]:
            issues[IdentifierIssue.CHEMBL_NEVER_ISSUED] = "never issued by ChEMBL (likely a typo)"
        elif row["chembl_entity_type"] != "COMPOUND":
            issues[IdentifierIssue.CHEMBL_NOT_A_COMPOUND] = (
                f"ChEMBL {row['chembl_entity_type'].lower()} ID, not a compound"
            )
        elif row["chembl_status"] != "ACTIVE":
            last_active = row["chembl_last_active"]
            last = f" (last in ChEMBL {int(last_active)})" if pd.notna(last_active) else ""
            issues[IdentifierIssue.CHEMBL_RETIRED] = f"retired ChEMBL compound{last}"
        else:
            issues[IdentifierIssue.CHEMBL_NOT_RETURNED] = (
                "active ChEMBL compound that the molecule API does not return"
            )
    if db == "ChEBI" and row.get("chebi_id_type") == "SECONDARY_ID":
        issues[IdentifierIssue.CHEBI_SECONDARY_ID] = (
            f"secondary ChEBI ID (primary {row['chebi_primary_id']})"
        )

    # a child form is the right link if it has the GtoPdb structure, so the parent is suggested only if they differ
    parent = row.get("chembl_parent_id")
    if (
        db == "ChEMBL"
        and isinstance(parent, str)
        and parent != xref_id
        and first_difference not in ("identical", "not compared")
    ):
        issues[IdentifierIssue.CHEMBL_CHILD_FORM] = f"ChEMBL salt/child form (parent {parent})"
    return issues


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### InChI layers

    InChI layers and stereo, read from standard InChI strings (no RDKit).
    """)
    return


@app.function
def inchi_parts(inchi):
    """The layers of a standard InChI as written after 'InChI=1S/', as (key, text) pairs: the formula, then each
    layer keyed by its letter, with the isotopic sublayers as 'ih', 'ib', 'it', 'im' and 'is'. Raises ValueError for
    a non-standard InChI, whose /f and /r layers repeat the others."""
    if not inchi.startswith("InChI=1S/"):
        raise ValueError(f"not a standard InChI: {inchi}")
    parts = inchi.split("/")
    pairs = [("formula", parts[1])]
    section = ""  # "" main layers, "i" isotopic sublayers
    for part in parts[2:]:
        if part[:1] == "i":
            section = "i"
            pairs.append(("i", part))
        else:
            pairs.append((section + part[0], part))
    return pairs


@app.function
def parse_standard_inchi(inchi):
    """Split a standard InChI into {formula, c, h, q, p, b, t, m, s, i}, with the isotopic sublayers as 'ih', 'ib',
    'it', 'im' and 'is' (`inchi_parts`), each without its letter. Raises ValueError for a non-standard InChI."""
    return {key: text if key == "formula" else text[1:] for key, text in inchi_parts(inchi)}


@app.function
def share_skeleton(layers_a, layers_b):
    """Whether two parsed InChIs have the same atoms connected the same way, hydrogens aside (so numbered alike)."""
    formulas = {
        re.sub(r"H(?![a-z])\d*", "", layers["formula"]) for layers in (layers_a, layers_b)
    }
    return layers_a.get("c") == layers_b.get("c") and len(formulas) == 1


@app.function
def strip_layers(inchi, drop):
    """A standard InChI without the layers in `drop`, from the isotopic layer too ('i' drops all of it)."""
    layers = parse_standard_inchi(inchi)
    kept = ["InChI=1S", layers["formula"]] + [
        key + layers[key]
        for key in ("c", "h", "q", "p", "b", "t", "m", "s")
        if key in layers and key not in drop
    ]
    if "i" in layers and "i" not in drop:
        kept += ["i" + layers["i"]] + [
            key + layers["i" + key]
            for key in ("h", "b", "t", "m", "s")
            if "i" + key in layers and key not in drop
        ]
    return "/".join(kept)


@app.function
def stereo_layers_of(inchi):
    """The stereo layers (/b, /t, /m, /s) of an InChI as written, or "none"."""
    layers = parse_standard_inchi(inchi)
    return " ".join(f"/{key}{layers[key]}" for key in STEREO_LAYERS if key in layers) or "none"


@app.function
def layer_components(layer):
    """Per-component descriptors of a /b or /t layer, with 'n*' repeats expanded ('2*2-' -> ['2-', '2-'])."""
    components = []
    for part in layer.split(";"):
        repeat = re.match(r"(\d+)\*", part)
        components += [part[repeat.end() :] if repeat else part] * (
            int(repeat.group(1)) if repeat else 1
        )
    return components


@app.function
def stereo_items(layers, key):
    """(component, centre or bond, parity) for each item of a parsed InChI's /t or /b layer, e.g. (0, '3', '-') for
    '/t3-'."""
    for component, descriptors in enumerate(layer_components(layers.get(key, ""))):
        for item in filter(None, descriptors.split(",")):
            yield component, item[:-1], item[-1]


@app.function
def defined_centres(layers):
    """{(component, atom): parity} for defined tetrahedral centres, with the /m inversion applied. /m has one
    character per component: '0' as drawn, '1' inverted, '.' no flag (e.g. '/m.0')."""
    flip = {"+": "-", "-": "+"}
    inverted = layers.get("m", "")
    return {
        (component, atom): flip[parity] if inverted[component : component + 1] == "1" else parity
        for component, atom, parity in stereo_items(layers, "t")
        if parity in "+-"
    }


@app.function
def defined_double_bonds(layers):
    """{(component, bond): parity} for double bonds with defined E/Z geometry."""
    return {(component, bond): parity for component, bond, parity in stereo_items(layers, "b") if parity in "+-"}


@app.function
def stereo_listing(layers):
    """{(layer, component, centre or bond): whether defined} for every centre ('t') and double bond ('b') an InChI
    lists. InChI lists a component's undefined ('?') items only if it defines another of the same kind there."""
    return {
        (key, component, element): parity != "?"
        for key in ("t", "b")
        for component, element, parity in stereo_items(layers, key)
    }


@app.function
def stereo_definitions(gtp_layers, xref_layers):
    """How many stereocentres and double bonds each record defines, out of those either lists, e.g.
    'stereocentres defined: GtoPdb 3 of 3, linked record 2 of 3'."""
    listings = [stereo_listing(layers) for layers in (gtp_layers, xref_layers)]
    counts = []
    for kind, key in (("stereocentres", "t"), ("double bonds", "b")):
        listed = {element for listing in listings for element in listing if element[0] == key}
        if listed:
            gtp_count, xref_count = (
                sum(defined for (layer, *_), defined in listing.items() if layer == key) for listing in listings
            )
            counts.append(
                f"{kind} defined: GtoPdb {gtp_count} of {len(listed)}, linked record {xref_count} of {len(listed)}"
            )
    return "; ".join(counts)


@app.function
def formula_components(formula):
    """Number of components in an InChI formula layer: 2 for '2C3H7NO2', 3 for '2C20H24N2O2.H2O4S'."""
    return sum(int(re.match(r"\d*", part).group() or 1) for part in formula.split("."))


@app.function
def stereo_difference(gtp_layers, xref_layers):
    """How two InChIs of the same constitution differ in stereo: (`stereo_relation`, detail).

    Parities are compared only when both records list and define the same centres and double bonds, because next to
    an undefined centre InChI can write one configuration with either parity. Then a differing double bond makes an
    E/Z isomer, some centres inverted an epimer or diastereomer, and all inverted an enantiomer ("enantiomer or
    diastereomer" if both leave some stereo undefined). If no parity differs, the layers were misread: ValueError.

    Otherwise a record defines more stereo if it defines all the other does, written the same way; any other
    difference is `STEREO_DEFINED_DIFFERENTLY`. The detail counts what each defines (`stereo_definitions`).
    """
    gtp_centres, xref_centres = defined_centres(gtp_layers), defined_centres(xref_layers)
    gtp_bonds, xref_bonds = defined_double_bonds(gtp_layers), defined_double_bonds(xref_layers)
    if not (xref_centres or xref_bonds):
        return "linked record has no stereo", "linked record has no stereo"
    if not (gtp_centres or gtp_bonds):
        return "GtoPdb structure has no stereo", "GtoPdb structure has no stereo"
    gtp_listing, xref_listing = stereo_listing(gtp_layers), stereo_listing(xref_layers)
    if gtp_listing == xref_listing:
        opposite = [key for key in gtp_centres if gtp_centres[key] != xref_centres[key]]
        if gtp_bonds != xref_bonds:
            return "E/Z isomer", "E/Z isomer"
        if opposite and len(opposite) == len(gtp_centres):
            # same listings, so both leave the same stereo undefined
            relation = "enantiomer" if all(gtp_listing.values()) else "enantiomer or diastereomer"
            return relation, relation
        if opposite:
            return "epimer or diastereomer", "epimer or diastereomer"
    else:
        gtp_defined = {("centre", *key): parity for key, parity in gtp_centres.items()} | {
            ("bond", *key): parity for key, parity in gtp_bonds.items()
        }
        xref_defined = {("centre", *key): parity for key, parity in xref_centres.items()} | {
            ("bond", *key): parity for key, parity in xref_bonds.items()
        }
        shared = gtp_defined.keys() & xref_defined.keys()
        inverted = sum(gtp_defined[key] != xref_defined[key] for key in shared)
        if not inverted and gtp_defined.keys() > xref_defined.keys():
            relation = "GtoPdb defines more stereo"
        elif not inverted and xref_defined.keys() > gtp_defined.keys():
            relation = "linked record defines more stereo"
        else:
            relation = STEREO_DEFINED_DIFFERENTLY
        counts = stereo_definitions(gtp_layers, xref_layers)
        if inverted:
            counts += f"; {inverted} of the {len(shared)} defined on both written inverted"
        return relation, f"{relation} ({counts})"
    written = [
        " ".join(f"/{key}{layers[key]}" for key in STEREO_LAYERS if key in layers)
        for layers in (gtp_layers, xref_layers)
    ]
    raise ValueError(
        "stereo layers differ but no defined centre or double bond does: "
        f"GtoPdb {written[0]}, linked record {written[1]}"
    )


@app.function
def isotope_difference(gtp_layers, xref_layers):
    """Which record carries isotope labels, for two InChIs whose /i layers differ."""
    if "i" in gtp_layers and "i" not in xref_layers:
        return "GtoPdb labelled, linked record unlabelled"
    if "i" in xref_layers and "i" not in gtp_layers:
        return "linked record labelled, GtoPdb unlabelled"
    return "different labelling"


@app.function
def isotope_layer(layers):
    """A parsed InChI's /i layer with its isotopic /h sublayer, or "none"."""
    if "i" not in layers:
        return "none"
    return f"/i{layers['i']}" + (f"/h{layers['ih']}" if "ih" in layers else "")


@app.function
def labelled_stereo(layers):
    """A parsed InChI's stereo layers including stereo that isotope labels create: each isotopic sublayer where
    given (InChI writes one only where it differs), otherwise the main layer."""
    return {
        key: layers.get("i" + key, layers.get(key))
        for key in STEREO_LAYERS
        if "i" + key in layers or key in layers
    }


@app.function
def comparable_stereo(layers_a, layers_b):
    """The stereo layers to compare for two parsed InChIs of one constitution: `labelled_stereo` when both carry
    the same labels, otherwise the main layers (a labelled /m read against an unlabelled meso /t inverts every
    centre)."""
    if isotope_layer(layers_a) == isotope_layer(layers_b):
        return labelled_stereo(layers_a), labelled_stereo(layers_b)
    return tuple(
        {key: layers[key] for key in STEREO_LAYERS if key in layers} for layers in (layers_a, layers_b)
    )


@app.function
def compare_layers(layers_a, layers_b):
    """(shared, first differing) layers of two parsed InChIs, read in `LAYERS` order: the layers before the first
    difference that either record has, and the first layer that differs (None if none does). Stereo includes the
    stereo that isotope labels create when both carry the same labels (`comparable_stereo`)."""
    stereo_a, stereo_b = comparable_stereo(layers_a, layers_b)
    values = {
        "formula": (layers_a["formula"], layers_b["formula"]),
        "connections": (layers_a.get("c"), layers_b.get("c")),
        "hydrogens": (layers_a.get("h"), layers_b.get("h")),
        "charge": (layers_a.get("q"), layers_b.get("q")),
        "protons": (layers_a.get("p"), layers_b.get("p")),
        "stereo": (stereo_a or None, stereo_b or None),
        "isotopes": tuple(isotope_layer(layers) if "i" in layers else None for layers in (layers_a, layers_b)),
    }
    shared = []
    for layer in LAYERS:
        value_a, value_b = values[layer]
        if value_a != value_b:
            return shared, layer
        if value_a is not None:
            shared.append(layer)
    return shared, None


@app.function
def first_difference_keys(layers_a, layers_b):
    """The keys (as `parse_standard_inchi` names them) of the first layer that differs (`compare_layers`) whose values
    differ between two parsed InChIs, or an empty set if none does. Stereo includes the isotopic stereo sublayers
    when both carry the same labels, as `comparable_stereo` compares them."""
    _, first = compare_layers(layers_a, layers_b)
    keys = {
        None: (),
        "formula": ("formula",),
        "connections": ("c",),
        "hydrogens": ("h",),
        "charge": ("q",),
        "protons": ("p",),
        "stereo": STEREO_LAYERS
        + (tuple("i" + key for key in STEREO_LAYERS) if isotope_layer(layers_a) == isotope_layer(layers_b) else ()),
        "isotopes": ("i", "ih"),  # as `isotope_layer` writes them
    }[first]
    return {key for key in keys if layers_a.get(key) != layers_b.get(key)}


@app.function
def net_charge(layers):
    """Net charge of a parsed InChI: its /q charges plus the /p protons."""
    charge = 0
    for component in layers.get("q", "").split(";"):
        count, _, value = component.rpartition("*")
        if value:
            charge += int(count or 1) * int(value)
    return charge + int(layers.get("p") or 0)


@app.function
def charge_difference(gtp_layers, xref_layers):
    """Each record's net charge, and how many more or fewer protons GtoPdb has when the formulas are equal."""

    def signed(charge):
        return f"{charge:+d}" if charge else "0"

    gtp_charge, xref_charge = net_charge(gtp_layers), net_charge(xref_layers)
    details = []
    if gtp_charge != xref_charge:
        details.append(
            f"net charge GtoPdb {signed(gtp_charge)}, linked record {signed(xref_charge)}"
        )
    elif gtp_layers.get("q") != xref_layers.get("q"):
        details.append("same net charge, different component charge layers")
    protons = int(gtp_layers.get("p") or 0) - int(xref_layers.get("p") or 0)
    if protons and gtp_layers["formula"] == xref_layers["formula"]:
        details.append(
            f"GtoPdb has {abs(protons)} {'proton' if abs(protons) == 1 else 'protons'} "
            f"{'more' if protons > 0 else 'fewer'}"
        )
    return "; ".join(details)


@app.function
def stereo_difference_keys(layers_a, layers_b):
    """(centres, double bonds) whose stereo differs between two InChIs of one constitution: opposite parity, or
    defined on one side only."""
    centres_a, centres_b = defined_centres(layers_a), defined_centres(layers_b)
    bonds_a, bonds_b = defined_double_bonds(layers_a), defined_double_bonds(layers_b)
    return (
        {
            key
            for key in centres_a.keys() | centres_b.keys()
            if centres_a.get(key) != centres_b.get(key)
        },
        {
            key
            for key in bonds_a.keys() | bonds_b.keys()
            if bonds_a.get(key) != bonds_b.get(key)
        },
    )


@app.function
def fixed_hydrogens(layers):
    """{(component, atom): hydrogen count} from a parsed InChI's /h layer, without the mobile hydrogens (in
    parentheses), whose position is not defined."""
    counts = {}
    for component, part in enumerate(layer_components(layers.get("h", ""))):
        for atoms, n in re.findall(r"([\d,-]+)H(\d*)", re.sub(r"\([^)]*\)", "", part)):
            for item in filter(None, atoms.split(",")):
                start, _, end = item.partition("-")
                for atom in range(int(start), int(end or start) + 1):
                    counts[component, atom] = int(n or 1)
    return counts


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### RDKit helpers

    RDKit structure helpers.
    """)
    return


@app.function
def mol_from_inchi(inchi):
    """RDKit molecule from an InChI, or None when the InChI is missing or unparsable."""
    return Chem.MolFromInchi(inchi) if isinstance(inchi, str) and inchi else None


@app.function
def heavy_and_hydrogen_count(smiles):
    """Number of atoms including hydrogens in a SMILES, or None when it cannot be parsed."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    return Chem.AddHs(mol).GetNumAtoms() if mol is not None else None


@app.function
def inchi_matching_key(smiles, inchikey):
    """Standard InChI built from `smiles`, kept only if it matches `inchikey`."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    inchi = Chem.MolToInchi(mol) if mol is not None else None
    return inchi if inchi and rdinchi.InchiToInchiKey(inchi) == inchikey else None


@app.function
def inorganic_smiles(name):
    """SMILES of an inorganic ligand named by its formula (a monatomic ion, or one listed below), otherwise None."""
    molecules = {
        "ClO4-": "[O-]Cl(=O)(=O)=O",
        "CO": "[C-]#[O+]",
        "H2O2": "OO",
        "H2S": "S",
        "HCl": "Cl",
        "NaHCO3": "[Na+].OC([O-])=O",
        "NH4Cl": "[NH4+].[Cl-]",
        "NO": "[N]=O",
        "O3": "[O-][O+]=O",
        "SCN-": "[S-]C#N",
    }
    ion = re.fullmatch(r"([A-Z][a-z]?)(\d?)([+-])", name)
    return molecules.get(name) or (f"[{ion[1]}{ion[3]}{ion[2]}]" if ion else None)


@app.function
def has_metal(inchi):
    """Whether a standard InChI is a metal complex (has an element outside SALT_AND_NONMETAL_ELEMENTS), which
    standard InChI writes as its ligands and disconnected metal ions."""
    elements = set(re.findall(r"[A-Z][a-z]?", parse_standard_inchi(inchi)["formula"]))
    return not elements <= SALT_AND_NONMETAL_ELEMENTS


@app.function
def parent_inchi(inchi, canonical_tautomer=False, keep_stereo=False):
    """InChI of the neutralised, unlabelled largest organic fragment, without stereo unless `keep_stereo`,
    optionally as the canonical tautomer. None for a metal complex, whose largest fragment would be a ligand."""
    mol = mol_from_inchi(inchi)
    if mol is None or has_metal(inchi):
        return None
    mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    if not keep_stereo:
        Chem.RemoveStereochemistry(mol)
    if canonical_tautomer:
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(1000)
        mol = enumerator.Canonicalize(mol)
    return Chem.MolToInchi(mol) or None


@app.function
def missing_stereo(more_inchi, less_inchi):
    """Kinds of stereo `more_inchi` defines and `less_inchi` (same constitution) does not: 'R/S centre', 'E/Z' or
    'ring cis/trans' (a centre with two symmetry-equivalent neighbours)."""
    more, less = comparable_stereo(parse_standard_inchi(more_inchi), parse_standard_inchi(less_inchi))
    centres = defined_centres(more).keys() - defined_centres(less).keys()
    kinds = (
        {"E/Z"}
        if defined_double_bonds(more).keys() - defined_double_bonds(less).keys()
        else set()
    )
    mol = mol_from_inchi(more_inchi)
    if centres and mol is not None:
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False))
        components = Chem.GetMolFrags(mol)  # RDKit keeps the InChI component and atom order
        for component, atom_number in centres:
            atom = mol.GetAtomWithIdx(components[component][int(atom_number) - 1])
            neighbour_ranks = [ranks[neighbour.GetIdx()] for neighbour in atom.GetNeighbors()]
            symmetric = len(set(neighbour_ranks)) < len(neighbour_ranks)
            kinds.add("ring cis/trans" if symmetric else "R/S centre")
    return " + ".join(sorted(kinds))


@app.function
def morgan_similarity(inchi_a, inchi_b):
    """Morgan radius-2 Tanimoto similarity of two InChIs, or None if either cannot be parsed."""
    mol_a, mol_b = mol_from_inchi(inchi_a), mol_from_inchi(inchi_b)
    if mol_a is None or mol_b is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return round(
        DataStructs.TanimotoSimilarity(
            generator.GetFingerprint(mol_a), generator.GetFingerprint(mol_b)
        ),
        3,
    )


@app.function
def radical_count(smiles_or_inchi):
    """Number of radical electrons in a SMILES or InChI (0 if it cannot be parsed)."""
    mol = (
        mol_from_inchi(smiles_or_inchi)
        if smiles_or_inchi.startswith("InChI=")
        else Chem.MolFromSmiles(smiles_or_inchi)
    )
    return (
        sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()) if mol is not None else 0
    )


@app.function
def isotopes_written(smiles, name_html):
    """Isotopes in a ligand's SMILES or name, such as ["125I", "3H"]. In a name, an isotope is a superscript mass
    before an element symbol (<sup>125</sup>I); other superscripts are sequence positions (Thr<sup>4</sup>)."""
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    found = (
        {
            f"{atom.GetIsotope()}{atom.GetSymbol()}"
            for atom in mol.GetAtoms()
            if atom.GetIsotope()
        }
        if mol is not None
        else set()
    )
    found |= {
        mass + symbol
        for mass, symbol in re.findall(
            r"(?:^|[\[,(\s\d-])<sup>(\d+)</sup>([A-Z][a-z]?)(?![a-z])", name_html
        )
    }
    return sorted(found)


@app.function
def is_radioisotope(isotope):
    """Ground-state radioactivity for recognised nuclides, or None for an unrecognised label."""
    mass, symbol = re.fullmatch(r"(\d+)([A-Z][a-z]?)", isotope).groups()
    table = Chem.GetPeriodicTable()
    try:
        atomic_number = table.GetAtomicNumber(symbol)
    except RuntimeError:  # not an element symbol
        return None
    if not table.GetMassForIsotope(atomic_number, int(mass)):
        return None
    return isotope not in STABLE_NUCLIDES


@app.function
def ring_systems(mol):
    """Atom index sets of the fused ring systems of `mol`."""
    systems = []
    for ring in mol.GetRingInfo().AtomRings():
        merged = set(ring)
        for system in [system for system in systems if system & merged]:
            merged |= system
            systems.remove(system)
        systems.append(merged)
    return systems


@app.function
def fragment_smiles(mol, atoms):
    """SMILES of the given atoms of `mol`, with * at each bond to the rest of the molecule."""
    cut = [
        bond.GetIdx()
        for bond in mol.GetBonds()
        if (bond.GetBeginAtomIdx() in atoms) != (bond.GetEndAtomIdx() in atoms)
    ]
    pieces = (
        Chem.FragmentOnBonds(mol, cut, dummyLabels=[(0, 0)] * len(cut))
        if cut
        else Chem.Mol(mol)
    )
    return ".".join(
        sorted(
            Chem.MolFragmentToSmiles(pieces, fragment, isomericSmiles=False)
            for fragment in Chem.GetMolFrags(pieces)
            if atoms & set(fragment)
        )
    )


@app.function
def common_substructure(inchi_a, inchi_b, max_atoms=60):
    """(molecules, SMARTS pattern of their largest common substructure), or None. Structures above `max_atoms`
    heavy atoms are too slow to search and give None."""
    mols = [mol_from_inchi(inchi_a), mol_from_inchi(inchi_b)]
    if None in mols or max(mol.GetNumAtoms() for mol in mols) > max_atoms:
        return None
    mcs = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, completeRingsOnly=True, timeout=10)
    if not mcs.numAtoms or mcs.canceled:
        return None
    return mols, Chem.MolFromSmarts(mcs.smartsString)


@app.function
def connectivity_difference(inchi_a, inchi_b):
    """The parts of two same-formula structures outside their largest common substructure, as SMILES with * at the
    attachment points, or None. Parts that look alike are widened to their ring system to show where they attach."""
    found = common_substructure(inchi_a, inchi_b)
    if found is None:
        return None
    mols, pattern = found
    differing = [
        {atom.GetIdx() for atom in mol.GetAtoms()} - set(mol.GetSubstructMatch(pattern))
        for mol in mols
    ]
    parts = [fragment_smiles(mol, atoms) for mol, atoms in zip(mols, differing)]
    if parts[0] == parts[1]:
        differing = [
            atoms.union(
                *[
                    system
                    for system in ring_systems(mol)
                    if any(mol.GetBondBetweenAtoms(i, j) for i in system for j in atoms)
                ]
            )
            for mol, atoms in zip(mols, differing)
        ]
        parts = [fragment_smiles(mol, atoms) for mol, atoms in zip(mols, differing)]
    if parts[0] == parts[1] or not all(parts):
        return None
    return f"GtoPdb {parts[0]}; linked record {parts[1]}"


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### Structure drawings

    Structure drawings: single molecules, and pairs laid out alike with their differences highlighted.
    """)
    return


@app.function
def radical_figure(smiles, legend, scale=15, max_width=600):
    """(HTML, width) of a SMILES drawing with radical atoms in orange, scaled down to fit `max_width`."""
    mol = Chem.MolFromSmiles(smiles)
    radicals = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons()}
    orange = dict.fromkeys(radicals, (1.0, 0.85, 0.5))
    caption = html.escape(legend)
    figure, width, _ = molecule_figure(mol, caption, scale, radicals, atom_colours=orange)
    if width > max_width:
        figure, width, _ = molecule_figure(
            mol, caption, scale * max_width / width, radicals, atom_colours=orange
        )
    return figure, width


@app.function
def drawable_mol(inchi, missing="no structure"):
    """(molecule, note) for drawing an InChI. The note says why there is no molecule, or that a valence error was
    ignored so that the structure can still be drawn."""
    if not isinstance(inchi, str) or not inchi:
        return None, missing
    mol = Chem.MolFromInchi(inchi)
    if mol is not None:
        return mol, ""
    mol = Chem.MolFromInchi(inchi, sanitize=False)
    if mol is None:
        return None, "InChI cannot be parsed by RDKit"
    problems = "; ".join(problem.Message() for problem in Chem.DetectChemistryProblems(mol))
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
        )
    except ValueError as error:  # RDKit's sanitisation errors derive from ValueError
        return None, f"InChI cannot be drawn by RDKit: {problems or error}"
    mol.UpdatePropertyCache(strict=False)
    return mol, f"drawn without the valence check: {problems}"


@app.function
def molecule_figure(
    mol,
    legend,
    scale=20,
    highlight_atoms=(),
    highlight_bonds=(),
    atom_colours=None,
    bond_colours=None,
    stereo_labels=False,
    min_width=120,
    note="",
):
    """(HTML, width, height) of `mol` on a fitted canvas at `scale` pixels per ångström, captioned with `legend`, which
    is HTML (RDKit's own legend is unreadably small on a fitted canvas), and `note` in red below it. The colour maps
    give highlighted atoms and bonds an RGB colour. `width` is the drawing's; a longer caption wraps at the width of
    the element holding the figure."""
    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    options = drawer.drawOptions()
    options.scalingFactor = scale
    if stereo_labels:
        options.addStereoAnnotation = True
        options.annotationFontScale = 0.7  # R/S and E/Z labels readable at table-cell size
    drawer.DrawMolecule(
        mol,
        highlightAtoms=sorted(highlight_atoms),
        highlightBonds=sorted(highlight_bonds),
        highlightAtomColors=atom_colours or {},
        highlightBondColors=bond_colours or {},
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    width = max(drawer.Width(), min_width)
    if note:
        legend += f'<div style="color: #8b1a1a">{html.escape(note)}</div>'
    figure = (
        f'<figure style="margin: 0; min-width: {width}px; display: flex; flex-direction: column; align-items: center">'
        f'<div style="flex: 1; display: flex; align-items: center">{svg[svg.index("<svg") :]}</div>'
        f'<figcaption style="font-size: 12px; line-height: 1.3; text-align: center">{legend}</figcaption></figure>'
    )
    return figure, width, drawer.Height() + 16


@app.function
def largest_common_substructure(mols, timeout=30, unseeded_atoms=60):
    """Query molecule of the largest common substructure of two molecules, or None.

    The search is seeded with the largest part of the first molecule whose Morgan environments the second shares,
    which keeps it fast for peptides. Pairs up to `unseeded_atoms` heavy atoms, and pairs without a seed, are also
    searched unseeded, keeping the larger result. Rings match only whole rings and fused systems of the same shape.
    """
    seed = ""
    try:
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        environments = []
        for mol in mols:
            output = rdFingerprintGenerator.AdditionalOutput()
            output.AllocateAtomToBits()
            generator.GetSparseCountFingerprint(mol, additionalOutput=output)
            environments.append([bits[-1] for bits in output.GetAtomToBits()])
        own_counts = pd.Series(environments[0]).value_counts()
        other_counts = pd.Series(environments[1]).value_counts()
        # atoms of the first molecule whose environment the second has at least as often
        shared = {
            i
            for i, bit in enumerate(environments[0])
            if own_counts[bit] <= other_counts.get(bit, 0)
        }
        cut = [
            bond.GetIdx()
            for bond in mols[0].GetBonds()
            if (bond.GetBeginAtomIdx() in shared) != (bond.GetEndAtomIdx() in shared)
        ]
        parts = Chem.GetMolFrags(
            Chem.FragmentOnBonds(mols[0], cut, addDummies=False) if cut else mols[0]
        )
        seeds = [set(part) for part in parts if set(part) <= shared]
        if seeds:
            seed = Chem.MolFragmentToSmarts(mols[0], atomsToUse=sorted(max(seeds, key=len)))
    except (RuntimeError, ValueError):  # molecules drawn without the valence check may not fingerprint
        seed = ""
    params = rdFMCS.MCSParameters()
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareOrder
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    params.BondCompareParameters.MatchFusedRings = True
    params.BondCompareParameters.MatchFusedRingsStrict = True
    params.Timeout = timeout
    found = []
    if seed:
        params.InitialSeed = seed
        found.append(rdFMCS.FindMCS(mols, params))
    if not seed or max(mol.GetNumAtoms() for mol in mols) <= unseeded_atoms:
        params.InitialSeed = ""
        found.append(rdFMCS.FindMCS(mols, params))
    mcs = max(found, key=lambda result: (result.numAtoms, result.numBonds))
    return Chem.MolFromSmarts(mcs.smartsString) if mcs.numAtoms else None


@app.function
def structure_pair_svg(
    inchi_a,
    inchi_b,
    legends=("GtoPdb", "linked record"),
    missing=("no structure", "no structure"),
    timeout=30,
    scale=22.5,
    max_width=600,
    max_height=1500,
    highlight_above=0.25,
):
    """(HTML, width) of two structures, stacked and laid out alike, with their differences highlighted, each captioned
    with its legend (HTML). A side that cannot be drawn shows the reason (`missing` for an absent InChI), and one drawn
    without the valence check says so under its drawing.

    - same constitution: stereocentres and double bonds defined differently, with R/S and E/Z labels
    - same skeleton: atoms whose fixed hydrogen count differs, and the bonds between them
    - otherwise: atoms outside the largest common substructure (widened to ring systems) in orange, and the bonds
      attaching them in red, unless the similarity is at most `highlight_above`

    Atoms whose formal charge differs are highlighted wherever the atoms correspond. Both structures are drawn at
    `scale` pixels per ångström, reduced to fit `max_width` × `max_height`.
    """
    mols, notes = zip(
        drawable_mol(inchi_a, missing[0]), drawable_mol(inchi_b, missing[1]), strict=True
    )
    mols, notes = list(mols), list(notes)
    drawn = [i for i, mol in enumerate(mols) if mol is not None]
    text_width, gap = 220, 16
    atoms, bonds = [set(), set()], [set(), set()]
    colours = [({}, {}), ({}, {})]  # per structure: atom and bond colours, where not RDKit's default
    pattern = None
    if len(drawn) == 2:
        layers = [parse_standard_inchi(inchi_a), parse_standard_inchi(inchi_b)]
        same_skeleton = share_skeleton(*layers)
        if all(layers[0].get(key) == layers[1].get(key) for key in ("formula", "c", "h")):
            centres, double_bonds = stereo_difference_keys(*comparable_stereo(*layers))
            for i, mol in enumerate(mols):
                components = Chem.GetMolFrags(mol)  # RDKit keeps the InChI component and atom order
                for component, atom in centres:
                    atoms[i].add(components[component][int(atom) - 1])
                for component, bond in double_bonds:
                    ends = [components[component][int(n) - 1] for n in bond.split("-")]
                    atoms[i].update(ends)
                    bonds[i].add(mol.GetBondBetweenAtoms(*ends).GetIdx())
        elif same_skeleton:
            hydrogens = [fixed_hydrogens(layer) for layer in layers]
            differing = {
                key
                for key in hydrogens[0].keys() | hydrogens[1].keys()
                if hydrogens[0].get(key, 0) != hydrogens[1].get(key, 0)
            }
            for i, mol in enumerate(mols):
                components = Chem.GetMolFrags(mol)
                atoms[i] = {
                    components[component][atom - 1]
                    for component, atom in differing
                    if component < len(components) and atom <= len(components[component])
                }
                bonds[i] = {
                    bond.GetIdx()
                    for bond in mol.GetBonds()
                    if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} <= atoms[i]
                }
        else:
            pattern = largest_common_substructure(mols, timeout)
            if pattern is not None:
                differing, attachment = (1.0, 0.85, 0.5), (0.85, 0.15, 0.15)  # orange, red
                for i, mol in enumerate(mols):
                    outside = set(range(mol.GetNumAtoms())) - set(
                        mol.GetSubstructMatch(pattern)
                    )
                    outside = outside.union(
                        *[system for system in ring_systems(mol) if system & outside]
                    )
                    atom_colours, bond_colours = colours[i]
                    for bond in mol.GetBonds():
                        ends = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
                        if ends <= outside:
                            bond_colours[bond.GetIdx()] = differing
                        elif ends & outside:  # joins the differing part to the common one
                            bond_colours[bond.GetIdx()] = attachment
                            atom_colours.update(dict.fromkeys(ends, attachment))
                    atom_colours.update(
                        {atom: differing for atom in outside if atom not in atom_colours}
                    )
                    atoms[i], bonds[i] = set(atom_colours), set(bond_colours)
        # corresponding atoms: by InChI number for the same skeleton, otherwise by the common substructure
        if same_skeleton and mols[0].GetNumAtoms() == mols[1].GetNumAtoms():
            components = [Chem.GetMolFrags(mol) for mol in mols]
            atom_map = [(a, b) for pair in zip(*components) for a, b in zip(*pair)]
        elif pattern is not None:
            atom_map = list(zip(*(mol.GetSubstructMatch(pattern) for mol in mols)))
        else:
            atom_map = []
        for a, b in atom_map:
            if (
                mols[0].GetAtomWithIdx(a).GetFormalCharge()
                != mols[1].GetAtomWithIdx(b).GetFormalCharge()
            ):
                atoms[0].add(a)
                atoms[1].add(b)
        if not same_skeleton and (morgan_similarity(inchi_a, inchi_b) or 0) <= highlight_above:
            # evidently different compounds: aligned, but nothing highlighted
            atoms, bonds = [set(), set()], [set(), set()]
            colours = [({}, {}), ({}, {})]
        for mol in mols:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)  # CIP labels for the annotations
        rdDepictor.Compute2DCoords(mols[0])
        if same_skeleton and mols[0].GetNumAtoms() == mols[1].GetNumAtoms():
            # atom for atom, whatever the bond orders
            rdDepictor.GenerateDepictionMatching2DStructure(mols[1], mols[0], atom_map)
        else:
            # an independent layout if the pattern does not match
            rdDepictor.GenerateDepictionMatching2DStructure(
                mols[1], mols[0], refPatt=pattern, acceptFailure=True
            )
    # one scale for both, reduced until the pair fits
    for _ in range(5):
        figures = {
            i: molecule_figure(
                mols[i],
                legends[i],
                scale,
                atoms[i],
                bonds[i],
                stereo_labels=len(drawn) == 2,
                atom_colours=colours[i][0],
                bond_colours=colours[i][1],
                note=notes[i],
            )
            for i in drawn
        }
        width = max([w for _, w, _ in figures.values()] + [text_width] * (2 - len(drawn)))
        height = sum(h for _, _, h in figures.values()) + gap
        shrink = min(max_width / width, max_height / height if height else 1)
        if shrink >= 1:
            break
        scale *= shrink
    # the missing side is a box with the legend and the reason
    panels = [
        figures[i][0]
        if i in figures
        else f'<div style="box-sizing: border-box; width: {text_width}px; padding: 8px; font-size: 13px; '
        f'color: #52514e; align-self: center"><div style="font-weight: 600">{legends[i]}</div>{html.escape(notes[i])}</div>'
        for i in range(2)
    ]
    body = f'<div style="display: flex; flex-direction: column; gap: {gap}px; align-items: center">{"".join(panels)}</div>'
    return body, width


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### Tables

    marimo tables with GtoPdb, ChEMBL, ChEBI and PubChem IDs shown as links.
    """)
    return


@app.function
def anchor(url, text):
    """Link to `url` that opens in a new tab, styled to stand out in a table cell."""
    return (
        f'<a href="{html.escape(url, quote=True)}" target="_blank" rel="noopener" '
        f'style="color: #1a5fb4; text-decoration: underline;">{html.escape(str(text))}</a>'
    )


@app.function
def link_cell(inner_html):
    """Table cell holding links: wraps like text (marimo renders HTML cells on one line otherwise)."""
    return mo.Html(
        f'<div style="white-space: normal; overflow-wrap: anywhere">{inner_html}</div>'
    )


@app.function
def xref_anchors(db, ids):
    """Links to the records that the |-separated `ids` name in `db`; IDs that are not usable stay plain text."""
    parts = []
    for raw_id in ids.split("|"):
        xref_id = normalise_xref(db, raw_id)
        parts.append(
            anchor(XREF_URLS[db].format(xref_id), raw_id) if xref_id else html.escape(raw_id)
        )
    return " | ".join(parts)


@app.function
def linked_ids(text):
    """`text` as a table cell, with the ChEMBL, ChEBI, PubChem and GtoPdb IDs it names (`ACTION_IDS`) linked to their
    records."""
    if not isinstance(text, str):
        return text

    def link(match):
        chembl, chebi, cid_words, cid, ligand_words, ligands = match.groups()
        if chembl:
            return anchor(XREF_URLS["ChEMBL"].format(chembl), chembl)
        if chebi:
            return anchor(XREF_URLS["ChEBI"].format(chebi), chebi)
        if cid:
            return cid_words + anchor(XREF_URLS["PubChem"].format(cid), cid)
        return ligand_words + ", ".join(
            anchor(GTP_LIGAND_URL.format(i), i) for i in ligands.split(", ")
        )

    return link_cell(ACTION_IDS.sub(link, html.escape(text)))


@app.function
def drawing_caption(source, ids_html, name):
    """Caption of a structure drawing: the source and its linked ID, then the name, cut to two lines (the whole name
    shows on hover)."""
    caption = f"<div>{source} {ids_html}</div>" if ids_html else f"<div>{source}</div>"
    if isinstance(name, str) and name:
        clamp = (
            "display: -webkit-box; -webkit-box-orient: vertical; -webkit-line-clamp: 2; overflow: hidden; "
            "overflow-wrap: anywhere"
        )
        caption += f'<div title="{html.escape(name, quote=True)}" style="{clamp}">{html.escape(name)}</div>'
    return caption


@app.function
def inchi_pair_cell(inchis, labels):
    """Table cell with two InChIs, the first above the second, each after its label, and the layers of the first
    difference (`first_difference_keys`) highlighted in both."""
    try:
        highlight = first_difference_keys(*(parse_standard_inchi(inchi) for inchi in inchis))
    except (AttributeError, ValueError):  # a missing or non-standard InChI: nothing to compare layer by layer
        highlight = set()
    lines = []
    for inchi, label in zip(inchis, labels):
        if not (isinstance(inchi, str) and inchi):
            text = '<span style="font-style: italic">no InChI</span>'
        else:
            try:
                parts = inchi_parts(inchi)
            except ValueError:
                text = html.escape(inchi)
            else:
                text = "InChI=1S" + "".join(
                    "/"
                    + (
                        f'<mark style="background: #ffd43b; color: #1f1f1f">{html.escape(part)}</mark>'
                        if key in highlight
                        else html.escape(part)
                    )
                    for key, part in parts
                )
        lines.append(
            f'<div><span style="font-weight: 600">{html.escape(str(label))}</span> '
            f'<span style="font-family: ui-monospace, monospace; font-size: 12px">{text}</span></div>'
        )
    return link_cell('<div style="margin: 4px 0; opacity: 0.6">vs</div>'.join(lines))


@app.function
def header_with_breaks(name, width, char_width=6.5):
    """`name` with line breaks after underscores or at spaces to fit a column `width` pixels wide (marimo wraps
    headers only at line breaks)."""
    fits = max(1, int((width - 16) / char_width))
    lines = [""]
    for part in re.findall(r"[^_ ]*[_ ]|[^_ ]+$", name):
        if lines[-1] and len((lines[-1] + part).rstrip(" ")) > fits:
            lines.append("")
        lines[-1] += part
    return "\n".join(line.rstrip(" ") for line in lines)


@app.function
def content_width(values, name, char_width=7, low=60, high=220):
    """Width in pixels for a text column not in `COLUMN_WIDTHS`: enough for its 90th-percentile text and its longest
    header word, between `low` and `high`."""
    lengths = values.dropna().astype(str).str.len()
    text = lengths.quantile(0.9) * char_width + 24 if len(lengths) else 0
    word = max(map(len, re.split(r"[_ ]", name))) * 6.5 + 16
    return int(min(high, max(low, text, word)))


@app.function
def linked_table(table, **kwargs):
    """`mo.ui.table` of `table` with its IDs shown as links.

    If it has `gtp_inchi` and `xref_inchi`, a first column, `structures`, draws both, each captioned with its linked ID
    and name, and a last shown column, `inchis`, gives both InChIs with the layers that differ first highlighted
    (`inchi_pair_cell`). The columns these give (`SHOWN_IN_DRAWINGS`) are hidden. The IDs in `suggested_action` are
    links, and the columns of the records it names (`NAMED_IN_ACTION`) are hidden. Hidden columns come last and keep
    their plain values, so the table's search still finds a row by them; they stay in the download, and can be shown
    from the table's column panel.

    A column's width comes from `column_widths`, then `COLUMN_WIDTHS`, then for other text columns their content
    (`content_width`). Text columns wrap within it, and their headers get line breaks to fit. Numeric columns do not
    wrap: marimo's 200 px minimum for wrapped cells would push them out of narrow columns. Other `mo.ui.table` options
    are passed on."""
    columns = list(table.columns)
    hidden = list(kwargs.pop("hidden_columns", []))
    format_mapping = dict(kwargs.pop("format_mapping", {}))
    linked, widths = {}, {}
    if {"gtp_inchi", "xref_inchi"} <= set(columns):
        nothing = [None] * len(table)
        gtp_legends = [
            drawing_caption("GtoPdb", anchor(GTP_LIGAND_URL.format(i), i) if pd.notna(i) else "", name)
            for i, name in zip(table.get("ligand_id", nothing), table.get("name", nothing))
        ]
        if {"xref_db", "xref_id_as_written"} <= set(columns):
            xref_legends = [
                drawing_caption(db, xref_anchors(db, ids) if isinstance(ids, str) else "", name)
                for db, ids, name in zip(
                    table.xref_db, table.xref_id_as_written, table.get("xref_name", nothing)
                )
            ]
        else:
            xref_legends = ["linked record"] * len(table)
        # why a side has no InChI: the `difference` if the table has it, otherwise a generic reason
        differences = table.difference if "difference" in table else [""] * len(table)
        reasons = [
            (
                "GtoPdb gives no InChI for this ligand",
                difference.lower()
                if difference
                in (
                    "Linked record has no structure",
                    "Identifier not resolvable",
                    "Identifier not usable",
                )
                else "the linked record has no structure",
            )
            for difference in differences
        ]
        drawings = [
            structure_pair_svg(g, x, (a, b), missing)
            for g, x, a, b, missing in zip(
                table.gtp_inchi, table.xref_inchi, gtp_legends, xref_legends, reasons
            )
        ]
        linked["structures"] = [mo.Html(drawing) for drawing, _ in drawings]
        # as wide as the widest drawing, or the table squeezes them all
        widths["structures"] = max((width for _, width in drawings), default=300) + 16
        columns.insert(0, "structures")
        # both InChIs in one column, after the others shown
        xref_labels = table.xref_db if "xref_db" in table else ["linked record"] * len(table)
        linked["inchis"] = [
            inchi_pair_cell((g, x), ("GtoPdb", db))
            for g, x, db in zip(table.gtp_inchi, table.xref_inchi, xref_labels)
        ]
        columns.append("inchis")
        hidden += [c for c in SHOWN_IN_DRAWINGS if c in columns]
    if "suggested_action" in columns:
        format_mapping.setdefault("suggested_action", linked_ids)
        hidden += [c for c in NAMED_IN_ACTION if c in columns]
    hidden = list(dict.fromkeys(hidden))
    # hidden columns keep their plain values, which the search reads (it skips HTML cells)
    if "ligand_id" in table and "ligand_id" not in hidden:
        linked["ligand_id"] = [
            link_cell(anchor(GTP_LIGAND_URL.format(i), i)) for i in table.ligand_id
        ]
    for column in ("ligand_ids", "other_gtp_ligands_with_xref_structure"):
        if column in table and column not in hidden:
            linked[column] = [
                link_cell(
                    ", ".join(
                        anchor(GTP_LIGAND_URL.format(i.strip()), i.strip())
                        for i in ids.split(",")
                    )
                )
                if isinstance(ids, str) and ids
                else ids
                for ids in table[column]
            ]
    pair = table.pair.str.split(" – ") if "pair" in table else None
    databases = {
        "xref_id_as_written": table.get("xref_db"),
        "xref_id_normalised": table.get("xref_db"),
        "xref_id_with_gtp_structure": table.get("xref_db"),
        "pubchem_cid": "PubChem",
        "pubchem_cid_with_structure": "PubChem",
        "chembl_parent_id": "ChEMBL",
        "xref_id_a": pair.str[0] if pair is not None else None,
        "xref_id_b": pair.str[1] if pair is not None else None,
    }
    for column, db in databases.items():
        if column not in table or column in hidden or db is None:
            continue
        dbs = [db] * len(table) if isinstance(db, str) else list(db)
        linked[column] = [
            link_cell(xref_anchors(d, ids)) if isinstance(ids, str) and ids else ids
            for d, ids in zip(dbs, table[column])
        ]
    columns = [c for c in columns if c not in hidden] + hidden
    shown = table.assign(**linked)[columns].reset_index(drop=True)
    text_columns = [c for c in columns if not pd.api.types.is_numeric_dtype(shown[c])]
    overrides = kwargs.pop("column_widths", {})
    for c in columns:
        if c in overrides:
            widths[c] = overrides[c]
        elif c in COLUMN_WIDTHS:
            widths[c] = COLUMN_WIDTHS[c]
        elif c in text_columns and c not in widths:
            widths[c] = content_width(table[c], c)
    headers = {c: header_with_breaks(c, widths[c]) if c in widths else c for c in columns}

    def cell_style(_row_id, column, _value):
        return {} if column == headers.get("structures") else TABLE_CELL_STYLE

    return mo.ui.table(
        shown.rename(columns=headers),
        selection=None,
        column_widths={headers[c]: w for c, w in widths.items()},
        wrapped_columns=[headers[c] for c in kwargs.pop("wrapped_columns", text_columns)],
        hidden_columns=[headers[c] for c in hidden] or None,
        format_mapping={headers.get(c, c): f for c, f in format_mapping.items()} or None,
        style_cell=kwargs.pop("style_cell", cell_style),
        **kwargs,
    )


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### Structure comparison and suggested actions

    How a linked record's structure compares with the GtoPdb ligand's, layer by layer, and the suggested action.
    """)
    return


@app.class_definition
class Action(StrEnum):
    """What a suggested action asks the curator to do (`action`)."""

    FIX_ID = "fix ID"
    REPLACE_LINK = "replace link"
    # a record with exactly the GtoPdb structure exists: link it, or change the GtoPdb structure if the link is right
    REPLACE_OR_CHANGE = "replace link or change GtoPdb structure"
    REMOVE_LINK = "remove link"
    CHECK_GTP_STRUCTURE = "check GtoPdb structure"
    CHECK_WHICH = "check which structure is intended"
    REVIEW = "review"
    KEEP_LINK = "keep link"


@app.function
def classify_inchi_pair(gtp_inchi, xref_inchi):
    """How a linked record's standard InChI compares with the GtoPdb ligand's, as `structure_comparison` fields.

    The layers are read in InChI order (`LAYERS`) up to the first that differs: `first_difference` names it and
    `difference_detail` says how it differs. For stereo, `stereo_relation` gives the relation and, where the stereo is
    defined on one side only, `more_stereo_in` and `missing_stereo` the record that defines more and what the other
    lacks. Raises ValueError for a non-standard InChI."""
    # parsed first so that a non-standard InChI is refused even when both are the same
    gtp_layers, xref_layers = parse_standard_inchi(gtp_inchi), parse_standard_inchi(xref_inchi)
    if gtp_inchi == xref_inchi:
        return structure_comparison("identical", similarity=1.0)
    shared, first = compare_layers(gtp_layers, xref_layers)
    if first is None:
        raise ValueError(f"the InChIs differ but none of their layers does: {gtp_inchi}, {xref_inchi}")
    similarity = morgan_similarity(gtp_inchi, xref_inchi)
    if first == "stereo":
        relation, detail = stereo_difference(*comparable_stereo(gtp_layers, xref_layers))
        more_stereo_in, missing = MORE_STEREO_IN.get(relation), None
        if more_stereo_in == "GtoPdb":
            missing = missing_stereo(gtp_inchi, xref_inchi)
        elif more_stereo_in == "linked record":
            missing = missing_stereo(xref_inchi, gtp_inchi)
        return structure_comparison(first, shared, similarity, detail, relation, more_stereo_in, missing)
    if first == "isotopes":
        detail = isotope_difference(gtp_layers, xref_layers)
    elif first in ("charge", "protons"):
        detail = charge_difference(gtp_layers, xref_layers)
    else:
        detail = constitution_difference(first, gtp_inchi, xref_inchi, gtp_layers, xref_layers)
    return structure_comparison(first, shared, similarity, detail)


@app.function
def constitution_difference(first, gtp_inchi, xref_inchi, gtp_layers, xref_layers):
    """How two records whose formula, connections or hydrogens differ first compare. Formula: both formulas, their
    numbers of components, and whether the records have the same atoms and connections. Connections: the atoms
    connected differently. Formula and hydrogens: how their standardised parents compare (`standardised_parents`)."""
    parts = []
    if first == "formula":
        # a metal that standard InChI disconnects reads back with unpaired electrons, so radicals are only reported
        # without one
        metal = has_metal(gtp_inchi) or has_metal(xref_inchi)
        radical = " (radical in structure)" if not metal and radical_count(gtp_inchi) + radical_count(xref_inchi) else ""
        parts.append(f"{gtp_layers['formula']} vs {xref_layers['formula']}{radical}")
        components = [formula_components(layers["formula"]) for layers in (gtp_layers, xref_layers)]
        if components[0] != components[1]:
            parts.append(f"GtoPdb {components[0]} component(s), linked record {components[1]}")
        if has_metal(gtp_inchi) and has_metal(xref_inchi):
            # InChI disconnects metals, so neither parents nor fingerprints describe the complexes
            parts.append("metal complexes, compared by formula only")
            return "; ".join(parts)
        if share_skeleton(gtp_layers, xref_layers):
            # hydrogens missing on one side: a double bond left out, or a radical
            parts.append("same atoms and connections")
    if first == "connections":
        parts.append(
            connectivity_difference(gtp_inchi, xref_inchi)
            or "the atoms connected differently could not be isolated automatically"
        )
    else:
        same_parent, parents = standardised_parents(gtp_inchi, xref_inchi)
        if same_parent or first == "hydrogens":
            parts.append(parents)
    return "; ".join(filter(None, parts))


@app.function
def standardised_parents(gtp_inchi, xref_inchi):
    """(same parent, description) of the two records' parents after standardisation (`parent_inchi`): the same, the
    same as RDKit canonical tautomers ("apart from stereo" where the parents' stereo differs), or different. (None, "")
    when a record has no parent (a metal complex)."""
    inchis = (gtp_inchi, xref_inchi)
    for canonical_tautomer, same in (
        (False, "same parent after standardisation"),
        (True, "same parent as RDKit canonical tautomers"),
    ):
        parents = [parent_inchi(inchi, canonical_tautomer) for inchi in inchis]
        if None in parents:
            return None, ""
        if parents[0] == parents[1]:
            with_stereo = [parent_inchi(inchi, canonical_tautomer, keep_stereo=True) for inchi in inchis]
            return True, same if with_stereo[0] == with_stereo[1] else f"{same}, apart from stereo"
    return False, "different parents, also as RDKit canonical tautomers"


@app.function
def structure_comparison(
    first_difference,
    shared=(),
    similarity=None,
    detail="",
    stereo_relation="",
    more_stereo_in=None,
    missing=None,
):
    """The structure comparison fields of one cross-reference."""
    return {
        "first_difference": first_difference,
        "difference": difference_label(first_difference, shared),
        "difference_detail": detail,
        "stereo_relation": stereo_relation,
        "more_stereo_in": more_stereo_in,
        "missing_stereo": missing,
        "tanimoto_similarity": similarity,
    }


@app.function
def assess_xref(row):
    """`classify_inchi_pair` and `identifier_issues` for one cross-reference (a row of `xref_table`). The structures
    are compared whenever the ID is usable and both exist; otherwise `difference` says why not.

    `metal_differs` flags a structure with a metal (`has_metal`) whose linked record does not have exactly the GtoPdb
    InChI, or has none: standard InChI removes the bonds to the metal, so its layers cannot say how two such structures
    differ. A GtoPdb ligand without a structure, or an identifier that cannot be resolved, is not flagged."""
    if not isinstance(row["xref_id_normalised"], str):
        not_compared = "identifier not usable"
    elif not isinstance(row["gtp_inchi"], str):
        not_compared = "GtoPdb ligand has no structure"
    elif not isinstance(row["xref_inchi"], str):
        not_compared = (
            "linked record has no structure" if row["xref_found"] else "identifier not resolvable"
        )
    else:
        not_compared = None
    if not_compared:
        comparison = structure_comparison("not compared") | {"difference": capitalised(not_compared)}
    else:
        comparison = classify_inchi_pair(row["gtp_inchi"], row["xref_inchi"])
    issues = identifier_issues(row, comparison["first_difference"])
    metal = any(isinstance(inchi, str) and has_metal(inchi) for inchi in (row["gtp_inchi"], row["xref_inchi"]))
    return {
        **comparison,
        "metal_differs": metal
        and (comparison["first_difference"] in LAYERS or not_compared == "linked record has no structure"),
        "identifier_issue": "; ".join(issues.values()),
        "identifier_issue_kinds": "|".join(issues),
    }


@app.function
def difference_label(first_difference, shared):
    """The `shared` layers and the first difference as a sentence: "Identical", "Differs in formula", or e.g.
    "Shares formula, connections and hydrogens; differs in stereo"."""
    if first_difference not in LAYERS:
        return capitalised(first_difference)
    differs = f"differs in {first_difference}"
    if not shared:
        return capitalised(differs)
    listed = f"{', '.join(shared[:-1])} and {shared[-1]}" if len(shared) > 1 else shared[0]
    return f"Shares {listed}; {differs}"


@app.function
def capitalised(text):
    """`text` with its first letter in upper case and the rest as written (unlike `str.capitalize`)."""
    return text[:1].upper() + text[1:]


@app.class_definition
class Step(NamedTuple):
    """One step of a suggested action: what it asks (`Action`), the sentence that says so, and any ID to link."""

    action: str
    text: str
    xref_id: str | None = None


@app.function
def suggested_action(row):
    """Curator's next steps for one cross-reference (a row of the notebook's `curation` table), as {`action`: the
    `Action` of each step, |-separated; `suggested_xref_id`: the ID to link instead, if any; `suggested_action`: the
    steps and notes as sentences}.

    The identifier steps (`identifier_steps`) come first and the structure step (`structure_step`) after, so that
    neither hides the other. A structure step that removes the link leaves out the identifier steps, and one that names
    a record with exactly the GtoPdb structure leaves out the identifier steps that suggest another record. Notes
    (`action_notes`) follow; they never change the action."""
    structure = structure_step(row)
    identifier = identifier_steps(row)
    if structure and structure.action == Action.REMOVE_LINK:
        identifier = []
    elif structure and structure.action == Action.REPLACE_OR_CHANGE:
        identifier = [step for step in identifier if step.action not in (Action.REPLACE_LINK, Action.REMOVE_LINK)]
    steps = [*identifier, *([structure] if structure else [])]
    return {
        "action": "|".join(dict.fromkeys(step.action for step in steps)),
        # the structure step's record before an identifier step's
        "suggested_xref_id": next((step.xref_id for step in reversed(steps) if step.xref_id), None),
        "suggested_action": ". ".join([step.text for step in steps] + action_notes(row)),
    }


@app.function
def identifier_steps(row):
    """One step for each identifier issue (`identifier_issue_kinds`), naming a record with exactly the GtoPdb
    structure where one was found."""
    issues = set(row["identifier_issue_kinds"].split("|"))
    exact, has = exact_records(row)
    exact_ids = row["xref_id_with_gtp_structure"] if exact else None
    hint = f"; {exact} {has} the GtoPdb structure" if exact else ""
    steps = []
    text_to_fix = issues & {IdentifierIssue.INVISIBLE_CHARACTER, IdentifierIssue.PREFIX_CASE}
    normalised = row["xref_id_normalised"]
    if text_to_fix and isinstance(normalised, str):
        named = records(row["xref_db"], normalised, row["xref_name"])
        steps.append(Step(Action.FIX_ID, f"Fix the identifier text: {named}", normalised))
    elif text_to_fix or IdentifierIssue.NOT_VALID in issues:
        steps.append(Step(Action.FIX_ID, f"Fix the identifier text{hint}", exact_ids))
    if IdentifierIssue.CHEMBL_NEVER_ISSUED in issues:
        steps.append(Step(Action.FIX_ID, f"Fix the identifier: it was never issued (likely a typo){hint}", exact_ids))
    if IdentifierIssue.CHEMBL_NOT_A_COMPOUND in issues:
        text = f"Replace: this ChEMBL ID belongs to another kind of ChEMBL record, not a compound{hint}"
        steps.append(Step(Action.REPLACE_LINK, text, exact_ids))
    if IdentifierIssue.CHEMBL_RETIRED in issues:
        if exact:
            text = f"Replace the retired ChEMBL compound ID with {exact}, which {has} the GtoPdb structure"
            steps.append(Step(Action.REPLACE_LINK, text, exact_ids))
        elif isinstance(row["gtp_inchi"], str):
            text = "Replace or remove the retired ChEMBL compound ID: no current ChEMBL compound has exactly the GtoPdb structure"
            steps.append(Step(Action.REMOVE_LINK, text))
        else:
            steps.append(Step(Action.REPLACE_LINK, "Replace the retired ChEMBL compound ID with the current one"))
    if issues & {IdentifierIssue.CHEBI_NOT_FOUND, IdentifierIssue.PUBCHEM_NO_STRUCTURE}:
        steps.append(Step(Action.REMOVE_LINK, "Remove or replace the obsolete identifier"))
    if IdentifierIssue.CHEBI_SECONDARY_ID in issues:
        primary = row.get("chebi_primary_id")
        primary = primary if isinstance(primary, str) else None
        text = "Replace with the primary ChEBI ID"
        steps.append(Step(Action.REPLACE_LINK, f"{text} {primary}" if primary else text, primary))
    if IdentifierIssue.CHEMBL_CHILD_FORM in issues:
        parent = row.get("chembl_parent_id")
        parent = parent if isinstance(parent, str) else None
        text = "Consider linking the ChEMBL parent compound"
        if parent:
            text = f"{text} {records('ChEMBL', parent, row.get('chembl_parent_name'))}"
        steps.append(Step(Action.REPLACE_LINK, text, parent))
    if IdentifierIssue.SEVERAL_IDS_IN_FIELD in issues:
        count = row.get("ids_in_field", 2)
        which = "both IDs" if count <= 2 else f"all {count} IDs"
        steps.append(Step(Action.FIX_ID, f"Check whether {which} are intended"))
    return steps


@app.function
def structure_step(row):
    """The step for the first layer that differs (`first_difference`), or None if the structures do not differ.

    Whatever the layer, a record in the linked database with exactly the GtoPdb structure decides it: GtoPdb already
    linking it makes this link redundant, and otherwise it is linked instead, unless the current link is the intended
    compound and the GtoPdb structure is the one to change. Without one, a link with another GtoPdb ligand's structure
    probably belongs to that ligand (except a labelled ligand's link to the unlabelled compound, which is often a GtoPdb
    ligand too), and linked databases that agree with each other but not with GtoPdb send the curator to the GtoPdb
    structure. Otherwise the layer sets the step (`layer_step`).

    A metal-containing ligand whose linked record has no structure (`metal_differs`) is compared by hand."""
    first = row["first_difference"]
    if row.get("metal_differs") and first not in LAYERS:
        return Step(Action.REVIEW, "Compare the structures by hand: the linked record has no structure to compare with")
    if first not in LAYERS:
        return None
    exact, has = exact_records(row)
    if exact and row.get("exact_record_already_linked"):
        text = f"Remove the cross-reference: GtoPdb also links {exact}, which {has} the GtoPdb structure"
        return Step(Action.REMOVE_LINK, text)
    if exact:
        text = (
            f"Replace with {exact}, which {has} exactly the GtoPdb structure, or change the GtoPdb structure if the "
            f"current link is the intended compound (the current link {difference(row)})"
        )
        return Step(Action.REPLACE_OR_CHANGE, text, row["xref_id_with_gtp_structure"])
    others = other_ligands(row)
    if others and not (first == "isotopes" and unlabelled_link(row)):
        return Step(Action.REMOVE_LINK, f"Remove the cross-reference: the linked record has the structure of {others} instead")
    if row["external_consensus_differs"]:
        text = "Check the GtoPdb structure: all linked databases agree with each other but not with GtoPdb"
        return Step(Action.CHECK_GTP_STRUCTURE, text)
    return layer_step(row, others)


@app.function
def layer_step(row, others):
    """The step for the first layer that differs, where the linked database has no record with exactly the GtoPdb
    structure and no other evidence settles it:

    - formula, connections and charge: check which structure is intended
    - hydrogens: keep the link where the records have the same parent after standardisation (`same_parent`), and
      otherwise check which structure is intended
    - protons: keep the link
    - stereo: check the GtoPdb stereochemistry for an opposite configuration, review stereo defined differently on each
      side, keep a link that leaves stereo undefined, and check the GtoPdb structure where the link defines more
    - isotopes: keep a labelled ligand's link to the unlabelled compound, check the GtoPdb structure of a radioactive
      ligand drawn without its labels, and otherwise check which structure is intended

    Whatever the layer, a structure with a metal (`metal_differs`) sends the curator to check which structure is
    intended, since standard InChI removes the bonds to the metal. A kept link always says why: how it differs, and
    whether the lookup found an exact record."""
    first, db = row["first_difference"], row["xref_db"]
    absent = f"; {db} has no record with exactly the GtoPdb structure" if searched(row) else ""
    if row.get("metal_differs"):
        text = (
            f"Check which structure is intended ({difference(row)}): standard InChI removes the bonds to metals, so it "
            f"cannot tell a different compound from the same one drawn differently{absent}"
        )
        return Step(Action.CHECK_WHICH, text)
    check_which = Step(Action.CHECK_WHICH, f"Check which structure is intended ({difference(row)}){absent}")
    keep = Step(Action.KEEP_LINK, f"Keep the link ({difference(row)}){absent}")
    if first == "hydrogens":
        return keep if same_parent(row["difference_detail"]) else check_which
    if first == "protons":
        return keep
    if first == "stereo":
        relation, missing = row["stereo_relation"], row["missing_stereo"]
        if relation in STEREO_CONFLICT_RELATIONS:
            # a radioactive ligand drawn without its isotope was not searched for, so only the conflict itself is stated
            reason = (
                f"{db} has no record with the GtoPdb stereochemistry"
                if searched(row)
                else "the linked record has the opposite configuration"
            )
            text = f"Check stereochemistry: {reason}, so review the GtoPdb structure"
            return Step(Action.CHECK_GTP_STRUCTURE, text)
        if relation == STEREO_DEFINED_DIFFERENTLY:
            text = "Review: the records define different stereocentres, so the InChIs cannot tell whether they are the same stereoisomer"
            return Step(Action.REVIEW, text)
        if row["more_stereo_in"] == "linked record":
            text = f"Check the GtoPdb structure: the linked record defines stereo that GtoPdb leaves undefined ({missing})"
            return Step(Action.CHECK_GTP_STRUCTURE, text)
        undefined = f"Keep the link (the linked record leaves stereo undefined: {missing})"
        removed = row.get("pubchem_standardisation_removes")
        if isinstance(removed, str):
            text = f"{undefined}; PubChem's standardisation removes this stereo ({removed}) from the GtoPdb structure and reproduces the linked InChI"
            return Step(Action.KEEP_LINK, text)
        no_record = f"; {db} has no record with the GtoPdb stereochemistry" if searched(row) else ""
        return Step(Action.KEEP_LINK, f"{undefined}{no_record}")
    if first == "isotopes":
        if unlabelled_link(row):
            which = f", {others}" if others else ""
            text = f"Keep the link (the linked record is the unlabelled compound{which}); {db} has no record with the GtoPdb isotope labels"
            return Step(Action.KEEP_LINK, text)
        if row["radioactive"] and not labelled(row["gtp_inchi"]):
            text = "Check the GtoPdb structure: this radioactive ligand is drawn without the isotope labels that the linked record carries"
            return Step(Action.CHECK_GTP_STRUCTURE, text)
    return check_which


@app.function
def action_notes(row):
    """Notes on a structure difference, which inform the action without changing it: where a record with exactly the
    GtoPdb structure is named, that the current link has another GtoPdb ligand's structure, or that the linked
    databases agree with each other but not with GtoPdb; that the only record agreeing with GtoPdb is the PubChem
    compound its own deposit is merged into; that a radioactive ligand drawn without its isotope was not searched for;
    and that the linked record of a labelled ligand is unlabelled (where the isotopes are not the layer that differs,
    whose step says so)."""
    first = row["first_difference"]
    if first not in LAYERS:
        return []
    notes = []
    exact = isinstance(row.get("xref_id_with_gtp_structure"), str) and not row.get("exact_record_already_linked")
    others = other_ligands(row)
    if exact and others:
        notes.append(f"Note: the current link has the structure of {others}")
    if exact and row["external_consensus_differs"]:
        notes.append("Note: all linked databases agree with each other but not with GtoPdb")
    if row.get("only_gtp_deposit_agrees"):
        notes.append(
            "Note: the only record that agrees with the GtoPdb structure is the PubChem compound that GtoPdb's own "
            "deposit is merged into"
        )
    if first != "isotopes":
        if not searched(row):
            notes.append(
                "Note: GtoPdb draws this radioactive ligand without its isotope labels, so the linked database was not "
                "searched for a record with the labelled structure"
            )
        elif unlabelled_link(row):
            notes.append("Note: the linked record is unlabelled")
    return notes


@app.function
def same_parent(detail):
    """Whether a `difference_detail` reports the same parent after standardisation (`standardised_parents`): directly,
    as RDKit canonical tautomers, or apart from stereo."""
    return isinstance(detail, str) and "same parent" in detail


@app.function
def difference(row):
    """How the current link differs, as an action says it: 'differs in formula: C13H18O2 vs C14H20O2'."""
    first, detail = row["first_difference"], row["difference_detail"]
    return f"differs in {first}: {detail}" if detail else f"differs in {first}"


@app.function
def other_ligands(row):
    """The other GtoPdb ligands with the linked record's structure, as an action names them, or None."""
    others = row.get("other_gtp_ligands_with_xref_structure")
    if not (isinstance(others, str) and others):
        return None
    return f"GtoPdb ligand{'s' if ',' in others else ''} {others}"


@app.function
def labelled(inchi):
    """Whether a standard InChI carries isotope labels."""
    return isinstance(inchi, str) and "/i" in inchi


@app.function
def unlabelled_link(row):
    """Whether the GtoPdb structure carries isotope labels and the linked record none."""
    return labelled(row["gtp_inchi"]) and not labelled(row["xref_inchi"])


@app.function
def searched(row):
    """Whether the linked database was searched for a record with exactly the GtoPdb structure: not for a radioactive
    ligand drawn without its isotope, whose InChIKey is the unlabelled compound's."""
    return not (row["radioactive"] and not labelled(row["gtp_inchi"]))


@app.function
def exact_records(row):
    """(the records in the linked database with exactly the GtoPdb structure, as an action names them, and the verb
    agreeing with their number), or (None, None) if none was found."""
    ids = row.get("xref_id_with_gtp_structure")
    if not isinstance(ids, str):
        return None, None
    has = {1: "has", 2: "both have"}.get(len(ids.split("|")), "all have")
    return records(row["xref_db"], ids, row.get("xref_name_with_gtp_structure")), has


@app.function
def records(db, ids, names):
    """The |-separated `ids` of `db`, each named by `record_description` from the |-separated `names`, joined with
    'or'."""
    names = names.split("|") if isinstance(names, str) else []
    return " or ".join(
        record_description(db, xref_id, names[i] if i < len(names) else None)
        for i, xref_id in enumerate(ids.split("|"))
    )


@app.function
def record_description(db, xref_id, name):
    """A record as a suggested action names it: '(name, ID)', or the ID alone if unnamed (or titled 'CID n').
    PubChem IDs are written 'PubChem CID n'."""
    label = f"PubChem CID {xref_id}" if db == "PubChem" else xref_id
    if isinstance(name, str) and name and name != f"CID {xref_id}":
        return f"({name}, {label})"
    return label


@app.function
def inorganic_structure_action(linked_cids, cids_with_structure, cid_names):
    """Curator's next step for an inorganic ligand without a structure: add it, and link a PubChem compound that has
    it if no linked one does."""
    add = "Add this SMILES, InChI and InChIKey"
    if not cids_with_structure:
        return f"{add} (PubChem has no compound with this structure)"
    if set(linked_cids) & set(cids_with_structure):
        return add
    records = " or ".join(
        record_description("PubChem", cid, cid_names.get(cid)) for cid in cids_with_structure
    )
    has = "has" if len(cids_with_structure) == 1 else "have"
    if linked_cids:
        return f"{add}, and replace the PubChem link with {records}, which {has} this structure"
    return f"{add}, and link {records}"


@app.function
def smiles_inchi_problem(smiles, inchi):
    """(InChI from a GtoPdb SMILES, how it disagrees with GtoPdb's own InChI or "", suggested action). The InChI is
    compared as the GtoPdb structure and the SMILES as the linked record."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        unsanitised = Chem.MolFromSmiles(smiles, sanitize=False)
        if unsanitised is None:
            return None, "SMILES cannot be parsed", "Fix the SMILES"
        try:
            Chem.SanitizeMol(unsanitised)
        except Chem.MolSanitizeException as error:  # e.g. RDKit names the atom whose bonds exceed its valence
            return (
                None,
                f"invalid SMILES: {error}",
                "Fix the SMILES: it is not a valid structure",
            )
        return None, "SMILES not accepted by RDKit", "Fix the SMILES"
    smiles_inchi = Chem.MolToInchi(mol)
    if not smiles_inchi:
        return None, "no InChI can be generated from the SMILES", "Fix the SMILES"
    if smiles_inchi == inchi:
        return smiles_inchi, "", ""
    result = classify_inchi_pair(inchi, smiles_inchi)
    if result["more_stereo_in"] == "GtoPdb":
        return (
            smiles_inchi,
            f"the InChI defines stereo that the SMILES lacks ({result['missing_stereo']})",
            "Add the missing stereo to the SMILES",
        )
    if result["more_stereo_in"] == "linked record":
        return (
            smiles_inchi,
            f"the SMILES defines stereo that the InChI lacks ({result['missing_stereo']})",
            "Regenerate the InChI from the SMILES",
        )
    if result["first_difference"] == "isotopes":
        isotopes = [isotope_layer(parse_standard_inchi(text)) for text in (inchi, smiles_inchi)]
        return (
            smiles_inchi,
            f"isotope layers differ: InChI {isotopes[0]}, SMILES {isotopes[1]}",
            "Make the SMILES and InChI describe the same isotopes",
        )
    detail = (
        result["difference_detail"]
        .replace("GtoPdb", "InChI")
        .replace("linked record", "SMILES")
    )
    return (
        smiles_inchi,
        f"{result['difference']}: {detail}" if detail else result["difference"],
        "Make the SMILES and InChI describe the same structure",
    )


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### UniChem checks

    UniChem checks: how UniChem holds GtoPdb's structures.
    """)
    return


@app.function
def unichem_status(gtp_inchi, unichem_inchi):
    """How the structure UniChem holds for a GtoPdb ligand compares with GtoPdb's own InChI."""
    if not isinstance(unichem_inchi, str):
        return "not in UniChem"
    if unichem_inchi == gtp_inchi:
        return "identical"
    if unichem_inchi == strip_layers(gtp_inchi, STEREO_LAYERS):
        return "stereo layers missing in UniChem"
    return "other difference"


@app.function
def unichem_structure_source(gtp_inchi, smiles_inchi, unichem_inchi):
    """For a GtoPdb ligand whose SMILES and InChI disagree, which of the two structures UniChem holds, compared without
    their stereo layers. `smiles_inchi` is the InChI generated from the SMILES, None for an invalid SMILES."""
    if not isinstance(unichem_inchi, str):
        return "not in UniChem"
    held = strip_layers(unichem_inchi, STEREO_LAYERS)
    from_inchi = held == strip_layers(gtp_inchi, STEREO_LAYERS)
    if not isinstance(smiles_inchi, str):
        return "InChI (the SMILES is invalid)" if from_inchi else "neither (the SMILES is invalid)"
    from_smiles = held == strip_layers(smiles_inchi, STEREO_LAYERS)
    if from_inchi and from_smiles:
        return "either (the two differ only in stereo)"
    return "SMILES" if from_smiles else "InChI" if from_inchi else "neither"


if __name__ == "__main__":
    app.run()
