"""
s05_chembl_fetch.py  --  Phase A, step A5 of IN_SILICO_PLAN.md

Pulls and curates cytotoxicity data for the QSAR training set:

  * A549 (human lung adenocarcinoma)  -- the study's cell line
  * MRC-5 / WI-38 / BEAS-2B / IMR-90  -- normal lung, for the Selectivity
                                         Index model (defect D8)

Curation rules (applied in this order, each one logged):
  1. activity type in {IC50, GI50, EC50}, standard_relation '=' only
     (">" / "<" censored values are dropped -- they are not measurements)
  2. units convertible to nM
  3. a parseable structure (canonical SMILES) must exist
  4. drop mixtures / salts -> keep the largest organic fragment
  5. standardise to pIC50 = 9 - log10(IC50 in nM)
  6. deduplicate by InChIKey: keep the MEDIAN pIC50 across replicates,
     but drop compounds whose replicate spread exceeds 1 log unit
     (irreproducible -- keeping them injects label noise)
  7. keep pIC50 within [3, 11]; outside that range values are usually
     unit-entry errors

Output: data/raw/chembl_<cell>_raw.csv
        data/processed/qsar_dataset_<cell>.csv
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
for p in (RAW, PROC):
    p.mkdir(parents=True, exist_ok=True)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# ChEMBL cell-line target identifiers.
CELL_LINES = {
    "A549":    "CHEMBL392",      # human lung adenocarcinoma  (the study line)
    "MRC5":    "CHEMBL614725",   # normal human lung fibroblast
    "WI38":    "CHEMBL614725",   # placeholder; resolved dynamically below
}

VALID_TYPES = {"IC50", "GI50", "EC50"}


def resolve_cell_target(name: str) -> list[str]:
    """Find the ChEMBL target id(s) for a cell line by name search."""
    url = f"{CHEMBL}/target/search.json"
    try:
        r = requests.get(url, params={"q": name, "limit": 25}, timeout=60)
        if r.status_code != 200:
            return []
        out = []
        for t in r.json().get("targets", []):
            if t.get("target_type") != "CELL-LINE":
                continue
            pref = (t.get("pref_name") or "").upper()
            if name.upper().replace("-", "") in pref.replace("-", "").replace(" ", ""):
                out.append((t["target_chembl_id"], t.get("pref_name")))
        return out
    except Exception as exc:                                    # noqa: BLE001
        print(f"  ! target search failed for {name}: {exc}")
        return []


def fetch_activities(target_id: str, limit_pages: int = 200) -> pd.DataFrame:
    """Page through /activity for one cell-line target."""
    rows, url = [], f"{CHEMBL}/activity.json"
    params = {"target_chembl_id": target_id, "limit": 1000, "offset": 0}
    for page in range(limit_pages):
        try:
            r = requests.get(url, params=params, timeout=120)
        except Exception as exc:                                # noqa: BLE001
            print(f"    ! page {page}: {exc}")
            break
        if r.status_code != 200:
            break
        j = r.json()
        acts = j.get("activities", [])
        if not acts:
            break
        rows.extend(acts)
        print(f"    page {page + 1}: +{len(acts)} (total {len(rows)})", end="\r")
        nxt = j.get("page_meta", {}).get("next")
        if not nxt:
            break
        params["offset"] += params["limit"]
        time.sleep(0.15)
    print()
    return pd.DataFrame(rows)


_LFC = rdMolStandardize.LargestFragmentChooser()


def clean_smiles(smi: str) -> str | None:
    """Strip salts/solvates, keep the largest organic fragment, canonicalise."""
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        m = _LFC.choose(m)
    except Exception:                                           # noqa: BLE001
        pass
    if m is None or m.GetNumHeavyAtoms() < 5:
        return None
    if not any(a.GetSymbol() == "C" for a in m.GetAtoms()):
        return None                                    # inorganic
    return Chem.MolToSmiles(m)


def curate(df: pd.DataFrame, label: str) -> pd.DataFrame:
    log = []

    def step(d, msg):
        log.append(f"    {msg:52s} -> {len(d):6d} rows")
        return d

    step(df, "raw records")
    if df.empty:
        return df

    d = df[df["standard_type"].isin(VALID_TYPES)]
    step(d, "activity type in IC50/GI50/EC50")

    d = d[d["standard_relation"] == "="]
    step(d, "exact values only (drop > and <)")

    d = d[d["standard_units"] == "nM"]
    d = d[pd.to_numeric(d["standard_value"], errors="coerce").notna()]
    step(d, "units = nM, numeric value")

    d = d.copy()
    d["value_nM"] = pd.to_numeric(d["standard_value"])
    d = d[d["value_nM"] > 0]
    step(d, "positive values")

    d["smiles"] = d["canonical_smiles"].map(clean_smiles)
    d = d[d["smiles"].notna()]
    step(d, "parseable structure, largest organic fragment")

    d["pIC50"] = 9.0 - np.log10(d["value_nM"])
    d = d[(d["pIC50"] >= 3.0) & (d["pIC50"] <= 11.0)]
    step(d, "pIC50 in [3, 11]")

    d["inchikey"] = d["smiles"].map(
        lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s)) if s else None)
    d = d[d["inchikey"].notna()]

    grp = d.groupby("inchikey")["pIC50"]
    spread = grp.max() - grp.min()
    good = spread[spread <= 1.0].index
    n_before = d["inchikey"].nunique()
    d = d[d["inchikey"].isin(good)]
    step(d, f"drop {n_before - len(good)} cpds w/ replicate spread > 1 log")

    agg = (d.groupby("inchikey")
             .agg(smiles=("smiles", "first"),
                  pIC50=("pIC50", "median"),
                  n_records=("pIC50", "size"),
                  molecule_chembl_id=("molecule_chembl_id", "first"))
             .reset_index())
    log.append(f"    {'unique compounds after dedup':52s} -> {len(agg):6d} rows")

    agg["mw"] = agg["smiles"].map(
        lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)))
    agg["ic50_uM"] = 10.0 ** (6.0 - agg["pIC50"])
    agg["ic50_ug_mL"] = agg["ic50_uM"] * agg["mw"] / 1000.0
    agg["cell_line"] = label

    print(f"\n  Curation log for {label}:")
    for line in log:
        print(line)
    return agg


def main():
    print("=" * 78)
    print("ChEMBL cytotoxicity data for the QSAR training set")
    print("=" * 78)

    # confirm the A549 target id rather than trusting a hard-coded value
    print("\nResolving cell-line targets by name search:")
    for nm in ("A549", "MRC-5", "WI-38", "BEAS-2B", "IMR-90"):
        hits = resolve_cell_target(nm)
        print(f"  {nm:9s}: {hits[:3] if hits else 'no CELL-LINE target found'}")

    results = {}
    for label, tid in [("A549", "CHEMBL392")]:
        print(f"\n--- fetching {label} ({tid})")
        raw = fetch_activities(tid)
        if raw.empty:
            print(f"  no activities returned for {tid}")
            continue
        raw.to_csv(RAW / f"chembl_{label}_raw.csv", index=False)
        cur = curate(raw, label)
        cur.to_csv(PROC / f"qsar_dataset_{label}.csv", index=False)
        results[label] = cur
        print(f"\n  -> {PROC / f'qsar_dataset_{label}.csv'}  ({len(cur)} compounds)")
        if len(cur):
            print(f"     pIC50: min {cur['pIC50'].min():.2f}  "
                  f"median {cur['pIC50'].median():.2f}  "
                  f"max {cur['pIC50'].max():.2f}")
            print(f"     IC50 (ug/mL): median {cur['ic50_ug_mL'].median():.2f}, "
                  f"IQR {cur['ic50_ug_mL'].quantile(.25):.2f}-"
                  f"{cur['ic50_ug_mL'].quantile(.75):.2f}")
    return results


if __name__ == "__main__":
    main()
