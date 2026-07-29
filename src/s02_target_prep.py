"""
s02_target_prep.py  --  Phase A, step A3 of IN_SILICO_PLAN.md

The plan listed candidate PDB IDs marked "verify".  This script does the
verifying against the RCSB Data API, and it applies per-target acceptance
criteria rather than a blanket rule.

A naive "highest resolution wins" rule was tried first and produced
biologically WRONG picks, all of which the criteria below now catch:

  * ERK2  -> 4QTB is ERK1 (MAPK3), a different kinase.       [gene check]
  * COX-2 -> 3NT1 / 1CX2 are murine, not human.              [organism check]
  * KRAS  -> 4OBE contains only GDP, the natural nucleotide;
             no inhibitor to define an allosteric docking box. [cofactor check]
  * TopoII-> 1ZXM is the ATPase domain (AMP-PNP), the wrong
             site entirely for doxorubicin, which intercalates
             at the DNA cleavage complex.                     [site check]
  * Tubulin / mTOR -> genuinely only solved at ~3.5 A; a blanket
             3.0 A cutoff wrongly discarded the whole target. [per-target res]

Acceptance requires: gene-symbol match, allowed organism, resolution within
the target-specific limit, and at least one DRUG-LIKE co-crystallised ligand
(cofactors and nucleotides do not count) so the Gate G1 redocking check
(RMSD <= 2.0 A) can actually be run.

Output: data/processed/target_panel_verified.csv
        data/processed/target_panel_selected.csv
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
PROC.mkdir(parents=True, exist_ok=True)

DATA_API = "https://data.rcsb.org/rest/v1/core"
CACHE = RAW / "rcsb_cache.json"

# Crystallisation additives, ions, buffers -- never a docking-box ligand.
JUNK = {
    "HOH", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "MPD", "DMS", "ACT",
    "CL", "NA", "MG", "ZN", "CA", "K", "MN", "FE", "NI", "CD", "IOD", "BR",
    "TRS", "EPE", "MES", "IMD", "FMT", "ACE", "NH4", "CO3", "NO3", "AZI",
    "BME", "DTT", "CIT", "TLA", "MLI", "SCN", "PGE", "1PE", "P6G", "OLC",
    "2PE", "BOG", "NAG", "MAN", "BMA", "FUC", "GAL", "SIA", "XYP",
}

# Cofactors / nucleotides: real molecules, but they mark the natural
# cofactor site, not a druggable inhibitor pocket. A structure whose ONLY
# ligand is one of these cannot anchor an inhibitor docking box.
COFACTORS = {
    "FAD", "FMN", "NAD", "NAP", "NDP", "NAI", "SAM", "SAH", "COA", "PLP",
    "GDP", "GTP", "GNP", "GSP", "ADP", "ATP", "ANP", "ACP", "AMP", "APC",
    "HEM", "HEC", "SF4", "FES", "MGF", "ALF", "BEF",
}


def _t(name, genes, candidates, *, organisms=("homo sapiens",),
       max_res=3.0, note=""):
    return {"name": name, "genes": {g.upper() for g in genes},
            "organisms": {o.lower() for o in organisms},
            "candidates": candidates, "max_res": max_res, "note": note}


TARGETS = [
    _t("EGFR (kinase domain)", ["EGFR", "ERBB1"], ["1M17", "4HJO", "1XKK"]),
    _t("PI3Ka (PIK3CA)", ["PIK3CA"], ["4JPS", "4L23", "5DXT"]),
    _t("AKT1", ["AKT1"], ["4EJN", "3O96", "4GV1"]),
    _t("mTOR (kinase)", ["MTOR", "FRAP1"], ["4JSV", "4JT6", "4JT5", "4DRI", "4DRH"],
       max_res=3.8, note="mTOR kinase domain is only solved at ~3.0-3.6 A"),
    _t("ERK2 (MAPK1)", ["MAPK1", "ERK2", "PRKM1"], ["2OJG", "6G54", "4QTB", "5NHO"],
       note="4QTB is ERK1/MAPK3 -- must be rejected by the gene check"),
    _t("KRAS (switch-II pocket)", ["KRAS", "KRAS2"], ["6OIM", "7RPZ", "6GJ8", "4OBE"],
       note="needs a covalent/allosteric inhibitor, not just GDP"),
    _t("Bcl-2", ["BCL2"], ["4LVT", "6O0K", "6QGG"]),
    _t("Bcl-xL", ["BCL2L1"], ["3QKD", "4QVX", "2YXJ"]),
    _t("Caspase-3", ["CASP3", "CPP32"], ["1NME", "2XYG", "4QTX", "3KJF"]),
    _t("COX-2 (PTGS2)", ["PTGS2", "COX2", "COX-2"], ["5KIR", "5F19", "3NT1", "1CX2"],
       note="3NT1/1CX2 are murine"),
    _t("VEGFR2 (KDR)", ["KDR", "FLK1"], ["4ASD", "3VHE", "2OH4"]),
    _t("Tubulin (colchicine site)", ["TUBB", "TUBA1B", "STMN4", "TTL"],
       ["1SA0", "4O2B", "5LYJ", "402B"],
       organisms=("bos taurus", "rattus norvegicus", "gallus gallus", "ovis aries",
                  "sus scrofa", "homo sapiens"),
       max_res=4.0,
       note="tubulin is solved from bovine/ovine brain; >99% identical to human"),
    _t("Topoisomerase II (DNA cleavage complex)", ["TOP2A", "TOP2B", "TOP2"],
       ["5GWK", "3QX3", "4G0U"],
       organisms=("homo sapiens", "synthetic construct"),
       max_res=3.6,
       note="doxorubicin acts at the DNA cleavage complex, NOT the ATPase "
            "domain -- 1ZXM (AMP-PNP) is the wrong site"),
    _t("KEAP1 (Kelch domain)", ["KEAP1", "INRF2", "KLHL19"],
       ["4L7B", "4XMB", "5FNU", "4IQK"]),
    _t("Thioredoxin reductase 1", ["TXNRD1"], ["2ZZ0", "3QFA", "3QFB"],
       note="only FAD is present -- flagged as cofactor-only"),
    _t("NQO1", ["NQO1"], ["2F1O", "1D4A", "1H69"]),
]


def _cache_load() -> dict:
    return json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}


def _cache_save(c: dict) -> None:
    CACHE.write_text(json.dumps(c, indent=1), encoding="utf-8")


def _get(url: str, cache: dict, pause: float = 0.12):
    if url in cache:
        return cache[url]
    try:
        r = requests.get(url, timeout=30)
        time.sleep(pause)
        val = r.json() if r.status_code == 200 else None
    except Exception:                                           # noqa: BLE001
        val = None
    cache[url] = val
    return val


def inspect_entry(pdb_id: str, cache: dict) -> dict | None:
    e = _get(f"{DATA_API}/entry/{pdb_id}", cache)
    if not e:
        return None

    info = e.get("rcsb_entry_info", {})
    res = info.get("resolution_combined")
    resolution = float(res[0]) if isinstance(res, list) and res else None
    method = ", ".join(m.get("method", "") for m in e.get("exptl", []))
    title = e.get("struct", {}).get("title", "")
    ids = e.get("rcsb_entry_container_identifiers", {})

    names, genes, organisms = [], [], []
    for pe in ids.get("polymer_entity_ids", []) or []:
        p = _get(f"{DATA_API}/polymer_entity/{pdb_id}/{pe}", cache)
        if not p:
            continue
        nm = p.get("rcsb_polymer_entity", {}).get("pdbx_description")
        if nm:
            names.append(nm)
        for g in p.get("rcsb_entity_source_organism", []) or []:
            if g.get("scientific_name"):
                organisms.append(g["scientific_name"])
            for n in g.get("rcsb_gene_name", []) or []:
                if n.get("value"):
                    genes.append(n["value"])

    druglike, cofac = [], []
    for ne in ids.get("non_polymer_entity_ids", []) or []:
        n = _get(f"{DATA_API}/nonpolymer_entity/{pdb_id}/{ne}", cache)
        if not n:
            continue
        cid = n.get("pdbx_entity_nonpoly", {}).get("comp_id")
        if not cid or cid in JUNK:
            continue
        (cofac if cid in COFACTORS else druglike).append(cid)

    return {
        "pdb_id": pdb_id,
        "title": title[:90],
        "method": method,
        "resolution_A": resolution,
        "protein": "; ".join(dict.fromkeys(names))[:90],
        "genes": ",".join(dict.fromkeys(genes))[:44],
        "organism": "; ".join(dict.fromkeys(organisms))[:44],
        "druglike_ligands": ",".join(druglike),
        "cofactor_ligands": ",".join(cofac),
        "n_druglike": len(druglike),
        "contains_dna": any("DNA" in (nm or "").upper() for nm in names),
    }


def judge(rec: dict | None, spec: dict) -> tuple[str, str]:
    if rec is None:
        return "REJECT", "entry not retrievable"

    problems = []

    # --- gene identity: the check that catches ERK1-for-ERK2 ---
    found = {g.strip().upper() for g in rec["genes"].split(",") if g.strip()}
    # normalise things like "COX-2 PGHS-B"
    found |= {tok for g in found for tok in g.replace("-", "").split()}
    if found and not (found & {g.replace("-", "") for g in spec["genes"]} or
                      found & spec["genes"]):
        problems.append(f"gene mismatch: got {sorted(found)[:3]}, "
                        f"want {sorted(spec['genes'])}")

    # --- organism (case-insensitive: RCSB returns 'HOMO SAPIENS' sometimes) ---
    orgs = {o.strip().lower() for o in rec["organism"].split(";") if o.strip()}
    if orgs and not (orgs & spec["organisms"]):
        problems.append(f"organism {sorted(orgs)[:2]} not allowed")

    # --- resolution, per target ---
    if rec["resolution_A"] is None:
        problems.append("no resolution (NMR/other)")
    elif rec["resolution_A"] > spec["max_res"]:
        problems.append(f"resolution {rec['resolution_A']:.2f} > {spec['max_res']} A")

    # --- a drug-like ligand is required for the G1 redocking check ---
    if rec["n_druglike"] == 0:
        if rec["cofactor_ligands"]:
            problems.append(f"cofactor-only ligand ({rec['cofactor_ligands']}) "
                            f"-- no inhibitor pocket to redock")
        else:
            problems.append("no co-crystallised ligand -> G1 redocking impossible")

    return ("ACCEPT" if not problems else "REJECT"), "; ".join(problems) or "ok"


def main() -> pd.DataFrame:
    cache = _cache_load()
    rows = []
    print("Verifying candidate structures against the RCSB Data API")
    print("(per-target gene / organism / resolution / ligand criteria)\n")

    for spec in TARGETS:
        print(f"--- {spec['name']}")
        if spec["note"]:
            print(f"    note: {spec['note']}")
        for pid in spec["candidates"]:
            rec = inspect_entry(pid, cache)
            v, why = judge(rec, spec)
            base = {"target": spec["name"], "pdb_id": pid,
                    "verdict": v, "reason": why}
            if rec is None:
                rows.append(base)
                print(f"    {pid}  {v:6s}  {why}")
                continue
            rows.append({**rec, **base})
            resd = f"{rec['resolution_A']:.2f}A" if rec["resolution_A"] else "  n/a"
            print(f"    {pid}  {v:6s} {resd:>7s}  "
                  f"lig=[{rec['druglike_ligands'] or '-'}]  {why if v=='REJECT' else ''}")
        print()

    _cache_save(cache)
    df = pd.DataFrame(rows)
    df.to_csv(PROC / "target_panel_verified.csv", index=False)

    acc = df[df["verdict"] == "ACCEPT"].copy()
    best = (acc.sort_values(["target", "resolution_A"])
               .groupby("target", as_index=False).first())

    print("=" * 92)
    print("SELECTED PANEL -- best ACCEPTED structure per target")
    print("=" * 92)
    print(best[["target", "pdb_id", "resolution_A", "genes",
                "druglike_ligands"]].to_string(index=False))
    best.to_csv(PROC / "target_panel_selected.csv", index=False)

    missing = sorted({s["name"] for s in TARGETS} - set(best["target"]))
    print(f"\nAccepted {len(best)} / {len(TARGETS)} targets.")
    if missing:
        print("NO ACCEPTED STRUCTURE:")
        for m in missing:
            sub = df[df["target"] == m]
            print(f"  - {m}")
            for _, r in sub.iterrows():
                print(f"      {r['pdb_id']}: {r['reason']}")
        print("  -> widen candidates, relax the target-specific limit, or drop"
              " the target and record that decision in the manuscript.")
    return df


if __name__ == "__main__":
    main()
