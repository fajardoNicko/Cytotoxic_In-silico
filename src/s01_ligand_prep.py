"""
s01_ligand_prep.py  --  Phase A, step A2 of IN_SILICO_PLAN.md

Builds the candidate-compound library for Moringa oleifera (bark-weighted),
plus the reference drugs, by pulling canonical structures from PubChem and
computing RDKit descriptors / drug-likeness locally.

This runs BEFORE the extract arrives.  When the GC-MS/LC-MS peak table comes
back from MSU-IIT it is reconciled against this library by PubChem CID
(see s05_reconcile.py), and anything new is appended.

Output: data/processed/ligand_library.csv
        data/structures/ligands.sdf   (3D, MMFF-minimised)
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
STRUCT = ROOT / "data" / "structures"
RAW = ROOT / "data" / "raw"
for p in (PROC, STRUCT, RAW):
    p.mkdir(parents=True, exist_ok=True)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CACHE = RAW / "pubchem_cache.json"


# ---------------------------------------------------------------------------
# the library.  `part` records where in the plant the compound is reported;
# bark/stembark entries are the ones most relevant to this study.
# ---------------------------------------------------------------------------

MORINGA_COMPOUNDS = [
    # --- glucosinolates / isothiocyanates: the signature Moringa chemotype ---
    ("moringin",                     "bark/seed", "isothiocyanate"),
    ("glucomoringin",                "bark/seed", "glucosinolate"),
    ("benzyl isothiocyanate",        "bark/root", "isothiocyanate"),
    ("4-hydroxybenzyl isothiocyanate", "root",    "isothiocyanate"),
    ("sinigrin",                     "seed",      "glucosinolate"),
    ("niazimicin",                   "bark",      "thiocarbamate"),
    ("niazirin",                     "bark",      "nitrile glycoside"),
    # --- alkaloids reported from bark ---
    ("moringine",                    "bark",      "alkaloid"),
    ("moringinine",                  "bark",      "alkaloid"),
    ("spirochin",                    "root/bark", "alkaloid"),
    # --- flavonoids ---
    ("quercetin",                    "leaf/bark", "flavonol"),
    ("kaempferol",                   "leaf/bark", "flavonol"),
    ("rutin",                        "leaf",      "flavonol glycoside"),
    ("isoquercitrin",                "leaf",      "flavonol glycoside"),
    ("apigenin",                     "leaf",      "flavone"),
    ("luteolin",                     "leaf",      "flavone"),
    ("catechin",                     "bark",      "flavan-3-ol"),
    ("epicatechin",                  "bark",      "flavan-3-ol"),
    # --- phenolic acids / tannin precursors (bark is tannin-rich) ---
    ("gallic acid",                  "bark",      "phenolic acid"),
    ("ellagic acid",                 "bark",      "phenolic acid"),
    ("chlorogenic acid",             "leaf/bark", "phenolic acid"),
    ("caffeic acid",                 "bark",      "phenolic acid"),
    ("ferulic acid",                 "bark",      "phenolic acid"),
    ("p-coumaric acid",              "bark",      "phenolic acid"),
    ("syringic acid",                "bark",      "phenolic acid"),
    ("vanillin",                     "bark",      "phenolic aldehyde"),
    ("protocatechuic acid",          "bark",      "phenolic acid"),
    # --- sterols / triterpenes ---
    ("beta-sitosterol",              "bark",      "phytosterol"),
    ("stigmasterol",                 "bark",      "phytosterol"),
    ("campesterol",                  "bark",      "phytosterol"),
    ("lupeol",                       "bark",      "triterpene"),
    ("oleanolic acid",               "bark",      "triterpene"),
    ("ursolic acid",                 "bark",      "triterpene"),
    ("beta-amyrin",                  "bark",      "triterpene"),
    # --- fatty acids commonly dominating GC-MS of ethanolic extracts ---
    ("palmitic acid",                "all",       "fatty acid"),
    ("oleic acid",                   "all",       "fatty acid"),
    ("linoleic acid",                "all",       "fatty acid"),
    ("stearic acid",                 "all",       "fatty acid"),
    ("octacosanoic acid",            "bark",      "fatty acid"),
    # --- misc reported constituents ---
    ("4-hydroxymellein",             "bark",      "isocoumarin"),
    ("vanillic acid",                "bark",      "phenolic acid"),
    ("scopoletin",                   "bark",      "coumarin"),
    ("squalene",                     "all",       "terpene"),
    ("phytol",                       "leaf",      "diterpene alcohol"),
]

REFERENCE_DRUGS = [
    ("doxorubicin",  "reference", "positive control (manuscript)"),
    ("paclitaxel",   "reference", "positive control (named once, p.4)"),
    ("cisplatin",    "reference", "NSCLC standard of care"),
    ("erlotinib",    "reference", "EGFR TKI benchmark"),
    ("etoposide",    "reference", "topoisomerase II benchmark"),
]


# ---------------------------------------------------------------------------
# PubChem retrieval, with an on-disk cache so reruns are offline & reproducible
# ---------------------------------------------------------------------------

def _load_cache() -> dict:
    if CACHE.exists():
        return json.loads(CACHE.read_text(encoding="utf-8"))
    return {}


def _save_cache(c: dict) -> None:
    CACHE.write_text(json.dumps(c, indent=1), encoding="utf-8")


def fetch_pubchem(name: str, cache: dict, pause: float = 0.22) -> dict | None:
    """Resolve a compound name to CID + canonical SMILES + formula + MW."""
    if name in cache:
        return cache[name]

    # PubChem renamed these properties: the current API returns "SMILES"
    # (isomeric) and "ConnectivitySMILES" (flat). The old CanonicalSMILES /
    # IsomericSMILES names are still accepted in the REQUEST but are absent
    # from the RESPONSE, so both are requested and all spellings are tried.
    props = ("MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,"
             "CanonicalSMILES,IsomericSMILES,IUPACName,InChIKey")
    url = f"{PUBCHEM}/compound/name/{requests.utils.quote(name)}/property/{props}/JSON"
    try:
        r = requests.get(url, timeout=30)
        time.sleep(pause)                      # PubChem asks for <= 5 req/s
        if r.status_code != 200:
            cache[name] = None
            return None
        rec = r.json()["PropertyTable"]["Properties"][0]
        smiles = (rec.get("SMILES") or rec.get("IsomericSMILES")
                  or rec.get("ConnectivitySMILES") or rec.get("CanonicalSMILES"))
        out = {
            "pubchem_cid": rec.get("CID"),
            "formula": rec.get("MolecularFormula"),
            "mw": rec.get("MolecularWeight"),
            "smiles": smiles,
            "inchikey": rec.get("InChIKey"),
            "iupac": rec.get("IUPACName"),
        }
        cache[name] = out
        return out
    except Exception as exc:                                    # noqa: BLE001
        print(f"    ! {name}: {exc}")
        cache[name] = None
        return None


# ---------------------------------------------------------------------------
# descriptors
# ---------------------------------------------------------------------------

def describe(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd = rdMolDescriptors.CalcNumHBD(m)
    hba = rdMolDescriptors.CalcNumHBA(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)

    lipinski_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_ok = (rotb <= 10) and (tpsa <= 140)

    return {
        "rdkit_mw": round(mw, 3),
        "logp": round(logp, 3),
        "hbd": hbd, "hba": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotb,
        "rings": rdMolDescriptors.CalcNumRings(m),
        "heavy_atoms": m.GetNumHeavyAtoms(),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
        "qed": round(QED.qed(m), 3),
        "lipinski_violations": lipinski_viol,
        "lipinski_pass": lipinski_viol <= 1,
        "veber_pass": veber_ok,
    }


def embed_3d(smiles: str, name: str):
    """3D conformer, MMFF94-minimised -- the docking-ready structure."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 20260729
    if AllChem.EmbedMolecule(m, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(m, ps) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    except Exception:                                           # noqa: BLE001
        pass
    m.SetProp("_Name", name)
    return m


def main() -> pd.DataFrame:
    cache = _load_cache()
    rows, mols, failed = [], [], []

    allc = ([(n, p, c, "moringa") for n, p, c in MORINGA_COMPOUNDS]
            + [(n, p, c, "reference") for n, p, c in REFERENCE_DRUGS])

    print(f"Resolving {len(allc)} compounds against PubChem "
          f"({len(cache)} already cached)...")
    for name, part, cls, src in allc:
        rec = fetch_pubchem(name, cache)
        if not rec or not rec.get("smiles"):
            failed.append(name)
            print(f"  MISS  {name}")
            continue
        d = describe(rec["smiles"])
        if not d:
            failed.append(name)
            print(f"  BADSMILES  {name}")
            continue
        rows.append({"name": name, "plant_part": part, "chem_class": cls,
                     "source": src, **rec, **d})
        m = embed_3d(rec["smiles"], name)
        if m is not None:
            mols.append(m)
        print(f"  ok    {name:32s} CID {str(rec['pubchem_cid']):>10s}  "
              f"MW {d['rdkit_mw']:7.2f}  logP {d['logp']:6.2f}  "
              f"QED {d['qed']:.2f}")

    _save_cache(cache)

    df = pd.DataFrame(rows)
    out_csv = PROC / "ligand_library.csv"
    df.to_csv(out_csv, index=False)

    sdf = STRUCT / "ligands.sdf"
    w = Chem.SDWriter(str(sdf))
    for m in mols:
        w.write(m)
    w.close()

    print("\n" + "=" * 74)
    print(f"Resolved   : {len(df)} / {len(allc)}")
    print(f"3D embedded: {len(mols)}")
    if failed:
        print(f"Failed     : {failed}")
    print(f"CSV        : {out_csv}")
    print(f"SDF        : {sdf}")

    if len(df):
        mo = df[df["source"] == "moringa"]
        print("\nDrug-likeness of the Moringa library:")
        print(f"  Lipinski pass (<=1 violation): "
              f"{int(mo['lipinski_pass'].sum())}/{len(mo)}")
        print(f"  Veber pass                   : "
              f"{int(mo['veber_pass'].sum())}/{len(mo)}")
        print(f"  median MW {mo['rdkit_mw'].median():.1f}, "
              f"median logP {mo['logp'].median():.2f}, "
              f"median QED {mo['qed'].median():.3f}")
        print("\nBy chemical class:")
        print(mo.groupby("chem_class").size().sort_values(ascending=False).to_string())
    return df


if __name__ == "__main__":
    main()
