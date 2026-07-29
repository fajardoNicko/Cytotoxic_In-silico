"""
s06_qsar_train.py  --  Phase A, steps A4 + A6 of IN_SILICO_PLAN.md

Trains the A549 cytotoxicity QSAR model and runs the two gates that decide
whether the rest of the study is worth doing:

  Gate G2  model validity
           5-fold CV Q^2 >= 0.5, external test R^2 >= 0.6,
           y-randomisation Q^2_rand < 0.2

  Gate G4  positive-control anchor
           doxorubicin's predicted A549 IC50 must land within 3-fold of
           published values.  CRITICALLY, doxorubicin and every reference
           drug are REMOVED from the training set before this test -- other-
           wise the model is simply recalling a memorised label and the
           "anchor" proves nothing.

Two splits are reported:
  * random split      -- optimistic, the number most papers quote
  * scaffold split    -- honest: test compounds share no Bemis-Murcko
                         scaffold with training compounds, which is the
                         situation the Moringa constituents are actually in

Output: results/tables/qsar_model_card.csv
        results/tables/qsar_moringa_predictions.csv
        data/processed/qsar_model.pkl
"""

from __future__ import annotations

import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SEED = 20260729
# 1024 bits rather than 2048: on ~7k compounds the extra bits did not improve
# scaffold-split R^2 but roughly doubled the cost of every forest fit, and the
# validation protocol below needs dozens of fits.
NBITS = 1024
MORGAN = GetMorganGenerator(radius=2, fpSize=NBITS)

# Published doxorubicin IC50 against A549 spans roughly 0.1-10 uM depending on
# exposure time (24/48/72 h) and readout. The anchor is judged against this
# window, not a single literature point.
DOX_LIT_UM = (0.1, 10.0)


# ---------------------------------------------------------------------------
# featurisation
# ---------------------------------------------------------------------------

DESC_FUNCS = [
    ("MolWt", Descriptors.MolWt), ("LogP", Crippen.MolLogP),
    ("TPSA", rdMolDescriptors.CalcTPSA), ("HBD", rdMolDescriptors.CalcNumHBD),
    ("HBA", rdMolDescriptors.CalcNumHBA),
    ("RotB", rdMolDescriptors.CalcNumRotatableBonds),
    ("Rings", rdMolDescriptors.CalcNumRings),
    ("AromRings", rdMolDescriptors.CalcNumAromaticRings),
    ("FracCSP3", rdMolDescriptors.CalcFractionCSP3),
    ("HeavyAtoms", lambda m: m.GetNumHeavyAtoms()),
    ("MolMR", Crippen.MolMR), ("BalabanJ", Descriptors.BalabanJ),
    ("BertzCT", Descriptors.BertzCT), ("HallKierAlpha", Descriptors.HallKierAlpha),
    ("NumHeteroatoms", rdMolDescriptors.CalcNumHeteroatoms),
]


def featurise(smiles_list):
    fps, descs, keep = [], [], []
    for i, smi in enumerate(smiles_list):
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if m is None:
            continue
        fp = MORGAN.GetFingerprint(m)
        arr = np.zeros(NBITS, dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        row = []
        for _, f in DESC_FUNCS:
            try:
                v = float(f(m))
            except Exception:                                   # noqa: BLE001
                v = 0.0
            row.append(v if np.isfinite(v) else 0.0)
        fps.append(arr)
        descs.append(row)
        keep.append(i)
    X = np.hstack([np.array(fps, dtype=np.float32),
                   np.array(descs, dtype=np.float32)])
    return X, keep


def morgan_bits(smiles_list):
    out = []
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        out.append(MORGAN.GetFingerprint(m) if m else None)
    return out


# ---------------------------------------------------------------------------
# splits
# ---------------------------------------------------------------------------

def scaffold_split(smiles, test_frac=0.2, seed=SEED):
    """Bemis-Murcko scaffold split: no scaffold appears in both sets."""
    groups = defaultdict(list)
    for i, smi in enumerate(smiles):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        try:
            sc = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
        except Exception:                                       # noqa: BLE001
            sc = ""
        groups[sc].append(i)
    sets = sorted(groups.values(), key=len, reverse=True)
    rng = np.random.default_rng(seed)
    n_test = int(len(smiles) * test_frac)
    test, train = [], []
    for s in sets:                       # big scaffolds to train, singletons to test
        (test if len(test) < n_test else train).extend(s)
    rng.shuffle(test)
    rng.shuffle(train)
    return np.array(train), np.array(test)


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def metrics(y, yhat) -> dict:
    return {"R2": r2_score(y, yhat),
            "RMSE": float(np.sqrt(mean_squared_error(y, yhat))),
            "MAE": float(mean_absolute_error(y, yhat)),
            "n": len(y)}


def applicability(train_fps, query_fps):
    """Max Tanimoto similarity of each query to the training set."""
    out = []
    tf = [f for f in train_fps if f is not None]
    for q in query_fps:
        if q is None:
            out.append(np.nan)
            continue
        out.append(max(DataStructs.BulkTanimotoSimilarity(q, tf)))
    return np.array(out)


def main():
    ds = PROC / "qsar_dataset_A549.csv"
    if not ds.exists():
        raise SystemExit(f"missing {ds} -- run s05_chembl_fetch.py first")
    df = pd.read_csv(ds).dropna(subset=["smiles", "pIC50"]).reset_index(drop=True)
    print(f"Loaded {len(df)} curated A549 compounds")

    lib_path = PROC / "ligand_library.csv"
    lib = pd.read_csv(lib_path)
    lib = lib[lib["smiles"].notna()].reset_index(drop=True)
    lib["inchikey_calc"] = lib["smiles"].map(
        lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
        if Chem.MolFromSmiles(s) else None)

    # ---- leakage control: remove anything we intend to PREDICT ----
    predict_keys = set(lib["inchikey_calc"].dropna())
    overlap = df[df["inchikey"].isin(predict_keys)]
    print(f"\nLeakage control: {len(overlap)} of the {len(predict_keys)} library "
          f"compounds already have A549 labels in ChEMBL.")
    if len(overlap):
        merged = overlap.merge(lib, left_on="inchikey", right_on="inchikey_calc")
        for _, r in merged.iterrows():
            print(f"   held out: {r['name']:22s} measured pIC50 = {r['pIC50']:.2f} "
                  f"(IC50 {r['ic50_ug_mL']:.2f} ug/mL)")
    train_df = df[~df["inchikey"].isin(predict_keys)].reset_index(drop=True)
    print(f"Training pool after removing them: {len(train_df)} compounds")

    X, keep = featurise(train_df["smiles"].tolist())
    y = train_df["pIC50"].to_numpy()[keep]
    smi = train_df["smiles"].to_numpy()[keep]
    print(f"Feature matrix: {X.shape}")

    rows = []
    models = {
        "RandomForest": RandomForestRegressor(
            n_estimators=200, min_samples_leaf=3, max_features="sqrt",
            n_jobs=-1, random_state=SEED),
        "Ridge": Ridge(alpha=10.0, random_state=SEED),
    }
    try:
        from xgboost import XGBRegressor
        models["XGBoost"] = XGBRegressor(
            n_estimators=350, max_depth=6, learning_rate=0.08,
            subsample=0.85, colsample_bytree=0.5, n_jobs=-1,
            random_state=SEED, tree_method="hist")
    except ImportError:
        print("  (xgboost unavailable, skipping)", flush=True)

    idx = np.random.default_rng(SEED).permutation(len(y))
    splits = {
        "random": (idx[:int(.8 * len(idx))], idx[int(.8 * len(idx)):]),
        "scaffold": scaffold_split(smi),
    }

    for split_name, (tr, te) in splits.items():
        print(f"\n{'=' * 74}\nSPLIT: {split_name}  (train {len(tr)}, test {len(te)})",
              flush=True)
        for name, mdl in models.items():
            m = mdl.__class__(**mdl.get_params())
            m.fit(X[tr], y[tr])
            te_m = metrics(y[te], m.predict(X[te]))

            cv_pred = cross_val_predict(
                mdl.__class__(**mdl.get_params()), X[tr], y[tr],
                cv=KFold(5, shuffle=True, random_state=SEED), n_jobs=1)
            q2 = r2_score(y[tr], cv_pred)

            rows.append({"split": split_name, "model": name, "cv_Q2": q2,
                         "test_R2": te_m["R2"], "test_RMSE": te_m["RMSE"],
                         "test_MAE": te_m["MAE"], "y_rand_Q2": np.nan,
                         "n_train": len(tr), "n_test": len(te), "gate_G2": None})
            print(f"  {name:14s} Q2(cv)={q2:6.3f}  R2(test)={te_m['R2']:6.3f}  "
                  f"RMSE={te_m['RMSE']:.3f}  MAE={te_m['MAE']:.3f}", flush=True)

    # y-randomisation is a check on the FEATURES + protocol, not on each model
    # variant, so it is run once for the strongest model on the harder split.
    card0 = pd.DataFrame(rows)
    sc = card0[card0["split"] == "scaffold"].sort_values("test_R2")
    best_name = sc.iloc[-1]["model"]
    tr, te = splits["scaffold"]
    print(f"\ny-randomisation ({best_name}, scaffold split, 3-fold)...", flush=True)
    yr = np.random.default_rng(SEED).permutation(y[tr])
    yr_pred = cross_val_predict(
        models[best_name].__class__(**models[best_name].get_params()),
        X[tr], yr, cv=KFold(3, shuffle=True, random_state=SEED), n_jobs=1)
    q2_rand = r2_score(yr, yr_pred)
    print(f"  Q2(y-randomised) = {q2_rand:+.4f}  "
          f"(must be < 0.2; near 0 means the model is not fitting noise)",
          flush=True)

    card0["y_rand_Q2"] = q2_rand
    card0["gate_G2"] = ((card0["cv_Q2"] >= 0.5) & (card0["test_R2"] >= 0.6)
                        & (q2_rand < 0.2))
    rows = card0.to_dict("records")
    for r in rows:
        print(f"  G2 {r['split']:9s} {r['model']:14s}: "
              f"{'PASS' if r['gate_G2'] else 'FAIL'}", flush=True)

    card = pd.DataFrame(rows)
    card.to_csv(TAB / "qsar_model_card.csv", index=False)

    # ---- final model on everything, then predict the library ----
    best = card[card["split"] == "scaffold"].sort_values("test_R2").iloc[-1]
    print(f"\nBest by scaffold-split test R2: {best['model']}")
    final = models[best["model"]].__class__(**models[best["model"]].get_params())
    final.fit(X, y)

    with open(PROC / "qsar_model.pkl", "wb") as fh:
        pickle.dump({"model": final, "nbits": NBITS,
                     "desc_names": [n for n, _ in DESC_FUNCS], "seed": SEED}, fh)

    Xq, keepq = featurise(lib["smiles"].tolist())
    libq = lib.iloc[keepq].reset_index(drop=True)
    pred = final.predict(Xq)

    train_fps = morgan_bits(smi.tolist())
    sim = applicability(train_fps, morgan_bits(libq["smiles"].tolist()))

    libq["pred_pIC50"] = pred
    libq["pred_IC50_uM"] = 10.0 ** (6.0 - pred)
    libq["pred_IC50_ug_mL"] = libq["pred_IC50_uM"] * libq["rdkit_mw"] / 1000.0
    libq["max_tanimoto_to_train"] = sim
    libq["in_applicability_domain"] = sim >= 0.30

    out = libq[["name", "plant_part", "chem_class", "source", "pubchem_cid",
                "rdkit_mw", "pred_pIC50", "pred_IC50_uM", "pred_IC50_ug_mL",
                "max_tanimoto_to_train", "in_applicability_domain"]]
    out.to_csv(TAB / "qsar_moringa_predictions.csv", index=False)

    # ---- Gate G4: the doxorubicin anchor ----
    print(f"\n{'=' * 74}\nGATE G4 -- positive-control anchor (doxorubicin held out of training)")
    dox = out[out["name"] == "doxorubicin"]
    if len(dox):
        d = dox.iloc[0]
        um = d["pred_IC50_uM"]
        inside = DOX_LIT_UM[0] <= um <= DOX_LIT_UM[1]
        print(f"  predicted doxorubicin IC50 = {um:.3f} uM "
              f"({d['pred_IC50_ug_mL']:.3f} ug/mL)")
        print(f"  published A549 window      = {DOX_LIT_UM[0]}-{DOX_LIT_UM[1]} uM")
        print(f"  applicability (max Tanimoto to train) = "
              f"{d['max_tanimoto_to_train']:.3f}")
        print(f"  G4: {'PASS' if inside else 'FAIL'}")

    print(f"\n{'=' * 74}\nTop predicted Moringa constituents (in applicability domain):")
    mo = out[(out["source"] == "moringa") & out["in_applicability_domain"]]
    print(mo.sort_values("pred_IC50_ug_mL").head(15).to_string(
        index=False, float_format=lambda x: f"{x:9.3f}"))

    n_out = int((~out[out["source"] == "moringa"]["in_applicability_domain"]).sum())
    print(f"\n{n_out} Moringa compounds fall OUTSIDE the applicability domain "
          f"(max Tanimoto < 0.30); their predictions are flagged unreliable.")
    print(f"\nWrote {TAB / 'qsar_model_card.csv'}")
    print(f"Wrote {TAB / 'qsar_moringa_predictions.csv'}")


if __name__ == "__main__":
    main()
