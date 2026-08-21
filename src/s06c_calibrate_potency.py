"""
s06c_calibrate_potency.py  --  turn raw QSAR output into a usable potency table

The reality check (s06b) established three things:
  * Gate G2 FAILED on the scaffold split (test R^2 = 0.562 < 0.60)
  * on 5 Moringa constituents with measured A549 values, the typical error is
    ~3-fold and the model is systematically OVER-potent (mean log10 error
    -0.328, i.e. it calls compounds ~2.1x more potent than measured)
  * predictions are shrunk toward the training mean (SD ratio 0.370)

n = 5 is far too thin to calibrate on. This script estimates bias and spread
from the scaffold-split TEST SET instead (n ~ 1400), and specifically from the
subset of it that lies inside the physicochemical envelope of the Moringa
library -- i.e. the compounds most relevant to this study.

It then produces a calibrated potency table with an explicit per-compound
uncertainty, using this precedence:

    1. MEASURED A549 value from ChEMBL, where one exists   (uncertainty: small)
    2. bias-corrected QSAR prediction                      (uncertainty: envelope RMSE)
    3. flagged unusable if outside the applicability domain

The over-potency bias matters for the dose-range decision: correcting it moves
the predicted extract IC50 UPWARD, toward and past the 200 ppm ceiling. An
uncorrected model would have understated the risk in defect D2.

Output: results/tables/calibrated_potency.csv
        results/tables/calibration_report.txt
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from sklearn.metrics import r2_score

from s06_qsar_train import featurise, scaffold_split, SEED

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

# Properties used to define "chemically like the Moringa library".
ENVELOPE_COLS = ["rdkit_mw", "logp", "tpsa", "fraction_csp3"]
DESC_INDEX = {"MolWt": 0, "LogP": 1, "TPSA": 2, "FracCSP3": 8}


def main():
    out_lines: list[str] = []

    def say(s=""):
        print(s, flush=True)
        out_lines.append(s)

    df = pd.read_csv(PROC / "qsar_dataset_A549.csv").dropna(
        subset=["smiles", "pIC50"]).reset_index(drop=True)
    lib = pd.read_csv(PROC / "ligand_library.csv")
    lib = lib[lib["smiles"].notna()].reset_index(drop=True)
    lib["inchikey_calc"] = lib["smiles"].map(
        lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
        if Chem.MolFromSmiles(s) else None)
    pred = pd.read_csv(TAB / "qsar_moringa_predictions.csv")

    # reproduce exactly the training pool used in s06 (library held out)
    predict_keys = set(lib["inchikey_calc"].dropna())
    train_df = df[~df["inchikey"].isin(predict_keys)].reset_index(drop=True)

    X, keep = featurise(train_df["smiles"].tolist())
    y = train_df["pIC50"].to_numpy()[keep]
    smi = train_df["smiles"].to_numpy()[keep]

    tr, te = scaffold_split(smi)
    say("=" * 78)
    say("CALIBRATION FROM THE SCAFFOLD-SPLIT TEST SET")
    say("=" * 78)
    say(f"train {len(tr)}   test {len(te)}")

    from xgboost import XGBRegressor
    mdl = XGBRegressor(n_estimators=350, max_depth=6, learning_rate=0.08,
                       subsample=0.85, colsample_bytree=0.5, n_jobs=4,
                       random_state=SEED, tree_method="hist")
    mdl.fit(X[tr], y[tr])
    yhat = mdl.predict(X[te])
    resid = y[te] - yhat            # positive => model UNDER-estimates pIC50

    say(f"\nAll test compounds (n={len(te)}):")
    say(f"  R^2                        {r2_score(y[te], yhat):+.3f}")
    say(f"  bias  mean(meas - pred)    {resid.mean():+.3f} log units")
    say(f"  RMSE                       {np.sqrt((resid ** 2).mean()):.3f} log units")
    say(f"  |error| median             {np.median(np.abs(resid)):.3f} "
        f"({10 ** np.median(np.abs(resid)):.1f}-fold)")

    # --- restrict to the Moringa physicochemical envelope ---
    mo = lib[lib["source"] == "moringa"]
    lo = {c: float(mo[c].quantile(0.05)) for c in ENVELOPE_COLS}
    hi = {c: float(mo[c].quantile(0.95)) for c in ENVELOPE_COLS}
    say("\nMoringa physicochemical envelope (5th-95th pct):")
    for c in ENVELOPE_COLS:
        say(f"  {c:15s} {lo[c]:8.2f} .. {hi[c]:8.2f}")

    # descriptors sit in the tail of X, after NBITS fingerprint columns
    n_desc = len(DESC_INDEX) and X.shape[1]
    from s06_qsar_train import NBITS
    d_mw = X[te, NBITS + DESC_INDEX["MolWt"]]
    d_lp = X[te, NBITS + DESC_INDEX["LogP"]]
    d_tp = X[te, NBITS + DESC_INDEX["TPSA"]]
    d_cs = X[te, NBITS + DESC_INDEX["FracCSP3"]]

    inside = ((d_mw >= lo["rdkit_mw"]) & (d_mw <= hi["rdkit_mw"]) &
              (d_lp >= lo["logp"]) & (d_lp <= hi["logp"]) &
              (d_tp >= lo["tpsa"]) & (d_tp <= hi["tpsa"]) &
              (d_cs >= lo["fraction_csp3"]) & (d_cs <= hi["fraction_csp3"]))

    r_in = resid[inside]
    say(f"\nTest compounds INSIDE the Moringa envelope (n={int(inside.sum())}):")
    if inside.sum() < 30:
        say("  too few for a stable estimate -- falling back to all test compounds")
        bias, sigma = float(resid.mean()), float(np.sqrt((resid ** 2).mean()))
        basis = f"all scaffold-test compounds (n={len(te)})"
    else:
        say(f"  R^2                        "
            f"{r2_score(y[te][inside], yhat[inside]):+.3f}")
        say(f"  bias  mean(meas - pred)    {r_in.mean():+.3f} log units")
        say(f"  RMSE                       {np.sqrt((r_in ** 2).mean()):.3f} log units")
        say(f"  |error| median             {np.median(np.abs(r_in)):.3f} "
            f"({10 ** np.median(np.abs(r_in)):.1f}-fold)")
        bias, sigma = float(r_in.mean()), float(np.sqrt((r_in ** 2).mean()))
        basis = f"scaffold-test compounds inside the Moringa envelope (n={int(inside.sum())})"

    say(f"\nCALIBRATION ADOPTED (basis: {basis})")
    say(f"  pIC50_calibrated = pIC50_predicted + ({bias:+.3f})")
    say(f"  1 sigma uncertainty = {sigma:.3f} log units "
        f"(= x/{10 ** sigma:.1f} to x{10 ** sigma:.1f})")
    direction = ("UPWARD (model was over-potent)" if bias < 0
                 else "DOWNWARD (model was under-potent)")
    say(f"  net effect on IC50 estimates: {direction}")

    # corroboration from the 5 real measurements
    m = (pred.merge(lib[["name", "inchikey_calc"]], on="name", how="left")
             .merge(df[["inchikey", "pIC50"]], left_on="inchikey_calc",
                    right_on="inchikey", how="inner"))
    qsar_usable = True
    if len(m):
        np_bias = float((m["pIC50"] - m["pred_pIC50"]).mean())
        say(f"\nCorroboration -- {len(m)} measured Moringa constituents give "
            f"bias {np_bias:+.3f}")
        say(f"  envelope estimate {bias:+.3f}")
        sign_conflict = np.sign(np_bias) != np.sign(bias)
        if sign_conflict:
            say("  !! SIGN CONFLICT: the two estimates disagree on the DIRECTION")
            say("     of the bias. A correction derived from the envelope would")
            say("     push potency the opposite way to what the only real")
            say("     natural-product measurements indicate.")
        if sigma > 0.7:
            say(f"  !! 1 sigma = {sigma:.2f} log units (~{10 ** sigma:.0f}-fold). The")
            say("     uncertainty is an order of magnitude larger than the bias")
            say("     correction, so correcting the bias changes nothing material.")
        if sign_conflict or sigma > 0.7:
            qsar_usable = False

    say("")
    say("=" * 78)
    if qsar_usable:
        say("VERDICT: QSAR potency values may be used, with the band above.")
    else:
        say("VERDICT: QSAR IS NOT USABLE FOR PER-COMPOUND POTENCY.")
        say("=" * 78)
        say("""  Do NOT sum these values into an extract-level IC50. With ~11-fold
  1-sigma error per component and an unresolved bias direction, a bottom-up
  mixture calculation would produce a confident-looking number with no
  information in it.

  Use instead:
    * MEASURED values for the constituents that have them (below).
    * A LITERATURE PRIOR on crude M. oleifera extract IC50 vs A549 for the
      extract-level prediction and the dose-range recommendation.
    * The QSAR for RANKING/TRIAGE only -- deciding which constituents are
      worth docking -- never for absolute potency.""")
    say("=" * 78)

    # ---------------- build the calibrated table ----------------
    cal = pred.copy()
    measured = df[["inchikey", "pIC50", "ic50_ug_mL"]].rename(
        columns={"pIC50": "measured_pIC50", "ic50_ug_mL": "measured_IC50_ug_mL"})
    cal = (cal.merge(lib[["name", "inchikey_calc"]], on="name", how="left")
              .merge(measured, left_on="inchikey_calc", right_on="inchikey",
                     how="left"))

    cal["pIC50_calibrated"] = cal["pred_pIC50"] + bias
    cal["potency_source"] = "QSAR (bias-corrected)"
    cal["sigma_log"] = sigma

    has_meas = cal["measured_pIC50"].notna()
    cal.loc[has_meas, "pIC50_calibrated"] = cal.loc[has_meas, "measured_pIC50"]
    cal.loc[has_meas, "potency_source"] = "measured (ChEMBL)"
    cal.loc[has_meas, "sigma_log"] = 0.30      # inter-lab spread on a single value

    out_of_domain = ~cal["in_applicability_domain"].fillna(False)
    cal.loc[out_of_domain & ~has_meas, "potency_source"] = "UNUSABLE (out of domain)"

    cal["IC50_uM_calibrated"] = 10.0 ** (6.0 - cal["pIC50_calibrated"])
    cal["pred_IC50_ug_mL_calibrated"] = (cal["IC50_uM_calibrated"]
                                        * cal["rdkit_mw"] / 1000.0)
    # weak/strong ends of a 1-sigma band, in ug/mL
    cal["IC50_ug_mL_weak"] = cal["pred_IC50_ug_mL_calibrated"] * 10 ** cal["sigma_log"]
    cal["IC50_ug_mL_strong"] = cal["pred_IC50_ug_mL_calibrated"] / 10 ** cal["sigma_log"]

    keep_cols = ["name", "plant_part", "chem_class", "source", "pubchem_cid",
                 "rdkit_mw", "potency_source", "pred_pIC50", "pIC50_calibrated",
                 "sigma_log", "pred_IC50_ug_mL", "pred_IC50_ug_mL_calibrated",
                 "IC50_ug_mL_strong", "IC50_ug_mL_weak",
                 "max_tanimoto_to_train", "in_applicability_domain"]
    cal[keep_cols].to_csv(TAB / "calibrated_potency.csv", index=False)

    say("\n" + "=" * 78)
    say("CALIBRATED POTENCY -- Moringa constituents, most potent first")
    say("=" * 78)
    show = cal[cal["source"] == "moringa"].sort_values(
        "pred_IC50_ug_mL_calibrated")[
        ["name", "chem_class", "potency_source", "pred_IC50_ug_mL",
         "pred_IC50_ug_mL_calibrated", "IC50_ug_mL_strong", "IC50_ug_mL_weak"]]
    show = show.rename(columns={"pred_IC50_ug_mL": "raw",
                                "pred_IC50_ug_mL_calibrated": "calibrated",
                                "IC50_ug_mL_strong": "1sig_strong",
                                "IC50_ug_mL_weak": "1sig_weak"})
    say(show.head(20).to_string(index=False, float_format=lambda x: f"{x:9.2f}"))

    n_meas = int(has_meas.sum())
    n_unus = int((cal["potency_source"] == "UNUSABLE (out of domain)").sum())
    say(f"\n{n_meas} compounds use a measured value; "
        f"{n_unus} are unusable (out of domain).")
    say(f"\nWrote {TAB / 'calibrated_potency.csv'}")

    (TAB / "calibration_report.txt").write_text("\n".join(out_lines),
                                                encoding="utf-8")


if __name__ == "__main__":
    main()
