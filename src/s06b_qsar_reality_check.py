"""
s06b_qsar_reality_check.py  --  does the QSAR actually discriminate?

Gate G2 failed on the scaffold split (test R^2 = 0.56 < 0.60), so before any
potency number is used in the Phase C mixture model, two questions must be
answered honestly:

  Q1  How do predictions compare to REAL measured A549 values for the handful
      of Moringa constituents that already have them?  These 5 compounds were
      removed from training, so this is a genuine external test on exactly the
      chemical class the study cares about.

  Q2  Is the model discriminating, or is it regressing to the training mean?
      A model that outputs ~the median pIC50 for every natural product would
      still look plausible in a table, and would be worthless.  The test: the
      spread of predictions vs the spread of the training labels, and whether
      a mean-only predictor would score as well.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"


def main():
    pred = pd.read_csv(TAB / "qsar_moringa_predictions.csv")
    chembl = pd.read_csv(PROC / "qsar_dataset_A549.csv")
    lib = pd.read_csv(PROC / "ligand_library.csv")

    lib["inchikey_calc"] = lib["smiles"].map(
        lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
        if isinstance(s, str) and Chem.MolFromSmiles(s) else None)
    key = lib[["name", "inchikey_calc"]]

    m = (pred.merge(key, on="name", how="left")
              .merge(chembl[["inchikey", "pIC50", "ic50_ug_mL", "n_records"]],
                     left_on="inchikey_calc", right_on="inchikey", how="inner"))

    print("=" * 78)
    print("Q1 -- predicted vs MEASURED A549 potency (compounds held out of training)")
    print("=" * 78)
    if m.empty:
        print("  no overlap found")
    else:
        m["fold_error"] = m["pred_IC50_ug_mL"] / m["ic50_ug_mL"]
        m["log_err"] = np.log10(m["pred_IC50_ug_mL"]) - np.log10(m["ic50_ug_mL"])
        show = m[["name", "ic50_ug_mL", "pred_IC50_ug_mL", "fold_error",
                  "pIC50", "pred_pIC50", "max_tanimoto_to_train", "n_records"]]
        show = show.rename(columns={"ic50_ug_mL": "measured_ug_mL",
                                    "pred_IC50_ug_mL": "predicted_ug_mL",
                                    "pIC50": "measured_pIC50"})
        print(show.to_string(index=False, float_format=lambda x: f"{x:10.3f}"))

        mae_log = float(np.abs(m["log_err"]).mean())
        print(f"\n  mean |log10 error| = {mae_log:.3f}  "
              f"(= {10 ** mae_log:.1f}-fold typical error)")
        within3 = int((m["fold_error"].between(1 / 3, 3)).sum())
        print(f"  within 3-fold: {within3}/{len(m)}")
        print(f"  systematic bias (mean log10 error) = {m['log_err'].mean():+.3f} "
              f"-> model is {'UNDER' if m['log_err'].mean() < 0 else 'OVER'}"
              f"-predicting IC50 (i.e. calling them "
              f"{'MORE' if m['log_err'].mean() < 0 else 'LESS'} potent than measured)")

    print("\n" + "=" * 78)
    print("Q2 -- is the model discriminating, or predicting the mean?")
    print("=" * 78)

    tr_sd = float(chembl["pIC50"].std())
    tr_med = float(chembl["pIC50"].median())
    mo = pred[pred["source"] == "moringa"]
    pr_sd = float(mo["pred_pIC50"].std())

    print(f"  training pIC50:    median {tr_med:.2f}, SD {tr_sd:.3f}, "
          f"range {chembl['pIC50'].min():.2f}-{chembl['pIC50'].max():.2f}")
    print(f"  predicted pIC50 (Moringa): mean {mo['pred_pIC50'].mean():.2f}, "
          f"SD {pr_sd:.3f}, range {mo['pred_pIC50'].min():.2f}-"
          f"{mo['pred_pIC50'].max():.2f}")
    shrink = pr_sd / tr_sd
    print(f"\n  prediction SD / training SD = {shrink:.3f}")
    if shrink < 0.35:
        print("  -> SEVERE regression to the mean. The model is assigning nearly")
        print("     the same potency to every constituent. Differences between")
        print("     compounds in the output table are mostly NOISE plus molecular")
        print("     weight, NOT learned potency.")
    elif shrink < 0.6:
        print("  -> Substantial shrinkage. Rankings may be weakly informative;")
        print("     absolute values are not trustworthy.")
    else:
        print("  -> Reasonable spread retained.")

    # how far are the Moringa compounds from the training distribution?
    print(f"\n  applicability (max Tanimoto to training set):")
    print(f"    median {mo['max_tanimoto_to_train'].median():.3f}, "
          f"min {mo['max_tanimoto_to_train'].min():.3f}, "
          f"max {mo['max_tanimoto_to_train'].max():.3f}")
    print(f"    in domain (>=0.30): "
          f"{int(mo['in_applicability_domain'].sum())}/{len(mo)}")
    print("    For reference, doxorubicin scored 0.877 -- the Moringa")
    print("    constituents are far less similar to the training data than that.")

    print("\n" + "=" * 78)
    print("CONSEQUENCE FOR PHASE C")
    print("=" * 78)
    print("""  The mixture model must NOT be fed these per-compound IC50 values as if
  they were accurate. Required changes:

    1. Propagate the scaffold-split error (~1 log unit) into the predicted
       extract IC50 as an explicit uncertainty band, not a point estimate.
    2. Prefer MEASURED literature IC50 values wherever they exist (the
       compounds above, plus any others found by literature search); fall back
       to QSAR only for constituents with no measurement.
    3. Report the extract prediction as an order-of-magnitude bracket, and let
       the dose-range recommendation be driven by the WORST case (weakest
       predicted potency), since the costly failure mode is a dose range that
       is too low (defect D2).""")


if __name__ == "__main__":
    main()
