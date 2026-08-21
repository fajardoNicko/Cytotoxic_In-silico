"""
s11_observed_data.py  --  Phase F, step 1 of IN_SILICO_PLAN.md

Ingests the real MTT absorbances returned by the MSU-IIT CCCBA laboratory and
rebuilds percent viability from raw optical density, rather than trusting the
pre-computed percentages in the report.

Source documents
  data/raw/mtt_observed_absorbance.csv   transcribed from the laboratory's
                                         raw-data spreadsheet, 6 wells per
                                         concentration (2 trials x 3 wells).
                                         Red-font cells in the source are
                                         recorded as lab_grubbs_outlier.
  data/raw/mtt_observed_reported.csv     Table 1 of the laboratory report,
                                         the percentages the lab itself
                                         calculated.

Why rebuild from raw OD
  The prediction was sealed against a quantity the lab computes with an
  undocumented step, described in the report only as "data were corrected
  prior to calculation". Recomputing from OD is the only way to know whether
  a predicted-versus-observed gap is biology or arithmetic. It also lets the
  outlier rule be audited, which matters here because Grubbs' test at n = 3
  is close to degenerate.

Outputs
  results/tables/F1_observed_viability.csv    per-well viability, all rules
  results/tables/F1_dose_summary.csv          per-dose mean and SD
  results/tables/F1_reconciliation.csv        our recomputation vs the lab's
  results/tables/F1_grubbs_audit.csv          our outlier flags vs the lab's
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05


# ---------------------------------------------------------------------------
# outlier rule
# ---------------------------------------------------------------------------

def grubbs_flags(x: np.ndarray, alpha: float = ALPHA) -> np.ndarray:
    """Two-sided Grubbs test, single outlier, applied once.

    Returns a boolean mask of flagged points. At n = 3 the attainable range
    of G is narrow, from 1.0000 for three equally spaced values to 1.1547
    when two coincide, while the critical value at alpha = 0.05 is 1.1543.
    The test therefore rejects almost nothing at three wells. This is a
    property of the test at that sample size, and it is what makes the
    laboratory's flag count worth auditing.
    """
    x = np.asarray(x, float)
    n = x.size
    out = np.zeros(n, bool)
    if n < 3:
        return out
    sd = x.std(ddof=1)
    if sd == 0:
        return out
    g = np.abs(x - x.mean()) / sd
    i = int(np.argmax(g))
    t = sps.t.ppf(1 - alpha / (2 * n), n - 2)
    g_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t ** 2 / (n - 2 + t ** 2))
    if g[i] > g_crit:
        out[i] = True
    return out


def grubbs_critical(n: int, alpha: float = ALPHA) -> float:
    t = sps.t.ppf(1 - alpha / (2 * n), n - 2)
    return ((n - 1) / np.sqrt(n)) * np.sqrt(t ** 2 / (n - 2 + t ** 2))


# ---------------------------------------------------------------------------
# viability
# ---------------------------------------------------------------------------

def viability(abs_sample, abs_blank, abs_control):
    """The laboratory's formula, identical to the manuscript's."""
    return (np.asarray(abs_sample, float) - abs_blank) / (abs_control - abs_blank) * 100.0


def per_trial_refs(raw: pd.DataFrame, trial: int, drop_outliers: bool) -> tuple[float, float]:
    """Blank and negative-control means for one trial."""
    d = raw[raw["trial"] == trial]
    if drop_outliers:
        d = d[~d["lab_grubbs_outlier"]]
    blank = d.loc[d["sample"] == "BLANK", "abs_570"].mean()
    ctrl = d.loc[d["sample"] == "NC", "abs_570"].mean()
    return float(blank), float(ctrl)


def build(raw: pd.DataFrame, drop_outliers: bool, label: str) -> pd.DataFrame:
    rows = []
    for trial in sorted(raw["trial"].unique()):
        blank, ctrl = per_trial_refs(raw, trial, drop_outliers)
        d = raw[(raw["trial"] == trial) & (~raw["sample"].isin(["BLANK"]))]
        if drop_outliers:
            d = d[~d["lab_grubbs_outlier"]]
        for _, w in d.iterrows():
            rows.append({
                "rule": label,
                "sample": w["sample"],
                "conc_ppm": w["conc_ppm"],
                "trial": trial,
                "replicate": w["replicate"],
                "abs_570": w["abs_570"],
                "abs_blank_used": blank,
                "abs_control_used": ctrl,
                "viability_pct": float(viability(w["abs_570"], blank, ctrl)),
            })
    out = pd.DataFrame(rows)
    out["cytotoxicity_pct"] = 100.0 - out["viability_pct"]
    return out


def dose_summary(viab: pd.DataFrame) -> pd.DataFrame:
    """Per-dose mean and SD, and the trial-mean SD the lab reports."""
    rows = []
    for (rule, sample, conc), g in viab.groupby(["rule", "sample", "conc_ppm"], sort=False):
        trial_means = g.groupby("trial")["viability_pct"].mean()
        rows.append({
            "rule": rule, "sample": sample, "conc_ppm": conc,
            "n_wells": len(g),
            "mean_viability": g["viability_pct"].mean(),
            "sd_wells": g["viability_pct"].std(ddof=1),
            "mean_of_trial_means": trial_means.mean(),
            "sd_of_trial_means": trial_means.std(ddof=1) if trial_means.size > 1 else np.nan,
            "mean_inhibition": 100.0 - g["viability_pct"].mean(),
        })
    return pd.DataFrame(rows).sort_values(["rule", "sample", "conc_ppm"]).reset_index(drop=True)


# ---------------------------------------------------------------------------

def main() -> None:
    raw = pd.read_csv(RAW / "mtt_observed_absorbance.csv")
    raw["lab_grubbs_outlier"] = raw["lab_grubbs_outlier"].astype(str).str.upper().eq("TRUE")
    rep = pd.read_csv(RAW / "mtt_observed_reported.csv")

    print("=" * 78)
    print("PHASE F STEP 1 -- observed MTT data, rebuilt from raw absorbance")
    print("=" * 78)
    print(f"\nWells transcribed        : {len(raw)}")
    print(f"Flagged by the lab       : {int(raw['lab_grubbs_outlier'].sum())} "
          f"({100*raw['lab_grubbs_outlier'].mean():.1f}%)")

    # --- outlier audit -----------------------------------------------------
    audit = []
    for (sample, conc), g in raw.groupby(["sample", "conc_ppm"], sort=False):
        for trial, gt in g.groupby("trial"):
            mine = grubbs_flags(gt["abs_570"].to_numpy())
            audit.append({
                "sample": sample, "conc_ppm": conc, "trial": trial,
                "n": len(gt),
                "lab_flagged": int(gt["lab_grubbs_outlier"].sum()),
                "grubbs_n3_flagged": int(mine.sum()),
                "values": ", ".join(f"{v:.3f}" for v in gt["abs_570"]),
            })
        pooled = grubbs_flags(g["abs_570"].to_numpy())
        audit.append({
            "sample": sample, "conc_ppm": conc, "trial": "pooled(6)",
            "n": len(g),
            "lab_flagged": int(g["lab_grubbs_outlier"].sum()),
            "grubbs_n3_flagged": int(pooled.sum()),
            "values": ", ".join(f"{v:.3f}" for v in g["abs_570"]),
        })
    audit = pd.DataFrame(audit)
    audit.to_csv(TAB / "F1_grubbs_audit.csv", index=False)

    # A one-standard-deviation rule on the pooled six wells, for comparison.
    raw["sd1_rule"] = False
    for (_s, _c), g in raw.groupby(["sample", "conc_ppm"], sort=False):
        v = g["abs_570"]
        raw.loc[g.index, "sd1_rule"] = (np.abs(v - v.mean()) > v.std(ddof=1)).values

    n3 = audit[audit["trial"] != "pooled(6)"]
    pooled6 = audit[audit["trial"] == "pooled(6)"]
    print(f"\nOutlier rule audit")
    print(f"  Grubbs critical G, n = 3, alpha = 0.05 : {grubbs_critical(3):.4f}")
    print(f"  Attainable G at n = 3                  : 1.0000 to {2/np.sqrt(3):.4f}")
    print(f"  Grubbs critical G, n = 6, alpha = 0.05 : {grubbs_critical(6):.4f}")
    print("  At three wells the attainable range of G barely reaches the")
    print("  critical value, so Grubbs rejects almost nothing at this design.")
    print(f"  Grubbs per trial, n = 3   : {int(n3['grubbs_n3_flagged'].sum()):3d} wells flagged")
    print(f"  Grubbs pooled, n = 6      : {int(pooled6['grubbs_n3_flagged'].sum()):3d} wells flagged")
    print(f"  Beyond 1 SD, pooled n = 6 : {int(raw['sd1_rule'].sum()):3d} wells flagged")
    print(f"  Laboratory                : {int(raw['lab_grubbs_outlier'].sum()):3d} wells flagged")
    agree = float((raw["lab_grubbs_outlier"] == raw["sd1_rule"]).mean())
    print(f"  The lab's exclusions are far more aggressive than a Grubbs test")
    print(f"  at either grouping. A plain one-standard-deviation rule on the")
    print(f"  pooled six wells gives a similar count and agrees on {100*agree:.0f}% of")
    print(f"  wells, but not on all of them, so the exact rule used could not")
    print(f"  be identified from the report. Roughly a third of the data was")
    print(f"  excluded before the reported percentages were computed.")

    # --- viability, three rules -------------------------------------------
    keep = build(raw, drop_outliers=False, label="all_wells")
    drop = build(raw, drop_outliers=True, label="lab_outliers_removed")
    viab = pd.concat([keep, drop], ignore_index=True)
    viab.to_csv(TAB / "F1_observed_viability.csv", index=False)

    summ = dose_summary(viab)
    summ.to_csv(TAB / "F1_dose_summary.csv", index=False)

    # --- reconciliation against the lab's own percentages ------------------
    recon = []
    for _, r in rep[rep["sample"].isin(["MBEE", "DOX"])].iterrows():
        for rule in ("all_wells", "lab_outliers_removed"):
            g = viab[(viab["rule"] == rule) & (viab["sample"] == r["sample"])
                     & (viab["conc_ppm"] == r["conc_ppm"])]
            if g.empty:
                continue
            tm = g.groupby("trial")["viability_pct"].mean()
            ours = float(tm.mean())
            recon.append({
                "sample": r["sample"], "conc_ppm": r["conc_ppm"], "rule": rule,
                "lab_reported": r["viab_mean"],
                "recomputed": ours,
                "difference_pp": ours - r["viab_mean"],
            })
    recon = pd.DataFrame(recon)
    recon.to_csv(TAB / "F1_reconciliation.csv", index=False)

    best = recon.groupby("rule")["difference_pp"].apply(lambda s: np.sqrt((s ** 2).mean()))
    print("\nReconciliation with the laboratory's reported percentages")
    print("  RMS difference in percentage points, by outlier rule:")
    for rule, v in best.items():
        print(f"    {rule:24s} {v:6.2f}")
    print("  Neither rule reproduces the report exactly. The lab states that")
    print("  'data were corrected prior to calculation' without defining the")
    print("  correction, so an exact match was not reachable. The laboratory's")
    print("  own percentages are therefore carried as primary for validation,")
    print("  and this recomputation is reported alongside as an audit.")

    mb = summ[(summ["rule"] == "lab_outliers_removed") & (summ["sample"] == "MBEE")]
    print("\nMBEE dose response, outliers removed, rebuilt from OD:")
    print(mb[["conc_ppm", "n_wells", "mean_viability", "sd_wells"]]
          .to_string(index=False, float_format=lambda x: f"{x:9.2f}"))

    print(f"\nWrote F1_observed_viability.csv, F1_dose_summary.csv, "
          f"F1_reconciliation.csv, F1_grubbs_audit.csv -> {TAB}")


if __name__ == "__main__":
    main()
