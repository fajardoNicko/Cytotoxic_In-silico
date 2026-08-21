"""
s10_validation.py  --  Phase F of IN_SILICO_PLAN.md

Scores the sealed pre-registration against the observed MTT result.

The prediction in results/prediction_registry/prediction_v1.json was frozen
and hashed on 2026-08-20, before any absorbance existed. Nothing in this
script may alter that file. Every comparison below is against the numbers as
sealed.

Checks, one per sealed prediction:
  P1  extract IC50, fold error against the observed value        (Gate G5)
  P2  percent viability at each planned dose, RMSE and MAE       (Gate G5)
  P3  the 12.5 to 200 ppm series fails to bracket 50% inhibition
  P4  dose rank order, Spearman rho                              (Gate G5)
  P5  Tukey pattern against the negative control
  P8  log-linear IC50 lands outside the tested range

Also computes Lin's concordance correlation coefficient and a Bland-Altman
analysis, refits both IC50 estimators on the full observed series, and runs
the paper's statistical chain on the real data through the same stats mirror
used for the simulation.

Outputs
  results/tables/F2_validation_metrics.csv
  results/tables/F2_prediction_scorecard.csv
  results/tables/F3_statistics.txt
  results/figures/F1_predicted_vs_observed.png
  results/figures/F2_bland_altman.png
  results/figures/F3_ic50_regression.png
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import stats_mirror as sm
from mtt_model import loglinear_ic50, fourpl_ic50, nci_activity_class

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
REG = ROOT / "results" / "prediction_registry"
for p in (TAB, FIG):
    p.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
PLANNED_DOSES = [12.5, 25.0, 50.0, 100.0, 200.0]

# The laboratory's own log-linear solution, Table 2 of the report.
LAB_IC50 = 911.84
LAB_SLOPE, LAB_INTERCEPT, LAB_R2 = -42.437, 175.61, 0.9996
LAB_FIT_DOSES = [400.0, 100.0, 50.0]
LAB_DOX_IC50 = 28.59
TOP_DOSE = 400.0


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def lins_ccc(x, y):
    """Lin's concordance correlation coefficient."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    sx, sy = x.var(ddof=0), y.var(ddof=0)
    sxy = ((x - x.mean()) * (y - y.mean())).mean()
    return float(2 * sxy / (sx + sy + (x.mean() - y.mean()) ** 2))


def fold_error(pred, obs):
    return float(pred / obs)


def within_fold(pred, obs, k=3.0):
    fe = fold_error(pred, obs)
    return bool(1.0 / k <= fe <= k)


# ---------------------------------------------------------------------------

def main() -> None:
    pred = json.loads((REG / "prediction_v1.json").read_text(encoding="utf-8"))
    rep = pd.read_csv(RAW / "mtt_observed_reported.csv")
    viab = pd.read_csv(TAB / "F1_observed_viability.csv")

    mbee = rep[rep["sample"] == "MBEE"].sort_values("conc_ppm").reset_index(drop=True)
    obs_doses = mbee["conc_ppm"].to_numpy(float)
    obs_viab = mbee["viab_mean"].to_numpy(float)

    print("=" * 78)
    print("PHASE F -- sealed prediction versus observed MTT result")
    print("=" * 78)
    print(f"\nSealed   : {pred['sealed_utc']}  commit {pred['git_commit'][:12]}")
    print(f"Basis    : {pred['basis']['method']}")
    print(f"Observed : MSU-IIT CCCBA laboratory, released August 2026")

    rows, score = [], []

    # ---- P1 IC50 -----------------------------------------------------------
    p_sm = pred["PREDICTION_1_extract_ic50_ug_per_mL"]["point_solvent_matched"]
    p_pool = pred["PREDICTION_1_extract_ic50_ug_per_mL"]["point_pooled_optimistic"]
    lo68, hi68 = pred["PREDICTION_1_extract_ic50_ug_per_mL"]["interval_68pct_solvent_matched"]
    lo95, hi95 = pred["PREDICTION_1_extract_ic50_ug_per_mL"]["interval_95pct_solvent_matched"]

    fe_sm = fold_error(p_sm, LAB_IC50)
    fe_pool = fold_error(p_pool, LAB_IC50)
    p1_pass = within_fold(p_sm, LAB_IC50)

    print("\n" + "-" * 78)
    print("P1  EXTRACT IC50")
    print("-" * 78)
    print(f"  predicted, solvent-matched : {p_sm:10.1f} ug/mL")
    print(f"  68% interval               : {lo68:.1f} to {hi68:.1f}")
    print(f"  95% interval               : {lo95:.1f} to {hi95:.1f}")
    print(f"  observed, lab log-linear   : {LAB_IC50:10.2f} ug/mL")
    print(f"  fold error                 : {fe_sm:10.3f}   ({'PASS' if p1_pass else 'FAIL'} at 3-fold)")
    print(f"  observed inside 68% band   : {lo68 <= LAB_IC50 <= hi68}")
    print(f"  pooled prior for reference : {p_pool:.1f} ug/mL, fold error {fe_pool:.3f}")
    print(f"  NCI class, predicted       : {nci_activity_class(p_sm)}")
    print(f"  NCI class, observed        : {nci_activity_class(LAB_IC50)}")

    rows += [
        {"metric": "predicted_ic50_solvent_matched_ug_ml", "value": p_sm},
        {"metric": "predicted_ic50_pooled_ug_ml", "value": p_pool},
        {"metric": "observed_ic50_lab_loglinear_ug_ml", "value": LAB_IC50},
        {"metric": "ic50_fold_error_solvent_matched", "value": fe_sm},
        {"metric": "ic50_fold_error_pooled", "value": fe_pool},
        {"metric": "observed_within_68pct_interval", "value": int(lo68 <= LAB_IC50 <= hi68)},
        {"metric": "observed_within_95pct_interval", "value": int(lo95 <= LAB_IC50 <= hi95)},
    ]
    score.append({
        "prediction": "P1 extract IC50",
        "sealed": f"{p_sm:.1f} ug/mL (68% CI {lo68:.0f} to {hi68:.0f})",
        "observed": f"{LAB_IC50:.2f} ug/mL",
        "criterion": "within 3-fold",
        "result": f"fold error {fe_sm:.3f}",
        "verdict": "PASS" if p1_pass else "FAIL",
    })

    # ---- P2 viability ------------------------------------------------------
    pv = pred["PREDICTION_2_viability_at_planned_doses"]["solvent_matched_prior"]
    pv_pool = pred["PREDICTION_2_viability_at_planned_doses"]["pooled_prior"]
    pred_v = np.array([pv[str(d)] for d in PLANNED_DOSES], float)
    pred_vp = np.array([pv_pool[str(d)] for d in PLANNED_DOSES], float)
    obs_v = np.array([float(mbee.loc[mbee["conc_ppm"] == d, "viab_mean"].iloc[0])
                      for d in PLANNED_DOSES])

    resid = pred_v - obs_v
    rmse = float(np.sqrt((resid ** 2).mean()))
    mae = float(np.abs(resid).mean())
    bias = float(resid.mean())
    ccc = lins_ccc(pred_v, obs_v)
    rmse_pool = float(np.sqrt(((pred_vp - obs_v) ** 2).mean()))
    p2_pass = rmse <= pred["VALIDATION_CRITERIA"]["viability_rmse_threshold_pct_points"]

    print("\n" + "-" * 78)
    print("P2  PERCENT VIABILITY AT THE PLANNED DOSES")
    print("-" * 78)
    print(f"  {'dose':>8} {'predicted':>11} {'observed':>10} {'residual':>10}")
    for d, p_, o_, r_ in zip(PLANNED_DOSES, pred_v, obs_v, resid):
        print(f"  {d:8.1f} {p_:11.1f} {o_:10.2f} {r_:+10.2f}")
    print(f"\n  RMSE  {rmse:6.2f} percentage points   ({'PASS' if p2_pass else 'FAIL'} at 15)")
    print(f"  MAE   {mae:6.2f}   mean signed bias {bias:+6.2f}")
    print(f"  Lin's CCC {ccc:.4f}")
    print(f"  RMSE under the pooled prior, for reference: {rmse_pool:.2f}")

    rows += [
        {"metric": "viability_rmse_pp", "value": rmse},
        {"metric": "viability_mae_pp", "value": mae},
        {"metric": "viability_mean_bias_pp", "value": bias},
        {"metric": "lins_ccc", "value": ccc},
        {"metric": "viability_rmse_pooled_prior_pp", "value": rmse_pool},
    ]
    score.append({
        "prediction": "P2 viability at 5 doses",
        "sealed": ", ".join(f"{v:.1f}" for v in pred_v),
        "observed": ", ".join(f"{v:.1f}" for v in obs_v),
        "criterion": "RMSE at or below 15 pp",
        "result": f"RMSE {rmse:.2f} pp",
        "verdict": "PASS" if p2_pass else "FAIL",
    })

    # ---- P3 dose range -----------------------------------------------------
    v200 = float(mbee.loc[mbee["conc_ppm"] == 200, "viab_mean"].iloc[0])
    v400 = float(mbee.loc[mbee["conc_ppm"] == 400, "viab_mean"].iloc[0])
    brackets = v200 < 50.0
    p3_pass = not brackets
    max_inhib = 100.0 - obs_viab.min()

    print("\n" + "-" * 78)
    print("P3  DOSE RANGE ADEQUACY")
    print("-" * 78)
    print(f"  sealed claim : the 12.5 to 200 ppm series will not bracket 50% inhibition")
    print(f"  falsified if : viability at 200 ppm is below 50%")
    print(f"  observed     : {v200:.2f}% viability at 200 ppm, {v400:.2f}% at 400 ppm")
    print(f"  the lab added 400 ppm and still did not reach 50%. Maximum")
    print(f"  inhibition anywhere in the series was {max_inhib:.1f}%.")
    print(f"  verdict      : {'CONFIRMED' if p3_pass else 'FALSIFIED'}")

    rows += [
        {"metric": "observed_viability_at_200ppm", "value": v200},
        {"metric": "observed_viability_at_400ppm", "value": v400},
        {"metric": "max_observed_inhibition_pct", "value": max_inhib},
        {"metric": "series_brackets_ic50", "value": int(brackets)},
    ]
    score.append({
        "prediction": "P3 dose range fails",
        "sealed": "series will not bracket 50% inhibition",
        "observed": f"{v200:.1f}% viability at 200 ppm, {v400:.1f}% at 400 ppm",
        "criterion": "falsified if viability at 200 ppm below 50%",
        "result": f"max inhibition {max_inhib:.1f}%",
        "verdict": "CONFIRMED" if p3_pass else "FALSIFIED",
    })

    # ---- P4 rank order -----------------------------------------------------
    rho, rho_p = sps.spearmanr(obs_doses, obs_viab)
    p4_pass = abs(rho) >= pred["VALIDATION_CRITERIA"]["spearman_rho_threshold"]

    print("\n" + "-" * 78)
    print("P4  DOSE RANK ORDER")
    print("-" * 78)
    print(f"  sealed   : monotonic decrease, expected Spearman rho = -1.0")
    print(f"  observed : rho = {rho:+.4f}, p = {rho_p:.2e} over {len(obs_doses)} doses")
    print(f"  verdict  : {'PASS' if p4_pass else 'FAIL'}")

    rows += [{"metric": "spearman_rho_dose_vs_viability", "value": float(rho)},
             {"metric": "spearman_p", "value": float(rho_p)}]
    score.append({
        "prediction": "P4 dose rank order",
        "sealed": "monotonic decrease, rho = -1.0",
        "observed": f"rho = {rho:+.4f}",
        "criterion": "absolute rho at or above 0.9",
        "result": f"rho = {rho:+.4f}, p = {rho_p:.1e}",
        "verdict": "PASS" if p4_pass else "FAIL",
    })

    # ---- P8 estimator behaviour -------------------------------------------
    ours_ll = loglinear_ic50(obs_doses, obs_viab)
    try:
        ours_4pl = fourpl_ic50(obs_doses, obs_viab)
        fourpl_val = ours_4pl.absolute_ic50
        fourpl_r2 = ours_4pl.r_squared
    except Exception:
        fourpl_val, fourpl_r2 = float("nan"), float("nan")

    extrapolated = LAB_IC50 > TOP_DOSE
    ratio = LAB_IC50 / TOP_DOSE

    print("\n" + "-" * 78)
    print("P8  ESTIMATOR BEHAVIOUR OUTSIDE THE TESTED RANGE")
    print("-" * 78)
    print(f"  lab log-linear, 3 points {LAB_FIT_DOSES} : {LAB_IC50:9.2f} ug/mL, R2 = {LAB_R2}")
    print(f"  our log-linear, all 6 doses              : {ours_ll.ic50:9.2f} ug/mL, R2 = {ours_ll.r_squared:.4f}")
    print(f"  our 4PL, all 6 doses                     : {fourpl_val:9.2f} ug/mL")
    print(f"  top dose actually tested                 : {TOP_DOSE:9.2f} ug/mL")
    print(f"  the reported IC50 sits {ratio:.2f}x beyond the highest dose tested,")
    print(f"  so it is an extrapolation, not a determination, while the")
    print(f"  regression still reports R2 = {LAB_R2}. This is the failure mode")
    print(f"  sealed as P8: {'CONFIRMED' if extrapolated else 'not observed'}")

    rows += [
        {"metric": "our_loglinear_ic50_all_doses", "value": float(ours_ll.ic50)},
        {"metric": "our_loglinear_r2", "value": float(ours_ll.r_squared)},
        {"metric": "our_4pl_absolute_ic50", "value": float(fourpl_val)},
        {"metric": "lab_ic50_over_top_dose_ratio", "value": float(ratio)},
        {"metric": "lab_ic50_is_extrapolation", "value": int(extrapolated)},
        {"metric": "observed_doxorubicin_ic50_ug_ml", "value": LAB_DOX_IC50},
    ]
    score.append({
        "prediction": "P8 estimator leaves range",
        "sealed": "log-linear IC50 will fall outside the tested range with high R2",
        "observed": f"IC50 {LAB_IC50:.2f} vs top dose {TOP_DOSE:.0f}, R2 = {LAB_R2}",
        "criterion": "IC50 above the highest dose tested",
        "result": f"{ratio:.2f}x beyond the top dose",
        "verdict": "CONFIRMED" if extrapolated else "NOT OBSERVED",
    })

    # ---- P5 statistics on the real data -----------------------------------
    print("\n" + "-" * 78)
    print("P5  THE PAPER'S STATISTICAL CHAIN ON THE REAL DATA")
    print("-" * 78)

    d = viab[(viab["rule"] == "lab_outliers_removed")
             & (viab["sample"].isin(["MBEE", "NC"]))].copy()
    d["group"] = np.where(d["sample"] == "NC", "negative_control",
                          d["conc_ppm"].map(lambda c: f"{c:g}ppm"))
    report = sm.run_analysis(d, "group", "viability_pct")
    text = sm.format_report(report, "OBSERVED MTT DATA, MBEE VERSUS NEGATIVE CONTROL")
    (TAB / "F3_statistics.txt").write_text(text, encoding="utf-8")
    print(text[:2200])

    ph = report.posthoc.table if report.posthoc is not None else pd.DataFrame()
    sig_vs_ctrl = {}
    if not ph.empty:
        cols = ph.columns.tolist()
        gcol = [c for c in cols if "group" in c.lower()]
        pcol = [c for c in cols if c.lower().startswith("p")]
        if len(gcol) >= 2 and pcol:
            for _, r in ph.iterrows():
                pair = {str(r[gcol[0]]), str(r[gcol[1]])}
                if "negative_control" in pair:
                    other = (pair - {"negative_control"}).pop()
                    sig_vs_ctrl[other] = bool(float(r[pcol[0]]) < ALPHA)

    p5_sealed = pred["PREDICTION_5_statistics_under_planned_design"]["predicted_tukey_vs_control_significant"]
    print("\n  Sealed Tukey pattern versus what the data gave:")
    p5_hits, p5_total = 0, 0
    for dose_key, claim in p5_sealed.items():
        got = sig_vs_ctrl.get(dose_key)
        if isinstance(claim, bool) and got is not None:
            p5_total += 1
            ok = (claim == got)
            p5_hits += ok
            print(f"    {dose_key:>8}  sealed {str(claim):>5}   observed {str(got):>5}   {'match' if ok else 'MISS'}")
        else:
            print(f"    {dose_key:>8}  sealed {str(claim):>9}   observed {str(got):>5}")

    rows += [
        {"metric": "anova_p", "value": float(report.anova.p_value) if report.anova else np.nan},
        {"metric": "anova_F", "value": float(report.anova.f_stat) if report.anova else np.nan},
        {"metric": "tukey_pattern_matches", "value": p5_hits},
        {"metric": "tukey_pattern_checked", "value": p5_total},
    ]
    score.append({
        "prediction": "P5 Tukey pattern",
        "sealed": "no dose from 12.5 to 100 ppm separates from control, 200 ppm borderline",
        "observed": ", ".join(f"{k}={v}" for k, v in sig_vs_ctrl.items()) or "post hoc not reached",
        "criterion": "match on the doses with a definite sealed claim",
        "result": f"{p5_hits} of {p5_total} matched",
        "verdict": "PASS" if p5_total and p5_hits == p5_total else ("PARTIAL" if p5_hits else "FAIL"),
    })

    # ---- Bland-Altman ------------------------------------------------------
    mean_pv = (pred_v + obs_v) / 2
    diff_pv = pred_v - obs_v
    ba_bias = float(diff_pv.mean())
    ba_sd = float(diff_pv.std(ddof=1))
    loa = (ba_bias - 1.96 * ba_sd, ba_bias + 1.96 * ba_sd)
    rows += [{"metric": "bland_altman_bias_pp", "value": ba_bias},
             {"metric": "bland_altman_loa_lower", "value": loa[0]},
             {"metric": "bland_altman_loa_upper", "value": loa[1]}]

    # ---- overall gate ------------------------------------------------------
    g5 = p1_pass and p2_pass and p4_pass
    print("\n" + "=" * 78)
    print(f"GATE G5  {'PASSED' if g5 else 'FAILED'}")
    print("=" * 78)
    print(f"  IC50 within 3-fold        : {p1_pass}  (fold error {fe_sm:.3f})")
    print(f"  viability RMSE <= 15 pp   : {p2_pass}  ({rmse:.2f} pp)")
    print(f"  Spearman rho >= 0.9       : {p4_pass}  ({abs(rho):.4f})")
    print(f"  dose-range claim P3       : {'CONFIRMED' if p3_pass else 'FALSIFIED'}")
    print(f"  estimator claim P8        : {'CONFIRMED' if extrapolated else 'NOT OBSERVED'}")
    rows.append({"metric": "gate_G5_passed", "value": int(g5)})

    pd.DataFrame(rows).to_csv(TAB / "F2_validation_metrics.csv", index=False)
    pd.DataFrame(score).to_csv(TAB / "F2_prediction_scorecard.csv", index=False)

    # ---- figures -----------------------------------------------------------
    make_figures(obs_doses, obs_viab, pred_v, pred_vp, obs_v, mean_pv, diff_pv,
                 ba_bias, loa, ours_ll, rep)

    print(f"\nWrote F2_validation_metrics.csv, F2_prediction_scorecard.csv, "
          f"F3_statistics.txt -> {TAB}")
    print(f"Wrote F1, F2, F3 figures -> {FIG}")


def make_figures(obs_doses, obs_viab, pred_v, pred_vp, obs_v,
                 mean_pv, diff_pv, ba_bias, loa, ours_ll, rep):
    # F1 predicted versus observed dose-response
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.axhline(50, ls=":", c="grey", lw=1.2)
    ax.axhline(100, ls="-", c="lightgrey", lw=1)
    ax.plot(PLANNED_DOSES, pred_v, "s--", color="#2166ac", lw=1.8, ms=7,
            label="predicted, solvent-matched prior (sealed)")
    ax.plot(PLANNED_DOSES, pred_vp, "^--", color="#67a9cf", lw=1.4, ms=6,
            label="predicted, pooled prior (sealed)")
    ax.errorbar(obs_doses, obs_viab,
                yerr=rep[rep["sample"] == "MBEE"].sort_values("conc_ppm")["viab_sd"],
                fmt="o-", color="#b2182b", lw=2, ms=8, capsize=4,
                label="observed, MSU-IIT MTT assay")
    ax.set_xscale("log")
    ax.set_xticks(list(obs_doses))
    ax.set_xticklabels([f"{d:g}" for d in obs_doses])
    ax.minorticks_off()
    ax.set_xlabel("concentration (µg/mL)")
    ax.set_ylabel("% cell viability")
    ax.set_title("Sealed prediction versus observed dose response, A549")
    ax.annotate("50% inhibition never reached", xy=(15, 52), fontsize=8, color="grey")
    ax.legend(fontsize=8.5, loc="lower left")
    ax.set_ylim(0, 150)
    fig.tight_layout()
    fig.savefig(FIG / "F1_predicted_vs_observed.png", dpi=170)
    plt.close(fig)

    # F2 Bland-Altman
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.scatter(mean_pv, diff_pv, s=70, color="#2166ac", zorder=3)
    for m, dd, lab in zip(mean_pv, diff_pv, PLANNED_DOSES):
        ax.annotate(f"{lab:g}", (m, dd), textcoords="offset points",
                    xytext=(7, 4), fontsize=8)
    ax.axhline(ba_bias, color="#b2182b", lw=1.6, label=f"bias {ba_bias:+.2f} pp")
    ax.axhline(loa[0], color="grey", ls="--", lw=1.2,
               label=f"95% limits {loa[0]:+.1f} to {loa[1]:+.1f}")
    ax.axhline(loa[1], color="grey", ls="--", lw=1.2)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("mean of predicted and observed % viability")
    ax.set_ylabel("predicted minus observed (pp)")
    ax.set_title("Bland-Altman, sealed prediction versus observation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "F2_bland_altman.png", dpi=170)
    plt.close(fig)

    # F3 the IC50 regression and its extrapolation
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    lx = np.log10(obs_doses)
    ax.scatter(lx, obs_viab, s=70, color="#b2182b", zorder=3, label="observed doses")
    xs = np.linspace(np.log10(10), np.log10(2000), 200)
    ax.plot(xs, LAB_SLOPE * xs + LAB_INTERCEPT, color="#2166ac", lw=1.8,
            label=f"lab fit on 3 points, R² = {LAB_R2}")
    ax.axhline(50, ls=":", c="grey", lw=1.2)
    ax.axvline(np.log10(TOP_DOSE), color="black", ls="--", lw=1.2)
    ax.axvspan(np.log10(TOP_DOSE), np.log10(2000), color="red", alpha=0.07)
    ax.annotate("beyond the highest dose tested\nIC50 here is extrapolated",
                xy=(np.log10(520), 20), fontsize=8.5, color="darkred")
    ax.scatter([np.log10(LAB_IC50)], [50], marker="*", s=260, color="darkred",
               zorder=4, label=f"reported IC50 = {LAB_IC50:.0f} µg/mL")
    ax.axvline(np.log10(911.84), color="darkred", ls=":", lw=1)
    ax.set_xlabel("log₁₀ concentration (µg/mL)")
    ax.set_ylabel("% cell viability")
    ax.set_title("Where the reported IC50 comes from")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_ylim(0, 140)
    fig.tight_layout()
    fig.savefig(FIG / "F3_ic50_regression.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
