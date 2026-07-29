"""
s08_monte_carlo.py  --  Phase D-4/D-5 of IN_SILICO_PLAN.md

Runs many thousands of virtual repetitions of the MTT experiment to answer
questions that cannot be answered after the fact:

  E1  Dose-range adequacy (defect D2)
      Over a sweep of plausible true IC50 values, how often does the
      12.5-200 ppm series actually bracket 50% inhibition, and how badly does
      the manuscript's log-linear method mis-estimate IC50 when it does not?

  E2  Statistical power (defect D9)
      With n technical replicates x m biological runs, what is the probability
      that ANOVA is significant and that Tukey flags each individual dose
      against the negative control?

  E3  Vehicle-control bias (defect D1)
      How much does normalising to a medium-only control instead of a DMSO
      vehicle control distort the reported IC50?

  E4  MTT interference bias (defect D6)
      How much does direct tetrazolium reduction by polyphenols inflate
      apparent viability, and does the cell-free control recover the truth?

  E5  Edge effect
      Cost of using the perimeter wells of the plate.

Outputs: results/tables/*.csv and results/figures/*.png
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mtt_model import loglinear_ic50, fourpl_ic50, hill_viability
from virtual_plate import (AssayConfig, simulate_experiment, compute_viability,
                           dose_response_summary, true_absolute_ic50,
                           build_layout, DOSES_PLAN)

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05


# ---------------------------------------------------------------------------
# fast analysis path (avoids the full report machinery inside the MC loop)
# ---------------------------------------------------------------------------

def anova_and_tukey_vs_control(viab: pd.DataFrame, doses, control="negative_control"):
    """One-way ANOVA p, plus Tukey p for each dose vs the control group."""
    groups = [control] + [f"{d:g}ppm" for d in doses]
    data = {g: viab.loc[viab["group"] == g, "viability_pct"].to_numpy() for g in groups}

    k = len(groups)
    n_tot = sum(v.size for v in data.values())
    grand = np.concatenate(list(data.values()))
    G = grand.sum()
    corr = G ** 2 / n_tot
    ss_treat = sum(v.sum() ** 2 / v.size for v in data.values()) - corr
    ss_total = (grand ** 2).sum() - corr
    ss_error = ss_total - ss_treat
    df_t, df_e = k - 1, n_tot - k
    mse = ss_error / df_e
    F = (ss_treat / df_t) / mse
    p_anova = float(sps.f.sf(F, df_t, df_e))

    ctrl = data[control]
    qs, ns = [], []
    for d in doses:
        v = data[f"{d:g}ppm"]
        se = np.sqrt((mse / 2.0) * (1.0 / ctrl.size + 1.0 / v.size))
        qs.append(abs(ctrl.mean() - v.mean()) / se)
    p_tukey = sps.studentized_range.sf(np.array(qs), k, df_e)
    return p_anova, np.asarray(p_tukey, float), mse


def one_experiment(cfg: AssayConfig, normalize_to="negative_control",
                   subtract_interference=False, do_4pl=True):
    raw = simulate_experiment(cfg)
    viab = compute_viability(raw, normalize_to=normalize_to,
                             subtract_interference=subtract_interference)
    s = dose_response_summary(viab)
    ll = loglinear_ic50(s["conc_ppm"], s["mean_viability"])
    fp = fourpl_ic50(s["conc_ppm"], s["mean_viability"]) if do_4pl else None

    top_v = float(s["mean_viability"].iloc[-1])
    bot_v = float(s["mean_viability"].iloc[0])
    bracketed = bool(top_v < 50.0 <= bot_v)

    return {"viab": viab, "summary": s, "ll": ll, "fp": fp,
            "bracketed": bracketed, "top_viability": top_v}


# ---------------------------------------------------------------------------
# E1 -- dose-range adequacy
# ---------------------------------------------------------------------------

def experiment_1(n_iter=1500, seed0=100000):
    print("\n[E1] Dose-range adequacy across plausible true IC50 values")
    midpoints = [25, 50, 75, 100, 150, 200, 300, 450, 700, 1000]
    rows = []
    for t in midpoints:
        # The estimand is the ABSOLUTE IC50 (viability = 50%), which is what
        # the manuscript defines -- not the 4PL midpoint parameter t.
        truth = true_absolute_ic50(AssayConfig(true_ic50=t))
        ll_est, fp_est, brack, extrap = [], [], [], []
        for i in range(n_iter):
            cfg = AssayConfig(true_ic50=t, seed=seed0 + i)
            r = one_experiment(cfg)
            ll_est.append(r["ll"].ic50)
            fp_est.append(r["fp"].absolute_ic50 if r["fp"] else np.nan)
            brack.append(r["bracketed"])
            extrap.append(r["ll"].extrapolated)
        ll_est = np.array(ll_est, float)
        fp_est = np.array(fp_est, float)
        rows.append({
            "midpoint_param": t,
            "true_absolute_ic50": truth,
            "p_bracketed": float(np.mean(brack)),
            "p_extrapolated": float(np.mean(extrap)),
            "loglin_median": float(np.nanmedian(ll_est)),
            "loglin_bias_pct": float(100 * (np.nanmedian(ll_est) - truth) / truth),
            "loglin_cv_pct": float(100 * np.nanstd(ll_est) / np.nanmean(ll_est)),
            "fourpl_median": float(np.nanmedian(fp_est)),
            "fourpl_bias_pct": float(100 * (np.nanmedian(fp_est) - truth) / truth),
            "fourpl_cv_pct": float(100 * np.nanstd(fp_est) / np.nanmean(fp_est)),
        })
        print(f"   true abs IC50 {truth:7.1f}: bracketed {rows[-1]['p_bracketed']*100:5.1f}%  "
              f"log-lin bias {rows[-1]['loglin_bias_pct']:+9.1f}%  "
              f"4PL bias {rows[-1]['fourpl_bias_pct']:+7.1f}%")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E1_dose_range_adequacy.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].plot(df["true_absolute_ic50"], 100 * df["p_bracketed"], "o-", color="#2166ac")
    ax[0].axvline(200, ls=":", c="red", lw=1.4)
    ax[0].axhline(80, ls="--", c="grey", lw=1)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("true absolute IC$_{50}$ (µg/mL)")
    ax[0].set_ylabel("% of experiments that bracket 50% inhibition")
    ax[0].set_title("Can the 12.5–200 ppm series determine IC$_{50}$?")
    ax[0].annotate("top dose = 200 ppm", xy=(205, 50), color="darkred", fontsize=9)

    ax[1].axhline(0, c="grey", lw=1)
    ax[1].plot(df["true_absolute_ic50"], df["loglin_bias_pct"], "o-",
               label="log-linear (manuscript)", color="#b2182b")
    ax[1].plot(df["true_absolute_ic50"], df["fourpl_bias_pct"], "s-",
               label="4PL (recommended)", color="#2166ac")
    ax[1].axvline(200, ls=":", c="red", lw=1.4)
    ax[1].set_xscale("log")
    ax[1].set_yscale("symlog", linthresh=10)
    ax[1].set_xlabel("true absolute IC$_{50}$ (µg/mL)")
    ax[1].set_ylabel("median bias in estimated IC$_{50}$ (%)")
    ax[1].set_title("IC$_{50}$ estimator bias by method")
    ax[1].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "E1_dose_range_adequacy.png", dpi=160)
    plt.close(fig)
    return df


# ---------------------------------------------------------------------------
# E2 -- statistical power
# ---------------------------------------------------------------------------

def experiment_2(n_iter=1500, seed0=200000):
    print("\n[E2] Statistical power vs replicate structure")
    rows = []
    for true_ic50 in (150.0, 300.0):
        for n_tech in (3, 4, 6, 8):
            for n_bio in (1, 3):
                sig_anova = 0
                sig_dose = np.zeros(len(DOSES_PLAN))
                nonnormal = 0
                for i in range(n_iter):
                    cfg = AssayConfig(true_ic50=true_ic50, n_tech=n_tech,
                                      n_bio=n_bio, seed=seed0 + i)
                    raw = simulate_experiment(cfg)
                    viab = compute_viability(raw)
                    p_a, p_t, _ = anova_and_tukey_vs_control(viab, DOSES_PLAN)
                    sig_anova += p_a < ALPHA
                    sig_dose += (p_t < ALPHA)
                    y = viab.loc[viab["group"] == "negative_control",
                                 "viability_pct"].to_numpy()
                    if y.size >= 3 and sps.shapiro(y)[1] <= ALPHA:
                        nonnormal += 1
                needs_edge = bool(build_layout(
                    AssayConfig(n_tech=n_tech))["perimeter_required"].iloc[0])
                r = {"true_ic50": true_ic50, "n_tech": n_tech, "n_bio": n_bio,
                     "wells_per_group": n_tech * n_bio,
                     "perimeter_required": needs_edge,
                     "power_anova": sig_anova / n_iter,
                     "p_shapiro_reject": nonnormal / n_iter}
                for d, s in zip(DOSES_PLAN, sig_dose / n_iter):
                    r[f"power_{d:g}ppm"] = float(s)
                rows.append(r)
                print(f"   IC50={true_ic50:5.0f} n_tech={n_tech} n_bio={n_bio} "
                      f"(N={r['wells_per_group']:2d}/group"
                      f"{', EDGE WELLS NEEDED' if needs_edge else ''}): "
                      f"ANOVA {r['power_anova']*100:5.1f}%  "
                      f"12.5ppm {r['power_12.5ppm']*100:5.1f}%  "
                      f"25ppm {r['power_25ppm']*100:5.1f}%  "
                      f"50ppm {r['power_50ppm']*100:5.1f}%")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E2_power.csv", index=False)

    sub = df[(df["true_ic50"] == 150.0) & (df["n_bio"] == 3)]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    x = np.arange(len(DOSES_PLAN))
    width = 0.26
    for j, (_, r) in enumerate(sub.iterrows()):
        vals = [r[f"power_{d:g}ppm"] * 100 for d in DOSES_PLAN]
        ax.bar(x + (j - 1) * width, vals, width,
               label=f"n={int(r['n_tech'])} tech × 3 runs (N={int(r['wells_per_group'])})")
    ax.axhline(80, ls="--", c="red", lw=1.2)
    ax.text(len(x) - 0.5, 82, "80% power", color="red", fontsize=9, ha="right")
    ax.set_xticks(x, [f"{d:g}" for d in DOSES_PLAN])
    ax.set_xlabel("extract concentration (ppm)")
    ax.set_ylabel("power: Tukey vs negative control (%)")
    ax.set_title("Power to detect each dose (true IC$_{50}$ = 150 µg/mL)")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "E2_power.png", dpi=160)
    plt.close(fig)
    return df


# ---------------------------------------------------------------------------
# E3 -- vehicle-control bias
# ---------------------------------------------------------------------------

def experiment_3(n_iter=1000, seed0=300000):
    print("\n[E3] Vehicle-control bias (defect D1)")
    rows = []
    # The estimand is the extract's own absolute IC50, free of solvent effect;
    # that is the vehicle-normalised truth, and it does not depend on dmso_v.
    truth = true_absolute_ic50(AssayConfig(true_ic50=150.0),
                               normalize_to="vehicle_control")
    print(f"   estimand (extract-only absolute IC50) = {truth:.1f} µg/mL")
    for dmso_v in (100.0, 96.0, 90.0, 80.0):
        for norm in ("negative_control", "vehicle_control"):
            est = []
            for i in range(n_iter):
                cfg = AssayConfig(true_ic50=150.0, dmso_viability=dmso_v,
                                  seed=seed0 + i)
                r = one_experiment(cfg, normalize_to=norm)
                est.append(r["fp"].absolute_ic50)
            est = np.array(est, float)
            rows.append({"dmso_viability_pct": dmso_v, "normalized_to": norm,
                         "true_ic50": truth,
                         "median_ic50": float(np.nanmedian(est)),
                         "bias_pct": float(100 * (np.nanmedian(est) - truth) / truth)})
            print(f"   DMSO leaves {dmso_v:5.1f}% viability, normalise to "
                  f"{norm:17s}: IC50 {rows[-1]['median_ic50']:7.1f} "
                  f"({rows[-1]['bias_pct']:+6.1f}%)")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E3_vehicle_bias.csv", index=False)
    return df


# ---------------------------------------------------------------------------
# E4 -- MTT interference
# ---------------------------------------------------------------------------

def experiment_4(n_iter=1000, seed0=400000):
    print("\n[E4] MTT interference by polyphenols (defect D6)")
    rows = []
    truth = true_absolute_ic50(AssayConfig(true_ic50=150.0))
    print(f"   estimand (absolute IC50) = {truth:.1f} µg/mL")
    for k in (0.0, 0.0005, 0.0010, 0.0020):
        for sub in (False, True):
            est, top = [], []
            for i in range(n_iter):
                cfg = AssayConfig(true_ic50=150.0, interference_k=k, seed=seed0 + i)
                r = one_experiment(cfg, subtract_interference=sub)
                est.append(r["fp"].absolute_ic50)
                top.append(r["top_viability"])
            est = np.array(est, float)
            rows.append({"interference_k": k, "cell_free_control_used": sub,
                         "true_ic50": truth,
                         "median_ic50": float(np.nanmedian(est)),
                         "bias_pct": float(100 * (np.nanmedian(est) - truth) / truth),
                         "apparent_viability_200ppm": float(np.nanmedian(top))})
            print(f"   k={k:.4f} OD/(µg/mL), cell-free control "
                  f"{'ON ' if sub else 'OFF'}: IC50 {rows[-1]['median_ic50']:7.1f} "
                  f"({rows[-1]['bias_pct']:+7.1f}%)  "
                  f"apparent viability @200ppm = {rows[-1]['apparent_viability_200ppm']:5.1f}%")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E4_mtt_interference.csv", index=False)
    return df


# ---------------------------------------------------------------------------
# E5 -- edge effect
# ---------------------------------------------------------------------------

def experiment_5(n_iter=1000, seed0=500000):
    print("\n[E5] Plate edge effect")
    rows = []
    truth = true_absolute_ic50(AssayConfig(true_ic50=150.0))
    for perim in (False, True):
        est = []
        for i in range(n_iter):
            cfg = AssayConfig(true_ic50=150.0, use_perimeter=perim, seed=seed0 + i)
            r = one_experiment(cfg)
            est.append(r["fp"].absolute_ic50)
        est = np.array(est, float)
        rows.append({"perimeter_wells_used": perim, "true_ic50": truth,
                     "median_ic50": float(np.nanmedian(est)),
                     "bias_pct": float(100 * (np.nanmedian(est) - truth) / truth),
                     "cv_pct": float(100 * np.nanstd(est) / np.nanmean(est))})
        print(f"   perimeter {'USED   ' if perim else 'EMPTY  '}: "
              f"IC50 {rows[-1]['median_ic50']:7.1f} ({rows[-1]['bias_pct']:+6.1f}%)  "
              f"CV {rows[-1]['cv_pct']:5.1f}%")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E5_edge_effect.csv", index=False)
    return df


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    t0 = time.time()
    print("=" * 78)
    print(f"MONTE CARLO -- {n} iterations per condition")
    print("=" * 78)
    experiment_1(n)
    experiment_2(n)
    experiment_3(max(n // 2, 200))
    experiment_4(max(n // 2, 200))
    experiment_5(max(n // 2, 200))
    print(f"\nDone in {time.time()-t0:.1f} s. Tables -> {TAB}, figures -> {FIG}")
