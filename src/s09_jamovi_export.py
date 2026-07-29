"""
s09_jamovi_export.py  --  Gate G3(b) of IN_SILICO_PLAN.md

Generates a FIXED reference dataset (seeded, reproducible), runs it through
stats_mirror.py, and writes:

  results/tables/G3_jamovi_input.csv       <- open this in jamovi
  docs/G3_JAMOVI_CHECK.md                  <- expected values + click path

The manuscript's analysis will be done in jamovi. This proves that the
simulation and the real analysis compute the same numbers, so the Monte Carlo
power estimates transfer to the actual experiment.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import stats_mirror as sm
from virtual_plate import AssayConfig, simulate_experiment, compute_viability, DOSES_PLAN

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
DOCS = ROOT / "docs"
TAB.mkdir(parents=True, exist_ok=True)
DOCS.mkdir(parents=True, exist_ok=True)

SEED = 424242


def main():
    cfg = AssayConfig(true_ic50=150.0, n_tech=3, n_bio=3, seed=SEED)
    raw = simulate_experiment(cfg)
    viab = compute_viability(raw, normalize_to="negative_control")

    keep = ["negative_control", "vehicle_control", "doxorubicin"] + \
           [f"{d:g}ppm" for d in DOSES_PLAN]
    df = viab[viab["group"].isin(keep)][["group", "viability_pct"]].copy()
    df["group"] = pd.Categorical(df["group"], categories=keep, ordered=True)
    df = df.sort_values("group").reset_index(drop=True)

    csv = TAB / "G3_jamovi_input.csv"
    df.to_csv(csv, index=False)

    # Which branch the data-driven rule takes on THIS dataset. It is not always
    # the parametric one: with 9 wells per group, Shapiro-Wilk rejects often
    # enough that the manuscript's fallback path gets exercised for real.
    auto = sm.run_analysis(df, "group", "viability_pct")

    # Both branches are verified, because the manuscript may end up on either.
    rep = sm.run_analysis(df, "group", "viability_pct", force_branch="parametric")
    npar = sm.run_analysis(df, "group", "viability_pct", force_branch="nonparametric")

    txt = (sm.format_report(rep, f"G3 REFERENCE -- PARAMETRIC (seed={SEED})")
           + "\n\n" + sm.format_report(npar, "G3 REFERENCE -- NON-PARAMETRIC"))
    (TAB / "G3_python_reference.txt").write_text(txt, encoding="utf-8")
    print(txt)

    a = rep.anova
    ph = rep.posthoc.table

    lines = [
        "# Gate G3(b) — jamovi cross-check",
        "",
        "Gate G3 requires the simulation's statistics and the manuscript's",
        "actual analysis software to agree **to 3 decimal places**. Part (a)",
        "(agreement with statsmodels/scipy) is automated in",
        "`src/test_stats_mirror.py` and passes. Part (b) is this manual check.",
        "",
        "## 1. Open the data",
        "",
        f"Open `results/tables/G3_jamovi_input.csv` in jamovi "
        f"({len(df)} rows, 2 columns).",
        "Set `group` to **Nominal** and `viability_pct` to **Continuous**.",
        "",
        "## 2. Run these analyses",
        "",
        "**ANOVA → One-Way ANOVA**",
        "- Dependent: `viability_pct`, Grouping: `group`",
        "- Variances: tick **Assume equal (Fisher's)**",
        "- Additional Statistics: tick **Descriptives table**",
        "- Assumption Checks: tick **Normality test** and **Homogeneity test**",
        "",
        "**ANOVA → ANOVA** (for the post hoc)",
        "- Dependent: `viability_pct`, Fixed Factors: `group`",
        "- Post Hoc Tests: move `group` across, correction **Tukey**",
        "",
        "## 3. Expected values",
        "",
        "### One-way ANOVA",
        "",
        "| Quantity | Expected |",
        "|---|---|",
        f"| F | **{a.F:.3f}** |",
        f"| df (between, within) | {a.df_treat}, {a.df_error} |",
        f"| p | {a.p:.3g} |",
        f"| Sum of Squares (group) | {a.ss_treat:.3f} |",
        f"| Sum of Squares (residual) | {a.ss_error:.3f} |",
        f"| Mean Square (group) = MST | {a.mst:.3f} |",
        f"| Mean Square (residual) = MSE | {a.mse:.3f} |",
        f"| η² | {a.eta_sq:.4f} |",
        "",
        "> jamovi's *One-Way ANOVA* panel defaults to **Welch's**. Switch to",
        "> **Fisher's (assume equal variances)** or the F will not match — Welch",
        "> is a different test, not a rounding difference.",
        "",
        "### Assumption checks",
        "",
        "| Test | Expected |",
        "|---|---|",
        f"| Shapiro–Wilk on residuals, W | {rep.normality.residuals_W:.3f} |",
        f"| Shapiro–Wilk on residuals, p | {rep.normality.residuals_p:.3f} |",
        f"| Levene's (center = median), F | {rep.levene.W:.3f} |",
        f"| Levene's, p | {rep.levene.p:.3f} |",
        "",
        "> jamovi's homogeneity test uses `car::leveneTest`, which centres on the",
        "> **median** (Brown–Forsythe). `stats_mirror.levene()` matches that by",
        "> default. Centring on the mean gives a different number.",
        "",
        "### Descriptives (mean ± SD)",
        "",
        "| Group | n | Mean | SD |",
        "|---|---|---|---|",
    ]
    for _, r in rep.descriptives.iterrows():
        lines.append(f"| {r['group']} | {int(r['n'])} | {r['mean']:.3f} | {r['sd']:.3f} |")

    lines += [
        "",
        f"### Tukey HSD — critical value",
        "",
        f"Balanced design, so the manuscript's closed form applies:",
        "",
        f"    HSD = q(0.05, k={a.k}, df={a.df_error}) × √(MSE/n) = "
        f"**{rep.posthoc.hsd_critical:.4f}**",
        "",
        "Any absolute mean difference above this is significant.",
        "",
        "### Tukey HSD — comparisons vs. negative control",
        "",
        "| Comparison | Mean diff | p (Tukey) |",
        "|---|---|---|",
    ]
    for _, r in ph[(ph["group_1"] == "negative_control") |
                   (ph["group_2"] == "negative_control")].iterrows():
        other = r["group_2"] if r["group_1"] == "negative_control" else r["group_1"]
        lines.append(f"| negative_control vs {other} | {r['mean_diff']:.3f} | "
                     f"{r['p_tukey']:.4f} |")

    kw = npar.kruskal
    lines += [
        "",
        "## 5. Non-parametric branch — also verify",
        "",
        f"On this dataset the manuscript's own decision rule actually selects "
        f"**{auto.branch}** "
        f"(Shapiro–Wilk on residuals p = {auto.normality.residuals_p:.4f}).",
        "This is not a defect — with 9 wells per group the normality test",
        "rejects a noticeable fraction of the time, so the fallback path in the",
        "research plan is a live possibility and must be verified too.",
        "",
        "**ANOVA → One-Way ANOVA (Non-parametric) / Kruskal-Wallis**",
        "",
        "| Quantity | Expected |",
        "|---|---|",
        f"| χ² (H) | **{kw['H']:.3f}** |",
        f"| df | {kw['df']} |",
        f"| p | {kw['p']:.3g} |",
        f"| ε² | {kw['epsilon_sq']:.4f} |",
        "",
        "> jamovi's non-parametric post hoc is **DSCF** (Dwass-Steel-Critchlow-",
        "> Fligner), whereas `stats_mirror` implements **Dunn**. These are",
        "> different procedures and will NOT match — that is expected, not a",
        "> failure. Compare the Kruskal-Wallis H and p only, and state in the",
        "> manuscript which post hoc was used.",
        "",
        "## 4. Pass criterion",
        "",
        "Every value above must match jamovi to **3 decimal places**.",
        "Sign conventions on mean differences may be reversed (jamovi reports",
        "group2 − group1); magnitudes and p-values must agree.",
        "",
        "If anything disagrees, do **not** proceed to Phase F — the Monte Carlo",
        "power estimates would not describe the analysis actually being run.",
        "",
        "Record the outcome, the jamovi version, and the date here:",
        "",
        "| Date | jamovi version | Result | Checked by |",
        "|---|---|---|---|",
        "|  |  |  |  |",
    ]

    (DOCS / "G3_JAMOVI_CHECK.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {csv}")
    print(f"Wrote {DOCS / 'G3_JAMOVI_CHECK.md'}")


if __name__ == "__main__":
    main()
