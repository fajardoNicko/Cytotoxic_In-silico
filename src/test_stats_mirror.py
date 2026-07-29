"""
test_stats_mirror.py  --  independent cross-validation of stats_mirror.py

Gate G3 has two halves:
  (a) agreement with an independent Python implementation (statsmodels / scipy)
      -- this file, automated;
  (b) agreement with jamovi to 3 dp -- manual, see docs/G3_JAMOVI_CHECK.md.

Run:  python src/test_stats_mirror.py
"""
from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from scipy import stats as sps

import stats_mirror as sm

TOL = 1e-6
failures: list[str] = []


def check(label: str, mine: float, theirs: float, tol: float = TOL) -> None:
    ok = np.isclose(mine, theirs, rtol=tol, atol=tol)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:44s} "
          f"mine={mine:>14.9f}  ref={theirs:>14.9f}  d={abs(mine-theirs):.2e}")
    if not ok:
        failures.append(label)


def make_data(seed: int, n: int, balanced: bool = True) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    spec = [("negative_control", 100.0, 5.0), ("12.5ppm", 93.0, 5.0),
            ("25ppm", 85.0, 5.5), ("50ppm", 71.0, 6.0),
            ("100ppm", 55.0, 6.5), ("200ppm", 36.0, 5.0),
            ("doxorubicin", 22.0, 4.0)]
    frames = []
    for i, (name, mu, sd) in enumerate(spec):
        k = n if balanced else n + (i % 3)
        frames.append(pd.DataFrame({"group": name,
                                    "viability_pct": rng.normal(mu, sd, k)}))
    return pd.concat(frames, ignore_index=True)


def test_anova(df):
    print("\n[1] one-way ANOVA vs scipy.stats.f_oneway + statsmodels OLS")
    import statsmodels.api as sa
    from statsmodels.formula.api import ols

    a = sm.one_way_anova(df, "group", "viability_pct")
    samples = [s["viability_pct"].to_numpy() for _, s in df.groupby("group", sort=False)]
    F_sp, p_sp = sps.f_oneway(*samples)
    check("F statistic (vs scipy)", a.F, F_sp)
    check("p value (vs scipy)", a.p, p_sp)

    tbl = sa.stats.anova_lm(ols("viability_pct ~ C(group)", data=df).fit(), typ=2)
    check("SS_treat (vs statsmodels)", a.ss_treat, tbl.loc["C(group)", "sum_sq"])
    check("SS_error (vs statsmodels)", a.ss_error, tbl.loc["Residual", "sum_sq"])
    check("df_error  (vs statsmodels)", a.df_error, tbl.loc["Residual", "df"])
    check("MSE       (vs statsmodels)", a.mse,
          tbl.loc["Residual", "sum_sq"] / tbl.loc["Residual", "df"])
    # identity that must hold exactly
    check("SS_total = SS_treat + SS_error", a.ss_total, a.ss_treat + a.ss_error)


def test_tukey(df):
    print("\n[2] Tukey HSD vs statsmodels pairwise_tukeyhsd")
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    a = sm.one_way_anova(df, "group", "viability_pct")
    ph = sm.tukey_hsd(df, "group", "viability_pct", a.mse, a.df_error)

    ref = pairwise_tukeyhsd(df["viability_pct"], df["group"], alpha=0.05)
    # NOTE: ref.summary() rounds p-adj to 4 dp for display. Compare against the
    # unrounded attributes instead, otherwise every p < 1e-4 spuriously "fails".
    pairs = [(ref.groupsunique[i], ref.groupsunique[j])
             for i, j in zip(*np.triu_indices(len(ref.groupsunique), k=1))]
    ref_map = {
        frozenset(pair): (float(md), float(p), float(lo), float(hi))
        for pair, md, p, (lo, hi) in zip(
            pairs, ref.meandiffs, ref.pvalues, ref.confint)
    }

    n_ok = 0
    for _, r in ph.table.iterrows():
        key = frozenset((r["group_1"], r["group_2"]))
        md, padj, lo, hi = ref_map[key]
        # statsmodels orients the difference as group2 - group1
        same_sign = np.isclose(abs(r["mean_diff"]), abs(md), atol=1e-8)
        p_ok = np.isclose(r["p_tukey"], padj, atol=1e-6)
        ci_ok = np.isclose(sorted([abs(r["ci_low"]), abs(r["ci_high"])]),
                           sorted([abs(lo), abs(hi)]), atol=1e-6).all()
        if same_sign and p_ok and ci_ok:
            n_ok += 1
        else:
            failures.append(f"tukey {tuple(key)}")
            print(f"  [FAIL] {tuple(key)} mine p={r['p_tukey']:.8f} ref p={padj:.8f}")
    print(f"  [{'PASS' if n_ok == len(ph.table) else 'FAIL'}] "
          f"all {len(ph.table)} pairwise comparisons match "
          f"(diff, p-adj, CI): {n_ok}/{len(ph.table)}")

    # the manuscript's closed form must equal the per-comparison margin
    if ph.hsd_critical is not None:
        margin = (ph.table["ci_high"] - ph.table["mean_diff"]).to_numpy()
        check("HSD = q*sqrt(MSE/n) equals CI margin",
              ph.hsd_critical, float(margin[0]))


def test_kruskal_dunn(df):
    print("\n[3] Kruskal-Wallis vs scipy; Dunn internal consistency")
    kw = sm.kruskal_wallis(df, "group", "viability_pct")
    samples = [s["viability_pct"].to_numpy() for _, s in df.groupby("group", sort=False)]
    H, p = sps.kruskal(*samples)
    check("H statistic", kw["H"], H)
    check("p value", kw["p"], p)

    d = sm.dunn(df, "group", "viability_pct", adjust="bonferroni")
    m = len(d.table)
    raw = d.table["p_unadj"].to_numpy()
    exp = np.minimum(raw * m, 1.0)
    check("Bonferroni adjustment", float(np.abs(d.table["p_adj"] - exp).max()), 0.0)

    dh = sm.dunn(df, "group", "viability_pct", adjust="holm")
    mono = np.all(np.diff(dh.table.sort_values("p_unadj")["p_adj"].to_numpy()) >= -1e-12)
    print(f"  [{'PASS' if mono else 'FAIL'}] Holm adjusted p-values are monotone")
    if not mono:
        failures.append("holm monotone")


def test_shapiro_levene(df):
    print("\n[4] Shapiro-Wilk / Levene vs scipy")
    n = sm.shapiro_wilk(df, "group", "viability_pct")
    first = df["group"].iloc[0]
    y = df.loc[df["group"] == first, "viability_pct"].to_numpy()
    W, p = sps.shapiro(y)
    check(f"Shapiro W ({first})", float(n.per_group.iloc[0]["W"]), W)
    check(f"Shapiro p ({first})", float(n.per_group.iloc[0]["p"]), p)

    lv = sm.levene(df, "group", "viability_pct")
    samples = [s["viability_pct"].to_numpy() for _, s in df.groupby("group", sort=False)]
    W2, p2 = sps.levene(*samples, center="median")
    check("Levene W", lv.W, W2)
    check("Levene p", lv.p, p2)


def test_unbalanced():
    print("\n[5] unbalanced design (Tukey-Kramer) vs statsmodels")
    df = make_data(7, 4, balanced=False)
    test_anova(df)
    test_tukey(df)


def test_hand_worked():
    """A tiny fixed dataset whose ANOVA can be verified by hand."""
    print("\n[6] hand-checkable fixed dataset")
    df = pd.DataFrame({
        "group": ["A"] * 3 + ["B"] * 3 + ["C"] * 3,
        "viability_pct": [100.0, 98.0, 102.0, 80.0, 82.0, 78.0, 50.0, 52.0, 48.0],
    })
    a = sm.one_way_anova(df, "group", "viability_pct")
    # T_A=300, T_B=240, T_C=150, G=690, n=9, correction=690^2/9=52900
    # SS_treat = (300^2+240^2+150^2)/3 - 52900 = (90000+57600+22500)/3 - 52900
    #          = 170100/3 - 52900 = 56700 - 52900 = 3800
    # each group has deviations -2,0,+2 -> SS_error = 3*8 = 24
    check("SS_treat (hand)", a.ss_treat, 3800.0)
    check("SS_error (hand)", a.ss_error, 24.0)
    check("grand total G (hand)", a.grand_total, 690.0)
    check("MST (hand)", a.mst, 1900.0)
    check("MSE (hand)", a.mse, 4.0)
    check("F   (hand)", a.F, 475.0)


if __name__ == "__main__":
    print("=" * 78)
    print("GATE G3 (a) -- stats_mirror.py vs independent implementations")
    print("=" * 78)

    df = make_data(20260729, 9)
    test_shapiro_levene(df)
    test_anova(df)
    test_tukey(df)
    test_kruskal_dunn(df)
    test_unbalanced()
    test_hand_worked()

    print("\n" + "=" * 78)
    if failures:
        print(f"RESULT: {len(failures)} FAILURE(S): {failures}")
        sys.exit(1)
    print("RESULT: ALL CHECKS PASSED -- Gate G3(a) satisfied.")
    print("Next: run docs/G3_JAMOVI_CHECK.md for the jamovi half.")
