"""
stats_mirror.py  --  Phase A, step A7 of IN_SILICO_PLAN.md

Reimplements, in explicit form, the exact statistical chain specified in the
research plan (pp. 12-14):

    Shapiro-Wilk  ->  one-way ANOVA   (normal)      ->  Tukey HSD
                  ->  Kruskal-Wallis  (non-normal)  ->  Dunn

The ANOVA is written out using the same sum-of-squares decomposition printed
in the paper, so every intermediate quantity (T_i, G, SS, MST, MSE, F) can be
checked by hand against the manuscript:

    SS_treat = sum(T_i^2 / n_i) - G^2 / n
    SS_error = sum(Y_ij^2)      - sum(T_i^2 / n_i)
    MST      = SS_treat / (k - 1)
    MSE      = SS_error / (n - k)
    F        = MST / MSE
    HSD      = q(alpha, k, df_error) * sqrt(MSE / n)

Gate G3: running the same CSV through this module and through jamovi must
agree to 3 decimal places.  See `verify_against_jamovi()` at the bottom for
the export helper that produces the jamovi-ready file.

No statsmodels dependency: Tukey/Games-Howell use scipy's studentized range
distribution directly, which is the same distribution jamovi (via R) uses.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from itertools import combinations
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

ALPHA = 0.05


# ---------------------------------------------------------------------------
# result containers
# ---------------------------------------------------------------------------

@dataclass
class NormalityResult:
    per_group: pd.DataFrame          # group, n, W, p, normal
    residuals_W: float               # Shapiro-Wilk on pooled residuals
    residuals_p: float
    all_groups_normal: bool
    residuals_normal: bool

    @property
    def use_parametric(self) -> bool:
        """Decision rule of the research plan: parametric iff normality holds.

        We require BOTH the per-group test and the pooled-residual test to
        pass.  The residual test is the one jamovi reports for ANOVA; the
        per-group test is what the manuscript describes.  Requiring both is
        the conservative reading.
        """
        return self.all_groups_normal and self.residuals_normal


@dataclass
class AnovaResult:
    k: int
    n_total: int
    group_totals: dict
    grand_total: float
    ss_treat: float
    ss_error: float
    ss_total: float
    df_treat: int
    df_error: int
    mst: float
    mse: float
    F: float
    p: float
    eta_sq: float
    omega_sq: float
    significant: bool


@dataclass
class LeveneResult:
    W: float
    p: float
    homogeneous: bool
    center: str = "median"


@dataclass
class PostHocResult:
    method: str
    table: pd.DataFrame
    hsd_critical: float | None = None   # only for Tukey with balanced design


@dataclass
class AnalysisReport:
    normality: NormalityResult
    levene: LeveneResult
    anova: AnovaResult | None
    kruskal: dict | None
    posthoc: PostHocResult
    descriptives: pd.DataFrame
    branch: str                       # "ANOVA+Tukey" or "Kruskal-Wallis+Dunn"


# ---------------------------------------------------------------------------
# descriptives
# ---------------------------------------------------------------------------

def descriptives(df: pd.DataFrame, group: str, value: str) -> pd.DataFrame:
    """mean +/- SD per group, as the manuscript requires for reporting."""
    g = df.groupby(group, sort=False)[value]
    out = pd.DataFrame({
        "n": g.count(),
        "mean": g.mean(),
        "sd": g.std(ddof=1),
        "sem": g.std(ddof=1) / np.sqrt(g.count()),
        "min": g.min(),
        "max": g.max(),
    })
    out["cv_pct"] = 100.0 * out["sd"] / out["mean"]
    return out.reset_index()


# ---------------------------------------------------------------------------
# assumption checks
# ---------------------------------------------------------------------------

def shapiro_wilk(df: pd.DataFrame, group: str, value: str,
                 alpha: float = ALPHA) -> NormalityResult:
    rows = []
    resid = []
    for name, sub in df.groupby(group, sort=False):
        y = sub[value].to_numpy(dtype=float)
        resid.append(y - y.mean())
        if len(y) < 3:
            rows.append({"group": name, "n": len(y), "W": np.nan,
                         "p": np.nan, "normal": True})
            continue
        W, p = stats.shapiro(y)
        rows.append({"group": name, "n": len(y), "W": W, "p": p,
                     "normal": bool(p > alpha)})
    per_group = pd.DataFrame(rows)

    pooled = np.concatenate(resid)
    rW, rp = stats.shapiro(pooled) if len(pooled) >= 3 else (np.nan, np.nan)

    return NormalityResult(
        per_group=per_group,
        residuals_W=float(rW),
        residuals_p=float(rp),
        all_groups_normal=bool(per_group["normal"].all()),
        residuals_normal=bool(rp > alpha),
    )


def levene(df: pd.DataFrame, group: str, value: str,
           alpha: float = ALPHA, center: str = "median") -> LeveneResult:
    """Homogeneity of variance.

    center='median' is the Brown-Forsythe variant, which is what jamovi
    reports (it calls car::leveneTest, default center = median).
    """
    samples = [s[value].to_numpy(dtype=float) for _, s in df.groupby(group, sort=False)]
    W, p = stats.levene(*samples, center=center)
    return LeveneResult(W=float(W), p=float(p),
                        homogeneous=bool(p > alpha), center=center)


# ---------------------------------------------------------------------------
# one-way ANOVA, written out exactly as in the manuscript
# ---------------------------------------------------------------------------

def one_way_anova(df: pd.DataFrame, group: str, value: str,
                  alpha: float = ALPHA) -> AnovaResult:
    groups = list(df[group].unique())
    k = len(groups)

    T = {}                      # T_i : total of observations in group i
    n_i = {}
    sum_sq = 0.0                # sum over all Y_ij^2
    for name in groups:
        y = df.loc[df[group] == name, value].to_numpy(dtype=float)
        T[name] = float(y.sum())
        n_i[name] = int(y.size)
        sum_sq += float((y ** 2).sum())

    G = float(sum(T.values()))              # grand total
    n = int(sum(n_i.values()))
    correction = G ** 2 / n                 # G^2 / n

    ss_treat = sum(T[g] ** 2 / n_i[g] for g in groups) - correction
    ss_total = sum_sq - correction
    ss_error = ss_total - ss_treat

    df_treat = k - 1
    df_error = n - k

    mst = ss_treat / df_treat
    mse = ss_error / df_error
    F = mst / mse
    p = float(stats.f.sf(F, df_treat, df_error))

    eta_sq = ss_treat / ss_total
    omega_sq = (ss_treat - df_treat * mse) / (ss_total + mse)

    return AnovaResult(
        k=k, n_total=n, group_totals=T, grand_total=G,
        ss_treat=ss_treat, ss_error=ss_error, ss_total=ss_total,
        df_treat=df_treat, df_error=df_error,
        mst=mst, mse=mse, F=F, p=p,
        eta_sq=eta_sq, omega_sq=omega_sq,
        significant=bool(p < alpha),
    )


# ---------------------------------------------------------------------------
# post hoc
# ---------------------------------------------------------------------------

def tukey_hsd(df: pd.DataFrame, group: str, value: str,
              mse: float, df_error: int, alpha: float = ALPHA) -> PostHocResult:
    """Tukey-Kramer HSD.

    For a balanced design this reduces to the manuscript's formula
        HSD = q(alpha, k, df_error) * sqrt(MSE / n)
    which is reported as `hsd_critical`.  For unbalanced designs the
    Tukey-Kramer standard error is used per comparison.
    """
    stats_by_group = {
        name: (sub[value].mean(), len(sub))
        for name, sub in df.groupby(group, sort=False)
    }
    names = list(stats_by_group)
    k = len(names)

    q_crit = float(stats.studentized_range.ppf(1 - alpha, k, df_error))

    ns = {n for _, n in stats_by_group.values()}
    hsd_crit = q_crit * np.sqrt(mse / ns.pop()) if len(ns) == 1 else None

    rows = []
    for a, b in combinations(names, 2):
        (m_a, n_a), (m_b, n_b) = stats_by_group[a], stats_by_group[b]
        diff = m_a - m_b
        se = np.sqrt((mse / 2.0) * (1.0 / n_a + 1.0 / n_b))
        q_obs = abs(diff) / se
        p_adj = float(stats.studentized_range.sf(q_obs, k, df_error))
        margin = q_crit * se
        rows.append({
            "group_1": a, "group_2": b,
            "mean_diff": diff, "se": se, "q": q_obs,
            "p_tukey": p_adj,
            "ci_low": diff - margin, "ci_high": diff + margin,
            "significant": bool(p_adj < alpha),
        })
    return PostHocResult("Tukey HSD", pd.DataFrame(rows), hsd_crit)


def games_howell(df: pd.DataFrame, group: str, value: str,
                 alpha: float = ALPHA) -> PostHocResult:
    """Post hoc for unequal variances (use when Levene fails)."""
    g = {name: sub[value].to_numpy(dtype=float)
         for name, sub in df.groupby(group, sort=False)}
    names = list(g)
    k = len(names)

    rows = []
    for a, b in combinations(names, 2):
        xa, xb = g[a], g[b]
        na, nb = xa.size, xb.size
        va, vb = xa.var(ddof=1), xb.var(ddof=1)
        diff = xa.mean() - xb.mean()
        se = np.sqrt(va / na + vb / nb)
        t = abs(diff) / se
        # Welch-Satterthwaite
        dfree = (va / na + vb / nb) ** 2 / (
            (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
        p_adj = float(stats.studentized_range.sf(t * np.sqrt(2), k, dfree))
        q_crit = float(stats.studentized_range.ppf(1 - alpha, k, dfree))
        margin = q_crit / np.sqrt(2) * se
        rows.append({
            "group_1": a, "group_2": b, "mean_diff": diff, "se": se,
            "t": t, "df": dfree, "p_games_howell": p_adj,
            "ci_low": diff - margin, "ci_high": diff + margin,
            "significant": bool(p_adj < alpha),
        })
    return PostHocResult("Games-Howell", pd.DataFrame(rows))


def kruskal_wallis(df: pd.DataFrame, group: str, value: str,
                   alpha: float = ALPHA) -> dict:
    samples = [s[value].to_numpy(dtype=float) for _, s in df.groupby(group, sort=False)]
    H, p = stats.kruskal(*samples)
    k = len(samples)
    n = sum(s.size for s in samples)
    return {"H": float(H), "df": k - 1, "p": float(p),
            "epsilon_sq": float((H - k + 1) / (n - k)),
            "significant": bool(p < alpha)}


def dunn(df: pd.DataFrame, group: str, value: str,
         alpha: float = ALPHA, adjust: str = "holm") -> PostHocResult:
    """Dunn's test with tie correction; Holm or Bonferroni adjustment."""
    d = df[[group, value]].dropna().copy()
    d["_rank"] = stats.rankdata(d[value].to_numpy(dtype=float))
    N = len(d)

    grp = d.groupby(group, sort=False)["_rank"]
    mean_rank = grp.mean().to_dict()
    n_i = grp.count().to_dict()
    names = list(mean_rank)
    k = len(names)

    # tie correction: sum(t^3 - t) over tied groups
    _, counts = np.unique(d[value].to_numpy(dtype=float), return_counts=True)
    ties = float(np.sum(counts ** 3 - counts))
    sigma_base = (N * (N + 1) / 12.0) - ties / (12.0 * (N - 1))

    raw = []
    for a, b in combinations(names, 2):
        diff = mean_rank[a] - mean_rank[b]
        se = np.sqrt(sigma_base * (1.0 / n_i[a] + 1.0 / n_i[b]))
        z = diff / se
        raw.append({"group_1": a, "group_2": b, "mean_rank_diff": diff,
                    "z": z, "p_unadj": float(2 * stats.norm.sf(abs(z)))})

    out = pd.DataFrame(raw)
    m = len(out)
    if adjust == "bonferroni":
        out["p_adj"] = np.minimum(out["p_unadj"] * m, 1.0)
    elif adjust == "holm":
        order = np.argsort(out["p_unadj"].to_numpy())
        adj = np.empty(m)
        running = 0.0
        for rank, idx in enumerate(order):
            val = (m - rank) * out["p_unadj"].to_numpy()[idx]
            running = max(running, min(val, 1.0))
            adj[idx] = running
        out["p_adj"] = adj
    else:
        out["p_adj"] = out["p_unadj"]
    out["significant"] = out["p_adj"] < alpha
    return PostHocResult(f"Dunn ({adjust})", out)


# ---------------------------------------------------------------------------
# the full chain
# ---------------------------------------------------------------------------

def run_analysis(df: pd.DataFrame, group: str = "group",
                 value: str = "viability_pct",
                 alpha: float = ALPHA,
                 force_branch: str | None = None) -> AnalysisReport:
    """Execute the manuscript's decision tree end to end.

    force_branch: "parametric" | "nonparametric" | None (follow the data).
    """
    desc = descriptives(df, group, value)
    norm = shapiro_wilk(df, group, value, alpha)
    lev = levene(df, group, value, alpha)

    parametric = norm.use_parametric if force_branch is None \
        else (force_branch == "parametric")

    if parametric:
        aov = one_way_anova(df, group, value, alpha)
        if lev.homogeneous:
            ph = tukey_hsd(df, group, value, aov.mse, aov.df_error, alpha)
        else:
            # Variance heterogeneity: Tukey's assumption is violated.
            ph = games_howell(df, group, value, alpha)
        return AnalysisReport(norm, lev, aov, None, ph, desc, "ANOVA+" + ph.method)

    kw = kruskal_wallis(df, group, value, alpha)
    ph = dunn(df, group, value, alpha)
    return AnalysisReport(norm, lev, None, kw, ph, desc, "Kruskal-Wallis+Dunn")


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

def format_report(rep: AnalysisReport, title: str = "") -> str:
    L = []
    w = 78
    if title:
        L += ["=" * w, title, "=" * w]

    L += ["", "DESCRIPTIVES (mean +/- SD)", "-" * w,
          rep.descriptives.to_string(index=False,
                                     float_format=lambda x: f"{x:9.3f}")]

    L += ["", "ASSUMPTION CHECKS", "-" * w, "Shapiro-Wilk, per group:",
          rep.normality.per_group.to_string(index=False,
                                            float_format=lambda x: f"{x:9.4f}"),
          f"Shapiro-Wilk on pooled residuals: W = {rep.normality.residuals_W:.4f}, "
          f"p = {rep.normality.residuals_p:.4f}",
          f"Levene (center={rep.levene.center}): W = {rep.levene.W:.4f}, "
          f"p = {rep.levene.p:.4f}  -> variances "
          f"{'homogeneous' if rep.levene.homogeneous else 'NOT homogeneous'}",
          f"Branch taken: {rep.branch}"]

    if rep.anova is not None:
        a = rep.anova
        L += ["", "ONE-WAY ANOVA", "-" * w,
              f"  Grand total G          = {a.grand_total:.4f}",
              f"  N                      = {a.n_total}   k = {a.k}",
              f"  SS_treatments          = {a.ss_treat:.4f}   df = {a.df_treat}",
              f"  SS_error               = {a.ss_error:.4f}   df = {a.df_error}",
              f"  SS_total               = {a.ss_total:.4f}",
              f"  MST                    = {a.mst:.4f}",
              f"  MSE                    = {a.mse:.4f}",
              f"  F = MST/MSE            = {a.F:.4f}",
              f"  p                      = {a.p:.6g}",
              f"  eta^2 = {a.eta_sq:.4f}    omega^2 = {a.omega_sq:.4f}",
              f"  -> {'SIGNIFICANT' if a.significant else 'not significant'} at alpha=0.05"]

    if rep.kruskal is not None:
        k = rep.kruskal
        L += ["", "KRUSKAL-WALLIS", "-" * w,
              f"  H = {k['H']:.4f}   df = {k['df']}   p = {k['p']:.6g}",
              f"  epsilon^2 = {k['epsilon_sq']:.4f}",
              f"  -> {'SIGNIFICANT' if k['significant'] else 'not significant'}"]

    L += ["", f"POST HOC: {rep.posthoc.method}", "-" * w]
    if rep.posthoc.hsd_critical is not None:
        L.append(f"  Critical HSD = q * sqrt(MSE/n) = {rep.posthoc.hsd_critical:.4f}")
        L.append("  (any |mean difference| exceeding this is significant)")
    L.append(rep.posthoc.table.to_string(index=False,
                                         float_format=lambda x: f"{x:9.4f}"))
    return "\n".join(L)


def verify_against_jamovi(df: pd.DataFrame, path: str,
                          group: str = "group", value: str = "viability_pct") -> str:
    """Write the long-format CSV to open in jamovi for the Gate G3 check.

    In jamovi:
      ANOVA > One-Way ANOVA (Fisher's, assume equal variances) for F and p
      ANOVA > One-Way ANOVA > Normality test (Shapiro-Wilk), Homogeneity test
      ANOVA > ANOVA > Post Hoc Tests > Tukey
      ANOVA > One-Way ANOVA (Non-parametric) > Kruskal-Wallis + DSCF/Dunn
    Every reported value must match this module to 3 decimal places.
    """
    df[[group, value]].to_csv(path, index=False)
    return path


if __name__ == "__main__":
    # self-test on a fixed synthetic dataset with a known answer
    rng = np.random.default_rng(20260729)
    frames = []
    for name, mu in [("negative_control", 100.0), ("12.5ppm", 92.0),
                     ("25ppm", 84.0), ("50ppm", 70.0),
                     ("100ppm", 52.0), ("200ppm", 33.0),
                     ("doxorubicin", 21.0)]:
        frames.append(pd.DataFrame({
            "group": name,
            "viability_pct": rng.normal(mu, 5.0, 9),
        }))
    demo = pd.concat(frames, ignore_index=True)
    rep = run_analysis(demo)
    print(format_report(rep, "SELF-TEST  --  synthetic dose-response"))
