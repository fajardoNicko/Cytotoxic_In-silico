"""
s07b_literature_prior.py  --  evidence-based prior for the extract IC50,
                              and the dose-range decision (Phase C5)

Why this exists instead of the bottom-up mixture model: s06c established that
the QSAR cannot supply per-compound potency for this chemical class (1 sigma
~11-fold, and the bias direction is unresolved). Summing such numbers would
manufacture false precision. A literature prior on CRUDE EXTRACT potency is
both more reliable and directly comparable to what the assay measures.

Key observation from the literature: the extraction SOLVENT dominates. Against
the same A549 line, aqueous and alkaloid extracts land near 160 ug/mL while an
ETHANOLIC extract lands near 1063 ug/mL -- a ~6.7-fold difference. This study
uses an ethanolic extract, so the solvent-matched anchor is the relevant one.

Bark-specific quantitative A549 data could not be found, which is consistent
with the manuscript's own premise that bark is underutilised. The prior is
therefore built from leaf/ethanolic analogues with deliberately widened
uncertainty, and this is recorded as a limitation.

Output: results/tables/literature_prior.csv
        results/prediction_registry/DOSE_RANGE_MEMO.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mtt_model import hill_viability, loglinear_ic50, fourpl_ic50, nci_activity_class

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
REG = ROOT / "results" / "prediction_registry"
for p in (TAB, FIG, REG):
    p.mkdir(parents=True, exist_ok=True)

DOSES_PLAN = [12.5, 25.0, 50.0, 100.0, 200.0]

# ---------------------------------------------------------------------------
# Published IC50 values, M. oleifera extracts vs A549 (and a normal-cell pair)
# ---------------------------------------------------------------------------
LITERATURE = [
    dict(extract="aqueous leaf", solvent="aqueous", cell="A549",
         ic50_ug_ml=166.7, assay="MTT",
         source="Tiloke et al. 2013, BMC Complement Altern Med 13:226",
         url="https://link.springer.com/article/10.1186/1472-6882-13-226"),
    dict(extract="alkaloid extract (MOAE)", solvent="alkaloid fraction",
         cell="A549", ic50_ug_ml=158.67, assay="MTT",
         source="Xie et al. 2021, Evid Based Complement Alternat Med 5591687",
         url="https://onlinelibrary.wiley.com/doi/10.1155/2021/5591687"),
    dict(extract="ethanolic leaf", solvent="ethanol", cell="A549",
         ic50_ug_ml=1062.87, assay="MTS",
         source="Trends in Sciences (2022), MOE vs A549",
         url="https://tis.wu.ac.th/index.php/tis/article/view/3202"),
    dict(extract="ethanolic leaf", solvent="ethanol", cell="MCF-12A (normal)",
         ic50_ug_ml=1424.04, assay="MTS",
         source="Trends in Sciences (2022), normal-cell comparator",
         url="https://tis.wu.ac.th/index.php/tis/article/view/3202"),
]


def build_priors() -> dict:
    df = pd.DataFrame(LITERATURE)
    df.to_csv(TAB / "literature_prior.csv", index=False)

    a549 = df[df["cell"] == "A549"]
    logs = np.log10(a549["ic50_ug_ml"].to_numpy())
    pooled_mu, pooled_sd = float(logs.mean()), float(logs.std(ddof=1))

    eth = a549[a549["solvent"] == "ethanol"]
    eth_mu = float(np.log10(eth["ic50_ug_ml"].to_numpy()).mean())
    # only one ethanolic datapoint, so borrow the between-study spread
    eth_sd = pooled_sd

    return {
        "table": df,
        "pooled": (pooled_mu, pooled_sd),
        "solvent_matched": (eth_mu, eth_sd),
    }


def p_exceeds(mu, sd, threshold) -> float:
    return float(sps.norm.sf(np.log10(threshold), mu, sd))


def top_dose_for_coverage(mu, sd, coverage=0.95) -> float:
    return float(10 ** (mu + sps.norm.ppf(coverage) * sd))


def main():
    pri = build_priors()
    df = pri["table"]
    p_mu, p_sd = pri["pooled"]
    e_mu, e_sd = pri["solvent_matched"]

    L: list[str] = []

    def say(s=""):
        print(s, flush=True)
        L.append(s)

    say("=" * 78)
    say("LITERATURE PRIOR -- M. oleifera extract IC50 vs A549")
    say("=" * 78)
    say(df[["extract", "solvent", "cell", "ic50_ug_ml", "assay", "source"]]
        .to_string(index=False))

    si = 1424.04 / 1062.87
    say(f"\nSelectivity Index from the ethanolic study "
        f"(MCF-12A / A549) = {si:.2f}")
    say("  -> barely selective; worth stating plainly in the discussion.")

    say("\n" + "-" * 78)
    say("PRIOR DISTRIBUTIONS (log10 IC50, ug/mL)")
    say("-" * 78)
    say(f"  pooled over all A549 extracts : GM {10 ** p_mu:8.1f}  "
        f"1sigma {10 ** (p_mu - p_sd):.0f} - {10 ** (p_mu + p_sd):.0f}")
    say(f"  solvent-matched (ethanol only): GM {10 ** e_mu:8.1f}  "
        f"1sigma {10 ** (e_mu - e_sd):.0f} - {10 ** (e_mu + e_sd):.0f}")
    say("  This study uses an ETHANOLIC extract -> the solvent-matched prior is")
    say("  the decision-relevant one; the pooled prior is the optimistic bound.")

    top = max(DOSES_PLAN)
    say("\n" + "-" * 78)
    say(f"RISK THAT THE PLANNED RANGE FAILS (top dose = {top:g} ppm)")
    say("-" * 78)
    pp = p_exceeds(p_mu, p_sd, top)
    pe = p_exceeds(e_mu, e_sd, top)
    say(f"  P(IC50 > {top:g} ppm) under pooled prior          = {pp * 100:5.1f}%")
    say(f"  P(IC50 > {top:g} ppm) under solvent-matched prior = {pe * 100:5.1f}%")
    say("")
    say("  If IC50 exceeds the top dose, viability never crosses 50%, the")
    say("  log-linear solution becomes an extrapolation, and Specific Question 3")
    say("  (IC50) cannot be answered. Monte Carlo E1 showed 0% of experiments")
    say("  bracket 50% inhibition once the true IC50 passes ~214 ug/mL.")

    say("\n" + "-" * 78)
    say("WHAT THE PLANNED EXPERIMENT WOULD ACTUALLY SEE")
    say("-" * 78)
    for label, mu in (("solvent-matched", e_mu), ("pooled", p_mu)):
        ic50 = 10 ** mu
        v = hill_viability(np.array(DOSES_PLAN), ic50, hill=1.2, v_min=8.0)
        say(f"  prior = {label:16s} (IC50 {ic50:7.1f} ug/mL)")
        say(f"    dose (ppm)      : " + "  ".join(f"{d:7.1f}" for d in DOSES_PLAN))
        say(f"    %viability      : " + "  ".join(f"{x:7.1f}" for x in v))
        say(f"    %inhibition     : " + "  ".join(f"{100 - x:7.1f}" for x in v))
        say(f"    max inhibition at {top:g} ppm = {100 - v[-1]:.1f}%")
        say("")

    say("-" * 78)
    say("RECOMMENDED DOSE SERIES")
    say("-" * 78)
    for label, mu, sd in (("solvent-matched", e_mu, e_sd), ("pooled", p_mu, p_sd)):
        need = top_dose_for_coverage(mu, sd, 0.95)
        say(f"  {label:16s}: top dose ~{need:7.0f} ppm for 95% coverage")

    # a practical 2-fold series that keeps the manuscript's anchors
    series = [50, 100, 200, 400, 800, 1600]
    say(f"\n  PRACTICAL RECOMMENDATION (2-fold series, 6 points):")
    say(f"    {', '.join(str(s) for s in series)} ppm")
    say("    Rationale: keeps 50/100/200 ppm so the original design is nested")
    say("    inside the new one, and extends far enough that the solvent-matched")
    say("    prior is bracketed. Drop 12.5 and 25 ppm -- both priors predict")
    say("    <5% inhibition there, so those wells buy no information.")
    say("")
    say("    CONSTRAINTS TO CHECK BEFORE ADOPTING:")
    say("      * extract solubility at 1600 ug/mL in the stock solvent")
    say("      * final DMSO <= 0.5% v/v in the top-dose well")
    say("      * if solubility caps below the IC50, report 'IC50 > [max soluble")
    say("        dose]' as a VALID result -- do not extrapolate the regression.")

    # ---- what the manuscript's own method would report if run as planned ----
    say("\n" + "-" * 78)
    say("IF THE PLANNED RANGE IS USED ANYWAY (solvent-matched prior)")
    say("-" * 78)
    ic50 = 10 ** e_mu
    v = hill_viability(np.array(DOSES_PLAN), ic50, hill=1.2, v_min=8.0)
    ll = loglinear_ic50(DOSES_PLAN, v)
    fp = fourpl_ic50(DOSES_PLAN, v)
    say(f"  true IC50                    = {ic50:9.1f} ug/mL")
    say(f"  log-linear estimate (paper)  = {ll.ic50:9.1f} ug/mL "
        f"({ll.ic50 / ic50:.1f}x the truth)")
    say(f"  4PL estimate                 = {fp.ic50:9.1f} ug/mL")
    say(f"  regression R^2               = {ll.r_squared:.4f}  <-- looks fine!")
    say("  The R^2 is high even though the IC50 is badly wrong. A good-looking")
    say("  regression is NOT evidence that the dose range was adequate.")
    say(f"  NCI class at the true value  : {nci_activity_class(ic50)}")

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    x = np.logspace(np.log10(1), np.log10(4000), 400)
    for label, mu, c in (("solvent-matched (ethanol)", e_mu, "#b2182b"),
                         ("pooled (all solvents)", p_mu, "#2166ac")):
        ax[0].plot(x, hill_viability(x, 10 ** mu, 1.2, 8.0), color=c, label=label)
    ax[0].axhline(50, ls="--", c="grey", lw=1)
    ax[0].axvspan(min(DOSES_PLAN), max(DOSES_PLAN), color="orange", alpha=.18,
                  label="planned range 12.5–200 ppm")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("extract concentration (µg/mL)")
    ax[0].set_ylabel("% viability")
    ax[0].set_title("Predicted dose–response vs the planned range")
    ax[0].legend(fontsize=8)

    lo, hi = 0.5, 4.0
    xs = np.linspace(lo, hi, 500)
    ax[1].plot(10 ** xs, sps.norm.pdf(xs, e_mu, e_sd), color="#b2182b",
               label="solvent-matched prior")
    ax[1].plot(10 ** xs, sps.norm.pdf(xs, p_mu, p_sd), color="#2166ac",
               label="pooled prior")
    ax[1].axvline(max(DOSES_PLAN), color="black", ls=":", lw=1.6)
    ax[1].annotate("200 ppm ceiling", xy=(210, 0.05), fontsize=8, rotation=90)
    for _, r in df[df["cell"] == "A549"].iterrows():
        ax[1].axvline(r["ic50_ug_ml"], color="grey", lw=0.8, alpha=.7)
    ax[1].set_xscale("log")
    ax[1].set_xlabel("IC$_{50}$ (µg/mL)")
    ax[1].set_ylabel("prior density (per log$_{10}$ unit)")
    ax[1].set_title("Prior on extract IC$_{50}$ (grey = published values)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "C5_dose_range_decision.png", dpi=160)
    plt.close(fig)

    memo = ["# Dose-range recommendation — ACTION REQUIRED before the MTT assay",
            "",
            "**To:** MSU-IIT partner institution / research adviser",
            "**From:** in silico arm",
            "**Status:** advisory, evidence-based, issued before any wet-lab data exists",
            "",
            "```",
            *L,
            "```",
            "",
            "## Sources",
            ""]
    for _, r in df.iterrows():
        memo.append(f"- {r['extract']} vs {r['cell']}: "
                    f"IC₅₀ = {r['ic50_ug_ml']} µg/mL ({r['assay']}) — "
                    f"{r['source']}. <{r['url']}>")
    memo += ["",
             "## Limitations of this advisory",
             "",
             "- No quantitative A549 data for *M. oleifera* **bark** could be located;",
             "  the prior is built from **leaf** extracts. Bark differs in tannin,",
             "  alkaloid and isothiocyanate content, so the true value may sit outside",
             "  the prior. This is the study's novelty and also its main uncertainty.",
             "- One of the four anchors used an **MTS** rather than MTT readout.",
             "- The single ethanolic datapoint borrows its spread from the",
             "  between-study variance of the other anchors.",
             "- This advisory does **not** depend on the QSAR model, which was found",
             "  unusable for per-compound potency (~11-fold 1σ error).",
             ]
    (REG / "DOSE_RANGE_MEMO.md").write_text("\n".join(memo), encoding="utf-8")

    say(f"\nWrote {TAB / 'literature_prior.csv'}")
    say(f"Wrote {REG / 'DOSE_RANGE_MEMO.md'}")
    say(f"Wrote {FIG / 'C5_dose_range_decision.png'}")


if __name__ == "__main__":
    main()
