"""
s07_mixture_model.py  --  Phase C of IN_SILICO_PLAN.md

Turns per-compound potency into an extract-level IC50 in ug/mL (= ppm), which
is the quantity the wet lab will actually measure.

Two mixture models, deliberately bracketing the truth:

  Concentration addition (Loewe).  Components act on the same target/pathway.
  Component i is present at p_i * C when the extract is dosed at C, so the
  mixture reaches 50% effect when sum(p_i * C / IC50_i) = 1:

        IC50_mix = 1 / sum( p_i / IC50_i )

  Bliss independence.  Components act independently; survival fractions
  multiply:

        V_mix(C) = prod_i V_i(p_i * C)

  solved numerically for V_mix = 50%.

Two normalisations of p_i, because MS never identifies 100% of the mass:

  CONSERVATIVE  p_i over TOTAL extract mass; unidentified mass assumed inert
                -> weaker (higher) IC50.  This is the defensible headline.
  OPTIMISTIC    p_i renormalised over the IDENTIFIED fraction only; assumes
                the unknown mass is as potent as the known mass
                -> stronger (lower) IC50.

Report the interval, not a point estimate.  If MS identifies only half the
mass the interval will be wide -- that is honest, and hiding it would not be.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from mtt_model import hill_viability, nci_activity_class

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
PROC = ROOT / "data" / "processed"
TAB.mkdir(parents=True, exist_ok=True)


@dataclass
class MixtureResult:
    ic50_ca_conservative: float
    ic50_ca_optimistic: float
    ic50_bliss_conservative: float
    ic50_bliss_optimistic: float
    identified_mass_fraction: float
    n_components: int
    n_in_domain: int
    driver_table: pd.DataFrame

    @property
    def interval(self) -> tuple[float, float]:
        vals = [self.ic50_ca_conservative, self.ic50_ca_optimistic,
                self.ic50_bliss_conservative, self.ic50_bliss_optimistic]
        vals = [v for v in vals if np.isfinite(v)]
        return (min(vals), max(vals)) if vals else (np.nan, np.nan)

    def summary(self) -> str:
        lo, hi = self.interval
        L = [
            f"components used           : {self.n_components} "
            f"({self.n_in_domain} inside QSAR applicability domain)",
            f"identified mass fraction  : {self.identified_mass_fraction * 100:.1f}%",
            "",
            f"  concentration addition, conservative : "
            f"{self.ic50_ca_conservative:10.1f} ug/mL",
            f"  concentration addition, optimistic   : "
            f"{self.ic50_ca_optimistic:10.1f} ug/mL",
            f"  Bliss independence,     conservative : "
            f"{self.ic50_bliss_conservative:10.1f} ug/mL",
            f"  Bliss independence,     optimistic   : "
            f"{self.ic50_bliss_optimistic:10.1f} ug/mL",
            "",
            f"PREDICTED EXTRACT IC50 INTERVAL : {lo:.1f} - {hi:.1f} ug/mL (ppm)",
            f"NCI class (conservative CA)     : "
            f"{nci_activity_class(self.ic50_ca_conservative)}",
        ]
        return "\n".join(L)


def concentration_addition(mass_fractions, ic50s) -> float:
    p = np.asarray(mass_fractions, float)
    ic = np.asarray(ic50s, float)
    ok = np.isfinite(p) & np.isfinite(ic) & (ic > 0) & (p > 0)
    if not ok.any():
        return np.nan
    denom = float((p[ok] / ic[ok]).sum())
    return 1.0 / denom if denom > 0 else np.inf


def bliss_independence(mass_fractions, ic50s, hills=None, v_min=0.0,
                       hi_bound=1e7) -> float:
    p = np.asarray(mass_fractions, float)
    ic = np.asarray(ic50s, float)
    ok = np.isfinite(p) & np.isfinite(ic) & (ic > 0) & (p > 0)
    p, ic = p[ok], ic[ok]
    if p.size == 0:
        return np.nan
    h = np.ones_like(p) if hills is None else np.asarray(hills, float)[ok]

    def v_mix(C):
        v = hill_viability(p * C, ic, h, v_min, 100.0) / 100.0
        return float(np.prod(v)) * 100.0

    if v_mix(hi_bound) > 50.0:
        return np.inf
    if v_mix(1e-9) < 50.0:
        return np.nan
    return float(brentq(lambda C: v_mix(C) - 50.0, 1e-9, hi_bound, xtol=1e-6))


def predict_extract_ic50(composition: pd.DataFrame,
                         potency: pd.DataFrame,
                         key: str = "pubchem_cid",
                         require_in_domain: bool = False) -> MixtureResult:
    """
    composition: needs `key` and `pct_area` (relative peak area, % of total)
    potency:     needs `key`, `pred_IC50_ug_mL`, `in_applicability_domain`
    """
    m = composition.merge(potency, on=key, how="left", suffixes=("", "_p"))

    total_area = float(composition["pct_area"].sum())
    m["mass_fraction_total"] = m["pct_area"] / 100.0     # % of whole extract

    usable = m["pred_IC50_ug_mL"].notna() & (m["pred_IC50_ug_mL"] > 0)
    if require_in_domain:
        usable &= m["in_applicability_domain"].fillna(False)
    u = m[usable].copy()

    identified = float(u["mass_fraction_total"].sum())
    u["mass_fraction_identified"] = (u["mass_fraction_total"] / identified
                                     if identified > 0 else np.nan)

    ic = u["pred_IC50_ug_mL"].to_numpy()
    ca_cons = concentration_addition(u["mass_fraction_total"], ic)
    ca_opt = concentration_addition(u["mass_fraction_identified"], ic)
    bl_cons = bliss_independence(u["mass_fraction_total"], ic)
    bl_opt = bliss_independence(u["mass_fraction_identified"], ic)

    # which components actually drive the result under CA?
    u["ca_contribution"] = u["mass_fraction_total"] / ic
    tot = u["ca_contribution"].sum()
    u["pct_of_effect"] = 100.0 * u["ca_contribution"] / tot if tot > 0 else np.nan
    drivers = u.sort_values("pct_of_effect", ascending=False)[
        [c for c in ["name", "pct_area", "pred_IC50_ug_mL",
                     "in_applicability_domain", "pct_of_effect"] if c in u.columns]]

    return MixtureResult(
        ic50_ca_conservative=ca_cons, ic50_ca_optimistic=ca_opt,
        ic50_bliss_conservative=bl_cons, ic50_bliss_optimistic=bl_opt,
        identified_mass_fraction=identified,
        n_components=int(len(u)),
        n_in_domain=int(u["in_applicability_domain"].fillna(False).sum())
        if "in_applicability_domain" in u else 0,
        driver_table=drivers,
    )


def dose_range_recommendation(ic50_lo: float, ic50_hi: float,
                              doses=(12.5, 25, 50, 100, 200)) -> str:
    """Phase C5: is the manuscript's 12.5-200 ppm series adequate?"""
    lo, hi = min(doses), max(doses)
    L = [f"Planned series      : {', '.join(f'{d:g}' for d in doses)} ppm",
         f"Predicted IC50      : {ic50_lo:.1f} - {ic50_hi:.1f} ppm", ""]
    if ic50_hi > hi:
        top = 2 ** np.ceil(np.log2(ic50_hi * 2 / hi)) * hi
        L += ["VERDICT: RANGE TOO LOW -- ACTION REQUIRED.",
              f"  The predicted IC50 may exceed the top dose ({hi:g} ppm).",
              "  If it does, viability never crosses 50%, the log-linear",
              "  solution becomes an EXTRAPOLATION, and the IC50 objective",
              "  (Specific Question 3) cannot be answered.",
              f"  RECOMMEND extending the series to ~{top:g} ppm, e.g. "
              f"{', '.join(f'{hi*2**i:g}' for i in range(1, 4))} ppm.",
              "  Check extract solubility and keep DMSO <= 0.5% v/v at the top dose."]
    elif ic50_lo < lo:
        L += ["VERDICT: RANGE TOO HIGH.",
              f"  Predicted IC50 is below the lowest dose ({lo:g} ppm);",
              "  add lower doses (e.g. 1.5625, 3.125, 6.25 ppm) to define the top",
              "  of the curve."]
    else:
        L += ["VERDICT: RANGE ADEQUATE.",
              "  The predicted IC50 falls inside the planned series; the design",
              "  should bracket 50% inhibition. Keep the interference and",
              "  vehicle controls (defects D6, D1) regardless."]
    return "\n".join(L)


def _demo():
    """End-to-end demonstration on a SYNTHETIC composition.

    This proves the pipeline runs before MSU-IIT's data arrives. The numbers
    are NOT a prediction -- the real peak table replaces the composition here.
    """
    pred_path = TAB / "qsar_moringa_predictions.csv"
    if not pred_path.exists():
        print(f"(need {pred_path}; run s06_qsar_train.py first)")
        return
    pot = pd.read_csv(pred_path)
    mo = pot[pot["source"] == "moringa"].copy()

    rng = np.random.default_rng(20260729)
    n = min(18, len(mo))
    pick = mo.sample(n, random_state=20260729).copy()
    w = rng.dirichlet(np.ones(n) * 0.6) * 62.0     # 62% of mass identified
    pick["pct_area"] = w

    comp = pick[["name", "pubchem_cid", "pct_area"]]
    res = predict_extract_ic50(comp, mo, key="pubchem_cid")

    print("=" * 74)
    print("PHASE C DEMONSTRATION -- synthetic composition, real QSAR potencies")
    print("=" * 74)
    print("!! The composition below is SYNTHETIC. Replace with the MSU-IIT")
    print("!! peak table. These numbers are a pipeline test, not a prediction.")
    print()
    print(res.summary())
    print("\nTop drivers of the predicted effect (concentration addition):")
    print(res.driver_table.head(10).to_string(
        index=False, float_format=lambda x: f"{x:9.3f}"))
    lo, hi = res.interval
    print("\n" + "-" * 74)
    print(dose_range_recommendation(lo, hi))


if __name__ == "__main__":
    _demo()
