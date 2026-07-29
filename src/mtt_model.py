"""
mtt_model.py  --  dose-response mathematics for the MTT assay

Implements, side by side:
  (1) the research plan's own formulas, verbatim -- so the in silico study
      speaks the manuscript's language; and
  (2) a 4-parameter logistic (4PL) nonlinear fit -- the defensible standard,
      used as the sensitivity analysis required by defect D5.

Manuscript formulas (pp. 10-11):
    %Viability    = [(Abs_sample - Abs_blank) / (Abs_control - Abs_blank)] * 100
    %Cytotoxicity = 100 - %Viability
    IC50          : regress %viability on log10(concentration), solve y = mx + b
                    at y = 50  ->  IC50 = 10 ** ((50 - b) / m)

The log-linear IC50 is only valid when the solution lies INSIDE the tested
concentration range.  `LogLinearIC50.extrapolated` flags when it does not --
this is defect D2 of the plan, and it is the single most likely way the real
experiment fails to answer Specific Question 3.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# the manuscript's formulas
# ---------------------------------------------------------------------------

def percent_viability(abs_sample, abs_blank, abs_control):
    """%Viability, exactly as written on p. 10 of the research plan."""
    return (np.asarray(abs_sample, float) - abs_blank) / (abs_control - abs_blank) * 100.0


def percent_cytotoxicity(viability_pct):
    """%Cytotoxicity = 100 - %Viability  (p. 11)."""
    return 100.0 - np.asarray(viability_pct, float)


# ---------------------------------------------------------------------------
# dose-response models
# ---------------------------------------------------------------------------

def hill_viability(conc, ic50, hill=1.0, v_min=0.0, v_max=100.0):
    """4PL / Hill curve for %viability as a function of concentration.

    v_min is the bottom plateau (residual viability at infinite dose);
    a crude plant extract typically bottoms out at 5-20%, not 0.
    """
    conc = np.asarray(conc, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return v_min + (v_max - v_min) / (1.0 + (conc / ic50) ** hill)


@dataclass
class LogLinearIC50:
    """IC50 by the manuscript's method: linear fit of %viability on log10(C)."""
    slope: float
    intercept: float
    r_squared: float
    p_value: float
    std_err: float
    log_ic50: float
    ic50: float
    conc_min: float
    conc_max: float
    extrapolated: bool
    n_points: int
    note: str = ""

    @property
    def equation(self) -> str:
        return f"y = {self.slope:.4f}x + {self.intercept:.4f}   (R^2 = {self.r_squared:.4f})"


def loglinear_ic50(concentrations, viability_pct) -> LogLinearIC50:
    """The research plan's IC50 method (p. 11).

    Regress %viability (y) on log10(concentration) (x), then solve for y = 50.
    Concentrations must be > 0; the untreated control (C = 0) is excluded, as
    log10(0) is undefined -- this is standard and matches how the manuscript's
    source study is done.
    """
    c = np.asarray(concentrations, float)
    v = np.asarray(viability_pct, float)
    keep = c > 0
    c, v = c[keep], v[keep]
    if c.size < 3:
        raise ValueError("need >= 3 non-zero concentrations for the regression")

    x = np.log10(c)
    reg = stats.linregress(x, v)

    note = ""
    if abs(reg.slope) < 1e-12:
        log_ic50, ic50 = np.nan, np.nan
        note = "slope ~ 0: no dose-response, IC50 undefined"
    else:
        log_ic50 = (50.0 - reg.intercept) / reg.slope
        ic50 = float(10.0 ** log_ic50)

    cmin, cmax = float(c.min()), float(c.max())
    extrap = not (cmin <= ic50 <= cmax) if np.isfinite(ic50) else True
    if extrap and np.isfinite(ic50):
        where = "ABOVE" if ic50 > cmax else "BELOW"
        note = (f"EXTRAPOLATION: IC50 = {ic50:.1f} lies {where} the tested range "
                f"[{cmin:g}, {cmax:g}]. Not a valid determination.")

    return LogLinearIC50(
        slope=float(reg.slope), intercept=float(reg.intercept),
        r_squared=float(reg.rvalue ** 2), p_value=float(reg.pvalue),
        std_err=float(reg.stderr), log_ic50=float(log_ic50), ic50=ic50,
        conc_min=cmin, conc_max=cmax, extrapolated=bool(extrap),
        n_points=int(c.size), note=note,
    )


@dataclass
class FourPLFit:
    ic50: float                       # 4PL midpoint parameter ("relative IC50")
    hill: float
    v_min: float
    v_max: float
    r_squared: float
    ic50_ci: tuple[float, float] | None
    converged: bool
    extrapolated: bool
    note: str = ""

    @property
    def absolute_ic50(self) -> float:
        """Concentration at which %viability actually equals 50.

        This -- not the midpoint parameter -- is what the research plan means
        by IC50 ("the concentration required to inhibit 50% of cell viability",
        p. 11).  The two coincide only when the bottom plateau is 0 and the top
        is 100.  With a residual-viability plateau of v_min they differ by a
        factor of [ (v_max - 50) / (50 - v_min) ] ** (1/hill).

        Returns NaN when the curve never crosses 50% (v_min >= 50), which is
        itself the correct answer: no IC50 exists.
        """
        if not np.isfinite(self.ic50) or self.v_min >= 50.0 or self.v_max <= 50.0:
            return np.nan
        ratio = (self.v_max - 50.0) / (50.0 - self.v_min)
        return float(self.ic50 * ratio ** (1.0 / self.hill))


def fourpl_ic50(concentrations, viability_pct,
                fix_top: float | None = 100.0) -> FourPLFit:
    """4PL nonlinear fit -- the sensitivity analysis demanded by defect D5.

    fix_top=100 anchors the upper plateau at the normalised control (100%),
    which is the right constraint when viability is already control-normalised
    and there is no supra-control data to define a free top.
    """
    c = np.asarray(concentrations, float)
    v = np.asarray(viability_pct, float)
    keep = c > 0
    c, v = c[keep], v[keep]

    cmin, cmax = float(c.min()), float(c.max())
    p0_ic50 = float(np.sqrt(cmin * cmax))

    try:
        # The bottom plateau may legitimately sit well above 50% when the
        # extract is weak and the dose range never reaches half-inhibition
        # (defect D2). Allow that, rather than letting the fit fail.
        if fix_top is None:
            def f(x, ic50, hill, vmin, vmax):
                return hill_viability(x, ic50, hill, vmin, vmax)
            p0 = [p0_ic50, 1.0, max(v.min() - 5, 0.0), 100.0]
            bounds = ([cmin * 1e-3, 0.1, 0.0, 50.0], [cmax * 1e3, 10.0, 95.0, 150.0])
        else:
            def f(x, ic50, hill, vmin):
                return hill_viability(x, ic50, hill, vmin, fix_top)
            p0 = [p0_ic50, 1.0, max(v.min() - 5, 0.0)]
            bounds = ([cmin * 1e-3, 0.1, 0.0], [cmax * 1e3, 10.0, 95.0])

        # curve_fit rejects a p0 outside the bounds; clip it in.
        lo, hi = np.asarray(bounds[0], float), np.asarray(bounds[1], float)
        p0 = np.clip(np.asarray(p0, float), lo, hi)

        popt, pcov = curve_fit(f, c, v, p0=p0, bounds=bounds, maxfev=20000)
        pred = f(c, *popt)
        ss_res = float(((v - pred) ** 2).sum())
        ss_tot = float(((v - v.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

        ic50 = float(popt[0])
        hill = float(popt[1])
        vmin = float(popt[2])
        vmax = float(popt[3]) if fix_top is None else float(fix_top)

        se = float(np.sqrt(np.diag(pcov))[0]) if np.all(np.isfinite(pcov)) else np.nan
        ci = (ic50 - 1.96 * se, ic50 + 1.96 * se) if np.isfinite(se) else None

        extrap = not (cmin <= ic50 <= cmax)
        note = ""
        if extrap:
            where = "ABOVE" if ic50 > cmax else "BELOW"
            note = (f"EXTRAPOLATION: 4PL IC50 = {ic50:.1f} lies {where} the tested "
                    f"range [{cmin:g}, {cmax:g}].")
        if vmin > v.min() + 1e-6 and v.min() > 50:
            note += " Curve never reaches 50% viability in the tested range."

        return FourPLFit(ic50, hill, vmin, vmax, r2, ci, True, extrap, note)

    except Exception as exc:                                    # noqa: BLE001
        return FourPLFit(np.nan, np.nan, np.nan, np.nan, np.nan, None,
                         False, True, f"fit failed: {exc}")


# ---------------------------------------------------------------------------
# unit conversion helpers (QSAR gives uM; the assay is dosed in ppm = ug/mL)
# ---------------------------------------------------------------------------

def um_to_ugml(ic50_um, mw_g_per_mol):
    """IC50 [uM] -> IC50 [ug/mL].  ug/mL = uM * MW / 1000."""
    return np.asarray(ic50_um, float) * np.asarray(mw_g_per_mol, float) / 1000.0


def ugml_to_um(ic50_ugml, mw_g_per_mol):
    return np.asarray(ic50_ugml, float) * 1000.0 / np.asarray(mw_g_per_mol, float)


def pic50_to_ic50_um(pic50):
    """pIC50 = -log10(IC50 in M)  ->  IC50 in uM."""
    return 10.0 ** (6.0 - np.asarray(pic50, float))


def ic50_um_to_pic50(ic50_um):
    return 6.0 - np.log10(np.asarray(ic50_um, float))


# ---------------------------------------------------------------------------
# NCI crude-extract activity classification
# ---------------------------------------------------------------------------

def nci_activity_class(ic50_ugml: float) -> str:
    """Common crude-extract cytotoxicity bands (US NCI / Geran criteria).

    Thresholds vary between sources; the widely cited NCI cut-off for a crude
    extract to be considered 'active' is IC50 <= 20-30 ug/mL. Cite whichever
    source the manuscript adopts and keep it consistent.
    """
    if not np.isfinite(ic50_ugml):
        return "indeterminate"
    if ic50_ugml <= 20:
        return "active (IC50 <= 20 ug/mL)"
    if ic50_ugml <= 100:
        return "moderately active (20-100 ug/mL)"
    if ic50_ugml <= 500:
        return "weakly active (100-500 ug/mL)"
    if ic50_ugml <= 1000:
        return "very weak (500-1000 ug/mL)"
    return "inactive (> 1000 ug/mL)"


if __name__ == "__main__":
    print("=" * 74)
    print("mtt_model.py self-test")
    print("=" * 74)

    doses = np.array([12.5, 25, 50, 100, 200])

    for label, true_ic50 in [("IC50 inside range", 80.0),
                             ("IC50 above range (defect D2)", 450.0)]:
        v = hill_viability(doses, true_ic50, hill=1.2, v_min=8.0)
        ll = loglinear_ic50(doses, v)
        fp = fourpl_ic50(doses, v)
        print(f"\n--- {label}: true IC50 = {true_ic50} ug/mL ---")
        print(f"  viability at {doses.tolist()} = "
              f"{np.round(v, 1).tolist()}")
        print(f"  log-linear (manuscript): IC50 = {ll.ic50:8.2f}   {ll.equation}")
        if ll.note:
            print(f"      !! {ll.note}")
        print(f"  4PL       (sensitivity): IC50 = {fp.ic50:8.2f}   "
              f"hill = {fp.hill:.2f}  R^2 = {fp.r_squared:.4f}")
        if fp.note:
            print(f"      !! {fp.note.strip()}")
        err = abs(ll.ic50 - fp.ic50) / fp.ic50 * 100
        print(f"  disagreement between methods: {err:.1f}%  "
              f"(> 20% is itself a finding, per D5)")
        print(f"  NCI class (4PL): {nci_activity_class(fp.ic50)}")
