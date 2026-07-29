"""
virtual_plate.py  --  Phase D of IN_SILICO_PLAN.md

Simulates the 96-well MTT experiment described in the research plan, well by
well, including the error sources that actually determine whether the real
experiment succeeds:

  * pipetting / seeding variation      (well-level CV)
  * plate-to-plate offset              (between independent biological runs)
  * edge evaporation                   (perimeter wells read high)
  * direct MTT reduction by polyphenols (defect D6) -- extract reduces the
    tetrazolium with NO cells present, inflating apparent viability
  * solvent (DMSO) cytotoxicity        (defect D1) -- and, critically, the
    BIAS introduced by normalising to a medium-only control instead of a
    vehicle control

Absorbance model for one well:

    A = A_blank
        + (A_ctrl - A_blank) * V(C)/100 * (1 + eps_well) * plate_gain   [cells only]
        + k_interference * C                                            [extract present]
        + edge_bias                                                     [perimeter only]

where V(C) is the 4PL viability curve from mtt_model.hill_viability.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from mtt_model import (hill_viability, percent_viability, loglinear_ic50,
                       fourpl_ic50)

ROWS = "ABCDEFGH"
NCOLS = 12
DOSES_PLAN = [12.5, 25.0, 50.0, 100.0, 200.0]      # the manuscript's series


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

@dataclass
class AssayConfig:
    """Everything the wet lab does not specify (defect D7) is pinned here."""

    # --- biology (the quantity Phase C predicts) ---
    true_ic50: float = 150.0          # ug/mL, extract vs A549
    hill: float = 1.2                 # Hill slope
    v_min: float = 8.0                # bottom plateau, % viability

    # --- doses ---
    doses: list[float] = field(default_factory=lambda: list(DOSES_PLAN))

    # --- positive control ---
    dox_ic50: float = 1.0             # ug/mL, doxorubicin vs A549 (48 h)
    dox_dose: float = 1.0             # single reference dose
    dox_hill: float = 1.0
    dox_v_min: float = 5.0

    # --- replication ---
    n_tech: int = 3                   # technical replicate wells per group
    n_bio: int = 3                    # independent biological runs

    # --- signal ---
    abs_blank: float = 0.05           # medium + MTT, no cells
    abs_control: float = 1.00         # untreated A549, 48 h

    # --- noise ---
    well_cv: float = 0.06             # 6% well-to-well
    plate_cv: float = 0.08            # 8% run-to-run gain
    read_noise: float = 0.005         # plate-reader OD noise
    edge_bias: float = 0.08           # extra OD on perimeter wells
    use_perimeter: bool = False       # True = naive layout, False = recommended

    # --- interference (defect D6) ---
    interference_k: float = 0.0       # OD units per ug/mL of extract, no cells

    # --- vehicle (defect D1) ---
    dmso_pct: float = 0.5             # final % v/v in every treated well
    dmso_viability: float = 96.0      # % viability of vehicle control vs medium

    seed: int | None = None


# ---------------------------------------------------------------------------
# plate layout
# ---------------------------------------------------------------------------

def _well_names(use_perimeter: bool) -> list[tuple[str, int, bool]]:
    out = []
    for r in ROWS:
        for c in range(1, NCOLS + 1):
            edge = (r in "AH") or (c in (1, NCOLS))
            if edge and not use_perimeter:
                continue
            out.append((r, c, edge))
    return out


def build_layout(cfg: AssayConfig) -> pd.DataFrame:
    """Assign groups to physical wells.

    Groups, in the order they are laid out:
      blank                 medium + MTT, no cells           (n_tech)
      negative_control      cells + medium                   (2 * n_tech)
      vehicle_control       cells + medium + DMSO            (2 * n_tech)  [D1]
      doxorubicin           cells + DMSO + dox               (2 * n_tech)
      <dose> ppm            cells + DMSO + extract           (2 * n_tech each)
      interference_<dose>   medium + MTT + extract, NO cells (1 each)      [D6]
    """
    spec: list[tuple[str, float, bool, bool]] = []   # group, conc, has_cells, has_dmso

    spec += [("blank", 0.0, False, False)] * max(cfg.n_tech, 3)
    spec += [("negative_control", 0.0, True, False)] * cfg.n_tech
    spec += [("vehicle_control", 0.0, True, True)] * cfg.n_tech
    spec += [("doxorubicin", 0.0, True, True)] * cfg.n_tech
    for d in cfg.doses:
        spec += [(f"{d:g}ppm", d, True, True)] * cfg.n_tech
    for d in cfg.doses:
        spec += [(f"interference_{d:g}ppm", d, False, True)]

    slots = _well_names(cfg.use_perimeter)
    used_perimeter = cfg.use_perimeter
    if len(spec) > len(slots) and not cfg.use_perimeter:
        # Real constraint, not a coding limit: the full control set (blank,
        # negative, vehicle, doxorubicin, 5 interference wells) plus 5 doses
        # stops fitting the 60 interior wells at n_tech >= 8. The lab's only
        # single-plate option then is to use the evaporation-prone perimeter.
        slots = _well_names(True)
        used_perimeter = True

    if len(spec) > len(slots):
        raise ValueError(
            f"layout needs {len(spec)} wells but a 96-well plate has "
            f"{len(slots)}. n_tech={cfg.n_tech} with {len(cfg.doses)} doses "
            f"does not fit one plate; split across plates or reduce n_tech.")

    rows = []
    for (grp, conc, cells, dmso), (r, c, edge) in zip(spec, slots):
        rows.append({"well": f"{r}{c:02d}", "row": r, "col": c, "is_edge": edge,
                     "group": grp, "conc_ppm": conc,
                     "has_cells": cells, "has_dmso": dmso,
                     "perimeter_required": used_perimeter})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# simulation
# ---------------------------------------------------------------------------

def _true_viability(cfg: AssayConfig, group: str, conc: float,
                    has_dmso: bool) -> float:
    """Ground-truth % viability, before any measurement error."""
    if group == "doxorubicin":
        v = hill_viability(cfg.dox_dose, cfg.dox_ic50, cfg.dox_hill, cfg.dox_v_min)
    elif conc > 0:
        v = hill_viability(conc, cfg.true_ic50, cfg.hill, cfg.v_min)
    else:
        v = 100.0
    # solvent cytotoxicity multiplies on top of the extract effect (D1)
    if has_dmso:
        v *= cfg.dmso_viability / 100.0
    return float(v)


def true_absolute_ic50(cfg: AssayConfig,
                       normalize_to: str = "negative_control") -> float:
    """Ground-truth concentration at which measured %viability equals 50.

    This is the estimand: the quantity the manuscript's IC50 is trying to
    recover.  It is NOT cfg.true_ic50 -- that is the 4PL midpoint parameter.
    They differ because (a) the bottom plateau v_min is above 0, and (b) when
    normalising to a medium-only negative control, the DMSO cytotoxicity of
    the vehicle is folded into the apparent extract effect (defect D1).

    Solved numerically so the DMSO factor is handled exactly.
    """
    from scipy.optimize import brentq

    dmso_factor = 1.0
    if normalize_to == "negative_control":
        dmso_factor = cfg.dmso_viability / 100.0

    def f(c):
        return hill_viability(c, cfg.true_ic50, cfg.hill, cfg.v_min) * dmso_factor - 50.0

    lo, hi = 1e-6, 1e9
    if f(lo) < 0 or f(hi) > 0:
        return float("nan")          # curve never crosses 50%
    return float(brentq(f, lo, hi, xtol=1e-9, rtol=1e-12))


def simulate_plate(cfg: AssayConfig, run_id: int, rng) -> pd.DataFrame:
    """One independent biological run = one plate."""
    lay = build_layout(cfg).copy()
    plate_gain = float(rng.normal(1.0, cfg.plate_cv))

    abs_vals, true_v = [], []
    for _, w in lay.iterrows():
        v = _true_viability(cfg, w["group"], w["conc_ppm"], w["has_dmso"])
        true_v.append(v)

        a = cfg.abs_blank
        if w["has_cells"]:
            eps = float(rng.normal(0.0, cfg.well_cv))
            a += (cfg.abs_control - cfg.abs_blank) * (v / 100.0) * (1 + eps) * plate_gain
        # direct MTT reduction by the extract itself -- happens with or without cells
        a += cfg.interference_k * w["conc_ppm"]
        if w["is_edge"]:
            a += cfg.edge_bias
        a += float(rng.normal(0.0, cfg.read_noise))
        abs_vals.append(max(a, 0.0))

    lay["run"] = run_id
    lay["plate_gain"] = plate_gain
    lay["true_viability"] = true_v
    lay["abs_570"] = abs_vals
    return lay


def simulate_experiment(cfg: AssayConfig) -> pd.DataFrame:
    rng = np.random.default_rng(cfg.seed)
    return pd.concat([simulate_plate(cfg, i + 1, rng) for i in range(cfg.n_bio)],
                     ignore_index=True)


# ---------------------------------------------------------------------------
# analysis -- exactly the manuscript's pipeline
# ---------------------------------------------------------------------------

def compute_viability(raw: pd.DataFrame,
                      normalize_to: str = "negative_control",
                      subtract_interference: bool = False) -> pd.DataFrame:
    """Apply the research plan's %viability formula (p. 10), per run.

    normalize_to:
        "negative_control" -- what the manuscript specifies (medium only)
        "vehicle_control"  -- the methodologically correct choice (defect D1)

    subtract_interference: use the cell-free extract wells to remove direct
    MTT reduction before normalising (defect D6).
    """
    out = []
    for run, sub in raw.groupby("run"):
        blank = sub.loc[sub["group"] == "blank", "abs_570"].mean()
        ctrl = sub.loc[sub["group"] == normalize_to, "abs_570"].mean()

        interf = {}
        if subtract_interference:
            for _, w in sub[sub["group"].str.startswith("interference_")].iterrows():
                interf[w["conc_ppm"]] = w["abs_570"] - blank

        s = sub[~sub["group"].str.startswith("interference_")
                & (sub["group"] != "blank")].copy()
        adj = s["abs_570"] - s["conc_ppm"].map(lambda c: interf.get(c, 0.0))
        s["viability_pct"] = percent_viability(adj, blank, ctrl)
        s["cytotoxicity_pct"] = 100.0 - s["viability_pct"]
        s["abs_blank_used"] = blank
        s["abs_control_used"] = ctrl
        out.append(s)
    return pd.concat(out, ignore_index=True)


def dose_response_summary(viab: pd.DataFrame) -> pd.DataFrame:
    """mean +/- SD per dose, pooled across runs -- the manuscript's Table."""
    d = viab[viab["conc_ppm"] > 0]
    d = d[~d["group"].isin(["doxorubicin"])]
    g = d.groupby("conc_ppm")["viability_pct"]
    return pd.DataFrame({"conc_ppm": g.mean().index,
                         "n": g.count().values,
                         "mean_viability": g.mean().values,
                         "sd": g.std(ddof=1).values,
                         "mean_inhibition": 100 - g.mean().values}).reset_index(drop=True)


def estimate_ic50(viab: pd.DataFrame):
    """Both IC50 methods on the simulated data."""
    s = dose_response_summary(viab)
    ll = loglinear_ic50(s["conc_ppm"], s["mean_viability"])
    fp = fourpl_ic50(s["conc_ppm"], s["mean_viability"])
    return ll, fp, s


if __name__ == "__main__":
    import stats_mirror as sm

    print("=" * 78)
    print("virtual_plate.py self-test -- one full simulated experiment")
    print("=" * 78)

    cfg = AssayConfig(true_ic50=150.0, seed=20260729)
    lay = build_layout(cfg)
    print(f"\nLayout: {len(lay)} wells used of "
          f"{'96' if cfg.use_perimeter else '60 interior'} available "
          f"(perimeter {'used' if cfg.use_perimeter else 'left empty'})")
    print(lay.groupby("group", sort=False).size().to_string())

    raw = simulate_experiment(cfg)
    print(f"\nSimulated {cfg.n_bio} runs x {len(lay)} wells = {len(raw)} wells")
    print(f"Plate gains: {raw.groupby('run')['plate_gain'].first().round(4).tolist()}")

    viab = compute_viability(raw, normalize_to="negative_control")
    print("\nDose-response (mean +/- SD, %viability):")
    print(dose_response_summary(viab).to_string(
        index=False, float_format=lambda x: f"{x:8.2f}"))

    ll, fp, _ = estimate_ic50(viab)
    print(f"\nlog-linear IC50 (manuscript) = {ll.ic50:8.2f} ug/mL   {ll.equation}")
    if ll.note:
        print(f"   !! {ll.note}")
    print(f"4PL IC50        (sensitivity) = {fp.ic50:8.2f} ug/mL   R^2 = {fp.r_squared:.4f}")
    print(f"true IC50 used to generate    = {cfg.true_ic50:8.2f} ug/mL")

    stat_df = viab[viab["group"].isin(
        ["negative_control", "vehicle_control", "doxorubicin"]
        + [f"{d:g}ppm" for d in cfg.doses])]
    rep = sm.run_analysis(stat_df, "group", "viability_pct")
    print("\n" + sm.format_report(rep, "STATISTICS ON THE SIMULATED EXPERIMENT"))
