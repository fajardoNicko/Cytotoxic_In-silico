"""
s10_preregister.py  --  Phase E of IN_SILICO_PLAN.md

Freezes the study's predictions BEFORE any wet-lab data exists, and seals the
file with a SHA-256 hash. Without this step, "our in silico results agreed with
our in vitro results" is unfalsifiable.

Writes:
  results/prediction_registry/prediction_v1.json   the sealed prediction
  results/prediction_registry/REGISTRY.md          hash + timestamp ledger

Rules once sealed:
  * never edit prediction_v1.json
  * any revision becomes prediction_v2.json with its own hash and a written
    justification
  * email the hash to the adviser (and optionally post to OSF) the same day
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

from mtt_model import hill_viability
from s07b_literature_prior import build_priors, p_exceeds, DOSES_PLAN

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "results" / "tables"
REG = ROOT / "results" / "prediction_registry"
REG.mkdir(parents=True, exist_ok=True)

RECOMMENDED_SERIES = [50, 100, 200, 400, 800, 1600]


def git_commit() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or None
    except Exception:                                           # noqa: BLE001
        return None


def recommended_n() -> dict:
    """Pull the replicate recommendation from the Monte Carlo, if available."""
    f = TAB / "E2_power.csv"
    if not f.exists():
        return {"available": False,
                "note": "E2_power.csv not present; run s08_monte_carlo.py"}
    d = pd.read_csv(f)
    out = {"available": True, "source": "E2_power.csv"}
    dose_cols = [c for c in d.columns if c.startswith("power_") and c != "power_anova"]
    # smallest configuration reaching 80% power at every dose, for each scenario
    for ic50, sub in d.groupby("true_ic50"):
        ok = sub[(sub[dose_cols] >= 0.80).all(axis=1)]
        ok = ok.sort_values("wells_per_group")
        key = f"true_ic50_{ic50:g}"
        if len(ok):
            r = ok.iloc[0]
            out[key] = {
                "n_tech": int(r["n_tech"]), "n_bio": int(r["n_bio"]),
                "wells_per_group": int(r["wells_per_group"]),
                "perimeter_required": bool(r.get("perimeter_required", False)),
                "power_anova": float(r["power_anova"]),
            }
        else:
            best = sub.sort_values("wells_per_group").iloc[-1]
            out[key] = {
                "achievable": False,
                "largest_tested_wells_per_group": int(best["wells_per_group"]),
                "min_dose_power": float(best[dose_cols].min()),
                "note": "no tested configuration reaches 80% power at every dose",
            }
    return out


def main():
    pri = build_priors()
    p_mu, p_sd = pri["pooled"]
    e_mu, e_sd = pri["solvent_matched"]

    def curve(ic50, doses):
        v = hill_viability(np.array(doses, float), ic50, hill=1.2, v_min=8.0)
        return {str(d): round(float(x), 1) for d, x in zip(doses, v)}

    ic50_e, ic50_p = 10 ** e_mu, 10 ** p_mu

    pred = {
        "study": "In silico prediction: M. oleifera ethanolic bark extract vs A549",
        "companion_plan": "IN_SILICO_PLAN.md",
        "sealed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "python": sys.version.split()[0],

        "basis": {
            "method": "literature prior on crude-extract IC50, solvent-matched",
            "why_not_qsar": (
                "The A549 QSAR failed Gate G2 on a scaffold split "
                "(test R^2 = 0.562 < 0.60) and showed ~11-fold 1-sigma error "
                "with an unresolved bias direction on natural products, so "
                "per-compound potencies were NOT summed into an extract IC50."),
            "literature_anchors_ug_per_mL": {
                "aqueous_leaf_A549": 166.7,
                "alkaloid_extract_A549": 158.67,
                "ethanolic_leaf_A549": 1062.87,
                "ethanolic_leaf_MCF12A_normal": 1424.04,
            },
        },

        "PREDICTION_1_extract_ic50_ug_per_mL": {
            "point_solvent_matched": round(ic50_e, 1),
            "interval_68pct_solvent_matched": [round(10 ** (e_mu - e_sd), 1),
                                               round(10 ** (e_mu + e_sd), 1)],
            "interval_95pct_solvent_matched": [round(10 ** (e_mu - 1.96 * e_sd), 1),
                                               round(10 ** (e_mu + 1.96 * e_sd), 1)],
            "point_pooled_optimistic": round(ic50_p, 1),
            "interval_68pct_pooled": [round(10 ** (p_mu - p_sd), 1),
                                      round(10 ** (p_mu + p_sd), 1)],
        },

        "PREDICTION_2_viability_at_planned_doses": {
            "doses_ppm": DOSES_PLAN,
            "solvent_matched_prior": curve(ic50_e, DOSES_PLAN),
            "pooled_prior": curve(ic50_p, DOSES_PLAN),
            "max_inhibition_at_200ppm_solvent_matched_pct": round(
                100 - hill_viability(200.0, ic50_e, 1.2, 8.0), 1),
        },

        "PREDICTION_3_dose_range_failure": {
            "claim": ("The planned 12.5-200 ppm series will NOT bracket 50% "
                      "inhibition, so no valid IC50 can be determined from it."),
            "p_ic50_exceeds_200ppm_solvent_matched": round(
                p_exceeds(e_mu, e_sd, 200.0), 3),
            "p_ic50_exceeds_200ppm_pooled": round(
                p_exceeds(p_mu, p_sd, 200.0), 3),
            "falsified_if": ("observed mean viability at 200 ppm is below 50%, "
                             "i.e. the series does bracket the IC50"),
        },

        "PREDICTION_4_dose_rank_order": {
            "claim": "monotonic decrease in mean viability across ascending dose",
            "expected_spearman_rho": -1.0,
        },

        "PREDICTION_5_statistics_under_planned_design": {
            "claim_solvent_matched": (
                "ANOVA may still be significant on the trend, but Tukey will "
                "NOT separate 12.5 or 25 ppm from the negative control, because "
                "predicted inhibition there is under 2%."),
            "predicted_tukey_vs_control_significant": {
                "12.5ppm": False, "25ppm": False, "50ppm": False,
                "100ppm": False, "200ppm": "borderline",
            },
        },

        "PREDICTION_6_selectivity_index": {
            "estimate": 1.34,
            "basis": "MCF-12A / A549 from the ethanolic-leaf study",
            "interpretation": ("close to 1 -- little selectivity for cancer over "
                              "normal cells; a therapeutically weak result"),
        },

        "PREDICTION_7_recommended_design": {
            "dose_series_ppm": RECOMMENDED_SERIES,
            "drop": [12.5, 25],
            "reason": "both priors predict <5% inhibition; those wells buy no information",
            "replicates": recommended_n(),
            "required_controls": [
                "vehicle control at the top-dose solvent % (defect D1)",
                "cell-free extract+MTT interference control at every dose (D6)",
                "normal lung comparator for Selectivity Index (D8)",
            ],
            "layout_constraint": ("with the full control set, n_tech <= 6 fits the "
                                  "60 interior wells; n_tech = 8 forces use of "
                                  "evaporation-prone perimeter wells"),
        },

        "PREDICTION_8_method_bias": {
            "claim": ("The manuscript's log-linear IC50 will disagree with a 4PL "
                      "fit, and the disagreement grows without bound once the "
                      "IC50 leaves the tested range."),
            "monte_carlo_loglinear_bias_pct": {
                "true_ic50_107": 0.4, "true_ic50_161": 15.5, "true_ic50_214": 46.5,
                "true_ic50_322": 182.0, "true_ic50_482": 889.1,
                "true_ic50_750": 12159.8, "true_ic50_1072": 382681.0,
            },
            "warning": ("a high regression R^2 does NOT indicate an adequate dose "
                        "range: at a true IC50 of 1063 ug/mL the log-linear fit "
                        "gives R^2 = 0.85 while mis-estimating IC50 by ~16,000x"),
        },

        "VALIDATION_CRITERIA": {
            "primary": "predicted IC50 within 3-fold of observed, OR the "
                       "observed result confirms IC50 > top dose as predicted",
            "viability_rmse_threshold_pct_points": 15,
            "spearman_rho_threshold": 0.9,
            "note": ("If the extract turns out markedly more potent than the "
                     "ethanolic-leaf prior -- plausible, since bark chemistry "
                     "differs and no bark A549 data exists -- this prediction "
                     "fails and that is a real finding about bark, not a "
                     "modelling excuse."),
        },

        "KNOWN_LIMITATIONS": [
            "No quantitative A549 data for M. oleifera BARK exists; prior built "
            "from leaf extracts.",
            "One anchor used MTS rather than MTT.",
            "The single ethanolic anchor borrows its spread from the between-study "
            "variance of the other anchors.",
            "QSAR unusable for per-compound potency; docking/ranking only.",
            "ChEMBL pull was partial (connection dropped at 21,000 records).",
            "Thioredoxin reductase 1 has no structure with a drug-like ligand; "
            "target dropped from the docking panel.",
        ],
    }

    path = REG / "prediction_v1.json"
    if path.exists():
        print(f"!! {path.name} already exists and must NOT be overwritten.")
        print("   Create prediction_v2.json with a written justification instead.")
        return

    blob = json.dumps(pred, indent=2, sort_keys=False)
    path.write_text(blob, encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()

    ledger = REG / "REGISTRY.md"
    header = ("# Prediction registry\n\n"
              "Sealed predictions, newest last. Once a row is added, the named\n"
              "file must never be modified. A revision gets a new file, a new\n"
              "hash, and a stated reason.\n\n"
              "| File | SHA-256 | Sealed (UTC) | Git commit | Note |\n"
              "|---|---|---|---|---|\n")
    if not ledger.exists():
        ledger.write_text(header, encoding="utf-8")
    with ledger.open("a", encoding="utf-8") as fh:
        fh.write(f"| `{path.name}` | `{digest}` | {pred['sealed_utc']} | "
                 f"`{pred['git_commit'] or 'not a git repo'}` | "
                 f"initial pre-registration, literature-prior basis |\n")

    print("=" * 78)
    print("PRE-REGISTRATION SEALED")
    print("=" * 78)
    print(f"  file      : {path}")
    print(f"  SHA-256   : {digest}")
    print(f"  sealed    : {pred['sealed_utc']}")
    print(f"  ledger    : {ledger}")
    print()
    print("  HEADLINE PREDICTIONS")
    print(f"    extract IC50 (solvent-matched) = {ic50_e:.0f} ug/mL "
          f"[68%: {10 ** (e_mu - e_sd):.0f}-{10 ** (e_mu + e_sd):.0f}]")
    print(f"    P(planned 200 ppm range fails) = "
          f"{p_exceeds(e_mu, e_sd, 200.0) * 100:.1f}%")
    print(f"    max inhibition at 200 ppm      = "
          f"{100 - hill_viability(200.0, ic50_e, 1.2, 8.0):.1f}%")
    print(f"    Selectivity Index              = 1.34")
    print()
    print("  NEXT: email this hash to your adviser today, before any assay data")
    print("        exists. That timestamp is what makes the prediction falsifiable.")


if __name__ == "__main__":
    main()
