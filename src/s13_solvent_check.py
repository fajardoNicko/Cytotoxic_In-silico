"""
s13_solvent_check.py  --  Phase F follow-up, defect D1 quantified on real data

The original plan listed "no vehicle control" as defect D1 and called it
critical. The laboratory report shows why.

From the report's own Sample Preparation section:
    2 mg of extract dissolved in 200 uL of neat DMSO  ->  10,000 ug/mL stock
    working concentrations prepared by diluting that stock with culture medium

If the working doses came straight off that stock, the final solvent fraction
in each well is fixed by the dose:

    DMSO %v/v = 100 * dose / 10,000

which puts 4% DMSO in the top-dose well. A549 is routinely held at or below
0.5% v/v, and above roughly 1% the solvent itself kills.

The problem this creates is not that DMSO explains the result. It is that
DMSO and extract concentration are perfectly collinear by construction, so no
statistical treatment of this dataset can separate them. Only a vehicle
control run at the same solvent percentage can, and the raw data contains
only NC and Blank columns.

Outputs
  results/tables/F4_solvent_confound.csv
  results/figures/F4_solvent_confound.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

STOCK_UG_ML = 10_000.0        # 2 mg in 200 uL neat DMSO
A549_ROUTINE_LIMIT = 0.5      # % v/v, standard working ceiling
A549_TOXIC_ABOVE = 1.0        # % v/v, commonly cytotoxic beyond this


def main() -> None:
    rep = pd.read_csv(RAW / "mtt_observed_reported.csv")
    mb = rep[rep["sample"] == "MBEE"].sort_values("conc_ppm").reset_index(drop=True)
    mb["dmso_pct_if_direct"] = 100.0 * mb["conc_ppm"] / STOCK_UG_ML
    mb["over_routine_limit"] = mb["dmso_pct_if_direct"] > A549_ROUTINE_LIMIT
    mb["over_toxic_threshold"] = mb["dmso_pct_if_direct"] > A549_TOXIC_ABOVE
    mb["inhibition_pct"] = 100.0 - mb["viab_mean"]

    print("=" * 78)
    print("DEFECT D1 ON REAL DATA -- solvent and dose are perfectly collinear")
    print("=" * 78)
    print(f"\nStock as described in the report : {STOCK_UG_ML:,.0f} ug/mL in neat DMSO")
    print("Implied final solvent fraction, if working doses were diluted")
    print("directly from that stock:\n")
    print(f"  {'dose':>8} {'DMSO %':>8} {'viability':>10} {'inhibition':>11}  flag")
    for _, r in mb.iterrows():
        flag = ("ABOVE the cytotoxic threshold" if r["over_toxic_threshold"]
                else "above routine limit" if r["over_routine_limit"] else "")
        print(f"  {r['conc_ppm']:8.1f} {r['dmso_pct_if_direct']:8.3f} "
              f"{r['viab_mean']:10.2f} {r['inhibition_pct']:11.2f}  {flag}")

    lo = mb[~mb["over_routine_limit"]]
    hi = mb[mb["over_routine_limit"]]
    print(f"\n  At or below {A549_ROUTINE_LIMIT}% DMSO ({', '.join(f'{d:g}' for d in lo['conc_ppm'])} ug/mL):")
    print(f"    viability {lo['viab_mean'].min():.2f} to {lo['viab_mean'].max():.2f}%, "
          f"every value ABOVE the untreated control")
    print(f"  Above {A549_ROUTINE_LIMIT}% DMSO ({', '.join(f'{d:g}' for d in hi['conc_ppm'])} ug/mL):")
    print(f"    viability {hi['viab_mean'].min():.2f} to {hi['viab_mean'].max():.2f}%, "
          f"every value BELOW it")

    r_, p_ = sps.pearsonr(mb["dmso_pct_if_direct"], mb["viab_mean"])
    print(f"\n  Viability against DMSO %: Pearson r = {r_:.4f}, p = {p_:.4f}")
    print("  This correlation is NOT evidence that DMSO caused the effect.")
    print("  Solvent fraction is an exact linear function of dose here, so the")
    print("  two variables carry identical information and the correlation with")
    print("  dose is necessarily the same number. That is precisely the problem.")
    print("  The design cannot attribute the observed inhibition to the extract,")
    print("  to the solvent, or to any mixture of the two.")

    print("\n  What would resolve it:")
    print("   1. Confirm whether an intermediate dilution in medium was used.")
    print("      If it was, the final DMSO may be far below these figures and")
    print("      the concern largely disappears. The report does not say.")
    print("   2. State the final % v/v DMSO in the top-dose well.")
    print("   3. Supply the vehicle-control absorbances. The report's viability")
    print("      formula names an 'Ave Vehicle Absorbance' term, but the raw")
    print("      data sheet carries only NC and Blank columns.")

    mb.to_csv(TAB / "F4_solvent_confound.csv", index=False)

    # figure
    fig, ax = plt.subplots(figsize=(7.4, 4.7))
    colors = ["#2166ac" if not f else "#b2182b" for f in mb["over_routine_limit"]]
    ax.bar([f"{d:g}" for d in mb["conc_ppm"]], mb["viab_mean"],
           yerr=mb["viab_sd"], color=colors, capsize=4, width=0.62)
    ax.axhline(100, color="black", lw=1.1)
    for i, r in mb.iterrows():
        ax.annotate(f"{r['dmso_pct_if_direct']:.3g}%", (i, 4),
                    ha="center", fontsize=8.5, color="white", weight="bold")
    ax.set_xlabel("extract concentration (µg/mL), DMSO % v/v shown inside each bar")
    ax.set_ylabel("% cell viability")
    ax.set_title("Every dose that inhibited also exceeded the A549 solvent limit")
    ax.set_ylim(0, 150)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color="#2166ac", label=f"DMSO at or below {A549_ROUTINE_LIMIT}% v/v"),
        Patch(color="#b2182b", label=f"DMSO above {A549_ROUTINE_LIMIT}% v/v")],
        fontsize=8.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "F4_solvent_confound.png", dpi=170)
    plt.close(fig)

    print(f"\nWrote F4_solvent_confound.csv -> {TAB}")
    print(f"Wrote F4_solvent_confound.png -> {FIG}")


if __name__ == "__main__":
    main()
