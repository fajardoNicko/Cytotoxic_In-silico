# In Silico Study — *Moringa oleifera* Bark Extract vs A549

Computational companion to the research plan *"In Vitro Cytotoxic Activity of Malunggay's (Moringa oleifera) Ethanolic Bark Extract against A549 Lung Cancer Cells."*

The full design rationale is in **[IN_SILICO_PLAN.md](IN_SILICO_PLAN.md)**. This file is the operational entry point: what exists, how to run it, and what is still blocked.

---

## Status

**Phase A (extract-agnostic) — built and running.** Everything here works without the extract data.

| Component | File | State |
|---|---|---|
| Statistics mirror (Shapiro–Wilk → ANOVA/KW → Tukey/Dunn) | `src/stats_mirror.py` | ✅ Gate G3(a) passed |
| Cross-validation vs statsmodels/scipy | `src/test_stats_mirror.py` | ✅ all checks pass |
| jamovi verification package | `src/s09_jamovi_export.py` | ✅ generated, **needs manual jamovi run** |
| Dose–response maths (manuscript + 4PL) | `src/mtt_model.py` | ✅ |
| Virtual 96-well plate simulator | `src/virtual_plate.py` | ✅ |
| Monte Carlo (power, D1/D2/D6 bias) | `src/s08_monte_carlo.py` | ⚠️ E1 done; **stopped mid-E2** |
| Ligand library (47 compounds, 3D) | `src/s01_ligand_prep.py` | ✅ |
| Target panel verification (RCSB) | `src/s02_target_prep.py` | ✅ 15/16 accepted |
| ChEMBL A549 curation (6,972 compounds) | `src/s05_chembl_fetch.py` | ✅ |
| QSAR training + Gates G2/G4 | `src/s06_qsar_train.py` | ⚠️ **G4 pass, G2 FAIL (scaffold)** |
| QSAR reality check | `src/s06b_qsar_reality_check.py` | ⏳ **written, not yet run** |
| Mixture model (Phase C) | `src/s07_mixture_model.py` | ✅ code ready, needs composition |

**Blocked on MSU-IIT:** everything downstream of the GC-MS/LC-MS peak table. See **[docs/MSU_IIT_DATA_REQUEST.md](docs/MSU_IIT_DATA_REQUEST.md)** — send this now.

---

## ▶ Resume here (next session)

Work stopped 2026-07-29. Nothing is left running.

1. **Run the QSAR reality check first — it gates everything downstream.**
   ```bash
   cd src && python -u s06b_qsar_reality_check.py
   ```
   Gate G2 **failed** on the scaffold split (XGBoost test R² = 0.562 vs the
   pre-declared ≥ 0.60), while passing on the random split (0.712). Since the
   *Moringa* constituents are novel scaffolds, the scaffold number is the one
   that counts. The predictions also look compressed — every constituent lands
   at pIC₅₀ 4.6–5.5, and oleic acid at 0.99 µg/mL is not a credible potency for
   a fatty acid. This script tests whether the model is simply predicting the
   training mean, and compares against the 5 constituents that have **real
   measured** A549 values (held out of training): benzyl isothiocyanate 3.06,
   ursolic acid 15.86, apigenin 27.29, stigmasterol 40.53, oleanolic acid
   43.94 µg/mL. **Do not use the per-compound IC₅₀ values in Phase C until
   this is resolved.**

2. **Finish the Monte Carlo** (E1 completed; stopped during E2):
   ```bash
   python -u s08_monte_carlo.py 1200
   ```
   ~1 h single-core. Use `python -u`; without it the output buffers and you
   cannot watch progress. Don't run it alongside the QSAR — they contend for
   every core.

3. **Do the manual jamovi check** — `docs/G3_JAMOVI_CHECK.md`, ~15 min. This is
   the only outstanding piece of Gate G3.

4. **Send `docs/MSU_IIT_DATA_REQUEST.md`** to the partner institution. Phases
   B–D are blocked on the GC-MS/LC-MS peak table, and the dose-range advisory
   is only useful if it reaches them before they run the assay.

Still unwritten: `s03_docking.py`, `s04_admet.py`, `s10_validation.py`.

## Setup

```bash
pip install numpy scipy pandas scikit-learn matplotlib rdkit statsmodels xgboost seaborn requests
```

Tested on Python 3.13, Windows. No compiled dependencies beyond wheels.

## Running

Scripts are ordered by their Phase-A step number and are safe to re-run — network calls are cached to `data/raw/*_cache.json`.

```bash
cd src

python test_stats_mirror.py     # Gate G3(a): stats vs independent implementations
python mtt_model.py             # self-test: IC50 estimators
python virtual_plate.py         # self-test: one simulated experiment

python s01_ligand_prep.py       # PubChem -> ligand_library.csv + ligands.sdf
python s02_target_prep.py       # RCSB   -> verified/selected target panel
python s05_chembl_fetch.py      # ChEMBL -> qsar_dataset_A549.csv
python s06_qsar_train.py        # Gates G2 + G4, Moringa potency predictions
python s08_monte_carlo.py 1200  # power + bias analysis (slow; ~1 h)
python s09_jamovi_export.py     # Gate G3(b) package
python s07_mixture_model.py     # Phase C demo (synthetic composition)
```

## Layout

```
data/raw/          ChEMBL pulls, API caches, extract_composition_TEMPLATE.csv
data/processed/    curated datasets, ligand library, verified target panel
data/structures/   3D ligand SDF
src/               pipeline (see table above)
results/tables/    CSV outputs
results/figures/   PNG figures
results/prediction_registry/   sealed pre-registration (Phase E)
docs/              data request, jamovi check protocol
```

---

## What this study will and will not claim

**Will:** predicted binding of constituents to A549-relevant targets; predicted ADMET; a predicted extract IC₅₀ with an explicit uncertainty interval; a quantitatively justified assay design; a prospectively sealed, falsifiable prediction.

**Will not:** that docking scores prove mechanism; that a predicted IC₅₀ substitutes for measurement; that in vitro cytotoxicity implies clinical efficacy.

---

## Key caveats recorded so far

- **ChEMBL pull is partial.** The connection dropped at 21,000 records; curation ran on what was retrieved. Re-run `s05_chembl_fetch.py` to extend it — the script is resumable by re-fetching.
- **Peak area ≈ mass fraction** is assumed by the mixture model unless MSU-IIT supplies response factors. This is Limitation #1.
- **Thioredoxin reductase 1 has no accepted structure** — every candidate had only FAD bound, giving no inhibitor pocket to anchor a docking box or run the G1 redocking check. Either widen the candidate list or drop the target and say so.
- **Some auto-selected structures still need a human look**: the Bcl-2 pick lists both `BCL2` and `BCL2L1` genes, the mTOR pick is an FRB–rapamycin–FKBP complex rather than the kinase domain, and the tubulin pick needs confirming as colchicine-site. Resolution alone is not a sufficient selection rule.
