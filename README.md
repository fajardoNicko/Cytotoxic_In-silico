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
| Monte Carlo (power, D1/D2/D6 bias) | `src/s08_monte_carlo.py` | ✅ E1 to E5 complete (2026-08-20) |
| Ligand library (47 compounds, 3D) | `src/s01_ligand_prep.py` | ✅ |
| Target panel verification (RCSB) | `src/s02_target_prep.py` | ✅ 15/16 accepted |
| ChEMBL A549 curation (6,972 compounds) | `src/s05_chembl_fetch.py` | ✅ |
| QSAR training + Gates G2/G4 | `src/s06_qsar_train.py` | ⚠️ **G4 pass, G2 FAIL (scaffold)** |
| QSAR reality check | `src/s06b_qsar_reality_check.py` | ✅ run — model is compressed toward the mean |
| Mixture model (Phase C) | `src/s07_mixture_model.py` | ⏸️ retained, **no longer the primary route** |
| Potency calibration / QSAR verdict | `src/s06c_calibrate_potency.py` | ✅ **QSAR ruled unusable for potency** |
| Literature prior + dose-range advisory | `src/s07b_literature_prior.py` | ✅ memo issued |
| Pre-registration + SHA-256 seal | `src/s10_preregister.py` | ✅ **sealed 2026-08-20** |
| **Full methods write-up** | `docs/METHODS.md` / `.docx` | ✅ |

**⚠ Composition data is not coming.** MSU-IIT confirmed that relative peak area
and PubChem CID are **not** included with the MTT deliverable — there will be no
quantitative GC-MS/LC-MS table. Phases B–D are **not** blocked: the potency
prediction was already rebuilt on a literature prior, and Phase B proceeds on the
literature-derived 47-compound library with every result labelled as
literature-based rather than batch-verified. See the amendment in
**[IN_SILICO_PLAN.md](IN_SILICO_PLAN.md)** §3 and
**[docs/MSU_IIT_DATA_REQUEST.md](docs/MSU_IIT_DATA_REQUEST.md)**.

**Still needed from MSU-IIT:** (a) a qualitative compound-name list if any
characterisation exists, (b) extract yield / solvent / stock mg/mL, and
(c) **raw per-well absorbance** after the assay — without (c) the sealed
prediction cannot be scored at all.

---

## ▶ Resume here (next session)

1. **Do the manual jamovi check** — `docs/G3_JAMOVI_CHECK.md`, ~15 min. This is
   the only outstanding piece of Gate G3, and the only task that cannot be
   automated.

2. **Email the SHA-256 hash** from `results/prediction_registry/REGISTRY.md` to
   the adviser. The prediction was sealed 2026-08-20 as
   `7f73b66d7d33f191046f842738c3abb1395fffa482a9313f8d7be87aedfd76c2`. The
   timestamp is what makes it falsifiable, so send it before any assay data
   exists. Do not edit `prediction_v1.json` afterwards.

3. **Send the dose-range memo** — `results/prediction_registry/DOSE_RANGE_MEMO.md`.
   Time-critical: it is only useful if it reaches MSU-IIT *before* they run the
   assay. Recommended series **50/100/200/400/800/1600 ppm**; P(true IC₅₀ >
   200 ppm) is 65–94% depending on prior.

4. **Write `s03_docking.py`** against the literature library — this is now the
   main remaining build, and it is unblocked.

Gate G2 **failed** on the scaffold split (XGBoost test R² = 0.562 vs the
pre-declared ≥ 0.60) while passing on the random split (0.712). Since the
*Moringa* constituents are novel scaffolds, the scaffold number governs.
`s06c_calibrate_potency.py` resolved the consequence: per-compound error is
~11-fold at 1σ inside the *Moringa* physicochemical envelope, and the bias
direction is unresolved (envelope says +0.213 log, the five real measured
constituents say −0.328 log — opposite signs). **The QSAR is therefore used for
ranking/triage only, never for absolute potency**, and Phase C runs on the
literature prior. Gate G4 passed: doxorubicin predicted at 0.634 µM against a
published 0.1–10 µM window, held out of training.

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
