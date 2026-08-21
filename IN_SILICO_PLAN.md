# In Silico Replication Plan
## *Moringa oleifera* Ethanolic Bark Extract vs. A549 Lung Cancer Cells

**Companion to:** "In Vitro Cytotoxic Activity of Malunggay's (*Moringa oleifera*) Ethanolic Bark Extract against A549 Lung Cancer Cells" (Dalistan, Roque, Soriano, Salvilla)
**Purpose:** A computational study that (a) predicts the wet-lab outcome *before* it is run, (b) mirrors the exact assay, formulas and statistics of the research plan, and (c) is falsifiable against the real MTT data when it arrives.
**Status:** Plan — no code written yet.

---

## 0. Ground truth extracted from the paper

Everything below is what the in silico model must reproduce *exactly*. Deviating from any of these breaks the "replication" claim.

| Element | Value in the paper |
|---|---|
| Cell line | A549 human lung adenocarcinoma |
| Test article | Ethanolic bark extract, *M. oleifera* |
| Extraction | Soxhlet, 70% EtOH, 1:10 w/v (25 g powder : 250 mL), 3 cycles, 4–6 h, ×4 runs → 1 L; rotary-evaporated to crude extract by MSU-IIT |
| Drying | 70 °C oven, 1 h 30 min; ground to coarse powder |
| Doses | **12.5, 25, 50, 100, 200 ppm** (= µg/mL), Setups A–E |
| Negative control | "McCoy's medium", untreated |
| Positive control | **Doxorubicin** |
| Assay | MTT proliferation assay, OD read at **570 nm** |
| Viability | `%V = [(Abs_sample − Abs_blank) / (Abs_control − Abs_blank)] × 100` |
| Inhibition | `%Cytotoxicity = 100 − %Viability` |
| IC₅₀ | Linear regression of **log₁₀(concentration) vs %viability**, solve `y = mx + b` at y = 50 |
| Design | True experimental, **posttest-only control group** |
| Statistics | Shapiro–Wilk → one-way ANOVA (or Kruskal–Wallis if non-normal) → **Tukey HSD** post hoc, α = 0.05, in **jamovi** |
| Reporting | mean ± SD |

**A549 genotype the model must respect** (this is what makes the in silico work biologically specific rather than generic):
- `KRAS` **G12S** homozygous → constitutive RAS/MAPK signalling
- `STK11`/`LKB1` **null** → loss of AMPK-mediated growth restraint
- `TP53` **wild type** → apoptosis via p53 axis is intact and druggable
- `CDKN2A` deleted → unrestrained CDK4/6–Rb
- `KEAP1` **G333C** → constitutive **NRF2** activation → intrinsic chemoresistance and high antioxidant capacity

That last one matters a lot: *Moringa* isothiocyanates are textbook NRF2 **activators**. Modelling this gives a real, defensible mechanistic hypothesis — the extract may hit an already-NRF2-saturated cell, which would predict *weaker*-than-expected potency. This is a genuine scientific contribution, not decoration.

---

## 1. Defects found in the research plan (fix these; the in silico study addresses several directly)

These are ordered by how badly they threaten the wet-lab result.

**D1 — No vehicle control. (Critical.)**
A crude rotary-evaporated extract must be reconstituted in DMSO (or EtOH) before dosing. Without a vehicle-only control at the *same* final solvent %, any cytotoxicity observed is confounded by solvent. A549 tolerates ≤0.5% v/v DMSO; above ~1% you get solvent kill. **Add a vehicle control group.** The in silico plate simulation will include and quantify this.

**D2 — 200 ppm ceiling is likely too low. (Critical.)**
Crude plant extracts on A549 commonly show IC₅₀ in the 100–1000 µg/mL range; published *M. oleifera* leaf/bark extract IC₅₀ values against A549 frequently exceed 200 µg/mL. If the true IC₅₀ > 200 ppm, the dose–response never crosses 50%, and the `y = mx + b` solution becomes an **extrapolation outside the tested range — statistically invalid**, and Specific Question 3 becomes unanswerable. Phase C of this plan exists specifically to predict this before you spend the samples, and to recommend an extended series (e.g. add 400, 800, 1600 ppm) if needed.

**D3 — "McCoy's medium" is almost certainly a copy-paste artifact.**
McCoy's 5A is the standard medium for **HT-29 colorectal** cells. ATCC's recommended medium for A549 (CCL-185) is **F-12K**; RPMI-1640 or DMEM + 10% FBS are also standard. Either switch to F-12K/RPMI, or have MSU-IIT confirm and document validated A549 growth in McCoy's 5A with a growth curve.

**D4 — Template artifacts from the source (watermelon/MCF-7) paper.**
- p.4: "compared with **Paclitaxel 3**" contradicts doxorubicin everywhere else.
- p.4: "varying concentrations of the **peel solution**" — should be bark extract.
- p.3: "stop or slow the growth of **breast cancer** cells."
- p.6: "safer, more affordable **breast cancer** therapies."
- p.9: "**Molinga** oleifera" typo.
- p.7 Table 1 header: "cValue" — should be "IC₅₀ Value"; "blank (McCoy **2**)".
- References list Abd-Rabou et al. as **2017**, body cites **2023**.
- Two identical Bhadresha et al. (2022a/2022b) entries.
Fix all before submission — reviewers read these as evidence the methods weren't authored deliberately.

**D5 — Log-linear IC₅₀ regression is the weakest defensible method.**
Dose–response is sigmoidal, not linear. Fitting a straight line to log(C) vs %V is only acceptable over the near-linear 20–80% region. **Report the paper's log-linear IC₅₀ as primary (for protocol fidelity) and a 4-parameter logistic (4PL) nonlinear fit as a sensitivity analysis.** If the two disagree by >20%, that itself is a finding.

**D6 — Polyphenol/MTT chemical interference not controlled.**
Phenolic-rich extracts can reduce MTT tetrazolium directly, with no cells present, inflating apparent viability at high doses. **Add a cell-free control: extract + medium + MTT at every concentration.** Modelled explicitly in Phase D.

**D7 — Unspecified assay parameters.**
Seeding density, incubation time (24/48/72 h), exposure time, replicate structure (n wells, n independent runs), MTT concentration and solubilisation agent are all absent. The in silico design will lock these to defensible values and hand you a written protocol.

**D8 — No selectivity control.**
Cytotoxicity to A549 alone says nothing about therapeutic value. Add a normal lung line (MRC-5, WI-38 or BEAS-2B) and report **Selectivity Index = IC₅₀(normal) / IC₅₀(A549)**. If the budget can't take it, Phase B predicts SI computationally and you report it as a limitation with a modelled estimate.

**D9 — n = 3 may be underpowered.**
Never stated in the plan. Phase D's Monte Carlo answers "what n do I need to detect the effect I expect?" quantitatively.

---

## 2. Objectives and success criteria

**Primary objective.** Produce a computational model that predicts, *before the wet lab runs*, the %viability at 12.5/25/50/100/200 ppm and the IC₅₀ of the extract against A549 — with stated uncertainty — and then test that prediction against the real data.

**Secondary objectives.**
1. Identify which *M. oleifera* bark constituents drive the cytotoxicity and via which A549-relevant targets.
2. Predict whether the dose range is adequate (D2) and what n is required (D9).
3. Predict selectivity vs. normal lung cells (D8).
4. Rank the extract against doxorubicin on the same in silico scale.

**Success criteria (pre-declared, so the study can fail honestly).**
| Gate | Criterion |
|---|---|
| G1 Docking validity | Redocking of each co-crystallised ligand reproduces the crystal pose at **RMSD ≤ 2.0 Å** |
| G2 QSAR validity | 5-fold CV **Q² ≥ 0.5**, external test **R² ≥ 0.6**, y-randomisation Q²_rand < 0.2 |
| G3 Stats fidelity | Python and jamovi agree to **3 decimal places** on the same CSV |
| G4 Positive-control anchor | Predicted doxorubicin A549 IC₅₀ within **3-fold** of published values |
| G5 Final validation | Predicted extract IC₅₀ within **3-fold** of the observed wet-lab IC₅₀; %viability RMSE ≤ 15 percentage points across the 5 doses; Spearman ρ ≥ 0.9 for dose rank order |

G4 is the honesty check that runs *now*: if the pipeline can't retrodict a drug whose answer is already known, its extract prediction is worthless.

---

## 3. Data intake specification — what to demand from MSU-IIT

> ### ⚠ AMENDMENT (2026-08-20) — the composition table is not coming
>
> MSU-IIT confirmed that **`pct_area` (relative peak area) and `pubchem_cid`
> are not included** with the MTT deliverable. There will be no quantitative
> GC-MS/LC-MS composition table for the tested batch.
>
> **Consequences, bounded explicitly:**
>
> | Lost | Unaffected |
> |---|---|
> | Composition vector `pᵢ` → any bottom-up mixture calculation (CA or Bliss). *Already retired on independent grounds in §6, so the practical loss is nil.* | Literature-prior potency estimate and the dose-range advisory — never depended on composition |
> | Batch-verified constituent identity. Docking/ADMET/network results now describe **literature-reported** constituents, not **batch-verified** ones | The entire virtual MTT assay, Monte Carlo power analysis and control-set evaluation |
> | Calibration of the MTT-interference coefficient from measured total phenolic content → E4 sweeps a plausible range instead | Statistics mirror, pre-registration, and Phase F validation |
>
> **Phases B–D are therefore NOT blocked.** Phase B proceeds on the
> literature-derived 47-compound library (§4 A2), with every result labelled
> as literature-based. This becomes Limitation #1 in the manuscript.
>
> **Revised minimum ask** (see `docs/MSU_IIT_DATA_REQUEST.md` amendment):
> (a) a *qualitative* compound-name list if any characterisation exists —
> names and match scores only, we can resolve identifiers ourselves;
> (b) extract yield, reconstitution solvent, stock mg/mL, final solvent % v/v;
> (c) **raw per-well absorbance** after the assay — now the single most
> important item, since without it the sealed prediction cannot be scored.

Send this list to the partner institution **now**, so the chemistry arrives in a usable form. Item 1 is withdrawn as blocking per the amendment above; items 2, 5, 6 and 8 remain live.

**Mandatory**
1. **GC-MS and/or LC-MS/MS peak table** — for each peak: compound name, retention time, % **relative peak area**, molecular formula, and **CAS number or PubChem CID**. Ask for the library-match score (NIST similarity) too; drop hits below 80%.
2. **Extract yield** (% w/w of dry bark) and the exact reconstitution solvent + stock concentration (mg/mL).
3. **Total phenolic content** (mg GAE/g extract) and **total flavonoid content** (mg QE/g) — these drive the MTT-interference term (D6).

**Strongly requested**
4. Phytochemical screening (alkaloids, saponins, tannins, glycosides, terpenoids, steroids) — qualitative +/− is fine.
5. DPPH or ABTS radical-scavenging IC₅₀ — feeds the NRF2/redox arm.
6. Moisture content and solvent-residue check on the crude extract.
7. BPI plant-identification certificate reference number.

**Deliver as:** one CSV, `data/raw/extract_composition.csv`, columns:
`compound_name, pubchem_cid, cas, formula, mw, rt_min, pct_area, match_score, method`

If the peak table is a PDF or an image, transcribe it into this CSV manually and have a second person verify — a transcription error here silently corrupts every downstream number.

---

## 4. Phase A — Extract-agnostic build (start immediately, no extract data needed)

Everything here can be done today and is the majority of the work.

**A1. Environment and repository.**
```
Lung_Cancer_In_Silico/
├── data/
│   ├── raw/            # extract_composition.csv, ChEMBL/NCI-60 downloads
│   ├── structures/     # cleaned PDBs, ligand SDF/PDBQT
│   └── processed/
├── src/                       # (as built -- module names avoid a leading
│   ├── stats_mirror.py        #  digit so they remain importable)
│   ├── test_stats_mirror.py
│   ├── mtt_model.py
│   ├── virtual_plate.py
│   ├── s01_ligand_prep.py
│   ├── s02_target_prep.py
│   ├── s03_docking.py         # pending
│   ├── s04_admet.py           # pending
│   ├── s05_chembl_fetch.py
│   ├── s06_qsar_train.py
│   ├── s07_mixture_model.py
│   ├── s08_monte_carlo.py
│   ├── s09_jamovi_export.py
│   └── s10_validation.py      # pending (Phase F)
├── results/
│   ├── figures/
│   ├── tables/
│   └── prediction_registry/    # SHA-256-sealed pre-registration
├── notebooks/
├── IN_SILICO_PLAN.md
└── environment.yml
```
Stack: Python 3.11, RDKit, Open Babel, AutoDock Vina (or smina), PDBFixer/OpenMM, scikit-learn + XGBoost, Mordred, SciPy, statsmodels, pingouin, pandas, matplotlib/seaborn. Version-pin everything in `environment.yml`; record the exact commit for every result.
**No-install fallback** (if compute is a constraint): CB-Dock2 / SwissDock for docking, SwissADME + ADMETlab 3.0 + ProTox-3.0 for ADMET, SwissTargetPrediction + STRING + ShinyGO for network pharmacology, WebGro for short MD. The plan works either way; note which route each result came from.

**A2. Candidate-compound library (built from literature, refined later by the real MS data).**
Assemble every *M. oleifera* constituent with a reported structure, prioritising bark-associated ones. Starting set: 4-(α-L-rhamnopyranosyloxy)benzyl isothiocyanate (moringin), glucomoringin, niazimicin, niaziminin, benzyl isothiocyanate, benzyl glucosinolate, moringine/moringinine, quercetin, kaempferol, rutin, chlorogenic acid, gallic acid, ellagic acid, catechin, epicatechin, vanillin, β-sitosterol, stigmasterol, campesterol, lupeol, oleanolic acid, ursolic acid, hexadecanoic acid, oleic acid, 4-hydroxymellein, octacosanoic acid, apigenin, luteolin, ferulic acid, caffeic acid, syringic acid, p-coumaric acid. Bark is tannin/alkaloid/saponin-rich, so include condensed tannin monomers.
Retrieve canonical SMILES from PubChem, generate 3D conformers (RDKit ETKDGv3 + MMFF94 minimisation), protonate at pH 7.4, export SDF + PDBQT.

**A3. Target panel (A549-specific, not generic).**
For each target, pull the highest-resolution **holo** structure, confirm the co-crystallised ligand, and validate by redocking (Gate G1). Candidate PDB IDs below are starting points — **verify each in the PDB before use**, do not trust them blind.

| Rationale | Gene / protein | Candidate structure |
|---|---|---|
| A549 driver | KRAS G12S / G12D switch-II pocket | verify current G12S/G12D holo entry |
| Proliferation | EGFR kinase domain | 1M17 (erlotinib) — verify |
| Survival | PI3Kα (PIK3CA) | 4JPS — verify |
| Survival | AKT1 | 4EJN — verify |
| Survival | mTOR kinase | verify |
| MAPK | ERK2 / MAPK1 | verify |
| Apoptosis | Bcl-2 | 4LVT — verify |
| Apoptosis | Bcl-xL | verify |
| Apoptosis | Caspase-3 | verify |
| Inflammation | COX-2 (PTGS2) | verify |
| Angiogenesis | VEGFR2 (KDR) | verify |
| Mitosis | β-tubulin colchicine site | 1SA0 — verify |
| **Doxorubicin anchor** | Topoisomerase IIα DNA cleavage complex | verify |
| **A549 chemoresistance** | KEAP1 Kelch domain (NRF2 pocket) | verify |
| Redox | Thioredoxin reductase (TXNRD1) | verify |
| Metabolic | LKB1/AMPK axis (context only) | — |

**A4. Positive-control anchoring (Gate G4).**
Run doxorubicin — and, because the paper once says paclitaxel, run paclitaxel too — through the entire pipeline (docking, ADMET, QSAR). Compare the QSAR-predicted A549 IC₅₀ against published values. **If G4 fails, stop and fix the pipeline before touching the extract.** This is the single most important checkpoint in the study.

**A5. QSAR training set.**
Pull all A549 cytotoxicity records from **ChEMBL** (cell-line target = A549; endpoints IC₅₀/GI₅₀/EC₅₀) and the **NCI-60 / NCI DTP** A549 panel. Curate: single-compound records only, standardise units to µM → pIC₅₀, drop qualified values (">", "<"), deduplicate by InChIKey (median of replicates), remove records with no structure. Expect a few thousand usable rows. Build a parallel **normal-lung** set (MRC-5 / WI-38 / BEAS-2B / IMR-90) for the Selectivity Index model (D8). Hold out a random 20% as an external test set and never touch it until the end.

**A6. QSAR model.**
Features: RDKit descriptors + Morgan fingerprints (r = 2, 2048 bits) + selected Mordred descriptors. Models: Random Forest and XGBoost, plus a Ridge baseline. Tune by nested 5-fold CV. Report Q², external R², RMSE, y-randomisation. Define the **applicability domain** by Tanimoto similarity to the training set (flag any *Moringa* compound with max similarity < 0.3 as out-of-domain — its prediction is a guess and must be labelled as such). Gate G2.

**A7. Statistics mirror (Gate G3).**
Implement the paper's exact analysis chain in Python: Shapiro–Wilk, Levene, one-way ANOVA (with the MST/MSE decomposition written out as in pp. 12–13), Kruskal–Wallis fallback, Tukey HSD (and Games-Howell if variances are unequal; Dunn if KW). Then take one simulated CSV, run it in **jamovi** and in Python, and require agreement to 3 dp. This proves your simulation speaks the same statistical language as the real analysis.

**Phase A deliverables:** validated docking setup, validated QSAR model, doxorubicin anchor result, verified stats mirror. All independent of the extract.

---

## 5. Phase B — Compound-level analysis (starts when the MS data arrives)

**B1. Reconcile.** Match the MS peak table to the A2 library by PubChem CID. Compounds present in MS but absent from A2 → add them. Compounds in A2 but absent from MS → drop (or keep flagged as "literature-reported, not detected in this batch").

**B2. Molecular docking.** Dock every detected compound + doxorubicin against every target in A3. AutoDock Vina, exhaustiveness 16–32, 9 poses, grid box centred on the co-crystallised ligand with ≥5 Å padding. Triplicate runs with different seeds; report mean ± SD binding affinity (kcal/mol). Scale: ~40 compounds × ~14 targets × 3 seeds ≈ 1700 runs ≈ 15–30 h on a laptop — run it overnight.
Output: compound × target affinity heatmap; interaction fingerprints (H-bonds, hydrophobic, π-stacking) for the top 10 complexes via PLIP.

**B3. ADMET and toxicity.** SwissADME / ADMETlab 3.0 / pkCSM for Lipinski, Veber, TPSA, logP, GI absorption, BBB, CYP inhibition, P-gp substrate; ProTox-3.0 for predicted oral LD₅₀, toxicity class, hepatotoxicity, carcinogenicity, cytotoxicity. This gives your discussion a *safety* dimension that the wet lab cannot supply — a real differentiator.

**B4. Network pharmacology.** SwissTargetPrediction / STITCH → predicted targets per compound. Intersect with lung-adenocarcinoma gene sets (GeneCards, OMIM, DisGeNET, TCGA-LUAD DEGs). Build the compound–target–pathway network in Cytoscape; STRING PPI on the intersection; hub genes by degree/betweenness (cytoHubba); KEGG + GO enrichment (DAVID or ShinyGO, BH-corrected q < 0.05). Expect NSCLC, PI3K-AKT, MAPK, apoptosis, and (given KEAP1) the NRF2/oxidative-stress pathway to surface.

**B5. Molecular dynamics (optional but high-value).** Take the single best compound–target complex plus the doxorubicin reference. GROMACS, CHARMM36 or AMBER ff14SB, TIP3P water, 0.15 M NaCl, 100 ns production. Analyse RMSD, RMSF, Rg, SASA, H-bond occupancy, PCA; compute **MM-PBSA/MM-GBSA** binding free energy with per-residue decomposition. If GROMACS is unavailable, WebGro gives a 20–50 ns run for free. Skip this rather than do it badly — a 5 ns run proves nothing.

**B6. Per-compound potency prediction.** Run the QSAR model (A6) over all detected compounds → predicted pIC₅₀ → IC₅₀ in µM → convert to **µg/mL** via `IC50_µg/mL = IC50_µM × MW / 1000`. Attach the applicability-domain flag to every prediction. Also predict normal-lung IC₅₀ for the Selectivity Index.

---

## 6. Phase C — Extract-level potency prediction (the bridge to ppm)

> ### ⚠ AMENDMENT (2026-07-30) — the bottom-up route below was tested and rejected
>
> Phase C as originally written sums per-compound QSAR potencies into an extract
> IC₅₀. That was implemented, then **abandoned on evidence**:
>
> - Gate G2 **failed** on a scaffold split (XGBoost test R² = 0.562 vs the
>   pre-declared ≥ 0.60), while passing on a random split (0.712). The *Moringa*
>   constituents are novel scaffolds, so the scaffold figure governs.
> - Per-compound error is **~11-fold at 1σ** (scaffold-test RMSE = 1.041 log units
>   inside the *Moringa* physicochemical envelope, n = 754).
> - The **bias direction is unresolved**: the envelope estimate says the model is
>   under-potent (+0.213 log), the five real *Moringa* measurements say over-potent
>   (−0.328 log). Opposite signs.
>
> Summing components with 11-fold errors and an unknown bias sign would have
> produced a confident-looking number carrying no information. **Phase C is
> therefore driven by a literature prior on crude-extract potency**
> (`s07b_literature_prior.py`), which is both better anchored and directly
> comparable to what the MTT assay measures.
>
> The QSAR is retained for **ranking/triage only** — choosing which constituents
> are worth docking — never for absolute potency. The concentration-addition and
> Bliss machinery in `s07_mixture_model.py` is kept, and becomes usable if
> measured IC₅₀ values are obtained for the major constituents.
>
> Gate G4 (doxorubicin anchor) **passed**: 0.634 µM predicted against a published
> 0.1–10 µM window, with doxorubicin held out of training.

This is the step that turns per-molecule chemistry into the number the wet lab will actually measure. Do not skip the assumption bookkeeping here — it is where this kind of study usually falls apart.

**C1. Composition vector.** From the MS table, convert % relative peak area to mass fraction `pᵢ`. State the assumption explicitly: *relative peak area ≈ relative mass*, which is only approximate (it ignores per-compound detector response factors). If MSU-IIT can supply response factors or an internal standard, use them. Otherwise this is Limitation #1 and must be written into the paper.

**C2. Concentration addition (Loewe) — the primary estimate.**
For total extract concentration `C` (µg/mL), component *i* is present at `pᵢ·C`. The mixture reaches 50% effect when `Σ(pᵢ·C / IC50ᵢ) = 1`, so:

```
IC50_extract = 1 / Σ (pᵢ / IC50ᵢ)          [pᵢ = mass fraction, IC50ᵢ in µg/mL]
```

Compute this **twice**, to bound the answer:
- **Conservative bound:** `pᵢ` normalised over the *total* extract mass; all unidentified mass treated as inert → higher (weaker) IC₅₀.
- **Optimistic bound:** `pᵢ` renormalised over the *identified* fraction only → lower (stronger) IC₅₀.
Report the interval. If MS identifies only 55% of the mass, this interval will be wide — say so plainly rather than reporting a false point estimate.

**C3. Bliss independence — the secondary estimate.**
Assume independent action: `V_mix(C) = Π_i V_i(pᵢ·C)`, each `V_i` a Hill function from the QSAR IC₅₀ and an assumed Hill slope. Solve numerically for `V_mix = 50%`. CA and Bliss bracket most real mixture behaviour; if they agree within 2-fold, confidence is decent.

**C4. Literature reality check.** Independently collect published IC₅₀ values for *M. oleifera* extracts against A549 (and other NSCLC lines) — bark if available, leaf/root/seed otherwise. If your model's prediction sits far outside the published envelope, the model is wrong, not the literature. Adjust or explain.

**C5. Dose-range adequacy call (resolves D2).**
Compare predicted IC₅₀ to the 12.5–200 ppm window.
- Predicted IC₅₀ within 25–150 ppm → range is good, proceed.
- Predicted IC₅₀ > 200 ppm → **issue a formal recommendation to extend the series** (add 400, 800, 1600 ppm; or replace 12.5 ppm with a higher anchor). Deliver this *before* MSU-IIT runs the assay. This alone justifies the entire in silico study.
- Predicted IC₅₀ < 12.5 ppm → recommend adding lower doses (1.5625, 3.125, 6.25 ppm).

**C6. NCI activity classification.** Apply the standard crude-extract criterion (NCI: active if IC₅₀ ≤ 20–30 µg/mL) to the predicted value and state the expected verdict up front.

---

## 7. Phase D — Virtual MTT assay (the actual "replication")

This is where you simulate the experiment that will physically happen, well by well.

**D-1. Plate layout.** Standard 96-well, simulated explicitly:

| Group | Wells | Content |
|---|---|---|
| Blank | 3–6 | Medium + MTT, **no cells** |
| Negative control | 6 | Cells + complete medium (F-12K or RPMI + 10% FBS) |
| **Vehicle control** | 6 | Cells + medium + DMSO at the highest final % used (**fixes D1**) |
| Positive control | 6 | Cells + doxorubicin (dose series, or single reference dose) |
| Extract A–E | 6 each | 12.5 / 25 / 50 / 100 / 200 ppm |
| **Interference control** | 5 | Extract at each dose + medium + MTT, **no cells** (**fixes D6**) |

Locked assay parameters (fixes D7): seeding 5 × 10³ – 1 × 10⁴ cells/well, 24 h attachment, 48 h treatment, MTT 0.5 mg/mL for 4 h, formazan solubilised in DMSO, read at **570 nm** with 630 nm reference. **n = 3 technical replicates × 3 independent biological runs** as the starting design — then let D-4 tell you if that's enough.

**D-2. Generative model.** True viability from a 4-parameter logistic:
```
V(C) = V_min + (V_max − V_min) / (1 + (C / IC50)^H)
```
with `V_max = 100`, `V_min` = predicted plateau (5–20% typical for crude extracts), `IC50` from Phase C, Hill slope `H` sampled 0.8–1.8.

Simulated absorbance per well:
```
A_well = A_blank + (A_ctrl − A_blank) · V(C)/100 · (1 + ε_well) + δ_edge + δ_plate + k_interf·C
```
- `A_blank ≈ 0.05`, `A_ctrl ≈ 0.8–1.2` (A549, 48 h, proper seeding)
- `ε_well ~ N(0, 0.06)` — pipetting + seeding variance (well CV 5–8%)
- `δ_plate ~ N(0, 0.08)` — plate-to-plate offset across the 3 biological runs
- `δ_edge` — positive bias on perimeter wells from evaporation (or design the layout to leave the perimeter unused, which is the better fix and should be the recommendation)
- `k_interf·C` — direct MTT reduction by polyphenols, calibrated from the total phenolic content (item 3 of §3); recovered and subtracted by the interference control wells

**D-3. Run the paper's pipeline on the simulated data.** Apply `%V = [(Abs_sample − Abs_blank)/(Abs_control − Abs_blank)] × 100`, then `%Cytotoxicity = 100 − %V`, then the log-linear IC₅₀ *and* the 4PL IC₅₀ (D5). Produce the exact tables and figures the paper will need: mean ± SD per group, dose–response curve, log-linear regression plot with the `y = mx + b` equation and R² displayed.

**D-4. Monte Carlo — 10,000 virtual experiments.** Re-run D-2/D-3 ten thousand times with fresh noise. Report:
- Distribution of the **estimated** IC₅₀ vs the **true** input IC₅₀ → bias and 95% coverage interval
- Probability that ANOVA returns p < 0.05
- Probability that Tukey HSD flags each individual dose vs. negative control
- **Statistical power at n = 3, 4, 6, 8** → a concrete recommended n (**fixes D9**)
- Probability that Shapiro–Wilk rejects normality, i.e. how often the real experiment will end up on the Kruskal–Wallis branch
- How often the log-linear IC₅₀ is an out-of-range extrapolation (**quantifies D2**)

**D-5. Sensitivity analysis.** Sweep true IC₅₀ across a 10-fold band, well CV from 3% to 15%, Hill slope 0.8–1.8, and unidentified-mass-fraction 0–60%. Tornado plot of which assumption most affects the conclusion. This is what makes the study bulletproof — you will know exactly which input, if wrong, breaks the result.

---

## 8. Phase E — Pre-registration (do this before any wet-lab result exists)

1. Freeze `results/prediction_registry/prediction_v1.json` containing: predicted IC₅₀ (point + interval, both CA and Bliss), predicted %viability ± SD at each of the 5 doses, predicted dose rank order, predicted ANOVA verdict, predicted Tukey pattern, predicted Selectivity Index, and the recommended n.
2. Compute its **SHA-256 hash**, record it in `results/prediction_registry/REGISTRY.md` with a timestamp, and email the hash to your adviser (and optionally post it to OSF).
3. Do not modify the file afterwards. Any revision becomes `prediction_v2.json` with its own hash and a written justification.

A sealed prediction is the difference between a modelling exercise and a genuine test. Without it, "our in silico results agreed with our in vitro results" is unfalsifiable.

---

## 9. Phase F — Validation against the real MTT data

When MSU-IIT returns the absorbance readings:

1. Run the real raw absorbances through the *same* `src/09_stats_mirror.py` used for the simulation.
2. Compare predicted vs observed:
   - **Fold error** on IC₅₀ (`FE = predicted / observed`); success = within 3-fold (G5)
   - **RMSE** and **MAE** on %viability across the 5 doses
   - **Lin's concordance correlation coefficient**
   - **Spearman ρ** on dose rank order
   - **Bland–Altman** plot of predicted vs observed %viability
3. Overlay predicted and observed dose–response curves in one figure. This will be the paper's headline figure.
4. If the prediction fails, diagnose honestly against the D-5 sensitivity analysis: was it the composition vector, the QSAR applicability domain, the additivity assumption, MTT interference, or NRF2-mediated resistance? **A well-diagnosed failed prediction is a stronger paper than an unexamined success.**
5. Report the docking/MD mechanism as a *hypothesis* consistent (or not) with the observed potency — never as proof of mechanism. In silico binding affinity is not evidence of cellular activity.

---

## 10. Timeline

| Week | Work | Blocked on extract data? |
|---|---|---|
| 1 | Repo, environment, target panel, structure prep, redocking validation (G1) | No |
| 2 | Ligand library, ChEMBL/NCI-60 curation, QSAR training (G2) | No |
| 3 | Doxorubicin anchor (G4), stats mirror + jamovi cross-check (G3) | No |
| 4 | Virtual MTT engine + Monte Carlo on *literature-prior* IC₅₀; **deliver preliminary dose-range recommendation to MSU-IIT** | No |
| 5 | ← MS data arrives → reconcile composition, full docking sweep | **Yes** |
| 6 | ADMET, network pharmacology, per-compound QSAR | Yes |
| 7 | Mixture model, final IC₅₀ prediction, **seal pre-registration** | Yes |
| 8 | MD (optional), figures, manuscript sections | Yes |
| 9+ | ← MTT results arrive → validation, Bland–Altman, discussion | Yes |

Weeks 1–4 are ~60% of the total effort and are unblocked right now. Week 4's dose-range recommendation must reach MSU-IIT **before** they run the assay to be useful.

---

## 11. Risk register

| # | Risk | Impact | Mitigation |
|---|---|---|---|
| R1 | MS identifies <50% of extract mass | Wide prediction interval | Report CA bounds honestly; request higher-resolution LC-MS/MS |
| R2 | Detected compounds fall outside QSAR applicability domain | Predictions unreliable | Flag every out-of-domain prediction; fall back to docking-rank + literature IC₅₀ |
| R3 | True IC₅₀ > 200 ppm | Wet-lab IC₅₀ unobtainable | Phase C5 recommendation delivered pre-assay |
| R4 | Additivity assumption wrong (real synergy/antagonism) | Prediction off | CA + Bliss bracket; report both; discuss as limitation |
| R5 | Peak area ≠ mass fraction | Systematic composition error | Request response factors / internal standard; sensitivity-test in D-5 |
| R6 | No holo structure for a target | Can't dock reliably | Substitute the closest homolog, or drop the target and say so |
| R7 | MS data never arrives / arrives late | Phases B–D stall | Phase A + Monte Carlo on literature priors still yields a complete, publishable methods-and-power study |
| R8 | Compute too limited for docking sweep | Timeline slip | Web-server fallbacks (CB-Dock2, SwissDock); reduce target panel to the top 8 |
| R9 | jamovi and Python disagree | Stats claim undermined | G3 checkpoint catches it in Week 3, not at submission |

---

## 12. Deliverables

1. **Reproducible repository** with pinned environment and a one-command rerun.
2. **Target-panel table** with PDB IDs, resolution, co-crystal ligand, and redocking RMSD (G1 evidence).
3. **Compound × target docking heatmap** + top-10 interaction diagrams.
4. **ADMET/toxicity table** for all detected compounds.
5. **Network pharmacology figures**: compound–target–pathway network, STRING PPI, KEGG/GO enrichment bar plots.
6. **QSAR model card**: training set size, Q², external R², y-randomisation, applicability-domain coverage.
7. **Predicted dose–response curve** at 12.5–200 ppm with uncertainty band, and predicted IC₅₀ with CA/Bliss interval.
8. **Monte Carlo power report** with the recommended n and the probability of a successful IC₅₀ determination.
9. **Sealed pre-registration** file + SHA-256 hash + timestamp.
10. **Dose-range recommendation memo** to MSU-IIT (delivered by Week 4).
11. **Improved wet-lab protocol** incorporating fixes D1, D3, D6, D7, D8.
12. **Validation report** comparing in silico vs in vitro, with Bland–Altman and an honest diagnosis.

---

## 13. What this study can and cannot claim

**Can claim:** predicted binding of specific constituents to A549-relevant targets; predicted ADMET/toxicity profile; predicted extract potency with stated uncertainty; a quantitatively justified assay design; a prospectively sealed, testable prediction.

**Cannot claim:** that docking scores prove mechanism; that predicted IC₅₀ substitutes for measurement; that in vitro cytotoxicity implies clinical efficacy. Write these limitations into the paper explicitly — every competent reviewer will look for them, and stating them first is a strength.
