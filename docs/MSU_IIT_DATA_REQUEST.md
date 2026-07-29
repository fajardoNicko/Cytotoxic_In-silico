# Data request — MSU-IIT partner institution

**Study:** In vitro cytotoxic activity of *Moringa oleifera* ethanolic bark extract against A549 lung cancer cells
**Requested by:** Dalistan, Roque, Soriano, Salvilla
**Purpose:** A parallel in silico model predicts the assay outcome *before* the MTT run. Items 1–3 block the entire prediction; the rest improve its accuracy.

---

## 1. Chemical composition of the crude extract — **REQUIRED, blocking**

A **GC-MS and/or LC-MS/MS peak table** of the rotary-evaporated crude extract. For every peak:

| Field | Meaning |
|---|---|
| `compound_name` | Library-assigned name |
| `pubchem_cid` **or** `cas` | Identifier — *at least one is essential*; a name alone is often ambiguous |
| `formula`, `mw` | Molecular formula and weight |
| `rt_min` | Retention time (min) |
| `pct_area` | **Relative peak area, % of total** |
| `match_score` | Library (e.g. NIST) similarity score — we discard hits < 80 |
| `method` | `GC-MS` or `LC-MS/MS` |

Use the template at `data/raw/extract_composition_TEMPLATE.csv`. **CSV or Excel, not PDF** — a PDF has to be transcribed by hand, which risks silent errors.

> **Please also state whether relative peak area may be treated as relative mass.** If per-compound response factors or an internal standard were used, send those too. Without them we must record "peak area ≈ mass" as a stated limitation.

## 2. Extract handling — **REQUIRED, blocking**

- Extract **yield** (% w/w of dry bark)
- **Reconstitution solvent** and **stock concentration** (mg/mL)
- **Final solvent % v/v** in the highest-dose well (needed for the vehicle control)
- Storage conditions and time from evaporation to assay

## 3. Bulk phytochemistry — **REQUIRED, blocking**

- **Total phenolic content** (mg gallic acid equivalent / g extract)
- **Total flavonoid content** (mg quercetin equivalent / g extract)

These calibrate the predicted MTT interference (see §6).

## 4. Strongly requested

- Qualitative phytochemical screening (alkaloids, saponins, tannins, glycosides, terpenoids, steroids)
- DPPH and/or ABTS radical-scavenging IC₅₀
- Moisture content; residual-solvent check on the crude extract
- BPI plant-identification certificate reference number

---

## 5. Assay parameters we need confirmed *before* the run

The research plan does not specify these, and the in silico power analysis needs them fixed:

| Parameter | Please confirm |
|---|---|
| Culture medium | The plan says **McCoy's 5A** — that is the standard medium for **HT-29 colorectal** cells. ATCC specifies **F-12K** for A549 (CCL-185); RPMI-1640 or DMEM + 10% FBS are also standard. Please confirm which will be used, and if McCoy's, whether A549 growth in it has been validated. |
| Seeding density | cells/well (we assume 5×10³–1×10⁴) |
| Attachment time | h before treatment (we assume 24) |
| Treatment exposure | h (we assume 48) |
| MTT concentration / incubation | mg/mL and h (we assume 0.5 mg/mL, 4 h) |
| Formazan solubiliser | DMSO / acidified isopropanol / SDS-HCl |
| Read wavelength | 570 nm confirmed; is a 630 nm reference subtracted? |
| Replicates | technical wells per group × independent runs |
| Doxorubicin | dose or dose series, and vendor/lot |

## 6. Three control groups we ask you to add

Each addresses a specific way the experiment can produce an uninterpretable result.

1. **Vehicle control** — cells + medium + solvent at the *same* final % v/v as the treated wells. Without it, solvent cytotoxicity is indistinguishable from extract cytotoxicity. Our simulation shows that normalising to a medium-only control instead biases IC₅₀ by **−6% to −39%** as solvent toxicity rises from 4% to 20%.
2. **Cell-free interference control** — extract + medium + MTT at every concentration, **no cells**. Polyphenol-rich extracts reduce tetrazolium directly, inflating apparent viability.
3. **Normal-cell comparator** (if budget allows) — MRC-5, WI-38 or BEAS-2B, to give a Selectivity Index. If not possible, we will report a modelled estimate and flag it as a limitation.

## 7. Dose range — advisory to follow

Published IC₅₀ values for crude *Moringa* extracts against A549 are frequently **above 200 µg/mL**. If the true IC₅₀ exceeds the planned top dose, viability never crosses 50%, and IC₅₀ can only be obtained by extrapolation — which is not a valid determination.

We will send a **dose-range recommendation** once the composition data (§1) arrives. Please hold the final dose series until then if scheduling allows. If not, we recommend adding **400 and 800 ppm** to the existing 12.5–200 ppm series as insurance, subject to extract solubility and keeping DMSO ≤ 0.5% v/v.

---

## 8. Raw data format we need back after the assay

Please return **raw absorbance values per well**, not pre-computed viability percentages, as a CSV:

```
run, plate, well, group, conc_ppm, abs_570, abs_630
```

Raw OD lets us apply the manuscript's formula ourselves, verify the blank subtraction, quantify plate and edge effects, and compare against the sealed prediction. Pre-averaged percentages make all of that impossible.
