# Dose-range recommendation — ACTION REQUIRED before the MTT assay

**To:** MSU-IIT partner institution / research adviser
**From:** in silico arm
**Status:** advisory, evidence-based, issued before any wet-lab data exists

```
==============================================================================
LITERATURE PRIOR -- M. oleifera extract IC50 vs A549
==============================================================================
                extract           solvent             cell  ic50_ug_ml assay                                                      source
           aqueous leaf           aqueous             A549      166.70   MTT        Tiloke et al. 2013, BMC Complement Altern Med 13:226
alkaloid extract (MOAE) alkaloid fraction             A549      158.67   MTT Xie et al. 2021, Evid Based Complement Alternat Med 5591687
         ethanolic leaf           ethanol             A549     1062.87   MTS                      Trends in Sciences (2022), MOE vs A549
         ethanolic leaf           ethanol MCF-12A (normal)     1424.04   MTS           Trends in Sciences (2022), normal-cell comparator

Selectivity Index from the ethanolic study (MCF-12A / A549) = 1.34
  -> barely selective; worth stating plainly in the discussion.

------------------------------------------------------------------------------
PRIOR DISTRIBUTIONS (log10 IC50, ug/mL)
------------------------------------------------------------------------------
  pooled over all A549 extracts : GM    304.1  1sigma 103 - 899
  solvent-matched (ethanol only): GM   1062.9  1sigma 359 - 3143
  This study uses an ETHANOLIC extract -> the solvent-matched prior is
  the decision-relevant one; the pooled prior is the optimistic bound.

------------------------------------------------------------------------------
RISK THAT THE PLANNED RANGE FAILS (top dose = 200 ppm)
------------------------------------------------------------------------------
  P(IC50 > 200 ppm) under pooled prior          =  65.0%
  P(IC50 > 200 ppm) under solvent-matched prior =  93.8%

  If IC50 exceeds the top dose, viability never crosses 50%, the
  log-linear solution becomes an extrapolation, and Specific Question 3
  (IC50) cannot be answered. Monte Carlo E1 showed 0% of experiments
  bracket 50% inhibition once the true IC50 passes ~214 ug/mL.

------------------------------------------------------------------------------
WHAT THE PLANNED EXPERIMENT WOULD ACTUALLY SEE
------------------------------------------------------------------------------
  prior = solvent-matched  (IC50  1062.9 ug/mL)
    dose (ppm)      :    12.5     25.0     50.0    100.0    200.0
    %viability      :    99.6     99.0     97.7     94.9     89.1
    %inhibition     :     0.4      1.0      2.3      5.1     10.9
    max inhibition at 200 ppm = 10.9%

  prior = pooled           (IC50   304.1 ug/mL)
    dose (ppm)      :    12.5     25.0     50.0    100.0    200.0
    %viability      :    98.0     95.6     90.5     80.8     65.3
    %inhibition     :     2.0      4.4      9.5     19.2     34.7
    max inhibition at 200 ppm = 34.7%

------------------------------------------------------------------------------
RECOMMENDED DOSE SERIES
------------------------------------------------------------------------------
  solvent-matched : top dose ~   6323 ppm for 95% coverage
  pooled          : top dose ~   1809 ppm for 95% coverage

  PRACTICAL RECOMMENDATION (2-fold series, 6 points):
    50, 100, 200, 400, 800, 1600 ppm
    Rationale: keeps 50/100/200 ppm so the original design is nested
    inside the new one, and extends far enough that the solvent-matched
    prior is bracketed. Drop 12.5 and 25 ppm -- both priors predict
    <5% inhibition there, so those wells buy no information.

    CONSTRAINTS TO CHECK BEFORE ADOPTING:
      * extract solubility at 1600 ug/mL in the stock solvent
      * final DMSO <= 0.5% v/v in the top-dose well
      * if solubility caps below the IC50, report 'IC50 > [max soluble
        dose]' as a VALID result -- do not extrapolate the regression.

------------------------------------------------------------------------------
IF THE PLANNED RANGE IS USED ANYWAY (solvent-matched prior)
------------------------------------------------------------------------------
  true IC50                    =    1062.9 ug/mL
  log-linear estimate (paper)  = 17113009.6 ug/mL (16100.8x the truth)
  4PL estimate                 =    1062.9 ug/mL
  regression R^2               = 0.8519  <-- looks fine!
  The R^2 is high even though the IC50 is badly wrong. A good-looking
  regression is NOT evidence that the dose range was adequate.
  NCI class at the true value  : inactive (> 1000 ug/mL)
```

## Sources

- aqueous leaf vs A549: IC₅₀ = 166.7 µg/mL (MTT) — Tiloke et al. 2013, BMC Complement Altern Med 13:226. <https://link.springer.com/article/10.1186/1472-6882-13-226>
- alkaloid extract (MOAE) vs A549: IC₅₀ = 158.67 µg/mL (MTT) — Xie et al. 2021, Evid Based Complement Alternat Med 5591687. <https://onlinelibrary.wiley.com/doi/10.1155/2021/5591687>
- ethanolic leaf vs A549: IC₅₀ = 1062.87 µg/mL (MTS) — Trends in Sciences (2022), MOE vs A549. <https://tis.wu.ac.th/index.php/tis/article/view/3202>
- ethanolic leaf vs MCF-12A (normal): IC₅₀ = 1424.04 µg/mL (MTS) — Trends in Sciences (2022), normal-cell comparator. <https://tis.wu.ac.th/index.php/tis/article/view/3202>

## Limitations of this advisory

- No quantitative A549 data for *M. oleifera* **bark** could be located;
  the prior is built from **leaf** extracts. Bark differs in tannin,
  alkaloid and isothiocyanate content, so the true value may sit outside
  the prior. This is the study's novelty and also its main uncertainty.
- One of the four anchors used an **MTS** rather than MTT readout.
- The single ethanolic datapoint borrows its spread from the
  between-study variance of the other anchors.
- This advisory does **not** depend on the QSAR model, which was found
  unusable for per-compound potency (~11-fold 1σ error).