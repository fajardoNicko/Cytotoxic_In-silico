# Gate G3(b) — jamovi cross-check

Gate G3 requires the simulation's statistics and the manuscript's
actual analysis software to agree **to 3 decimal places**. Part (a)
(agreement with statsmodels/scipy) is automated in
`src/test_stats_mirror.py` and passes. Part (b) is this manual check.

## 1. Open the data

Open `results/tables/G3_jamovi_input.csv` in jamovi (72 rows, 2 columns).
Set `group` to **Nominal** and `viability_pct` to **Continuous**.

## 2. Run these analyses

**ANOVA → One-Way ANOVA**
- Dependent: `viability_pct`, Grouping: `group`
- Variances: tick **Assume equal (Fisher's)**
- Additional Statistics: tick **Descriptives table**
- Assumption Checks: tick **Normality test** and **Homogeneity test**

**ANOVA → ANOVA** (for the post hoc)
- Dependent: `viability_pct`, Fixed Factors: `group`
- Post Hoc Tests: move `group` across, correction **Tukey**

## 3. Expected values

### One-way ANOVA

| Quantity | Expected |
|---|---|
| F | **235.652** |
| df (between, within) | 7, 64 |
| p | 3.73e-43 |
| Sum of Squares (group) | 27680.669 |
| Sum of Squares (residual) | 1073.960 |
| Mean Square (group) = MST | 3954.381 |
| Mean Square (residual) = MSE | 16.781 |
| η² | 0.9627 |

> jamovi's *One-Way ANOVA* panel defaults to **Welch's**. Switch to
> **Fisher's (assume equal variances)** or the F will not match — Welch
> is a different test, not a rounding difference.

### Assumption checks

| Test | Expected |
|---|---|
| Shapiro–Wilk on residuals, W | 0.974 |
| Shapiro–Wilk on residuals, p | 0.147 |
| Levene's (center = median), F | 2.136 |
| Levene's, p | 0.052 |

> jamovi's homogeneity test uses `car::leveneTest`, which centres on the
> **median** (Brown–Forsythe). `stats_mirror.levene()` matches that by
> default. Centring on the mean gives a different number.

### Descriptives (mean ± SD)

| Group | n | Mean | SD |
|---|---|---|---|
| negative_control | 9 | 100.000 | 3.858 |
| vehicle_control | 9 | 95.746 | 6.460 |
| doxorubicin | 9 | 50.388 | 2.596 |
| 12.5ppm | 9 | 90.100 | 4.542 |
| 25ppm | 9 | 86.544 | 4.708 |
| 50ppm | 9 | 76.817 | 4.172 |
| 100ppm | 9 | 63.768 | 2.669 |
| 200ppm | 9 | 44.480 | 1.887 |

### Tukey HSD — critical value

Balanced design, so the manuscript's closed form applies:

    HSD = q(0.05, k=8, df=64) × √(MSE/n) = **6.0507**

Any absolute mean difference above this is significant.

### Tukey HSD — comparisons vs. negative control

| Comparison | Mean diff | p (Tukey) |
|---|---|---|
| negative_control vs vehicle_control | 4.254 | 0.3640 |
| negative_control vs doxorubicin | 49.612 | 0.0000 |
| negative_control vs 12.5ppm | 9.900 | 0.0001 |
| negative_control vs 25ppm | 13.456 | 0.0000 |
| negative_control vs 50ppm | 23.183 | 0.0000 |
| negative_control vs 100ppm | 36.232 | 0.0000 |
| negative_control vs 200ppm | 55.520 | 0.0000 |

## 5. Non-parametric branch — also verify

On this dataset the manuscript's own decision rule actually selects **Kruskal-Wallis+Dunn** (Shapiro–Wilk on residuals p = 0.1472).
This is not a defect — with 9 wells per group the normality test
rejects a noticeable fraction of the time, so the fallback path in the
research plan is a live possibility and must be verified too.

**ANOVA → One-Way ANOVA (Non-parametric) / Kruskal-Wallis**

| Quantity | Expected |
|---|---|
| χ² (H) | **66.171** |
| df | 7 |
| p | 8.74e-12 |
| ε² | 0.9245 |

> jamovi's non-parametric post hoc is **DSCF** (Dwass-Steel-Critchlow-
> Fligner), whereas `stats_mirror` implements **Dunn**. These are
> different procedures and will NOT match — that is expected, not a
> failure. Compare the Kruskal-Wallis H and p only, and state in the
> manuscript which post hoc was used.

## 4. Pass criterion

Every value above must match jamovi to **3 decimal places**.
Sign conventions on mean differences may be reversed (jamovi reports
group2 − group1); magnitudes and p-values must agree.

If anything disagrees, do **not** proceed to Phase F — the Monte Carlo
power estimates would not describe the analysis actually being run.

Record the outcome, the jamovi version, and the date here:

| Date | jamovi version | Result | Checked by |
|---|---|---|---|
|  |  |  |  |