# Technical appendix to the CCCBA report query

**Re:** IC50 Values of Malunggay (*Moringa oleifera*) Bark Ethanolic Extract Against A549 Cancer Cell Lines Using MTT Assay, CCCBA Laboratory, MSU-IIT, released August 2026

Every figure below is derived from the laboratory's own report and its raw
absorbance file. Nothing here uses external data. All of it is reproducible
from `src/s11_observed_data.py`, `src/s10_validation.py` and
`src/s13_solvent_check.py` in the accompanying repository.

Findings are ordered by how much they affect the reported result.

---

## 1. The two independent trials disagree systematically, not randomly

This is the most consequential finding in the report and it is not mentioned
anywhere in the document.

| Dose (ug/mL) | Trial 1 | Trial 2 | Difference |
|---|---|---|---|
| 12.5 | 105.01 | 135.59 | +30.58 |
| 25 | 98.76 | 115.88 | +17.12 |
| 50 | 97.53 | 108.92 | +11.39 |
| 100 | 80.64 | 101.67 | +21.03 |
| 200 | 70.88 | 96.23 | +25.35 |
| 400 | 55.37 | 74.71 | +19.34 |

Trial 2 reads higher at **every single dose**, with a mean offset of **+20.80
percentage points**. A paired t test against zero gives t = 7.671, p = 0.0006.
Six out of six in the same direction has a probability of 1 in 32 under chance
alone, and the magnitude makes it stronger than that.

This is a systematic between run shift, not measurement scatter. Two
consequences follow.

First, the reported mean is the midpoint of two results that do not agree.
Trial 1 alone would put the extract near 55 percent viability at the top dose.
Trial 2 alone would put it near 75 percent. The reported 65.04 describes
neither run.

Second, the standard deviations in Table 1 are not measuring precision. They
are measuring the trial offset. An SD computed from two numbers carries one
degree of freedom and is reported to two decimal places throughout the table.

The offset is not a plate gain effect. If it were, the controls would move with
it. Mean negative control absorbance fell from 1.332 in trial 1 to 1.221 in
trial 2, and mean doxorubicin absorbance fell from 0.936 to 0.828. Only the
extract wells rose, by 43 percent on the raw optical density scale.

**Question for the laboratory.** What accounts for the systematic elevation of
the MBEE wells in trial 2 relative to their own controls, and should the two
trials be pooled at all.

---

## 2. Point selection for the IC50 regression

Table 2 states that the MBEE fit used 400, 100 and 50 ug/mL, omitting 200, 25
and 12.5. No criterion is given for the omission.

There are 20 possible ways to choose 3 concentrations from 6. Fitting the same
log linear model to each gives IC50 values from **912 to 3827 ug/mL**, with a
median of 1699. The subset the laboratory used returns **911.7, the lowest of
all 20**. Using all six concentrations returns 1482.

The same pattern appears in the doxorubicin fit. Table 2 records 50, 12.5 and
6.2 ug/mL, omitting 25 and 3.1. Of the 10 possible 3 point subsets, theirs
returns the **second lowest** IC50 at 28.56 against a range of 27.74 to 66.00.
Using all five returns 34.29.

There is a straightforward and innocent explanation. Both chosen subsets are
among the most linear available. The doxorubicin subset has the highest R
squared of all ten at 0.9976, and the MBEE subset has the second highest of all
twenty at 0.9996. Selecting the points that fall closest to a straight line is
a defensible instinct.

The difficulty is that it is still selection on the outcome. Choosing which
points to fit after seeing how well they fit biases the estimate, and here it
moved the reported potency to the most favourable value available in one case
and the second most favourable in the other. Omitting 200 ug/mL, which sits
between two included points, has no principled justification we can identify.
Omitting 25 ug/mL from the doxorubicin fit removes the concentration closest to
the 50 percent crossing, which is the single most informative point in the
series.

**Question for the laboratory.** What criterion determined which
concentrations entered each regression, and was it fixed before the fits were
run.

---

## 3. The MBEE IC50 lies outside the tested range

Highest concentration tested was 400 ug/mL. Viability there was 65.04 percent.
The response never crosses 50 percent at any concentration, so no IC50 exists
within the data.

The reported 911.84 ug/mL is obtained by extending the fitted line **2.28 times
beyond the highest measured point**. A four parameter logistic fit on the same
data returns no solution at all, which is the correct behaviour when a curve has
no midpoint inside the range. The straight line returns a number regardless.

The correct statement is that the IC50 exceeds 400 ug/mL.

---

## 4. R squared is being asked to do work it cannot do

Three points fitted with a two parameter linear model leave one residual degree
of freedom. An R squared near unity in that situation indicates that three
points lie close to a line, which is close to guaranteed, and says nothing about
whether the concentration range was adequate.

Concretely, the subset (25, 100, 200) gives R squared of 0.9998, higher than the
subset actually used, and an IC50 of 3696 ug/mL, four times larger. Fit quality
and estimate validity are moving independently here.

---

## 5. The outlier criterion cannot be Grubbs

The raw data file states that red values are outliers based on the Grubbs test.
**27 of 78 wells are marked, which is 34.6 percent of the dataset.**

A two sided Grubbs test at alpha 0.05 cannot produce that at this design.

Applied per trial at n equals 3, the test flags **one well across all 78**. The
reason is structural. At three observations the test statistic can only range
from 1.0000 for three equally spaced values to 1.1547 when two coincide, while
the critical value is 1.1543. The test has almost no power to reject at n equals
3.

Applied to the pooled six wells per concentration, at n equals 6 with a critical
value of 1.8871, it flags **zero wells**.

A plain rule of discarding anything more than one standard deviation from the
group mean gives 26 flags, close to the 27 observed, and agrees with the
laboratory on 76 percent of individual wells. That rule would discard roughly a
third of any normally distributed dataset by construction.

**Question for the laboratory.** What criterion and significance level were
actually applied.

---

## 6. The blank wells show a reproducible position effect, not outliers

| Trial | R1 | R2 | R3 |
|---|---|---|---|
| 1 | 0.068 | 0.065 | **0.166** |
| 2 | 0.071 | 0.067 | **0.157** |

Blank wells contain medium and dye with no cells, so they should differ only by
pipetting error. The third replicate reads roughly **2.4 times higher than the
other two, in both trials independently**.

An artifact that reproduces across two separate runs in the same position is
systematic. It points to a plate position effect, a dispensing issue, or
carryover affecting that column. Labelling it as a statistical outlier records
the symptom and discards the evidence.

This matters beyond the blanks. Blank absorbance is subtracted from every well
in the viability calculation, so a systematic error in the blank propagates into
every reported percentage.

**Question for the laboratory.** Was a plate position or dispensing check
performed, and are the flagged wells in the sample plates also concentrated in
particular positions.

---

## 7. Solvent and dose cannot be separated in this design

The Sample Preparation section describes a 10,000 ug/mL stock made by dissolving
2 mg of extract in 200 uL of neat DMSO, with working concentrations prepared by
diluting that stock into culture medium.

If the doses were prepared directly from that stock, the final solvent fraction
is fixed by the dose.

| Dose (ug/mL) | DMSO % v/v | Viability | Above A549 limit |
|---|---|---|---|
| 12.5 | 0.125 | 120.30 | no |
| 25 | 0.25 | 107.32 | no |
| 50 | 0.5 | 103.22 | no |
| 100 | 1.0 | 91.16 | **yes** |
| 200 | 2.0 | 83.56 | **yes** |
| 400 | 4.0 | 65.04 | **yes** |

A549 is routinely maintained at or below 0.5 percent v/v and the solvent alone
commonly becomes cytotoxic above about 1 percent. The top dose would sit at 4
percent, eight times the routine ceiling.

The three doses that inhibited are exactly the three that would exceed the
limit. The three at or below the limit all read above the untreated control.

This is not a claim that DMSO produced the effect. Solvent fraction is an exact
linear function of dose in this design, so the two variables carry identical
information and **no analysis of this dataset can separate them**. Only a
vehicle control at matched solvent percentage can.

**Question for the laboratory.** Was an intermediate dilution in medium
performed before dosing, and what was the final percentage by volume of DMSO in
the top dose well.

---

## 8. The vehicle control is named in the formula but absent from the data

The Percent Cell Viability formula on the Raw Data page normalises to "Ave
Vehicle Absorbance". The Analysis section states that doxorubicin showed
decreased viability "relative to the vehicle control".

The raw data file contains two control columns, NC and Blank. There is no
vehicle column.

Either vehicle wells were run and were not supplied, or normalisation was to the
untreated negative control and the formula as printed is incorrect. Both are
fixable, but they are different corrections.

---

## 9. Treatment exposure time is not stated

The report gives the seeding density at 1.8 x 10^4 cells per well, an adhesion
period of at least 24 hours, a 4 hour dye incubation, and a 1 hour solubilisation
stand. It does not state how long cells were exposed to the extract before the
dye was added.

The 4 hour figure is the formazan development step, timed from the addition of
dye solution to wells that already contained treated cells. It is not the
exposure.

Exposure duration is a required method parameter. IC50 values are exposure
dependent and cannot be compared to published figures without it.

We note that doxorubicin returned an IC50 of 28.59 ug/mL, which is roughly 53
micromolar. Published values against A549 sit between 0.1 and 10 micromolar
depending on exposure and readout. A short exposure would account for that
offset.

---

## 10. The untreated control is reported with zero variance

Table 1 lists the untreated control as 100.00 in both trials with SD 0.00.

That is true by construction, since every well is normalised to the control
mean. But the underlying control wells do vary. Across the six negative control
wells the raw absorbances are 1.242, 1.382, 1.373, 1.166, 1.274 and 1.222, a
coefficient of variation of 6.7 percent.

Presenting the reference group with no variance understates the uncertainty on
every other row in the table, because that 6.7 percent propagates into all of
them.

---

## 11. Table 2 typographical error

The coefficient of determination for MBEE is printed as "0.0.9996".

---

## 12. No dose response statistics were performed

The report supplies means, standard deviations and a regression. It contains no
test of whether any concentration differs from the untreated control.

We ran the analysis ourselves on the supplied raw data. Because the 12.5 ug/mL
group failed a Shapiro-Wilk normality check, the analysis moves to the
Kruskal-Wallis branch, which returned H = 12.85 on 6 degrees of freedom with p =
0.0455. Dunn post hoc comparisons with Holm correction separated **no individual
concentration** from the negative control.

We raise this not as a defect in the deliverable, since the commission was for
IC50 determination, but because the result is relevant to how the data can be
described. No concentration in this dataset can be called significantly
cytotoxic on its own.

---

## Summary

| # | Finding | Effect on the reported result |
|---|---|---|
| 1 | Trials disagree systematically by 20.8 pp, p = 0.0006 | Mean describes neither run, SDs measure the offset |
| 2 | 3 of 6 points fitted, no criterion, returns lowest of 20 possible IC50 values | Reported potency is the most favourable available |
| 3 | IC50 sits 2.28x beyond the highest dose tested | Value is extrapolated, not determined |
| 4 | R squared on 3 points with 1 df | Does not validate the range |
| 5 | 27 of 78 wells excluded, criterion cannot be Grubbs | A third of the data removed by an unstated rule |
| 6 | Blank R3 elevated 2.4x in both trials | Systematic artifact recorded as random outlier |
| 7 | Solvent fraction collinear with dose, up to 4 percent DMSO | Extract and solvent effects inseparable |
| 8 | Vehicle control in the formula, absent from the data | Normalisation basis unclear |
| 9 | Exposure time not stated | Result not comparable to literature |
| 10 | Control reported with zero variance | Understates uncertainty on every row |
| 11 | "0.0.9996" | Typographical |
| 12 | No inferential statistics | No concentration is significant vs control |

---

## A note on how to read this

Items 3, 4, 9, 10, 11 and 12 are presentation and completeness issues. They are
straightforward to correct and none of them means the bench work was wrong.

Items 1, 2, 5, 6 and 7 are different. Each one affects what the numbers mean,
and together they mean the reported IC50 of 911.84 ug/mL should not be treated
as a determination.

The assay itself appears sound. Doxorubicin produced a clean dose dependent kill
and the negative controls behaved. The concerns here are with how the data was
selected, reported and described, not with whether the experiment was capable of
producing a usable answer.
