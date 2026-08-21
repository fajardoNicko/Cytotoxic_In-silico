# Request for clarification and a revised report

**To:** Cell Culture and Cell-Based Assay (CCCBA) Laboratory, MSU-IIT
**Attention:** Francis N. Limbag, Research Assistant, and Mylene M. Uy, DSc, RCh, Primary Investigator
**From:** Reme Allew Alistan et al., Valenzuela City School of Mathematics and Science
**Re:** IC50 Values of Malunggay (*Moringa oleifera*) Bark Ethanolic Extract Against A549 Cancer Cell Lines Using MTT Assay, released August 2026

---

Good day.

Thank you for conducting the assay and for supplying the raw absorbance file
along with the report. Having the per well data has been genuinely useful.

We are preparing this result for a research paper that will be defended and
reviewed, so we need the reported figures to be internally consistent and the
method fully specified. On review we found several points that we cannot
resolve from the document as issued. We are requesting clarification on each,
and a revised report where a correction is warranted.

## A. Points where the report appears inconsistent

**A1. The closing sentence of the Analysis section.**
It reads "Overall, MBEE demonstrated inhibitory activity against A549 cells at
<50 ug/ml." Table 1 reports viability of 103.22 percent at 50 ug/ml, 107.32
percent at 25 ug/ml, and 120.30 percent at 12.5 ug/ml, all above the untreated
control. Expressed as inhibition those are minus 3.22, minus 7.32 and minus
20.30 percent. We cannot reconcile the sentence with the table.

Two corrections would each make it consistent. If the intended statement was
"greater than 50 ug/ml" that matches the data, since inhibition appears at 100,
200 and 400 ug/ml. If the intended statement was that inhibition did not exceed
50 percent, that is also correct, since the maximum was 34.96 percent. Please
confirm which was meant so we quote it accurately.

**A2. Table 2 prints the coefficient of determination for MBEE as "0.0.9996".**
We assume 0.9996 is intended. Please confirm.

## B. Points where the method is not fully specified

**B1. Treatment exposure time is not stated anywhere in the report.**
The document gives the seeding density of 1.8 x 10^4 cells per well, an
adhesion period of at least 24 hours, a 4 hour dye incubation, and a 1 hour
solubilisation stand. It does not state how long the cells were exposed to the
extract before the dye solution was added. This is required for any comparison
against published IC50 values, which are exposure dependent. Please supply the
exposure duration.

**B2. Final DMSO concentration in each well is not stated.**
The Sample Preparation section describes a 10,000 ug/ml stock prepared by
dissolving 2 mg of sample in 200 uL of neat DMSO, with working concentrations
then prepared by diluting that stock with culture medium. If the working doses
were prepared directly from that stock, the resulting solvent fraction would be
4 percent v/v at 400 ug/ml, 2 percent at 200, and 1 percent at 100. A549 is
normally maintained at or below 0.5 percent v/v.

We note that the three doses showing inhibition are the same three that would
exceed that limit, and the three at or below it all read above the untreated
control. We are not asserting that solvent caused the effect. We are asking
whether an intermediate dilution in medium was performed, because if it was,
the actual solvent fraction may be far lower and this concern does not arise.
Please state the final percentage by volume of DMSO in the top dose well.

**B3. Vehicle control.**
The Percent Cell Viability formula on the Raw Data page normalises to "Ave
Vehicle Absorbance", and the Analysis section refers to doxorubicin showing
decreased viability "relative to the vehicle control". The raw data file
contains only NC and Blank columns. Please confirm whether vehicle only wells
were run, and if so supply their absorbances. If the normalisation was in fact
to the untreated negative control, please correct the formula in the report.

**B4. Outlier criterion.**
The raw data file states that values in red are outliers based on the Grubbs
test, and 27 of 78 wells are marked. We were unable to reproduce that. Applied
per trial at n equals 3, a two sided Grubbs test at alpha 0.05 flags one well
across the entire dataset, because at three observations the attainable range
of the test statistic barely reaches its own critical value. Applied to the
pooled six wells per concentration it flags none. Please state the criterion
and the significance level actually used.

## C. Points about the reported IC50

**C1. The MBEE IC50 lies outside the tested range.**
The highest concentration tested was 400 ug/ml, where viability was 65.04
percent. The response never reaches 50 percent, so the reported value of 911.84
ug/ml is obtained by extending the fitted line 2.28 times beyond the highest
measured point. We intend to report this as "IC50 greater than 400 ug/ml" with
the extrapolated figure shown separately and labelled as an extrapolation.
Please advise if you disagree with that presentation.

**C2. Three of six concentrations were used in the regression.**
Table 2 records the concentrations used for MBEE as 400, 100 and 50 ug/ml,
omitting 200, 25 and 12.5. No basis is given. Refitting the same log linear
model on all six concentrations returns approximately 1482 ug/ml rather than
911.84. Please state the criterion for point selection.

**C3. On the coefficient of determination.**
With three points and a two parameter linear model there is one residual degree
of freedom, so a high R squared is close to guaranteed and does not indicate
that the concentration range was adequate. We raise this only because the value
may be read as validating the IC50, and we would prefer the report not to be
misinterpreted on that point.

## D. What we are requesting

1. Written answers to A1, A2, B1, B2, B3, B4, C2.
2. Vehicle control absorbances, if those wells were run.
3. A revised report correcting the Analysis sentence and the Table 2 typographical
   error, and stating the exposure time, the final DMSO percentage, and the
   outlier criterion in the Experimental Section.

We would like to acknowledge the laboratory in our paper and cite the report,
so we would rather resolve these now than have them raised during review.

Thank you for your time and assistance.

---

*Supporting analysis for every point above, including the raw data audit, is
available on request.*
