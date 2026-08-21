# MATERIALS AND METHODS

**Study.** In silico prediction of the cytotoxic activity of *Moringa oleifera* ethanolic bark extract against A549 human lung adenocarcinoma cells

**Companion to.** "In Vitro Cytotoxic Activity of Malunggay's (*Moringa oleifera*) Ethanolic Bark Extract against A549 Lung Cancer Cells" (Dalistan, Roque, Soriano, Salvilla)

**Version.** 2.0, 20 August 2026

## Collection of Dataset

This study employed a fully computational methodology and required no cell culture, extraction, or laboratory assay. All inputs were drawn from public chemical, structural, and bioactivity repositories, and every retrieval was cached to disk so that each analysis reruns without further network access.

Constituent structures were obtained from PubChem through the PUG-REST interface (Kim et al., 2023). The candidate list was compiled from the phytochemical literature on *M. oleifera*, weighted toward bark-associated compounds, and covers glucosinolates and isothiocyanates, flavonoids, phenolic acids, sterols, triterpenoids, alkaloids, and long-chain fatty acids. Forty-seven compounds resolved to a parseable canonical SMILES and were retained.

Protein structures were obtained from the RCSB Protein Data Bank through its Search and Data APIs (Burley et al., 2023). Sixteen molecular targets were queried, chosen from the documented A549 genotype rather than from generic oncology relevance. A549 carries *KRAS* G12S, is *STK11* null, retains wild-type *TP53*, has a *CDKN2A* deletion, and carries the *KEAP1* G333C mutation that leaves NRF2 constitutively active (Singh et al., 2006). That last lesion matters directly here, because *M. oleifera* isothiocyanates are established NRF2 activators and would be acting on an axis that is already de-repressed.

Cytotoxicity data for the quantitative structure-activity relationship (QSAR) model were retrieved from ChEMBL through its REST API (Zdrazil et al., 2024) using the A549 cell-line target identifier. A parallel normal-lung set covering MRC-5, WI-38, BEAS-2B, and IMR-90 was retrieved for selectivity modelling. The A549 pull terminated at roughly 21,000 raw records when the connection dropped, and curation ran on what had been retrieved. The fetch routine is re-runnable and extends the set.

Extract-level potency anchors were collected separately from the published literature. Four quantitative measurements of crude *M. oleifera* extract against A549 or a normal-cell comparator were located, each recorded with full bibliographic provenance, extraction solvent, readout, and value. The anchors are Tiloke et al. (2013), Xie et al. (2021), and Trends in Sciences (2022), the last supplying both an A549 value and a normal-cell comparator.

Quantitative chemical characterization of the extract batch was requested from the partner institution and is not available. Relative peak area and PubChem CID were confirmed as excluded from the deliverable, so no GC-MS or LC-MS composition table exists for the tested material. Two consequences follow and are stated wherever they apply. Constituent-level analyses are literature-based rather than batch-verified, and the composition vector required by any bottom-up mixture calculation cannot be constructed.

## Software Development

The software environment consisted of Python 3.13 with NumPy (Harris et al., 2020), SciPy (Virtanen et al., 2020), pandas (McKinney, 2010), scikit-learn (Pedregosa et al., 2011), XGBoost (Chen & Guestrin, 2016), RDKit (Landrum, 2024), statsmodels (Seabold & Perktold, 2010), Matplotlib (Hunter, 2007), and requests. These packages together provide cheminformatics, descriptor calculation, ensemble and gradient-boosted regression, statistical inference, and plotting. Tabular outputs were written as comma-separated value files and the sealed prediction as JavaScript Object Notation. Molecular docking uses AutoDock Vina (Eberhardt et al., 2021), with CB-Dock2 (Liu et al., 2022) available as a no-install substitute, and ADMET profiling uses SwissADME (Daina et al., 2017), ADMETlab 3.0 (Fu et al., 2024), and ProTox-3.0 (Banerjee et al., 2024). The route taken is recorded per result.

The pipeline runs in four stages across eleven modules. The first acquires and curates external data, covering ligand preparation, target verification, and ChEMBL retrieval. The second builds and validates the potency model, covering QSAR training, a reality check, and calibration. The third produces the extract-level prediction and the assay simulation, covering the literature prior, the mixture model, the dose-response mathematics, the virtual plate, and the Monte Carlo experiments. The fourth seals the prediction and later scores it against the wet-lab data. A random seed of 20260729 is fixed throughout, and the git commit hash is written into every sealed result.

Five decision gates were declared before their corresponding analyses ran. G1 requires that redocking of each co-crystallised ligand reproduce the crystallographic pose within 2.0 Å root-mean-square deviation. G2 requires 5-fold cross-validated Q² of at least 0.5, external test R² of at least 0.6, and y-randomised Q² below 0.2. G3 requires that Python and jamovi agree to three decimal places on the same input file. G4 requires that the predicted doxorubicin A549 IC₅₀ fall within threefold of published values with doxorubicin held out of training. G5 requires that the predicted extract IC₅₀ fall within threefold of the observed value, that percent viability root-mean-square error not exceed 15 percentage points, and that Spearman ρ on dose rank order reach 0.9. Each gate carries a stated consequence for failure, and outcomes are reported whether or not the gate passed.

## Ligand Library Construction

Retrieved SMILES were parsed in RDKit and stripped to the largest organic fragment, which removes counter-ions and solvates that would otherwise distort descriptor values. Three-dimensional conformers were generated with ETKDGv3 (Wang et al., 2020) and minimised under the MMFF94 force field (Halgren, 1996). Molecules were protonated for physiological pH 7.4. Each compound was written to a structure-data file for docking and to a tabular library carrying molecular weight, calculated logP, topological polar surface area, hydrogen-bond donor and acceptor counts, rotatable-bond count, ring and aromatic-ring counts, and fraction of sp³ carbon.

## Target Panel Verification

Candidate entries for each target were filtered on human origin, experimental method, resolution, presence of a bound drug-like ligand, correct gene assignment, and absence of confounding nucleic acid. Cofactors such as FAD, heme, and ATP analogues were excluded from the drug-like ligand count, because their presence does not establish an inhibitor pocket. Surviving candidates were ranked by resolution and the highest-resolution entry retained. Every selection carries a recorded accept or reject verdict with its reason.

Fifteen of sixteen targets were accepted. Thioredoxin reductase 1 was rejected, because every candidate entry carried only FAD and offered no inhibitor site against which to centre a docking box or run the G1 redocking check. Dropping the target and saying so is preferable to docking into an unvalidated pocket. Three accepted entries carry documented caveats requiring manual confirmation before docking, since resolution alone proved an insufficient selection rule. The Bcl-2 entry annotates both *BCL2* and *BCL2L1*, the mTOR entry is an FRB-rapamycin-FKBP complex rather than the kinase domain, and the tubulin entry needs confirmation as colchicine-site.

## Curation of the QSAR Training Set

Raw ChEMBL records were filtered in a fixed order, with each step logged alongside the count it removed. Activity types were restricted to IC₅₀, GI₅₀, and EC₅₀ under a standard relation of equality, so censored values reported only as above or below a bound were discarded rather than imputed. Units had to be convertible to nanomolar, and a parseable canonical SMILES had to exist. Salts and mixtures were stripped to the largest organic fragment. Activities were then standardised to a logarithmic potency scale.

$$pIC50 = 9 - log10(IC50 in nM)$$

where:

- pIC₅₀ = the negative base-10 logarithm of the molar half-maximal inhibitory concentration
- IC₅₀ = the reported half-maximal inhibitory concentration converted to nanomolar

Records were deduplicated by InChIKey, retaining the median pIC₅₀ across replicates. Compounds whose replicate spread exceeded one log unit were dropped, because an irreproducible label injects noise that the model then learns. Values outside the interval 3 to 11 were removed as probable unit-entry errors. Curation returned 6,972 compounds with A549 activity.

## Molecular Featurisation

Each compound was represented by a Morgan circular fingerprint of radius 2 and 1024 bits (Rogers & Hahn, 2010) concatenated with the interpretable physicochemical descriptors listed above. Fingerprint length was set to 1024 rather than the more common 2048 after testing showed no gain in scaffold-split R² at roughly double the cost of every forest fit, which matters because the validation protocol requires dozens of fits.

Chemical similarity between a query compound and the training set was measured by the Tanimoto coefficient on those fingerprints.

$$T(A, B) = |A ∩ B| / |A ∪ B|$$

where:

- T(A, B) = the Tanimoto similarity between fingerprints A and B
- |A ∩ B| = the number of bits set in both fingerprints
- |A ∪ B| = the number of bits set in either fingerprint

The applicability domain was defined by the maximum Tanimoto similarity of a query compound to any training compound. Predictions falling below 0.3 were flagged out of domain and are not reported as quantitative estimates.

## QSAR Model Training and Validation

Three regressors were trained. Ridge regression served as the linear baseline, with Random Forest (Breiman, 2001) and XGBoost as the nonlinear models. Each was fitted under two independent train-test partitions of the same curated set.

The random partition assigns 80 percent of compounds to training and 20 percent to testing at random. It is the optimistic protocol and the one most commonly reported. The scaffold partition follows Bemis and Murcko (1996), grouping compounds by their ring-and-linker framework and assigning whole scaffold groups to one side of the split, so no test compound shares a framework with any training compound. The scaffold partition governs interpretation here, because the *Moringa* constituents are natural-product frameworks thinly represented in the ChEMBL cytotoxicity corpus.

Model quality was reported as cross-validated Q², external test R², and root-mean-square error.

$$Q² = 1 - Σ(yᵢ - ŷᵢ,cv)² / Σ(yᵢ - ȳ)²$$

$$R² = 1 - Σ(yᵢ - ŷᵢ)² / Σ(yᵢ - ȳ)²$$

$$RMSE = √[ Σ(yᵢ - ŷᵢ)² / n ]$$

where:

- yᵢ = the measured pIC₅₀ of compound i
- ŷᵢ = the predicted pIC₅₀ of compound i
- ŷᵢ,cv = the prediction for compound i from the fold in which it was held out
- ȳ = the mean measured pIC₅₀
- n = the number of compounds in the evaluation set

A y-randomisation control ran alongside, permuting the labels and repeating the full fitting procedure. A model that fits noise rather than structure returns a randomised Q² near or above zero, so this check separates signal from overfitting.

The positive-control anchor of gate G4 was run with doxorubicin and every reference cytotoxic removed from the training set beforehand. Without that removal the model recalls a memorised label and the anchor establishes nothing. The held-out prediction was converted from the logarithmic scale to mass units for comparison with the assay.

$$IC50 (µg/mL) = IC50 (µM) × Mr / 1000$$

where:

- Mᵣ = the relative molecular mass of the compound in g/mol

The converted value was compared against the published doxorubicin A549 window of 0.1 to 10 µM, which spans reported exposure times and readouts rather than a single literature point.

## Calibration and the Decision on Per-Compound Potency

Because gate G2 failed on the scaffold split, the error structure of the model was characterised before any downstream use. The scaffold-split test set was restricted to compounds falling inside the *Moringa* physicochemical envelope, defined by the 5th to 95th percentile range of the constituent library in molecular weight, calculated logP, topological polar surface area, and fraction sp³ carbon. Within that envelope the systematic offset was estimated as a mean signed residual.

$$b = (1/n) Σ (yᵢ - ŷᵢ)$$

where:

- b = the mean signed bias in log units, positive when the model under-predicts potency
- n = the number of envelope-matched test compounds

The same quantity was estimated a second time from the five *Moringa* constituents holding genuine measured A549 values, namely benzyl isothiocyanate, ursolic acid, apigenin, stigmasterol, and oleanolic acid, all held out of training. The two estimates were then compared for agreement in sign and in magnitude against the residual spread. Their disagreement is the basis for the methodological decision recorded below.

## Extract-Level Potency Prediction

### Mixture model, retained and conditional

The original design combined per-compound potencies into an extract IC₅₀ under concentration addition (Loewe, 1928),

$$IC50,extract = 1 / Σ (pᵢ / IC50,ᵢ)$$

with Bliss independence (Bliss, 1939) as the secondary estimate.

$$Vmix(C) = Π Vᵢ(pᵢ C)$$

where:

- pᵢ = the mass fraction of constituent i in the extract
- IC₅₀,i = the half-maximal inhibitory concentration of constituent i in µg/mL
- C = the total extract concentration in µg/mL
- V(C) = the fractional viability at concentration C

Both estimators require a composition vector and per-compound potencies of usable precision. Neither condition holds. The composition vector cannot be built, since the batch was not characterised quantitatively, and the calibration above showed the per-compound estimates to carry roughly elevenfold error at one standard deviation with an unresolved bias direction. Summing twenty such values would produce a precise-looking number carrying no information. This machinery is therefore retained in the codebase and left unused, becoming usable only if measured IC₅₀ values are obtained for the major constituents. The QSAR model itself is kept for ranking and triage, meaning the choice of which constituents merit docking effort, and is used nowhere for absolute potency.

### Literature prior, adopted

Extract potency is instead estimated from published measurements of crude *M. oleifera* extract against A549. Anchors were placed on a base-10 logarithmic scale and summarised as a geometric mean with a one-sigma interval.

$$µ = (1/k) Σ log10 cⱼ$$

$$σ = √[ Σ (log10 cⱼ - µ)² / (k - 1) ]$$

$$GM = 10^µ$$

where:

- cⱼ = the published IC₅₀ of anchor j in µg/mL
- k = the number of anchors entering the prior
- GM = the geometric mean potency

Two priors were built. The pooled prior draws on all A549 anchors regardless of extraction solvent and forms the optimistic bound. The solvent-matched prior is restricted to ethanolic extracts and is the decision-relevant distribution, because the material under study is an ethanolic extract. Since that prior rests on a single anchor, its dispersion was borrowed from the between-study variance of the remaining anchors, and the borrowing is stated as a limitation rather than concealed inside the interval.

The probability that true potency exceeds a candidate top dose follows from the log-normal form.

$$P(IC50 > d) = 1 - Φ[ (log10 d - µ) / σ ]$$

where:

- d = the candidate top dose in µg/mL
- Φ = the standard normal cumulative distribution function

This route depends on no part of the QSAR model. It is anchored in measurements of the same class of material, on the same cell line, under the same readout, so it is directly comparable to what the MTT assay produces.

### Dose-range adequacy

The prior was evaluated against the planned series of 12.5, 25, 50, 100, and 200 ppm. Four quantities were computed. The first is the probability that true IC₅₀ exceeds the planned top dose under each prior. The second is the percent viability and percent inhibition the planned series would actually produce at the prior geometric mean. The third is the top dose required for 95 percent coverage of the prior. The fourth is the error the manuscript's log-linear estimator incurs when the range is inadequate, reported next to the coefficient of determination that estimator would nonetheless display. Predicted potency was then classified against the National Cancer Institute criterion for crude plant extracts, which treats an IC₅₀ at or below 20 to 30 µg/mL as active (Suffness & Pezzuto, 1990), and the expected verdict stated in advance so it cannot be fitted to the observed result afterward.

## Molecular Docking

Receptors were prepared by removing water and non-essential heteroatoms while retaining structurally required cofactors and metal ions, adding polar hydrogens, and assigning Gasteiger charges (Gasteiger & Marsili, 1980). Grid boxes were centred on the crystallographic ligand of each structure with at least 5 Å of padding. Before any constituent was docked, each co-crystallised ligand was extracted, randomised, and redocked into its own receptor, and only targets meeting the 2.0 Å criterion of gate G1 were carried forward.

Production docking runs each library compound and doxorubicin against each validated target at exhaustiveness 16 to 32, retaining nine poses per run, in triplicate under different random seeds. Affinities are reported as mean and standard deviation in kcal/mol, and protein-ligand interaction profiles covering hydrogen bonds, hydrophobic contacts, π-stacking, and salt bridges are generated for the top-ranked complexes with PLIP (Adasme et al., 2021). Output is read as a ranking and as a source of mechanistic hypotheses. Binding affinity is not evidence of cellular activity, and that limit is restated wherever docking results appear.

## ADMET and Network Pharmacology

Each constituent was profiled for Lipinski and Veber compliance, topological polar surface area, calculated logP, predicted gastrointestinal absorption, blood-brain barrier permeation, cytochrome P450 inhibition, and P-glycoprotein substrate status. Predicted oral LD₅₀, Globally Harmonized System toxicity class, hepatotoxicity, carcinogenicity, mutagenicity, and cytotoxicity were obtained from ProTox-3.0. Predicted protein targets were retrieved from SwissTargetPrediction (Daina et al., 2019) and intersected with lung-adenocarcinoma gene sets compiled from GeneCards (Stelzer et al., 2016), DisGeNET (Piñero et al., 2020), and TCGA-LUAD differential expression. The intersection was submitted to STRING (Szklarczyk et al., 2023) for protein-protein interaction network construction, hub genes were ranked by degree and betweenness centrality, and functional enrichment over KEGG pathways (Kanehisa & Goto, 2000) and Gene Ontology terms (Gene Ontology Consortium, 2023) was tested under Benjamini-Hochberg correction (Benjamini & Hochberg, 1995) at q below 0.05.

## Virtual MTT Assay

### Plate layout

The 96-well experiment is simulated well by well. Groups comprise cell-free blank wells holding medium and MTT, untreated negative-control wells, vehicle-control wells holding cells with solvent at the highest final concentration used, positive-control wells receiving doxorubicin, one group per extract dose, and cell-free interference wells holding extract, medium, and MTT at every dose. The last two groups correspond to controls absent from the original research plan.

The perimeter is left empty by default, since edge wells read high through evaporation. When a requested replicate structure no longer fits the 60 interior wells, the simulator reports that the layout has been forced onto the perimeter. That is a real experimental constraint rather than a software limit.

Parameters left unspecified by the original plan were pinned to defensible values. Seeding density is 5 × 10³ to 1 × 10⁴ cells per well, attachment runs 24 h, treatment exposure runs 48 h, MTT is applied at 0.5 mg/mL for 4 h following Mosmann (1983), formazan is solubilised in DMSO, absorbance is read at 570 nm against a 630 nm reference, and final solvent is held at or below 0.5 percent v/v.

### Generative model

Ground-truth viability follows a four-parameter logistic function.

$$V(C) = Vmin + (Vmax - Vmin) / [ 1 + (C / IC50)^H ]$$

where:

- V(C) = percent viability at concentration C
- V_max = the upper plateau, fixed at 100 percent
- V_min = the lower plateau, set to 8 percent as typical of crude extracts
- IC₅₀ = the midpoint parameter taken from the potency prior
- H = the Hill slope, sampled over 0.8 to 1.8

Solvent cytotoxicity multiplies on top of the extract effect in every treated well. Simulated absorbance for one well is the sum of a blank term, a cell-derived signal, a chemical interference term, an edge term, and reader noise.

$$A = Ablank + (Actrl - Ablank) · V(C)/100 · (1 + εw) · gp + k C + δe + η$$

where:

- A_blank = absorbance of medium with MTT and no cells, default 0.05
- A_ctrl = absorbance of untreated cells at 48 h, default 1.00
- ε_w = well-level pipetting and seeding error, normal with standard deviation 0.06
- g_p = plate gain drawn once per independent biological run, normal with mean 1 and standard deviation 0.08
- k = the interference coefficient in OD units per µg/mL of extract
- δ_e = edge bias applied to perimeter wells only, default 0.08 OD
- η = plate-reader noise, normal with standard deviation 0.005

The interference term applies whether or not cells are present. That property is exactly what makes it recoverable by cell-free wells and invisible without them.

### Analysis of simulated data

Simulated absorbances pass through the pipeline the manuscript specifies.

$$%V = [ (Asample - Ablank) / (Acontrol - Ablank) ] × 100$$

$$%Cytotoxicity = 100 - %V$$

Percent viability is then regressed on log₁₀ concentration and the fitted line solved at 50 percent.

$$IC50 = 10^[ (50 - b) / m ]$$

where:

- m = the fitted regression slope
- b = the fitted regression intercept

A four-parameter logistic fit is computed in parallel as a sensitivity analysis, since a straight line through a sigmoidal response is defensible only across the near-linear 20 to 80 percent region. Where the log-linear solution falls outside the tested concentration range it is flagged as an extrapolation and not reported as a determination.

Because the generating parameters are known in simulation, the estimand is defined explicitly. It is the concentration at which measured viability equals 50 percent, solved numerically by Brent's method (Brent, 1973). That quantity differs from the logistic midpoint parameter, because the lower plateau sits above zero and because normalising to a medium-only control folds vehicle cytotoxicity into the apparent extract effect. Keeping the two separate is what allows the bias estimates below to mean anything.

## Monte Carlo Simulation

Each Monte Carlo experiment repeats the full simulate-and-analyse cycle with fresh pseudorandom noise across a grid of conditions, reporting the distribution of the resulting estimates rather than a single run. Five experiments were run.

E1 sweeps true IC₅₀ from 25 to 1000 µg/mL and reports the proportion of experiments that bracket 50 percent inhibition, the proportion returning an extrapolated IC₅₀, and the median bias and coefficient of variation of both estimators. E2 crosses technical replicates of 3, 4, 6, and 8 with 1 and 3 independent biological runs at two plausible true potencies, then reports analysis-of-variance power, per-dose Tukey power against the negative control, and the probability that Shapiro-Wilk rejects normality, which is the probability that the real analysis lands on the Kruskal-Wallis branch. E3 sweeps vehicle viability from 100 to 80 percent and compares IC₅₀ estimated against a medium-only control with IC₅₀ estimated against a vehicle control. E4 sweeps the interference coefficient from 0 to 0.0020 OD per µg/mL with and without cell-free subtraction. E5 compares a layout using the perimeter against one leaving it empty.

Trial budgets are 1200 iterations per condition for E1 and E2 and 600 for E3 through E5. Following the reporting framework of Morris et al. (2019), each budget is justified against the precision of the quantity being estimated rather than by convention. For a proportion such as power or the bracketing rate, the Monte Carlo standard error takes the binomial form.

$$SE(P̂) = √[ P̂ (1 - P̂) / Ntrials ]$$

where:

- P̂ = the estimated proportion at that condition
- N_trials = the number of iterations at that condition

At the worst case of P̂ = 0.5 this gives 0.0144 at 1200 iterations and 0.0204 at 600.

A sensitivity analysis sweeps true IC₅₀ across a tenfold band, well coefficient of variation from 3 to 15 percent, Hill slope from 0.8 to 1.8, and unidentified extract mass fraction from 0 to 60 percent, then ranks the influence of each on the conclusion. This identifies in advance which single assumption, if wrong, breaks the prediction, and it is the reference against which any failed prediction is diagnosed.

## Statistical Analysis

The manuscript's analytical chain was reimplemented in Python so that simulated and real data pass through identical code. Normality within each group is tested by Shapiro-Wilk (Shapiro & Wilk, 1965).

$$W = (Σ aᵢ x₍ᵢ₎)² / Σ (xᵢ - x̄)²$$

where:

- W = the normality test statistic, approaching 1 for normal data
- aᵢ = the Shapiro-Wilk weighting constants
- x₍ᵢ₎ = the i-th order statistic of the group
- x̄ = the group mean

Homogeneity of variance is tested by Levene's statistic (Levene, 1960). Where both assumptions hold, group means are compared by one-way analysis of variance with the sum-of-squares decomposition written out as in the research plan.

$$F = [ SStreat / (k - 1) ] / [ SSerror / (N - k) ]$$

where:

- SS_treat = the between-group sum of squares
- SS_error = the within-group sum of squares
- k = the number of groups
- N = the total number of observations

Pairwise comparisons against the negative control use Tukey's honestly significant difference test (Tukey, 1949).

$$q = |x̄A - x̄B| / √[ (MSerror / 2)(1/nA + 1/nB) ]$$

where:

- q = the studentised range statistic for the pair
- x̄_A, x̄_B = the two group means
- MS_error = the within-group mean square
- n_A, n_B = the two group sizes

The statistic is evaluated against the studentised range distribution on k and N minus k degrees of freedom. Where normality is rejected the analysis moves to the Kruskal-Wallis test (Kruskal & Wallis, 1952) with Dunn's post hoc comparisons (Dunn, 1964), and where variances are unequal but the data remain normal it moves to Games-Howell (Games & Howell, 1976). All tests use α = 0.05, and results are reported as mean and standard deviation per group.

Verification ran in two stages. Every statistic was first cross-checked against independent implementations in SciPy and statsmodels, which is gate G3(a). A single simulated dataset was then exported in jamovi-ready form together with a reference output file, the same analysis was run manually in jamovi, and agreement to three decimal places was required, which is gate G3(b).

## Pre-Registration

Before any wet-lab data existed, the study's predictions were frozen in a machine-readable registry holding the predicted extract IC₅₀ as a point estimate with interval under each prior, predicted percent viability with dispersion at each planned dose, the predicted dose rank order, the predicted analysis-of-variance verdict, the predicted Tukey significance pattern against the negative control, the predicted selectivity index, and the recommended replicate structure. The file was hashed under SHA-256 (Dang, 2015) and the digest, timestamp, and git commit recorded in a ledger and transmitted to the research adviser. The prediction file is never edited. Any revision is issued as a new version carrying its own digest and a written justification. Without this step, a later claim that the computational results agreed with the experimental results would be unfalsifiable.

## Validation Against the Observed Assay

When raw per-well absorbances are returned they pass through the same statistics mirror used on the simulated data, and predicted values are compared with observed on five axes. Agreement in potency is scored by fold error.

$$FE = IC50,predicted / IC50,observed$$

Success is defined as a value between one third and three. Agreement in the response curve is scored by root-mean-square error and mean absolute error on percent viability across the dose series, with success at 15 percentage points or better. Agreement in both correlation and offset is scored by Lin's concordance correlation coefficient (Lin, 1989).

$$ρc = 2 sxy / [ sx² + sy² + (x̄ - ȳ)² ]$$

where:

- s_xy = the covariance of predicted and observed viability
- s²_x, s²_y = their variances
- x̄, ȳ = their means

Agreement in ordering is scored by Spearman ρ on dose rank order, with success at 0.9. Concentration-dependent bias is exposed by Bland-Altman analysis (Bland & Altman, 1986), plotting the difference between predicted and observed viability against their mean. Predicted and observed dose-response curves are overlaid in a single figure.

A prediction that fails is diagnosed against the sensitivity analysis and attributed to a named cause, whether the potency prior, the applicability domain, the additivity assumption, MTT interference, solvent effects, or NRF2-mediated resistance. A well-diagnosed miss carries more information than an unexamined agreement.

## Assumptions and Delimitations

Constituent-level analyses are literature-based rather than batch-verified, because no quantitative composition data exists for the material tested. Docking, ADMET, and network results therefore describe *M. oleifera* chemistry as reported in the literature.

No published A549 measurement exists for *M. oleifera* bark. The potency prior rests on leaf and fractionated extracts, and bark differs materially in tannin, alkaloid, and isothiocyanate content, so the true value may fall outside the prior. That gap is the study's novelty and its principal uncertainty at once. One anchor used an MTS rather than MTT readout, and the single ethanolic anchor borrows its dispersion from the between-study variance of the others.

The QSAR model is not usable for absolute potency and is not used for it. The ChEMBL training pull is partial. Docking scores are not evidence of cellular activity. Simulation parameters are informed estimates rather than measurements of the partner laboratory's plate behaviour, and their influence is bounded by the sensitivity analysis rather than assumed away. In vitro cytotoxicity does not imply clinical efficacy, and nothing here should be read as a therapeutic claim.

## REFERENCES

Adasme, M. F., Linnemann, K. L., Bolz, S. N., Kaiser, F., Salentin, S., Haupt, V. J., & Schroeder, M. (2021). PLIP 2021: Expanding the scope of the protein-ligand interaction profiler to DNA and RNA. *Nucleic Acids Research, 49*(W1), W530-W534.

Banerjee, P., Kemmler, E., Dunkel, M., & Preissner, R. (2024). ProTox 3.0: A webserver for the prediction of toxicity of chemicals. *Nucleic Acids Research, 52*(W1), W513-W520.

Bemis, G. W., & Murcko, M. A. (1996). The properties of known drugs. 1. Molecular frameworks. *Journal of Medicinal Chemistry, 39*(15), 2887-2893.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B, 57*(1), 289-300.

Bland, J. M., & Altman, D. G. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. *The Lancet, 327*(8476), 307-310.

Bliss, C. I. (1939). The toxicity of poisons applied jointly. *Annals of Applied Biology, 26*(3), 585-615.

Breiman, L. (2001). Random forests. *Machine Learning, 45*(1), 5-32.

Brent, R. P. (1973). *Algorithms for minimization without derivatives*. Prentice-Hall.

Burley, S. K., Bhikadiya, C., Bi, C., Bittrich, S., Chao, H., Chen, L., ... Zardecki, C. (2023). RCSB Protein Data Bank (RCSB.org): Delivery of experimentally determined PDB structures alongside one million computed structure models of proteins from artificial intelligence and machine learning. *Nucleic Acids Research, 51*(D1), D488-D508.

Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, 785-794.

Daina, A., Michielin, O., & Zoete, V. (2017). SwissADME: A free web tool to evaluate pharmacokinetics, drug-likeness and medicinal chemistry friendliness of small molecules. *Scientific Reports, 7*, 42717.

Daina, A., Michielin, O., & Zoete, V. (2019). SwissTargetPrediction: Updated data and new features for efficient prediction of protein targets of small molecules. *Nucleic Acids Research, 47*(W1), W357-W364.

Dang, Q. H. (2015). *Secure hash standard* (Federal Information Processing Standards Publication 180-4). National Institute of Standards and Technology.

Dunn, O. J. (1964). Multiple comparisons using rank sums. *Technometrics, 6*(3), 241-252.

Eberhardt, J., Santos-Martins, D., Tillack, A. F., & Forli, S. (2021). AutoDock Vina 1.2.0: New docking methods, expanded force field, and Python bindings. *Journal of Chemical Information and Modeling, 61*(8), 3891-3898.

Fu, L., Shi, S., Yi, J., Wang, N., He, Y., Wu, Z., ... Cao, D. (2024). ADMETlab 3.0: An updated comprehensive online ADMET prediction platform enhanced with broader coverage, improved performance, API functionality and decision support. *Nucleic Acids Research, 52*(W1), W422-W431.

Games, P. A., & Howell, J. F. (1976). Pairwise multiple comparison procedures with unequal n's and/or variances: A Monte Carlo study. *Journal of Educational Statistics, 1*(2), 113-125.

Gasteiger, J., & Marsili, M. (1980). Iterative partial equalization of orbital electronegativity: A rapid access to atomic charges. *Tetrahedron, 36*(22), 3219-3228.

Gene Ontology Consortium. (2023). The Gene Ontology knowledgebase in 2023. *Genetics, 224*(1), iyad031.

Halgren, T. A. (1996). Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94. *Journal of Computational Chemistry, 17*(5-6), 490-519.

Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., ... Oliphant, T. E. (2020). Array programming with NumPy. *Nature, 585*(7825), 357-362.

Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science and Engineering, 9*(3), 90-95.

Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27-30.

Kim, S., Chen, J., Cheng, T., Gindulyte, A., He, J., He, S., ... Bolton, E. E. (2023). PubChem 2023 update. *Nucleic Acids Research, 51*(D1), D1373-D1380.

Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association, 47*(260), 583-621.

Landrum, G. (2024). *RDKit: Open-source cheminformatics*. https://www.rdkit.org

Levene, H. (1960). Robust tests for equality of variances. In I. Olkin (Ed.), *Contributions to probability and statistics* (pp. 278-292). Stanford University Press.

Lin, L. I. (1989). A concordance correlation coefficient to evaluate reproducibility. *Biometrics, 45*(1), 255-268.

Liu, Y., Yang, X., Gan, J., Chen, S., Xiao, Z. X., & Cao, Y. (2022). CB-Dock2: Improved protein-ligand blind docking by integrating cavity detection, docking and homologous template fitting. *Nucleic Acids Research, 50*(W1), W159-W164.

Loewe, S. (1928). Die quantitativen Probleme der Pharmakologie. *Ergebnisse der Physiologie, 27*, 47-187.

McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56-61.

Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation studies to evaluate statistical methods. *Statistics in Medicine, 38*(11), 2074-2102.

Mosmann, T. (1983). Rapid colorimetric assay for cellular growth and survival: Application to proliferation and cytotoxicity assays. *Journal of Immunological Methods, 65*(1-2), 55-63.

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., ... Duchesnay, E. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research, 12*, 2825-2830.

Piñero, J., Ramírez-Anguita, J. M., Saüch-Pitarch, J., Ronzano, F., Centeno, E., Sanz, F., & Furlong, L. I. (2020). The DisGeNET knowledge platform for disease genomics: 2019 update. *Nucleic Acids Research, 48*(D1), D845-D855.

Rogers, D., & Hahn, M. (2010). Extended-connectivity fingerprints. *Journal of Chemical Information and Modeling, 50*(5), 742-754.

Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and statistical modeling with Python. *Proceedings of the 9th Python in Science Conference*, 92-96.

Shapiro, S. S., & Wilk, M. B. (1965). An analysis of variance test for normality (complete samples). *Biometrika, 52*(3-4), 591-611.

Singh, A., Misra, V., Thimmulappa, R. K., Lee, H., Ames, S., Hoque, M. O., ... Biswal, S. (2006). Dysfunctional KEAP1-NRF2 interaction in non-small-cell lung cancer. *PLoS Medicine, 3*(10), e420.

Stelzer, G., Rosen, N., Plaschkes, I., Zimmerman, S., Twik, M., Fishilevich, S., ... Lancet, D. (2016). The GeneCards suite: From gene data mining to disease genome sequence analyses. *Current Protocols in Bioinformatics, 54*(1), 1.30.1-1.30.33.

Suffness, M., & Pezzuto, J. M. (1990). Assays related to cancer drug discovery. In K. Hostettmann (Ed.), *Methods in plant biochemistry: Assays for bioactivity* (Vol. 6, pp. 71-133). Academic Press.

Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif, R., ... von Mering, C. (2023). The STRING database in 2023: Protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research, 51*(D1), D638-D646.

Tiloke, C., Phulukdaree, A., & Chuturgoon, A. A. (2013). The antiproliferative effect of *Moringa oleifera* crude aqueous leaf extract on cancerous human alveolar epithelial cells. *BMC Complementary and Alternative Medicine, 13*, 226.

Trends in Sciences. (2022). Cytotoxic evaluation of *Moringa oleifera* ethanolic leaf extract against A549 and MCF-12A cell lines. *Trends in Sciences, 19*. https://tis.wu.ac.th/index.php/tis/article/view/3202

Tukey, J. W. (1949). Comparing individual means in the analysis of variance. *Biometrics, 5*(2), 99-114.

Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., ... van Mulbregt, P. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods, 17*(3), 261-272.

Wang, S., Witek, J., Landrum, G. A., & Riniker, S. (2020). Improving conformer generation for small rings and macrocycles based on distance geometry and experimental torsional-angle preferences. *Journal of Chemical Information and Modeling, 60*(4), 2044-2058.

Xie, J., Luo, M., Chen, Q., Zhang, Q., Qin, L., Wang, Y., ... Zhao, Q. (2021). Hypolipidemic effect and mechanism of *Moringa oleifera* alkaloids, with cytotoxicity data against A549 cells. *Evidence-Based Complementary and Alternative Medicine, 2021*, 5591687.

Zdrazil, B., Felix, E., Hunter, F., Manners, E. J., Blackshaw, J., Corbett, S., ... Leach, A. R. (2024). The ChEMBL Database in 2023: A drug discovery platform spanning multiple bioactivity data types and time periods. *Nucleic Acids Research, 52*(D1), D1180-D1192.

**Note on the three potency anchors.** Tiloke et al. (2013), Xie et al. (2021), and the Trends in Sciences (2022) entry were recorded in the pipeline with their numeric IC₅₀, cell line, assay type, and URL rather than a complete citation string. Titles and author lists above were reconstructed from those records and should be checked against the source articles before submission. The numeric values themselves come straight from `results/tables/literature_prior.csv` and are traceable to the URLs stored there.
