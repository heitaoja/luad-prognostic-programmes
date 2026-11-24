# Statistical and measurement methods

The analysis specification is [analysis_plan_v2.md](analysis_plan_v2.md). The main manuscript includes the clinical interpretation and study limitations.

## Specimen and molecular features

The unit of inference is the patient. The primary TCGA sample rule takes the first eligible primary vial in barcode order. The reverse choice provides an alternate specimen sensitivity. The assay is an expression-based readout of the bulk tumour specimen.

The association analysis keeps the feature with the highest cohort mean log expression when several features map to one symbol. For model transport, aggregate duplicates by the within-specimen mean of log expression before ranking. Platform identifiers and mapping files define the shared gene universe. This step uses no outcomes, and duplicate aggregation uses no statistics from other participants. The expression scale, gene membership and per-set coverage are audited in outputs.

For a sample with N genes, let r_i be the descending expression midrank and v_i = N - r_i + 1. For a set of k measured genes, the score is the difference between the mean v_i of genes inside the set weighted by v_i^0.25 and the unweighted mean v_i outside the set, divided by N. For untied expression this is the normalized integral of the weighted rank walk. Tied genes receive the same contribution. In the association analysis, each component is z-scored within cohort and the equally weighted component average is z-scored again. Age is centered and scaled within cohort.

## Association and specificity

Cox models use overall survival time in months and a continuous programme score, with age, sex and stage as covariates. The complete-case count is documented by cohort. A secondary set of models assesses stage-stratified baseline hazards, stage I–II restriction, a 3-df spline and four-programme mutual adjustment. Schoenfeld tests describe model assumptions. Programme HRs from linear models describe a summary tendency when the spline suggests nonlinearity.

Synthesis uses REML to estimate between-cohort variance, inverse-variance weighting, and a modified Knapp–Hartung interval with variance multiplier bounded below by one and k–1 degrees of freedom. The normal-approximation DerSimonian–Laird synthesis is a sensitivity analysis. Holm adjustment applies to the four programme tests in each synthesis scope. The original external cohorts and the expanded external set have separate estimates.

Random programmes use a one-to-one mapping from genes in the observed component sets to other measured genes in the same average-rank decile. Reusing a replacement across component sets preserves overlap. The upper-tail empirical probability uses the plus-one correction. TCGA survival permutations shuffle time and status together and serve as a pipeline check. The alternative mean expression score averages gene-wise z-scores over the union of component genes.

## Transport of fitted models

The main transport gene universe is restricted to symbols available across TCGA and the four GEO series. TCGA defines the mean and standard deviation of each component, the composite reference, and the MKI67 percentile reference. Frozen features are combined with stage categories I–IV, age per decade centred at 65 and sex with female as the reference. TCGA alone determines model coefficients and baseline hazard. Test-cohort calibration slopes describe how the spread of the fixed predictions relates to survival.

Discrimination is measured by Harrell C. Paired percentile intervals use 1,000 patient resamples in each test cohort, with the fitted models fixed. For 36-month risk, observed mortality is estimated by Kaplan–Meier. The Brier score uses inverse censoring weights from the test cohort: G(t-) for observed deaths by the horizon and G(36) for patients followed beyond 36 months. Patients censored by the horizon carry zero weight in this score. Calibration plots group fixed predicted risks in quartiles and show Greenwood-based Kaplan–Meier intervals. The underlying C and Brier measures and score contrasts are available as numerical tables.

A platform match is an assay-design condition. Reproducing these fixed models in an additional cohort requires the archived gene list and the same preprocessing and reference values.

## RNA-seq evaluation

CPTAC estimates use cohort-standardised programme scores and Cox adjustment for stage I versus II–IV, age and sex. For transport, start from the primary gene list and remove the two symbols listed in the CPTAC assay audit. Repeat the TCGA fitting procedure with the 17,593 covered genes, and apply the resulting TCGA parameters to CPTAC. The sensitivity analysis uses the same five model formulae, risk horizon, metrics and bootstrap design as the GEO analyses. The RNA-seq coefficients and score references are available with the associated tables. [Audit](../results/tables/05b_rnaseq_audit.csv) · [Coefficients and 36-month baseline hazards](../results/tables/05b_rnaseq_frozen_coefficients.csv) · [Component and programme references](../results/tables/05b_rnaseq_score_reference.csv)
