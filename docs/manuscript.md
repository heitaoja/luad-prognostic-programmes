# Cross-platform reproducibility, redundancy and risk transport of Hallmark transcriptional programmes in lung adenocarcinoma

## Abstract

**Background.** We assessed the reproducibility, redundancy and predictive transport of transcriptional programmes associated with overall survival in lung adenocarcinoma.

**Methods.** We analysed 1,443 patients across six cohorts, using TCGA-LUAD for discovery and training. We scored four Hallmark programmes with tie-aware ssGSEA and pooled covariate-adjusted associations using two-stage REML with modified Knapp–Hartung inference (REML-mKH) and one-stage stratified individual participant data (IPD) Cox models. We assessed purity, smoking, driver status and scoring scale, benchmarked established signatures and evaluated frozen models using Harrell's C, Uno's inverse probability of censoring weighted (IPCW) C and 36-month baseline recalibration.

**Results.** Proliferation HRs from REML-mKH at $k=3,4,5$ were 1.32, 1.34 and 1.35 (five-cohort 95% CI 1.15–1.59, Holm-adjusted $P=0.027$). For proliferation and hypoxia, cohort-stratified IPD HRs were 1.34 and 1.36 at $k=4$ ($N=830$, adjusted $P=2.3\times10^{-6}$ and $1.1\times10^{-5}$), and 1.35 and 1.40 at $k=5$ ($N=935$, both adjusted $P=4.4\times10^{-7}$). Purity adjustment increased the hypoxia HR from 1.34 to 1.54. Proliferation estimates changed modestly after smoking and driver adjustment. Sample centring increased its TCGA mean gene-wise z-score HR from 1.08 to 1.24. The composite correlated strongly with CCP-31 and meta-PCNA ($\rho=0.927$–$0.965$), with no clear discrimination advantage over either rank-based signature. In CPTAC, the composite increased Uno's C by 0.071 relative to clinical predictors. In GSE31210, baseline recalibration reduced mean predicted mortality from 25.6% to 8.2% (observed 8.0%) and Brier error from 0.0935 to 0.0675.

**Conclusions.** Proliferation was the most reproducible programme, with substantial overlap among signatures. Absolute risk benefited from local baseline updating. Discrimination gains were generally modest, and recalibrated risks were evaluated in the same data used for updating.

## Introduction

Lung adenocarcinoma (LUAD) has substantial variation in survival within pathological stage. Tumour transcriptomes offer information on cellular state that may partly explain this variation and have long been evaluated as a complement to clinicopathological prognostication. The Director's Challenge compared expression-based survival models under blinded multisite validation and established a demanding standard for reproducibility across institutions. Kratz and colleagues subsequently translated an expression model into an assay applicable to fixed tissue and evaluated it internationally in resected non-squamous lung cancer. Wistuba and colleagues evaluated a 31-gene cell cycle progression (CCP) signature in early-stage LUAD, applying a cell-cycle measurement principle previously studied in prostate cancer by Cuzick and colleagues. These studies motivate a benchmark that examines the signal shared by different signatures while preserving the distinction between an expression score and a clinically validated assay. (1–4)

A recurrent difficulty in interpreting prognostic signatures is that many are correlated with the same biological process. Venet and colleagues showed in breast cancer that many arbitrary gene sets were associated with outcome and identified a core set of transcripts strongly correlated with proliferating cell nuclear antigen (meta-PCNA). Ringnér and Staaf found that a proliferation-related transcriptional programme underlay much of the agreement among published LUAD expression phenotypes and risk predictors. Adding genes or combining pathways could stabilise measurement of proliferation across platforms, yet an association with outcome can persist even when the resulting score contains little information beyond an individual marker or a simpler pathway. A comparison therefore needs both conditional survival models and paired external prediction contrasts against credible molecular comparators. (5, 6)

Bulk tumour expression also conflates variation in the malignant cells with variation in the proportions and states of surrounding cells. Systematic analyses using independent purity estimates have shown how composition influences transcriptional contrasts across cancers. DNA-derived estimates such as ABSOLUTE offer a complementary measure for assessing that dependence in the same specimen. Purity may affect immune and stromal signatures directly through changes in cell abundance and may influence hypoxia-related transcripts expressed in multiple compartments. Smoking and somatic alterations are additional sources of variation in LUAD and can associate with both tumour transcriptional state and outcome. Comparing coefficients before and after adjustment on identical patients provides one way to assess sensitivity to these measured factors. (7, 8)

Cross-platform analyses introduce a separate measurement problem. Expression matrices can differ in gene coverage, probe annotation, dynamic range and the role of sequencing depth in the expression scale. Scores that standardise each gene across patients retain sample-level offsets unless preprocessing addresses them. Within-sample ranks are insensitive to a common monotone transformation on a fixed gene set, making them an attractive basis for scoring data from different assays. Gene-dependent measurement differences and ties in low-count profiles can still affect the score. Comparing a rank-based score with expression averages before and after sample centring helps separate sensitivity to a global offset from changes in how the programme itself is summarised. (9)

A marker can have a reproducible hazard association while adding little to a clinical prediction model or producing poorly calibrated absolute risks in a new setting. REMARK addresses the transparent reporting of marker associations, whereas TRIPOD and TRIPOD+AI emphasise fit specification, external performance and reporting of prediction models. External validation should examine discrimination and calibration separately, with attention to censoring when comparing risk rankings and predicting at a fixed horizon. We accordingly designed a retrospective benchmark around four Hallmark programmes. We contrasted cohort-specific and pooled associations, tested sensitivity to biological and measurement factors, and evaluated how much externally transported clinical predictions gained from the proliferation composite and simpler or established signatures. (10–16)

## Methods

### Study design, cohorts and endpoint

We analysed primary LUAD tumour profiles from TCGA-LUAD, Gene Expression Omnibus (GEO) series GSE31210, GSE37745, GSE50081 and GSE72094, and CPTAC-LUAD. TCGA served as the discovery and model development cohort. GSE31210, GSE37745 and GSE50081 formed the original external scope. GSE72094 extended that scope to four microarray cohorts, and CPTAC added independent RNA-seq data to the five-cohort external scope. All external cases were held out from TCGA model fitting. The four programme choices followed exploratory analysis in TCGA. The original array findings were available during the revision, with the driver, IPD, recalibration and published-signature analyses documented as subsequent extensions. Cohort inclusion and expression input versions were audited against public source records and a pinned file manifest. (17–23) [Analysis specification](analysis_plan_v2.md) · [Data provenance](data_provenance.md)

We matched TCGA primary-tumour specimens to survival and clinical records using exact sample identifiers and then the 12-character participant barcode for unmatched records. Repeated vials were checked for agreement of survival information. The first eligible vial in barcode order represented each patient and the last eligible vial was used for a sensitivity analysis. For GSE31210 we applied the source exclusions for incomplete resection or adjuvant treatment. GSE37745 and GSE50081 were restricted to adenocarcinoma. GSE72094 patient identifiers and CPTAC case identifiers were checked within their sources, and CPTAC was linked across expression, clinical and survival files by case. The endpoint was all-cause death, with censoring at the last available follow-up. Days were divided by 30.4375 and years multiplied by 12. Patients with valid positive overall survival (OS) and known vital status entered the cohort audit. Complete-case models additionally required the specified covariates. We estimated median potential follow-up by reverse Kaplan–Meier, treating censoring times as events for that calculation. (24)

### Measurement and scoring

TCGA expression consisted of Xena-mirrored GDC STAR counts on the $\log_2(\mathrm{count}+1)$ scale. CPTAC used upper-quartile normalised, log-transformed RSEM estimates from LinkedOmics. GEO expression was analysed on the log scale with platform-specific probe annotation. We removed Ensembl version suffixes before symbol mapping and used pinned annotations for the GPL570 and HuRSTA arrays. In association analyses, duplicate features mapping to a single symbol were resolved using the feature with the highest cohort mean log expression. For model transport the duplicate log-expression features were averaged within each specimen before ranking, avoiding a target-cohort estimate in that part of preprocessing. Each platform's measured genes and gene-set coverage were audited. [Assay coverage](../results/tables/02_gene_coverage.csv) · [CPTAC coverage](../results/tables/05b_rnaseq_coverage.csv)

Four programmes combined 14 gene sets from the MSigDB Hallmark collection. Proliferation combined E2F targets, the G2/M checkpoint and MYC targets V1 and V2. Hypoxia combined hypoxia and glycolysis. EMT/stroma combined epithelial–mesenchymal transition, TGF-β signalling and angiogenesis. Immune combined inflammatory response, interferon-α and interferon-γ responses, allograft rejection and complement. For $G$ genes in a sample, we calculated descending expression midranks $r_g$ and used $v_g=G-r_g+1$. The single-sample gene set enrichment analysis (ssGSEA) statistic for a component was the difference between the mean $v_g$ within the set weighted by $v_g^{0.25}$ and its unweighted mean outside the set, divided by $G$. Midranks assigned identical contributions to tied genes. Within each cohort we standardised the component scores, averaged them with equal weights and standardised the average. This defines each association HR per cohort-specific SD for a predeclared programme. (25–28)

We also averaged gene-wise z-scores across the union of measured genes in each programme, excluding genes with zero variance. The original variant standardised each gene across specimens directly. The sample-centred variant first subtracted the median across all measured genes in each specimen, then applied the same gene-wise and final score standardisations. If a count can be approximated by a sample library factor times underlying abundance, then at sufficiently large counts the log profile contains an approximately additive term $\log_2(\mathrm{library\ factor})$. Gene-wise standardisation centres across patients, leaving the sample term in the transformed values. Subtracting the within-sample median removes any common additive term exactly on the log scale. This is a diagnostic for sample-level offsets, including possible depth variation. Offsets can also arise from composition or batch, and the additive approximation is weaker at low counts because of the pseudocount. The within-sample rank statistic has exact invariance to strictly increasing transformations of the matrix column used for scoring. (9)

### Clinical association and external synthesis

We fitted a separate Cox proportional hazards model for each programme, with OS time in months and adjustment for pathological stage (I–IV), cohort-standardised age and sex. The core CPTAC cohort-specific model used stage I versus II–IV to limit degrees of freedom with 23 deaths. We used four-stage coding for the one-stage pooled models and for fixed prediction. For the TCGA and GEO models we described proportional hazards (PH) using score-specific and global Schoenfeld tests, assessed the linear score term against a natural spline with three degrees of freedom and fitted sensitivities with a stage-stratified baseline or restriction to stages I–II. We jointly entered the four programmes to examine mutually adjusted associations. Exploratory within-stage TCGA estimates were adjusted for age and sex. Median-split Kaplan–Meier curves were used for visualisation only. (29)

In two-stage synthesis we first estimated clinical-adjusted log HRs and their standard errors by cohort. We then used inverse-variance random-effects weighting with restricted maximum likelihood (REML) estimates of between-cohort variance. Modified Knapp–Hartung inference bounded the variance multiplier below by one and used a $t$ reference distribution with $k-1$ degrees of freedom, where $k$ is the number of cohorts. Normal-approximation DerSimonian–Laird (DL) estimates and array leave-one-cohort-out summaries assessed sensitivity. One-stage IPD models pooled complete cases, used the within-cohort standardised programme and age values, and estimated a common programme coefficient from a Cox model with cohort-specific baseline hazards and fixed terms for stage, age and sex. A second IPD specification stratified the baseline by cohort and stage and included age and sex as covariates. Holm adjustment applied across the four programmes within each external scope and synthesis method. All external scopes excluded TCGA. (29–32) [Array synthesis](../results/tables/04_meta_analysis.csv) · [IPD and five-cohort synthesis](../results/tables/04b_ipd_and_pan_external_synthesis.csv)

### Purity, smoking, driver alterations and specificity

We linked DNA-derived ABSOLUTE purity from the Aran pan-cancer release to TCGA by exact specimen barcode. We examined Spearman correlations and paired Cox fits, preserving the same patients for the clinical-plus-score and clinical-plus-score-and-purity models. For smoking and driver sensitivity analyses, each pair likewise used the complete cases for the extended model. TCGA provided standardised pack-years. GSE31210 provided ever versus never smoking and a three-level alteration term (EGFR, KRAS or ALK, and EGFR/KRAS/ALK triple-negative). GSE50081 and GSE72094 provided ever versus never smoking. In a separate model GSE72094 additionally supported joint adjustment for recorded EGFR, KRAS, TP53 and STK11 mutation status. The GSE31210 and GSE72094 driver subgroups were exploratory, retaining the original cohort score scaling. For an event-free stage or sex category, we applied ridge penalisation to the stage and sex design columns ($\theta=1$, unscaled), leaving the programme and age terms unpenalised. (7, 8, 33) [Purity models](../results/tables/03c_purity_sensitivity.csv) · [Smoking and driver models](../results/tables/03d_driver_smoking_adjustment.csv)

The random-programme control asked whether the observed tumour programme produced a larger HR than gene sets with comparable structure. Each observed gene was replaced one-to-one by a measured gene outside the programme union, matched on its average-rank decile. The same replacement was assigned wherever the original gene recurred, preserving component sizes and overlap. We rescored and refitted proliferation and hypoxia in TCGA (500 draws per programme) and in GSE31210 and GSE72094 (200 each). The upper-tail empirical probability was $(1+b)/(1+B)$ for $b$ null HRs at least as large as observed and $B$ draws. As a separate pipeline control, we permuted TCGA OS times and event indicators together relative to the predictors 1,000 times per programme. A GSE31210 recurrence-endpoint analysis, winsorisation and stage restrictions assessed additional sensitivities.

### Molecular benchmark and frozen transport

The association benchmark compared the composite with log expression of MKI67 and the E2F and G2/M Hallmark scores. For each comparator we estimated clinical-adjusted marker HRs and the conditional HR for the composite on identical patients. We implemented CCP-31 from the Wistuba LUAD publication and meta-PCNA from Venet and colleagues' 131-gene list, resolving previous HGNC symbols with explicit aliases. We audited coverage and calculated Spearman correlations with the composite. We fitted clinical-adjusted and mutually conditional Cox models in both directions for the published signatures, analysing ranked scores and sample-centred mean gene-wise z-scores in parallel. We pooled these associations in the three external scopes using REML-mKH, DL and cohort-stratified IPD, reporting nominal inference for these secondary comparisons. Public expression implementations allow common-platform comparison of the gene sets. Reproduction of a commercial PCR assay would require its original assay procedures and decision thresholds. (3–5) [Signature definitions](../config/published_signatures.csv) · [Published benchmark](../results/tables/03e_published_signature_benchmark.csv)

We defined the assay universe for transport by intersecting annotated symbols across TCGA and the four GEO platforms, yielding 17,595 genes. For CPTAC, we removed two unavailable symbols and performed a separate TCGA fit on the 17,593 compatible genes before evaluating CPTAC. All patients in a transport comparison were scored on the corresponding fixed universe. TCGA defined the means and SDs for every component and the composite, the reference for the MKI67 within-sample percentile rank, Cox coefficients and the reference cumulative baseline hazard. The clinical model used stage I and female sex as reference categories and age in decades centred at 65 years. We added the composite, MKI67, E2F or G2/M one at a time. For CCP-31 and meta-PCNA we also fitted TCGA models with rank scores and, separately, a within-sample mean of measured log-expression values standardised to the TCGA distribution. All references and regression coefficients were held fixed during external evaluation.

We calculated Harrell's C and Uno's inverse probability of censoring weighted (IPCW) C from each frozen linear predictor. We used the survival concordance implementation with time weights proportional to $n(t)/\widehat G(t)^2$ for Uno's C over the available follow-up. Paired 95% percentile intervals for C differences used 1,000 resamples of the same test patients for all models in each comparison, with no refitting of the TCGA models. We compared each molecular model with the fixed clinical model and used paired Harrell differences to compare the composite directly with the molecular comparators. For 36-month assessment we required at least 10 deaths by the horizon and at least 10 patients at risk at the horizon. The observed risk was one minus Kaplan–Meier survival. The IPCW Brier score was the mean weighted squared prediction error, using inverse censoring survival just before an observed death by 36 months or at 36 months for those followed beyond it. Censoring on or before the horizon contributed zero weight. (14, 34, 35)

For baseline recalibration we fitted a Cox model in each eligible target cohort with the TCGA linear predictor $\eta_i$ as an offset with coefficient fixed at one. We evaluated the target cumulative baseline hazard $\widehat H_{0,\mathrm{recal}}(36)$ at $\eta=0$ and calculated individual mortality as $1-\exp[-\widehat H_{0,\mathrm{recal}}(36)\exp(\eta_i)]$. Frozen risks used the TCGA hazard in the same expression. We reused target censoring weights for the two Brier scores and identical risk groups for the frozen and recalibrated quartile plots of clinical and composite predictions. The plots show mean predicted risk and Kaplan–Meier observed risk with Greenwood-based intervals. A separate diagnostic estimated the calibration slope by freely regressing OS on $\eta_i$. Prediction used only the fixed slope. Recalibration performance describes an update estimated and evaluated on the same target outcome data. (13, 15, 16)

## Results

### 1. Cohort composition and technical coverage

The six cohorts contributed 1,443 eligible OS records with 476 deaths, and 1,420 patients with 469 deaths entered the complete-case association analyses (Table 1). In TCGA, 515 of 528 primary specimens had usable OS. Twelve participants had repeated vials among these records, leaving 503 distinct patients and 182 deaths after linkage and selection. In CPTAC, 105 of 110 source cases had usable OS and all 105 had complete clinical covariates. Median potential follow-up ranged from 22.8 months in CPTAC to 126.4 months in GSE37745. This variation, together with the low death count in GSE31210 relative to other array cohorts, motivated separate presentation of cohort estimates and risk calibration. [TCGA participant flow](../results/tables/01_tcga_participant_flow.csv) · [Cohort audit](../results/tables/01_cohort_audit.csv) · [CPTAC audit](../results/tables/05b_rnaseq_audit.csv)

**Table 1.** Eligible patients, complete cases and reverse Kaplan–Meier median follow-up.

| Cohort | Eligible OS | Deaths | Complete cases | Complete-case deaths | Median follow-up, months |
| --- | --- | --- | --- | --- | --- |
| TCGA-LUAD | 503 | 182 | 485 | 177 | 29.0 |
| GSE31210 | 204 | 30 | 204 | 30 | 61.0 |
| GSE37745 | 106 | 77 | 106 | 77 | 126.4 |
| GSE50081 | 127 | 51 | 127 | 51 | 65.9 |
| GSE72094 | 398 | 113 | 393 | 111 | 29.9 |
| CPTAC-LUAD | 105 | 23 | 105 | 23 | 22.8 |
| Total | 1443 | 476 | 1420 | 469 | — |

Every Hallmark set retained at least 96.5% measured coverage on both shared transport universes. CPTAC lacked ACTL10 and NPBWR1 from the array-compatible list. The cohort-standardised programme distributions differed in shape despite sharing a mean and scale by design (Figure 1). Switching to the last eligible TCGA vial yielded a proliferation HR of 1.30 (95% CI 1.11–1.53) compared with 1.31 (1.12–1.54) under the primary rule on the same 485 patients. [TCGA specimen sensitivity](../results/tables/03c_specimen_sensitivity.csv) · [CPTAC coverage](../results/tables/05b_rnaseq_coverage.csv)

### 2. Programme associations and cross-platform synthesis

In TCGA, clinical-adjusted HRs were 1.31 (95% CI 1.12–1.54) for proliferation and 1.42 (1.22–1.66) for hypoxia. Both estimates exceeded one in every external cohort. In CPTAC they were 1.57 (1.00–2.48) and 2.02 (1.24–3.27), respectively, with uncertainty reflecting 23 deaths. The fitted CPTAC lower interval limit for proliferation was 1.002 before rounding to two decimals. EMT/stroma and immune estimates varied in direction. The TCGA median-split curves and exploratory stage-specific results are shown in Figure 2. Proliferation estimates for stages I, II and III were 1.29, 1.10 and 1.68, respectively, with a stage interaction likelihood ratio $P=0.310$. The global PH test for the proliferation model in TCGA gave $P=0.039$ and the score-specific test gave $P=0.081$. Stratifying the baseline by stage yielded an HR of 1.31. [Cohort estimates](../results/tables/04_cohort_estimates.csv) · [CPTAC estimates](../results/tables/05b_rnaseq_association.csv) · [Diagnostics](../results/tables/03c_model_diagnostics.csv) · [Stage interaction](../results/tables/03_stage_interaction.csv)

Proliferation HRs were similar across the external synthesis scopes (Table 2), including the original three-cohort REML-mKH estimate of 1.32 (95% CI 0.93–1.86) and the expanded four-array estimate of 1.34 (1.10–1.62). Holm-adjusted $P$ values were 0.297 and 0.069, respectively. Including CPTAC yielded a five-cohort REML-mKH HR of 1.35 (1.15–1.59), with Holm-adjusted $P=0.027$ and estimated $I^2=0\%$. Hypoxia had a five-cohort REML-mKH HR of 1.42 (1.04–1.94), with Holm-adjusted $P=0.103$ and $I^2=57.6\%$. The one-stage cohort-stratified estimates supported associations for both scores, with similar coefficients after stratifying by cohort and stage. In exploratory models containing all four scores, the four-array pooled REML-mKH HRs were 1.23 (0.96–1.57) for proliferation and 1.25 (0.94–1.65) for hypoxia. Neither EMT/stroma nor immune showed a reproducible pooled association. [Original meta-analysis](../results/tables/04_meta_analysis.csv) · [IPD and five-cohort synthesis](../results/tables/04b_ipd_and_pan_external_synthesis.csv) · [Joint models](../results/tables/03c_marker_synthesis.csv)

**Table 2.** External HRs per SD (95% CI) by evidence scope and model. Holm-adjusted $P$ values appear in brackets.

| Cohorts (patients/deaths) | Synthesis | Proliferation | Hypoxia | EMT/stroma | Immune |
| --- | --- | --- | --- | --- | --- |
| 3 (437/158) | Two-stage REML-mKH | 1.32 (0.93–1.86) [0.297] | 1.37 (0.66–2.82) [0.613] | 1.02 (0.52–2.00) [1] | 1.09 (0.57–2.08) [1] |
| 3 (437/158) | One-stage IPD, cohort | 1.34 (1.15–1.56) [0.000903] | 1.33 (1.12–1.58) [0.00379] | 0.99 (0.83–1.18) [1] | 1.05 (0.90–1.24) [1] |
| 3 (437/158) | One-stage IPD, cohort and stage | 1.31 (1.12–1.53) [0.00351] | 1.31 (1.10–1.56) [0.00689] | 0.97 (0.82–1.16) [1] | 1.05 (0.89–1.22) [1] |
| 4 (830/269) | Two-stage REML-mKH | 1.34 (1.10–1.62) [0.0688] | 1.35 (0.95–1.90) [0.212] | 1.03 (0.74–1.43) [1] | 1.03 (0.73–1.46) [1] |
| 4 (830/269) | One-stage IPD, cohort | 1.34 (1.20–1.51) [2.3e-06] | 1.36 (1.19–1.54) [1.09e-05] | 1.02 (0.90–1.16) [1] | 0.99 (0.88–1.12) [1] |
| 4 (830/269) | One-stage IPD, cohort and stage | 1.32 (1.18–1.49) [1.35e-05] | 1.34 (1.18–1.52) [3.03e-05] | 1.02 (0.90–1.15) [1] | 1.00 (0.89–1.13) [1] |
| 5 (935/292) | Two-stage REML-mKH | 1.35 (1.15–1.59) [0.0271] | 1.42 (1.04–1.94) [0.103] | 1.09 (0.81–1.46) [0.97] | 1.02 (0.79–1.31) [0.97] |
| 5 (935/292) | Two-stage DL | 1.35 (1.20–1.51) [1.09e-06] | 1.42 (1.16–1.73) [0.00174] | 1.08 (0.88–1.33) [0.883] | 1.02 (0.86–1.21) [0.883] |
| 5 (935/292) | One-stage IPD, cohort | 1.35 (1.21–1.51) [4.4e-07] | 1.40 (1.23–1.58) [4.4e-07] | 1.05 (0.93–1.19) [0.839] | 1.00 (0.89–1.12) [0.993] |
| 5 (935/292) | One-stage IPD, cohort and stage | 1.34 (1.19–1.50) [2.14e-06] | 1.38 (1.22–1.56) [2e-06] | 1.05 (0.93–1.19) [0.83] | 1.00 (0.89–1.13) [0.955] |

### 3. Specificity and sensitivity to the expression scale

Against the structure-matched null in TCGA, the upper-tail empirical $P$ value was 0.0020 for both proliferation and hypoxia. In GSE72094 both were 0.0050. In GSE31210 the corresponding values were 0.0945 and 0.0199. Thus the comparatively small Japanese cohort gave more specific evidence for hypoxia than for proliferation under the random-programme comparison, despite a positive proliferation HR. The paired-outcome permutation in TCGA gave an empirical $P$ of 0.0010 for both scores. Using the exploratory GSE31210 recurrence annotation, with 54 events, gave HRs of 1.68 (95% CI 1.32–2.13) for proliferation and 1.83 (1.37–2.43) for hypoxia. [Structure-matched draws](../results/tables/06_random_set_null.csv) · [Outcome permutation](../results/tables/06_permutation_test.csv) · [Sensitivity analysis](../results/tables/06_sensitivity.csv)

Scoring scale made a substantial difference in TCGA. The untreated gene-wise z-score mean gave a proliferation HR of 1.08 (95% CI 0.92–1.25), rising to 1.24 (1.06–1.45) with sample centring and 1.31 (1.12–1.54) with ssGSEA. Hypoxia gave HRs of 1.03 (0.89–1.20), 1.21 (1.04–1.41) and 1.42 (1.22–1.66), respectively. On the GSE31210 array, the proliferation HR changed from 1.43 (1.02–2.00) to 1.33 (0.96–1.84) after centring. The gain in TCGA was therefore specific to the measurement setting examined and the choice of log-expression scale. Both scoring variants are shown alongside the other model sensitivities in Figure 3. [Scoring sensitivity](../results/tables/06_sensitivity.csv)

### 4. Purity, smoking and driver profiles

Within the 323 TCGA cases with exact-specimen ABSOLUTE estimates and 120 deaths, proliferation had a weak correlation with purity ($\rho=+0.05$), and its adjusted HR changed from 1.33 to 1.34 (Table 3). Hypoxia had a negative correlation ($\rho=-0.29$) and an HR increase from 1.34 to 1.54. Immune had the largest absolute correlation ($\rho=-0.70$) and its HR increased from 0.98 to 1.26 with a 95% CI of 0.96–1.64. EMT/stroma also shifted upwards with purity adjustment. These within-sample comparisons separate the observed coefficient changes from the loss of patients without a linked DNA measurement. [Purity sensitivity](../results/tables/03c_purity_sensitivity.csv)

**Table 3.** Clinical-adjusted HRs (95% CI) before and after ABSOLUTE purity adjustment on matched TCGA specimens.

| Programme | Spearman ρ | Clinical | Clinical + purity |
| --- | --- | --- | --- |
| Proliferation | 0.05 | 1.33 (1.10–1.62) | 1.34 (1.10–1.62) |
| Hypoxia | −0.29 | 1.34 (1.11–1.61) | 1.54 (1.27–1.88) |
| EMT/stroma | −0.50 | 1.05 (0.88–1.27) | 1.22 (0.97–1.52) |
| Immune | −0.70 | 0.98 (0.82–1.17) | 1.26 (0.96–1.64) |

In TCGA the proliferation and hypoxia estimates changed minimally after pack-years adjustment among 333 patients. In GSE31210, adding smoking and the grouped alteration variable moved the proliferation HR from 1.42 to 1.38 (95% CI 0.98–1.95) and the hypoxia HR from 1.96 to 2.03 (1.37–3.00). Separately adjusting for smoking or the four driver indicators in GSE72094 changed both estimates modestly (Table 4). Among the 116 EGFR-mutant GSE31210 cases with 12 deaths, the proliferation estimate was 2.10 (1.10–4.03). Among the 41 EGFR-mutant GSE72094 cases with four deaths, the ridge-assisted estimate was 1.22 (0.36–4.14). STK11-mutant GSE72094 estimates were 0.84 (0.48–1.47) and 0.89 (0.52–1.51) for proliferation and hypoxia versus 1.43 and 1.41 in wild-type cases. These subgroup intervals show the limited precision available to investigate effect modification. [Driver and smoking analysis](../results/tables/03d_driver_smoking_adjustment.csv)

**Table 4.** Matched complete-case comparison of programme HRs with clinical adjustment and with additional smoking or driver adjustment.

| Cohort | Added covariates | Patients | Deaths | Proliferation (95% CI) | Hypoxia (95% CI) |
| --- | --- | --- | --- | --- | --- |
| TCGA-LUAD | Pack-years | 333 | 118 | 1.33 (1.09–1.62) → 1.34 (1.10–1.63) | 1.37 (1.14–1.65) → 1.37 (1.14–1.65) |
| GSE31210 | Smoking and alteration group | 204 | 30 | 1.42 (1.01–1.99) → 1.38 (0.98–1.95) | 1.96 (1.33–2.90) → 2.03 (1.37–3.00) |
| GSE50081 | Smoking | 115 | 46 | 1.20 (0.92–1.57) → 1.19 (0.90–1.57) | 1.25 (0.91–1.72) → 1.25 (0.91–1.72) |
| GSE72094 | Smoking | 328 | 89 | 1.46 (1.19–1.78) → 1.46 (1.19–1.79) | 1.39 (1.12–1.73) → 1.39 (1.12–1.73) |
| GSE72094 | EGFR, KRAS, TP53 and STK11 | 393 | 111 | 1.36 (1.13–1.62) → 1.38 (1.13–1.68) | 1.36 (1.12–1.65) → 1.30 (1.07–1.58) |

### 5. Redundancy with simple markers and published signatures

Spearman correlations with the composite across TCGA and the four GEO cohorts ranged from 0.69 to 0.85 for MKI67, from 0.97 to 0.98 for E2F and from 0.94 to 0.97 for G2/M. The simple proliferation measures all had positive clinical-adjusted associations in the four-array REML-mKH synthesis. The composite HR conditional on G2/M was 0.89 (95% CI 0.37–2.11), and its estimates conditional on E2F and MKI67 were also imprecise (Table 5). Across all six cohorts the rank-based CCP-31 and meta-PCNA scores correlated with the composite at 0.93–0.96 and 0.94–0.96, respectively, to two decimals. Each association platform measured all 31 CCP genes and 130 or 131 meta-PCNA genes. The fixed transport spaces retained 31 and 129, respectively. [Simple marker correlations](../results/tables/03c_marker_correlations.csv) · [Published benchmark](../results/tables/03e_published_signature_benchmark.csv) · [Transport coverage](../results/tables/05d_published_signature_transport.csv)

**Table 5.** Clinical-adjusted and conditional molecular associations. The four-array REML-mKH comparison uses the same cohorts for all scores. The published rank-based scores also have five-cohort estimates.

| Score | Cohorts | Clinical + score | Composite given score | Score given composite |
| --- | --- | --- | --- | --- |
| Composite | 4 | 1.34 (1.10–1.62) | — | — |
| MKI67 | 4 | 1.28 (1.04–1.56) | 1.30 (0.96–1.76) | — |
| E2F | 4 | 1.32 (1.09–1.60) | 1.26 (0.44–3.57) | — |
| G2/M | 4 | 1.36 (1.12–1.65) | 0.89 (0.37–2.11) | — |
| CCP-31 | 4 | 1.42 (1.14–1.77) | 0.65 (0.16–2.66) | 2.31 (0.43–12.41) |
| meta-PCNA | 4 | 1.37 (1.12–1.67) | 0.76 (0.14–4.10) | 1.89 (0.29–12.39) |
| CCP-31 | 5 | 1.45 (1.20–1.74) | 0.60 (0.21–1.74) | 2.61 (0.73–9.37) |
| meta-PCNA | 5 | 1.38 (1.17–1.64) | 0.78 (0.23–2.57) | 1.89 (0.50–7.17) |

Under five-cohort IPD, the clinical-adjusted HRs were 1.45 (95% CI 1.29–1.64) for CCP-31 and 1.39 (1.24–1.56) for meta-PCNA. The composite HR was 0.75 (0.54–1.06) conditional on CCP-31 and 0.96 (0.65–1.40) conditional on meta-PCNA. The reciprocal IPD estimate for CCP-31 conditional on the composite was 1.91 (1.34–2.73), whereas the corresponding REML-mKH interval spanned one. The reciprocal result is therefore sensitive to inference and model assumptions. In all ten external head-to-head contrasts between the composite and a rank-based published score, the paired Harrell difference interval included zero. One of ten exploratory mean-log-expression contrasts favoured the composite at the nominal 95% level, against meta-PCNA in GSE72094 (difference 0.020, 95% CI 0.005–0.036). The direct intervals versus the simple markers also spanned zero. Figure 4 relates these conditional and prediction contrasts to comparisons with clinical adjustment alone. [Published associations](../results/tables/03e_published_signature_benchmark.csv) · [GEO marker contrasts](../results/tables/05_simple_marker_contrasts.csv) · [CPTAC marker contrasts](../results/tables/05b_rnaseq_marker_contrasts.csv) · [Published transport](../results/tables/05d_published_signature_transport.csv)

### 6. Frozen discrimination, absolute risk and baseline recalibration

Adding the composite to the frozen clinical model gave positive point differences in Harrell's and Uno's C in all five external cohorts (Table 6). Both paired intervals excluded zero only in GSE72094. In CPTAC, Harrell's C rose by 0.010 (95% CI −0.071 to 0.104) and Uno's C by 0.071 (−0.023 to 0.192). The difference between concordance estimators was largest there, where follow-up and the number of deaths were limited. Among the arrays, Uno increments ranged from 0.016 in GSE37745 to 0.068 in GSE31210. These external estimates reflect fixed TCGA predictions and their sampling variability within each target cohort. [GEO transport](../results/tables/05_incremental_value.csv) · [CPTAC transport](../results/tables/05b_rnaseq_transport.csv)

**Table 6.** External discrimination of the frozen clinical and clinical-plus-composite models.

| External cohort | Harrell C, clinical | Harrell C, composite | Harrell ΔC (95% CI) | Uno C, clinical | Uno C, composite | Uno ΔC (95% CI) |
| --- | --- | --- | --- | --- | --- | --- |
| GSE31210 | 0.688 | 0.743 | 0.054 (−0.012 to 0.119) | 0.676 | 0.743 | 0.068 (−0.042 to 0.176) |
| GSE37745 | 0.601 | 0.618 | 0.017 (−0.038 to 0.069) | 0.607 | 0.623 | 0.016 (−0.037 to 0.065) |
| GSE50081 | 0.628 | 0.664 | 0.036 (−0.007 to 0.084) | 0.625 | 0.683 | 0.059 (−0.025 to 0.148) |
| GSE72094 | 0.640 | 0.674 | 0.034 (0.008 to 0.062) | 0.601 | 0.653 | 0.051 (0.006 to 0.102) |
| CPTAC-LUAD | 0.740 | 0.750 | 0.010 (−0.071 to 0.104) | 0.694 | 0.765 | 0.071 (−0.023 to 0.192) |

The largest frozen overprediction occurred in GSE31210. Its observed 36-month mortality was 8.0%, compared with 28.6% from the frozen clinical model and 25.6% from the frozen composite model. Baseline recalibration brought both to 8.2% and reduced composite IPCW Brier error from 0.0935 to 0.0675 (Table 7). The composite Brier score also fell in the remaining four cohorts, with a small change in GSE72094 where mean risk was already close to the observed rate. In CPTAC the clinical model's recalibrated mean was closer to the observed rate while its Brier score increased slightly. Calibration across the risk distribution was evaluated separately in Figure 5. In GSE31210 the frozen composite calibration slope was 1.73, indicating that a corrected mean can coexist with differences in the spread of relative risks. [GEO recalibration curves](../results/tables/05c_recalibration.csv) · [CPTAC recalibration curves](../results/tables/05b_rnaseq_recalibration.csv)

**Table 7.** Mean 36-month mortality predictions, Kaplan–Meier observed mortality and IPCW Brier scores for the clinical and composite models.

| External cohort | Model | Observed KM | Frozen mean risk | Recalibrated mean risk | Frozen Brier | Recalibrated Brier |
| --- | --- | --- | --- | --- | --- | --- |
| GSE31210 | Clinical | 8.0% | 28.6% | 8.2% | 0.1047 | 0.0683 |
| GSE31210 | Clinical + composite | 8.0% | 25.6% | 8.2% | 0.0935 | 0.0675 |
| GSE37745 | Clinical | 39.6% | 34.6% | 38.5% | 0.2393 | 0.2379 |
| GSE37745 | Clinical + composite | 39.6% | 37.9% | 38.5% | 0.2309 | 0.2308 |
| GSE50081 | Clinical | 24.5% | 32.3% | 24.6% | 0.1802 | 0.1748 |
| GSE50081 | Clinical + composite | 24.5% | 30.8% | 24.8% | 0.1736 | 0.1711 |
| GSE72094 | Clinical | 36.6% | 36.8% | 36.5% | 0.2108 | 0.2108 |
| GSE72094 | Clinical + composite | 36.6% | 37.2% | 36.6% | 0.2064 | 0.2063 |
| CPTAC-LUAD | Clinical | 34.4% | 37.6% | 33.8% | 0.1796 | 0.1798 |
| CPTAC-LUAD | Clinical + composite | 34.4% | 31.8% | 36.0% | 0.1722 | 0.1679 |

## Discussion

Proliferation showed the strongest cross-platform reproducibility of the four programmes. The movement in its REML-mKH inference between three, four and five external cohorts was accompanied by little movement in the pooled HR. With few cohorts, the modified procedure retains a variance safeguard and uses a $t$ distribution with two or three degrees of freedom in the array scopes. This produces wider intervals than a normal approximation even when the estimated heterogeneity is zero. The five-cohort result increased precision while extending the evidence across assay types. The IPD analyses estimated a common score slope with cohort-specific baseline risks, using individual risk sets and harmonised covariates. Their more precise proliferation and hypoxia results answer a complementary question about a shared association. The patients also contribute to the two-stage synthesis, and the IPD estimates rely on a common slope even when cohort-specific hypoxia effects vary materially. The larger hypoxia heterogeneity and its remaining multiplicity-adjusted uncertainty in the two-stage analysis matter when extrapolating to another cohort. (29–32)

The purity analysis suggests that tissue composition can mask prognostic associations in bulk expression. The proliferation score showed little monotonic relation to ABSOLUTE purity ($\rho=+0.05$), and its hazard association barely changed with adjustment. This is consistent with a largely tumour-intrinsic proliferative component in this dataset, although proliferating non-malignant cells can contribute to expression. Hypoxia and glycolysis combine metabolic and stress responses from both malignant and non-malignant cells. Their negative correlation with purity and the shift from HR 1.34 to 1.54 suggest that comparing tumours with similar cellularity exposes a stronger adverse association. Immune expression showed the strongest relation to purity ($\rho=-0.70$), consistent with immune-cell abundance contributing to the bulk signal. Its HR rose from 0.98 to 1.26, but the adjusted interval remained compatible with a small or absent association. ABSOLUTE is a proxy for the malignant fraction. Resolving individual stromal and immune populations or the causal role of purity would require additional measurements and designs. (7, 8)

The smoking and driver sensitivities address differences in clinical context between the Japanese GSE31210 series, enriched for EGFR-mutant disease and never smokers, and the US GSE72094 series, with many ever smokers and a reported spectrum of KRAS, TP53 and STK11 alterations. The TCGA pack-years analysis likewise retained the proliferation association among patients with documented exposure. The modest paired coefficient changes support a contribution from proliferation beyond the measured smoking and driver terms in the available samples. GSE31210's adjusted interval spanning one and GSE72094's smaller hypoxia estimate with driver adjustment qualify the strength of that inference. The discrepant point estimates in the EGFR-mutant strata and the attenuated STK11-mutant estimates could reflect biological differences, small numbers of events, or a combination of both. In particular, the four deaths in the GSE72094 EGFR-mutant stratum required stabilisation of nuisance terms. Reliable estimation of interactions requires more events. The available clinical data also leave room for differences in treatment, exposure measurement and comorbidity. (17, 19, 22)

The scoring sensitivity explains how a change in preprocessing can obscure an otherwise reproducible association. If two log-count profiles differ by a common additive shift, the gene-wise z-score calculation sends that shift into each programme gene's standardised value with an inverse-SD weight. The mean across the set can then largely measure a shared sample offset, particularly when sequencing counts have received no library-size normalisation upstream of the log transform. In TCGA, proliferation attenuated to HR 1.08 under the uncentred expression mean. Removing the column median raised it to 1.24, partway towards the rank-based HR of 1.31. A column shift cannot alter the order of genes, which explains this aspect of the robustness of ssGSEA. Quantifying the contribution of sequencing depth would require a direct analysis of library size and other sample-level factors. The median also removes biological global shifts, and gene-specific assay effects can alter ranks. Variation across platforms supports publishing the exact scoring universe, transforms and reference parameters with the results. (9, 26–28)

The comparison with established markers favours parsimony when the clinical objective is risk ranking. The composite correlated strongly with E2F, G2/M, CCP-31 and meta-PCNA, with correlations reaching 0.98. MKI67 was less correlated and is a distinct, simpler measurement. Strong correlation makes it difficult to attribute prognostic information to one component in a conditional model. The wide conditional intervals, the reciprocal CCP-31 IPD association and the paired external discrimination results collectively suggest shared information with some model-dependent differences. The data give little evidence that combining the four Hallmarks produces a reproducible improvement over the ranked published signatures or a simple pathway score. Equivalence would require a prespecified margin and adequate precision. The performance of an implemented CCP assay would depend on its analytical platform, reference genes and clinical use population. (3–6)

Discrimination can transport across assay platforms when the ordering induced by the frozen predictor is preserved sufficiently well, while absolute risk depends on the target baseline hazard and covariate mix. An offset update preserves the ranking and therefore leaves both concordance indices unchanged, yet can materially change absolute prediction error, as shown by GSE31210. The modest and uncertain C increments elsewhere, the early censoring and event counts in CPTAC, and the variation in cohort calibration slopes constrain the predictive interpretation. The results also depend on available stage, age and sex data, inconsistent historical follow-up origins and absence of uniform treatment information. The cohort-specific censoring weights assume that the marginal censoring estimate is adequate for evaluation. Recalibrated probabilities were evaluated using the same outcome data used for the baseline update, and would require independent evaluation before clinical use. These findings support reporting discrimination together with calibration and comparing molecular complexity on the basis of performance measured in the target population. (10–16, 34, 35)

## Data and code availability

The analysis code, versioned input manifest, public source accessions, gene definitions and model references are provided with this study. The pipeline regenerates the summary tables and figures from the pinned inputs. Patient-level derived expression scores are produced locally from those sources. [Input manifest](../data/checksums.csv) · [Methods](methods.md) · [Reporting crosswalk](reporting_crosswalk.md) · [Results workbook](../results/tables/project_tables.xlsx)

## Figure Legends

### Figure 1. Cohort composition and score distributions

**A,** Bars report OS-eligible and clinical complete-case counts for six cohorts. Points display reverse Kaplan–Meier median follow-up in months. **B,** Density curves compare within-cohort standardised scores for proliferation, hypoxia, EMT/stroma and immune across the five external cohorts and TCGA-LUAD. Each colour denotes a cohort. [Figure 1 PDF](../results/figures/Figure1_Cohort_and_Scores.pdf)

### Figure 2. Programme associations and external synthesis

**A,** Points and 95% confidence intervals show clinical-adjusted HRs per within-cohort SD in six cohorts, two-stage REML-mKH summaries for three, four and five external cohorts, and cohort-stratified one-stage IPD Cox summaries for four and five. The individual CPTAC association uses stage I versus II–IV. The vertical line indicates HR 1. **B,** TCGA Kaplan–Meier curves compare patients above versus at or below each programme median for proliferation and hypoxia. Crosses indicate censoring. **C,** Exploratory TCGA HRs by stages I, II and III adjust for age and sex. [Figure 2 PDF](../results/figures/Figure2_Association_and_IPD_Synthesis.pdf)

### Figure 3. Specificity, sample composition and sensitivity analyses

**A,** Histograms show Cox HRs for random programmes matched on component size, overlap and gene rank decile. Vertical lines show the observed proliferation or hypoxia HR in TCGA, GSE31210 and GSE72094. **B,** Paired TCGA estimates compare clinical adjustment with clinical and ABSOLUTE adjustment on the same specimens. Labels give Spearman correlations between score and purity. **C,** Forest plots show matched clinical and extended smoking or driver models and exploratory proliferation and hypoxia estimates within driver strata. The sparse EGFR-mutant GSE72094 fit penalises stage and sex. **D,** The sensitivity heatmap shows HRs under alternative endpoints, covariates and scoring methods, including both mean gene-wise z-score variants. [Figure 3 PDF](../results/figures/Figure3_Specificity_Purity_and_Drivers.pdf)

### Figure 4. Redundancy and external discrimination

**A,** Four-array REML-mKH estimates show clinical-adjusted associations for the composite, MKI67, E2F, G2/M, CCP-31 and meta-PCNA. Conditional rows show the composite given each comparator and the published scores given the composite. **B,** Harrell and Uno IPCW concordance differences are relative to the frozen clinical model in the four GEO cohorts and CPTAC. Published signatures use ssGSEA. Intervals are paired patient bootstrap 95% intervals. [Figure 4 PDF](../results/figures/Figure4_Redundancy_Benchmark.pdf)

### Figure 5. Frozen and baseline-recalibrated 36-month mortality

**A,** In each of the five test cohorts, points show mean predicted and Kaplan–Meier observed mortality in identical risk groups for the clinical and composite models. Intervals describe uncertainty in observed mortality. The diagonal marks perfect agreement. **B,** Paired points compare overall mean predicted mortality, IPCW Brier scores and their reductions after updating the baseline hazard. The observed mortality marker is common to both models in a cohort. Positive reductions indicate lower Brier error. [Figure 5 PDF](../results/figures/Figure5_Calibration_and_Recalibration.pdf)

## Supplementary Figure Legends

### Figure S1. External RNA-seq evaluation
**A,** Clinical-adjusted programme HRs in CPTAC. **B,** Paired bootstrap Harrell C contrasts in CPTAC for molecular models trained on the compatible TCGA gene universe. [figS1](../results/figures/figS1_rnaseq_validation.pdf)

### Figure S2. Cohort size and follow-up

Eligible and complete-case counts and reverse Kaplan–Meier follow-up for TCGA and the four array cohorts. [fig1](../results/figures/fig1_cohort_overview.pdf)

### Figure S3. Distributions of programme scores

Densities for the four cohort-standardised programmes in TCGA and the four array cohorts. [fig2](../results/figures/fig2_axis_distribution.pdf)

### Figure S4. Individual and pooled associations in the array analysis

Forest plots show TCGA and GEO HRs and the three- and four-cohort external array syntheses with 95% confidence intervals. [fig3](../results/figures/fig3_forest_all_cohorts.pdf)

### Figure S5. TCGA survival by score median

Kaplan–Meier curves contrast above-median and at-or-below-median groups for visualisation. [fig4](../results/figures/fig4_km_tcga.pdf)

### Figure S6. TCGA estimates by pathological stage

Within-stage estimates and 95% confidence intervals compare exploratory programme associations adjusted for age and sex. [fig5](../results/figures/fig5_stage_stratified.pdf)

### Figure S7. Incremental external Harrell concordance

Test-cohort differences from the frozen clinical model are shown with paired 95% patient bootstrap intervals. [fig6](../results/figures/fig6_incremental.pdf)

### Figure S8. Structure-matched random programmes

Histograms contrast observed programme HRs with null HRs generated from measured genes of similar rank abundance. [fig7](../results/figures/fig7_negative_control.pdf)

### Figure S9. Endpoint, scoring and adjustment sensitivities

Cells summarise HRs obtained under alternative scoring, clinical adjustment, stage restriction and recurrence analyses. [fig8](../results/figures/fig8_sensitivity.pdf)

### Figure S10. Comparison with simpler proliferation measures

Forest plots compare clinical-adjusted and conditional HRs for the composite and simpler markers across array cohorts. [fig9](../results/figures/fig9_marker_benchmark.pdf)

### Figure S11. Frozen risk calibration

Quartile means of frozen 36-month risks are plotted against Kaplan–Meier observed mortality in each array cohort. [fig10](../results/figures/fig10_calibration.pdf)

### Figure S12. Purity sensitivity

Matched TCGA estimates show clinical and purity-adjusted programme and marker HRs with 95% confidence intervals. [fig11](../results/figures/fig11_purity_sensitivity.pdf)

## References

1. Shedden K, Taylor JM, Enkemann SA, et al. Gene expression-based survival prediction in lung adenocarcinoma: a multi-site, blinded validation study. *Nat Med.* 2008, 14, 822–827. [DOI](https://doi.org/10.1038/nm.1790)
2. Kratz JR, He J, Van Den Eeden SK, et al. A practical molecular assay to predict survival in resected non-squamous, non-small-cell lung cancer: development and international validation studies. *Lancet.* 2012, 379, 823–832. [PMID](https://pubmed.ncbi.nlm.nih.gov/22285053/)
3. Wistuba II, Behrens C, Lombardi F, et al. Validation of a proliferation-based expression signature as prognostic marker in early stage lung adenocarcinoma. *Clin Cancer Res.* 2013, 19, 6261–6271. [DOI](https://doi.org/10.1158/1078-0432.CCR-13-0596)
4. Cuzick J, Swanson GP, Fisher G, et al. Prognostic value of an RNA expression signature derived from cell cycle proliferation genes in patients with prostate cancer: a retrospective study. *Lancet Oncol.* 2011, 12, 245–255. [DOI](https://doi.org/10.1016/S1470-2045(10)70295-3)
5. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol.* 2011, 7, e1002240. [DOI](https://doi.org/10.1371/journal.pcbi.1002240)
6. Ringnér M, Staaf J. Consensus of gene expression phenotypes and prognostic risk predictors in primary lung adenocarcinoma. *Oncotarget.* 2016, 7, 52957–52973. [DOI](https://doi.org/10.18632/oncotarget.10641)
7. Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nat Commun.* 2015, 6, 8971. [DOI](https://doi.org/10.1038/ncomms9971)
8. Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA alterations in human cancer. *Nat Biotechnol.* 2012, 30, 413–421. [DOI](https://doi.org/10.1038/nbt.2203)
9. Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol.* 2010, 11, R25. [DOI](https://doi.org/10.1186/gb-2010-11-3-r25)
10. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). *J Natl Cancer Inst.* 2005, 97, 1180–1184. [DOI](https://doi.org/10.1093/jnci/dji237)
11. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. *Ann Intern Med.* 2015, 162, 55–63. [DOI](https://doi.org/10.7326/M14-0697)
12. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ.* 2024, 385, e078378. [DOI](https://doi.org/10.1136/bmj-2023-078378)
13. Riley RD, Archer L, Snell KIE, et al. Evaluation of clinical prediction models (part 2): how to undertake an external validation study. *BMJ.* 2024, 384, e074820. [DOI](https://doi.org/10.1136/bmj-2023-074820)
14. Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med.* 2011, 30, 1105–1117. [DOI](https://doi.org/10.1002/sim.4154)
15. Steyerberg EW. *Clinical Prediction Models: A Practical Approach to Development, Validation, and Updating*. 2nd ed. Springer, 2019. [DOI](https://doi.org/10.1007/978-3-030-16399-0)
16. Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW. Calibration: the Achilles heel of predictive analytics. *BMC Med.* 2019, 17, 230. [DOI](https://doi.org/10.1186/s12916-019-1466-7)
17. The Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. *Nature.* 2014, 511, 543–550. [DOI](https://doi.org/10.1038/nature13385)
18. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell.* 2018, 173, 400–416.e11. [DOI](https://doi.org/10.1016/j.cell.2018.02.052)
19. Okayama H, Kohno T, Ishii Y, et al. Identification of genes upregulated in ALK-positive and EGFR/KRAS/ALK-negative lung adenocarcinomas. *Cancer Res.* 2012, 72, 100–111. [DOI](https://doi.org/10.1158/0008-5472.CAN-11-1403)
20. Botling J, Edlund K, Lohr M, et al. Biomarker discovery in non-small cell lung cancer: integrating gene expression profiling, meta-analysis, and tissue microarray validation. *Clin Cancer Res.* 2013, 19, 194–204. [DOI](https://doi.org/10.1158/1078-0432.CCR-12-1139)
21. Der SD, Sykes J, Pintilie M, et al. Validation of a histology-independent prognostic gene signature for early-stage, non-small-cell lung cancer including stage IA patients. *J Thorac Oncol.* 2014, 9, 59–64. [DOI](https://doi.org/10.1097/JTO.0000000000000042)
22. Schabath MB, Welsh EA, Fulp WJ, et al. Differential association of STK11 and TP53 with KRAS mutation-associated gene expression, proliferation and immune surveillance in lung adenocarcinoma. *Oncogene.* 2016, 35, 3209–3216. [DOI](https://doi.org/10.1038/onc.2015.375)
23. Gillette MA, Satpathy S, Cao S, et al. Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. *Cell.* 2020, 182, 200–225.e35. [DOI](https://doi.org/10.1016/j.cell.2020.06.013)
24. Schemper M, Smith TL. A note on quantifying follow-up in studies of failure time. *Control Clin Trials.* 1996, 17, 343–346. [DOI](https://doi.org/10.1016/0197-2456(96)00075-X)
25. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database Hallmark gene set collection. *Cell Syst.* 2015, 1, 417–425. [DOI](https://doi.org/10.1016/j.cels.2015.12.004)
26. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA.* 2005, 102, 15545–15550. [DOI](https://doi.org/10.1073/pnas.0506580102)
27. Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature.* 2009, 462, 108–112. [DOI](https://doi.org/10.1038/nature08460)
28. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics.* 2013, 14, 7. [DOI](https://doi.org/10.1186/1471-2105-14-7)
29. Cox DR. Regression models and life-tables. *J R Stat Soc Series B.* 1972, 34, 187–202. [DOI](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)
30. Röver C, Knapp G, Friede T. Hartung–Knapp–Sidik–Jonkman approach and its modification for random-effects meta-analysis with few studies. *BMC Med Res Methodol.* 2015, 15, 99. [DOI](https://doi.org/10.1186/s12874-015-0091-1)
31. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Control Clin Trials.* 1986, 7, 177–188. [DOI](https://doi.org/10.1016/0197-2456(86)90046-2)
32. Riley RD, Lambert PC, Abo-Zaid G. Meta-analysis of individual participant data: rationale, conduct, and reporting. *BMJ.* 2010, 340, c221. [DOI](https://doi.org/10.1136/bmj.c221)
33. Heinze G, Schemper M. A solution to the problem of monotone likelihood in Cox regression. *Biometrics.* 2001, 57, 114–119. [DOI](https://doi.org/10.1111/j.0006-341X.2001.00114.x)
34. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Stat Med.* 1996, 15, 361–387. [PMID](https://pubmed.ncbi.nlm.nih.gov/8668867/)
35. Graf E, Schmoor C, Sauerbrei W, Schumacher M. Assessment and comparison of prognostic classification schemes for survival data. *Stat Med.* 1999, 18, 2529–2545. [PMID](https://pubmed.ncbi.nlm.nih.gov/10474158/)
