# Cross-platform reproducibility, molecular redundancy, and risk transport of Hallmark transcriptional programmes in lung adenocarcinoma

**Running title:** Hallmark transcriptional programmes in lung adenocarcinoma

**Authors:** [To be completed by authors]

**Affiliations:** [To be completed by authors]

**Corresponding author:** [To be completed by authors]

**Keywords:** Lung adenocarcinoma; Transcriptional signatures; Prognostic benchmarking; External validation; Calibration; Tumour purity; ssGSEA

---

## Abstract

**Background:** Although numerous transcriptional prognostic signatures have been proposed for lung adenocarcinoma (LUAD), their cross-platform reproducibility, sensitivity to tumour composition and expression scaling, overlap with cell-cycle markers, and out-of-sample calibration remain incompletely resolved.

**Methods:** We evaluated 1,443 patients with primary LUAD across six cohorts (1,420 complete clinical cases; 469 deaths), reserving TCGA-LUAD for discovery and model development while holding out four Gene Expression Omnibus (GEO) microarray series (GSE31210, GSE37745, GSE50081, and GSE72094) and the CPTAC-LUAD RNA-seq cohort for external evaluation. Four transcriptional programmes (proliferation, hypoxia/glycolysis, EMT/stroma, and immune) built from 14 MSigDB Hallmark sets were scored using tie-aware single-sample gene set enrichment analysis (ssGSEA) and pooled with two-stage restricted maximum likelihood–modified Knapp–Hartung (REML-mKH) meta-analysis and one-stage stratified individual participant data (IPD) Cox models. We quantified sensitivity to DNA-derived ABSOLUTE purity, smoking history, somatic drivers, and sample-level expression centring, benchmarked redundancy against simple markers and published signatures (CCP-31 and meta-PCNA), and evaluated frozen TCGA-derived models using Harrell's C, Uno's inverse probability of censoring weighted (IPCW) C, and 36-month baseline hazard recalibration.

**Results:** Proliferation was consistently adverse across external cohorts, with REML-mKH hazard ratios (HRs) per standard deviation of 1.32 (*k* = 3), 1.34 (*k* = 4), and 1.35 (*k* = 5; 95% CI 1.15–1.59, Holm-adjusted *P* = 0.027, *I*² = 0%). In cohort-stratified IPD models, proliferation and hypoxia yielded HRs of 1.34 and 1.36 at *k* = 4 (*N* = 830; adjusted *P* = 2.3 × 10⁻⁶ and 1.1 × 10⁻⁵) and 1.35 and 1.40 at *k* = 5 (*N* = 935; both adjusted *P* = 4.4 × 10⁻⁷), whereas EMT/stroma and immune programmes showed no reproducible pooled association. Adjusting for ABSOLUTE purity in matched TCGA specimens left the proliferation HR unchanged (1.33 to 1.34) but increased the hypoxia HR from 1.34 to 1.54. Proliferation estimates were largely stable after smoking and driver-mutation adjustments, whereas sample-median centring recovered the TCGA gene-wise z-score proliferation HR from 1.08 to 1.24. The proliferation composite correlated strongly with rank-based CCP-31 and meta-PCNA (Spearman ρ = 0.927–0.965) and showed no consistent discrimination advantage over either comparator in external testing. Adding the composite to the frozen clinical model improved Uno's C by 0.071 in CPTAC (95% CI −0.023 to 0.192) and by 0.051 in GSE72094 (95% CI 0.006 to 0.102). In GSE31210, baseline recalibration reduced 36-month mean predicted mortality from 25.6% to 8.2% (observed 8.0%) and IPCW Brier error from 0.0935 to 0.0675.

**Conclusions:** Proliferation is the most reproducible Hallmark transcriptional programme in LUAD, sharing nearly identical prognostic information with established proliferation signatures. While relative risk ranking transports across microarray and RNA-seq platforms, accurate absolute risk prediction requires local baseline hazard updating, which here was evaluated in the same target data used for recalibration.

---

## Introduction

Patients with lung adenocarcinoma (LUAD) experience markedly heterogeneous survival even within the same pathological stage, motivating molecular predictors that capture tumour cell state beyond anatomical staging. Early multi-institutional efforts such as the Director's Challenge demonstrated that gene-expression models can inform LUAD prognosis when subjected to blinded external validation [1], and subsequent assays adapted expression signatures to formalin-fixed clinical specimens in resected non-squamous lung cancer [2]. Parallel work translated a 31-gene cell-cycle progression (CCP-31) signature—originally developed for prostate cancer [4]—into early-stage LUAD [3]. However, establishing whether newly proposed pathway scores contribute reproducible, non-redundant clinical value requires evaluating them head-to-head against both standard clinicopathological factors and credible molecular comparators, while distinguishing exploratory expression scores from locked diagnostic assays.

Interpreting multi-gene prognostic signatures is complicated by pervasive biological redundancy and cellular composition effects. In breast cancer, Venet and colleagues demonstrated that most randomly selected gene sets correlate with survival simply by tracking a core mitotic programme represented by proliferating cell nuclear antigen (meta-PCNA) [5], and Ringnér and Staaf similarly showed that proliferation dominates consensus survival signals across published LUAD expression phenotypes [6]. Combining multiple pathways into a composite score may stabilise measurement across platforms, yet an elevated hazard ratio alone does not prove incremental value over a single proliferation gene or pathway unless tested in conditional and paired external prediction models. At the same time, bulk tumour profiles mix malignant cell states with stromal and immune infiltrates. Pan-cancer analyses pairing DNA-derived ABSOLUTE purity estimates [8] with bulk transcriptomes have revealed how tumour cellularity systematically shifts transcriptional associations [7]. Together with geographical differences in smoking exposure and somatic driver alterations, tissue composition can mask or confound programme associations unless evaluated in matched-patient sensitivity models.

Cross-platform transport introduces further statistical and measurement hurdles. Sequencing and microarray cohorts differ in gene coverage, dynamic range, and the influence of library depth on log-transformed expression matrices. Because gene-wise z-score standardisation centres genes across patients rather than within specimens, unnormalised sample-level offsets remain embedded in each gene's standardised value [9]; within-sample rank statistics are invariant to common strictly monotone transformations, though ties and gene-specific assay effects can still alter rankings. Furthermore, as emphasised by REMARK [10], TRIPOD [11], and TRIPOD+AI [12] guidelines as well as contemporary external validation frameworks [13, 15, 16], a marker with a stable hazard ratio may still fail to improve censoring-adjusted discrimination [14] or may miscalibrate absolute survival probabilities when baseline hazard differs across clinical settings.

To examine these issues within a single reproducible framework, we conducted a retrospective cross-platform benchmark built around four predefined Hallmark transcriptional programmes in 1,443 patients across six LUAD cohorts. Using TCGA-LUAD strictly for discovery and model fitting and holding out four GEO microarray cohorts and the CPTAC-LUAD RNA-seq cohort, we evaluated (1) cross-platform reproducibility under two-stage random-effects and one-stage individual participant data (IPD) syntheses, (2) sensitivity to ABSOLUTE purity, smoking history, somatic drivers, and expression-scale centring, (3) molecular redundancy against simple markers (`MKI67`, `E2F`, `G2/M`) and published signatures (`CCP-31`, `meta-PCNA`), and (4) frozen risk-model discrimination and 36-month baseline hazard recalibration.

---

## Methods

### Study design, cohorts, and clinical endpoint

We analysed primary tumour transcriptomes and clinical follow-up from six public LUAD cohorts: TCGA-LUAD [17, 18], four Gene Expression Omnibus (GEO) microarray series (GSE31210 [19], GSE37745 [20], GSE50081 [21], and GSE72094 [22]), and CPTAC-LUAD [23]. TCGA-LUAD served exclusively as the discovery and model-training cohort, whereas all five external cohorts—spanning three primary Affymetrix GPL570 series (GSE31210, GSE37745, GSE50081), an extended Rosetta/Merck HuRSTA series (GSE72094), and an independent RNA-seq cohort (CPTAC-LUAD)—were strictly excluded from model fitting and reference derivation. Because the four Hallmark programmes were chosen after an initial Hallmark scan in TCGA-LUAD, external syntheses provide the primary test of prognostic reproducibility. All input file snapshots, byte sizes, and MD5 checksums were pinned prior to analysis.

For TCGA-LUAD, primary tumour specimens (barcodes ending in `-01A`, `-01B`, or `-01C`) were linked to survival and clinical records by exact sample identifier and, where unmatched, by 12-character participant barcode [18]. After confirming survival concordance across repeated vials from the same participant, we retained the first eligible vial in barcode order as the primary participant unit and tested the last eligible vial in sensitivity analyses. In external cohorts, we retained source-defined exclusions for incomplete resection or adjuvant therapy in GSE31210 [19], restricted GSE37745 and GSE50081 to adenocarcinoma histology [20, 21], verified unique patient identifiers in GSE72094 [22], and linked CPTAC-LUAD expression, clinical, and survival tables by unique case identifier [23]. Overall survival (OS) was defined as months from study baseline to death from any cause, with censoring at last contact (days divided by 30.4375; years multiplied by 12). Patients with positive OS time and known vital status entered the cohort audit, and multivariable Cox models used complete cases for the specified clinical covariates. Median potential follow-up was estimated using the reverse Kaplan–Meier method [24].

### Expression preprocessing and transcriptional programme scoring

TCGA-LUAD expression comprised UCSC Xena-mirrored GDC STAR counts on the log₂(count + 1) scale, whereas CPTAC-LUAD comprised upper-quartile-normalised, log₂-transformed RSEM estimates from LinkedOmics (pan-cancer freeze v1.2). GEO cohorts were analysed on log₂-transformed expression matrices mapped using `hgu133plus2.db` for GPL570 and the pinned original-centre HuRSTA annotation for GPL15048 (retaining single-gene probes without conflicting symbols). TCGA Ensembl identifiers were mapped to HGNC symbols using `org.Hs.eg.db` after stripping version suffixes. For cohort-specific association analyses, duplicated symbols were resolved by retaining the feature with the highest cohort mean log expression; for frozen cross-cohort model transport, duplicate features mapping to the same symbol were averaged within each specimen prior to ranking so that preprocessing remained completely independent of target-cohort summary statistics.

Four transcriptional programmes were constructed from 14 gene sets in the MSigDB Hallmark 2020 collection [25]: proliferation (`E2F_TARGETS`, `G2M_CHECKPOINT`, `MYC_TARGETS_V1`, `MYC_TARGETS_V2`), hypoxia (`HYPOXIA`, `GLYCOLYSIS`), EMT/stroma (`EPITHELIAL_MESENCHYMAL_TRANSITION`, `TGF_BETA_SIGNALING`, `ANGIOGENESIS`), and immune (`INFLAMMATORY_RESPONSE`, `INTERFERON_ALPHA_RESPONSE`, `INTERFERON_GAMMA_RESPONSE`, `ALLOGRAFT_REJECTION`, `COMPLEMENT`). For a sample measuring *G* genes with descending expression midranks *r_g*, weights *v_g* = *G* − *r_g* + 1 were computed so that tied expression values received identical contributions [26–28]. Each component ssGSEA score was calculated as the difference between the mean of *v_g* weighted by *v_g*⁰·²⁵ inside the set and the unweighted mean of *v_g* outside the set, divided by *G*. Within each cohort, component scores were standardised to zero mean and unit variance, averaged with equal weights across member sets, and standardised again to express hazard ratios (HRs) per within-cohort standard deviation (SD).

To dissect expression-scale effects, we compared ssGSEA against the mean of gene-wise z-scores across the union of non-constant measured genes in each programme under two variants: direct cross-patient standardisation and a sample-centred variant that subtracted each specimen's median expression across all measured genes before gene-wise standardisation. Under an approximate multiplicative library-size model at moderate-to-high counts, log₂ profiles carry a sample-wide additive shift [9] that persists after gene-wise centring across patients; subtracting the within-sample median removes any common additive shift on the log scale (though it can also remove genuine global shifts and degrades at low counts owing to the pseudocount), whereas within-sample midranks are strictly invariant to any strictly increasing column transformation.

### Clinical-adjusted survival modelling and external evidence synthesis

Separate Cox proportional hazards models [29] were fitted for each programme with OS in months, adjusting for pathological stage (I–IV), cohort-standardised age, and sex. Because CPTAC-LUAD contained 23 deaths, its standalone cohort model collapsed stage into stage I versus stages II–IV to limit degrees of freedom; four-level stage coding was retained in all pooled IPD and transport models. Model assumptions and specification sensitivity were examined via score-specific and global Schoenfeld residual tests, 3-degree-of-freedom natural splines, stage-stratified baseline hazards, restriction to stages I–II, joint four-programme mutual adjustment, alternative TCGA vial selection, and exploratory within-stage TCGA models adjusted for age and sex. Median-split Kaplan–Meier curves were generated solely for visual illustration; all statistical inference used continuous scores without cutoff optimisation.

External evidence was synthesised across three predefined scopes excluding TCGA-LUAD: the original GPL570 array cohorts (*k* = 3; *N* = 437, 158 deaths), the expanded microarray cohorts (*k* = 4; *N* = 830, 269 deaths), and all five external cohorts including CPTAC-LUAD (*k* = 5; *N* = 935, 292 deaths). In two-stage random-effects meta-analyses, cohort-specific log HRs were pooled using restricted maximum likelihood (REML) heterogeneity estimation with modified Knapp–Hartung (mKH) confidence intervals (variance multiplier bounded below by 1; *t* reference distribution with *k* − 1 degrees of freedom) [30], supplemented by DerSimonian–Laird (DL) normal-approximation estimates [31] and leave-one-cohort-out analyses. In complementary one-stage IPD Cox models [32], complete cases were pooled across cohorts using within-cohort standardised scores and age alongside four-level stage and sex, fitting both cohort-stratified (`strata(cohort)`) and cohort-by-stage-stratified (`strata(cohort, stage)`) baseline hazards. Within each external scope and synthesis method, *P* values across the four programmes were adjusted using the Holm procedure.

### Tumour purity, clinical confounding, and specificity controls

DNA-derived ABSOLUTE tumour purity estimates from Aran et al. [7, 8] (pinned via `TCGAbiolinks`) were matched to TCGA-LUAD using exact 16-character specimen barcodes. Spearman correlations between purity and programme scores were evaluated alongside paired Cox fits in the exact same matched specimens before and after adding purity. Identical matched-patient comparisons were conducted for smoking and somatic driver covariates: standardised pack-years in TCGA-LUAD; ever- versus never-smoking and grouped alteration status (`EGFR`, `KRAS`/`ALK`, or triple-negative) in GSE31210; smoking status in GSE50081; and smoking status as well as a separate four-driver model (`EGFR`, `KRAS`, `TP53`, and `STK11`) in GSE72094. Exploratory driver-stratified subgroup analyses in GSE31210 and GSE72094 retained cohort score scaling; when a stage or sex stratum contained zero events, ridge penalisation (θ = 1, unscaled) was applied strictly to the stage and sex nuisance columns while leaving programme and age coefficients unpenalised [33].

Specificity was evaluated against structure-matched random programmes generated by replacing each gene in a programme's union one-to-one with a non-member measured gene from the same mean-rank expression decile, reusing identical replacements across overlapping member sets to preserve set sizes and inter-set overlap. We drew 500 matched null programmes in TCGA-LUAD and 200 each in GSE31210 and GSE72094 for proliferation and hypoxia, computing upper-tail empirical *P* values as (1 + *b*)/(1 + *B*). As a pipeline control in TCGA-LUAD, survival time and event indicators were permuted together relative to predictors 1,000 times per programme. Additional sensitivity analyses tested 1%/99% score winsorisation and the exploratory relapse-free survival endpoint in GSE31210.

### Molecular redundancy benchmarking

We benchmarked the 4-Hallmark proliferation composite against `MKI67` expression, single `E2F` and `G2/M` Hallmark scores, the 31-gene Wistuba LUAD `CCP-31` signature [3, 4], and the 131-gene Venet `meta-PCNA` signature [5] after resolving historical HGNC symbols using explicit alias mappings. For each comparator, we audited gene coverage, computed Spearman correlations with the composite, and fitted clinical-adjusted as well as bidirectional mutually conditional Cox models across cohorts using both tie-aware rank scoring and sample-centred mean gene-wise z-scores. Pooled estimates across the three external scopes used REML-mKH, DL, and cohort-stratified IPD models with nominal inference. In external transport comparisons, `CCP-31` and `meta-PCNA` were evaluated both as rank scores and as TCGA-standardised within-specimen mean log-expression scores. These comparisons evaluate shared gene-set information on common expression matrices rather than proprietary PCR assay implementations.

### Frozen risk-model transport and baseline hazard recalibration

To ensure zero target-cohort leakage during external prediction, we defined a fixed transport universe by intersecting annotated symbols across TCGA-LUAD and the four GEO platforms (17,595 genes; minimum Hallmark coverage 96.5%). Because `ACTL10` and `NPBWR1` were unavailable in CPTAC-LUAD, a separate TCGA reference and model fit were derived on the 17,593 compatible genes (minimum Hallmark coverage 96.5%) prior to evaluating CPTAC-LUAD. All component means and SDs, `MKI67` within-sample percentile-rank references, Cox regression coefficients, and 36-month cumulative baseline hazards were estimated solely in TCGA-LUAD (clinical reference: stage I, female sex, age per decade centred at 65 years) and held fixed during external evaluation.

Discrimination of each frozen linear predictor η_i was quantified using Harrell's concordance index [34] and Uno's IPCW concordance index (time weights proportional to *n*(*t*)/*Ĝ*(*t*)²) [14], with 95% percentile confidence intervals for ΔC relative to the frozen clinical model and for head-to-head molecular contrasts derived from 1,000 paired bootstrap resamples of target-cohort patients without model refitting. Absolute 36-month mortality calibration and IPCW Brier scores [35] were assessed in external cohorts meeting pre-specified stability criteria (≥ 10 deaths and ≥ 10 patients at risk at 36 months). Calibration slopes were estimated diagnostically by regressing target OS on η_i. For baseline hazard recalibration ("recalibration-in-the-large") [13, 15, 16], we fitted a target-cohort Cox model with η_i as a fixed offset (slope fixed at 1), extracted the updated 36-month cumulative baseline hazard *Ĥ*₀,recal(36) at η = 0, and computed updated individual 36-month risks as 1 − exp[−*Ĥ*₀,recal(36) exp(η_i)]. Both frozen and recalibrated quartile calibration curves used identical risk-group membership and test-cohort censoring weights; accordingly, recalibrated absolute risks reflect in-sample target updating rather than independent post-update validation.

### Statistical software

All analyses were conducted in R (v4.6.2) using `survival` (v3.8.12), `data.table` (v1.18.6.1), `ggplot2` (v4.0.3), `patchwork` (v1.3.2), `openxlsx` (v4.2.9), `AnnotationDbi` (v1.74.0), `org.Hs.eg.db` (v3.23.1), and `hgu133plus2.db` (v3.13.0), with random seed `20261124` and two-sided hypothesis tests at α = 0.05.

---

## Results

### Cohort characteristics and platform coverage

Across the six cohorts, 1,443 patients had valid OS data (476 deaths), of whom 1,420 complete clinical cases (469 deaths) entered covariate-adjusted analyses (Table 1; Supplementary Fig. S2). In TCGA-LUAD, 515 of 528 primary tumour specimens had usable OS; resolving 12 participants with repeated vials yielded 503 unique patients and 182 deaths, with 485 complete cases and 177 deaths (Supplementary Table S1). All 105 eligible CPTAC-LUAD cases (out of 110 source cases; 23 deaths) had complete clinical covariates (Supplementary Table S11). Reverse Kaplan–Meier median potential follow-up varied substantially across cohorts, from 22.8 months in CPTAC-LUAD and 29.0 months in TCGA-LUAD to 126.4 months in GSE37745 (Table 1; Fig. 1A). Together with the low 36-month event rate in GSE31210 (16 deaths by 36 months among 204 patients), this heterogeneity motivated both censoring-weighted discrimination metrics and explicit baseline calibration assessment.

**Table 1.** Cohort characteristics, complete-case counts, and reverse Kaplan–Meier median potential follow-up.

| Cohort | Role | Platform | Eligible OS | Deaths | Complete cases | Complete-case deaths | Median follow-up (months) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| TCGA-LUAD | Discovery / training | RNA-seq (STAR) | 503 | 182 | 485 | 177 | 29.0 |
| GSE31210 | External (GPL570) | Microarray | 204 | 30 | 204 | 30 | 61.0 |
| GSE37745 | External (GPL570) | Microarray | 106 | 77 | 106 | 77 | 126.4 |
| GSE50081 | External (GPL570) | Microarray | 127 | 51 | 127 | 51 | 65.9 |
| GSE72094 | External (HuRSTA) | Microarray | 398 | 113 | 393 | 111 | 29.9 |
| CPTAC-LUAD | External (RNA-seq) | RNA-seq (RSEM) | 105 | 23 | 105 | 23 | 22.8 |
| **Total** | — | — | **1,443** | **476** | **1,420** | **469** | — |

All 14 Hallmark gene sets achieved ≥ 96.5% coverage across both fixed transport universes (17,595 and 17,593 genes; Supplementary Table S1). Although within-cohort standardisation anchored all programme scores to mean 0 and SD 1, score distributions differed in skewness and tail behaviour across platforms (Fig. 1B; Supplementary Fig. S3). In TCGA-LUAD, switching from the first to the last eligible vial produced virtually identical results on the same 485 patients (proliferation HR 1.30, 95% CI 1.11–1.53 vs 1.31, 95% CI 1.12–1.54; Supplementary Tables S7–S8).

### Cross-platform programme associations and evidence synthesis

In TCGA-LUAD, clinical-adjusted HRs were 1.31 (95% CI 1.12–1.54, *P* = 0.00077) for proliferation and 1.42 (95% CI 1.22–1.66, *P* = 9.7 × 10⁻⁶) for hypoxia (Fig. 2A; Supplementary Fig. S5). Both programmes yielded point estimates exceeding 1.0 in every external cohort, including CPTAC-LUAD (proliferation HR 1.57, 95% CI 1.00–2.48, *P* = 0.049, exact lower bound 1.002; hypoxia HR 2.02, 95% CI 1.24–3.27, *P* = 0.0045; Supplementary Fig. S1A and Supplementary Table S12). By contrast, EMT/stroma and immune estimates varied in direction across cohorts (Fig. 2A; Supplementary Fig. S4). In exploratory within-stage TCGA models adjusted for age and sex, proliferation HRs were 1.29 (stage I), 1.10 (stage II), and 1.68 (stage III), with no evidence of stage interaction (likelihood ratio *P* = 0.310; Fig. 2C and Supplementary Fig. S6). Although the TCGA proliferation model showed global Schoenfeld non-proportionality (*P* = 0.039; score-specific *P* = 0.081), stratifying the baseline hazard by stage yielded an identical HR of 1.31 (Supplementary Table S3).

Across external evidence syntheses (Table 2; Fig. 2A), two-stage REML-mKH HRs for proliferation were remarkably consistent across *k* = 3 (1.32, 95% CI 0.93–1.86; Holm *P* = 0.297), *k* = 4 (1.34, 95% CI 1.10–1.62; Holm *P* = 0.069), and *k* = 5 (1.35, 95% CI 1.15–1.59; Holm *P* = 0.027), with zero estimated between-cohort heterogeneity (*I*² = 0%). Hypoxia yielded a five-cohort REML-mKH HR of 1.42 (95% CI 1.04–1.94; Holm *P* = 0.103) alongside substantial heterogeneity (*I*² = 57.6%). One-stage IPD Cox models strongly supported both proliferation (cohort-stratified *k* = 5 HR 1.35, 95% CI 1.21–1.51, Holm *P* = 4.4 × 10⁻⁷) and hypoxia (HR 1.40, 95% CI 1.23–1.58, Holm *P* = 4.4 × 10⁻⁷), with nearly identical estimates under cohort-by-stage stratification (Table 2; Supplementary Table S18). In exploratory four-programme joint models across the four array cohorts, pooled REML-mKH HRs attenuated to 1.23 (95% CI 0.96–1.57) for proliferation and 1.25 (95% CI 0.94–1.65) for hypoxia. Neither EMT/stroma nor immune achieved significance in any pooled external model.

**Table 2.** External hazard ratios per SD (95% CI) across evidence scopes and synthesis models, with Holm-adjusted *P* values in brackets.

| External cohorts (patients/deaths) | Synthesis model | Proliferation | Hypoxia | EMT/stroma | Immune |
| --- | --- | --- | --- | --- | --- |
| 3 array (437/158) | Two-stage REML-mKH | 1.32 (0.93–1.86) [0.30] | 1.37 (0.66–2.82) [0.61] | 1.02 (0.52–2.00) [1.00] | 1.09 (0.57–2.08) [1.00] |
| 3 array (437/158) | One-stage IPD, cohort | 1.34 (1.15–1.56) [9.0 × 10⁻⁴] | 1.33 (1.12–1.58) [0.0038] | 0.99 (0.83–1.18) [1.00] | 1.05 (0.90–1.24) [1.00] |
| 3 array (437/158) | One-stage IPD, cohort × stage | 1.31 (1.12–1.53) [0.0035] | 1.31 (1.10–1.56) [0.0069] | 0.97 (0.82–1.16) [1.00] | 1.05 (0.89–1.22) [1.00] |
| 4 array (830/269) | Two-stage REML-mKH | 1.34 (1.10–1.62) [0.069] | 1.35 (0.95–1.90) [0.21] | 1.03 (0.74–1.43) [1.00] | 1.03 (0.73–1.46) [1.00] |
| 4 array (830/269) | One-stage IPD, cohort | 1.34 (1.20–1.51) [2.3 × 10⁻⁶] | 1.36 (1.19–1.54) [1.1 × 10⁻⁵] | 1.02 (0.90–1.16) [1.00] | 0.99 (0.88–1.12) [1.00] |
| 4 array (830/269) | One-stage IPD, cohort × stage | 1.32 (1.18–1.49) [1.3 × 10⁻⁵] | 1.34 (1.18–1.52) [3.0 × 10⁻⁵] | 1.02 (0.90–1.15) [1.00] | 1.00 (0.89–1.13) [1.00] |
| 5 all external (935/292) | Two-stage REML-mKH | 1.35 (1.15–1.59) [0.027] | 1.42 (1.04–1.94) [0.10] | 1.09 (0.81–1.46) [0.97] | 1.02 (0.79–1.31) [0.97] |
| 5 all external (935/292) | Two-stage DL (normal) | 1.35 (1.20–1.51) [1.1 × 10⁻⁶] | 1.42 (1.16–1.73) [0.0017] | 1.08 (0.88–1.33) [0.88] | 1.02 (0.86–1.21) [0.88] |
| 5 all external (935/292) | One-stage IPD, cohort | 1.35 (1.21–1.51) [4.4 × 10⁻⁷] | 1.40 (1.23–1.58) [4.4 × 10⁻⁷] | 1.05 (0.93–1.19) [0.84] | 1.00 (0.89–1.12) [0.99] |
| 5 all external (935/292) | One-stage IPD, cohort × stage | 1.34 (1.19–1.50) [2.1 × 10⁻⁶] | 1.38 (1.22–1.56) [2.0 × 10⁻⁶] | 1.05 (0.93–1.19) [0.83] | 1.00 (0.89–1.13) [0.96] |

### Programme specificity and sensitivity to expression scaling

Against structure-matched random gene sets preserving gene counts, abundance deciles, and component overlap, observed HRs exceeded null expectations in TCGA-LUAD (empirical *P* = 0.0020 for both proliferation and hypoxia across 500 draws) and GSE72094 (empirical *P* = 0.0050 for both across 200 draws; Fig. 3A and Supplementary Fig. S8). In the smaller GSE31210 cohort (30 deaths), empirical *P* values were 0.0945 for proliferation and 0.0199 for hypoxia, indicating stronger specificity above random gene sets for hypoxia than for proliferation in that series. Paired outcome permutations in TCGA yielded empirical *P* = 0.0010 for both programmes. In the exploratory GSE31210 relapse endpoint (54 events), HRs strengthened to 1.68 (95% CI 1.32–2.13) for proliferation and 1.83 (95% CI 1.37–2.43) for hypoxia (Supplementary Table S5).

Choice of expression-scoring scale produced pronounced differences in TCGA-LUAD (Fig. 3D; Supplementary Fig. S9). Using uncentred mean gene-wise z-scores attenuated the TCGA proliferation HR to 1.08 (95% CI 0.92–1.25, *P* = 0.35) and hypoxia HR to 1.03 (95% CI 0.89–1.20, *P* = 0.67). Subtracting each specimen's median expression prior to gene-wise standardisation restored these estimates to 1.24 (95% CI 1.06–1.45, *P* = 0.0071) and 1.21 (95% CI 1.04–1.41, *P* = 0.012), approaching tie-aware ssGSEA (1.31 and 1.42). Consistent upward recovery with sample centring was observed in GSE50081 (proliferation HR 1.11 to 1.25; hypoxia HR 1.06 to 1.43) and GSE72094 (proliferation HR 1.23 to 1.38), whereas in GSE31210 centring slightly reduced the mean gene-wise z-score proliferation HR from 1.43 (95% CI 1.02–2.00) to 1.33 (95% CI 0.96–1.84).

### Influence of tumour purity, smoking history, and driver alterations

Among 323 TCGA-LUAD patients (120 deaths) with exact-specimen ABSOLUTE DNA purity estimates, proliferation was essentially uncorrelated with purity (ρ = +0.05) and its HR remained invariant after purity adjustment (1.33 to 1.34; Table 3 and Fig. 3B). Conversely, hypoxia, EMT/stroma, and immune scores all correlated negatively with tumour purity (ρ = −0.29, −0.50, and −0.70, respectively). Controlling for purity increased the hypoxia HR from 1.34 (95% CI 1.11–1.61) to 1.54 (95% CI 1.27–1.88) and shifted immune (0.98 to 1.26, 95% CI 0.96–1.64) and EMT/stroma (1.05 to 1.22, 95% CI 0.97–1.52) upwards (Supplementary Fig. S12).

**Table 3.** Paired TCGA-LUAD Cox models before and after ABSOLUTE purity adjustment (*N* = 323, 120 deaths).

| Programme | Spearman ρ with purity | Clinical (95% CI) | Clinical + purity (95% CI) |
| --- | ---: | --- | --- |
| Proliferation | +0.05 | 1.33 (1.10–1.62) | 1.34 (1.10–1.62) |
| Hypoxia | −0.29 | 1.34 (1.11–1.61) | 1.54 (1.27–1.88) |
| EMT/stroma | −0.50 | 1.05 (0.88–1.27) | 1.22 (0.97–1.52) |
| Immune | −0.70 | 0.98 (0.82–1.17) | 1.26 (0.96–1.64) |

In matched complete-case comparisons across TCGA-LUAD, GSE31210, GSE50081, and GSE72094, adding smoking history (pack-years or ever/never status) or somatic driver terms produced modest shifts in proliferation and hypoxia point estimates (Table 4; Fig. 3C; Supplementary Table S16). Jointly adjusting for `EGFR`, `KRAS`, `TP53`, and `STK11` status in GSE72094 (*N* = 393, 111 deaths) shifted proliferation from 1.36 to 1.38 (95% CI 1.13–1.68) and hypoxia from 1.36 to 1.30 (95% CI 1.07–1.58). In exploratory driver-stratified models, proliferation yielded HR 2.10 (95% CI 1.10–4.03) among 116 `EGFR`-mutant GSE31210 cases (12 deaths) versus 1.22 (95% CI 0.36–4.14) among 41 `EGFR`-mutant GSE72094 cases (4 deaths, ridge-stabilised), while `STK11`-mutant GSE72094 cases (64 patients, 18 deaths) showed wide intervals overlapping unity (proliferation HR 0.84, 95% CI 0.48–1.47; hypoxia HR 0.89, 95% CI 0.52–1.51 vs 1.43 and 1.41 in `STK11` wild-type).

**Table 4.** Matched complete-case sensitivity to additional smoking and driver-mutation adjustments.

| Cohort | Added covariates | Patients | Deaths | Proliferation HR (95% CI) | Hypoxia HR (95% CI) |
| --- | --- | ---: | ---: | --- | --- |
| TCGA-LUAD | Pack-years | 333 | 118 | 1.33 (1.09–1.62) → 1.34 (1.10–1.63) | 1.37 (1.14–1.65) → 1.37 (1.14–1.65) |
| GSE31210 | Smoking + alteration group | 204 | 30 | 1.42 (1.01–1.99) → 1.38 (0.98–1.95) | 1.96 (1.33–2.90) → 2.03 (1.37–3.00) |
| GSE50081 | Smoking | 115 | 46 | 1.20 (0.92–1.57) → 1.19 (0.90–1.57) | 1.25 (0.91–1.72) → 1.25 (0.91–1.72) |
| GSE72094 | Smoking | 328 | 89 | 1.46 (1.19–1.78) → 1.46 (1.19–1.79) | 1.39 (1.12–1.73) → 1.39 (1.12–1.73) |
| GSE72094 | `EGFR`, `KRAS`, `TP53`, `STK11` | 393 | 111 | 1.36 (1.13–1.62) → 1.38 (1.13–1.68) | 1.36 (1.12–1.65) → 1.30 (1.07–1.58) |

### Redundancy with simple markers and published proliferation signatures

Across cohorts, the 4-Hallmark proliferation composite correlated tightly with `E2F` (Spearman ρ = 0.97–0.98) and `G2/M` (ρ = 0.94–0.97), moderately with `MKI67` (ρ = 0.69–0.85), and strongly with rank-based `CCP-31` (ρ = 0.93–0.96) and `meta-PCNA` (ρ = 0.94–0.96) across all six cohorts (Supplementary Table S10 and Supplementary Table S17). Every comparator showed positive clinical-adjusted survival associations in four-array REML-mKH synthesis (Table 5; Fig. 4A; Supplementary Fig. S10). However, collinearity severely limited precision in mutually conditional models: conditional on `G2/M`, `CCP-31`, or `meta-PCNA`, the composite HR crossed unity in both four- and five-cohort REML-mKH syntheses (Table 5). Under five-cohort cohort-stratified IPD, clinical-adjusted HRs were 1.45 (95% CI 1.29–1.64) for `CCP-31` and 1.39 (95% CI 1.24–1.56) for `meta-PCNA`; conditional on the composite, `CCP-31` retained HR 1.91 (95% CI 1.34–2.73) under IPD whereas its REML-mKH interval spanned 1 (HR 2.61, 95% CI 0.73–9.37; Supplementary Table S17).

**Table 5.** Clinical-adjusted and mutually conditional four-array and five-cohort REML-mKH associations.

| Comparator score | External cohorts | Clinical + score HR (95% CI) | Composite given score HR (95% CI) | Score given composite HR (95% CI) |
| --- | ---: | --- | --- | --- |
| Composite | 4 | 1.34 (1.10–1.62) | — | — |
| `MKI67` | 4 | 1.28 (1.04–1.56) | 1.30 (0.96–1.76) | — |
| `E2F` | 4 | 1.32 (1.09–1.60) | 1.26 (0.44–3.57) | — |
| `G2/M` | 4 | 1.36 (1.12–1.65) | 0.89 (0.37–2.11) | — |
| `CCP-31` | 4 | 1.42 (1.14–1.77) | 0.65 (0.16–2.66) | 2.31 (0.43–12.41) |
| `meta-PCNA` | 4 | 1.37 (1.12–1.67) | 0.76 (0.14–4.10) | 1.89 (0.29–12.39) |
| `CCP-31` | 5 | 1.45 (1.20–1.74) | 0.60 (0.21–1.74) | 2.61 (0.73–9.37) |
| `meta-PCNA` | 5 | 1.38 (1.17–1.64) | 0.78 (0.23–2.57) | 1.89 (0.50–7.17) |

In external frozen transport across all five held-out cohorts, paired Harrell ΔC 95% bootstrap intervals comparing the composite head-to-head against rank-based `CCP-31`, rank-based `meta-PCNA`, `MKI67`, `E2F`, and `G2/M` all included zero (Fig. 4B; Supplementary Tables S6, S14, S20). Only one of ten exploratory contrasts using mean log-expression scoring favoured the composite nominally (`meta-PCNA` in GSE72094: ΔC = 0.020, 95% CI 0.005–0.036; Supplementary Table S20).

### External discrimination, calibration, and baseline recalibration

Applying the frozen TCGA models across all five external cohorts yielded positive point gains in both Harrell's C and Uno's IPCW C relative to the frozen clinical model in every test cohort (Table 6; Fig. 4B; Supplementary Figs. S1B and S7). Paired bootstrap 95% intervals excluded zero for both estimators in GSE72094 (Harrell ΔC = +0.034, 95% CI 0.008–0.062; Uno ΔC = +0.051, 95% CI 0.006–0.102). In CPTAC-LUAD, where follow-up was shortest (22.8 months) and deaths were fewest (23 events), Harrell ΔC was +0.010 (95% CI −0.071 to 0.104) whereas Uno's IPCW ΔC reached +0.071 (95% CI −0.023 to 0.192).

**Table 6.** External discrimination of frozen TCGA clinical and clinical-plus-composite models.

| External cohort | Harrell C, clinical | Harrell C, composite | Harrell ΔC (95% CI) | Uno C, clinical | Uno C, composite | Uno ΔC (95% CI) |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| GSE31210 | 0.688 | 0.743 | +0.054 (−0.012 to 0.119) | 0.676 | 0.743 | +0.068 (−0.042 to 0.176) |
| GSE37745 | 0.601 | 0.618 | +0.017 (−0.038 to 0.069) | 0.607 | 0.623 | +0.016 (−0.037 to 0.065) |
| GSE50081 | 0.628 | 0.664 | +0.036 (−0.007 to 0.084) | 0.625 | 0.683 | +0.059 (−0.025 to 0.148) |
| GSE72094 | 0.640 | 0.674 | +0.034 (0.008 to 0.062) | 0.601 | 0.653 | +0.051 (0.006 to 0.102) |
| CPTAC-LUAD | 0.740 | 0.750 | +0.010 (−0.071 to 0.104) | 0.694 | 0.765 | +0.071 (−0.023 to 0.192) |

All five external cohorts satisfied stability criteria for 36-month absolute risk evaluation (Supplementary Tables S4, S5, S9, S13, S15). Direct transport of the TCGA 36-month baseline hazard substantially overpredicted mortality in GSE31210, where Kaplan–Meier observed mortality was 8.0% versus 28.6% predicted by the frozen clinical model and 25.6% by the frozen composite model (Table 7; Fig. 5; Supplementary Fig. S11). Baseline hazard recalibration aligned both mean predictions to 8.2% and lowered the composite IPCW Brier score from 0.0935 to 0.0675 (and from 0.1047 to 0.0683 for the clinical model; Supplementary Table S19). Smaller Brier reductions followed recalibration in GSE50081 (0.1736 to 0.1711) and CPTAC-LUAD (0.1722 to 0.1679; Supplementary Table S21), whereas GSE72094 was already well calibrated before updating (37.2% vs 36.6% observed). Diagnostic calibration slopes varied across cohorts (e.g., 1.73 for the composite in GSE31210 vs 0.66 in GSE37745), indicating that correcting mean risk does not remove differences in relative-risk spread across risk quartiles (Fig. 5A).

**Table 7.** Frozen versus baseline-recalibrated 36-month predicted mortality and IPCW Brier scores.

| External cohort | Model | Observed KM | Frozen mean risk | Recalibrated mean risk | Frozen Brier | Recalibrated Brier |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
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

---

## Discussion

Across six independent RNA-seq and microarray cohorts totalling 1,443 patients, proliferation emerged as the most reproducible Hallmark transcriptional programme associated with overall survival in LUAD. Why two-stage REML-mKH inference shifted from Holm-adjusted *P* = 0.297 (*k* = 3) and 0.069 (*k* = 4) to 0.027 (*k* = 5) despite near-identical pooled HRs (1.32, 1.34, and 1.35) reflects the small-sample conservatism of modified Knapp–Hartung *t* intervals when *k* is small, even at *I*² = 0% [30]. Adding CPTAC-LUAD improved degrees of freedom while extending replication to an external RNA-seq platform. One-stage stratified IPD Cox models [32] provided highly consistent point estimates (1.35 for proliferation, 1.40 for hypoxia at *k* = 5), though readers should note that IPD and two-stage syntheses reuse the same external patients and that IPD imposes a common score slope even where cohort-specific hypoxia effects diverge (*I*² = 57.6% in two-stage pooling).

Decomposing biological and technical confounders clarified why transcriptional signatures behave differently across cohorts. Proliferation showed negligible correlation with DNA-derived ABSOLUTE purity (ρ = +0.05) and remained stable after purity adjustment, consistent with a predominantly tumour-intrinsic signal [7, 8]. By contrast, hypoxia/glycolysis correlated inversely with purity (ρ = −0.29) and strengthened markedly after adjustment (HR 1.34 to 1.54), suggesting that varying stromal dilution obscures hypoxic stress signals in bulk profiles. Similarly, immune expression tracked non-malignant cellularity closely (ρ = −0.70), and its null unadjusted association (HR 0.98) shifted only modestly and imprecisely after purity adjustment (HR 1.26, 95% CI 0.96–1.64). Across cohorts with markedly different clinical backgrounds—such as the predominantly never-smoking, `EGFR`-mutant Japanese GSE31210 series [19] versus the smoking- and `KRAS`/`STK11`-enriched US GSE72094 series [22] and TCGA pack-years subset [17]—matched-case adjustment for smoking and somatic drivers altered proliferation estimates only slightly, indicating that proliferative state captures prognostic variation beyond recorded clinical and mutational covariates.

Our scoring-scale comparison also demonstrates a practical preprocessing pitfall in RNA-seq survival studies. On TCGA log₂(STAR count + 1) profiles lacking upstream library-size normalisation, standardising each gene across patients leaves any sample-wide additive log offset intact in every standardised gene value [9]. Averaging uncentred gene-wise z-scores attenuated the TCGA proliferation HR from 1.31 to 1.08 and hypoxia from 1.42 to 1.03, whereas subtracting each specimen's median expression before gene standardisation recovered HRs to 1.24 and 1.21. Because within-sample midranks are mathematically invariant to strictly monotone column shifts, rank-based ssGSEA [26–28] avoids this artefact naturally, reinforcing the need to specify exact scoring universes, normalisation transforms, and frozen references alongside any expression predictor.

Benchmarking against simpler and established comparators strongly favours parsimony for clinical risk ranking. Correlations between the four-Hallmark proliferation composite and `E2F`, `G2/M`, `CCP-31` [3, 4], and `meta-PCNA` [5] reached 0.93–0.98, echoing earlier observations that diverse prognostic signatures converge onto a shared cell-cycle axis [5, 6]. In conditional models, collinearity yielded wide intervals, and across all ten external head-to-head frozen transport tests, paired Harrell ΔC intervals against rank-based `CCP-31` or `meta-PCNA` included zero. Combining multiple proliferation Hallmarks therefore offers little reproducible gain over a single cell-cycle pathway or published gene set on public expression matrices, though translation of a commercial PCR assay depends on platform chemistry and clinical thresholds rather than gene lists alone [2, 3].

Finally, separating discrimination from calibration [10–16] highlights that relative risk ordering transfers more readily across assays than absolute survival probability. Using within-specimen duplicate averaging on a fixed shared gene universe ensured strict zero-leakage transport, yielding positive point increments in Harrell's C [34] and Uno's IPCW C [14] across all five test cohorts, including +0.071 Uno ΔC in CPTAC-LUAD where early censoring magnified the divergence from Harrell's C. Yet direct transport of TCGA's 36-month baseline hazard drastically overpredicted mortality in low-risk early-stage GSE31210 (25.6% predicted vs 8.0% observed), and simple recalibration-in-the-large [13, 15] corrected mean predicted risk (8.2%) and IPCW Brier error (0.0935 to 0.0675) [35] without altering risk ranking or resolving slope differences across risk quartiles [16].

Several limitations qualify these conclusions. First, the four Hallmark programmes were selected after exploratory screening in TCGA-LUAD, and our observational cohorts consist largely of surgically resected tumours with heterogeneous staging editions, follow-up origins, incomplete treatment records, and limited power within rare driver strata (notably four deaths in GSE72094 `EGFR`-mutant cases requiring ridge stabilisation [33]). Public identifiers cannot completely exclude patient overlap across repositories, and transport to unresectable or contemporary immunotherapy-treated populations remains untested. Second, ABSOLUTE purity provides only a bulk malignant-fraction proxy rather than cell-type-specific expression, and formal equivalence margins were not pre-specified for head-to-head molecular comparisons. Third, our 36-month baseline hazard updates were estimated and evaluated within the same target cohorts; independent prospective evaluation is required before clinical deployment of recalibrated absolute risks.

---

## Conclusions

Among four predefined Hallmark transcriptional programmes in LUAD, proliferation demonstrates the strongest cross-platform reproducibility across microarray and RNA-seq cohorts, remaining robust to tumour purity, smoking, and driver-mutation adjustments. However, the composite proliferation score exhibits near-complete redundancy with single cell-cycle pathways (`E2F`, `G2/M`) and established signatures (`CCP-31`, `meta-PCNA`) and yields modest discrimination gains over clinical staging alone. While relative risk ranking transports reliably under zero-leakage frozen scoring, accurate absolute mortality estimation requires local baseline hazard recalibration validated in independent target populations.

---

## Declarations

### Ethics approval and consent to participate

Not applicable. All analyses were performed exclusively on publicly available, de-identified human transcriptomic and clinical datasets; ethics approval and informed consent were obtained by the original contributing studies (TCGA-LUAD, GSE31210, GSE37745, GSE50081, GSE72094, and CPTAC-LUAD).

### Consent for publication

Not applicable.

### Availability of data and materials

All datasets analysed in this study are publicly available from their source repositories: TCGA-LUAD via the NCI Genomic Data Commons / UCSC Xena (`TCGA-LUAD.star_counts.tsv.gz`, `TCGA-LUAD.survival.tsv.gz`, `TCGA-LUAD.clinical.tsv.gz`), NCBI Gene Expression Omnibus under accessions `GSE31210`, `GSE37745`, `GSE50081`, and `GSE72094` (with `GPL15048` HuRSTA probe annotation), LinkedOmics CPTAC-LUAD pan-cancer data freeze v1.2, MSigDB Hallmark 2020 gene sets via Enrichr, and pan-cancer ABSOLUTE tumour purity estimates from Aran et al. via Bioconductor `TCGAbiolinks`. Pinned input checksums, analysis scripts, and reproducible output tables are provided in the accompanying project repository (`data/checksums.csv`, `scripts/`, `R/`, and `results/tables/project_tables.xlsx`).

### Competing interests

The authors declare that they have no competing interests.

### Funding

[To be completed by authors]

### Authors' contributions

[To be completed by authors]

### Acknowledgements

The authors thank the TCGA Research Network, CPTAC Consortium, LinkedOmics, and the contributing investigators of GSE31210, GSE37745, GSE50081, and GSE72094 for making primary genomic and clinical data publicly accessible.

---

## Figure Legends

### Figure 1. Cohort composition and programme score distributions

**A,** Eligible overall survival (OS) counts and complete clinical case counts across the six LUAD cohorts, paired with reverse Kaplan–Meier median potential follow-up in months. **B,** Density distributions of within-cohort standardised ssGSEA scores for proliferation, hypoxia/glycolysis, EMT/stroma, and immune programmes across TCGA-LUAD and the five external cohorts.

### Figure 2. Programme survival associations and cross-platform synthesis

**A,** Clinical-adjusted hazard ratios (HRs) per within-cohort SD with 95% confidence intervals across all six cohorts, two-stage REML-mKH syntheses for three, four, and five external cohorts, and cohort-stratified one-stage IPD Cox syntheses for four and five external cohorts. Individual CPTAC estimates adjust for stage I versus II–IV, age, and sex; dashed vertical line indicates HR = 1. **B,** Kaplan–Meier overall survival curves in TCGA-LUAD stratified at the cohort median for proliferation and hypoxia (visualisation only; crosses denote censoring). **C,** Exploratory within-stage TCGA HRs (stages I, II, and III) adjusted for age and sex.

### Figure 3. Programme specificity, tumour purity, clinical confounders, and scoring-scale sensitivity

**A,** Null HR distributions from structure-matched random gene sets preserving component size, abundance decile, and gene overlap; vertical lines mark observed HRs in TCGA-LUAD, GSE31210, and GSE72094. **B,** Paired TCGA estimates comparing clinical adjustment against clinical plus ABSOLUTE tumour purity adjustment in identical matched specimens (*N* = 323), annotated with Spearman correlations between score and purity. **C,** Matched complete-case comparisons with additional smoking or somatic driver covariates alongside exploratory driver-stratified subgroup HRs in GSE31210 and GSE72094 (ridge-stabilised stage/sex terms in sparse `EGFR`-mutant GSE72094). **D,** Sensitivity heatmap across alternative endpoints, clinical adjustments, winsorisation, and gene-wise z-score variants (uncentred vs sample-median-centred).

### Figure 4. Molecular redundancy and external concordance contrasts

**A,** Four-array REML-mKH clinical-adjusted and mutually conditional HRs comparing the proliferation composite against `MKI67`, `E2F`, `G2/M`, `CCP-31`, and `meta-PCNA`. **B,** Paired bootstrap 95% confidence intervals for Harrell's ΔC and Uno's IPCW ΔC relative to the frozen clinical model across the four GEO cohorts and CPTAC-LUAD (`CCP-31` and `meta-PCNA` scored via ssGSEA).

### Figure 5. Frozen versus baseline-recalibrated 36-month mortality prediction

**A,** Quartile calibration plots comparing mean predicted 36-month mortality against Kaplan–Meier observed mortality (with 95% Greenwood intervals) before and after offset-Cox baseline hazard recalibration for the clinical and composite models across all five external cohorts. **B,** Cohort-level shifts in 36-month mean predicted mortality and IPCW Brier scores (with positive Brier reductions denoting lower prediction error after baseline recalibration).

---

## Supplementary Information

### Supplementary Figure Legends

**Figure S1. Independent RNA-seq evaluation in CPTAC-LUAD.** **A,** Clinical-adjusted HRs for the four Hallmark programmes in CPTAC-LUAD (*N* = 105, 23 deaths). **B,** Paired bootstrap Harrell's ΔC intervals relative to the frozen clinical model derived on the 17,593-gene compatible TCGA universe.

**Figure S2. Cohort sizes, complete cases, and potential follow-up.** Comparison of OS-eligible participants, complete clinical cases, and reverse Kaplan–Meier median follow-up across cohorts.

**Figure S3. Cohort-standardised programme distributions.** Density plots comparing the four transcriptional programme distributions across TCGA-LUAD and external array cohorts.

**Figure S4. Individual cohort and array meta-analysis forest plots.** Cohort-level HRs alongside three- and four-cohort REML-mKH external syntheses with 95% confidence intervals.

**Figure S5. TCGA-LUAD Kaplan–Meier curves by programme median.** Exploratory visualisation comparing patients above versus at or below the median score for each programme.

**Figure S6. Stage-stratified exploratory associations in TCGA-LUAD.** Within-stage HRs adjusted for age and sex across stages I, II, and III.

**Figure S7. External Harrell's C increments over the clinical model.** Paired 1,000-resample bootstrap intervals for frozen TCGA models evaluated in GEO microarray cohorts.

**Figure S8. Structure-matched random-programme null distributions.** Empirical null HR histograms for proliferation and hypoxia matched on gene-set size, rank decile, and overlap.

**Figure S9. Comprehensive sensitivity matrix.** HRs across alternative survival endpoints, covariate specifications, stage restrictions, and scoring methods.

**Figure S10. Benchmark against simpler proliferation measures.** Forest plots of clinical-adjusted and conditional HRs comparing the composite with `MKI67`, `E2F`, and `G2/M`.

**Figure S11. Frozen 36-month quartile calibration in microarray cohorts.** Predicted versus observed 36-month mortality across quartiles of frozen risk.

**Figure S12. DNA-derived tumour purity sensitivity in TCGA-LUAD.** Matched-specimen HRs before and after ABSOLUTE purity adjustment across programmes and simple proliferation markers (*N* = 323).

### Supplementary Tables

All supporting summary tables are compiled in `results/tables/project_tables.xlsx` alongside machine-readable CSV files in `results/tables/`: **Tables T1–T9** (main summary extracts corresponding to cohort audits, adjusted HRs, meta-analyses, transport, sensitivities, null controls, stage stratification, marker benchmarks, and purity) and **Supplementary Tables S1–S21** (gene-set coverage, leave-one-out meta-analyses, Schoenfeld/spline diagnostics, frozen TCGA coefficients and score references, simple marker contrasts, TCGA vial sensitivity, paired specimen concordance, held-out calibration quartiles, marker correlations, CPTAC-LUAD audit/association/transport/recalibration tables, smoking and driver adjustments, published signature associations and transport, IPD and five-cohort syntheses, and GEO quartile recalibration).

---

## References

1. Shedden K, Taylor JM, Enkemann SA, et al. Gene expression-based survival prediction in lung adenocarcinoma: a multi-site, blinded validation study. *Nat Med.* 2008;14:822–7. [https://doi.org/10.1038/nm.1790](https://doi.org/10.1038/nm.1790)
2. Kratz JR, He J, Van Den Eeden SK, et al. A practical molecular assay to predict survival in resected non-squamous, non-small-cell lung cancer: development and international validation studies. *Lancet.* 2012;379:823–32. [https://pubmed.ncbi.nlm.nih.gov/22285053/](https://pubmed.ncbi.nlm.nih.gov/22285053/)
3. Wistuba II, Behrens C, Lombardi F, et al. Validation of a proliferation-based expression signature as prognostic marker in early stage lung adenocarcinoma. *Clin Cancer Res.* 2013;19:6261–71. [https://doi.org/10.1158/1078-0432.CCR-13-0596](https://doi.org/10.1158/1078-0432.CCR-13-0596)
4. Cuzick J, Swanson GP, Fisher G, et al. Prognostic value of an RNA expression signature derived from cell cycle proliferation genes in patients with prostate cancer: a retrospective study. *Lancet Oncol.* 2011;12:245–55. [https://doi.org/10.1016/S1470-2045(10)70295-3](https://doi.org/10.1016/S1470-2045(10)70295-3)
5. Venet D, Dumont JE, Detours V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput Biol.* 2011;7:e1002240. [https://doi.org/10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240)
6. Ringnér M, Staaf J. Consensus of gene expression phenotypes and prognostic risk predictors in primary lung adenocarcinoma. *Oncotarget.* 2016;7:52957–73. [https://doi.org/10.18632/oncotarget.10641](https://doi.org/10.18632/oncotarget.10641)
7. Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nat Commun.* 2015;6:8971. [https://doi.org/10.1038/ncomms9971](https://doi.org/10.1038/ncomms9971)
8. Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA alterations in human cancer. *Nat Biotechnol.* 2012;30:413–21. [https://doi.org/10.1038/nbt.2203](https://doi.org/10.1038/nbt.2203)
9. Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol.* 2010;11:R25. [https://doi.org/10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)
10. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). *J Natl Cancer Inst.* 2005;97:1180–4. [https://doi.org/10.1093/jnci/dji237](https://doi.org/10.1093/jnci/dji237)
11. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a multivariable prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. *Ann Intern Med.* 2015;162:55–63. [https://doi.org/10.7326/M14-0697](https://doi.org/10.7326/M14-0697)
12. Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ.* 2024;385:e078378. [https://doi.org/10.1136/bmj-2023-078378](https://doi.org/10.1136/bmj-2023-078378)
13. Riley RD, Archer L, Snell KIE, et al. Evaluation of clinical prediction models (part 2): how to undertake an external validation study. *BMJ.* 2024;384:e074820. [https://doi.org/10.1136/bmj-2023-074820](https://doi.org/10.1136/bmj-2023-074820)
14. Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med.* 2011;30:1105–17. [https://doi.org/10.1002/sim.4154](https://doi.org/10.1002/sim.4154)
15. Steyerberg EW. *Clinical Prediction Models: A Practical Approach to Development, Validation, and Updating.* 2nd ed. Cham: Springer; 2019. [https://doi.org/10.1007/978-3-030-16399-0](https://doi.org/10.1007/978-3-030-16399-0)
16. Van Calster B, McLernon DJ, van Smeden M, Wynants L, Steyerberg EW. Calibration: the Achilles heel of predictive analytics. *BMC Med.* 2019;17:230. [https://doi.org/10.1186/s12916-019-1466-7](https://doi.org/10.1186/s12916-019-1466-7)
17. Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. *Nature.* 2014;511:543–50. [https://doi.org/10.1038/nature13385](https://doi.org/10.1038/nature13385)
18. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell.* 2018;173:400–416.e11. [https://doi.org/10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052)
19. Okayama H, Kohno T, Ishii Y, et al. Identification of genes upregulated in ALK-positive and EGFR/KRAS/ALK-negative lung adenocarcinomas. *Cancer Res.* 2012;72:100–11. [https://doi.org/10.1158/0008-5472.CAN-11-1403](https://doi.org/10.1158/0008-5472.CAN-11-1403)
20. Botling J, Edlund K, Lohr M, et al. Biomarker discovery in non-small cell lung cancer: integrating gene expression profiling, meta-analysis, and tissue microarray validation. *Clin Cancer Res.* 2013;19:194–204. [https://doi.org/10.1158/1078-0432.CCR-12-1139](https://doi.org/10.1158/1078-0432.CCR-12-1139)
21. Der SD, Sykes J, Pintilie M, et al. Validation of a histology-independent prognostic gene signature for early-stage, non-small-cell lung cancer including stage IA patients. *J Thorac Oncol.* 2014;9:59–64. [https://doi.org/10.1097/JTO.0000000000000042](https://doi.org/10.1097/JTO.0000000000000042)
22. Schabath MB, Welsh EA, Fulp WJ, et al. Differential association of STK11 and TP53 with KRAS mutation-associated gene expression, proliferation and immune surveillance in lung adenocarcinoma. *Oncogene.* 2016;35:3209–16. [https://doi.org/10.1038/onc.2015.375](https://doi.org/10.1038/onc.2015.375)
23. Gillette MA, Satpathy S, Cao S, et al. Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. *Cell.* 2020;182:200–225.e35. [https://doi.org/10.1016/j.cell.2020.06.013](https://doi.org/10.1016/j.cell.2020.06.013)
24. Schemper M, Smith TL. A note on quantifying follow-up in studies of failure time. *Control Clin Trials.* 1996;17:343–6. [https://doi.org/10.1016/0197-2456(96)00075-X](https://doi.org/10.1016/0197-2456(96)00075-X)
25. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database Hallmark gene set collection. *Cell Syst.* 2015;1:417–25. [https://doi.org/10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
26. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci U S A.* 2005;102:15545–50. [https://doi.org/10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)
27. Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature.* 2009;462:108–12. [https://doi.org/10.1038/nature08460](https://doi.org/10.1038/nature08460)
28. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics.* 2013;14:7. [https://doi.org/10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)
29. Cox DR. Regression models and life-tables. *J R Stat Soc Series B.* 1972;34:187–202. [https://doi.org/10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)
30. Röver C, Knapp G, Friede T. Hartung–Knapp–Sidik–Jonkman approach and its modification for random-effects meta-analysis with few studies. *BMC Med Res Methodol.* 2015;15:99. [https://doi.org/10.1186/s12874-015-0091-1](https://doi.org/10.1186/s12874-015-0091-1)
31. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Control Clin Trials.* 1986;7:177–88. [https://doi.org/10.1016/0197-2456(86)90046-2](https://doi.org/10.1016/0197-2456(86)90046-2)
32. Riley RD, Lambert PC, Abo-Zaid G. Meta-analysis of individual participant data: rationale, conduct, and reporting. *BMJ.* 2010;340:c221. [https://doi.org/10.1136/bmj.c221](https://doi.org/10.1136/bmj.c221)
33. Heinze G, Schemper M. A solution to the problem of monotone likelihood in Cox regression. *Biometrics.* 2001;57:114–9. [https://doi.org/10.1111/j.0006-341X.2001.00114.x](https://doi.org/10.1111/j.0006-341X.2001.00114.x)
34. Harrell FE Jr, Lee KL, Mark DB. Multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing errors. *Stat Med.* 1996;15:361–87. [https://pubmed.ncbi.nlm.nih.gov/8668867/](https://pubmed.ncbi.nlm.nih.gov/8668867/)
35. Graf E, Schmoor C, Sauerbrei W, Schumacher M. Assessment and comparison of prognostic classification schemes for survival data. *Stat Med.* 1999;18:2529–45. [https://pubmed.ncbi.nlm.nih.gov/10474158/](https://pubmed.ncbi.nlm.nih.gov/10474158/)
