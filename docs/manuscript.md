# Cross-cohort prognostic association and model transport of Hallmark transcriptional programmes in lung adenocarcinoma

## Abstract

**Background.** Many lung adenocarcinoma prognostic signatures capture proliferation. The added information from a multi-pathway proliferation score can be assessed by comparison with simpler markers and by applying the same fitted model to new cohorts.

**Methods.** We analysed overall survival in TCGA-LUAD, four Gene Expression Omnibus array cohorts and the CPTAC-LUAD RNA-seq cohort. We adjusted cohort-specific hazard ratios (HRs) per standard deviation for stage, age and sex and pooled array estimates using REML with modified Knapp–Hartung intervals. We compared the proliferation composite with MKI67, E2F and G2/M, examined DNA-derived tumour purity and evaluated TCGA-fitted models externally.

**Results.** The six cohorts contained 1,443 eligible participants. The four array cohorts contributed 830 complete cases and 269 deaths to the synthesis. Pooled proliferation HRs were 1.32 (95% CI 0.93–1.86) in the original three cohorts and 1.34 (1.10–1.62, Holm-adjusted *P* = .069) in the four-cohort synthesis. In 323 TCGA patients, the proliferation HR was 1.33 before and 1.34 after purity adjustment. In CPTAC (105 patients, 23 deaths), the adjusted proliferation HR was 1.575 (1.002–2.476). The composite improved Harrell C by 0.017–0.054 across GEO cohorts and by 0.010 (−0.071 to 0.104) in CPTAC using a TCGA model fitted to compatible genes. All paired discrimination intervals comparing the composite with simpler markers included zero. Absolute risk calibration varied across cohorts.

**Conclusions.** The strongest evidence was for a directionally reproducible association of proliferation with survival. The incremental performance of the frozen score varied, and the extra complexity of the composite had an uncertain advantage over simpler proliferation measures.

## Introduction

Stage, age and sex are widely available prognostic variables in lung adenocarcinoma. Published expression models have demonstrated associations with survival across datasets, including in the Director's Challenge lung adenocarcinoma study. Proliferation-related expression explains part of the agreement among previously published molecular risk classifiers. Cross-cohort differences in treatment, follow-up, case mix and measurement technology can influence the performance of a fixed prediction model. [Shedden et al.](https://doi.org/10.1038/nm.1790) · [Ringnér and Staaf](https://pubmed.ncbi.nlm.nih.gov/27437773/) · [Riley et al.](https://www.bmj.com/content/384/bmj-2023-074820)

We compared four defined hallmark programmes using the same outcome and clinical adjustment across cohorts. We examined whether a composite proliferation measure retains an association after adjustment for tumour purity, the other programmes or simpler proliferation markers. In a separate analysis, we evaluated the ability of fixed TCGA models to discriminate and calibrate risk in external cohorts. The analysis paired biological programme comparison with explicit tests of assay portability and score redundancy.

## Methods

### Design, cohorts and outcome

This retrospective study used five public cohorts in the main analysis and CPTAC-LUAD in a separate RNA-seq evaluation. TCGA-LUAD served for discovery and model development. GSE31210, GSE37745 and GSE50081 were the original replication sets. GSE72094 expanded the external comparison after we recovered a public, versioned mapping for the custom HuRSTA array. We restricted GSE37745 and GSE50081 to adenocarcinoma and applied the source study's prognosis exclusions to GSE31210. In TCGA we used the first available primary 01A–01C vial in barcode order for each participant with a valid survival time. We checked concordance of survival records across vials and repeated the analysis using the last eligible vial. For CPTAC we linked the LinkedOmics v1.2 tumour RNA-seq, clinical and survival files by case identifier. Overall survival was death from any cause with censoring at last recorded follow-up. We converted time to months and used Kaplan–Meier to estimate median survival and reverse Kaplan–Meier to estimate potential follow-up. [Cohort sources and eligibility](data_provenance.md) · [Gillette et al.](https://doi.org/10.1016/j.cell.2020.06.013)

The four programmes were selected after TCGA exploration during the original project. The revised analysis specification documents corrections and planned comparisons. The expanded analysis preserves the original three-cohort synthesis and reports the estimate after adding GSE72094. [Analysis specification](analysis_plan_v2.md)

### Gene measurement and scoring

We mapped Ensembl identifiers on RNA-seq and probe sets on microarrays to gene symbols, using pinned inputs and software versions. For programme association, the retained feature for each symbol was the one with highest mean expression within its series. For transport, duplicate features were averaged within each specimen on the log-expression scale. The primary transport gene universe was the annotation-based intersection of TCGA, GPL570 and GPL15048. The source GMT was parsed with its description column intact, and each of the 14 member sets was required to reach at least 80% coverage.

The four programmes were formed from 14 Hallmark sets. Proliferation aggregated E2F targets, G2/M checkpoint and MYC targets V1 and V2. The other programmes grouped hypoxia and glycolysis, EMT and stromal activity, and immune and inflammatory transcription. We implemented a rank-integral ssGSEA statistic with exponent 0.25 and mean ranks for ties. For association, we standardised each component within cohort, averaged the components and standardised the resulting score. High expression gets a higher rank contribution. For fixed model transport, both rounds of centring and scaling came from TCGA, along with the reference for MKI67 rank percentile. [Liberzon et al.](https://doi.org/10.1016/j.cels.2015.12.004) · [Barbie et al.](https://doi.org/10.1038/nature08460) · [Implementation details](methods.md)

### Association, synthesis and controls

Separate Cox models included the programme score as a continuous variable and adjusted for pathological stage (I–IV), age and sex in complete cases. We report cohort-specific HRs and 95% confidence intervals. The primary revised synthesis used REML random effects with modified Knapp–Hartung confidence intervals, bounding the variance multiplier below by one. We used Holm adjustment for the four programme tests in each scope. DerSimonian–Laird normal intervals and leave-one-cohort-out synthesis assessed sensitivity. Schoenfeld tests, score splines, stage-stratified baseline hazards and stage I–II restriction assessed model assumptions. We additionally fitted a model with all four programmes.

We compared the composite with MKI67 gene expression and the E2F and G2/M set scores. Conditional models added the composite to each simple measure and the clinical variables on the same cases. For tumour purity, exact specimen matches linked TCGA RNA with copy-number-derived ABSOLUTE values in the Aran pan-cancer data. We used the same subset and score scale for estimates with and without purity adjustment. In TCGA, GSE31210 and GSE72094, random programmes were matched for the number of genes, deciles of mean within-sample rank and overlap between component sets. TCGA time/status permutation and the GSE31210 source relapse endpoint were supporting analyses. [Röver et al.](https://doi.org/10.1186/s12874-015-0091-1) · [Aran et al.](https://doi.org/10.1038/ncomms9971)

### Frozen model transport

We fitted clinical-only and four clinical-plus-molecular Cox models in TCGA and held the score reference, coefficients and baseline cumulative hazards fixed across the four array cohorts. We applied the same gene list, sample-level averaging of duplicate features and ranking across the training and test cohorts. Stage used I as reference, age was entered per decade centred at 65, and sex used female as reference. We assessed Harrell C and paired changes using 1,000 external patient bootstrap resamples with predictions held fixed. Absolute performance at 36 months was assessed using Kaplan–Meier observed mortality, inverse censoring-weighted Brier score and risk-quartile calibration curves. Calibration slopes were estimated as diagnostics. [Riley et al.](https://www.bmj.com/content/384/bmj-2023-074820)

In CPTAC, the association models grouped stage as I versus II–IV to limit the number of coefficients given the event count. For model transport we checked gene availability against the 17,595-gene reference. We fitted a separate set of models in TCGA using the subset shared with CPTAC. The same model formulae and evaluation methods were used. CPTAC supplied the test outcomes.

## Results

### Cohorts and assay coverage

TCGA contained 528 primary tumour specimens, of which 515 linked to valid survival records. Twelve participants had two eligible vials. The primary analysis included 503 distinct participants (182 deaths) and 485 complete cases (177 deaths). Among 442 GSE72094 source samples, 398 had usable OS and 393 also had complete clinical variables. The common measurement universe contained 17,595 genes and covered 96.5–100% of the genes in individual Hallmark sets. Table 1 lists eligibility and follow-up for the five main cohorts. [Participant flow](../results/tables/01_tcga_participant_flow.csv) · [Coverage table](../results/tables/02_gene_coverage.csv)

| Cohort | Patients with OS | Deaths | Model complete cases | Complete-case deaths | Median follow-up, months |
|---|---:|---:|---:|---:|---:|
| TCGA-LUAD | 503 | 182 | 485 | 177 | 29.0 |
| GSE31210 | 204 | 30 | 204 | 30 | 61.0 |
| GSE37745 | 106 | 77 | 106 | 77 | 126.4 |
| GSE50081 | 127 | 51 | 127 | 51 | 65.9 |
| GSE72094 | 398 | 113 | 393 | 111 | 29.9 |
| Main analysis total | 1,338 | 453 | 1,315 | 446 | |

### Associations in the original and expanded external cohorts

The adjusted proliferation HR was 1.31 (1.12–1.54) in TCGA. External proliferation HRs ranged from 1.28 to 1.42. In GSE72094 the HR was 1.36 (1.13–1.62). Across the four external cohorts, the REML estimate was 1.34 (1.10–1.62; nominal *P* = .017; Holm-adjusted *P* = .069). The between-cohort variance estimate was zero under REML. With three original cohorts, the interval was 0.93–1.86. Leave-one-cohort-out expanded estimates ranged from 1.32 to 1.35. Omitting GSE72094 gave the widest proliferation interval. Normal-approximation synthesis yielded narrower intervals and stronger nominal evidence. [Cohort estimates](../results/tables/04_cohort_estimates.csv) · [Synthesis](../results/tables/04_meta_analysis.csv) · [Leave-one-out results](../results/tables/04_leave_one_out.csv)

In the expanded synthesis, hypoxia/glycolysis had HR 1.35 (0.95–1.90; I² 49.5%), EMT/stroma 1.03 (0.74–1.43; I² 56.1%), and immune/inflammation 1.03 (0.73–1.46; I² 62.9%). The extended joint four-programme model gave proliferation HR 1.23 (0.96–1.57). The proliferation HRs under stage-specific baseline hazards were close to the main estimates. Score-specific Schoenfeld tests gave nominal *P* values above .05 for proliferation in every cohort. Global tests indicated potential model departures in TCGA and GSE31210, while spline comparisons indicated possible nonlinearity in several GEO cohorts. The linear estimates provide a consistent summary with these qualifications. [Model diagnostics](../results/tables/03c_model_diagnostics.csv) · [Joint estimates](../results/tables/03c_marker_synthesis.csv)

### Redundancy and robustness

G2/M, E2F and MKI67 also had positive pooled estimates. Expanded HRs were 1.36 (1.12–1.65) for G2/M, 1.32 (1.09–1.60) for E2F and 1.28 (1.04–1.56) for MKI67. The composite conditional on MKI67 had HR 1.30 (0.96–1.76). Conditioning on either of its E2F or G2/M components gave wide intervals. Spearman correlations between the composite and G2/M were 0.94–0.97 across cohorts, and correlations with E2F were 0.97–0.98. [Marker benchmark](../results/tables/03c_marker_synthesis.csv) · [Score correlations](../results/tables/03c_marker_correlations.csv)

Among 323 TCGA patients with paired ABSOLUTE data (120 deaths), adjusting for purity shifted the proliferation HR from 1.33 (1.10–1.62) to 1.34 (1.10–1.62). The Spearman correlation between the score and purity was 0.05. The corresponding correlation was −0.70 for the immune/inflammation score. The hypoxia/glycolysis HR increased from 1.34 (1.11–1.61) to 1.54 (1.27–1.88) in the matched subset after purity adjustment. The TCGA proliferation HR under the alternate vial rule was 1.30 (1.11–1.53). [Purity sensitivity](../results/tables/03c_purity_sensitivity.csv) · [Alternate specimens](../results/tables/03c_specimen_sensitivity.csv)

The empirical upper-tail probabilities for proliferation against structure-matched random programmes were .002 in TCGA, .095 in GSE31210 and .005 in GSE72094. Replacing ssGSEA with average gene-wise z-scores attenuated the TCGA proliferation HR to 1.08 (0.92–1.25), while the alternative scoring result was stronger in some GEO cohorts. In GSE31210 the exploratory relapse model gave HR 1.68 (1.32–2.13) for proliferation. [Matched random controls](../results/tables/06_random_set_null.csv) · [Sensitivity analyses](../results/tables/06_sensitivity.csv)

### Performance after transporting frozen models

The composite raised fixed-model discrimination in all four external cohorts. The gain varied from 0.017 to 0.054 and the paired interval excluded zero in GSE72094. Paired C-index contrasts for the composite against MKI67, E2F and G2/M included zero in every cohort. Prediction error at 36 months, measured with the IPCW Brier score, was lower by 0.004–0.011 for the composite than for the clinical model. Calibration varied. For example, the composite predicted mean 36-month mortality of 25.6% in GSE31210, where the Kaplan–Meier estimate was 8.0%. [Frozen model performance](../results/tables/05_incremental_value.csv) · [Paired marker contrasts](../results/tables/05_simple_marker_contrasts.csv)

| Test cohort | Clinical C | Composite C | Change in C (95% CI) | Predicted vs observed mortality at 36 months |
|---|---:|---:|---:|---:|
| GSE31210 | .688 | .743 | +.054 (−.012 to .119) | 25.6% vs 8.0% |
| GSE37745 | .601 | .618 | +.017 (−.038 to .069) | 37.9% vs 39.6% |
| GSE50081 | .628 | .664 | +.036 (−.007 to .084) | 30.8% vs 24.5% |
| GSE72094 | .640 | .674 | +.034 (.008 to .062) | 37.2% vs 36.6% |

### Independent RNA-seq evaluation

Of 110 CPTAC cases linked to both clinical and survival records, 105 had usable OS and complete covariates. There were 23 deaths and a median potential follow-up of 22.8 months. The adjusted HR was 1.575 (95% CI 1.002–2.476) for proliferation and 2.02 (1.24–3.27) for hypoxia/glycolysis. The EMT/stroma and immune/inflammation intervals included one. These estimates are shown with nominal intervals in Supplementary Figure S1. [CPTAC sample audit](../results/tables/05b_rnaseq_audit.csv) · [CPTAC associations](../results/tables/05b_rnaseq_association.csv)

CPTAC covered 17,593 of the 17,595 symbols in the main transport reference, with the same coverage of the 14 Hallmark sets. The clinical and composite model C-indexes on the compatible gene list were 0.740 and 0.750, respectively. The paired change was 0.010 (95% CI −0.071 to 0.104), and all paired comparisons with simpler markers included zero. At 36 months, 22 deaths had occurred and 27 patients remained at risk. Mean composite-predicted mortality was 31.8% and Kaplan–Meier mortality was 34.4%. [CPTAC transport](../results/tables/05b_rnaseq_transport.csv) · [CPTAC contrasts](../results/tables/05b_rnaseq_marker_contrasts.csv)

## Discussion

The proliferation score had the narrowest range of HR point estimates across the four GEO cohorts, and its estimate in the matched TCGA subset changed little after DNA-derived purity adjustment. CPTAC RNA-seq gave a compatible proliferation estimate with a wide interval. Hypoxia/glycolysis also had an elevated HR in CPTAC, a cohort with relatively few events and a distinct case mix. The original three-cohort synthesis was imprecise, and the expanded proliferation estimate had a Holm-adjusted *P* value of .069.

The G2/M set alone had a pooled association of comparable magnitude and paired transported discrimination close to the composite. The proliferation estimates in mutually adjusted and simpler-marker models were imprecise. The small gain for the composite in CPTAC and the broad paired intervals across assay settings suggest that simpler measures may capture much of the same signal. Absolute risk from the clinical-plus-composite model varied by setting. In particular, the calibration gap in GSE31210 was large compared with the observed event rate.

Compared with studies developing many candidate risk signatures, this analysis focused on four interpretable programme scores and a compact benchmark. The control procedure preserves abundance and observed gene-set overlap, two properties of the molecular measurements. Purity and alternate specimens address distinct sources of bulk-tissue and specimen variation. A frozen evaluation then confronts the same biological hypotheses with a prediction transport task. Together these analyses characterise where the incremental evidence is strongest.

The main limitations are the retrospective selection of the four programmes in TCGA, limited independent cohort counts, differences in follow-up origin and stage measurement, and incomplete treatment and smoking covariates. Public identifiers offer limited ability to check person-level overlap across repositories. The main external synthesis used arrays. The CPTAC analysis used a modified gene reference and had 23 deaths. A larger RNA-seq cohort would improve the precision of independent estimates. Purity data represent a subset of TCGA, and bulk immune/inflammation scores combine expression from multiple cell types. The discrimination bootstrap intervals quantify test-cohort sampling uncertainty conditional on the fitted models. Clinical utility would require prospective evaluation.

## Conclusions

Bulk transcriptional proliferation showed a directionally consistent association with overall survival across four external array cohorts and the CPTAC RNA-seq cohort. The incremental discrimination of TCGA-derived models varied by cohort. The advantage of combining multiple proliferation gene sets over using simpler markers remains uncertain.

## Data and code availability

The accession and exact source snapshots are in [data provenance](data_provenance.md). The analysis can be replayed from the manifest with the R pipeline. Summary results, gene coverage, frozen coefficients and figures are supplied with the repository. Author affiliations, institutional declarations, funding and competing interests require confirmation by the submitting authors.

## Figure legends

1. **Cohort overview.** Participant counts, complete clinical cases and reverse Kaplan–Meier follow-up by cohort.
2. **Programme distributions.** Within-cohort distributions of the four rank-derived programme scores.
3. **External association.** Adjusted HRs per SD and modified Knapp–Hartung intervals for the original and expanded external syntheses. TCGA is shown as discovery evidence.
4. **TCGA survival.** Overall survival curves for above-median and at-or-below-median groups. Primary analysis uses continuous scores.
5. **Stage exploration.** Stage-specific TCGA estimates adjusted for age and sex.
6. **Transported discrimination.** Paired changes in Harrell C for four molecular alternatives to the fixed clinical model.
7. **Matched random programme comparison.** Empirical null distribution of adjusted HRs and the observed score effect.
8. **Sensitivity estimates.** Point estimates across alternate adjustment, scoring and outcome definitions. Intervals are available in the corresponding table.
9. **Simpler-marker benchmark.** External synthesis of the composite, simple measures and conditional composite models.
10. **Absolute risk calibration.** Mean predicted and Kaplan–Meier observed 36-month mortality by predicted-risk quartile.
11. **Purity sensitivity.** Adjusted estimates before and after inclusion of ABSOLUTE estimates in identical TCGA specimens.

**Supplementary Figure S1. Independent RNA-seq evaluation.** Panel A shows adjusted HRs per SD in CPTAC-LUAD with nominal 95% confidence intervals. Panel B shows paired C-index changes relative to the fixed clinical model. The five models were trained in TCGA using the gene list compatible with CPTAC. [Supplementary Figure S1](../results/figures/figS1_rnaseq_validation.pdf)

## Selected references

1. Shedden K, Taylor JMG, Enkemann SA, et al. Gene expression-based survival prediction in lung adenocarcinoma: a multi-site, blinded validation study. *Nat Med*. 2008;14:822–827. [doi:10.1038/nm.1790](https://doi.org/10.1038/nm.1790)
2. Ringnér M, Staaf J. Consensus of gene expression phenotypes and prognostic risk predictors in primary lung adenocarcinoma. *Oncotarget*. 2016. [PMID:27437773](https://pubmed.ncbi.nlm.nih.gov/27437773/)
3. Schabath MB, Welsh EA, Fulp WJ, et al. Differential association of STK11 and TP53 with KRAS mutation-associated gene expression, proliferation and immune surveillance in lung adenocarcinoma. *Oncogene*. 2016;35:3209–3216. [doi:10.1038/onc.2015.375](https://doi.org/10.1038/onc.2015.375)
4. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database hallmark gene set collection. *Cell Syst*. 2015;1:417–425. [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
5. Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*. 2009;462:108–112. [doi:10.1038/nature08460](https://doi.org/10.1038/nature08460)
6. Röver C, Knapp G, Friede T. Hartung–Knapp–Sidik–Jonkman approach and its modification for random-effects meta-analysis with few studies. *BMC Med Res Methodol*. 2015;15:99. [doi:10.1186/s12874-015-0091-1](https://doi.org/10.1186/s12874-015-0091-1)
7. Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nat Commun*. 2015;6:8971. [doi:10.1038/ncomms9971](https://doi.org/10.1038/ncomms9971)
8. Riley RD, Archer L, Snell KIE, et al. Evaluation of clinical prediction models (part 2): how to undertake an external validation study. *BMJ*. 2024;384:e074820. [doi:10.1136/bmj-2023-074820](https://doi.org/10.1136/bmj-2023-074820)
9. McShane LM, Altman DG, Sauerbrei W, et al. Reporting recommendations for tumor marker prognostic studies (REMARK). *J Natl Cancer Inst*. 2005;97:1180–1184. [doi:10.1093/jnci/dji237](https://doi.org/10.1093/jnci/dji237)
10. Gillette MA, Satpathy S, Cao S, et al. Proteogenomic characterization reveals therapeutic vulnerabilities in lung adenocarcinoma. *Cell*. 2020;182:200–225.e35. [doi:10.1016/j.cell.2020.06.013](https://doi.org/10.1016/j.cell.2020.06.013)
