# Reporting crosswalk

The manuscript follows the reporting domains of REMARK for tumour markers and TRIPOD/TRIPOD+AI for prognostic models. The external model evaluation separates discrimination, calibration and model updating. [REMARK](https://doi.org/10.1093/jnci/dji237) · [TRIPOD](https://doi.org/10.7326/M14-0697) · [TRIPOD+AI](https://doi.org/10.1136/bmj-2023-078378) · [External validation guidance](https://doi.org/10.1136/bmj-2023-074820)

## Manuscript, results and figures

| Reporting domain | Manuscript | Source outputs |
|---|---|---|
| Cohort selection, missingness, participant counts and follow-up | Results 1, Table 1, Figure 1A | [Data provenance](data_provenance.md) · [Cohort audit](../results/tables/01_cohort_audit.csv) · [TCGA participant flow](../results/tables/01_tcga_participant_flow.csv) · [CPTAC audit](../results/tables/05b_rnaseq_audit.csv) |
| Programme composition, platform coverage and sample ranks | Methods and Results 1, Figure 1B | [Programme definitions](../results/tables/02_axis_definition.csv) · [Hallmark coverage](../results/tables/02_gene_coverage.csv) · [Transport gene list](../results/tables/02_transport_gene_universe.csv) |
| Clinical-adjusted effects, external synthesis and stage-specific estimates | Results 2, Table 2, Figure 2 | [Cohort HRs](../results/tables/04_cohort_estimates.csv) · [CPTAC HRs](../results/tables/05b_rnaseq_association.csv) · [Array meta-analysis](../results/tables/04_meta_analysis.csv) · [IPD and five-cohort synthesis](../results/tables/04b_ipd_and_pan_external_synthesis.csv) · [Stage-specific HRs](../results/tables/03_stage_stratified.csv) |
| Specificity, modelling and scoring sensitivity | Results 2–3, Figure 3A and D | [Diagnostics](../results/tables/03c_model_diagnostics.csv) · [Random-set null](../results/tables/06_random_set_null.csv) · [Outcome permutations](../results/tables/06_permutation_test.csv) · [Scoring and endpoint sensitivity](../results/tables/06_sensitivity.csv) |
| Matched purity, smoking and driver adjustment | Results 4, Tables 3–4, Figure 3B and C | [Purity](../results/tables/03c_purity_sensitivity.csv) · [Smoking and drivers](../results/tables/03d_driver_smoking_adjustment.csv) |
| Comparator definitions, conditional effects and head-to-head performance | Results 5, Table 5, Figure 4 | [Gene lists](../config/published_signatures.csv) · [Simple marker synthesis](../results/tables/03c_marker_synthesis.csv) · [Published benchmark](../results/tables/03e_published_signature_benchmark.csv) · [Published transport](../results/tables/05d_published_signature_transport.csv) |
| Frozen modelling and separate CPTAC-compatible fit | Methods, Results 6 | [GEO score reference](../results/tables/05_frozen_score_reference.csv) · [GEO coefficients](../results/tables/05_frozen_coefficients.csv) · [CPTAC-compatible score reference](../results/tables/05b_rnaseq_score_reference.csv) · [CPTAC-compatible coefficients](../results/tables/05b_rnaseq_frozen_coefficients.csv) |
| External discrimination, 36-month risk and recalibration | Results 6, Tables 6–7, Figures 4B and 5 | [GEO evaluation](../results/tables/05_incremental_value.csv) · [CPTAC evaluation](../results/tables/05b_rnaseq_transport.csv) · [GEO recalibration](../results/tables/05c_recalibration.csv) · [CPTAC recalibration](../results/tables/05b_rnaseq_recalibration.csv) |
| Design history, interpretation and limitations | Methods and Discussion | [Analysis specification](analysis_plan_v2.md) · [Limitations](limitations.md) · [Changelog](changelog.md) |

The five main figures appear as PDF and 300 dpi PNG in the figures directory. The supplementary numbering reserves Figure S1 for the standalone CPTAC evaluation. Figures S2–S12 correspond in order to the retained fig1–fig11 exports. [Main and supplementary legends](manuscript.md#figure-legends) · [Figure directory](../results/figures/)

## Full release workbook

The 30-sheet workbook contains the original summaries and six worksheets for the extension analyses. The CSVs are the machine-readable source of the worksheet values. [Excel workbook](../results/tables/project_tables.xlsx)

| Worksheet | Source |
|---|---|
| S16 Driver and smoking | [03d driver and smoking adjustment](../results/tables/03d_driver_smoking_adjustment.csv) |
| S17 Published signatures | [03e published signature benchmark](../results/tables/03e_published_signature_benchmark.csv) |
| S18 IPD and external meta | [04b IPD and external synthesis](../results/tables/04b_ipd_and_pan_external_synthesis.csv) |
| S19 GEO recalibration | [05c GEO quartile recalibration](../results/tables/05c_recalibration.csv) |
| S20 Signature transport | [05d published signature transport](../results/tables/05d_published_signature_transport.csv) |
| S21 CPTAC recalibration | [05b CPTAC quartile recalibration](../results/tables/05b_rnaseq_recalibration.csv) |
