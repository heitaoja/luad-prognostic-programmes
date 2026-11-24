# Reporting crosswalk

This map identifies reporting material for revision against REMARK and against guidance for external model validation. [REMARK](https://doi.org/10.1093/jnci/dji237) · [External validation](https://www.bmj.com/content/384/bmj-2023-074820)

| Reporting topic | Location |
|---|---|
| Marker rationale, outcome and hypotheses | Manuscript, Introduction and Methods; analysis_plan_v2.md |
| Sources, eligibility, specimen selection and missingness | data_provenance.md; Table 1; TCGA participant flow |
| Assay processing, coverage and scoring | methods.md; 02_gene_coverage.csv; fixed transport gene universe |
| Clinical adjustment and estimates with intervals | 04_cohort_estimates.csv; 04_meta_analysis.csv |
| Sensitivity to assumptions, alternative measures and purity | 03c_model_diagnostics.csv; 03c_marker_synthesis.csv; 03c_purity_sensitivity.csv |
| Frozen preprocessing and model definition | 05_frozen_score_reference.csv; 05_frozen_coefficients.csv; code and environment versions |
| External test-cohort discrimination and calibration | 05_incremental_value.csv; 05_calibration.csv; Figures 6 and 10 |
| Independent RNA-seq evaluation and assay compatibility | CPTAC entry in data_provenance.md, 05b_rnaseq_*.csv, and Supplementary Figure S1 |
| Limits of interpretation and implications | Manuscript, Discussion |
