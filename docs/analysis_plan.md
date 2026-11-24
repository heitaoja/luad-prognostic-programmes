# Pre-specified Analysis Plan (v1.0, Archived)

> **Note:** This document archives the initial pre-specified protocol finalised on 2026-11-06 prior to analysing the external replication cohorts. The active specification incorporating cross-platform transport, GSE72094, CPTAC-LUAD, IPD synthesis, published-signature benchmarking, and baseline recalibration is documented in [analysis_plan_v2.md](analysis_plan_v2.md) and implemented in [R/lib_score.R](../R/lib_score.R).

## Scope and limitations of pre-specification

- **Locked before external evaluation:** Programme scoring rules, primary Cox model specification, clinical covariates, continuous exposure coding, survival endpoint, random-effects pooling, negative controls, and sensitivity variants were fixed before external cohort evaluation.
- **Discovery selection bias:** The four transcriptional programmes were chosen from the MSigDB Hallmark collection after an initial scan in TCGA-LUAD. TCGA effect sizes are therefore subject to discovery selection bias, and external cohorts provide the primary replication evidence (see [changelog.md](changelog.md)).

## 1. Initial cohorts

| Cohort | Role | Inclusion criteria |
|---|---|---|
| TCGA-LUAD | Discovery | Primary tumour (`-01`), valid overall survival (OS) time and vital status |
| GSE31210 | External validation | Primary lung tumour, excluding cases flagged by source authors for incomplete resection or adjuvant therapy |
| GSE37745 | External validation | Adenocarcinoma histology within the NSCLC series, valid OS |
| GSE50081 | External validation | Adenocarcinoma histology within the NSCLC series, valid OS |

At v1.0, GSE72094 was deferred because the custom HuRSTA platform lacked a standard Bioconductor annotation package.

## 2. Predeclared exposures

| Programme | MSigDB Hallmark component gene sets |
|---|---|
| `proliferation` | `E2F_TARGETS`, `G2M_CHECKPOINT`, `MYC_TARGETS_V1`, `MYC_TARGETS_V2` |
| `hypoxia` | `HYPOXIA`, `GLYCOLYSIS` |
| `emt_stroma` | `EPITHELIAL_MESENCHYMAL_TRANSITION`, `TGF_BETA_SIGNALING`, `ANGIOGENESIS` |
| `immune` | `INFLAMMATORY_RESPONSE`, `INTERFERON_GAMMA_RESPONSE`, `INTERFERON_ALPHA_RESPONSE`, `ALLOGRAFT_REJECTION`, `COMPLEMENT` |

**Scoring rule:**
1. Compute per-sample ssGSEA enrichment ($\alpha = 0.25$, normalised by gene count $G$) for each component gene set.
2. Standardise each component score to mean 0 and SD 1 within each cohort.
3. Average component z-scores with equal weights and standardise the composite score within each cohort.

## 3. Primary endpoint, model and inference

- **Endpoint:** Overall survival (months; all-cause death).
- **Exposure:** Continuous programme score (HR per within-cohort SD). No dichotomisation or optimal cutoff search is permitted for inference; median splits are restricted to Kaplan–Meier visualisation.
- **Model:** `coxph(Surv(os_time, os_event) ~ programme + stage + age + sex)`.
- **Synthesis:** Random-effects meta-analysis across external cohorts with heterogeneity ($I^2$) reported alongside pooled effects and negative results.

## 4. Secondary, exploratory and negative-control analyses

- **Secondary:** Incremental concordance beyond clinical covariates, likelihood ratio tests, and clinical-adjusted comparisons against simpler proliferation markers.
- **Exploratory:** Stage-stratified TCGA models (adjusted for age and sex), four-programme mutually adjusted models, and calibration assessment.
- **Negative controls:** Structure-matched random gene sets (500 draws in TCGA, 200 in external cohorts) and 1,000 paired survival-outcome permutations in TCGA.
- **Sensitivities:** Stage-only adjustment, restriction to stages I–II, alternative scoring methods (mean gene-wise z-scores), 1%/99% winsorisation, and relapse-free survival in GSE31210.
