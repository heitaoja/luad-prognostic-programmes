# Interpretation and limitations

## Study design and populations

This is an observational benchmark of overall survival in cohorts assembled largely from resected lung adenocarcinomas. The programme choices followed TCGA exploration, and the later sensitivities and comparisons were designed with earlier results available. External evidence should be read alongside the distinction between the initial array scope, the extended array scope and the later five-cohort synthesis. Public patient identifiers cannot fully resolve overlap between repositories. Results apply to the observed populations and assays, with uncertain transport to contemporary treatment settings or unresectable disease.

Stage editions, follow-up origins and censoring differ by source. Baseline clinical adjustment used stage, age and sex. Treatment, comorbidity, histologic details and exposure histories were incompletely recorded or inconsistently measured. Additional adjustments are conditional on available fields and on the corresponding complete cases. The GSE31210 recurrence endpoint uses its source-specific relapse annotation. Its interpretation requires the same endpoint definition and independent replication.

## Association, composition and scoring

The two-stage random-effects synthesis has few independent cohorts, leaving heterogeneity uncertain. A zero estimate of $I^2$ does not establish identical effects. The one-stage IPD results estimate a common score slope conditional on the stratification scheme. The two approaches reuse the same external patients. Their intervals answer different questions, particularly for hypoxia, where cohort estimates vary. Core CPTAC individual associations use a coarser stage adjustment because only 23 deaths were available. Some score relationships varied in spline diagnostics, and global PH diagnostics indicate that a single Cox HR may represent a time-averaged association.

ABSOLUTE is available for a subset of exactly matched TCGA specimens and measures malignant fraction at the bulk level. Changes with adjustment can reflect differences in tumour composition, correlation among covariates or other properties of the Cox model. These analyses lack direct measures of cell-type-specific transcription, oxygenation and causal pathways. Driver subgroups have wide intervals, especially GSE72094 EGFR-mutant cases with four deaths. Ridge stabilises the stage and sex terms in that fit, and subgroup differences remain exploratory.

The TCGA gene-wise z-score sensitivity is consistent with attenuation from a common log-expression offset. Sample centring can remove technical and biological shifts together, and the analysis does not directly attribute the offset to sequencing depth. Within-sample ranks are invariant to common strictly monotone transforms, given the same measured genes and ties. Annotation, platform-dependent gene effects, dynamic range and zero counts can still change the ordering. The fixed gene lists and TCGA reference values are part of the model specification.

## Molecular benchmark and prediction

High correlation between the composite, individual pathways and published proliferation signatures limits precision in conditional models. Head-to-head prediction intervals describe the observed uncertainty with the same patients. These secondary comparisons lack prespecified equivalence margins. The public expression implementations of CCP-31 and meta-PCNA share a benchmark pipeline with the Hallmark scores. Performance of an original assay would also depend on its chemistry, reference genes, score construction and clinical thresholds.

Harrell's and Uno's C measure discrimination across the observed follow-up. The latter depends on a reliable estimate of censoring survival, particularly where follow-up is sparse. The 36-month IPCW Brier estimates use marginal test-cohort censoring weights, and observed quartile mortality relies on adequate follow-up within each group. Bootstrap intervals condition on the fixed TCGA fit. CPTAC supports an external RNA-seq check but has limited precision for both discrimination and calibration. Model performance here does not quantify benefit from altering treatment.

The offset update estimates baseline hazard and assesses its performance using the same target outcomes. It preserves patient ordering and can correct mean mortality without correcting calibration slopes or quartile residuals. The reported changes are descriptive and potentially optimistic for new samples. Independent validation would be needed to assess updated absolute risks for a clinical population.

[Manuscript discussion](manuscript.md#discussion) · [Analysis specification](analysis_plan_v2.md) · [Data provenance](data_provenance.md)
