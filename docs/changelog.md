# Study Revision Log

This file records protocol deviations and major revisions relative to the initial pre-specified analysis plan ([analysis_plan.md](analysis_plan.md); current specification in [analysis_plan_v2.md](analysis_plan_v2.md)).

## v1.0 (2026-11-06) — Initial pre-specified array study

- **Scope:** Discovery and training in TCGA-LUAD, with external replication in three GPL570 microarray cohorts (GSE31210, GSE37745, GSE50081).
- **Programme selection:** Initial exploratory work in TCGA considered a dendritic-cell maturation hypothesis and a two-programme proliferation–immune balance score (`z(proliferation) - z(immune)`). In TCGA, the immune score showed no univariate/adjusted association (adjusted HR 0.99, *P* = 0.91) and no proliferation-by-immune interaction (*P* = 0.92). The primary exposure was therefore reframed prior to reading external cohorts as a comparison of four predefined Hallmark programmes (proliferation, hypoxia, EMT/stroma, and immune). Because the four programmes were selected after scanning Hallmark sets in TCGA, TCGA effect estimates are treated as discovery-stage and external pooled estimates provide the primary replication evidence.
- **Excluded candidate:** GSE72094 (Rosetta/Merck HuRSTA custom array) was initially omitted because standard Bioconductor annotation packages did not map its probes to HGNC gene symbols.
- **Exploratory analyses:** Within-stage TCGA models and four-programme mutually adjusted models were designated as exploratory.

## v2.0 (2026-11-24) — Cross-platform benchmark release

- **Participant unit and scoring:** Corrected GMT description parsing for empty columns, adopted midranks for expression ties in ssGSEA, moved TCGA from specimen level to participant level with an alternate-vial sensitivity check, and separated Kaplan–Meier median overall survival from reverse Kaplan–Meier median potential follow-up.
- **Expanded external scope and evidence synthesis:** Restored GSE72094 using the pinned original-centre HuRSTA annotation mapping; added an independent RNA-seq cohort (CPTAC-LUAD) using LinkedOmics pan-cancer freeze v1.2; retained identifiable three- and four-cohort microarray syntheses alongside a five-cohort REML–modified Knapp–Hartung (mKH) synthesis and one-stage stratified IPD Cox models (*k* = 3, 4, 5) with Holm multiplicity adjustment.
- **Confounding and expression-scale sensitivity:** Added exact-specimen DNA ABSOLUTE tumour purity adjustment in TCGA, matched-complete-case smoking and driver-mutation (`EGFR`, `KRAS`, `TP53`, `STK11`) models and exploratory driver subgroup fits, and sample-median-centred gene-wise z-score comparisons to quantify sample-level offset attenuation in RNA-seq counts.
- **Molecular redundancy benchmark:** Added clinical-adjusted, mutually conditional, and external transport comparisons against simple proliferation markers (`MKI67`, `E2F`, `G2/M`) and published signatures (`CCP-31` and `meta-PCNA`) with explicit HGNC alias resolution.
- **Frozen transport and baseline recalibration:** Replaced within-cohort incremental comparisons and target-mean probe selection with a fixed shared gene universe (17,595 genes for TCGA/GEO; 17,593 compatible genes for CPTAC) using within-specimen duplicate averaging, fixed TCGA parameters, paired bootstrap Harrell's C and Uno's IPCW C contrasts, 36-month IPCW Brier scores, and offset-Cox baseline hazard recalibration.
- **Outputs and documentation:** Assembled five multipanel main figures (`Figure1`–`Figure5`), twelve supplementary figures (`FigureS1`–`FigureS12`), a 30-sheet Excel evidence workbook, and accompanying documentation ([manuscript.md](manuscript.md), [reporting_crosswalk.md](reporting_crosswalk.md)).
