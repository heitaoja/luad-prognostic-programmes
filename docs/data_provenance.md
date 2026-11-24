# Data provenance and harmonisation

The exact input snapshots, byte sizes, MD5 checksums and download locations are recorded in [data/checksums.csv](../data/checksums.csv). Running the pipeline uses these snapshots. If a remote file changes, retain the original verified snapshot or document and analyse the updated version separately.

## Cohorts

| Cohort | Cohort definition | Clinical fields used for overall survival | Time unit |
|---|---|---|---|
| TCGA-LUAD | Primary tumour barcode 01A, 01B or 01C. Use one specimen per participant. Match exact specimen identifiers first and use the participant barcode for unmatched clinical links. | Xena-mirrored GDC STAR counts, survival and clinical tables. Take OS and OS.time and pathological AJCC stage. | Days |
| GSE31210 | Primary lung tumour with usable OS. Exclude samples flagged for incomplete resection or adjuvant treatment in the series. | Days before death or censoring, death, pathological stage, age, gender. For the secondary relapse endpoint use days before relapse or censoring and the reported relapse indicator. | Days |
| GSE37745 | Adenocarcinoma within the original NSCLC series. | Days to determined death status, dead, tumour stage, age, sex. | Days |
| GSE50081 | Adenocarcinoma within the original resected NSCLC series. | Survival time, vital status, stage, age, sex. The series records time in years. | Years |
| GSE72094 | Public series of resected lung adenocarcinomas. Verify that patient IDs are unique. | Survival time in days, vital status, stage, age at diagnosis, gender. Missing or nonpositive time is ineligible. | Days |
| CPTAC-LUAD | Histologically confirmed primary lung adenocarcinomas sampled at resection. Match tumour expression, clinical and survival files by unique case ID. | LinkedOmics CPTAC pan-cancer data freeze v1.2. Use OS_days and OS_event, with age, sex and AJCC pathological stage from the clinical table. | Days |

Days were divided by 30.4375 and years multiplied by 12. Overall survival events are deaths of any cause. The recorded start of follow-up may differ across sources. The recruitment sources are Japan for GSE31210, Sweden for GSE37745, Toronto for GSE50081 and Moffitt for GSE72094. The four GEO cohorts are combined in the expanded synthesis. TCGA provides discovery and training evidence, and CPTAC provides a separate RNA-seq evaluation. CPTAC collected surgically resected tumours at multiple tissue source sites. The original study confirms LUAD histology and AJCC seventh-edition staging. Public case identifiers offer limited ability to check person-level overlap across repositories.

- [TCGA-LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD) · [GDC TCGA barcode guide](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)
- [GSE31210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210) · [GSE37745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37745) · [GSE50081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50081) · [GSE72094](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72094)
- [CPTAC-LUAD source study](https://doi.org/10.1016/j.cell.2020.06.013) · [LinkedOmics CPTAC pan-cancer LUAD downloads](https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/)

## Annotation, scores and purity

The TCGA Ensembl identifiers were mapped by org.Hs.eg.db after removing version suffixes. GPL570 was mapped with hgu133plus2.db. For GPL15048, the public updated HuRSTA annotation from the original study centre provides probe IDs, the number of mapped genes and symbols. Retain single-gene probes and discard conflicting symbol labels for a probe. The source table is pinned to commit 091f06eb26f0d73824412a9d41d7adbd572a21ae. Its Artistic-2.0 licence belongs to the upstream package. A copy is downloaded on demand into the ignored raw-data folder. [HuRSTA annotation documentation](https://steveneschrich.github.io/Creating-HuRSTA-Packages/) · [GPL15048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15048)

Hallmark definitions come from the Enrichr MSigDB_Hallmark_2020 release. The cached GMT retains the original empty description column. Participant-level derived scores can be regenerated in the local, ignored processed-data folder. [Hallmark collection](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H) · [Enrichr](https://maayanlab.cloud/Enrichr/)

ABSOLUTE estimates are obtained from the Aran et al. pan-cancer data distributed by Bioconductor's TCGAbiolinks. The exact data object is pinned to Git commit f43e902d183671b165ff1eec94271c7e9b9719c9. The source encodes decimal fractions with commas. Convert those fractions and use exact 16-character LUAD specimen barcode matches before fitting paired models. [Aran et al.](https://doi.org/10.1038/ncomms9971) · [TCGAbiolinks data documentation](https://github.com/BioinformaticsFMRP/TCGAbiolinks/blob/master/R/TumorPurityInfo.R)

The CPTAC tumour matrix contains upper-quartile-normalised, log2-transformed RSEM expression estimates. The case IDs matched 110 clinical records and 110 survival records. Five cases had absent or incomplete overall survival. The 105 retained cases had complete clinical covariates and 23 deaths. ACTL10 and NPBWR1 were unavailable after mapping among the 17,595 genes in the original transport reference. The RNA-seq sensitivity used the remaining 17,593, retaining at least 96.5% of genes in every Hallmark set. [CPTAC eligibility and coverage](../results/tables/05b_rnaseq_audit.csv)

## Additional candidate cohorts

GSE68465 remains a candidate for additional replication after auditing the platform mapping and harmonising stage and survival definitions. [Director's Challenge](https://doi.org/10.1038/nm.1790)

## Reproduction

Run `Rscript scripts/fetch_inputs.R` followed by `Rscript scripts/run_pipeline.R` from the repository root. The fetcher verifies existing inputs and stages downloads to temporary files until they match the manifest. The R package versions and platform used for the release are recorded in [sessionInfo.txt](sessionInfo.txt). The full run writes participant-level intermediates under the ignored data directory and summary tables and figures under results.
