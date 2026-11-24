# Independent RNA-seq evaluation on the public CPTAC-LUAD outcome snapshot.
source("R/lib_inputs.R")
source("R/lib_score.R")
source("R/lib_surv.R")
source("R/lib_transport.R")
id <- "CPTAC-LUAD"
tab <- "results/tables"
raw <- "data/raw/cptac"
inputs <- input_manifest()
inputs <- inputs[startsWith(inputs$path, paste0(raw, "/")), ]
stopifnot(nrow(inputs) == 3L)
verify_input_files(inputs)

rna <- fread(file.path(raw, "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"))
meta <- fread(file.path(raw, "LUAD_meta.txt"))[case_id != "data_type"]
surv <- fread(file.path(raw, "LUAD_survival.txt"))
cases <- names(rna)[-1L]
stopifnot(!anyDuplicated(cases), !anyNA(meta$case_id), !anyDuplicated(meta$case_id),
          !anyNA(surv$case_id), !anyDuplicated(surv$case_id),
          setequal(cases, meta$case_id), setequal(cases, surv$case_id),
          all(is.na(surv$OS_event) | surv$OS_event %in% 0:1))
meta <- meta[match(cases, case_id)]
surv <- surv[match(cases, case_id)]
stopifnot(identical(meta$case_id, cases), identical(surv$case_id, cases))
clin <- data.frame(sample = cases, os_time = surv$OS_days / 30.4375, os_event = surv$OS_event,
                   age = as.numeric(meta$Age), sex = tolower(meta$Sex),
                   stage = factor(sub("^Stage ", "", meta$Stage), levels = c("I", "II", "III", "IV")))
eligible <- with(clin, is.finite(os_time) & os_time > 0 & !is.na(os_event))
clin <- clin[eligible, ]
clin$sex <- factor(clin$sex, levels = c("female", "male"))
complete <- complete.cases(clin[, c("stage", "age", "sex")])
stopifnot(any(complete), all(is.na(clin$age) | clin$age > 0))

mapped <- collapse_gene_expression(as.matrix(rna[, ..cases]), ensembl_symbols(rna$idx),
                                   transport = TRUE)
sets <- programme_sets()
original <- fread(file.path(tab, "02_transport_gene_universe.csv"))$gene_symbol
covered <- intersect(original, rownames(mapped$transport_expr))
missing <- setdiff(original, covered)
coverage <- rbindlist(list(
  gene_set_coverage(rownames(mapped$expr), sets, id, "cohort"),
  gene_set_coverage(covered, sets, id, "rna_seq_transport")
))
check_gene_set_coverage(coverage)

expr <- mapped$expr[, clin$sample, drop = FALSE]
score <- programme_scores(ssgsea_fast(expr, sets))$scores
d <- cbind(clin, score)
d$stage_group <- factor(ifelse(d$stage == "I", "I", "II–IV"), levels = c("I", "II–IV"))
d$age_z <- z(d$age)
association <- rbindlist(lapply(names(PROGRAMMES), function(axis)
  cox_row(d, axis, c("stage_group", "age_z", "sex"), "stage_I_vs_II_IV_age_sex", id)))

tcga <- readRDS("data/processed/transport_components.rds")[["TCGA-LUAD"]]
tcga_expr <- readRDS("data/processed/transport_expression.rds")[["TCGA-LUAD"]]
if (length(missing)) {
  tcga <- transport_components(tcga_expr, tcga$clin, covered, sets)
}
test <- transport_components(mapped$transport_expr[, clin$sample, drop = FALSE], clin, covered, sets)
frozen <- fit_frozen_transport(tcga)
set.seed(20261124)
performance <- evaluate_frozen_transport(frozen, test, id)

followup <- survfit(Surv(os_time, 1L - os_event) ~ 1, data = clin)
audit <- data.frame(
  cohort = id, source_cases = length(cases), with_os = nrow(clin), excluded_invalid_os = sum(!eligible),
  deaths = sum(clin$os_event), complete_cases = sum(complete), complete_case_deaths = sum(clin$os_event[complete]),
  stage_known = sum(!is.na(clin$stage)), age_known = sum(!is.na(clin$age)), sex_known = sum(!is.na(clin$sex)),
  median_followup_months = unname(summary(followup)$table["median"]),
  original_genes_required = length(original), original_genes_measured = length(covered),
  missing_symbols = paste(missing, collapse = ";"),
  transport_reference = if (length(missing)) "TCGA_fit_on_compatible_gene_subset" else "original_TCGA_fit"
)
outputs <- list(audit = audit, coverage = coverage, association = association,
                transport = performance$evaluation, marker_contrasts = performance$contrasts,
                calibration = performance$calibration, frozen_coefficients = transport_coefficients(frozen),
                score_reference = transport_score_reference(frozen), recalibration = performance$recalibration)
for (name in names(outputs))
  fwrite(outputs[[name]], file.path(tab, paste0("05b_rnaseq_", name, ".csv")))
print(audit)
print(association[, .(exposure, n, events, HR, lo, hi)])
print(performance$evaluation[, .(model, c_index, delta_c, delta_lo, delta_hi, brier_36)])

source("R/lib_external_synthesis.R")
write_external_synthesis(readRDS("data/processed/scores.rds"), tab, cptac = d)

source("R/lib_published_signatures.R")
signatures <- read_published_signatures()
published <- readRDS(PUBLISHED_ASSOCIATION_CACHE)
published[[id]] <- published_association_cohort(expr, d, id, signatures)
write_published_benchmark(published, tab)

path <- file.path(tab, "05d_published_signature_transport.csv")
arrays <- fread(path)[cohort != id]
components <- list(`TCGA-LUAD` = tcga, `CPTAC-LUAD` = test)
expression <- list(`TCGA-LUAD` = tcga_expr, `CPTAC-LUAD` = mapped$transport_expr)
cptac_published <- published_transport_results(components, expression, covered, signatures)
fwrite(rbindlist(list(arrays, cptac_published), use.names = TRUE), path)
