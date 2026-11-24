# Check scientific invariants and the release deliverables.
source("R/lib_inputs.R")
source("R/lib_score.R")
source("R/lib_surv.R")
source("R/lib_release_tables.R")
verify_input_files()
sets <- programme_sets()
stopifnot(length(sets) == 14L, length(sets$E2F_TARGETS) == 200L,
          "AURKA" %in% sets$E2F_TARGETS, "PGK1" %in% sets$HYPOXIA)

# Regression coverage for empty GMT descriptions and order-invariant ties.
example <- matrix(c(3, 3, 2, 2, 1, 1, 0, 0), ncol = 1L,
                  dimnames = list(LETTERS[1:8], "sample"))
gene_sets <- list(example = c("A", "C", "E", "G", "H"))
direct <- ssgsea_fast(example, gene_sets)
stopifnot(isTRUE(all.equal(direct, ssgsea_fast(example[8:1, , drop = FALSE], gene_sets),
                           tolerance = 1e-12)),
          isTRUE(all.equal(direct, ssgsea_fast(log2(example + 1), gene_sets), tolerance = 1e-12)))
m <- meta_reml_mkh(rep(log(1.3), 3), rep(0.13, 3))
stopifnot(m$tau2 == 0, isTRUE(all.equal(m$pooled_re, 1.3)),
          isTRUE(all.equal(log(m$hi / m$pooled_re), qt(0.975, 2) * 0.13 / sqrt(3))))

full <- Sys.getenv("LUAD_CORE_ONLY") != "1"
paths <- c("01_cohort_audit", "01_tcga_participant_flow", "02_axis_definition", "02_gene_coverage",
           "02_paired_specimen_scores", "02_transport_gene_universe", "02_transport_coverage",
           "03c_marker_benchmark", "03c_marker_synthesis", "03c_marker_correlations", "03c_model_diagnostics",
           "03c_specimen_sensitivity", "03c_purity_sensitivity", "03d_driver_smoking_adjustment",
           "03e_published_signature_benchmark", "04_cohort_estimates",
           "04_meta_analysis", "04_leave_one_out", "04b_ipd_and_pan_external_synthesis",
           "05_incremental_value", "05_frozen_coefficients",
           "05_frozen_score_reference", "05_simple_marker_contrasts", "05_calibration", "05c_recalibration",
           "05d_published_signature_transport",
           "06_random_set_null", "06_permutation_test", "06_sensitivity",
           "table1_cohorts", "table3_meta_analysis", "table4_delta_cindex", "table8_marker_benchmark")
if (full) {
  paths <- c(paths, paste0("05b_rnaseq_", c("audit", "coverage", "association", "transport",
             "marker_contrasts", "calibration", "frozen_coefficients", "score_reference", "recalibration")))
}
read_table <- function(x) fread(file.path("results/tables", paste0(x, ".csv")))
verify_files <- function(files) {
  if (any(!file.exists(files))) stop("Missing outputs: ", paste(files[!file.exists(files)], collapse = ", "))
  if (any(file.size(files) <= 0)) stop("Empty output file")
}
verify_files(file.path("results/tables", paste0(paths, ".csv")))
verify_files(c("results/tables/project_tables.xlsx", "README.md",
               file.path("docs", c("manuscript.md", "methods.md", "analysis_plan_v2.md",
                                   "limitations.md", "reporting_crosswalk.md", "changelog.md"))))
verify_files("config/published_signatures.csv")
main_figures <- c("Figure1_Cohort_and_Scores", "Figure2_Association_and_IPD_Synthesis",
                  "Figure3_Specificity_Purity_and_Drivers", "Figure4_Redundancy_Benchmark",
                  "Figure5_Calibration_and_Recalibration")
supp_figures <- c("FigureS2_Cohort_Overview", "FigureS3_Score_Distributions",
                  "FigureS4_Array_Forest_Plots", "FigureS5_TCGA_Kaplan_Meier",
                  "FigureS6_TCGA_Stage_Stratified", "FigureS7_External_Harrell_C",
                  "FigureS8_Random_Programme_Null", "FigureS9_Sensitivity_Heatmap",
                  "FigureS10_Simple_Marker_Benchmark", "FigureS11_Frozen_Calibration")
if (full) supp_figures <- c("FigureS1_CPTAC_RNAseq_Validation", supp_figures, "FigureS12_Purity_Sensitivity")
verify_files(file.path("results/figures/main", unlist(lapply(main_figures, paste0, c(".pdf", ".png")))))
verify_files(file.path("results/figures/supplementary", unlist(lapply(supp_figures, paste0, c(".pdf", ".png")))))

workbook_path <- "results/tables/project_tables.xlsx"
sheet_names <- openxlsx::getSheetNames(workbook_path)
evidence <- release_evidence_sheets(full)
core_sheets <- if (full) 24L else 19L
stopifnot(length(sheet_names) == core_sheets + length(evidence),
          !anyDuplicated(sheet_names), all(nchar(sheet_names) <= 31L),
          identical(tail(sheet_names, length(evidence)), names(evidence)))
for (sheet in names(evidence)) {
  csv <- read_table(evidence[[sheet]])
  actual <- openxlsx::read.xlsx(workbook_path, sheet = sheet, skipEmptyRows = FALSE,
                                skipEmptyCols = FALSE, check.names = FALSE)
  stopifnot(nrow(actual) == nrow(csv), identical(names(actual), names(csv)))
  for (column in names(csv)) {
    expected <- csv[[column]]
    observed <- actual[[column]]
    if (all(is.na(expected))) {
      stopifnot(all(is.na(observed)))
      next
    }
    if (is.character(expected)) {
      # Blank CSV fields and blank workbook cells both represent missing text.
      expected[is.na(expected)] <- ""
      observed[is.na(observed)] <- ""
    }
    if (!isTRUE(all.equal(observed, expected, tolerance = 1e-12))) {
      stop("Workbook sheet differs from ", evidence[[sheet]], ".csv: ", sheet, " / ", column, call. = FALSE)
    }
  }
}

aud <- read_table("01_cohort_audit")
flow <- read_table("01_tcga_participant_flow")
p <- readRDS("data/processed/scores.rds")
est <- read_table("04_cohort_estimates")
expected_cohorts <- if (full) 5L else 4L
stopifnot(length(p) == expected_cohorts,
          all(vapply(p, function(d) !anyDuplicated(d$sample), logical(1L))),
          !anyDuplicated(p[["TCGA-LUAD"]]$patient),
          nrow(p[["TCGA-LUAD"]]) == flow$participants_with_os,
          flow$specimens_with_os > flow$participants_with_os,
          sum(est$role == "discovery") == 4L,
          nrow(est) == 4L * nrow(aud),
          all(est$n == aud$n_complete[match(est$cohort, aud$cohort)]),
          all(est$events == aud$events_complete[match(est$cohort, aud$cohort)]))
coverage <- read_table("02_gene_coverage")
stopifnot(all(coverage$n_found / coverage$n_in_set >= .8))
meta <- read_table("04_meta_analysis")
stopifnot(all(meta[scope == "validation", k] == 3L))
if (full) {
  stopifnot(nrow(meta[scope == "extended"]) == 8L, all(meta[scope == "extended", k] == 4L),
            all(meta[scope == "extended", n_total] == sum(aud[cohort != "TCGA-LUAD", n_complete])))
  purity <- read_table("03c_purity_sensitivity")
  stopifnot(nrow(purity) == 14L, length(unique(purity$n)) == 1L,
            length(unique(purity$events)) == 1L, all(purity$n >= 50L))
  cptac <- read_table("05b_rnaseq_audit")
  cptac_coverage <- read_table("05b_rnaseq_coverage")
  cptac_association <- read_table("05b_rnaseq_association")
  cptac_transport <- read_table("05b_rnaseq_transport")
  required <- read_table("02_transport_gene_universe")$gene_symbol
  omitted <- unlist(strsplit(na.omit(as.character(cptac$missing_symbols)), ";", fixed = TRUE),
                    use.names = FALSE)
  stopifnot(nrow(cptac) == 1L, cptac$source_cases == cptac$with_os + cptac$excluded_invalid_os,
            cptac$original_genes_required == length(required), all(omitted %in% required),
            cptac$original_genes_measured == length(required) - length(omitted),
            nrow(cptac_coverage) == 28L, all(cptac_coverage$n_found / cptac_coverage$n_in_set >= .8),
            identical(cptac_association$exposure, names(PROGRAMMES)),
            all(cptac_association$n == cptac$complete_cases),
            all(cptac_association$events == cptac$complete_case_deaths),
            nrow(cptac_transport) == 5L, all(cptac_transport$n == cptac$complete_cases),
            all(cptac_transport$events == cptac$complete_case_deaths),
            all(cptac_transport$c_index >= 0 & cptac_transport$c_index <= 1),
            nrow(read_table("05b_rnaseq_marker_contrasts")) == 3L,
            !"CPTAC-LUAD" %in% est$cohort)
}
transport <- read_table("05_incremental_value")
stopifnot(nrow(transport) == 5L * (nrow(aud) - 1L), !"TCGA-LUAD" %in% transport$cohort,
          all(transport[model == "clinical", delta_c] == 0),
          all(transport$c_index >= 0 & transport$c_index <= 1),
          all(transport$predicted_36 >= 0 & transport$predicted_36 <= 1))
stopifnot(all(read_table("06_random_set_null")$empirical_p > 0))

valid_effects <- function(rows) {
  all(is.finite(rows$HR) & rows$HR > 0 & is.finite(rows$lo) & rows$lo > 0 &
        rows$lo <= rows$HR & is.finite(rows$hi) & rows$hi >= rows$HR &
        is.finite(rows$p) & rows$p >= 0 & rows$p <= 1 &
        rows$events > 0L & rows$events <= rows$n)
}
driver_smoking <- read_table("03d_driver_smoking_adjustment")
stopifnot(identical(names(driver_smoking), c("cohort", "analysis_type", "subgroup_or_model", "exposure",
                                            "covariates", "n", "events", "HR", "lo", "hi", "p")),
          nrow(driver_smoking) == (if (full) 60L else 28L),
          valid_effects(driver_smoking),
          !anyDuplicated(driver_smoking[, .(cohort, analysis_type, subgroup_or_model, exposure)]))
adjustment_pairs <- driver_smoking[analysis_type != "driver_subgroup",
  .(matched = .N == 2L && setequal(subgroup_or_model, c("base_clinical", "extended_adjusted")) &&
      uniqueN(n) == 1L && uniqueN(events) == 1L), by = .(cohort, analysis_type, exposure)]
stopifnot(nrow(adjustment_pairs) == (if (full) 20L else 12L),
          all(adjustment_pairs$matched), all(adjustment_pairs$exposure %in% names(PROGRAMMES)))

score_sensitivity <- read_table("06_sensitivity")[variant %in% c(
  "Mean gene-wise z-score", "Mean gene-wise z-score (sample-centred)")]
score_pairs <- score_sensitivity[, .(matched = .N == 2L && uniqueN(variant) == 2L &&
  uniqueN(n) == 1L && uniqueN(events) == 1L), by = .(cohort, axis)]
stopifnot(nrow(score_pairs) == 2L * expected_cohorts, all(score_pairs$matched),
          setequal(score_pairs$cohort, names(p)), all(score_pairs$axis %in% c("proliferation", "hypoxia")),
          valid_effects(score_sensitivity))

synthesis <- read_table("04b_ipd_and_pan_external_synthesis")
expected_scopes <- data.table(scope = c("validation_array_k3", "extended_array_k4", "all_external_k5"),
                              k = c(3L, 4L, 5L), n_total = c(437L, 830L, 935L),
                              events_total = c(158L, 269L, 292L))
if (!full) expected_scopes <- expected_scopes[1L]
ipd_methods <- c("IPD_strata_cohort", "IPD_strata_cohort_stage")
meta_methods <- c("REML_mKH", "DL_normal")
expected_rows <- expected_scopes[, {
  methods <- if (.BY$k == 5L) c(ipd_methods, meta_methods) else ipd_methods
  CJ(method = methods, axis = names(PROGRAMMES))
}, by = .(scope, k, n_total, events_total)]
keys <- names(expected_rows)
stopifnot(nrow(synthesis) == nrow(expected_rows), fsetequal(synthesis[, ..keys], expected_rows),
          valid_effects(synthesis[, .(HR = pooled_HR, lo, hi, p, n = n_total, events = events_total)]),
          all(synthesis[method == "REML_mKH", df] == 4L),
          all(is.na(synthesis[method != "REML_mKH", df])),
          all(is.finite(synthesis$p_holm) & synthesis$p_holm >= 0 & synthesis$p_holm <= 1))
holm <- synthesis[, .(matches = all(abs(p_holm - p.adjust(p, method = "holm")) <= 1e-12)),
                  by = .(scope, method)]
pan_meta <- synthesis[method %in% meta_methods]
stopifnot(all(holm$matches), all(is.finite(pan_meta$I2) & pan_meta$I2 >= 0 & pan_meta$I2 <= 100),
          all(is.finite(pan_meta$tau2) & pan_meta$tau2 >= 0), all(is.finite(pan_meta$Q) & pan_meta$Q >= 0))

verify_transport_extension <- function(evaluation, calibration, recalibration) {
  recal_fields <- c("recal_H0_36", "recal_predicted_36", "recal_brier_36", "recal_brier_gain_36")
  appended <- c("uno_c", "delta_uno_c", "delta_uno_lo", "delta_uno_hi", recal_fields)
  probability <- function(x) all(is.finite(x) & x >= 0 & x <= 1)
  clinical <- evaluation[model == "clinical"]
  baseline <- match(evaluation$cohort, clinical$cohort)
  models <- c("clinical", "composite", "MKI67", "E2F", "G2M")
  model_checks <- evaluation[, .(complete = setequal(model, models)), by = cohort]
  stopifnot(identical(tail(names(evaluation), length(appended)), appended),
            !anyDuplicated(evaluation[, .(cohort, model)]), all(model_checks$complete),
            probability(evaluation$uno_c),
            all(abs(evaluation$delta_uno_c - (evaluation$uno_c - clinical$uno_c[baseline])) <= 1e-12),
            all(evaluation$delta_uno_lo >= -1 & evaluation$delta_uno_hi <= 1 &
                  evaluation$delta_uno_lo <= evaluation$delta_uno_hi),
            all(clinical$delta_uno_lo == 0 & clinical$delta_uno_hi == 0))
  eligible <- evaluation$events_36 >= 10L & evaluation$n_at_risk_36 >= 10L
  active <- evaluation[eligible]
  gain <- clinical$recal_brier_36[baseline[eligible]] - active$recal_brier_36
  stopifnot(all(is.na(evaluation[!eligible, ..recal_fields])),
            all(is.finite(active$recal_H0_36) & active$recal_H0_36 > 0),
            probability(active$recal_predicted_36), probability(active$recal_brier_36),
            all(abs(active$recal_brier_gain_36 - gain) <= 1e-12))

  fields <- names(calibration)
  shared <- setdiff(fields, "predicted")
  frozen <- recalibration[calibration_type == "frozen", ..fields]
  updated <- recalibration[calibration_type == "baseline_recalibrated", ..fields]
  stopifnot(identical(names(recalibration), c(fields, "calibration_type")),
            nrow(recalibration) == 2L * nrow(calibration), identical(frozen, calibration),
            identical(updated[, ..shared], calibration[, ..shared]),
            probability(updated$predicted),
            all(updated$n > 0L & updated$at_risk >= 0L & updated$at_risk <= updated$n))
  at_horizon <- updated[at_risk > 0L]
  stopifnot(probability(at_horizon$observed), probability(at_horizon$observed_lo),
            probability(at_horizon$observed_hi),
            all(at_horizon$observed_lo <= at_horizon$observed & at_horizon$observed <= at_horizon$observed_hi))
  totals <- updated[, .(quartiles = .N == 4L && setequal(group, 1:4),
                        n = sum(n), at_risk = sum(at_risk), predicted = weighted.mean(predicted, n)),
                    by = .(cohort, model)]
  expected <- active[model %in% c("clinical", "composite")]
  stopifnot(all(totals$quartiles),
            fsetequal(totals[, .(cohort, model, n, at_risk)], expected[, .(cohort, model, n, at_risk = n_at_risk_36)]),
            all(abs(totals$predicted - expected[totals, on = .(cohort, model), recal_predicted_36]) <= 1e-12))
}
verify_transport_extension(transport, read_table("05_calibration"), read_table("05c_recalibration"))
if (full) {
  verify_transport_extension(cptac_transport, read_table("05b_rnaseq_calibration"),
                             read_table("05b_rnaseq_recalibration"))
}
gse31210 <- transport[cohort == "GSE31210"]
stopifnot(nrow(gse31210) == 5L, all(gse31210$recal_brier_36 < gse31210$brier_36),
          all(abs(gse31210$recal_predicted_36 - gse31210$observed_36) < 0.01),
          all(abs(gse31210$recal_predicted_36 - gse31210$observed_36) <
                abs(gse31210$predicted_36 - gse31210$observed_36)))

source("R/lib_published_signatures.R")
signatures <- read_published_signatures()
stopifnot(all(c("CDC2", "CENPM", "SKA1") %in% signatures$CCP_31), "PCNA" %in% signatures$meta_PCNA,
          identical(match_published_signatures(c("SKA3", "SKA1"), list(test = "C18orf24"))$test, "SKA1"))
verify_published_coverage <- function(rows) {
  minimum <- ifelse(rows$signature_id == "meta_PCNA", 0.85, 0.8)
  stopifnot(all(rows$n_in_set == lengths(signatures)[rows$signature_id]),
            all(rows$n_found <= rows$n_in_set), all(rows$coverage >= minimum & rows$coverage <= 1),
            all(abs(rows$coverage - rows$n_found / rows$n_in_set) <= 1e-12))
}
benchmark <- read_table("03e_published_signature_benchmark")
model_keys <- rbindlist(lapply(names(signatures), function(id)
  data.table(signature_id = id,
             model = c("clinical_plus_marker", paste0("composite_given_", id), paste0(id, "_given_composite")),
             exposure = c(id, "proliferation", id))))
model_keys <- model_keys[, .(scoring = PUBLISHED_ASSOCIATION_SCORES), by = .(signature_id, model, exposure)]
key <- names(model_keys)
cohort_rows <- benchmark[scope == "cohort"]
pooled_rows <- benchmark[scope != "cohort"]
cohort_counts <- aud[, .(cohort, n = n_complete, events = events_complete)]
if (full) {
  cohort_counts <- rbind(cohort_counts, cptac[, .(cohort, n = complete_cases, events = complete_case_deaths)])
}
pool_counts <- expected_scopes[, .(scope = c("validation", "extended", "all_external_k5")[k - 2L],
                                    k, n = n_total, events = events_total)]
pool_counts <- pool_counts[, .(method = c("REML_mKH", "DL_normal", "IPD_strata_cohort")),
                            by = .(scope, k, n, events)]
groups <- benchmark[, .(complete = .N == nrow(model_keys) && fsetequal(.SD, model_keys)),
                     by = .(cohort, scope, method), .SDcols = key]
stopifnot(valid_effects(benchmark), all(groups$complete), all(cohort_rows$method == "Cox" & cohort_rows$k == 1L),
          fsetequal(unique(cohort_rows[, .(cohort, n, events)]), cohort_counts),
          fsetequal(unique(pooled_rows[, .(scope, k, n, events, method)]), pool_counts),
          nrow(benchmark) == nrow(model_keys) * (nrow(cohort_counts) + nrow(pool_counts)),
          all(is.finite(cohort_rows$spearman) & abs(cohort_rows$spearman) <= 1))
verify_published_coverage(cohort_rows)

published_transport <- read_table("05d_published_signature_transport")
external <- cohort_counts[cohort != "TCGA-LUAD"]
expected <- CJ(cohort = external$cohort, signature_id = names(signatures), scoring = PUBLISHED_TRANSPORT_SCORES)
key <- names(expected)
verify_published_coverage(published_transport)
sharing <- nrow(read_table("02_transport_gene_universe"))
universes <- rep(sharing, nrow(published_transport))
original_transport <- transport
if (full) {
  stopifnot(sharing == 17595L, cptac$original_genes_measured == 17593L)
  universes[published_transport$cohort == "CPTAC-LUAD"] <- cptac$original_genes_measured
  original_transport <- rbind(transport, cptac_transport)
}
clinical <- original_transport[model == "clinical"][published_transport, on = "cohort", .(c_index, uno_c)]
composite <- original_transport[model == "composite"][published_transport, on = "cohort", c_index]
probabilities <- published_transport[, .(c_index, uno_c, brier_36, recal_brier_36)]
stopifnot(nrow(published_transport) == nrow(expected), fsetequal(published_transport[, ..key], expected),
          all(published_transport$model == paste(published_transport$signature_id, published_transport$scoring, sep = "_")),
          all(published_transport$n_genes_universe == universes),
          fsetequal(unique(published_transport[, .(cohort, n, events)]), external),
          all(vapply(probabilities, function(x) all(is.finite(x) & x >= 0 & x <= 1), logical(1))),
          all(abs(published_transport$clinical_c_index - clinical$c_index) <= 1e-12),
          all(abs(published_transport$composite_c_index - composite) <= 1e-12),
          all(abs(published_transport$delta_c - (published_transport$c_index - clinical$c_index)) <= 1e-12),
          all(abs(published_transport$delta_uno_c - (published_transport$uno_c - clinical$uno_c)) <= 1e-12),
          all(abs(published_transport$composite_delta_c - (composite - published_transport$c_index)) <= 1e-12),
          all(published_transport$composite_delta_lo >= -1 & published_transport$composite_delta_hi <= 1 &
                published_transport$composite_delta_lo <= published_transport$composite_delta_hi))
cat("PASS: pinned inputs, scoring, synthesis, participant flow, model transport and release outputs\n")
