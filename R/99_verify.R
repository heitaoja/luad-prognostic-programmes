# Check scientific invariants and the release deliverables.
source("R/lib_inputs.R")
source("R/lib_score.R")
source("R/lib_surv.R")
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
           "03c_specimen_sensitivity", "03c_purity_sensitivity", "04_cohort_estimates",
           "04_meta_analysis", "04_leave_one_out", "05_incremental_value", "05_frozen_coefficients",
           "05_frozen_score_reference", "05_simple_marker_contrasts", "05_calibration",
           "06_random_set_null", "06_permutation_test", "06_sensitivity",
           "table1_cohorts", "table3_meta_analysis", "table4_delta_cindex", "table8_marker_benchmark")
if (full) {
  paths <- c(paths, paste0("05b_rnaseq_", c("audit", "coverage", "association", "transport",
             "marker_contrasts", "calibration", "frozen_coefficients", "score_reference")))
}
read_table <- function(x) fread(file.path("results/tables", paste0(x, ".csv")))
verify_files <- function(files) {
  if (any(!file.exists(files))) stop("Missing outputs: ", paste(files[!file.exists(files)], collapse = ", "))
  if (any(file.size(files) <= 0)) stop("Empty output file")
}
verify_files(file.path("results/tables", paste0(paths, ".csv")))
verify_files(c("results/tables/project_tables.xlsx", "docs/manuscript.md", "docs/analysis_plan_v2.md"))
figures <- c("fig1_cohort_overview", "fig2_axis_distribution", "fig3_forest_all_cohorts",
             "fig4_km_tcga", "fig5_stage_stratified", "fig6_incremental", "fig7_negative_control",
             "fig8_sensitivity", "fig9_marker_benchmark", "fig10_calibration")
if (full) figures <- c(figures, "fig11_purity_sensitivity", "figS1_rnaseq_validation")
verify_files(file.path("results/figures", unlist(lapply(figures, paste0, c(".pdf", ".png")))))

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
cat("PASS: pinned inputs, scoring, synthesis, participant flow, model transport and release outputs\n")
