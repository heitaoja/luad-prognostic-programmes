# =============================================================================
# 99_verify.R  --  self-check: are all deliverables present and are the headline
#                  numbers still what the README claims?
# -----------------------------------------------------------------------------
# Run after ./run_all.sh:   Rscript R/99_verify.R
# Exit status 0 = all checks passed, 1 = something is missing or has drifted.
# =============================================================================
PROJ <- normalizePath(".", winslash = "/")
TAB  <- file.path(PROJ, "results", "tables")
FIG  <- file.path(PROJ, "results", "figures")
DOC  <- file.path(PROJ, "docs")

fails <- character(0)
check <- function(ok, label, detail = "") {
  cat(sprintf("%s  %-52s %s\n", if (ok) "PASS" else "FAIL", label, detail))
  if (!ok) fails <<- c(fails, label)
  invisible(ok)
}
near <- function(x, target, tol) abs(x - target) <= tol

# ---- 1. files ---------------------------------------------------------------
need_tab <- c("00_params.csv", "01_cohort_audit.csv", "02_axis_definition.csv",
              "03_cox_adjusted.csv", "03_cindex_ladder.csv", "03_stage_stratified.csv",
              "03b_mutual_adjustment.csv", "04_cohort_estimates.csv",
              "04_meta_analysis.csv", "05_incremental_value.csv", "05_dca_curves.csv",
              "06_random_set_null.csv", "06_permutation_test.csv", "06_sensitivity.csv",
              "table1_cohorts.csv", "table3_meta_analysis.csv", "project_tables.xlsx")
for (f in need_tab)
  check(file.exists(file.path(TAB, f)) && file.size(file.path(TAB, f)) > 200,
        paste0("results/tables/", f))

need_fig <- c("fig1_cohort_overview", "fig2_axis_distribution", "fig3_forest_all_cohorts",
              "fig4_km_tcga", "fig5_stage_stratified", "fig6_incremental",
              "fig7_negative_control", "fig8_sensitivity")
for (f in need_fig)
  check(all(file.exists(file.path(FIG, paste0(f, c(".png", ".pdf"))))),
        paste0("results/figures/", f, " (.png + .pdf)"))

for (f in c("analysis_plan.md", "methods.md", "limitations.md", "data_provenance.md",
            "changelog.md", "sessionInfo.txt", "slides_LUAD_prognostic_programmes.pptx"))
  check(file.exists(file.path(DOC, f)), paste0("docs/", f))

check(dir.exists(file.path(PROJ, "data", "processed")) &&
      file.exists(file.path(PROJ, "data", "processed", "cohorts.rds")),
      "data/processed/cohorts.rds")

# ---- 2. cohort sizes --------------------------------------------------------
aud <- utils::read.csv(file.path(TAB, "01_cohort_audit.csv"))
get_n <- function(co) sum(aud$n_samples[aud$cohort %in% co])
check(get_n("TCGA-LUAD") == 515, "TCGA-LUAD has 515 samples", get_n("TCGA-LUAD"))
check(sum(get_n(c("GSE31210", "GSE37745", "GSE50081"))) == 437,
      "validation cohorts total 437 samples",
      sum(get_n(c("GSE31210", "GSE37745", "GSE50081"))))
check(sum(aud$n_events) == 343, "total deaths across cohorts = 343", sum(aud$n_events))

# ---- 3. headline numbers ----------------------------------------------------
meta <- utils::read.csv(file.path(TAB, "04_meta_analysis.csv"))
val <- function(a) meta$pooled_HR[meta$axis == a & meta$scope == "validation"]
i2  <- function(a) meta$I2[meta$axis == a & meta$scope == "validation"]
check(near(val("proliferation"), 1.317, 0.02), "proliferation pooled HR ~ 1.32",
      round(val("proliferation"), 3))
check(i2("proliferation") == 0, "proliferation I2 = 0%", i2("proliferation"))
check(val("immune") < 1.25, "immune pooled HR not protective/large",
      round(val("immune"), 3))

inc <- utils::read.csv(file.path(TAB, "05_incremental_value.csv"))
dc <- function(co, a) inc$delta_c[inc$cohort == co & inc$axis == a]
check(near(dc("TCGA-LUAD", "proliferation"), 0.019, 0.005),
      "TCGA delta C-index (proliferation) ~ 0.019", round(dc("TCGA-LUAD", "proliferation"), 3))
check(near(dc("TCGA-LUAD", "hypoxia"), 0.033, 0.005),
      "TCGA delta C-index (hypoxia) ~ 0.033", round(dc("TCGA-LUAD", "hypoxia"), 3))

perm <- utils::read.csv(file.path(TAB, "06_permutation_test.csv"))
check(all(perm$empirical_p < 0.01), "permutation p < 0.01 for both axes",
      paste(signif(perm$empirical_p, 2), collapse = " / "))

nullt <- utils::read.csv(file.path(TAB, "06_random_set_null.csv"))
t0 <- nullt[nullt$cohort == "TCGA-LUAD", ]
check(all(t0$observed_HR > t0$null_q99), "observed HR above the 99th percentile of random sets",
      paste(sprintf("%.2f>%.2f", t0$observed_HR, t0$null_q99), collapse = " "))

sens <- utils::read.csv(file.path(TAB, "06_sensitivity.csv"))
check(all(tapply(sens$HR > 1, sens$axis, all)),
      "all sensitivity variants keep HR > 1")

est <- utils::read.csv(file.path(TAB, "04_cohort_estimates.csv"))
check(nrow(est) == 16, "16 cohort x axis estimates (4 axes x 4 cohorts)", nrow(est))

# ---- 4. summary -------------------------------------------------------------
cat("\n")
if (length(fails)) {
  cat("FAILED:", length(fails), "check(s)\n")
  for (f in fails) cat("  -", f, "\n")
  quit(status = 1)
}
cat("ALL CHECKS PASSED\n")
