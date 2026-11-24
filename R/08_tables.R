# Assemble rounded reporting tables and the full evidence workbook.
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
tab <- "results/tables"
read_table <- function(name) fread(file.path(tab, paste0(name, ".csv")))
write_table <- function(d, name) fwrite(d, file.path(tab, paste0(name, ".csv")))
ci <- function(x, lo, hi, digits = 2L) {
  form <- paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)")
  sprintf(form, x, lo, hi)
}

aud <- read_table("01_cohort_audit")
t1 <- aud[, .(cohort, participants = n_samples, deaths = n_events, complete_cases = n_complete,
              deaths_complete = events_complete, median_os_mo, median_followup_mo,
              events_36m, at_risk_36m, stage_known = n_stage, age_known = n_age, sex_known = n_sex)]
t2 <- dcast(read_table("04_cohort_estimates")[, .(programme = exposure, cohort,
                                                   adjusted_HR = ci(HR, lo, hi))],
            programme ~ cohort, value.var = "adjusted_HR")
t3 <- read_table("04_meta_analysis")[, .(programme = axis, scope, method, k,
                  participants = n_total, deaths = events_total, HR_95_CI = ci(pooled_HR, lo, hi),
                  p = signif(p, 4), p_holm = signif(p_holm, 4), I2_pct = round(I2, 1), tau2)]
t4 <- read_table("05_incremental_value")[, .(cohort, model, n, events,
                  Harrell_C = round(c_index, 4), delta_C_95_CI = ci(delta_c, delta_lo, delta_hi, 3L),
                  predicted_36 = round(predicted_36, 3), observed_36 = round(observed_36, 3),
                  calibration_slope = ci(calibration_slope, slope_lo, slope_hi),
                  Brier_36 = round(brier_36, 4))]
t5 <- read_table("06_sensitivity")[, .(cohort, axis, variant, n, events, HR_95_CI = ci(HR, lo, hi), p)]
t6 <- read_table("06_random_set_null")[, .(cohort, axis, B, observed_HR = round(observed_HR, 3),
                                      null_q99 = round(null_q99, 3), empirical_p)]
t7 <- dcast(read_table("03_stage_stratified")[, .(axis = exposure, stage, HR_95_CI = ci(HR, lo, hi))],
            axis ~ stage, value.var = "HR_95_CI")
t8 <- read_table("03c_marker_synthesis")[, .(scope, model, exposure, k, HR_95_CI = ci(pooled_HR, lo, hi), p)]
t9 <- read_table("03c_purity_sensitivity")[, .(exposure, model, n, events, HR_95_CI = ci(HR, lo, hi),
                                            spearman = round(correlation_with_purity, 3))]
sheet_names <- c("T1 Cohorts", "T2 Adjusted HR", "T3 External synthesis", "T4 Frozen transport",
                  "T5 Sensitivity", "T6 Matched null", "T7 By stage", "T8 Marker benchmark", "T9 DNA purity")
tables <- setNames(list(t1, t2, t3, t4, t5, t6, t7, t8, t9), sheet_names)
outputs <- c("table1_cohorts", "table2_cox_by_cohort", "table3_meta_analysis", "table4_delta_cindex",
             "table5_sensitivity", "table6_negative_controls", "table7_stage_stratified",
             "table8_marker_benchmark", "table9_purity_sensitivity")
for (i in seq_along(outputs)) write_table(tables[[i]], outputs[i])
tables[["S1 Gene-set coverage"]] <- read_table("02_gene_coverage")
tables[["S2 Leave-one-out"]] <- read_table("04_leave_one_out")
tables[["S3 Model diagnostics"]] <- read_table("03c_model_diagnostics")
tables[["S4 Frozen coefficients"]] <- read_table("05_frozen_coefficients")
tables[["S5 Frozen score reference"]] <- read_table("05_frozen_score_reference")
tables[["S6 Marker contrasts"]] <- read_table("05_simple_marker_contrasts")
tables[["S7 Specimen sensitivity"]] <- read_table("03c_specimen_sensitivity")
tables[["S8 Paired specimens"]] <- read_table("02_paired_specimen_scores")
tables[["S9 Held-out calibration"]] <- read_table("05_calibration")
tables[["S10 Marker correlations"]] <- read_table("03c_marker_correlations")
if (Sys.getenv("LUAD_CORE_ONLY") != "1") {
  tables[["S11 CPTAC cohort"]] <- read_table("05b_rnaseq_audit")
  tables[["S12 CPTAC adjusted HR"]] <- read_table("05b_rnaseq_association")
  tables[["S13 CPTAC transport"]] <- read_table("05b_rnaseq_transport")
  tables[["S14 CPTAC marker contrasts"]] <- read_table("05b_rnaseq_marker_contrasts")
  tables[["S15 CPTAC calibration"]] <- read_table("05b_rnaseq_calibration")
}
wb <- createWorkbook()
header <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom")
for (name in names(tables)) {
  data <- tables[[name]]
  addWorksheet(wb, name)
  writeData(wb, name, data, headerStyle = header, withFilter = nrow(data) > 0L)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}
saveWorkbook(wb, file.path(tab, "project_tables.xlsx"), overwrite = TRUE)
cat("Wrote", length(outputs), "reporting tables and", length(tables), "workbook sheets\n")
