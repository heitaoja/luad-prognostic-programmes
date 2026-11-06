# =============================================================================
# 08_tables.R  --  assemble the manuscript tables (CSV + one Excel workbook)
# -----------------------------------------------------------------------------
# Output: results/tables/table1_cohorts.csv ... table5_sensitivity.csv
#         results/tables/project_tables.xlsx   (all sheets, formatted)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
PROJ <- normalizePath(".", winslash = "/")
TAB  <- file.path(PROJ, "results", "tables")

aud  <- fread(file.path(TAB, "01_cohort_audit.csv"))
est  <- fread(file.path(TAB, "04_cohort_estimates.csv"))
meta <- fread(file.path(TAB, "04_meta_analysis.csv"))
inc  <- fread(file.path(TAB, "05_incremental_value.csv"))
sens <- fread(file.path(TAB, "06_sensitivity.csv"))
null <- fread(file.path(TAB, "06_random_set_null.csv"))
perm <- fread(file.path(TAB, "06_permutation_test.csv"))
strat<- fread(file.path(TAB, "03_stage_stratified.csv"))

AXIS_LAB <- c(proliferation = "Proliferation", hypoxia = "Hypoxia/glycolysis",
              emt_stroma = "EMT/stroma", immune = "Immune infiltration")

# ---- Table 1: cohort characteristics ---------------------------------------
t1 <- aud[, .(Cohort = cohort, Samples = n_samples,
              `Deaths` = n_events,
              `Median OS (months)` = median_os_mo,
              `Stage known (%)` = pct_stage_known,
              `Age known` = n_age, `Sex known` = n_sex,
              `Raw time unit` = raw_time_unit,
              `Expression range (log2)` = expr_range)]
fwrite(t1, file.path(TAB, "table1_cohorts.csv"))

# ---- Table 2: adjusted HR by cohort ----------------------------------------
t2 <- dcast(est[, .(axis = exposure, cohort, HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi), p = p)],
            axis ~ cohort, value.var = c("HR_CI", "p"))
fwrite(t2, file.path(TAB, "table2_cox_by_cohort.csv"))

# ---- Table 3: meta-analysis -------------------------------------------------
t3 <- dcast(meta[, .(axis, scope, HR_CI = sprintf("%.2f (%.2f-%.2f)", pooled_HR, lo, hi))],
            axis ~ scope, value.var = "HR_CI")
setnames(t3, "all", "all_cohorts_pooled")
setnames(t3, "validation", "validation_pooled")
t3 <- merge(t3, meta[scope == "validation", .(axis, k_validation = k, I2_percent = I2,
                                              tau2, p)], by = "axis")
t3 <- t3[, .(axis, validation_pooled, all_cohorts_pooled, k_validation, I2_percent, tau2, p)]
t3[, axis := AXIS_LAB[axis]]
fwrite(t3, file.path(TAB, "table3_meta_analysis.csv"))

# ---- Table 4: incremental value ---------------------------------------------
t4 <- dcast(inc[, .(axis, cohort, dC = sprintf("%+.3f (%.3f, %.3f)", delta_c, delta_lo, delta_hi),
                    lr_p = signif(lr_p, 2), nri = round(nri, 2))],
            axis ~ cohort, value.var = "dC")
t4b <- inc[, .(axis, cohort, lr_p = signif(lr_p, 3), NRI = round(nri, 2))]
fwrite(t4, file.path(TAB, "table4_delta_cindex.csv"))
fwrite(t4b, file.path(TAB, "table4b_lr_and_nri.csv"))

# ---- Table 5: sensitivity ----------------------------------------------------
t5 <- dcast(sens, axis + variant ~ cohort, value.var = "HR")
fwrite(t5, file.path(TAB, "table5_sensitivity.csv"))

# ---- Table 6: negative controls ---------------------------------------------
t6 <- merge(null[, .(cohort, axis, B, observed_HR, null_q95, null_q99, empirical_p)],
            perm[, .(axis, perm_empirical_p = empirical_p)], by = "axis", all.x = TRUE)
fwrite(t6, file.path(TAB, "table6_negative_controls.csv"))

# ---- Table 7: stage-stratified ----------------------------------------------
t7 <- dcast(strat[, .(HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi)), by = .(exposure, stage)],
            exposure ~ stage, value.var = "HR_CI")
fwrite(t7, file.path(TAB, "table7_stage_stratified.csv"))

# ---- Excel workbook ----------------------------------------------------------
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", border = "bottom")
add_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
add_sheet("T1 cohorts", t1)
add_sheet("T2 HR by cohort", t2)
add_sheet("T3 meta", t3)
add_sheet("T4 delta C-index", t4)
add_sheet("T4b LR and NRI", t4b)
add_sheet("T5 sensitivity", t5)
add_sheet("T6 negative controls", t6)
add_sheet("T7 by stage", t7)
saveWorkbook(wb, file.path(TAB, "project_tables.xlsx"), overwrite = TRUE)

cat("tables written:\n")
for (f in c("table1_cohorts", "table2_cox_by_cohort", "table3_meta_analysis",
            "table4_delta_cindex", "table5_sensitivity", "table6_negative_controls",
            "table7_stage_stratified")) cat("  ", f, ".csv\n")
cat("   project_tables.xlsx\n\n")
cat("=== T1 ===\n"); print(t1, row.names = FALSE)
cat("\n=== T3 ===\n"); print(t3, row.names = FALSE)
cat("\nDONE 08\n")
