# =============================================================================
# 04_external_meta.R  --  external validation and random-effects meta-analysis
# -----------------------------------------------------------------------------
# Same adjusted model as TCGA (stage + age + sex) in each independent cohort,
# then DerSimonian-Laird pooling of the log-HR per axis.
# Output: results/tables/04_*.csv
# =============================================================================
suppressPackageStartupMessages({ library(survival); library(data.table) })
source("R/lib_surv.R")
PROJ <- normalizePath(".", winslash = "/")
OUT  <- file.path(PROJ, "results", "tables")
AXES <- c("proliferation", "hypoxia", "emt_stroma", "immune")
COV  <- c("stage", "age_z", "sex")

s <- readRDS(file.path(PROJ, "data", "processed", "scores.rds"))
COH  <- c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081")
ROLE <- c("discovery", "validation", "validation", "validation")

res <- list()
for (i in seq_along(COH)) {
  id <- COH[i]
  d <- s[[id]]
  d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
  d$age_z <- z(d$age)
  cat(sprintf("\n== %-10s %-10s n=%3d events=%3d\n", id, ROLE[i], nrow(d), sum(d$os_event)))
  for (a in AXES) {
    r <- cox_row(d, a, COV, "adjusted", id)
    r$role <- ROLE[i]
    r$median_fu <- round(median(d$os_time, na.rm = TRUE), 1)
    res[[length(res) + 1]] <- r
  }
}
res <- rbindlist(res)
res[, se := (log(hi) - log(lo)) / (2 * 1.96)]
fwrite(res, file.path(OUT, "04_cohort_estimates.csv"))

cat("\n===== adjusted HR per 1 SD, by cohort =====\n")
for (a in AXES) {
  cat("\n--", a, "--\n")
  sub <- res[exposure == a]
  for (j in seq_len(nrow(sub)))
    cat(sprintf("  %-10s %-10s HR %.2f (%.2f-%.2f) p=%.2g  C=%.3f\n", sub$cohort[j],
        sub$role[j], sub$HR[j], sub$lo[j], sub$hi[j], sub$p[j], sub$cindex[j]))
}

# ---- meta-analysis: validation cohorts only, and all cohorts -----------------
meta_rows <- list()
for (a in AXES) {
  for (scope in c("validation", "all")) {
    sub <- if (scope == "validation") res[exposure == a & role == "validation"] else res[exposure == a]
    m <- meta_dl(log(sub$HR), sub$se, sub$cohort)
    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      axis = a, scope = scope, k = m$k,
      pooled_HR = m$pooled_re, lo = m$lo, hi = m$hi, p = m$p,
      I2 = round(100 * m$I2), tau2 = round(m$tau2, 4), Q = round(m$Q, 2),
      n_total = sum(sub$n), events_total = sum(sub$events))
  }
}
meta_tab <- rbindlist(meta_rows)
fwrite(meta_tab, file.path(OUT, "04_meta_analysis.csv"))
cat("\n===== random-effects meta-analysis =====\n")
print(meta_tab, row.names = FALSE, digits = 3)

# ---- consistency: direction of effect across cohorts ------------------------
cons <- res[, .(n_cohorts = .N,
                n_consistent = sum(HR > 1),
                min_HR = round(min(HR), 2), max_HR = round(max(HR), 2)), by = exposure]
fwrite(cons, file.path(OUT, "04_direction_consistency.csv"))
cat("\n===== direction consistency (HR > 1) =====\n")
print(cons, row.names = FALSE)

saveRDS(list(res = res, meta = meta_tab), file.path(PROJ, "results", "objects", "04_meta.rds"))
cat("\nDONE 04\n")
