# Clinical-adjusted association estimates and small-k external synthesis.
source("R/lib_surv.R")
source("R/lib_score.R")
scores <- readRDS("data/processed/scores.rds")
audit <- fread("results/tables/01_cohort_audit.csv")
tab <- "results/tables"
core <- c("GSE31210", "GSE37745", "GSE50081")
axes <- names(PROGRAMMES)

estimates <- rbindlist(lapply(names(scores), function(id) {
  d <- scores[[id]]
  d$age_z <- z(d$age)
  rows <- rbindlist(lapply(axes, function(a)
    cox_row(d, a, c("stage", "age_z", "sex"), "adjusted", id)))
  rows$role <- if (id == "TCGA-LUAD") "discovery" else if (id %in% core) "validation" else "extension"
  rows$median_fu <- audit[cohort == id, median_followup_mo]
  rows
}))
estimates[, se := (log(hi) - log(lo)) / (2 * qnorm(0.975))]
fwrite(estimates, file.path(tab, "04_cohort_estimates.csv"))

external <- estimates[role != "discovery"]
scopes <- list(validation = external[cohort %in% core])
if (any(external$role == "extension")) scopes$extended <- external

pooled <- rbindlist(lapply(names(scopes), function(s) {
  scope_rows <- scopes[[s]]
  rbindlist(lapply(c("REML_mKH", "DL_normal"), function(m) {
    out <- rbindlist(lapply(axes, function(a)
      cbind(axis = a, scope = s, pool_estimates(scope_rows[exposure == a], m))))
    out[, p_holm := p.adjust(p, method = "holm")]
    out
  }))
}))
fwrite(pooled, file.path(tab, "04_meta_analysis.csv"))

loo <- rbindlist(lapply(names(scopes), function(s) {
  scope_rows <- scopes[[s]]
  rbindlist(lapply(unique(scope_rows$cohort), function(omit)
    rbindlist(lapply(axes, function(a)
      cbind(axis = a, scope = s, omitted = omit,
            pool_estimates(scope_rows[exposure == a & cohort != omit]))))))
}))
fwrite(loo, file.path(tab, "04_leave_one_out.csv"))
direction <- estimates[, .(n_cohorts = .N, n_positive = sum(HR > 1),
                           min_HR = min(HR), max_HR = max(HR)), by = exposure]
fwrite(direction, file.path(tab, "04_direction_consistency.csv"))
print(pooled[, .(axis, scope, method, k, pooled_HR, lo, hi, p_holm, I2)], digits = 3)

if (Sys.getenv("LUAD_CORE_ONLY") == "1") {
  source("R/lib_external_synthesis.R")
  write_external_synthesis(scores, tab)
}
