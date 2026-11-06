# =============================================================================
# 05_incremental.R  --  what each axis adds on top of stage
# -----------------------------------------------------------------------------
# In every cohort: C-index with bootstrap CI for the delta, likelihood-ratio
# test, continuous NRI at 36 months, decision-curve analysis and calibration.
# Output: results/tables/05_*.csv
# =============================================================================
suppressPackageStartupMessages({ library(survival); library(data.table) })
source("R/lib_surv.R")
PROJ <- normalizePath(".", winslash = "/")
OUT  <- file.path(PROJ, "results", "tables")
T_HORIZON <- 36          # months
AXES <- c("proliferation", "hypoxia")
COH  <- c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081")

s <- readRDS(file.path(PROJ, "data", "processed", "scores.rds"))

c_mat <- function(fit, d) {
  # C-index computed on the very same complete-case set for both models
  unname(summary(coxph(fit, data = d))$concordance[1])
}

res <- list(); dca_list <- list(); cal_list <- list()
for (id in COH) {
  d <- s[[id]]
  d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
  d$age_z <- z(d$age)
  vars <- c("os_time", "os_event", "stage", "age_z", "sex", AXES)
  d <- d[stats::complete.cases(d[, vars]), ]
  cat(sprintf("\n== %-10s n=%3d events=%3d\n", id, nrow(d), sum(d$os_event)))

  f_clin <- Surv(os_time, os_event) ~ stage + age_z + sex
  d$c_clin <- c_mat(f_clin, d)

  for (a in AXES) {
    f_add <- as.formula(paste("Surv(os_time, os_event) ~", a, "+ stage + age_z + sex"))
    fit0 <- coxph(f_clin, data = d); fit1 <- coxph(f_add, data = d)
    c0 <- unname(summary(fit0)$concordance[1]); c1 <- unname(summary(fit1)$concordance[1])

    set.seed(7); B <- 200
    boot <- replicate(B, {
      b <- d[sample(nrow(d), replace = TRUE), ]
      tryCatch(unname(summary(coxph(f_add, data = b))$concordance[1]) -
               unname(summary(coxph(f_clin, data = b))$concordance[1]),
               error = function(e) NA_real_)
    })
    lr <- anova(fit0, fit1)
    r <- data.frame(cohort = id, axis = a, n = nrow(d), events = sum(d$os_event),
                    c_clin = c0, c_add = c1, delta_c = c1 - c0,
                    delta_lo = quantile(boot, .025, na.rm = TRUE),
                    delta_hi = quantile(boot, .975, na.rm = TRUE),
                    lr_chisq = lr$Chisq[2], lr_p = lr$`Pr(>|Chi|)`[2],
                    aic_clin = AIC(fit0), aic_add = AIC(fit1))

    if (sum(d$os_event) >= 30) {
      r0 <- tryCatch(risk_at(fit0, d, T_HORIZON), error = function(e) NULL)
      r1 <- tryCatch(risk_at(fit1, d, T_HORIZON), error = function(e) NULL)
      if (!is.null(r0) && !is.null(r1)) {
        nri <- nri_continuous(r0, r1, d$os_time, d$os_event, T_HORIZON)
        r$nri <- nri["nri"]; r$nri_event <- nri["nri_event"]; r$nri_nonevent <- nri["nri_nonevent"]
        if (id == "TCGA-LUAD") {
          dc <- dca_surv(d, r1, d$os_time, d$os_event, T_HORIZON,
                         thresholds = seq(0.05, 0.60, 0.01))
          dc$model <- paste0("stage+clin+", a)
          dca_list[[a]] <- dc
          dc0 <- dca_surv(d, r0, d$os_time, d$os_event, T_HORIZON,
                          thresholds = seq(0.05, 0.60, 0.01))
          dc0$model <- "stage+clin"
          dca_list[[paste0(a, "_base")]] <- dc0
        }
        # calibration at the horizon: deciles of predicted risk
        grp <- cut(r1, breaks = quantile(r1, seq(0, 1, 0.2)), include.lowest = TRUE)
        cal <- data.frame(cohort = id, axis = a, grp = levels(grp),
                          predicted = tapply(r1, grp, mean, na.rm = TRUE))
        km <- survfit(Surv(os_time, os_event) ~ grp, data = d)
        s36 <- summary(km, times = T_HORIZON, extend = TRUE)
        cal$observed <- 1 - s36$surv
        cal_list[[paste(id, a)]] <- cal
      }
    }
    res[[length(res) + 1]] <- r
    cat(sprintf("   %-14s C %.3f -> %.3f  dC=%+.3f (%.3f, %.3f)  LR p=%.2g\n",
                a, c0, c1, c1 - c0, r$delta_lo, r$delta_hi, r$lr_p))
  }
}

inc <- rbindlist(res, fill = TRUE)
fwrite(inc, file.path(OUT, "05_incremental_value.csv"))
if (length(dca_list)) fwrite(rbindlist(dca_list, idcol = "key"), file.path(OUT, "05_dca_curves.csv"))
if (length(cal_list)) fwrite(rbindlist(cal_list, fill = TRUE), file.path(OUT, "05_calibration.csv"))

cat("\n===== incremental value over stage + age + sex =====\n")
print(inc[, .(cohort, axis, delta_c = round(delta_c, 3),
              ci = sprintf("%.3f..%.3f", delta_lo, delta_hi),
              lr_p = signif(lr_p, 2), nri = round(nri, 3))], row.names = FALSE)

# pooled delta C-index across validation cohorts
cat("\n===== pooled delta C-index (validation cohorts) =====\n")
pool <- inc[cohort != "TCGA-LUAD", .(mean_delta = round(mean(delta_c), 4),
                                     min_delta = round(min(delta_c), 4),
                                     max_delta = round(max(delta_c), 4)), by = axis]
print(pool, row.names = FALSE)
fwrite(pool, file.path(OUT, "05_delta_c_pooled.csv"))

saveRDS(list(inc = inc, dca = dca_list, cal = cal_list),
        file.path(PROJ, "results", "objects", "05_incremental.rds"))
cat("\nDONE 05\n")
