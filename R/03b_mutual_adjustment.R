# =============================================================================
# 03b_mutual_adjustment.R  --  all four axes in one model (mutual adjustment)
# -----------------------------------------------------------------------------
# Addresses the obvious reviewer question: are the axes just proxies for each
# other (and for tumour purity / stromal content)?
# =============================================================================
suppressPackageStartupMessages({ library(survival); library(data.table) })
source("R/lib_surv.R")
PROJ <- normalizePath(".", winslash = "/"); TAB <- file.path(PROJ, "results", "tables")
AXES <- c("proliferation", "hypoxia", "emt_stroma", "immune")
s <- readRDS(file.path(PROJ, "data", "processed", "scores.rds"))
out <- list()
for (id in names(s)) {
  d <- s[[id]]; d$stage <- factor(d$stage, levels = c("I","II","III","IV")); d$age_z <- z(d$age)
  f <- as.formula(paste("Surv(os_time, os_event) ~",
                        paste(c(AXES, "stage", "age_z", "sex"), collapse = " + ")))
  fit <- coxph(f, data = d)
  s2 <- summary(fit)$coef; ci <- confint(fit)
  for (a in AXES) out[[length(out)+1]] <- data.frame(
    cohort = id, axis = a, n = fit$n, events = fit$nevent,
    HR = exp(s2[a, "coef"]), lo = exp(ci[a, 1]), hi = exp(ci[a, 2]),
    p = s2[a, "Pr(>|z|)"], cindex = unname(summary(fit)$concordance[1]))
}
res <- rbindlist(out)
fwrite(res, file.path(TAB, "03b_mutual_adjustment.csv"))
cat("=== all four axes + stage + age + sex in one model ===\n")
for (id in unique(res$cohort)) {
  cat("\n--", id, "--\n"); print(res[cohort == id, .(axis, HR = round(HR, 2),
      CI = sprintf("%.2f-%.2f", lo, hi), p = signif(p, 2))], row.names = FALSE)
}
cat("\nDONE 03b\n")
