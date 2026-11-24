# Primary TCGA-LUAD survival models, concordance ladder and stage-stratified analyses.
suppressPackageStartupMessages({ library(survival); library(data.table) })
source("R/lib_surv.R")
PROJ <- normalizePath(".", winslash = "/")
OUT  <- file.path(PROJ, "results", "tables")
AXES <- c("proliferation", "hypoxia", "emt_stroma", "immune")

s <- readRDS(file.path(PROJ, "data", "processed", "scores.rds"))
d <- s[["TCGA-LUAD"]]
d$stage  <- factor(d$stage, levels = c("I", "II", "III", "IV"))
d$age_z  <- z(d$age)
COV      <- c("stage", "age_z", "sex")

cat("=== TCGA-LUAD: n =", nrow(d), " events =", sum(d$os_event), "===\n")
cat("stage:", paste(names(table(d$stage, useNA = "ifany")),
                    table(d$stage, useNA = "ifany"), sep = "=", collapse = " "), "\n")

# Cox models
rows_uni <- rbindlist(lapply(AXES, function(a)
  cox_row(d, a, character(0), "univariable", "TCGA-LUAD")))
rows_adj <- rbindlist(lapply(AXES, function(a)
  cox_row(d, a, COV, "adjusted", "TCGA-LUAD")))
fwrite(rows_uni, file.path(OUT, "03_cox_univariable.csv"))
fwrite(rows_adj, file.path(OUT, "03_cox_adjusted.csv"))

cat("\n--- univariable (per 1 SD) ---\n")
for (i in seq_len(nrow(rows_uni)))
  cat(sprintf("%-14s HR %.2f (%.2f-%.2f) p=%.2g\n", rows_uni$exposure[i],
      rows_uni$HR[i], rows_uni$lo[i], rows_uni$hi[i], rows_uni$p[i]))
cat("\n--- adjusted for stage + age + sex ---\n")
for (i in seq_len(nrow(rows_adj)))
  cat(sprintf("%-14s HR %.2f (%.2f-%.2f) p=%.2g  C=%.3f\n", rows_adj$exposure[i],
      rows_adj$HR[i], rows_adj$lo[i], rows_adj$hi[i], rows_adj$p[i], rows_adj$cindex[i]))

# Adjusted stage effects from the model with the lowest p-value
best <- rows_adj$exposure[which.min(rows_adj$p)]
fit_best <- coxph(as.formula(paste("Surv(os_time, os_event) ~", best,
                                  "+", paste(COV, collapse = " + "))), data = d)
st <- data.frame(term = names(coef(fit_best)), HR = exp(coef(fit_best)),
                 lo = exp(confint(fit_best))[, 1], hi = exp(confint(fit_best))[, 2],
                 p = summary(fit_best)$coef[, "Pr(>|z|)"])
fwrite(st, file.path(OUT, "03_adjusted_model_terms.csv"))
cat("\n--- covariates in the model with the strongest axis (", best, ") ---\n")
print(st, row.names = FALSE, digits = 3)

# Concordance ladder
dc <- d[complete.cases(d[, c("os_time", "os_event", COV, AXES)]), ]
clinical_formula <- Surv(os_time, os_event) ~ stage + age_z + sex
fit_clin <- coxph(clinical_formula, data = dc)
c_clin <- unname(summary(fit_clin)$concordance[1L])
models <- c("clinical", "stage", AXES)
lad <- rbindlist(lapply(models, function(a) {
  f <- switch(a, clinical = clinical_formula, stage = Surv(os_time, os_event) ~ stage,
              reformulate(c(COV, a), response = "Surv(os_time, os_event)"))
  fit <- coxph(f, data = dc)
  cindex <- unname(summary(fit)$concordance[1L])
  p <- if (a %in% AXES) anova(fit_clin, fit)[["Pr(>|Chi|)"]][2L] else NA_real_
  data.frame(model = if (a %in% AXES) paste("clinical +", a) else a,
             variable = a, n = nrow(dc), cindex = cindex,
             delta = cindex - c_clin, LR_p = p, AIC = AIC(fit))
}))
fwrite(lad, file.path(OUT, "03_cindex_ladder.csv"))
cat("\n--- C-index ladder (TCGA) ---\n")
print(lad[, c("model", "cindex", "delta", "LR_p", "AIC")], row.names = FALSE, digits = 3)

# Stage-stratified models
strat <- rbindlist(lapply(AXES, function(a) rbindlist(lapply(levels(d$stage)[1:3], function(st) {
  x <- d[which(d$stage == st), ]
  if (nrow(x) < 40 || sum(x$os_event) < 10) return(NULL)
  x$age_z <- z(x$age)
  r <- cox_row(x, a, c("age_z", "sex"), paste0("stage ", st), "TCGA-LUAD")
  r$stage <- st; r
}))))
fwrite(strat, file.path(OUT, "03_stage_stratified.csv"))
cat("\n--- stage-stratified (adjusted for age, sex) ---\n")
for (i in seq_len(nrow(strat)))
  cat(sprintf("%-14s %-8s n=%3d ev=%2d HR %.2f (%.2f-%.2f) p=%.2g\n", strat$exposure[i],
      strat$stage[i], strat$n[i], strat$events[i], strat$HR[i], strat$lo[i],
      strat$hi[i], strat$p[i]))

# Exploratory stage interaction tests
interactions <- rbindlist(lapply(AXES, function(a) {
  additive <- coxph(reformulate(c(a, COV), response = "Surv(os_time, os_event)"), data = dc)
  interaction <- coxph(reformulate(c(paste0(a, " * stage"), "age_z", "sex"),
                                   response = "Surv(os_time, os_event)"), data = dc)
  data.frame(axis = a, omnibus_LR_p = anova(additive, interaction)[["Pr(>|Chi|)"]][2L])
}))
fwrite(interactions, file.path(OUT, "03_stage_interaction.csv"))

# Secondary proliferation-immune balance model
fb <- cox_row(d, "balance", COV, "adjusted", "TCGA-LUAD")
fpi <- coxph(Surv(os_time, os_event) ~ proliferation * immune + stage + age_z + sex, data = d)
int_p <- summary(fpi)$coef["proliferation:immune", "Pr(>|z|)"]
bal <- data.frame(
  analysis = c("balance = z(proliferation) - z(immune), adjusted",
               "proliferation x immune interaction p"),
  estimate = c(fb$HR, NA), lo = c(fb$lo, NA), hi = c(fb$hi, NA),
  p = c(fb$p, int_p))
fwrite(bal, file.path(OUT, "03_balance_secondary.csv"))
cat("\n--- secondary: balance model ---\n")
cat(sprintf("balance adj HR %.2f (%.2f-%.2f) p=%.3g | proliferation x immune p=%.2f\n",
            fb$HR, fb$lo, fb$hi, fb$p, int_p))

saveRDS(list(d = d, lad = lad, rows_uni = rows_uni, rows_adj = rows_adj,
             strat = strat, balance = bal), file.path(PROJ, "results", "objects", "03_tcga.rds"))
cat("\nDONE 03\n")
