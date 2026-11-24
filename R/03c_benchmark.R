# Compare the programme with simple markers on identical patients and examine
# the effects of modelling assumptions, specimen choice and tumour purity.
source("R/lib_surv.R")
source("R/lib_score.R")
scores <- readRDS("data/processed/scores.rds")
alt <- readRDS("data/processed/tcga_alternate_scores.rds")
tab <- "results/tables"
axes <- names(PROGRAMMES)
markers <- c("MKI67", "e2f", "g2m")
clinical <- c("stage", "age_z", "sex")
complete_data <- function(d, features = axes) {
  d$age_z <- z(d$age)
  d <- d[complete.cases(d[, c("os_time", "os_event", clinical, features)]), ]
  d$stage <- droplevels(d$stage)
  d
}

bench <- rbindlist(lapply(names(scores), function(id) {
  d <- complete_data(scores[[id]], c("proliferation", markers))
  simple <- rbindlist(lapply(c("proliferation", markers), function(feature)
    cox_row(d, feature, clinical, "clinical_plus_marker", id)))
  added <- rbindlist(lapply(markers, function(feature)
    cox_row(d, "proliferation", c(clinical, feature),
            paste0("composite_given_", feature), id)))
  rbindlist(list(simple, added))
}))
fwrite(bench, file.path(tab, "03c_marker_benchmark.csv"))
correlations <- rbindlist(lapply(names(scores), function(id) {
  d <- complete_data(scores[[id]], c("proliferation", markers))
  rho <- cor(d$proliferation, d[, markers], method = "spearman")
  data.frame(cohort = id, marker = markers, n = nrow(d), spearman = as.numeric(rho))
}))
fwrite(correlations, file.path(tab, "03c_marker_correlations.csv"))

scopes <- list(validation = setdiff(names(scores), c("TCGA-LUAD", "GSE72094")))
if ("GSE72094" %in% names(scores)) scopes$extended <- setdiff(names(scores), "TCGA-LUAD")
joint <- fread(file.path(tab, "03b_mutual_adjustment.csv"))
comparisons <- unique(bench[, .(model, exposure)])
pool <- rbindlist(lapply(names(scopes), function(scope) {
  ids <- scopes[[scope]]
  baseline <- rbindlist(lapply(seq_len(nrow(comparisons)), function(i) {
    key <- comparisons[i]
    sub <- bench[cohort %in% ids & model == key$model & exposure == key$exposure]
    cbind(scope = scope, key, pool_estimates(sub))
  }))
  mutual <- rbindlist(lapply(axes, function(axis_name)
    cbind(scope = scope, model = "four_programmes_joint", exposure = axis_name,
          pool_estimates(joint[cohort %in% ids & axis == axis_name]))))
  rbindlist(list(baseline, mutual))
}))
fwrite(pool, file.path(tab, "03c_marker_synthesis.csv"))

checks <- rbindlist(lapply(names(scores), function(id) {
  d <- complete_data(scores[[id]])
  early <- d[d$stage %in% c("I", "II"), ]
  rbindlist(lapply(axes, function(a) {
    f <- reformulate(c(a, clinical), response = "Surv(os_time, os_event)")
    fit <- coxph(f, data = d, x = TRUE)
    ph <- cox.zph(fit)$table
    flexible <- coxph(update(f, paste(". ~ . -", a, "+ splines::ns(", a, ", df = 3)")),
                      data = d)
    stratified_fit <- coxph(reformulate(c(a, "age_z", "sex", "strata(stage)"),
                                      response = "Surv(os_time, os_event)"), data = d)
    strat_ci <- exp(confint(stratified_fit)[a, ])
    restricted <- cox_row(early, a, clinical, "stage_I_II", id)
    data.frame(cohort = id, axis = a, score_ph_p = ph[a, "p"], global_ph_p = ph["GLOBAL", "p"],
               spline_vs_linear_p = anova(fit, flexible)[["Pr(>|Chi|)"]][2L],
               stage_stratified_HR = unname(exp(coef(stratified_fit)[a])),
               stage_stratified_lo = strat_ci[1L], stage_stratified_hi = strat_ci[2L],
               early_n = restricted$n, early_events = restricted$events,
               early_HR = restricted$HR, early_lo = restricted$lo, early_hi = restricted$hi)
  }))
}))
fwrite(checks, file.path(tab, "03c_model_diagnostics.csv"))

primary <- complete_data(scores[["TCGA-LUAD"]])
alternate <- complete_data(alt)
common <- intersect(primary$patient, alternate$patient)
specimen <- rbindlist(lapply(axes, function(a)
  rbindlist(list(
    cox_row(primary[primary$patient %in% common, ], a, clinical, "prefer_first_vial", "TCGA-LUAD"),
    cox_row(alternate[alternate$patient %in% common, ], a, clinical, "prefer_last_vial", "TCGA-LUAD")
  ))))
fwrite(specimen, file.path(tab, "03c_specimen_sensitivity.csv"))

purity_file <- "data/raw/tcga/TCGA_purity_Aran.rda"
if (!file.exists(purity_file) && Sys.getenv("LUAD_CORE_ONLY") != "1")
  stop("Missing DNA-based purity input specified in docs/data_provenance.md")
purity_out <- data.frame(cohort = character(), model = character(), exposure = character(),
                         n = integer(), events = integer(), HR = numeric(), lo = numeric(),
                         hi = numeric(), p = numeric(), cindex = numeric(), aic = numeric(),
                         correlation_with_purity = numeric())
if (file.exists(purity_file) && Sys.getenv("LUAD_CORE_ONLY") != "1") {
  source_data <- new.env(parent = emptyenv())
  load(purity_file, envir = source_data)
  purity <- source_data$Tumor.purity
  purity <- purity[as.character(purity$Cancer.type) == "LUAD", ]
  key <- substr(as.character(purity$Sample.ID), 1L, 16L)
  # This upstream data snapshot stores decimal fractions with commas.
  value <- as.numeric(sub(",", ".", as.character(purity$ABSOLUTE), fixed = TRUE))
  index <- is.finite(value)
  map <- unique(data.frame(sample = key[index], purity = value[index]))
  stopifnot(!anyDuplicated(map$sample), all(map$purity >= 0 & map$purity <= 1))
  primary$purity <- map$purity[match(primary$sample, map$sample)]
  matched <- primary[!is.na(primary$purity), ]
  if (nrow(matched) < 50L || sum(matched$os_event) < 20L)
    stop("Insufficient exact specimen matches with ABSOLUTE purity")
  purity_out <- rbindlist(lapply(c(axes, markers), function(a) {
    rows <- rbindlist(list(cox_row(matched, a, clinical, "clinical", "TCGA-LUAD"),
                         cox_row(matched, a, c(clinical, "purity"), "clinical_and_purity", "TCGA-LUAD")))
    rows$correlation_with_purity <- cor(matched[[a]], matched$purity, method = "spearman")
    rows
  }))
}
fwrite(purity_out, file.path(tab, "03c_purity_sensitivity.csv"))
cat("Matched clinical benchmark, model diagnostics and alternate specimen analyses written\n")
