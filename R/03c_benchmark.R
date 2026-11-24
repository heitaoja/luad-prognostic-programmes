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

adjustment_row <- function(d, a, covars, id, analysis_type, label, penalise = FALSE) {
  terms <- paste(covars, collapse = " + ")
  if (penalise) {
    stage_sex <- model.matrix(~ stage + sex, data = d)[, -1L, drop = FALSE]
    f <- reformulate(c(a, "age_z", "ridge(stage_sex, theta = 1, scale = FALSE)"),
                     response = "Surv(os_time, os_event)")
    fit <- coxph(f, data = d)
    ci <- exp(confint(fit)[a, ])
    row <- data.frame(n = fit$n, events = fit$nevent, HR = unname(exp(coef(fit)[a])),
                      lo = ci[1L], hi = ci[2L], p = summary(fit)$coefficients[a, "p"])
    terms <- paste0(terms, " (ridge on stage and sex, theta = 1)")
  } else {
    row <- cox_row(d, a, covars)
  }
  data.frame(cohort = id, analysis_type = analysis_type, subgroup_or_model = label,
             exposure = a, covariates = terms, row[c("n", "events", "HR", "lo", "hi", "p")])
}

paired_adjustment <- function(d, id, analysis_type, extra) {
  d <- complete_data(d, c(axes, extra))
  rbindlist(lapply(axes, function(a) {
    rbind(adjustment_row(d, a, clinical, id, analysis_type, "base_clinical"),
          adjustment_row(d, a, c(clinical, extra), id, analysis_type, "extended_adjusted"))
  }))
}

subgroup_axes <- c("proliferation", "hypoxia")
driver_subgroup <- function(d, id, field, status, label) {
  d <- complete_data(d, subgroup_axes)
  d <- droplevels(d[d[[field]] %in% status, ])
  deaths <- d$os_event == 1L
  # Event-free stage or sex categories can cause separation, as in EGFR-mutant GSE72094.
  penalise <- !all(d$stage %in% d$stage[deaths]) || !all(d$sex %in% d$sex[deaths])
  rbindlist(lapply(subgroup_axes, function(a)
    adjustment_row(d, a, clinical, id, "driver_subgroup", label, penalise)))
}

# Any source label outside the known categories remains missing.
recode_factor <- function(x, map) {
  factor(unname(map[tolower(trimws(x))]), levels = unique(unname(map)))
}
smoking_labels <- c(never = "Never", "never-smoker" = "Never",
                    ever = "Ever", "ever-smoker" = "Ever", current = "Ever",
                    ex = "Ever", "ex-smoker" = "Ever")
driver_labels <- c("egfr mutation +" = "EGFR", "egfr mutation (+)" = "EGFR",
                   "kras mutation +" = "KRAS/ALK", "kras mutation (+)" = "KRAS/ALK",
                   "alk-fusion +" = "KRAS/ALK", "eml4-alk fusion (+)" = "KRAS/ALK",
                   "egfr/kras/alk -" = "Triple-negative", "triple negative" = "Triple-negative")

tcga <- scores[["TCGA-LUAD"]]
tcga$pack_years_z <- z(tcga$pack_years)
gse31210 <- scores[["GSE31210"]]
gse31210$smoking_ever <- recode_factor(gse31210$smoking, smoking_labels)
gse31210$driver_group <- recode_factor(gse31210$driver, driver_labels)
gse50081 <- scores[["GSE50081"]]
gse50081$smoking_ever <- recode_factor(gse50081$smoking, smoking_labels)

adjustments <- list(
  paired_adjustment(tcga, "TCGA-LUAD", "pack_years_adjustment", "pack_years_z"),
  paired_adjustment(gse31210, "GSE31210", "smoking_driver_adjustment", c("smoking_ever", "driver_group")),
  paired_adjustment(gse50081, "GSE50081", "smoking_adjustment", "smoking_ever")
)
subgroups <- list(
  driver_subgroup(gse31210, "GSE31210", "driver_group", "EGFR", "EGFR mutant"),
  driver_subgroup(gse31210, "GSE31210", "driver_group", "Triple-negative", "Wild-type / triple-negative")
)

if ("GSE72094" %in% names(scores)) {
  gse72094 <- scores[["GSE72094"]]
  gse72094$smoking_ever <- recode_factor(gse72094$smoking, smoking_labels)
  drivers <- c("egfr", "kras", "tp53", "stk11")
  mutation_labels <- c(wt = "Negative", negative = "Negative", mut = "Positive", positive = "Positive")
  gse72094[drivers] <- lapply(gse72094[drivers], recode_factor, map = mutation_labels)
  adjustments <- c(adjustments, list(
    paired_adjustment(gse72094, "GSE72094", "smoking_adjustment", "smoking_ever"),
    paired_adjustment(gse72094, "GSE72094", "driver_adjustment", drivers)
  ))
  status_labels <- c(Positive = "mut", Negative = "wt")
  subgroups <- c(subgroups, lapply(drivers, function(gene) {
    rbindlist(lapply(names(status_labels), function(status) {
      label <- paste(toupper(gene), status_labels[[status]])
      driver_subgroup(gse72094, "GSE72094", gene, status, label)
    }))
  }))
}

driver_smoking <- rbindlist(c(adjustments, subgroups))
fwrite(driver_smoking, file.path(tab, "03d_driver_smoking_adjustment.csv"))
cat("Wrote", nrow(driver_smoking), "smoking and driver adjustment or subgroup estimates\n")

source("R/lib_external_synthesis.R")
source("R/lib_published_signatures.R")
signatures <- read_published_signatures()
cohorts <- readRDS("data/processed/cohorts.rds")
published <- setNames(lapply(names(scores), function(id)
  published_association_cohort(cohorts[[id]]$expr, scores[[id]], id, signatures)), names(scores))
write_published_benchmark(published, tab)
saveRDS(published, PUBLISHED_ASSOCIATION_CACHE)
