# External IPD models and the five-cohort synthesis of cohort estimates.
external_ipd_cohort <- function(d, id, axes) {
  d$cohort <- id
  d$age_z <- z(d$age)
  d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
  d$sex <- factor(d$sex, levels = c("female", "male"))
  d <- d[, c("cohort", "os_time", "os_event", "stage", "age_z", "sex", axes)]
  d[complete.cases(d), ]
}

fit_external_ipd <- function(d, scope, axes) {
  terms <- list(IPD_strata_cohort = c("stage", "age_z", "sex", "strata(cohort)"),
                IPD_strata_cohort_stage = c("age_z", "sex", "strata(cohort, stage)"))
  k <- uniqueN(d$cohort)
  rbindlist(lapply(names(terms), function(method) {
    rbindlist(lapply(axes, function(a) {
      f <- reformulate(c(a, terms[[method]]), response = "Surv(os_time, os_event)")
      fit <- coxph(f, data = d, na.action = na.fail)
      s <- summary(fit)
      data.frame(axis = a, scope = scope, method = method, k = k,
                 pooled_HR = s$conf.int[a, "exp(coef)"],
                 lo = s$conf.int[a, "lower .95"], hi = s$conf.int[a, "upper .95"],
                 p = s$coefficients[a, "Pr(>|z|)"], I2 = NA_real_, tau2 = NA_real_, Q = NA_real_,
                 n_total = fit$n, events_total = fit$nevent, df = NA_integer_)
    }))
  }))
}

pool_pan_external <- function(tab, ids, axes) {
  fields <- c("cohort", "exposure", "n", "events", "HR", "lo", "hi")
  rows <- rbindlist(list(
    fread(file.path(tab, "04_cohort_estimates.csv"), select = fields),
    fread(file.path(tab, "05b_rnaseq_association.csv"), select = fields)
  ))[cohort %in% ids]
  stopifnot(nrow(rows) == length(ids) * length(axes), setequal(rows$cohort, ids),
            setequal(rows$exposure, axes), uniqueN(rows, by = c("cohort", "exposure")) == nrow(rows))
  rbindlist(lapply(c("REML_mKH", "DL_normal"), function(m) {
    out <- rbindlist(lapply(axes, function(a)
      cbind(axis = a, scope = "all_external_k5", pool_estimates(rows[exposure == a], m))))
    out[, df := if (m == "REML_mKH") k - 1L else NA_integer_]
    out
  }))
}

write_external_synthesis <- function(scores, tab, cptac = NULL) {
  axes <- names(PROGRAMMES)
  scopes <- list(validation_array_k3 = c("GSE31210", "GSE37745", "GSE50081"))
  if (!is.null(cptac)) {
    scores[["CPTAC-LUAD"]] <- cptac
    scopes$extended_array_k4 <- c(scopes$validation_array_k3, "GSE72094")
    scopes$all_external_k5 <- c(scopes$extended_array_k4, "CPTAC-LUAD")
  }
  ids <- unique(unlist(scopes, use.names = FALSE))
  stopifnot(all(ids %in% names(scores)))
  patients <- setNames(lapply(ids, function(id) external_ipd_cohort(scores[[id]], id, axes)), ids)
  out <- rbindlist(lapply(names(scopes), function(scope)
    fit_external_ipd(rbindlist(patients[scopes[[scope]]]), scope, axes)))
  if (!is.null(cptac)) {
    out <- rbindlist(list(out, pool_pan_external(tab, scopes$all_external_k5, axes)))
  }
  out[, p_holm := p.adjust(p, method = "holm"), by = .(scope, method)]
  fwrite(out, file.path(tab, "04b_ipd_and_pan_external_synthesis.csv"))
  cat("Wrote", nrow(out), "external IPD and meta-analysis estimates\n")
}
