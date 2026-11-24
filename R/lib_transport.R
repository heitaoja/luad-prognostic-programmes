# TCGA fitting and external evaluation shared by the array and RNA-seq analyses.
TRANSPORT_MARKERS <- c(composite = "proliferation", MKI67 = "MKI67", E2F = "e2f", G2M = "g2m")
TRANSPORT_CLINICAL <- c("stage", "age_10", "sex")

transport_design <- function(data, reference, mki67_reference) {
  d <- data$clin
  stopifnot(identical(d$sample, rownames(data$scores)))
  molecular <- programme_scores(data$scores, reference)$scores
  molecular$MKI67 <- standardise(data$mki67_quantile, mki67_reference)
  d <- cbind(d, molecular)
  d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
  d$sex <- factor(d$sex, levels = c("female", "male"))
  d$age_10 <- (d$age - 65) / 10
  fields <- c("os_time", "os_event", TRANSPORT_CLINICAL, TRANSPORT_MARKERS)
  d[complete.cases(d[, fields]), ]
}

fit_frozen_transport <- function(train, horizon = 36) {
  reference <- programme_scores(train$scores)$reference
  mki67_reference <- z_parameters(train$mki67_quantile)
  training <- transport_design(train, reference, mki67_reference)
  fits <- lapply(c(clinical = "", TRANSPORT_MARKERS), function(feature) {
    formula <- reformulate(c(TRANSPORT_CLINICAL, feature[nzchar(feature)]),
                           response = "Surv(os_time, os_event)")
    coxph(formula, data = training, singular.ok = FALSE, x = TRUE)
  })
  baseline <- vapply(fits, function(fit) {
    hazard <- basehaz(fit, centered = FALSE)
    c(0, hazard$hazard)[1L + findInterval(horizon, hazard$time)]
  }, numeric(1))
  list(fits = fits, reference = reference, mki67_reference = mki67_reference,
       horizon = horizon, baseline = baseline)
}

transport_coefficients <- function(model) {
  rbindlist(lapply(names(model$fits), function(name) {
    beta <- coef(model$fits[[name]])
    data.frame(model = name, term = names(beta), coefficient = unname(beta),
               H0_36 = model$baseline[[name]])
  }))
}

transport_score_reference <- function(model) {
  reference <- model$reference
  rbindlist(list(
    data.frame(component = names(reference$components),
               do.call(rbind, reference$components), row.names = NULL),
    data.frame(component = names(reference$programmes),
               do.call(rbind, reference$programmes), row.names = NULL),
    data.frame(component = "MKI67_rank_percentile", t(model$mki67_reference))
  ))
}

c_value <- function(time, status, lp) {
  survival::concordance(Surv(time, status) ~ lp, reverse = TRUE)$concordance
}

ipcw_weights <- function(d, horizon) {
  fit <- survfit(Surv(os_time, 1L - os_event) ~ 1, data = d)
  surv_at <- function(time, left = FALSE) {
    i <- if (left) sum(fit$time < time) else sum(fit$time <= time)
    c(1, fit$surv)[i + 1L]
  }
  event <- d$os_event == 1L & d$os_time <= horizon
  control <- d$os_time > horizon
  weights <- rep(0, nrow(d))
  censoring_survival <- surv_at(horizon)
  if (censoring_survival <= 0) stop("Censoring survival is zero at the horizon")
  weights[event] <- 1 / vapply(d$os_time[event], surv_at, numeric(1), left = TRUE)
  weights[control] <- 1 / censoring_survival
  list(event = event, weights = weights)
}

evaluate_frozen_transport <- function(model, data, cohort, B = 1000L) {
  d <- transport_design(data, model$reference, model$mki67_reference)
  fits <- model$fits
  horizon <- model$horizon
  lp <- vapply(fits, predict, numeric(nrow(d)), newdata = d, type = "lp", reference = "zero")
  risk <- -expm1(-sweep(exp(lp), 2L, model$baseline, "*"))
  c_index <- apply(lp, 2L, function(x) c_value(d$os_time, d$os_event, x))
  boot <- replicate(B, {
    i <- sample.int(nrow(d), replace = TRUE)
    apply(lp[i, , drop = FALSE], 2L, function(x) c_value(d$os_time[i], d$os_event[i], x))
  })
  events_36 <- sum(d$os_event == 1L & d$os_time <= horizon)
  n_at_risk <- sum(d$os_time >= horizon)
  early_censored <- sum(d$os_event == 0L & d$os_time <= horizon)
  observed <- NA_real_
  brier <- setNames(rep(NA_real_, length(fits)), names(fits))
  if (events_36 >= 10L && n_at_risk >= 10L) {
    km <- survfit(Surv(os_time, os_event) ~ 1, data = d)
    observed <- 1 - summary(km, times = horizon)$surv
    weight <- ipcw_weights(d, horizon)
    brier <- colMeans(weight$weights * (risk - weight$event)^2)
  }
  c_interval <- function(x) quantile(x, c(.025, .975), na.rm = TRUE)
  evaluation <- rbindlist(lapply(names(fits), function(m) {
    c_ci <- c_interval(boot[m, ])
    delta_ci <- c_interval(boot[m, ] - boot["clinical", ])
    d$lp <- lp[, m]
    slope <- coxph(Surv(os_time, os_event) ~ lp, data = d)
    slope_ci <- confint(slope)["lp", ]
    data.frame(cohort = cohort, model = m, n = nrow(d), events = sum(d$os_event),
               horizon_months = horizon, events_36 = events_36, n_at_risk_36 = n_at_risk,
               censored_by_36 = early_censored, c_index = c_index[m],
               c_lo = c_ci[1L], c_hi = c_ci[2L],
               delta_c = c_index[m] - c_index["clinical"],
               delta_lo = delta_ci[1L], delta_hi = delta_ci[2L],
               calibration_slope = coef(slope)[["lp"]],
               slope_lo = slope_ci[1L], slope_hi = slope_ci[2L],
               predicted_36 = mean(risk[, m]), observed_36 = observed,
               brier_36 = brier[[m]], brier_gain_36 = brier[["clinical"]] - brier[[m]])
  }))
  contrasts <- rbindlist(lapply(c("MKI67", "E2F", "G2M"), function(m) {
    ci <- c_interval(boot["composite", ] - boot[m, ])
    data.frame(cohort = cohort, comparison = m, delta_c = c_index["composite"] - c_index[m],
               delta_lo = ci[1L], delta_hi = ci[2L])
  }))
  calibration <- data.table(cohort = character(), model = character(), group = integer(),
                            n = integer(), at_risk = integer(), predicted = numeric(),
                            observed = numeric(), observed_lo = numeric(), observed_hi = numeric())
  if (is.finite(observed)) {
    calibration <- rbindlist(lapply(c("clinical", "composite"), function(m) {
      breaks <- unique(quantile(risk[, m], seq(0, 1, length.out = 5)))
      group <- if (length(breaks) == 1L) rep(1L, nrow(d)) else
        cut(risk[, m], breaks, include.lowest = TRUE, labels = FALSE)
      rbindlist(lapply(sort(unique(group)), function(g) {
        keep <- group == g
        sf <- survfit(Surv(os_time, os_event) ~ 1, data = d[keep, ])
        surv <- summary(sf, times = horizon)
        data.frame(cohort = cohort, model = m, group = g, n = sum(keep),
                   at_risk = sum(keep & d$os_time >= horizon), predicted = mean(risk[keep, m]),
                   observed = if (length(surv$surv)) 1 - surv$surv else NA_real_,
                   observed_lo = if (length(surv$upper)) 1 - surv$upper else NA_real_,
                   observed_hi = if (length(surv$lower)) 1 - surv$lower else NA_real_)
      }))
    }))
  }
  list(evaluation = evaluation, contrasts = contrasts, calibration = calibration)
}
