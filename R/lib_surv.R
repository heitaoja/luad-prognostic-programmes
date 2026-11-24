# Survival modelling, random-effects meta-analysis and Kaplan-Meier helpers.

suppressPackageStartupMessages({
  library(survival); library(data.table); library(ggplot2)
})

z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(ifelse(is.na(x), NA_real_, 0))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

# Cox regression helpers
# Fit Surv ~ exposure + covariates and return one tidy row.
cox_row <- function(data, exposure, covars = character(0), label = "model",
                    cohort = NA_character_, stratify = character(0)) {
  required <- unique(c("os_time", "os_event", exposure, covars, stratify))
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing model variables: ", paste(missing, collapse = ", "))
  strata_term <- if (length(stratify)) paste0("strata(", paste(stratify, collapse = ", "), ")")
  f <- as.formula(paste("Surv(os_time, os_event) ~", paste(c(exposure, covars, strata_term), collapse = " + ")))
  d <- data[stats::complete.cases(data[, required]), ]
  fit <- coxph(f, data = d)
  s <- summary(fit)
  ci <- suppressMessages(confint(fit))
  data.frame(cohort = cohort, model = label, exposure = exposure,
             n = nrow(d), events = sum(d$os_event),
             HR = exp(coef(fit))[1], lo = exp(ci[1, 1]), hi = exp(ci[1, 2]),
             p = s$coef[1, "Pr(>|z|)"],
             cindex = unname(s$concordance[1]),
             aic = stats::AIC(fit),
             stringsAsFactors = FALSE)
}

# Random-effects meta-analysis helpers
# DerSimonian-Laird random effects on the log-HR scale.
meta_dl <- function(est, se, labels = NULL) {
  k <- length(est); w <- 1 / se^2
  fe <- sum(w * est) / sum(w)
  q  <- sum(w * (est - fe)^2)
  c_ <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - (k - 1)) / c_)
  wr <- 1 / (se^2 + tau2)
  re <- sum(wr * est) / sum(wr); se_re <- sqrt(1 / sum(wr))
  I2 <- ifelse(q > 0, max(0, (q - (k - 1)) / q), 0)
  list(k = k, pooled_fe = exp(fe), pooled_re = exp(re),
       lo = exp(re - 1.96 * se_re), hi = exp(re + 1.96 * se_re),
       p = 2 * pnorm(-abs(re / se_re)), tau2 = tau2, I2 = I2, Q = q)
}

# REML estimate of between-study variance and a modified Knapp–Hartung interval.
# The Hartung–Knapp variance multiplier is bounded below by one.
meta_reml_mkh <- function(est, se) {
  stopifnot(length(est) == length(se), length(est) >= 2L,
            all(is.finite(est)), all(is.finite(se) & se > 0))
  inverse_variance <- function(tau2) 1 / (se^2 + tau2)
  loss <- function(tau2) {
    w <- inverse_variance(tau2)
    center <- weighted.mean(est, w)
    sum(log(se^2 + tau2)) + log(sum(w)) + sum(w * (est - center)^2)
  }
  bound <- max(1, 10 * stats::var(est), 10 * max(se^2))
  minimum <- optimize(loss, interval = c(0, bound), tol = 1e-10)
  tau2 <- if (loss(0) <= minimum$objective) 0 else minimum$minimum
  w <- inverse_variance(tau2)
  center <- weighted.mean(est, w)
  df <- length(est) - 1L
  se_mkh <- sqrt(max(1, sum(w * (est - center)^2) / df) / sum(w))
  ci <- center + c(-1, 1) * qt(0.975, df) * se_mkh
  w0 <- inverse_variance(0)
  typical_var <- df * sum(w0) / (sum(w0)^2 - sum(w0^2))
  list(k = length(est), pooled_re = exp(center), lo = exp(ci[1L]), hi = exp(ci[2L]),
       p = 2 * pt(-abs(center / se_mkh), df), I2 = tau2 / (tau2 + typical_var),
       tau2 = tau2, Q = sum(w0 * (est - weighted.mean(est, w0))^2))
}

pool_estimates <- function(rows, method = "REML_mKH") {
  se <- (log(rows$hi) - log(rows$lo)) / (2 * qnorm(0.975))
  fit <- switch(method, REML_mKH = meta_reml_mkh(log(rows$HR), se),
                DL_normal = meta_dl(log(rows$HR), se),
                stop("Unknown synthesis method: ", method))
  data.frame(method = method, k = fit$k, pooled_HR = fit$pooled_re, lo = fit$lo,
             hi = fit$hi, p = fit$p, I2 = 100 * fit$I2, tau2 = fit$tau2, Q = fit$Q,
             n_total = sum(rows$n), events_total = sum(rows$events))
}

# Kaplan-Meier data preparation
km_df <- function(data, group, time = "os_time", event = "os_event") {
  d <- data.frame(time = data[[time]], event = data[[event]], group = data[[group]])
  d <- d[stats::complete.cases(d), ]
  fit <- survfit(Surv(time, event) ~ group, data = d)
  s <- summary(fit)$table
  df <- data.frame(time = fit$time, surv = fit$surv,
                   n.censor = fit$n.censor,
                   group = rep(names(fit$strata), fit$strata))
  list(df = df, table = s, data = d)
}
