# =============================================================================
# lib_surv.R  --  shared survival / validation helpers
# -----------------------------------------------------------------------------
# Sourced by 03-07.  Everything here is deliberately explicit (no black-box
# wrappers) so that each number in the paper can be traced to a line of code.
# =============================================================================

suppressPackageStartupMessages({
  library(survival); library(data.table); library(ggplot2)
})

z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

# ---------------------------------------------------------------- Cox models --
# Fit Surv ~ exposure + covariates and return one tidy row.
cox_row <- function(data, exposure, covars = character(0), label = "model",
                    cohort = NA_character_) {
  covars <- covars[covars %in% names(data)]
  covars <- covars[vapply(covars, function(v) sum(!is.na(data[[v]])) > 0, logical(1))]
  f <- as.formula(paste("Surv(os_time, os_event) ~", paste(c(exposure, covars), collapse = " + ")))
  d <- data[stats::complete.cases(data[, c("os_time", "os_event", exposure, covars)]), ]
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

# C-index of a fitted model, with an optional bootstrap CI
cindex_ci <- function(data, exposure, covars = character(0), B = 200) {
  covars <- covars[covars %in% names(data)]
  f <- as.formula(paste("Surv(os_time, os_event) ~", paste(c(exposure, covars), collapse = " + ")))
  d <- data[stats::complete.cases(data[, c("os_time", "os_event", exposure, covars)]), ]
  fit <- coxph(f, data = d)
  est <- unname(summary(fit)$concordance[1])
  if (B < 1) return(c(cindex = est, lo = NA, hi = NA))
  set.seed(11)
  boot <- replicate(B, {
    b <- d[sample(nrow(d), replace = TRUE), ]
    tryCatch(unname(summary(coxph(f, data = b))$concordance[1]), error = function(e) NA_real_)
  })
  c(cindex = est, lo = quantile(boot, .025, na.rm = TRUE), hi = quantile(boot, .975, na.rm = TRUE))
}

# ------------------------------------------------------------------- meta ----
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

# ------------------------------------------------------- decision curve ------
# Survival DCA at horizon t with IPCW (Vickers).  Returns net benefit by
# threshold for "treat all", "treat none" and the model.
dca_surv <- function(data, risk, time, event, t, thresholds = seq(0.01, 0.5, 0.01)) {
  d <- data.frame(risk = risk, time = time, event = event)
  d <- d[stats::complete.cases(d), ]
  # censoring distribution (KM of censoring)
  cf <- survfit(Surv(time, 1 - event) ~ 1, data = d)
  Ghat <- stepfun(cf$time, c(1, cf$surv))
  g_at <- function(tm) pmax(Ghat(tm), 1e-6)

  ev_t  <- d$event == 1 & d$time <= t
  ok_ne <- d$time > t
  w <- rep(0, nrow(d))
  w[ev_t]  <- 1 / g_at(d$time[ev_t])
  w[ok_ne] <- 1 / g_at(t)

  out <- data.frame(threshold = thresholds, nb_model = NA_real_,
                    nb_all = NA_real_, nb_none = 0)
  n <- nrow(d)
  for (i in seq_along(thresholds)) {
    pt <- thresholds[i]
    pos <- d$risk >= pt
    tp <- sum(w[ev_t & pos]); fp <- sum(w[ok_ne & pos])
    out$nb_model[i] <- (tp - fp * pt / (1 - pt)) / n
    tp_all <- sum(w[ev_t]); fp_all <- sum(w[ok_ne])
    out$nb_all[i] <- (tp_all - fp_all * pt / (1 - pt)) / n
  }
  out
}

# --------------------------------------------------------- classification -----
# continuous NRI at horizon t (event vs no event by t), censored before t dropped
nri_continuous <- function(risk_old, risk_new, time, event, t) {
  ev <- event == 1 & time <= t
  ne <- time > t
  keep <- ev | ne
  up  <- risk_new[keep] > risk_old[keep]
  dn  <- risk_new[keep] < risk_old[keep]
  e <- ev[keep]
  nri_e  <- mean(up[e])  - mean(dn[e])
  nri_ne <- mean(dn[!e]) - mean(up[!e])
  c(nri = nri_e + nri_ne, nri_event = nri_e, nri_nonevent = nri_ne,
    n = sum(keep), events = sum(e))
}

# predicted absolute risk at horizon t from a coxph model
risk_at <- function(fit, newdata, t) {
  sf <- survfit(fit, newdata = newdata)
  s <- summary(sf, times = t, extend = TRUE)$surv
  1 - as.numeric(s)
}

# ------------------------------------------------------------- KM plotting ---
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

logrank_p <- function(data, group) {
  sd <- survdiff(as.formula(paste("Surv(os_time, os_event) ~", group)), data = data)
  1 - pchisq(sd$chisq, length(sd$n) - 1)
}
