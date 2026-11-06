# =============================================================================
# 06_negative_control.R  --  negative controls and sensitivity analyses
# -----------------------------------------------------------------------------
# (1) random gene sets matched to the real sets' sizes -> null distribution of
#     the adjusted hazard ratio (is 1.31 "just any 800 genes"?)
# (2) survival-label permutation -> HR distribution centred on 1
# (3) sensitivity matrix: alternative scoring, stage IV excluded, stage-only
#     adjustment, winsorised scores, RFS instead of OS (GSE31210)
# Output: results/tables/06_*.csv
# =============================================================================
suppressPackageStartupMessages({ library(survival); library(data.table) })
source("R/lib_surv.R"); source("R/lib_score.R")
PROJ <- normalizePath(".", winslash = "/")
OUT  <- file.path(PROJ, "results", "tables")
S    <- readRDS(file.path(PROJ, "data", "processed", "scores.rds"))
CO   <- readRDS(file.path(PROJ, "data", "processed", "cohorts.rds"))
AXES <- list(
  proliferation = c("E2F_TARGETS", "G2_M_CHECKPOINT", "MYC_TARGETS_V1", "MYC_TARGETS_V2"),
  hypoxia       = c("HYPOXIA", "GLYCOLYSIS"))
COV <- c("stage", "age_z", "sex")

read_gmt <- function(f) {
  ln <- readLines(f); ln <- ln[nzchar(ln)]
  sets <- lapply(strsplit(ln, "\t"), function(v) { v <- v[nzchar(v)]; unique(v[-c(1, 2)]) })
  names(sets) <- gsub("[^A-Z0-9]+", "_", toupper(vapply(strsplit(ln, "\t"), `[`, character(1), 1)))
  sets
}
HALL <- read_gmt(file.path(PROJ, "data", "genesets", "Hallmark_Enrichr.gmt"))

ranks_of <- function(expr) {
  N <- nrow(expr)
  apply(expr, 2, function(x) { o <- order(x, decreasing = TRUE)
                               r <- integer(N); r[o] <- seq_len(N); r })
}

prep <- function(id) {
  d <- S[[id]]
  d$stage <- factor(d$stage, levels = c("I", "II", "III", "IV"))
  d$age_z <- z(d$age)
  d <- d[stats::complete.cases(d[, c("os_time", "os_event", "stage", "age_z", "sex")]), ]
  e <- CO[[id]]$expr[, d$sample, drop = FALSE]
  list(d = d, expr = e, ranks = ranks_of(e))
}

safe_adj_hr <- function(...) tryCatch(suppressWarnings(adj_hr(...)),
                                        error = function(e) c(HR = NA_real_, p = NA_real_))
adj_hr <- function(d, score, covars = COV) {
  if (length(score) != nrow(d)) stop("length mismatch: ", length(score), " vs ", nrow(d))
  if (sum(!is.na(score)) < 20) stop("score has too few non-missing values")
  d$score <- score
  covars <- covars[covars %in% names(d)]
  f <- as.formula(paste("Surv(os_time, os_event) ~ score +", paste(covars, collapse = " + ")))
  m <- coxph(f, data = d)
  c(HR = unname(exp(coef(m))[1]), p = unname(summary(m)$coef[1, "Pr(>|z|)"]))
}

set.seed(2026)
null_rows <- list()
for (id in c("TCGA-LUAD", "GSE31210")) {
  P <- prep(id)
  message("null distribution for ", id, " (", nrow(P$d), " samples)")
  for (ax in names(AXES)) {
    sizes <- vapply(AXES[[ax]], function(s) sum(rownames(P$expr) %in% HALL[[s]]), integer(1))
    B <- if (id == "TCGA-LUAD") 500 else 200
    hr <- numeric(B)
    for (b in seq_len(B)) {
      idx <- lapply(sizes, function(k) sample(nrow(P$expr), k))
      sc  <- ssgsea_from_ranks(P$ranks, idx)
      hr[b] <- adj_hr(P$d, z(rowMeans(apply(sc, 2, z))))[["HR"]]
    }
    obs <- adj_hr(P$d, z(P$d[[ax]]))[["HR"]]
    null_rows[[paste(id, ax)]] <- data.frame(
      cohort = id, axis = ax, B = B, observed_HR = obs,
      null_median = median(hr), null_q95 = quantile(hr, .95), null_q99 = quantile(hr, .99),
      empirical_p = (sum(hr >= obs) + 1) / (B + 1),
      n_random_above = sum(hr >= obs))
    cat(sprintf("  %-14s observed %.2f | null median %.2f, q95 %.2f, q99 %.2f | p=%.4f\n",
                ax, obs, median(hr), quantile(hr, .95), quantile(hr, .99),
                (sum(hr >= obs) + 1) / (B + 1)))
    fwrite(data.frame(cohort = id, axis = ax, HR = hr),
           file.path(OUT, paste0("06_null_hr_", gsub("[^A-Za-z0-9]", "_", id), "_", ax, ".csv")))
  }
}
null_tab <- rbindlist(null_rows)
fwrite(null_tab, file.path(OUT, "06_random_set_null.csv"))

# ---- label permutation ------------------------------------------------------
perm_rows <- list(); perm_draws <- list()
P <- prep("TCGA-LUAD")
for (ax in names(AXES)) {
  obs <- adj_hr(P$d, z(P$d[[ax]]))[["HR"]]
  B <- 1000; hr <- numeric(B)
  for (b in seq_len(B)) {
    dd <- P$d
    dd[, c("os_time", "os_event")] <- dd[sample(nrow(dd)), c("os_time", "os_event")]
    hr[b] <- adj_hr(dd, z(dd[[ax]]))[["HR"]]
  }
  perm_draws[[ax]] <- data.frame(axis = ax, HR = hr)
  perm_rows[[ax]] <- data.frame(axis = ax, observed_HR = obs, B = B,
                                perm_median = median(hr),
                                perm_q025 = quantile(hr, .025), perm_q975 = quantile(hr, .975),
                                empirical_p = (sum(hr >= obs) + 1) / (B + 1))
  cat(sprintf("permutation %-14s observed %.2f | perm 95%% range %.2f-%.2f | p=%.4f\n",
              ax, obs, quantile(hr, .025), quantile(hr, .975),
              (sum(hr >= obs) + 1) / (B + 1)))
}
perm_tab <- rbindlist(perm_rows)
if (length(perm_draws)) fwrite(rbindlist(perm_draws), file.path(OUT, "06_perm_hr.csv"))
fwrite(perm_tab, file.path(OUT, "06_permutation_test.csv"))
saveRDS(list(null = null_tab), file.path(PROJ, "results", "objects", "06_null_tab.rds"))


# ---- sensitivity matrix -----------------------------------------------------
sens <- list()
add <- function(id, ax, variant, val) {
  sens[[length(sens) + 1]] <<- data.frame(cohort = id, axis = ax, variant = variant,
    HR = round(val[["HR"]], 3), p = signif(val[["p"]], 3))
}
mean_z_score <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  if (!length(g)) return(rep(NA_real_, ncol(expr)))
  m <- apply(expr[g, , drop = FALSE], 1, z)   # samples x genes
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  z(rowMeans(m))                              # one value per sample
}

for (id in c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081")) {
  P <- prep(id); d <- P$d
  for (ax in names(AXES)) {
    add(id, ax, "primary (adjusted, all stages)", safe_adj_hr(d, z(d[[ax]])))
    add(id, ax, "stage-only adjustment", safe_adj_hr(d, z(d[[ax]]), c("stage")))
    no4 <- which(d$stage != "IV")
    if (length(no4) > 50) add(id, ax, "stage IV excluded", safe_adj_hr(d[no4, ], z(d[[ax]])[no4]))
    add(id, ax, "mean z-score instead of ssGSEA",
        safe_adj_hr(d, mean_z_score(P$expr[, d$sample, drop = FALSE], unlist(HALL[AXES[[ax]]]))))
    q <- quantile(d[[ax]], c(.01, .99))
    w <- pmin(pmax(d[[ax]], q[1]), q[2])
    add(id, ax, "winsorised at 1%/99%", safe_adj_hr(d, w))
  }
  if (!is.null(d$rfs_time)) {
    dd <- d[!is.na(d$rfs_time) & !is.na(d$rfs_event) & d$rfs_time > 0, ]
    for (ax in names(AXES)) {
      dd$os_time <- dd$rfs_time; dd$os_event <- dd$rfs_event
      add(id, ax, "recurrence-free survival", safe_adj_hr(dd, z(dd[[ax]])))
    }
  }
}
sens <- rbindlist(sens)
fwrite(sens, file.path(OUT, "06_sensitivity.csv"))
cat("\n===== sensitivity: HR range per axis/variant (TCGA) =====\n")
print(sens[cohort == "TCGA-LUAD"], row.names = FALSE)
cat("\n===== sensitivity: does the direction hold everywhere? =====\n")
print(sens[, .(n = .N, n_HR_above_1 = sum(HR > 1)), by = axis], row.names = FALSE)

saveRDS(list(null = null_tab, perm = perm_tab, sens = sens),
        file.path(PROJ, "results", "objects", "06_negative.rds"))
cat("\nDONE 06\n")
