# =============================================================================
# 02_score.R  --  ssGSEA scoring of four pre-specified tumour programmes
# -----------------------------------------------------------------------------
# Input : data/processed/cohorts.rds, data/genesets/Hallmark_Enrichr.gmt
# Output: data/processed/scores.rds   (cohort -> clin + 4 axis scores + balance)
#         results/tables/02_gene_coverage.csv
#         results/tables/02_axis_definition.csv
#
# Pre-specified axes (docs/analysis_plan.md).  Each axis score is the mean of the
# within-cohort z-scored ssGSEA values of its member gene sets, then z-scored
# again, so every axis has mean 0 / SD 1 inside every cohort.
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })
PROJ <- normalizePath(".", winslash = "/")
PDIR <- file.path(PROJ, "data", "processed")
OUT  <- file.path(PROJ, "results", "tables")
dir.create(OUT, TRUE, showWarnings = FALSE)
cohorts <- readRDS(file.path(PDIR, "cohorts.rds"))

read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  sets <- lapply(strsplit(ln, "\t"), function(v) { v <- v[nzchar(v)]; unique(v[-c(1, 2)]) })
  names(sets) <- gsub("[^A-Z0-9]+", "_", toupper(vapply(strsplit(ln, "\t"),
                                                        `[`, character(1), 1)))
  sets
}
HALL <- read_gmt(file.path(PROJ, "data", "genesets", "Hallmark_Enrichr.gmt"))
message("Hallmark gene sets loaded: ", length(HALL))

AXES <- list(
  proliferation = c("E2F_TARGETS", "G2_M_CHECKPOINT", "MYC_TARGETS_V1", "MYC_TARGETS_V2"),
  hypoxia       = c("HYPOXIA", "GLYCOLYSIS"),
  emt_stroma    = c("EPITHELIAL_MESENCHYMAL_TRANSITION", "TGF_BETA_SIGNALING",
                    "ANGIOGENESIS"),
  immune        = c("INFLAMMATORY_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
                    "INTERFERON_ALPHA_RESPONSE", "ALLOGRAFT_REJECTION", "COMPLEMENT")
)
AXES <- lapply(AXES, function(x) x[x %in% names(HALL)])
stopifnot(all(vapply(AXES, length, integer(1)) >= 2))
ALL_SETS <- unique(unlist(AXES))
message("axes: ", paste(names(AXES), vapply(AXES, length, integer(1)), sep = "=",
                        collapse = ", "))

# ------------------------------------------------------------------ ssGSEA ---
# ssGSEA (Barbie 2009) as in GSVA: rank genes per sample, weighted walk with
# |rank|^alpha, score = sum of the running sum, normalised by the gene count.
score_sets <- function(expr, sets, alpha = 0.25, min_genes = 5) {
  genes <- rownames(expr)
  idx <- lapply(sets, function(g) which(genes %in% g))
  keep <- vapply(idx, length, integer(1)) >= min_genes
  if (!any(keep)) stop("no gene set with >= ", min_genes, " genes in this cohort")
  idx <- idx[keep]
  ns <- ncol(expr); N <- nrow(expr)
  res <- matrix(NA_real_, ns, length(idx), dimnames = list(colnames(expr), names(idx)))
  for (j in seq_len(ns)) {
    o <- order(expr[, j], decreasing = TRUE)
    hm <- vapply(idx, function(i) { v <- logical(N); v[i] <- TRUE; v[o] }, logical(N))
    for (k in seq_along(idx)) {
      hit <- hm[, k]; kk <- sum(hit); w <- numeric(N)
      w[hit] <- (N - which(hit) + 1L)^alpha; w <- w / sum(w)
      st <- rep(-1 / (N - kk), N); st[hit] <- w[hit]
      res[j, k] <- sum(cumsum(st)) / N
    }
  }
  res
}
zscore <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

# --------------------------------------------------------------- per cohort --
out <- list(); cov_tab <- list(); def_tab <- list()
for (id in names(cohorts)) {
  e <- cohorts[[id]]$expr
  message("\n== ", id, " (", ncol(e), " samples x ", nrow(e), " genes)")

  sc <- score_sets(e, HALL[ALL_SETS])
  scz <- apply(sc, 2, zscore)

  cov_tab[[id]] <- data.frame(
    cohort = id, set = colnames(sc),
    n_in_set = vapply(ALL_SETS, function(s) length(HALL[[s]]), integer(1))[colnames(sc)],
    n_found = vapply(colnames(sc), function(s) sum(rownames(e) %in% HALL[[s]]), integer(1)),
    stringsAsFactors = FALSE)

  d <- cohorts[[id]]$clin
  for (a in names(AXES)) {
    sets_a <- intersect(AXES[[a]], colnames(scz))
    d[[a]] <- zscore(rowMeans(scz[, sets_a, drop = FALSE]))
  }
  # the "balance" hypothesis that was pre-specified first, kept as a secondary
  # exploratory analysis (see docs/changelog.md)
  d$balance <- zscore(scale(d$proliferation)[, 1] - scale(d$immune)[, 1])
  d$cohort <- id
  out[[id]] <- d

  def_tab[[id]] <- data.frame(cohort = id, axis = names(AXES),
      member_sets = vapply(AXES, paste, character(1), collapse = "; "),
      stringsAsFactors = FALSE)

  # sanity: every real set must sit far in the tail of a matched random-set null
  set.seed(1)
  null_sets <- lapply(1:50, function(i) sample(rownames(e), 200))
  names(null_sets) <- paste0("NULL", 1:50)
  nr <- score_sets(e, null_sets)
  z_null <- (colMeans(sc) - mean(colMeans(nr))) / sd(colMeans(nr))
  cat(sprintf("   null z of real sets: %s\n", paste(sprintf("%+.1f", z_null), collapse = " ")))
  if ("MKI67" %in% rownames(e))
    cat(sprintf("   cor(proliferation, MKI67) = %.2f\n",
        cor(d$proliferation, as.numeric(e["MKI67", d$sample]), use = "complete.obs")))
}

saveRDS(out, file.path(PDIR, "scores.rds"))
fwrite(do.call(rbind, cov_tab), file.path(OUT, "02_gene_coverage.csv"))
fwrite(def_tab[[1]], file.path(OUT, "02_axis_definition.csv"))
for (id in names(out)) {
  fwrite(out[[id]][, c("sample", "cohort", "os_time", "os_event", "age", "sex", "stage",
                       names(AXES), "balance")],
         file.path(OUT, paste0("02_scores_", gsub("[^A-Za-z0-9]", "_", id), ".csv")))
}

cat("\n===== axis score summary (mean / sd within cohort) =====\n")
for (id in names(out)) {
  d <- out[[id]]
  cat(sprintf("%-10s n=%3d | %s\n", id, nrow(d),
      paste(sprintf("%s=%+.2f/%.2f", names(AXES),
                    sapply(names(AXES), function(a) mean(d[[a]])),
                    sapply(names(AXES), function(a) sd(d[[a]]))), collapse = "  ")))
}
cat("\ncorrelation matrix of axes (TCGA):\n")
print(round(cor(as.matrix(out[["TCGA-LUAD"]][, names(AXES)])), 2))
cat("\nsaved:", file.path(PDIR, "scores.rds"), "\n")
