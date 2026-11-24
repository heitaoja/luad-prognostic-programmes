# Outcome permutations, structure-matched random programmes, and robustness.
source("R/lib_surv.R")
source("R/lib_score.R")
scores <- readRDS("data/processed/scores.rds")
cohorts <- readRDS("data/processed/cohorts.rds")
tab <- "results/tables"
sets <- programme_sets()
axes <- PROGRAMMES[c("proliferation", "hypoxia")]
covariates <- c("stage", "age_z", "sex")
set.seed(20261124)

prepare <- function(id) {
  d <- scores[[id]]
  d$age_z <- z(d$age)
  expr <- cohorts[[id]]$expr[, d$sample, drop = FALSE]
  list(d = d, expr = expr)
}

fit_hr <- function(d, score, covars = covariates) {
  d$score <- score
  row <- cox_row(d, "score", covars)
  unlist(row[c("n", "events", "HR", "lo", "hi", "p")], use.names = TRUE)
}

# Replace each unique target gene exactly once. Reusing the replacement where
# the target gene occurred preserves all component sizes and overlaps.
null_drawer <- function(gene_names, abundance_bin, components) {
  members <- lapply(components, function(g) which(gene_names %in% g))
  target <- unique(unlist(members, use.names = FALSE))
  by_bin <- split(target, abundance_bin[target])
  available <- setdiff(seq_along(gene_names), target)
  pools <- split(available, abundance_bin[available])
  if (any(lengths(by_bin) > lengths(pools[names(by_bin)])))
    stop("Abundance matching has too few control genes")
  idx <- lapply(members, match, table = target)
  target_pos <- lapply(by_bin, match, table = target)
  function() {
    replacement <- integer(length(target))
    for (bin in names(by_bin)) {
      pool <- pools[[bin]]
      replacement[target_pos[[bin]]] <- pool[sample.int(length(pool), length(by_bin[[bin]]))]
    }
    lapply(idx, function(m) replacement[m])
  }
}

null_summary <- list()
null_cohorts <- intersect(c("TCGA-LUAD", "GSE31210", "GSE72094"), names(scores))
for (id in null_cohorts) {
  p <- prepare(id)
  ranks <- ranks_of(p$expr)
  abundance <- rowMeans(ranks)
  bin <- pmin(10L, ceiling(10 * rank(abundance, ties.method = "average") / length(abundance)))
  B <- if (id == "TCGA-LUAD") 500L else 200L
  for (a in names(axes)) {
    draw <- null_drawer(rownames(p$expr), bin, sets[axes[[a]]])
    hr <- replicate(B, {
      sc <- ssgsea_from_ranks(ranks, draw())
      combined <- standardise(rowMeans(apply(sc, 2L, standardise)))
      fit_hr(p$d, combined)[["HR"]]
    })
    obs <- fit_hr(p$d, p$d[[a]])[["HR"]]
    null_summary[[length(null_summary) + 1L]] <- data.frame(
      cohort = id, axis = a, B = B, observed_HR = obs, null_median = median(hr),
      null_q95 = quantile(hr, .95), null_q99 = quantile(hr, .99),
      empirical_p = (sum(hr >= obs) + 1) / (B + 1), n_random_above = sum(hr >= obs))
    fwrite(data.frame(cohort = id, axis = a, HR = hr),
           file.path(tab, paste0("06_null_hr_", gsub("-", "_", id), "_", a, ".csv")))
    cat(id, a, "matched-null upper-tail probability", signif((sum(hr >= obs) + 1)/(B + 1), 3), "\n")
  }
}
fwrite(rbindlist(null_summary), file.path(tab, "06_random_set_null.csv"))

d <- prepare("TCGA-LUAD")$d
d <- d[complete.cases(d[, c("os_time", "os_event", covariates)]), ]
perms <- lapply(names(axes), function(a) {
  obs <- fit_hr(d, d[[a]])[["HR"]]
  B <- 1000L
  hr <- replicate(B, {
    t <- d
    t[, c("os_time", "os_event")] <- d[sample.int(nrow(d)), c("os_time", "os_event")]
    fit_hr(t, t[[a]])[["HR"]]
  })
  list(summary = data.frame(axis = a, observed_HR = obs, B = B, perm_median = median(hr),
                            perm_q025 = quantile(hr, .025), perm_q975 = quantile(hr, .975),
                            empirical_p = (sum(hr >= obs) + 1) / (B + 1)),
       draws = data.frame(axis = a, HR = hr))
})
fwrite(rbindlist(lapply(perms, function(d) d$summary)), file.path(tab, "06_permutation_test.csv"))
fwrite(rbindlist(lapply(perms, function(d) d$draws)), file.path(tab, "06_perm_hr.csv"))

sens <- list()
for (id in names(scores)) {
  p <- prepare(id)
  d <- p$d
  centred_expr <- sweep(p$expr, 2L, apply(p$expr, 2L, median), "-")
  for (a in names(axes)) {
    add <- function(variant, result) {
      sens[[length(sens) + 1L]] <<- cbind(data.frame(cohort = id, axis = a, variant = variant),
                                         as.data.frame(as.list(result)))
    }
    add("Primary clinical adjustment", fit_hr(d, d[[a]]))
    add("Stage only", fit_hr(d, d[[a]], "stage"))
    early <- !is.na(d$stage) & d$stage != "IV"
    add("Exclude stage IV", fit_hr(d[early, ], d[[a]][early]))
    genes <- intersect(rownames(p$expr), unlist(sets[axes[[a]]]))
    e <- t(p$expr[genes, , drop = FALSE])
    keep <- apply(e, 2L, sd) > 0
    add("Mean gene-wise z-score", fit_hr(d, standardise(rowMeans(apply(e[, keep, drop = FALSE], 2L, standardise)))))
    centred <- t(centred_expr[genes, , drop = FALSE])
    keep_centred <- apply(centred, 2L, sd) > 0
    centred_score <- standardise(rowMeans(apply(centred[, keep_centred, drop = FALSE], 2L, standardise)))
    add("Mean gene-wise z-score (sample-centred)", fit_hr(d, centred_score))
    lim <- quantile(d[[a]], c(.01, .99))
    add("Winsorise at 1% and 99%", fit_hr(d, standardise(pmin(pmax(d[[a]], lim[1L]), lim[2L]))))
    if (!is.null(d$rfs_time)) {
      keep_rfs <- !is.na(d$rfs_time) & d$rfs_time > 0 & !is.na(d$rfs_event)
      r <- d[keep_rfs, ]
      r$os_time <- r$rfs_time
      r$os_event <- r$rfs_event
      add("Recurrence endpoint (exploratory)", fit_hr(r, r[[a]]))
    }
  }
}
fwrite(rbindlist(sens), file.path(tab, "06_sensitivity.csv"))
