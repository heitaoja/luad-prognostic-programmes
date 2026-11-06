`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# lib_score.R  --  ssGSEA with a closed-form running sum (fast path)
# -----------------------------------------------------------------------------
# The classic implementation walks the ranked list (O(N) per set).  Because the
# enrichment score is a sum over positions it has a closed form:
#
#   ES = sum_j step_j * (N - j + 1)
#      = sum_hits w_j (N-j+1)  -  (1/(N-k)) * sum_nonhits (N-j+1)
#
# with w_j = (N-j+1)^alpha / Z.  This makes scoring O(k) per set, which is what
# makes the 1000-random-set null in 06_negative_control.R affordable.
# Identical to the walk version up to floating point (checked in 06).
# =============================================================================

ssgsea_fast <- function(expr, sets, alpha = 0.25, min_genes = 5) {
  sets <- lapply(sets, function(g) unique(g))
  genes <- rownames(expr)
  pos   <- lapply(sets, function(g) match(g, genes))
  pos   <- lapply(pos, function(p) p[!is.na(p)])
  keep  <- vapply(pos, length, integer(1)) >= min_genes
  pos   <- pos[keep]
  N <- nrow(expr); ns <- ncol(expr)
  out <- matrix(NA_real_, ns, length(pos), dimnames = list(colnames(expr), names(pos)))
  tot <- N * (N + 1) / 2

  for (j in seq_len(ns)) {
    o <- order(expr[, j], decreasing = TRUE)
    rank_of_gene <- integer(N); rank_of_gene[o] <- seq_len(N)   # 1 = highest
    for (k in seq_along(pos)) {
      r  <- rank_of_gene[pos[[k]]]          # positions in the ranked list
      kk <- length(r)
      v  <- (N - r + 1)                     # rank weight base
      w  <- v^alpha / sum(v^alpha)
      out[j, k] <- sum(w * v) - (tot - sum(v)) / (N - kk)
    }
  }
  out / N
}

# ---------------------------------------------------------- rank-based path --
# Precompute each sample's gene ranking once, then score any number of gene sets
# without re-ordering.  Used for the random-set null distributions.
ranks_of <- function(expr) {
  N <- nrow(expr)
  apply(expr, 2, function(x) { o <- order(x, decreasing = TRUE)
                               r <- integer(N); r[o] <- seq_len(N); r })
}

ssgsea_from_ranks <- function(rank_mat, idx_list, alpha = 0.25, min_genes = 5) {
  N <- nrow(rank_mat); tot <- N * (N + 1) / 2
  out <- matrix(NA_real_, ncol(rank_mat), length(idx_list),
                dimnames = list(colnames(rank_mat), names(idx_list) %||% seq_along(idx_list)))
  for (k in seq_along(idx_list)) {
    idx <- idx_list[[k]]
    if (length(idx) < min_genes) next
    R  <- rank_mat[idx, , drop = FALSE]
    v  <- N - R + 1
    va <- v^alpha
    w  <- va / rep(colSums(va), each = nrow(va))
    out[, k] <- (colSums(w * v) - (tot - colSums(v)) / (N - length(idx))) / N
  }
  out
}
