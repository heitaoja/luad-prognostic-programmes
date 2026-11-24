# Shared gene-set definitions and a tie-aware, sample-rank enrichment score.
PROGRAMMES <- list(
  proliferation = c("E2F_TARGETS", "G2_M_CHECKPOINT", "MYC_TARGETS_V1", "MYC_TARGETS_V2"),
  hypoxia = c("HYPOXIA", "GLYCOLYSIS"),
  emt_stroma = c("EPITHELIAL_MESENCHYMAL_TRANSITION", "TGF_BETA_SIGNALING", "ANGIOGENESIS"),
  immune = c("INFLAMMATORY_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
             "INTERFERON_ALPHA_RESPONSE", "ALLOGRAFT_REJECTION", "COMPLEMENT")
)

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  fields <- strsplit(lines[nzchar(lines)], "\t", fixed = TRUE)
  if (any(lengths(fields) < 3L)) stop("Malformed GMT row in ", path)
  ids <- gsub("[^A-Z0-9]+", "_", toupper(vapply(fields, function(row) row[1L], "")))
  if (anyDuplicated(ids)) stop("Duplicate gene set names in ", path)
  sets <- lapply(fields, function(row) {
    genes <- trimws(row[-c(1L, 2L)])
    unique(genes[nzchar(genes)])
  })
  setNames(sets, ids)
}

programme_sets <- function(path = "data/genesets/Hallmark_Enrichr.gmt") {
  sets <- read_gmt(path)
  required <- unlist(PROGRAMMES, use.names = FALSE)
  missing <- setdiff(required, names(sets))
  if (length(missing)) stop("Missing programme sets: ", paste(missing, collapse = ", "))
  sets[required]
}

gene_set_coverage <- function(symbols, sets, cohort, universe) {
  data.frame(cohort = cohort, universe = universe, set = names(sets),
             n_in_set = lengths(sets),
             n_found = vapply(sets, function(genes) sum(genes %in% symbols), integer(1)))
}

check_gene_set_coverage <- function(coverage) {
  bad <- coverage$n_found / coverage$n_in_set < 0.8
  if (any(bad))
    stop("Less than 80% measured: ",
         paste(paste(coverage$cohort[bad], coverage$set[bad]), collapse = ", "))
}

z_parameters <- function(x) {
  spread <- stats::sd(x)
  if (!is.finite(spread) || spread <= 0) stop("Cannot standardise a constant or incomplete score")
  c(center = mean(x), scale = spread)
}

standardise <- function(x, parameters = z_parameters(x)) {
  (x - parameters[["center"]]) / parameters[["scale"]]
}

programme_scores <- function(components, reference = NULL) {
  if (is.null(reference))
    reference <- list(components = lapply(as.data.frame(components), z_parameters))
  z <- vapply(names(reference$components), function(s)
    standardise(components[, s], reference$components[[s]]), numeric(nrow(components)))
  axes <- vapply(PROGRAMMES, function(sets)
    rowMeans(z[, sets, drop = FALSE]), numeric(nrow(components)))
  if (is.null(reference$programmes))
    reference$programmes <- lapply(as.data.frame(axes), z_parameters)
  for (axis in names(PROGRAMMES))
    axes[, axis] <- standardise(axes[, axis], reference$programmes[[axis]])
  scores <- cbind(axes, e2f = z[, "E2F_TARGETS"], g2m = z[, "G2_M_CHECKPOINT"])
  list(scores = as.data.frame(scores), reference = reference)
}

ranks_of <- function(expr) {
  if (nrow(expr) < 2L || any(!is.finite(expr))) stop("Non-finite or empty expression matrix")
  apply(expr, 2L, function(x) rank(-x, ties.method = "average"))
}

ssgsea_fast <- function(expr, sets, alpha = 0.25, min_genes = 5L) {
  idx <- lapply(sets, function(genes) which(rownames(expr) %in% genes))
  ssgsea_from_ranks(ranks_of(expr), idx, alpha, min_genes)
}

transport_components <- function(expr, clin, genes, sets) {
  stopifnot(identical(colnames(expr), clin$sample), !anyDuplicated(genes),
            all(genes %in% rownames(expr)), "MKI67" %in% genes)
  ranks <- ranks_of(expr[genes, , drop = FALSE])
  indices <- lapply(sets, function(members) which(genes %in% members))
  list(clin = clin, scores = ssgsea_from_ranks(ranks, indices),
       mki67_quantile = (nrow(ranks) - ranks[match("MKI67", genes), ]) / (nrow(ranks) - 1))
}

# The integral of the rank walk is the sum of each step times the number
# of positions it contributes to. Midranks give tied genes equal contributions.
ssgsea_from_ranks <- function(rank_mat, idx_list, alpha = 0.25, min_genes = 5L) {
  n <- nrow(rank_mat)
  total_ranks <- n * (n + 1) / 2
  labels <- names(idx_list)
  if (is.null(labels)) labels <- as.character(seq_along(idx_list))
  sizes <- lengths(idx_list)
  invalid <- sizes < min_genes | sizes >= n
  if (any(invalid)) stop("Invalid measured gene set sizes: ", paste(labels[invalid], collapse = ", "))
  out <- matrix(NA_real_, ncol(rank_mat), length(idx_list),
                dimnames = list(colnames(rank_mat), labels))
  for (k in seq_along(idx_list)) {
    idx <- idx_list[[k]]
    v <- n - rank_mat[idx, , drop = FALSE] + 1
    weights <- v ^ alpha
    hit <- colSums(weights * v) / colSums(weights)
    miss <- (total_ranks - colSums(v)) / (n - length(idx))
    out[, k] <- (hit - miss) / n
  }
  out
}
