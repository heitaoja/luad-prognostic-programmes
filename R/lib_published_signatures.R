# Benchmark published gene sets with the same association and transport designs.
PUBLISHED_SIGNATURE_SIZES <- c(CCP_31 = 31L, meta_PCNA = 131L)
PUBLISHED_ASSOCIATION_SCORES <- c("ssGSEA", "sample_centred_mean_z")
PUBLISHED_TRANSPORT_SCORES <- c("ssGSEA", "mean_log_z")
PUBLISHED_ASSOCIATION_CACHE <- "data/processed/published_signature_association.rds"

read_published_signatures <- function() {
  rows <- fread("config/published_signatures.csv")
  ids <- names(PUBLISHED_SIGNATURE_SIZES)
  stopifnot(identical(names(rows), c("signature_id", "display_name", "citation", "gene_symbol")),
            !anyNA(rows), setequal(rows$signature_id, ids), all(nzchar(rows$gene_symbol)),
            all(nzchar(rows$display_name)), all(nzchar(rows$citation)),
            !anyDuplicated(rows[, .(signature_id, gene_symbol = toupper(gene_symbol))]))
  sets <- setNames(lapply(ids, function(id) rows[signature_id == id, gene_symbol]), ids)
  stopifnot(identical(lengths(sets), PUBLISHED_SIGNATURE_SIZES))
  sets
}

# HGNC previous symbols and unambiguous aliases. Venet Table S1 Entrez 146909
# identifies KIF18B. C18orf24 and C13orf3 identify different SKA subunits.
PUBLISHED_GENE_ALIASES <- list(
  AHSP = "ERAF", BUB1B = "BUBR1", CDC45 = c("CDC45L", "CDC45L2"), CDK1 = "CDC2",
  CENPF = "HCP-1", CENPU = "MLF1IP", DDX39A = "DDX39", `H3-3A` = c("H3F3A", "H3F3"),
  HJURP = "DKFZP762E1312", KIF18B = "LOC146909", MIS18A = c("C21ORF45", "C21ORF46"),
  NSD2 = "WHSC1", ORC6 = "ORC6L", PCLAF = "KIAA0101", SKA1 = "C18ORF24", SKA3 = "C13ORF3",
  SRSF2 = "SFRS2", TSPO2 = "BZRPL1"
)

match_published_signatures <- function(genes, signatures) {
  aliases <- PUBLISHED_GENE_ALIASES
  approved <- setNames(rep(names(aliases), lengths(aliases)), unlist(aliases, use.names = FALSE))
  symbols <- toupper(genes)
  measured <- lapply(signatures, function(set) {
    set <- toupper(set)
    old <- set %in% names(approved)
    set[old] <- approved[set[old]]
    indices <- vapply(set, function(symbol) {
      hits <- match(c(symbol, aliases[[symbol]]), symbols)
      hits[!is.na(hits)][1L]
    }, integer(1))
    found <- genes[na.omit(indices)]
    stopifnot(!anyDuplicated(found))
    found
  })
  minimum <- ifelse(names(signatures) == "meta_PCNA", 0.85, 0.8)
  stopifnot(all(lengths(measured) / lengths(signatures) >= minimum))
  measured
}

published_coverage <- function(signatures, measured) {
  data.table(signature_id = names(signatures), n_in_set = lengths(signatures),
             n_found = lengths(measured), coverage = lengths(measured) / lengths(signatures))
}

published_association_cohort <- function(expr, d, id, signatures) {
  expr <- expr[, d$sample, drop = FALSE]
  measured <- match_published_signatures(rownames(expr), signatures)
  medians <- apply(expr, 2L, median)
  centred_scores <- vapply(measured, function(genes) {
    e <- t(sweep(expr[genes, , drop = FALSE], 2L, medians, "-"))
    keep <- apply(e, 2L, sd) > 0
    rowMeans(apply(e[, keep, drop = FALSE], 2L, standardise))
  }, numeric(ncol(expr)))
  raw <- list(ssGSEA = ssgsea_fast(expr, measured), sample_centred_mean_z = centred_scores)
  patients <- lapply(raw, function(matrix) {
    scores <- lapply(as.data.frame(matrix), standardise)
    external_ipd_cohort(cbind(d, scores), id, c("proliferation", names(signatures)))
  })
  list(patients = patients, coverage = published_coverage(signatures, measured))
}

published_cox_rows <- function(d, id = NA_character_, stratify = character(0)) {
  clinical <- c("stage", "age_z", "sex")
  rbindlist(lapply(names(PUBLISHED_SIGNATURE_SIZES), function(signature) {
    rows <- rbindlist(list(
      cox_row(d, signature, clinical, "clinical_plus_marker", id, stratify),
      cox_row(d, "proliferation", c(clinical, signature), paste0("composite_given_", signature), id, stratify),
      cox_row(d, signature, c(clinical, "proliferation"), paste0(signature, "_given_composite"), id, stratify)
    ))
    rows[, signature_id := signature]
    rows[, .(cohort, signature_id, model, exposure, n, events, HR, lo, hi, p)]
  }))
}

write_published_benchmark <- function(cohorts, tab) {
  scopes <- list(validation = c("GSE31210", "GSE37745", "GSE50081"))
  if ("GSE72094" %in% names(cohorts)) scopes$extended <- c(scopes$validation, "GSE72094")
  if ("CPTAC-LUAD" %in% names(cohorts)) scopes$all_external_k5 <- c(scopes$extended, "CPTAC-LUAD")
  out <- rbindlist(lapply(PUBLISHED_ASSOCIATION_SCORES, function(score_type) {
    patients <- lapply(cohorts, function(x) x$patients[[score_type]])
    individual <- rbindlist(lapply(names(patients), function(id) {
      d <- patients[[id]]
      rows <- published_cox_rows(d, id)
      coverage <- cohorts[[id]]$coverage
      coverage$spearman <- vapply(coverage$signature_id, function(s)
        cor(d$proliferation, d[[s]], method = "spearman"), numeric(1))
      rows <- merge(rows, coverage, by = "signature_id", sort = FALSE)
      rows[, `:=`(scope = "cohort", method = "Cox", k = 1L)]
      rows
    }))
    pooled <- rbindlist(lapply(names(scopes), function(scope_name) {
      ids <- scopes[[scope_name]]
      rows <- individual[cohort %in% ids]
      two_stage <- rbindlist(lapply(c("REML_mKH", "DL_normal"), function(m) {
        rows[, {
          est <- pool_estimates(.SD, m)
          .(method = m, k = est$k, n = est$n_total, events = est$events_total,
            HR = est$pooled_HR, lo = est$lo, hi = est$hi, p = est$p)
        }, by = .(signature_id, model, exposure)]
      }))
      ipd <- published_cox_rows(do.call(rbind, patients[ids]), stratify = "cohort")
      ipd[, `:=`(method = "IPD_strata_cohort", k = length(ids))]
      combined <- rbindlist(list(two_stage, ipd), use.names = TRUE, fill = TRUE)
      combined[, scope := scope_name]
      combined
    }), use.names = TRUE)
    rows <- rbindlist(list(individual, pooled), use.names = TRUE, fill = TRUE)
    rows[, scoring := score_type]
    rows
  }))
  setcolorder(out, c("cohort", "scope", "method", "scoring", "signature_id", "model",
                     "exposure", "n_in_set", "n_found", "coverage", "spearman", "k",
                     "n", "events", "HR", "lo", "hi", "p"))
  fwrite(out, file.path(tab, "03e_published_signature_benchmark.csv"))
  cat("Wrote", nrow(out), "published signature association estimates\n")
}

published_transport_input <- function(components, expression, universe, measured) {
  stopifnot(all(universe %in% rownames(expression)), all(components$clin$sample %in% colnames(expression)))
  expr <- expression[universe, components$clin$sample, drop = FALSE]
  ranks <- ssgsea_fast(expr, measured)
  colnames(ranks) <- paste0(names(measured), "_ssGSEA")
  # Standardise the within-sample mean log expression using the TCGA reference at fitting.
  means <- vapply(measured, function(genes) colMeans(expr[genes, , drop = FALSE]), numeric(ncol(expr)))
  dimnames(means) <- list(colnames(expr), paste0(names(measured), "_mean_log_z"))
  components$extra_scores <- cbind(ranks, means)
  components
}

published_transport_results <- function(components, expression, universe, signatures) {
  measured <- match_published_signatures(universe, signatures)
  variants <- CJ(signature_id = names(signatures), scoring = PUBLISHED_TRANSPORT_SCORES, sorted = FALSE)
  variants[, model := paste(signature_id, scoring, sep = "_")]
  inputs <- setNames(lapply(names(components), function(id)
    published_transport_input(components[[id]], expression[[id]], universe, measured)), names(components))
  markers <- c(TRANSPORT_MARKERS["composite"], setNames(variants$model, variants$model))
  frozen <- fit_frozen_transport(inputs[["TCGA-LUAD"]], markers = markers)
  coverage <- published_coverage(signatures, measured)
  variants <- merge(variants, coverage, by = "signature_id", sort = FALSE)
  set.seed(20261124)
  rbindlist(lapply(setdiff(names(inputs), "TCGA-LUAD"), function(id) {
    test <- evaluate_frozen_transport(frozen, inputs[[id]], id, comparisons = variants$model)
    metrics <- test$evaluation
    contrasts <- test$contrasts[match(variants$model, test$contrasts$comparison)]
    published <- metrics[match(variants$model, metrics$model)]
    published[, model := NULL]
    rows <- cbind(variants, n_genes_universe = length(universe), published,
                  clinical_c_index = metrics[model == "clinical", c_index],
                  composite_c_index = metrics[model == "composite", c_index],
                  contrasts[, .(composite_delta_c = delta_c, composite_delta_lo = delta_lo,
                                composite_delta_hi = delta_hi)])
    setcolorder(rows, c("cohort", "signature_id", "scoring", "model", "n_genes_universe"))
    rows
  }))
}
