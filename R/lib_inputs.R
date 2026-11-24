# Input manifest verification and gene-symbol mapping helpers.
input_manifest <- function() {
  inputs <- read.csv("data/checksums.csv", colClasses = c(
    path = "character", md5 = "character", bytes = "numeric", scope = "character", source = "character"
  ))
  if (Sys.getenv("LUAD_CORE_ONLY") == "1") inputs <- inputs[inputs$scope == "core", ]
  inputs
}

verify_input_files <- function(inputs = input_manifest()) {
  missing <- inputs$path[!file.exists(inputs$path)]
  if (length(missing))
    stop("Missing pinned inputs. Run Rscript scripts/fetch_inputs.R. Missing: ", paste(missing, collapse = ", "))
  actual <- unname(tools::md5sum(inputs$path))
  bad <- is.na(actual) | actual != inputs$md5 | file.size(inputs$path) != inputs$bytes
  if (any(bad))
    stop("Input checksum or size differs: ", paste(inputs$path[bad], collapse = ", "))
  invisible(inputs)
}

ensembl_symbols <- function(ids) {
  ids <- sub("\\..*$", "", ids)
  symbols <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(ids), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  as.character(symbols[ids])
}

collapse_gene_expression <- function(mat, symbols, transport = FALSE) {
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (!nrow(mat)) stop("No expression features mapped to gene symbols")

  if (transport) {
    sample_means <- rowsum(mat, group = symbols, reorder = FALSE)
    probe_counts <- tabulate(match(symbols, rownames(sample_means)), nbins = nrow(sample_means))
    sample_means <- sample_means / probe_counts
  }

  sorted <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  picked <- sorted[!duplicated(symbols[sorted])]
  expr <- mat[picked, , drop = FALSE]
  rownames(expr) <- symbols[picked]
  if (transport) list(expr = expr, transport_expr = sample_means) else expr
}
