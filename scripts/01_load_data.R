# Load, harmonise and audit all expression and clinical cohorts.
# Inputs: raw TCGA, GEO and annotation files (log2 scale, survival in months).
# Outputs: data/processed/cohorts.rds, transport_expression.rds, tcga_alternate.rds, 01_cohort_audit.csv.

suppressPackageStartupMessages({
  library(data.table); library(AnnotationDbi); library(org.Hs.eg.db); library(survival)
})
source("R/lib_inputs.R")

PROJ <- normalizePath(".", winslash = "/")
OUT  <- file.path(PROJ, "results", "tables");    dir.create(OUT, TRUE, showWarnings = FALSE)
PDIR <- file.path(PROJ, "data", "processed");    dir.create(PDIR, TRUE, showWarnings = FALSE)
RAW  <- file.path(PROJ, "data", "raw")
DAY_PER_MONTH <- 30.4375

# Helpers
`%||%` <- function(a, b) if (is.null(a) || !length(a) || all(is.na(a))) b else a

read_geo <- function(file) {
  con <- gzfile(file, "r"); i <- 0L; hdr <- character(0)
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE); i <- i + 1L
    if (!length(ln) || grepl("^!series_matrix_table_begin", ln) || i > 600L) break
    hdr <- c(hdr, ln)
  }
  close(con)
  tab <- suppressWarnings(fread(file, skip = i, header = TRUE, sep = "\t",
               quote = "\"", data.table = FALSE, check.names = FALSE,
               showProgress = FALSE))
  tab <- tab[!grepl("^!", tab[[1]]), , drop = FALSE]      # drop table_end footer
  ids <- tab[[1]]
  mat <- as.matrix(tab[, -1, drop = FALSE]); rownames(mat) <- ids
  storage.mode(mat) <- "double"

  ch_lines <- grep("^!Sample_characteristics_ch1", hdr, value = TRUE)
  gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", hdr, value = TRUE)[1],
                                "\t")[[1]])[-1]
  rec <- lapply(seq_along(gsm), function(k) list())
  for (ln in ch_lines) {
    vals <- gsub('"', "", strsplit(ln, "\t")[[1]])[-1]
    for (k in seq_along(vals)) {
      kv <- strsplit(vals[k], ": ", fixed = TRUE)[[1]]
      if (length(kv) >= 2) {
        key <- trimws(kv[1]); val <- trimws(paste(kv[-1], collapse = ": "))
        if (nzchar(val) && val != "NA") rec[[k]][[key]] <- val
      }
    }
  }
  keys <- unique(unlist(lapply(rec, names)))
  pheno <- as.data.frame(setNames(lapply(keys, function(k)
    vapply(rec, function(r) r[[k]] %||% NA_character_, character(1))), keys),
    stringsAsFactors = FALSE, check.names = FALSE)
  pheno$gsm <- gsm
  list(expr = mat, pheno = pheno)
}

to_symbol <- function(mat, platform, transport = FALSE) {
  if (platform == "ensembl") {
    sym <- ensembl_symbols(rownames(mat))
  } else if (platform == "hgu133plus2") {
    sym <- suppressMessages(AnnotationDbi::mapIds(hgu133plus2.db::hgu133plus2.db,
             keys = rownames(mat), keytype = "PROBEID", column = "SYMBOL",
             multiVals = "first"))
    sym <- as.character(sym)
  } else if (platform == "hursta") {
    annotation <- fread(file.path(RAW, "geo", "GPL15048_hursta_annotation.txt"),
                        select = c("ProbeID", "NumGenes", "Symbol"), na.strings = "---")
    annotation <- unique(annotation[NumGenes == 1L & !is.na(Symbol)])
    unambiguous <- annotation[, .N, by = ProbeID][N == 1L, ProbeID]
    annotation <- annotation[ProbeID %in% unambiguous]
    sym <- annotation$Symbol[match(rownames(mat), annotation$ProbeID)]
  } else stop("unknown platform: ", platform)

  collapse_gene_expression(mat, sym, transport)
}

# "Stage IA" / "IA" / "1a" / "2b" / "IV" -> I, II, III, IV
stage_from <- function(x) {
  s <- toupper(gsub("[^0-9IVXivx]", "", as.character(x)))
  out <- rep(NA_character_, length(s))
  # longest match first; only fill entries that are still NA so a later,
  # shorter pattern cannot overwrite a longer one ("IIIA" must stay III)
  fill <- function(pattern, value) {
    idx <- is.na(out) & grepl(pattern, s)
    out[idx] <<- value
  }
  fill("^IV", "IV");  fill("^III", "III"); fill("^II", "II")
  fill("^I($|[^VI])", "I")
  fill("^4", "IV");   fill("^3", "III");   fill("^2", "II");   fill("^1", "I")
  factor(out, levels = c("I", "II", "III", "IV"))
}

sex_from <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[grepl("^m", s)] <- "male"; out[grepl("^f", s)] <- "female"
  out[s %in% c("1", "男")] <- "male"; out[s %in% c("2", "女")] <- "female"
  out
}

# "yes/no" or "dead/alive" -> 1/0
event_from <- function(x, pos = c("yes", "dead", "deceased", "1", "relapsed")) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(s))
  out[s %in% pos] <- 1L
  out[s %in% c("no", "alive", "0", "not relapsed", "censored")] <- 0L
  out
}

num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", as.character(x))))

audit <- list()
audit_add <- function(id, expr, clin, unit) {
  km_median <- function(event) {
    fit <- survfit(Surv(clin$os_time, event) ~ 1)
    round(unname(summary(fit)$table["median"]), 1)
  }
  complete <- complete.cases(clin[, c("stage", "age", "sex")])
  stopifnot(identical(colnames(expr), clin$sample), !anyDuplicated(clin$sample),
            all(clin$os_event %in% 0:1), all(is.finite(clin$os_time) & clin$os_time > 0))
  audit[[id]] <<- data.frame(
    cohort = id, n_samples = ncol(expr), n_genes = nrow(expr),
    n_os = sum(!is.na(clin$os_time) & !is.na(clin$os_event)),
    n_events = sum(clin$os_event == 1, na.rm = TRUE),
    median_os_mo = km_median(clin$os_event),
    median_followup_mo = km_median(1L - clin$os_event),
    n_complete = sum(complete), events_complete = sum(clin$os_event[complete]),
    events_36m = sum(clin$os_event == 1L & clin$os_time <= 36),
    at_risk_36m = sum(clin$os_time >= 36),
    n_stage = sum(!is.na(clin$stage)), pct_stage_known = round(100 * mean(!is.na(clin$stage))),
    n_age = sum(!is.na(clin$age)), n_sex = sum(!is.na(clin$sex)),
    raw_time_unit = unit,
    expr_range = sprintf("%.1f..%.1f", min(expr), max(expr)),
    stringsAsFactors = FALSE)
}

cohorts <- list()

# TCGA-LUAD
message("== TCGA-LUAD")
tc <- file.path(RAW, "tcga")
ex <- fread(file.path(tc, "TCGA-LUAD.star_counts.tsv.gz"), data.table = FALSE, showProgress = FALSE)
rownames(ex) <- ex[[1]]; ex <- as.matrix(ex[, -1, drop = FALSE]); storage.mode(ex) <- "double"
ex <- ex[, grepl("-01[ABC]?$", colnames(ex)), drop = FALSE]        # primary tumours
n_primary <- ncol(ex)
surv <- fread(file.path(tc, "TCGA-LUAD.survival.tsv.gz"), data.table = FALSE)
cl   <- fread(file.path(tc, "TCGA-LUAD.clinical.tsv.gz"), data.table = FALSE)

match_tcga <- function(sample, source) {
  exact <- match(sample, source)
  missing <- is.na(exact)
  priority <- order(!grepl("-01[ABC]?$", source), source)
  exact[missing] <- priority[match(substr(sample[missing], 1L, 12L),
                                   substr(source[priority], 1L, 12L))]
  exact
}
i_s <- match_tcga(colnames(ex), surv$sample)
i_c <- match_tcga(colnames(ex), cl$sample)

tcga_clin <- data.frame(
  sample   = colnames(ex),
  patient  = substr(colnames(ex), 1L, 12L),
  os_time  = as.numeric(surv$OS.time[i_s]) / DAY_PER_MONTH,   # raw unit: days
  os_event = as.integer(surv$OS[i_s]),
  age      = num(cl$age_at_index.demographic[i_c]),
  sex      = sex_from(cl$gender.demographic[i_c]),
  stage    = stage_from(cl$ajcc_pathologic_stage.diagnoses[i_c]),
  pack_years = num(cl$pack_years_smoked.exposures[i_c]),
  stringsAsFactors = FALSE)
tcga_clin <- tcga_clin[!is.na(tcga_clin$os_time) & !is.na(tcga_clin$os_event) &
                       tcga_clin$os_time > 0, ]
linked <- tcga_clin[order(tcga_clin$sample), ]
repeat_os <- unique(linked[, c("patient", "os_time", "os_event")])
if (anyDuplicated(repeat_os$patient))
  stop("Repeated TCGA samples disagree on patient-level survival")
tcga_clin <- linked[!duplicated(linked$patient), ]
tcga_alt <- linked[!duplicated(linked$patient, fromLast = TRUE), ]
tcga_alt <- tcga_alt[match(tcga_clin$patient, tcga_alt$patient), ]
tcga <- to_symbol(ex[, tcga_clin$sample, drop = FALSE], "ensembl", transport = TRUE)
tcga_alt_expr <- to_symbol(ex[, tcga_alt$sample, drop = FALSE], "ensembl")
fwrite(data.frame(primary_specimens = n_primary, specimens_with_os = nrow(linked),
                  participants_with_os = nrow(tcga_clin),
                  participants_with_repeated_specimens = sum(table(linked$patient) > 1L),
                  deaths = sum(tcga_clin$os_event)), file.path(OUT, "01_tcga_participant_flow.csv"))
rm(ex)
cohorts[["TCGA-LUAD"]] <- c(tcga, list(clin = tcga_clin,
  meta = list(source = "GDC Data Portal (via gdc-hub mirror)",
              platform = "Illumina RNA-seq, STAR counts", scale = "log2(count+1)",
              role = "discovery / training", raw_time_unit = "days")))
audit_add("TCGA-LUAD", tcga$expr, tcga_clin, "days")

# GSE31210
message("== GSE31210")
g <- read_geo(file.path(RAW, "geo", "GSE31210_series_matrix.txt.gz"))
mapped <- to_symbol(log2(g$expr + 1), "hgu133plus2", transport = TRUE)
p <- g$pheno
excl <- p[["exclude for prognosis analysis due to incomplete resection or adjuvant therapy"]]
clin <- data.frame(
  sample   = p$gsm,
  os_time  = num(p[["days before death/censor"]]) / DAY_PER_MONTH,
  os_event = event_from(p$death),
  rfs_time = num(p[["days before relapse/censor"]]) / DAY_PER_MONTH,
  rfs_event= event_from(p$relapse, c("yes", "dead", "relapsed")),
  age      = num(p[["age (years)"]]),
  sex      = sex_from(p$gender),
  stage    = stage_from(p[["pathological stage"]]),
  smoking  = p[["smoking status"]],
  driver   = p[["gene alteration status"]],
  stringsAsFactors = FALSE)
keep <- grepl("primary lung tumor", p$tissue, ignore.case = TRUE) &
        !is.na(clin$os_time) & !is.na(clin$os_event) & clin$os_time > 0 &
        !(grepl("^exclude", excl %||% "none", ignore.case = TRUE))
clin <- clin[keep, ]
mapped <- lapply(mapped, function(x) x[, clin$sample, drop = FALSE])
cohorts[["GSE31210"]] <- c(mapped, list(clin = clin,
  meta = list(source = "GEO", platform = "Affymetrix HG-U133 Plus 2.0 (GPL570)",
              scale = "log2(intensity+1)", role = "external validation",
              raw_time_unit = "days")))
audit_add("GSE31210", mapped$expr, clin, "days")

# GSE37745
message("== GSE37745")
g <- read_geo(file.path(RAW, "geo", "GSE37745_series_matrix.txt.gz"))
mapped <- to_symbol(g$expr, "hgu133plus2", transport = TRUE)
p <- g$pheno
clin <- data.frame(
  sample   = p$gsm,
  os_time  = num(p[["days to determined death status"]]) / DAY_PER_MONTH,
  os_event = event_from(p$dead),
  age      = num(p$age),
  sex      = sex_from(p$gender),
  stage    = stage_from(p[["tumor stage"]]),
  histology= p$histology,
  stringsAsFactors = FALSE)
keep <- !is.na(clin$os_time) & !is.na(clin$os_event) & clin$os_time > 0 &
        grepl("^adeno", clin$histology, ignore.case = TRUE)
clin <- clin[keep, ]
mapped <- lapply(mapped, function(x) x[, clin$sample, drop = FALSE])
cohorts[["GSE37745"]] <- c(mapped, list(clin = clin,
  meta = list(source = "GEO", platform = "Affymetrix HG-U133 Plus 2.0 (GPL570)",
              scale = "log2", role = "external validation", raw_time_unit = "days",
              selection = "histology == adeno")))
audit_add("GSE37745", mapped$expr, clin, "days")

# GSE50081
message("== GSE50081")
g <- read_geo(file.path(RAW, "geo", "GSE50081_series_matrix.txt.gz"))
mapped <- to_symbol(g$expr, "hgu133plus2", transport = TRUE)
p <- g$pheno
clin <- data.frame(
  sample   = p$gsm,
  os_time  = num(p[["survival time"]]) * 12,          # raw unit: YEARS
  os_event = event_from(p$status),
  age      = num(p$age),
  sex      = sex_from(p$Sex),
  stage    = stage_from(p$Stage),
  histology= p$histology,
  smoking  = p$smoking,
  stringsAsFactors = FALSE)
keep <- !is.na(clin$os_time) & !is.na(clin$os_event) & clin$os_time > 0 &
        grepl("^adenocarcinoma$", trimws(tolower(clin$histology)))
clin <- clin[keep, ]
mapped <- lapply(mapped, function(x) x[, clin$sample, drop = FALSE])
cohorts[["GSE50081"]] <- c(mapped, list(clin = clin,
  meta = list(source = "GEO", platform = "Affymetrix HG-U133 Plus 2.0 (GPL570)",
              scale = "log2", role = "external validation", raw_time_unit = "years",
              selection = "histology == adenocarcinoma")))
audit_add("GSE50081", mapped$expr, clin, "years")

# GSE72094
extra_files <- file.path(RAW, "geo", c("GSE72094_series_matrix.txt.gz",
                                     "GPL15048_hursta_annotation.txt"))
if (Sys.getenv("LUAD_CORE_ONLY") != "1") {
  if (!all(file.exists(extra_files)))
    stop("GSE72094 requires both files in docs/data_provenance.md. Set LUAD_CORE_ONLY=1 for a core-only replay.")
  message("== GSE72094")
  g <- read_geo(extra_files[1L])
  p <- g$pheno
  clin <- data.frame(sample = p$gsm, patient = p$patient_id,
                     os_time = num(p$survival_time_in_days) / DAY_PER_MONTH,
                     os_event = event_from(p$vital_status),
                     age = num(p$age_at_diagnosis), sex = sex_from(p$gender),
                     stage = stage_from(p$Stage), smoking = p$smoking_status,
                     kras = p$kras_status, egfr = p$egfr_status,
                     stk11 = p$stk11_status, tp53 = p$tp53_status,
                     stringsAsFactors = FALSE)
  stopifnot(!anyNA(clin$patient), !anyDuplicated(clin$patient))
  keep <- !is.na(clin$os_time) & !is.na(clin$os_event) & clin$os_time > 0
  clin <- clin[keep, ]
  mapped <- to_symbol(g$expr[, clin$sample, drop = FALSE], "hursta", transport = TRUE)
  stopifnot("MKI67" %in% rownames(mapped$expr))
  cohorts[["GSE72094"]] <- c(mapped, list(clin = clin,
    meta = list(source = "GEO (Schabath et al.)", platform = "HuRSTA (GPL15048)",
                scale = "log2", role = "external extension", raw_time_unit = "days")))
  audit_add("GSE72094", mapped$expr, clin, "days")
}

# Save processed objects and cohort audit table
aud <- do.call(rbind, audit); rownames(aud) <- NULL
fwrite(aud, file.path(OUT, "01_cohort_audit.csv"))
transport_expression <- lapply(cohorts, function(x) x$transport_expr)
saveRDS(transport_expression, file.path(PDIR, "transport_expression.rds"))
cohorts <- lapply(cohorts, function(x) x[c("expr", "clin", "meta")])
saveRDS(cohorts, file.path(PDIR, "cohorts.rds"))
saveRDS(list(expr = tcga_alt_expr, clin = tcga_alt), file.path(PDIR, "tcga_alternate.rds"))

cat("\n================ cohort audit ================\n")
print(aud, row.names = FALSE)
cat("\ngenes per cohort:", paste(sapply(cohorts, function(z) nrow(z$expr)), collapse = ", "), "\n")
cat("saved:", file.path(PDIR, "cohorts.rds"), "\n")
