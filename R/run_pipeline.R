# Cross-platform entry point. Run from the repository root.
steps <- c("00_setup", "01_load_data", "02_score", "03_tcga_primary",
           "03b_mutual_adjustment", "03c_benchmark", "04_external_meta",
           "05_incremental", "05b_rnaseq_validation", "06_negative_control",
           "07_figures", "08_tables", "99_verify")
selected <- commandArgs(trailingOnly = TRUE)
if (length(selected)) {
  hits <- lapply(selected, function(prefix) steps[startsWith(steps, prefix)])
  if (any(lengths(hits) == 0L)) stop("Unknown step prefix: ", paste(selected[lengths(hits) == 0L], collapse = ", "))
  steps <- steps[steps %in% unlist(hits)]
}
if (Sys.getenv("LUAD_CORE_ONLY") == "1") steps <- setdiff(steps, "05b_rnaseq_validation")
for (d in c("data/processed", "results/logs", "results/objects", "results/tables", "results/figures"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
start <- Sys.time()
for (step in steps) {
  log <- file.path("results/logs", paste0(step, ".log"))
  cat("Running", step, "\n")
  status <- system2(rscript, c("--vanilla", shQuote(file.path("R", paste0(step, ".R")))),
                    stdout = log, stderr = log)
  if (status != 0L) {
    cat(tail(readLines(log, warn = FALSE), 30), sep = "\n")
    stop("Step failed: ", step, ". See ", log, call. = FALSE)
  }
}
cat("Completed", length(steps), "steps in", round(as.numeric(difftime(Sys.time(), start, units = "secs"))), "seconds\n")
