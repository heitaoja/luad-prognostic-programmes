# Check package dependencies, verify pinned inputs and record session metadata.
source("R/lib_inputs.R")
source("R/lib_score.R")
packages <- c("survival", "data.table", "ggplot2", "patchwork", "openxlsx",
              "org.Hs.eg.db", "hgu133plus2.db", "AnnotationDbi", "R.utils")
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(available)) stop("Missing R packages: ", paste(packages[!available], collapse = ", "))
for (p in c("data/processed", "results/objects", "results/tables",
            "results/figures/main", "results/figures/supplementary", "results/logs"))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
inputs <- verify_input_files()
if (!capabilities("cairo")) stop("Cairo graphics is required for the PDF figures")
variants <- c(
  profile = if (Sys.getenv("LUAD_CORE_ONLY") == "1") "core" else "full",
  outcome = "overall survival",
  unit = "month",
  score = "tie-aware sample rank ssGSEA",
  alpha = 0.25,
  coverage_min = 0.8,
  synthesis = "REML with modified Knapp-Hartung",
  multiplicity = "Holm across the four programme tests",
  validation_bootstrap = 1000L,
  prediction_horizon_months = 36,
  random_sets_tcga = 500L,
  random_sets_geo = 200L,
  paired_label_permutations = 1000L,
  seed = 20261124
)
parameters <- c(variants, vapply(PROGRAMMES, paste, "", collapse = "; "))
write.csv(data.frame(parameter = names(parameters), value = unname(parameters)),
          "results/tables/00_params.csv", row.names = FALSE)
versions <- vapply(packages, function(name) as.character(packageVersion(name)), "")
writeLines(c(R.version.string, paste("Platform", R.version$platform),
             "", paste(names(versions), versions, sep = " "), "",
             trimws(capture.output(sessionInfo()), which = "right")), "docs/sessionInfo.txt")
cat("Verified", nrow(inputs), "pinned inputs and recorded", length(versions), "package versions\n")
