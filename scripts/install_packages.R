# Install CRAN and Bioconductor dependencies required by the pipeline.
cran_pkgs <- c("survival", "data.table", "ggplot2", "patchwork", "openxlsx", "R.utils", "BiocManager")
missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("AnnotationDbi", "org.Hs.eg.db", "hgu133plus2.db")
missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
