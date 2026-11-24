# Score association cohorts and prepare sample-independent transport features.
suppressPackageStartupMessages(library(data.table))
source("R/lib_score.R")
cohorts <- readRDS("data/processed/cohorts.rds")
sets <- programme_sets()
tab <- "results/tables"

associate <- function(cohort) {
  components <- ssgsea_fast(cohort$expr, sets)
  fit <- programme_scores(components)
  scores <- fit$scores
  scores$MKI67 <- standardise(cohort$expr["MKI67", ])
  scores$balance <- standardise(scores$proliferation - scores$immune)
  stopifnot(identical(rownames(components), cohort$clin$sample))
  list(data = cbind(cohort$clin, scores), reference = fit$reference)
}

association <- lapply(cohorts, associate)
scores <- lapply(names(cohorts), function(id) {
  d <- association[[id]]$data
  d$cohort <- id
  d
})
names(scores) <- names(cohorts)
raw_coverage <- lapply(names(cohorts), function(id)
  gene_set_coverage(rownames(cohorts[[id]]$expr), sets, id, "cohort"))
check_gene_set_coverage(rbindlist(raw_coverage))

alt <- readRDS("data/processed/tcga_alternate.rds")
alt_data <- associate(alt)$data
base <- scores[["TCGA-LUAD"]]
changed <- base$sample != alt$clin$sample
paired <- programme_scores(ssgsea_fast(alt$expr[, changed, drop = FALSE], sets),
                           association[["TCGA-LUAD"]]$reference)$scores
fwrite(rbindlist(lapply(names(PROGRAMMES), function(axis) {
  comparable <- paired[[axis]]
  primary <- base[[axis]][changed]
  data.frame(axis = axis, paired_samples = sum(changed),
             spearman = cor(primary, comparable, method = "spearman"),
             median_absolute_score_difference = median(abs(primary - comparable)))
})), file.path(tab, "02_paired_specimen_scores.csv"))
saveRDS(alt_data, "data/processed/tcga_alternate_scores.rds")

# Use the intersection of annotated measurement platforms, without outcome
# information. Every participant is ranked on this identical gene universe.
transport_expression <- readRDS("data/processed/transport_expression.rds")
# radix sorting is locale-independent, so the gene-universe file is byte-identical
# across platforms and collation settings (default sort() depends on LC_COLLATE)
common <- sort(Reduce(intersect, lapply(transport_expression, rownames)), method = "radix")
stopifnot("MKI67" %in% common)
fwrite(data.frame(gene_symbol = common), file.path(tab, "02_transport_gene_universe.csv"))
transport <- list()
shared_coverage <- list()
compatibility <- list()
for (id in names(cohorts)) {
  x <- cohorts[[id]]
  expr <- transport_expression[[id]]
  present <- common %in% rownames(expr)
  compatibility[[id]] <- data.frame(cohort = id, n_genes_required = length(common),
                                    n_genes_present = sum(present),
                                    transport_eligible = all(present))
  if (!all(present)) next
  shared_coverage[[id]] <- gene_set_coverage(common, sets, id, "fixed_transport")
  check_gene_set_coverage(shared_coverage[[id]])
  transport[[id]] <- transport_components(expr, x$clin, common, sets)
}
fwrite(rbindlist(compatibility), file.path(tab, "02_transport_coverage.csv"))
fwrite(rbindlist(c(raw_coverage, shared_coverage)), file.path(tab, "02_gene_coverage.csv"))
fwrite(data.frame(axis = names(PROGRAMMES),
                  member_sets = vapply(PROGRAMMES, paste, "", collapse = "; ")),
       file.path(tab, "02_axis_definition.csv"))
patient_csv <- "data/processed/patient_scores"
dir.create(patient_csv, showWarnings = FALSE)
for (id in names(scores)) {
  d <- scores[[id]]
  fields <- c("sample", "cohort", "os_time", "os_event", "age", "sex", "stage",
              names(PROGRAMMES), "MKI67", "e2f", "g2m", "balance")
  fwrite(d[, fields], file.path(patient_csv, paste0("02_scores_", gsub("[^A-Za-z0-9]", "_", id), ".csv")))
}
saveRDS(scores, "data/processed/scores.rds")
saveRDS(transport, "data/processed/transport_components.rds")
cat("Scored", length(cohorts), "cohorts,", length(common), "genes in the fixed transport universe\n")
