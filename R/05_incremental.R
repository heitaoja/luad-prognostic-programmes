# Fit in TCGA and evaluate the same models in each external array cohort.
source("R/lib_surv.R")
source("R/lib_score.R")
source("R/lib_transport.R")
input <- readRDS("data/processed/transport_components.rds")
tab <- "results/tables"
frozen <- fit_frozen_transport(input[["TCGA-LUAD"]])
fwrite(transport_coefficients(frozen), file.path(tab, "05_frozen_coefficients.csv"))
fwrite(transport_score_reference(frozen), file.path(tab, "05_frozen_score_reference.csv"))

set.seed(20261124)
external <- setdiff(names(input), "TCGA-LUAD")
results <- lapply(external, function(id) evaluate_frozen_transport(frozen, input[[id]], id))
evaluation <- rbindlist(lapply(results, function(x) x$evaluation))
fwrite(evaluation, file.path(tab, "05_incremental_value.csv"))
fwrite(rbindlist(lapply(results, function(x) x$contrasts)),
       file.path(tab, "05_simple_marker_contrasts.csv"))
fwrite(rbindlist(lapply(results, function(x) x$calibration)), file.path(tab, "05_calibration.csv"))
print(evaluation[model %in% c("clinical", "composite"),
                 .(cohort, model, c_index, delta_c, delta_lo, delta_hi, brier_36)])
