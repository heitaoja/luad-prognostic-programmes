# =============================================================================
# 00_setup.R  --  project design constants and shared helpers
# -----------------------------------------------------------------------------
# Single place to read the design of this project.  every threshold, cohort and
# gene-set definition used by the pipeline is recorded here and written to
# results/tables/00_params.csv.
# =============================================================================

# ---- paths --------------------------------------------------------------
PROJ    <- normalizePath(".", winslash = "/")
RAW_DIR <- file.path(PROJ, "data", "raw")
PDIR    <- file.path(PROJ, "data", "processed")
CFG_DIR <- file.path(PROJ, "config")
TAB_DIR <- file.path(PROJ, "results", "tables")
FIG_DIR <- file.path(PROJ, "results", "figures")
OBJ_DIR <- file.path(PROJ, "results", "objects")
LOG_DIR <- file.path(PROJ, "results", "logs")
DOC_DIR <- file.path(PROJ, "docs")
for (d in c(PDIR, TAB_DIR, FIG_DIR, OBJ_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEED <- 20261101
set.seed(SEED)

# ---- project design -----------------------------------------------------
PARAMS <- list(
  seed = SEED,

  cohorts = list(
    discovery  = "TCGA-LUAD",
    validation = c("GSE31210", "GSE37745", "GSE50081"),
    excluded   = "GSE72094 (custom HuRSTA array, no annotation package available)"
  ),

  # four pre-specified Hallmark programmes (docs/analysis_plan.md)
  axes = list(
    proliferation = c("E2F_TARGETS", "G2_M_CHECKPOINT", "MYC_TARGETS_V1", "MYC_TARGETS_V2"),
    hypoxia       = c("HYPOXIA", "GLYCOLYSIS"),
    emt_stroma    = c("EPITHELIAL_MESENCHYMAL_TRANSITION", "TGF_BETA_SIGNALING",
                      "ANGIOGENESIS"),
    immune        = c("INFLAMMATORY_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
                      "INTERFERON_ALPHA_RESPONSE", "ALLOGRAFT_REJECTION", "COMPLEMENT")
  ),

  scoring = list(
    method    = "ssgsea",      # closed-form implementation in lib_score.R
    alpha     = 0.25,          # GSVA default
    min_genes = 5L,
    normalize = "z_within_cohort",
    aggregate = "mean_of_z"    # mean of member-set z scores, then z again
  ),

  model = list(
    endpoint   = "OS",
    time_unit  = "month",
    covariates = c("stage", "age_z", "sex"),
    exposure   = "continuous_per_sd",
    cutoff     = "median",     # KM visualisation only
    infer_on   = "continuous"
  ),

  inference = list(
    meta_model = "random_effects_DerSimonian-Laird",
    ci_level   = 0.95
  ),

  negative_control = list(
    n_random_sets_tcga = 500L,
    n_random_sets_geo  = 200L,
    n_permutations     = 1000L,
    null_set_size      = 200L
  ),

  incremental = list(
    horizon_months = 36,
    bootstrap_B    = 200L,
    dca_lower      = 0.05,
    dca_upper      = 0.60
  ),

  figures = list(width_mm = 180, dpi = 300),

  assert = list(max_na_pct = 0.30, min_cohort_n = 50L, epv_min = 10L)
)

# ---- helpers ------------------------------------------------------------
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", sep = "", file = file.path(LOG_DIR, "pipeline.log"), append = TRUE)
  invisible(msg)
}

assert_nonempty <- function(x, what) {
  if (is.null(x) || NROW(x) == 0L) stop("EMPTY RESULT: ", what, call. = FALSE)
  invisible(x)
}

assert_cols <- function(x, cols, what = deparse(substitute(x))) {
  miss <- setdiff(cols, colnames(x))
  if (length(miss))
    stop("MISSING COLUMNS in ", what, ": ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(x)
}

write_params <- function() {
  p <- unlist(PARAMS)
  utils::write.csv(data.frame(parameter = names(p), value = unname(p)),
                   file.path(TAB_DIR, "00_params.csv"), row.names = FALSE)
  invisible(p)
}

log_msg("design loaded: 4 programmes x ",
        length(PARAMS$cohorts$validation) + 1L, " cohorts (",
        PARAMS$cohorts$discovery, " + ",
        paste(PARAMS$cohorts$validation, collapse = ", "), ")")
write_params()
