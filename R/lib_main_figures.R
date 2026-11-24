# Multipanel journal figures assembled from the release tables.
source("R/lib_score.R")
source("R/lib_inputs.R")

MAIN_INK <- "#223447"
MAIN_BLUE <- "#265B80"
MAIN_TEAL <- "#23766E"
MAIN_RUST <- "#A64F36"
MAIN_SLATE <- "#626F80"
MAIN_PURPLE <- "#73568B"
MAIN_GRID <- "#E3E9EF"
MAIN_AXES <- c(proliferation = "Proliferation", hypoxia = "Hypoxia / glycolysis",
               emt_stroma = "EMT / stroma", immune = "Immune / inflammation")
MAIN_COHORT_COLOURS <- c(`TCGA-LUAD` = MAIN_BLUE, GSE31210 = MAIN_TEAL, GSE37745 = MAIN_PURPLE,
                         GSE50081 = "#95701D", GSE72094 = "#AA5360", `CPTAC-LUAD` = "#46576A")
MAIN_MARKERS <- c(composite = "Composite", MKI67 = "MKI67", E2F = "E2F", G2M = "G2/M",
                  CCP_31 = "CCP-31", meta_PCNA = "meta-PCNA")

main_theme <- function() {
  theme_minimal(base_size = 9) +
    theme(text = element_text(colour = MAIN_INK), axis.text = element_text(colour = MAIN_INK, size = 8),
          axis.title = element_text(size = 8.5), panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = MAIN_GRID, linewidth = .25),
          panel.grid.major.y = element_blank(), panel.spacing = grid::unit(2.6, "mm"),
          strip.background = element_rect(fill = "#EDF2F6", colour = NA),
          strip.text = element_text(face = "bold", size = 8.2, margin = margin(4, 3, 4, 3)),
          plot.title.position = "plot", plot.title = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(colour = MAIN_SLATE, size = 8, margin = margin(b = 4)),
          plot.caption = element_text(colour = MAIN_SLATE, size = 7.4, hjust = 0),
          plot.margin = margin(6, 7, 5, 6), legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 8), legend.key.width = grid::unit(4, "mm"))
}

main_axis <- function(axis) factor(unname(MAIN_AXES[as.character(axis)]), levels = unname(MAIN_AXES))

save_main_figure <- function(plot, name, width, height) {
  tagged <- plot + plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(5, 6, 5, 6)))
  tagged <- tagged & theme(plot.tag = element_text(face = "bold", size = 11.5, colour = MAIN_INK))
  save_figure(tagged, name, width, height)
}

main_forest <- function(plot, position = "identity", breaks = c(.25, .5, 1, 2, 4, 8, 16)) {
  plot + geom_vline(xintercept = 1, colour = "#8796A3", linetype = 2, linewidth = .36) +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y", position = position, linewidth = .42) +
    geom_point(position = position, size = 1.7) +
    scale_x_log10(breaks = breaks, labels = function(x) sprintf("%g", x),
                  expand = expansion(mult = c(.07, .09))) +
    labs(x = "Adjusted hazard ratio per SD (95% CI)", y = NULL) + main_theme()
}

cptac_programme_scores <- function() {
  path <- "data/raw/cptac"
  rna <- fread(file.path(path, "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"))
  survival <- fread(file.path(path, "LUAD_survival.txt"))
  cases <- names(rna)[-1L]
  stopifnot(!anyDuplicated(cases), !anyDuplicated(survival$case_id), setequal(cases, survival$case_id))
  survival <- survival[match(cases, case_id)]
  eligible <- with(survival, is.finite(OS_days) & OS_days > 0 & !is.na(OS_event))
  # As in 05b, collapse duplicate symbols before retaining cases with valid survival.
  expr <- collapse_gene_expression(as.matrix(rna[, ..cases]), ensembl_symbols(rna$idx))
  stopifnot(sum(eligible) == read_table("05b_rnaseq_audit")$with_os)
  scores <- programme_scores(ssgsea_fast(expr[, eligible, drop = FALSE], programme_sets()))$scores
  cbind(cohort = "CPTAC-LUAD", scores)
}

figure_cohorts_and_scores <- function(audit, scores) {
  audit$cohort <- factor(audit$cohort, levels = rev(audit$cohort))
  measures <- c("Participants with survival data", "Median follow-up (months)")
  counts <- copy(audit)[, measure := factor(measures[1L], levels = measures)]
  followup <- copy(audit)[, measure := factor(measures[2L], levels = measures)]
  panel_a <- ggplot(mapping = aes(y = cohort)) +
    geom_col(data = counts, aes(x = participants), width = .56, fill = "#D9E4ED") +
    geom_col(data = counts, aes(x = complete), width = .56, fill = MAIN_BLUE) +
    geom_text(data = counts, aes(x = participants, label = paste(complete, participants, sep = " / ")),
              hjust = -.13, colour = MAIN_INK, size = 2.9) +
    geom_segment(data = followup, aes(x = 0, xend = followup, yend = cohort), colour = "#C8D7E0", linewidth = .6) +
    geom_point(data = followup, aes(x = followup), colour = MAIN_TEAL, size = 2) +
    geom_text(data = followup, aes(x = followup, label = sprintf("%.1f", followup)),
              hjust = -.4, colour = MAIN_INK, size = 2.9) +
    facet_wrap(~ measure, scales = "free_x", nrow = 1) +
    scale_x_continuous(expand = expansion(mult = c(.015, .24))) +
    labs(title = "Cohort size and follow-up",
         subtitle = "Counts show complete cases / participants. Follow-up uses reverse Kaplan-Meier.",
         x = NULL, y = NULL) + main_theme()

  fields <- c("cohort", names(MAIN_AXES))
  density <- rbindlist(lapply(scores, function(d)
    melt(as.data.table(d)[, ..fields], id.vars = "cohort", variable.name = "programme", value.name = "score")))
  density[, `:=`(cohort = factor(cohort, levels = rev(levels(audit$cohort))), programme = main_axis(programme))]
  panel_b <- ggplot(density, aes(score, colour = cohort)) +
    geom_density(linewidth = .56, adjust = 1.05, key_glyph = "path") +
    facet_wrap(~ programme, ncol = 2) + scale_colour_manual(values = MAIN_COHORT_COLOURS) +
    guides(colour = guide_legend(nrow = 1)) +
    labs(title = "Hallmark programme scores across cohorts", x = "Cohort-standardised score (SD)",
         y = "Density") + main_theme()
  save_main_figure(panel_a / panel_b + plot_layout(heights = c(1, 1.45)), "Figure1_Cohort_and_Scores", 198, 173)
}

figure_association_and_ipd <- function(ids, km, full) {
  single <- read_table("04_cohort_estimates")[,
    .(label = cohort, programme = exposure, HR, lo, hi, series = fifelse(role == "discovery", "Discovery", "External cohort"))]
  if (full) {
    rna <- read_table("05b_rnaseq_association")[,
      .(label = cohort, programme = exposure, HR, lo, hi, series = "External cohort")]
    single <- rbind(single, rna)
  }
  arrays <- read_table("04_meta_analysis")[method == "REML_mKH",
    .(label = sprintf("REML-mKH (k=%d)", k), programme = axis, HR = pooled_HR, lo, hi, series = "Two-stage")]
  synthesis <- read_table("04b_ipd_and_pan_external_synthesis")
  pan <- synthesis[method == "REML_mKH",
    .(label = sprintf("REML-mKH (k=%d)", k), programme = axis, HR = pooled_HR, lo, hi, series = "Two-stage")]
  ipd <- synthesis[method == "IPD_strata_cohort" & (k >= 4L | !full),
    .(label = sprintf("Stratified IPD (k=%d)", k), programme = axis, HR = pooled_HR, lo, hi, series = "One-stage")]
  forest <- rbindlist(list(single, arrays, pan, ipd))
  labels <- c(ids, unique(c(arrays$label, pan$label, ipd$label)))
  forest[, `:=`(label = factor(label, levels = rev(labels)), programme = main_axis(programme),
                series = factor(series, levels = c("Discovery", "External cohort", "Two-stage", "One-stage")))]
  note <- "IPD Cox models stratify baseline risk by cohort."
  if (full) note <- paste(note, "CPTAC cohort estimates use stage I vs II-IV.")
  panel_a <- main_forest(ggplot(forest, aes(HR, label, colour = series, shape = series)) +
                           geom_hline(yintercept = length(labels) - length(ids) + .5, colour = MAIN_GRID),
                         breaks = c(.5, 1, 2, 3)) +
    facet_wrap(~ programme, ncol = 4) +
    scale_colour_manual(values = c(Discovery = MAIN_SLATE, `External cohort` = MAIN_BLUE,
                                   `Two-stage` = MAIN_RUST, `One-stage` = MAIN_TEAL)) +
    scale_shape_manual(values = c(Discovery = 1, `External cohort` = 16, `Two-stage` = 15, `One-stage` = 18)) +
    labs(title = "Programme associations and external synthesis", subtitle = note)

  km <- copy(km)
  km[, axis := factor(axis, levels = unname(MAIN_AXES)[1:2])]
  panel_b <- ggplot(km, aes(time, surv, colour = group, linetype = group)) +
    geom_step(linewidth = .55) +
    geom_point(data = km[n.censor > 0], shape = 3, size = .6, alpha = .42) +
    facet_wrap(~ axis, ncol = 2) + coord_cartesian(xlim = c(0, 96), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(0, 24, 48, 72, 96)) +
    scale_y_continuous(labels = scales::label_percent(), breaks = c(0, .25, .5, .75, 1)) +
    scale_colour_manual(values = c(`At or below median` = MAIN_BLUE, `Above median` = MAIN_RUST)) +
    scale_linetype_manual(values = c(`At or below median` = "solid", `Above median` = "dashed")) +
    labs(title = "TCGA survival by median score", x = "Months from baseline", y = "Overall survival") +
    main_theme() + theme(panel.grid.major.y = element_line(colour = MAIN_GRID, linewidth = .25))

  stages <- c("I", "II", "III")
  stage <- read_table("03_stage_stratified")[stage %in% stages]
  stage[, `:=`(programme = main_axis(exposure), stage = factor(stage, levels = rev(stages)))]
  panel_c <- main_forest(ggplot(stage, aes(HR, stage)), breaks = c(.5, .75, 1, 1.5, 2, 3)) +
    facet_wrap(~ programme, ncol = 2) + labs(title = "TCGA estimates by stage", subtitle = "Exploratory models",
                                           x = "Adjusted HR per SD (95% CI)", y = "Stage")
  lower <- (panel_b | panel_c) + plot_layout(widths = c(1, 1.07))
  save_main_figure(panel_a / lower + plot_layout(heights = c(1.08, 1)),
                   "Figure2_Association_and_IPD_Synthesis", 205, 178)
}

main_purity_panel <- function() {
  purity <- read_table("03c_purity_sensitivity")[exposure %in% names(MAIN_AXES)]
  if (!nrow(purity)) {
    return(ggplot() + theme_void() + labs(title = "ABSOLUTE purity adjustment",
                                         subtitle = "ABSOLUTE estimates are available in the full analysis.") +
             theme(plot.title = element_text(face = "bold", size = 10.5, colour = MAIN_INK),
                   plot.subtitle = element_text(size = 8, colour = MAIN_SLATE),
                   plot.margin = margin(6, 7, 5, 6)))
  }
  labels <- unique(purity[, .(exposure, correlation_with_purity)])
  labels[, label := sprintf("%s\n\u03c1 = %+.2f", MAIN_AXES[exposure], correlation_with_purity)]
  labels <- setNames(labels$label, labels$exposure)
  purity[, `:=`(label = factor(labels[exposure], levels = rev(unname(labels[names(MAIN_AXES)]))),
                adjustment = factor(model, levels = c("clinical", "clinical_and_purity"),
                                    labels = c("Clinical", "Clinical + purity")))]
  hypoxia <- purity[exposure == "hypoxia", HR[order(adjustment)]]
  position <- position_dodge(width = .45, orientation = "y")
  main_forest(ggplot(purity, aes(HR, label, colour = adjustment, shape = adjustment)), position,
              breaks = c(.75, 1, 1.25, 1.5, 2)) +
    scale_colour_manual(values = c(Clinical = MAIN_SLATE, `Clinical + purity` = MAIN_TEAL)) +
    scale_shape_manual(values = c(Clinical = 1, `Clinical + purity` = 16)) +
    labs(title = "ABSOLUTE purity adjustment",
         subtitle = sprintf("Matched TCGA cases (n = %d)\nHypoxia HR %.2f to %.2f",
                            purity$n[1L], hypoxia[1L], hypoxia[2L]),
         caption = "\u03c1 is Spearman correlation with purity.", x = "Adjusted HR per SD (95% CI)")
}

main_drivers_panel <- function() {
  rows <- read_table("03d_driver_smoking_adjustment")[exposure %in% names(MAIN_AXES)[1:2]]
  descriptions <- c(pack_years_adjustment = "Smoking pack-years",
                    smoking_driver_adjustment = "Smoking and driver groups",
                    smoking_adjustment = "Smoking history",
                    driver_adjustment = "EGFR, KRAS, TP53, STK11")
  rows[, description := fifelse(analysis_type == "driver_subgroup", subgroup_or_model,
                                unname(descriptions[analysis_type]))]
  rows[cohort == "GSE31210" & description == "Wild-type / triple-negative", description := "Triple-negative"]
  sparse <- rows$cohort == "GSE72094" & rows$description == "EGFR mut"
  rows[sparse, description := "EGFR mut*"]
  rows[, label := paste(cohort, description, sep = " | ")]
  labels <- rev(unique(rows$label))
  models <- c(base_clinical = "Clinical", extended_adjusted = "Extended", driver_subgroup = "Driver subgroup")
  rows[, series := subgroup_or_model]
  rows[analysis_type == "driver_subgroup", series := "driver_subgroup"]
  rows[, `:=`(programme = main_axis(exposure), label = factor(label, levels = labels),
               series = factor(models[series], levels = unname(models)))]
  position <- position_dodge(width = .55, orientation = "y")
  divider <- uniqueN(rows[analysis_type == "driver_subgroup", label]) + .5
  note <- if (any(sparse)) "* EGFR mutant subgroup used ridge-penalised stage and sex terms (4 deaths)."
  main_forest(ggplot(rows, aes(HR, label, colour = series, shape = series)) +
                geom_hline(yintercept = divider, colour = MAIN_GRID), position) +
    facet_wrap(~ programme, ncol = 2, scales = "free_x") +
    scale_colour_manual(values = c(Clinical = MAIN_SLATE, Extended = MAIN_BLUE, `Driver subgroup` = MAIN_RUST)) +
    scale_shape_manual(values = c(Clinical = 1, Extended = 16, `Driver subgroup` = 18)) +
    labs(title = "Smoking and driver mutation sensitivity",
         subtitle = "Paired adjustments use the same complete cases. Subgroup models include clinical variables.",
         caption = note)
}

main_sensitivity_panel <- function(ids) {
  variants <- c(`Primary clinical adjustment` = "Clinical adjustment", `Stage only` = "Stage only",
                `Exclude stage IV` = "Exclude stage IV", `Mean gene-wise z-score` = "Mean gene-wise z-score",
                `Mean gene-wise z-score (sample-centred)` = "Mean gene-wise z-score\n(sample-centred)",
                `Winsorise at 1% and 99%` = "Winsorise at 1% and 99%",
                `Recurrence endpoint (exploratory)` = "Recurrence (exploratory)")
  sens <- read_table("06_sensitivity")
  sens[, `:=`(cohort = factor(cohort, levels = ids), programme = main_axis(axis),
               variant_label = factor(variants[variant], levels = rev(unname(variants))))]
  bound <- max(abs(log2(sens$HR)))
  ggplot(sens, aes(cohort, variant_label, fill = log2(HR))) +
    geom_tile(colour = "white", linewidth = .5) +
    geom_tile(data = sens[startsWith(variant, "Mean gene-wise z-score")],
              colour = "#61748A", fill = NA, linewidth = .30) +
    geom_text(aes(label = sprintf("%.2f", HR)), size = 2.6, colour = MAIN_INK) +
    facet_wrap(~ programme, ncol = 2) +
    scale_x_discrete(labels = function(x) sub("^TCGA-LUAD$", "TCGA\nLUAD", x)) +
    scale_fill_gradient2(low = "#9FBED5", mid = "#FAFAF8", high = "#D1937D", midpoint = 0,
                         limits = c(-bound, bound), breaks = log2(c(.5, 1, 1.5, 2)),
                         labels = c("0.5", "1", "1.5", "2")) +
    labs(title = "Scoring and endpoint sensitivity", x = NULL, y = NULL, fill = "Hazard ratio") +
    main_theme() + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6.9),
                         legend.title = element_text(size = 8), legend.key.width = grid::unit(8, "mm"),
                         legend.key.height = grid::unit(2.2, "mm"))
}

figure_specificity_and_purity <- function(ids, draws, null) {
  draws <- copy(draws)
  null <- copy(null)
  draws[, `:=`(programme = main_axis(axis), cohort = factor(cohort, levels = ids))]
  null[, `:=`(programme = main_axis(axis), cohort = factor(cohort, levels = ids))]
  panel_a <- ggplot(draws, aes(HR)) +
    geom_histogram(bins = 24, fill = "#8EADC4", colour = "white", linewidth = .2) +
    geom_vline(data = null, aes(xintercept = observed_HR), colour = MAIN_RUST, linewidth = .62) +
    facet_grid(cohort ~ programme, scales = "free") +
    scale_x_continuous(breaks = scales::breaks_extended(n = 4)) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 3)) +
    labs(title = "Matched random gene sets", subtitle = "Rust lines show observed hazard ratios.",
         x = "Adjusted hazard ratio", y = "Random gene sets") +
    main_theme() + theme(strip.text.y = element_text(angle = 0, size = 7.4),
                         strip.text.x = element_text(size = 7.7))
  panel_b <- main_purity_panel()
  panel_c <- main_drivers_panel()
  panel_d <- main_sensitivity_panel(setdiff(ids, "CPTAC-LUAD"))
  # Give each plot its own axis margins and keep tags beside the wrapper titles.
  wrap_panel <- function(panel) {
    title <- panel$labels$title
    subtitle <- panel$labels$subtitle
    body <- panel + labs(title = NULL, subtitle = NULL) + theme(plot.margin = margin(3, 2, 2, 2))
    wrap_elements(panel = body) + labs(title = title, subtitle = subtitle) + main_theme()
  }
  top <- (wrap_panel(panel_a) | wrap_panel(panel_b)) + plot_layout(widths = c(1.10, 1))
  layout <- top / wrap_panel(panel_c) / wrap_panel(panel_d)
  save_main_figure(layout + plot_layout(heights = c(1.00, 1.30, .95)),
                   "Figure3_Specificity_Purity_and_Drivers", 210, 277)
}

main_benchmark_panel <- function() {
  simple <- read_table("03c_marker_synthesis")
  published <- read_table("03e_published_signature_benchmark")
  common_scope <- if ("extended" %in% intersect(simple$scope, published$scope)) "extended" else "validation"
  simple_labels <- c(proliferation = "Composite", MKI67 = "MKI67", e2f = "E2F", g2m = "G2/M")
  clinical_label <- "Clinical + score"
  simple <- simple[scope == common_scope & method == "REML_mKH" & model != "four_programmes_joint"]
  simple[, `:=`(marker = unname(simple_labels[fifelse(model == "clinical_plus_marker", exposure,
                                                      sub("^composite_given_", "", model))]),
                contrast = fifelse(model == "clinical_plus_marker", clinical_label, "Composite given marker"))]
  published <- published[scope == common_scope & method == "REML_mKH" & scoring == "ssGSEA"]
  published[, `:=`(marker = unname(MAIN_MARKERS[signature_id]), contrast = fcase(
    model == "clinical_plus_marker", clinical_label, startsWith(model, "composite_given_"), "Composite given marker",
    default = "Marker given composite"))]
  rows <- rbind(simple[, .(marker, contrast, HR = pooled_HR, lo, hi)], published[, .(marker, contrast, HR, lo, hi)])
  rows[, `:=`(marker = factor(marker, levels = rev(unname(MAIN_MARKERS))),
               contrast = factor(contrast, levels = c(clinical_label, "Composite given marker", "Marker given composite")))]
  position <- position_dodge(width = .62, orientation = "y")
  main_forest(ggplot(rows, aes(HR, marker, colour = contrast, shape = contrast)), position,
              breaks = c(.25, .5, .75, 1, 1.5, 2, 3, 4, 6)) +
    scale_colour_manual(values = c(`Clinical + score` = MAIN_BLUE, `Composite given marker` = MAIN_RUST,
                                   `Marker given composite` = MAIN_TEAL)) +
    scale_shape_manual(values = c(`Clinical + score` = 16, `Composite given marker` = 17,
                                  `Marker given composite` = 15)) +
    labs(title = "External associations and conditional effects",
         subtitle = sprintf("REML-mKH on the same %d GEO cohorts. Reverse conditional estimates are available for published signatures.",
                            simple$k[1L]))
}

figure_redundancy <- function(ids, transport) {
  published <- read_table("05d_published_signature_transport")[scoring == "ssGSEA"]
  published[, model := signature_id]
  fields <- c("cohort", "model", "delta_c", "delta_lo", "delta_hi", "delta_uno_c", "delta_uno_lo", "delta_uno_hi")
  scores <- rbind(transport[model != "clinical", ..fields], published[, ..fields])
  comparison <- rbind(scores[, .(cohort, model, index = "Harrell C", delta = delta_c, lo = delta_lo, hi = delta_hi)],
                      scores[, .(cohort, model, index = "Uno IPCW C", delta = delta_uno_c,
                                  lo = delta_uno_lo, hi = delta_uno_hi)])
  comparison[, `:=`(cohort = factor(cohort, levels = ids),
                     model = factor(MAIN_MARKERS[model], levels = rev(unname(MAIN_MARKERS))))]
  position <- position_dodge(width = .50, orientation = "y")
  panel_b <- ggplot(comparison, aes(delta, model, colour = index, shape = index)) +
    geom_vline(xintercept = 0, colour = "#8796A3", linetype = 2, linewidth = .36) +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y", position = position, linewidth = .42) +
    geom_point(position = position, size = 1.65) +
    facet_wrap(~ cohort, ncol = 3) +
    scale_colour_manual(values = c(`Harrell C` = MAIN_BLUE, `Uno IPCW C` = MAIN_RUST)) +
    scale_shape_manual(values = c(`Harrell C` = 16, `Uno IPCW C` = 17)) +
    scale_x_continuous(breaks = seq(-.15, .30, by = .05), labels = scales::label_number(accuracy = .01),
                       expand = expansion(mult = .08)) +
    labs(title = "Discrimination in external cohorts",
         subtitle = "Published signatures use ssGSEA. Intervals use paired patient bootstrap.",
         x = "Change in C versus the frozen clinical model (95% CI)", y = NULL) + main_theme()
  save_main_figure(main_benchmark_panel() / panel_b + plot_layout(heights = c(.88, 1.23)),
                   "Figure4_Redundancy_Benchmark", 205, 184)
}

figure_recalibration <- function(ids, transport, full) {
  external <- setdiff(ids, "TCGA-LUAD")
  model_labels <- c(clinical = "Clinical", composite = "Composite")
  stages <- c(frozen = "Frozen", baseline_recalibrated = "Baseline recalibrated")
  phase_colours <- c(Frozen = MAIN_SLATE, `Baseline recalibrated` = MAIN_TEAL, `Observed KM` = MAIN_INK)
  phase_shapes <- c(Frozen = 1, `Baseline recalibrated` = 16, `Observed KM` = 4)
  quartiles <- read_table("05c_recalibration")
  if (full) quartiles <- rbind(quartiles, read_table("05b_rnaseq_recalibration"))
  quartiles[, `:=`(cohort = factor(cohort, levels = external),
                    model = factor(model_labels[model], levels = unname(model_labels)),
                    phase = factor(stages[calibration_type], levels = unname(stages)))]
  max_risk <- min(1, ceiling(10 * max(quartiles$predicted, quartiles$observed_hi, na.rm = TRUE)) / 10)
  panel_a <- ggplot(quartiles, aes(predicted, observed, colour = phase, shape = phase)) +
    geom_abline(slope = 1, intercept = 0, colour = "#94A2AF", linewidth = .38, linetype = 2) +
    geom_linerange(aes(ymin = observed_lo, ymax = observed_hi), linewidth = .33, alpha = .51, na.rm = TRUE) +
    geom_line(aes(group = phase, linetype = phase), linewidth = .49, na.rm = TRUE) +
    geom_point(size = 1.55, na.rm = TRUE) +
    facet_grid(model ~ cohort) + coord_equal(xlim = c(0, max_risk), ylim = c(0, max_risk)) +
    scale_x_continuous(labels = scales::label_percent(), breaks = c(0, .25, .5, .75, 1)) +
    scale_y_continuous(labels = scales::label_percent(), breaks = c(0, .25, .5, .75, 1)) +
    scale_colour_manual(values = phase_colours) + scale_shape_manual(values = phase_shapes) +
    scale_linetype_manual(values = c(Frozen = "dashed", `Baseline recalibrated` = "solid")) +
    labs(title = "Quartile calibration at 36 months", subtitle = "Points include Kaplan-Meier 95% confidence intervals.",
         x = "Mean predicted mortality", y = "Observed mortality") + main_theme() +
    theme(panel.spacing = grid::unit(2, "mm"), strip.text.y = element_text(angle = 0, size = 7.7),
          strip.text.x = element_text(size = 7.6), axis.text = element_text(size = 7.2),
          panel.grid.major.y = element_line(colour = MAIN_GRID, linewidth = .25))

  eligible <- transport[model %in% names(model_labels) & is.finite(recal_predicted_36)]
  measures <- c("Mean mortality risk", "IPCW Brier score", "Reduction in Brier score")
  pairs <- rbindlist(list(
    eligible[, .(cohort, model, measure = measures[1L], Frozen = predicted_36,
                 `Baseline recalibrated` = recal_predicted_36, `Observed KM` = observed_36)],
    eligible[, .(cohort, model, measure = measures[2L], Frozen = brier_36,
                 `Baseline recalibrated` = recal_brier_36, `Observed KM` = NA_real_)],
    eligible[, .(cohort, model, measure = measures[3L], Frozen = 0,
                 `Baseline recalibrated` = brier_36 - recal_brier_36, `Observed KM` = NA_real_)]
  ))
  pairs[, `:=`(cohort = factor(cohort, levels = rev(external)), model = factor(model_labels[model], levels = unname(model_labels)),
                measure = factor(measure, levels = measures))]
  points <- melt(pairs, id.vars = c("cohort", "model", "measure"), variable.name = "phase", value.name = "value",
                 variable.factor = TRUE, na.rm = TRUE)
  reductions <- pairs[measure == measures[3L]]
  reductions[, label_x := 1.10 * max(abs(`Baseline recalibrated`))]
  panel_b <- ggplot(mapping = aes(y = cohort)) +
    geom_segment(data = pairs, aes(x = Frozen, xend = `Baseline recalibrated`, yend = cohort),
                 colour = "#9CAEBB", linewidth = .68) +
    geom_point(data = points, aes(x = value, colour = phase, shape = phase), size = 2) +
    geom_text(data = reductions, aes(x = label_x, label = sprintf("%+.4f", `Baseline recalibrated`)),
              hjust = 0, colour = MAIN_INK, size = 2.45) +
    facet_grid(model ~ measure, scales = "free_x") +
    scale_x_continuous(labels = scales::label_number(accuracy = .01), breaks = scales::breaks_extended(n = 5),
                       expand = expansion(mult = c(.08, .24))) +
    scale_colour_manual(values = phase_colours) + scale_shape_manual(values = phase_shapes) +
    labs(title = "Mean mortality and prediction error at 36 months",
         subtitle = "Positive Brier reductions indicate lower prediction error after recalibration.",
         x = "Risk or score (proportion)", y = NULL) + main_theme() +
    theme(strip.text.y = element_text(angle = 0, size = 8))
  save_main_figure(panel_a / panel_b + plot_layout(heights = c(1.13, 1)),
                   "Figure5_Calibration_and_Recalibration", 205, 195)
}

write_main_figures <- function(scores, km, draws, null, full) {
  audit <- read_table("01_cohort_audit")[, .(cohort, participants = n_os, complete = n_complete,
                                              followup = median_followup_mo)]
  transport <- read_table("05_incremental_value")
  if (full) {
    rna <- read_table("05b_rnaseq_audit")[, .(cohort, participants = with_os, complete = complete_cases,
                                               followup = median_followup_months)]
    audit <- rbind(audit, rna)
    scores[["CPTAC-LUAD"]] <- cptac_programme_scores()
    transport <- rbind(transport, read_table("05b_rnaseq_transport"))
  }
  ids <- audit$cohort
  figure_cohorts_and_scores(audit, scores)
  figure_association_and_ipd(ids, km, full)
  figure_specificity_and_purity(ids, draws, null)
  figure_redundancy(ids, transport)
  figure_recalibration(ids, transport, full)
}
