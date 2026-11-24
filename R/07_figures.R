# Publication figures from the revised analyses.
source("R/lib_surv.R")
suppressPackageStartupMessages(library(patchwork))
fig <- "results/figures"
tab <- "results/tables"
read_table <- function(name) fread(file.path(tab, paste0(name, ".csv")))
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
cohorts <- read_table("01_cohort_audit")
ids <- cohorts$cohort
axes <- c(proliferation = "Proliferation", hypoxia = "Hypoxia / glycolysis",
          emt_stroma = "EMT / stroma", immune = "Immune / inflammation")
marker <- c(proliferation = "Composite", MKI67 = "MKI67", e2f = "E2F", g2m = "G2/M")
models <- c(clinical = "Clinical", composite = "Composite", MKI67 = "MKI67",
            E2F = "E2F", G2M = "G2/M")
ink <- "#203045"
blue <- "#285D86"
red <- "#9F472C"
palette <- c("#285D86", "#9F472C", "#2F7359", "#7B5A95", "#526274")
base <- theme_minimal(base_size = 10) +
  theme(text = element_text(colour = ink), panel.grid.minor = element_blank(),
        panel.spacing.x = grid::unit(8, "mm"),
        plot.margin = margin(6, 12, 6, 6, unit = "pt"),
        strip.text = element_text(face = "bold"),
        plot.title.position = "plot", plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")
save_figure <- function(p, name, width = 180, height = 108) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = width, height = height,
         units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = width, height = height,
         units = "mm", dpi = 300, bg = "white")
}
cohort_order <- function(x) factor(x, levels = rev(ids))
axis_order <- function(x) factor(unname(axes[x]), levels = unname(axes))
hr_plot <- function(p) {
  p + geom_vline(xintercept = 1, linetype = 2, colour = "#768290") +
    geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = .44) +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
                  labels = function(x) formatC(x, format = "fg"),
                  expand = expansion(mult = c(.07, .08))) +
    labs(x = "Adjusted hazard ratio per SD (95% CI)", y = NULL) + base
}

cohorts$cohort <- cohort_order(cohorts$cohort)
p_count <- ggplot(cohorts, aes(n_samples, cohort)) +
  geom_col(fill = "#DCE7F0") +
  geom_col(aes(x = n_complete), fill = blue, width = .38) +
  geom_text(aes(label = n_samples), hjust = -.24, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, .13))) +
  labs(x = "Participants (dark bars are complete cases)", y = NULL) + base
p_fu <- ggplot(cohorts, aes(median_followup_mo, cohort)) +
  geom_point(colour = blue, size = 2.4) +
  geom_text(aes(label = sprintf("%.1f", median_followup_mo)), hjust = -.28, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(.05, .16))) +
  labs(x = "Reverse Kaplan-Meier follow-up (months)", y = NULL) + base
save_figure((p_count | p_fu) + plot_annotation(title = "Cohort size and follow-up"),
            "fig1_cohort_overview", 205, 82)

scores <- readRDS("data/processed/scores.rds")
d <- rbindlist(lapply(scores, function(d) {
  x <- as.data.table(d)[, c("cohort", names(axes)), with = FALSE]
  melt(x, id.vars = "cohort", variable.name = "axis")
}))
d[, axis := axis_order(axis)]
p_dist <- ggplot(d, aes(value, colour = cohort)) +
  geom_density(linewidth = .65) + facet_wrap(~ axis, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = setNames(palette[seq_along(ids)], ids)) +
  labs(title = "Distribution of cohort-standardised programme scores",
       x = "Score (SD)", y = "Density", colour = "Cohort") + base
save_figure(p_dist, "fig2_axis_distribution", 180, 113)

est <- read_table("04_cohort_estimates")
meta <- read_table("04_meta_analysis")[method == "REML_mKH"]
meta[, cohort := fifelse(scope == "validation", "Original external synthesis", "Expanded external synthesis")]
fp <- rbindlist(list(
  est[, .(cohort, axis = exposure, HR, lo, hi, source = fifelse(role == "discovery", "Discovery", "External"))],
  meta[, .(cohort, axis, HR = pooled_HR, lo, hi, source = "External synthesis")]
))
fp[, cohort := factor(cohort, levels = rev(c(ids, "Original external synthesis",
                                            "Expanded external synthesis")))]
fp[, axis := axis_order(axis)]
p_assoc <- hr_plot(ggplot(fp, aes(HR, cohort, colour = source, shape = source))) +
  facet_wrap(~ axis, ncol = 2) +
  scale_colour_manual(values = c(Discovery = "#586573", External = blue, "External synthesis" = red)) +
  labs(title = "Programme associations and external synthesis", colour = NULL, shape = NULL)
save_figure(p_assoc, "fig3_forest_all_cohorts", 208, 147)

km <- rbindlist(lapply(c("proliferation", "hypoxia"), function(a) {
  d <- scores[["TCGA-LUAD"]]
  d$group <- factor(ifelse(d[[a]] > median(d[[a]]), "Above median", "At or below median"))
  result <- km_df(d, "group")$df
  result$axis <- axes[a]
  result$group <- sub("^group=", "", result$group)
  rbind(data.frame(time = 0, surv = 1, n.censor = 0,
                   group = levels(d$group), axis = unname(axes[a])), result)
}))
p_km <- ggplot(km, aes(time, surv, colour = group)) +
  geom_step(linewidth = .65) +
  geom_point(data = km[n.censor > 0], shape = 3, size = .75, alpha = .5) +
  facet_wrap(~ axis) + coord_cartesian(xlim = c(0, 100), ylim = c(0, 1)) +
  scale_colour_manual(values = c(red, blue)) +
  labs(title = "TCGA-LUAD survival by median score",
       subtitle = "Median groups are used for visualisation",
       x = "Follow-up (months)", y = "Overall survival probability", colour = NULL) + base
save_figure(p_km, "fig4_km_tcga", 178, 86)

stage <- read_table("03_stage_stratified")
stage[, axis := axis_order(exposure)]
p_stage <- hr_plot(ggplot(stage, aes(HR, stage, colour = stage))) +
  facet_wrap(~ axis, ncol = 2) + theme(legend.position = "none") +
  labs(title = "TCGA-LUAD stage-specific estimates", subtitle = "Exploratory models")
save_figure(p_stage, "fig5_stage_stratified", 179, 98)

transport <- read_table("05_incremental_value")
transport[, cohort := factor(cohort, levels = ids)]
transport[, model_label := factor(models[model], levels = rev(unname(models)))]
d <- transport[model != "clinical"]
p_transport <- ggplot(d, aes(delta_c, model_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#768290") +
  geom_pointrange(aes(xmin = delta_lo, xmax = delta_hi, colour = model_label), linewidth = .45) +
  facet_wrap(~ cohort, ncol = 2) +
  scale_colour_manual(values = setNames(palette, unname(models))) +
  labs(title = "Discrimination after transport from TCGA-LUAD",
       x = "Change in Harrell C versus the fixed clinical model (95% CI)", y = NULL) +
  base + theme(legend.position = "none")
save_figure(p_transport, "fig6_incremental", 184, 119)

null <- read_table("06_random_set_null")
draws <- rbindlist(lapply(seq_len(nrow(null)), function(i)
  read_table(paste0("06_null_hr_", gsub("-", "_", null$cohort[i]), "_", null$axis[i]))))
p_null <- ggplot(draws, aes(HR)) +
  geom_histogram(bins = 30, fill = "#9FBDD4", colour = "white", linewidth = .2) +
  geom_vline(data = null, aes(xintercept = observed_HR), colour = red, linewidth = .7) +
  facet_grid(cohort ~ axis, scales = "free_y", labeller = labeller(axis = as_labeller(axes))) +
  labs(title = "Observed effects against matched random programmes",
       subtitle = "Null programmes preserve component gene counts, rank abundance and overlap",
       x = "Adjusted hazard ratio", y = "Random draws") + base
save_figure(p_null, "fig7_negative_control", 188, 132)

sens <- read_table("06_sensitivity")
sens[, cohort := cohort_order(cohort)]
sens[, axis := axis_order(axis)]
p_sens <- ggplot(sens, aes(variant, cohort, fill = log2(HR))) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", HR)), size = 2.9) +
  facet_wrap(~ axis, ncol = 1) +
  scale_fill_gradient2(low = "#C8D9EB", mid = "#F7F7F5", high = "#CB7757", midpoint = 0) +
  labs(title = "Sensitivity to endpoint, scoring and clinical adjustment",
       x = NULL, y = NULL, fill = "log2 HR") +
  base + theme(axis.text.x = element_text(angle = 18, hjust = 1))
save_figure(p_sens, "fig8_sensitivity", 195, 120)

bench <- read_table("03c_marker_synthesis")
selected_scope <- if ("extended" %in% bench$scope) "extended" else "validation"
b <- bench[scope == selected_scope & model != "four_programmes_joint"]
b[, label := fifelse(model == "clinical_plus_marker", marker[exposure],
                     paste0("Composite + adjustment for ",
                            marker[sub("^composite_given_", "", model)]))]
b[, contrast := fifelse(model == "clinical_plus_marker", "With clinical variables",
                         "Composite conditional on a simpler score")]
b[, label := factor(label, levels = rev(unique(label)))]
b[, HR := pooled_HR]
p_bench <- hr_plot(ggplot(b, aes(HR, label, colour = contrast))) +
  scale_colour_manual(values = c("With clinical variables" = blue,
                                 "Composite conditional on a simpler score" = red)) +
  labs(title = "External synthesis against simpler proliferation measures", colour = NULL)
save_figure(p_bench, "fig9_marker_benchmark", 191, 97)

cal <- read_table("05_calibration")
cal[, model := factor(models[model], levels = unname(models)[1:2])]
p_cal <- ggplot(cal, aes(predicted, observed, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#768290") +
  geom_linerange(aes(ymin = observed_lo, ymax = observed_hi), alpha = .57, linewidth = .48) +
  geom_point(size = 2.2) + facet_wrap(~ cohort, ncol = 2) +
  scale_colour_manual(values = unname(c("#586573", blue))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1)) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) +
  labs(title = "Calibration of fixed 36-month risks",
       subtitle = "Points are quartiles of predicted risk within test cohorts",
       x = "Mean predicted mortality", y = "Kaplan-Meier mortality", colour = NULL) + base
save_figure(p_cal, "fig10_calibration", 175, 154)

purity <- read_table("03c_purity_sensitivity")
if (nrow(purity)) {
  purity <- purity[exposure %in% names(axes)]
  purity[, axis := axis_order(exposure)]
  purity[, adjustment := factor(model, levels = c("clinical", "clinical_and_purity"),
                                labels = c("Clinical", "Clinical + ABSOLUTE purity"))]
  p_purity <- hr_plot(ggplot(purity, aes(HR, axis, colour = adjustment)) +
                      geom_blank()) +
    scale_colour_manual(values = c("#586573", blue)) +
    facet_wrap(~ adjustment) +
    labs(title = "DNA-derived tumour purity sensitivity",
         subtitle = paste0("Identical TCGA subsets in both models (n = ", unique(purity$n), ")"),
         colour = NULL) + theme(legend.position = "none")
  save_figure(p_purity, "fig11_purity_sensitivity", 188, 90)
}

if (Sys.getenv("LUAD_CORE_ONLY") != "1") {
  rna_audit <- read_table("05b_rnaseq_audit")
  association <- read_table("05b_rnaseq_association")
  association[, axis := factor(axes[exposure], levels = rev(unname(axes)))]
  p_rna_hr <- hr_plot(ggplot(association, aes(HR, axis))) +
    labs(title = "Programme associations")
  rna_transport <- read_table("05b_rnaseq_transport")[model != "clinical"]
  rna_transport[, model_label := factor(models[model], levels = rev(unname(models[-1L])))]
  p_rna_c <- ggplot(rna_transport, aes(delta_c, model_label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "#768290") +
    geom_pointrange(aes(xmin = delta_lo, xmax = delta_hi), colour = blue, linewidth = .45) +
    labs(title = "Increment over the clinical model", x = "Change in Harrell C (95% CI)", y = NULL) +
    base
  subtitle <- sprintf("%d participants and %d deaths", rna_audit$complete_cases, rna_audit$complete_case_deaths)
  save_figure((p_rna_hr | p_rna_c) + plot_annotation(
    title = "Independent RNA-seq evaluation in CPTAC-LUAD", subtitle = subtitle, tag_levels = "A"),
    "figS1_rnaseq_validation", 204, 89)
}

source("R/lib_main_figures.R")
write_main_figures(scores, km, draws, null, full = Sys.getenv("LUAD_CORE_ONLY") != "1")
