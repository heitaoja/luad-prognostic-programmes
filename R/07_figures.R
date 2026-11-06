# =============================================================================
# 07_figures.R  --  main figures
# -----------------------------------------------------------------------------
# All figures read only results/tables/*.csv (plus the score table), so they can
# be rebuilt without rerunning the analysis.
# Output: results/figures/fig1..fig8 (.pdf + .png, 300 dpi)
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(survival); library(patchwork)
})
PROJ <- normalizePath(".", winslash = "/")
TAB  <- file.path(PROJ, "results", "tables")
FIG  <- file.path(PROJ, "results", "figures")
dir.create(FIG, TRUE, showWarnings = FALSE)
source("R/lib_surv.R")

AXES <- c("proliferation", "hypoxia", "emt_stroma", "immune")
AXIS_LAB <- c(proliferation = "Proliferation", hypoxia = "Hypoxia / glycolysis",
              emt_stroma = "EMT / stroma", immune = "Immune infiltration")
COL <- c(proliferation = "#C0392B", hypoxia = "#8E44AD", emt_stroma = "#16A085",
         immune = "#2980B9")
COL_COH <- c(`TCGA-LUAD` = "#111111", GSE31210 = "#E67E22", GSE37745 = "#2980B9",
             GSE50081 = "#27AE60")

theme_pub <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = NA),
          plot.title = element_text(face = "bold", size = base + 1),
          legend.key = element_blank())
}
theme_set(theme_pub())

save_fig <- function(p, name, w = 180, h = 120) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h,
         units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h,
         units = "mm", dpi = 300)
  cat("  saved", name, "\n")
}

aud   <- fread(file.path(TAB, "01_cohort_audit.csv"))
est   <- fread(file.path(TAB, "04_cohort_estimates.csv"))
meta  <- fread(file.path(TAB, "04_meta_analysis.csv"))
inc   <- fread(file.path(TAB, "05_incremental_value.csv"))
dca   <- fread(file.path(TAB, "05_dca_curves.csv"))
nullt <- fread(file.path(TAB, "06_random_set_null.csv"))
perm  <- fread(file.path(TAB, "06_permutation_test.csv"))
sens  <- fread(file.path(TAB, "06_sensitivity.csv"))
strat <- fread(file.path(TAB, "03_stage_stratified.csv"))
lad   <- fread(file.path(TAB, "03_cindex_ladder.csv"))
scores <- readRDS(file.path(PROJ, "data", "processed", "scores.rds"))

# ---------------------------------------------------------------- Fig 1 ------
# cohort overview: n, events, median follow-up
aud[, cohort := factor(cohort, levels = rev(c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081")))]
a1 <- melt(as.data.table(aud)[, .(cohort, Samples = n_samples, Deaths = n_events)],
           id.vars = "cohort")
p1a <- ggplot(a1, aes(value, cohort, fill = variable)) +
  geom_col(position = position_dodge(width = .7), width = .62) +
  geom_text(aes(label = value), position = position_dodge(width = .7),
            hjust = -0.15, size = 2.7) +
  scale_fill_manual(values = c(Samples = "grey55", Deaths = "#C0392B"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .18))) +
  labs(x = "Patients", y = NULL, title = "Cohorts") + theme(legend.position = "top")

p1b <- ggplot(aud, aes(x = cohort, y = median_os_mo)) +
  geom_col(fill = "#34495E", width = .6) +
  geom_text(aes(label = sprintf("%.0f mo", median_os_mo)), vjust = -0.4, size = 2.7) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = "Median OS (months)", title = "Follow-up") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

txt <- data.table(
  x = 1, y = c(4, 3, 2, 1),
  lab = c("Discovery: TCGA-LUAD
(n = 515, 185 deaths)",
          "4 pre-specified Hallmark programmes
-> ssGSEA score per sample",
          "Validation: 3 independent GEO cohorts
(n = 437, 158 deaths)",
          "Incremental value + negative controls
+ sensitivity (21 variants)"))
p1c <- ggplot(txt) +
  geom_text(aes(x, y, label = lab), hjust = 0, size = 2.9, lineheight = 1.05) +
  geom_point(aes(x = 0.85, y), colour = "#C0392B", size = 2.6) +
  xlim(0.7, 8.4) + ylim(0.4, 4.6) +
  labs(title = "Analysis plan") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 10))

p1 <- (p1a | p1b | p1c) + plot_layout(widths = c(1.05, .95, 1.6)) +
  plot_annotation(tag_levels = "A")
save_fig(p1, "fig1_cohort_overview", 250, 85)

# ---------------------------------------------------------------- Fig 2 ------
# axis score distributions and their correlations (TCGA)
sc <- rbindlist(lapply(names(scores), function(id) {
  d <- as.data.table(scores[[id]])[, c("cohort", AXES), with = FALSE]
  melt(d, id.vars = "cohort")
}))
sc[, cohort := factor(cohort, levels = c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081"))]
p2a <- ggplot(sc, aes(value, fill = cohort)) +
  geom_density(alpha = .45, colour = NA) +
  facet_wrap(~variable, scales = "free", nrow = 1,
             labeller = as_labeller(AXIS_LAB)) +
  scale_fill_manual(values = COL_COH, name = NULL) +
  labs(x = "Score (z within cohort)", y = "Density") +
  theme(legend.position = "top", axis.text.x = element_text(size = 7))

cm <- cor(as.matrix(scores[["TCGA-LUAD"]][, AXES]))
cmd <- melt(as.data.table(cm, keep.rownames = "from"), id.vars = "from")
cmd[, variable := factor(variable, levels = AXES)]
p2b <- ggplot(cmd, aes(from, variable, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.8) +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-.7, .7), name = "r") +
  scale_x_discrete(labels = AXIS_LAB) + scale_y_discrete(labels = AXIS_LAB) +
  labs(x = NULL, y = NULL, title = "Correlation between axes (TCGA)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- (p2a / p2b) + plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "A")
save_fig(p2, "fig2_axis_distribution", 200, 150)

# ---------------------------------------------------------------- Fig 3 ------
# forest plot: adjusted HR per SD, by cohort, with the pooled estimate as a row
fp <- est[, .(cohort, axis = exposure, HR, lo, hi, role)]
pooled <- meta[scope == "all", .(cohort = "Pooled (random effects)", axis = axis,
                                 HR = pooled_HR, lo = lo, hi = hi, role = "pooled")]
fp <- rbind(fp, pooled, fill = TRUE)
fp[, ylev := factor(cohort, levels = rev(c("Pooled (random effects)", "GSE50081",
                                           "GSE37745", "GSE31210", "TCGA-LUAD")))]
fp[, axis := factor(axis, levels = AXES)]
fp[, is_pool := cohort == "Pooled (random effects)"]

p3 <- ggplot(fp, aes(x = HR, y = ylev, colour = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .55) +
  geom_point(aes(shape = role), size = 2.3) +
  facet_wrap(~axis, nrow = 1, labeller = as_labeller(AXIS_LAB)) +
  scale_x_log10(breaks = c(.8, 1, 1.25, 1.5, 2)) +
  scale_colour_manual(values = c(COL_COH, `Pooled (random effects)` = "#C0392B"),
                      guide = "none") +
  scale_shape_manual(values = c(discovery = 16, validation = 17, pooled = 23),
                     name = NULL) +
  labs(x = "Adjusted hazard ratio per 1 SD (95% CI)", y = NULL,
       title = "Prognostic value of each programme across cohorts") +
  theme(legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey93"),
        axis.text.x = element_text(size = 7))
save_fig(p3, "fig3_forest_all_cohorts", 230, 110)

# ---------------------------------------------------------------- Fig 4 ------
# KM in TCGA for the two axes that survive
km_panel <- function(axis) {
  d <- scores[["TCGA-LUAD"]]
  d$grp <- factor(ifelse(d[[axis]] > median(d[[axis]], na.rm = TRUE), "High", "Low"),
                  levels = c("Low", "High"))
  km <- km_df(d, "grp")
  km$df$group <- factor(sub("^[^=]+=", "", km$df$group), levels = c("Low", "High"))
  p <- logrank_p(d, "grp")
  plab <- if (p < 0.001) "log-rank p < 0.001" else sprintf("log-rank p = %.3f", p)
  ggplot(km$df, aes(time, surv, colour = group, fill = group)) +
    geom_step(linewidth = .7) +
    annotate("text", x = 0, y = 0.08, hjust = 0, size = 2.9, label = plab) +
    scale_colour_manual(values = c(Low = "#2980B9", High = "#C0392B"), name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Months", y = "Overall survival",
         title = AXIS_LAB[[axis]]) +
    theme(legend.position = c(.78, .86), legend.background = element_blank(),
          legend.key = element_blank(), legend.title = element_blank())
}
p4 <- (km_panel("proliferation") | km_panel("hypoxia")) +
  plot_annotation(tag_levels = "A", title = "TCGA-LUAD, median split (visualisation only)")
save_fig(p4, "fig4_km_tcga", 200, 95)

# ---------------------------------------------------------------- Fig 5 ------
# stage-stratified hazard ratios
st <- strat[exposure %in% c("proliferation", "hypoxia")]
st[, axis := factor(exposure, levels = c("proliferation", "hypoxia"))]
st[, stage := factor(stage, levels = c("I", "II", "III"))]
p5 <- ggplot(st, aes(HR, stage, colour = axis)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .5) +
  geom_point(size = 2.2) +
  scale_x_log10(breaks = c(.7, 1, 1.5, 2, 2.5)) +
  scale_y_discrete(limits = rev(c("I", "II", "III"))) +
  scale_colour_manual(values = COL[c("proliferation", "hypoxia")],
                      labels = AXIS_LAB[c("proliferation", "hypoxia")], name = NULL) +
  facet_wrap(~axis, labeller = as_labeller(AXIS_LAB)) +
  labs(x = "Adjusted HR per 1 SD (95% CI)", y = "AJCC stage",
       title = "Stage-stratified effects in TCGA-LUAD") +
  theme(legend.position = "none")
save_fig(p5, "fig5_stage_stratified", 170, 95)

# ---------------------------------------------------------------- Fig 6 ------
# incremental value: delta C-index and decision curves
i6 <- inc[axis %in% c("proliferation", "hypoxia")]
i6[, axis := factor(axis, levels = c("proliferation", "hypoxia"))]
i6[, cohort := factor(cohort, levels = c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081"))]
p6a <- ggplot(i6, aes(delta_c, cohort, colour = axis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45") +
  geom_linerange(aes(xmin = delta_lo, xmax = delta_hi), linewidth = .5) +
  geom_point(size = 2.2, position = position_dodge(width = .5)) +
  scale_colour_manual(values = COL[c("proliferation", "hypoxia")],
                      labels = AXIS_LAB[c("proliferation", "hypoxia")], name = NULL) +
  labs(x = expression(Delta * "C-index (95% CI)"), y = NULL,
       title = "Added discrimination beyond stage + age + sex") +
  theme(legend.position = "top")

d6 <- dca[model %in% c("stage+clin", "stage+clin+hypoxia", "stage+clin+proliferation")]
d6[, model := factor(model, levels = c("stage+clin", "stage+clin+proliferation", "stage+clin+hypoxia"),
                     labels = c("Stage + clinical", "+ proliferation", "+ hypoxia"))]
p6b <- ggplot(d6, aes(threshold, nb_model, colour = model)) +
  geom_line(linewidth = .7) +
  geom_line(aes(y = nb_all), linetype = 3, colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey70") +
  coord_cartesian(xlim = c(.05, .55), ylim = c(-0.02, 0.22)) +
  scale_colour_manual(values = c("grey20", COL[["proliferation"]], COL[["hypoxia"]]),
                      name = NULL) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "Decision curves, 36-month mortality (TCGA)") +
  theme(legend.position = "top")
p6 <- (p6a | p6b) + plot_layout(widths = c(1, 1.2)) + plot_annotation(tag_levels = "A")
save_fig(p6, "fig6_incremental", 240, 105)

# ---------------------------------------------------------------- Fig 7 ------
# negative controls
nl <- rbindlist(lapply(c("proliferation", "hypoxia"), function(a) {
  f <- file.path(TAB, paste0("06_null_hr_TCGA_LUAD_", a, ".csv"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f); x[, axis := a]; x
}))
obs <- nullt[cohort == "TCGA-LUAD"]
p7a <- ggplot(nl, aes(HR, fill = axis)) +
  geom_histogram(bins = 40, alpha = .55, colour = NA, position = "identity") +
  geom_vline(data = obs, aes(xintercept = observed_HR, colour = axis),
             linetype = 2, linewidth = .8) +
  facet_wrap(~axis, labeller = as_labeller(AXIS_LAB), nrow = 1) +
  scale_fill_manual(values = COL[c("proliferation", "hypoxia")], guide = "none") +
  scale_colour_manual(values = COL[c("proliferation", "hypoxia")], guide = "none") +
  labs(x = "Adjusted HR of a random gene set", y = "Count",
       title = sprintf("Null distribution (%s random sets of matched size); dashed = observed",
                       unique(nl[, .N, by = axis]$N[1])))

pdraws_f <- file.path(TAB, "06_perm_hr.csv")
p7b <- if (file.exists(pdraws_f)) {
  pd <- fread(pdraws_f)
  pd[, axis := factor(axis, levels = c("proliferation", "hypoxia"))]
  perm[, axis := factor(axis, levels = c("proliferation", "hypoxia"))]
  ggplot(pd, aes(HR)) +
    geom_histogram(bins = 40, fill = "grey60", colour = NA) +
    geom_vline(data = perm, aes(xintercept = observed_HR, colour = axis),
               linetype = 2, linewidth = .8, show.legend = FALSE) +
    geom_vline(xintercept = 1, linetype = 3, colour = "grey30") +
    facet_wrap(~axis, labeller = as_labeller(AXIS_LAB), nrow = 1) +
    scale_colour_manual(values = COL[c("proliferation", "hypoxia")]) +
    labs(x = "HR after shuffling the survival labels", y = "Count",
         title = "Label permutation (1000 shuffles); dashed = observed, dotted = 1")
} else {
  ggplot() + annotate("text", x = 1, y = 1, label = "permutation draws unavailable")
}
p7 <- (p7a / p7b) + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "A")
save_fig(p7, "fig7_negative_control", 200, 150)

# ---------------------------------------------------------------- Fig 8 ------
# sensitivity heatmap
s8 <- sens[axis %in% c("proliferation", "hypoxia") & !is.na(HR)]
s8[, variant := factor(variant, levels = unique(variant))]
s8[, cohort := factor(cohort, levels = rev(c("TCGA-LUAD", "GSE31210", "GSE37745", "GSE50081")))]
s8[, lab := sprintf("%.2f%s", HR, ifelse(p < .05, "*", ""))]
p8 <- ggplot(s8, aes(variant, cohort, fill = HR)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = lab), size = 2.5) +
  facet_wrap(~axis, labeller = as_labeller(AXIS_LAB), nrow = 1) +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                       midpoint = 1, limits = c(0.7, 2), name = "HR") +
  labs(x = NULL, y = NULL,
       title = "Sensitivity analysis: adjusted HR across 21 variants (* p < 0.05)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7))
save_fig(p8, "fig8_sensitivity", 220, 130)

cat("\nall figures written to", FIG, "\n")
