# scripts/figures/figure3_tet_mechanism.R
# =============================================================================
# FIGURE 3: TET-Impediment, Not DNMT3A Recruitment  (paper section R2)
#
# Adjudicates the mechanism of BAP1-KO gene-body methylation remodeling:
# the coordinated 5mC-up / 5hmC-down phenotype is driven by impaired
# TET-mediated active demethylation (substrate-limited), NOT by H2AK119ub
# directly recruiting DNMT3A. At neuronal loci the mC-for-hmC exchange is
# stoichiometric (slope ~= -1, dehydroxymethylase-like), distinct from the
# non-stoichiometric (TET-inhibition + de-novo) signature genome-wide.
#
# Panels (per approved master plan, Figure 3 table):
#   A  Sec23  Dose-response bar      baseline-5hmC decile vs median hmC delta
#             NEW computation: bin all genes into deciles by WT 5hmC
#             (mean_mod_group1 from pre-loaded hmc_dmr), median hmc_diff/bin.
#             Annotate Model A AUC = 0.762 (baseline_hmc_model_comparison.tsv).
#   B  Sec24  Model comparison       TET 0.793 vs DNMT3A 0.696 pointrange
#             (dnmt3a_model_comparison.tsv; DeLong p = 9.43e-49).
#   C  Sec22  Chromatin-state lollipop  median delta_ratio per state
#             (demethylation_ratio_by_chromatin_state.tsv; p_adj/sig_label).
#   D  Sec78  Stoichiometry scatter  Δ5mC vs Δ5hmC over the SAME diffbind
#             universe Section 78 fit (diffbind_gene_level_all_marks.tsv,
#             20,915 genes; axes x=Δ5hmC, y=Δ5mC to match section_78 78b)
#             + 2 regression lines (all genes -0.959, neuronal -0.995)
#             from 78_stoichiometry_slopes.tsv; annotate differs-from-(-1).
#   E  Sec26  Attenuation bar        BAP1-KO vs TET-KO median delta_ratio
#             (tet_ko_comparison_summary.tsv; 3.3% absolute attenuation).
#
# Layout: design = "AABB\nCCDD\n#EE#", heights = c(1, 1, 0.6)
# Dimensions: 180 x 200 mm
#
# Data sources (all confirmed present in plots/visualizations/tables/):
#   demethylation_ratio_all_genes.tsv          (20,915 genes; panel A)
#   baseline_hmc_model_comparison.tsv          (3 models;     panel A AUC)
#   dnmt3a_model_comparison.tsv                 (5 models;     panel B)
#   demethylation_ratio_by_chromatin_state.tsv (7 states;     panel C)
#   diffbind_gene_level_all_marks.tsv           (20,915 genes; panel D scatter)
#   78_stoichiometry_slopes.tsv                 (7 groups;     panel D slopes)
#   72_neuronal_gene_set_go_derived.tsv         (5,614 genes;  panel D neuronal)
#   tet_ko_comparison_summary.tsv               (37 metrics;   panel E)
# =============================================================================

# Universal contract: source shared config (from downstream/) then figure config.
source("scripts/viz_sections/_shared_config.R")
source("scripts/figures/_figure_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

cat("\n========================================\n")
cat("FIGURE 3: TET-Impediment, Not DNMT3A Recruitment\n")
cat("========================================\n\n")

# DeLong p-value (DNMT3A recruitment vs TET impediment, shared models).
# Exact value lives in dnmt3a_exclusive_model_comparison.tsv (delong_p column)
# and is documented in docs/results/section_22_26_tet_dnmt3a_mechanism.md.
DELONG_P_SHARED <- 9.43e-49

# =============================================================================
# PANEL A: baseline-5hmC dose-response (NEW computation)
# Show the substrate-availability thesis as raw data: genes carrying more WT
# 5hmC (more oxidized TET substrate) show the largest hmC loss in BAP1-KO.
# =============================================================================
cat("Panel A: baseline-5hmC dose-response (decile binning)...\n")

ratio_genes <- read_table_tsv("demethylation_ratio_all_genes.tsv")
stopifnot(all(c("gene", "hmc_diff") %in% colnames(ratio_genes)))

# WT 5hmC baseline = mean_mod_group1 (control) from the pre-loaded, gene-deduped
# hmc_dmr object (per plan: "add WT 5hmC from hmc_dmr$mean_mod_group1").
hmc_baseline <- hmc_dmr[, c("gene", "mean_mod_group1")]
colnames(hmc_baseline) <- c("gene", "wt_hmc")

panelA_data <- ratio_genes %>%
  dplyr::select(gene, hmc_diff) %>%
  dplyr::inner_join(hmc_baseline, by = "gene") %>%
  dplyr::filter(is.finite(hmc_diff), is.finite(wt_hmc))

if (nrow(panelA_data) == 0) {
  stop("Panel A: no genes after joining ratio table to hmc_dmr WT 5hmC baseline.")
}

# Decile bins by WT 5hmC (1 = lowest substrate, 10 = highest substrate).
panelA_data$decile <- dplyr::ntile(panelA_data$wt_hmc, 10)

panelA_summary <- panelA_data %>%
  dplyr::group_by(decile) %>%
  dplyr::summarise(
    n              = dplyr::n(),
    median_hmc_diff = median(hmc_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(decile)

cat(sprintf("  Panel A: %d genes binned into %d deciles\n",
            nrow(panelA_data), nrow(panelA_summary)))
cat(sprintf("  Decile 1 median hmc_diff = %+.4f; Decile 10 = %+.4f\n",
            panelA_summary$median_hmc_diff[1],
            panelA_summary$median_hmc_diff[nrow(panelA_summary)]))

# AUC of the baseline-5hmC predictor (Model A) for on-panel annotation.
baseline_models <- read_table_tsv("baseline_hmc_model_comparison.tsv")
modelA_auc <- baseline_models$auc[grepl("Baseline 5hmC", baseline_models$model)]
stopifnot(length(modelA_auc) == 1)
cat(sprintf("  Model A (baseline 5hmC) AUC = %.3f\n", modelA_auc))

panelA <- ggplot(panelA_summary,
                 aes(x = factor(decile), y = median_hmc_diff)) +
  geom_col(fill = PUB_COLORS$hmC, width = 0.8) +
  geom_line(aes(group = 1), colour = "grey30", linewidth = 0.4) +
  geom_point(colour = "grey20", size = 1) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  annotate("text", x = 1, y = max(panelA_summary$median_hmc_diff) * 0.5,
           label = stat_label("AUC", modelA_auc, "%.3f"),
           hjust = 0, vjust = 1, size = 2.6, fontface = "italic") +
  labs(
    title = "Baseline 5hmC predicts 5hmC loss",
    x     = "WT 5hmC decile (low → high substrate)",
    y     = "Median Δ5hmC (mut − ctrl)"
  ) +
  theme_pub()

# =============================================================================
# PANEL B: TET-impediment vs DNMT3A-recruitment model comparison
# =============================================================================
cat("Panel B: model comparison (TET vs DNMT3A)...\n")

dnmt3a_models <- read_table_tsv("dnmt3a_model_comparison.tsv")
stopifnot(all(c("model", "auc", "auc_ci_lower", "auc_ci_upper") %in%
              colnames(dnmt3a_models)))

# Order the four interpretable models worst -> best for a clean ladder.
model_order <- c("K119ub only", "DNMT3A recruitment", "TET impediment", "Full")
panelB_data <- dnmt3a_models %>%
  dplyr::filter(model %in% model_order) %>%
  dplyr::mutate(model = factor(model, levels = model_order))

tet_auc    <- panelB_data$auc[panelB_data$model == "TET impediment"]
dnmt3a_auc <- panelB_data$auc[panelB_data$model == "DNMT3A recruitment"]
stopifnot(length(tet_auc) == 1, length(dnmt3a_auc) == 1)
cat(sprintf("  TET impediment AUC = %.3f; DNMT3A recruitment AUC = %.3f\n",
            tet_auc, dnmt3a_auc))

# Highlight the two competing mechanistic models; Full/K119ub as context.
model_palette <- c(
  "TET impediment"     = PUB_COLORS$mC,
  "DNMT3A recruitment" = "#7570B3",
  "Full"               = "grey40",
  "K119ub only"        = "grey70"
)

panelB <- ggplot(panelB_data,
                 aes(x = auc, y = model, colour = model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             linewidth = 0.3, colour = "grey50") +
  geom_pointrange(aes(xmin = auc_ci_lower, xmax = auc_ci_upper),
                  linewidth = 0.6, size = 0.5) +
  scale_colour_manual(values = model_palette, guide = "none") +
  annotate("text",
           x = min(panelB_data$auc_ci_lower),
           y = 0.6,
           label = sprintf("DeLong %s", stat_label("p", DELONG_P_SHARED, "%.2e")),
           hjust = 0, vjust = 0, size = 2.6, fontface = "italic") +
  coord_cartesian(xlim = c(min(panelB_data$auc_ci_lower) - 0.02, 1.0),
                  clip = "off") +
  labs(
    title = "TET model out-discriminates DNMT3A",
    x     = "AUC (predicting mC hypermethylation)",
    y     = NULL
  ) +
  theme_pub()

# =============================================================================
# PANEL C: demethylation-ratio impairment by chromatin state (lollipop)
# =============================================================================
cat("Panel C: chromatin-state demethylation-ratio lollipop...\n")

state_ratio <- read_table_tsv("demethylation_ratio_by_chromatin_state.tsv")
stopifnot(all(c("chromatin_state", "median_delta", "sig_label") %in%
              colnames(state_ratio)))

# Order states by median delta_ratio so the most-impaired (most negative)
# state sits at the TOP of the lollipop. ggplot places the first factor level
# at the bottom of a discrete y-axis, so we sort descending (least negative
# first) and assign that as the level order -> most negative becomes the last
# level -> top of the panel.
panelC_data <- state_ratio %>%
  dplyr::arrange(dplyr::desc(median_delta)) %>%
  dplyr::mutate(chromatin_state = factor(chromatin_state,
                                         levels = chromatin_state))

ap_delta <- state_ratio$median_delta[state_ratio$chromatin_state == "Active_Promoter"]
cat(sprintf("  Active_Promoter median delta_ratio = %+.4f (strongest impairment)\n",
            ap_delta))

# Significance label position: just past the point, on the side of its sign.
panelC_data$label_x <- panelC_data$median_delta +
  ifelse(panelC_data$median_delta < 0, -0.0015, 0.0015)

panelC <- ggplot(panelC_data,
                 aes(x = median_delta, y = chromatin_state,
                     colour = chromatin_state)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
  geom_segment(aes(x = 0, xend = median_delta,
                   y = chromatin_state, yend = chromatin_state),
               linewidth = 0.6) +
  geom_point(size = 2) +
  geom_text(aes(x = label_x, label = sig_label),
            hjust = ifelse(panelC_data$median_delta < 0, 1, 0),
            size = 2.6, colour = "black", fontface = "bold") +
  scale_colour_manual(values = CHROMATIN_STATE_COLORS, guide = "none") +
  labs(
    title = "TET impairment is strongest at active chromatin",
    x     = "Median Δ demethylation ratio (mut − ctrl)",
    y     = NULL
  ) +
  theme_pub()

# =============================================================================
# PANEL D: stoichiometry scatter (Δ5hmC vs Δ5mC) + regression lines
#
# Reproduces the exact universe + axis orientation that Section 78 used to fit
# the authoritative slopes in 78_stoichiometry_slopes.tsv. Section 78 builds its
# working table from diffbind_gene_level_all_marks.tsv (20,915 genes) and fits
#   lm(mc_diff ~ hmc_diff)
# i.e. the slope is Δ5mC per unit Δ5hmC, with hmc_diff on x and mc_diff on y
# (section_78 78b: aes(x = hmc_diff, y = mc_diff), labs x="delta-5hmC",
# y="delta-5mC"). We therefore (a) scatter the SAME 20,915-gene diffbind
# universe — NOT coordinated_changes.tsv, which is a smaller 8,371-gene subset
# the slope was NOT fit on — and (b) keep x = Δ5hmC / y = Δ5mC so a line drawn
# with the TSV slope is geometrically correct. The neuronal subset resolves to
# 4,118 genes, matching the "Neuronal (broad)" slope's n exactly.
# =============================================================================
cat("Panel D: stoichiometry scatter with regression slopes...\n")

# Section-78 source table (20,915 genes); columns mc_diff (Δ5mC) + hmc_diff
# (Δ5hmC) are the same per-gene deltas used for the OLS fits.
stoich_genes <- read_table_tsv("diffbind_gene_level_all_marks.tsv")
stopifnot(all(c("gene", "mc_diff", "hmc_diff") %in% colnames(stoich_genes)))

# Broad GO-derived neuronal gene set — the membership used to fit the
# "Neuronal (broad)" slope in 78_stoichiometry_slopes.tsv (n = 4,118 within
# this diffbind universe).
neuronal_set <- read_table_tsv("72_neuronal_gene_set_go_derived.tsv")
stopifnot("gene" %in% colnames(neuronal_set))
neuronal_genes <- unique(neuronal_set$gene)

panelD_data <- stoich_genes %>%
  dplyr::select(gene, mc_diff, hmc_diff) %>%
  dplyr::filter(is.finite(mc_diff), is.finite(hmc_diff)) %>%
  dplyr::mutate(
    gene_class = ifelse(gene %in% neuronal_genes, "Neuronal", "Other")
  )
panelD_data$gene_class <- factor(panelD_data$gene_class,
                                 levels = c("Other", "Neuronal"))
cat(sprintf("  Panel D scatter: %d genes (%d neuronal)\n",
            nrow(panelD_data), sum(panelD_data$gene_class == "Neuronal")))

# Slopes from the stoichiometry TSV (the authoritative OLS fits).
slopes <- read_table_tsv("78_stoichiometry_slopes.tsv")
stopifnot(all(c("group", "slope", "differs_from_neg1") %in% colnames(slopes)))

slope_all  <- slopes$slope[slopes$group == "All genes"]
slope_neur <- slopes$slope[slopes$group == "Neuronal (broad)"]
diff_all   <- slopes$differs_from_neg1[slopes$group == "All genes"]
diff_neur  <- slopes$differs_from_neg1[slopes$group == "Neuronal (broad)"]
stopifnot(length(slope_all) == 1, length(slope_neur) == 1)
cat(sprintf("  Slope all genes = %.3f (differs from -1: %s)\n",
            slope_all, diff_all))
cat(sprintf("  Slope neuronal  = %.3f (differs from -1: %s)\n",
            slope_neur, diff_neur))

# Anchor each fitted line at its set's data centroid so the displayed slope
# (from the TSV; Δ5mC per unit Δ5hmC) passes through the cloud's centre of
# mass. x = Δ5hmC (hmc_diff), y = Δ5mC (mc_diff), matching Section 78.
mean_hmc_all <- mean(panelD_data$hmc_diff)
mean_mc_all  <- mean(panelD_data$mc_diff)
neur_sub <- panelD_data[panelD_data$gene_class == "Neuronal", ]
mean_hmc_neur <- mean(neur_sub$hmc_diff)
mean_mc_neur  <- mean(neur_sub$mc_diff)

int_all  <- mean_mc_all  - slope_all  * mean_hmc_all
int_neur <- mean_mc_neur - slope_neur * mean_hmc_neur

# Build segment endpoints spanning the x (Δ5hmC) range for each line.
x_rng <- range(panelD_data$hmc_diff)
line_df <- data.frame(
  group = factor(c("All genes", "Neuronal (broad)"),
                 levels = c("All genes", "Neuronal (broad)")),
  x     = c(x_rng[1], x_rng[1]),
  xend  = c(x_rng[2], x_rng[2]),
  y     = c(int_all + slope_all * x_rng[1],
            int_neur + slope_neur * x_rng[1]),
  yend  = c(int_all + slope_all * x_rng[2],
            int_neur + slope_neur * x_rng[2])
)

scatter_palette <- c("Other" = "grey75", "Neuronal" = PUB_COLORS$hmC)

slope_annot <- sprintf(
  "All genes: slope = %.3f (≠ −1)\nNeuronal: slope = %.3f (≈ −1)",
  slope_all, slope_neur
)

panelD <- ggplot(panelD_data, aes(x = hmc_diff, y = mc_diff)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_abline(slope = -1, intercept = 0, linetype = "dotted",
              linewidth = 0.4, colour = "#E69F00") +
  geom_point(aes(colour = gene_class), size = 0.3, alpha = 0.4) +
  geom_segment(
    data = line_df,
    aes(x = x, xend = xend, y = y, yend = yend, linetype = group),
    colour = "black", linewidth = 0.6, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = scatter_palette, name = "Gene class") +
  scale_linetype_manual(values = c("All genes" = "solid",
                                   "Neuronal (broad)" = "dashed"),
                        name = "OLS slope") +
  annotate("text", x = max(panelD_data$hmc_diff), y = max(panelD_data$mc_diff),
           label = slope_annot, hjust = 1, vjust = 1,
           size = 2.4, fontface = "italic") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(
    title = "Stoichiometric mC-for-hmC exchange at neuronal genes",
    x     = "Δ5hmC (mut − ctrl)",
    y     = "Δ5mC (mut − ctrl)"
  ) +
  theme_pub()

# =============================================================================
# PANEL E: BAP1-KO vs TET-KO demethylation-ratio attenuation bar
# "dimmer not switch": BAP1-KO reproduces only ~3.3% of the TET-KO shift.
# =============================================================================
cat("Panel E: TET-KO attenuation bar...\n")

tet_ko <- read_table_tsv("tet_ko_comparison_summary.tsv")
stopifnot(all(c("metric", "value") %in% colnames(tet_ko)))

get_metric <- function(name) {
  v <- tet_ko$value[tet_ko$metric == name]
  if (length(v) != 1) stop("Panel E: metric not found uniquely: ", name)
  as.numeric(v)
}

bap1_delta <- get_metric("BAP1-KO median delta")
tetko_delta <- get_metric("TET-KO median delta")
abs_atten  <- get_metric("Absolute attenuation (BAP1/TET)")
cat(sprintf("  BAP1-KO median delta = %.4f; TET-KO median delta = %.4f\n",
            bap1_delta, tetko_delta))
cat(sprintf("  Absolute attenuation (BAP1/TET) = %.4f (%.1f%%)\n",
            abs_atten, abs_atten * 100))

panelE_data <- data.frame(
  perturbation = factor(c("BAP1-KO", "TET-KO"),
                        levels = c("BAP1-KO", "TET-KO")),
  median_delta = c(bap1_delta, tetko_delta)
)

panelE_palette <- c("BAP1-KO" = PUB_COLORS$Mutant, "TET-KO" = "grey45")

panelE <- ggplot(panelE_data,
                 aes(x = perturbation, y = median_delta, fill = perturbation)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  geom_text(aes(label = sprintf("%.3f", median_delta)),
            vjust = 1.3, size = 2.6) +
  scale_fill_manual(values = panelE_palette, guide = "none") +
  annotate("text", x = 1.5, y = tetko_delta * 0.45,
           label = sprintf("BAP1-KO = %.1f%% of TET-KO\n(\"dimmer, not switch\")",
                           abs_atten * 100),
           size = 2.5, fontface = "italic") +
  labs(
    title = "BAP1-KO partially attenuates TET pathway",
    x     = NULL,
    y     = "Median Δ demethylation ratio"
  ) +
  theme_pub()

# =============================================================================
# COMPOSE: patchwork layout per master plan
#   design = "AABB\nCCDD\n#EE#", heights = c(1, 1, 0.6)
# =============================================================================
cat("\nComposing Figure 3...\n")

layout_design <- "AABB\nCCDD\n#EE#"

# Single plot_annotation() call carries BOTH the panel-letter tags and the
# figure title, avoiding any field-merge ambiguity from chaining two calls.
figure3 <- panelA + panelB + panelC + panelD + panelE +
  plot_layout(design = layout_design, heights = c(1, 1, 0.6)) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 3. TET-impediment, not DNMT3A recruitment, drives BAP1-KO methylation remodeling",
    theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))
  )

# =============================================================================
# SAVE: composite + load-bearing sub-panels (PDF + SVG + JPG)
# =============================================================================
cat("Saving Figure 3 (composite + sub-panels)...\n")

save_figure(figure3, "figure3_tet_mechanism", width_mm = 180, height_mm = 200)

# Sub-panels saved individually for flexible manuscript assembly.
save_figure(panelB, "figure3_panelB_model_comparison",
            width_mm = 90, height_mm = 70)
save_figure(panelD, "figure3_panelD_stoichiometry",
            width_mm = 100, height_mm = 90)

cat("\nFigure 3 complete.\n")
cat("Output: ", file.path(FIGURE_DIR, "figure3_tet_mechanism"), "\n")
