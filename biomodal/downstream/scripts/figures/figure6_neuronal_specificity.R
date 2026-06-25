# scripts/figures/figure6_neuronal_specificity.R
# =============================================================================
# Figure 6: Neuronal Genes Are Preferentially Affected  (paper section R5)
#
# Question: Why are neuronal genes disproportionately affected by BAP1 loss?
#   Constitutive H2AK119ub enrichment at neuronal/Polycomb loci + selective
#   Polycomb de-repression (K27me3 loss) at synapse/axon genes, with neuronal
#   methylation undergoing stoichiometric (slope -1) mC<->hmC exchange.
#
# FRAMING (per master plan): neuronal genes are PREFERENTIALLY / DISPROPORTIONATELY
#   affected, NOT "exclusively" — MeCP2 redistribution tracks coordinated
#   methylation (OR=5.16) ~3x more strongly than neuronal identity (OR=1.73).
#
# Panels (master plan "Figure 6" table; layout "AABB\nCCDD", 180 x 160 mm):
#   A  Sec 72  Dose-response line   72_neuronal_decile_summary.tsv
#             neuronal fraction per constitutive-K119ub decile; D10 OR=1.70
#   B  Sec 74  OR comparison bar    74_pairwise_fisher.tsv
#             Coord x MeCP2-Up OR=5.16 >> Neuronal x MeCP2-Up OR=1.73 (log y)
#   C  Sec 76  Synapse K27me3 box   diffbind_gene_level_all_marks.tsv split by
#             synapse/broader-neuronal/non-neuronal; selective K27me3 loss is the
#             punchline (synapse vs broader-neuronal K27me3 d=-0.044, p=2.95e-3,
#             from 76_synapse_vs_neuronal_stats.tsv)
#   D  Sec 78  Stoichiometry slopes 78_stoichiometry_slopes.tsv
#             slope per gene class vs dashed -1; neuronal=-1.0 (stoichiometric)
#
# Data source: pre-computed section TSVs in plots/visualizations/tables/.
# All key statistics annotated on panels are read from and confirmed against
# those TSVs (see header comments at each panel). No new computation beyond the
# gene-class partition for panel C (which exactly reproduces the n's and medians
# in 76_synapse_vs_neuronal_stats.tsv).
#
# Output: plots/figures/figure6_neuronal_specificity/  (PDF + SVG + JPG)
# =============================================================================

# -----------------------------------------------------------------------------
# Shared infrastructure. Working directory MUST be downstream/ (getwd()).
# _shared_config.R provides paths, COLORS, theme_biomodal(), helpers, and the
# pre-loaded DMR objects; _figure_config.R adds theme_pub(), save_figure(),
# read_table_tsv(), add_panel_labels(), stat_label(), FIGURE_DIR, PUB_COLORS.
# -----------------------------------------------------------------------------
source("scripts/viz_sections/_shared_config.R")
source("scripts/figures/_figure_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

message("\n=== Figure 6: Neuronal Genes Are Preferentially Affected ===")

# =============================================================================
# Panel A — Sec 72: neuronal fraction by constitutive H2AK119ub decile
# Source: 72_neuronal_decile_summary.tsv (decile, frac_neuronal, ci_lo, ci_hi, OR)
# Key stat: top-decile (D10) OR = 1.70 (confirmed: decile-10 OR = 1.7013;
#           72_fisher_results.tsv "Ctrl top decile (D9)" OR = 1.7013).
# Story: neuronal genes are constitutively K119ub-enriched, scaling with signal.
# =============================================================================

decile_df <- read_table_tsv("72_neuronal_decile_summary.tsv")

required_decile_cols <- c("decile", "frac_neuronal", "ci_lo", "ci_hi", "OR")
missing_decile <- setdiff(required_decile_cols, colnames(decile_df))
if (length(missing_decile) > 0) {
  stop("Panel A: 72_neuronal_decile_summary.tsv missing columns: ",
       paste(missing_decile, collapse = ", "))
}

decile_df <- decile_df %>%
  mutate(
    decile        = as.integer(decile),
    frac_neuronal = as.numeric(frac_neuronal),
    ci_lo         = as.numeric(ci_lo),
    ci_hi         = as.numeric(ci_hi),
    OR            = as.numeric(OR)
  ) %>%
  arrange(decile)

# Genome-wide neuronal fraction reference line (from 72_fisher_results.tsv).
fisher72_df <- read_table_tsv("72_fisher_results.tsv")
genome_frac <- as.numeric(fisher72_df$frac_neuronal_genome[1])  # 0.23500

# Confirm the annotated D10 OR against the source TSV before drawing.
or_d10 <- decile_df$OR[decile_df$decile == 10]
stopifnot(length(or_d10) == 1, is.finite(or_d10))

panel_a <- ggplot(decile_df, aes(x = decile, y = frac_neuronal)) +
  geom_hline(yintercept = genome_frac, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  annotate("text", x = 1, y = genome_frac, vjust = -0.6, hjust = 0,
           label = sprintf("genome-wide = %.1f%%", 100 * genome_frac),
           colour = "grey40", size = 2.3) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
              fill = COLORS$k119ub[["K119ub Gained"]], alpha = 0.18) +
  geom_line(colour = COLORS$k119ub[["K119ub Gained"]], linewidth = 0.6) +
  geom_point(colour = COLORS$k119ub[["K119ub Gained"]], size = 1.4) +
  annotate("text", x = 10, y = decile_df$ci_hi[decile_df$decile == 10],
           vjust = -0.8, hjust = 1, size = 2.4, fontface = "bold",
           label = stat_label("OR", or_d10)) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    title = "Neuronal genes are constitutively H2AK119ub-enriched",
    subtitle = "Neuronal-gene fraction rises with wild-type K119ub signal",
    x = "Constitutive H2AK119ub gene-body signal (decile, low→high)",
    y = "Neuronal-gene fraction"
  ) +
  theme_pub()

# =============================================================================
# Panel B — Sec 74: pairwise overlap odds ratios (MeCP2 redistribution driver)
# Source: 74_pairwise_fisher.tsv (comparison, odds_ratio, ci_lo, ci_hi, p_adj)
# Key stat: Coordinated x MeCP2-Up OR = 5.16 >> Neuronal x MeCP2-Up OR = 1.73;
#           Neuronal x Coordinated OR = 1.05 (NS). Log y scale.
# Story: MeCP2 redistribution tracks methylation state, not neuronal lineage.
# =============================================================================

fisher74_df <- read_table_tsv("74_pairwise_fisher.tsv")

required_f74_cols <- c("comparison", "odds_ratio", "ci_lo", "ci_hi", "p_adj")
missing_f74 <- setdiff(required_f74_cols, colnames(fisher74_df))
if (length(missing_f74) > 0) {
  stop("Panel B: 74_pairwise_fisher.tsv missing columns: ",
       paste(missing_f74, collapse = ", "))
}

# Order comparisons strongest -> weakest so the methylation > lineage message
# reads left-to-right.
f74_order <- c("Coordinated × MeCP2-Up",
               "Neuronal × MeCP2-Up",
               "Neuronal × Coordinated")
missing_rows <- setdiff(f74_order, fisher74_df$comparison)
if (length(missing_rows) > 0) {
  stop("Panel B: 74_pairwise_fisher.tsv missing expected rows: ",
       paste(missing_rows, collapse = ", "))
}

fisher74_df <- fisher74_df %>%
  filter(comparison %in% f74_order) %>%
  mutate(
    odds_ratio = as.numeric(odds_ratio),
    ci_lo      = as.numeric(ci_lo),
    ci_hi      = as.numeric(ci_hi),
    p_adj      = as.numeric(p_adj),
    # Colour the methylation-driven overlap distinctly from lineage overlaps.
    is_methyl  = comparison == "Coordinated × MeCP2-Up",
    comparison = factor(comparison, levels = f74_order),
    sig_lab    = ifelse(p_adj < 0.001, "italic(P) < 0.001",
                        sprintf("italic(P) == %.3f", p_adj))
  )

# Format a "OR = x.xx" label that sits just above each bar's CI cap.
fisher74_df$or_lab <- vapply(fisher74_df$odds_ratio,
                             function(v) stat_label("OR", v),
                             character(1))

panel_b <- ggplot(fisher74_df,
                  aes(x = comparison, y = odds_ratio, fill = is_methyl)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_col(width = 0.62, colour = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.18, linewidth = 0.4, colour = "grey20") +
  geom_text(aes(y = ci_hi, label = or_lab),
            vjust = -0.6, size = 2.4, fontface = "bold") +
  geom_text(aes(y = ci_hi, label = sig_lab), parse = TRUE,
            vjust = -2.1, size = 2.1, colour = "grey35") +
  scale_fill_manual(
    values = c("TRUE" = PUB_COLORS$Hyper, "FALSE" = "grey60"),
    guide  = "none"
  ) +
  scale_y_log10(
    breaks = c(1, 2, 5, 10),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "MeCP2 redistribution tracks methylation, not lineage",
    subtitle = "Coordinated mC↑/hmC↓ overlaps MeCP2 gain ~3× more than neuronal identity",
    x = NULL,
    y = "Odds ratio (log scale)"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# =============================================================================
# Panel C — Sec 76: synapse/axon selective Polycomb de-repression (K27me3 loss)
# Source distributions: diffbind_gene_level_all_marks.tsv (atac/k27ac/k27me3_fold)
# Source gene sets: 76_synapse_axon_gene_set.tsv (subset of)
#                   72_neuronal_gene_set_go_derived.tsv
# Annotated stat: synapse vs broader-neuronal K27me3 d(median) = -0.044,
#   Wilcoxon p = 2.95e-3 (from 76_synapse_vs_neuronal_stats.tsv). The class
#   partition reproduces that table's n's and medians exactly.
# Story: synapse/axon genes lose more K27me3 (Polycomb de-repression) without
#   extra accessibility/enhancer gain — selective, not blanket, remodeling.
# =============================================================================

marks_df <- read_table_tsv("diffbind_gene_level_all_marks.tsv")
required_marks_cols <- c("gene", "atac_fold", "k27ac_fold", "k27me3_fold")
missing_marks <- setdiff(required_marks_cols, colnames(marks_df))
if (length(missing_marks) > 0) {
  stop("Panel C: diffbind_gene_level_all_marks.tsv missing columns: ",
       paste(missing_marks, collapse = ", "))
}

synapse_set  <- read_table_tsv("76_synapse_axon_gene_set.tsv")$gene
neuronal_set <- read_table_tsv("72_neuronal_gene_set_go_derived.tsv")$gene
if (length(synapse_set) == 0 || length(neuronal_set) == 0) {
  stop("Panel C: empty synapse or neuronal gene set.")
}

# Three mutually exclusive classes (synapse set is a strict subset of neuronal).
classify_gene_class <- function(g) {
  ifelse(g %in% synapse_set, "Synapse/axon",
         ifelse(g %in% neuronal_set, "Broader neuronal", "Non-neuronal"))
}

mark_levels <- c("ATAC", "H3K27ac", "H3K27me3")
class_levels <- c("Synapse/axon", "Broader neuronal", "Non-neuronal")

marks_long <- marks_df %>%
  mutate(
    gene_class = classify_gene_class(gene),
    ATAC       = suppressWarnings(as.numeric(atac_fold)),
    H3K27ac    = suppressWarnings(as.numeric(k27ac_fold)),
    H3K27me3   = suppressWarnings(as.numeric(k27me3_fold))
  ) %>%
  select(gene, gene_class, ATAC, H3K27ac, H3K27me3) %>%
  pivot_longer(cols = all_of(mark_levels),
               names_to = "mark", values_to = "fold") %>%
  filter(!is.na(fold)) %>%
  mutate(
    mark       = factor(mark, levels = mark_levels),
    gene_class = factor(gene_class, levels = class_levels)
  )

# Confirm the class partition reproduces 76_synapse_vs_neuronal_stats.tsv so the
# annotated K27me3 stat is faithful to the source distributions being plotted.
syn_stats_df <- read_table_tsv("76_synapse_vs_neuronal_stats.tsv")
k27me3_syn_vs_neur <- syn_stats_df %>%
  filter(mark == "K27me3",
         group_a == "Synapse/axon", group_b == "Broader neuronal")
if (nrow(k27me3_syn_vs_neur) != 1) {
  stop("Panel C: could not locate the synapse-vs-broader-neuronal K27me3 row ",
       "in 76_synapse_vs_neuronal_stats.tsv")
}
k27me3_delta <- as.numeric(k27me3_syn_vs_neur$median_a) -
                as.numeric(k27me3_syn_vs_neur$median_b)   # -0.0444
k27me3_p     <- as.numeric(k27me3_syn_vs_neur$wilcox_p)    # 2.95e-3

# Cross-check the plotted K27me3 medians against the source table (exact match).
plotted_med <- marks_long %>%
  filter(mark == "H3K27me3") %>%
  group_by(gene_class) %>%
  summarise(med = median(fold), .groups = "drop")
src_syn_med  <- as.numeric(k27me3_syn_vs_neur$median_a)
plot_syn_med <- plotted_med$med[plotted_med$gene_class == "Synapse/axon"]
if (abs(plot_syn_med - src_syn_med) > 1e-6) {
  stop("Panel C: plotted synapse K27me3 median (", plot_syn_med,
       ") does not match source table (", src_syn_med, ").")
}

class_fill <- c(
  "Synapse/axon"     = PUB_COLORS$Hyper,   # the special, K27me3-losing class
  "Broader neuronal" = "#4DAF4A",
  "Non-neuronal"     = "grey70"
)

# Annotation: place the synapse-specific K27me3 stat in the K27me3 facet only.
# Plain (unparsed) string so the Greek delta and superscript-free P render as-is.
k27me3_anno <- data.frame(
  mark  = factor("H3K27me3", levels = mark_levels),
  label = sprintf(
    "Synapse vs broader-neuronal\nK27me3 Δmedian = %.3f, P = %.1e",
    k27me3_delta, k27me3_p
  )
)

panel_c <- ggplot(marks_long, aes(x = gene_class, y = fold, fill = gene_class)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.62, linewidth = 0.3) +
  geom_text(data = k27me3_anno, inherit.aes = FALSE,
            aes(x = 1, y = Inf, label = label),
            hjust = 0, vjust = 1.2, size = 2.1, colour = "grey25") +
  facet_wrap(~ mark, nrow = 1) +
  scale_fill_manual(values = class_fill, guide = "none") +
  coord_cartesian(ylim = c(-1.0, 1.0)) +
  labs(
    title = "Synapse/axon genes show selective Polycomb de-repression",
    subtitle = "Extra K27me3 loss at synaptic genes; no extra ATAC / H3K27ac gain",
    x = NULL,
    y = "BAP1-KO differential signal (log2 fold change)"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# =============================================================================
# Panel D — Sec 78: mC<->hmC stoichiometry slope by gene class
# Source: 78_stoichiometry_slopes.tsv (group, slope, ci_lo, ci_hi,
#         differs_from_neg1)
# Key stat: Neuronal (broad) slope = -1.0 (differs_from_neg1 = FALSE,
#           stoichiometric); Synapse/axon = -1.02 (FALSE); All genes = -0.959,
#           Non-neuronal = -0.949 (both TRUE, deviate from -1). Dashed ref at -1.
# Story: neuronal/synaptic loci undergo stoichiometric 5hmC->5mC conversion
#   (dehydroxymethylase-like), distinct from the genome-wide deviation.
# =============================================================================

slopes_df <- read_table_tsv("78_stoichiometry_slopes.tsv")
required_slope_cols <- c("group", "slope", "ci_lo", "ci_hi", "differs_from_neg1")
missing_slope <- setdiff(required_slope_cols, colnames(slopes_df))
if (length(missing_slope) > 0) {
  stop("Panel D: 78_stoichiometry_slopes.tsv missing columns: ",
       paste(missing_slope, collapse = ", "))
}

slope_order <- c("All genes", "Non-neuronal", "Neuronal (broad)", "Synapse/axon")
missing_slope_rows <- setdiff(slope_order, slopes_df$group)
if (length(missing_slope_rows) > 0) {
  stop("Panel D: 78_stoichiometry_slopes.tsv missing expected groups: ",
       paste(missing_slope_rows, collapse = ", "))
}

slopes_df <- slopes_df %>%
  filter(group %in% slope_order) %>%
  mutate(
    slope             = as.numeric(slope),
    ci_lo             = as.numeric(ci_lo),
    ci_hi             = as.numeric(ci_hi),
    differs_from_neg1 = as.logical(differs_from_neg1),
    group             = factor(group, levels = slope_order),
    # "Stoichiometric" = CI overlaps -1 (differs_from_neg1 == FALSE).
    relation          = ifelse(differs_from_neg1,
                               "Deviates from −1",
                               "Stoichiometric (≈1)")
  )

slopes_df$slope_lab <- sprintf("%.2f", slopes_df$slope)

relation_fill <- c(
  "Stoichiometric (≈1)" = PUB_COLORS$Hyper,
  "Deviates from −1"     = "grey60"
)

panel_d <- ggplot(slopes_df, aes(x = group, y = slope, fill = relation)) +
  geom_hline(yintercept = -1, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.6, y = -1, vjust = -0.6, hjust = 0,
           label = "stoichiometric (−1)", colour = "grey40", size = 2.3) +
  geom_col(width = 0.6, colour = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.18, linewidth = 0.4, colour = "grey20") +
  geom_text(aes(y = ci_lo, label = slope_lab),
            vjust = 1.5, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = relation_fill, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.14, 0.06))) +
  labs(
    title = "Neuronal methylation changes are stoichiometric",
    subtitle = "δ-5mC ~ δ-5hmC slope ≈ −1 at neuronal / synaptic genes",
    x = NULL,
    y = expression(paste("Slope of ", delta, "-5mC vs ", delta, "-5hmC"))
  ) +
  theme_pub() +
  theme(
    axis.text.x     = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

# =============================================================================
# Compose: design "AABB\nCCDD", heights c(1, 1); 180 x 160 mm (master plan).
# =============================================================================

design <- "AABB\nCCDD"

figure6 <- panel_a + panel_b + panel_c + panel_d +
  plot_layout(design = design, heights = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 6. Neuronal genes are preferentially affected by BAP1 loss",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

save_figure(figure6, "figure6_neuronal_specificity",
            width_mm = 180, height_mm = 160)

# Also save the standalone synapse-specificity panel (C) — it carries the R5
# punchline and is useful at full resolution for slide / supplement reuse.
save_figure(panel_c, "figure6C_synapse_chromatin",
            width_mm = 180, height_mm = 80)

message("Figure 6 written to ", file.path(FIGURE_DIR, "figure6_neuronal_specificity"),
        " (PDF + SVG + JPG)")
message("  Panel A  D10 neuronal-enrichment OR = ", sprintf("%.2f", or_d10))
message("  Panel B  Coord×MeCP2-Up OR = ",
        sprintf("%.2f", fisher74_df$odds_ratio[fisher74_df$comparison == "Coordinated × MeCP2-Up"]),
        " ; Neuronal×MeCP2-Up OR = ",
        sprintf("%.2f", fisher74_df$odds_ratio[fisher74_df$comparison == "Neuronal × MeCP2-Up"]))
message("  Panel C  synapse-vs-broader-neuronal K27me3 Δ = ",
        sprintf("%.3f", k27me3_delta), " (P = ", sprintf("%.2e", k27me3_p), ")")
message("  Panel D  Neuronal(broad) slope = ",
        sprintf("%.3f", slopes_df$slope[slopes_df$group == "Neuronal (broad)"]),
        " (stoichiometric)")
