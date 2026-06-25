# scripts/viz_sections/section_82_neuronal_specificity.R
# Section 82: Neuronal Genes Are Preferentially Affected
#
# Consolidated section drawing from sections 72, 74, 76, 78.
# Neuronal genes are PREFERENTIALLY / DISPROPORTIONATELY affected, NOT
# "exclusively" -- MeCP2 redistribution tracks coordinated methylation
# (OR=5.16) ~3x more strongly than neuronal identity (OR=1.73).
#
# Panels:
#   A  Sec 72  K119ub decile dose-response (72_neuronal_decile_summary.tsv)
#   B  Sec 74  OR comparison bar (74_pairwise_fisher.tsv)
#   C  Sec 76  Synapse K27me3 boxplot (diffbind_gene_level_all_marks.tsv)
#   D  Sec 78  Stoichiometry slopes (78_stoichiometry_slopes.tsv)
#
# Output: plots/visualizations/82_neuronal_specificity/
# =============================================================================

source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

SEC82_DIR <- file.path(OUTPUT_DIR, "82_neuronal_specificity")
dir.create(SEC82_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC82_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

message("\n=== Section 82: Neuronal Genes Are Preferentially Affected ===")

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
  dplyr::mutate(
    decile        = as.integer(decile),
    frac_neuronal = as.numeric(frac_neuronal),
    ci_lo         = as.numeric(ci_lo),
    ci_hi         = as.numeric(ci_hi),
    OR            = as.numeric(OR)
  ) %>%
  dplyr::arrange(decile)

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
    x = "Constitutive H2AK119ub gene-body signal (decile, low->high)",
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
  dplyr::filter(comparison %in% f74_order) %>%
  dplyr::mutate(
    odds_ratio = as.numeric(odds_ratio),
    ci_lo      = as.numeric(ci_lo),
    ci_hi      = as.numeric(ci_hi),
    p_adj      = as.numeric(p_adj),
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
    values = c("TRUE" = COLORS$direction[["Hypermethylated"]], "FALSE" = "grey60"),
    guide  = "none"
  ) +
  scale_y_log10(
    breaks = c(1, 2, 5, 10),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "MeCP2 redistribution tracks methylation, not lineage",
    subtitle = "Coordinated mC-up/hmC-down overlaps MeCP2 gain ~3x more than neuronal identity",
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
  dplyr::mutate(
    gene_class = classify_gene_class(gene),
    ATAC       = suppressWarnings(as.numeric(atac_fold)),
    H3K27ac    = suppressWarnings(as.numeric(k27ac_fold)),
    H3K27me3   = suppressWarnings(as.numeric(k27me3_fold))
  ) %>%
  dplyr::select(gene, gene_class, ATAC, H3K27ac, H3K27me3) %>%
  tidyr::pivot_longer(cols = tidyr::all_of(mark_levels),
                      names_to = "mark", values_to = "fold") %>%
  dplyr::filter(!is.na(fold)) %>%
  dplyr::mutate(
    mark       = factor(mark, levels = mark_levels),
    gene_class = factor(gene_class, levels = class_levels)
  )

# Confirm the class partition reproduces 76_synapse_vs_neuronal_stats.tsv so the
# annotated K27me3 stat is faithful to the source distributions being plotted.
syn_stats_df <- read_table_tsv("76_synapse_vs_neuronal_stats.tsv")
k27me3_syn_vs_neur <- syn_stats_df %>%
  dplyr::filter(mark == "K27me3",
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
  dplyr::filter(mark == "H3K27me3") %>%
  dplyr::group_by(gene_class) %>%
  dplyr::summarise(med = median(fold), .groups = "drop")
src_syn_med  <- as.numeric(k27me3_syn_vs_neur$median_a)
plot_syn_med <- plotted_med$med[plotted_med$gene_class == "Synapse/axon"]
if (abs(plot_syn_med - src_syn_med) > 1e-6) {
  stop("Panel C: plotted synapse K27me3 median (", plot_syn_med,
       ") does not match source table (", src_syn_med, ").")
}

class_fill <- c(
  "Synapse/axon"     = COLORS$direction[["Hypermethylated"]],   # the special, K27me3-losing class
  "Broader neuronal" = "#4DAF4A",
  "Non-neuronal"     = "grey70"
)

# Annotation: place the synapse-specific K27me3 stat in the K27me3 facet only.
# Plain (unparsed) string so the Greek delta and superscript-free P render as-is.
k27me3_anno <- data.frame(
  mark  = factor("H3K27me3", levels = mark_levels),
  label = sprintf(
    "Synapse vs broader-neuronal\nK27me3 delta-median = %.3f, P = %.1e",
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
  dplyr::filter(group %in% slope_order) %>%
  dplyr::mutate(
    slope             = as.numeric(slope),
    ci_lo             = as.numeric(ci_lo),
    ci_hi             = as.numeric(ci_hi),
    differs_from_neg1 = as.logical(differs_from_neg1),
    group             = factor(group, levels = slope_order),
    relation          = ifelse(differs_from_neg1,
                               "Deviates from -1",
                               "Stoichiometric (~-1)")
  )

slopes_df$slope_lab <- sprintf("%.2f", slopes_df$slope)

relation_fill <- c(
  "Stoichiometric (~-1)" = COLORS$direction[["Hypermethylated"]],
  "Deviates from -1"     = "grey60"
)

panel_d <- ggplot(slopes_df, aes(x = group, y = slope, fill = relation)) +
  geom_hline(yintercept = -1, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.6, y = -1, vjust = -0.6, hjust = 0,
           label = "stoichiometric (-1)", colour = "grey40", size = 2.3) +
  geom_col(width = 0.6, colour = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.18, linewidth = 0.4, colour = "grey20") +
  geom_text(aes(y = ci_lo, label = slope_lab),
            vjust = 1.5, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = relation_fill, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.14, 0.06))) +
  labs(
    title = "Neuronal methylation changes are stoichiometric",
    subtitle = "d-5mC ~ d-5hmC slope ~ -1 at neuronal / synaptic genes",
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

sec82_composite <- panel_a + panel_b + panel_c + panel_d +
  plot_layout(design = design, heights = c(1, 1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Section 82: Neuronal genes are preferentially affected by BAP1 loss",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

save_plot(panel_a, "82a_k119ub_decile", w = 8, h = 6)
save_plot(panel_b, "82b_fisher_or_comparison", w = 8, h = 6)
save_plot(panel_c, "82c_synapse_chromatin", w = 14, h = 6)
save_plot(panel_d, "82d_stoichiometry_slopes", w = 8, h = 6)
save_plot(sec82_composite, "82_composite", w = 14, h = 10)

message("\nSection 82 written to ", SEC82_DIR, " (PDF + SVG + JPG)")
message("  Panel A  D10 neuronal-enrichment OR = ", sprintf("%.2f", or_d10))
message("  Panel B  Coord x MeCP2-Up OR = ",
        sprintf("%.2f", fisher74_df$odds_ratio[fisher74_df$comparison == "Coordinated × MeCP2-Up"]),
        " ; Neuronal x MeCP2-Up OR = ",
        sprintf("%.2f", fisher74_df$odds_ratio[fisher74_df$comparison == "Neuronal × MeCP2-Up"]))
message("  Panel C  synapse-vs-broader-neuronal K27me3 delta = ",
        sprintf("%.3f", k27me3_delta), " (P = ", sprintf("%.2e", k27me3_p), ")")
message("  Panel D  Neuronal(broad) slope = ",
        sprintf("%.3f", slopes_df$slope[slopes_df$group == "Neuronal (broad)"]),
        " (stoichiometric)")
