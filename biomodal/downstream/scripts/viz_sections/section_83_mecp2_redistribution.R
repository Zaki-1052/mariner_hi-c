# scripts/viz_sections/section_83_mecp2_redistribution.R
# =============================================================================
# FIGURE 4 — MeCP2 Redistribution and Developmental Amplification (paper R1)
#
# Question: How does MeCP2 binding change in BAP1-KO cerebellum? It redistributes
# to distal-intergenic sites (vacating gene bodies), preferentially overlaps
# hypermethylated DMRs, and the normal developmental MeCP2 ramp is amplified in
# the mutant. The new panel F asks whether the mutant-unique developmental
# (aging) MeCP2 gain co-localizes with the coordinated mC up / hmC down DMR set.
#
# MeCP2 data is CUT&RUN (Epicypher), NOT ChIP — labels say "signal"/"binding"/
# "mark", never "ChIP".
#
# Panels (per approved master plan, design = "AABB\nCCDD\nEEFF"):
#   A | Sec 11 | MeCP2 adult volcano       | MECP2_FILES$annotated (202,650 peaks)
#               7,686 UP / 1,200 DOWN sig peaks; UP orange / DOWN purple.
#   B | Sec 75 | Peak-annotation fraction  | 75_peak_annotation_distribution.tsv
#               UP 51.7% distal-intergenic vs 2.2% promoter (vs DOWN 61.8% intron).
#   C | Sec 11 | DMR x MeCP2 2x2 overlap   | mecp2_dmr_overlap_summary.tsv
#               Fisher OR = 5.13, p = 1.27e-33 (hyper DMRs -> MeCP2-Up peaks).
#   D | Sec 77 | Aging peak counts         | 77_aging_peak_summary.tsv
#               Mutant 23,117 vs Control 10,930 aging-UP peaks (2.1x).
#   E | Sec 77 | Shared-peak fold scatter  | 77_shared_peak_fold_comparison.tsv
#               7,305 shared loci; median fold 2.24 (mut) vs 1.83 (ctrl), +22%.
#   F | NEW    | Aging x methylation Fisher| re-derived from MeCP2 aging DiffBind
#               (MECP2_FILES$ctrl_aging / $mut_aging) x mecp2_coordinated_genes.tsv
#               over the 20,915-gene universe. Saved to
#               plots/figures/tables/fig4f_aging_methylation_overlap.tsv.
#
# Dimensions: 180 x 220 mm.
#
# Sourcing contract (run from downstream/, getwd() == BASE_DIR):
#   1. scripts/viz_sections/_shared_config.R  (data objects, COLORS, helpers,
#      ChIPseeker/TxDb, save_multiformat_ggplot via multi_format_output.R)
#   2. scripts/viz_sections/_figure_config.R   (theme_pub, save_figure,
#      read_table_tsv, add_panel_labels, stat_label, FIGURE_DIR/_TABLES_DIR)
# =============================================================================

source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

SEC83_DIR <- file.path(OUTPUT_DIR, "83_mecp2_redistribution")
dir.create(SEC83_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC83_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  # ChIPseeker is NOT loaded by _shared_config.R (only GenomicRanges, IRanges,
  # org.Mm.eg.db, and TxDb.Mmusculus.UCSC.mm10.knownGene are). Panel F's aging
  # x methylation analysis annotates peaks to genes via annotatePeak(), so load
  # it here — mirroring section_77_mecp2_aging_trajectory.R.
  library(ChIPseeker)
})

cat("\n=============================================================\n")
cat("FIGURE 4 — MeCP2 Redistribution and Developmental Amplification\n")
cat("=============================================================\n\n")

# Master font size for panel-level on-plot annotations (geom_text size unit ~pt/2.8)
ANNO_SIZE <- 2.4
LABEL_SIZE <- 2.2

# =============================================================================
# PANEL A — MeCP2 adult volcano (Section 11)
# Source: MECP2_FILES$annotated (REPO_ROOT/peaks/mecp2/MeCP2_annotated.txt)
# 202,650 peaks; columns seqnames/start/end/Fold/p.value/FDR/annotation/SYMBOL.
# Fold > 0 in this file = higher in mutant (Conc_mut > Conc_ctrl); UP = MeCP2 Up.
# Verified: 7,686 sig UP, 1,200 sig DOWN at FDR < 0.05.
# =============================================================================

cat("--- Panel A: MeCP2 adult volcano ---\n")

mecp2_anno_path <- MECP2_FILES$annotated
if (!file.exists(mecp2_anno_path)) {
  stop("Panel A: MeCP2 annotated peak file not found: ", mecp2_anno_path)
}

# load_diffbind_flex standardizes seqnames -> Chr and validates Fold/FDR.
mecp2_peaks <- load_diffbind_flex(mecp2_anno_path, mark_name = "MeCP2 adult", fdr_threshold = Q_THRESHOLD)

# Pull SYMBOL through for top-gene labeling (load_diffbind_flex preserves all cols).
if (!"SYMBOL" %in% colnames(mecp2_peaks)) {
  stop("Panel A: expected SYMBOL column in annotated MeCP2 peaks.")
}

mecp2_peaks <- mecp2_peaks %>%
  mutate(
    neg_log10_fdr = -log10(pmax(FDR, 1e-300)),
    direction = case_when(
      FDR < Q_THRESHOLD & Fold > 0 ~ "MeCP2 Up",
      FDR < Q_THRESHOLD & Fold < 0 ~ "MeCP2 Down",
      TRUE ~ "Not Significant"
    )
  )
mecp2_peaks$direction <- factor(
  mecp2_peaks$direction,
  levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant")
)

n_up_peaks   <- sum(mecp2_peaks$direction == "MeCP2 Up")
n_down_peaks <- sum(mecp2_peaks$direction == "MeCP2 Down")
cat(sprintf("  MeCP2 adult peaks: %d UP / %d DOWN (FDR < %.2f)\n",
            n_up_peaks, n_down_peaks, Q_THRESHOLD))
stopifnot(n_up_peaks == 7686, n_down_peaks == 1200)

# Label the strongest significant peaks (largest |Fold| among significant peaks
# with a valid gene symbol) plus any of the paper KEY_GENES that are significant.
sig_peaks <- mecp2_peaks %>%
  filter(direction != "Not Significant",
         !is.na(SYMBOL), SYMBOL != "")

top_by_fold <- sig_peaks %>%
  arrange(desc(abs(Fold))) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  slice_head(n = 12)

key_gene_peaks <- sig_peaks %>%
  filter(SYMBOL %in% KEY_GENES) %>%
  arrange(FDR) %>%
  distinct(SYMBOL, .keep_all = TRUE)

label_peaks <- bind_rows(top_by_fold, key_gene_peaks) %>%
  distinct(SYMBOL, .keep_all = TRUE)

fdr_line_y <- -log10(Q_THRESHOLD)

p_A <- ggplot(mecp2_peaks, aes(x = Fold, y = neg_log10_fdr, colour = direction)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_hline(yintercept = fdr_line_y, linetype = "dashed",
             colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "grey60", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = label_peaks,
    aes(label = SYMBOL),
    size = LABEL_SIZE, colour = "black",
    max.overlaps = 20, segment.color = NA,
    box.padding = 0.3, show.legend = FALSE
  ) +
  scale_colour_manual(values = COLORS$mecp2, name = NULL) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("%s UP\n%s DOWN",
                           format(n_up_peaks, big.mark = ","),
                           format(n_down_peaks, big.mark = ",")),
           hjust = 1.05, vjust = 1.2, size = ANNO_SIZE, lineheight = 0.95) +
  labs(
    title = "Adult MeCP2 binding (CUT&RUN)",
    x = expression(log[2]~"fold change (mutant / control)"),
    y = expression(-log[10]~"FDR")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  theme_pub()

# =============================================================================
# PANEL B — Peak-annotation distribution (Section 75)
# Source: 75_peak_annotation_distribution.tsv (peak_status, annotation_broad, n, pct).
# UP peaks 51.7% Distal Intergenic, 41.0% Intron, 2.2% Promoter;
# DOWN peaks 61.8% Intron, 19.5% Distal Intergenic, 8.0% Promoter.
# =============================================================================

cat("--- Panel B: MeCP2 peak-annotation distribution ---\n")

peak_anno <- read_table_tsv("75_peak_annotation_distribution.tsv")
stopifnot(all(c("peak_status", "annotation_broad", "n", "pct") %in% colnames(peak_anno)))

# Confirm the load-bearing values match the plan before annotating.
up_distal <- peak_anno %>%
  filter(peak_status == "UP", annotation_broad == "Distal Intergenic") %>% pull(pct)
up_promoter <- peak_anno %>%
  filter(peak_status == "UP", annotation_broad == "Promoter") %>% pull(pct)
cat(sprintf("  UP distal-intergenic = %.1f%%, UP promoter = %.1f%%\n",
            up_distal, up_promoter))
stopifnot(round(up_distal, 1) == 51.7, round(up_promoter, 1) == 2.2)

# Order annotation categories from most-genic (promoter) to most-distal so the
# stacked bar reads top-down as "MeCP2 gain is distal".
anno_levels <- c("Promoter", "Exon/UTR", "Intron", "Downstream", "Distal Intergenic")
peak_anno <- peak_anno %>%
  mutate(
    annotation_broad = factor(annotation_broad, levels = anno_levels),
    peak_status = factor(peak_status, levels = c("UP", "DOWN"),
                         labels = c("MeCP2 Up", "MeCP2 Down")),
    frac = pct / 100
  )

anno_palette <- c(
  "Promoter"          = "#1B9E77",
  "Exon/UTR"          = "#66A61E",
  "Intron"            = "#7570B3",
  "Downstream"        = "#E7298A",
  "Distal Intergenic" = "#D95F02"
)

p_B <- ggplot(peak_anno, aes(x = peak_status, y = frac, fill = annotation_broad)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = anno_palette, name = "Genomic feature",
                    drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  # "Distal Intergenic" is the last factor level, so geom_col stacks it at the
  # bottom of the bar (y in [0, up_distal/100]); centre the label on that segment.
  annotate("text", x = "MeCP2 Up", y = (up_distal / 100) / 2,
           label = sprintf("%.0f%%\ndistal", up_distal),
           size = ANNO_SIZE, colour = "white", fontface = "bold",
           lineheight = 0.9) +
  labs(
    title = "MeCP2 gains land at distal sites",
    x = NULL, y = "Fraction of significant peaks"
  ) +
  theme_pub()

# =============================================================================
# PANEL C — DMR x MeCP2 2x2 overlap (Section 11)
# Source: mecp2_dmr_overlap_summary.tsv
# (DMR_Direction, MeCP2_Direction, Overlap_Count, Total_DMRs, Percentage,
#  Fisher_OR, Fisher_pvalue).
# Hyper DMRs preferentially overlap MeCP2-Up peaks: OR = 5.13, p = 1.27e-33.
# =============================================================================

cat("--- Panel C: DMR x MeCP2 overlap ---\n")

dmr_overlap <- read_table_tsv("mecp2_dmr_overlap_summary.tsv")
stopifnot(all(c("DMR_Direction", "MeCP2_Direction", "Percentage",
                "Fisher_OR", "Fisher_pvalue") %in% colnames(dmr_overlap)))

fisher_or  <- unique(dmr_overlap$Fisher_OR)
fisher_p   <- unique(dmr_overlap$Fisher_pvalue)
stopifnot(length(fisher_or) == 1, length(fisher_p) == 1)
cat(sprintf("  DMR x MeCP2 Fisher OR = %.2f, p = %.2e\n", fisher_or, fisher_p))
stopifnot(round(fisher_or, 2) == 5.13)

dmr_overlap <- dmr_overlap %>%
  mutate(
    DMR_Direction = factor(DMR_Direction,
                           levels = c("Hypermethylated", "Hypomethylated")),
    MeCP2_Direction = factor(MeCP2_Direction,
                             levels = c("MeCP2 Up", "MeCP2 Down")),
    frac = Percentage / 100
  )

p_C <- ggplot(dmr_overlap,
              aes(x = DMR_Direction, y = frac, fill = MeCP2_Direction)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)),
            position = position_dodge(width = 0.75), vjust = -0.3,
            size = ANNO_SIZE) +
  scale_fill_manual(values = COLORS$mecp2[c("MeCP2 Up", "MeCP2 Down")],
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  annotate("text", x = 1.5, y = max(dmr_overlap$frac) * 1.12,
           label = sprintf("Fisher %s\n%s",
                           stat_label("OR", fisher_or),
                           stat_label("p", fisher_p, "%.1e")),
           size = ANNO_SIZE, lineheight = 0.95) +
  labs(
    title = "Hyper DMRs overlap MeCP2 gains",
    x = "DMR direction (5mC)",
    y = "Peaks overlapping DMRs"
  ) +
  theme_pub()

# =============================================================================
# PANEL D — Aging peak counts (Section 77)
# Source: 77_aging_peak_summary.tsv (genotype, status, count).
# Aging-UP: Control 10,930 vs Mutant 23,117 (2.1x); DOWN 2,822 vs 10,646.
# (NS is dropped from the bar — it dwarfs the differential peaks.)
# =============================================================================

cat("--- Panel D: MeCP2 aging peak counts ---\n")

aging_summary <- read_table_tsv("77_aging_peak_summary.tsv")
stopifnot(all(c("genotype", "status", "count") %in% colnames(aging_summary)))

aging_diff <- aging_summary %>%
  filter(status %in% c("UP", "DOWN")) %>%
  mutate(
    genotype = factor(genotype, levels = c("Control", "Mutant")),
    status   = factor(status, levels = c("UP", "DOWN"),
                      labels = c("Age-increased", "Age-decreased"))
  )

ctrl_up <- aging_summary %>% filter(genotype == "Control", status == "UP") %>% pull(count)
mut_up  <- aging_summary %>% filter(genotype == "Mutant",  status == "UP") %>% pull(count)
up_ratio <- mut_up / ctrl_up
cat(sprintf("  Aging-UP: ctrl = %d, mut = %d (%.2fx)\n", ctrl_up, mut_up, up_ratio))
stopifnot(ctrl_up == 10930, mut_up == 23117)

p_D <- ggplot(aging_diff, aes(x = status, y = count, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = format(count, big.mark = ",")),
            position = position_dodge(width = 0.7), vjust = -0.3,
            size = ANNO_SIZE) +
  scale_fill_manual(values = COLORS$condition[c("Control", "Mutant")],
                    name = NULL) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  annotate("text", x = 1, y = mut_up * 1.10,
           label = sprintf("%.1fx more\nin mutant", up_ratio),
           size = ANNO_SIZE, lineheight = 0.95) +
  labs(
    title = "Developmental MeCP2 ramp is amplified",
    subtitle = "Young -> adult DiffBind, per genotype",
    x = "Age-related binding change",
    y = "Significant peaks"
  ) +
  theme_pub()

# =============================================================================
# PANEL E — Shared-peak aging-fold scatter (Section 77)
# Source: 77_shared_peak_fold_comparison.tsv (ctrl_fold, mut_fold, gene, is_neuronal).
# 7,305 loci where both genotypes age-gain MeCP2; median fold 1.83 (ctrl) vs
# 2.24 (mut) — mutant climbs higher (+22%) at shared loci. 2,266 are neuronal.
# =============================================================================

cat("--- Panel E: shared aging-peak fold scatter ---\n")

shared_fold <- read_table_tsv("77_shared_peak_fold_comparison.tsv")
stopifnot(all(c("ctrl_fold", "mut_fold", "gene", "is_neuronal") %in% colnames(shared_fold)))

n_shared <- nrow(shared_fold)
med_ctrl <- median(shared_fold$ctrl_fold, na.rm = TRUE)
med_mut  <- median(shared_fold$mut_fold,  na.rm = TRUE)
n_neuronal <- sum(shared_fold$is_neuronal == TRUE | shared_fold$is_neuronal == "TRUE",
                  na.rm = TRUE)
cat(sprintf("  Shared loci = %d; median fold ctrl = %.3f, mut = %.3f; neuronal = %d\n",
            n_shared, med_ctrl, med_mut, n_neuronal))
stopifnot(n_shared == 7305, round(med_ctrl, 2) == 1.83, round(med_mut, 2) == 2.24)

shared_fold <- shared_fold %>%
  mutate(neuronal_lab = ifelse(is_neuronal == TRUE | is_neuronal == "TRUE",
                               "Neuronal gene", "Other gene"))

fold_lim <- range(c(shared_fold$ctrl_fold, shared_fold$mut_fold), na.rm = TRUE)

p_E <- ggplot(shared_fold, aes(x = ctrl_fold, y = mut_fold, colour = neuronal_lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.3) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = med_mut, linetype = "dotted",
             colour = COLORS$condition[["Mutant"]], linewidth = 0.3) +
  geom_vline(xintercept = med_ctrl, linetype = "dotted",
             colour = COLORS$condition[["Control"]], linewidth = 0.3) +
  scale_colour_manual(values = c("Neuronal gene" = "#1B9E77",
                                 "Other gene" = "grey60"), name = NULL) +
  coord_equal(xlim = fold_lim, ylim = fold_lim) +
  annotate("text", x = fold_lim[1], y = fold_lim[2],
           label = sprintf("%s shared loci\nmedian mut %.2f vs ctrl %.2f",
                           format(n_shared, big.mark = ","), med_mut, med_ctrl),
           hjust = 0, vjust = 1, size = ANNO_SIZE, lineheight = 0.95) +
  labs(
    title = "Mutant ages higher at shared loci",
    x = expression("Control aging fold ("*log[2]*" adult/young)"),
    y = expression("Mutant aging fold ("*log[2]*" adult/young)")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  theme_pub()

# =============================================================================
# COMPOSE — 5 panels (A-E). Panel F (aging x methylation overlap) is now
# standalone section 84.
# =============================================================================

cat("\n--- Composing Section 83 ---\n")

layout_design <- "AABB\nCCDD\n#EE#"

sec83_composite <- p_A + p_B + p_C + p_D + p_E +
  patchwork::plot_layout(design = layout_design, heights = c(1.2, 1, 0.8))

sec83_composite <- sec83_composite +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Section 83: MeCP2 redistribution and developmental amplification",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

save_plot(p_A, "83a_mecp2_volcano",     w = 8, h = 7)
save_plot(p_B, "83b_peak_annotation",   w = 8, h = 6)
save_plot(p_C, "83c_dmr_overlap",       w = 8, h = 6)
save_plot(p_D, "83d_aging_peak_counts", w = 8, h = 6)
save_plot(p_E, "83e_shared_fold",       w = 8, h = 7)
save_plot(sec83_composite, "83_composite", w = 14, h = 14)

cat("\nSection 83 complete.\n")
