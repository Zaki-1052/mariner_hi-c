# biomodal/downstream/scripts/viz_sections/section_39_h3k36me2_boundary_analysis.R
# H3K36me2 Polycomb Boundary Analysis (NSD/DNMT3A axis)
#
# H3K36me2 is deposited by NSD1/NSD2, broadly distributed across intergenic and
# genic regions. It recruits DNMT3A and acts as a barrier preventing PRC2/H3K27me3
# from spreading into active domains. This section tests whether that barrier is
# compromised in BAP1-KO, and whether me2 changes correlate with methylation changes.
#
# Data: DiffBind results (403,123 peaks, NO gene annotations — requires ChIPseeker)
# Figures: 39a-39g (7 total)
#
# Usage:
#   cd downstream/
#   Rscript scripts/viz_sections/section_39_h3k36me2_boundary_analysis.R

source("scripts/viz_sections/_shared_config.R")
suppressPackageStartupMessages(library(ChIPseeker))

SECTION_DIR <- "39_h3k36me2_boundary"
FDR_THRESHOLD <- 0.05

dir.create(file.path(OUTPUT_DIR, SECTION_DIR), recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("================================================================================\n")
cat("SECTION 39: H3K36me2 Polycomb Boundary Analysis (NSD/DNMT3A Axis)\n")
cat("================================================================================\n\n")

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading H3K36me2 DiffBind data...\n")
me2_all <- load_diffbind_flex(H3K36ME2_FILES$master, "H3K36me2 (all peaks)")

# Load demethylation ratio table
RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
stopifnot(file.exists(RATIO_TABLE))
ratio_data <- read.table(RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Ratio table: %d genes loaded\n", nrow(ratio_data)))

# Load differential BED files for coordinate-level overlap
cat("Loading H3K36me2 differential BEDs...\n")
me2_up_gr <- load_chip_peaks(H3K36ME2_FILES$up, "H3K36me2 Up")
me2_down_gr <- load_chip_peaks(H3K36ME2_FILES$down, "H3K36me2 Down")

# Load H3K27me3 peaks for boundary analysis
cat("Loading H3K27me3 peaks for boundary analysis...\n")
k27me3_gr <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")

# Load H3K27me3 DiffBind results for cross-mark comparison
cat("Loading H3K27me3 DiffBind data...\n")
k27me3_diffbind <- load_diffbind_flex(DIFFBIND_FILES$k27me3, "H3K27me3 DiffBind")

# =============================================================================
# FILTER TO SIGNIFICANT AND ANNOTATE VIA CHIPSEEKER
# =============================================================================

cat("\nAnnotating significant H3K36me2 peaks with ChIPseeker...\n")

me2_sig <- me2_all %>% dplyr::filter(FDR < FDR_THRESHOLD)
cat(sprintf("  %d significant peaks to annotate\n", nrow(me2_sig)))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

me2_sig_gr <- GRanges(
  seqnames = me2_sig$Chr,
  ranges = IRanges(start = me2_sig$Start, end = me2_sig$End),
  Fold = me2_sig$Fold,
  FDR = me2_sig$FDR,
  Conc = me2_sig$Conc
)

me2_anno <- suppressMessages(annotatePeak(
  me2_sig_gr, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"
))
me2_anno_df <- as.data.frame(me2_anno)
cat(sprintf("  Annotated %d peaks\n", nrow(me2_anno_df)))

# Gene-level aggregation (nearest-to-TSS)
me2_gene <- me2_anno_df %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    me2_fold = Fold[which.min(abs(distanceToTSS))],
    me2_fdr  = FDR[which.min(abs(distanceToTSS))],
    me2_n_peaks = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::rename(gene = SYMBOL)

cat(sprintf("  %d unique genes with significant H3K36me2 peaks\n", nrow(me2_gene)))

# Merge with methylation data
merged <- ratio_data %>%
  left_join(me2_gene, by = "gene")

n_with_me2 <- sum(!is.na(merged$me2_fold))
cat(sprintf("  %d genes with both methylation and me2 data (%.1f%%)\n",
            n_with_me2, 100 * n_with_me2 / nrow(merged)))

# =============================================================================
# FIGURE 39a: H3K36me2 DiffBind Overview (Volcano + Annotation Bar)
# =============================================================================

cat("\n--- FIGURE 39a: H3K36me2 DiffBind Overview ---\n\n")

me2_all$sig <- me2_all$FDR < FDR_THRESHOLD
me2_all$neg_log10_fdr <- -log10(pmax(me2_all$FDR, 1e-300))
me2_all$neg_log10_fdr[is.infinite(me2_all$neg_log10_fdr)] <- 300

me2_all$color_group <- case_when(
  me2_all$sig & me2_all$Fold > 0 ~ "Up (Gained)",
  me2_all$sig & me2_all$Fold < 0 ~ "Down (Lost)",
  TRUE ~ "Not Significant"
)

n_up <- sum(me2_all$color_group == "Up (Gained)")
n_down <- sum(me2_all$color_group == "Down (Lost)")
n_ns <- sum(me2_all$color_group == "Not Significant")

# Subsample non-significant heavily (403K peaks)
set.seed(42)
ns_subset <- me2_all %>% dplyr::filter(color_group == "Not Significant") %>% dplyr::sample_frac(0.05)
sig_peaks <- me2_all %>% dplyr::filter(color_group != "Not Significant")
plot_data_volcano <- bind_rows(ns_subset, sig_peaks)

p_volcano <- ggplot(plot_data_volcano, aes(x = Fold, y = neg_log10_fdr, color = color_group)) +
  geom_point(alpha = 0.4, size = 0.5) +
  scale_color_manual(
    values = c("Up (Gained)" = "#E6AB02", "Down (Lost)" = "#66A61E", "Not Significant" = "grey70"),
    name = "Direction"
  ) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "H3K36me2 DiffBind: All Peaks",
    subtitle = sprintf("%s up, %s down, %s NS (FDR<%.2f) | 5%% NS subsampled",
                       format(n_up, big.mark = ","), format(n_down, big.mark = ","),
                       format(n_ns, big.mark = ","), FDR_THRESHOLD),
    x = "Log2 Fold Change (Mutant / Control)",
    y = "-log10(FDR)"
  ) +
  theme_biomodal()

# Annotation bar chart for significant peaks (up vs down)
anno_df <- me2_anno_df %>%
  mutate(
    direction = ifelse(Fold > 0, "Gained", "Lost"),
    anno_simple = case_when(
      grepl("Promoter", annotation) ~ "Promoter",
      grepl("Exon", annotation)     ~ "Exon",
      grepl("3' UTR|5' UTR", annotation) ~ "UTR",
      grepl("Intron", annotation)   ~ "Intron",
      grepl("Downstream", annotation) ~ "Downstream",
      grepl("Intergenic", annotation) ~ "Intergenic",
      TRUE ~ "Other"
    )
  )

anno_bar_data <- anno_df %>%
  count(direction, anno_simple) %>%
  group_by(direction) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

anno_colors <- c("Promoter" = "#E41A1C", "Exon" = "#377EB8", "Intron" = "#4DAF4A",
                 "UTR" = "#984EA3", "Intergenic" = "#A65628", "Downstream" = "#F781BF",
                 "Other" = "#999999")

p_anno_bar <- ggplot(anno_bar_data, aes(x = direction, y = pct, fill = anno_simple)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = anno_colors, name = "Annotation") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Annotation: Up vs Down", x = "", y = "% of Peaks") +
  theme_biomodal() +
  theme(axis.text.x = element_text(face = "bold", size = 11))

p_39a <- p_volcano + p_anno_bar + plot_layout(widths = c(2, 1))

save_multiformat_ggplot(p_39a, file.path(OUTPUT_DIR, SECTION_DIR, "39a_h3k36me2_volcano_annotation"),
                        width = 14, height = 7)

# =============================================================================
# FIGURE 39b: O/E Heatmap — me2 Direction x mC Direction
# =============================================================================

cat("\n--- FIGURE 39b: me2 x mC Direction O/E Heatmap ---\n\n")

# Reuse the build_oe_heatmap from conceptual pattern
build_oe_heatmap <- function(direction_a, direction_b, levels_a, levels_b,
                             title, x_label, y_label, predicted_pairs = list()) {
  quad_table <- table(factor(direction_a, levels = levels_a),
                      factor(direction_b, levels = levels_b))
  fisher_result <- fisher.test(quad_table)
  total_n <- sum(quad_table)
  row_totals <- rowSums(quad_table)
  col_totals <- colSums(quad_table)
  expected <- outer(row_totals, col_totals) / total_n
  enrichment <- quad_table / expected

  heatmap_df <- as.data.frame(as.table(quad_table))
  colnames(heatmap_df) <- c("row_dir", "col_dir", "Observed")
  expected_df <- as.data.frame(as.table(round(expected, 1)))
  colnames(expected_df) <- c("row_dir", "col_dir", "Expected")
  enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
  colnames(enrichment_df) <- c("row_dir", "col_dir", "Enrichment")

  heatmap_data <- heatmap_df %>%
    left_join(expected_df, by = c("row_dir", "col_dir")) %>%
    left_join(enrichment_df, by = c("row_dir", "col_dir")) %>%
    mutate(
      label = sprintf("Obs: %d\nExp: %.0f\nO/E: %.2f", Observed, Expected, Enrichment),
      row_dir = factor(row_dir, levels = levels_a),
      col_dir = factor(col_dir, levels = levels_b),
      is_predicted = FALSE
    )

  for (pair in predicted_pairs) {
    heatmap_data$is_predicted[heatmap_data$row_dir == pair[1] &
                              heatmap_data$col_dir == pair[2]] <- TRUE
  }

  sub_text <- sprintf("Fisher OR = %.2f, p = %.2e (n = %d genes)",
                      fisher_result$estimate, fisher_result$p.value, total_n)
  if (length(predicted_pairs) > 0) {
    sub_text <- paste0(sub_text, " | Black borders = predicted")
  }

  p <- ggplot(heatmap_data, aes(x = col_dir, y = row_dir, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2)

  if (any(heatmap_data$is_predicted)) {
    p <- p + geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
                       aes(x = col_dir, y = row_dir),
                       fill = NA, color = "black", linewidth = 1.5, linetype = "solid")
  }

  p <- p +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "O/E",
                         limits = c(min(heatmap_data$Enrichment) * 0.9,
                                    max(heatmap_data$Enrichment) * 1.1)) +
    labs(title = title, subtitle = sub_text, x = x_label, y = y_label) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 11, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  list(plot = p, data = heatmap_data, fisher = fisher_result, table = quad_table)
}

oe_mc_df <- merged %>%
  dplyr::filter(mc_sig & !is.na(me2_fold)) %>%
  mutate(
    mc_direction = ifelse(mc_diff > 0, "mC Up", "mC Down"),
    me2_direction = ifelse(me2_fold > 0, "me2 Gained", "me2 Lost")
  )

cat(sprintf("  Genes with both sig mC DMR and sig me2: %d\n", nrow(oe_mc_df)))

if (nrow(oe_mc_df) >= 4) {
  result_39b <- build_oe_heatmap(
    oe_mc_df$mc_direction, oe_mc_df$me2_direction,
    c("mC Up", "mC Down"), c("me2 Gained", "me2 Lost"),
    title = "H3K36me2 Direction x 5mC Direction (Gene Level)",
    x_label = "H3K36me2 Direction", y_label = "5mC DMR Direction",
    predicted_pairs = list(c("mC Up", "me2 Gained"))
  )

  save_multiformat_ggplot(result_39b$plot,
    file.path(OUTPUT_DIR, SECTION_DIR, "39b_me2_x_mc_oe_heatmap"),
    width = 8, height = 6)

  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_39b$fisher$estimate, result_39b$fisher$p.value))
} else {
  cat("  WARNING: Too few genes for mC x me2 O/E analysis\n")
}

# =============================================================================
# FIGURE 39c: H3K36me2 at H3K27me3 Boundary Zones
# =============================================================================

cat("\n--- FIGURE 39c: me2 Fold at H3K27me3 Boundaries ---\n\n")

# Define H3K27me3 boundary zones (flanking regions ±5kb from peak edges)
BOUNDARY_FLANK <- 5000

k27me3_boundaries <- c(
  flank(k27me3_gr, width = BOUNDARY_FLANK, start = TRUE),
  flank(k27me3_gr, width = BOUNDARY_FLANK, start = FALSE)
)
k27me3_boundaries <- reduce(k27me3_boundaries)

# All significant me2 peaks as GRanges
me2_sig_peaks <- GRanges(
  seqnames = me2_sig$Chr,
  ranges = IRanges(start = me2_sig$Start, end = me2_sig$End),
  Fold = me2_sig$Fold
)

# Classify peaks as boundary vs non-boundary
at_boundary <- countOverlaps(me2_sig_peaks, k27me3_boundaries) > 0

boundary_df <- data.frame(
  fold = me2_sig$Fold,
  location = ifelse(at_boundary, "At H3K27me3\nBoundary", "Away from\nH3K27me3")
)

n_boundary <- sum(at_boundary)
n_away <- sum(!at_boundary)
cat(sprintf("  me2 peaks at boundary: %d (%.1f%%)\n", n_boundary, 100 * n_boundary / nrow(me2_sig)))

wilcox_boundary <- wilcox.test(
  boundary_df$fold[boundary_df$location == "At H3K27me3\nBoundary"],
  boundary_df$fold[boundary_df$location == "Away from\nH3K27me3"]
)

p_39c <- ggplot(boundary_df, aes(x = location, y = fold, fill = location)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("At H3K27me3\nBoundary" = "#B2182B", "Away from\nH3K27me3" = "#999999")) +
  annotate("text", x = 1.5, y = max(boundary_df$fold, na.rm = TRUE) * 0.95,
           label = sprintf("Wilcoxon p = %.2e\nBoundary n=%d, Away n=%d",
                           wilcox_boundary$p.value, n_boundary, n_away),
           size = 3.5, fontface = "italic", lineheight = 0.9) +
  labs(
    title = "H3K36me2 Fold Change: H3K27me3 Boundary vs Non-Boundary",
    subtitle = sprintf("Boundary = within %dkb of H3K27me3 peak edges", BOUNDARY_FLANK / 1000),
    x = "", y = "H3K36me2 Log2 Fold Change"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_39c, file.path(OUTPUT_DIR, SECTION_DIR, "39c_me2_at_k27me3_boundary"),
                        width = 8, height = 7)

cat(sprintf("  Wilcoxon p = %.2e\n", wilcox_boundary$p.value))
cat(sprintf("  Median fold — Boundary: %+.3f, Away: %+.3f\n",
            median(boundary_df$fold[boundary_df$location == "At H3K27me3\nBoundary"], na.rm = TRUE),
            median(boundary_df$fold[boundary_df$location == "Away from\nH3K27me3"], na.rm = TRUE)))

# =============================================================================
# FIGURE 39d: O/E Heatmap — me2 Direction x K27me3 Direction (Mutual Exclusivity)
# =============================================================================

cat("\n--- FIGURE 39d: me2 x K27me3 Direction O/E Heatmap ---\n\n")

# Annotate K27me3 DiffBind peaks to genes (same approach)
k27me3_valid <- k27me3_diffbind %>%
  dplyr::filter(!is.na(k27me3_diffbind$Chr) & !is.na(k27me3_diffbind$Start) & !is.na(k27me3_diffbind$End))

k27me3_gr_db <- GRanges(
  seqnames = k27me3_valid$Chr,
  ranges = IRanges(start = k27me3_valid$Start, end = k27me3_valid$End),
  Fold = k27me3_valid$Fold,
  FDR = k27me3_valid$FDR
)

k27me3_anno <- suppressMessages(annotatePeak(
  k27me3_gr_db, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"
))
k27me3_anno_df <- as.data.frame(k27me3_anno)

k27me3_gene <- k27me3_anno_df %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    k27me3_fold = Fold[which.min(abs(distanceToTSS))],
    k27me3_fdr  = FDR[which.min(abs(distanceToTSS))],
    .groups = "drop"
  ) %>%
  dplyr::rename(gene = SYMBOL)

# Cross-mark analysis: genes with both sig me2 and sig K27me3
cross_df <- me2_gene %>%
  inner_join(k27me3_gene, by = "gene") %>%
  dplyr::filter(me2_fdr < FDR_THRESHOLD & k27me3_fdr < FDR_THRESHOLD) %>%
  mutate(
    me2_direction = ifelse(me2_fold > 0, "me2 Gained", "me2 Lost"),
    k27me3_direction = ifelse(k27me3_fold > 0, "K27me3 Gained", "K27me3 Lost")
  )

cat(sprintf("  Genes with both sig me2 and sig K27me3: %d\n", nrow(cross_df)))

if (nrow(cross_df) >= 4) {
  # Predicted: me2 Loss + K27me3 Gain (PRC2 barrier breach)
  result_39d <- build_oe_heatmap(
    cross_df$me2_direction, cross_df$k27me3_direction,
    c("me2 Gained", "me2 Lost"), c("K27me3 Gained", "K27me3 Lost"),
    title = "H3K36me2 x H3K27me3 Direction (Mutual Exclusivity Test)",
    x_label = "H3K27me3 Direction", y_label = "H3K36me2 Direction",
    predicted_pairs = list(c("me2 Lost", "K27me3 Gained"))
  )

  save_multiformat_ggplot(result_39d$plot,
    file.path(OUTPUT_DIR, SECTION_DIR, "39d_me2_x_k27me3_oe_heatmap"),
    width = 8, height = 6)

  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_39d$fisher$estimate, result_39d$fisher$p.value))
} else {
  cat("  WARNING: Too few genes for me2 x K27me3 cross-mark analysis\n")
}

# =============================================================================
# FIGURE 39e: Scatter — me2 Fold vs mC/hmC Diff
# =============================================================================

cat("\n--- FIGURE 39e: me2 Fold vs Methylation Change Scatters ---\n\n")

scatter_df <- merged %>%
  dplyr::filter(!is.na(me2_fold))

scatter_df$is_key <- scatter_df$gene %in% KEY_GENES

# mC scatter
cor_mc <- cor.test(scatter_df$me2_fold, scatter_df$mc_diff, method = "spearman")

p_scatter_mc <- ggplot(scatter_df, aes(x = me2_fold, y = mc_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.3, size = 0.8, color = "grey50") +
  geom_smooth(method = "lm", color = "#E6AB02", fill = "#E6AB02", alpha = 0.15, linewidth = 0.8) +
  {if (any(scatter_df$is_key))
    ggrepel::geom_text_repel(data = scatter_df %>% dplyr::filter(is_key),
                             aes(label = gene), size = 2.5, fontface = "italic",
                             max.overlaps = 15, segment.size = 0.2)
  } +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3.5,
           label = sprintf("rho = %+.3f\np = %.2e\nn = %d",
                           cor_mc$estimate, cor_mc$p.value, nrow(scatter_df))) +
  labs(title = "H3K36me2 Fold vs 5mC Change", x = "H3K36me2 Log2FC", y = "5mC Difference") +
  theme_biomodal() +
  theme(legend.position = "none")

# hmC scatter
cor_hmc <- cor.test(scatter_df$me2_fold, scatter_df$hmc_diff, method = "spearman")

p_scatter_hmc <- ggplot(scatter_df, aes(x = me2_fold, y = hmc_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.3, size = 0.8, color = "grey50") +
  geom_smooth(method = "lm", color = "#66A61E", fill = "#66A61E", alpha = 0.15, linewidth = 0.8) +
  {if (any(scatter_df$is_key))
    ggrepel::geom_text_repel(data = scatter_df %>% dplyr::filter(is_key),
                             aes(label = gene), size = 2.5, fontface = "italic",
                             max.overlaps = 15, segment.size = 0.2)
  } +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3.5,
           label = sprintf("rho = %+.3f\np = %.2e\nn = %d",
                           cor_hmc$estimate, cor_hmc$p.value, nrow(scatter_df))) +
  labs(title = "H3K36me2 Fold vs 5hmC Change", x = "H3K36me2 Log2FC", y = "5hmC Difference") +
  theme_biomodal() +
  theme(legend.position = "none")

p_39e <- p_scatter_mc + p_scatter_hmc +
  plot_annotation(title = "Gene-Level H3K36me2 Fold vs Methylation Changes",
                  subtitle = "Spearman correlation | Key genes labeled",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                                plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")))

save_multiformat_ggplot(p_39e, file.path(OUTPUT_DIR, SECTION_DIR, "39e_me2_vs_methylation_scatter"),
                        width = 14, height = 7)

cat(sprintf("  mC correlation:  rho = %+.3f, p = %.2e\n", cor_mc$estimate, cor_mc$p.value))
cat(sprintf("  hmC correlation: rho = %+.3f, p = %.2e\n", cor_hmc$estimate, cor_hmc$p.value))

# =============================================================================
# FIGURE 39f: me2 Status by Chromatin State at DMRs
# =============================================================================

cat("\n--- FIGURE 39f: me2 Status by Chromatin State ---\n\n")

state_mc_df <- merged %>%
  dplyr::filter(mc_sig & !is.na(chromatin_state)) %>%
  mutate(
    me2_status = case_when(
      is.na(me2_fold)     ~ "No Sig me2",
      me2_fold > 0        ~ "me2 Gained",
      me2_fold < 0        ~ "me2 Lost",
      TRUE                ~ "me2 Unchanged"
    ),
    chromatin_state = factor(chromatin_state, levels = CHROMATIN_STATE_ORDER)
  )

state_bar_data <- state_mc_df %>%
  count(chromatin_state, me2_status) %>%
  group_by(chromatin_state) %>%
  mutate(pct = 100 * n / sum(n), total = sum(n)) %>%
  ungroup()

me2_status_colors <- c("me2 Gained" = "#E6AB02", "me2 Lost" = "#66A61E",
                        "me2 Unchanged" = "#BDBDBD", "No Sig me2" = "#636363")

p_39f <- ggplot(state_bar_data, aes(x = chromatin_state, y = pct, fill = me2_status)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  geom_text(data = state_bar_data %>% group_by(chromatin_state) %>%
              summarise(total = total[1], .groups = "drop"),
            aes(x = chromatin_state, y = 105, label = sprintf("n=%d", total)),
            inherit.aes = FALSE, size = 2.8, color = "grey40") +
  scale_fill_manual(values = me2_status_colors, name = "H3K36me2 Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "H3K36me2 Status at Significant mC DMRs by Chromatin State",
    subtitle = "Fraction of DMR genes with me2 Gained / Lost / No significant change",
    x = "Chromatin State", y = "% of Genes"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10))

save_multiformat_ggplot(p_39f, file.path(OUTPUT_DIR, SECTION_DIR, "39f_me2_by_chromatin_state"),
                        width = 11, height = 7)

# =============================================================================
# FIGURE 39g: me2 Fold in Genic vs Intergenic Compartments
# =============================================================================

cat("\n--- FIGURE 39g: me2 Fold by Genomic Compartment ---\n\n")

compartment_df <- me2_anno_df %>%
  mutate(
    compartment = case_when(
      grepl("Intergenic|Downstream", annotation) ~ "Intergenic",
      grepl("Promoter|Exon|Intron|UTR", annotation) ~ "Genic",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(compartment %in% c("Genic", "Intergenic"))

wilcox_compartment <- wilcox.test(
  compartment_df$Fold[compartment_df$compartment == "Genic"],
  compartment_df$Fold[compartment_df$compartment == "Intergenic"]
)

comp_counts <- compartment_df %>% count(compartment)

p_39g <- ggplot(compartment_df, aes(x = compartment, y = Fold, fill = compartment)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Genic" = "#E6AB02", "Intergenic" = "#66A61E")) +
  annotate("text", x = 1.5, y = max(compartment_df$Fold, na.rm = TRUE) * 0.95,
           label = sprintf("Wilcoxon p = %.2e\nGenic n=%d, Intergenic n=%d",
                           wilcox_compartment$p.value,
                           comp_counts$n[comp_counts$compartment == "Genic"],
                           comp_counts$n[comp_counts$compartment == "Intergenic"]),
           size = 3.5, fontface = "italic", lineheight = 0.9) +
  labs(
    title = "H3K36me2 Fold Change: Genic vs Intergenic Peaks",
    subtitle = "Genic = promoter/exon/intron/UTR | Intergenic = distal intergenic/downstream",
    x = "", y = "H3K36me2 Log2 Fold Change"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_39g, file.path(OUTPUT_DIR, SECTION_DIR, "39g_me2_genic_vs_intergenic"),
                        width = 8, height = 7)

cat(sprintf("  Wilcoxon p = %.2e\n", wilcox_compartment$p.value))

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting Tables ---\n\n")

# Gene-level summary
export_gene <- me2_gene %>%
  left_join(ratio_data %>% dplyr::select(gene, mc_diff, mc_q, mc_sig, hmc_diff, hmc_q, hmc_sig,
                                          chromatin_state, dmr_status),
            by = "gene")
write.table(export_gene, file.path(TABLES_DIR, "h3k36me2_gene_level_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved h3k36me2_gene_level_summary.tsv (%d rows)\n", nrow(export_gene)))

# Boundary analysis summary
boundary_summary <- data.frame(
  location = c("At H3K27me3 Boundary", "Away from H3K27me3"),
  n_peaks = c(n_boundary, n_away),
  median_fold = c(
    median(boundary_df$fold[boundary_df$location == "At H3K27me3\nBoundary"], na.rm = TRUE),
    median(boundary_df$fold[boundary_df$location == "Away from\nH3K27me3"], na.rm = TRUE)
  ),
  wilcoxon_p = wilcox_boundary$p.value
)
write.table(boundary_summary, file.path(TABLES_DIR, "h3k36me2_k27me3_boundary_analysis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved h3k36me2_k27me3_boundary_analysis.tsv\n"))

cat("\n================================================================================\n")
cat("SECTION 39 COMPLETE\n")
cat("================================================================================\n")
