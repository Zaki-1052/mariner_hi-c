# biomodal/downstream/scripts/viz_sections/section_31_mecp2_loop_anchor_integration.R
# Section 31: MeCP2 x Hi-C Loop Anchor Integration
# Standalone script - sources shared config for all dependencies and data
#
# Tests whether MeCP2 differential binding associates with Hi-C loop anchor
# positions. Connects the methylation reader (MeCP2) to 3D architecture.
#
# Biological prediction (Mellen et al. 2017):
#   MeCP2 binds 5mC with HIGH affinity, 5hmC with LOW affinity.
#   BAP1 loss -> TET blocked -> 5mC up / 5hmC down at gene bodies.
#   Section 27 showed: lost-loop anchors are 43% hypermethylated.
#   Therefore: lost-loop anchors with more mC should show increased MeCP2.
#
# Sub-analyses:
#   31a: Direct genomic overlap — MeCP2 peaks at loop anchors
#   31b: MeCP2 direction x loop direction heatmap (gene-level)
#   31c: MeCP2 fold change by loop direction (violin/boxplot)
#   31d: Loop logFC vs MeCP2 fold change scatter
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_31_mecp2_loop_anchor_integration.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# MeCP2-SPECIFIC HELPER FUNCTIONS (self-contained for section independence)
# =============================================================================

#' Load MeCP2 annotated differential binding data
load_mecp2_annotated <- function(filepath) {
  stopifnot("MeCP2 annotated file not found" = file.exists(filepath))

  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   fill = TRUE, quote = "")

  numeric_cols <- c("start", "end", "width", "Conc", "Conc_mut", "Conc_ctrl",
                    "Fold", "p.value", "FDR", "geneStart", "geneEnd",
                    "geneLength", "distanceToTSS")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }

  cat(sprintf("  Loaded %d MeCP2 peaks\n", nrow(df)))
  cat(sprintf("  Significant (FDR < 0.05): %d up, %d down\n",
              sum(df$FDR < 0.05 & df$Fold > 0, na.rm = TRUE),
              sum(df$FDR < 0.05 & df$Fold < 0, na.rm = TRUE)))

  return(df)
}

#' Aggregate MeCP2 peaks per gene
aggregate_mecp2_by_gene <- function(mecp2_df) {
  gene_summary <- mecp2_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      mean_fold = mean(Fold, na.rm = TRUE),
      max_fold = max(Fold, na.rm = TRUE),
      min_fold = min(Fold, na.rm = TRUE),
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      nearest_fdr = FDR[which.min(abs(distanceToTSS))],
      min_distance_tss = min(abs(distanceToTSS), na.rm = TRUE),
      n_peaks = n(),
      any_sig = any(FDR < 0.05),
      .groups = "drop"
    )

  return(gene_summary)
}

# =============================================================================
# CONFIGURATION AND INPUT VALIDATION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 31: MeCP2 x HI-C LOOP ANCHOR INTEGRATION\n")
cat("================================================================================\n\n")

# Validate all inputs — no fallbacks
stopifnot("MeCP2 annotated file not found" = file.exists(MECP2_FILES$annotated))
stopifnot("MeCP2 up peaks not found" = file.exists(MECP2_FILES$up))
stopifnot("MeCP2 down peaks not found" = file.exists(MECP2_FILES$down))
stopifnot("Loop annotation file not found" = file.exists(LOOP_FILES$late))

ANCHOR_GENE_TABLE <- file.path(TABLES_DIR, "anchor_gene_associations_great.tsv")
stopifnot("anchor_gene_associations_great.tsv not found -- run Section 27 first" =
            file.exists(ANCHOR_GENE_TABLE))

# Section output directory
SECTION_DIR <- file.path(OUTPUT_DIR, "31_mecp2_loop_anchor_integration")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# Loop direction display labels and colors
LOOP_DIRECTION_COLORS <- c("Gained" = "#4575b4", "Lost" = "#d73027", "Background" = "grey70")

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# DATA LOADING
# =============================================================================

# Load MeCP2 data
cat("Loading MeCP2 differential binding data...\n")
mecp2_annotated <- load_mecp2_annotated(MECP2_FILES$annotated)
mecp2_up_gr <- load_chip_peaks(MECP2_FILES$up, "MeCP2 Up")
mecp2_down_gr <- load_chip_peaks(MECP2_FILES$down, "MeCP2 Down")

# Aggregate to gene level
cat("\nAggregating MeCP2 peaks per gene...\n")
mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated)
cat(sprintf("  %d unique genes with MeCP2 peaks\n", nrow(mecp2_gene)))

# Load loop data
cat("\nLoading Hi-C loop annotations...\n")
loops <- read.table(LOOP_FILES$late, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d differential loops\n", nrow(loops)))
cat(sprintf("    up_in_mutant (Gained):   %d\n", sum(loops$direction == "up_in_mutant")))
cat(sprintf("    down_in_mutant (Lost):   %d\n", sum(loops$direction == "down_in_mutant")))

# Load pre-computed anchor-gene GREAT associations
cat("\nLoading anchor-gene GREAT associations...\n")
anchor_genes <- read.table(ANCHOR_GENE_TABLE, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE)
cat(sprintf("  %d anchor-gene associations (%d unique anchors, %d unique genes)\n",
            nrow(anchor_genes), length(unique(anchor_genes$anchor_id)),
            length(unique(anchor_genes$symbol))))

# =============================================================================
# FIGURE 31a: Direct Genomic Overlap — MeCP2 Peaks at Loop Anchors
# =============================================================================
cat("\nCreating Figure 31a: MeCP2 peak overlap at loop anchors by loop direction...\n")

# Build anchor GRanges
anchor1_gr <- GRanges(seqnames = loops$chr1,
                      ranges = IRanges(start = loops$start1, end = loops$end1))
anchor2_gr <- GRanges(seqnames = loops$chr2,
                      ranges = IRanges(start = loops$start2, end = loops$end2))

# Compute MeCP2 overlap at each anchor
loops$anchor1_mecp2_up   <- countOverlaps(anchor1_gr, mecp2_up_gr) > 0
loops$anchor1_mecp2_down <- countOverlaps(anchor1_gr, mecp2_down_gr) > 0
loops$anchor2_mecp2_up   <- countOverlaps(anchor2_gr, mecp2_up_gr) > 0
loops$anchor2_mecp2_down <- countOverlaps(anchor2_gr, mecp2_down_gr) > 0

# Pool anchors: any anchor with MeCP2 overlap
loops$any_mecp2_up   <- loops$anchor1_mecp2_up | loops$anchor2_mecp2_up
loops$any_mecp2_down <- loops$anchor1_mecp2_down | loops$anchor2_mecp2_down

cat(sprintf("\n  Anchor MeCP2 overlap summary:\n"))
cat(sprintf("    Loops with any anchor MeCP2 Up:   %d/%d (%.1f%%)\n",
            sum(loops$any_mecp2_up), nrow(loops),
            100 * sum(loops$any_mecp2_up) / nrow(loops)))
cat(sprintf("    Loops with any anchor MeCP2 Down: %d/%d (%.1f%%)\n",
            sum(loops$any_mecp2_down), nrow(loops),
            100 * sum(loops$any_mecp2_down) / nrow(loops)))

# Stratify by loop direction
up_loops <- loops[loops$direction == "up_in_mutant", ]
dn_loops <- loops[loops$direction == "down_in_mutant", ]

# Compute overlap rates
overlap_stats <- data.frame(
  loop_direction = rep(c("Gained", "Lost"), each = 2),
  mecp2_direction = rep(c("MeCP2 Up", "MeCP2 Down"), 2),
  count = c(
    sum(up_loops$any_mecp2_up), sum(up_loops$any_mecp2_down),
    sum(dn_loops$any_mecp2_up), sum(dn_loops$any_mecp2_down)
  ),
  total = rep(c(nrow(up_loops), nrow(dn_loops)), each = 2),
  stringsAsFactors = FALSE
)
overlap_stats$pct <- 100 * overlap_stats$count / overlap_stats$total

overlap_stats$loop_direction <- factor(overlap_stats$loop_direction,
                                       levels = c("Gained", "Lost"))
overlap_stats$mecp2_direction <- factor(overlap_stats$mecp2_direction,
                                         levels = c("MeCP2 Up", "MeCP2 Down"))

for (i in 1:nrow(overlap_stats)) {
  cat(sprintf("  %s loops, %s: %d/%d (%.1f%%)\n",
              overlap_stats$loop_direction[i], overlap_stats$mecp2_direction[i],
              overlap_stats$count[i], overlap_stats$total[i],
              overlap_stats$pct[i]))
}

# Fisher's test: gained loops vs MeCP2-up enrichment
fisher_up <- fisher.test(matrix(
  c(sum(up_loops$any_mecp2_up), nrow(up_loops) - sum(up_loops$any_mecp2_up),
    sum(dn_loops$any_mecp2_up), nrow(dn_loops) - sum(dn_loops$any_mecp2_up)),
  nrow = 2
))
# Fisher's test: lost loops vs MeCP2-down enrichment
fisher_dn <- fisher.test(matrix(
  c(sum(dn_loops$any_mecp2_down), nrow(dn_loops) - sum(dn_loops$any_mecp2_down),
    sum(up_loops$any_mecp2_down), nrow(up_loops) - sum(up_loops$any_mecp2_down)),
  nrow = 2
))

cat(sprintf("  Fisher (gained loops vs MeCP2 Up): OR=%.2f, %s\n",
            fisher_up$estimate, fmt_p(fisher_up$p.value)))
cat(sprintf("  Fisher (lost loops vs MeCP2 Down): OR=%.2f, %s\n",
            fisher_dn$estimate, fmt_p(fisher_dn$p.value)))

p_31a <- ggplot(overlap_stats, aes(x = loop_direction, y = pct, fill = mecp2_direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, count)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = COLORS$mecp2[c("MeCP2 Up", "MeCP2 Down")],
                    name = "MeCP2 Binding Change") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "MeCP2 Peak Overlap at Hi-C Loop Anchors",
    subtitle = sprintf(
      "Fisher (gained/MeCP2 Up): OR=%.2f, %s | Fisher (lost/MeCP2 Down): OR=%.2f, %s",
      fisher_up$estimate, fmt_p(fisher_up$p.value),
      fisher_dn$estimate, fmt_p(fisher_dn$p.value)
    ),
    x = "Loop Direction (Mutant vs Control)", y = "% of Loops with Anchor MeCP2 Overlap"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_31a, file.path(SECTION_DIR, "31a_mecp2_peak_overlap_at_loop_anchors"),
                        width = 10, height = 8)

# Save per-loop MeCP2 overlap table
loop_mecp2_export <- loops %>%
  dplyr::select(loop_id, chr1, start1, end1, chr2, start2, end2,
                logFC, FDR, direction, loop_type_extended,
                anchor1_mecp2_up, anchor1_mecp2_down,
                anchor2_mecp2_up, anchor2_mecp2_down,
                any_mecp2_up, any_mecp2_down)

write.table(loop_mecp2_export, file.path(TABLES_DIR, "mecp2_loop_anchor_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: mecp2_loop_anchor_overlap.tsv\n")

# =============================================================================
# FIGURE 31b: MeCP2 Direction x Loop Direction Heatmap (Gene-Level)
# =============================================================================
cat("\nCreating Figure 31b: MeCP2 direction x loop direction heatmap...\n")

# Join anchor genes to MeCP2 gene-level data
anchor_mecp2 <- anchor_genes %>%
  inner_join(mecp2_gene, by = c("symbol" = "SYMBOL"))

cat(sprintf("  %d anchor-gene associations matched to MeCP2 data\n", nrow(anchor_mecp2)))
cat(sprintf("  %d unique genes at loop anchors with MeCP2 data\n",
            length(unique(anchor_mecp2$symbol))))

# Define MeCP2 status per gene (use nearest-to-TSS peak)
anchor_mecp2 <- anchor_mecp2 %>%
  mutate(
    mecp2_status = case_when(
      nearest_fdr < 0.05 & nearest_fold > 0 ~ "MeCP2 Up",
      nearest_fdr < 0.05 & nearest_fold < 0 ~ "MeCP2 Down",
      TRUE ~ "Not Significant"
    ),
    loop_dir_label = ifelse(loop_direction == "up_in_mutant", "Gained", "Lost")
  )

# Deduplicate: each gene counted once per loop direction
# (a gene can appear at multiple anchors — keep unique gene x loop_direction)
anchor_mecp2_unique <- anchor_mecp2 %>%
  distinct(symbol, loop_dir_label, .keep_all = TRUE)

# Build 2x2 table: loop direction x MeCP2 sig-up status
# Binary: MeCP2 sig up vs other (including sig down and NS)
heatmap_input <- anchor_mecp2_unique %>%
  mutate(mecp2_binary = ifelse(mecp2_status == "MeCP2 Up", "MeCP2 Sig Up", "Other"))

quad_table <- table(heatmap_input$loop_dir_label, heatmap_input$mecp2_binary)
cat("  Contingency table (loop direction x MeCP2 sig up):\n")
print(quad_table)

fisher_31b <- fisher.test(quad_table)
cat(sprintf("  Fisher's exact test: OR = %.2f, %s\n",
            fisher_31b$estimate, fmt_p(fisher_31b$p.value)))

# Compute expected and enrichment
total_n <- sum(quad_table)
row_totals <- rowSums(quad_table)
col_totals <- colSums(quad_table)
expected <- outer(row_totals, col_totals) / total_n
enrichment <- quad_table / expected

# Build heatmap data
heatmap_df <- as.data.frame(as.table(quad_table))
colnames(heatmap_df) <- c("loop_direction", "mecp2_status", "Observed")

expected_df <- as.data.frame(as.table(round(expected, 1)))
colnames(expected_df) <- c("loop_direction", "mecp2_status", "Expected")

enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
colnames(enrichment_df) <- c("loop_direction", "mecp2_status", "Enrichment")

heatmap_data <- heatmap_df %>%
  left_join(expected_df, by = c("loop_direction", "mecp2_status")) %>%
  left_join(enrichment_df, by = c("loop_direction", "mecp2_status")) %>%
  mutate(
    label = sprintf("Obs: %d\nExp: %.0f\nEnr: %.2f", Observed, Expected, Enrichment),
    loop_direction = factor(loop_direction, levels = c("Gained", "Lost")),
    mecp2_status = factor(mecp2_status, levels = c("MeCP2 Sig Up", "Other"))
  )

p_31b <- ggplot(heatmap_data, aes(x = mecp2_status, y = loop_direction, fill = Enrichment)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), size = 4.5, lineheight = 1.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                       limits = c(min(heatmap_data$Enrichment) * 0.9,
                                  max(heatmap_data$Enrichment) * 1.1)) +
  labs(
    title = "MeCP2 Significant Gain x Loop Direction (Gene-Level)",
    subtitle = sprintf("Fisher's exact test: OR = %.2f, %s | GREAT-style anchor-gene mapping",
                       fisher_31b$estimate, fmt_p(fisher_31b$p.value)),
    x = "MeCP2 Binding Status", y = "Loop Direction"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_multiformat_ggplot(p_31b, file.path(SECTION_DIR, "31b_mecp2_loop_direction_fisher_heatmap"),
                        width = 9, height = 7)

# Save Fisher's test results
fisher_31b_export <- heatmap_data %>%
  dplyr::select(loop_direction, mecp2_status, Observed, Expected, Enrichment) %>%
  mutate(fisher_OR = fisher_31b$estimate,
         fisher_pvalue = fisher_31b$p.value)

write.table(fisher_31b_export, file.path(TABLES_DIR, "mecp2_loop_direction_fisher.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: mecp2_loop_direction_fisher.tsv\n")

# =============================================================================
# FIGURE 31c: MeCP2 Fold Change by Loop Direction (Violin/Boxplot)
# =============================================================================
cat("\nCreating Figure 31c: MeCP2 fold change by loop direction...\n")

# Get unique genes per loop direction with MeCP2 data
gained_genes <- anchor_mecp2_unique %>%
  dplyr::filter(loop_dir_label == "Gained") %>%
  dplyr::select(symbol, nearest_fold) %>%
  mutate(group = "Gained\n(Loop Up)")

lost_genes <- anchor_mecp2_unique %>%
  dplyr::filter(loop_dir_label == "Lost") %>%
  dplyr::select(symbol, nearest_fold) %>%
  mutate(group = "Lost\n(Loop Down)")

# Background: genes with MeCP2 data NOT at any loop anchor
anchor_gene_symbols <- unique(anchor_genes$symbol)
background_genes <- mecp2_gene %>%
  dplyr::filter(!SYMBOL %in% anchor_gene_symbols) %>%
  dplyr::select(symbol = SYMBOL, nearest_fold) %>%
  mutate(group = "Background\n(No Loop)")

violin_data <- rbind(gained_genes, lost_genes, background_genes)
violin_data$group <- factor(violin_data$group,
                             levels = c("Gained\n(Loop Up)", "Lost\n(Loop Down)",
                                        "Background\n(No Loop)"))

# Wilcoxon tests
wilcox_gl <- wilcox.test(gained_genes$nearest_fold, lost_genes$nearest_fold)
wilcox_gb <- wilcox.test(gained_genes$nearest_fold, background_genes$nearest_fold)
wilcox_lb <- wilcox.test(lost_genes$nearest_fold, background_genes$nearest_fold)

cat(sprintf("  Median MeCP2 fold — Gained: %.4f, Lost: %.4f, Background: %.4f\n",
            median(gained_genes$nearest_fold, na.rm = TRUE),
            median(lost_genes$nearest_fold, na.rm = TRUE),
            median(background_genes$nearest_fold, na.rm = TRUE)))
cat(sprintf("  Wilcoxon Gained vs Lost: %s\n", fmt_p(wilcox_gl$p.value)))
cat(sprintf("  Wilcoxon Gained vs Background: %s\n", fmt_p(wilcox_gb$p.value)))
cat(sprintf("  Wilcoxon Lost vs Background: %s\n", fmt_p(wilcox_lb$p.value)))

# Group colors
violin_colors <- c("Gained\n(Loop Up)" = "#4575b4",
                    "Lost\n(Loop Down)" = "#d73027",
                    "Background\n(No Loop)" = "grey70")

# Build significance annotation text
sig_text <- sprintf("Gained vs Lost: %s\nGained vs Bg: %s\nLost vs Bg: %s",
                    fmt_p(wilcox_gl$p.value),
                    fmt_p(wilcox_gb$p.value),
                    fmt_p(wilcox_lb$p.value))

max_fold_val <- quantile(abs(violin_data$nearest_fold), 0.99, na.rm = TRUE)
y_upper <- max_fold_val * 1.2

p_31c <- ggplot(violin_data, aes(x = group, y = nearest_fold, fill = group)) +
  geom_violin(alpha = 0.5, scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = violin_colors) +
  coord_cartesian(ylim = c(-y_upper, y_upper)) +
  annotate("text", x = 2, y = y_upper * 0.9,
           label = sig_text, size = 3, hjust = 0.5, lineheight = 1.2) +
  labs(
    title = "MeCP2 Fold Change at Loop Anchor Genes",
    subtitle = sprintf("Gene-level MeCP2 binding change by loop direction | n: Gained=%d, Lost=%d, Bg=%d",
                       nrow(gained_genes), nrow(lost_genes), nrow(background_genes)),
    x = "", y = "MeCP2 log2 Fold Change\n(Mutant/Control)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_31c, file.path(SECTION_DIR, "31c_mecp2_fold_by_loop_direction"),
                        width = 9, height = 8)

# Save gene-level table
fold_export <- violin_data %>%
  dplyr::select(gene = symbol, mecp2_nearest_fold = nearest_fold, loop_group = group)

write.table(fold_export, file.path(TABLES_DIR, "mecp2_fold_at_loop_anchor_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: mecp2_fold_at_loop_anchor_genes.tsv\n")

# =============================================================================
# FIGURE 31d: Loop logFC vs MeCP2 Fold Change Scatter
# =============================================================================
cat("\nCreating Figure 31d: Loop logFC vs MeCP2 fold change scatter...\n")

# For each loop, compute mean MeCP2 fold across anchor genes
loop_mecp2_summary <- anchor_mecp2 %>%
  group_by(anchor_id) %>%
  summarise(
    mean_mecp2_fold = mean(nearest_fold, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  )

# Extract loop_id from anchor_id (e.g., "loop_8162_A1" -> "loop_8162")
loop_mecp2_summary <- loop_mecp2_summary %>%
  mutate(loop_id = sub("_A[12]$", "", anchor_id))

# Average across both anchors per loop
loop_mecp2_per_loop <- loop_mecp2_summary %>%
  group_by(loop_id) %>%
  summarise(
    mean_mecp2_fold = mean(mean_mecp2_fold, na.rm = TRUE),
    total_anchor_genes = sum(n_genes),
    .groups = "drop"
  )

# Join to loop data
scatter_loop <- loops %>%
  inner_join(loop_mecp2_per_loop, by = "loop_id") %>%
  mutate(loop_dir_label = ifelse(direction == "up_in_mutant", "Gained", "Lost"))

cat(sprintf("  %d loops with MeCP2 gene-level data\n", nrow(scatter_loop)))

# Spearman correlation
cor_loop_mecp2 <- cor.test(scatter_loop$logFC, scatter_loop$mean_mecp2_fold,
                            method = "spearman")
cat(sprintf("  Spearman rho = %.3f, %s\n",
            cor_loop_mecp2$estimate, fmt_p(cor_loop_mecp2$p.value)))

# Quadrant counts
q1 <- sum(scatter_loop$logFC > 0 & scatter_loop$mean_mecp2_fold > 0)
q2 <- sum(scatter_loop$logFC < 0 & scatter_loop$mean_mecp2_fold > 0)
q3 <- sum(scatter_loop$logFC < 0 & scatter_loop$mean_mecp2_fold < 0)
q4 <- sum(scatter_loop$logFC > 0 & scatter_loop$mean_mecp2_fold < 0)

# Get top genes for labeling (highest |mean_mecp2_fold| per loop)
top_label_loops <- scatter_loop %>%
  arrange(desc(abs(mean_mecp2_fold))) %>%
  head(15)

# Find representative gene name for each labeled loop
loop_gene_labels <- anchor_mecp2 %>%
  mutate(loop_id = sub("_A[12]$", "", anchor_id)) %>%
  dplyr::filter(loop_id %in% top_label_loops$loop_id) %>%
  group_by(loop_id) %>%
  slice_min(order_by = min_distance_tss, n = 1) %>%
  ungroup() %>%
  dplyr::select(loop_id, label_gene = symbol)

scatter_loop <- scatter_loop %>%
  left_join(loop_gene_labels, by = "loop_id") %>%
  mutate(label_gene = ifelse(is.na(label_gene), "", label_gene))

p_31d <- ggplot(scatter_loop, aes(x = logFC, y = mean_mecp2_fold, color = loop_dir_label)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE, alpha = 0.12) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text_repel(aes(label = label_gene), size = 2.8, max.overlaps = 12,
                  fontface = "italic", color = "grey20",
                  segment.color = "grey60", segment.size = 0.3) +
  scale_color_manual(values = c("Gained" = "#4575b4", "Lost" = "#d73027"),
                     name = "Loop Direction") +
  # Quadrant annotations
  annotate("text",
           x = max(scatter_loop$logFC, na.rm = TRUE) * 0.7,
           y = max(scatter_loop$mean_mecp2_fold, na.rm = TRUE) * 0.9,
           label = sprintf("Loop\u2191 MeCP2\u2191\nn=%d", q1),
           size = 3, color = "#4575b4", fontface = "bold") +
  annotate("text",
           x = min(scatter_loop$logFC, na.rm = TRUE) * 0.7,
           y = min(scatter_loop$mean_mecp2_fold, na.rm = TRUE) * 0.9,
           label = sprintf("Loop\u2193 MeCP2\u2193\nn=%d", q3),
           size = 3, color = "#d73027", fontface = "bold") +
  annotate("text",
           x = max(scatter_loop$logFC, na.rm = TRUE) * 0.7,
           y = min(scatter_loop$mean_mecp2_fold, na.rm = TRUE) * 0.9,
           label = sprintf("Loop\u2191 MeCP2\u2193\nn=%d", q4),
           size = 3, color = "grey50") +
  annotate("text",
           x = min(scatter_loop$logFC, na.rm = TRUE) * 0.7,
           y = max(scatter_loop$mean_mecp2_fold, na.rm = TRUE) * 0.9,
           label = sprintf("Loop\u2193 MeCP2\u2191\nn=%d", q2),
           size = 3, color = "grey50") +
  labs(
    title = "Loop Strength Change vs MeCP2 Binding Change",
    subtitle = sprintf("Spearman \u03C1 = %.3f, %s | Gene-level MeCP2 fold averaged across anchor genes",
                       cor_loop_mecp2$estimate, fmt_p(cor_loop_mecp2$p.value)),
    x = "Loop log2 Fold Change (Mutant/Control)",
    y = "Mean MeCP2 log2 Fold Change\n(Anchor Genes)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_31d, file.path(SECTION_DIR, "31d_loop_logfc_vs_mecp2_fold_scatter"),
                        width = 11, height = 9)

# Save scatter data table
scatter_export <- scatter_loop %>%
  dplyr::select(loop_id, chr1, start1, end1, chr2, start2, end2,
                logFC, FDR, direction, loop_type_extended,
                mean_mecp2_fold, total_anchor_genes, label_gene)

write.table(scatter_export, file.path(TABLES_DIR, "mecp2_loop_gene_level_scatter_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: mecp2_loop_gene_level_scatter_data.tsv\n")

# =============================================================================
# SECTION 31 SUMMARY
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("SECTION 31 SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("Total differential loops analyzed: %d\n", nrow(loops)))
cat(sprintf("Loops with MeCP2 anchor overlap: %d up, %d down\n",
            sum(loops$any_mecp2_up), sum(loops$any_mecp2_down)))
cat(sprintf("Anchor genes matched to MeCP2: %d unique genes\n",
            length(unique(anchor_mecp2$symbol))))
cat(sprintf("\n31a - Direct overlap:\n"))
cat(sprintf("  Gained loops/MeCP2 Up: Fisher OR=%.2f, %s\n",
            fisher_up$estimate, fmt_p(fisher_up$p.value)))
cat(sprintf("  Lost loops/MeCP2 Down: Fisher OR=%.2f, %s\n",
            fisher_dn$estimate, fmt_p(fisher_dn$p.value)))
cat(sprintf("\n31b - Gene-level direction association:\n"))
cat(sprintf("  Fisher OR=%.2f, %s\n", fisher_31b$estimate, fmt_p(fisher_31b$p.value)))
cat(sprintf("\n31c - MeCP2 fold by loop direction:\n"))
cat(sprintf("  Gained vs Lost: %s\n", fmt_p(wilcox_gl$p.value)))
cat(sprintf("  Gained vs Background: %s\n", fmt_p(wilcox_gb$p.value)))
cat(sprintf("  Lost vs Background: %s\n", fmt_p(wilcox_lb$p.value)))
cat(sprintf("\n31d - Loop logFC vs MeCP2 fold:\n"))
cat(sprintf("  Spearman rho=%.3f, %s\n",
            cor_loop_mecp2$estimate, fmt_p(cor_loop_mecp2$p.value)))

cat("\nSection 31 complete.\n\n")
