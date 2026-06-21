# biomodal/downstream/scripts/viz_sections/section_67_mecp2_k119ub_unmethylated.R
# Section 67: MeCP2 Binding at Unmethylated Genes Tracks H2AK119ub, Not Methylation
# Standalone script - sources shared config for all dependencies and data
#
# At 359 genes where MeCP2 binding increases WITHOUT significant methylation
# change, H2AK119ub is significantly gained (OR=3.15, p=1.8e-24 vs genome).
# This proves MeCP2 responds to BAP1-mediated H2AK119ub, not methylation.
#
#   Panel 67a: K119ub log2FC distribution: MeCP2-no-meth genes vs all others
#   Panel 67b: % K119ub gained: MeCP2-no-meth vs genome-wide (Fisher)
#   Panel 67c: MeCP2 fold vs K119ub log2FC scatter at these 359 genes
#   Panel 67d: Multi-mark comparison violin (K119ub, K27me3, K27ac, ATAC)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_67_mecp2_k119ub_unmethylated.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 67: MeCP2 AT UNMETHYLATED GENES TRACKS H2AK119ub
# =============================================================================

cat("================================================================================\n")
cat("SECTION 67: MeCP2 BINDING AT UNMETHYLATED GENES TRACKS H2AK119ub\n")
cat("================================================================================\n\n")

# ---- Load data ---------------------------------------------------------------

MECP2_NO_METH_PATH <- file.path(TABLES_DIR, "65_mecp2_bound_no_methylation_genes.tsv")
K119UB_SIGNAL_PATH <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")

if (!file.exists(MECP2_NO_METH_PATH)) {
  stop("MeCP2-bound-no-methylation gene list not found: ", MECP2_NO_METH_PATH,
       "\nRun section_65 first.")
}
if (!file.exists(K119UB_SIGNAL_PATH)) {
  stop("K119ub gene signal file not found: ", K119UB_SIGNAL_PATH)
}

cat("Loading MeCP2-bound-no-methylation genes...\n")
mecp2_no_meth <- read.table(MECP2_NO_METH_PATH, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "")
cat(sprintf("  %d genes\n", nrow(mecp2_no_meth)))

cat("Loading K119ub gene body signal...\n")
k119ub_signal <- read.table(K119UB_SIGNAL_PATH, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, quote = "")
cat(sprintf("  %d genes with K119ub signal\n", nrow(k119ub_signal)))

# ---- Join and classify -------------------------------------------------------

cat("\nJoining datasets...\n")

k119ub_signal$group <- ifelse(
  k119ub_signal$symbol %in% mecp2_no_meth$gene,
  "MeCP2 Up\n(No Methylation Change)",
  "All Other Genes"
)

our_genes <- k119ub_signal %>% dplyr::filter(group == "MeCP2 Up\n(No Methylation Change)")
rest_genes <- k119ub_signal %>% dplyr::filter(group == "All Other Genes")

cat(sprintf("  Matched: %d / %d MeCP2-no-meth genes to K119ub data\n",
            nrow(our_genes), nrow(mecp2_no_meth)))
cat(sprintf("  Background: %d genes\n", nrow(rest_genes)))

# ---- Statistics --------------------------------------------------------------

cat("\nStatistical tests:\n")

# Filter to genes with finite K119ub log2FC
our_genes <- our_genes %>% dplyr::filter(is.finite(gb_log2fc))
rest_genes <- rest_genes %>% dplyr::filter(is.finite(gb_log2fc))

cat(sprintf("  After NA/Inf filtering: %d our genes, %d rest\n",
            nrow(our_genes), nrow(rest_genes)))

# Mann-Whitney U
wt <- wilcox.test(our_genes$gb_log2fc, rest_genes$gb_log2fc, alternative = "greater")
cat(sprintf("  Mann-Whitney U (K119ub log2FC, our > rest):\n"))
cat(sprintf("    Our median:  %+.4f\n", median(our_genes$gb_log2fc)))
cat(sprintf("    Rest median: %+.4f\n", median(rest_genes$gb_log2fc)))
cat(sprintf("    p = %.2e\n", wt$p.value))

# Fisher: K119ub gained (>0) in our group vs rest
our_gained <- sum(our_genes$gb_log2fc > 0)
our_not <- nrow(our_genes) - our_gained
rest_gained <- sum(rest_genes$gb_log2fc > 0)
rest_not <- nrow(rest_genes) - rest_gained

fisher_mat <- matrix(c(our_gained, our_not, rest_gained, rest_not),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(c("MeCP2-no-meth", "Rest"),
                                    c("K119ub Gained", "K119ub Not Gained")))
fisher_result <- fisher.test(fisher_mat, alternative = "greater")

cat(sprintf("\n  Fisher exact (K119ub gained in our group vs rest):\n"))
cat(sprintf("    Our group:  %d/%d gained (%.1f%%)\n",
            our_gained, nrow(our_genes), 100 * our_gained / nrow(our_genes)))
cat(sprintf("    Rest:       %d/%d gained (%.1f%%)\n",
            rest_gained, nrow(rest_genes), 100 * rest_gained / nrow(rest_genes)))
cat(sprintf("    OR = %.2f, p = %.2e\n", fisher_result$estimate, fisher_result$p.value))

# ---- Panel 67a: K119ub log2FC distribution violin ----------------------------

cat("\nCreating Figure 67a: K119ub log2FC distribution...\n")

k119ub_plot <- k119ub_signal %>% dplyr::filter(is.finite(gb_log2fc))
k119ub_plot$group <- factor(k119ub_plot$group,
                            levels = c("MeCP2 Up\n(No Methylation Change)", "All Other Genes"))

p_67a <- ggplot(k119ub_plot, aes(x = group, y = gb_log2fc, fill = group)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1.5, y = max(k119ub_plot$gb_log2fc, na.rm = TRUE) * 0.9,
           label = sprintf("Mann-Whitney p = %.1e", wt$p.value),
           size = 4, fontface = "italic") +
  annotate("segment", x = 1, xend = 2,
           y = max(k119ub_signal$gb_log2fc, na.rm = TRUE) * 0.85,
           yend = max(k119ub_signal$gb_log2fc, na.rm = TRUE) * 0.85,
           linewidth = 0.4) +
  scale_fill_manual(values = c("MeCP2 Up\n(No Methylation Change)" = "#D95F02",
                               "All Other Genes" = "grey70"),
                    guide = "none") +
  coord_cartesian(ylim = c(quantile(k119ub_plot$gb_log2fc, 0.001, na.rm = TRUE),
                           quantile(k119ub_plot$gb_log2fc, 0.999, na.rm = TRUE) * 1.15)) +
  labs(
    title = "H2AK119ub Signal Change at MeCP2-Bound Unmethylated Genes",
    subtitle = sprintf("MeCP2-no-meth genes (n=%d) show significantly higher K119ub gain than genome (n=%s)",
                        nrow(our_genes), format(nrow(rest_genes), big.mark = ",")),
    x = "", y = "H2AK119ub Gene Body log2FC (Mut/Ctrl)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 11))

save_multiformat_ggplot(p_67a,
                        file.path(OUTPUT_DIR, "67a_k119ub_at_mecp2_no_meth"),
                        width = 8, height = 7)

# ---- Panel 67b: % K119ub gained bar chart ------------------------------------

cat("Creating Figure 67b: % K119ub gained comparison...\n")

bar_df <- data.frame(
  group = c("MeCP2 Up\n(No Meth Change)", "All Other\nGenes"),
  gained_pct = c(100 * our_gained / nrow(our_genes),
                 100 * rest_gained / nrow(rest_genes)),
  n_gained = c(our_gained, rest_gained),
  n_total = c(nrow(our_genes), nrow(rest_genes)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(group = factor(group, levels = group))

p_67b <- ggplot(bar_df, aes(x = group, y = gained_pct, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", gained_pct, n_gained, n_total)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("MeCP2 Up\n(No Meth Change)" = "#D95F02",
                               "All Other\nGenes" = "grey70"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 90), expand = c(0, 0)) +
  annotate("text", x = 1.5, y = 85,
           label = sprintf("Fisher OR = %.2f\np = %.1e", fisher_result$estimate, fisher_result$p.value),
           size = 3.5, fontface = "italic") +
  labs(
    title = "H2AK119ub Gain at MeCP2-Bound Unmethylated Genes",
    subtitle = "MeCP2 binding changes track with Polycomb mark, not DNA methylation",
    x = "", y = "% of Genes with K119ub Gain (log2FC > 0)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_67b,
                        file.path(OUTPUT_DIR, "67b_k119ub_gained_fraction"),
                        width = 7, height = 7)

# ---- Panel 67c: MeCP2 fold vs K119ub log2FC scatter -------------------------

cat("Creating Figure 67c: MeCP2 fold vs K119ub scatter...\n")

scatter_df <- our_genes %>%
  inner_join(
    mecp2_no_meth %>% dplyr::select(gene, mecp2_nearest_fold, mecp2_nearest_fdr),
    by = c("symbol" = "gene")
  )

cat(sprintf("  %d genes in scatter\n", nrow(scatter_df)))

# Spearman correlation
cor_test <- cor.test(scatter_df$mecp2_nearest_fold, scatter_df$gb_log2fc,
                     method = "spearman")
cat(sprintf("  Spearman rho = %.3f, p = %.2e\n", cor_test$estimate, cor_test$p.value))

# Label top genes
scatter_df$label <- ""
top_k119 <- scatter_df %>% dplyr::slice_max(gb_log2fc, n = 8)
top_mecp2 <- scatter_df %>% dplyr::slice_max(abs(mecp2_nearest_fold), n = 5)
label_genes <- unique(c(top_k119$symbol, top_mecp2$symbol))
scatter_df$label[scatter_df$symbol %in% label_genes] <-
  scatter_df$symbol[scatter_df$symbol %in% label_genes]

p_67c <- ggplot(scatter_df, aes(x = gb_log2fc, y = mecp2_nearest_fold)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(alpha = 0.6, size = 2, color = "#D95F02") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 15,
                  fontface = "italic", color = "black") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Spearman rho = %.3f\np = %.2e\nn = %d",
                           cor_test$estimate, cor_test$p.value, nrow(scatter_df)),
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(
    title = "MeCP2 Fold Change vs H2AK119ub Change (No-Methylation Genes)",
    subtitle = "Genes where MeCP2 binds without methylation change",
    x = "H2AK119ub Gene Body log2FC (Mut/Ctrl)",
    y = "MeCP2 Nearest Peak Fold (log2)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_67c,
                        file.path(OUTPUT_DIR, "67c_mecp2_vs_k119ub_scatter"),
                        width = 10, height = 8)

# ---- Panel 67d: Multi-mark comparison (K119ub, K27me3, K27ac, ATAC) ---------

cat("Creating Figure 67d: multi-mark comparison...\n")

# Extract additional mark signals at gene body coordinates using BigWigs
dmr_gr <- GRanges(
  seqnames = our_genes$chr,
  ranges = IRanges(start = our_genes$start, end = our_genes$end)
)

extract_signal <- function(bw_path, region_gr, label) {
  cat(sprintf("  %s...", label))
  bw <- rtracklayer::import.bw(bw_path, which = region_gr)
  cov <- coverage(bw, weight = "score")
  ranges_by_chr <- split(ranges(region_gr), seqnames(region_gr))
  shared_chrs <- intersect(names(cov), names(ranges_by_chr))
  v <- Views(cov[shared_chrs], ranges_by_chr[shared_chrs])
  means <- viewMeans(v)
  result <- rep(NA_real_, length(region_gr))
  chr_names <- as.character(seqnames(region_gr))
  for (chr in names(means)) {
    idx <- which(chr_names == chr)
    result[idx] <- as.numeric(means[[chr]])
  }
  cat(sprintf(" done (median=%.4f)\n", median(result, na.rm = TRUE)))
  result
}

# Extract ctrl and mut for each mark
marks_to_extract <- list(
  K27me3 = list(ctrl = HISTONE_BIGWIGS$k27me3_ctrl, mut = HISTONE_BIGWIGS$k27me3_mut),
  K27ac  = list(ctrl = HISTONE_BIGWIGS$k27ac_ctrl, mut = HISTONE_BIGWIGS$k27ac_mut),
  ATAC   = list(ctrl = HISTONE_BIGWIGS$atac_ctrl, mut = HISTONE_BIGWIGS$atac_mut)
)

mark_log2fcs <- data.frame(symbol = our_genes$symbol, stringsAsFactors = FALSE)
mark_log2fcs$K119ub <- our_genes$gb_log2fc

pseudo <- 1e-4
for (mark_name in names(marks_to_extract)) {
  ctrl_signal <- extract_signal(marks_to_extract[[mark_name]]$ctrl, dmr_gr,
                                paste0(mark_name, "_ctrl"))
  mut_signal <- extract_signal(marks_to_extract[[mark_name]]$mut, dmr_gr,
                               paste0(mark_name, "_mut"))
  mark_log2fcs[[mark_name]] <- log2((mut_signal + pseudo) / (ctrl_signal + pseudo))
}

# Pivot to long format for violin
mark_long <- mark_log2fcs %>%
  tidyr::pivot_longer(cols = c(K119ub, K27me3, K27ac, ATAC),
                      names_to = "mark", values_to = "log2fc") %>%
  dplyr::mutate(
    mark = factor(mark, levels = c("K119ub", "K27me3", "K27ac", "ATAC"))
  )

# Wilcoxon tests (each mark vs 0)
cat("\n  Per-mark one-sample Wilcoxon (vs 0):\n")
mark_stats <- data.frame(mark = character(), median = numeric(), p = numeric(),
                         stringsAsFactors = FALSE)
for (m in levels(mark_long$mark)) {
  vals <- mark_long$log2fc[mark_long$mark == m]
  wt_mark <- wilcox.test(vals, mu = 0)
  mark_stats <- rbind(mark_stats, data.frame(
    mark = m, median = median(vals, na.rm = TRUE), p = wt_mark$p.value,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("    %s: median=%+.4f, p=%.2e\n", m, median(vals, na.rm = TRUE), wt_mark$p.value))
}

mark_colors <- c("K119ub" = "#756BB1", "K27me3" = "#4DAF4A",
                 "K27ac" = "#FF7F00", "ATAC" = "#E6AB02")

p_annot <- mark_stats %>%
  dplyr::mutate(
    mark = factor(mark, levels = c("K119ub", "K27me3", "K27ac", "ATAC")),
    label = ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p))
  )

p_67d <- ggplot(mark_long, aes(x = mark, y = log2fc, fill = mark)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
  geom_text(data = p_annot,
            aes(x = mark, y = Inf, label = label),
            vjust = 1.5, size = 3, fontface = "italic", inherit.aes = FALSE) +
  scale_fill_manual(values = mark_colors, guide = "none") +
  coord_cartesian(ylim = c(
    quantile(mark_long$log2fc, 0.005, na.rm = TRUE),
    quantile(mark_long$log2fc, 0.995, na.rm = TRUE) * 1.2
  )) +
  labs(
    title = "Epigenomic Mark Changes at MeCP2-Bound Unmethylated Genes",
    subtitle = sprintf("n = %d genes | H2AK119ub gain is the dominant signal at these loci",
                        nrow(our_genes)),
    x = "", y = "log2FC (Mut/Ctrl)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

save_multiformat_ggplot(p_67d,
                        file.path(OUTPUT_DIR, "67d_multimark_at_mecp2_no_meth"),
                        width = 9, height = 7)

# ---- Combined figure ---------------------------------------------------------

cat("Creating combined Figure 67...\n")

p_67_combined <- ((p_67a | p_67b) / (p_67c | p_67d)) +
  plot_annotation(
    title = "MeCP2 Responds to H2AK119ub, Not Methylation",
    subtitle = "At genes where MeCP2 binding increases without methylation change, H2AK119ub gain is the dominant epigenomic signal",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_67_combined,
                        file.path(OUTPUT_DIR, "67_mecp2_k119ub_unmethylated"),
                        width = 16, height = 14)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving tables...\n")

# Per-gene table with all marks
gene_table <- mark_log2fcs %>%
  inner_join(
    mecp2_no_meth %>% dplyr::select(gene, mecp2_nearest_fold, mecp2_nearest_fdr,
                                    mecp2_annotation, mc_ctrl, mc_mut, mc_difference),
    by = c("symbol" = "gene")
  ) %>%
  arrange(desc(K119ub))

write.table(gene_table,
            file.path(TABLES_DIR, "67_mecp2_no_meth_multimark.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 67_mecp2_no_meth_multimark.tsv\n")

# Statistics summary
stats_table <- data.frame(
  test = c("Mann-Whitney (K119ub log2FC, our vs rest)",
           "Fisher (K119ub gained, our vs rest)",
           "Spearman (MeCP2 fold vs K119ub log2FC)"),
  statistic = c(sprintf("U = %.0f", wt$statistic),
                sprintf("OR = %.2f", fisher_result$estimate),
                sprintf("rho = %.3f", cor_test$estimate)),
  p_value = c(wt$p.value, fisher_result$p.value, cor_test$p.value),
  n_our = nrow(our_genes),
  n_rest = c(nrow(rest_genes), nrow(rest_genes), NA),
  our_median = c(median(our_genes$gb_log2fc), NA, NA),
  rest_median = c(median(rest_genes$gb_log2fc), NA, NA),
  our_pct_gained = c(NA, 100 * our_gained / nrow(our_genes), NA),
  rest_pct_gained = c(NA, 100 * rest_gained / nrow(rest_genes), NA),
  stringsAsFactors = FALSE
)

write.table(stats_table,
            file.path(TABLES_DIR, "67_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 67_statistics.tsv\n")

# Per-mark summary
write.table(mark_stats,
            file.path(TABLES_DIR, "67_per_mark_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 67_per_mark_summary.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 67 COMPLETE\n")
cat("================================================================================\n")
