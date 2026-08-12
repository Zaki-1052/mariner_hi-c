# biomodal/downstream/scripts/viz_sections/section_74_geneset_overlap_methylation.R
# Section 74: Gene Set Overlap Counts & Neuronal Methylation Levels
#
# Creates a detailed 3-way Venn (Neuronal, Coordinated, MeCP2-Up) with pairwise
# Fisher enrichment tests, plus methylation level plots for the neuronal gene
# subset showing total (mC+hmC), 5mC, and 5hmC in ctrl vs mut.
#
# Panels:
#   74a: 3-way Venn diagram with gene counts + pairwise Fisher ORs
#   74b: Total methylation (mC+hmC) ctrl vs mut at neuronal genes
#   74c: 5mC levels ctrl vs mut at neuronal genes
#   74d: 5hmC levels ctrl vs mut at neuronal genes
#   74e: Composite
#
# Input:
#   tables/59_quadrant_master.tsv, tables/mecp2_coordinated_genes.tsv,
#   tables/72_neuronal_gene_set_go_derived.tsv
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_74_geneset_overlap_methylation.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(patchwork)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 74: GENE SET OVERLAP & NEURONAL METHYLATION LEVELS\n")
cat("  Detailed overlap counts + methylation characterization\n")
cat("================================================================================\n\n")

SEC74_DIR <- file.path(OUTPUT_DIR, "74_geneset_overlap_methylation")
dir.create(SEC74_DIR, recursive = TRUE, showWarnings = FALSE)

QUADRANT_PATH    <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
COORDINATED_PATH <- file.path(TABLES_DIR, "mecp2_coordinated_genes.tsv")
NEURONAL_PATH    <- file.path(TABLES_DIR, "72_neuronal_gene_set_go_derived.tsv")

stopifnot(
  "59_quadrant_master.tsv not found" = file.exists(QUADRANT_PATH),
  "mecp2_coordinated_genes.tsv not found" = file.exists(COORDINATED_PATH),
  "72_neuronal_gene_set_go_derived.tsv not found (run section 72 first)" = file.exists(NEURONAL_PATH)
)

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC74_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n\n")

qm <- read.table(QUADRANT_PATH, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Quadrant master:    %d genes\n", nrow(qm)))

coord <- read.table(COORDINATED_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Coordinated genes:  %d genes\n", nrow(coord)))

neuronal_genes <- read.table(NEURONAL_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")$gene
cat(sprintf("  Neuronal genes:     %d genes (GO-derived, section 72)\n",
            length(neuronal_genes)))

# --- Define gene sets within the quadrant master universe ---

universe <- qm$gene

mecp2_up_genes <- qm$gene[!is.na(qm$mecp2_nearest_fdr) &
                           qm$mecp2_nearest_fdr < Q_THRESHOLD &
                           !is.na(qm$mecp2_mean_fold) &
                           qm$mecp2_mean_fold > 0]
coordinated_in_univ <- intersect(coord$gene, universe)
neuronal_in_univ    <- intersect(neuronal_genes, universe)

cat(sprintf("\n  Universe:           %d genes\n", length(universe)))
cat(sprintf("  Neuronal in univ:   %d (%.1f%%)\n",
            length(neuronal_in_univ),
            100 * length(neuronal_in_univ) / length(universe)))
cat(sprintf("  Coordinated in univ:%d (%.1f%%)\n",
            length(coordinated_in_univ),
            100 * length(coordinated_in_univ) / length(universe)))
cat(sprintf("  MeCP2 Up:           %d (%.1f%%)\n",
            length(mecp2_up_genes),
            100 * length(mecp2_up_genes) / length(universe)))

# --- Flag genes ---

qm$is_neuronal    <- qm$gene %in% neuronal_in_univ
qm$is_coordinated <- qm$gene %in% coordinated_in_univ
qm$is_mecp2_up    <- qm$gene %in% mecp2_up_genes

# =============================================================================
# COMPUTE EXCLUSIVE PARTITION COUNTS
# =============================================================================

cat("\n--- Computing exclusive overlap partitions ---\n\n")

n <- nrow(qm)
N <- qm$is_neuronal
C <- qm$is_coordinated
M <- qm$is_mecp2_up

partitions <- data.frame(
  partition = c("Neuronal only", "Coordinated only", "MeCP2-Up only",
                "Neuronal + Coordinated", "Neuronal + MeCP2-Up",
                "Coordinated + MeCP2-Up", "All three", "None"),
  count = c(
    sum(N & !C & !M),
    sum(!N & C & !M),
    sum(!N & !C & M),
    sum(N & C & !M),
    sum(N & !C & M),
    sum(!N & C & M),
    sum(N & C & M),
    sum(!N & !C & !M)
  ),
  stringsAsFactors = FALSE
)
partitions$pct <- round(100 * partitions$count / n, 2)

for (i in seq_len(nrow(partitions))) {
  cat(sprintf("  %-26s %6d  (%.2f%%)\n",
              partitions$partition[i], partitions$count[i], partitions$pct[i]))
}
cat(sprintf("  %-26s %6d\n", "TOTAL", sum(partitions$count)))

write.table(partitions,
            file.path(TABLES_DIR, "74_geneset_overlap_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 74_geneset_overlap_counts.tsv\n")

# --- Pairwise Fisher enrichment tests ---

cat("\n--- Pairwise Fisher tests (enrichment of overlap) ---\n\n")

fisher_2x2 <- function(set_a, set_b, n_total, name) {
  a_and_b   <- sum(set_a & set_b)
  a_not_b   <- sum(set_a & !set_b)
  b_not_a   <- sum(!set_a & set_b)
  neither   <- sum(!set_a & !set_b)
  ft <- fisher.test(matrix(c(a_and_b, b_not_a, a_not_b, neither), nrow = 2))
  data.frame(
    comparison = name,
    overlap = a_and_b,
    set_a_size = sum(set_a),
    set_b_size = sum(set_b),
    odds_ratio = as.numeric(ft$estimate),
    ci_lo = ft$conf.int[1],
    ci_hi = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

fisher_results <- rbind(
  fisher_2x2(N, C, n, "Neuronal × Coordinated"),
  fisher_2x2(N, M, n, "Neuronal × MeCP2-Up"),
  fisher_2x2(C, M, n, "Coordinated × MeCP2-Up")
)
fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "BH")

for (i in seq_len(nrow(fisher_results))) {
  r <- fisher_results[i, ]
  cat(sprintf("  %s: overlap=%d, OR=%.2f [%.2f, %.2f], %s (adj %s)\n",
              r$comparison, r$overlap, r$odds_ratio, r$ci_lo, r$ci_hi,
              fmt_p(r$p_value), fmt_p(r$p_adj)))
}

write.table(fisher_results,
            file.path(TABLES_DIR, "74_pairwise_fisher.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 74_pairwise_fisher.tsv\n")

# =============================================================================
# 74a: 3-WAY VENN DIAGRAM
# =============================================================================

cat("\n--- 74a: Venn diagram ---\n")

venn_list <- list(
  "Neuronal" = neuronal_in_univ,
  "Coordinated\n(mC↑/hmC↓)" = coordinated_in_univ,
  "MeCP2 Up" = mecp2_up_genes
)

fisher_label <- paste0(
  sprintf("Neuronal × Coord: OR=%.1f, %s", fisher_results$odds_ratio[1], fmt_p(fisher_results$p_adj[1])),
  "\n",
  sprintf("Neuronal × MeCP2: OR=%.1f, %s", fisher_results$odds_ratio[2], fmt_p(fisher_results$p_adj[2])),
  "\n",
  sprintf("Coord × MeCP2: OR=%.1f, %s", fisher_results$odds_ratio[3], fmt_p(fisher_results$p_adj[3]))
)

p_74a <- ggVennDiagram(venn_list, label = "count", label_alpha = 0,
                        edge_size = 0.8) +
  scale_fill_gradient(low = "white", high = "#756BB1", guide = "none") +
  labs(
    title = "Gene set overlap: Neuronal, Coordinated, MeCP2-Up",
    subtitle = sprintf("Universe: %s genes from quadrant master", format(length(universe), big.mark = ",")),
    caption = fisher_label
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.caption = element_text(hjust = 0, size = 9, color = "grey30")
  )

save_plot(p_74a, "74a_geneset_venn", w = 9, h = 8)

# =============================================================================
# 74b-d: METHYLATION LEVELS AT NEURONAL GENES
# =============================================================================

cat("\n--- 74b-d: Methylation levels at neuronal genes ---\n")

neur_df <- qm[qm$is_neuronal, ]

has_mc  <- !is.na(neur_df$mc_ctrl) & !is.na(neur_df$mc_mut)
has_hmc <- !is.na(neur_df$hmc_ctrl) & !is.na(neur_df$hmc_mut)
has_both <- has_mc & has_hmc
neur_meth <- neur_df[has_both, ]

neur_meth$total_ctrl <- neur_meth$mc_ctrl + neur_meth$hmc_ctrl
neur_meth$total_mut  <- neur_meth$mc_mut + neur_meth$hmc_mut

cat(sprintf("  Neuronal genes with methylation data: %d / %d\n",
            nrow(neur_meth), nrow(neur_df)))

# --- Save neuronal methylation table ---

write.table(
  neur_meth[, c("gene", "chr", "start", "end",
                "mc_ctrl", "mc_mut", "hmc_ctrl", "hmc_mut",
                "total_ctrl", "total_mut",
                "is_coordinated", "is_mecp2_up")],
  file.path(TABLES_DIR, "74_neuronal_methylation_levels.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("  Saved: 74_neuronal_methylation_levels.tsv\n")

# --- Helper: paired methylation violin ---

make_meth_violin <- function(ctrl_vals, mut_vals, mark_label, color_ctrl, color_mut) {
  long <- data.frame(
    condition = rep(c("Control", "Mutant"), each = length(ctrl_vals)),
    value = c(ctrl_vals, mut_vals),
    stringsAsFactors = FALSE
  )
  long$condition <- factor(long$condition, levels = c("Control", "Mutant"))

  wt <- wilcox.test(mut_vals, ctrl_vals, paired = TRUE)
  delta <- median(mut_vals - ctrl_vals)
  med_c <- median(ctrl_vals)
  med_m <- median(mut_vals)

  cat(sprintf("  %s: ctrl median=%.5f, mut median=%.5f, delta=%.5f, %s\n",
              mark_label, med_c, med_m, delta, fmt_p(wt$p.value)))

  ggplot(long, aes(x = condition, y = value, fill = condition)) +
    geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = c("Control" = unname(color_ctrl),
                                  "Mutant" = unname(color_mut))) +
    labs(
      title = mark_label,
      subtitle = sprintf("n=%s | Δ median=%.5f | %s",
                         format(length(ctrl_vals), big.mark = ","),
                         delta, fmt_p(wt$p.value)),
      x = NULL, y = "Methylation fraction"
    ) +
    theme_biomodal() +
    theme(plot.subtitle = element_text(size = 9))
}

# 74b: Total methylation
p_74b <- make_meth_violin(
  neur_meth$total_ctrl, neur_meth$total_mut,
  "Total methylation (mC + hmC)\nat neuronal genes",
  COLORS$condition["Control"], COLORS$condition["Mutant"]
)
save_plot(p_74b, "74b_neuronal_total_methylation", w = 6, h = 7)

# 74c: 5mC
p_74c <- make_meth_violin(
  neur_meth$mc_ctrl, neur_meth$mc_mut,
  "5mC levels\nat neuronal genes",
  COLORS$condition["Control"], COLORS$condition["Mutant"]
)
save_plot(p_74c, "74c_neuronal_5mc_levels", w = 6, h = 7)

# 74d: 5hmC
p_74d <- make_meth_violin(
  neur_meth$hmc_ctrl, neur_meth$hmc_mut,
  "5hmC levels\nat neuronal genes",
  COLORS$condition["Control"], COLORS$condition["Mutant"]
)
save_plot(p_74d, "74d_neuronal_5hmc_levels", w = 6, h = 7)

# =============================================================================
# 74e: COMPOSITE
# =============================================================================

cat("\n--- 74e: Composite ---\n")

p_74e <- (p_74a | (p_74b / p_74c / p_74d)) +
  plot_annotation(
    title = "Section 74: Gene Set Overlap & Neuronal Methylation Characterization",
    subtitle = sprintf("Neuronal: %s genes (GO-derived) | Coordinated: %s | MeCP2-Up: %d",
                        format(length(neuronal_in_univ), big.mark = ","),
                        format(length(coordinated_in_univ), big.mark = ","),
                        length(mecp2_up_genes)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

save_plot(p_74e, "74_composite", w = 18, h = 16)

cat("\n================================================================================\n")
cat("SECTION 74 COMPLETE\n")
cat("================================================================================\n")
