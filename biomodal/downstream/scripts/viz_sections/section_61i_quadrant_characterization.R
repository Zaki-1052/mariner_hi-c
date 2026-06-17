# biomodal/downstream/scripts/viz_sections/section_61i_quadrant_characterization.R
# Section 61i: Characterize MeCP2-Up + K119ub-Up quadrant genes
#
# Chromatin state composition + neuronal vs non-neuronal comparison
# for the 72 top-right quadrant genes from section 61h.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_61i_quadrant_characterization.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(patchwork)
})

cat("================================================================================\n")
cat("SECTION 61i: QUADRANT GENE CHARACTERIZATION\n")
cat("  Chromatin states + neuronal vs non-neuronal within 72 top-right genes\n")
cat("================================================================================\n\n")

SEC61I_DIR <- file.path(OUTPUT_DIR, "61i_quadrant_characterization")
dir.create(SEC61I_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC61I_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# --- Load data ---

quad_path <- file.path(TABLES_DIR, "61h_mecp2_up_k119ub_up_genes.tsv")
neur_path <- file.path(TABLES_DIR, "61_neuronal_gene_set.tsv")
db_path   <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
qm_path   <- file.path(TABLES_DIR, "59_quadrant_master.tsv")

for (f in c(quad_path, neur_path, db_path, qm_path)) {
  if (!file.exists(f)) stop("Missing: ", f)
}

quad_genes <- read.table(quad_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "")
neur_genes <- read.table(neur_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "")
db <- read.table(db_path, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")
qm <- read.table(qm_path, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")

cat(sprintf("  Quadrant genes: %d\n", nrow(quad_genes)))
cat(sprintf("  Neuronal gene set: %d\n", nrow(neur_genes)))

# --- Merge to get all marks ---

quad_full <- merge(quad_genes, db[, c("gene", "atac_fold", "k27ac_fold",
                                       "k27me3_fold", "k119ub_fold",
                                       "total_ctrl", "total_mut", "delta_ratio",
                                       "dmr_status")],
                   by = "gene", all.x = TRUE)

quad_full$is_neuronal <- quad_full$gene %in% neur_genes$gene
quad_full$gene_class <- ifelse(quad_full$is_neuronal, "Neuronal", "Non-neuronal")

n_neur <- sum(quad_full$is_neuronal)
n_nonneur <- sum(!quad_full$is_neuronal)
cat(sprintf("\n  Neuronal:     %d / %d (%.0f%%)\n", n_neur, nrow(quad_full),
            100 * n_neur / nrow(quad_full)))
cat(sprintf("  Non-neuronal: %d / %d (%.0f%%)\n", n_nonneur, nrow(quad_full),
            100 * n_nonneur / nrow(quad_full)))

# Print neuronal genes
cat(sprintf("\n  Neuronal genes in quadrant: %s\n",
            paste(sort(quad_full$gene[quad_full$is_neuronal]), collapse = ", ")))

# =============================================================================
# CHROMATIN STATE DISTRIBUTION
# =============================================================================

cat("\n--- Chromatin state distribution ---\n")

state_ct <- table(quad_full$chromatin_state)
cat("  Quadrant gene chromatin states:\n")
for (s in names(sort(state_ct, decreasing = TRUE))) {
  cat(sprintf("    %-22s %d (%.0f%%)\n", s, state_ct[s],
              100 * state_ct[s] / nrow(quad_full)))
}

# Compare to genome-wide MeCP2-Up distribution
qm$mecp2_up <- !is.na(qm$mecp2_nearest_fdr) &
  qm$mecp2_nearest_fdr < Q_THRESHOLD & qm$mecp2_mean_fold > 0
genome_states <- table(qm$chromatin_state[qm$mecp2_up])

cat("\n  Genome-wide MeCP2-Up chromatin states (for comparison):\n")
for (s in names(sort(genome_states, decreasing = TRUE))) {
  cat(sprintf("    %-22s %d (%.0f%%)\n", s, genome_states[s],
              100 * genome_states[s] / sum(genome_states)))
}

# Chromatin state bar
state_df <- data.frame(
  state = names(state_ct),
  count = as.numeric(state_ct),
  stringsAsFactors = FALSE
)
state_df$frac <- state_df$count / sum(state_df$count)
state_df <- state_df[order(-state_df$count), ]
state_df$state <- factor(state_df$state, levels = state_df$state)

p_state <- ggplot(state_df, aes(x = state, y = count, fill = state)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", count, 100 * frac)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS[levels(state_df$state)],
                    na.value = "grey70") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Chromatin states of MeCP2-Up + K119ub-Up genes (n=72)",
       x = NULL, y = "Count") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# =============================================================================
# NEURONAL vs NON-NEURONAL COMPARISON
# =============================================================================

cat("\n--- Neuronal vs Non-neuronal comparison ---\n")

metrics <- c("mecp2_mean_fold", "gb_log2fc", "mc_diff", "hmc_diff",
             "k27me3_fold", "k27ac_fold", "k119ub_fold", "atac_fold", "delta_ratio")
metric_labels <- c("MeCP2 fold", "K119ub log2FC", "5mC diff", "5hmC diff",
                   "K27me3 fold", "K27ac fold", "K119ub fold (DiffBind)",
                   "ATAC fold", "Delta ratio")

comp_results <- data.frame()
for (i in seq_along(metrics)) {
  col <- metrics[i]
  if (!col %in% colnames(quad_full)) next

  neur_vals <- quad_full[[col]][quad_full$is_neuronal]
  nonneur_vals <- quad_full[[col]][!quad_full$is_neuronal]

  neur_vals <- neur_vals[!is.na(neur_vals)]
  nonneur_vals <- nonneur_vals[!is.na(nonneur_vals)]

  if (length(neur_vals) < 3 || length(nonneur_vals) < 3) next

  wt <- tryCatch(wilcox.test(neur_vals, nonneur_vals), error = function(e) NULL)

  comp_results <- rbind(comp_results, data.frame(
    metric = metric_labels[i],
    n_neuronal = length(neur_vals),
    n_nonneuronal = length(nonneur_vals),
    mean_neuronal = mean(neur_vals),
    mean_nonneuronal = mean(nonneur_vals),
    median_neuronal = median(neur_vals),
    median_nonneuronal = median(nonneur_vals),
    wilcox_p = if (!is.null(wt)) wt$p.value else NA,
    stringsAsFactors = FALSE
  ))
}

if (nrow(comp_results) > 0) {
  comp_results$wilcox_q <- p.adjust(comp_results$wilcox_p, method = "BH")
  comp_results$sig <- sapply(comp_results$wilcox_q, function(q) {
    if (is.na(q)) return("")
    if (q < 0.001) "***" else if (q < 0.01) "**" else if (q < 0.05) "*" else "ns"
  })

  cat("\n  Metric comparison (Wilcoxon, BH-corrected):\n")
  for (i in seq_len(nrow(comp_results))) {
    r <- comp_results[i, ]
    cat(sprintf("    %-22s Neur=%.4f  NonNeur=%.4f  q=%.3f %s\n",
                r$metric, r$mean_neuronal, r$mean_nonneuronal,
                r$wilcox_q, r$sig))
  }
}

# Violin comparison for key metrics
plot_metrics <- c("mecp2_mean_fold", "gb_log2fc", "mc_diff", "hmc_diff")
plot_labels <- c("MeCP2 fold change", "K119ub log2FC (BigWig)", "5mC difference", "5hmC difference")

violin_panels <- list()
for (i in seq_along(plot_metrics)) {
  col <- plot_metrics[i]
  if (!col %in% colnames(quad_full)) next

  vdf <- quad_full[!is.na(quad_full[[col]]), c("gene", "gene_class", col)]
  colnames(vdf)[3] <- "value"

  p <- ggplot(vdf, aes(x = gene_class, y = value, fill = gene_class)) +
    geom_boxplot(alpha = 0.6, outlier.size = 1) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = c("Neuronal" = "#D73027", "Non-neuronal" = "#4575B4")) +
    labs(title = plot_labels[i], x = NULL, y = NULL) +
    theme_biomodal(base_size = 10) +
    theme(legend.position = "none")

  violin_panels[[length(violin_panels) + 1]] <- p
}

p_comp <- wrap_plots(violin_panels, ncol = 2) +
  plot_annotation(
    title = "Neuronal vs Non-neuronal within MeCP2-Up + K119ub-Up quadrant (n=72)",
    subtitle = sprintf("Neuronal: %d genes | Non-neuronal: %d genes",
                       n_neur, n_nonneur)
  )

# =============================================================================
# PER-GENE STRIP CHART
# =============================================================================

cat("\n--- Per-gene multi-mark strip ---\n")

strip_metrics <- c("mecp2_mean_fold", "gb_log2fc", "mc_diff", "hmc_diff",
                   "k27me3_fold", "k27ac_fold")
strip_labels <- c("MeCP2", "K119ub", "5mC", "5hmC", "K27me3", "K27ac")

strip_long <- do.call(rbind, lapply(seq_along(strip_metrics), function(i) {
  col <- strip_metrics[i]
  if (!col %in% colnames(quad_full)) return(NULL)
  data.frame(
    gene = quad_full$gene,
    gene_class = quad_full$gene_class,
    mark = strip_labels[i],
    value = quad_full[[col]],
    stringsAsFactors = FALSE
  )
}))
strip_long <- strip_long[!is.na(strip_long$value), ]
strip_long$mark <- factor(strip_long$mark, levels = strip_labels)

# Order genes by MeCP2 fold
gene_order <- quad_full$gene[order(quad_full$mecp2_mean_fold, decreasing = TRUE)]
strip_long$gene <- factor(strip_long$gene, levels = gene_order)

p_strip <- ggplot(strip_long, aes(x = gene, y = value, color = gene_class)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~ mark, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Neuronal" = "#D73027", "Non-neuronal" = "#4575B4"),
                     name = "Gene class") +
  labs(title = "Per-gene mark profiles: MeCP2-Up + K119ub-Up quadrant",
       x = NULL, y = "Change (log2FC or difference)") +
  theme_biomodal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        strip.text = element_text(face = "bold"))

# =============================================================================
# SAVE
# =============================================================================

save_plot(p_state, "61i_chromatin_state_bar", w = 10, h = 6)
save_plot(p_comp, "61i_neuronal_vs_nonneuronal", w = 10, h = 9)
save_plot(p_strip, "61i_pergene_strip", w = 16, h = 14)

# Composite
p_composite <- (p_state | p_comp) +
  plot_annotation(
    title = "Section 61i: MeCP2-Up + K119ub-Up Quadrant Characterization",
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )
save_plot(p_composite, "61i_composite", w = 20, h = 9)

write.table(comp_results, file.path(TABLES_DIR, "61i_neuronal_vs_nonneuronal_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61i_neuronal_vs_nonneuronal_comparison.tsv\n")

write.table(quad_full, file.path(TABLES_DIR, "61i_quadrant_genes_full.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61i_quadrant_genes_full.tsv\n")

cat("\n================================================================================\n")
cat("SECTION 61i COMPLETE\n")
cat(sprintf("  72 quadrant genes: %d neuronal, %d non-neuronal\n", n_neur, n_nonneur))
cat("================================================================================\n")
