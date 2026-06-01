# biomodal/downstream/scripts/viz_sections/section_49_homer_k119ub_k27ac.R
# Section 49: HOMER Motif Enrichment — H2AK119ub + H3K27ac differential sites
# Matrix-style dot plots (compareCluster layout) + TF family bar chart.

source("scripts/viz_sections/_shared_config.R")
source("../../scripts/utils/multi_format_output.R")

cat("================================================================================\n")
cat("SECTION 49: HOMER MOTIF ENRICHMENT (K119ub + K27ac)\n")
cat("================================================================================\n\n")

HOMER_DIR <- file.path(BASE_DIR, "results/homer_motif_enrichment")
SECTION_DIR <- file.path(OUTPUT_DIR, "49_homer_k119ub_k27ac")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N <- 20

# =============================================================================
# Comparison definitions
# =============================================================================

COMPARISONS <- list(
  B1 = list(dir = "B1_k119ub_gained_vs_lost",   label = "Gained\nvs Lost"),
  B2 = list(dir = "B2_k119ub_lost_vs_gained",   label = "Lost\nvs Gained"),
  B3 = list(dir = "B3_diffbind_gained_vs_lost",  label = "DiffBind\nGained vs Lost"),
  B4 = list(dir = "B4_diffbind_lost_vs_gained",  label = "DiffBind\nLost vs Gained"),
  C1 = list(dir = "C1_k27ac_gained_vs_lost",     label = "Gained\nvs Lost"),
  C2 = list(dir = "C2_k27ac_lost_vs_gained",     label = "Lost\nvs Gained"),
  C3 = list(dir = "C3_k27ac_gained_vs_genome",   label = "Gained\nvs Genome"),
  C4 = list(dir = "C4_k27ac_lost_vs_genome",     label = "Lost\nvs Genome")
)

K119UB_IDS <- c("B1", "B2", "B3", "B4")
K27AC_IDS  <- c("C1", "C2", "C3", "C4")

# =============================================================================
# Load all knownResults.txt
# =============================================================================

cat("Loading HOMER knownResults...\n")

load_known_results <- function(comp_id, comp_info) {
  fpath <- file.path(HOMER_DIR, comp_info$dir, "knownResults.txt")
  if (!file.exists(fpath)) {
    cat(sprintf("  WARNING: %s not found\n", fpath))
    return(NULL)
  }

  df <- read.delim(fpath, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE)
  colnames(df) <- c("motif_name", "consensus", "pvalue", "log_pvalue",
                     "qvalue", "target_count", "target_pct",
                     "bg_count", "bg_pct")

  df$target_pct <- as.numeric(gsub("%", "", df$target_pct))
  df$bg_pct     <- as.numeric(gsub("%", "", df$bg_pct))
  df$pvalue     <- as.numeric(df$pvalue)
  df$qvalue     <- as.numeric(df$qvalue)

  df$tf_name <- trimws(sub("\\(.*", "", df$motif_name))
  df$family  <- sub(".*?\\(([^)]+)\\).*", "\\1", df$motif_name)
  df$family[!grepl("\\(", df$motif_name)] <- "Other"

  df$fold_enrichment <- ifelse(df$bg_pct > 0.05, df$target_pct / df$bg_pct, NA_real_)
  df$neg_log10_p <- -log10(pmax(df$pvalue, 1e-300))

  df$comparison <- comp_id
  df$comp_label <- comp_info$label
  df
}

all_results <- bind_rows(lapply(names(COMPARISONS), function(id) {
  load_known_results(id, COMPARISONS[[id]])
}))

cat(sprintf("  Loaded %d motif-comparison pairs across %d comparisons\n",
            nrow(all_results), n_distinct(all_results$comparison)))

sig_counts <- all_results %>%
  filter(qvalue < 0.05) %>%
  count(comparison, comp_label, name = "n_significant")
cat("\nSignificant motifs (q < 0.05) per comparison:\n")
for (i in seq_len(nrow(sig_counts))) {
  cat(sprintf("  %s: %d\n", sig_counts$comparison[i], sig_counts$n_significant[i]))
}
cat("\n")

# =============================================================================
# Helper: build matrix dot plot data
# =============================================================================

build_matrix_dotplot_data <- function(comp_ids, top_n = TOP_N) {
  subset <- all_results %>%
    filter(comparison %in% comp_ids, qvalue < 0.05) %>%
    group_by(comparison, tf_name) %>%
    slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()

  top_per_comp <- subset %>%
    group_by(comparison) %>%
    slice_min(order_by = pvalue, n = top_n, with_ties = FALSE) %>%
    ungroup()

  union_tfs <- unique(top_per_comp$tf_name)

  global_order <- subset %>%
    filter(tf_name %in% union_tfs) %>%
    group_by(tf_name) %>%
    summarise(best_p = min(pvalue), .groups = "drop") %>%
    arrange(best_p) %>%
    pull(tf_name)

  subset %>%
    filter(tf_name %in% union_tfs) %>%
    mutate(
      tf_name = factor(tf_name, levels = rev(global_order)),
      comparison = factor(comparison, levels = comp_ids),
      comp_label = factor(comp_label,
                          levels = sapply(comp_ids, function(x) COMPARISONS[[x]]$label))
    )
}

# =============================================================================
# Plot 1a: K119ub matrix dot plot (B1-B4)
# =============================================================================

cat("Plot 1a: H2AK119ub matrix dot plot (B1-B4)...\n")

k119_data <- build_matrix_dotplot_data(K119UB_IDS)
cat(sprintf("  %d unique TFs in union set\n", n_distinct(k119_data$tf_name)))

p_k119 <- ggplot(k119_data, aes(x = comp_label, y = tf_name)) +
  geom_point(aes(size = target_pct, fill = neg_log10_p),
             shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C",
                      name = expression(-log[10](p))) +
  scale_size_continuous(range = c(2, 8), name = "% Target\nwith Motif") +
  labs(
    x = NULL, y = NULL,
    title = "H2AK119ub: TF Motif Enrichment at Differential Sites",
    subtitle = sprintf("Top %d motifs per comparison (q < 0.05), HOMER known motifs", TOP_N)
  ) +
  theme_biomodal() +
  theme(
    axis.text.x = element_text(size = 9, lineheight = 0.9),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey92"),
    legend.position = "right"
  )

save_multiformat_ggplot(p_k119, file.path(SECTION_DIR, "homer_k119ub_matrix_dotplot"),
                        width = 11, height = 12)
cat("  Saved: homer_k119ub_matrix_dotplot\n\n")

# =============================================================================
# Plot 1b: K27ac matrix dot plot (C1-C4)
# =============================================================================

cat("Plot 1b: H3K27ac matrix dot plot (C1-C4)...\n")

k27_data <- build_matrix_dotplot_data(K27AC_IDS)
cat(sprintf("  %d unique TFs in union set\n", n_distinct(k27_data$tf_name)))

p_k27 <- ggplot(k27_data, aes(x = comp_label, y = tf_name)) +
  geom_point(aes(size = target_pct, fill = neg_log10_p),
             shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C",
                      name = expression(-log[10](p))) +
  scale_size_continuous(range = c(2, 8), name = "% Target\nwith Motif") +
  labs(
    x = NULL, y = NULL,
    title = "H3K27ac: TF Motif Enrichment at Differential Sites",
    subtitle = sprintf("Top %d motifs per comparison (q < 0.05), HOMER known motifs", TOP_N)
  ) +
  theme_biomodal() +
  theme(
    axis.text.x = element_text(size = 9, lineheight = 0.9),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey92"),
    legend.position = "right"
  )

save_multiformat_ggplot(p_k27, file.path(SECTION_DIR, "homer_k27ac_matrix_dotplot"),
                        width = 11, height = 12)
cat("  Saved: homer_k27ac_matrix_dotplot\n\n")

# =============================================================================
# Plot 2: TF family bar chart (all 8 comparisons, faceted by mark)
# =============================================================================

cat("Plot 2: TF family bar chart (all 8 comparisons)...\n")

family_data <- all_results %>%
  filter(qvalue < 0.05) %>%
  mutate(mark = ifelse(grepl("^B", comparison), "H2AK119ub", "H3K27ac")) %>%
  count(comparison, comp_label, mark, family, name = "n_motifs")

top_families <- family_data %>%
  group_by(family) %>%
  summarise(total = sum(n_motifs), .groups = "drop") %>%
  slice_max(order_by = total, n = 10) %>%
  pull(family)

BAR_LABELS <- c(
  B1 = "Gained vs Lost", B2 = "Lost vs Gained",
  B3 = "DiffBind Gained vs Lost", B4 = "DiffBind Lost vs Gained",
  C1 = "Gained vs Lost", C2 = "Lost vs Gained",
  C3 = "Gained vs Genome", C4 = "Lost vs Genome"
)

family_plot_data <- family_data %>%
  filter(family %in% top_families) %>%
  mutate(
    bar_label = BAR_LABELS[comparison],
    mark = factor(mark, levels = c("H2AK119ub", "H3K27ac"))
  )

mark_palettes <- setNames(
  c("#7B3294", "#C2A5CF", "#E7298A", "#F1B6DA",
    "#1B7837", "#7FBC41", "#D95F02", "#FDAE61"),
  BAR_LABELS
)

p_family <- ggplot(family_plot_data,
                   aes(x = reorder(family, -n_motifs), y = n_motifs,
                       fill = bar_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = mark_palettes, name = "Comparison") +
  facet_wrap(~ mark, ncol = 2) +
  labs(
    x = "TF Family",
    y = "Significant Motifs (q < 0.05)",
    title = "TF Family Enrichment at Differential Chromatin Sites",
    subtitle = "H2AK119ub (left) and H3K27ac (right), all comparisons"
  ) +
  theme_biomodal() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
    strip.text = element_text(size = 12, face = "bold")
  )

save_multiformat_ggplot(p_family, file.path(SECTION_DIR, "homer_family_barchart_all"),
                        width = 16, height = 8)
cat("  Saved: homer_family_barchart_all\n\n")

# =============================================================================
# Summary table exports
# =============================================================================

cat("Exporting summary tables...\n")

summary_all <- all_results %>%
  filter(qvalue < 0.05) %>%
  group_by(comparison) %>%
  slice_min(order_by = pvalue, n = 25, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(comparison, comp_label, tf_name, family, consensus,
                pvalue, qvalue, target_pct, bg_pct, fold_enrichment) %>%
  arrange(comparison, pvalue)

write_tsv(summary_all, file.path(SECTION_DIR, "homer_top25_per_comparison.tsv"))
cat(sprintf("  Saved: homer_top25_per_comparison.tsv (%d rows)\n", nrow(summary_all)))

sig_all <- all_results %>%
  filter(qvalue < 0.05) %>%
  dplyr::select(comparison, comp_label, tf_name, family, consensus,
                pvalue, qvalue, target_pct, bg_pct, fold_enrichment) %>%
  arrange(comparison, pvalue)

write_tsv(sig_all, file.path(SECTION_DIR, "homer_all_significant_motifs.tsv"))
cat(sprintf("  Saved: homer_all_significant_motifs.tsv (%d rows)\n", nrow(sig_all)))

cat("\n================================================================================\n")
cat("Section 49 complete.\n")
cat(sprintf("Output directory: %s\n", SECTION_DIR))
cat("================================================================================\n")
