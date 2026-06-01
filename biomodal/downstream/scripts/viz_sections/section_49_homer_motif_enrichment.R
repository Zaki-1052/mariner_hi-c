# biomodal/downstream/scripts/viz_sections/section_49_homer_motif_enrichment.R
# Section 49: HOMER Motif Enrichment Visualization (K119ub + K27ac)
# Visualizes TF motif enrichment at differential H2AK119ub and H3K27ac sites.

source("scripts/viz_sections/_shared_config.R")
source("../../scripts/utils/multi_format_output.R")

if (!requireNamespace("ggwordcloud", quietly = TRUE)) {
  install.packages("ggwordcloud", repos = "https://cloud.r-project.org")
}
library(ggwordcloud)

cat("================================================================================\n")
cat("SECTION 49: HOMER MOTIF ENRICHMENT (K119ub + K27ac)\n")
cat("================================================================================\n\n")

# =============================================================================
# Configuration
# =============================================================================

HOMER_DIR <- file.path(BASE_DIR, "results/homer_motif_enrichment")
SECTION_DIR <- file.path(OUTPUT_DIR, "49_homer_motif_enrichment")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

COMPARISONS <- list(
  B1 = list(dir = "B1_k119ub_gained_vs_lost",     label = "K119ub Gained vs Lost"),
  B2 = list(dir = "B2_k119ub_lost_vs_gained",     label = "K119ub Lost vs Gained"),
  B3 = list(dir = "B3_diffbind_gained_vs_lost",   label = "DiffBind K119ub Gained vs Lost"),
  B4 = list(dir = "B4_diffbind_lost_vs_gained",   label = "DiffBind K119ub Lost vs Gained"),
  C1 = list(dir = "C1_k27ac_gained_vs_lost",      label = "K27ac Gained vs Lost"),
  C2 = list(dir = "C2_k27ac_lost_vs_gained",      label = "K27ac Lost vs Gained"),
  C3 = list(dir = "C3_k27ac_gained_vs_genome",    label = "K27ac Gained vs Genome"),
  C4 = list(dir = "C4_k27ac_lost_vs_genome",      label = "K27ac Lost vs Genome")
)

MOTIFS_OF_INTEREST <- c(
  "REST", "YY1", "SUZ12", "EZH2",
  "NeuroD1", "Mef2c", "Mef2a", "CREB",
  "Atoh1", "Olig2", "Sp1", "Ronin",
  "NeuroG2", "Ascl1", "TCF4", "Mesp1"
)

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

  df$fold_enrichment <- ifelse(df$bg_pct > 0, df$target_pct / df$bg_pct, NA_real_)
  df$neg_log10_q <- -log10(pmax(df$qvalue, 1e-300))
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
  cat(sprintf("  %s (%s): %d\n", sig_counts$comparison[i],
              sig_counts$comp_label[i], sig_counts$n_significant[i]))
}
cat("\n")

# =============================================================================
# Plot 1: Dot plot — Top 15 motifs per key comparison
# =============================================================================

cat("Plot 1: Dot plot (top 15 motifs, 4 key comparisons)...\n")

KEY_COMPS <- c("B1", "B3", "C1", "C2")

dot_data <- all_results %>%
  filter(comparison %in% KEY_COMPS, qvalue < 0.05) %>%
  group_by(comparison, tf_name) %>%
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
  group_by(comparison) %>%
  slice_min(order_by = pvalue, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    comp_label = factor(comp_label,
                        levels = sapply(KEY_COMPS, function(x) COMPARISONS[[x]]$label))
  )

dot_data <- dot_data %>%
  mutate(rank_label = paste0(comparison, ": ", tf_name)) %>%
  arrange(comparison, desc(pvalue)) %>%
  mutate(rank_label = factor(rank_label, levels = unique(rank_label)))

p_dot <- ggplot(dot_data, aes(x = neg_log10_p, y = rank_label)) +
  scale_y_discrete(labels = function(x) sub("^[^:]+: ", "", x)) +
  geom_point(aes(size = target_pct, color = fold_enrichment)) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Fold Enrichment\n(target/bg)") +
  scale_size_continuous(range = c(2, 7), name = "% Target\nwith Motif") +
  facet_wrap(~ comp_label, scales = "free_y", ncol = 2) +
  labs(
    x = expression(-log[10](p-value)),
    y = NULL,
    title = "Top Enriched TF Motifs at Differential Chromatin Sites",
    subtitle = "HOMER findMotifsGenome.pl - known motifs, q < 0.05, ranked by p-value"
  ) +
  theme_biomodal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

save_multiformat_ggplot(p_dot, file.path(SECTION_DIR, "homer_dotplot_top15"),
                        width = 14, height = 12)
cat("  Saved: homer_dotplot_top15\n\n")

# =============================================================================
# Plot 2: TF family bar chart
# =============================================================================

cat("Plot 2: TF family summary bar chart...\n")

family_data <- all_results %>%
  filter(comparison %in% KEY_COMPS, qvalue < 0.05) %>%
  count(comparison, comp_label, family, name = "n_motifs") %>%
  group_by(comparison) %>%
  slice_max(order_by = n_motifs, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(comp_label = factor(comp_label,
                             levels = sapply(KEY_COMPS, function(x) COMPARISONS[[x]]$label)))

top_families <- family_data %>%
  group_by(family) %>%
  summarise(total = sum(n_motifs)) %>%
  slice_max(order_by = total, n = 10) %>%
  pull(family)

family_data_top <- family_data %>% filter(family %in% top_families)

p_family <- ggplot(family_data_top, aes(x = reorder(family, -n_motifs), y = n_motifs,
                                         fill = comp_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = setNames(
    c("#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A"),
    sapply(KEY_COMPS, function(x) COMPARISONS[[x]]$label)
  ), name = "Comparison") +
  labs(
    x = "TF Family",
    y = "Number of Significant Motifs (q < 0.05)",
    title = "TF Family Enrichment at Differential Chromatin Sites",
    subtitle = "bHLH transcription factors dominate H2AK119ub-gained regions"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10))

save_multiformat_ggplot(p_family, file.path(SECTION_DIR, "homer_family_barchart"),
                        width = 12, height = 7)
cat("  Saved: homer_family_barchart\n\n")

# =============================================================================
# Plot 3: Word clouds (B1, B3, C2)
# =============================================================================

cat("Plot 3: Word clouds (B1, B3, C2)...\n")

CLOUD_COMPS <- c("B1", "B3", "C2")

cloud_data <- all_results %>%
  filter(comparison %in% CLOUD_COMPS, qvalue < 0.05) %>%
  group_by(comparison, tf_name) %>%
  slice_min(order_by = qvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    word_size = pmin(neg_log10_p, 200),
    comp_label = factor(comp_label,
                        levels = sapply(CLOUD_COMPS, function(x) COMPARISONS[[x]]$label))
  )

p_cloud <- ggplot(cloud_data, aes(label = tf_name, size = word_size,
                                   color = word_size)) +
  geom_text_wordcloud_area(rm_outside = TRUE, eccentricity = 0.65,
                            shape = "circle") +
  scale_size_area(max_size = 14) +
  scale_color_gradient(low = "#4393C3", high = "#D6604D") +
  facet_wrap(~ comp_label, ncol = 3) +
  labs(
    title = "TF Motif Enrichment Word Clouds",
    subtitle = "Word size proportional to -log10(p-value)"
  ) +
  theme_biomodal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "grey70", fill = NA)
  )

save_multiformat_ggplot(p_cloud, file.path(SECTION_DIR, "homer_wordclouds"),
                        width = 18, height = 7)
cat("  Saved: homer_wordclouds\n\n")

# =============================================================================
# Plot 4: Motifs-of-interest heatmap across all B/C comparisons
# =============================================================================

cat("Plot 4: Motifs-of-interest heatmap...\n")

moi_data <- all_results %>%
  filter(tf_name %in% MOTIFS_OF_INTEREST) %>%
  group_by(comparison, comp_label, tf_name) %>%
  slice_min(order_by = qvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

moi_matrix <- moi_data %>%
  mutate(neg_log10_q_capped = pmin(neg_log10_q, 100)) %>%
  dplyr::select(tf_name, comparison, neg_log10_q_capped) %>%
  pivot_wider(names_from = comparison, values_from = neg_log10_q_capped,
              values_fill = 0) %>%
  column_to_rownames("tf_name")

comp_order <- c("B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4")
comp_order <- intersect(comp_order, colnames(moi_matrix))
moi_matrix <- moi_matrix[, comp_order, drop = FALSE]

col_labels <- sapply(comp_order, function(x) COMPARISONS[[x]]$label)

sig_mat <- moi_data %>%
  mutate(is_sig = ifelse(qvalue < 0.05, "*", "")) %>%
  dplyr::select(tf_name, comparison, is_sig) %>%
  pivot_wider(names_from = comparison, values_from = is_sig, values_fill = "") %>%
  column_to_rownames("tf_name")
sig_mat <- sig_mat[rownames(moi_matrix), comp_order, drop = FALSE]

color_breaks <- seq(0, 100, length.out = 101)
color_palette <- colorRampPalette(c("white", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"))(100)

col_groups <- data.frame(
  Mark = c(rep("H2AK119ub", sum(grepl("^B", comp_order))),
           rep("H3K27ac", sum(grepl("^C", comp_order)))),
  row.names = comp_order
)

mark_colors <- list(Mark = c("H2AK119ub" = "#984EA3", "H3K27ac" = "#FF7F00"))

save_multiformat_pheatmap(
  quote(pheatmap::pheatmap(
    as.matrix(moi_matrix),
    display_numbers = as.matrix(sig_mat),
    number_color = "black",
    fontsize_number = 14,
    color = color_palette,
    breaks = color_breaks,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_col = col_labels,
    annotation_col = col_groups,
    annotation_colors = mark_colors,
    main = "Motifs of Interest - Enrichment Across K119ub & K27ac Comparisons",
    fontsize = 10,
    fontsize_row = 11,
    fontsize_col = 8,
    angle_col = 45,
    border_color = "grey90"
  )),
  file.path(SECTION_DIR, "homer_motifs_of_interest_heatmap"),
  width = 14, height = 8
)
cat("  Saved: homer_motifs_of_interest_heatmap\n\n")

# =============================================================================
# Summary table export
# =============================================================================

cat("Exporting summary tables...\n")

summary_table <- all_results %>%
  filter(comparison %in% KEY_COMPS, qvalue < 0.05) %>%
  group_by(comparison) %>%
  slice_min(order_by = qvalue, n = 25, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(comparison, comp_label, tf_name, family, consensus, pvalue, qvalue,
         target_pct, bg_pct, fold_enrichment) %>%
  arrange(comparison, qvalue)

write_tsv(summary_table, file.path(SECTION_DIR, "homer_top25_significant_motifs.tsv"))
cat(sprintf("  Saved: homer_top25_significant_motifs.tsv (%d rows)\n", nrow(summary_table)))

moi_table <- moi_data %>%
  dplyr::select(comparison, comp_label, tf_name, family, pvalue, qvalue,
         target_pct, bg_pct, fold_enrichment) %>%
  arrange(tf_name, comparison)

write_tsv(moi_table, file.path(SECTION_DIR, "homer_motifs_of_interest_all.tsv"))
cat(sprintf("  Saved: homer_motifs_of_interest_all.tsv (%d rows)\n", nrow(moi_table)))

# =============================================================================
# Done
# =============================================================================

cat("\n================================================================================\n")
cat("Section 49 complete.\n")
cat(sprintf("Output directory: %s\n", SECTION_DIR))
cat("================================================================================\n")
