# abc/scripts/step11b_homer_motif_visualization.R
# HOMER motif enrichment visualization for BAP1-KO enhancer subsets.
#
# Visualizes known motif enrichment results from findMotifsGenome.pl across
# 4 comparisons of enhancer subsets: Activity_Lost, K119ub_Only, Activity_Gain
# vs Stable, plus T3_high vs T1_low (K119ub dose dependence).
#
# Inputs (relative to abc/ working directory):
#   results/enhancer_subset_analysis/homer_results/homer_results/*/knownResults.txt
#
# Outputs:
#   results/enhancer_subset_analysis/homer_motif_viz/  -- 12 figures + summary TSV
#
# Usage:
#   cd abc && Rscript scripts/step11b_homer_motif_visualization.R

# =============================================================================
# CONFIGURATION
# =============================================================================

HOMER_BASE <- "abc/results/enhancer_subset_analysis/homer_results/homer_results"  # Note: repo-relative path, not bundled in data/

COMPARISONS <- c(
  "Activity_Lost_vs_Stable",
  "K119ub_Only_vs_Stable",
  "Activity_Gain_vs_Stable",
  "T3_high_vs_T1_low"
)

COMPARISON_LABELS <- c(
  Activity_Lost_vs_Stable = "Activity Lost vs Stable",
  K119ub_Only_vs_Stable   = "K119ub-Only vs Stable",
  Activity_Gain_vs_Stable = "Activity Gain vs Stable",
  T3_high_vs_T1_low       = "T3 High vs T1 Low (K119ub Dose)"
)

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

TF_FAMILY_COLORS <- c(
  bHLH     = "#E41A1C",
  ETS      = "#377EB8",
  Zf       = "#4DAF4A",
  Homeobox = "#984EA3",
  bZIP     = "#FF7F00",
  MADS     = "#A6761D",
  Forkhead = "#A65628",
  `T-box`  = "#F781BF",
  NR       = "#666666",
  HTH      = "#66C2A5",
  CTF      = "#FC8D62",
  Other    = "#BBBBBB"
)
FAMILY_ORDER <- names(TF_FAMILY_COLORS)

OUTPUT_DIR <- "data/plots/supplemental"  # Original: results/enhancer_subset_analysis/homer_motif_viz

cat("================================================================================\n")
cat("STEP 11B: HOMER MOTIF ENRICHMENT VISUALIZATION\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
homer_files <- setNames(
  file.path(HOMER_BASE, COMPARISONS, "knownResults.txt"),
  COMPARISONS
)
for (nm in names(homer_files)) {
  stopifnot(file.exists(homer_files[[nm]]))
  cat(sprintf("  [OK] %s\n", nm))
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutput directory: %s\n\n", OUTPUT_DIR))

# =============================================================================
# HELPERS
# =============================================================================

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

UTIL_PATH <- "data/scripts/_shared/multi_format_output.R"  # Original: UTIL_PATH <- "../scripts/utils/multi_format_output.R"
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

save_heatmap <- function(expr, name, w = 8, h = 10) {
  save_multiformat_pheatmap(expr, base_path = file.path(OUTPUT_DIR, name),
                            width = w, height = h, dpi = 300,
                            verbose = TRUE, use_subfolders = TRUE)
}

fmt_p <- function(p) {
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# TF FAMILY ASSIGNMENT
# =============================================================================

assign_family_group <- function(fam) {
  fam <- trimws(fam)
  if (fam %in% c("bHLH", "HLH")) return("bHLH")
  if (fam == "ETS") return("ETS")
  if (fam %in% c("Zf", "Zf?")) return("Zf")
  if (fam == "Homeobox") return("Homeobox")
  if (fam == "bZIP") return("bZIP")
  if (fam == "MADS") return("MADS")
  if (fam %in% c("Forkhead", "forkhead")) return("Forkhead")
  if (fam == "T-box") return("T-box")
  if (fam == "NR") return("NR")
  if (fam == "HTH") return("HTH")
  if (fam %in% c("CTF", "CCAAT")) return("CTF")
  # Compound families — assign by primary component
  if (grepl("bHLH", fam)) return("bHLH")
  if (grepl("ETS", fam)) return("ETS")
  if (grepl("Homeobox|POU", fam)) return("Homeobox")
  if (grepl("[Ff]orkhead", fam)) return("Forkhead")
  if (grepl("bZIP", fam)) return("bZIP")
  if (grepl("Zf", fam)) return("Zf")
  if (grepl("T-box", fam)) return("T-box")
  if (grepl("HTH", fam)) return("HTH")
  if (grepl("CTF", fam)) return("CTF")
  if (grepl("NR", fam)) return("NR")
  "Other"
}

# =============================================================================
# PARSE HOMER RESULTS
# =============================================================================

parse_homer_results <- function(path, comparison_name) {
  hdr <- readLines(path, n = 1)
  cols <- strsplit(hdr, "\t")[[1]]
  n_target_total <- as.integer(sub(".*\\(of (\\d+)\\).*", "\\1",
                                   grep("Target", cols, value = TRUE)[1]))
  n_bg_total <- as.integer(sub(".*\\(of (\\d+)\\).*", "\\1",
                                grep("Background", cols, value = TRUE)[1]))

  raw <- read.delim(path, header = FALSE, skip = 1, stringsAsFactors = FALSE,
                    quote = "")
  colnames(raw) <- c("motif_name", "consensus", "pvalue", "log_pvalue", "qvalue",
                      "n_target", "pct_target", "n_background", "pct_background")

  raw$pct_target <- as.numeric(sub("%", "", raw$pct_target))
  raw$pct_background <- as.numeric(sub("%", "", raw$pct_background))

  raw$tf_name <- sub("^([^(]+)\\(.*", "\\1", raw$motif_name)
  raw$family_raw <- sub("^[^(]+\\(([^)]+)\\).*", "\\1", raw$motif_name)

  raw$fold_enrichment <- raw$pct_target / pmax(raw$pct_background, 0.01)
  raw$neg_log10_p <- pmin(-raw$log_pvalue / log(10), 300)
  raw$comparison <- comparison_name
  raw$n_target_total <- n_target_total
  raw$n_background_total <- n_bg_total
  raw$family_group <- sapply(raw$family_raw, assign_family_group, USE.NAMES = FALSE)
  raw$family_group <- factor(raw$family_group, levels = FAMILY_ORDER)

  raw
}

cat("Parsing HOMER results...\n")
all_results <- lapply(names(homer_files), function(comp) {
  cat(sprintf("  Parsing %s...\n", comp))
  parse_homer_results(homer_files[[comp]], comp)
})
all_df <- do.call(rbind, all_results)
cat(sprintf("  Total motif-comparison pairs: %d\n", nrow(all_df)))

# Deduplicated view: one entry per TF per comparison (keep most significant)
all_dedup <- all_df %>%
  group_by(comparison, tf_name) %>%
  slice_max(neg_log10_p, n = 1, with_ties = FALSE) %>%
  ungroup()
cat(sprintf("  Unique TF-comparison pairs: %d\n\n", nrow(all_dedup)))

# =============================================================================
# SECTION A: TOP ENRICHED MOTIFS PER COMPARISON (lollipop plots)
# =============================================================================

cat("=== Section A: Top enriched motifs per comparison ===\n")

for (i in seq_along(COMPARISONS)) {
  comp <- COMPARISONS[i]
  top20 <- all_dedup %>%
    filter(comparison == comp) %>%
    arrange(desc(neg_log10_p)) %>%
    head(20)

  top20$tf_name <- factor(top20$tf_name, levels = rev(top20$tf_name))

  n_tgt <- top20$n_target_total[1]
  n_bg  <- top20$n_background_total[1]
  n_sig <- sum(all_dedup$comparison == comp & all_dedup$qvalue < 0.05)

  p <- ggplot(top20, aes(x = neg_log10_p, y = tf_name)) +
    geom_segment(aes(x = 0, xend = neg_log10_p, y = tf_name, yend = tf_name),
                 color = "grey70", linewidth = 0.5) +
    geom_point(aes(color = family_group), size = 3) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = TF_FAMILY_COLORS, name = "TF Family", drop = TRUE) +
    labs(
      title = sprintf("Top 20 Enriched Motifs: %s", COMPARISON_LABELS[comp]),
      subtitle = sprintf("Target: n=%s | Background: n=%s | Significant (q<0.05): %d",
                         format(n_tgt, big.mark = ","),
                         format(n_bg, big.mark = ","), n_sig),
      x = expression(-log[10](p-value)),
      y = NULL
    ) +
    theme_pub +
    theme(legend.position = "right")

  plot_name <- sprintf("%02d_top_motifs_%s", i, comp)
  save_plot(p, plot_name, w = 8, h = 7)
}

# =============================================================================
# SECTION B: POLYCOMB-RELATED MOTIF FOCUS (K119ub_Only_vs_Stable)
# =============================================================================

cat("\n=== Section B: Polycomb-related motif focus ===\n")

k119_df <- all_df %>% filter(comparison == "K119ub_Only_vs_Stable")
k119_dedup <- all_dedup %>% filter(comparison == "K119ub_Only_vs_Stable")

# --- Plot 05: Volcano-style (fold enrichment vs significance) ---

k119_dedup <- k119_dedup %>%
  mutate(
    sig_class = case_when(
      tf_name == "YY1"       ~ "YY1",
      tf_name == "REST-NRSF" ~ "REST",
      qvalue < 0.05          ~ "Significant (q<0.05)",
      TRUE                   ~ "Not significant"
    ),
    label = case_when(
      tf_name == "YY1"       ~ "YY1",
      tf_name == "REST-NRSF" ~ "REST",
      TRUE                   ~ NA_character_
    )
  )

sig_colors <- c(
  "YY1"                   = "#E41A1C",
  "REST"                  = "#377EB8",
  "Significant (q<0.05)"  = "grey30",
  "Not significant"       = "grey80"
)

label_df <- k119_dedup %>% filter(!is.na(label))

p05 <- ggplot(k119_dedup, aes(x = fold_enrichment, y = neg_log10_p)) +
  geom_point(aes(color = sig_class), size = 1.8, alpha = 0.7) +
  geom_text_repel(
    data = label_df, aes(label = label, color = sig_class),
    size = 4, fontface = "bold",
    box.padding = 0.5, max.overlaps = 20, show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = sig_colors, name = NULL) +
  labs(
    title = "Motif Enrichment in K119ub-Only Enhancers",
    subtitle = paste0(
      "PRC1 subunits (RING1B, BMI1) and PRC2 (SUZ12, EZH2) lack DNA-binding motifs\n",
      "and are absent from HOMER's database -- they are recruited by TFs like YY1"
    ),
    x = "Fold Enrichment (target / background)",
    y = expression(-log[10](p-value))
  ) +
  theme_pub

save_plot(p05, "05_volcano_K119ub_Only_vs_Stable", w = 8, h = 7)

# --- Plot 06: YY1 and REST target vs background bar chart ---

yy1_row <- k119_dedup %>% filter(tf_name == "YY1")
rest_row <- k119_dedup %>% filter(tf_name == "REST-NRSF")

focus_df <- bind_rows(
  data.frame(tf = "YY1", group = "K119ub-Only (target)",
             pct = yy1_row$pct_target, stringsAsFactors = FALSE),
  data.frame(tf = "YY1", group = "Stable (background)",
             pct = yy1_row$pct_background, stringsAsFactors = FALSE),
  data.frame(tf = "REST", group = "K119ub-Only (target)",
             pct = rest_row$pct_target, stringsAsFactors = FALSE),
  data.frame(tf = "REST", group = "Stable (background)",
             pct = rest_row$pct_background, stringsAsFactors = FALSE)
)

qval_labels <- data.frame(
  tf = c("YY1", "REST"),
  label = c(sprintf("q = %.3f", yy1_row$qvalue),
            sprintf("q = %.3f", rest_row$qvalue)),
  y_pos = c(max(yy1_row$pct_target, yy1_row$pct_background) + 0.3,
            max(rest_row$pct_target, rest_row$pct_background) + 0.05),
  stringsAsFactors = FALSE
)

p06 <- ggplot(focus_df, aes(x = tf, y = pct, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(data = qval_labels, aes(x = tf, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.5, fontface = "italic") +
  scale_fill_manual(
    values = c("K119ub-Only (target)" = "#B2182B", "Stable (background)" = "#D1E5F0"),
    name = NULL
  ) +
  labs(
    title = "PRC1-Recruiting TF Motifs in K119ub-Only Enhancers",
    subtitle = "YY1 recruits PRC1 to chromatin; REST (silencer TF) is NOT enriched",
    x = NULL,
    y = "% of Sequences with Motif"
  ) +
  theme_pub

save_plot(p06, "06_yy1_rest_barplot", w = 6, h = 5)

# =============================================================================
# SECTION C: CROSS-COMPARISON HEATMAP
# =============================================================================

cat("\n=== Section C: Cross-comparison heatmap ===\n")

# Top 25 TFs per comparison (deduplicated)
top_per_comp <- all_dedup %>%
  group_by(comparison) %>%
  arrange(desc(neg_log10_p)) %>%
  slice_head(n = 25) %>%
  ungroup()

top_tfs <- unique(top_per_comp$tf_name)
cat(sprintf("  Union of top TFs across comparisons: %d\n", length(top_tfs)))

# Build fold enrichment matrix
mat_df <- all_dedup %>%
  filter(tf_name %in% top_tfs) %>%
  select(tf_name, comparison, fold_enrichment) %>%
  pivot_wider(names_from = comparison, values_from = fold_enrichment,
              values_fill = 1.0)

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$tf_name
colnames(mat) <- COMPARISON_LABELS[colnames(mat)]

heatmap_colors <- colorRampPalette(c("white", "#FEE0D2", "#FC9272", "#DE2D26"))(50)

save_heatmap(
  quote(pheatmap(
    mat,
    color = heatmap_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Motif Fold Enrichment Across Comparisons",
    fontsize_row = 7,
    fontsize_col = 9,
    border_color = "grey90"
  )),
  "07_cross_comparison_heatmap", w = 9, h = 12
)

# =============================================================================
# SECTION D: T3_HIGH VS T1_LOW DOSE DEPENDENCE
# =============================================================================

cat("\n=== Section D: Dose-dependence analysis ===\n")

# --- Plot 08: Paired dot plot ---

top15_k119 <- all_dedup %>%
  filter(comparison == "K119ub_Only_vs_Stable") %>%
  arrange(desc(neg_log10_p)) %>%
  head(15) %>%
  pull(tf_name)

dose_compare <- all_dedup %>%
  filter(comparison %in% c("K119ub_Only_vs_Stable", "T3_high_vs_T1_low"),
         tf_name %in% top15_k119)

k119_order <- dose_compare %>%
  filter(comparison == "K119ub_Only_vs_Stable") %>%
  arrange(neg_log10_p) %>%
  pull(tf_name)

dose_compare <- dose_compare %>%
  mutate(
    tf_name = factor(tf_name, levels = k119_order),
    comparison = factor(
      comparison,
      levels = c("K119ub_Only_vs_Stable", "T3_high_vs_T1_low"),
      labels = c("K119ub-Only vs Stable", "T3 High vs T1 Low")
    )
  )

p08 <- ggplot(dose_compare, aes(x = neg_log10_p, y = tf_name)) +
  geom_segment(aes(x = 0, xend = neg_log10_p, y = tf_name, yend = tf_name),
               color = "grey70", linewidth = 0.4) +
  geom_point(aes(color = family_group), size = 3) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  facet_wrap(~comparison, scales = "free_x") +
  scale_color_manual(values = TF_FAMILY_COLORS, name = "TF Family", drop = TRUE) +
  labs(
    title = "K119ub-Enriched Motifs Are Not Dose-Dependent",
    subtitle = "Same motifs enriched in K119ub-Only show no differential enrichment between K119ub dose tertiles",
    x = expression(-log[10](p-value)),
    y = NULL
  ) +
  theme_pub +
  theme(legend.position = "right")

save_plot(p08, "08_dose_comparison_dotplot", w = 11, h = 6)

# --- Plot 09: Family-level enrichment in T3_high_vs_T1_low ---

t3_family <- all_dedup %>%
  filter(comparison == "T3_high_vs_T1_low") %>%
  group_by(family_group) %>%
  arrange(desc(neg_log10_p)) %>%
  slice_head(n = 5) %>%
  summarize(mean_fold = mean(fold_enrichment), .groups = "drop") %>%
  arrange(desc(mean_fold))

t3_family$family_group <- factor(t3_family$family_group,
                                  levels = rev(t3_family$family_group))

p09 <- ggplot(t3_family, aes(x = mean_fold, y = family_group)) +
  geom_col(aes(fill = family_group), width = 0.7, show.legend = FALSE) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = TF_FAMILY_COLORS) +
  annotate("text", x = max(t3_family$mean_fold) * 0.8, y = nrow(t3_family) * 0.5,
           label = "No significant motifs (q<0.05)\nbetween K119ub dose tertiles",
           size = 3.5, fontface = "italic", color = "grey40") +
  labs(
    title = "TF Family Enrichment: T3 High vs T1 Low (K119ub Dose)",
    subtitle = "Mean fold enrichment of top 5 motifs per family",
    x = "Mean Fold Enrichment",
    y = NULL
  ) +
  theme_pub

save_plot(p09, "09_dose_family_barplot", w = 7, h = 5)

# =============================================================================
# SECTION E: TF FAMILY SUMMARY ACROSS COMPARISONS
# =============================================================================

cat("\n=== Section E: TF family summary ===\n")

# Mean -log10(p) of top 5 motifs per family per comparison
family_enrich <- all_dedup %>%
  group_by(comparison, family_group) %>%
  arrange(desc(neg_log10_p)) %>%
  slice_head(n = 5) %>%
  summarize(mean_neg_log10_p = mean(neg_log10_p), .groups = "drop")

# Count significant motifs per family per comparison
family_sig <- all_dedup %>%
  group_by(comparison, family_group) %>%
  summarize(n_sig = sum(qvalue < 0.05), .groups = "drop")

family_summary <- left_join(family_enrich, family_sig,
                            by = c("comparison", "family_group"))

# Clean comparison labels for display
family_summary <- family_summary %>%
  mutate(comp_label = COMPARISON_LABELS[comparison])

# --- Plot 10: Family-level heatmap (geom_tile) ---

p10 <- ggplot(family_summary,
              aes(x = comp_label, y = family_group, fill = mean_neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_sig), size = 3, color = "grey20") +
  scale_fill_gradient(low = "white", high = "#DE2D26",
                      name = expression(Mean~-log[10](p))) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    title = "TF Family Enrichment Across Comparisons",
    subtitle = "Fill: mean -log10(p) of top 5 motifs | Text: number significant (q<0.05)",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 8),
    legend.position = "right",
    panel.grid = element_blank()
  )

save_plot(p10, "10_family_summary_heatmap", w = 9, h = 7)

# --- Plot 11: Stacked bar chart of significant motif counts ---

sig_counts <- family_summary %>%
  filter(n_sig > 0)

p11 <- ggplot(sig_counts, aes(x = comp_label, y = n_sig, fill = family_group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = TF_FAMILY_COLORS, name = "TF Family") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    title = "Significant Motifs (q<0.05) by TF Family",
    subtitle = "T3 High vs T1 Low has zero significant motifs",
    x = NULL,
    y = "Number of Significant Motifs"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(size = 9))

save_plot(p11, "11_family_sig_counts", w = 8, h = 6)

# =============================================================================
# SECTION F: NARRATIVE SUMMARY COMPOSITE
# =============================================================================

cat("\n=== Section F: Narrative summary composite ===\n")

# Panel A: K119ub_Only lollipop (regenerate for consistent sizing)
top20_k119 <- all_dedup %>%
  filter(comparison == "K119ub_Only_vs_Stable") %>%
  arrange(desc(neg_log10_p)) %>%
  head(20) %>%
  mutate(tf_name = factor(tf_name, levels = rev(tf_name)))

panel_a <- ggplot(top20_k119, aes(x = neg_log10_p, y = tf_name)) +
  geom_segment(aes(x = 0, xend = neg_log10_p, y = tf_name, yend = tf_name),
               color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = family_group), size = 2.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = TF_FAMILY_COLORS, name = "TF Family", drop = TRUE) +
  labs(title = "A) K119ub-Only: Top Motifs",
       x = expression(-log[10](p)), y = NULL) +
  theme_pub +
  theme(legend.position = "none", plot.title = element_text(size = 10))

# Panel B: YY1/REST volcano (reuse p05 with adjusted theme)
k119_vol <- k119_dedup %>%
  mutate(
    sig_class = case_when(
      tf_name == "YY1"       ~ "YY1",
      tf_name == "REST-NRSF" ~ "REST",
      qvalue < 0.05          ~ "Significant",
      TRUE                   ~ "NS"
    ),
    label = case_when(
      tf_name == "YY1"       ~ "YY1",
      tf_name == "REST-NRSF" ~ "REST",
      TRUE                   ~ NA_character_
    )
  )

vol_labels <- k119_vol %>% filter(!is.na(label))

panel_b <- ggplot(k119_vol, aes(x = fold_enrichment, y = neg_log10_p)) +
  geom_point(aes(color = sig_class), size = 1.2, alpha = 0.6) +
  geom_text_repel(data = vol_labels, aes(label = label, color = sig_class),
                  size = 3, fontface = "bold", show.legend = FALSE,
                  box.padding = 0.4, max.overlaps = 20) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c(YY1 = "#E41A1C", REST = "#377EB8",
               Significant = "grey30", NS = "grey80"),
    name = NULL
  ) +
  labs(title = "B) YY1/REST in K119ub-Only",
       x = "Fold Enrichment", y = expression(-log[10](p))) +
  theme_pub +
  theme(legend.position = "none", plot.title = element_text(size = 10))

# Panel C: Cross-comparison heatmap (ggplot geom_tile version)
heatmap_df <- all_dedup %>%
  filter(tf_name %in% top_tfs) %>%
  select(tf_name, comparison, fold_enrichment) %>%
  mutate(comp_label = COMPARISON_LABELS[comparison])

# Order TFs by hierarchical clustering on fold enrichment
fe_mat <- mat_df %>% select(-tf_name) %>% as.matrix()
rownames(fe_mat) <- mat_df$tf_name
tf_clust <- hclust(dist(fe_mat))
tf_order <- tf_clust$labels[tf_clust$order]

heatmap_df$tf_name <- factor(heatmap_df$tf_name, levels = tf_order)

panel_c <- ggplot(heatmap_df, aes(x = comp_label, y = tf_name, fill = fold_enrichment)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#DE2D26", name = "Fold\nEnrich.") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(title = "C) Cross-Comparison Heatmap", x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 7),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(size = 10)
  )

# Panel D: Family summary heatmap
panel_d <- ggplot(family_summary,
                  aes(x = comp_label, y = family_group, fill = mean_neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_sig), size = 2.5, color = "grey20") +
  scale_fill_gradient(low = "white", high = "#DE2D26",
                      name = expression(atop(Mean, -log[10](p)))) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(title = "D) TF Family Summary", x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 7),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(size = 10)
  )

composite <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    title = "HOMER Motif Enrichment: Sequence Features of K119ub-Targeted Enhancers",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(composite, "12_narrative_composite", w = 16, h = 14)

# =============================================================================
# SUMMARY STATISTICS TABLE
# =============================================================================

cat("\n=== Writing summary statistics ===\n")

summary_stats <- all_dedup %>%
  group_by(comparison) %>%
  summarize(
    n_motifs = n(),
    n_sig_q05 = sum(qvalue < 0.05),
    top_motif = tf_name[which.max(neg_log10_p)],
    top_pvalue = pvalue[which.max(neg_log10_p)],
    top_neg_log10_p = max(neg_log10_p),
    .groups = "drop"
  )

# Dominant family: family with most significant motifs per comparison
dominant_fam <- all_dedup %>%
  filter(qvalue < 0.05) %>%
  group_by(comparison, family_group) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(comparison) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  select(comparison, dominant_family = family_group) %>%
  ungroup()

summary_stats <- left_join(summary_stats, dominant_fam, by = "comparison")
summary_stats$dominant_family <- as.character(summary_stats$dominant_family)
summary_stats$dominant_family[is.na(summary_stats$dominant_family)] <- "None (no sig. motifs)"

summary_stats$comparison_label <- COMPARISON_LABELS[summary_stats$comparison]

write.table(summary_stats,
            file.path("data/tsvs/supplemental", "homer_motif_summary_stats.tsv"),  # Original: file.path(OUTPUT_DIR, "homer_motif_summary_stats.tsv")
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: homer_motif_summary_stats.tsv\n")

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 11B COMPLETE\n")
cat(sprintf("Output: %s/\n", OUTPUT_DIR))
cat(sprintf("  Figures: 12 (each in PDF/SVG/JPG)\n"))
cat(sprintf("  Summary: homer_motif_summary_stats.tsv\n"))
cat("================================================================================\n")
