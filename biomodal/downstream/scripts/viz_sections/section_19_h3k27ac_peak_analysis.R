# biomodal/downstream/scripts/viz_sections/section_19_h3k27ac_peak_analysis.R
# Section 19: H3K27ac Peak Analysis at DMR Genes
#
# Comprehensive analysis of condition-specific H3K27ac peaks at DMR genes,
# parallel to Sections 16-18 (K119ub). Uses ctrl/mut intersect peak files
# rather than continuous bigwig signal. Eight figures:
#   19a: Complete H3K27ac status breakdown (100% stacked bars)
#   19b: Conditional direction among genes with peaks (grouped bars + Fisher)
#   19c: Per-gene H3K27ac effect size distribution (strip + box)
#   19d: Methylation change vs H3K27ac peak change scatter (Spearman)
#   19e: Gene waterfall ranked by H3K27ac net peak change
#   19f: Condition-specific H3K27ac overlap at DMRs (grouped bars + Fisher)
#   19g: 2x2 O/E heatmaps (mC x H3K27ac, hmC x H3K27ac)
#   19h: Comparative O/E dot plot across 4 chromatin marks
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_19_h3k27ac_peak_analysis.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages(library(ChIPseeker))

# =============================================================================
# HELPER FUNCTIONS (self-contained for section independence)
# =============================================================================

annotate_peaks_to_genes <- function(gr, txdb) {
  anno <- suppressMessages(annotatePeak(
    gr, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"
  ))
  as.data.frame(anno)
}

aggregate_peaks_by_gene <- function(anno_df, label) {
  anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(!!paste0("n_", label) := n(), .groups = "drop")
}

classify_h3k27ac_status <- function(gene_list, h3k27ac_df) {
  gene_df <- data.frame(gene = gene_list, stringsAsFactors = FALSE) %>%
    left_join(h3k27ac_df, by = c("gene" = "SYMBOL")) %>%
    mutate(
      n_h3k27ac_ctrl = replace_na(n_h3k27ac_ctrl, 0),
      n_h3k27ac_mut = replace_na(n_h3k27ac_mut, 0),
      net_h3k27ac = n_h3k27ac_mut - n_h3k27ac_ctrl,
      h3k27ac_status = case_when(
        n_h3k27ac_ctrl == 0 & n_h3k27ac_mut == 0 ~ "No Peaks",
        net_h3k27ac > 0 ~ "Gained",
        net_h3k27ac < 0 ~ "Lost",
        TRUE ~ "Equal Peaks"
      )
    )
  gene_df
}

# Cliff's delta (non-parametric effect size)
cliffs_delta <- function(x, y) {
  n_x <- length(x)
  n_y <- length(y)
  more <- sum(outer(x, y, ">"))
  less <- sum(outer(x, y, "<"))
  (more - less) / (n_x * n_y)
}

# Load MeCP2 annotated data from HOMER-style file (for 19h)
load_mecp2_annotated <- function(filepath) {
  stopifnot(file.exists(filepath))
  df <- read.table(filepath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE, quote = "",
                   comment.char = "")
  numeric_cols <- c("Fold", "FDR", "p.value", "distanceToTSS",
                    "geneStart", "geneEnd", "geneLength")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) df[[col]] <- as.numeric(df[[col]])
  }
  df
}

# Aggregate MeCP2 peaks per gene — nearest peak to TSS (for 19h)
aggregate_mecp2_by_gene <- function(mecp2_df) {
  mecp2_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      mean_fold = mean(Fold, na.rm = TRUE),
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      nearest_fdr = FDR[which.min(abs(distanceToTSS))],
      n_peaks = n(),
      .groups = "drop"
    )
}

# Derive differential peaks from condition-specific peak sets (for 19h)
derive_differential_peaks <- function(ctrl_gr, mut_gr) {
  ctrl_in_mut <- countOverlaps(ctrl_gr, mut_gr) > 0
  mut_in_ctrl <- countOverlaps(mut_gr, ctrl_gr) > 0
  list(
    gained = mut_gr[!mut_in_ctrl],
    lost = ctrl_gr[!ctrl_in_mut],
    shared = ctrl_gr[ctrl_in_mut]
  )
}

# Build 2x2 O/E heatmap for methylation direction x chromatin mark direction
# Adapted from section 15. predicted_pairs = list() disables border highlighting.
build_2x2_heatmap <- function(met_direction, mark_direction,
                              met_levels, mark_levels,
                              predicted_pairs = list(),
                              title, subtitle_extra = "",
                              x_label, y_label) {
  quad_table <- table(
    factor(met_direction, levels = met_levels),
    factor(mark_direction, levels = mark_levels)
  )
  fisher_result <- fisher.test(quad_table)
  total_n <- sum(quad_table)
  row_totals <- rowSums(quad_table)
  col_totals <- colSums(quad_table)
  expected <- outer(row_totals, col_totals) / total_n
  enrichment <- quad_table / expected

  heatmap_df <- as.data.frame(as.table(quad_table))
  colnames(heatmap_df) <- c("met_dir", "mark_dir", "Observed")
  expected_df <- as.data.frame(as.table(round(expected, 1)))
  colnames(expected_df) <- c("met_dir", "mark_dir", "Expected")
  enrichment_df <- as.data.frame(as.table(round(enrichment, 2)))
  colnames(enrichment_df) <- c("met_dir", "mark_dir", "Enrichment")

  heatmap_data <- heatmap_df %>%
    left_join(expected_df, by = c("met_dir", "mark_dir")) %>%
    left_join(enrichment_df, by = c("met_dir", "mark_dir")) %>%
    mutate(
      label = sprintf("Obs: %d\nExp: %.0f\nO/E: %.2f", Observed, Expected, Enrichment),
      met_dir = factor(met_dir, levels = met_levels),
      mark_dir = factor(mark_dir, levels = mark_levels),
      is_predicted = FALSE
    )

  for (pair in predicted_pairs) {
    heatmap_data$is_predicted[heatmap_data$met_dir == pair[1] &
                              heatmap_data$mark_dir == pair[2]] <- TRUE
  }

  if (length(predicted_pairs) > 0) {
    sub_text <- sprintf("Fisher OR = %.2f, p = %.2e | Black borders = predicted%s",
                        fisher_result$estimate, fisher_result$p.value, subtitle_extra)
  } else {
    sub_text <- sprintf("Fisher OR = %.2f, p = %.2e%s",
                        fisher_result$estimate, fisher_result$p.value, subtitle_extra)
  }

  p <- ggplot(heatmap_data, aes(x = mark_dir, y = met_dir, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2)

  if (any(heatmap_data$is_predicted)) {
    p <- p + geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
                       aes(x = mark_dir, y = met_dir),
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

# Extract O/E for enriched quadrants (O/E > 1) from a build_2x2_heatmap result
extract_enriched_oe <- function(result, perspective, mark_name) {
  if (is.null(result)) return(NULL)
  enriched <- result$data %>% dplyr::filter(Enrichment > 1)
  if (nrow(enriched) == 0) return(NULL)
  data.frame(
    perspective = perspective,
    mark = mark_name,
    quadrant = paste(as.character(enriched$met_dir), "+",
                     as.character(enriched$mark_dir)),
    observed = as.integer(enriched$Observed),
    expected = as.numeric(enriched$Expected),
    oe_ratio = as.numeric(enriched$Enrichment),
    fisher_or = result$fisher$estimate,
    fisher_p = result$fisher$p.value,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# =============================================================================
# SECTION 19: H3K27AC PEAK ANALYSIS AT DMR GENES
# =============================================================================

cat("================================================================================\n")
cat("SECTION 19: H3K27AC PEAK ANALYSIS AT DMR GENES\n")
cat("================================================================================\n\n")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# =============================================================================
# DATA PREPARATION
# =============================================================================

# -----------------------------------------------------------------------
# Derive methylation gene groups
# -----------------------------------------------------------------------
cat("Preparing methylation gene groups...\n")

mc_sig <- mc_dmr %>% dplyr::filter(significant)
hmc_sig <- hmc_dmr %>% dplyr::filter(significant)

mc_up_genes   <- mc_sig$gene[mc_sig$mod_difference > 0]
mc_down_genes <- mc_sig$gene[mc_sig$mod_difference < 0]
hmc_down_genes <- hmc_sig$gene[hmc_sig$mod_difference < 0]
hmc_up_genes   <- hmc_sig$gene[hmc_sig$mod_difference > 0]

all_dmr_genes <- unique(c(mc_up_genes, mc_down_genes, hmc_down_genes, hmc_up_genes))

cat(sprintf("  mC Up:    %d genes\n", length(mc_up_genes)))
cat(sprintf("  mC Down:  %d genes\n", length(mc_down_genes)))
cat(sprintf("  hmC Down: %d genes\n", length(hmc_down_genes)))
cat(sprintf("  hmC Up:   %d genes\n", length(hmc_up_genes)))
cat(sprintf("  All DMR:  %d unique genes\n", length(all_dmr_genes)))

# -----------------------------------------------------------------------
# Load H3K27ac peaks and annotate to genes
# -----------------------------------------------------------------------
cat("\nLoading H3K27ac condition-specific peaks...\n")

h3k27ac_ctrl_gr <- load_chip_peaks(H3K27AC_FILES$ctrl, "H3K27ac Control")
h3k27ac_mut_gr  <- load_chip_peaks(H3K27AC_FILES$mut, "H3K27ac Mutant")

if (is.null(h3k27ac_ctrl_gr) || is.null(h3k27ac_mut_gr)) {
  stop("H3K27ac peak files not found. Cannot proceed with Section 19.")
}

cat("  Annotating H3K27ac peaks to genes via ChIPseeker...\n")
h3k27ac_ctrl_gene <- aggregate_peaks_by_gene(
  annotate_peaks_to_genes(h3k27ac_ctrl_gr, txdb), "h3k27ac_ctrl")
h3k27ac_mut_gene <- aggregate_peaks_by_gene(
  annotate_peaks_to_genes(h3k27ac_mut_gr, txdb), "h3k27ac_mut")

# Build per-gene H3K27ac summary
h3k27ac_gene <- full_join(h3k27ac_ctrl_gene, h3k27ac_mut_gene, by = "SYMBOL") %>%
  mutate(
    n_h3k27ac_ctrl = replace_na(n_h3k27ac_ctrl, 0),
    n_h3k27ac_mut = replace_na(n_h3k27ac_mut, 0),
    net_h3k27ac = n_h3k27ac_mut - n_h3k27ac_ctrl
  )

cat(sprintf("  %d unique genes with any H3K27ac peak\n", nrow(h3k27ac_gene)))
cat(sprintf("  Ctrl-only genes: %d\n", sum(h3k27ac_gene$n_h3k27ac_ctrl > 0 & h3k27ac_gene$n_h3k27ac_mut == 0)))
cat(sprintf("  Mut-only genes:  %d\n", sum(h3k27ac_gene$n_h3k27ac_ctrl == 0 & h3k27ac_gene$n_h3k27ac_mut > 0)))
cat(sprintf("  Both conditions: %d\n", sum(h3k27ac_gene$n_h3k27ac_ctrl > 0 & h3k27ac_gene$n_h3k27ac_mut > 0)))

# -----------------------------------------------------------------------
# Classify every DMR gene into 4 categories
# -----------------------------------------------------------------------
cat("\nClassifying H3K27ac status for all DMR genes...\n")

mc_up_classified   <- classify_h3k27ac_status(mc_up_genes, h3k27ac_gene) %>%
  mutate(met_group = "mC Up")
mc_down_classified <- classify_h3k27ac_status(mc_down_genes, h3k27ac_gene) %>%
  mutate(met_group = "mC Down")
hmc_down_classified <- classify_h3k27ac_status(hmc_down_genes, h3k27ac_gene) %>%
  mutate(met_group = "hmC Down")
hmc_up_classified  <- classify_h3k27ac_status(hmc_up_genes, h3k27ac_gene) %>%
  mutate(met_group = "hmC Up")

all_classified <- bind_rows(mc_up_classified, mc_down_classified,
                            hmc_down_classified, hmc_up_classified)

# Also classify ALL DMR genes for background rate
all_dmr_classified <- classify_h3k27ac_status(all_dmr_genes, h3k27ac_gene) %>%
  mutate(met_group = "All DMR Genes")

# Print overview
cat("\n  H3K27ac status breakdown:\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  sub <- all_classified %>% dplyr::filter(met_group == grp)
  total <- nrow(sub)
  no_peaks <- sum(sub$h3k27ac_status == "No Peaks")
  with_peaks <- total - no_peaks
  cat(sprintf("    %-10s: %d total, %d (%.1f%%) with peaks, %d (%.1f%%) no peaks\n",
              grp, total, with_peaks, 100 * with_peaks / total,
              no_peaks, 100 * no_peaks / total))
}

# =============================================================================
# FIGURE 19a: Complete H3K27ac Status Breakdown (100% Stacked Bars)
# =============================================================================

cat("\n--- FIGURE 19a: Complete H3K27ac Status Breakdown ---\n\n")

compute_status_summary <- function(classified_df) {
  classified_df %>%
    group_by(met_group, h3k27ac_status) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(met_group) %>%
    mutate(
      total = sum(n),
      pct = 100 * n / total
    ) %>%
    ungroup()
}

group_summary <- compute_status_summary(all_classified)
background_summary <- compute_status_summary(all_dmr_classified)

status_order <- c("Gained", "Equal Peaks", "No Peaks", "Lost")
status_colors_19 <- c(
  "Gained"      = "#FF7F00",  # orange
  "Equal Peaks" = "#BDBDBD",  # light grey
  "No Peaks"    = "#636363",  # dark grey
  "Lost"        = "#1F78B4"   # blue
)

group_summary$met_group <- factor(group_summary$met_group,
                                   levels = c("mC Up", "mC Down", "hmC Down", "hmC Up"))
group_summary$h3k27ac_status <- factor(group_summary$h3k27ac_status, levels = status_order)

group_summary$label <- sprintf("n=%d\n(%.0f%%)", group_summary$n, group_summary$pct)

p_19a <- ggplot(group_summary,
                aes(x = met_group, y = pct, fill = h3k27ac_status)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2.8, lineheight = 0.9, color = "white", fontface = "bold") +
  scale_fill_manual(values = status_colors_19, name = "H3K27ac Status",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "H3K27ac Status Across ALL DMR Genes",
    subtitle = "Including genes with no H3K27ac peaks in either condition",
    x = "Methylation Group (Significant DMR Genes)",
    y = "Percentage of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_19a,
                        file.path(OUTPUT_DIR, "19a_h3k27ac_status_breakdown"),
                        width = 9, height = 7)

cat("  Per-group breakdown:\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  sub <- group_summary %>%
    dplyr::filter(met_group == grp) %>%
    arrange(match(h3k27ac_status, status_order))
  cat(sprintf("    %s (n=%d):\n", grp, sub$total[1]))
  for (i in 1:nrow(sub)) {
    cat(sprintf("      %-12s: %4d (%5.1f%%)\n",
                as.character(sub$h3k27ac_status[i]), sub$n[i], sub$pct[i]))
  }
}

# =============================================================================
# FIGURE 19b: H3K27ac Direction Among Genes WITH Peaks
# =============================================================================

cat("\n--- FIGURE 19b: H3K27ac Conditional Direction ---\n\n")

conditional_summary <- all_classified %>%
  dplyr::filter(h3k27ac_status != "No Peaks") %>%
  group_by(met_group) %>%
  summarise(
    n_with_peaks = n(),
    n_gained = sum(h3k27ac_status == "Gained"),
    n_equal = sum(h3k27ac_status == "Equal Peaks"),
    n_lost = sum(h3k27ac_status == "Lost"),
    pct_gained = 100 * n_gained / n(),
    pct_equal = 100 * n_equal / n(),
    pct_lost = 100 * n_lost / n(),
    .groups = "drop"
  )

bg <- all_dmr_classified %>%
  dplyr::filter(h3k27ac_status != "No Peaks") %>%
  summarise(
    met_group = "All DMR Genes",
    n_with_peaks = n(),
    n_gained = sum(h3k27ac_status == "Gained"),
    n_equal = sum(h3k27ac_status == "Equal Peaks"),
    n_lost = sum(h3k27ac_status == "Lost"),
    pct_gained = 100 * n_gained / n(),
    pct_equal = 100 * n_equal / n(),
    pct_lost = 100 * n_lost / n()
  )

conditional_summary <- bind_rows(conditional_summary, bg)

# Pivot for grouped bar chart
conditional_long <- conditional_summary %>%
  dplyr::select(met_group, pct_gained, pct_equal, pct_lost) %>%
  pivot_longer(cols = c(pct_gained, pct_equal, pct_lost),
               names_to = "direction", values_to = "pct") %>%
  mutate(
    direction = recode(direction,
                       "pct_gained" = "Gained",
                       "pct_equal" = "Equal Peaks",
                       "pct_lost" = "Lost"),
    direction = factor(direction, levels = c("Gained", "Equal Peaks", "Lost"))
  )

# Merge counts for annotation
conditional_long <- conditional_long %>%
  left_join(
    conditional_summary %>%
      dplyr::select(met_group, n_with_peaks, n_gained, n_equal, n_lost) %>%
      pivot_longer(cols = c(n_gained, n_equal, n_lost),
                   names_to = "dir_count", values_to = "n_count") %>%
      mutate(direction = recode(dir_count,
                                "n_gained" = "Gained",
                                "n_equal" = "Equal Peaks",
                                "n_lost" = "Lost")) %>%
      dplyr::select(met_group, direction, n_count),
    by = c("met_group", "direction")
  )

conditional_long$met_group <- factor(
  conditional_long$met_group,
  levels = c("mC Up", "mC Down", "hmC Down", "hmC Up", "All DMR Genes")
)

direction_colors_19 <- c(
  "Gained"      = "#FF7F00",
  "Equal Peaks" = "#BDBDBD",
  "Lost"        = "#1F78B4"
)

bg_gain_rate <- bg$pct_gained
bg_loss_rate <- bg$pct_lost

# Fisher exact tests: each group's gain rate vs background
fisher_results <- list()
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  grp_row <- conditional_summary %>% dplyr::filter(met_group == grp)
  bg_row <- conditional_summary %>% dplyr::filter(met_group == "All DMR Genes")

  mat <- matrix(c(
    grp_row$n_gained, grp_row$n_with_peaks - grp_row$n_gained,
    bg_row$n_gained, bg_row$n_with_peaks - bg_row$n_gained
  ), nrow = 2, byrow = TRUE)

  ft <- fisher.test(mat)
  fisher_results[[grp]] <- ft
}

p_19b <- ggplot(conditional_long,
                aes(x = met_group, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d)", pct, n_count)),
            position = position_dodge(width = 0.75),
            vjust = -0.2, size = 2.6, lineheight = 0.85) +
  # Background gain rate reference line
  geom_hline(yintercept = bg_gain_rate, linetype = "dashed",
             color = "#FF7F00", linewidth = 0.5) +
  annotate("text", x = 0.55, y = bg_gain_rate + 1.5,
           label = sprintf("BG gain: %.1f%%", bg_gain_rate),
           size = 2.5, color = "#FF7F00", hjust = 0, fontface = "italic") +
  # Background loss rate reference line
  geom_hline(yintercept = bg_loss_rate, linetype = "dashed",
             color = "#1F78B4", linewidth = 0.5) +
  annotate("text", x = 0.55, y = bg_loss_rate + 1.5,
           label = sprintf("BG loss: %.1f%%", bg_loss_rate),
           size = 2.5, color = "#1F78B4", hjust = 0, fontface = "italic") +
  # Fisher p-value annotations
  annotate("text",
           x = 1:4,
           y = rep(max(conditional_long$pct, na.rm = TRUE) * 1.12, 4),
           label = sapply(c("mC Up", "mC Down", "hmC Down", "hmC Up"), function(g) {
             p <- fisher_results[[g]]$p.value
             if (p < 0.001) sprintf("p=%.1e", p)
             else sprintf("p=%.3f", p)
           }),
           size = 2.5, fontface = "italic", color = "grey30") +
  scale_fill_manual(values = direction_colors_19, name = "H3K27ac Direction") +
  scale_y_continuous(limits = c(0, max(conditional_long$pct, na.rm = TRUE) * 1.25),
                     expand = c(0, 0)) +
  labs(
    title = "H3K27ac Direction Among Genes WITH Peaks Only",
    subtitle = "Dashed lines = background gain/loss rate across all DMR genes with peaks",
    x = "Methylation Group",
    y = "% of Genes with H3K27ac Peaks"
  ) +
  theme_biomodal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 0.7))

save_multiformat_ggplot(p_19b,
                        file.path(OUTPUT_DIR, "19b_h3k27ac_conditional_direction"),
                        width = 10, height = 7)

cat("  Conditional direction (genes with peaks):\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up", "All DMR Genes")) {
  row <- conditional_summary %>% dplyr::filter(met_group == grp)
  p_str <- ""
  if (grp %in% names(fisher_results)) {
    p <- fisher_results[[grp]]$p.value
    p_str <- sprintf(", Fisher vs BG p=%.2e", p)
  }
  cat(sprintf("    %-14s: %d genes with peaks, %.1f%% gained, %.1f%% equal, %.1f%% lost%s\n",
              grp, row$n_with_peaks, row$pct_gained, row$pct_equal, row$pct_lost, p_str))
}

cat(sprintf("\n  Background gain rate: %.1f%%\n", bg_gain_rate))
cat(sprintf("  Background loss rate: %.1f%%\n", bg_loss_rate))

# =============================================================================
# FIGURE 19c: Per-Gene H3K27ac Effect Size Distribution
# =============================================================================

cat("\n--- FIGURE 19c: Per-Gene H3K27ac Effect Size ---\n\n")

effect_df <- all_classified %>%
  dplyr::filter(h3k27ac_status != "No Peaks") %>%
  dplyr::select(gene, met_group, net_h3k27ac, n_h3k27ac_ctrl, n_h3k27ac_mut)

effect_df$met_group <- factor(effect_df$met_group,
                               levels = c("mC Up", "mC Down", "hmC Down", "hmC Up"))

# Compute medians for annotation
medians <- effect_df %>%
  group_by(met_group) %>%
  summarise(
    median_net = median(net_h3k27ac),
    iqr_low = quantile(net_h3k27ac, 0.25),
    iqr_high = quantile(net_h3k27ac, 0.75),
    n = n(),
    .groups = "drop"
  )

# Identify extreme genes for labeling
extreme_genes <- effect_df %>%
  group_by(met_group) %>%
  slice_max(abs(net_h3k27ac), n = 3, with_ties = FALSE) %>%
  dplyr::filter(abs(net_h3k27ac) >= 5) %>%
  ungroup()

# Color by methylation type
effect_df$met_type <- ifelse(grepl("^mC", effect_df$met_group), "5mC", "5hmC")

p_19c <- ggplot(effect_df, aes(x = met_group, y = net_h3k27ac, color = met_type)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_jitter(width = 0.25, alpha = 0.3, size = 1, shape = 16) +
  geom_boxplot(aes(fill = met_type), width = 0.4, alpha = 0.15,
               outlier.shape = NA, color = "black", linewidth = 0.4) +
  # Annotate median + IQR
  geom_text(data = medians,
            aes(x = met_group, y = iqr_high + 2,
                label = sprintf("med=%+.0f\nn=%d", median_net, n)),
            inherit.aes = FALSE,
            size = 2.8, lineheight = 0.85, color = "grey20") +
  # Label extreme genes
  {if (nrow(extreme_genes) > 0)
    ggrepel::geom_text_repel(
      data = extreme_genes,
      aes(label = gene),
      size = 2.3, fontface = "italic", color = "grey30",
      max.overlaps = 10, nudge_x = 0.3,
      segment.size = 0.2, segment.color = "grey60"
    )
  } +
  scale_color_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8"),
                     guide = "none") +
  scale_fill_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8"),
                    guide = "none") +
  labs(
    title = "Per-Gene H3K27ac Change (Mutant - Control Peak Count)",
    subtitle = "Only genes with at least one H3K27ac peak in either condition",
    x = "Methylation Group",
    y = "Net H3K27ac Peaks (Mutant - Control)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_19c,
                        file.path(OUTPUT_DIR, "19c_h3k27ac_effect_sizes"),
                        width = 9, height = 7)

cat("  Effect size distributions (genes with peaks):\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  sub <- effect_df %>% dplyr::filter(met_group == grp)
  cat(sprintf("    %-10s: n=%d, median=%+.0f, IQR=[%+.0f, %+.0f], range=[%+.0f, %+.0f]\n",
              grp, nrow(sub),
              median(sub$net_h3k27ac),
              quantile(sub$net_h3k27ac, 0.25),
              quantile(sub$net_h3k27ac, 0.75),
              min(sub$net_h3k27ac),
              max(sub$net_h3k27ac)))
}

if (nrow(extreme_genes) > 0) {
  cat("\n  Extreme genes (|net| >= 5):\n")
  for (i in 1:nrow(extreme_genes)) {
    cat(sprintf("    %s (%s): net = %+d (ctrl=%d, mut=%d)\n",
                extreme_genes$gene[i], extreme_genes$met_group[i],
                extreme_genes$net_h3k27ac[i],
                extreme_genes$n_h3k27ac_ctrl[i],
                extreme_genes$n_h3k27ac_mut[i]))
  }
}

# =============================================================================
# FIGURE 19d: Methylation Change vs H3K27ac Peak Change Scatter
# =============================================================================

cat("\n--- FIGURE 19d: Methylation vs H3K27ac Scatter ---\n\n")

# Build scatter data: merge DMR mod_difference with net H3K27ac peak change
# Only genes with at least one peak in either condition
build_scatter <- function(dmr_df, modality_label, h3k27ac_df) {
  dmr_df %>%
    dplyr::filter(significant) %>%
    dplyr::select(gene, mod_difference) %>%
    inner_join(
      h3k27ac_df %>%
        dplyr::select(SYMBOL, net_h3k27ac),
      by = c("gene" = "SYMBOL")
    ) %>%
    mutate(modality = modality_label)
}

scatter_mc  <- build_scatter(mc_dmr, "5mC", h3k27ac_gene)
scatter_hmc <- build_scatter(hmc_dmr, "5hmC", h3k27ac_gene)
scatter_all <- bind_rows(scatter_mc, scatter_hmc)
scatter_all$modality <- factor(scatter_all$modality, levels = c("5mC", "5hmC"))

cat(sprintf("  5mC genes with H3K27ac peaks: %d\n", nrow(scatter_mc)))
cat(sprintf("  5hmC genes with H3K27ac peaks: %d\n", nrow(scatter_hmc)))

# Spearman correlations per modality
cor_results <- scatter_all %>%
  group_by(modality) %>%
  summarise(
    n = n(),
    rho = cor(mod_difference, net_h3k27ac, method = "spearman"),
    p_value = cor.test(mod_difference, net_h3k27ac, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("rho = %+.3f\np = %s\nn = %d",
                          rho,
                          ifelse(p_value < 0.001, sprintf("%.1e", p_value),
                                 sprintf("%.3f", p_value)),
                          n))

# Label KEY_GENES
scatter_all$is_key <- scatter_all$gene %in% KEY_GENES

p_19d <- ggplot(scatter_all, aes(x = mod_difference, y = net_h3k27ac)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.15, size = 0.8, color = "grey40") +
  geom_smooth(method = "lm", color = "#FF7F00", fill = "#FF7F00",
              alpha = 0.15, linewidth = 0.8) +
  # Key gene labels
  {if (any(scatter_all$is_key))
    ggrepel::geom_text_repel(
      data = scatter_all %>% dplyr::filter(is_key),
      aes(label = gene),
      size = 2.5, fontface = "italic", color = "black",
      max.overlaps = 15, nudge_y = 0.5,
      segment.size = 0.2, segment.color = "grey50"
    )
  } +
  # Correlation annotation
  geom_text(data = cor_results,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.3, size = 3.2, color = "black",
            lineheight = 0.9) +
  facet_wrap(~ modality, scales = "free_x") +
  labs(
    title = "Methylation Change vs H3K27ac Peak Change",
    subtitle = "Gene body mod_difference vs net H3K27ac peaks (mut - ctrl); Spearman correlation",
    x = "Methylation Difference (Mutant - Control)",
    y = "Net H3K27ac Peaks (Mutant - Control)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_19d,
                        file.path(OUTPUT_DIR, "19d_methylation_vs_h3k27ac_scatter"),
                        width = 14, height = 7)

cat("  Spearman correlations:\n")
for (i in 1:nrow(cor_results)) {
  cat(sprintf("    %s: rho = %+.4f, p = %s, n = %d\n",
              cor_results$modality[i], cor_results$rho[i],
              ifelse(cor_results$p_value[i] < 0.001,
                     sprintf("%.2e", cor_results$p_value[i]),
                     sprintf("%.4f", cor_results$p_value[i])),
              cor_results$n[i]))
}

# =============================================================================
# FIGURE 19e: Gene Waterfall Ranked by H3K27ac Change
# =============================================================================

cat("\n--- FIGURE 19e: H3K27ac Waterfall ---\n\n")

# All genes with any H3K27ac peak, ranked by net peak change
waterfall_base <- h3k27ac_gene %>%
  dplyr::select(SYMBOL, net_h3k27ac) %>%
  arrange(desc(net_h3k27ac)) %>%
  mutate(rank = row_number())

# Tag DMR group membership (priority: mC Up > hmC Down > mC Down > hmC Up)
waterfall_base$dmr_group <- "None"
waterfall_base$dmr_group[waterfall_base$SYMBOL %in% hmc_up_genes]   <- "hmC Up"
waterfall_base$dmr_group[waterfall_base$SYMBOL %in% mc_down_genes]  <- "mC Down"
waterfall_base$dmr_group[waterfall_base$SYMBOL %in% hmc_down_genes] <- "hmC Down"
waterfall_base$dmr_group[waterfall_base$SYMBOL %in% mc_up_genes]    <- "mC Up"

waterfall_base$dmr_group <- factor(waterfall_base$dmr_group,
  levels = c("mC Up", "mC Down", "hmC Down", "hmC Up", "None"))

waterfall_colors <- c(
  "mC Up"    = "#D7191C",
  "mC Down"  = "#2C7BB6",
  "hmC Down" = "#377EB8",
  "hmC Up"   = "#E41A1C",
  "None"     = "grey85"
)

# Compute where DMR genes cluster
n_total <- nrow(waterfall_base)
cat(sprintf("  Total genes with H3K27ac peaks: %d\n", n_total))
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  ranks <- waterfall_base$rank[waterfall_base$dmr_group == grp]
  if (length(ranks) > 0) {
    cat(sprintf("  %s: n=%d, median rank=%.0f / %d (percentile=%.1f%%)\n",
                grp, length(ranks), median(ranks), n_total,
                100 * median(ranks) / n_total))
  }
}

# Plot: draw "None" first (background), then DMR groups on top
waterfall_dmr <- waterfall_base %>% dplyr::filter(dmr_group != "None")
waterfall_none <- waterfall_base %>% dplyr::filter(dmr_group == "None")

p_19e <- ggplot() +
  # Background genes
  geom_point(data = waterfall_none,
             aes(x = rank, y = net_h3k27ac),
             color = "grey85", size = 0.3, alpha = 0.5) +
  # DMR genes on top
  geom_point(data = waterfall_dmr,
             aes(x = rank, y = net_h3k27ac, color = dmr_group),
             size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = waterfall_colors, name = "DMR Group",
                     breaks = c("mC Up", "mC Down", "hmC Down", "hmC Up")) +
  labs(
    title = "Gene Waterfall: All Genes Ranked by H3K27ac Peak Change",
    subtitle = "Grey = non-DMR genes; colored = significant DMR genes",
    x = sprintf("Gene Rank (1 = highest H3K27ac gain, %d = highest loss)", n_total),
    y = "Net H3K27ac Peaks (Mutant - Control)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_19e,
                        file.path(OUTPUT_DIR, "19e_h3k27ac_waterfall"),
                        width = 12, height = 6)

# =============================================================================
# FIGURE 19f: Condition-Specific H3K27ac Overlap at DMRs
# =============================================================================

cat("\n--- FIGURE 19f: Condition-Specific H3K27ac Overlap at DMRs ---\n\n")

# Get significant mC DMRs split by direction
mc_hyper <- mc_dmr %>% dplyr::filter(significant & mod_difference > 0)
mc_hypo  <- mc_dmr %>% dplyr::filter(significant & mod_difference < 0)

hyper_gr <- dmr_to_granges(mc_hyper)
hypo_gr  <- dmr_to_granges(mc_hypo)
n_hyper <- length(hyper_gr)
n_hypo  <- length(hypo_gr)

cat(sprintf("  Hypermethylated DMRs: %d\n", n_hyper))
cat(sprintf("  Hypomethylated DMRs: %d\n", n_hypo))

# Count overlaps: % of hyper/hypo DMRs overlapping H3K27ac peaks in ctrl vs mut
hyper_ctrl_ov <- sum(countOverlaps(hyper_gr, h3k27ac_ctrl_gr) > 0)
hyper_mut_ov  <- sum(countOverlaps(hyper_gr, h3k27ac_mut_gr) > 0)
hypo_ctrl_ov  <- sum(countOverlaps(hypo_gr, h3k27ac_ctrl_gr) > 0)
hypo_mut_ov   <- sum(countOverlaps(hypo_gr, h3k27ac_mut_gr) > 0)

cat(sprintf("  Hypermethylated DMRs: %.1f%% ctrl H3K27ac, %.1f%% mut H3K27ac\n",
            100 * hyper_ctrl_ov / n_hyper, 100 * hyper_mut_ov / n_hyper))
cat(sprintf("  Hypomethylated DMRs: %.1f%% ctrl H3K27ac, %.1f%% mut H3K27ac\n",
            100 * hypo_ctrl_ov / n_hypo, 100 * hypo_mut_ov / n_hypo))

# Fisher's exact test on 2x2 contingency
fisher_19f_mat <- matrix(c(hyper_ctrl_ov, hyper_mut_ov, hypo_ctrl_ov, hypo_mut_ov),
                         nrow = 2, byrow = TRUE,
                         dimnames = list(c("Hypermethylated", "Hypomethylated"),
                                         c("Control H3K27ac", "Mutant H3K27ac")))
fisher_19f <- fisher.test(fisher_19f_mat)

cat(sprintf("  Fisher's exact test: OR = %.2f, p = %.2e\n",
            fisher_19f$estimate, fisher_19f$p.value))

overlap_19f_df <- data.frame(
  DMR_Direction = rep(c("Hypermethylated", "Hypomethylated"), each = 2),
  Genotype = rep(c("Control", "Mutant"), 2),
  Count = c(hyper_ctrl_ov, hyper_mut_ov, hypo_ctrl_ov, hypo_mut_ov),
  Total = rep(c(n_hyper, n_hypo), each = 2)
) %>%
  mutate(Percentage = 100 * Count / Total)

overlap_19f_df$DMR_Direction <- factor(overlap_19f_df$DMR_Direction,
                                        levels = c("Hypermethylated", "Hypomethylated"))
overlap_19f_df$Genotype <- factor(overlap_19f_df$Genotype, levels = c("Control", "Mutant"))

p_19f <- ggplot(overlap_19f_df, aes(x = DMR_Direction, y = Percentage, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, Count)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = COLORS$condition, name = "H3K27ac Peaks") +
  scale_y_continuous(limits = c(0, max(overlap_19f_df$Percentage) * 1.3), expand = c(0, 0)) +
  labs(
    title = "H3K27ac Peak Overlap at Differentially Methylated Regions",
    subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e",
                       fisher_19f$estimate, fisher_19f$p.value),
    x = "DMR Direction (5mC)", y = "% of DMRs Overlapping H3K27ac Peaks"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_19f, file.path(OUTPUT_DIR, "19f_h3k27ac_condition_overlap"),
                        width = 8, height = 7)

# =============================================================================
# FIGURE 19g: 2x2 O/E Heatmaps (Methylation x H3K27ac Direction)
# =============================================================================

cat("\n--- FIGURE 19g: H3K27ac O/E Heatmaps ---\n\n")

# mC Direction x H3K27ac Direction (gene-level, net_h3k27ac != 0)
mc_h3k27ac_oe <- mc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(h3k27ac_gene %>% dplyr::select(SYMBOL, net_h3k27ac),
            by = c("gene" = "SYMBOL")) %>%
  mutate(net_h3k27ac = replace_na(net_h3k27ac, 0)) %>%
  dplyr::filter(net_h3k27ac != 0) %>%
  mutate(
    mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
    h3k27ac_direction = ifelse(net_h3k27ac > 0, "H3K27ac Gained", "H3K27ac Lost")
  )

cat(sprintf("  mC DMR genes with non-zero net H3K27ac: %d\n", nrow(mc_h3k27ac_oe)))

result_19g_mc <- build_2x2_heatmap(
  met_direction = mc_h3k27ac_oe$mc_direction,
  mark_direction = mc_h3k27ac_oe$h3k27ac_direction,
  met_levels = c("mC Up", "mC Down"),
  mark_levels = c("H3K27ac Gained", "H3K27ac Lost"),
  predicted_pairs = list(),
  title = "5mC Direction x H3K27ac Direction",
  x_label = "H3K27ac Direction", y_label = "5mC DMR Direction"
)

cat("  mC x H3K27ac contingency table:\n")
print(result_19g_mc$table)
cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
            result_19g_mc$fisher$estimate, result_19g_mc$fisher$p.value))

# hmC Direction x H3K27ac Direction
hmc_h3k27ac_oe <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(h3k27ac_gene %>% dplyr::select(SYMBOL, net_h3k27ac),
            by = c("gene" = "SYMBOL")) %>%
  mutate(net_h3k27ac = replace_na(net_h3k27ac, 0)) %>%
  dplyr::filter(net_h3k27ac != 0) %>%
  mutate(
    hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
    h3k27ac_direction = ifelse(net_h3k27ac > 0, "H3K27ac Gained", "H3K27ac Lost")
  )

cat(sprintf("  hmC DMR genes with non-zero net H3K27ac: %d\n", nrow(hmc_h3k27ac_oe)))

result_19g_hmc <- build_2x2_heatmap(
  met_direction = hmc_h3k27ac_oe$hmc_direction,
  mark_direction = hmc_h3k27ac_oe$h3k27ac_direction,
  met_levels = c("hmC Down", "hmC Up"),
  mark_levels = c("H3K27ac Gained", "H3K27ac Lost"),
  predicted_pairs = list(),
  title = "5hmC Direction x H3K27ac Direction",
  x_label = "H3K27ac Direction", y_label = "5hmC DMR Direction"
)

cat("  hmC x H3K27ac contingency table:\n")
print(result_19g_hmc$table)
cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
            result_19g_hmc$fisher$estimate, result_19g_hmc$fisher$p.value))

# Combine via patchwork
p_19g <- result_19g_mc$plot + result_19g_hmc$plot +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "O/E Enrichment: Methylation Direction x H3K27ac Direction",
    subtitle = "No pre-specified predictions for H3K27ac",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
    )
  )

save_multiformat_ggplot(p_19g, file.path(OUTPUT_DIR, "19g_h3k27ac_oe_heatmaps"),
                        width = 16, height = 8)

# =============================================================================
# FIGURE 19h: Comparative O/E Dot Plot Across 4 Chromatin Marks
# =============================================================================

cat("\n--- FIGURE 19h: Comparative O/E Dot Plot ---\n\n")

# --- Re-compute MeCP2 O/E (independence from section 15) ---
cat("  Re-computing MeCP2 O/E...\n")
mecp2_annotated_19h <- load_mecp2_annotated(MECP2_FILES$annotated)
mecp2_gene_19h <- aggregate_mecp2_by_gene(mecp2_annotated_19h)

# mC x MeCP2
mc_mecp2_19h <- mc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(mecp2_gene_19h %>% dplyr::select(SYMBOL, nearest_fold),
            by = c("gene" = "SYMBOL")) %>%
  dplyr::filter(!is.na(nearest_fold)) %>%
  mutate(
    mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
    mecp2_direction = ifelse(nearest_fold > 0, "MeCP2 Up", "MeCP2 Down")
  )

mc_mecp2_result_19h <- build_2x2_heatmap(
  mc_mecp2_19h$mc_direction, mc_mecp2_19h$mecp2_direction,
  c("mC Up", "mC Down"), c("MeCP2 Up", "MeCP2 Down"),
  list(c("mC Up", "MeCP2 Up"), c("mC Down", "MeCP2 Down")),
  title = "", x_label = "", y_label = ""
)

# hmC x MeCP2
hmc_mecp2_19h <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(mecp2_gene_19h %>% dplyr::select(SYMBOL, nearest_fold),
            by = c("gene" = "SYMBOL")) %>%
  dplyr::filter(!is.na(nearest_fold)) %>%
  mutate(
    hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
    mecp2_direction = ifelse(nearest_fold > 0, "MeCP2 Up", "MeCP2 Down")
  )

hmc_mecp2_result_19h <- build_2x2_heatmap(
  hmc_mecp2_19h$hmc_direction, hmc_mecp2_19h$mecp2_direction,
  c("hmC Down", "hmC Up"), c("MeCP2 Up", "MeCP2 Down"),
  list(c("hmC Down", "MeCP2 Up"), c("hmC Up", "MeCP2 Down")),
  title = "", x_label = "", y_label = ""
)

# --- Re-compute ATAC O/E ---
cat("  Re-computing ATAC O/E...\n")
atac_up_gr_19h   <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down_gr_19h <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")

atac_up_anno_19h   <- annotate_peaks_to_genes(atac_up_gr_19h, txdb)
atac_down_anno_19h <- annotate_peaks_to_genes(atac_down_gr_19h, txdb)

atac_up_gene_19h   <- aggregate_peaks_by_gene(atac_up_anno_19h, "atac_up")
atac_down_gene_19h <- aggregate_peaks_by_gene(atac_down_anno_19h, "atac_down")

atac_gene_19h <- full_join(atac_up_gene_19h, atac_down_gene_19h, by = "SYMBOL") %>%
  mutate(
    n_atac_up = replace_na(n_atac_up, 0),
    n_atac_down = replace_na(n_atac_down, 0),
    net_atac = n_atac_up - n_atac_down
  )

# mC x ATAC
mc_atac_19h <- mc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(atac_gene_19h %>% dplyr::select(SYMBOL, net_atac),
            by = c("gene" = "SYMBOL")) %>%
  mutate(net_atac = replace_na(net_atac, 0)) %>%
  dplyr::filter(net_atac != 0) %>%
  mutate(
    mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
    atac_direction = ifelse(net_atac > 0, "ATAC Up", "ATAC Down")
  )

mc_atac_result_19h <- build_2x2_heatmap(
  mc_atac_19h$mc_direction, mc_atac_19h$atac_direction,
  c("mC Up", "mC Down"), c("ATAC Up", "ATAC Down"),
  list(c("mC Up", "ATAC Down"), c("mC Down", "ATAC Up")),
  title = "", x_label = "", y_label = ""
)

# hmC x ATAC
hmc_atac_19h <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(atac_gene_19h %>% dplyr::select(SYMBOL, net_atac),
            by = c("gene" = "SYMBOL")) %>%
  mutate(net_atac = replace_na(net_atac, 0)) %>%
  dplyr::filter(net_atac != 0) %>%
  mutate(
    hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
    atac_direction = ifelse(net_atac > 0, "ATAC Up", "ATAC Down")
  )

hmc_atac_result_19h <- build_2x2_heatmap(
  hmc_atac_19h$hmc_direction, hmc_atac_19h$atac_direction,
  c("hmC Down", "hmC Up"), c("ATAC Up", "ATAC Down"),
  list(c("hmC Down", "ATAC Down"), c("hmC Up", "ATAC Up")),
  title = "", x_label = "", y_label = ""
)

# --- Re-compute K119ub O/E ---
cat("  Re-computing K119ub O/E...\n")
k119ub_ctrl_gr_19h <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Control")
k119ub_mut_gr_19h  <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mutant")

k119ub_ctrl_anno_19h <- annotate_peaks_to_genes(k119ub_ctrl_gr_19h, txdb)
k119ub_mut_anno_19h  <- annotate_peaks_to_genes(k119ub_mut_gr_19h, txdb)

k119ub_ctrl_gene_19h <- aggregate_peaks_by_gene(k119ub_ctrl_anno_19h, "k119ub_ctrl")
k119ub_mut_gene_19h  <- aggregate_peaks_by_gene(k119ub_mut_anno_19h, "k119ub_mut")

k119ub_gene_19h <- full_join(k119ub_ctrl_gene_19h, k119ub_mut_gene_19h, by = "SYMBOL") %>%
  mutate(
    n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
    n_k119ub_mut = replace_na(n_k119ub_mut, 0),
    net_k119ub = n_k119ub_mut - n_k119ub_ctrl
  )

# mC x K119ub
mc_k119ub_19h <- mc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(k119ub_gene_19h %>% dplyr::select(SYMBOL, net_k119ub),
            by = c("gene" = "SYMBOL")) %>%
  mutate(net_k119ub = replace_na(net_k119ub, 0)) %>%
  dplyr::filter(net_k119ub != 0) %>%
  mutate(
    mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
    k119ub_direction = ifelse(net_k119ub > 0, "K119ub Gained", "K119ub Lost")
  )

mc_k119ub_result_19h <- build_2x2_heatmap(
  mc_k119ub_19h$mc_direction, mc_k119ub_19h$k119ub_direction,
  c("mC Up", "mC Down"), c("K119ub Gained", "K119ub Lost"),
  list(c("mC Up", "K119ub Gained"), c("mC Down", "K119ub Lost")),
  title = "", x_label = "", y_label = ""
)

# hmC x K119ub
hmc_k119ub_19h <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  left_join(k119ub_gene_19h %>% dplyr::select(SYMBOL, net_k119ub),
            by = c("gene" = "SYMBOL")) %>%
  mutate(net_k119ub = replace_na(net_k119ub, 0)) %>%
  dplyr::filter(net_k119ub != 0) %>%
  mutate(
    hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
    k119ub_direction = ifelse(net_k119ub > 0, "K119ub Gained", "K119ub Lost")
  )

hmc_k119ub_result_19h <- build_2x2_heatmap(
  hmc_k119ub_19h$hmc_direction, hmc_k119ub_19h$k119ub_direction,
  c("hmC Down", "hmC Up"), c("K119ub Gained", "K119ub Lost"),
  list(c("hmC Down", "K119ub Gained"), c("hmC Up", "K119ub Lost")),
  title = "", x_label = "", y_label = ""
)

# --- H3K27ac O/E already computed in 19g (result_19g_mc, result_19g_hmc) ---

# --- Collect O/E values across all 4 marks x 2 perspectives ---
cat("  Collecting O/E comparison data...\n")

comparison_rows_19h <- list(
  extract_enriched_oe(mc_mecp2_result_19h, "5mC", "MeCP2"),
  extract_enriched_oe(hmc_mecp2_result_19h, "5hmC", "MeCP2"),
  extract_enriched_oe(mc_atac_result_19h, "5mC", "ATAC"),
  extract_enriched_oe(hmc_atac_result_19h, "5hmC", "ATAC"),
  extract_enriched_oe(mc_k119ub_result_19h, "5mC", "K119ub"),
  extract_enriched_oe(hmc_k119ub_result_19h, "5hmC", "K119ub"),
  extract_enriched_oe(result_19g_mc, "5mC", "H3K27ac"),
  extract_enriched_oe(result_19g_hmc, "5hmC", "H3K27ac")
)

comparison_df_19h <- do.call(rbind, comparison_rows_19h)

cat("  O/E comparison summary:\n")
for (i in 1:nrow(comparison_df_19h)) {
  cat(sprintf("    %s | %s | %s: O/E = %.2f (Obs=%d, Exp=%.0f)\n",
              comparison_df_19h$perspective[i], comparison_df_19h$mark[i],
              comparison_df_19h$quadrant[i],
              comparison_df_19h$oe_ratio[i], comparison_df_19h$observed[i],
              comparison_df_19h$expected[i]))
}

# Format labels for plot
comparison_df_19h$p_label <- ifelse(comparison_df_19h$fisher_p < 0.001,
                                    sprintf("p=%.1e", comparison_df_19h$fisher_p),
                                    sprintf("p=%.3f", comparison_df_19h$fisher_p))
comparison_df_19h$short_quadrant <- gsub(
  "hmC Down \\+ |hmC Up \\+ |mC Up \\+ |mC Down \\+ ", "", comparison_df_19h$quadrant)
comparison_df_19h$mark <- factor(comparison_df_19h$mark,
                                 levels = c("MeCP2", "ATAC", "K119ub", "H3K27ac"))
comparison_df_19h$perspective <- factor(comparison_df_19h$perspective,
                                        levels = c("5mC", "5hmC"))

p_19h <- ggplot(comparison_df_19h,
                aes(x = mark, y = oe_ratio, color = perspective, shape = perspective)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sprintf("%.2f", oe_ratio)),
            position = position_dodge(width = 0.5),
            vjust = -1.2, size = 3.5, show.legend = FALSE) +
  geom_text(aes(label = short_quadrant),
            position = position_dodge(width = 0.5),
            vjust = 2.2, size = 2.5, show.legend = FALSE, color = "grey40") +
  scale_color_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8"),
                     name = "DMR Perspective") +
  scale_shape_manual(values = c("5mC" = 16, "5hmC" = 17),
                     name = "DMR Perspective") +
  scale_y_continuous(limits = c(min(comparison_df_19h$oe_ratio, 1) * 0.85,
                                max(comparison_df_19h$oe_ratio, 1) * 1.2)) +
  facet_wrap(~ mark, scales = "free_x", nrow = 1) +
  labs(
    title = "Chromatin Mark O/E Enrichment Across 4 Marks",
    subtitle = "Each point = enriched quadrant from 2x2 integration | Dashed line = null (O/E=1)",
    x = "", y = "Observed / Expected Enrichment"
  ) +
  theme_biomodal() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1.5, "lines")
  )

save_multiformat_ggplot(p_19h, file.path(OUTPUT_DIR, "19h_chromatin_mark_oe_comparison"),
                        width = 14, height = 8)

# Save comparison table
write.table(comparison_df_19h,
            file.path(TABLES_DIR, "h3k27ac_all_marks_oe_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: h3k27ac_all_marks_oe_comparison.tsv\n")

# =============================================================================
# STATISTICAL TESTS (console output)
# =============================================================================

cat("\n--- STATISTICAL TESTS ---\n\n")

# Build net change data for stats (genes with peaks)
net_data <- all_classified %>%
  dplyr::filter(h3k27ac_status != "No Peaks") %>%
  dplyr::select(gene, met_group, net_h3k27ac)

bg_net <- all_dmr_classified %>%
  dplyr::filter(h3k27ac_status != "No Peaks") %>%
  pull(net_h3k27ac)

cat("1. Fisher's exact: each group's gain rate vs background (all DMR genes)\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  p <- fisher_results[[grp]]$p.value
  or <- fisher_results[[grp]]$estimate
  cat(sprintf("   %-10s: OR=%.2f, p=%s\n", grp, or,
              ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.4f", p))))
}

cat("\n2. Wilcoxon rank-sum: each group's net change vs All DMR Genes\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  grp_net <- net_data$net_h3k27ac[net_data$met_group == grp]
  w <- wilcox.test(grp_net, bg_net)
  cat(sprintf("   %-10s: W=%.0f, p=%s\n", grp, w$statistic,
              ifelse(w$p.value < 0.001, sprintf("%.2e", w$p.value),
                     sprintf("%.4f", w$p.value))))
}

cat("\n3. Wilcoxon one-sample: each group's net change != 0\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up", "All DMR Genes")) {
  if (grp == "All DMR Genes") {
    grp_net <- bg_net
  } else {
    grp_net <- net_data$net_h3k27ac[net_data$met_group == grp]
  }
  w <- wilcox.test(grp_net, mu = 0)
  cat(sprintf("   %-14s: V=%.0f, p=%s, median=%+.0f\n", grp, w$statistic,
              ifelse(w$p.value < 0.001, sprintf("%.2e", w$p.value),
                     sprintf("%.4f", w$p.value)),
              median(grp_net)))
}

cat("\n4. Spearman correlation: mod_difference vs net_h3k27ac\n")
for (mod in c("5mC", "5hmC")) {
  sub <- scatter_all %>% dplyr::filter(modality == mod)
  ct <- cor.test(sub$mod_difference, sub$net_h3k27ac, method = "spearman")
  cat(sprintf("   %s: rho=%+.4f, S=%.0f, p=%s, n=%d\n",
              mod, ct$estimate, ct$statistic,
              ifelse(ct$p.value < 0.001, sprintf("%.2e", ct$p.value),
                     sprintf("%.4f", ct$p.value)),
              nrow(sub)))
}

cat("\n5. Cliff's delta (effect size): each group vs All DMR Genes\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  grp_net <- net_data$net_h3k27ac[net_data$met_group == grp]
  d <- cliffs_delta(grp_net, bg_net)
  magnitude <- if (abs(d) < 0.147) "negligible"
               else if (abs(d) < 0.33) "small"
               else if (abs(d) < 0.474) "medium"
               else "large"
  cat(sprintf("   %-10s: delta=%+.4f (%s)\n", grp, d, magnitude))
}

cat("\n6. Fisher's exact: condition-specific H3K27ac overlap at DMRs (19f)\n")
cat(sprintf("   OR = %.2f, p = %s\n", fisher_19f$estimate,
            ifelse(fisher_19f$p.value < 0.001, sprintf("%.2e", fisher_19f$p.value),
                   sprintf("%.4f", fisher_19f$p.value))))
cat(sprintf("   Hyper DMRs: %.1f%% ctrl, %.1f%% mut H3K27ac overlap\n",
            100 * hyper_ctrl_ov / n_hyper, 100 * hyper_mut_ov / n_hyper))
cat(sprintf("   Hypo DMRs:  %.1f%% ctrl, %.1f%% mut H3K27ac overlap\n",
            100 * hypo_ctrl_ov / n_hypo, 100 * hypo_mut_ov / n_hypo))

cat("\n7. Fisher's exact: 2x2 O/E heatmaps (19g)\n")
cat(sprintf("   mC x H3K27ac:  Fisher OR = %.2f, p = %s (n=%d genes)\n",
            result_19g_mc$fisher$estimate,
            ifelse(result_19g_mc$fisher$p.value < 0.001,
                   sprintf("%.2e", result_19g_mc$fisher$p.value),
                   sprintf("%.4f", result_19g_mc$fisher$p.value)),
            nrow(mc_h3k27ac_oe)))
cat(sprintf("   hmC x H3K27ac: Fisher OR = %.2f, p = %s (n=%d genes)\n",
            result_19g_hmc$fisher$estimate,
            ifelse(result_19g_hmc$fisher$p.value < 0.001,
                   sprintf("%.2e", result_19g_hmc$fisher$p.value),
                   sprintf("%.4f", result_19g_hmc$fisher$p.value)),
            nrow(hmc_h3k27ac_oe)))

cat("\n8. O/E comparison across 4 chromatin marks (19h)\n")
for (m in c("MeCP2", "ATAC", "K119ub", "H3K27ac")) {
  mc_rows <- comparison_df_19h[comparison_df_19h$perspective == "5mC" &
                                comparison_df_19h$mark == m, ]
  hmc_rows <- comparison_df_19h[comparison_df_19h$perspective == "5hmC" &
                                 comparison_df_19h$mark == m, ]
  if (nrow(mc_rows) > 0 && nrow(hmc_rows) > 0) {
    mc_mean <- mean(mc_rows$oe_ratio)
    hmc_mean <- mean(hmc_rows$oe_ratio)
    winner <- ifelse(hmc_mean > mc_mean, "hmC", "mC")
    cat(sprintf("   %s: mC mean O/E = %.2f, hmC mean O/E = %.2f -> %s stronger\n",
                m, mc_mean, hmc_mean, winner))
  }
}

# =============================================================================
# SAVE SUMMARY TABLE
# =============================================================================

cat("\n--- Saving summary table ---\n\n")

table_data <- all_classified %>%
  group_by(met_group) %>%
  summarise(
    total_genes = n(),
    genes_with_peaks = sum(h3k27ac_status != "No Peaks"),
    genes_no_peaks = sum(h3k27ac_status == "No Peaks"),
    pct_no_peaks = 100 * genes_no_peaks / total_genes,
    genes_gained = sum(h3k27ac_status == "Gained"),
    genes_equal = sum(h3k27ac_status == "Equal Peaks"),
    genes_lost = sum(h3k27ac_status == "Lost"),
    pct_gained_of_total = 100 * genes_gained / total_genes,
    pct_gained_of_peaks = ifelse(genes_with_peaks > 0,
                                  100 * genes_gained / genes_with_peaks, NA),
    median_net_change = median(net_h3k27ac[h3k27ac_status != "No Peaks"], na.rm = TRUE),
    .groups = "drop"
  )

# Add background row
bg_table <- all_dmr_classified %>%
  summarise(
    met_group = "All DMR Genes",
    total_genes = n(),
    genes_with_peaks = sum(h3k27ac_status != "No Peaks"),
    genes_no_peaks = sum(h3k27ac_status == "No Peaks"),
    pct_no_peaks = 100 * genes_no_peaks / total_genes,
    genes_gained = sum(h3k27ac_status == "Gained"),
    genes_equal = sum(h3k27ac_status == "Equal Peaks"),
    genes_lost = sum(h3k27ac_status == "Lost"),
    pct_gained_of_total = 100 * genes_gained / total_genes,
    pct_gained_of_peaks = ifelse(genes_with_peaks > 0,
                                  100 * genes_gained / genes_with_peaks, NA),
    median_net_change = median(net_h3k27ac[h3k27ac_status != "No Peaks"], na.rm = TRUE)
  )

table_data <- bind_rows(table_data, bg_table)

# Add Cliff's delta and Wilcoxon p vs background
table_data$cliffs_delta_vs_bg <- NA_real_
table_data$wilcox_p_vs_bg <- NA_real_
for (i in 1:nrow(table_data)) {
  grp <- as.character(table_data$met_group[i])
  if (grp == "All DMR Genes") next
  grp_net <- net_data$net_h3k27ac[net_data$met_group == grp]
  if (length(grp_net) > 0) {
    table_data$cliffs_delta_vs_bg[i] <- cliffs_delta(grp_net, bg_net)
    table_data$wilcox_p_vs_bg[i] <- wilcox.test(grp_net, bg_net)$p.value
  }
}

# Add Fisher p-values
table_data$fisher_gain_p <- NA_real_
for (i in 1:nrow(table_data)) {
  grp <- as.character(table_data$met_group[i])
  if (grp %in% names(fisher_results)) {
    table_data$fisher_gain_p[i] <- fisher_results[[grp]]$p.value
  }
}

write.table(table_data,
            file.path(TABLES_DIR, "h3k27ac_peak_analysis_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: h3k27ac_peak_analysis_summary.tsv\n")

# =============================================================================
# SECTION 19 SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 19 SUMMARY: H3K27AC PEAK ANALYSIS AT DMR GENES\n")
cat("================================================================================\n\n")

cat("Per methylation group:\n")
cat(sprintf("  %-14s %6s %8s %8s %8s %8s %8s %10s\n",
            "Group", "Total", "No Peak", "Gained", "Equal", "Lost",
            "%NoPeak", "%Gain|Peak"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (i in 1:nrow(table_data)) {
  r <- table_data[i, ]
  cat(sprintf("  %-14s %6d %8d %8d %8d %8d %7.1f%% %9.1f%%\n",
              r$met_group, r$total_genes, r$genes_no_peaks,
              r$genes_gained, r$genes_equal, r$genes_lost,
              r$pct_no_peaks,
              ifelse(is.na(r$pct_gained_of_peaks), 0, r$pct_gained_of_peaks)))
}

cat(sprintf("\nBackground gain rate (all DMR genes with peaks): %.1f%%\n", bg_gain_rate))
cat(sprintf("Background loss rate (all DMR genes with peaks): %.1f%%\n", bg_loss_rate))

cat("\nFigures generated:\n")
cat("  19a: H3K27ac status breakdown (100% stacked bars)\n")
cat("  19b: H3K27ac conditional direction (grouped bars + Fisher)\n")
cat("  19c: H3K27ac per-gene effect sizes (strip + box)\n")
cat("  19d: Methylation change vs H3K27ac scatter (Spearman)\n")
cat("  19e: Gene waterfall ranked by H3K27ac change\n")
cat("  19f: Condition-specific H3K27ac overlap at DMRs\n")
cat("  19g: 2x2 O/E heatmaps (mC x H3K27ac, hmC x H3K27ac)\n")
cat("  19h: Comparative O/E dot plot across 4 chromatin marks\n")

cat("\nCorrelations (mod_difference vs net_h3k27ac):\n")
for (i in 1:nrow(cor_results)) {
  cat(sprintf("  %s: rho = %+.4f\n", cor_results$modality[i], cor_results$rho[i]))
}

cat("\nCondition-specific overlap (19f):\n")
cat(sprintf("  Fisher OR = %.2f, p = %.2e\n", fisher_19f$estimate, fisher_19f$p.value))

cat("\nO/E heatmaps (19g):\n")
cat(sprintf("  mC x H3K27ac:  Fisher OR = %.2f, p = %.2e\n",
            result_19g_mc$fisher$estimate, result_19g_mc$fisher$p.value))
cat(sprintf("  hmC x H3K27ac: Fisher OR = %.2f, p = %.2e\n",
            result_19g_hmc$fisher$estimate, result_19g_hmc$fisher$p.value))

cat("\nComparative O/E (19h): 4 marks x 2 perspectives\n")
for (m in c("MeCP2", "ATAC", "K119ub", "H3K27ac")) {
  mc_rows <- comparison_df_19h[comparison_df_19h$perspective == "5mC" &
                                comparison_df_19h$mark == m, ]
  hmc_rows <- comparison_df_19h[comparison_df_19h$perspective == "5hmC" &
                                 comparison_df_19h$mark == m, ]
  if (nrow(mc_rows) > 0 && nrow(hmc_rows) > 0) {
    cat(sprintf("  %s: mC O/E = %.2f, hmC O/E = %.2f\n",
                m, mean(mc_rows$oe_ratio), mean(hmc_rows$oe_ratio)))
  }
}

cat("\nSection 19 complete.\n\n")
