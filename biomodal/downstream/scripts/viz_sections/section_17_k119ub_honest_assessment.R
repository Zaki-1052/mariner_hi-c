# biomodal/downstream/scripts/viz_sections/section_17_k119ub_honest_assessment.R
# Section 17: Honest K119ub Assessment at DMR Genes
#
# Sections 15 (O/E heatmaps) and 16c (concordance bars) both present the
# K119ub-methylation relationship in ways that can mislead:
#   - Section 16c: 75.5% of mC Up genes show K119ub Gained -- but this excludes
#     the majority of mC Up genes that have NO K119ub signal at all, and doesn't
#     show the background gain rate (~48%)
#   - Section 15: O/E = 1.08 for hmC Down -> K119ub Gained (barely above chance)
#
# This section presents the data honestly:
#   17a: Complete K119ub status breakdown (100% stacked bars including "No Peaks")
#   17b: Conditional direction among genes WITH peaks (with background reference)
#   17c: Per-gene K119ub effect size distribution (strip + box plot)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_17_k119ub_honest_assessment.R

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

# =============================================================================
# SECTION 17: HONEST K119UB ASSESSMENT AT DMR GENES
# =============================================================================

cat("================================================================================\n")
cat("SECTION 17: HONEST K119UB ASSESSMENT AT DMR GENES\n")
cat("================================================================================\n\n")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# -----------------------------------------------------------------------
# Prepare methylation gene groups (character vectors)
# -----------------------------------------------------------------------
cat("Preparing methylation gene groups...\n")

mc_sig <- mc_dmr %>% dplyr::filter(significant)
hmc_sig <- hmc_dmr %>% dplyr::filter(significant)

mc_up_genes <- mc_sig$gene[mc_sig$mod_difference > 0]
mc_down_genes <- mc_sig$gene[mc_sig$mod_difference < 0]
hmc_down_genes <- hmc_sig$gene[hmc_sig$mod_difference < 0]
hmc_up_genes <- hmc_sig$gene[hmc_sig$mod_difference > 0]

cat(sprintf("  mC Up:    %d genes\n", length(mc_up_genes)))
cat(sprintf("  mC Down:  %d genes\n", length(mc_down_genes)))
cat(sprintf("  hmC Down: %d genes\n", length(hmc_down_genes)))
cat(sprintf("  hmC Up:   %d genes\n", length(hmc_up_genes)))

# -----------------------------------------------------------------------
# Load K119ub peaks and annotate to genes
# -----------------------------------------------------------------------
cat("\nLoading K119ub condition-specific peaks...\n")

k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Control")
k119ub_mut_gr <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mutant")

if (is.null(k119ub_ctrl_gr) || is.null(k119ub_mut_gr)) {
  stop("K119ub peak files not found. Cannot proceed with Section 17.")
}

cat("  Annotating K119ub peaks to genes via ChIPseeker...\n")
k119ub_ctrl_gene <- aggregate_peaks_by_gene(
  annotate_peaks_to_genes(k119ub_ctrl_gr, txdb), "k119ub_ctrl")
k119ub_mut_gene <- aggregate_peaks_by_gene(
  annotate_peaks_to_genes(k119ub_mut_gr, txdb), "k119ub_mut")

# Build per-gene K119ub summary (all genes that have ANY peak in either condition)
k119ub_gene <- full_join(k119ub_ctrl_gene, k119ub_mut_gene, by = "SYMBOL") %>%
  mutate(
    n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
    n_k119ub_mut = replace_na(n_k119ub_mut, 0),
    net_k119ub = n_k119ub_mut - n_k119ub_ctrl
  )

cat(sprintf("  %d unique genes with any K119ub peak\n", nrow(k119ub_gene)))

# -----------------------------------------------------------------------
# Classify every DMR gene into 4 categories
# -----------------------------------------------------------------------
cat("\nClassifying K119ub status for all DMR genes...\n")

# All unique significant DMR genes across both modalities
all_dmr_genes <- unique(c(mc_up_genes, mc_down_genes, hmc_down_genes, hmc_up_genes))
cat(sprintf("  %d unique significant DMR genes total\n", length(all_dmr_genes)))

classify_k119ub_status <- function(gene_list, k119ub_df) {
  # For each gene: look up ctrl/mut peak counts
  gene_df <- data.frame(gene = gene_list, stringsAsFactors = FALSE) %>%
    left_join(k119ub_df, by = c("gene" = "SYMBOL")) %>%
    mutate(
      n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
      n_k119ub_mut = replace_na(n_k119ub_mut, 0),
      net_k119ub = n_k119ub_mut - n_k119ub_ctrl,
      k119ub_status = case_when(
        n_k119ub_ctrl == 0 & n_k119ub_mut == 0 ~ "No Peaks",
        net_k119ub > 0 ~ "Gained",
        net_k119ub < 0 ~ "Lost",
        TRUE ~ "Equal Peaks"
      )
    )
  gene_df
}

# Classify each methylation group
mc_up_classified <- classify_k119ub_status(mc_up_genes, k119ub_gene) %>%
  mutate(met_group = "mC Up")
mc_down_classified <- classify_k119ub_status(mc_down_genes, k119ub_gene) %>%
  mutate(met_group = "mC Down")
hmc_down_classified <- classify_k119ub_status(hmc_down_genes, k119ub_gene) %>%
  mutate(met_group = "hmC Down")
hmc_up_classified <- classify_k119ub_status(hmc_up_genes, k119ub_gene) %>%
  mutate(met_group = "hmC Up")

all_classified <- bind_rows(mc_up_classified, mc_down_classified,
                            hmc_down_classified, hmc_up_classified)

# Also classify ALL DMR genes for background rate
all_dmr_classified <- classify_k119ub_status(all_dmr_genes, k119ub_gene) %>%
  mutate(met_group = "All DMR Genes")

# Compute summary per group
compute_status_summary <- function(classified_df) {
  classified_df %>%
    group_by(met_group, k119ub_status) %>%
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

# Print overview
cat("\n  K119ub status breakdown:\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  sub <- group_summary %>% dplyr::filter(met_group == grp)
  total <- sub$total[1]
  no_peaks <- sub$n[sub$k119ub_status == "No Peaks"]
  if (length(no_peaks) == 0) no_peaks <- 0
  cat(sprintf("    %-10s: %d total, %d (%.1f%%) no peaks\n",
              grp, total, no_peaks, 100 * no_peaks / total))
}

# =============================================================================
# FIGURE 17a: Complete K119ub Status Breakdown (100% Stacked Bars)
# =============================================================================
cat("\n--- FIGURE 17a: Complete K119ub Status Breakdown ---\n\n")

# Order groups and status categories
status_order <- c("Gained", "Equal Peaks", "No Peaks", "Lost")
status_colors <- c(
  "Gained" = "#756BB1",     # COLORS$k119ub purple
  "Equal Peaks" = "#BDBDBD", # light grey
  "No Peaks" = "#636363",    # dark grey
  "Lost" = "#74C476"         # COLORS$k119ub green
)

group_summary$met_group <- factor(group_summary$met_group,
                                   levels = c("mC Up", "mC Down", "hmC Down", "hmC Up"))
group_summary$k119ub_status <- factor(group_summary$k119ub_status, levels = status_order)

# Label each segment with count and percentage
group_summary$label <- sprintf("n=%d\n(%.0f%%)", group_summary$n, group_summary$pct)

p_17a <- ggplot(group_summary,
                aes(x = met_group, y = pct, fill = k119ub_status)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 2.8, lineheight = 0.9, color = "white", fontface = "bold") +
  scale_fill_manual(values = status_colors, name = "K119ub Status",
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "H2AK119ub Status Across ALL DMR Genes",
    subtitle = "Including genes with no K119ub peaks (majority in every group)",
    x = "Methylation Group (Significant DMR Genes)",
    y = "Percentage of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_17a,
                        file.path(OUTPUT_DIR, "17a_k119ub_full_breakdown"),
                        width = 9, height = 7)

# Print figure 17a details
cat("  Per-group breakdown:\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  sub <- group_summary %>%
    dplyr::filter(met_group == grp) %>%
    arrange(match(k119ub_status, status_order))
  cat(sprintf("    %s (n=%d):\n", grp, sub$total[1]))
  for (i in 1:nrow(sub)) {
    cat(sprintf("      %-12s: %4d (%5.1f%%)\n",
                as.character(sub$k119ub_status[i]), sub$n[i], sub$pct[i]))
  }
}

# =============================================================================
# FIGURE 17b: Conditional Direction Among Genes WITH K119ub Peaks
# =============================================================================
cat("\n--- FIGURE 17b: Conditional Direction (Genes WITH K119ub Peaks) ---\n\n")

# Compute conditional rates: among genes with peaks, what % gained/equal/lost?
conditional_summary <- all_classified %>%
  dplyr::filter(k119ub_status != "No Peaks") %>%
  group_by(met_group) %>%
  summarise(
    n_with_peaks = n(),
    n_gained = sum(k119ub_status == "Gained"),
    n_equal = sum(k119ub_status == "Equal Peaks"),
    n_lost = sum(k119ub_status == "Lost"),
    pct_gained = 100 * n_gained / n(),
    pct_equal = 100 * n_equal / n(),
    pct_lost = 100 * n_lost / n(),
    .groups = "drop"
  )

# Background: all DMR genes with peaks
bg <- all_dmr_classified %>%
  dplyr::filter(k119ub_status != "No Peaks") %>%
  summarise(
    met_group = "All DMR Genes",
    n_with_peaks = n(),
    n_gained = sum(k119ub_status == "Gained"),
    n_equal = sum(k119ub_status == "Equal Peaks"),
    n_lost = sum(k119ub_status == "Lost"),
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

# Bar colors matching k119ub scheme
direction_colors <- c(
  "Gained" = "#756BB1",
  "Equal Peaks" = "#BDBDBD",
  "Lost" = "#74C476"
)

# Background gain and loss rates for reference lines
bg_gain_rate <- bg$pct_gained
bg_loss_rate <- bg$pct_lost

# Fisher exact tests: each group's gain rate vs background gain rate
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

p_17b <- ggplot(conditional_long,
                aes(x = met_group, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d)", pct, n_count)),
            position = position_dodge(width = 0.75),
            vjust = -0.2, size = 2.6, lineheight = 0.85) +
  # Background gain rate reference line
  geom_hline(yintercept = bg_gain_rate, linetype = "dashed",
             color = "#756BB1", linewidth = 0.5) +
  annotate("text", x = 0.55, y = bg_gain_rate + 1.5,
           label = sprintf("BG gain: %.1f%%", bg_gain_rate),
           size = 2.5, color = "#756BB1", hjust = 0, fontface = "italic") +
  # Background loss rate reference line
  geom_hline(yintercept = bg_loss_rate, linetype = "dashed",
             color = "#74C476", linewidth = 0.5) +
  annotate("text", x = 0.55, y = bg_loss_rate + 1.5,
           label = sprintf("BG loss: %.1f%%", bg_loss_rate),
           size = 2.5, color = "#74C476", hjust = 0, fontface = "italic") +
  # Fisher p-value annotations above each methylation group (gain vs background)
  annotate("text",
           x = 1:4,
           y = rep(max(conditional_long$pct, na.rm = TRUE) * 1.12, 4),
           label = sapply(c("mC Up", "mC Down", "hmC Down", "hmC Up"), function(g) {
             p <- fisher_results[[g]]$p.value
             if (p < 0.001) sprintf("p=%.1e", p)
             else sprintf("p=%.3f", p)
           }),
           size = 2.5, fontface = "italic", color = "grey30") +
  scale_fill_manual(values = direction_colors, name = "K119ub Direction") +
  scale_y_continuous(limits = c(0, max(conditional_long$pct, na.rm = TRUE) * 1.25),
                     expand = c(0, 0)) +
  labs(
    title = "K119ub Direction Among Genes WITH Peaks Only",
    subtitle = "Dashed lines = background gain/loss rate across all DMR genes with peaks",
    x = "Methylation Group",
    y = "% of Genes with K119ub Peaks"
  ) +
  theme_biomodal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 0.7))

save_multiformat_ggplot(p_17b,
                        file.path(OUTPUT_DIR, "17b_k119ub_conditional_direction"),
                        width = 10, height = 7)

# Print figure 17b details
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
for (grp in c("mC Up", "hmC Down")) {
  row <- conditional_summary %>% dplyr::filter(met_group == grp)
  diff_pp <- row$pct_gained - bg_gain_rate
  cat(sprintf("  %s gain rate: %.1f%% (%+.1f pp above background)\n",
              grp, row$pct_gained, diff_pp))
}

# =============================================================================
# FIGURE 17c: Per-Gene K119ub Effect Size Distribution
# =============================================================================
cat("\n--- FIGURE 17c: Per-Gene K119ub Effect Size ---\n\n")

# Only genes with at least one peak in either condition
effect_df <- all_classified %>%
  dplyr::filter(k119ub_status != "No Peaks") %>%
  dplyr::select(gene, met_group, net_k119ub, n_k119ub_ctrl, n_k119ub_mut)

effect_df$met_group <- factor(effect_df$met_group,
                               levels = c("mC Up", "mC Down", "hmC Down", "hmC Up"))

# Compute medians for annotation
medians <- effect_df %>%
  group_by(met_group) %>%
  summarise(
    median_net = median(net_k119ub),
    iqr_low = quantile(net_k119ub, 0.25),
    iqr_high = quantile(net_k119ub, 0.75),
    n = n(),
    .groups = "drop"
  )

# Identify extreme genes for labeling
extreme_genes <- effect_df %>%
  group_by(met_group) %>%
  slice_max(abs(net_k119ub), n = 3, with_ties = FALSE) %>%
  dplyr::filter(abs(net_k119ub) >= 5) %>%
  ungroup()

# Color by methylation type
effect_df$met_type <- ifelse(grepl("^mC", effect_df$met_group), "5mC", "5hmC")

p_17c <- ggplot(effect_df, aes(x = met_group, y = net_k119ub, color = met_type)) +
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
    title = "Per-Gene K119ub Change (Mutant - Control Peak Count)",
    subtitle = "Only genes with at least one K119ub peak; most changes are +/-1 peak",
    x = "Methylation Group",
    y = "Net K119ub Peaks (Mutant - Control)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_17c,
                        file.path(OUTPUT_DIR, "17c_k119ub_effect_sizes"),
                        width = 9, height = 7)

# Print figure 17c details
cat("  Effect size distributions (genes with peaks):\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  sub <- effect_df %>% dplyr::filter(met_group == grp)
  cat(sprintf("    %-10s: n=%d, median=%+.0f, IQR=[%+.0f, %+.0f], range=[%+.0f, %+.0f]\n",
              grp, nrow(sub),
              median(sub$net_k119ub),
              quantile(sub$net_k119ub, 0.25),
              quantile(sub$net_k119ub, 0.75),
              min(sub$net_k119ub),
              max(sub$net_k119ub)))
}

if (nrow(extreme_genes) > 0) {
  cat("\n  Extreme genes (|net| >= 5):\n")
  for (i in 1:nrow(extreme_genes)) {
    cat(sprintf("    %s (%s): net = %+d (ctrl=%d, mut=%d)\n",
                extreme_genes$gene[i], extreme_genes$met_group[i],
                extreme_genes$net_k119ub[i],
                extreme_genes$n_k119ub_ctrl[i],
                extreme_genes$n_k119ub_mut[i]))
  }
}

# =============================================================================
# SAVE TABLE
# =============================================================================
cat("\n--- Saving summary table ---\n\n")

# Build comprehensive summary table
table_data <- all_classified %>%
  group_by(met_group) %>%
  summarise(
    total_genes = n(),
    genes_with_peaks = sum(k119ub_status != "No Peaks"),
    genes_no_peaks = sum(k119ub_status == "No Peaks"),
    pct_no_peaks = 100 * genes_no_peaks / total_genes,
    genes_gained = sum(k119ub_status == "Gained"),
    genes_equal = sum(k119ub_status == "Equal Peaks"),
    genes_lost = sum(k119ub_status == "Lost"),
    pct_gained_of_total = 100 * genes_gained / total_genes,
    pct_gained_of_peaks = ifelse(genes_with_peaks > 0,
                                  100 * genes_gained / genes_with_peaks, NA),
    .groups = "drop"
  )

# Add background row
bg_table <- all_dmr_classified %>%
  summarise(
    met_group = "All DMR Genes",
    total_genes = n(),
    genes_with_peaks = sum(k119ub_status != "No Peaks"),
    genes_no_peaks = sum(k119ub_status == "No Peaks"),
    pct_no_peaks = 100 * genes_no_peaks / total_genes,
    genes_gained = sum(k119ub_status == "Gained"),
    genes_equal = sum(k119ub_status == "Equal Peaks"),
    genes_lost = sum(k119ub_status == "Lost"),
    pct_gained_of_total = 100 * genes_gained / total_genes,
    pct_gained_of_peaks = ifelse(genes_with_peaks > 0,
                                  100 * genes_gained / genes_with_peaks, NA)
  )

table_data <- bind_rows(table_data, bg_table)

write.table(table_data,
            file.path(TABLES_DIR, "k119ub_honest_breakdown.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: k119ub_honest_breakdown.tsv\n")

# =============================================================================
# SECTION 17 SUMMARY
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("SECTION 17 SUMMARY: HONEST K119UB ASSESSMENT\n")
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
              r$pct_no_peaks, r$pct_gained_of_peaks))
}

cat(sprintf("\nBackground gain rate (all DMR genes with peaks): %.1f%%\n", bg_gain_rate))

cat("\nGain enrichment above background:\n")
for (grp in c("mC Up", "mC Down", "hmC Down", "hmC Up")) {
  row <- conditional_summary %>% dplyr::filter(met_group == grp)
  diff_pp <- row$pct_gained - bg_gain_rate
  p <- fisher_results[[grp]]$p.value
  sig <- ifelse(p < 0.05, "*", "ns")
  cat(sprintf("  %-10s: %+.1f pp (%s, Fisher p=%.2e)\n",
              grp, diff_pp, sig, p))
}

cat("\nKey takeaways:\n")
cat("  1. The majority of genes in every group have NO K119ub peaks at all\n")
cat("  2. Gain enrichment at hyper sites is statistically modest above background\n")
cat("  3. Most per-gene changes are +/-1 peak (not large redistributions)\n")

cat("\nSection 17 complete.\n\n")
