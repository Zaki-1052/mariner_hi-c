# loops/scripts/se_abc_integration.R
#
# 3a. DEG-Centric Enhancer Contact Analysis (ABC Integration)
#
# Joins ABC pipeline contact strength metrics to the DEG vs. invariant gene
# framework, testing whether DEGs show more differential enhancer contact than
# matched invariant genes (3a), and sub-classifying by SE status and K27ac
# change at the enhancer (3b).
#
# Usage:
#   Rscript loops/scripts/se_abc_integration.R --timepoint late
#   Rscript loops/scripts/se_abc_integration.R --timepoint early
#   Rscript loops/scripts/se_abc_integration.R --timepoint both
#
# Output:
#   loops/output/superenhancer_analysis/3a_deg_abc_contact/{timepoint}/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

BASE_DIR <- getwd()

source(file.path(BASE_DIR, "loops/scripts/utils/multi_format_output.R"))
source(file.path(BASE_DIR, "loops/scripts/utils/se_utils.R"))

DEG_PADJ <- 0.05
DEG_LFC  <- 0.3
INVARIANT_RATIO <- 3

ABC_GENE_SUMMARY <- file.path(BASE_DIR, "abc/results/gene_level_summary.tsv")
ABC_PAIRS_FILE   <- file.path(BASE_DIR, "abc/results/delta_abc_all_pairs.tsv")
DIFFBIND_K27AC   <- file.path(BASE_DIR, "peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt")

SE_FILES <- list(
  P60    = file.path(BASE_DIR, "peaks/Superenhancers_P60.bed"),
  ENCODE = file.path(BASE_DIR, "peaks/Superenhancers_encode.bed")
)

TIMEPOINT_CONFIG <- list(
  late = list(
    rna_file = file.path(BASE_DIR, "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"),
    label = "Late (Adult)"
  ),
  early = list(
    rna_file = file.path(BASE_DIR, "tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx"),
    label = "Early (P12)"
  )
)

# =============================================================================
# 2. ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  timepoint <- "late"
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      timepoint <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(list(timepoint = timepoint))
}

parsed <- parse_args(args)

# =============================================================================
# 3. ABC LOADING HELPERS
# =============================================================================

load_abc_gene_summary <- function(file) {
  if (!file.exists(file)) stop(sprintf("ABC gene summary not found: %s", file))

  df <- read_tsv(file, show_col_types = FALSE)

  required <- c("TargetGene", "max_delta_abc", "sum_delta_abc",
                 "n_enhancers", "n_gained", "n_lost",
                 "sum_abc_wt", "sum_abc_ko")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing ABC columns: %s", paste(missing, collapse = ", ")))
  }

  df <- df %>%
    dplyr::rename(gene_symbol = TargetGene) %>%
    dplyr::select(gene_symbol, max_delta_abc, sum_delta_abc,
                  max_delta_unnorm = max_delta_unnorm,
                  n_enhancers, n_gained, n_lost,
                  sum_abc_wt, sum_abc_ko, sum_delta_unnorm)

  cat(sprintf("  Loaded ABC gene summary: %d genes\n", nrow(df)))
  return(df)
}

load_abc_pairs <- function(file) {
  if (!file.exists(file)) stop(sprintf("ABC pairs file not found: %s", file))

  df <- read_tsv(file, show_col_types = FALSE)

  required <- c("chr", "start", "end", "TargetGene", "delta_ABC",
                 "hic_contact_pl_scaled_adj_WT", "hic_contact_pl_scaled_adj_KO",
                 "class")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing ABC pair columns: %s", paste(missing, collapse = ", ")))
  }

  n_before <- nrow(df)
  df <- df %>% dplyr::filter(class != "promoter")
  cat(sprintf("  Loaded ABC pairs: %d total, %d distal (removed %d promoter)\n",
              n_before, nrow(df), n_before - nrow(df)))

  return(df)
}

compute_log2fc_contact <- function(contact_wt, contact_ko) {
  nonzero_wt <- contact_wt[contact_wt > 0]
  nonzero_ko <- contact_ko[contact_ko > 0]
  nonzero <- c(nonzero_wt, nonzero_ko)
  pseudo <- if (length(nonzero) > 0) min(nonzero) / 2 else 1e-6
  log2((contact_ko + pseudo) / (contact_wt + pseudo))
}

# =============================================================================
# 4. K27ac WITH FOLD VALUE (extend classify_k27ac_change to retain magnitude)
# =============================================================================

classify_k27ac_with_fold <- function(peaks_gr, diffbind_gr,
                                     lfc_thresh = 0.5, fdr_thresh = 0.05) {
  hits <- findOverlaps(peaks_gr, diffbind_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    return(tibble(
      peak_idx = seq_along(peaks_gr),
      k27ac_class = rep("no_diffbind_peak", length(peaks_gr)),
      k27ac_fold = NA_real_,
      k27ac_fdr = NA_real_
    ))
  }

  hit_df <- tibble(
    peak_idx = queryHits(hits),
    fold = diffbind_gr$fold[subjectHits(hits)],
    fdr = diffbind_gr$fdr[subjectHits(hits)]
  )

  best_hit <- hit_df %>%
    dplyr::group_by(peak_idx) %>%
    dplyr::slice_min(fdr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      k27ac_class = dplyr::case_when(
        fdr < fdr_thresh & fold > lfc_thresh  ~ "gained_k27ac",
        fdr < fdr_thresh & fold < -lfc_thresh ~ "lost_k27ac",
        TRUE                                   ~ "stable_k27ac"
      )
    ) %>%
    dplyr::rename(k27ac_fold = fold, k27ac_fdr = fdr)

  result <- tibble(peak_idx = seq_along(peaks_gr)) %>%
    dplyr::left_join(best_hit %>% dplyr::select(peak_idx, k27ac_class, k27ac_fold, k27ac_fdr),
                     by = "peak_idx") %>%
    dplyr::mutate(k27ac_class = ifelse(is.na(k27ac_class), "no_diffbind_peak", k27ac_class))

  cat(sprintf("  K27ac classification (with fold): %s\n",
              paste(names(table(result$k27ac_class)),
                    table(result$k27ac_class), sep = "=", collapse = ", ")))

  return(result)
}

# =============================================================================
# 5. PLOTTING FUNCTIONS
# =============================================================================

plot_delta_abc_violin <- function(gene_abc, metric, metric_label, output_dir) {
  comparisons <- list(
    c("DEG_up", "invariant"),
    c("DEG_down", "invariant"),
    c("DEG_up", "DEG_down")
  )

  p <- ggplot(gene_abc, aes(x = gene_class, y = .data[[metric]], fill = gene_class)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.format", tip.length = 0.02) +
    scale_fill_manual(values = DEG_CLASS_COLORS) +
    labs(title = metric_label, x = "Gene Class", y = metric_label) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  save_multiformat_ggplot(p, file.path(output_dir, tolower(gsub("[. ]", "_", metric))),
                          width = 6, height = 6)
  return(p)
}

plot_n_gained_lost_bar <- function(gene_abc, output_dir) {
  gained <- gene_abc %>%
    dplyr::group_by(gene_class) %>%
    dplyr::summarise(mean_count = mean(n_gained), se = sd(n_gained) / sqrt(n()), .groups = "drop") %>%
    dplyr::mutate(direction = "gained")
  lost <- gene_abc %>%
    dplyr::group_by(gene_class) %>%
    dplyr::summarise(mean_count = mean(n_lost), se = sd(n_lost) / sqrt(n()), .groups = "drop") %>%
    dplyr::mutate(direction = "lost")
  bar_df <- bind_rows(gained, lost)

  p <- ggplot(bar_df, aes(x = gene_class, y = mean_count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean_count - se, ymax = mean_count + se),
                  position = position_dodge(width = 0.7), width = 0.2) +
    scale_fill_manual(values = c("gained" = "#e0a730", "lost" = "#4575b4"),
                      labels = c("gained" = "Gained E-G Links", "lost" = "Lost E-G Links")) +
    labs(title = "ABC Enhancer-Gene Link Changes", x = "Gene Class",
         y = "Mean Number of Links per Gene", fill = "") +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(output_dir, "n_gained_lost_bar"), width = 7, height = 6)
  return(p)
}

plot_delta_abc_vs_lfc <- function(gene_abc, output_dir) {
  degs_only <- gene_abc %>% dplyr::filter(gene_class != "invariant")
  if (nrow(degs_only) < 5) return(NULL)

  rho <- cor.test(degs_only$sum_delta_abc, degs_only$log2FoldChange, method = "spearman")

  p <- ggplot(degs_only, aes(x = sum_delta_abc, y = log2FoldChange)) +
    geom_point(aes(color = gene_class), alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = DEG_CLASS_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.3f\np = %.2e", rho$estimate, rho$p.value),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(title = "ABC Contact Change vs Gene Expression",
         x = "Sum Delta-ABC", y = "log2FoldChange (KO/WT)", color = "Gene Class") +
    theme_classic(base_size = 14)

  save_multiformat_ggplot(p, file.path(output_dir, "delta_abc_vs_lfc_scatter"),
                          width = 7, height = 6)
  return(p)
}

plot_se_contact_violin <- function(pairs_df, output_dir) {
  plot_df <- pairs_df %>%
    dplyr::filter(gene_class != "invariant", !is.na(log2fc_contact), is.finite(log2fc_contact)) %>%
    dplyr::mutate(enhancer_type = ifelse(is_se, "Superenhancer", "Regular Enhancer"))

  if (nrow(plot_df) < 10) return(NULL)

  p <- ggplot(plot_df, aes(x = enhancer_type, y = log2fc_contact, fill = enhancer_type)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(~gene_class, nrow = 1) +
    scale_fill_manual(values = c("Superenhancer" = "#e0a730", "Regular Enhancer" = "#999999")) +
    coord_cartesian(ylim = c(-3, 3)) +
    labs(title = "Contact Frequency Change: SE vs Regular Enhancers",
         x = "", y = "log2FC Contact (KO/WT)") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", strip.text = element_text(size = 12))

  save_multiformat_ggplot(p, file.path(output_dir, "se_vs_regular_contact_violin"),
                          width = 10, height = 6)
  return(p)
}

plot_se_proportion_bar <- function(pairs_df, output_dir) {
  prop_df <- pairs_df %>%
    dplyr::count(gene_class, is_se) %>%
    dplyr::mutate(enhancer_type = ifelse(is_se, "SE", "Regular"))

  p <- ggplot(prop_df, aes(x = gene_class, y = n, fill = enhancer_type)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = c("SE" = "#e0a730", "Regular" = "#999999")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Proportion of Enhancer Contacts at Superenhancers",
         x = "Gene Class", y = "Proportion of E-G Pairs", fill = "Enhancer Type") +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  save_multiformat_ggplot(p, file.path(output_dir, "se_contact_proportion_bar"),
                          width = 7, height = 6)
  return(p)
}

plot_k27ac_contact_scatter <- function(pairs_df, output_dir) {
  plot_df <- pairs_df %>%
    dplyr::filter(!is.na(k27ac_fold), !is.na(log2fc_contact),
                  is.finite(k27ac_fold), is.finite(log2fc_contact),
                  k27ac_class != "no_diffbind_peak",
                  gene_class != "invariant")

  if (nrow(plot_df) < 10) return(NULL)

  rho <- cor.test(plot_df$k27ac_fold, plot_df$log2fc_contact, method = "spearman")

  p <- ggplot(plot_df, aes(x = k27ac_fold, y = log2fc_contact)) +
    geom_point(aes(color = gene_class), alpha = 0.4, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = DEG_CLASS_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.3f\np = %.2e", rho$estimate, rho$p.value),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(title = "K27ac Change vs Contact Frequency Change",
         x = "K27ac DiffBind Fold Change", y = "log2FC Contact (KO/WT)",
         color = "Gene Class") +
    theme_classic(base_size = 14)

  save_multiformat_ggplot(p, file.path(output_dir, "k27ac_contact_scatter"),
                          width = 7, height = 6)
  return(p)
}

plot_k27ac_class_violin <- function(pairs_df, output_dir) {
  plot_df <- pairs_df %>%
    dplyr::filter(k27ac_class %in% c("gained_k27ac", "stable_k27ac", "lost_k27ac"),
                  !is.na(log2fc_contact), is.finite(log2fc_contact))

  if (nrow(plot_df) < 10) return(NULL)

  comparisons <- list(
    c("gained_k27ac", "lost_k27ac"),
    c("gained_k27ac", "stable_k27ac"),
    c("lost_k27ac", "stable_k27ac")
  )

  p <- ggplot(plot_df, aes(x = k27ac_class, y = log2fc_contact, fill = k27ac_class)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.format", tip.length = 0.02) +
    scale_fill_manual(values = K27AC_CLASS_COLORS) +
    scale_x_discrete(labels = c("gained_k27ac" = "Gained K27ac",
                                "lost_k27ac" = "Lost K27ac",
                                "stable_k27ac" = "Stable K27ac")) +
    coord_cartesian(ylim = c(-3, 3)) +
    labs(title = "Contact Frequency Change by K27ac Status",
         x = "K27ac Change at Enhancer", y = "log2FC Contact (KO/WT)") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  save_multiformat_ggplot(p, file.path(output_dir, "k27ac_class_contact_violin"),
                          width = 7, height = 6)
  return(p)
}

plot_se_k27ac_subclass_bar <- function(pairs_df, output_dir) {
  plot_df <- pairs_df %>%
    dplyr::filter(is_se, gene_class != "invariant",
                  k27ac_class %in% c("gained_k27ac", "lost_k27ac", "stable_k27ac")) %>%
    dplyr::count(gene_class, k27ac_class)

  if (nrow(plot_df) == 0) return(NULL)

  p <- ggplot(plot_df, aes(x = gene_class, y = n, fill = k27ac_class)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = K27AC_CLASS_COLORS,
                      labels = c("gained_k27ac" = "Gained K27ac",
                                 "lost_k27ac" = "Lost K27ac",
                                 "stable_k27ac" = "Stable K27ac")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = c("DEG_up" = "Upregulated DEGs",
                                "DEG_down" = "Downregulated DEGs")) +
    labs(title = "K27ac Status at SE-DEG Contacts",
         x = "", y = "Proportion of SE-E-G Pairs", fill = "K27ac Change") +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  save_multiformat_ggplot(p, file.path(output_dir, "se_k27ac_subclass_bar"),
                          width = 7, height = 6)
  return(p)
}

# =============================================================================
# 6. STATISTICS HELPERS
# =============================================================================

wilcox_line <- function(label, vals_a, vals_b, name_a, name_b) {
  if (length(vals_a) < 3 || length(vals_b) < 3) {
    return(sprintf("  %s: insufficient data (%s n=%d, %s n=%d)",
                   label, name_a, length(vals_a), name_b, length(vals_b)))
  }
  wt <- wilcox.test(vals_a, vals_b)
  sprintf("  %s vs %s: median %.4f vs %.4f, W = %.0f, p = %.2e",
          name_a, name_b, median(vals_a), median(vals_b), wt$statistic, wt$p.value)
}

# =============================================================================
# 7. MAIN ANALYSIS
# =============================================================================

run_analysis <- function(timepoint) {
  config <- TIMEPOINT_CONFIG[[timepoint]]
  cat(sprintf("\n========== 3a. DEG-ABC Contact Analysis: %s ==========\n\n", config$label))

  output_dir <- file.path(BASE_DIR, "loops/output/superenhancer_analysis/3a_deg_abc_contact", timepoint)
  tables_dir <- file.path(output_dir, "tables")
  stats_dir  <- file.path(output_dir, "statistics")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

  stats_lines <- c(
    "DEG-ABC Contact Strength Analysis",
    sprintf("Timepoint: %s", config$label),
    sprintf("Date: %s", Sys.time()),
    sprintf("Thresholds: padj < %g, |LFC| > %g, invariant ratio %d:1\n",
            DEG_PADJ, DEG_LFC, INVARIANT_RATIO)
  )

  # --- Step 1: Load RNA-seq, build gene sets ---
  cat("Step 1: Loading RNA-seq and building gene sets...\n")
  rna <- load_rna_for_se(config$rna_file, padj_thresh = DEG_PADJ, lfc_thresh = DEG_LFC)
  degs <- rna$degs
  all_genes <- rna$all_genes

  all_symbols <- unique(c(degs$gene_symbol, all_genes$gene_symbol))
  gene_coords <- get_gene_coordinates(all_symbols)
  invariants <- select_invariant_genes(degs, all_genes, gene_coords, n_per_deg = INVARIANT_RATIO)

  coord_cols <- gene_coords %>% dplyr::select(gene_symbol, chr, tss, strand)
  deg_tbl <- degs %>%
    dplyr::select(gene_symbol, log2FoldChange, padj, baseMean, deg_direction) %>%
    dplyr::inner_join(coord_cols, by = "gene_symbol") %>%
    dplyr::mutate(gene_class = deg_direction)
  inv_tbl <- invariants %>%
    dplyr::select(gene_symbol, log2FoldChange, padj, baseMean, deg_direction) %>%
    dplyr::inner_join(coord_cols, by = "gene_symbol") %>%
    dplyr::mutate(gene_class = "invariant")
  gene_tbl <- bind_rows(deg_tbl, inv_tbl) %>%
    dplyr::filter(!is.na(chr), !is.na(tss)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  n_up   <- sum(gene_tbl$gene_class == "DEG_up")
  n_down <- sum(gene_tbl$gene_class == "DEG_down")
  n_inv  <- sum(gene_tbl$gene_class == "invariant")
  cat(sprintf("  Gene table: %d DEG_up, %d DEG_down, %d invariant\n", n_up, n_down, n_inv))
  stats_lines <- c(stats_lines, sprintf("Genes: %d DEG_up, %d DEG_down, %d invariant\n",
                                         n_up, n_down, n_inv))

  # --- Step 2: Gene-level ABC join ---
  cat("\nStep 2: Joining gene sets to ABC gene summary...\n")
  abc_genes <- load_abc_gene_summary(ABC_GENE_SUMMARY)

  gene_abc <- gene_tbl %>%
    dplyr::select(gene_symbol, gene_class, log2FoldChange, padj, baseMean) %>%
    dplyr::inner_join(abc_genes, by = "gene_symbol")

  cat(sprintf("  Matched %d / %d genes to ABC summary (%.1f%%)\n",
              nrow(gene_abc), nrow(gene_tbl),
              100 * nrow(gene_abc) / nrow(gene_tbl)))
  stats_lines <- c(stats_lines,
    sprintf("Gene-ABC join: %d / %d matched (%.1f%%)",
            nrow(gene_abc), nrow(gene_tbl), 100 * nrow(gene_abc) / nrow(gene_tbl)),
    sprintf("  DEG_up matched: %d", sum(gene_abc$gene_class == "DEG_up")),
    sprintf("  DEG_down matched: %d", sum(gene_abc$gene_class == "DEG_down")),
    sprintf("  invariant matched: %d\n", sum(gene_abc$gene_class == "invariant"))
  )

  # --- Step 3: Gene-level plots and stats (TODO 3a) ---
  cat("\nStep 3: Gene-level ABC analysis (3a)...\n")
  stats_lines <- c(stats_lines, "=== 3a: DEG vs Invariant Delta-ABC ===")

  p1 <- plot_delta_abc_violin(gene_abc, "sum_delta_abc",
                              "Sum Delta-ABC per Gene", output_dir)
  p2 <- plot_delta_abc_violin(gene_abc, "max_delta_abc",
                              "Max Delta-ABC per Gene", output_dir)
  p3 <- plot_n_gained_lost_bar(gene_abc, output_dir)
  p4 <- plot_delta_abc_vs_lfc(gene_abc, output_dir)

  for (metric in c("sum_delta_abc", "max_delta_abc")) {
    stats_lines <- c(stats_lines, sprintf("\n--- %s ---", metric))
    up_vals   <- gene_abc %>% dplyr::filter(gene_class == "DEG_up") %>% dplyr::pull(!!sym(metric))
    down_vals <- gene_abc %>% dplyr::filter(gene_class == "DEG_down") %>% dplyr::pull(!!sym(metric))
    inv_vals  <- gene_abc %>% dplyr::filter(gene_class == "invariant") %>% dplyr::pull(!!sym(metric))

    stats_lines <- c(stats_lines,
      wilcox_line(metric, up_vals, inv_vals, "DEG_up", "invariant"),
      wilcox_line(metric, down_vals, inv_vals, "DEG_down", "invariant"),
      wilcox_line(metric, up_vals, down_vals, "DEG_up", "DEG_down")
    )

    if (length(up_vals) >= 3 && length(down_vals) >= 3 && length(inv_vals) >= 3) {
      kt <- kruskal.test(list(DEG_up = up_vals, DEG_down = down_vals, invariant = inv_vals))
      stats_lines <- c(stats_lines,
        sprintf("  Kruskal-Wallis: chi-sq = %.2f, df = %d, p = %.2e",
                kt$statistic, kt$parameter, kt$p.value))
    }
  }

  degs_abc <- gene_abc %>% dplyr::filter(gene_class != "invariant")
  if (nrow(degs_abc) >= 5) {
    rho <- cor.test(degs_abc$sum_delta_abc, degs_abc$log2FoldChange, method = "spearman")
    stats_lines <- c(stats_lines,
      sprintf("\nSpearman: sum_delta_abc vs log2FC (DEGs only, n=%d)", nrow(degs_abc)),
      sprintf("  rho = %.4f, p = %.2e", rho$estimate, rho$p.value))
  }

  for (gc in c("DEG_up", "DEG_down", "invariant")) {
    gc_df <- gene_abc %>% dplyr::filter(gene_class == gc)
    stats_lines <- c(stats_lines,
      sprintf("  %s: n_gained>0 = %d/%d (%.1f%%), n_lost>0 = %d/%d (%.1f%%)",
              gc,
              sum(gc_df$n_gained > 0), nrow(gc_df), 100 * sum(gc_df$n_gained > 0) / max(1, nrow(gc_df)),
              sum(gc_df$n_lost > 0), nrow(gc_df), 100 * sum(gc_df$n_lost > 0) / max(1, nrow(gc_df))))
  }

  # --- Step 4: Enhancer-level ABC pairs join ---
  cat("\nStep 4: Loading enhancer-level ABC pairs...\n")
  abc_pairs <- load_abc_pairs(ABC_PAIRS_FILE)

  pairs_filtered <- abc_pairs %>%
    dplyr::filter(TargetGene %in% gene_tbl$gene_symbol) %>%
    dplyr::left_join(gene_tbl %>% dplyr::select(gene_symbol, gene_class),
                     by = c("TargetGene" = "gene_symbol"))

  pairs_filtered <- pairs_filtered %>%
    dplyr::mutate(
      log2fc_contact = compute_log2fc_contact(
        hic_contact_pl_scaled_adj_WT, hic_contact_pl_scaled_adj_KO
      )
    )

  cat(sprintf("  Filtered to %d E-G pairs for %d genes\n",
              nrow(pairs_filtered), length(unique(pairs_filtered$TargetGene))))
  stats_lines <- c(stats_lines,
    sprintf("\n=== 3b: Enhancer-Level Analysis ==="),
    sprintf("Enhancer-gene pairs (distal): %d for %d genes",
            nrow(pairs_filtered), length(unique(pairs_filtered$TargetGene))))

  # --- Step 5: SE classification ---
  cat("\nStep 5: Classifying enhancers as SE vs regular...\n")
  p60_gr    <- load_se_bed(SE_FILES$P60, "P60")
  encode_gr <- load_se_bed(SE_FILES$ENCODE, "ENCODE")
  se_combined <- c(p60_gr, encode_gr)

  enh_gr <- GRanges(
    seqnames = pairs_filtered$chr,
    ranges = IRanges(start = pairs_filtered$start, end = pairs_filtered$end)
  )

  se_hits <- findOverlaps(enh_gr, se_combined, ignore.strand = TRUE)

  se_idx <- tibble(pair_idx = queryHits(se_hits),
                   se_source = se_combined$se_source[subjectHits(se_hits)]) %>%
    dplyr::distinct(pair_idx, .keep_all = TRUE)

  pairs_filtered <- pairs_filtered %>%
    dplyr::mutate(pair_idx = row_number()) %>%
    dplyr::left_join(se_idx, by = "pair_idx") %>%
    dplyr::mutate(
      is_se = !is.na(se_source),
      se_source = ifelse(is.na(se_source), "none", se_source)
    ) %>%
    dplyr::select(-pair_idx)

  n_se <- sum(pairs_filtered$is_se)
  cat(sprintf("  SE overlap: %d / %d E-G pairs (%.1f%%)\n",
              n_se, nrow(pairs_filtered), 100 * n_se / nrow(pairs_filtered)))
  stats_lines <- c(stats_lines,
    sprintf("SE enhancer overlap: %d / %d pairs (%.1f%%)",
            n_se, nrow(pairs_filtered), 100 * n_se / nrow(pairs_filtered)))

  # --- Step 6: K27ac DiffBind join ---
  cat("\nStep 6: K27ac DiffBind classification...\n")
  diffbind_gr <- load_k27ac_diffbind(DIFFBIND_K27AC)

  k27ac_result <- classify_k27ac_with_fold(enh_gr, diffbind_gr)

  pairs_filtered <- pairs_filtered %>%
    dplyr::mutate(pair_idx = row_number()) %>%
    dplyr::left_join(k27ac_result %>% dplyr::rename(pair_idx = peak_idx),
                     by = "pair_idx") %>%
    dplyr::select(-pair_idx)

  # --- Step 7: Enhancer-level plots and stats (TODO 3b) ---
  cat("\nStep 7: Enhancer-level plots and statistics (3b)...\n")

  p5 <- plot_se_contact_violin(pairs_filtered, output_dir)
  p6 <- plot_se_proportion_bar(pairs_filtered, output_dir)
  p7 <- plot_k27ac_contact_scatter(pairs_filtered, output_dir)
  p8 <- plot_k27ac_class_violin(pairs_filtered, output_dir)
  p9 <- plot_se_k27ac_subclass_bar(pairs_filtered, output_dir)

  # SE proportion Fisher tests
  for (gc in c("DEG_up", "DEG_down")) {
    gc_se  <- sum(pairs_filtered$gene_class == gc & pairs_filtered$is_se)
    gc_reg <- sum(pairs_filtered$gene_class == gc & !pairs_filtered$is_se)
    inv_se  <- sum(pairs_filtered$gene_class == "invariant" & pairs_filtered$is_se)
    inv_reg <- sum(pairs_filtered$gene_class == "invariant" & !pairs_filtered$is_se)
    tbl <- matrix(c(gc_se, gc_reg, inv_se, inv_reg), nrow = 2)
    if (all(tbl > 0)) {
      ft <- fisher.test(tbl)
      stats_lines <- c(stats_lines,
        sprintf("\nSE proportion Fisher (%s vs invariant):", gc),
        sprintf("  %s: %d SE / %d regular (%.1f%%)", gc, gc_se, gc_reg,
                100 * gc_se / max(1, gc_se + gc_reg)),
        sprintf("  invariant: %d SE / %d regular (%.1f%%)", inv_se, inv_reg,
                100 * inv_se / max(1, inv_se + inv_reg)),
        sprintf("  OR = %.3f [%.3f - %.3f], p = %.2e",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
    }
  }

  # SE vs regular contact change (DEGs only)
  deg_pairs <- pairs_filtered %>%
    dplyr::filter(gene_class != "invariant", !is.na(log2fc_contact), is.finite(log2fc_contact))
  se_contact  <- deg_pairs %>% dplyr::filter(is_se) %>% dplyr::pull(log2fc_contact)
  reg_contact <- deg_pairs %>% dplyr::filter(!is_se) %>% dplyr::pull(log2fc_contact)
  if (length(se_contact) >= 3 && length(reg_contact) >= 3) {
    stats_lines <- c(stats_lines,
      wilcox_line("log2fc_contact (DEGs only)", se_contact, reg_contact, "SE", "regular"))
  }

  # K27ac fold vs contact (DEGs only)
  deg_k27ac <- pairs_filtered %>%
    dplyr::filter(gene_class != "invariant", !is.na(k27ac_fold),
                  !is.na(log2fc_contact), is.finite(log2fc_contact), is.finite(k27ac_fold),
                  k27ac_class != "no_diffbind_peak")
  if (nrow(deg_k27ac) >= 10) {
    rho_k <- cor.test(deg_k27ac$k27ac_fold, deg_k27ac$log2fc_contact, method = "spearman")
    stats_lines <- c(stats_lines,
      sprintf("\nSpearman: K27ac Fold vs log2fc_contact (DEGs, n=%d)", nrow(deg_k27ac)),
      sprintf("  rho = %.4f, p = %.2e", rho_k$estimate, rho_k$p.value))
  }

  # K27ac concordance at SE contacts
  se_deg_pairs <- pairs_filtered %>%
    dplyr::filter(is_se, gene_class %in% c("DEG_up", "DEG_down"),
                  k27ac_class %in% c("gained_k27ac", "lost_k27ac")) %>%
    dplyr::mutate(
      concordant = (gene_class == "DEG_up" & k27ac_class == "gained_k27ac") |
                   (gene_class == "DEG_down" & k27ac_class == "lost_k27ac")
    )
  if (nrow(se_deg_pairs) > 0) {
    tbl <- table(
      deg_direction = se_deg_pairs$gene_class,
      k27ac_direction = se_deg_pairs$k27ac_class
    )
    if (nrow(tbl) == 2 && ncol(tbl) == 2) {
      ft <- fisher.test(tbl)
      stats_lines <- c(stats_lines,
        sprintf("\nK27ac concordance at SE contacts (Fisher):"),
        sprintf("  OR = %.3f [%.3f - %.3f], p = %.2e",
                ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value),
        sprintf("  Concordant: %d / %d (%.1f%%)",
                sum(se_deg_pairs$concordant), nrow(se_deg_pairs),
                100 * sum(se_deg_pairs$concordant) / nrow(se_deg_pairs)))
    }
  }

  # K27ac class pairwise Wilcoxon on log2fc_contact
  for (pair in list(c("gained_k27ac", "lost_k27ac"),
                    c("gained_k27ac", "stable_k27ac"),
                    c("lost_k27ac", "stable_k27ac"))) {
    a <- pairs_filtered %>%
      dplyr::filter(k27ac_class == pair[1], !is.na(log2fc_contact), is.finite(log2fc_contact)) %>%
      dplyr::pull(log2fc_contact)
    b <- pairs_filtered %>%
      dplyr::filter(k27ac_class == pair[2], !is.na(log2fc_contact), is.finite(log2fc_contact)) %>%
      dplyr::pull(log2fc_contact)
    stats_lines <- c(stats_lines,
      wilcox_line("log2fc_contact", a, b, pair[1], pair[2]))
  }

  # --- Step 8: Save tables ---
  cat("\nStep 8: Saving outputs...\n")
  write_tsv(gene_abc, file.path(tables_dir, "deg_abc_gene_level.tsv"))

  pairs_out <- pairs_filtered %>%
    dplyr::select(chr, start, end, TargetGene, gene_class, delta_ABC, log2fc_contact,
                  class, is_se, se_source, k27ac_class, k27ac_fold, k27ac_fdr,
                  hic_contact_pl_scaled_adj_WT, hic_contact_pl_scaled_adj_KO)
  write_tsv(pairs_out, file.path(tables_dir, "deg_abc_pairs_classified.tsv"))

  se_summary <- pairs_filtered %>%
    dplyr::group_by(gene_class, is_se, k27ac_class) %>%
    dplyr::summarise(
      n_pairs = n(),
      mean_delta_abc = mean(delta_ABC, na.rm = TRUE),
      median_delta_abc = median(delta_ABC, na.rm = TRUE),
      mean_log2fc_contact = mean(log2fc_contact, na.rm = TRUE),
      .groups = "drop"
    )
  write_tsv(se_summary, file.path(tables_dir, "se_contact_summary.tsv"))

  writeLines(stats_lines, file.path(stats_dir, "abc_contact_statistics.txt"))

  cat(sprintf("  Wrote 3 tables + statistics file\n"))
  cat(sprintf("\n========== Done: %s ==========\n", config$label))
}

# =============================================================================
# 8. DISPATCH
# =============================================================================

if (parsed$timepoint == "both") {
  timepoints <- c("late", "early")
} else {
  timepoints <- parsed$timepoint
}

for (tp in timepoints) {
  run_analysis(tp)
}
