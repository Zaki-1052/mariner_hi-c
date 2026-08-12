# scripts/se_deg_proximity.R
#
# 1a. SE-DEG Proximity Analysis
#
# For each DEG, identifies superenhancers within a 2 Mb window and compares
# against invariant genes matched by expression level and chromosome.
# Quantifies SE proximity, K27ac sub-classification, and prepares APA inputs.
# Runs two enhancer comparison arms: ChromHMM active enhancers and K27ac ChIP peaks.
#
# Usage:
#   Rscript scripts/se_deg_proximity.R --timepoint late
#   Rscript scripts/se_deg_proximity.R --timepoint early
#   Rscript scripts/se_deg_proximity.R --timepoint both
#
# Output:
#   output/superenhancer_analysis/1a_se_deg_proximity/{timepoint}/

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

SE_2MB_WINDOW <- 2e6
DEG_PADJ <- 0.05
DEG_LFC <- 0.3
INVARIANT_RATIO <- 3

TIMEPOINT_CONFIG <- list(
  late = list(
    rna_file = file.path(BASE_DIR, "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"),
    k27ac_chip = file.path(BASE_DIR, "peaks/beds/H3K27acCerebellumLate2.bed"),
    label = "Late (Adult)"
  ),
  early = list(
    rna_file = file.path(BASE_DIR, "tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx"),
    k27ac_chip = file.path(BASE_DIR, "peaks/beds/H3K27acCerebellumEarly2.bed"),
    label = "Early (P12)"
  )
)

SE_FILES <- list(
  P60    = file.path(BASE_DIR, "peaks/Superenhancers_P60.bed"),
  ENCODE = file.path(BASE_DIR, "peaks/Superenhancers_encode.bed")
)

CHROMHMM_ENH <- file.path(BASE_DIR, "peaks/251230-Challana-EmissionState12-activeenhancer.bed")
DIFFBIND_K27AC <- file.path(BASE_DIR, "peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt")

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
# 3. PROXIMITY FUNCTIONS
# =============================================================================

build_gene_windows <- function(gene_tbl, window_size = SE_2MB_WINDOW) {
  gr <- GRanges(
    seqnames = gene_tbl$chr,
    ranges = IRanges(
      start = pmax(1L, as.integer(gene_tbl$tss - window_size)),
      end = as.integer(gene_tbl$tss + window_size)
    ),
    gene_symbol = gene_tbl$gene_symbol,
    gene_class = gene_tbl$gene_class
  )
  cat(sprintf("  Built %d gene windows (±%g Mb)\n", length(gr), window_size / 1e6))
  return(gr)
}

compute_proximity <- function(gene_windows_gr, feature_gr, feature_label) {
  hits <- findOverlaps(gene_windows_gr, feature_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    cat(sprintf("  No %s found within gene windows\n", feature_label))
    return(tibble(
      gene_symbol = character(), gene_class = character(),
      feature_idx = integer(), feature_label = character(),
      distance_to_tss = numeric()
    ))
  }

  gene_tss <- (start(gene_windows_gr) + end(gene_windows_gr)) / 2
  feature_mid <- (start(feature_gr) + end(feature_gr)) / 2

  pairs <- tibble(
    gene_symbol = gene_windows_gr$gene_symbol[queryHits(hits)],
    gene_class = gene_windows_gr$gene_class[queryHits(hits)],
    feature_idx = subjectHits(hits),
    feature_label = feature_label,
    distance_to_tss = abs(gene_tss[queryHits(hits)] - feature_mid[subjectHits(hits)])
  )

  cat(sprintf("  %s proximity: %d gene-feature pairs (%d unique genes)\n",
              feature_label, nrow(pairs), length(unique(pairs$gene_symbol))))

  return(pairs)
}

summarise_proximity <- function(pairs_df) {
  pairs_df %>%
    dplyr::group_by(gene_symbol, gene_class, feature_label) %>%
    dplyr::summarise(
      n_features = n(),
      min_distance = min(distance_to_tss),
      median_distance = median(distance_to_tss),
      .groups = "drop"
    )
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

plot_se_count_violin <- function(summary_df, se_source, output_dir) {
  plot_df <- summary_df %>%
    dplyr::filter(feature_label == se_source)

  if (nrow(plot_df) < 3) {
    cat(sprintf("  Skipping violin for %s — too few data points\n", se_source))
    return(NULL)
  }

  comparisons <- list()
  classes <- unique(plot_df$gene_class)
  if ("DEG_up" %in% classes && "invariant" %in% classes)
    comparisons <- c(comparisons, list(c("DEG_up", "invariant")))
  if ("DEG_down" %in% classes && "invariant" %in% classes)
    comparisons <- c(comparisons, list(c("DEG_down", "invariant")))
  if ("DEG_up" %in% classes && "DEG_down" %in% classes)
    comparisons <- c(comparisons, list(c("DEG_up", "DEG_down")))

  p <- ggplot(plot_df, aes(x = gene_class, y = n_features, fill = gene_class)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    scale_fill_manual(values = DEG_CLASS_COLORS) +
    labs(
      title = sprintf("%s SEs Within 2 Mb of Genes", se_source),
      x = "Gene Class", y = "Number of SEs Within 2 Mb"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  if (length(comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                                 label = "p.format", tip.length = 0.02)
  }

  save_multiformat_ggplot(p, file.path(output_dir, sprintf("se_count_violin_%s", tolower(se_source))),
                          width = 6, height = 6)
  return(p)
}

plot_k27ac_stacked_bar <- function(pairs_with_k27ac, se_source, output_dir) {
  plot_df <- pairs_with_k27ac %>%
    dplyr::filter(feature_label == se_source, gene_class != "invariant") %>%
    dplyr::count(gene_class, k27ac_class)

  if (nrow(plot_df) == 0) return(NULL)

  p <- ggplot(plot_df, aes(x = gene_class, y = n, fill = k27ac_class)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = K27AC_CLASS_COLORS,
                      labels = c("gained_k27ac" = "Gained K27ac",
                                 "lost_k27ac" = "Lost K27ac",
                                 "stable_k27ac" = "Stable K27ac",
                                 "no_diffbind_peak" = "No DiffBind Peak")) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(labels = c("DEG_up" = "Upregulated DEGs",
                                "DEG_down" = "Downregulated DEGs")) +
    labs(
      title = sprintf("K27ac Status of %s SEs Near DEGs", se_source),
      x = "", y = "Proportion of DEG-SE Pairs", fill = "K27ac Change"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  save_multiformat_ggplot(p, file.path(output_dir, sprintf("k27ac_subclass_%s", tolower(se_source))),
                          width = 7, height = 6)
  return(p)
}

plot_se_vs_lfc_scatter <- function(summary_df, rna_df, se_source, output_dir) {
  plot_df <- summary_df %>%
    dplyr::filter(feature_label == se_source, gene_class != "invariant") %>%
    dplyr::inner_join(rna_df %>% dplyr::select(gene_symbol, log2FoldChange),
                      by = "gene_symbol")

  if (nrow(plot_df) < 5) return(NULL)

  rho <- cor.test(plot_df$n_features, plot_df$log2FoldChange, method = "spearman")

  p <- ggplot(plot_df, aes(x = n_features, y = log2FoldChange)) +
    geom_point(aes(color = gene_class), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = DEG_CLASS_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.3f\np = %.2e", rho$estimate, rho$p.value),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(
      title = sprintf("Gene Expression vs %s SE Proximity", se_source),
      x = sprintf("Number of %s SEs Within 2 Mb", se_source),
      y = "log2FoldChange (mut/ctrl)",
      color = "Gene Class"
    ) +
    theme_classic(base_size = 14)

  save_multiformat_ggplot(p, file.path(output_dir, sprintf("se_vs_lfc_%s", tolower(se_source))),
                          width = 7, height = 6)
  return(p)
}

plot_se_vs_enhancer_comparison <- function(summary_all, output_dir) {
  compare_df <- summary_all %>%
    dplyr::group_by(gene_symbol, gene_class, feature_label) %>%
    dplyr::summarise(n_features = sum(n_features), .groups = "drop") %>%
    dplyr::group_by(gene_class, feature_label) %>%
    dplyr::summarise(
      mean_count = mean(n_features),
      se_count = sd(n_features) / sqrt(n()),
      .groups = "drop"
    )

  if (nrow(compare_df) == 0) return(NULL)

  p <- ggplot(compare_df, aes(x = gene_class, y = mean_count, fill = feature_label)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count),
                  position = position_dodge(width = 0.7), width = 0.2) +
    labs(
      title = "SEs vs Regular Enhancers Near Genes",
      x = "Gene Class", y = "Mean Feature Count Within 2 Mb", fill = "Feature Type"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1))

  save_multiformat_ggplot(p, file.path(output_dir, "se_vs_enhancer_comparison"),
                          width = 8, height = 6)
  return(p)
}

# =============================================================================
# 5. APA PREPARATION
# =============================================================================

write_apa_inputs <- function(se_pairs, enh_pairs, invariant_pairs,
                             gene_coords, se_gr, timepoint, apa_dir) {
  dir.create(apa_dir, recursive = TRUE, showWarnings = FALSE)

  write_bedpe <- function(pairs_df, feature_gr, filename) {
    if (nrow(pairs_df) == 0) {
      cat(sprintf("  Skipping %s — no pairs\n", filename))
      return()
    }

    coords <- gene_coords %>% dplyr::select(gene_symbol, chr, tss)
    bedpe <- pairs_df %>%
      dplyr::inner_join(coords, by = "gene_symbol") %>%
      dplyr::mutate(
        chr1 = as.character(seqnames(feature_gr))[feature_idx],
        x1 = start(feature_gr)[feature_idx],
        x2 = end(feature_gr)[feature_idx],
        chr2 = chr,
        y1 = as.integer(tss - 5000),
        y2 = as.integer(tss + 5000)
      ) %>%
      dplyr::select(chr1, x1, x2, chr2, y1, y2)

    write_tsv(bedpe, file.path(apa_dir, filename), col_names = FALSE)
    cat(sprintf("  Wrote %s (%d pairs)\n", filename, nrow(bedpe)))
  }

  write_bedpe(se_pairs %>% dplyr::filter(feature_label == "P60"),
              se_gr[se_gr$se_source == "P60"],
              sprintf("%s_deg_se_pairs.bedpe", timepoint))
  write_bedpe(invariant_pairs %>% dplyr::filter(feature_label == "P60"),
              se_gr[se_gr$se_source == "P60"],
              sprintf("%s_invariant_se_pairs.bedpe", timepoint))

  launcher <- file.path(apa_dir, "apa_se_launcher.R")
  writeLines(c(
    "# APA launcher for SE-DEG contact frequency analysis",
    "# Run on HPC (requires .hic files and mariner::pullHicMatrices)",
    "#",
    "# Template based on scripts/apa_analysis.R",
    "# Input BEDPEs generated by se_deg_proximity.R",
    "#",
    sprintf("# Generated: %s", Sys.time()),
    "",
    "# TODO: adapt apa_analysis.R to load these BEDPEs as GInteractions",
    "# and run pullHicMatrices + APA aggregation per gene class"
  ), launcher)
  cat(sprintf("  Wrote APA launcher stub: %s\n", launcher))
}

# =============================================================================
# 6. MAIN ANALYSIS
# =============================================================================

run_analysis <- function(timepoint) {
  config <- TIMEPOINT_CONFIG[[timepoint]]
  cat(sprintf("\n========== 1a. SE-DEG Proximity: %s ==========\n\n", config$label))

  output_dir <- file.path(BASE_DIR, "loops/output/superenhancer_analysis/1a_se_deg_proximity", timepoint)
  tables_dir <- file.path(output_dir, "tables")
  stats_dir <- file.path(output_dir, "statistics")
  apa_dir <- file.path(BASE_DIR, "loops/output/superenhancer_analysis/1a_se_deg_proximity/apa_inputs")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Load RNA-seq ---
  cat("Step 1: Loading RNA-seq data...\n")
  rna <- load_rna_for_se(config$rna_file, padj_thresh = DEG_PADJ, lfc_thresh = DEG_LFC)
  degs <- rna$degs
  all_genes <- rna$all_genes

  # --- Gene coordinates ---
  cat("\nStep 2: Mapping gene coordinates...\n")
  all_symbols <- unique(c(degs$gene_symbol, all_genes$gene_symbol))
  gene_coords <- get_gene_coordinates(all_symbols)

  # --- Invariant genes ---
  cat("\nStep 3: Selecting invariant genes...\n")
  invariants <- select_invariant_genes(degs, all_genes, gene_coords,
                                        n_per_deg = INVARIANT_RATIO)

  # --- Build gene table ---
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

  cat(sprintf("\n  Gene table: %d DEG_up, %d DEG_down, %d invariant\n",
              sum(gene_tbl$gene_class == "DEG_up"),
              sum(gene_tbl$gene_class == "DEG_down"),
              sum(gene_tbl$gene_class == "invariant")))

  # --- Load features ---
  cat("\nStep 4: Loading SE and enhancer references...\n")
  se_gr <- load_all_ses(SE_FILES$P60, SE_FILES$ENCODE)
  p60_gr <- se_gr[se_gr$se_source == "P60"]
  encode_gr <- se_gr[se_gr$se_source == "ENCODE"]
  chromhmm_enh_gr <- load_enhancer_bed(CHROMHMM_ENH, "ChromHMM_ActiveEnhancer")
  k27ac_chip_gr <- load_enhancer_bed(config$k27ac_chip, "K27ac_ChIP")

  # --- Load K27ac DiffBind ---
  cat("\nStep 5: Loading K27ac DiffBind for SE classification...\n")
  diffbind_gr <- load_k27ac_diffbind(DIFFBIND_K27AC)

  # --- Proximity analysis ---
  cat("\nStep 6: Computing proximity...\n")
  gene_windows_gr <- build_gene_windows(gene_tbl)

  pairs_p60 <- compute_proximity(gene_windows_gr, p60_gr, "P60")
  pairs_encode <- compute_proximity(gene_windows_gr, encode_gr, "ENCODE")
  pairs_chromhmm <- compute_proximity(gene_windows_gr, chromhmm_enh_gr, "ChromHMM_ActiveEnhancer")
  pairs_k27ac <- compute_proximity(gene_windows_gr, k27ac_chip_gr, "K27ac_ChIP")

  all_pairs <- bind_rows(pairs_p60, pairs_encode, pairs_chromhmm, pairs_k27ac)

  # --- K27ac sub-classification on SEs ---
  cat("\nStep 7: K27ac sub-classification of SEs...\n")
  k27ac_p60 <- classify_k27ac_change(p60_gr, diffbind_gr)
  k27ac_encode <- classify_k27ac_change(encode_gr, diffbind_gr)

  pairs_p60 <- pairs_p60 %>%
    dplyr::left_join(k27ac_p60 %>% dplyr::rename(feature_idx = peak_idx),
                     by = "feature_idx")
  pairs_encode <- pairs_encode %>%
    dplyr::left_join(k27ac_encode %>% dplyr::rename(feature_idx = peak_idx),
                     by = "feature_idx")

  se_pairs_with_k27ac <- bind_rows(pairs_p60, pairs_encode)

  # --- Summary per gene ---
  cat("\nStep 8: Summarizing proximity per gene...\n")
  summary_all <- summarise_proximity(all_pairs)

  # --- Statistics ---
  cat("\nStep 9: Running statistical tests...\n")
  stats_lines <- c(
    sprintf("SE-DEG Proximity Analysis"),
    sprintf("Timepoint: %s", config$label),
    sprintf("Date: %s\n", Sys.time()),
    sprintf("Genes: %d DEG_up, %d DEG_down, %d invariant\n",
            sum(gene_tbl$gene_class == "DEG_up"),
            sum(gene_tbl$gene_class == "DEG_down"),
            sum(gene_tbl$gene_class == "invariant"))
  )

  for (fl in c("P60", "ENCODE", "ChromHMM_ActiveEnhancer", "K27ac_ChIP")) {
    fl_summary <- summary_all %>% dplyr::filter(feature_label == fl)
    if (nrow(fl_summary) == 0) next

    stats_lines <- c(stats_lines, sprintf("\n--- %s ---", fl))

    for (gc in c("DEG_up", "DEG_down")) {
      gc_vals <- fl_summary %>% dplyr::filter(gene_class == gc) %>% dplyr::pull(n_features)
      inv_vals <- fl_summary %>% dplyr::filter(gene_class == "invariant") %>% dplyr::pull(n_features)

      if (length(gc_vals) >= 3 && length(inv_vals) >= 3) {
        wt <- wilcox.test(gc_vals, inv_vals)
        stats_lines <- c(stats_lines,
          sprintf("  %s vs invariant: median %g vs %g, W = %g, p = %.2e",
                  gc, median(gc_vals), median(inv_vals), wt$statistic, wt$p.value))
      }
    }
  }

  # Concordance Fisher test
  if (nrow(se_pairs_with_k27ac) > 0) {
    deg_se <- se_pairs_with_k27ac %>%
      dplyr::filter(gene_class %in% c("DEG_up", "DEG_down"),
                    k27ac_class %in% c("gained_k27ac", "lost_k27ac")) %>%
      dplyr::mutate(
        concordant = (gene_class == "DEG_up" & k27ac_class == "gained_k27ac") |
                     (gene_class == "DEG_down" & k27ac_class == "lost_k27ac")
      )

    if (nrow(deg_se) > 0) {
      tbl <- table(
        deg_direction = deg_se$gene_class,
        k27ac_direction = deg_se$k27ac_class
      )
      if (nrow(tbl) == 2 && ncol(tbl) == 2) {
        ft <- fisher.test(tbl)
        stats_lines <- c(stats_lines,
          sprintf("\nDEG-K27ac Concordance Fisher Test:"),
          sprintf("  OR = %.3f [%.3f - %.3f], p = %.2e",
                  ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value),
          sprintf("  Concordant pairs: %d / %d (%.1f%%)",
                  sum(deg_se$concordant), nrow(deg_se),
                  100 * sum(deg_se$concordant) / nrow(deg_se)))
      }
    }
  }

  writeLines(stats_lines, file.path(stats_dir, "se_proximity_statistics.txt"))
  cat(sprintf("  Wrote statistics file\n"))

  # --- Save tables ---
  cat("\nStep 10: Saving tables...\n")
  write_tsv(all_pairs, file.path(tables_dir, "gene_feature_proximity.tsv"))
  write_tsv(summary_all, file.path(tables_dir, "proximity_summary_by_gene.tsv"))
  write_tsv(se_pairs_with_k27ac, file.path(tables_dir, "se_pairs_with_k27ac.tsv"))

  gene_summary <- gene_tbl %>%
    dplyr::select(gene_symbol, gene_class, log2FoldChange, padj, baseMean, chr, tss)
  write_tsv(gene_summary, file.path(tables_dir, "gene_table.tsv"))
  cat(sprintf("  Wrote 4 TSV tables\n"))

  # --- Plots ---
  cat("\nStep 11: Generating plots...\n")
  plot_se_count_violin(summary_all, "P60", output_dir)
  plot_se_count_violin(summary_all, "ENCODE", output_dir)
  plot_k27ac_stacked_bar(se_pairs_with_k27ac, "P60", output_dir)
  plot_k27ac_stacked_bar(se_pairs_with_k27ac, "ENCODE", output_dir)
  plot_se_vs_lfc_scatter(summary_all, degs, "P60", output_dir)
  plot_se_vs_lfc_scatter(summary_all, degs, "ENCODE", output_dir)
  plot_se_vs_enhancer_comparison(summary_all, output_dir)

  # --- APA prep ---
  cat("\nStep 12: Writing APA inputs...\n")
  deg_se_pairs <- se_pairs_with_k27ac %>% dplyr::filter(gene_class %in% c("DEG_up", "DEG_down"))
  inv_se_pairs <- se_pairs_with_k27ac %>% dplyr::filter(gene_class == "invariant")
  write_apa_inputs(deg_se_pairs, pairs_chromhmm, inv_se_pairs,
                   gene_coords, se_gr, timepoint, apa_dir)

  cat(sprintf("\n========== Done: %s ==========\n", config$label))
}

# =============================================================================
# 7. DISPATCH
# =============================================================================

if (parsed$timepoint == "both") {
  timepoints <- c("late", "early")
} else {
  timepoints <- parsed$timepoint
}

for (tp in timepoints) {
  run_analysis(tp)
}
