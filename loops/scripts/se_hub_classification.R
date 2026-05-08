# scripts/se_hub_classification.R
#
# 1c. SE Hub Gene Classification
#
# Identifies stripes whose anchors overlap superenhancers ("SE-hub stripes"),
# extracts the body genes contacted by those stripes, and tests whether genes
# within the same SE hub are coordinately regulated.
#
# SE hub model: the SE at the stripe anchor modifies multiple promoters
# along the stripe body ("enhancers modifying multiple promoters").
#
# Usage:
#   Rscript scripts/se_hub_classification.R --timepoint late
#   Rscript scripts/se_hub_classification.R --timepoint early
#   Rscript scripts/se_hub_classification.R --timepoint both
#
# Output:
#   output/superenhancer_analysis/1c_se_hub_genes/{timepoint}/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

BASE_DIR <- getwd()

source(file.path(BASE_DIR, "loops/scripts/utils/multi_format_output.R"))
source(file.path(BASE_DIR, "loops/scripts/utils/se_utils.R"))

DEG_PADJ <- 0.05
DEG_LFC <- 0.3
N_PERMUTATIONS <- 10000
PERM_SEED <- 42

TIMEPOINT_CONFIG <- list(
  late = list(
    rna_file = file.path(BASE_DIR, "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"),
    stripe_files = list(
      allsig = file.path(BASE_DIR, "stripes/stripenn/outputs/250402/250402_stripes_allsig.bedpe")
    ),
    label = "Late (Adult)"
  ),
  early = list(
    rna_file = file.path(BASE_DIR, "tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx"),
    stripe_files = list(
      allsig = file.path(BASE_DIR, "stripes/stripenn/outputs/250831/250831_stripes_allsig.bedpe"),
      concordant = file.path(BASE_DIR, "stripes/stripenn/outputs/250831/250831_stripes_concordant.bedpe")
    ),
    label = "Early (P12)"
  )
)

SE_FILES <- list(
  P60    = file.path(BASE_DIR, "peaks/Superenhancers_P60.bed"),
  ENCODE = file.path(BASE_DIR, "peaks/Superenhancers_encode.bed")
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
# 3. CORE FUNCTIONS
# =============================================================================

load_stripes_for_hub <- function(stripe_file) {
  if (!file.exists(stripe_file)) stop(sprintf("Stripe file not found: %s", stripe_file))

  stripes <- read_tsv(stripe_file, show_col_types = FALSE)

  required_cols <- c("chr1", "x1", "x2", "direction", "logFC", "FDR",
                     "body_genes", "body_gene_count")
  missing <- setdiff(required_cols, colnames(stripes))
  if (length(missing) > 0) {
    stop(sprintf("Missing stripe columns: %s", paste(missing, collapse = ", ")))
  }

  stripes <- stripes %>%
    dplyr::filter(chr1 %in% VALID_CHROMS) %>%
    dplyr::mutate(stripe_id = paste0("stripe_", row_number()))

  cat(sprintf("  Loaded %d stripes from %s\n", nrow(stripes), basename(stripe_file)))
  return(stripes)
}

identify_se_hubs <- function(stripes, p60_gr, encode_gr) {
  anchor_gr <- GRanges(
    seqnames = stripes$chr1,
    ranges = IRanges(start = stripes$x1, end = stripes$x2)
  )

  hits_p60 <- findOverlaps(anchor_gr, p60_gr, ignore.strand = TRUE)
  hits_encode <- findOverlaps(anchor_gr, encode_gr, ignore.strand = TRUE)

  stripes$se_hub_P60 <- FALSE
  stripes$se_hub_ENCODE <- FALSE
  stripes$se_hub_P60[unique(queryHits(hits_p60))] <- TRUE
  stripes$se_hub_ENCODE[unique(queryHits(hits_encode))] <- TRUE
  stripes$se_hub_any <- stripes$se_hub_P60 | stripes$se_hub_ENCODE

  cat(sprintf("  SE-hub stripes: %d P60, %d ENCODE, %d any (%d total stripes)\n",
              sum(stripes$se_hub_P60), sum(stripes$se_hub_ENCODE),
              sum(stripes$se_hub_any), nrow(stripes)))

  return(stripes)
}

extract_hub_genes <- function(stripes) {
  hub_stripes <- stripes %>%
    dplyr::filter(se_hub_any == TRUE,
                  !is.na(body_genes), body_genes != "")

  if (nrow(hub_stripes) == 0) {
    cat("  No SE-hub stripes with body genes found\n")
    return(tibble())
  }

  hub_gene_assign <- hub_stripes %>%
    dplyr::select(stripe_id, direction, logFC, FDR, se_hub_P60, se_hub_ENCODE,
                  body_genes, body_gene_count) %>%
    dplyr::mutate(gene_list = strsplit(body_genes, ",\\s*")) %>%
    tidyr::unnest(gene_list) %>%
    dplyr::rename(gene_symbol = gene_list) %>%
    dplyr::filter(gene_symbol != "")

  cat(sprintf("  Extracted %d gene-hub assignments from %d SE-hub stripes\n",
              nrow(hub_gene_assign), nrow(hub_stripes)))
  cat(sprintf("    - Unique genes: %d\n", length(unique(hub_gene_assign$gene_symbol))))

  return(hub_gene_assign)
}

compute_hub_coordination <- function(hub_gene_assign, rna_df) {
  joined <- hub_gene_assign %>%
    dplyr::left_join(
      rna_df %>% dplyr::select(gene_symbol, log2FoldChange, padj, baseMean),
      by = "gene_symbol"
    )

  hub_summary <- joined %>%
    dplyr::group_by(stripe_id, direction, logFC) %>%
    dplyr::summarise(
      hub_size = n_distinct(gene_symbol),
      n_with_rna = sum(!is.na(log2FoldChange)),
      n_degs = sum(!is.na(padj) & padj < DEG_PADJ & abs(log2FoldChange) > DEG_LFC, na.rm = TRUE),
      pct_degs = ifelse(n_with_rna > 0, n_degs / n_with_rna, NA_real_),
      mean_gene_lfc = mean(log2FoldChange, na.rm = TRUE),
      median_gene_lfc = median(log2FoldChange, na.rm = TRUE),
      concordant_n = sum(sign(log2FoldChange) == sign(logFC[1]), na.rm = TRUE),
      concordant_pct = ifelse(n_with_rna > 0,
                              concordant_n / n_with_rna,
                              NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::rename(logFC_stripe = logFC)

  cat(sprintf("  Hub coordination: %d hubs, median size = %d, median concordance = %.1f%%\n",
              nrow(hub_summary), median(hub_summary$hub_size),
              100 * median(hub_summary$concordant_pct, na.rm = TRUE)))

  return(list(hub_summary = hub_summary, joined = joined))
}

run_permutation_test <- function(hub_gene_assign, rna_df, n_perm = N_PERMUTATIONS,
                                 seed = PERM_SEED) {
  joined <- hub_gene_assign %>%
    dplyr::left_join(
      rna_df %>% dplyr::select(gene_symbol, log2FoldChange),
      by = "gene_symbol"
    ) %>%
    dplyr::filter(!is.na(log2FoldChange))

  if (nrow(joined) < 10) {
    cat("  Too few gene-hub pairs for permutation test\n")
    return(NULL)
  }

  observed_concordance <- joined %>%
    dplyr::group_by(stripe_id, logFC) %>%
    dplyr::summarise(
      concordant_pct = mean(sign(log2FoldChange) == sign(logFC[1])),
      .groups = "drop"
    ) %>%
    dplyr::pull(concordant_pct) %>%
    mean(na.rm = TRUE)

  set.seed(seed)
  perm_concordances <- numeric(n_perm)
  all_lfcs <- joined$log2FoldChange

  for (i in seq_len(n_perm)) {
    shuffled <- joined %>%
      dplyr::mutate(log2FoldChange = sample(all_lfcs))

    perm_concordances[i] <- shuffled %>%
      dplyr::group_by(stripe_id, logFC) %>%
      dplyr::summarise(
        concordant_pct = mean(sign(log2FoldChange) == sign(logFC[1])),
        .groups = "drop"
      ) %>%
      dplyr::pull(concordant_pct) %>%
      mean(na.rm = TRUE)
  }

  p_value <- mean(perm_concordances >= observed_concordance)
  z_score <- (observed_concordance - mean(perm_concordances)) / sd(perm_concordances)

  cat(sprintf("  Permutation test (n=%d): observed = %.3f, null mean = %.3f, z = %.2f, p = %.4f\n",
              n_perm, observed_concordance, mean(perm_concordances), z_score, p_value))

  return(list(
    observed = observed_concordance,
    null_mean = mean(perm_concordances),
    null_sd = sd(perm_concordances),
    z_score = z_score,
    p_value = p_value,
    perm_distribution = perm_concordances
  ))
}

find_multi_hub_genes <- function(hub_gene_assign, rna_df) {
  multi_hub <- hub_gene_assign %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::filter(n_distinct(stripe_id) >= 2) %>%
    dplyr::summarise(
      n_hubs = n_distinct(stripe_id),
      hub_directions = paste(sort(unique(direction)), collapse = ","),
      direction_mixed = n_distinct(direction) > 1,
      mean_stripe_lfc = mean(logFC),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      rna_df %>% dplyr::select(gene_symbol, log2FoldChange, padj, baseMean),
      by = "gene_symbol"
    ) %>%
    dplyr::arrange(desc(n_hubs))

  cat(sprintf("  Multi-hub genes: %d genes in >= 2 SE-hub stripes (max = %d hubs)\n",
              nrow(multi_hub), max(multi_hub$n_hubs, 0)))

  return(multi_hub)
}

# =============================================================================
# 4. PLOTTING FUNCTIONS
# =============================================================================

plot_hub_size_histogram <- function(stripes, output_dir) {
  diff_stripes <- stripes %>%
    dplyr::filter(direction %in% c("gained", "lost"),
                  !is.na(body_gene_count), body_gene_count > 0)

  if (nrow(diff_stripes) < 3) return(NULL)

  plot_df <- diff_stripes %>%
    dplyr::mutate(hub_status = ifelse(se_hub_any, "SE-hub", "Non-SE-hub"))

  p <- ggplot(plot_df, aes(x = body_gene_count, fill = hub_status)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
    facet_wrap(~direction, ncol = 1, scales = "free_y",
               labeller = labeller(direction = c("gained" = "Gained Stripes",
                                                 "lost" = "Lost Stripes"))) +
    scale_fill_manual(values = c("SE-hub" = "#e0a730", "Non-SE-hub" = "#cccccc")) +
    labs(
      title = "Hub Size: SE-Hub vs Non-SE-Hub Stripes",
      x = "Number of Body Genes", y = "Count", fill = ""
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(output_dir, "hub_size_histogram"),
                          width = 7, height = 8)
  return(p)
}

plot_deg_enrichment_violin <- function(hub_summary, stripes, rna_df, output_dir) {
  non_hub_stripes <- stripes %>%
    dplyr::filter(!se_hub_any, direction %in% c("gained", "lost"),
                  !is.na(body_genes), body_genes != "") %>%
    dplyr::mutate(gene_list = strsplit(body_genes, ",\\s*")) %>%
    tidyr::unnest(gene_list) %>%
    dplyr::rename(gene_symbol = gene_list) %>%
    dplyr::filter(gene_symbol != "") %>%
    dplyr::left_join(
      rna_df %>% dplyr::select(gene_symbol, log2FoldChange, padj),
      by = "gene_symbol"
    ) %>%
    dplyr::mutate(stripe_id = paste0("stripe_", row_number())) %>%
    dplyr::group_by(stripe_id, direction) %>%
    dplyr::summarise(
      n_with_rna = sum(!is.na(log2FoldChange)),
      n_degs = sum(!is.na(padj) & padj < DEG_PADJ & abs(log2FoldChange) > DEG_LFC, na.rm = TRUE),
      pct_degs = ifelse(n_with_rna > 0, n_degs / n_with_rna, NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(hub_status = "Non-SE-hub")

  se_hub_df <- hub_summary %>%
    dplyr::filter(direction %in% c("gained", "lost")) %>%
    dplyr::mutate(hub_status = "SE-hub") %>%
    dplyr::select(stripe_id, direction, pct_degs, hub_status)

  plot_df <- bind_rows(
    se_hub_df,
    non_hub_stripes %>% dplyr::select(stripe_id, direction, pct_degs, hub_status)
  ) %>%
    dplyr::filter(!is.na(pct_degs))

  if (nrow(plot_df) < 3) return(NULL)

  p <- ggplot(plot_df, aes(x = hub_status, y = pct_degs, fill = hub_status)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.5) +
    facet_wrap(~direction, labeller = labeller(
      direction = c("gained" = "Gained Stripes", "lost" = "Lost Stripes"))) +
    scale_fill_manual(values = c("SE-hub" = "#e0a730", "Non-SE-hub" = "#cccccc")) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       comparisons = list(c("SE-hub", "Non-SE-hub"))) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "DEG Enrichment in SE-Hub vs Non-SE-Hub Stripes",
      x = "", y = "Fraction of Body Genes That Are DEGs"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  save_multiformat_ggplot(p, file.path(output_dir, "deg_enrichment_violin"),
                          width = 8, height = 6)
  return(p)
}

plot_coordination_scatter <- function(hub_summary, output_dir) {
  plot_df <- hub_summary %>%
    dplyr::filter(direction %in% c("gained", "lost"), n_with_rna >= 2)

  if (nrow(plot_df) < 5) return(NULL)

  rho <- cor.test(plot_df$logFC_stripe, plot_df$mean_gene_lfc, method = "spearman")

  p <- ggplot(plot_df, aes(x = logFC_stripe, y = mean_gene_lfc, color = direction)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
    scale_color_manual(values = DIRECTION_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.3f\np = %.2e", rho$estimate, rho$p.value),
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(
      title = "Stripe-Gene Coordination in SE Hubs",
      x = "Stripe logFC", y = "Mean Body Gene log2FC",
      color = "Stripe Direction"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(output_dir, "coordination_scatter"),
                          width = 7, height = 6)
  return(p)
}

plot_concordance_bar <- function(hub_summary, output_dir) {
  plot_df <- hub_summary %>%
    dplyr::filter(direction %in% c("gained", "lost"), !is.na(concordant_pct)) %>%
    dplyr::mutate(concordant_class = ifelse(concordant_pct > 0.5, "Concordant", "Discordant")) %>%
    dplyr::count(direction, concordant_class)

  if (nrow(plot_df) == 0) return(NULL)

  p <- ggplot(plot_df, aes(x = direction, y = n, fill = concordant_class)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_fill_manual(values = c("Concordant" = "#2ca02c", "Discordant" = "#d62728")) +
    scale_x_discrete(labels = c("gained" = "Gained", "lost" = "Lost")) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Within-Hub Regulation Concordance",
      subtitle = ">50% body genes matching stripe direction",
      x = "Stripe Direction", y = "Proportion of SE-Hub Stripes", fill = ""
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")

  save_multiformat_ggplot(p, file.path(output_dir, "concordance_bar"),
                          width = 5, height = 6)
  return(p)
}

plot_multi_hub_dotplot <- function(multi_hub, output_dir, n_top = 20) {
  if (nrow(multi_hub) == 0) return(NULL)

  plot_df <- multi_hub %>%
    dplyr::filter(!is.na(log2FoldChange)) %>%
    dplyr::slice_max(n_hubs, n = n_top, with_ties = FALSE) %>%
    dplyr::mutate(gene_symbol = fct_reorder(gene_symbol, n_hubs))

  p <- ggplot(plot_df, aes(x = n_hubs, y = gene_symbol)) +
    geom_point(aes(color = log2FoldChange, size = n_hubs)) +
    scale_color_gradient2(low = "#4575b4", mid = "grey80", high = "#d73027", midpoint = 0) +
    scale_size_continuous(range = c(3, 8)) +
    labs(
      title = "Genes in Multiple SE-Hub Stripes",
      x = "Number of SE-Hub Stripes",
      y = "",
      color = "Gene log2FC",
      size = "Hub Count"
    ) +
    theme_classic(base_size = 14) +
    theme(axis.text.y = element_text(size = 10))

  save_multiformat_ggplot(p, file.path(output_dir, "multi_hub_dotplot"),
                          width = 8, height = max(4, 0.35 * nrow(plot_df)))
  return(p)
}

# =============================================================================
# 5. MAIN ANALYSIS
# =============================================================================

run_analysis <- function(timepoint, stripe_set_name, stripe_file) {
  config <- TIMEPOINT_CONFIG[[timepoint]]
  cat(sprintf("\n========== 1c. SE Hub Classification: %s (%s) ==========\n\n",
              config$label, stripe_set_name))

  output_dir <- file.path(BASE_DIR,
    "loops/output/superenhancer_analysis/1c_se_hub_genes", timepoint, stripe_set_name)
  tables_dir <- file.path(output_dir, "tables")
  stats_dir <- file.path(output_dir, "statistics")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Load data ---
  cat("Step 1: Loading stripes...\n")
  stripes <- load_stripes_for_hub(stripe_file)

  cat("\nStep 2: Loading superenhancers...\n")
  p60_gr <- load_se_bed(SE_FILES$P60, "P60")
  encode_gr <- load_se_bed(SE_FILES$ENCODE, "ENCODE")

  cat("\nStep 3: Identifying SE-hub stripes...\n")
  stripes <- identify_se_hubs(stripes, p60_gr, encode_gr)

  cat("\nStep 4: Extracting hub genes...\n")
  hub_gene_assign <- extract_hub_genes(stripes)

  if (nrow(hub_gene_assign) == 0) {
    cat("  No hub genes found — skipping remaining analysis\n")
    writeLines("No SE-hub stripes with body genes found for this timepoint/stripe set.",
               file.path(stats_dir, "se_hub_statistics.txt"))
    return(invisible(NULL))
  }

  cat("\nStep 5: Loading RNA-seq...\n")
  rna <- load_rna_for_se(config$rna_file, padj_thresh = DEG_PADJ, lfc_thresh = DEG_LFC)
  rna_all <- rna$all_genes

  # --- Coordination analysis ---
  cat("\nStep 6: Computing hub coordination...\n")
  coord_result <- compute_hub_coordination(hub_gene_assign, rna_all)
  hub_summary <- coord_result$hub_summary
  hub_joined <- coord_result$joined

  # --- Permutation test ---
  cat("\nStep 7: Running permutation test...\n")
  perm_result <- run_permutation_test(hub_gene_assign, rna_all)

  # --- SE-hub vs non-SE-hub ---
  cat("\nStep 8: SE-hub vs non-SE-hub comparison...\n")
  non_hub_with_genes <- stripes %>%
    dplyr::filter(!se_hub_any, !is.na(body_genes), body_genes != "",
                  direction %in% c("gained", "lost"))
  if (nrow(non_hub_with_genes) > 0 && nrow(hub_summary) > 0) {
    non_hub_genes <- non_hub_with_genes %>%
      dplyr::mutate(stripe_id = paste0("nonhub_", row_number()),
                    gene_list = strsplit(body_genes, ",\\s*")) %>%
      tidyr::unnest(gene_list) %>%
      dplyr::rename(gene_symbol = gene_list) %>%
      dplyr::filter(gene_symbol != "") %>%
      dplyr::left_join(rna_all %>% dplyr::select(gene_symbol, log2FoldChange, padj),
                       by = "gene_symbol") %>%
      dplyr::group_by(stripe_id) %>%
      dplyr::summarise(
        n_with_rna = sum(!is.na(log2FoldChange)),
        n_degs = sum(!is.na(padj) & padj < DEG_PADJ & abs(log2FoldChange) > DEG_LFC, na.rm = TRUE),
        pct_degs = ifelse(n_with_rna > 0, n_degs / n_with_rna, NA_real_),
        .groups = "drop"
      )

    hub_pct <- hub_summary$pct_degs[!is.na(hub_summary$pct_degs)]
    nonhub_pct <- non_hub_genes$pct_degs[!is.na(non_hub_genes$pct_degs)]

    if (length(hub_pct) >= 3 && length(nonhub_pct) >= 3) {
      wt <- wilcox.test(hub_pct, nonhub_pct)
      cat(sprintf("  SE-hub vs non-SE-hub DEG enrichment: median %.3f vs %.3f, p = %.2e\n",
                  median(hub_pct), median(nonhub_pct), wt$p.value))
    }
  }

  # --- Multi-hub genes ---
  cat("\nStep 9: Finding multi-hub genes...\n")
  multi_hub <- find_multi_hub_genes(hub_gene_assign, rna_all)

  # --- Save tables ---
  cat("\nStep 10: Saving tables...\n")
  write_tsv(stripes %>% dplyr::select(stripe_id, chr1, x1, x2, direction, logFC, FDR,
                                       body_genes, body_gene_count,
                                       se_hub_P60, se_hub_ENCODE, se_hub_any),
            file.path(tables_dir, "se_hub_stripes.tsv"))
  write_tsv(hub_joined, file.path(tables_dir, "hub_gene_assignments.tsv"))
  write_tsv(hub_summary, file.path(tables_dir, "hub_summary_stats.tsv"))
  if (nrow(multi_hub) > 0) {
    write_tsv(multi_hub, file.path(tables_dir, "multi_hub_genes.tsv"))
  }
  cat(sprintf("  Wrote %d TSV tables\n", 3 + as.integer(nrow(multi_hub) > 0)))

  # --- Statistics ---
  cat("\nStep 11: Writing statistics...\n")
  stats_lines <- c(
    "SE Hub Gene Classification Analysis",
    sprintf("Timepoint: %s", config$label),
    sprintf("Stripe set: %s", stripe_set_name),
    sprintf("Date: %s\n", Sys.time()),
    sprintf("Total stripes: %d", nrow(stripes)),
    sprintf("SE-hub stripes (any): %d (%.1f%%)",
            sum(stripes$se_hub_any), 100 * sum(stripes$se_hub_any) / nrow(stripes)),
    sprintf("  P60: %d, ENCODE: %d", sum(stripes$se_hub_P60), sum(stripes$se_hub_ENCODE)),
    sprintf("Hub gene assignments: %d", nrow(hub_gene_assign)),
    sprintf("Unique hub genes: %d\n", length(unique(hub_gene_assign$gene_symbol))),
    "--- Hub Coordination ---",
    sprintf("Median hub size: %d genes", median(hub_summary$hub_size)),
    sprintf("Median concordance: %.1f%%",
            100 * median(hub_summary$concordant_pct, na.rm = TRUE)),
    sprintf("Mean body gene lfc (gained hubs): %.3f",
            mean(hub_summary$mean_gene_lfc[hub_summary$direction == "gained"], na.rm = TRUE)),
    sprintf("Mean body gene lfc (lost hubs): %.3f",
            mean(hub_summary$mean_gene_lfc[hub_summary$direction == "lost"], na.rm = TRUE))
  )

  if (!is.null(perm_result)) {
    stats_lines <- c(stats_lines,
      sprintf("\n--- Permutation Test (n=%d) ---", N_PERMUTATIONS),
      sprintf("Observed mean concordance: %.4f", perm_result$observed),
      sprintf("Null mean: %.4f (sd: %.4f)", perm_result$null_mean, perm_result$null_sd),
      sprintf("Z-score: %.2f", perm_result$z_score),
      sprintf("Permutation p-value: %.4f", perm_result$p_value)
    )
  }

  if (nrow(multi_hub) > 0) {
    stats_lines <- c(stats_lines,
      sprintf("\n--- Multi-Hub Genes ---"),
      sprintf("Genes in >= 2 hubs: %d", nrow(multi_hub)),
      sprintf("Max hubs per gene: %d", max(multi_hub$n_hubs)),
      sprintf("Direction-mixed genes: %d / %d",
              sum(multi_hub$direction_mixed), nrow(multi_hub))
    )
  }

  writeLines(stats_lines, file.path(stats_dir, "se_hub_statistics.txt"))
  cat("  Wrote statistics file\n")

  # --- Plots ---
  cat("\nStep 12: Generating plots...\n")
  plot_hub_size_histogram(stripes, output_dir)
  plot_deg_enrichment_violin(hub_summary, stripes, rna_all, output_dir)
  plot_coordination_scatter(hub_summary, output_dir)
  plot_concordance_bar(hub_summary, output_dir)
  plot_multi_hub_dotplot(multi_hub, output_dir)

  cat(sprintf("\n========== Done: %s / %s ==========\n", config$label, stripe_set_name))
  return(invisible(hub_summary))
}

# =============================================================================
# 6. DISPATCH
# =============================================================================

if (parsed$timepoint == "both") {
  timepoints <- c("late", "early")
} else {
  timepoints <- parsed$timepoint
}

for (tp in timepoints) {
  config <- TIMEPOINT_CONFIG[[tp]]
  for (ss_name in names(config$stripe_files)) {
    run_analysis(tp, ss_name, config$stripe_files[[ss_name]])
  }
}
