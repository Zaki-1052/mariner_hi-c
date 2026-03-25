# biomodal/downstream/scripts/viz_sections/section_33_multi_mark_diffbind_integration.R
# Section 33: Multi-Mark DiffBind Quantitative Integration
# Standalone script - sources shared config for all dependencies and data
#
# Integrates quantitative differential binding (DiffBind) results for 4 chromatin
# marks (ATAC, H3K27ac, H3K27me3, H2AK119ub) with methylation changes to test
# the BAP1 mechanistic cascade:
#   K119ub gain -> K27me3 gain -> accessibility loss -> K27ac loss -> methylation change
#
# Upgrades binary peak-based integration (Section 19h) to quantitative fold-change
# analysis. Each mark's DiffBind output provides per-peak log2FC, FDR, and summit
# positions, enabling continuous correlation with gene-level methylation changes.
#
# Figures:
#   33a: Per-mark volcano plots (2x2 grid)
#   33b: Cross-mark Spearman correlation heatmap (7x7: 4 marks + mC + hmC + ratio)
#   33c: Quantitative O/E dot plot (upgrade from 19h binary)
#   33d: Methylation vs mark fold-change scatters (4x2 grid)
#   33e: Multivariate logistic regression forest plot (4-mark model)
#   33f: Convergence analysis (concordance count distribution + intersection bars)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_33_multi_mark_diffbind_integration.R

# =============================================================================
# SETUP
# =============================================================================

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(pROC)
  library(pheatmap)
})

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 33: MULTI-MARK DIFFBIND QUANTITATIVE INTEGRATION\n")
cat("===========================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

RATIO_TABLE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
FEATURE_MATRIX <- file.path(TABLES_DIR, "dnmt3a_feature_matrix.tsv")

SECTION_DIR <- "33_multi_mark_diffbind"
FDR_THRESHOLD <- 0.05

# Mark metadata: display name, colors, predicted direction for hypermethylation
MARK_META <- list(
  atac   = list(display = "ATAC-seq",    color = "#E6AB02", pred_dir = "down"),
  k27ac  = list(display = "H3K27ac",     color = "#FF7F00", pred_dir = "down"),
  k27me3 = list(display = "H3K27me3",    color = "#B2182B", pred_dir = "up"),
  k119ub = list(display = "H2AK119ub",   color = "#756BB1", pred_dir = "up")
)
MARK_ORDER <- c("atac", "k27ac", "k27me3", "k119ub")

# O/E direction labels per mark
OE_LABELS <- list(
  atac   = list(up = "ATAC Up",       down = "ATAC Down"),
  k27ac  = list(up = "K27ac Gained",  down = "K27ac Lost"),
  k27me3 = list(up = "K27me3 Gained", down = "K27me3 Lost"),
  k119ub = list(up = "K119ub Gained", down = "K119ub Lost")
)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# =============================================================================
# INPUT VALIDATION
# =============================================================================

cat("--- Input Validation ---\n\n")

for (mark in MARK_ORDER) {
  fp <- DIFFBIND_FILES[[mark]]
  if (!file.exists(fp)) {
    stop("DiffBind file missing for ", mark, ": ", fp)
  }
  cat(sprintf("  [OK] %s: %s\n", MARK_META[[mark]]$display, basename(fp)))
}

stopifnot(
  "demethylation_ratio_all_genes.tsv not found" = file.exists(RATIO_TABLE)
)
cat(sprintf("  [OK] Ratio table: %s\n", basename(RATIO_TABLE)))

stopifnot(
  "mc_dmr not loaded from shared config" = exists("mc_dmr"),
  "hmc_dmr not loaded from shared config" = exists("hmc_dmr")
)
cat("  [OK] mc_dmr and hmc_dmr loaded from shared config\n")

has_feature_matrix <- file.exists(FEATURE_MATRIX)
if (has_feature_matrix) {
  cat(sprintf("  [OK] Feature matrix: %s\n", basename(FEATURE_MATRIX)))
} else {
  cat(sprintf("  [--] Feature matrix not found (extended model will be skipped)\n"))
}

cat("\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load a DiffBind results file and validate its schema
#' @param filepath Path to the DiffBind TSV
#' @param mark_name Display name for logging
#' @return data.frame with DiffBind results
load_diffbind <- function(filepath, mark_name) {
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  required_cols <- c("Summit_Chr", "Summit_Start", "Summit_End", "Fold", "FDR", "p.value")
  stopifnot(all(required_cols %in% colnames(df)))

  n_up   <- sum(df$FDR < FDR_THRESHOLD & df$Fold > 0, na.rm = TRUE)
  n_down <- sum(df$FDR < FDR_THRESHOLD & df$Fold < 0, na.rm = TRUE)
  cat(sprintf("  %s: %d peaks (%d sig up, %d sig down at FDR<%.2f)\n",
              mark_name, nrow(df), n_up, n_down, FDR_THRESHOLD))
  df
}

#' Annotate DiffBind peaks to genes via ChIPseeker and aggregate per gene
#' @param diffbind_df DiffBind data.frame
#' @param mark_prefix Column name prefix for the mark
#' @param txdb_obj TxDb object
#' @return data.frame with gene-level fold, fdr, n_peaks
diffbind_to_gene <- function(diffbind_df, mark_prefix, txdb_obj) {
  # Filter rows with valid summit coordinates
  valid <- !is.na(diffbind_df$Summit_Chr) &
           !is.na(diffbind_df$Summit_Start) &
           !is.na(diffbind_df$Summit_End)
  if (sum(!valid) > 0) {
    cat(sprintf("    Filtered %d peaks with NA summit coordinates\n", sum(!valid)))
  }
  diffbind_df <- diffbind_df[valid, ]

  gr <- GRanges(
    seqnames = diffbind_df$Summit_Chr,
    ranges = IRanges(start = diffbind_df$Summit_Start,
                     end = diffbind_df$Summit_End),
    Fold = diffbind_df$Fold,
    FDR = diffbind_df$FDR
  )

  anno <- suppressMessages(annotatePeak(
    gr, TxDb = txdb_obj, annoDb = "org.Mm.eg.db", level = "gene"
  ))
  anno_df <- as.data.frame(anno)

  gene_summary <- anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      !!paste0(mark_prefix, "_fold") := Fold[which.min(abs(distanceToTSS))],
      !!paste0(mark_prefix, "_fdr")  := FDR[which.min(abs(distanceToTSS))],
      !!paste0(mark_prefix, "_n_peaks") := dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(gene = SYMBOL)

  cat(sprintf("    -> %d unique genes annotated\n", nrow(gene_summary)))
  gene_summary
}

#' Build 2x2 O/E heatmap for methylation direction x chromatin mark direction
#' Adapted from section_19:108-180 for independence.
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

#' Extract O/E for enriched quadrants (O/E > 1) from a build_2x2_heatmap result
#' Adapted from section_19:183-199.
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
# STEP 1: LOAD & ANNOTATE DIFFBIND DATA
# =============================================================================

cat("--- Step 1: Loading & annotating DiffBind data ---\n\n")

gene_tables <- list()
raw_diffbind <- list()

for (mark in MARK_ORDER) {
  cat(sprintf("Loading %s...\n", MARK_META[[mark]]$display))
  raw_diffbind[[mark]] <- load_diffbind(DIFFBIND_FILES[[mark]], MARK_META[[mark]]$display)
  gene_tables[[mark]] <- diffbind_to_gene(raw_diffbind[[mark]], mark, txdb)
  cat("\n")
}

# =============================================================================
# STEP 2: BUILD MULTI-MARK GENE PROFILE
# =============================================================================

cat("--- Step 2: Building multi-mark gene profile ---\n\n")

ratio_data <- read.table(RATIO_TABLE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Ratio table: %d genes loaded\n", nrow(ratio_data)))

multi_mark <- ratio_data
for (mark in MARK_ORDER) {
  multi_mark <- multi_mark %>%
    left_join(gene_tables[[mark]], by = "gene")
}

cat(sprintf("  Multi-mark profile: %d genes total\n", nrow(multi_mark)))
for (mark in MARK_ORDER) {
  fold_col <- paste0(mark, "_fold")
  n_with <- sum(!is.na(multi_mark[[fold_col]]))
  cat(sprintf("    %s coverage: %d genes (%.1f%%)\n",
              MARK_META[[mark]]$display, n_with, 100 * n_with / nrow(multi_mark)))
}

# Complete cases (all 4 marks)
fold_cols <- paste0(MARK_ORDER, "_fold")
complete_mask <- complete.cases(multi_mark[, fold_cols])
n_complete <- sum(complete_mask)
cat(sprintf("  Complete cases (all 4 marks): %d genes (%.1f%%)\n",
            n_complete, 100 * n_complete / nrow(multi_mark)))

cat("\n")

# =============================================================================
# FIGURE 33a: PER-MARK VOLCANO PLOTS
# =============================================================================

cat("--- FIGURE 33a: Per-Mark Volcano Plots ---\n\n")

volcano_panels <- list()
for (mark in MARK_ORDER) {
  df <- raw_diffbind[[mark]] %>%
    mutate(
      neg_log10_fdr = -log10(pmax(FDR, 1e-300)),
      sig_status = case_when(
        FDR < FDR_THRESHOLD & Fold > 0 ~ "Up",
        FDR < FDR_THRESHOLD & Fold < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )

  n_up   <- sum(df$sig_status == "Up")
  n_down <- sum(df$sig_status == "Down")

  volcano_panels[[mark]] <- ggplot(df, aes(x = Fold, y = neg_log10_fdr, color = sig_status)) +
    geom_point(size = 0.5, alpha = 0.4) +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "grey40") +
    scale_color_manual(
      values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
      name = ""
    ) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("Up: %s\nDown: %s", format(n_up, big.mark = ","),
                             format(n_down, big.mark = ",")),
             size = 3.5, color = "grey30") +
    labs(
      title = MARK_META[[mark]]$display,
      x = "log2 Fold Change (mut/ctrl)",
      y = "-log10(FDR)"
    ) +
    theme_biomodal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
}

p_33a <- (volcano_panels[["atac"]] + volcano_panels[["k27ac"]]) /
         (volcano_panels[["k27me3"]] + volcano_panels[["k119ub"]]) +
  plot_annotation(
    title = "DiffBind Differential Analysis: 4 Chromatin Marks",
    subtitle = sprintf("FDR < %.2f threshold shown as dashed line", FDR_THRESHOLD),
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
    )
  )

save_multiformat_ggplot(p_33a, file.path(OUTPUT_DIR, SECTION_DIR, "33a_diffbind_volcano_plots"),
                        width = 14, height = 10)
cat("  Saved: 33a_diffbind_volcano_plots\n\n")

# =============================================================================
# FIGURE 33b: CROSS-MARK CORRELATION HEATMAP
# =============================================================================

cat("--- FIGURE 33b: Cross-Mark Correlation Heatmap ---\n\n")

cor_cols <- c(fold_cols, "mc_diff", "hmc_diff", "delta_ratio")
cor_labels <- c(
  "ATAC-seq fold", "H3K27ac fold", "H3K27me3 fold", "H2AK119ub fold",
  "mC difference", "hmC difference", "Delta ratio"
)

cor_data <- multi_mark[, cor_cols]

# Pairwise complete Spearman correlations
n_vars <- length(cor_cols)
cor_mat <- matrix(NA, n_vars, n_vars, dimnames = list(cor_labels, cor_labels))
pval_mat <- matrix(NA, n_vars, n_vars, dimnames = list(cor_labels, cor_labels))

for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    if (i == j) {
      cor_mat[i, j] <- 1.0
      pval_mat[i, j] <- 0
    } else {
      pair_complete <- complete.cases(cor_data[, c(i, j)])
      if (sum(pair_complete) >= 10) {
        test <- cor.test(cor_data[pair_complete, i], cor_data[pair_complete, j],
                         method = "spearman", exact = FALSE)
        cor_mat[i, j] <- test$estimate
        pval_mat[i, j] <- test$p.value
      }
    }
  }
}

# Significance stars for display
sig_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01, "**",
  ifelse(p < 0.05, "*", ""))))
}

display_mat <- matrix(
  sprintf("%.2f%s", cor_mat, sapply(pval_mat, sig_stars)),
  nrow = n_vars, dimnames = list(cor_labels, cor_labels)
)
diag(display_mat) <- ""

# pheatmap
hm_breaks <- seq(-1, 1, length.out = 101)
hm_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

hm_call <- bquote(pheatmap(
  .(cor_mat),
  display_numbers = .(display_mat),
  color = .(hm_colors),
  breaks = .(hm_breaks),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  fontsize = 11,
  fontsize_number = 9,
  main = "Cross-Mark Spearman Correlations (pairwise complete cases)",
  border_color = "grey80",
  na_col = "grey90",
  number_color = "black"
))

save_multiformat_pheatmap(
  hm_call,
  file.path(OUTPUT_DIR, SECTION_DIR, "33b_cross_mark_correlation_heatmap"),
  width = 10, height = 9
)

cat("  Saved: 33b_cross_mark_correlation_heatmap\n")

# Print key correlations
cat("\n  Key cross-mark correlations:\n")
key_pairs <- list(
  c("H2AK119ub fold", "H3K27me3 fold",  "K119ub <-> K27me3 (cascade step 1)"),
  c("H3K27me3 fold",  "H3K27ac fold",   "K27me3 <-> K27ac (opposition)"),
  c("ATAC-seq fold",  "H3K27ac fold",   "ATAC <-> K27ac (concordant openness)"),
  c("H2AK119ub fold", "mC difference",  "K119ub <-> mC (mechanistic link)"),
  c("H2AK119ub fold", "hmC difference", "K119ub <-> hmC (TET impediment)")
)
for (kp in key_pairs) {
  rho <- cor_mat[kp[1], kp[2]]
  pv  <- pval_mat[kp[1], kp[2]]
  if (!is.na(rho)) {
    cat(sprintf("    %s: rho=%.3f, p=%.2e\n", kp[3], rho, pv))
  } else {
    cat(sprintf("    %s: insufficient data\n", kp[3]))
  }
}
cat("\n")

# =============================================================================
# FIGURE 33c: QUANTITATIVE O/E DOT PLOT
# =============================================================================

cat("--- FIGURE 33c: Quantitative O/E Dot Plot ---\n\n")

oe_results <- list()

for (mark in MARK_ORDER) {
  fold_col <- paste0(mark, "_fold")
  fdr_col  <- paste0(mark, "_fdr")
  labels   <- OE_LABELS[[mark]]

  # mC perspective
  mc_oe_df <- mc_dmr %>%
    dplyr::filter(significant) %>%
    left_join(gene_tables[[mark]] %>% dplyr::select(gene, !!sym(fold_col), !!sym(fdr_col)),
              by = "gene") %>%
    dplyr::filter(!is.na(!!sym(fold_col)) & !!sym(fdr_col) < FDR_THRESHOLD) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      mark_direction = ifelse(!!sym(fold_col) > 0, labels$up, labels$down)
    )

  if (nrow(mc_oe_df) >= 4) {
    mc_result <- build_2x2_heatmap(
      mc_oe_df$mc_direction, mc_oe_df$mark_direction,
      c("mC Up", "mC Down"), c(labels$up, labels$down),
      predicted_pairs = list(),
      title = "", x_label = "", y_label = ""
    )
    oe_results <- append(oe_results, list(extract_enriched_oe(mc_result, "5mC", MARK_META[[mark]]$display)))
    cat(sprintf("  mC x %s: Fisher OR=%.2f, p=%.2e (N=%d)\n",
                MARK_META[[mark]]$display, mc_result$fisher$estimate,
                mc_result$fisher$p.value, nrow(mc_oe_df)))
  } else {
    cat(sprintf("  mC x %s: insufficient overlap (N=%d)\n", MARK_META[[mark]]$display, nrow(mc_oe_df)))
  }

  # hmC perspective
  hmc_oe_df <- hmc_dmr %>%
    dplyr::filter(significant) %>%
    left_join(gene_tables[[mark]] %>% dplyr::select(gene, !!sym(fold_col), !!sym(fdr_col)),
              by = "gene") %>%
    dplyr::filter(!is.na(!!sym(fold_col)) & !!sym(fdr_col) < FDR_THRESHOLD) %>%
    mutate(
      hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
      mark_direction = ifelse(!!sym(fold_col) > 0, labels$up, labels$down)
    )

  if (nrow(hmc_oe_df) >= 4) {
    hmc_result <- build_2x2_heatmap(
      hmc_oe_df$hmc_direction, hmc_oe_df$mark_direction,
      c("hmC Down", "hmC Up"), c(labels$up, labels$down),
      predicted_pairs = list(),
      title = "", x_label = "", y_label = ""
    )
    oe_results <- append(oe_results, list(extract_enriched_oe(hmc_result, "5hmC", MARK_META[[mark]]$display)))
    cat(sprintf("  hmC x %s: Fisher OR=%.2f, p=%.2e (N=%d)\n",
                MARK_META[[mark]]$display, hmc_result$fisher$estimate,
                hmc_result$fisher$p.value, nrow(hmc_oe_df)))
  } else {
    cat(sprintf("  hmC x %s: insufficient overlap (N=%d)\n", MARK_META[[mark]]$display, nrow(hmc_oe_df)))
  }
}

comparison_df <- do.call(rbind, oe_results)

if (!is.null(comparison_df) && nrow(comparison_df) > 0) {
  comparison_df$p_label <- ifelse(comparison_df$fisher_p < 0.001,
                                  sprintf("p=%.1e", comparison_df$fisher_p),
                                  sprintf("p=%.3f", comparison_df$fisher_p))
  comparison_df$short_quadrant <- gsub(
    "hmC Down \\+ |hmC Up \\+ |mC Up \\+ |mC Down \\+ ", "", comparison_df$quadrant)
  comparison_df$mark <- factor(comparison_df$mark,
                               levels = sapply(MARK_META, `[[`, "display"))
  comparison_df$perspective <- factor(comparison_df$perspective,
                                     levels = c("5mC", "5hmC"))

  p_33c <- ggplot(comparison_df,
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
    scale_y_continuous(limits = c(min(comparison_df$oe_ratio, 1) * 0.85,
                                  max(comparison_df$oe_ratio, 1) * 1.2)) +
    facet_wrap(~ mark, scales = "free_x", nrow = 1) +
    labs(
      title = "Quantitative O/E Enrichment: DiffBind FDR-Filtered Direction",
      subtitle = "Each point = enriched quadrant from 2x2 integration | Dashed line = null (O/E=1)",
      x = "", y = "Observed / Expected Enrichment"
    ) +
    theme_biomodal() +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1.5, "lines")
    )

  save_multiformat_ggplot(p_33c, file.path(OUTPUT_DIR, SECTION_DIR, "33c_quantitative_oe_dotplot"),
                          width = 14, height = 8)
  cat("\n  Saved: 33c_quantitative_oe_dotplot\n\n")
} else {
  cat("\n  WARNING: No O/E results to plot (insufficient overlaps)\n\n")
}

# =============================================================================
# FIGURE 33d: METHYLATION VS MARK FOLD SCATTERS
# =============================================================================

cat("--- FIGURE 33d: Methylation vs Mark Fold Scatters ---\n\n")

scatter_panels <- list()
panel_idx <- 0

for (met_type in c("mc_diff", "hmc_diff")) {
  met_label <- ifelse(met_type == "mc_diff", "mC Difference", "hmC Difference")

  for (mark in MARK_ORDER) {
    panel_idx <- panel_idx + 1
    fold_col <- paste0(mark, "_fold")
    fdr_col  <- paste0(mark, "_fdr")
    sig_col  <- ifelse(met_type == "mc_diff", "mc_sig", "hmc_sig")

    plot_df <- multi_mark %>%
      dplyr::filter(!is.na(!!sym(fold_col)) & !is.na(!!sym(met_type))) %>%
      mutate(
        sig_group = case_when(
          !!sym(sig_col) & !!sym(fdr_col) < FDR_THRESHOLD ~ "Both sig",
          !!sym(sig_col) | !!sym(fdr_col) < FDR_THRESHOLD ~ "One sig",
          TRUE ~ "Neither"
        )
      )

    # Spearman correlation
    rho_test <- cor.test(plot_df[[fold_col]], plot_df[[met_type]],
                         method = "spearman", exact = FALSE)

    scatter_panels[[panel_idx]] <- ggplot(plot_df, aes(x = !!sym(fold_col), y = !!sym(met_type))) +
      geom_point(aes(color = sig_group), size = 0.3, alpha = 0.3) +
      geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      scale_color_manual(
        values = c("Both sig" = "#E41A1C", "One sig" = "#FFA500", "Neither" = "grey70"),
        name = ""
      ) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
               label = sprintf("rho=%.3f\np=%.1e\nn=%s",
                               rho_test$estimate, rho_test$p.value,
                               format(nrow(plot_df), big.mark = ",")),
               size = 3, color = "grey20") +
      labs(
        title = MARK_META[[mark]]$display,
        x = "DiffBind log2FC",
        y = met_label
      ) +
      theme_biomodal() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  }
}

# 4x2 layout: top row = mC, bottom row = hmC
p_33d <- (scatter_panels[[1]] + scatter_panels[[2]] + scatter_panels[[3]] + scatter_panels[[4]]) /
         (scatter_panels[[5]] + scatter_panels[[6]] + scatter_panels[[7]] + scatter_panels[[8]]) +
  plot_annotation(
    title = "Methylation Change vs Chromatin Mark Fold-Change",
    subtitle = "Top: mC difference | Bottom: hmC difference | Loess smooth in black",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40")
    )
  )

save_multiformat_ggplot(p_33d, file.path(OUTPUT_DIR, SECTION_DIR, "33d_methylation_vs_mark_scatters"),
                        width = 16, height = 10)
cat("  Saved: 33d_methylation_vs_mark_scatters\n\n")

# =============================================================================
# FIGURE 33e: MULTIVARIATE LOGISTIC REGRESSION FOREST PLOT
# =============================================================================

cat("--- FIGURE 33e: Multivariate Logistic Regression ---\n\n")

# Build model data: need hyper_dmr outcome + 4 mark folds
if (has_feature_matrix) {
  feat_mat <- read.table(FEATURE_MATRIX, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  model_base <- feat_mat %>% dplyr::select(gene, hyper_dmr, baseline_mc, baseline_hmc)
} else {
  model_base <- multi_mark %>%
    mutate(hyper_dmr = as.integer(mc_sig & mc_diff > 0)) %>%
    dplyr::select(gene, hyper_dmr)
}

model_data <- model_base %>%
  left_join(gene_tables[["atac"]]   %>% dplyr::select(gene, atac_fold),   by = "gene") %>%
  left_join(gene_tables[["k27ac"]]  %>% dplyr::select(gene, k27ac_fold),  by = "gene") %>%
  left_join(gene_tables[["k27me3"]] %>% dplyr::select(gene, k27me3_fold), by = "gene") %>%
  left_join(gene_tables[["k119ub"]] %>% dplyr::select(gene, k119ub_fold), by = "gene")

# Complete cases for 4-mark model
model_data_cc <- model_data %>%
  dplyr::filter(complete.cases(hyper_dmr, atac_fold, k27ac_fold, k27me3_fold, k119ub_fold))
cat(sprintf("  Model data: %d genes (complete cases for 4-mark model)\n", nrow(model_data_cc)))
cat(sprintf("  Outcome: %d hyper-DMR (%.1f%%), %d non-hyper\n",
            sum(model_data_cc$hyper_dmr), 100 * mean(model_data_cc$hyper_dmr),
            sum(!model_data_cc$hyper_dmr)))

# McFadden pseudo-R2
null_model <- glm(hyper_dmr ~ 1, data = model_data_cc, family = binomial)
null_ll <- as.numeric(logLik(null_model))
mcfadden <- function(model) {
  1 - as.numeric(logLik(model)) / null_ll
}

# Extract OR with 95% CI (from section_24 pattern)
extract_or <- function(model) {
  cc <- confint.default(model)
  coefs <- coef(model)
  data.frame(
    term = names(coefs),
    estimate = coefs,
    or = exp(coefs),
    or_lower = exp(cc[, 1]),
    or_upper = exp(cc[, 2]),
    p_value = summary(model)$coefficients[, 4],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# --- Model 1: 4-mark DiffBind ---
cat("  Fitting Model 1: 4-mark DiffBind...\n")
model_4mark <- glm(
  hyper_dmr ~ atac_fold + k27ac_fold + k27me3_fold + k119ub_fold,
  data = model_data_cc, family = binomial
)

roc_4mark <- roc(model_data_cc$hyper_dmr, predict(model_4mark, type = "response"), quiet = TRUE)
ci_4mark  <- ci.auc(roc_4mark, quiet = TRUE)
cat(sprintf("    AUC = %.3f [%.3f, %.3f], McFadden R2 = %.4f\n",
            auc(roc_4mark), ci_4mark[1], ci_4mark[3], mcfadden(model_4mark)))

# --- Model 2: Extended (4 marks + baseline) if feature matrix available ---
if (has_feature_matrix && all(c("baseline_mc", "baseline_hmc") %in% names(model_data_cc))) {
  model_data_ext <- model_data_cc %>%
    dplyr::filter(!is.na(baseline_mc) & !is.na(baseline_hmc))

  if (nrow(model_data_ext) >= 50) {
    cat("  Fitting Model 2: Extended (4 marks + baseline mC/hmC)...\n")
    null_ext <- glm(hyper_dmr ~ 1, data = model_data_ext, family = binomial)
    null_ll_ext <- as.numeric(logLik(null_ext))

    model_extended <- glm(
      hyper_dmr ~ atac_fold + k27ac_fold + k27me3_fold + k119ub_fold +
        baseline_mc + baseline_hmc,
      data = model_data_ext, family = binomial
    )
    roc_extended <- roc(model_data_ext$hyper_dmr, predict(model_extended, type = "response"), quiet = TRUE)
    ci_extended  <- ci.auc(roc_extended, quiet = TRUE)
    cat(sprintf("    AUC = %.3f [%.3f, %.3f], McFadden R2 = %.4f (N=%d)\n",
                auc(roc_extended), ci_extended[1], ci_extended[3],
                1 - as.numeric(logLik(model_extended)) / null_ll_ext,
                nrow(model_data_ext)))
  }
}

# --- Forest plot for 4-mark model ---
or_df <- extract_or(model_4mark) %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(
    display_name = case_when(
      term == "atac_fold"   ~ "ATAC-seq",
      term == "k27ac_fold"  ~ "H3K27ac",
      term == "k27me3_fold" ~ "H3K27me3",
      term == "k119ub_fold" ~ "H2AK119ub",
      TRUE ~ term
    ),
    sig_label = ifelse(p_value < 0.001, "***",
                ifelse(p_value < 0.01, "**",
                ifelse(p_value < 0.05, "*", "ns"))),
    display_name = factor(display_name,
                          levels = c("ATAC-seq", "H3K27ac", "H3K27me3", "H2AK119ub"))
  )

p_33e <- ggplot(or_df, aes(x = or, y = display_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = or_lower, xmax = or_upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5, color = "#E41A1C") +
  geom_text(aes(label = sprintf("OR=%.2f %s", or, sig_label)),
            hjust = -0.2, size = 3.5) +
  scale_x_log10() +
  labs(
    title = "4-Mark DiffBind Logistic Regression: Predicting mC Hypermethylation",
    subtitle = sprintf("N=%s genes | AUC=%.3f [%.3f, %.3f] | OR per unit log2FC increase",
                        format(nrow(model_data_cc), big.mark = ","),
                        auc(roc_4mark), ci_4mark[1], ci_4mark[3]),
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme_biomodal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

save_multiformat_ggplot(p_33e, file.path(OUTPUT_DIR, SECTION_DIR, "33e_logistic_regression_forest"),
                        width = 10, height = 6)
cat("\n  Saved: 33e_logistic_regression_forest\n\n")

# =============================================================================
# FIGURE 33f: CONVERGENCE ANALYSIS
# =============================================================================

cat("--- FIGURE 33f: Convergence Analysis ---\n\n")

# Score concordance per gene (predicted BAP1 cascade direction)
convergence_df <- multi_mark %>%
  mutate(
    atac_concordant   = !is.na(atac_fold)   & atac_fdr < FDR_THRESHOLD   & atac_fold < 0,
    k27ac_concordant  = !is.na(k27ac_fold)  & k27ac_fdr < FDR_THRESHOLD  & k27ac_fold < 0,
    k27me3_concordant = !is.na(k27me3_fold) & k27me3_fdr < FDR_THRESHOLD & k27me3_fold > 0,
    k119ub_concordant = !is.na(k119ub_fold) & k119ub_fdr < FDR_THRESHOLD & k119ub_fold > 0,
    n_concordant = as.integer(atac_concordant) + as.integer(k27ac_concordant) +
                   as.integer(k27me3_concordant) + as.integer(k119ub_concordant),
    dmr_group = case_when(
      mc_sig & mc_diff > 0  ~ "mC Hyper",
      mc_sig & mc_diff <= 0 ~ "mC Hypo",
      TRUE                  ~ "Non-DMR"
    )
  )

cat("  Concordance distribution (all genes):\n")
conc_table <- table(convergence_df$n_concordant)
for (i in names(conc_table)) {
  cat(sprintf("    %s concordant marks: %d genes\n", i, conc_table[i]))
}

# 33f-i: Stacked bar — concordance by DMR group
conc_by_group <- convergence_df %>%
  dplyr::filter(dmr_group != "Non-DMR" | n_concordant > 0) %>%
  group_by(dmr_group, n_concordant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dmr_group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    dmr_group = factor(dmr_group, levels = c("mC Hyper", "mC Hypo", "Non-DMR")),
    n_concordant = factor(n_concordant)
  )

p_33f_i <- ggplot(conc_by_group, aes(x = dmr_group, y = pct, fill = n_concordant)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(pct >= 3, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_brewer(palette = "YlOrRd", name = "Concordant\nMarks", direction = 1) +
  labs(
    title = "BAP1 Cascade Concordance by DMR Status",
    subtitle = "Concordant = mark changes in predicted direction (FDR < 0.05)",
    x = "", y = "Percentage of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

# 33f-ii: Intersection bar chart for genes with 2+ concordant marks
multi_conc <- convergence_df %>%
  dplyr::filter(n_concordant >= 2) %>%
  mutate(
    combo = paste0(
      ifelse(atac_concordant,   "ATAC", ""),
      ifelse(k27ac_concordant,  "+K27ac", ""),
      ifelse(k27me3_concordant, "+K27me3", ""),
      ifelse(k119ub_concordant, "+K119ub", "")
    ),
    combo = gsub("^\\+", "", combo)
  )

if (nrow(multi_conc) > 0) {
  combo_counts <- multi_conc %>%
    count(combo, sort = TRUE) %>%
    slice_head(n = 15)

  p_33f_ii <- ggplot(combo_counts, aes(x = reorder(combo, n), y = n)) +
    geom_col(fill = "#E41A1C", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(
      title = "Mark Combination Frequency (2+ Concordant Marks)",
      subtitle = sprintf("N = %s genes with 2+ marks changing in predicted direction",
                          format(nrow(multi_conc), big.mark = ",")),
      x = "Mark Combination",
      y = "Number of Genes"
    ) +
    theme_biomodal()

  p_33f <- p_33f_i + p_33f_ii +
    plot_layout(widths = c(1, 1.3)) +
    plot_annotation(
      title = "Epigenomic Convergence at BAP1-Responsive Genes",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )
    )
} else {
  p_33f <- p_33f_i +
    plot_annotation(
      title = "Epigenomic Convergence at BAP1-Responsive Genes",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
      )
    )
}

save_multiformat_ggplot(p_33f, file.path(OUTPUT_DIR, SECTION_DIR, "33f_convergence_analysis"),
                        width = 14, height = 8)
cat("\n  Saved: 33f_convergence_analysis\n\n")

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("--- Exporting Tables ---\n\n")

# Table 1: Full multi-mark gene profile
write.table(multi_mark,
            file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: diffbind_gene_level_all_marks.tsv (%d genes)\n", nrow(multi_mark)))

# Table 2: Cross-mark correlation matrix
cor_export <- as.data.frame(cor_mat)
cor_export$variable <- rownames(cor_export)
cor_export <- cor_export[, c("variable", cor_labels)]
write.table(cor_export,
            file.path(TABLES_DIR, "diffbind_cross_mark_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: diffbind_cross_mark_correlations.tsv\n")

# Table 3: Quantitative O/E comparison
if (!is.null(comparison_df) && nrow(comparison_df) > 0) {
  write.table(comparison_df,
              file.path(TABLES_DIR, "diffbind_quantitative_oe_comparison.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: diffbind_quantitative_oe_comparison.tsv\n")
}

# Table 4: Logistic model coefficients
write.table(or_df,
            file.path(TABLES_DIR, "diffbind_logistic_model_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: diffbind_logistic_model_coefficients.tsv\n")

# Table 5: Per-gene convergence scores
convergence_export <- convergence_df %>%
  dplyr::select(gene, chromatin_state, mc_diff, hmc_diff, delta_ratio, mc_sig, hmc_sig,
                atac_fold, atac_fdr, atac_concordant,
                k27ac_fold, k27ac_fdr, k27ac_concordant,
                k27me3_fold, k27me3_fdr, k27me3_concordant,
                k119ub_fold, k119ub_fdr, k119ub_concordant,
                n_concordant, dmr_group)
write.table(convergence_export,
            file.path(TABLES_DIR, "diffbind_convergence_per_gene.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: diffbind_convergence_per_gene.tsv (%d genes)\n", nrow(convergence_export)))

# Table 6: Convergent genes (3+ marks)
convergent_3plus <- convergence_df %>%
  dplyr::filter(n_concordant >= 3) %>%
  dplyr::select(gene, chromatin_state, mc_diff, hmc_diff, delta_ratio,
                atac_fold, k27ac_fold, k27me3_fold, k119ub_fold, n_concordant, dmr_group) %>%
  arrange(desc(n_concordant), gene)
write.table(convergent_3plus,
            file.path(TABLES_DIR, "diffbind_convergent_genes_3plus.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: diffbind_convergent_genes_3plus.tsv (%d genes)\n", nrow(convergent_3plus)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 33 SUMMARY\n")
cat("===========================================================================\n\n")

cat("--- Per-Mark Statistics ---\n")
for (mark in MARK_ORDER) {
  df <- raw_diffbind[[mark]]
  n_up   <- sum(df$FDR < FDR_THRESHOLD & df$Fold > 0, na.rm = TRUE)
  n_down <- sum(df$FDR < FDR_THRESHOLD & df$Fold < 0, na.rm = TRUE)
  n_genes <- nrow(gene_tables[[mark]])
  cat(sprintf("  %-12s %6d peaks | %5d up, %5d down (FDR<%.2f) | %5d genes\n",
              paste0(MARK_META[[mark]]$display, ":"),
              nrow(df), n_up, n_down, FDR_THRESHOLD, n_genes))
}

cat(sprintf("\n--- Multi-Mark Profile: %d genes, %d complete cases ---\n",
            nrow(multi_mark), n_complete))

cat(sprintf("\n--- Logistic Model: AUC=%.3f [%.3f, %.3f], N=%d ---\n",
            auc(roc_4mark), ci_4mark[1], ci_4mark[3], nrow(model_data_cc)))
for (i in 1:nrow(or_df)) {
  cat(sprintf("  %-12s OR=%.3f [%.3f, %.3f] p=%.2e\n",
              or_df$display_name[i], or_df$or[i], or_df$or_lower[i],
              or_df$or_upper[i], or_df$p_value[i]))
}

cat(sprintf("\n--- Convergence: %d genes with 3+ concordant marks ---\n",
            nrow(convergent_3plus)))
if (nrow(convergent_3plus) > 0) {
  cat(sprintf("  Hyper-DMR: %d | Hypo-DMR: %d | Non-DMR: %d\n",
              sum(convergent_3plus$dmr_group == "mC Hyper"),
              sum(convergent_3plus$dmr_group == "mC Hypo"),
              sum(convergent_3plus$dmr_group == "Non-DMR")))
}

cat("\n--- Output Locations ---\n")
cat(sprintf("  Figures: %s/33_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/diffbind_*.tsv\n", TABLES_DIR))

cat("\n=== Section 33 Complete ===\n\n")
