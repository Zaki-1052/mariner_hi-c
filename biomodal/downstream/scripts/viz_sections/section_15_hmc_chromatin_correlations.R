# biomodal/downstream/scripts/viz_sections/section_15_hmc_chromatin_correlations.R
# Section 15: 5hmC DMR Correlations with Chromatin Marks
# Standalone script - sources shared config for all dependencies and data
#
# Tests the BAP1 mechanistic model using 5hmC DMR direction (complement to
# sections 11/12/14 which use 5mC direction). Since mC and hmC go in opposite
# directions in BAP1-KO (92.3% coordinated mC up/hmC down), the predicted
# enrichment quadrants flip:
#
#   MeCP2:  hmC Down + MeCP2 Up   (less hmC -> more MeCP2 binding)
#   ATAC:   hmC Down + ATAC Down  (hmC tracks with accessibility)
#   K119ub: hmC Down + K119ub Gained (K119ub blocks TET -> hmC depletes)
#
# Key question: Do hmC O/E enrichments match or exceed those from mC?
# The hmC Down group (4,361 genes) is much larger than hmC Up (536 genes).
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_15_hmc_chromatin_correlations.R

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# Load ChIPseeker for gene annotation (not in shared config)
suppressPackageStartupMessages(library(ChIPseeker))

# =============================================================================
# HELPER FUNCTIONS (adapted from sections 11, 12, 14 for independence)
# =============================================================================

#' Load MeCP2 annotated data from HOMER-style file
#' @param filepath Path to MeCP2 annotated text file
#' @return data.frame with MeCP2 binding data
load_mecp2_annotated <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("MeCP2 annotated file not found: ", filepath)
    return(NULL)
  }

  df <- read.table(filepath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE, quote = "",
                   comment.char = "")

  # Ensure numeric columns
  numeric_cols <- c("Fold", "FDR", "p.value", "distanceToTSS",
                    "geneStart", "geneEnd", "geneLength")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }

  return(df)
}

#' Aggregate MeCP2 peaks per gene (nearest peak to TSS)
#' @param mecp2_df MeCP2 annotated data.frame
#' @return data.frame with gene-level MeCP2 summary
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

#' Annotate peaks to nearest genes using ChIPseeker
#' @param gr GRanges of peaks
#' @param txdb TxDb object
#' @return data.frame with peak-to-gene mapping
annotate_peaks_to_genes <- function(gr, txdb) {
  anno <- suppressMessages(annotatePeak(
    gr, TxDb = txdb, annoDb = "org.Mm.eg.db", level = "gene"
  ))
  as.data.frame(anno)
}

#' Aggregate peaks per gene with a directional label
#' @param anno_df Annotated peak data.frame (from ChIPseeker)
#' @param label Character label for direction
#' @return data.frame with gene-level peak counts
aggregate_peaks_by_gene <- function(anno_df, label) {
  anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      !!paste0("n_", label) := n(),
      .groups = "drop"
    )
}

#' Derive differential (gained/lost) peaks from condition-specific sets
#' @param ctrl_gr GRanges of control peaks
#' @param mut_gr GRanges of mutant peaks
#' @return list with gained, lost, and shared GRanges
derive_differential_peaks <- function(ctrl_gr, mut_gr) {
  ctrl_in_mut <- countOverlaps(ctrl_gr, mut_gr) > 0
  mut_in_ctrl <- countOverlaps(mut_gr, ctrl_gr) > 0

  list(
    gained = mut_gr[!mut_in_ctrl],
    lost = ctrl_gr[!ctrl_in_mut],
    shared = ctrl_gr[ctrl_in_mut]
  )
}

#' Build 2x2 O/E heatmap for methylation direction x chromatin mark direction
#' @param met_direction Character vector of methylation direction per gene
#' @param mark_direction Character vector of chromatin mark direction per gene
#' @param met_levels Character vector of 2 levels for methylation axis
#' @param mark_levels Character vector of 2 levels for chromatin mark axis
#' @param predicted_pairs list of 2 vectors defining predicted quadrant pairs
#' @param title Character, plot title
#' @param subtitle_extra Character, appended to subtitle
#' @param x_label Character, x-axis label
#' @param y_label Character, y-axis label
#' @return list with ggplot object, heatmap data, fisher result
build_2x2_heatmap <- function(met_direction, mark_direction,
                              met_levels, mark_levels,
                              predicted_pairs,
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

  # Mark predicted quadrants
  for (pair in predicted_pairs) {
    heatmap_data$is_predicted[heatmap_data$met_dir == pair[1] &
                              heatmap_data$mark_dir == pair[2]] <- TRUE
  }

  p <- ggplot(heatmap_data, aes(x = mark_dir, y = met_dir, fill = Enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 4, lineheight = 1.2) +
    geom_tile(data = heatmap_data %>% dplyr::filter(is_predicted),
              aes(x = mark_dir, y = met_dir),
              fill = NA, color = "black", linewidth = 1.5, linetype = "solid") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                         limits = c(min(heatmap_data$Enrichment) * 0.9,
                                    max(heatmap_data$Enrichment) * 1.1)) +
    labs(
      title = title,
      subtitle = sprintf("Fisher's exact test: OR = %.2f, p = %.2e | Black borders = predicted quadrants%s",
                         fisher_result$estimate, fisher_result$p.value, subtitle_extra),
      x = x_label, y = y_label
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  list(plot = p, data = heatmap_data, fisher = fisher_result, table = quad_table)
}

# =============================================================================
# SECTION 15: 5hmC DMR CORRELATIONS WITH CHROMATIN MARKS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 15: 5hmC DMR CORRELATIONS WITH CHROMATIN MARKS\n")
cat("================================================================================\n\n")

# -----------------------------------------------------------------------
# Prepare hmC DMR gene sets
# -----------------------------------------------------------------------
cat("Preparing hmC DMR gene sets...\n")

hmc_sig <- hmc_dmr %>% dplyr::filter(significant)
hmc_down_genes <- hmc_sig %>% dplyr::filter(mod_difference < 0)
hmc_up_genes <- hmc_sig %>% dplyr::filter(mod_difference > 0)

cat(sprintf("  Total significant hmC DMRs: %d\n", nrow(hmc_sig)))
cat(sprintf("  hmC Down (decreased 5hmC): %d (%.1f%%)\n",
            nrow(hmc_down_genes), 100 * nrow(hmc_down_genes) / nrow(hmc_sig)))
cat(sprintf("  hmC Up (increased 5hmC):   %d (%.1f%%)\n",
            nrow(hmc_up_genes), 100 * nrow(hmc_up_genes) / nrow(hmc_sig)))

# Also prepare mC DMR gene sets (needed for comparison in 15d)
mc_sig <- mc_dmr %>% dplyr::filter(significant)
mc_up_genes <- mc_sig %>% dplyr::filter(mod_difference > 0)
mc_down_genes <- mc_sig %>% dplyr::filter(mod_difference < 0)

cat(sprintf("\n  Total significant mC DMRs: %d (for comparison)\n", nrow(mc_sig)))
cat(sprintf("  mC Up: %d, mC Down: %d\n", nrow(mc_up_genes), nrow(mc_down_genes)))

# -----------------------------------------------------------------------
# Load TxDb for ChIPseeker annotations
# -----------------------------------------------------------------------
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# =============================================================================
# FIGURE 15a: hmC Direction x MeCP2 Direction (2x2 Heatmap)
# =============================================================================
cat("\n--- FIGURE 15a: hmC Direction x MeCP2 Direction ---\n\n")

cat("Loading MeCP2 differential binding data...\n")
mecp2_annotated <- load_mecp2_annotated(MECP2_FILES$annotated)

if (is.null(mecp2_annotated)) {
  cat("  ERROR: MeCP2 data not found. Skipping 15a.\n")
  result_15a <- NULL
} else {

  mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated)
  cat(sprintf("  %d unique genes with MeCP2 peaks\n", nrow(mecp2_gene)))

  # Join hmC DMRs to MeCP2 gene-level data
  hmc_mecp2 <- hmc_sig %>%
    left_join(mecp2_gene %>% dplyr::select(SYMBOL, nearest_fold),
              by = c("gene" = "SYMBOL")) %>%
    dplyr::filter(!is.na(nearest_fold)) %>%
    mutate(
      hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
      mecp2_direction = ifelse(nearest_fold > 0, "MeCP2 Up", "MeCP2 Down")
    )

  cat(sprintf("  %d significant hmC DMR genes matched to MeCP2 data\n", nrow(hmc_mecp2)))

  result_15a <- build_2x2_heatmap(
    met_direction = hmc_mecp2$hmc_direction,
    mark_direction = hmc_mecp2$mecp2_direction,
    met_levels = c("hmC Down", "hmC Up"),
    mark_levels = c("MeCP2 Up", "MeCP2 Down"),
    predicted_pairs = list(c("hmC Down", "MeCP2 Up"), c("hmC Up", "MeCP2 Down")),
    title = "Integration: 5hmC Direction x MeCP2 Binding Direction",
    x_label = "MeCP2 Binding Direction (Mutant vs Control)",
    y_label = "5hmC DMR Direction"
  )

  save_multiformat_ggplot(result_15a$plot,
                          file.path(OUTPUT_DIR, "15a_hmc_mecp2_heatmap"),
                          width = 9, height = 7)

  cat("  Contingency table:\n")
  print(result_15a$table)
  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_15a$fisher$estimate, result_15a$fisher$p.value))

  # Save table
  write.table(result_15a$data %>% dplyr::select(-is_predicted),
              file.path(TABLES_DIR, "hmc_mecp2_integration.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: hmc_mecp2_integration.tsv\n")

  # --- Also compute mC x MeCP2 for comparison (15d) ---
  mc_mecp2 <- mc_sig %>%
    left_join(mecp2_gene %>% dplyr::select(SYMBOL, nearest_fold),
              by = c("gene" = "SYMBOL")) %>%
    dplyr::filter(!is.na(nearest_fold)) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      mecp2_direction = ifelse(nearest_fold > 0, "MeCP2 Up", "MeCP2 Down")
    )

  mc_mecp2_result <- build_2x2_heatmap(
    met_direction = mc_mecp2$mc_direction,
    mark_direction = mc_mecp2$mecp2_direction,
    met_levels = c("mC Up", "mC Down"),
    mark_levels = c("MeCP2 Up", "MeCP2 Down"),
    predicted_pairs = list(c("mC Up", "MeCP2 Up"), c("mC Down", "MeCP2 Down")),
    title = "(Reference) 5mC Direction x MeCP2 Direction",
    x_label = "MeCP2 Direction", y_label = "5mC DMR Direction"
  )
}

# =============================================================================
# FIGURE 15b: hmC Direction x ATAC Direction (2x2 Heatmap)
# =============================================================================
cat("\n--- FIGURE 15b: hmC Direction x ATAC Direction ---\n\n")

cat("Loading ATAC-seq differential peaks...\n")
atac_up_gr <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down_gr <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")

if (is.null(atac_up_gr) || is.null(atac_down_gr)) {
  cat("  ERROR: ATAC data not found. Skipping 15b.\n")
  result_15b <- NULL
} else {

  # Annotate ATAC peaks to genes
  cat("  Annotating ATAC peaks to genes...\n")
  atac_up_anno <- annotate_peaks_to_genes(atac_up_gr, txdb)
  atac_down_anno <- annotate_peaks_to_genes(atac_down_gr, txdb)

  atac_up_gene <- aggregate_peaks_by_gene(atac_up_anno, "atac_up")
  atac_down_gene <- aggregate_peaks_by_gene(atac_down_anno, "atac_down")

  atac_gene <- full_join(atac_up_gene, atac_down_gene, by = "SYMBOL") %>%
    mutate(
      n_atac_up = replace_na(n_atac_up, 0),
      n_atac_down = replace_na(n_atac_down, 0),
      net_atac = n_atac_up - n_atac_down
    )

  cat(sprintf("  %d genes with ATAC data\n", nrow(atac_gene)))

  # Join hmC DMRs to ATAC gene-level data (exclude net_atac == 0)
  hmc_atac <- hmc_sig %>%
    left_join(atac_gene %>% dplyr::select(SYMBOL, net_atac),
              by = c("gene" = "SYMBOL")) %>%
    mutate(net_atac = replace_na(net_atac, 0)) %>%
    dplyr::filter(net_atac != 0) %>%
    mutate(
      hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
      atac_direction = ifelse(net_atac > 0, "ATAC Up", "ATAC Down")
    )

  cat(sprintf("  %d significant hmC DMR genes with non-zero net ATAC\n", nrow(hmc_atac)))

  result_15b <- build_2x2_heatmap(
    met_direction = hmc_atac$hmc_direction,
    mark_direction = hmc_atac$atac_direction,
    met_levels = c("hmC Down", "hmC Up"),
    mark_levels = c("ATAC Up", "ATAC Down"),
    predicted_pairs = list(c("hmC Down", "ATAC Down"), c("hmC Up", "ATAC Up")),
    title = "Integration: 5hmC Direction x ATAC-seq Direction",
    x_label = "ATAC-seq Direction (Mutant vs Control)",
    y_label = "5hmC DMR Direction"
  )

  save_multiformat_ggplot(result_15b$plot,
                          file.path(OUTPUT_DIR, "15b_hmc_atac_heatmap"),
                          width = 9, height = 7)

  cat("  Contingency table:\n")
  print(result_15b$table)
  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_15b$fisher$estimate, result_15b$fisher$p.value))

  # Save table
  write.table(result_15b$data %>% dplyr::select(-is_predicted),
              file.path(TABLES_DIR, "hmc_atac_integration.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: hmc_atac_integration.tsv\n")

  # --- Also compute mC x ATAC for comparison (15d) ---
  mc_atac <- mc_sig %>%
    left_join(atac_gene %>% dplyr::select(SYMBOL, net_atac),
              by = c("gene" = "SYMBOL")) %>%
    mutate(net_atac = replace_na(net_atac, 0)) %>%
    dplyr::filter(net_atac != 0) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      atac_direction = ifelse(net_atac > 0, "ATAC Up", "ATAC Down")
    )

  mc_atac_result <- build_2x2_heatmap(
    met_direction = mc_atac$mc_direction,
    mark_direction = mc_atac$atac_direction,
    met_levels = c("mC Up", "mC Down"),
    mark_levels = c("ATAC Up", "ATAC Down"),
    predicted_pairs = list(c("mC Up", "ATAC Down"), c("mC Down", "ATAC Up")),
    title = "(Reference) 5mC Direction x ATAC Direction",
    x_label = "ATAC Direction", y_label = "5mC DMR Direction"
  )
}

# =============================================================================
# FIGURE 15c: hmC Direction x K119ub Direction (2x2 Heatmap)
# =============================================================================
cat("\n--- FIGURE 15c: hmC Direction x K119ub Direction ---\n\n")

cat("Loading H2AK119ub condition-specific peaks...\n")
k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Control")
k119ub_mut_gr <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mutant")

if (is.null(k119ub_ctrl_gr) || is.null(k119ub_mut_gr)) {
  cat("  ERROR: K119ub data not found. Skipping 15c.\n")
  result_15c <- NULL
} else {

  # Derive gained/lost peaks
  diff_peaks <- derive_differential_peaks(k119ub_ctrl_gr, k119ub_mut_gr)
  cat(sprintf("  K119ub Gained (mut-only): %d\n", length(diff_peaks$gained)))
  cat(sprintf("  K119ub Lost (ctrl-only):  %d\n", length(diff_peaks$lost)))

  # Annotate to genes
  cat("  Annotating K119ub peaks to genes...\n")
  k119ub_ctrl_anno <- annotate_peaks_to_genes(k119ub_ctrl_gr, txdb)
  k119ub_mut_anno <- annotate_peaks_to_genes(k119ub_mut_gr, txdb)

  k119ub_ctrl_gene <- aggregate_peaks_by_gene(k119ub_ctrl_anno, "k119ub_ctrl")
  k119ub_mut_gene <- aggregate_peaks_by_gene(k119ub_mut_anno, "k119ub_mut")

  k119ub_gene <- full_join(k119ub_ctrl_gene, k119ub_mut_gene, by = "SYMBOL") %>%
    mutate(
      n_k119ub_ctrl = replace_na(n_k119ub_ctrl, 0),
      n_k119ub_mut = replace_na(n_k119ub_mut, 0),
      net_k119ub = n_k119ub_mut - n_k119ub_ctrl
    )

  cat(sprintf("  %d genes with K119ub data\n", nrow(k119ub_gene)))

  # Join hmC DMRs to K119ub gene-level data (exclude net == 0)
  hmc_k119ub <- hmc_sig %>%
    left_join(k119ub_gene %>% dplyr::select(SYMBOL, net_k119ub),
              by = c("gene" = "SYMBOL")) %>%
    mutate(net_k119ub = replace_na(net_k119ub, 0)) %>%
    dplyr::filter(net_k119ub != 0) %>%
    mutate(
      hmc_direction = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
      k119ub_direction = ifelse(net_k119ub > 0, "K119ub Gained", "K119ub Lost")
    )

  cat(sprintf("  %d significant hmC DMR genes with non-zero net K119ub\n", nrow(hmc_k119ub)))

  result_15c <- build_2x2_heatmap(
    met_direction = hmc_k119ub$hmc_direction,
    mark_direction = hmc_k119ub$k119ub_direction,
    met_levels = c("hmC Down", "hmC Up"),
    mark_levels = c("K119ub Gained", "K119ub Lost"),
    predicted_pairs = list(c("hmC Down", "K119ub Gained"), c("hmC Up", "K119ub Lost")),
    title = "Integration: 5hmC Direction x H2AK119ub Direction",
    x_label = "H2AK119ub Direction (Mutant vs Control)",
    y_label = "5hmC DMR Direction"
  )

  save_multiformat_ggplot(result_15c$plot,
                          file.path(OUTPUT_DIR, "15c_hmc_k119ub_heatmap"),
                          width = 9, height = 7)

  cat("  Contingency table:\n")
  print(result_15c$table)
  cat(sprintf("  Fisher OR = %.2f, p = %.2e\n",
              result_15c$fisher$estimate, result_15c$fisher$p.value))

  # Save table
  write.table(result_15c$data %>% dplyr::select(-is_predicted),
              file.path(TABLES_DIR, "hmc_k119ub_integration.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: hmc_k119ub_integration.tsv\n")

  # --- Also compute mC x K119ub for comparison (15d) ---
  mc_k119ub <- mc_sig %>%
    left_join(k119ub_gene %>% dplyr::select(SYMBOL, net_k119ub),
              by = c("gene" = "SYMBOL")) %>%
    mutate(net_k119ub = replace_na(net_k119ub, 0)) %>%
    dplyr::filter(net_k119ub != 0) %>%
    mutate(
      mc_direction = ifelse(mod_difference > 0, "mC Up", "mC Down"),
      k119ub_direction = ifelse(net_k119ub > 0, "K119ub Gained", "K119ub Lost")
    )

  mc_k119ub_result <- build_2x2_heatmap(
    met_direction = mc_k119ub$mc_direction,
    mark_direction = mc_k119ub$k119ub_direction,
    met_levels = c("mC Up", "mC Down"),
    mark_levels = c("K119ub Gained", "K119ub Lost"),
    predicted_pairs = list(c("mC Up", "K119ub Gained"), c("mC Down", "K119ub Lost")),
    title = "(Reference) 5mC Direction x K119ub Direction",
    x_label = "K119ub Direction", y_label = "5mC DMR Direction"
  )
}

# =============================================================================
# FIGURE 15d: Comparative Summary - mC vs hmC Enrichment
# =============================================================================
cat("\n--- FIGURE 15d: Comparative Summary (mC vs hmC Enrichment) ---\n\n")

# Collect O/E values from all 6 analyses (3 marks x 2 perspectives)
comparison_rows <- list()

# Helper to extract O/E for predicted quadrants from a build_2x2_heatmap result
extract_predicted_oe <- function(result, perspective, mark_name) {
  if (is.null(result)) return(NULL)
  predicted <- result$data %>% dplyr::filter(is_predicted)
  rows <- lapply(1:nrow(predicted), function(i) {
    data.frame(
      perspective = perspective,
      mark = mark_name,
      quadrant = paste(predicted$met_dir[i], "+", predicted$mark_dir[i]),
      observed = predicted$Observed[i],
      expected = predicted$Expected[i],
      oe_ratio = predicted$Enrichment[i],
      fisher_or = result$fisher$estimate,
      fisher_p = result$fisher$p.value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# hmC perspective
comparison_rows[[1]] <- extract_predicted_oe(result_15a, "5hmC", "MeCP2")
comparison_rows[[2]] <- extract_predicted_oe(result_15b, "5hmC", "ATAC")
comparison_rows[[3]] <- extract_predicted_oe(result_15c, "5hmC", "K119ub")

# mC perspective (re-computed above)
if (exists("mc_mecp2_result")) {
  comparison_rows[[4]] <- extract_predicted_oe(mc_mecp2_result, "5mC", "MeCP2")
}
if (exists("mc_atac_result")) {
  comparison_rows[[5]] <- extract_predicted_oe(mc_atac_result, "5mC", "ATAC")
}
if (exists("mc_k119ub_result")) {
  comparison_rows[[6]] <- extract_predicted_oe(mc_k119ub_result, "5mC", "K119ub")
}

comparison_df <- do.call(rbind, comparison_rows)

if (is.null(comparison_df) || nrow(comparison_df) == 0) {
  cat("  ERROR: No comparison data available. Skipping 15d.\n")
} else {

  cat("  O/E comparison summary:\n")
  for (i in 1:nrow(comparison_df)) {
    cat(sprintf("    %s | %s | %s: O/E = %.2f (Obs=%d, Exp=%.0f)\n",
                comparison_df$perspective[i], comparison_df$mark[i],
                comparison_df$quadrant[i],
                comparison_df$oe_ratio[i], comparison_df$observed[i],
                comparison_df$expected[i]))
  }

  # Format p-value labels
  comparison_df$p_label <- ifelse(comparison_df$fisher_p < 0.001,
                                  sprintf("p=%.1e", comparison_df$fisher_p),
                                  sprintf("p=%.3f", comparison_df$fisher_p))

  # Create short quadrant labels for readability
  comparison_df$short_quadrant <- gsub("hmC Down \\+ |hmC Up \\+ |mC Up \\+ |mC Down \\+ ",
                                       "", comparison_df$quadrant)

  # Factor ordering
  comparison_df$mark <- factor(comparison_df$mark, levels = c("MeCP2", "ATAC", "K119ub"))
  comparison_df$perspective <- factor(comparison_df$perspective, levels = c("5mC", "5hmC"))

  # Grouped dot plot
  p_15d <- ggplot(comparison_df,
                  aes(x = mark, y = oe_ratio,
                      color = perspective, shape = perspective)) +
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
      title = "mC vs hmC: O/E Enrichment in Predicted Quadrants",
      subtitle = "Each point = one predicted quadrant from a 2x2 integration heatmap | Dashed line = null (O/E=1)",
      x = "", y = "Observed / Expected Enrichment"
    ) +
    theme_biomodal() +
    theme(
      legend.position = "top",
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1.5, "lines")
    )

  save_multiformat_ggplot(p_15d,
                          file.path(OUTPUT_DIR, "15d_mc_vs_hmc_enrichment_comparison"),
                          width = 12, height = 8)

  # Save comparison table
  write.table(comparison_df,
              file.path(TABLES_DIR, "hmc_vs_mc_enrichment_comparison.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: hmc_vs_mc_enrichment_comparison.tsv\n")
}

# =============================================================================
# SECTION 15 SUMMARY
# =============================================================================
cat("\n")
cat("================================================================================\n")
cat("SECTION 15 SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("Significant hmC DMRs: %d (Down: %d, Up: %d)\n",
            nrow(hmc_sig), nrow(hmc_down_genes), nrow(hmc_up_genes)))
cat(sprintf("Significant mC DMRs: %d (Up: %d, Down: %d)\n",
            nrow(mc_sig), nrow(mc_up_genes), nrow(mc_down_genes)))

if (!is.null(result_15a)) {
  cat(sprintf("\n15a hmC x MeCP2:  Fisher OR = %.2f, p = %.2e\n",
              result_15a$fisher$estimate, result_15a$fisher$p.value))
  pred_15a <- result_15a$data %>% dplyr::filter(is_predicted)
  for (i in 1:nrow(pred_15a)) {
    cat(sprintf("     %s + %s: O/E = %.2f (n=%d)\n",
                pred_15a$met_dir[i], pred_15a$mark_dir[i],
                pred_15a$Enrichment[i], pred_15a$Observed[i]))
  }
}

if (!is.null(result_15b)) {
  cat(sprintf("\n15b hmC x ATAC:   Fisher OR = %.2f, p = %.2e\n",
              result_15b$fisher$estimate, result_15b$fisher$p.value))
  pred_15b <- result_15b$data %>% dplyr::filter(is_predicted)
  for (i in 1:nrow(pred_15b)) {
    cat(sprintf("     %s + %s: O/E = %.2f (n=%d)\n",
                pred_15b$met_dir[i], pred_15b$mark_dir[i],
                pred_15b$Enrichment[i], pred_15b$Observed[i]))
  }
}

if (!is.null(result_15c)) {
  cat(sprintf("\n15c hmC x K119ub: Fisher OR = %.2f, p = %.2e\n",
              result_15c$fisher$estimate, result_15c$fisher$p.value))
  pred_15c <- result_15c$data %>% dplyr::filter(is_predicted)
  for (i in 1:nrow(pred_15c)) {
    cat(sprintf("     %s + %s: O/E = %.2f (n=%d)\n",
                pred_15c$met_dir[i], pred_15c$mark_dir[i],
                pred_15c$Enrichment[i], pred_15c$Observed[i]))
  }
}

if (!is.null(comparison_df) && nrow(comparison_df) > 0) {
  cat("\n--- mC vs hmC Comparison ---\n")
  for (m in c("MeCP2", "ATAC", "K119ub")) {
    mc_rows <- comparison_df[comparison_df$perspective == "5mC" & comparison_df$mark == m, ]
    hmc_rows <- comparison_df[comparison_df$perspective == "5hmC" & comparison_df$mark == m, ]
    if (nrow(mc_rows) > 0 && nrow(hmc_rows) > 0) {
      mc_mean <- mean(mc_rows$oe_ratio)
      hmc_mean <- mean(hmc_rows$oe_ratio)
      winner <- ifelse(hmc_mean > mc_mean, "hmC", "mC")
      cat(sprintf("  %s: mC mean O/E = %.2f, hmC mean O/E = %.2f -> %s stronger\n",
                  m, mc_mean, hmc_mean, winner))
    }
  }
}

cat("\nSection 15 complete.\n\n")
