# biomodal/downstream/scripts/viz_sections/section_13_atac_chromatin_and_loops.R
# Section 13: ATAC Chromatin State Enrichment + Loop Anchor Accessibility
# Standalone script - sources shared config for all dependencies and data
#
# Part A (13a-c): Are ATAC-down peaks enriched at Polycomb/H3K27me3 regions?
#   BAP1 loss -> H2AK119ub -> PRC1 -> chromatin compaction
#   Prediction: ATAC-down enriched in Polycomb/Repressed_Promoter
#
# Part B (13d-f): Do Hi-C loop anchors show coordinated accessibility changes?
#   Prediction: up_in_mutant loop anchors -> more ATAC-up overlap
#               down_in_mutant loop anchors -> more ATAC-down overlap
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_13_atac_chromatin_and_loops.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages(library(ChIPseeker))

# =============================================================================
# SECTION 13: ATAC CHROMATIN STATE ENRICHMENT + LOOP ANCHOR ACCESSIBILITY
# =============================================================================

cat("================================================================================\n")
cat("SECTION 13: ATAC CHROMATIN STATE ENRICHMENT + LOOP ANCHOR ACCESSIBILITY\n")
cat("================================================================================\n\n")

# -----------------------------------------------------------------------
# Load ATAC differential peaks
# -----------------------------------------------------------------------

cat("Loading ATAC-seq differential peaks...\n")
atac_up_gr <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down_gr <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")

if (is.null(atac_up_gr) || is.null(atac_down_gr)) {
  cat("  ERROR: ATAC differential peak files not found. Skipping Section 13.\n")
} else {

  # -----------------------------------------------------------------------
  # Load ChIP-seq peaks for chromatin state classification
  # -----------------------------------------------------------------------

  cat("\nLoading ChIP-seq peaks for chromatin state classification...\n")
  chip_peaks <- list(
    ctcf     = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
    h3k27ac  = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
    h3k27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
    h3k4me1  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
    h3k4me3  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
    bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
  )

  # -----------------------------------------------------------------------
  # Classify ATAC peaks into 7 chromatin states
  # -----------------------------------------------------------------------

  cat("\nClassifying ATAC peaks into chromatin states...\n")
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)
  tss_gr <- resize(genes_gr, width = 1, fix = "start")

  # ATAC Up peaks
  overlaps_up <- compute_chip_overlaps(atac_up_gr, chip_peaks)
  dist_up <- distanceToNearest(atac_up_gr, tss_gr)
  distance_up <- rep(NA_real_, length(atac_up_gr))
  distance_up[queryHits(dist_up)] <- mcols(dist_up)$distance
  state_up <- classify_chromatin_state(overlaps_up, distance_up, TSS_THRESHOLD)

  # ATAC Down peaks
  overlaps_down <- compute_chip_overlaps(atac_down_gr, chip_peaks)
  dist_down <- distanceToNearest(atac_down_gr, tss_gr)
  distance_down <- rep(NA_real_, length(atac_down_gr))
  distance_down[queryHits(dist_down)] <- mcols(dist_down)$distance
  state_down <- classify_chromatin_state(overlaps_down, distance_down, TSS_THRESHOLD)

  cat(sprintf("  ATAC Up:   %d peaks classified\n", length(state_up)))
  cat(sprintf("  ATAC Down: %d peaks classified\n", length(state_down)))

  # =====================================================================
  # PART A: ATAC-down Polycomb Enrichment (Figures 13a-c)
  # =====================================================================

  cat("\n--- PART A: ATAC Chromatin State Enrichment ---\n\n")

  # -----------------------------------------------------------------------
  # FIGURE 13a: Chromatin State Distribution of ATAC-up vs ATAC-down Peaks
  # -----------------------------------------------------------------------
  cat("Creating Figure 13a: Chromatin state distribution...\n")

  state_counts_up <- as.data.frame(table(state_up))
  colnames(state_counts_up) <- c("state", "count")
  state_counts_up$direction <- "ATAC Up"
  state_counts_up$pct <- 100 * state_counts_up$count / sum(state_counts_up$count)

  state_counts_down <- as.data.frame(table(state_down))
  colnames(state_counts_down) <- c("state", "count")
  state_counts_down$direction <- "ATAC Down"
  state_counts_down$pct <- 100 * state_counts_down$count / sum(state_counts_down$count)

  state_dist <- rbind(state_counts_up, state_counts_down)
  state_dist$state <- factor(state_dist$state, levels = CHROMATIN_STATE_ORDER)
  state_dist$direction <- factor(state_dist$direction, levels = c("ATAC Up", "ATAC Down"))

  # Print summary
  for (s in CHROMATIN_STATE_ORDER) {
    pct_up <- state_dist$pct[state_dist$state == s & state_dist$direction == "ATAC Up"]
    pct_dn <- state_dist$pct[state_dist$state == s & state_dist$direction == "ATAC Down"]
    if (length(pct_up) > 0 && length(pct_dn) > 0) {
      cat(sprintf("  %s: Up=%.1f%%, Down=%.1f%%\n", s, pct_up, pct_dn))
    }
  }

  p_13a <- ggplot(state_dist, aes(x = state, y = pct, fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
    scale_fill_manual(values = COLORS$atac[c("ATAC Up", "ATAC Down")],
                      name = "ATAC Direction") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Chromatin State Distribution of Differential ATAC Peaks",
      x = "Chromatin State", y = "% of Peaks"
    ) +
    theme_biomodal() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
      legend.position = "top"
    )

  save_multiformat_ggplot(p_13a, file.path(OUTPUT_DIR, "13a_atac_chromatin_state_distribution"),
                          width = 11, height = 8)

  # Save table
  state_table <- state_dist %>%
    dplyr::select(state, direction, count, pct) %>%
    arrange(state, direction)
  write.table(state_table, file.path(TABLES_DIR, "atac_chromatin_state_distribution.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: atac_chromatin_state_distribution.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 13b: Direct ChIP Mark Overlap Enrichment
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 13b: ChIP mark overlap enrichment...\n")

  mark_names <- c("CTCF", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "Bivalent")
  overlap_cols <- c("CTCF_overlap", "H3K27ac_overlap", "H3K27me3_overlap",
                    "H3K4me1_overlap", "H3K4me3_overlap", "Bivalent_overlap")

  fisher_results <- data.frame(
    mark = character(),
    atac_up_overlap = numeric(),
    atac_up_total = numeric(),
    atac_up_pct = numeric(),
    atac_down_overlap = numeric(),
    atac_down_total = numeric(),
    atac_down_pct = numeric(),
    odds_ratio = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(mark_names)) {
    up_ov <- sum(overlaps_up[[overlap_cols[i]]])
    up_no <- sum(!overlaps_up[[overlap_cols[i]]])
    dn_ov <- sum(overlaps_down[[overlap_cols[i]]])
    dn_no <- sum(!overlaps_down[[overlap_cols[i]]])

    ft <- fisher.test(matrix(c(dn_ov, dn_no, up_ov, up_no), nrow = 2))
    fisher_results <- rbind(fisher_results, data.frame(
      mark = mark_names[i],
      atac_up_overlap = up_ov,
      atac_up_total = length(atac_up_gr),
      atac_up_pct = 100 * up_ov / length(atac_up_gr),
      atac_down_overlap = dn_ov,
      atac_down_total = length(atac_down_gr),
      atac_down_pct = 100 * dn_ov / length(atac_down_gr),
      odds_ratio = as.numeric(ft$estimate),
      pvalue = ft$p.value,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %s: Up=%.1f%%, Down=%.1f%%, OR=%.2f (Down/Up), p=%.2e\n",
                mark_names[i],
                100 * up_ov / length(atac_up_gr),
                100 * dn_ov / length(atac_down_gr),
                as.numeric(ft$estimate), ft$p.value))
  }

  fisher_results$padj <- p.adjust(fisher_results$pvalue, method = "BH")
  fisher_results$sig_label <- ifelse(fisher_results$padj < 0.001, "***",
                               ifelse(fisher_results$padj < 0.01, "**",
                               ifelse(fisher_results$padj < 0.05, "*", "ns")))

  # Build grouped bar chart
  bar_df <- rbind(
    data.frame(mark = fisher_results$mark, direction = "ATAC Up",
               pct = fisher_results$atac_up_pct, stringsAsFactors = FALSE),
    data.frame(mark = fisher_results$mark, direction = "ATAC Down",
               pct = fisher_results$atac_down_pct, stringsAsFactors = FALSE)
  )
  bar_df$mark <- factor(bar_df$mark, levels = mark_names)
  bar_df$direction <- factor(bar_df$direction, levels = c("ATAC Up", "ATAC Down"))

  # Add significance labels at the max height of each mark pair
  sig_labels <- fisher_results %>%
    mutate(y_pos = pmax(atac_up_pct, atac_down_pct) + 2)

  p_13b <- ggplot(bar_df, aes(x = mark, y = pct, fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7),
             width = 0.6, color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
    geom_text(data = sig_labels, aes(x = mark, y = y_pos, label = sig_label),
              inherit.aes = FALSE, size = 4, fontface = "bold", color = "red") +
    scale_fill_manual(values = COLORS$atac[c("ATAC Up", "ATAC Down")],
                      name = "ATAC Direction") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title = "ChIP-seq Mark Overlap: ATAC Up vs ATAC Down Peaks",
      subtitle = "Fisher's exact test per mark (ATAC Down / ATAC Up), BH-adjusted | *p<0.05, **p<0.01, ***p<0.001",
      x = "ChIP-seq Mark", y = "% of Peaks Overlapping Mark"
    ) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p_13b, file.path(OUTPUT_DIR, "13b_atac_chip_overlap_enrichment"),
                          width = 10, height = 8)

  write.table(fisher_results, file.path(TABLES_DIR, "atac_chip_overlap_enrichment.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Saved: atac_chip_overlap_enrichment.tsv\n")

  # -----------------------------------------------------------------------
  # FIGURE 13c: Enrichment Heatmap - ATAC Direction x Chromatin State
  # -----------------------------------------------------------------------
  cat("\nCreating Figure 13c: ATAC direction x chromatin state enrichment heatmap...\n")

  # Build 2x7 contingency data
  n_up <- length(atac_up_gr)
  n_down <- length(atac_down_gr)
  n_total <- n_up + n_down

  heatmap_rows <- list()
  for (s in CHROMATIN_STATE_ORDER) {
    obs_up <- sum(state_up == s)
    obs_dn <- sum(state_down == s)

    # Expected under independence
    total_state <- obs_up + obs_dn
    exp_up <- total_state * (n_up / n_total)
    exp_dn <- total_state * (n_down / n_total)

    # Enrichment
    enrich_up <- ifelse(exp_up > 0, obs_up / exp_up, NA)
    enrich_dn <- ifelse(exp_dn > 0, obs_dn / exp_dn, NA)

    # Fisher test: (in_state vs not_in_state) x (up vs down)
    ft <- fisher.test(matrix(c(obs_dn, n_down - obs_dn, obs_up, n_up - obs_up), nrow = 2))

    heatmap_rows[[length(heatmap_rows) + 1]] <- data.frame(
      state = s, direction = "ATAC Up",
      observed = obs_up, expected = round(exp_up, 1),
      enrichment = round(enrich_up, 3),
      pvalue = ft$p.value, stringsAsFactors = FALSE
    )
    heatmap_rows[[length(heatmap_rows) + 1]] <- data.frame(
      state = s, direction = "ATAC Down",
      observed = obs_dn, expected = round(exp_dn, 1),
      enrichment = round(enrich_dn, 3),
      pvalue = ft$p.value, stringsAsFactors = FALSE
    )
  }
  heatmap_data <- do.call(rbind, heatmap_rows)
  heatmap_data$padj <- p.adjust(heatmap_data$pvalue, method = "BH")
  heatmap_data$sig <- heatmap_data$padj < 0.05

  heatmap_data$state <- factor(heatmap_data$state, levels = CHROMATIN_STATE_ORDER)
  heatmap_data$direction <- factor(heatmap_data$direction, levels = c("ATAC Up", "ATAC Down"))

  heatmap_data$label <- sprintf("O/E: %.2f\n(n=%d)", heatmap_data$enrichment, heatmap_data$observed)

  p_13c <- ggplot(heatmap_data, aes(x = state, y = direction, fill = enrichment)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = label), size = 3.5, lineheight = 1.1) +
    # Bold border for significant cells
    geom_tile(data = heatmap_data %>% dplyr::filter(sig),
              aes(x = state, y = direction),
              fill = NA, color = "black", linewidth = 1.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "Enrichment\n(Obs/Exp)",
                         limits = c(
                           min(heatmap_data$enrichment, na.rm = TRUE) * 0.9,
                           max(heatmap_data$enrichment, na.rm = TRUE) * 1.1
                         )) +
    labs(
      title = "ATAC Direction x Chromatin State Enrichment",
      subtitle = "Observed/Expected ratio | Black borders = significant (BH-adjusted p < 0.05)",
      x = "Chromatin State", y = "ATAC Direction"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  save_multiformat_ggplot(p_13c, file.path(OUTPUT_DIR, "13c_atac_chromatin_enrichment_heatmap"),
                          width = 12, height = 6)

  # =====================================================================
  # PART B: ATAC at Hi-C Loop Anchors (Figures 13d-f)
  # =====================================================================

  cat("\n--- PART B: ATAC at Hi-C Loop Anchors ---\n\n")

  # -----------------------------------------------------------------------
  # Load loop annotation data
  # -----------------------------------------------------------------------

  if (!file.exists(LOOP_FILES$late)) {
    cat("  ERROR: Loop annotation file not found. Skipping Part B.\n")
  } else {

    cat("Loading Hi-C loop annotations...\n")
    loops <- read.table(LOOP_FILES$late, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
    cat(sprintf("  Loaded %d differential loops\n", nrow(loops)))
    cat(sprintf("    up_in_mutant:   %d\n", sum(loops$direction == "up_in_mutant")))
    cat(sprintf("    down_in_mutant: %d\n", sum(loops$direction == "down_in_mutant")))

    # Extract anchor GRanges
    anchor1_gr <- GRanges(seqnames = loops$chr1,
                          ranges = IRanges(start = loops$start1, end = loops$end1))
    anchor2_gr <- GRanges(seqnames = loops$chr2,
                          ranges = IRanges(start = loops$start2, end = loops$end2))

    # Compute ATAC overlap at each anchor
    loops$anchor1_atac_up   <- countOverlaps(anchor1_gr, atac_up_gr) > 0
    loops$anchor1_atac_down <- countOverlaps(anchor1_gr, atac_down_gr) > 0
    loops$anchor2_atac_up   <- countOverlaps(anchor2_gr, atac_up_gr) > 0
    loops$anchor2_atac_down <- countOverlaps(anchor2_gr, atac_down_gr) > 0

    # Pool anchors: any anchor with ATAC overlap
    loops$any_atac_up   <- loops$anchor1_atac_up | loops$anchor2_atac_up
    loops$any_atac_down <- loops$anchor1_atac_down | loops$anchor2_atac_down

    cat(sprintf("\n  Anchor ATAC overlap summary:\n"))
    cat(sprintf("    Loops with any anchor ATAC Up:   %d/%d (%.1f%%)\n",
                sum(loops$any_atac_up), nrow(loops),
                100 * sum(loops$any_atac_up) / nrow(loops)))
    cat(sprintf("    Loops with any anchor ATAC Down: %d/%d (%.1f%%)\n",
                sum(loops$any_atac_down), nrow(loops),
                100 * sum(loops$any_atac_down) / nrow(loops)))

    # -----------------------------------------------------------------------
    # FIGURE 13d: ATAC Overlap at Loop Anchors by Loop Direction
    # -----------------------------------------------------------------------
    cat("\nCreating Figure 13d: ATAC overlap at loop anchors by loop direction...\n")

    up_loops <- loops[loops$direction == "up_in_mutant", ]
    dn_loops <- loops[loops$direction == "down_in_mutant", ]

    # Compute overlap rates
    overlap_stats <- data.frame(
      loop_direction = rep(c("Up in Mutant", "Down in Mutant"), each = 2),
      atac_direction = rep(c("ATAC Up", "ATAC Down"), 2),
      count = c(
        sum(up_loops$any_atac_up), sum(up_loops$any_atac_down),
        sum(dn_loops$any_atac_up), sum(dn_loops$any_atac_down)
      ),
      total = rep(c(nrow(up_loops), nrow(dn_loops)), each = 2),
      stringsAsFactors = FALSE
    )
    overlap_stats$pct <- 100 * overlap_stats$count / overlap_stats$total

    overlap_stats$loop_direction <- factor(overlap_stats$loop_direction,
                                           levels = c("Up in Mutant", "Down in Mutant"))
    overlap_stats$atac_direction <- factor(overlap_stats$atac_direction,
                                           levels = c("ATAC Up", "ATAC Down"))

    for (i in 1:nrow(overlap_stats)) {
      cat(sprintf("  %s loops, %s: %d/%d (%.1f%%)\n",
                  overlap_stats$loop_direction[i], overlap_stats$atac_direction[i],
                  overlap_stats$count[i], overlap_stats$total[i],
                  overlap_stats$pct[i]))
    }

    # Fisher's test: loop direction x anchor ATAC-up
    fisher_up <- fisher.test(matrix(
      c(sum(up_loops$any_atac_up), nrow(up_loops) - sum(up_loops$any_atac_up),
        sum(dn_loops$any_atac_up), nrow(dn_loops) - sum(dn_loops$any_atac_up)),
      nrow = 2
    ))
    # Fisher's test: loop direction x anchor ATAC-down
    fisher_dn <- fisher.test(matrix(
      c(sum(dn_loops$any_atac_down), nrow(dn_loops) - sum(dn_loops$any_atac_down),
        sum(up_loops$any_atac_down), nrow(up_loops) - sum(up_loops$any_atac_down)),
      nrow = 2
    ))

    cat(sprintf("  Fisher (up loops vs ATAC Up): OR=%.2f, p=%.2e\n",
                fisher_up$estimate, fisher_up$p.value))
    cat(sprintf("  Fisher (down loops vs ATAC Down): OR=%.2f, p=%.2e\n",
                fisher_dn$estimate, fisher_dn$p.value))

    p_13d <- ggplot(overlap_stats, aes(x = loop_direction, y = pct, fill = atac_direction)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.6, color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, count)),
                position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
      scale_fill_manual(values = COLORS$atac[c("ATAC Up", "ATAC Down")],
                        name = "Anchor ATAC Change") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(
        title = "ATAC-seq Overlap at Hi-C Loop Anchors",
        subtitle = sprintf(
          "Fisher (up loops/ATAC Up): OR=%.2f, p=%.2e | Fisher (down loops/ATAC Down): OR=%.2f, p=%.2e",
          fisher_up$estimate, fisher_up$p.value,
          fisher_dn$estimate, fisher_dn$p.value
        ),
        x = "Loop Direction (Mutant vs Control)", y = "% of Loops with Anchor ATAC Overlap"
      ) +
      theme_biomodal() +
      theme(legend.position = "top")

    save_multiformat_ggplot(p_13d, file.path(OUTPUT_DIR, "13d_loop_anchor_atac_overlap"),
                            width = 10, height = 8)

    # Save per-loop ATAC overlap table
    loop_export <- loops %>%
      dplyr::select(loop_id, chr1, start1, end1, chr2, start2, end2,
                    logFC, FDR, direction, loop_type_extended,
                    anchor1_type_extended, anchor2_type_extended,
                    anchor1_atac_up, anchor1_atac_down,
                    anchor2_atac_up, anchor2_atac_down,
                    any_atac_up, any_atac_down)

    write.table(loop_export, file.path(TABLES_DIR, "loop_anchor_atac_overlap.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Saved: loop_anchor_atac_overlap.tsv\n")

    # -----------------------------------------------------------------------
    # FIGURE 13e: ATAC Overlap by Anchor Chromatin State
    # -----------------------------------------------------------------------
    cat("\nCreating Figure 13e: ATAC overlap by anchor chromatin state...\n")

    # Use pre-computed anchor types from the loop TSV
    # Include CTCF_Site as the 8th category from the loop annotation pipeline
    anchor_state_order <- c(CHROMATIN_STATE_ORDER, "CTCF_Site")

    # Build per-anchor data (pool anchor1 + anchor2)
    anchor_data <- rbind(
      data.frame(
        anchor_type = loops$anchor1_type_extended,
        loop_direction = loops$direction,
        atac_up = loops$anchor1_atac_up,
        atac_down = loops$anchor1_atac_down,
        stringsAsFactors = FALSE
      ),
      data.frame(
        anchor_type = loops$anchor2_type_extended,
        loop_direction = loops$direction,
        atac_up = loops$anchor2_atac_up,
        atac_down = loops$anchor2_atac_down,
        stringsAsFactors = FALSE
      )
    )

    # Filter to valid anchor types
    anchor_data <- anchor_data[anchor_data$anchor_type %in% anchor_state_order, ]

    # Summarize ATAC overlap per anchor type
    anchor_summary <- anchor_data %>%
      group_by(anchor_type) %>%
      summarise(
        n_anchors = n(),
        n_atac_up = sum(atac_up),
        n_atac_down = sum(atac_down),
        pct_atac_up = 100 * sum(atac_up) / n(),
        pct_atac_down = 100 * sum(atac_down) / n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_anchors >= 20)

    cat("  Anchor ATAC overlap by chromatin state:\n")
    for (i in 1:nrow(anchor_summary)) {
      cat(sprintf("    %s (n=%d): ATAC Up=%.1f%%, ATAC Down=%.1f%%\n",
                  anchor_summary$anchor_type[i], anchor_summary$n_anchors[i],
                  anchor_summary$pct_atac_up[i], anchor_summary$pct_atac_down[i]))
    }

    # Reshape for grouped bar chart
    anchor_bar <- anchor_summary %>%
      tidyr::pivot_longer(cols = c(pct_atac_up, pct_atac_down),
                          names_to = "atac_dir", values_to = "pct") %>%
      mutate(
        atac_direction = ifelse(atac_dir == "pct_atac_up", "ATAC Up", "ATAC Down"),
        anchor_type = factor(anchor_type,
                             levels = intersect(anchor_state_order, anchor_summary$anchor_type))
      )

    # Define colors for anchor types (extend CHROMATIN_STATE_COLORS with CTCF_Site)
    anchor_type_colors <- c(CHROMATIN_STATE_COLORS, "CTCF_Site" = "#1b9e77")

    p_13e <- ggplot(anchor_bar, aes(x = anchor_type, y = pct, fill = atac_direction)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7),
               width = 0.6, color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f%%", pct)),
                position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
      scale_fill_manual(values = COLORS$atac[c("ATAC Up", "ATAC Down")],
                        name = "ATAC Direction") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(
        title = "ATAC Overlap at Loop Anchors by Chromatin State",
        x = "Anchor Chromatin State", y = "% of Anchors Overlapping ATAC Peaks"
      ) +
      theme_biomodal() +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
        legend.position = "top"
      )

    save_multiformat_ggplot(p_13e, file.path(OUTPUT_DIR, "13e_anchor_atac_by_chromatin_state"),
                            width = 12, height = 8)

    # -----------------------------------------------------------------------
    # FIGURE 13f: Loop Type Enrichment for ATAC Concordance
    # -----------------------------------------------------------------------
    cat("\nCreating Figure 13f: Loop type ATAC concordance...\n")

    # Define ATAC concordance: at least one anchor overlaps the expected ATAC direction
    # up_in_mutant -> expected ATAC Up at anchor
    # down_in_mutant -> expected ATAC Down at anchor
    loops$atac_concordant <- ifelse(
      loops$direction == "up_in_mutant",
      loops$any_atac_up,
      loops$any_atac_down
    )

    # Also track discordant: anchor has opposite ATAC direction
    loops$atac_discordant <- ifelse(
      loops$direction == "up_in_mutant",
      loops$any_atac_down & !loops$any_atac_up,
      loops$any_atac_up & !loops$any_atac_down
    )

    overall_concordance <- 100 * sum(loops$atac_concordant) / nrow(loops)
    cat(sprintf("  Overall ATAC concordance: %d/%d (%.1f%%)\n",
                sum(loops$atac_concordant), nrow(loops), overall_concordance))

    # Compute concordance per loop type (filter to types with >= 10 loops)
    concordance_by_type <- loops %>%
      group_by(loop_type_extended) %>%
      summarise(
        n_loops = n(),
        n_concordant = sum(atac_concordant),
        n_discordant = sum(atac_discordant),
        concordance_pct = 100 * sum(atac_concordant) / n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_loops >= 10) %>%
      arrange(desc(concordance_pct))

    cat(sprintf("  %d loop types with >= 10 loops\n", nrow(concordance_by_type)))
    cat("  Top 10 by concordance:\n")
    for (i in 1:min(10, nrow(concordance_by_type))) {
      cat(sprintf("    %s: %.1f%% (n=%d)\n",
                  concordance_by_type$loop_type_extended[i],
                  concordance_by_type$concordance_pct[i],
                  concordance_by_type$n_loops[i]))
    }

    # Bar chart showing concordance % per loop type
    concordance_by_type$loop_type_extended <- factor(
      concordance_by_type$loop_type_extended,
      levels = rev(concordance_by_type$loop_type_extended)
    )

    # Color by concordance rate
    p_13f <- ggplot(concordance_by_type,
                    aes(x = loop_type_extended, y = concordance_pct, fill = concordance_pct)) +
      geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.0f%% (n=%d)", concordance_pct, n_loops)),
                hjust = -0.05, size = 2.8, color = "grey30") +
      geom_hline(yintercept = overall_concordance, linetype = "dashed",
                 color = "red", linewidth = 0.6) +
      annotate("text", x = nrow(concordance_by_type) * 0.95, y = overall_concordance + 1.5,
               label = sprintf("Overall: %.1f%%", overall_concordance),
               size = 3.5, color = "red", fontface = "italic") +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = overall_concordance,
                           name = "Concordance %") +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
      labs(
        title = "ATAC-Loop Direction Concordance by Loop Type",
        subtitle = "Concordant = anchor has ATAC peak matching loop direction (up loop/ATAC up, down loop/ATAC down)",
        x = "Loop Type", y = "% of Loops with ATAC Concordance"
      ) +
      theme_biomodal() +
      theme(legend.position = "right")

    save_multiformat_ggplot(p_13f, file.path(OUTPUT_DIR, "13f_loop_atac_concordance_by_type"),
                            width = 13, height = max(6, nrow(concordance_by_type) * 0.35 + 2))

    # Save concordance table
    write.table(concordance_by_type, file.path(TABLES_DIR, "loop_atac_concordance_by_type.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Saved: loop_atac_concordance_by_type.tsv\n")

  }  # end loop file check

  # -----------------------------------------------------------------------
  # Print summary
  # -----------------------------------------------------------------------
  cat("\n")
  cat("================================================================================\n")
  cat("SECTION 13 SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("ATAC Up peaks: %d\n", length(atac_up_gr)))
  cat(sprintf("ATAC Down peaks: %d\n", length(atac_down_gr)))
  cat("\nPart A: Chromatin State Enrichment\n")
  for (s in c("Polycomb", "Repressed_Promoter", "Active_Enhancer", "Active_Promoter")) {
    pct_up <- state_counts_up$pct[state_counts_up$state == s]
    pct_dn <- state_counts_down$pct[state_counts_down$state == s]
    if (length(pct_up) > 0 && length(pct_dn) > 0) {
      cat(sprintf("  %s: ATAC Up=%.1f%% vs Down=%.1f%% (%s)\n",
                  s, pct_up, pct_dn,
                  ifelse(pct_dn > pct_up, "Down enriched", "Up enriched")))
    }
  }

  h3k27me3_row <- fisher_results[fisher_results$mark == "H3K27me3", ]
  if (nrow(h3k27me3_row) > 0) {
    cat(sprintf("\n  H3K27me3 enrichment in ATAC Down vs Up: OR=%.2f, p=%.2e\n",
                h3k27me3_row$odds_ratio, h3k27me3_row$pvalue))
  }

  if (file.exists(LOOP_FILES$late)) {
    cat(sprintf("\nPart B: Loop Anchor Accessibility\n"))
    cat(sprintf("  Differential loops: %d (up=%d, down=%d)\n",
                nrow(loops),
                sum(loops$direction == "up_in_mutant"),
                sum(loops$direction == "down_in_mutant")))
    cat(sprintf("  Overall ATAC concordance: %.1f%%\n", overall_concordance))
  }

}  # end ATAC file check

cat("\nSection 13 complete.\n\n")
