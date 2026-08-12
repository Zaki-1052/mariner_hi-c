# biomodal/downstream/scripts/viz_sections/section_92_noncg_method_comparison.R
# Section 92: Non-CG Methylation Method Comparison (WGBS vs EM-seq vs evoC + Ecker)
#
# Compares non-CG (CHG+CHH) methylation signal across four methods at 10kb
# genomic bins. WGBS, EM-seq, and evoC were run on the same 4 BAP1-KO/ctrl
# adult cerebellum samples. Ecker provides a deep-seq wildtype reference.
#
# evoC reports CHG and CHH separately; they are averaged per bin to produce
# a combined non-CG track comparable to the WGBS/EM-seq combined signal.
#
#   Panel 92a: Signal distribution (4 methods)
#   Panel 92b: Spatial concordance (6 pairwise hex-scatter)
#   Panel 92c: Coverage overlap (3-way Venn: WGBS/EM-seq/evoC)
#   Panel 92d: Ctrl vs Mut within each method (3 methods)
#   Panel 92e: Chromosomal distribution (3 methods)
#   Panel 92f: Gene body signal (4 methods incl. Ecker)
#   Panel 92g: Neuronal marker gene signal (3 methods)
#   Panel 92h: Composite
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_92_noncg_method_comparison.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 92 CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 92: NON-CG METHOD COMPARISON\n")
cat("================================================================================\n\n")

SEC92_DIR <- file.path(OUTPUT_DIR, "92_noncg_method_comparison")
dir.create(SEC92_DIR, recursive = TRUE, showWarnings = FALSE)

METHOD_COLORS <- c("WGBS" = "#4DAF4A", "EM-seq" = "#377EB8",
                   "evoC" = "#E41A1C", "Ecker" = "#984EA3")
AUTOSOMES <- paste0("chr", c(1:19, "X"))
BIN_WIDTH <- 10000

# ---- BigWig paths (WGBS/EM-seq are local to this section) --------------------

NONCG_BIGWIGS <- list(
  wgbs_all   = file.path(REPO_ROOT, "biomodal/wgbs/nonCG_all-wgbs.bw"),
  wgbs_ctrl  = file.path(REPO_ROOT, "biomodal/wgbs/nonCG_ctrl-wgbs.bw"),
  wgbs_mut   = file.path(REPO_ROOT, "biomodal/wgbs/nonCG_mut-wgbs.bw"),
  emseq_all  = file.path(REPO_ROOT, "biomodal/emseq/nonCG_all-emseq.bw"),
  emseq_ctrl = file.path(REPO_ROOT, "biomodal/emseq/nonCG_ctrl-emseq.bw"),
  emseq_mut  = file.path(REPO_ROOT, "biomodal/emseq/nonCG_mut-emseq.bw")
)

cat("Validating BigWig files...\n")
for (nm in names(NONCG_BIGWIGS)) {
  if (!file.exists(NONCG_BIGWIGS[[nm]])) {
    stop("Required BigWig not found: ", NONCG_BIGWIGS[[nm]])
  }
  cat(sprintf("  %s: OK\n", nm))
}

evoc_bw_paths <- list(
  chg_ctrl = METHYLATION_BIGWIGS$chg_mc_ctrl,
  chg_mut  = METHYLATION_BIGWIGS$chg_mc_mut,
  chh_ctrl = METHYLATION_BIGWIGS$chh_mc_ctrl,
  chh_mut  = METHYLATION_BIGWIGS$chh_mc_mut
)
for (nm in names(evoc_bw_paths)) {
  if (!file.exists(evoc_bw_paths[[nm]])) {
    stop("evoC BigWig not found: ", evoc_bw_paths[[nm]])
  }
  cat(sprintf("  evoc_%s: OK\n", nm))
}

ECKER_AVAILABLE <- file.exists(ECKER_BIGWIGS$ch)
if (ECKER_AVAILABLE) {
  cat(sprintf("  ecker_ch: OK (%.1f GB)\n", file.info(ECKER_BIGWIGS$ch)$size / 1e9))
} else {
  cat("  ecker_ch: NOT FOUND — Ecker panels will be skipped\n")
}
cat("\n")

# =============================================================================
# GENOME TILING
# =============================================================================

cat("Creating 10kb genome tiles...\n")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
sl <- seqlengths(txdb)
sl <- sl[names(sl) %in% AUTOSOMES]
tiles <- tileGenome(sl, tilewidth = BIN_WIDTH, cut.last.tile.in.chrom = TRUE)
cat(sprintf("  %d tiles across %d chromosomes\n", length(tiles), length(sl)))

# =============================================================================
# SIGNAL EXTRACTION
# =============================================================================

extract_binned_signal <- function(bw_path, tiles_gr, label) {
  cat(sprintf("  %s...", label))
  bw <- rtracklayer::import.bw(bw_path, which = tiles_gr)
  cov <- coverage(bw, weight = "score")
  ranges_by_chr <- split(ranges(tiles_gr), seqnames(tiles_gr))
  shared_chrs <- intersect(names(cov), names(ranges_by_chr))
  v <- Views(cov[shared_chrs], ranges_by_chr[shared_chrs])
  means <- viewMeans(v)
  result <- rep(NA_real_, length(tiles_gr))
  chr_names <- as.character(seqnames(tiles_gr))
  for (chr in names(means)) {
    idx <- which(chr_names == chr)
    result[idx] <- as.numeric(means[[chr]])
  }
  nonzero <- sum(result > 0, na.rm = TRUE)
  cat(sprintf(" median=%.6f, non-zero=%d/%d\n",
              median(result, na.rm = TRUE), nonzero, length(result)))
  result
}

cat("\n--- Extracting BAP1 method signals at 10kb bins ---\n")
tiles_df <- data.frame(
  chr   = as.character(seqnames(tiles)),
  start = start(tiles),
  end   = end(tiles),
  stringsAsFactors = FALSE
)

tiles_df$wgbs_all   <- extract_binned_signal(NONCG_BIGWIGS$wgbs_all, tiles, "WGBS all")
tiles_df$wgbs_ctrl  <- extract_binned_signal(NONCG_BIGWIGS$wgbs_ctrl, tiles, "WGBS ctrl")
tiles_df$wgbs_mut   <- extract_binned_signal(NONCG_BIGWIGS$wgbs_mut, tiles, "WGBS mut")
tiles_df$emseq_all  <- extract_binned_signal(NONCG_BIGWIGS$emseq_all, tiles, "EM-seq all")
tiles_df$emseq_ctrl <- extract_binned_signal(NONCG_BIGWIGS$emseq_ctrl, tiles, "EM-seq ctrl")
tiles_df$emseq_mut  <- extract_binned_signal(NONCG_BIGWIGS$emseq_mut, tiles, "EM-seq mut")

cat("\n--- Extracting evoC CHG/CHH signals ---\n")
evoc_chg_ctrl <- extract_binned_signal(evoc_bw_paths$chg_ctrl, tiles, "evoC CHG ctrl")
evoc_chg_mut  <- extract_binned_signal(evoc_bw_paths$chg_mut, tiles, "evoC CHG mut")
evoc_chh_ctrl <- extract_binned_signal(evoc_bw_paths$chh_ctrl, tiles, "evoC CHH ctrl")
evoc_chh_mut  <- extract_binned_signal(evoc_bw_paths$chh_mut, tiles, "evoC CHH mut")

# Average CHG + CHH to match WGBS/EM-seq combined non-CG.
# Not identical to a weighted average across all non-CG cytosines,
# but at <1% methylation the difference is negligible.
tiles_df$evoc_ctrl <- (evoc_chg_ctrl + evoc_chh_ctrl) / 2
tiles_df$evoc_mut  <- (evoc_chg_mut + evoc_chh_mut) / 2
tiles_df$evoc_all  <- (tiles_df$evoc_ctrl + tiles_df$evoc_mut) / 2

if (ECKER_AVAILABLE) {
  cat("\n--- Extracting Ecker CH signal (6.1 GB — loading per-chromosome) ---\n")
  ecker_result <- rep(NA_real_, length(tiles))
  chr_names <- as.character(seqnames(tiles))
  for (chr in AUTOSOMES) {
    chr_idx <- which(chr_names == chr)
    if (length(chr_idx) == 0) next
    cat(sprintf("  Ecker %s (%d bins)...", chr, length(chr_idx)))
    chr_tiles <- tiles[chr_idx]
    bw <- rtracklayer::import.bw(ECKER_BIGWIGS$ch, which = chr_tiles)
    if (length(bw) == 0) {
      cat(" no data\n")
      next
    }
    cov <- coverage(bw, weight = "score")
    chr_ranges <- ranges(chr_tiles)
    if (chr %in% names(cov)) {
      v <- Views(cov[[chr]], chr_ranges)
      ecker_result[chr_idx] <- as.numeric(viewMeans(v))
    }
    rm(bw, cov)
    gc(verbose = FALSE)
    nonzero <- sum(ecker_result[chr_idx] > 0, na.rm = TRUE)
    cat(sprintf(" non-zero=%d/%d\n", nonzero, length(chr_idx)))
  }
  tiles_df$ecker <- ecker_result
  cat(sprintf("  Ecker total: median=%.6f, non-zero=%d/%d\n",
              median(ecker_result, na.rm = TRUE),
              sum(ecker_result > 0, na.rm = TRUE), sum(!is.na(ecker_result))))
}

cat("\n--- Signal extraction summary ---\n")
for (method in c("wgbs_all", "emseq_all", "evoc_all")) {
  vals <- tiles_df[[method]]
  n_nonzero <- sum(vals > 0, na.rm = TRUE)
  cat(sprintf("  %-12s  detected: %d / %d (%.1f%%),  median(>0): %.6f\n",
              method, n_nonzero, sum(!is.na(vals)),
              100 * n_nonzero / sum(!is.na(vals)),
              median(vals[vals > 0], na.rm = TRUE)))
}
if (ECKER_AVAILABLE) {
  vals <- tiles_df$ecker
  n_nonzero <- sum(vals > 0, na.rm = TRUE)
  cat(sprintf("  %-12s  detected: %d / %d (%.1f%%),  median(>0): %.6f\n",
              "ecker", n_nonzero, sum(!is.na(vals)),
              100 * n_nonzero / sum(!is.na(vals)),
              median(vals[vals > 0], na.rm = TRUE)))
}
cat("\n")

# =============================================================================
# PANEL 92a: SIGNAL DISTRIBUTION COMPARISON
# =============================================================================

cat("================================================================================\n")
cat("Panel 92a: Signal Distribution Comparison\n")
cat("================================================================================\n\n")

dist_methods <- c("WGBS" = "wgbs_all", "EM-seq" = "emseq_all", "evoC" = "evoc_all")
if (ECKER_AVAILABLE) dist_methods <- c(dist_methods, "Ecker" = "ecker")

dist_long <- do.call(rbind, lapply(names(dist_methods), function(nm) {
  vals <- tiles_df[[dist_methods[nm]]]
  data.frame(method = nm, signal = vals[!is.na(vals)], stringsAsFactors = FALSE)
}))
dist_long$method <- factor(dist_long$method, levels = names(dist_methods))

dist_nonzero <- dist_long %>% dplyr::filter(signal > 0)

zero_fracs <- dist_long %>%
  group_by(method) %>%
  summarise(n_total = n(), n_zero = sum(signal == 0),
            pct_zero = 100 * n_zero / n_total, .groups = "drop")

cat("Zero-signal fractions:\n")
for (i in seq_len(nrow(zero_fracs))) {
  cat(sprintf("  %s: %.1f%% zero (%d / %d bins)\n",
              zero_fracs$method[i], zero_fracs$pct_zero[i],
              zero_fracs$n_zero[i], zero_fracs$n_total[i]))
}

dist_medians <- dist_nonzero %>%
  group_by(method) %>%
  summarise(med = median(signal), .groups = "drop")

p_92a <- ggplot(dist_nonzero, aes(x = signal, fill = method, color = method)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(data = dist_medians, aes(xintercept = med, color = method),
             linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = METHOD_COLORS) +
  scale_color_manual(values = METHOD_COLORS) +
  coord_cartesian(xlim = c(0, 0.05)) +
  labs(
    title = "Non-CG Methylation Signal Distribution (bins with signal > 0)",
    subtitle = paste(sapply(seq_len(nrow(zero_fracs)), function(i)
      sprintf("%s: %.0f%% zero", zero_fracs$method[i], zero_fracs$pct_zero[i])),
      collapse = "  |  "),
    x = "Non-CG methylation fraction",
    y = "Density",
    fill = "Method", color = "Method"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_92a,
  file.path(SEC92_DIR, "92a_signal_distribution"), width = 10, height = 7)

# =============================================================================
# PANEL 92b: SPATIAL CONCORDANCE (PAIRWISE SCATTER)
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92b: Spatial Concordance (Pairwise Hex-Scatter)\n")
cat("================================================================================\n\n")

pair_methods <- list(
  c("WGBS", "wgbs_all", "EM-seq", "emseq_all"),
  c("WGBS", "wgbs_all", "evoC", "evoc_all"),
  c("EM-seq", "emseq_all", "evoC", "evoc_all")
)
if (ECKER_AVAILABLE) {
  pair_methods <- c(pair_methods, list(
    c("WGBS", "wgbs_all", "Ecker", "ecker"),
    c("EM-seq", "emseq_all", "Ecker", "ecker"),
    c("evoC", "evoc_all", "Ecker", "ecker")
  ))
}

concordance_results <- list()
scatter_plots <- list()

for (pair in pair_methods) {
  nm_x <- pair[1]; col_x <- pair[2]
  nm_y <- pair[3]; col_y <- pair[4]
  label <- paste0(nm_x, " vs ", nm_y)

  x <- tiles_df[[col_x]]
  y <- tiles_df[[col_y]]
  both_pos <- !is.na(x) & !is.na(y) & (x > 0 | y > 0)
  n_shared <- sum(both_pos)

  if (n_shared < 10) {
    cat(sprintf("  %s: only %d shared bins — skipping\n", label, n_shared))
    concordance_results[[label]] <- data.frame(
      pair = label, n_shared = n_shared,
      spearman_rho = NA, pearson_r = NA, p_value = NA,
      stringsAsFactors = FALSE)
    next
  }

  sp <- cor.test(x[both_pos], y[both_pos], method = "spearman")
  pe <- cor.test(x[both_pos], y[both_pos], method = "pearson")

  concordance_results[[label]] <- data.frame(
    pair = label, n_shared = n_shared,
    spearman_rho = sp$estimate, pearson_r = pe$estimate,
    p_value = sp$p.value, stringsAsFactors = FALSE)

  cat(sprintf("  %s: n=%d, rho=%.3f, r=%.3f\n",
              label, n_shared, sp$estimate, pe$estimate))

  plot_df <- data.frame(x = x[both_pos], y = y[both_pos])
  max_val <- quantile(c(plot_df$x, plot_df$y), 0.99, na.rm = TRUE)

  scatter_plots[[label]] <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(option = "C", trans = "log10") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    labs(
      title = label,
      subtitle = sprintf("n = %s, rho = %.3f", format(n_shared, big.mark = ","), sp$estimate),
      x = paste(nm_x, "non-CG fraction"),
      y = paste(nm_y, "non-CG fraction")
    ) +
    theme_biomodal() +
    theme(legend.position = "right", legend.key.height = unit(12, "mm"))
}

concordance_df <- do.call(rbind, concordance_results)
rownames(concordance_df) <- NULL

if (length(scatter_plots) > 0) {
  ncols <- min(3, length(scatter_plots))
  nrows <- ceiling(length(scatter_plots) / ncols)
  p_92b <- wrap_plots(scatter_plots, ncol = ncols) +
    plot_annotation(title = "Non-CG Spatial Concordance at 10kb Bins")
  save_multiformat_ggplot(p_92b,
    file.path(SEC92_DIR, "92b_spatial_concordance"),
    width = 6 * ncols, height = 5.5 * nrows)
}

# =============================================================================
# PANEL 92c: COVERAGE OVERLAP (VENN)
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92c: Coverage Overlap (Venn Diagram)\n")
cat("================================================================================\n\n")

detected_sets <- list(
  "WGBS"   = which(tiles_df$wgbs_all > 0 & !is.na(tiles_df$wgbs_all)),
  "EM-seq" = which(tiles_df$emseq_all > 0 & !is.na(tiles_df$emseq_all)),
  "evoC"   = which(tiles_df$evoc_all > 0 & !is.na(tiles_df$evoc_all))
)

for (nm in names(detected_sets)) {
  cat(sprintf("  %s: %d bins with signal > 0\n", nm, length(detected_sets[[nm]])))
}

p_92c <- ggVennDiagram(detected_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#377EB8") +
  labs(
    title = "Non-CG Detection Overlap at 10kb Bins",
    subtitle = sprintf("Genome-wide: %s total bins", format(nrow(tiles_df), big.mark = ",")),
    fill = "Bin count"
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_92c,
  file.path(SEC92_DIR, "92c_coverage_overlap"), width = 8, height = 7)

# =============================================================================
# PANEL 92d: CTRL VS MUT WITHIN EACH METHOD
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92d: Ctrl vs Mut Within Each Method\n")
cat("================================================================================\n\n")

cvm_specs <- list(
  list(method = "WGBS",   ctrl = "wgbs_ctrl",  mut = "wgbs_mut"),
  list(method = "EM-seq", ctrl = "emseq_ctrl", mut = "emseq_mut"),
  list(method = "evoC",   ctrl = "evoc_ctrl",  mut = "evoc_mut")
)

cvm_results <- list()
cvm_long <- list()

for (sp in cvm_specs) {
  ctrl_vals <- tiles_df[[sp$ctrl]]
  mut_vals <- tiles_df[[sp$mut]]
  both_valid <- !is.na(ctrl_vals) & !is.na(mut_vals) & (ctrl_vals > 0 | mut_vals > 0)
  n_valid <- sum(both_valid)

  if (n_valid < 10) {
    cat(sprintf("  %s: only %d valid bins — skipping\n", sp$method, n_valid))
    cvm_results[[sp$method]] <- data.frame(
      method = sp$method, median_ctrl = NA, median_mut = NA,
      delta = NA, wilcox_p = NA, direction = "N/A", stringsAsFactors = FALSE)
    next
  }

  med_ctrl <- median(ctrl_vals[both_valid], na.rm = TRUE)
  med_mut <- median(mut_vals[both_valid], na.rm = TRUE)
  wt <- wilcox.test(ctrl_vals[both_valid], mut_vals[both_valid])
  direction <- ifelse(med_mut > med_ctrl, "mut > ctrl", "ctrl > mut")

  cvm_results[[sp$method]] <- data.frame(
    method = sp$method, median_ctrl = med_ctrl, median_mut = med_mut,
    delta = med_mut - med_ctrl, wilcox_p = wt$p.value,
    direction = direction, stringsAsFactors = FALSE)

  cat(sprintf("  %s: ctrl=%.6f, mut=%.6f, delta=%.6f, p=%.2e (%s)\n",
              sp$method, med_ctrl, med_mut, med_mut - med_ctrl,
              wt$p.value, direction))

  cvm_long[[sp$method]] <- rbind(
    data.frame(method = sp$method, condition = "Control",
               signal = ctrl_vals[both_valid], stringsAsFactors = FALSE),
    data.frame(method = sp$method, condition = "Mutant",
               signal = mut_vals[both_valid], stringsAsFactors = FALSE)
  )
}

cvm_df <- do.call(rbind, cvm_results)
rownames(cvm_df) <- NULL

directions <- cvm_df$direction[cvm_df$direction != "N/A"]
direction_consensus <- ifelse(length(unique(directions)) == 1,
  paste("All methods agree:", unique(directions)),
  "Methods disagree on direction")
cat(sprintf("\n  Direction consensus: %s\n", direction_consensus))

cvm_plot_df <- do.call(rbind, cvm_long)
cvm_plot_df$method <- factor(cvm_plot_df$method, levels = c("WGBS", "EM-seq", "evoC"))
cvm_plot_df$condition <- factor(cvm_plot_df$condition, levels = c("Control", "Mutant"))

cvm_annot <- cvm_df %>% dplyr::filter(direction != "N/A")

p_92d <- ggplot(cvm_plot_df, aes(x = signal, fill = condition, color = condition)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  facet_wrap(~method, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Control" = COLORS$condition["Control"],
                                "Mutant" = COLORS$condition["Mutant"])) +
  scale_color_manual(values = c("Control" = COLORS$condition["Control"],
                                 "Mutant" = COLORS$condition["Mutant"])) +
  coord_cartesian(xlim = c(0, 0.03)) +
  labs(
    title = "Ctrl vs Mut Non-CG Methylation by Method",
    subtitle = paste0(direction_consensus,
                      "\nNote: WGBS mut includes underpowered sample mut_3201 (4.3M reads)"),
    x = "Non-CG methylation fraction", y = "Density",
    fill = "Condition", color = "Condition"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_92d,
  file.path(SEC92_DIR, "92d_ctrl_vs_mut"), width = 14, height = 6)

# =============================================================================
# PANEL 92e: CHROMOSOMAL DISTRIBUTION
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92e: Chromosomal Distribution\n")
cat("================================================================================\n\n")

chr_methods <- c("WGBS" = "wgbs_all", "EM-seq" = "emseq_all", "evoC" = "evoc_all")
chr_long <- do.call(rbind, lapply(names(chr_methods), function(nm) {
  tiles_df %>%
    dplyr::filter(!is.na(.data[[chr_methods[nm]]])) %>%
    group_by(chr) %>%
    summarise(
      mean_signal = mean(.data[[chr_methods[nm]]], na.rm = TRUE),
      n_bins = n(),
      n_detected = sum(.data[[chr_methods[nm]]] > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(method = nm)
}))

chr_long$chr <- factor(chr_long$chr, levels = AUTOSOMES)
chr_long$method <- factor(chr_long$method, levels = names(chr_methods))

chr_wide <- chr_long %>%
  dplyr::select(chr, method, mean_signal) %>%
  tidyr::pivot_wider(names_from = method, values_from = mean_signal)

if (nrow(chr_wide) >= 5) {
  rho_wgbs_emseq <- cor(chr_wide$WGBS, chr_wide$`EM-seq`,
                         use = "complete.obs", method = "spearman")
  rho_wgbs_evoc <- cor(chr_wide$WGBS, chr_wide$evoC,
                        use = "complete.obs", method = "spearman")
  cat(sprintf("  Chromosome profile Spearman: WGBS vs EM-seq = %.3f, WGBS vs evoC = %.3f\n",
              rho_wgbs_emseq, rho_wgbs_evoc))
  chr_subtitle <- sprintf("Chr-level Spearman: WGBS-EM-seq = %.2f, WGBS-evoC = %.2f",
                           rho_wgbs_emseq, rho_wgbs_evoc)
} else {
  chr_subtitle <- ""
}

p_92e <- ggplot(chr_long, aes(x = chr, y = mean_signal, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = METHOD_COLORS) +
  labs(
    title = "Mean Non-CG Methylation by Chromosome",
    subtitle = chr_subtitle,
    x = "Chromosome", y = "Mean non-CG fraction",
    fill = "Method"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

save_multiformat_ggplot(p_92e,
  file.path(SEC92_DIR, "92e_chromosomal_distribution"), width = 14, height = 7)

# =============================================================================
# PANEL 92f: GENE BODY SIGNAL
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92f: Gene Body Signal Comparison\n")
cat("================================================================================\n\n")

gene_gr <- genes(txdb)
gene_gr <- gene_gr[seqnames(gene_gr) %in% AUTOSOMES]
cat(sprintf("  %d genes on autosomes + chrX\n", length(gene_gr)))

gene_methods <- c("WGBS" = "wgbs_all", "EM-seq" = "emseq_all", "evoC" = "evoc_all")

gene_signals <- list()
for (nm in names(gene_methods)) {
  bw_path <- if (nm == "WGBS") NONCG_BIGWIGS$wgbs_all
             else if (nm == "EM-seq") NONCG_BIGWIGS$emseq_all
             else NULL

  if (!is.null(bw_path)) {
    gene_signals[[nm]] <- extract_binned_signal(bw_path, gene_gr, paste(nm, "gene bodies"))
  } else {
    # evoC: extract CHG+CHH at gene bodies and average
    chg_sig <- extract_binned_signal(evoc_bw_paths$chg_ctrl, gene_gr, "evoC CHG ctrl genes")
    chh_sig <- extract_binned_signal(evoc_bw_paths$chh_ctrl, gene_gr, "evoC CHH ctrl genes")
    chg_mut_sig <- extract_binned_signal(evoc_bw_paths$chg_mut, gene_gr, "evoC CHG mut genes")
    chh_mut_sig <- extract_binned_signal(evoc_bw_paths$chh_mut, gene_gr, "evoC CHH mut genes")
    gene_signals[[nm]] <- ((chg_sig + chh_sig) / 2 + (chg_mut_sig + chh_mut_sig) / 2) / 2
  }
}

if (ECKER_AVAILABLE) {
  cat("  Ecker gene bodies (per-chromosome to manage memory)...\n")
  ecker_gene_result <- rep(NA_real_, length(gene_gr))
  gene_chr_names <- as.character(seqnames(gene_gr))
  for (chr in AUTOSOMES) {
    chr_idx <- which(gene_chr_names == chr)
    if (length(chr_idx) == 0) next
    chr_genes <- gene_gr[chr_idx]
    bw <- rtracklayer::import.bw(ECKER_BIGWIGS$ch, which = chr_genes)
    if (length(bw) > 0) {
      cov <- coverage(bw, weight = "score")
      if (chr %in% names(cov)) {
        v <- Views(cov[[chr]], ranges(chr_genes))
        ecker_gene_result[chr_idx] <- as.numeric(viewMeans(v))
      }
      rm(bw, cov)
    }
    gc(verbose = FALSE)
  }
  gene_signals[["Ecker"]] <- ecker_gene_result
  nonzero <- sum(ecker_gene_result > 0, na.rm = TRUE)
  cat(sprintf("  Ecker gene bodies: non-zero=%d/%d\n", nonzero, sum(!is.na(ecker_gene_result))))
}

gene_long <- do.call(rbind, lapply(names(gene_signals), function(nm) {
  vals <- gene_signals[[nm]]
  data.frame(method = nm, signal = vals[!is.na(vals)], stringsAsFactors = FALSE)
}))
gene_long$method <- factor(gene_long$method, levels = names(gene_signals))

gene_summary <- list()
for (nm in names(gene_signals)) {
  vals <- gene_signals[[nm]]
  vals_clean <- vals[!is.na(vals)]
  n_det <- sum(vals_clean > 0)
  gene_summary[[nm]] <- data.frame(
    method = nm, n_genes = length(vals_clean),
    n_detected = n_det, pct_detected = 100 * n_det / length(vals_clean),
    median_signal = median(vals_clean), mean_signal = mean(vals_clean),
    median_detected = median(vals_clean[vals_clean > 0]),
    stringsAsFactors = FALSE)
  cat(sprintf("  %s: %d genes, %d detected (%.1f%%), median(>0)=%.6f\n",
              nm, length(vals_clean), n_det,
              100 * n_det / length(vals_clean),
              median(vals_clean[vals_clean > 0], na.rm = TRUE)))
}
gene_summary_df <- do.call(rbind, gene_summary)
rownames(gene_summary_df) <- NULL

gene_nonzero <- gene_long %>% dplyr::filter(signal > 0)

gene_medians <- gene_nonzero %>%
  group_by(method) %>%
  summarise(med = median(signal), n = n(), .groups = "drop")

p_92f <- ggplot(gene_nonzero, aes(x = method, y = signal, fill = method)) +
  geom_violin(alpha = 0.5, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  geom_text(data = gene_medians,
            aes(x = method, y = med, label = sprintf("med=%.5f\nn=%s",
                med, format(n, big.mark = ","))),
            vjust = -0.5, size = 2.8, color = "black") +
  scale_fill_manual(values = METHOD_COLORS) +
  coord_cartesian(ylim = c(0, quantile(gene_nonzero$signal, 0.95))) +
  labs(
    title = "Non-CG Methylation at Gene Bodies (genes with signal > 0)",
    subtitle = paste(sapply(seq_len(nrow(gene_summary_df)), function(i)
      sprintf("%s: %.0f%% detected", gene_summary_df$method[i],
              gene_summary_df$pct_detected[i])), collapse = "  |  "),
    x = "Method", y = "Non-CG methylation fraction"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_92f,
  file.path(SEC92_DIR, "92f_gene_body_signal"), width = 10, height = 7)

# =============================================================================
# PANEL 92g: NEURONAL MARKER GENE SIGNAL
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92g: Neuronal Marker Gene Signal\n")
cat("================================================================================\n\n")

NEURONAL_GENES <- unique(c(KEY_GENES, "Nrxn1", "Nrxn3", "Grin2b", "Dlg2"))

gene_ids <- tryCatch({
  mapping <- AnnotationDbi::select(org.Mm.eg.db,
    keys = NEURONAL_GENES, keytype = "SYMBOL",
    columns = c("ENTREZID", "SYMBOL"))
  mapping <- mapping[!is.na(mapping$ENTREZID), ]
  mapping
}, error = function(e) {
  cat(sprintf("  Warning: gene mapping failed: %s\n", e$message))
  NULL
})

if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
  all_gene_gr <- genes(txdb)
  matched_genes <- all_gene_gr[names(all_gene_gr) %in% gene_ids$ENTREZID]
  matched_genes <- matched_genes[seqnames(matched_genes) %in% AUTOSOMES]

  id_to_sym <- setNames(gene_ids$SYMBOL, gene_ids$ENTREZID)

  neuro_results <- list()
  neuro_methods <- c("WGBS" = "wgbs_all", "EM-seq" = "emseq_all")

  for (nm in names(neuro_methods)) {
    bw_path <- NONCG_BIGWIGS[[neuro_methods[nm]]]
    sig <- extract_binned_signal(bw_path, matched_genes, paste(nm, "neuronal"))
    for (i in seq_along(sig)) {
      gid <- names(matched_genes)[i]
      sym <- ifelse(gid %in% names(id_to_sym), id_to_sym[gid], gid)
      neuro_results[[paste(nm, sym)]] <- data.frame(
        gene = sym, method = nm, signal = sig[i],
        detected = !is.na(sig[i]) & sig[i] > 0, stringsAsFactors = FALSE)
    }
  }

  # evoC: average CHG+CHH at neuronal genes
  chg_ctrl_sig <- extract_binned_signal(evoc_bw_paths$chg_ctrl, matched_genes, "evoC CHG ctrl neuro")
  chh_ctrl_sig <- extract_binned_signal(evoc_bw_paths$chh_ctrl, matched_genes, "evoC CHH ctrl neuro")
  chg_mut_sig <- extract_binned_signal(evoc_bw_paths$chg_mut, matched_genes, "evoC CHG mut neuro")
  chh_mut_sig <- extract_binned_signal(evoc_bw_paths$chh_mut, matched_genes, "evoC CHH mut neuro")
  evoc_neuro <- ((chg_ctrl_sig + chh_ctrl_sig) / 2 + (chg_mut_sig + chh_mut_sig) / 2) / 2

  for (i in seq_along(evoc_neuro)) {
    gid <- names(matched_genes)[i]
    sym <- ifelse(gid %in% names(id_to_sym), id_to_sym[gid], gid)
    neuro_results[[paste("evoC", sym)]] <- data.frame(
      gene = sym, method = "evoC", signal = evoc_neuro[i],
      detected = !is.na(evoc_neuro[i]) & evoc_neuro[i] > 0, stringsAsFactors = FALSE)
  }

  neuro_df <- do.call(rbind, neuro_results)
  rownames(neuro_df) <- NULL
  neuro_df$method <- factor(neuro_df$method, levels = c("WGBS", "EM-seq", "evoC"))

  gene_order <- neuro_df %>%
    dplyr::filter(method == "EM-seq") %>%
    arrange(desc(signal)) %>%
    pull(gene)
  neuro_df$gene <- factor(neuro_df$gene, levels = gene_order)

  neuro_df$signal_display <- ifelse(is.na(neuro_df$signal) | neuro_df$signal == 0,
                                     0, neuro_df$signal)

  p_92g <- ggplot(neuro_df, aes(x = gene, y = signal_display, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = METHOD_COLORS) +
    labs(
      title = "Non-CG Signal at Neuronal Marker Genes",
      subtitle = "Sorted by EM-seq signal. Zero bars = no detectable signal.",
      x = "Gene", y = "Non-CG methylation fraction",
      fill = "Method"
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  save_multiformat_ggplot(p_92g,
    file.path(SEC92_DIR, "92g_neuronal_gene_signal"), width = 12, height = 7)
} else {
  cat("  Skipping panel 92g — gene mapping unavailable\n")
  neuro_df <- data.frame()
}

# =============================================================================
# PANEL 92h: COMPOSITE
# =============================================================================

cat("\n================================================================================\n")
cat("Panel 92h: Composite Figure\n")
cat("================================================================================\n\n")

p_92h <- (p_92a + p_92c) / (p_92d + p_92e) +
  plot_annotation(
    title = "Section 92: Non-CG Methylation Method Comparison",
    subtitle = "WGBS vs EM-seq vs evoC (BAP1-KO cerebellum, mm10)",
    tag_levels = "A"
  )

save_multiformat_ggplot(p_92h,
  file.path(SEC92_DIR, "92h_composite"), width = 20, height = 14)

# =============================================================================
# SAVE SUMMARY TABLES
# =============================================================================

cat("\n================================================================================\n")
cat("Saving summary tables\n")
cat("================================================================================\n\n")

# Method signal summary
signal_summary <- list()
all_methods <- c("WGBS" = "wgbs_all", "EM-seq" = "emseq_all", "evoC" = "evoc_all")
if (ECKER_AVAILABLE) all_methods <- c(all_methods, "Ecker" = "ecker")

for (nm in names(all_methods)) {
  vals <- tiles_df[[all_methods[nm]]]
  vals_clean <- vals[!is.na(vals)]
  n_det <- sum(vals_clean > 0)
  signal_summary[[nm]] <- data.frame(
    method = nm, n_bins = length(vals_clean),
    n_detected = n_det, pct_detected = 100 * n_det / length(vals_clean),
    median = median(vals_clean), mean = mean(vals_clean),
    p95 = quantile(vals_clean, 0.95),
    median_detected = median(vals_clean[vals_clean > 0], na.rm = TRUE),
    stringsAsFactors = FALSE)
}
signal_summary_df <- do.call(rbind, signal_summary)
rownames(signal_summary_df) <- NULL
write.table(signal_summary_df, file.path(TABLES_DIR, "92_method_signal_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved 92_method_signal_summary.tsv\n")

write.table(concordance_df, file.path(TABLES_DIR, "92_pairwise_concordance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved 92_pairwise_concordance.tsv\n")

write.table(chr_long, file.path(TABLES_DIR, "92_chromosomal_signal.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved 92_chromosomal_signal.tsv\n")

write.table(cvm_df, file.path(TABLES_DIR, "92_ctrl_mut_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved 92_ctrl_mut_comparison.tsv\n")

write.table(gene_summary_df, file.path(TABLES_DIR, "92_gene_body_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved 92_gene_body_summary.tsv\n")

if (nrow(neuro_df) > 0) {
  write.table(neuro_df, file.path(TABLES_DIR, "92_neuronal_gene_signal.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved 92_neuronal_gene_signal.tsv\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 92 COMPLETE\n")
cat("================================================================================\n\n")

cat("Key findings:\n")
cat(sprintf("  Detection rates: WGBS %.1f%%, EM-seq %.1f%%, evoC %.1f%%\n",
            signal_summary_df$pct_detected[signal_summary_df$method == "WGBS"],
            signal_summary_df$pct_detected[signal_summary_df$method == "EM-seq"],
            signal_summary_df$pct_detected[signal_summary_df$method == "evoC"]))

wgbs_emseq_row <- concordance_df[concordance_df$pair == "WGBS vs EM-seq", ]
if (nrow(wgbs_emseq_row) > 0 && !is.na(wgbs_emseq_row$spearman_rho)) {
  cat(sprintf("  WGBS-EM-seq concordance: rho = %.3f (n = %s shared bins)\n",
              wgbs_emseq_row$spearman_rho,
              format(wgbs_emseq_row$n_shared, big.mark = ",")))
}

cat(sprintf("  Ctrl vs Mut direction: %s\n", direction_consensus))
cat(sprintf("  Gene body detection: WGBS %.1f%%, EM-seq %.1f%%, evoC %.1f%%\n",
            gene_summary_df$pct_detected[gene_summary_df$method == "WGBS"],
            gene_summary_df$pct_detected[gene_summary_df$method == "EM-seq"],
            gene_summary_df$pct_detected[gene_summary_df$method == "evoC"]))

cat(sprintf("\nOutputs: %s\n", SEC92_DIR))
cat(sprintf("Tables:  %s/92_*.tsv\n", TABLES_DIR))
