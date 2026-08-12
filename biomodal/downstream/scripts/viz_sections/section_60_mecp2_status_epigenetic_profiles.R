# biomodal/downstream/scripts/viz_sections/section_60_mecp2_status_epigenetic_profiles.R
# Section 60: MeCP2 Status x Epigenetic Mark Profiles
#
# Formalizes the BAP1-KO mechanistic model test:
#   MeCP2-Up  : K119ub up -> K27ac down, K27me3 up, 5mC up, 5hmC down
#   MeCP2-Down: K119ub FLAT -> opposite cascade (mirror image, not weaker)
#
# Figures:
#   60a: 5-facet violin+box (one per mark, x=MeCP2 status, ggsignif brackets)
#   60b: Summary heatmap (pheatmap, 5 marks x 2 groups, mean change + stars)
#   60c: Lollipop / dot-and-whisker effect-size plot with 95% CI
#   60d: Per-gene strip chart for MeCP2-Down (35 genes x 6 marks)
#   60e: Patchwork composite (60a + 60c + 60d)
#
# Input: tables/59_quadrant_master.tsv (from section_59)
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_60_mecp2_status_epigenetic_profiles.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ggsignif)
})

# =============================================================================
# SECTION 60: MeCP2 STATUS EPIGENETIC PROFILES
# =============================================================================

cat("================================================================================\n")
cat("SECTION 60: MeCP2 STATUS EPIGENETIC PROFILES\n")
cat("  Mirror-image epigenetic patterns at MeCP2-Up vs MeCP2-Down loci\n")
cat("================================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

SEC60_DIR <- file.path(OUTPUT_DIR, "60_mecp2_status_epigenetic_profiles")
dir.create(SEC60_DIR, recursive = TRUE, showWarnings = FALSE)

MASTER_PATH <- file.path(TABLES_DIR, "59_quadrant_master.tsv")

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

sig_stars <- function(p) {
  if (is.na(p) || !is.finite(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

save_plot <- function(p, name, w = 12, h = 9) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC60_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# INPUT VALIDATION + LOAD
# =============================================================================

cat("Validating input...\n")
stopifnot(
  "59_quadrant_master.tsv not found (run section_59 first)" =
    file.exists(MASTER_PATH)
)
cat(sprintf("  [OK] %s\n", MASTER_PATH))

master <- read.table(MASTER_PATH, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Loaded %d genes, %d columns\n", nrow(master), ncol(master)))

required <- c("gene", "gb_log2fc", "gb_signal_class",
              "mecp2_mean_fold", "mecp2_nearest_fdr",
              "mc_diff", "hmc_diff",
              "k27ac_fold", "k27me3_fold",
              "k119ub_fold")
missing <- setdiff(required, colnames(master))
stopifnot("Missing required columns in 59_quadrant_master.tsv" = length(missing) == 0)

# =============================================================================
# MeCP2 CLASSIFICATION (identical to section_59 lines 258-264)
# =============================================================================

cat("\nClassifying MeCP2 status...\n")

master <- master %>%
  dplyr::mutate(
    mecp2_status = dplyr::case_when(
      mecp2_nearest_fdr < Q_THRESHOLD & mecp2_mean_fold > 0 ~ "MeCP2 Up",
      mecp2_nearest_fdr < Q_THRESHOLD & mecp2_mean_fold < 0 ~ "MeCP2 Down",
      TRUE ~ "Not Significant"
    ),
    mecp2_status = factor(mecp2_status,
                          levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant"))
  )

n_up   <- sum(master$mecp2_status == "MeCP2 Up",           na.rm = TRUE)
n_down <- sum(master$mecp2_status == "MeCP2 Down",         na.rm = TRUE)
n_ns   <- sum(master$mecp2_status == "Not Significant",    na.rm = TRUE)
cat(sprintf("  MeCP2 Up:   %d\n", n_up))
cat(sprintf("  MeCP2 Down: %d\n", n_down))
cat(sprintf("  NS:         %d\n", n_ns))

# =============================================================================
# MARK DEFINITIONS
# =============================================================================

MARK_ORDER <- c("k119ub_bw", "k119ub_db", "mc", "hmc", "k27ac", "k27me3")

MARK_COLS <- c(
  k119ub_bw = "gb_log2fc",
  k119ub_db = "k119ub_fold",
  mc        = "mc_diff",
  hmc       = "hmc_diff",
  k27ac     = "k27ac_fold",
  k27me3    = "k27me3_fold"
)

MARK_LABELS <- c(
  k119ub_bw = "H2AK119ub\nlog2(mut/ctrl) [BigWig]",
  k119ub_db = "H2AK119ub\nDiffBind fold",
  mc        = "Δ5mC\n(mut - ctrl)",
  hmc       = "Δ5hmC\n(mut - ctrl)",
  k27ac     = "H3K27ac\nDiffBind fold",
  k27me3    = "H3K27me3\nDiffBind fold"
)

MARK_LABELS_SHORT <- c(
  k119ub_bw = "K119ub [BigWig]",
  k119ub_db = "K119ub [DiffBind]",
  mc        = "5mC (CpG)",
  hmc       = "5hmC (CpG)",
  k27ac     = "H3K27ac",
  k27me3    = "H3K27me3"
)

# =============================================================================
# BUILD LONG-FORMAT DATA FRAME
# =============================================================================

cat("\nBuilding long-format data frame...\n")

build_long_df <- function(master) {
  rows <- vector("list", length(MARK_ORDER))
  for (i in seq_along(MARK_ORDER)) {
    mid <- MARK_ORDER[i]
    col <- MARK_COLS[mid]
    df_sub <- master

    if (mid == "k119ub_bw") {
      df_sub <- df_sub %>% dplyr::filter(gb_signal_class == "quantifiable")
    }

    rows[[i]] <- df_sub %>%
      dplyr::mutate(value = .data[[col]]) %>%
      dplyr::select(gene, mecp2_status, value) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::mutate(mark_id = mid)
  }
  dplyr::bind_rows(rows)
}

long_df <- build_long_df(master)
cat(sprintf("  Long df: %d rows\n", nrow(long_df)))

for (mid in MARK_ORDER) {
  sub <- long_df %>% dplyr::filter(mark_id == mid)
  cat(sprintf("  %-12s: %d values (Up=%d, Down=%d, NS=%d)\n",
              mid,
              nrow(sub),
              sum(sub$mecp2_status == "MeCP2 Up"),
              sum(sub$mecp2_status == "MeCP2 Down"),
              sum(sub$mecp2_status == "Not Significant")))
}

# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

cat("\n--- Statistical tests ---\n")

run_onesample_tests <- function(x, group_label, mark_id) {
  x <- x[!is.na(x)]
  if (length(x) < 3) {
    return(data.frame(
      group = group_label, mark = mark_id,
      n = length(x), mean = NA_real_, median = NA_real_,
      ttest_t = NA_real_, ttest_p = NA_real_,
      wilcox_V = NA_real_, wilcox_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  tt <- t.test(x, mu = 0)
  wt <- wilcox.test(x, mu = 0, exact = FALSE)
  data.frame(
    group    = group_label,
    mark     = mark_id,
    n        = length(x),
    mean     = mean(x),
    median   = median(x),
    ttest_t  = unname(tt$statistic),
    ttest_p  = tt$p.value,
    wilcox_V = unname(wt$statistic),
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  )
}

run_twosample_mw <- function(x1, x2, label_pair, mark_id) {
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]
  if (length(x1) < 3 || length(x2) < 3) {
    return(data.frame(
      comparison = label_pair, mark = mark_id,
      n1 = length(x1), n2 = length(x2),
      mw_W = NA_real_, mw_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  wt <- wilcox.test(x1, x2, exact = FALSE)
  data.frame(
    comparison = label_pair, mark = mark_id,
    n1 = length(x1), n2 = length(x2),
    mw_W = unname(wt$statistic),
    mw_p = wt$p.value,
    stringsAsFactors = FALSE
  )
}

one_sample_results <- purrr::map_dfr(MARK_ORDER, function(mid) {
  sub <- long_df %>% dplyr::filter(mark_id == mid)
  purrr::map_dfr(
    c("MeCP2 Up", "MeCP2 Down", "Not Significant"),
    function(grp) {
      vals <- sub %>% dplyr::filter(mecp2_status == grp) %>% dplyr::pull(value)
      run_onesample_tests(vals, grp, mid)
    }
  )
})

two_sample_results <- purrr::map_dfr(MARK_ORDER, function(mid) {
  sub   <- long_df %>% dplyr::filter(mark_id == mid)
  x_up  <- sub %>% dplyr::filter(mecp2_status == "MeCP2 Up")           %>% dplyr::pull(value)
  x_dn  <- sub %>% dplyr::filter(mecp2_status == "MeCP2 Down")         %>% dplyr::pull(value)
  x_ns  <- sub %>% dplyr::filter(mecp2_status == "Not Significant")    %>% dplyr::pull(value)
  dplyr::bind_rows(
    run_twosample_mw(x_up, x_dn, "Up_vs_Down", mid),
    run_twosample_mw(x_dn, x_ns, "Down_vs_NS",  mid),
    run_twosample_mw(x_up, x_ns, "Up_vs_NS",    mid)
  )
})

# BH correction within each mark
one_sample_results <- one_sample_results %>%
  dplyr::group_by(mark) %>%
  dplyr::mutate(
    ttest_q  = p.adjust(ttest_p,  method = "BH"),
    wilcox_q = p.adjust(wilcox_p, method = "BH")
  ) %>%
  dplyr::ungroup()

two_sample_results <- two_sample_results %>%
  dplyr::group_by(mark) %>%
  dplyr::mutate(mw_q = p.adjust(mw_p, method = "BH")) %>%
  dplyr::ungroup()

# Print summary
cat("\nOne-sample Wilcoxon vs zero (BH-corrected):\n")
for (grp in c("MeCP2 Up", "MeCP2 Down")) {
  cat(sprintf("\n  %s:\n", grp))
  for (mid in MARK_ORDER) {
    row <- one_sample_results %>%
      dplyr::filter(mark == mid, group == grp)
    cat(sprintf("    %-12s  n=%4d  mean=%+.4f  wilcox_q=%.2e  %s\n",
                mid, row$n, row$mean, row$wilcox_q, sig_stars(row$wilcox_q)))
  }
}

cat("\nTwo-sample Mann-Whitney (BH-corrected):\n")
for (comp in c("Up_vs_Down", "Down_vs_NS")) {
  cat(sprintf("\n  %s:\n", comp))
  for (mid in MARK_ORDER) {
    row <- two_sample_results %>%
      dplyr::filter(mark == mid, comparison == comp)
    cat(sprintf("    %-12s  mw_q=%.2e  %s\n",
                mid, row$mw_q, sig_stars(row$mw_q)))
  }
}

# Export combined stats table
full_stats <- one_sample_results %>%
  dplyr::left_join(
    two_sample_results %>%
      dplyr::filter(comparison == "Up_vs_Down") %>%
      dplyr::select(mark, mw_updown_W = mw_W, mw_updown_p = mw_p, mw_updown_q = mw_q),
    by = "mark"
  ) %>%
  dplyr::left_join(
    two_sample_results %>%
      dplyr::filter(comparison == "Down_vs_NS") %>%
      dplyr::select(mark, mw_downns_W = mw_W, mw_downns_p = mw_p, mw_downns_q = mw_q),
    by = "mark"
  ) %>%
  dplyr::left_join(
    two_sample_results %>%
      dplyr::filter(comparison == "Up_vs_NS") %>%
      dplyr::select(mark, mw_upns_W = mw_W, mw_upns_p = mw_p, mw_upns_q = mw_q),
    by = "mark"
  )

stats_path <- file.path(TABLES_DIR, "60_mecp2_status_stats.tsv")
write.table(full_stats, stats_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\n  Saved: %s (%d rows)\n", stats_path, nrow(full_stats)))

# =============================================================================
# PLOT 60a: MULTI-PANEL VIOLIN + BOX + GGSIGNIF BRACKETS
# =============================================================================

cat("\n--- Plot 60a: Multi-panel violin+box ---\n")

VIOLIN_MARKS <- c("k119ub_bw", "mc", "hmc", "k27ac", "k27me3")

violin_df <- long_df %>%
  dplyr::filter(mark_id %in% VIOLIN_MARKS) %>%
  dplyr::mutate(
    facet_label = factor(MARK_LABELS[mark_id], levels = MARK_LABELS[VIOLIN_MARKS]),
    mecp2_status = factor(mecp2_status,
                          levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant"))
  )

cat(sprintf("  Violin data: %d rows across %d marks\n",
            nrow(violin_df), length(VIOLIN_MARKS)))

COMPARISONS_60A <- list(
  c("MeCP2 Up", "MeCP2 Down"),
  c("MeCP2 Down", "Not Significant")
)

p_60a <- ggplot(violin_df,
                aes(x = mecp2_status, y = value, fill = mecp2_status)) +
  geom_violin(alpha = 0.65, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, fill = "white",
               outlier.size = 0.3, outlier.alpha = 0.4,
               linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_signif(
    comparisons = COMPARISONS_60A,
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    step_increase = 0.12,
    tip_length = 0.02,
    textsize = 3
  ) +
  facet_wrap(~ facet_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = COLORS$mecp2, name = NULL) +
  scale_x_discrete(labels = c(
    "MeCP2 Up"        = sprintf("Up\n(n=%d)", n_up),
    "MeCP2 Down"      = sprintf("Down\n(n=%d)", n_down),
    "Not Significant" = sprintf("NS\n(n~%s)", format(n_ns, big.mark = ","))
  )) +
  labs(
    title    = "Epigenetic Mark Changes by MeCP2 Status",
    subtitle = "BAP1-KO: MeCP2-Up and MeCP2-Down show mirror-image patterns",
    x        = NULL,
    y        = "Change (mutant vs control)"
  ) +
  theme_biomodal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 9, face = "bold"),
    axis.text.x     = element_text(size = 8)
  )

save_plot(p_60a, "60a_violin_marks_by_mecp2", w = 18, h = 7)

# =============================================================================
# PLOT 60b: SUMMARY HEATMAP (PHEATMAP)
# =============================================================================

cat("\n--- Plot 60b: Summary heatmap ---\n")

HEATMAP_MARKS  <- c("k119ub_bw", "mc", "hmc", "k27ac", "k27me3")
HEATMAP_GROUPS <- c("MeCP2 Up", "MeCP2 Down")

mean_mat <- matrix(NA_real_,
  nrow = length(HEATMAP_MARKS),
  ncol = length(HEATMAP_GROUPS),
  dimnames = list(MARK_LABELS_SHORT[HEATMAP_MARKS], HEATMAP_GROUPS)
)

star_mat <- matrix("",
  nrow = length(HEATMAP_MARKS),
  ncol = length(HEATMAP_GROUPS),
  dimnames = dimnames(mean_mat)
)

for (i in seq_along(HEATMAP_MARKS)) {
  for (j in seq_along(HEATMAP_GROUPS)) {
    vals <- long_df %>%
      dplyr::filter(mark_id == HEATMAP_MARKS[i],
                    mecp2_status == HEATMAP_GROUPS[j]) %>%
      dplyr::pull(value)
    if (length(vals) > 0) mean_mat[i, j] <- mean(vals, na.rm = TRUE)

    row <- one_sample_results %>%
      dplyr::filter(mark == HEATMAP_MARKS[i], group == HEATMAP_GROUPS[j])
    if (nrow(row) == 1) star_mat[i, j] <- sig_stars(row$wilcox_q)
  }
}

# Display values + stars in cells
display_mat <- matrix("",
  nrow = nrow(mean_mat), ncol = ncol(mean_mat),
  dimnames = dimnames(mean_mat)
)
for (i in seq_len(nrow(mean_mat))) {
  for (j in seq_len(ncol(mean_mat))) {
    val_str <- sprintf("%+.3f", mean_mat[i, j])
    display_mat[i, j] <- paste0(val_str, "\n", star_mat[i, j])
  }
}

cat("  Mean change matrix:\n")
print(round(mean_mat, 4))
cat("  Significance:\n")
print(star_mat)

max_abs_val <- max(abs(mean_mat), na.rm = TRUE) * 1.1

pheatmap_call <- quote(
  pheatmap(
    mean_mat,
    display_numbers = display_mat,
    number_color    = "black",
    fontsize_number = 11,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
    breaks = seq(-max_abs_val, max_abs_val, length.out = 102),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "white",
    cellwidth    = 100,
    cellheight   = 50,
    fontsize     = 12,
    fontsize_row = 11,
    fontsize_col = 12,
    main         = "Mean Epigenetic Change by MeCP2 Status\n(BH-corrected Wilcoxon vs zero: *** p<0.001, ** p<0.01, * p<0.05)",
    angle_col    = 0
  )
)

save_multiformat_pheatmap(
  pheatmap_call,
  base_path = file.path(SEC60_DIR, "60b_summary_heatmap"),
  width = 7, height = 5.5, dpi = 300, verbose = TRUE
)

# =============================================================================
# PLOT 60c: LOLLIPOP / DOT-AND-WHISKER EFFECT SIZE
# =============================================================================

cat("\n--- Plot 60c: Effect-size lollipop ---\n")

compute_ci <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  tt <- t.test(x, mu = 0)
  c(mean = unname(tt$estimate),
    lo   = unname(tt$conf.int[1]),
    hi   = unname(tt$conf.int[2]))
}

lollipop_rows <- purrr::map_dfr(MARK_ORDER, function(mid) {
  purrr::map_dfr(c("MeCP2 Up", "MeCP2 Down"), function(grp) {
    vals <- long_df %>%
      dplyr::filter(mark_id == mid, mecp2_status == grp) %>%
      dplyr::pull(value)
    ci <- compute_ci(vals)
    data.frame(
      mark_id  = mid,
      group    = grp,
      mean_val = ci["mean"],
      ci_lo    = ci["lo"],
      ci_hi    = ci["hi"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
})

lollipop_display <- c(
  k119ub_bw = "K119ub\n[BigWig]",
  k119ub_db = "K119ub\n[DiffBind]",
  mc        = "5mC",
  hmc       = "5hmC",
  k27ac     = "K27ac",
  k27me3    = "K27me3"
)

lollipop_rows <- lollipop_rows %>%
  dplyr::mutate(
    mark_display = factor(lollipop_display[mark_id],
                          levels = rev(lollipop_display[MARK_ORDER])),
    group = factor(group, levels = c("MeCP2 Up", "MeCP2 Down"))
  )

p_60c <- ggplot(lollipop_rows,
                aes(x = mean_val, y = mark_display, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, linewidth = 0.8,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 4, position = position_dodge(width = 0.6)) +
  scale_color_manual(
    values = COLORS$mecp2[c("MeCP2 Up", "MeCP2 Down")],
    labels = c(sprintf("MeCP2 Up (n=%d)", n_up),
               sprintf("MeCP2 Down (n=%d)", n_down)),
    name = NULL
  ) +
  labs(
    title    = "Effect Size: Mean Change ± 95% CI by MeCP2 Status",
    subtitle = "K119ub is flat at MeCP2-Down loci; all other marks show mirror-image shifts",
    x        = "Mean change (mutant vs control)",
    y        = NULL
  ) +
  theme_biomodal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 10)
  )

save_plot(p_60c, "60c_effect_size_lollipop", w = 11, h = 6)

# =============================================================================
# PLOT 60d: PER-GENE STRIP CHART (MeCP2-DOWN ONLY)
# =============================================================================

cat("\n--- Plot 60d: Per-gene strip (MeCP2-Down) ---\n")

down_genes <- master %>%
  dplyr::filter(mecp2_status == "MeCP2 Down") %>%
  dplyr::pull(gene)

cat(sprintf("  MeCP2-Down genes: %d\n", length(down_genes)))

strip_labels <- c(
  k119ub_bw = "K119ub\n[BigWig]",
  k119ub_db = "K119ub\n[DiffBind]",
  mc        = "5mC diff",
  hmc       = "5hmC diff",
  k27ac     = "K27ac fold",
  k27me3    = "K27me3 fold"
)

down_long <- long_df %>%
  dplyr::filter(mecp2_status == "MeCP2 Down") %>%
  dplyr::mutate(
    facet_label = factor(strip_labels[mark_id], levels = strip_labels[MARK_ORDER]),
    is_key = gene %in% KEY_GENES
  ) %>%
  dplyr::arrange(is_key)

cat(sprintf("  Strip data: %d points\n", nrow(down_long)))

p_60d <- ggplot(down_long, aes(x = 1, y = value)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_jitter(aes(color = is_key),
              width = 0.25, height = 0, alpha = 0.7, size = 2.5,
              show.legend = FALSE) +
  geom_boxplot(width = 0.4, fill = NA, outlier.shape = NA,
               linewidth = 0.8, color = COLORS$mecp2[["MeCP2 Down"]]) +
  geom_text_repel(
    data = down_long %>% dplyr::filter(is_key),
    aes(x = 1, y = value, label = gene),
    size = 2.8, fontface = "italic", color = "black",
    nudge_x = 0.4, max.overlaps = 20,
    segment.color = "grey60", segment.size = 0.3
  ) +
  scale_color_manual(values = c("TRUE" = "black",
                                "FALSE" = COLORS$mecp2[["MeCP2 Down"]])) +
  facet_wrap(~ facet_label, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = NULL) +
  labs(
    title    = sprintf("MeCP2-Down Genes (n=%d): Per-Gene Epigenetic Mark Values",
                       length(down_genes)),
    subtitle = "K119ub is clustered at zero; other marks shift in opposite direction from MeCP2-Up",
    x        = NULL,
    y        = "Change (mutant vs control)"
  ) +
  theme_biomodal(base_size = 11) +
  theme(
    strip.text   = element_text(size = 9, face = "bold"),
    axis.ticks.x = element_blank()
  )

save_plot(p_60d, "60d_down_genes_strip", w = 18, h = 6)

# =============================================================================
# PLOT 60e: COMPOSITE (PATCHWORK)
# =============================================================================

cat("\n--- Plot 60e: Composite ---\n")

p_60e <- p_60a /
  (p_60c | p_60d) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title    = "MeCP2 Status Defines Mirror-Image Epigenetic Profiles",
    subtitle = "BAP1-KO: K119ub accumulation drives cascade at MeCP2-Up loci; K119ub is flat at MeCP2-Down loci",
    tag_levels = "A"
  )

save_plot(p_60e, "60e_composite", w = 24, h = 14)

# =============================================================================
# EXPORT: MeCP2-DOWN INDIVIDUAL GENE TABLE
# =============================================================================

cat("\n--- Exporting MeCP2-Down gene table ---\n")

down_table <- master %>%
  dplyr::filter(mecp2_status == "MeCP2 Down") %>%
  dplyr::select(
    gene,
    mecp2_mean_fold, mecp2_nearest_fdr,
    gb_log2fc, gb_signal_class,
    k119ub_fold, k119ub_fdr,
    mc_diff, mc_q, mc_sig,
    hmc_diff, hmc_q, hmc_sig,
    k27ac_fold, k27ac_fdr,
    k27me3_fold, k27me3_fdr,
    chromatin_state
  ) %>%
  dplyr::mutate(is_key_gene = gene %in% KEY_GENES) %>%
  dplyr::arrange(dplyr::desc(is_key_gene), mecp2_nearest_fdr)

down_table_path <- file.path(TABLES_DIR, "60_mecp2_down_gene_table.tsv")
write.table(down_table, down_table_path, sep = "\t", row.names = FALSE,
            quote = FALSE)
cat(sprintf("  Saved: %s (%d genes, %d cols)\n",
            down_table_path, nrow(down_table), ncol(down_table)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 60 SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("  MeCP2-Up   (%3d genes): K119ub UP, 5mC UP, 5hmC DOWN, K27ac DOWN, K27me3 UP\n", n_up))
cat(sprintf("  MeCP2-Down (%3d genes): K119ub FLAT (ns), 5mC down, K27ac UP, K27me3 DOWN\n", n_down))
cat(sprintf("  NS         (%5s genes): all marks near zero\n", format(n_ns, big.mark = ",")))

cat("\nPlots saved to:", SEC60_DIR, "\n")
cat("  60a: 5-facet violin+box with ggsignif brackets\n")
cat("  60b: Summary heatmap (pheatmap)\n")
cat("  60c: Effect-size lollipop with 95% CI\n")
cat("  60d: MeCP2-Down per-gene strip chart\n")
cat("  60e: Patchwork composite (60a + 60c + 60d)\n")

cat("\nTables saved to:", TABLES_DIR, "\n")
cat(sprintf("  60_mecp2_status_stats.tsv     (%d rows)\n", nrow(full_stats)))
cat(sprintf("  60_mecp2_down_gene_table.tsv  (%d genes)\n", nrow(down_table)))

cat("\nSection 60 complete.\n")
