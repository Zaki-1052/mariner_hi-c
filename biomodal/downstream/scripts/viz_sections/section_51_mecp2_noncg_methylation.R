# biomodal/downstream/scripts/viz_sections/section_51_mecp2_noncg_methylation.R
# Section 51: Non-CG Methylation at MeCP2 Binding Sites
#
# Question: At MeCP2 peaks — especially those in regions with low CG methylation
# (where MeCP2 might be reading mCH instead of mCG) — is there measurable non-CG
# methylation, and does it differ between BAP1-KO mutant and wildtype?
#
# 51c: Stratify MeCP2-bound genes by CG level, measure CHG/CHH within each stratum
# 51d: Compare non-CG methylation at MeCP2-bound vs non-MeCP2-bound genes
#
# Note: evoC cannot distinguish mC from hmC at non-CG contexts. Only mC is analyzed.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_51_mecp2_noncg_methylation.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 51: NON-CG METHYLATION AT MeCP2 BINDING SITES
# =============================================================================

cat("================================================================================\n")
cat("SECTION 51: NON-CG METHYLATION AT MeCP2 BINDING SITES\n")
cat("================================================================================\n\n")

SEC51_DIR <- file.path(OUTPUT_DIR, "51_mecp2_noncg")
dir.create(SEC51_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Helper ------------------------------------------------------------------

load_extract <- function(filepath, context_label) {
  stopifnot(file.exists(filepath))
  df <- read.table(gzfile(filepath), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  if (!grepl("^chr", df$Chromosome[1])) {
    df$Chromosome <- paste0("chr", df$Chromosome)
  }
  cat(sprintf("  %s extract: %d genes\n", context_label, nrow(df)))
  return(df)
}

# ---- Load MeCP2 peaks -------------------------------------------------------

cat("Loading MeCP2 peaks...\n")
mecp2_up_gr <- load_chip_peaks(MECP2_FILES$up, "MeCP2 Up")
mecp2_down_gr <- load_chip_peaks(MECP2_FILES$down, "MeCP2 Down")
mecp2_all_gr <- c(mecp2_up_gr, mecp2_down_gr)

# ---- Load per-sample extracts (CG + CHG + CHH) ------------------------------

cat("\nLoading per-sample extracts...\n")
cg_extract <- load_extract(EXTRACT_PATHS$mc_regional_frac, "CG")
chg_extract <- load_extract(CHG_EXTRACT_PATHS$mc_regional_frac, "CHG")
chh_extract <- load_extract(CHH_EXTRACT_PATHS$mc_regional_frac, "CHH")

CTRL_COLS <- grep("ctrl", colnames(chg_extract), value = TRUE)
MUT_COLS <- grep("mut", colnames(chg_extract), value = TRUE)
CG_CTRL_COLS <- grep("ctrl", colnames(cg_extract), value = TRUE)
CG_MUT_COLS <- grep("mut", colnames(cg_extract), value = TRUE)

cat(sprintf("  %d ctrl, %d mut samples\n", length(CTRL_COLS), length(MUT_COLS)))

# ---- Identify MeCP2-overlapping genes ---------------------------------------

cat("\nFinding genes overlapping MeCP2 peaks...\n")
gene_gr <- GRanges(
  seqnames = chg_extract$Chromosome,
  ranges = IRanges(start = chg_extract$Start, end = chg_extract$End)
)
mecp2_overlap <- countOverlaps(gene_gr, mecp2_all_gr) > 0
mecp2_up_overlap <- countOverlaps(gene_gr, mecp2_up_gr) > 0

cat(sprintf("  %d / %d genes overlap any MeCP2 peak\n",
            sum(mecp2_overlap), length(mecp2_overlap)))
cat(sprintf("  %d / %d genes overlap MeCP2-Up peaks\n",
            sum(mecp2_up_overlap), length(mecp2_up_overlap)))

# ---- Compute per-gene CG methylation (ctrl baseline) ------------------------

cat("\nComputing per-gene CG methylation for stratification...\n")
cg_ctrl_matrix <- as.matrix(cg_extract[, CG_CTRL_COLS])
cg_gene_mean <- rowMeans(cg_ctrl_matrix, na.rm = TRUE)

# Median split among MeCP2-overlapping genes
mecp2_cg_vals <- cg_gene_mean[mecp2_overlap]
cg_median <- median(mecp2_cg_vals, na.rm = TRUE)
cat(sprintf("  CG median at MeCP2-bound genes: %.4f (%.1f%%)\n", cg_median, cg_median * 100))

# =============================================================================
# FIGURE 51c: Non-CG methylation at MeCP2 peaks, stratified by CG level
# =============================================================================

cat("\n--- Figure 51c: Non-CG methylation at MeCP2 peaks, stratified by CG ---\n")

build_stratified_noncg <- function(extract_df, mecp2_mask, cg_gene_mean, cg_median,
                                   ctrl_cols, mut_cols, context_label) {
  sub <- extract_df[mecp2_mask, ]
  sub_cg <- cg_gene_mean[mecp2_mask]

  sub$cg_stratum <- ifelse(sub_cg >= cg_median, "CG-High", "CG-Low")

  long <- sub %>%
    pivot_longer(cols = all_of(c(ctrl_cols, mut_cols)),
                 names_to = "sample", values_to = "methylation") %>%
    mutate(
      condition = ifelse(grepl("ctrl", sample), "Control", "Mutant"),
      context = context_label
    )

  sample_means <- long %>%
    group_by(sample, condition, cg_stratum, context) %>%
    summarise(mean_meth = mean(methylation, na.rm = TRUE), .groups = "drop")

  return(sample_means)
}

strat_chg <- build_stratified_noncg(chg_extract, mecp2_overlap, cg_gene_mean,
                                    cg_median, CTRL_COLS, MUT_COLS, "CHG")
strat_chh <- build_stratified_noncg(chh_extract, mecp2_overlap, cg_gene_mean,
                                    cg_median, CTRL_COLS, MUT_COLS, "CHH")

strat_all <- bind_rows(strat_chg, strat_chh) %>%
  mutate(
    context = factor(context, levels = c("CHG", "CHH")),
    cg_stratum = factor(cg_stratum, levels = c("CG-Low", "CG-High"))
  )

# Statistical tests per stratum per context
cat("\n  Wilcoxon tests (ctrl vs mut) per stratum:\n")
test_results <- list()
for (ctx in c("CHG", "CHH")) {
  for (stratum in c("CG-Low", "CG-High")) {
    sub <- strat_all %>%
      dplyr::filter(context == ctx, cg_stratum == stratum)
    ctrl_vals <- sub$mean_meth[sub$condition == "Control"]
    mut_vals <- sub$mean_meth[sub$condition == "Mutant"]
    wt <- wilcox.test(ctrl_vals, mut_vals)
    key <- paste(ctx, stratum)
    test_results[[key]] <- wt
    cat(sprintf("    %s %s: p=%.3e (ctrl=%.2e, mut=%.2e)\n",
                ctx, stratum, wt$p.value, mean(ctrl_vals), mean(mut_vals)))
  }
}

n_mecp2 <- sum(mecp2_overlap)
n_cg_low <- sum(cg_gene_mean[mecp2_overlap] < cg_median)
n_cg_high <- n_mecp2 - n_cg_low

p_51c <- ggplot(strat_all, aes(x = condition, y = mean_meth, fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  facet_grid(context ~ cg_stratum, scales = "free_y") +
  scale_fill_manual(values = COLORS$condition, guide = "none") +
  labs(
    title = "Non-CG Methylation at MeCP2 Binding Sites, Stratified by CG Level",
    subtitle = sprintf(
      "CG-Low (<%0.f%% CG mC, n=%d genes) vs CG-High (n=%d genes) | CG median = %.1f%%",
      cg_median * 100, n_cg_low, n_cg_high, cg_median * 100
    ),
    x = "Condition", y = "Mean Regional Fraction (mC)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_51c, file.path(SEC51_DIR, "51c_noncg_at_mecp2_by_cg_stratum"),
                        width = 10, height = 7)

# =============================================================================
# FIGURE 51d: Non-CG methylation at MeCP2-bound vs non-MeCP2-bound genes
# =============================================================================

cat("\n--- Figure 51d: MeCP2-bound vs genome-wide non-CG methylation ---\n")

build_mecp2_vs_bg <- function(extract_df, mecp2_mask, ctrl_cols, mut_cols,
                              context_label) {
  extract_df$mecp2_group <- ifelse(mecp2_mask, "MeCP2-Bound", "Non-MeCP2")

  long <- extract_df %>%
    pivot_longer(cols = all_of(c(ctrl_cols, mut_cols)),
                 names_to = "sample", values_to = "methylation") %>%
    mutate(
      condition = ifelse(grepl("ctrl", sample), "Control", "Mutant"),
      context = context_label
    )

  sample_means <- long %>%
    group_by(sample, condition, mecp2_group, context) %>%
    summarise(mean_meth = mean(methylation, na.rm = TRUE), .groups = "drop")

  return(sample_means)
}

bg_chg <- build_mecp2_vs_bg(chg_extract, mecp2_overlap, CTRL_COLS, MUT_COLS, "CHG")
bg_chh <- build_mecp2_vs_bg(chh_extract, mecp2_overlap, CTRL_COLS, MUT_COLS, "CHH")

bg_all <- bind_rows(bg_chg, bg_chh) %>%
  mutate(
    context = factor(context, levels = c("CHG", "CHH")),
    mecp2_group = factor(mecp2_group, levels = c("MeCP2-Bound", "Non-MeCP2"))
  )

cat("\n  Wilcoxon tests (MeCP2-bound vs non-MeCP2, within ctrl):\n")
for (ctx in c("CHG", "CHH")) {
  sub <- bg_all %>% dplyr::filter(context == ctx, condition == "Control")
  bound_vals <- sub$mean_meth[sub$mecp2_group == "MeCP2-Bound"]
  bg_vals <- sub$mean_meth[sub$mecp2_group == "Non-MeCP2"]
  wt <- wilcox.test(bound_vals, bg_vals)
  cat(sprintf("    %s ctrl: MeCP2-bound=%.2e vs non-MeCP2=%.2e, p=%.3e\n",
              ctx, mean(bound_vals), mean(bg_vals), wt$p.value))
}

cat("\n  Wilcoxon tests (ctrl vs mut at MeCP2-bound genes):\n")
bg_tests <- list()
for (ctx in c("CHG", "CHH")) {
  sub <- bg_all %>% dplyr::filter(context == ctx, mecp2_group == "MeCP2-Bound")
  ctrl_vals <- sub$mean_meth[sub$condition == "Control"]
  mut_vals <- sub$mean_meth[sub$condition == "Mutant"]
  wt <- wilcox.test(ctrl_vals, mut_vals)
  bg_tests[[ctx]] <- wt
  cat(sprintf("    %s MeCP2-bound: ctrl=%.2e vs mut=%.2e, p=%.3e\n",
              ctx, mean(ctrl_vals), mean(mut_vals), wt$p.value))
}

p_51d <- ggplot(bg_all, aes(x = condition, y = mean_meth, fill = condition)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  facet_grid(context ~ mecp2_group, scales = "free_y") +
  scale_fill_manual(values = COLORS$condition, guide = "none") +
  labs(
    title = "Non-CG Methylation: MeCP2-Bound vs Non-MeCP2 Genes",
    subtitle = sprintf(
      "MeCP2-bound: %d genes | Non-MeCP2: %d genes | Is non-CG enriched at MeCP2 sites?",
      sum(mecp2_overlap), sum(!mecp2_overlap)
    ),
    x = "Condition", y = "Mean Regional Fraction (mC)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_51d, file.path(SEC51_DIR, "51d_mecp2_bound_vs_background_noncg"),
                        width = 10, height = 7)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n--- Writing summary table ---\n")

summary_rows <- list()

for (ctx in c("CHG", "CHH")) {
  for (stratum in c("CG-Low", "CG-High")) {
    key <- paste(ctx, stratum)
    wt <- test_results[[key]]
    sub <- strat_all %>% dplyr::filter(context == ctx, cg_stratum == stratum)
    ctrl_mean <- mean(sub$mean_meth[sub$condition == "Control"])
    mut_mean <- mean(sub$mean_meth[sub$condition == "Mutant"])
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      test = sprintf("Wilcoxon ctrl vs mut: %s at MeCP2 %s genes", ctx, stratum),
      p_value = wt$p.value,
      ctrl_mean = ctrl_mean,
      mut_mean = mut_mean,
      effect = mut_mean - ctrl_mean,
      n_genes = if (stratum == "CG-Low") n_cg_low else n_cg_high,
      stringsAsFactors = FALSE
    )
  }
}

for (ctx in c("CHG", "CHH")) {
  wt <- bg_tests[[ctx]]
  sub <- bg_all %>% dplyr::filter(context == ctx, mecp2_group == "MeCP2-Bound")
  ctrl_mean <- mean(sub$mean_meth[sub$condition == "Control"])
  mut_mean <- mean(sub$mean_meth[sub$condition == "Mutant"])
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    test = sprintf("Wilcoxon ctrl vs mut: %s at all MeCP2-bound genes", ctx),
    p_value = wt$p.value,
    ctrl_mean = ctrl_mean,
    mut_mean = mut_mean,
    effect = mut_mean - ctrl_mean,
    n_genes = sum(mecp2_overlap),
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
write.table(summary_df, file.path(TABLES_DIR, "mecp2_noncg_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: mecp2_noncg_summary.tsv\n")

cat("\n================================================================================\n")
cat("SECTION 51 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Figures saved to: %s\n", SEC51_DIR))
