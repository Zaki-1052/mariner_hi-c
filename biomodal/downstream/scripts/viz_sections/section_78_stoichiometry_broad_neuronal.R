# biomodal/downstream/scripts/viz_sections/section_78_stoichiometry_broad_neuronal.R
# Section 78: Stoichiometry & Mechanism with Unbiased Neuronal Gene Set
#
# Consolidates the key analyses from section 61 (and add-ons 61h-k) using the
# broader GO-derived neuronal gene set (section 72: synap|neuron|axon|dendrit|nervous,
# ~5,614 genes) and synapse/axon subset (section 76: synap|axon, ~2,321 genes)
# instead of section 61's DMR-enriched set (~1,149 genes from GO enrichment of
# significant DMRs — circular: genes selected for methylation changes test positive
# for methylation changes).
#
# Panels:
#   78a: Total methylation conservation — delta_total per gene group (forest plot)
#   78b: Stoichiometry scatter — mc_diff vs hmc_diff colored by neuronal status
#   78c: Stoichiometry slope forest — slopes per gene group vs reference -1
#   78d: BAP1-KO vs TET-KO phenocopy — delta_ratio correlation by gene group
#   78e: Narrow vs broad neuronal comparison (bias demonstration)
#   78g: MeCP2-Up + K119ub-Up quadrant characterization (replaces 61i)
#   78h: Composite
#
# Input:
#   tables/diffbind_gene_level_all_marks.tsv, tables/59_quadrant_master.tsv,
#   tables/72_neuronal_gene_set_go_derived.tsv, tables/76_synapse_axon_gene_set.tsv,
#   tables/61_neuronal_gene_set.tsv, tables/mecp2_coordinated_genes.tsv,
#   data/tet_ko_gene_signal.tsv
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_78_stoichiometry_broad_neuronal.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(patchwork)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 78: STOICHIOMETRY & MECHANISM — UNBIASED NEURONAL SET\n")
cat("  Redoing section 61 with GO-derived neuronal genes (not DMR-enriched)\n")
cat("================================================================================\n\n")

SEC78_DIR <- file.path(OUTPUT_DIR, "78_stoichiometry_broad_neuronal")
dir.create(SEC78_DIR, recursive = TRUE, showWarnings = FALSE)

DIFFBIND_PATH    <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
QUADRANT_PATH    <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
NEURONAL_BROAD   <- file.path(TABLES_DIR, "72_neuronal_gene_set_go_derived.tsv")
NEURONAL_NARROW  <- file.path(TABLES_DIR, "61_neuronal_gene_set.tsv")
SYNAPSE_PATH     <- file.path(TABLES_DIR, "76_synapse_axon_gene_set.tsv")
COORDINATED_PATH <- file.path(TABLES_DIR, "mecp2_coordinated_genes.tsv")
TETKO_PATH       <- file.path(BASE_DIR, "data/tet_ko_gene_signal.tsv")

for (f in c(DIFFBIND_PATH, QUADRANT_PATH, NEURONAL_BROAD, COORDINATED_PATH, TETKO_PATH)) {
  if (!file.exists(f)) stop(paste0("Required file not found: ", f))
}

# Synapse set may not exist if section 76 hasn't been run
has_synapse <- file.exists(SYNAPSE_PATH)
has_narrow  <- file.exists(NEURONAL_NARROW)

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

sig_stars <- function(q) {
  if (is.na(q) || !is.finite(q)) return("")
  if (q < 0.001) return("***")
  if (q < 0.01)  return("**")
  if (q < 0.05)  return("*")
  return("ns")
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC78_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n\n")

db <- read.table(DIFFBIND_PATH, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind all-marks: %d genes\n", nrow(db)))

qm <- read.table(QUADRANT_PATH, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Quadrant master:    %d genes\n", nrow(qm)))

neuronal_broad <- read.table(NEURONAL_BROAD, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")$gene
cat(sprintf("  Neuronal (broad):   %d genes\n", length(neuronal_broad)))

if (has_narrow) {
  neuronal_narrow <- read.table(NEURONAL_NARROW, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE, quote = "")$gene
  cat(sprintf("  Neuronal (narrow):  %d genes\n", length(neuronal_narrow)))
} else {
  neuronal_narrow <- character(0)
  cat("  Neuronal (narrow):  not found (section 61 not run)\n")
}

if (has_synapse) {
  synapse_genes <- read.table(SYNAPSE_PATH, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, quote = "")$gene
  cat(sprintf("  Synapse/axon:       %d genes\n", length(synapse_genes)))
} else {
  synapse_genes <- character(0)
  cat("  Synapse/axon:       not found (section 76 not run)\n")
}

coord_genes <- read.table(COORDINATED_PATH, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "")$gene
cat(sprintf("  Coordinated:        %d genes\n", length(coord_genes)))

tetko <- read.table(TETKO_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  TET-KO signal:      %d genes\n", nrow(tetko)))

# --- Build master table ---

mecp2_up_genes <- qm$gene[!is.na(qm$mecp2_nearest_fdr) &
                           qm$mecp2_nearest_fdr < Q_THRESHOLD &
                           !is.na(qm$mecp2_mean_fold) &
                           qm$mecp2_mean_fold > 0]

df <- db
df$is_neuronal_broad  <- df$gene %in% neuronal_broad
df$is_neuronal_narrow <- df$gene %in% neuronal_narrow
df$is_synapse         <- df$gene %in% synapse_genes
df$is_coordinated     <- df$gene %in% coord_genes
df$is_mecp2_up        <- df$gene %in% mecp2_up_genes

df$delta_total <- df$total_mut - df$total_ctrl

cat(sprintf("\n  Working table: %d genes\n", nrow(df)))
cat(sprintf("  Neuronal (broad):  %d (%.1f%%)\n",
            sum(df$is_neuronal_broad), 100 * mean(df$is_neuronal_broad)))
cat(sprintf("  Neuronal (narrow): %d (%.1f%%)\n",
            sum(df$is_neuronal_narrow), 100 * mean(df$is_neuronal_narrow)))
cat(sprintf("  Synapse/axon:      %d (%.1f%%)\n",
            sum(df$is_synapse), 100 * mean(df$is_synapse)))
cat(sprintf("  Coordinated:       %d (%.1f%%)\n",
            sum(df$is_coordinated), 100 * mean(df$is_coordinated)))
cat(sprintf("  MeCP2-Up:          %d\n", sum(df$is_mecp2_up)))

# =============================================================================
# 78a: TOTAL METHYLATION CONSERVATION — DELTA_TOTAL PER GENE GROUP
# =============================================================================

cat("\n--- 78a: Total methylation conservation test ---\n")

gene_groups <- list(
  "All genes"         = rep(TRUE, nrow(df)),
  "Non-neuronal"      = !df$is_neuronal_broad,
  "Neuronal (broad)"  = df$is_neuronal_broad,
  "Synapse/axon"      = df$is_synapse,
  "Coordinated"       = df$is_coordinated,
  "Neur + Coord"      = df$is_neuronal_broad & df$is_coordinated,
  "MeCP2-Up"          = df$is_mecp2_up
)

# Remove empty groups
gene_groups <- gene_groups[sapply(gene_groups, function(m) sum(m & !is.na(df$delta_total)) >= 3)]

stoich_results <- data.frame()
for (grp_name in names(gene_groups)) {
  mask <- gene_groups[[grp_name]] & !is.na(df$delta_total)
  vals <- df$delta_total[mask]
  n <- length(vals)
  if (n < 3) next

  wt <- wilcox.test(vals, mu = 0, conf.int = TRUE)
  stoich_results <- rbind(stoich_results, data.frame(
    group = grp_name, n = n,
    mean_delta = mean(vals),
    median_delta = median(vals),
    ci_lo = wt$conf.int[1],
    ci_hi = wt$conf.int[2],
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  ))
}
stoich_results$wilcox_q <- p.adjust(stoich_results$wilcox_p, method = "BH")
stoich_results$sig <- sapply(stoich_results$wilcox_q, sig_stars)

cat("\n  Total methylation change (Wilcoxon vs 0, BH-corrected):\n")
for (i in seq_len(nrow(stoich_results))) {
  r <- stoich_results[i, ]
  direction <- ifelse(r$mean_delta > 0, "UP", ifelse(r$mean_delta < 0, "DOWN", "FLAT"))
  cat(sprintf("    %-20s n=%5d  mean=%+.5f  median=%+.5f  %s  q=%s %s\n",
              r$group, r$n, r$mean_delta, r$median_delta, direction,
              fmt_p(r$wilcox_q), r$sig))
}

write.table(stoich_results,
            file.path(TABLES_DIR, "78_total_methylation_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Forest plot
stoich_results$group <- factor(stoich_results$group, levels = rev(stoich_results$group))

p_78a <- ggplot(stoich_results, aes(x = group, y = mean_delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n), color = "#4393C3") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, color = "#4393C3") +
  geom_text(aes(label = sprintf("%s (n=%s)", sig, format(n, big.mark = ","))),
            hjust = -0.2, size = 3.2) +
  coord_flip() +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(
    title = "Total methylation (mC+hmC) change in BAP1-KO",
    subtitle = "Mean delta (mut - ctrl) with Wilcoxon 95% CI | dashed = conservation",
    x = NULL, y = "Mean delta total methylation"
  ) +
  theme_biomodal()

save_plot(p_78a, "78a_total_methylation_forest", w = 11, h = 7)

# =============================================================================
# 78b: STOICHIOMETRY SCATTER (mc_diff vs hmc_diff)
# =============================================================================

cat("\n--- 78b: Stoichiometry scatter ---\n")

valid <- !is.na(df$mc_diff) & !is.na(df$hmc_diff)

# Global slope
global_lm <- lm(mc_diff ~ hmc_diff, data = df[valid, ])
global_slope <- coef(global_lm)["hmc_diff"]
global_ci <- confint(global_lm, "hmc_diff", level = 0.95)
global_r2 <- summary(global_lm)$r.squared

cat(sprintf("  Global slope: %.4f [%.4f, %.4f]  R2=%.4f\n",
            global_slope, global_ci[1], global_ci[2], global_r2))

scatter_df <- df[valid, c("gene", "mc_diff", "hmc_diff",
                           "is_neuronal_broad", "is_synapse")]
scatter_df$highlight <- dplyr::case_when(
  scatter_df$is_synapse ~ "Synapse/axon",
  scatter_df$is_neuronal_broad ~ "Broader neuronal",
  TRUE ~ "Non-neuronal"
)
scatter_df$highlight <- factor(scatter_df$highlight,
  levels = c("Synapse/axon", "Broader neuronal", "Non-neuronal"))

p_78b <- ggplot(scatter_df, aes(x = hmc_diff, y = mc_diff)) +
  geom_point(aes(color = highlight), alpha = 0.3, size = 0.6) +
  geom_smooth(method = "lm", color = "black", linewidth = 1, se = TRUE, alpha = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "#E69F00", linewidth = 1) +
  scale_color_manual(values = c(
    "Synapse/axon" = "#D95F02",
    "Broader neuronal" = "#756BB1",
    "Non-neuronal" = "grey80"
  ), name = NULL) +
  labs(
    title = "Stoichiometry: delta-5mC vs delta-5hmC",
    subtitle = sprintf("Black = OLS (slope=%.3f, R2=%.4f) | Orange = reference -1 (dehydroxymethylase)",
                       global_slope, global_r2),
    x = "delta-5hmC (mut - ctrl)", y = "delta-5mC (mut - ctrl)"
  ) +
  coord_fixed(ratio = 1) +
  theme_biomodal()

save_plot(p_78b, "78b_stoichiometry_scatter", w = 10, h = 9)

# =============================================================================
# 78c: STOICHIOMETRY SLOPE FOREST PLOT
# =============================================================================

cat("\n--- 78c: Stoichiometry slopes by gene group ---\n")

slope_groups <- list(
  "All genes"         = rep(TRUE, nrow(df)),
  "Non-neuronal"      = !df$is_neuronal_broad,
  "Neuronal (broad)"  = df$is_neuronal_broad,
  "Synapse/axon"      = df$is_synapse,
  "Coordinated"       = df$is_coordinated,
  "Neur + Coord"      = df$is_neuronal_broad & df$is_coordinated,
  "MeCP2-Up"          = df$is_mecp2_up
)

slope_results <- data.frame()
for (grp_name in names(slope_groups)) {
  mask <- slope_groups[[grp_name]] & valid
  n <- sum(mask)
  if (n < 10) next

  m <- lm(mc_diff ~ hmc_diff, data = df[mask, ])
  ci <- confint(m, "hmc_diff", level = 0.95)
  s <- summary(m)
  differs <- ci[1] > -1 | ci[2] < -1

  slope_results <- rbind(slope_results, data.frame(
    group = grp_name, n = n,
    slope = coef(m)["hmc_diff"],
    ci_lo = ci[1], ci_hi = ci[2],
    r_squared = s$r.squared,
    differs_from_neg1 = differs,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("    %-20s slope=%.4f [%.4f, %.4f]  R2=%.4f  n=%d  %s\n",
              grp_name, coef(m)["hmc_diff"], ci[1], ci[2],
              s$r.squared, n,
              ifelse(differs, "!= -1", "~ -1")))
}

write.table(slope_results,
            file.path(TABLES_DIR, "78_stoichiometry_slopes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

slope_results$group <- factor(slope_results$group, levels = rev(slope_results$group))

p_78c <- ggplot(slope_results, aes(x = group, y = slope)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "#E69F00", linewidth = 0.8) +
  geom_point(aes(color = differs_from_neg1), size = 4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, color = differs_from_neg1),
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sprintf("n=%s", format(n, big.mark = ","))),
            hjust = -0.3, size = 3) +
  scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4393C3"),
                      labels = c("TRUE" = "Differs from -1", "FALSE" = "Consistent with -1"),
                      name = NULL) +
  coord_flip() +
  labs(
    title = "Stoichiometry slopes: delta-5mC ~ delta-5hmC",
    subtitle = "Orange dashed = -1 (stoichiometric dehydroxymethylase) | Red = significantly differs",
    x = NULL, y = "OLS slope (95% CI)"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_plot(p_78c, "78c_stoichiometry_slope_forest", w = 11, h = 7)

# =============================================================================
# 78d: BAP1-KO vs TET-KO PHENOCOPY
# =============================================================================

cat("\n--- 78d: TET-KO phenocopy comparison ---\n")

tetko_merged <- merge(
  df[, c("gene", "delta_ratio", "is_neuronal_broad", "is_synapse", "is_coordinated")],
  tetko[, c("gene", "delta_ratio_tet")],
  by = "gene"
)
tetko_valid <- tetko_merged[!is.na(tetko_merged$delta_ratio) &
                             !is.na(tetko_merged$delta_ratio_tet), ]

cat(sprintf("  Genes with both BAP1-KO and TET-KO data: %d\n", nrow(tetko_valid)))

tetko_groups <- list(
  "All genes"        = rep(TRUE, nrow(tetko_valid)),
  "Non-neuronal"     = !tetko_valid$is_neuronal_broad,
  "Neuronal (broad)" = tetko_valid$is_neuronal_broad,
  "Synapse/axon"     = tetko_valid$is_synapse,
  "Coordinated"      = tetko_valid$is_coordinated
)

tetko_cors <- data.frame()
for (grp in names(tetko_groups)) {
  mask <- tetko_groups[[grp]]
  n <- sum(mask)
  if (n < 5) next

  ct <- cor.test(tetko_valid$delta_ratio[mask],
                 tetko_valid$delta_ratio_tet[mask],
                 method = "spearman")
  tetko_cors <- rbind(tetko_cors, data.frame(
    group = grp, n = n,
    rho = as.numeric(ct$estimate),
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("    %-20s rho=%.4f  %s  n=%d\n",
              grp, ct$estimate, fmt_p(ct$p.value), n))
}

write.table(tetko_cors,
            file.path(TABLES_DIR, "78_tetko_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

tetko_valid$highlight <- dplyr::case_when(
  tetko_valid$is_synapse ~ "Synapse/axon",
  tetko_valid$is_neuronal_broad ~ "Broader neuronal",
  TRUE ~ "Non-neuronal"
)

global_rho <- cor.test(tetko_valid$delta_ratio, tetko_valid$delta_ratio_tet,
                       method = "spearman")

p_78d <- ggplot(tetko_valid, aes(x = delta_ratio_tet, y = delta_ratio)) +
  geom_point(aes(color = highlight), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(
    "Synapse/axon" = "#D95F02",
    "Broader neuronal" = "#756BB1",
    "Non-neuronal" = "grey70"
  ), name = NULL) +
  labs(
    title = "BAP1-KO vs TET-KO: TET efficiency change",
    subtitle = sprintf("Spearman rho=%.3f, %s, n=%s",
                       global_rho$estimate, fmt_p(global_rho$p.value),
                       format(nrow(tetko_valid), big.mark = ",")),
    x = "TET-KO delta ratio (hmC/(mC+hmC) change)",
    y = "BAP1-KO delta ratio (hmC/(mC+hmC) change)"
  ) +
  theme_biomodal()

save_plot(p_78d, "78d_bap1_vs_tetko", w = 9, h = 8)

# =============================================================================
# 78e: NARROW vs BROAD NEURONAL COMPARISON (BIAS DEMONSTRATION)
# =============================================================================

cat("\n--- 78e: Narrow vs broad neuronal comparison ---\n")

if (has_narrow && sum(df$is_neuronal_narrow) >= 10) {

  # Build comparison data
  narrow_vals <- df$delta_total[df$is_neuronal_narrow & !is.na(df$delta_total)]
  broad_vals  <- df$delta_total[df$is_neuronal_broad & !is.na(df$delta_total)]
  broad_only  <- df$delta_total[df$is_neuronal_broad & !df$is_neuronal_narrow & !is.na(df$delta_total)]

  cat(sprintf("  Narrow neuronal set (DMR-enriched): n=%d, mean delta=%.5f, median=%.5f\n",
              length(narrow_vals), mean(narrow_vals), median(narrow_vals)))
  cat(sprintf("  Broad neuronal set (GO-derived):    n=%d, mean delta=%.5f, median=%.5f\n",
              length(broad_vals), mean(broad_vals), median(broad_vals)))
  cat(sprintf("  Broad-only (in broad, not narrow):  n=%d, mean delta=%.5f, median=%.5f\n",
              length(broad_only), mean(broad_only), median(broad_only)))

  wt_narrow <- wilcox.test(narrow_vals, mu = 0)
  wt_broad  <- wilcox.test(broad_vals, mu = 0)
  wt_diff   <- wilcox.test(narrow_vals, broad_only)

  cat(sprintf("  Narrow vs 0: %s | Broad vs 0: %s\n",
              fmt_p(wt_narrow$p.value), fmt_p(wt_broad$p.value)))
  cat(sprintf("  Narrow vs broad-only: %s\n", fmt_p(wt_diff$p.value)))

  comp_long <- rbind(
    data.frame(set = "Narrow\n(DMR-enriched, n=1,149)",
               delta_total = narrow_vals, stringsAsFactors = FALSE),
    data.frame(set = "Broad-only\n(GO-derived minus narrow)",
               delta_total = broad_only, stringsAsFactors = FALSE),
    data.frame(set = "Full broad\n(GO-derived, n=5,614)",
               delta_total = broad_vals, stringsAsFactors = FALSE)
  )
  comp_long$set <- factor(comp_long$set, levels = unique(comp_long$set))

  p_78e <- ggplot(comp_long, aes(x = set, y = delta_total, fill = set)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = c("#D73027", "#4393C3", "#756BB1")) +
    annotate("text",
      x = c(1, 2, 3),
      y = rep(max(comp_long$delta_total, na.rm = TRUE) * 0.95, 3),
      label = c(
        sprintf("mean=%+.5f\n%s", mean(narrow_vals), fmt_p(wt_narrow$p.value)),
        sprintf("mean=%+.5f", mean(broad_only)),
        sprintf("mean=%+.5f\n%s", mean(broad_vals), fmt_p(wt_broad$p.value))
      ),
      size = 3, vjust = 1
    ) +
    labs(
      title = "Bias demonstration: narrow vs broad neuronal gene set",
      subtitle = paste0(
        "Narrow set was derived from DMR enrichment (circular selection bias)\n",
        "Broad set is GO-derived from org.Mm.eg.db (unbiased by methylation data)"
      ),
      x = NULL, y = "Delta total methylation (mut - ctrl)"
    ) +
    theme_biomodal()

  save_plot(p_78e, "78e_narrow_vs_broad_bias", w = 10, h = 8)

} else {
  cat("  Narrow neuronal set not available — skipping bias comparison\n")
  p_78e <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Narrow neuronal set not available\n(run section 61 first)",
             size = 6, hjust = 0.5) +
    theme_void() +
    labs(title = "Narrow vs broad neuronal comparison")
  save_plot(p_78e, "78e_narrow_vs_broad_bias", w = 10, h = 8)
}

# =============================================================================
# 78g: QUADRANT CHARACTERIZATION (replaces 61i with broad neuronal set)
# =============================================================================

cat("\n--- 78g: MeCP2-Up + K119ub-Up quadrant characterization ---\n")

QUAD_PATH <- file.path(TABLES_DIR, "61h_mecp2_up_k119ub_up_genes.tsv")

if (file.exists(QUAD_PATH)) {
  quad_genes <- read.table(QUAD_PATH, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "")
  cat(sprintf("  Quadrant genes (MeCP2-Up + K119ub-Up): %d\n", nrow(quad_genes)))

  quad_full <- merge(quad_genes, db[, c("gene", "atac_fold", "k27ac_fold",
                                         "k27me3_fold", "k119ub_fold",
                                         "total_ctrl", "total_mut", "delta_ratio",
                                         "dmr_status")],
                     by = "gene", all.x = TRUE)

  quad_full$is_neuronal_broad  <- quad_full$gene %in% neuronal_broad
  quad_full$is_neuronal_narrow <- quad_full$gene %in% neuronal_narrow
  quad_full$is_synapse         <- quad_full$gene %in% synapse_genes

  n_broad  <- sum(quad_full$is_neuronal_broad)
  n_narrow <- sum(quad_full$is_neuronal_narrow)
  n_syn    <- sum(quad_full$is_synapse)

  cat(sprintf("  Neuronal (broad):  %d / %d (%.0f%%)\n", n_broad, nrow(quad_full),
              100 * n_broad / nrow(quad_full)))
  cat(sprintf("  Neuronal (narrow): %d / %d (%.0f%%)\n", n_narrow, nrow(quad_full),
              100 * n_narrow / nrow(quad_full)))
  cat(sprintf("  Synapse/axon:      %d / %d (%.0f%%)\n", n_syn, nrow(quad_full),
              100 * n_syn / nrow(quad_full)))

  cat(sprintf("\n  Broad-neuronal quadrant genes: %s\n",
              paste(sort(quad_full$gene[quad_full$is_neuronal_broad]), collapse = ", ")))

  # Neuronal vs non-neuronal comparison across marks
  quad_full$gene_class <- dplyr::case_when(
    quad_full$is_synapse ~ "Synapse/axon",
    quad_full$is_neuronal_broad ~ "Broader neuronal",
    TRUE ~ "Non-neuronal"
  )

  quad_metrics <- c("mecp2_mean_fold", "gb_log2fc", "mc_diff", "hmc_diff",
                    "k27me3_fold", "k27ac_fold", "k119ub_fold", "atac_fold")
  quad_labels  <- c("MeCP2 fold", "K119ub log2FC", "5mC diff", "5hmC diff",
                    "K27me3 fold", "K27ac fold", "K119ub DiffBind", "ATAC fold")

  quad_comp <- data.frame()
  for (i in seq_along(quad_metrics)) {
    col <- quad_metrics[i]
    if (!col %in% colnames(quad_full)) next
    for (cls in c("Synapse/axon", "Broader neuronal", "Non-neuronal")) {
      vals <- quad_full[[col]][quad_full$gene_class == cls]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) next
      quad_comp <- rbind(quad_comp, data.frame(
        metric = quad_labels[i], gene_class = cls,
        n = length(vals), mean = mean(vals), median = median(vals),
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("\n  Per-class mark means in quadrant genes:\n")
  for (i in seq_len(nrow(quad_comp))) {
    r <- quad_comp[i, ]
    cat(sprintf("    %-18s %-18s n=%2d  mean=%+.4f\n",
                r$metric, r$gene_class, r$n, r$mean))
  }

  write.table(quad_comp,
              file.path(TABLES_DIR, "78_quadrant_neuronal_comparison.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Boxplot+jitter for key marks
  plot_metrics <- c("mecp2_mean_fold", "gb_log2fc", "mc_diff", "hmc_diff")
  plot_labels_q <- c("MeCP2 fold change", "K119ub log2FC", "5mC difference", "5hmC difference")

  violin_panels <- list()
  for (i in seq_along(plot_metrics)) {
    col <- plot_metrics[i]
    if (!col %in% colnames(quad_full)) next
    vdf <- quad_full[!is.na(quad_full[[col]]), c("gene", "gene_class", col)]
    colnames(vdf)[3] <- "value"

    p <- ggplot(vdf, aes(x = gene_class, y = value, fill = gene_class)) +
      geom_boxplot(alpha = 0.6, outlier.size = 1) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
      scale_fill_manual(values = c("Synapse/axon" = "#D95F02",
                                    "Broader neuronal" = "#756BB1",
                                    "Non-neuronal" = "grey70")) +
      labs(title = plot_labels_q[i], x = NULL, y = NULL) +
      theme_biomodal(base_size = 10) +
      theme(legend.position = "none")

    violin_panels[[length(violin_panels) + 1]] <- p
  }

  p_78g <- wrap_plots(violin_panels, ncol = 2) +
    plot_annotation(
      title = sprintf("MeCP2-Up + K119ub-Up quadrant genes (n=%d) — broad neuronal set",
                       nrow(quad_full)),
      subtitle = sprintf("Synapse/axon: %d | Broader neuronal: %d | Non-neuronal: %d",
                          n_syn, n_broad - n_syn, nrow(quad_full) - n_broad)
    )

  save_plot(p_78g, "78g_quadrant_characterization", w = 10, h = 9)

} else {
  cat("  61h quadrant gene file not found — skipping 78g\n")
  p_78g <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Run section 61h first", size = 6) +
    theme_void() +
    labs(title = "Quadrant characterization")
  save_plot(p_78g, "78g_quadrant_characterization", w = 10, h = 9)
}

# =============================================================================
# 78h: COMPOSITE
# =============================================================================

cat("\n--- 78h: Composite ---\n")

p_78h <- (p_78a | p_78e) / (p_78b | p_78c) / (p_78d | p_78g) +
  plot_layout(heights = c(3, 3, 3)) +
  plot_annotation(
    title = "Section 78: Stoichiometry & Mechanism with Unbiased Neuronal Set",
    subtitle = sprintf("Broad neuronal: %s genes (GO-derived) | Narrow: %d (DMR-enriched) | Synapse/axon: %s",
                        format(sum(df$is_neuronal_broad), big.mark = ","),
                        sum(df$is_neuronal_narrow),
                        format(sum(df$is_synapse), big.mark = ",")),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

save_plot(p_78h, "78_composite", w = 18, h = 22)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 78 SUMMARY\n")
cat("================================================================================\n\n")

cat("  KEY COMPARISONS (narrow DMR-enriched vs broad GO-derived neuronal set):\n\n")

if (has_narrow && sum(df$is_neuronal_narrow) >= 10) {
  narrow_mean <- mean(df$delta_total[df$is_neuronal_narrow & !is.na(df$delta_total)])
  broad_mean  <- mean(df$delta_total[df$is_neuronal_broad & !is.na(df$delta_total)])
  cat(sprintf("    Narrow neuronal delta_total mean: %+.5f (%s)\n",
              narrow_mean, ifelse(narrow_mean > 0, "INCREASES", "DECREASES")))
  cat(sprintf("    Broad neuronal delta_total mean:  %+.5f (%s)\n",
              broad_mean, ifelse(broad_mean > 0, "INCREASES", "DECREASES")))

  if (sign(narrow_mean) != sign(broad_mean)) {
    cat("    >>> DIRECTION FLIPS between narrow and broad sets <<<\n")
    cat("    The narrow set's total methylation increase was driven by DMR selection bias.\n")
  }
}

narrow_slope <- slope_results$slope[slope_results$group == "Neuronal (broad)"]
if (length(narrow_slope) > 0) {
  cat(sprintf("\n    Broad neuronal stoichiometry slope: %.4f (ref: -1.0)\n", narrow_slope))
}

cat("\n")
cat("================================================================================\n")
cat("SECTION 78 COMPLETE\n")
cat("================================================================================\n")
