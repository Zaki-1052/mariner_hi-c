# biomodal/downstream/scripts/viz_sections/section_30_polycomb_target_enrichment.R
# Section 30: Polycomb Target Gene Systematic Enrichment (TODO #17)
# Tests whether classic Polycomb target genes (H3K27me3-marked, repressed)
# are the primary hypermethylation targets in BAP1-KO. The dual-mechanism
# model predicts they should NOT be — heterochromatin is inaccessible to
# DNMT3A, so normally active genes gaining ectopic H2AK119ub should be
# hypermethylated instead.
#
# Sub-analyses:
#   30a: Stacked bar — Polycomb vs Non-Polycomb hyper/hypo proportions
#   30b: Forest plot — Fisher's OR across all Polycomb definitions
#   30c: Violin — mC change magnitude, Polycomb vs Non-Polycomb by direction
#   30d: Bar chart — % hypermethylated per chromatin state (all 7)
#   30e: Violin — hmC change magnitude, same layout as 30c
#   30f: Composite summary panel (30a + 30b + 30d)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_30_polycomb_target_enrichment.R

source("scripts/viz_sections/_shared_config.R")

cat("\n")
cat("================================================================================\n")
cat("SECTION 30: POLYCOMB TARGET GENE SYSTEMATIC ENRICHMENT\n")
cat("================================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

RATIO_FILE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")

SECTION_DIR <- file.path(OUTPUT_DIR, "30_polycomb_enrichment")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# Helper: significance stars
sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# =============================================================================
# STEP 1: LOAD & VALIDATE DATA
# =============================================================================

cat("================================================================================\n")
cat("STEP 1: Load and validate data\n")
cat("================================================================================\n\n")

stopifnot("Ratio file not found" = file.exists(RATIO_FILE))

all_genes <- read.table(RATIO_FILE, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "")

# Validate required columns
required_cols <- c("gene", "chromatin_state", "mc_diff", "mc_sig",
                   "hmc_diff", "hmc_sig", "dmr_status")
missing_cols <- setdiff(required_cols, colnames(all_genes))
stopifnot("Missing required columns" = length(missing_cols) == 0)

cat(sprintf("  Loaded %d genes from %s\n", nrow(all_genes), basename(RATIO_FILE)))
cat(sprintf("  Chromatin state breakdown:\n"))
cs_counts <- table(all_genes$chromatin_state)
for (state in names(sort(cs_counts, decreasing = TRUE))) {
  cat(sprintf("    %s: %d\n", state, cs_counts[state]))
}

# =============================================================================
# STEP 2: DEFINE POLYCOMB GENE SETS (Task 17a)
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 2: Define Polycomb gene sets\n")
cat("================================================================================\n\n")

# --- Set 1: Chromatin State Polycomb (Repressed_Promoter + Polycomb + Bivalent_Promoter) ---
all_genes$is_polycomb_broad <- all_genes$chromatin_state %in%
  c("Repressed_Promoter", "Polycomb", "Bivalent_Promoter")

# --- Set 2: Strict Polycomb (excludes Bivalent) ---
all_genes$is_polycomb_strict <- all_genes$chromatin_state %in%
  c("Repressed_Promoter", "Polycomb")

# --- Set 3: Broad H3K27me3 overlap (gene body overlaps H3K27me3 peak) ---
cat("  Computing H3K27me3 gene body overlaps...\n")
h3k27me3_gr <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
stopifnot("H3K27me3 peaks not loaded" = !is.null(h3k27me3_gr))

gene_gr <- GRanges(
  seqnames = all_genes$chr,
  ranges = IRanges(start = all_genes$start, end = all_genes$end)
)
all_genes$h3k27me3_overlap <- countOverlaps(gene_gr, h3k27me3_gr) > 0

# --- Optional: published Polycomb gene lists ---
polycomb_list_dir <- file.path(BASE_DIR, "data/polycomb_gene_lists")
published_lists <- list()
if (dir.exists(polycomb_list_dir)) {
  list_files <- list.files(polycomb_list_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(list_files) > 0) {
    for (f in list_files) {
      list_name <- tools::file_path_sans_ext(basename(f))
      genes_in_list <- readLines(f)
      genes_in_list <- genes_in_list[genes_in_list != ""]
      matched <- all_genes$gene %in% genes_in_list
      published_lists[[list_name]] <- matched
      cat(sprintf("  Published list '%s': %d genes, %d matched in data\n",
                  list_name, length(genes_in_list), sum(matched)))
    }
  } else {
    cat("  No published Polycomb gene lists found in data/polycomb_gene_lists/\n")
  }
} else {
  cat("  Directory data/polycomb_gene_lists/ does not exist — skipping published lists\n")
}

# Print summary
cat(sprintf("\n  Polycomb gene set summary:\n"))
cat(sprintf("    Chromatin State Polycomb (Rep+Pol+Biv): %d genes\n",
            sum(all_genes$is_polycomb_broad)))
cat(sprintf("    Strict Polycomb (Rep+Pol only):          %d genes\n",
            sum(all_genes$is_polycomb_strict)))
cat(sprintf("    Broad H3K27me3 (gene body overlap):     %d genes\n",
            sum(all_genes$h3k27me3_overlap)))

# Overlap between sets
overlap_strict_h3k27me3 <- sum(all_genes$is_polycomb_strict & all_genes$h3k27me3_overlap)
cat(sprintf("    Strict ∩ H3K27me3: %d genes\n", overlap_strict_h3k27me3))
cat(sprintf("    Broad chromatin ∩ H3K27me3: %d genes\n",
            sum(all_genes$is_polycomb_broad & all_genes$h3k27me3_overlap)))

# =============================================================================
# STEP 3: DEFINE DMR CATEGORIES
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 3: Define DMR categories\n")
cat("================================================================================\n\n")

all_genes$mc_hyper <- all_genes$mc_sig & all_genes$mc_diff > 0
all_genes$mc_hypo  <- all_genes$mc_sig & all_genes$mc_diff < 0
all_genes$hmc_hyper <- all_genes$hmc_sig & all_genes$hmc_diff > 0
all_genes$hmc_hypo  <- all_genes$hmc_sig & all_genes$hmc_diff < 0

cat(sprintf("  mC hyper: %d | mC hypo: %d | mC not-sig: %d\n",
            sum(all_genes$mc_hyper), sum(all_genes$mc_hypo),
            sum(!all_genes$mc_sig)))
cat(sprintf("  hmC hyper: %d | hmC hypo: %d | hmC not-sig: %d\n",
            sum(all_genes$hmc_hyper), sum(all_genes$hmc_hypo),
            sum(!all_genes$hmc_sig)))
cat(sprintf("  Universe: %d genes\n", nrow(all_genes)))

# =============================================================================
# STEP 4: FISHER'S EXACT TESTS (Task 17b)
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 4: Fisher's exact tests — Polycomb x DMR enrichment\n")
cat("================================================================================\n\n")

# Reusable Fisher's test helper
run_polycomb_fisher <- function(universe_df, polycomb_flag, dmr_flag, test_name) {
  a <- sum(polycomb_flag & dmr_flag)
  b <- sum(!polycomb_flag & dmr_flag)
  c <- sum(polycomb_flag & !dmr_flag)
  d <- sum(!polycomb_flag & !dmr_flag)

  fisher_mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                       dimnames = list(c("DMR+", "DMR-"),
                                       c("Polycomb+", "Polycomb-")))
  fisher_res <- fisher.test(fisher_mat)

  pct_polycomb_dmr <- ifelse((a + c) > 0, 100 * a / (a + c), 0)
  pct_non_polycomb_dmr <- ifelse((b + d) > 0, 100 * b / (b + d), 0)

  cat(sprintf("  %s:\n", test_name))
  cat(sprintf("    OR = %.3f (95%% CI: %.3f-%.3f), %s\n",
              fisher_res$estimate, fisher_res$conf.int[1], fisher_res$conf.int[2],
              fmt_p(fisher_res$p.value)))
  cat(sprintf("    Polycomb+ DMR rate: %.1f%% | Non-Polycomb DMR rate: %.1f%%\n",
              pct_polycomb_dmr, pct_non_polycomb_dmr))
  cat(sprintf("    2x2: a=%d, b=%d, c=%d, d=%d\n\n", a, b, c, d))

  data.frame(
    test = test_name,
    n_universe = nrow(universe_df),
    n_polycomb = sum(polycomb_flag),
    n_dmr = sum(dmr_flag),
    a = a, b = b, c = c, d = d,
    odds_ratio = as.numeric(fisher_res$estimate),
    ci_lower = fisher_res$conf.int[1],
    ci_upper = fisher_res$conf.int[2],
    p_value = fisher_res$p.value,
    pct_polycomb_dmr = pct_polycomb_dmr,
    pct_non_polycomb_dmr = pct_non_polycomb_dmr,
    stringsAsFactors = FALSE
  )
}

# --- Core tests: per Polycomb definition × direction × mark ---
fisher_results <- list()

polycomb_defs <- list(
  "Chromatin State" = all_genes$is_polycomb_broad,
  "Strict (no Bivalent)" = all_genes$is_polycomb_strict,
  "H3K27me3 overlap" = all_genes$h3k27me3_overlap
)

for (def_name in names(polycomb_defs)) {
  pc_flag <- polycomb_defs[[def_name]]

  fisher_results <- c(fisher_results, list(
    run_polycomb_fisher(all_genes, pc_flag, all_genes$mc_hyper,
                        sprintf("%s × mC hyper", def_name)),
    run_polycomb_fisher(all_genes, pc_flag, all_genes$mc_hypo,
                        sprintf("%s × mC hypo", def_name)),
    run_polycomb_fisher(all_genes, pc_flag, all_genes$hmc_hyper,
                        sprintf("%s × hmC hyper", def_name)),
    run_polycomb_fisher(all_genes, pc_flag, all_genes$hmc_hypo,
                        sprintf("%s × hmC hypo", def_name))
  ))
}

# --- Published lists (if any) ---
for (list_name in names(published_lists)) {
  pc_flag <- published_lists[[list_name]]
  fisher_results <- c(fisher_results, list(
    run_polycomb_fisher(all_genes, pc_flag, all_genes$mc_hyper,
                        sprintf("Published:%s × mC hyper", list_name)),
    run_polycomb_fisher(all_genes, pc_flag, all_genes$mc_hypo,
                        sprintf("Published:%s × mC hypo", list_name))
  ))
}

# --- Per-chromatin-state tests (Task 17d) ---
cat("  --- Per-chromatin-state hypermethylation enrichment ---\n\n")

per_state_results <- list()
for (state in CHROMATIN_STATE_ORDER) {
  is_state <- all_genes$chromatin_state == state
  per_state_results <- c(per_state_results, list(
    run_polycomb_fisher(all_genes, is_state, all_genes$mc_hyper,
                        sprintf("State:%s × mC hyper", state))
  ))
}

# Combine all Fisher's results
fisher_all <- dplyr::bind_rows(c(fisher_results, per_state_results))

# BH FDR correction across all tests
fisher_all$q_value <- p.adjust(fisher_all$p_value, method = "BH")

cat(sprintf("  Total Fisher's tests: %d (BH-corrected)\n\n", nrow(fisher_all)))

# Save Fisher's results
write.table(fisher_all,
            file.path(TABLES_DIR, "polycomb_fisher_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: polycomb_fisher_tests.tsv (%d tests)\n", nrow(fisher_all)))

# =============================================================================
# STEP 5: METHYLATION MAGNITUDE COMPARISON (Task 17c)
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 5: Methylation magnitude comparison — Polycomb vs Non-Polycomb\n")
cat("================================================================================\n\n")

# Use the broad chromatin-state definition for magnitude analysis
magnitude_results <- list()

for (mark in c("mc", "hmc")) {
  diff_col <- paste0(mark, "_diff")
  sig_col <- paste0(mark, "_sig")

  for (direction in c("hyper", "hypo")) {
    if (direction == "hyper") {
      subset_idx <- all_genes[[sig_col]] & all_genes[[diff_col]] > 0
    } else {
      subset_idx <- all_genes[[sig_col]] & all_genes[[diff_col]] < 0
    }

    if (sum(subset_idx) < 5) {
      cat(sprintf("  %s %s: too few genes (%d), skipping\n", mark, direction, sum(subset_idx)))
      next
    }

    polycomb_vals <- abs(all_genes[[diff_col]][subset_idx & all_genes$is_polycomb_broad])
    non_polycomb_vals <- abs(all_genes[[diff_col]][subset_idx & !all_genes$is_polycomb_broad])

    if (length(polycomb_vals) < 3 || length(non_polycomb_vals) < 3) {
      cat(sprintf("  %s %s: insufficient Polycomb (%d) or Non-Polycomb (%d) genes\n",
                  mark, direction, length(polycomb_vals), length(non_polycomb_vals)))
      next
    }

    wilcox_res <- wilcox.test(polycomb_vals, non_polycomb_vals)

    cat(sprintf("  %s %s: Polycomb median=%.4f (n=%d) vs Non-Polycomb median=%.4f (n=%d), %s\n",
                toupper(mark), direction,
                median(polycomb_vals), length(polycomb_vals),
                median(non_polycomb_vals), length(non_polycomb_vals),
                fmt_p(wilcox_res$p.value)))

    magnitude_results <- c(magnitude_results, list(data.frame(
      mark = toupper(mark),
      direction = direction,
      polycomb_n = length(polycomb_vals),
      polycomb_median = median(polycomb_vals),
      polycomb_mean = mean(polycomb_vals),
      non_polycomb_n = length(non_polycomb_vals),
      non_polycomb_median = median(non_polycomb_vals),
      non_polycomb_mean = mean(non_polycomb_vals),
      wilcox_p = wilcox_res$p.value,
      wilcox_stat = wilcox_res$statistic,
      stringsAsFactors = FALSE
    )))
  }
}

magnitude_df <- dplyr::bind_rows(magnitude_results)
magnitude_df$wilcox_q <- p.adjust(magnitude_df$wilcox_p, method = "BH")

write.table(magnitude_df,
            file.path(TABLES_DIR, "polycomb_methylation_magnitude.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: polycomb_methylation_magnitude.tsv (%d comparisons)\n", nrow(magnitude_df)))

# =============================================================================
# STEP 6: VISUALIZATIONS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 6: Generating visualizations\n")
cat("================================================================================\n\n")

# --- 30a: Stacked bar — Polycomb vs Non-Polycomb DMR proportions ---
cat("  30a: Stacked bar — Polycomb vs Non-Polycomb...\n")

all_genes$polycomb_label <- ifelse(all_genes$is_polycomb_broad,
                                    "Polycomb\ntargets", "Non-Polycomb")
all_genes$mc_status <- case_when(
  all_genes$mc_hyper ~ "Hypermethylated",
  all_genes$mc_hypo  ~ "Hypomethylated",
  TRUE               ~ "Not significant"
)
all_genes$mc_status <- factor(all_genes$mc_status,
                               levels = c("Hypermethylated", "Hypomethylated", "Not significant"))

# Fisher's result for annotation
fisher_broad_hyper <- fisher_all %>%
  dplyr::filter(test == "Chromatin State × mC hyper")

bar_data_30a <- all_genes %>%
  dplyr::count(polycomb_label, mc_status) %>%
  dplyr::group_by(polycomb_label) %>%
  dplyr::mutate(pct = 100 * n / sum(n),
                total = sum(n)) %>%
  dplyr::ungroup()

p_30a <- ggplot(bar_data_30a, aes(x = polycomb_label, y = pct, fill = mc_status)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("Hypermethylated" = "#D7191C",
                                "Hypomethylated" = "#2C7BB6",
                                "Not significant" = "grey80")) +
  annotate("text", x = 1.5, y = 102,
           label = sprintf("Fisher's hyper: OR=%.2f, %s",
                            fisher_broad_hyper$odds_ratio,
                            fmt_p(fisher_broad_hyper$p_value)),
           size = 3.2, fontface = "italic") +
  labs(
    title = "mC DMR Status: Polycomb Targets vs Non-Polycomb",
    subtitle = sprintf("Polycomb = Repressed_Promoter + Polycomb + Bivalent (n=%d) | Non-Polycomb (n=%d)",
                       sum(all_genes$is_polycomb_broad),
                       sum(!all_genes$is_polycomb_broad)),
    x = NULL, y = "Percentage of genes",
    fill = "mC status"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_30a, file.path(SECTION_DIR, "30a_polycomb_vs_non_polycomb_stacked_bar"),
                        width = 8, height = 7)

# --- 30b: Forest plot — Fisher's ORs across definitions ---
cat("  30b: Forest plot — Fisher's ORs...\n")

# Select core tests for forest plot (exclude per-state tests)
core_tests <- fisher_all %>%
  dplyr::filter(!grepl("^State:", test))

core_tests$sig_label <- sapply(core_tests$q_value, sig_stars)

# Clean test names for display
core_tests$test_clean <- core_tests$test %>%
  gsub("Chromatin State", "ChrmSt", .) %>%
  gsub("Strict \\(no Bivalent\\)", "Strict", .) %>%
  gsub("H3K27me3 overlap", "H3K27me3", .) %>%
  gsub("Published:", "Pub:", .)

p_30b <- ggplot(core_tests, aes(x = odds_ratio, y = reorder(test_clean, odds_ratio))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, linewidth = 0.7) +
  geom_point(aes(color = ifelse(q_value < 0.05, "Significant", "Not significant")),
             size = 3) +
  geom_text(aes(label = sprintf("%.2f %s", odds_ratio, sig_label)),
            hjust = -0.15, size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not significant" = "grey50")) +
  labs(
    title = "Polycomb Target Enrichment in DMRs",
    subtitle = "Odds Ratios with 95% CI (BH-corrected significance)",
    x = "Odds Ratio (log scale)", y = NULL,
    color = "BH-adjusted"
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

save_multiformat_ggplot(p_30b, file.path(SECTION_DIR, "30b_fisher_forest_plot"),
                        width = 11, height = 8)

# --- 30c: Violin — mC change magnitude, Polycomb vs Non-Polycomb ---
cat("  30c: mC magnitude violin...\n")

mc_sig_genes <- all_genes %>%
  dplyr::filter(mc_sig) %>%
  dplyr::mutate(
    mc_direction = ifelse(mc_diff > 0, "Hypermethylated", "Hypomethylated"),
    polycomb_group = ifelse(is_polycomb_broad, "Polycomb", "Non-Polycomb"),
    abs_mc_diff = abs(mc_diff)
  )

mc_sig_genes$polycomb_group <- factor(mc_sig_genes$polycomb_group,
                                       levels = c("Non-Polycomb", "Polycomb"))
mc_sig_genes$mc_direction <- factor(mc_sig_genes$mc_direction,
                                     levels = c("Hypermethylated", "Hypomethylated"))

# Pairwise Wilcoxon for subtitle
wilcox_mc_hyper <- wilcox.test(
  abs_mc_diff ~ polycomb_group,
  data = mc_sig_genes %>% dplyr::filter(mc_direction == "Hypermethylated")
)
wilcox_mc_hypo <- wilcox.test(
  abs_mc_diff ~ polycomb_group,
  data = mc_sig_genes %>% dplyr::filter(mc_direction == "Hypomethylated")
)

p_30c <- ggplot(mc_sig_genes, aes(x = polycomb_group, y = abs_mc_diff, fill = polycomb_group)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  facet_wrap(~mc_direction) +
  scale_fill_manual(values = c("Non-Polycomb" = "#4575b4", "Polycomb" = "#756bb1")) +
  labs(
    title = "mC Change Magnitude: Polycomb vs Non-Polycomb Targets",
    subtitle = sprintf("Hyper Wilcoxon %s | Hypo Wilcoxon %s",
                       fmt_p(wilcox_mc_hyper$p.value),
                       fmt_p(wilcox_mc_hypo$p.value)),
    x = NULL, y = "|mC mod_difference|",
    fill = NULL
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_30c, file.path(SECTION_DIR, "30c_mc_magnitude_violin"),
                        width = 10, height = 7)

# --- 30d: Bar chart — % hypermethylated per chromatin state ---
cat("  30d: Per-chromatin-state hypermethylation rates...\n")

per_state_enrichment <- all_genes %>%
  dplyr::group_by(chromatin_state) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_mc_hyper = sum(mc_hyper),
    n_mc_hypo = sum(mc_hypo),
    n_mc_sig = sum(mc_sig),
    pct_hyper = 100 * sum(mc_hyper) / dplyr::n(),
    pct_hypo = 100 * sum(mc_hypo) / dplyr::n(),
    pct_sig = 100 * sum(mc_sig) / dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(chromatin_state = factor(chromatin_state, levels = CHROMATIN_STATE_ORDER))

# Add Fisher's test results (per-state tests from Step 4)
per_state_fisher <- fisher_all %>%
  dplyr::filter(grepl("^State:", test)) %>%
  dplyr::mutate(chromatin_state = gsub("State:(.*) × mC hyper", "\\1", test)) %>%
  dplyr::select(chromatin_state, or_state = odds_ratio, p_state = p_value, q_state = q_value)

per_state_enrichment <- per_state_enrichment %>%
  dplyr::left_join(per_state_fisher, by = "chromatin_state")

per_state_enrichment$sig_label <- sapply(per_state_enrichment$q_state, sig_stars)

# Genome-wide hypermethylation rate
genome_wide_hyper_rate <- 100 * sum(all_genes$mc_hyper) / nrow(all_genes)

p_30d <- ggplot(per_state_enrichment,
                aes(x = chromatin_state, y = pct_hyper, fill = chromatin_state)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = genome_wide_hyper_rate, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  annotate("text", x = 0.7, y = genome_wide_hyper_rate + 1.5,
           label = sprintf("Genome-wide: %.1f%%", genome_wide_hyper_rate),
           hjust = 0, size = 3, fontface = "italic") +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)\n%s", pct_hyper, n_total, sig_label)),
            vjust = -0.2, size = 2.8) +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
  labs(
    title = "Hypermethylation Rate by Chromatin State",
    subtitle = "Fisher's exact test vs all other states (BH-corrected: *** <0.001, ** <0.01, * <0.05)",
    x = NULL, y = "% genes hypermethylated (mC)",
    fill = "Chromatin state"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, max(per_state_enrichment$pct_hyper) * 1.25))

save_multiformat_ggplot(p_30d, file.path(SECTION_DIR, "30d_per_state_hypermethylation_rate"),
                        width = 10, height = 7)

# Save per-state enrichment table
write.table(per_state_enrichment,
            file.path(TABLES_DIR, "polycomb_per_chromatin_state_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: polycomb_per_chromatin_state_enrichment.tsv\n")

# --- 30e: Violin — hmC change magnitude ---
cat("  30e: hmC magnitude violin...\n")

hmc_sig_genes <- all_genes %>%
  dplyr::filter(hmc_sig) %>%
  dplyr::mutate(
    hmc_direction = ifelse(hmc_diff > 0, "hmC gain", "hmC loss"),
    polycomb_group = ifelse(is_polycomb_broad, "Polycomb", "Non-Polycomb"),
    abs_hmc_diff = abs(hmc_diff)
  )

hmc_sig_genes$polycomb_group <- factor(hmc_sig_genes$polycomb_group,
                                        levels = c("Non-Polycomb", "Polycomb"))
hmc_sig_genes$hmc_direction <- factor(hmc_sig_genes$hmc_direction,
                                       levels = c("hmC gain", "hmC loss"))

# Check if there are enough genes in each group for Wilcoxon tests
hmc_gain_counts <- table(hmc_sig_genes$polycomb_group[hmc_sig_genes$hmc_direction == "hmC gain"])
hmc_loss_counts <- table(hmc_sig_genes$polycomb_group[hmc_sig_genes$hmc_direction == "hmC loss"])

subtitle_parts <- c()
if (all(hmc_gain_counts >= 3)) {
  wilcox_hmc_gain <- wilcox.test(
    abs_hmc_diff ~ polycomb_group,
    data = hmc_sig_genes %>% dplyr::filter(hmc_direction == "hmC gain")
  )
  subtitle_parts <- c(subtitle_parts, sprintf("Gain Wilcoxon %s", fmt_p(wilcox_hmc_gain$p.value)))
} else {
  subtitle_parts <- c(subtitle_parts, "Gain: too few Polycomb genes")
}

if (all(hmc_loss_counts >= 3)) {
  wilcox_hmc_loss <- wilcox.test(
    abs_hmc_diff ~ polycomb_group,
    data = hmc_sig_genes %>% dplyr::filter(hmc_direction == "hmC loss")
  )
  subtitle_parts <- c(subtitle_parts, sprintf("Loss Wilcoxon %s", fmt_p(wilcox_hmc_loss$p.value)))
} else {
  subtitle_parts <- c(subtitle_parts, "Loss: too few Polycomb genes")
}

p_30e <- ggplot(hmc_sig_genes, aes(x = polycomb_group, y = abs_hmc_diff, fill = polycomb_group)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  facet_wrap(~hmc_direction) +
  scale_fill_manual(values = c("Non-Polycomb" = "#4575b4", "Polycomb" = "#756bb1")) +
  labs(
    title = "hmC Change Magnitude: Polycomb vs Non-Polycomb Targets",
    subtitle = paste(subtitle_parts, collapse = " | "),
    x = NULL, y = "|hmC mod_difference|",
    fill = NULL
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_30e, file.path(SECTION_DIR, "30e_hmc_magnitude_violin"),
                        width = 10, height = 7)

# --- 30f: Composite panel ---
cat("  30f: Composite panel...\n")

p_30f <- (p_30a | p_30b) / p_30d +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(
    title = "Polycomb Target Gene Enrichment in BAP1-KO Methylation Changes",
    subtitle = "Dual-mechanism prediction: Polycomb targets depleted from hypermethylation (heterochromatin inaccessible to DNMT3A)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic")
    )
  )

save_multiformat_ggplot(p_30f, file.path(SECTION_DIR, "30f_composite_polycomb_summary"),
                        width = 18, height = 14)

# =============================================================================
# STEP 7: EXPORT GENE CLASSIFICATION TABLE
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 7: Exporting gene classification table\n")
cat("================================================================================\n\n")

gene_classification <- all_genes %>%
  dplyr::select(gene, chr, start, end, chromatin_state,
                is_polycomb_broad, is_polycomb_strict, h3k27me3_overlap,
                mc_diff, mc_sig, mc_hyper, mc_hypo,
                hmc_diff, hmc_sig, hmc_hyper, hmc_hypo,
                dmr_status)

write.table(gene_classification,
            file.path(TABLES_DIR, "polycomb_gene_classification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: polycomb_gene_classification.tsv (%d genes)\n", nrow(gene_classification)))

# =============================================================================
# STEP 8: CONSOLE SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 30 SUMMARY: POLYCOMB TARGET GENE ENRICHMENT\n")
cat("================================================================================\n\n")

cat(sprintf("Universe: %d genes with chromatin state and DMR data\n\n", nrow(all_genes)))

cat("--- Polycomb Gene Set Sizes ---\n")
cat(sprintf("  Chromatin State (Rep+Pol+Biv): %d genes (%.1f%%)\n",
            sum(all_genes$is_polycomb_broad),
            100 * mean(all_genes$is_polycomb_broad)))
cat(sprintf("  Strict (Rep+Pol):              %d genes (%.1f%%)\n",
            sum(all_genes$is_polycomb_strict),
            100 * mean(all_genes$is_polycomb_strict)))
cat(sprintf("  H3K27me3 gene body overlap:    %d genes (%.1f%%)\n",
            sum(all_genes$h3k27me3_overlap),
            100 * mean(all_genes$h3k27me3_overlap)))

cat("\n--- Key Fisher's Tests (Polycomb × mC hyper) ---\n")
hyper_tests <- fisher_all %>%
  dplyr::filter(grepl("× mC hyper", test) & !grepl("^State:", test))
for (i in seq_len(nrow(hyper_tests))) {
  row <- hyper_tests[i, ]
  cat(sprintf("  %s: OR = %.3f (%.3f-%.3f), %s, q = %.2e\n",
              row$test, row$odds_ratio, row$ci_lower, row$ci_upper,
              fmt_p(row$p_value), row$q_value))
}

cat("\n--- Key Fisher's Tests (Polycomb × mC hypo) ---\n")
hypo_tests <- fisher_all %>%
  dplyr::filter(grepl("× mC hypo", test) & !grepl("^State:", test))
for (i in seq_len(nrow(hypo_tests))) {
  row <- hypo_tests[i, ]
  cat(sprintf("  %s: OR = %.3f (%.3f-%.3f), %s, q = %.2e\n",
              row$test, row$odds_ratio, row$ci_lower, row$ci_upper,
              fmt_p(row$p_value), row$q_value))
}

cat("\n--- Per-Chromatin-State Hypermethylation Rates ---\n")
cat(sprintf("  %-22s %8s %8s %8s %6s\n", "State", "N genes", "N hyper", "% hyper", "OR"))
for (i in seq_len(nrow(per_state_enrichment))) {
  row <- per_state_enrichment[i, ]
  cat(sprintf("  %-22s %8d %8d %7.1f%% %6.2f %s\n",
              as.character(row$chromatin_state), row$n_total, row$n_mc_hyper,
              row$pct_hyper, row$or_state, row$sig_label))
}
cat(sprintf("  %-22s %8d %8d %7.1f%%\n",
            "GENOME-WIDE", nrow(all_genes), sum(all_genes$mc_hyper),
            genome_wide_hyper_rate))

cat("\n--- Biological Interpretation ---\n")
# Check if Polycomb targets are depleted from hypermethylation
cs_hyper_or <- fisher_all %>%
  dplyr::filter(test == "Chromatin State × mC hyper")
cs_hypo_or <- fisher_all %>%
  dplyr::filter(test == "Chromatin State × mC hypo")

if (cs_hyper_or$odds_ratio < 1 && cs_hyper_or$q_value < 0.05) {
  cat("  CONFIRMED: Polycomb target genes are DEPLETED from hypermethylation\n")
  cat(sprintf("    OR = %.3f (q = %.2e) — consistent with dual-mechanism model\n",
              cs_hyper_or$odds_ratio, cs_hyper_or$q_value))
  cat("    Interpretation: heterochromatin is inaccessible to DNMT3A;\n")
  cat("    hypermethylation occurs at normally active genes gaining ectopic H2AK119ub\n")
} else if (cs_hyper_or$odds_ratio > 1 && cs_hyper_or$q_value < 0.05) {
  cat("  UNEXPECTED: Polycomb targets are ENRICHED in hypermethylation\n")
  cat(sprintf("    OR = %.3f (q = %.2e) — inconsistent with dual-mechanism model\n",
              cs_hyper_or$odds_ratio, cs_hyper_or$q_value))
  cat("    This would suggest DNMT3A can access Polycomb-marked chromatin,\n")
  cat("    or that Polycomb derepression precedes hypermethylation\n")
} else {
  cat("  INCONCLUSIVE: Polycomb targets show no significant enrichment/depletion\n")
  cat(sprintf("    OR = %.3f (q = %.2e)\n", cs_hyper_or$odds_ratio, cs_hyper_or$q_value))
}

if (cs_hypo_or$odds_ratio > 1 && cs_hypo_or$q_value < 0.05) {
  cat(sprintf("\n  Polycomb targets are ENRICHED in hypomethylation (OR = %.3f, q = %.2e)\n",
              cs_hypo_or$odds_ratio, cs_hypo_or$q_value))
  cat("    Consistent with Polycomb-mediated repression maintaining low mC baseline\n")
}

cat(sprintf("\n  Plots saved to: %s\n", SECTION_DIR))
cat(sprintf("  Tables saved to: %s\n", TABLES_DIR))
cat("\nSection 30 complete.\n")
