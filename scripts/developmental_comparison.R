# scripts/developmental_comparison.R
# Developmental Comparison: Does adult BAP1-KO resemble P13 wildtype?
#
# Tests the "blocked developmental remodeling" hypothesis by comparing
# adult-KO 3D genome structure to P13-WT across multiple structural layers.

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("================================================\n")
cat("Developmental Comparison: Adult-KO vs P13-WT\n")
cat("================================================\n\n")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(scales)
  library(RColorBrewer)
})

source("scripts/utils/multi_format_output.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILES <- list(
  early_all_loops = "outputs/250831-early_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe",
  late_all_loops  = "outputs/250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe",
  early_pc1_raw   = "tads/tad-pc-analysis/inputs/new/all_PC1.txt",
  late_pc1_raw    = "tads/tad-pc-analysis/inputs/late/diffPC/all_PC1.txt",
  early_tad       = "tads/results/early/final/tadcompare_final_filtered.tsv",
  late_tad        = "tads/results/late/final/tadcompare_final_filtered.tsv",
  early_stripes   = "stripes/stripenn/outputs/250831/cross_res_merged.tsv",
  late_stripes    = "stripes/stripenn/outputs/250402/cross_res_merged.tsv"
)

ANCHOR_TOLERANCE <- 10000
TAD_TOLERANCE    <- 20000
STRIPE_TOLERANCE <- 50000
N_PERM           <- 10000
SEED             <- 42

AUTOSOMES <- paste0("chr", 1:19)

COLORS <- c(
  early_WT = "#7570b3",
  early_KO = "#b3b0d8",
  late_WT  = "#1b9e77",
  late_KO  = "#66c2a5"
)

args <- commandArgs(trailingOnly = TRUE)
OUTPUT_DIR <- "output/developmental_comparison"
for (i in seq_along(args)) {
  if (args[i] == "--output" && i < length(args)) OUTPUT_DIR <- args[i + 1]
  if (args[i] == "--n-perm" && i < length(args)) N_PERM <- as.integer(args[i + 1])
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

check_inputs <- function(file_list) {
  missing <- character(0)
  for (name in names(file_list)) {
    if (!file.exists(file_list[[name]])) {
      missing <- c(missing, sprintf("  %s: %s", name, file_list[[name]]))
    }
  }
  if (length(missing) > 0) {
    stop("Missing input files:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  cat(sprintf("All %d input files verified.\n\n", length(file_list)))
}

load_pc1_matrix <- function(file_path, timepoint_label) {
  raw <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE)
  region_id <- raw[[1]]
  pc1_cols <- raw[, (ncol(raw) - 5):ncol(raw)]
  colnames(pc1_cols) <- paste0(timepoint_label, "_",
                               c("ctrl_M1", "ctrl_M2", "ctrl_M3",
                                 "mut_M1", "mut_M2", "mut_M3"))
  chr_col <- raw[[2]]
  result <- data.frame(region_id = region_id, chr = chr_col,
                       pc1_cols, stringsAsFactors = FALSE)
  result <- result[result$chr %in% AUTOSOMES, ]
  result$chr <- NULL
  result
}

load_bedpe <- function(file_path) {
  read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
}

build_loop_set <- function(bedpe_df, directions) {
  bedpe_df[bedpe_df$direction %in% directions, ]
}

compute_loop_distance <- function(bedpe_df) {
  abs(((bedpe_df$start2 + bedpe_df$end2) / 2) -
      ((bedpe_df$start1 + bedpe_df$end1) / 2))
}

find_matching_loops <- function(query_df, ref_df, tolerance_bp) {
  query_a1 <- GRanges(seqnames = query_df$chr1,
                      ranges = IRanges(start = query_df$start1,
                                       end = query_df$end1))
  query_a2 <- GRanges(seqnames = query_df$chr2,
                      ranges = IRanges(start = query_df$start2,
                                       end = query_df$end2))
  ref_a1 <- GRanges(seqnames = ref_df$chr1,
                    ranges = IRanges(start = ref_df$start1,
                                     end = ref_df$end1))
  ref_a2 <- GRanges(seqnames = ref_df$chr2,
                    ranges = IRanges(start = ref_df$start2,
                                     end = ref_df$end2))

  hits_a1 <- findOverlaps(query_a1, ref_a1, maxgap = tolerance_bp)

  if (length(hits_a1) == 0) return(0L)

  q_idx <- queryHits(hits_a1)
  s_idx <- subjectHits(hits_a1)

  q_a2_start <- start(query_a2[q_idx])
  q_a2_end   <- end(query_a2[q_idx])
  r_a2_start <- start(ref_a2[s_idx])
  r_a2_end   <- end(ref_a2[s_idx])

  gap_a2 <- pmax(q_a2_start - r_a2_end, r_a2_start - q_a2_end, 0L)
  both_match <- gap_a2 <= tolerance_bp

  length(unique(q_idx[both_match]))
}

jaccard_loops <- function(set_a_df, set_b_df, tolerance_bp) {
  n_ab <- find_matching_loops(set_a_df, set_b_df, tolerance_bp)
  n_ba <- find_matching_loops(set_b_df, set_a_df, tolerance_bp)
  n_intersection <- (n_ab + n_ba) / 2
  n_union <- nrow(set_a_df) + nrow(set_b_df) - n_intersection
  n_intersection / n_union
}

overlap_fraction <- function(query_df, ref_df, tolerance_bp) {
  n_match <- find_matching_loops(query_df, ref_df, tolerance_bp)
  n_match / nrow(query_df)
}

fisher_z_test <- function(r1, r2, n1, n2) {
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z <- (z1 - z2) / se
  p <- pnorm(abs(z), lower.tail = FALSE)
  list(z = z, p = p, r1 = r1, r2 = r2)
}

create_output_dirs <- function(base_dir) {
  dirs <- c("01_pc1_correlation", "02_loop_jaccard", "03_loop_overlap_gained",
            "04_distance_ecdf", "05_tad_boundary_concordance",
            "06_stripe_overlap", "tables")
  for (d in dirs) {
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
}

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Checking input files...\n")
check_inputs(INPUT_FILES)
create_output_dirs(OUTPUT_DIR)

# Collector for master summary stats
summary_stats <- list()


# =============================================================================
# MODULE 1: PC1 CORRELATION MATRIX
# =============================================================================

cat("=== Module 1: PC1 Correlation Matrix ===\n\n")

cat("Loading PC1 values...\n")
early_pc1 <- load_pc1_matrix(INPUT_FILES$early_pc1_raw, "early")
late_pc1  <- load_pc1_matrix(INPUT_FILES$late_pc1_raw, "late")

cat(sprintf("  Early: %d autosomal bins\n", nrow(early_pc1)))
cat(sprintf("  Late:  %d autosomal bins\n", nrow(late_pc1)))

pc1_joined <- inner_join(early_pc1, late_pc1, by = "region_id")
cat(sprintf("  Joined: %d shared bins\n\n", nrow(pc1_joined)))

pc1_mat <- as.matrix(pc1_joined[, -1])

cor_12x12 <- cor(pc1_mat, use = "pairwise.complete.obs")

group_avgs <- data.frame(
  early_WT = rowMeans(pc1_joined[, c("early_ctrl_M1", "early_ctrl_M2", "early_ctrl_M3")], na.rm = TRUE),
  early_KO = rowMeans(pc1_joined[, c("early_mut_M1", "early_mut_M2", "early_mut_M3")], na.rm = TRUE),
  late_WT  = rowMeans(pc1_joined[, c("late_ctrl_M1", "late_ctrl_M2", "late_ctrl_M3")], na.rm = TRUE),
  late_KO  = rowMeans(pc1_joined[, c("late_mut_M1", "late_mut_M2", "late_mut_M3")], na.rm = TRUE)
)
cor_4x4 <- cor(group_avgs, use = "pairwise.complete.obs")

n_bins <- nrow(pc1_joined)
r_lateKO_earlyWT <- cor_4x4["late_KO", "early_WT"]
r_lateWT_earlyWT <- cor_4x4["late_WT", "early_WT"]
r_lateKO_lateWT  <- cor_4x4["late_KO", "late_WT"]

fz <- fisher_z_test(r_lateKO_earlyWT, r_lateWT_earlyWT, n_bins, n_bins)

cat("PC1 Correlation Summary (group averages):\n")
cat(sprintf("  r(late_KO, early_WT) = %.4f\n", r_lateKO_earlyWT))
cat(sprintf("  r(late_WT, early_WT) = %.4f\n", r_lateWT_earlyWT))
cat(sprintf("  r(late_KO, late_WT)  = %.4f\n", r_lateKO_lateWT))
cat(sprintf("  Fisher Z-test (late_KO-earlyWT vs late_WT-earlyWT): z=%.3f, p=%.2e\n\n",
            fz$z, fz$p))

summary_stats$pc1_r_lateKO_earlyWT <- r_lateKO_earlyWT
summary_stats$pc1_r_lateWT_earlyWT <- r_lateWT_earlyWT
summary_stats$pc1_r_lateKO_lateWT  <- r_lateKO_lateWT
summary_stats$pc1_fisher_z <- fz$z
summary_stats$pc1_fisher_p <- fz$p
summary_stats$pc1_n_bins   <- n_bins

# -- PC1 replicate heatmap --
rep_labels <- gsub("(early|late)_(ctrl|mut)_(M[123])",
                   "\\1_\\2_\\3", colnames(cor_12x12))
colnames(cor_12x12) <- rep_labels
rownames(cor_12x12) <- rep_labels

ann_row <- data.frame(
  Timepoint = ifelse(grepl("^early", rep_labels), "Early (P13)", "Late (Adult)"),
  Condition = ifelse(grepl("_ctrl_", rep_labels), "Wildtype", "BAP1-KO"),
  row.names = rep_labels
)
ann_colors <- list(
  Timepoint = c("Early (P13)" = "#7570b3", "Late (Adult)" = "#1b9e77"),
  Condition = c("Wildtype" = "#4393c3", "BAP1-KO" = "#d6604d")
)

save_multiformat_pheatmap(
  quote(pheatmap(cor_12x12,
                 display_numbers = TRUE, number_format = "%.3f", fontsize_number = 7,
                 color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                 breaks = seq(0.5, 1, length.out = 101),
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 annotation_row = ann_row, annotation_col = ann_row,
                 annotation_colors = ann_colors,
                 main = "PC1 Correlation: 12 Replicates Across Timepoints")),
  file.path(OUTPUT_DIR, "01_pc1_correlation/pc1_correlation_replicate"),
  width = 12, height = 10
)

# -- PC1 4x4 summary heatmap --
save_multiformat_pheatmap(
  quote(pheatmap(cor_4x4,
                 display_numbers = TRUE, number_format = "%.4f", fontsize_number = 12,
                 color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                 breaks = seq(0.5, 1, length.out = 101),
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 main = "PC1 Correlation: Group Averages (4 Conditions)")),
  file.path(OUTPUT_DIR, "01_pc1_correlation/pc1_correlation_summary"),
  width = 7, height = 6
)

# -- Fisher Z test output --
writeLines(c(
  "Fisher Z-Transform Test: PC1 Correlation Comparison",
  "====================================================",
  sprintf("r(late_KO, early_WT) = %.6f", r_lateKO_earlyWT),
  sprintf("r(late_WT, early_WT) = %.6f", r_lateWT_earlyWT),
  sprintf("Difference           = %.6f", r_lateKO_earlyWT - r_lateWT_earlyWT),
  sprintf("n (shared bins)      = %d", n_bins),
  sprintf("Fisher Z-statistic   = %.4f", fz$z),
  sprintf("P-value (one-tailed) = %.4e", fz$p),
  "",
  "Interpretation:",
  if (fz$p < 0.05 && fz$z > 0) {
    "  Adult BAP1-KO PC1 is significantly MORE similar to P13-WT than Adult-WT is."
  } else if (fz$p < 0.05 && fz$z < 0) {
    "  Adult BAP1-KO PC1 is significantly LESS similar to P13-WT than Adult-WT is."
  } else {
    "  No significant difference in PC1 similarity to P13-WT between Adult-KO and Adult-WT."
  }
), file.path(OUTPUT_DIR, "01_pc1_correlation/pc1_fisher_z_test.txt"))

cat("Module 1 complete.\n\n")


# =============================================================================
# MODULE 2: LOOP JACCARD OVERLAP MATRIX
# =============================================================================

cat("=== Module 2: Loop Jaccard Overlap Matrix ===\n\n")

cat("Loading loop sets...\n")
early_bedpe <- load_bedpe(INPUT_FILES$early_all_loops)
late_bedpe  <- load_bedpe(INPUT_FILES$late_all_loops)

loop_sets <- list(
  early_WT = build_loop_set(early_bedpe, c("unchanged", "down_in_mutant")),
  early_KO = build_loop_set(early_bedpe, c("unchanged", "up_in_mutant")),
  late_WT  = build_loop_set(late_bedpe,  c("unchanged", "down_in_mutant")),
  late_KO  = build_loop_set(late_bedpe,  c("unchanged", "up_in_mutant"))
)

for (nm in names(loop_sets)) {
  cat(sprintf("  %s: %d loops\n", nm, nrow(loop_sets[[nm]])))
}
cat("\n")

cat("Computing pairwise Jaccard (10kb tolerance)...\n")
set_names <- names(loop_sets)
jaccard_mat <- matrix(NA, nrow = 4, ncol = 4,
                      dimnames = list(set_names, set_names))

pairs <- combn(set_names, 2)
for (i in seq_len(ncol(pairs))) {
  a <- pairs[1, i]
  b <- pairs[2, i]
  cat(sprintf("  %s vs %s... ", a, b))
  j <- jaccard_loops(loop_sets[[a]], loop_sets[[b]], ANCHOR_TOLERANCE)
  jaccard_mat[a, b] <- j
  jaccard_mat[b, a] <- j
  cat(sprintf("J = %.4f\n", j))
}
diag(jaccard_mat) <- 1.0

cat(sprintf("\nKey comparisons:\n"))
cat(sprintf("  J(late_KO, early_WT) = %.4f\n", jaccard_mat["late_KO", "early_WT"]))
cat(sprintf("  J(late_WT, early_WT) = %.4f\n", jaccard_mat["late_WT", "early_WT"]))
cat(sprintf("  J(late_KO, late_WT)  = %.4f\n", jaccard_mat["late_KO", "late_WT"]))

summary_stats$loop_J_lateKO_earlyWT <- jaccard_mat["late_KO", "early_WT"]
summary_stats$loop_J_lateWT_earlyWT <- jaccard_mat["late_WT", "early_WT"]
summary_stats$loop_J_lateKO_lateWT  <- jaccard_mat["late_KO", "late_WT"]

save_multiformat_pheatmap(
  quote(pheatmap(jaccard_mat,
                 display_numbers = TRUE, number_format = "%.4f", fontsize_number = 12,
                 color = colorRampPalette(c("white", "#fee08b", "#d73027"))(100),
                 breaks = seq(0, 1, length.out = 101),
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 main = "Loop Jaccard Similarity (10kb anchor tolerance)")),
  file.path(OUTPUT_DIR, "02_loop_jaccard/loop_jaccard_matrix"),
  width = 7, height = 6
)

jaccard_df <- as.data.frame(as.table(jaccard_mat))
colnames(jaccard_df) <- c("Set_A", "Set_B", "Jaccard")
write.table(jaccard_df, file.path(OUTPUT_DIR, "02_loop_jaccard/loop_jaccard_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nModule 2 complete.\n\n")


# =============================================================================
# MODULE 3: GAINED ADULT-KO LOOP OVERLAP WITH P13-WT
# =============================================================================

cat("=== Module 3: Gained Loop Overlap with P13-WT ===\n\n")

late_gained <- build_loop_set(late_bedpe, "up_in_mutant")
late_lost   <- build_loop_set(late_bedpe, "down_in_mutant")
late_bg     <- build_loop_set(late_bedpe, "unchanged")
early_wt    <- loop_sets$early_WT

cat(sprintf("  Late gained (up_in_mutant):   %d loops\n", nrow(late_gained)))
cat(sprintf("  Late lost (down_in_mutant):   %d loops\n", nrow(late_lost)))
cat(sprintf("  Late unchanged (background):  %d loops\n", nrow(late_bg)))
cat(sprintf("  Early WT (reference):         %d loops\n\n", nrow(early_wt)))

frac_gained <- overlap_fraction(late_gained, early_wt, ANCHOR_TOLERANCE)
frac_lost   <- overlap_fraction(late_lost, early_wt, ANCHOR_TOLERANCE)

cat(sprintf("Overlap with P13-WT:\n"))
cat(sprintf("  Late gained: %.1f%% (%d / %d)\n",
            frac_gained * 100, round(frac_gained * nrow(late_gained)), nrow(late_gained)))
cat(sprintf("  Late lost:   %.1f%% (%d / %d)\n\n",
            frac_lost * 100, round(frac_lost * nrow(late_lost)), nrow(late_lost)))

cat(sprintf("Running permutation test (%d iterations)...\n", N_PERM))
set.seed(SEED)
null_fracs <- numeric(N_PERM)
for (i in seq_len(N_PERM)) {
  if (i %% 100 == 0) cat(sprintf("  Iteration %d / %d\n", i, N_PERM))
  sampled_idx <- sample(nrow(late_bg), nrow(late_gained), replace = FALSE)
  sampled_loops <- late_bg[sampled_idx, ]
  null_fracs[i] <- overlap_fraction(sampled_loops, early_wt, ANCHOR_TOLERANCE)
}

z_score <- (frac_gained - mean(null_fracs)) / sd(null_fracs)
emp_p   <- (sum(null_fracs >= frac_gained) + 1) / (N_PERM + 1)

cat(sprintf("\nPermutation results:\n"))
cat(sprintf("  Observed overlap (gained): %.4f\n", frac_gained))
cat(sprintf("  Null mean:   %.4f\n", mean(null_fracs)))
cat(sprintf("  Null SD:     %.4f\n", sd(null_fracs)))
cat(sprintf("  Z-score:     %.3f\n", z_score))
cat(sprintf("  Empirical p: %.4f\n\n", emp_p))

summary_stats$gained_overlap_frac   <- frac_gained
summary_stats$lost_overlap_frac     <- frac_lost
summary_stats$perm_null_mean        <- mean(null_fracs)
summary_stats$perm_z_score          <- z_score
summary_stats$perm_empirical_p      <- emp_p

# -- Bar chart: gained vs lost overlap --
bar_df <- data.frame(
  category = c("Gained in adult KO", "Lost in adult KO", "Null (permutation mean)"),
  fraction = c(frac_gained, frac_lost, mean(null_fracs)),
  fill = c(COLORS["late_KO"], COLORS["late_WT"], "grey60")
)
bar_df$category <- factor(bar_df$category, levels = bar_df$category)

p_bar <- ggplot(bar_df, aes(x = category, y = fraction, fill = category)) +
  geom_col(width = 0.6) +
  geom_errorbar(data = bar_df[bar_df$category == "Null (permutation mean)", ],
                aes(ymin = mean(null_fracs) - sd(null_fracs),
                    ymax = mean(null_fracs) + sd(null_fracs)),
                width = 0.2, linewidth = 0.5) +
  scale_fill_manual(values = c("Gained in adult KO" = unname(COLORS["late_KO"]),
                                "Lost in adult KO" = unname(COLORS["late_WT"]),
                                "Null (permutation mean)" = "grey60")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Overlap of Adult Differential Loops with P13-WT Loops",
       subtitle = sprintf("Permutation z=%.2f, p=%.3f (n=%d iterations)",
                          z_score, emp_p, N_PERM),
       y = "Fraction overlapping P13-WT", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

save_multiformat_ggplot(p_bar,
                        file.path(OUTPUT_DIR, "03_loop_overlap_gained/gained_vs_early_wt_overlap"),
                        width = 8, height = 6)

# -- Permutation null distribution --
p_null <- ggplot(data.frame(frac = null_fracs), aes(x = frac)) +
  geom_histogram(bins = 40, fill = "grey70", color = "grey40") +
  geom_vline(xintercept = frac_gained, color = unname(COLORS["late_KO"]),
             linewidth = 1.2, linetype = "solid") +
  geom_vline(xintercept = frac_lost, color = unname(COLORS["late_WT"]),
             linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = frac_gained, y = Inf, label = sprintf("Gained (%.1f%%)", frac_gained * 100),
           vjust = 2, hjust = 1.1, color = unname(COLORS["late_KO"]), fontface = "bold", size = 4) +
  annotate("text", x = frac_lost, y = Inf, label = sprintf("Lost (%.1f%%)", frac_lost * 100),
           vjust = 4, hjust = -0.1, color = unname(COLORS["late_WT"]), fontface = "bold", size = 4) +
  scale_x_continuous(labels = percent_format()) +
  labs(title = "Permutation Null Distribution",
       subtitle = sprintf("Null: random %d loops from late-unchanged; z=%.2f",
                          nrow(late_gained), z_score),
       x = "Fraction overlapping P13-WT", y = "Count") +
  theme_minimal(base_size = 14)

save_multiformat_ggplot(p_null,
                        file.path(OUTPUT_DIR, "03_loop_overlap_gained/permutation_null_distribution"),
                        width = 8, height = 6)

write.table(
  data.frame(metric = c("observed_gained_frac", "observed_lost_frac",
                         "null_mean", "null_sd", "z_score", "empirical_p",
                         "n_gained", "n_lost", "n_perm", "tolerance_bp"),
             value = c(frac_gained, frac_lost, mean(null_fracs), sd(null_fracs),
                       z_score, emp_p, nrow(late_gained), nrow(late_lost),
                       N_PERM, ANCHOR_TOLERANCE)),
  file.path(OUTPUT_DIR, "03_loop_overlap_gained/permutation_stats.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Module 3 complete.\n\n")


# =============================================================================
# MODULE 4: LOOP DISTANCE eCDF OVERLAY
# =============================================================================

cat("=== Module 4: Loop Distance eCDF Overlay ===\n\n")

dist_data <- do.call(rbind, lapply(names(loop_sets), function(nm) {
  df <- loop_sets[[nm]]
  data.frame(
    group = nm,
    loop_distance = compute_loop_distance(df),
    stringsAsFactors = FALSE
  )
}))

dist_data$group <- factor(dist_data$group, levels = names(COLORS))

medians <- dist_data %>%
  group_by(group) %>%
  summarise(median_dist = median(loop_distance), .groups = "drop")

cat("Median loop distances:\n")
for (i in seq_len(nrow(medians))) {
  cat(sprintf("  %s: %.0f kb\n", medians$group[i], medians$median_dist[i] / 1000))
}

ks_tests <- list(
  lateKO_vs_earlyWT = ks.test(
    dist_data$loop_distance[dist_data$group == "late_KO"],
    dist_data$loop_distance[dist_data$group == "early_WT"]
  ),
  lateKO_vs_lateWT = ks.test(
    dist_data$loop_distance[dist_data$group == "late_KO"],
    dist_data$loop_distance[dist_data$group == "late_WT"]
  ),
  lateWT_vs_earlyWT = ks.test(
    dist_data$loop_distance[dist_data$group == "late_WT"],
    dist_data$loop_distance[dist_data$group == "early_WT"]
  )
)

cat("\nKS tests:\n")
for (nm in names(ks_tests)) {
  cat(sprintf("  %s: D=%.4f, p=%.2e\n", nm,
              ks_tests[[nm]]$statistic, ks_tests[[nm]]$p.value))
}

summary_stats$dist_D_lateKO_earlyWT <- ks_tests$lateKO_vs_earlyWT$statistic
summary_stats$dist_p_lateKO_earlyWT <- ks_tests$lateKO_vs_earlyWT$p.value
summary_stats$dist_D_lateKO_lateWT  <- ks_tests$lateKO_vs_lateWT$statistic
summary_stats$dist_p_lateKO_lateWT  <- ks_tests$lateKO_vs_lateWT$p.value
summary_stats$dist_D_lateWT_earlyWT <- ks_tests$lateWT_vs_earlyWT$statistic
summary_stats$dist_p_lateWT_earlyWT <- ks_tests$lateWT_vs_earlyWT$p.value

# -- eCDF plot --
p_ecdf <- ggplot(dist_data, aes(x = loop_distance, color = group)) +
  stat_ecdf(linewidth = 0.8) +
  scale_x_log10(labels = label_number(scale = 1e-6, suffix = " Mb"),
                breaks = c(1e5, 3e5, 1e6, 3e6, 1e7, 3e7)) +
  scale_color_manual(values = COLORS,
                     labels = c("Early WT (P13)", "Early KO (P13)",
                                "Late WT (Adult)", "Late KO (Adult)")) +
  geom_vline(data = medians, aes(xintercept = median_dist, color = group),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  labs(title = "Loop Distance eCDF: 4-Group Developmental Comparison",
       subtitle = sprintf("KS(late_KO vs early_WT): D=%.4f, p=%.1e | KS(late_KO vs late_WT): D=%.4f, p=%.1e",
                          ks_tests$lateKO_vs_earlyWT$statistic,
                          ks_tests$lateKO_vs_earlyWT$p.value,
                          ks_tests$lateKO_vs_lateWT$statistic,
                          ks_tests$lateKO_vs_lateWT$p.value),
       x = "Loop Distance", y = "Cumulative Fraction",
       color = "Condition") +
  theme_minimal(base_size = 14) +
  theme(legend.position = c(0.75, 0.3))

save_multiformat_ggplot(p_ecdf,
                        file.path(OUTPUT_DIR, "04_distance_ecdf/distance_ecdf_4groups"),
                        width = 10, height = 7)

ks_df <- data.frame(
  comparison = names(ks_tests),
  D_statistic = sapply(ks_tests, function(x) x$statistic),
  p_value = sapply(ks_tests, function(x) x$p.value),
  row.names = NULL
)
write.table(ks_df, file.path(OUTPUT_DIR, "04_distance_ecdf/distance_ks_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nModule 4 complete.\n\n")


# =============================================================================
# MODULE 5: TAD BOUNDARY DIRECTION CONCORDANCE
# =============================================================================

cat("=== Module 5: TAD Boundary Direction Concordance ===\n\n")

early_tad <- read.delim(INPUT_FILES$early_tad, stringsAsFactors = FALSE)
late_tad  <- read.delim(INPUT_FILES$late_tad, stringsAsFactors = FALSE)

early_diff <- early_tad %>%
  filter(Differential %in% c("Differential", "Shifted")) %>%
  filter(chr %in% AUTOSOMES)
late_diff <- late_tad %>%
  filter(Differential %in% c("Differential", "Shifted")) %>%
  filter(chr %in% AUTOSOMES)

cat(sprintf("  Early differential boundaries: %d\n", nrow(early_diff)))
cat(sprintf("  Late differential boundaries:  %d\n", nrow(late_diff)))

early_gr <- GRanges(seqnames = early_diff$chr,
                    ranges = IRanges(start = early_diff$Boundary,
                                     end = early_diff$Boundary))
late_gr <- GRanges(seqnames = late_diff$chr,
                   ranges = IRanges(start = late_diff$Boundary,
                                    end = late_diff$Boundary))

hits <- findOverlaps(early_gr, late_gr, maxgap = TAD_TOLERANCE)
cat(sprintf("  Shared differential boundaries (within %dkb): %d\n\n",
            TAD_TOLERANCE / 1000, length(hits)))

if (length(hits) > 0) {
  shared <- data.frame(
    early_enriched = early_diff$Enriched_In[queryHits(hits)],
    late_enriched  = late_diff$Enriched_In[subjectHits(hits)]
  )

  shared$early_dir <- ifelse(shared$early_enriched == "Matrix 1", "Ctrl-enriched", "Mut-enriched")
  shared$late_dir  <- ifelse(shared$late_enriched == "Matrix 1", "Ctrl-enriched", "Mut-enriched")

  cont_table <- table(Early = shared$early_dir, Late = shared$late_dir)
  fisher_res <- fisher.test(cont_table)

  cat("Contingency table (Early direction x Late direction):\n")
  print(cont_table)
  cat(sprintf("\nFisher's exact test: OR=%.3f, p=%.4e\n", fisher_res$estimate, fisher_res$p.value))

  summary_stats$tad_shared_n       <- length(hits)
  summary_stats$tad_fisher_OR      <- fisher_res$estimate
  summary_stats$tad_fisher_p       <- fisher_res$p.value

  # Concordance: late-KO-enriched that were early-ctrl-enriched
  n_lateKO_earlyCtrl <- cont_table["Ctrl-enriched", "Mut-enriched"]
  n_lateKO_total     <- sum(cont_table[, "Mut-enriched"])
  concordance_frac   <- n_lateKO_earlyCtrl / n_lateKO_total

  cat(sprintf("\nOf %d late-KO-enriched shared boundaries, %d (%.1f%%) were early-ctrl-enriched\n",
              n_lateKO_total, n_lateKO_earlyCtrl, concordance_frac * 100))
  cat("  (If blocking: expect >50%% — boundaries strong in WT at P13 become\n")
  cat("   KO-enriched in adult, meaning KO retained the immature state)\n")

  summary_stats$tad_concordance_frac <- concordance_frac

  # -- Contingency heatmap --
  cont_mat <- as.matrix(cont_table)
  save_multiformat_pheatmap(
    quote(pheatmap(cont_mat,
                   display_numbers = TRUE, number_format = "%.0f", fontsize_number = 16,
                   color = colorRampPalette(c("white", "#fc8d62"))(50),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   main = sprintf("TAD Boundary Direction: Early vs Late\n(Fisher OR=%.2f, p=%.2e)",
                                  fisher_res$estimate, fisher_res$p.value))),
    file.path(OUTPUT_DIR, "05_tad_boundary_concordance/boundary_direction_heatmap"),
    width = 7, height = 6
  )

  # -- Stacked bar of concordance --
  bar_data <- shared %>%
    count(early_dir, late_dir) %>%
    group_by(late_dir) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_stack <- ggplot(bar_data, aes(x = late_dir, y = prop, fill = early_dir)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, prop * 100)),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("Ctrl-enriched" = "#4393c3", "Mut-enriched" = "#d6604d")) +
    scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Shared Differential Boundaries: Direction Concordance",
         subtitle = sprintf("n=%d shared boundaries, Fisher OR=%.2f, p=%.2e",
                            length(hits), fisher_res$estimate, fisher_res$p.value),
         x = "Late (Adult) Enrichment", y = "Proportion",
         fill = "Early (P13) Enrichment") +
    theme_minimal(base_size = 14)

  save_multiformat_ggplot(p_stack,
                          file.path(OUTPUT_DIR, "05_tad_boundary_concordance/boundary_concordance_bar"),
                          width = 8, height = 6)

  write.table(cont_table,
              file.path(OUTPUT_DIR, "05_tad_boundary_concordance/boundary_contingency_table.tsv"),
              sep = "\t", quote = FALSE)

  writeLines(c(
    "TAD Boundary Direction Concordance",
    "===================================",
    sprintf("Shared differential boundaries: %d", length(hits)),
    sprintf("Tolerance: %d kb", TAD_TOLERANCE / 1000),
    "",
    "Contingency table:",
    capture.output(print(cont_table)),
    "",
    sprintf("Fisher's exact test: OR=%.4f, p=%.4e", fisher_res$estimate, fisher_res$p.value),
    sprintf("Late-KO-enriched that were early-ctrl-enriched: %d / %d (%.1f%%)",
            n_lateKO_earlyCtrl, n_lateKO_total, concordance_frac * 100)
  ), file.path(OUTPUT_DIR, "05_tad_boundary_concordance/fisher_exact_result.txt"))

} else {
  cat("WARNING: No shared differential boundaries found within tolerance.\n")
  summary_stats$tad_shared_n <- 0
}

cat("\nModule 5 complete.\n\n")


# =============================================================================
# MODULE 6: STRIPE FOOTPRINT COMPARISON (SUPPLEMENTARY)
# =============================================================================

cat("=== Module 6: Stripe Footprint Comparison ===\n\n")

early_stripes <- tryCatch(
  read.delim(INPUT_FILES$early_stripes, stringsAsFactors = FALSE),
  error = function(e) { cat("Could not load early stripes:", e$message, "\n"); NULL }
)
late_stripes <- tryCatch(
  read.delim(INPUT_FILES$late_stripes, stringsAsFactors = FALSE),
  error = function(e) { cat("Could not load late stripes:", e$message, "\n"); NULL }
)

if (!is.null(early_stripes) && !is.null(late_stripes)) {
  late_gained_stripes <- late_stripes %>%
    filter(source == "mutant_only" |
           (source == "shared" & pval_mut < 0.05 & pval_ctrl >= 0.05))

  early_wt_stripes <- early_stripes %>%
    filter(pval_ctrl < 0.1)

  cat(sprintf("  Late-gained stripes: %d\n", nrow(late_gained_stripes)))
  cat(sprintf("  Early-WT stripes (ctrl p<0.1): %d\n", nrow(early_wt_stripes)))

  if (nrow(late_gained_stripes) > 0 && nrow(early_wt_stripes) > 0) {
    late_anchor_gr <- GRanges(seqnames = late_gained_stripes$chr,
                              ranges = IRanges(start = late_gained_stripes$pos1,
                                               end = late_gained_stripes$pos2))
    early_anchor_gr <- GRanges(seqnames = early_wt_stripes$chr,
                               ranges = IRanges(start = early_wt_stripes$pos1,
                                                end = early_wt_stripes$pos2))

    stripe_hits <- findOverlaps(late_anchor_gr, early_anchor_gr, maxgap = STRIPE_TOLERANCE)
    n_overlap   <- length(unique(queryHits(stripe_hits)))
    frac_stripe <- n_overlap / nrow(late_gained_stripes)

    cat(sprintf("\n  Late-gained stripes overlapping early-WT anchor: %d / %d (%.1f%%)\n",
                n_overlap, nrow(late_gained_stripes), frac_stripe * 100))

    summary_stats$stripe_n_late_gained   <- nrow(late_gained_stripes)
    summary_stats$stripe_n_early_wt      <- nrow(early_wt_stripes)
    summary_stats$stripe_overlap_n       <- n_overlap
    summary_stats$stripe_overlap_frac    <- frac_stripe

    overlap_df <- data.frame(
      category = c("Overlap with\nP13-WT stripes", "No overlap"),
      count = c(n_overlap, nrow(late_gained_stripes) - n_overlap)
    )
    overlap_df$frac <- overlap_df$count / sum(overlap_df$count)

    p_stripe <- ggplot(overlap_df, aes(x = "", y = frac, fill = category)) +
      geom_col(width = 0.5) +
      geom_text(aes(label = sprintf("%d (%.0f%%)", count, frac * 100)),
                position = position_stack(vjust = 0.5), size = 5) +
      scale_fill_manual(values = c("Overlap with\nP13-WT stripes" = "#66c2a5",
                                    "No overlap" = "grey80")) +
      coord_flip() +
      labs(title = "Late-Gained Stripes: Overlap with P13-WT Stripe Anchors",
           subtitle = sprintf("n=%d late-gained stripes, %dkb tolerance",
                              nrow(late_gained_stripes), STRIPE_TOLERANCE / 1000),
           x = NULL, y = "Fraction", fill = NULL) +
      theme_minimal(base_size = 14) +
      theme(axis.text.y = element_blank())

    save_multiformat_ggplot(p_stripe,
                            file.path(OUTPUT_DIR, "06_stripe_overlap/stripe_anchor_overlap"),
                            width = 9, height = 4)
  } else {
    cat("Insufficient stripes for comparison.\n")
  }
} else {
  cat("Stripe data unavailable, skipping Module 6.\n")
}

cat("\nModule 6 complete.\n\n")


# =============================================================================
# MASTER SUMMARY TABLE
# =============================================================================

cat("=== Writing master summary ===\n\n")

stats_df <- data.frame(
  metric = names(summary_stats),
  value = sapply(summary_stats, function(x) sprintf("%.6g", x)),
  stringsAsFactors = FALSE,
  row.names = NULL
)

write.table(stats_df,
            file.path(OUTPUT_DIR, "tables/developmental_comparison_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -- Human-readable report --
report_lines <- c(
  "Developmental Comparison: Adult BAP1-KO vs P13 Wildtype",
  "========================================================",
  sprintf("Generated: %s", Sys.time()),
  "",
  "Question: Does adult BAP1-KO 3D genome resemble P13 wildtype,",
  "suggesting blocked developmental chromatin remodeling?",
  "",
  "--- PC1 Compartment Correlation ---",
  sprintf("r(late_KO, early_WT)  = %.4f", summary_stats$pc1_r_lateKO_earlyWT),
  sprintf("r(late_WT, early_WT)  = %.4f", summary_stats$pc1_r_lateWT_earlyWT),
  sprintf("r(late_KO, late_WT)   = %.4f", summary_stats$pc1_r_lateKO_lateWT),
  sprintf("Fisher Z: z=%.3f, p=%.2e", summary_stats$pc1_fisher_z, summary_stats$pc1_fisher_p),
  sprintf("Shared bins: %d", summary_stats$pc1_n_bins),
  "",
  "--- Loop Jaccard Similarity ---",
  sprintf("J(late_KO, early_WT)  = %.4f", summary_stats$loop_J_lateKO_earlyWT),
  sprintf("J(late_WT, early_WT)  = %.4f", summary_stats$loop_J_lateWT_earlyWT),
  sprintf("J(late_KO, late_WT)   = %.4f", summary_stats$loop_J_lateKO_lateWT),
  "",
  "--- Gained Loop Overlap with P13-WT ---",
  sprintf("Gained overlap:  %.1f%%", summary_stats$gained_overlap_frac * 100),
  sprintf("Lost overlap:    %.1f%%", summary_stats$lost_overlap_frac * 100),
  sprintf("Null mean:       %.1f%%", summary_stats$perm_null_mean * 100),
  sprintf("Permutation z:   %.3f", summary_stats$perm_z_score),
  sprintf("Permutation p:   %.4f", summary_stats$perm_empirical_p),
  "",
  "--- Loop Distance eCDF ---",
  sprintf("KS(late_KO vs early_WT): D=%.4f, p=%.2e",
          summary_stats$dist_D_lateKO_earlyWT, summary_stats$dist_p_lateKO_earlyWT),
  sprintf("KS(late_KO vs late_WT):  D=%.4f, p=%.2e",
          summary_stats$dist_D_lateKO_lateWT, summary_stats$dist_p_lateKO_lateWT),
  "",
  "--- TAD Boundary Direction Concordance ---"
)

if (!is.null(summary_stats$tad_shared_n) && summary_stats$tad_shared_n > 0) {
  report_lines <- c(report_lines,
    sprintf("Shared differential boundaries: %d", summary_stats$tad_shared_n),
    sprintf("Fisher OR: %.3f, p=%.2e", summary_stats$tad_fisher_OR, summary_stats$tad_fisher_p),
    sprintf("Late-KO-enriched from early-ctrl-enriched: %.1f%%",
            summary_stats$tad_concordance_frac * 100)
  )
} else {
  report_lines <- c(report_lines, "No shared differential boundaries found.")
}

if (!is.null(summary_stats$stripe_overlap_frac)) {
  report_lines <- c(report_lines,
    "",
    "--- Stripe Footprint (Supplementary) ---",
    sprintf("Late-gained stripes overlapping early-WT: %d / %d (%.1f%%)",
            summary_stats$stripe_overlap_n, summary_stats$stripe_n_late_gained,
            summary_stats$stripe_overlap_frac * 100)
  )
}

report_lines <- c(report_lines, "",
  "--- Interpretation ---",
  if (!is.null(summary_stats$pc1_fisher_p) && summary_stats$pc1_fisher_p < 0.05 &&
      summary_stats$pc1_fisher_z > 0) {
    "PC1: Adult KO is significantly more similar to P13-WT than Adult-WT is. SUPPORTS blocking."
  } else {
    "PC1: No significant shift toward P13-WT in Adult KO compartments."
  },
  if (!is.null(summary_stats$perm_z_score) && summary_stats$perm_z_score > 1.96) {
    "Loops: Gained adult-KO loops significantly overlap P13-WT positions. SUPPORTS blocking."
  } else {
    "Loops: Gained adult-KO loops do not specifically overlap P13-WT positions."
  },
  if (!is.null(summary_stats$tad_fisher_p) && summary_stats$tad_fisher_p < 0.05 &&
      !is.null(summary_stats$tad_concordance_frac) && summary_stats$tad_concordance_frac > 0.5) {
    "TADs: Boundary direction concordance supports a developmental blocking model."
  } else {
    "TADs: Boundary directions do not show clear developmental blocking pattern."
  }
)

writeLines(report_lines, file.path(OUTPUT_DIR, "developmental_comparison_report.txt"))

cat("Report written to: ", file.path(OUTPUT_DIR, "developmental_comparison_report.txt"), "\n")
cat("\n================================================\n")
cat("Developmental Comparison Analysis Complete\n")
cat("================================================\n\n")
