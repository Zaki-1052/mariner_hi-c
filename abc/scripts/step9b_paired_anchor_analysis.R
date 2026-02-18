# abc/scripts/step9b_paired_anchor_analysis.R
#
# Paired-anchor loop-ABC analysis: tests whether ABC-predicted enhancer-gene
# pairs are physically connected by the SAME differential Hi-C loop (enhancer
# at one anchor, gene TSS at the other).

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# === CONFIGURATION ===
LOOPS_FILE     <- "../characterized_loops.tsv"
ABC_FILE       <- "results/delta_abc_all_pairs.tsv"
TSS_FILE       <- "reference/mm10_tss.bed"
OUT_DIR        <- "results"
PLOT_DIR       <- "results/paired_anchor_plots"
ABC_DELTA_THRESH <- 0.01

# === VALIDATE INPUTS ===
stopifnot(file.exists(LOOPS_FILE))
stopifnot(file.exists(ABC_FILE))
stopifnot(file.exists(TSS_FILE))
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Step 9b: Paired-Anchor Loop-ABC Analysis ===\n\n")

# === STEP 1: LOAD DATA ===
cat("Loading data...\n")

loops <- fread(LOOPS_FILE)
cat(sprintf("  Loops: %d rows\n", nrow(loops)))

abc <- fread(ABC_FILE)
cat(sprintf("  ABC pairs: %d rows\n", nrow(abc)))

tss <- fread(TSS_FILE)
setnames(tss, c("chr", "start", "end", "name", "score", "strand", "ensembl_id", "gene_type"))
cat(sprintf("  TSS entries: %d rows\n", nrow(tss)))

# Merge TSS coordinates onto ABC pairs by gene name
abc_genes <- unique(abc$TargetGene)
tss_lookup <- tss[name %in% abc_genes, .(tss_chr = chr, tss_start = start, tss_end = end, name)]
missing_genes <- setdiff(abc_genes, tss_lookup$name)
if (length(missing_genes) > 0) {
  stop(sprintf("FATAL: %d ABC TargetGenes have no TSS entry: %s",
               length(missing_genes),
               paste(head(missing_genes, 10), collapse = ", ")))
}
cat(sprintf("  All %d ABC target genes matched to TSS\n", length(abc_genes)))

abc <- merge(abc, tss_lookup, by.x = "TargetGene", by.y = "name", all.x = TRUE)
stopifnot(!any(is.na(abc$tss_chr)))

# === STEP 2: BUILD GRanges OBJECTS ===
cat("\nBuilding GRanges objects...\n")

anchor1_gr <- GRanges(seqnames = loops$chr1,
                      ranges = IRanges(start = loops$start1, end = loops$end1))
anchor2_gr <- GRanges(seqnames = loops$chr2,
                      ranges = IRanges(start = loops$start2, end = loops$end2))

enh_gr <- GRanges(seqnames = abc$chr,
                  ranges = IRanges(start = abc$start, end = abc$end))
tss_gr <- GRanges(seqnames = abc$tss_chr,
                  ranges = IRanges(start = abc$tss_start, end = abc$tss_end))

cat(sprintf("  Anchor1/2: %d ranges each\n", length(anchor1_gr)))
cat(sprintf("  Enhancer GRanges: %d ranges\n", length(enh_gr)))
cat(sprintf("  TSS GRanges: %d ranges\n", length(tss_gr)))

# === STEP 3: CASE A — enhancer overlaps anchor1, TSS overlaps anchor2 ===
cat("\nFinding paired-anchor overlaps (Case A: enh→anchor1, TSS→anchor2)...\n")

hits_enh_a1 <- findOverlaps(enh_gr, anchor1_gr)
abc_idx_A <- queryHits(hits_enh_a1)
loop_idx_A <- subjectHits(hits_enh_a1)

# Pairwise check: does this ABC row's TSS overlap the same loop's anchor2?
pairwise_A <- (as.character(seqnames(tss_gr[abc_idx_A])) == as.character(seqnames(anchor2_gr[loop_idx_A]))) &
              (start(tss_gr[abc_idx_A]) <= end(anchor2_gr[loop_idx_A])) &
              (end(tss_gr[abc_idx_A]) >= start(anchor2_gr[loop_idx_A]))

matched_A <- data.table(abc_row = abc_idx_A[pairwise_A],
                        loop_row = loop_idx_A[pairwise_A],
                        case = "A")
cat(sprintf("  Case A: %d initial overlaps, %d paired matches\n",
            length(abc_idx_A), sum(pairwise_A)))

# === STEP 4: CASE B — enhancer overlaps anchor2, TSS overlaps anchor1 ===
cat("Finding paired-anchor overlaps (Case B: enh→anchor2, TSS→anchor1)...\n")

hits_enh_a2 <- findOverlaps(enh_gr, anchor2_gr)
abc_idx_B <- queryHits(hits_enh_a2)
loop_idx_B <- subjectHits(hits_enh_a2)

pairwise_B <- (as.character(seqnames(tss_gr[abc_idx_B])) == as.character(seqnames(anchor1_gr[loop_idx_B]))) &
              (start(tss_gr[abc_idx_B]) <= end(anchor1_gr[loop_idx_B])) &
              (end(tss_gr[abc_idx_B]) >= start(anchor1_gr[loop_idx_B]))

matched_B <- data.table(abc_row = abc_idx_B[pairwise_B],
                        loop_row = loop_idx_B[pairwise_B],
                        case = "B")
cat(sprintf("  Case B: %d initial overlaps, %d paired matches\n",
            length(abc_idx_B), sum(pairwise_B)))

# === STEP 5: COMBINE AND DEDUPLICATE ===
cat("\nCombining Cases A+B...\n")

matched_all <- rbind(matched_A, matched_B)
# Deduplicate by unique (abc_row, loop_row) — same pair found via both cases
matched_all <- unique(matched_all, by = c("abc_row", "loop_row"))
cat(sprintf("  Total unique (ABC pair, loop) matches: %d\n", nrow(matched_all)))

if (nrow(matched_all) == 0) {
  stop("FATAL: No paired-anchor matches found. Check coordinate systems.")
}

# Attach metadata from both tables
result <- data.table(
  # Loop info
  loop_id        = loops$loop_id[matched_all$loop_row],
  loop_chr1      = loops$chr1[matched_all$loop_row],
  loop_start1    = loops$start1[matched_all$loop_row],
  loop_end1      = loops$end1[matched_all$loop_row],
  loop_chr2      = loops$chr2[matched_all$loop_row],
  loop_start2    = loops$start2[matched_all$loop_row],
  loop_end2      = loops$end2[matched_all$loop_row],
  loop_logFC     = loops$logFC[matched_all$loop_row],
  loop_FDR       = loops$FDR[matched_all$loop_row],
  loop_direction = loops$direction[matched_all$loop_row],
  loop_type      = loops$loop_type[matched_all$loop_row],
  loop_significant = loops$significant[matched_all$loop_row],
  # ABC info
  enh_chr        = abc$chr[matched_all$abc_row],
  enh_start      = abc$start[matched_all$abc_row],
  enh_end        = abc$end[matched_all$abc_row],
  target_gene    = abc$TargetGene[matched_all$abc_row],
  tss_chr        = abc$tss_chr[matched_all$abc_row],
  tss_start      = abc$tss_start[matched_all$abc_row],
  tss_end        = abc$tss_end[matched_all$abc_row],
  abc_score_wt   = abc$ABC.Score_WT[matched_all$abc_row],
  abc_score_ko   = abc$ABC.Score_KO[matched_all$abc_row],
  delta_ABC      = abc$delta_ABC[matched_all$abc_row],
  delta_unnorm   = abc$delta_unnorm[matched_all$abc_row],
  distance       = abc$distance[matched_all$abc_row],
  match_case     = matched_all$case
)

# === STEP 6: ANALYSES ===

# --- 6a: Match statistics ---
cat("\n=== MATCH STATISTICS ===\n")
n_triplets     <- nrow(result)
n_unique_loops <- length(unique(result$loop_id))
n_unique_pairs <- nrow(unique(result[, .(enh_chr, enh_start, enh_end, target_gene)]))
n_unique_genes <- length(unique(result$target_gene))
n_unique_enh   <- nrow(unique(result[, .(enh_chr, enh_start, enh_end)]))

cat(sprintf("  Matched triplets (loop, enhancer, gene): %d\n", n_triplets))
cat(sprintf("  Unique loops with >= 1 match: %d / %d (%.1f%%)\n",
            n_unique_loops, nrow(loops), 100 * n_unique_loops / nrow(loops)))
cat(sprintf("  Unique E-G pairs confirmed by loops: %d / %d (%.2f%%)\n",
            n_unique_pairs, nrow(abc), 100 * n_unique_pairs / nrow(abc)))
cat(sprintf("  Unique genes in matched set: %d\n", n_unique_genes))
cat(sprintf("  Unique enhancers in matched set: %d\n", n_unique_enh))

# --- 6b: Directional concordance ---
cat("\n=== DIRECTIONAL CONCORDANCE ===\n")

# Define concordance: loop up + ΔABC positive = concordant
#                     loop down + ΔABC negative = concordant
# Filter to pairs with meaningful delta
result[, abc_direction := fifelse(delta_ABC > ABC_DELTA_THRESH, "gained",
                           fifelse(delta_ABC < -ABC_DELTA_THRESH, "lost", "unchanged"))]
result[, unnorm_direction := fifelse(delta_unnorm > 0, "gained",
                               fifelse(delta_unnorm < 0, "lost", "unchanged"))]

# ΔABC concordance (filtered to directional pairs)
directional_abc <- result[abc_direction != "unchanged"]
if (nrow(directional_abc) > 0) {
  directional_abc[, concordant := (loop_direction == "up_in_mutant" & abc_direction == "gained") |
                                  (loop_direction == "down_in_mutant" & abc_direction == "lost")]
  n_conc_abc <- sum(directional_abc$concordant)
  n_dir_abc  <- nrow(directional_abc)
  pct_abc    <- 100 * n_conc_abc / n_dir_abc
  binom_abc  <- binom.test(n_conc_abc, n_dir_abc, p = 0.5)

  cat(sprintf("  ΔABC concordance: %d / %d = %.1f%% (p = %.2e, 95%% CI: %.1f-%.1f%%)\n",
              n_conc_abc, n_dir_abc, pct_abc,
              binom_abc$p.value,
              100 * binom_abc$conf.int[1], 100 * binom_abc$conf.int[2]))
} else {
  cat("  No directional ABC pairs found (all |ΔABC| <= threshold)\n")
}

# Δ(A×C) unnormalized concordance (any nonzero)
directional_unnorm <- result[unnorm_direction != "unchanged"]
if (nrow(directional_unnorm) > 0) {
  directional_unnorm[, concordant := (loop_direction == "up_in_mutant" & unnorm_direction == "gained") |
                                     (loop_direction == "down_in_mutant" & unnorm_direction == "lost")]
  n_conc_un <- sum(directional_unnorm$concordant)
  n_dir_un  <- nrow(directional_unnorm)
  pct_un    <- 100 * n_conc_un / n_dir_un
  binom_un  <- binom.test(n_conc_un, n_dir_un, p = 0.5)

  cat(sprintf("  Δ(A×C) concordance: %d / %d = %.1f%% (p = %.2e, 95%% CI: %.1f-%.1f%%)\n",
              n_conc_un, n_dir_un, pct_un,
              binom_un$p.value,
              100 * binom_un$conf.int[1], 100 * binom_un$conf.int[2]))
}

cat(sprintf("  Reference: Step 9 independent-overlap concordance was 51.4%%\n"))

# --- 6c: Stratification by loop type ---
cat("\n=== STRATIFICATION BY LOOP TYPE ===\n")

loop_type_counts <- result[, .N, by = loop_type][order(-N)]
cat("  Loop type distribution in matched set:\n")
for (i in seq_len(nrow(loop_type_counts))) {
  cat(sprintf("    %s: %d (%.1f%%)\n",
              loop_type_counts$loop_type[i],
              loop_type_counts$N[i],
              100 * loop_type_counts$N[i] / nrow(result)))
}

# Concordance by loop type (using ΔABC)
if (nrow(directional_abc) > 0) {
  cat("\n  Concordance by loop type (ΔABC):\n")
  conc_by_type <- directional_abc[, .(
    concordant = sum(concordant),
    total = .N,
    pct = 100 * sum(concordant) / .N
  ), by = loop_type][order(-total)]

  for (i in seq_len(nrow(conc_by_type))) {
    cat(sprintf("    %s: %d/%d = %.1f%%\n",
                conc_by_type$loop_type[i],
                conc_by_type$concordant[i],
                conc_by_type$total[i],
                conc_by_type$pct[i]))
  }
}

# --- 6d: Effect size comparison ---
cat("\n=== EFFECT SIZE COMPARISON ===\n")

matched_delta <- abs(result$delta_ABC)
all_delta     <- abs(abc$delta_ABC)
wt <- wilcox.test(matched_delta, all_delta, alternative = "greater")
cat(sprintf("  Median |ΔABC| matched: %.4f vs all: %.4f\n",
            median(matched_delta), median(all_delta)))
cat(sprintf("  Wilcoxon p-value (matched > all): %.2e\n", wt$p.value))

matched_unnorm <- abs(result$delta_unnorm)
all_unnorm     <- abs(abc$delta_unnorm)
wt2 <- wilcox.test(matched_unnorm, all_unnorm, alternative = "greater")
cat(sprintf("  Median |Δ(A×C)| matched: %.6f vs all: %.6f\n",
            median(matched_unnorm), median(all_unnorm)))
cat(sprintf("  Wilcoxon p-value (matched > all): %.2e\n", wt2$p.value))

# === STEP 7: SAVE RESULTS ===
cat("\nSaving results...\n")

# Main results table
fwrite(result, file.path(OUT_DIR, "paired_anchor_matches.tsv"), sep = "\t")
cat(sprintf("  Wrote %s (%d rows)\n",
            file.path(OUT_DIR, "paired_anchor_matches.tsv"), nrow(result)))

# Summary text
summary_lines <- c(
  "=== Paired-Anchor Loop-ABC Analysis Summary ===",
  sprintf("Date: %s", Sys.time()),
  "",
  "--- Match Statistics ---",
  sprintf("Matched (loop, enhancer, gene) triplets: %d", n_triplets),
  sprintf("Unique loops with >= 1 match: %d / %d (%.1f%%)",
          n_unique_loops, nrow(loops), 100 * n_unique_loops / nrow(loops)),
  sprintf("Unique E-G pairs confirmed by loops: %d / %d (%.4f%%)",
          n_unique_pairs, nrow(abc), 100 * n_unique_pairs / nrow(abc)),
  sprintf("Unique genes: %d", n_unique_genes),
  sprintf("Unique enhancers: %d", n_unique_enh),
  ""
)

if (nrow(directional_abc) > 0) {
  summary_lines <- c(summary_lines,
    "--- Directional Concordance ---",
    sprintf("ΔABC concordance: %d / %d = %.1f%% (p = %.2e)",
            n_conc_abc, n_dir_abc, pct_abc, binom_abc$p.value),
    sprintf("Δ(A×C) concordance: %d / %d = %.1f%% (p = %.2e)",
            n_conc_un, n_dir_un, pct_un, binom_un$p.value),
    "Step 9 independent-overlap concordance: 51.4%",
    ""
  )
}

summary_lines <- c(summary_lines,
  "--- Effect Size ---",
  sprintf("Median |ΔABC| matched: %.4f vs all: %.4f (Wilcoxon p = %.2e)",
          median(matched_delta), median(all_delta), wt$p.value),
  sprintf("Median |Δ(A×C)| matched: %.6f vs all: %.6f (Wilcoxon p = %.2e)",
          median(matched_unnorm), median(all_unnorm), wt2$p.value)
)

writeLines(summary_lines, file.path(OUT_DIR, "paired_anchor_summary.txt"))
cat(sprintf("  Wrote %s\n", file.path(OUT_DIR, "paired_anchor_summary.txt")))

# === STEP 8: PLOTS ===
cat("\nGenerating plots...\n")

# --- Plot 1: Concordance comparison bar chart ---
if (nrow(directional_abc) > 0) {
  conc_df <- data.frame(
    method = c("Independent\noverlap (Step 9)", "Paired-anchor\n(dABC)", "Paired-anchor\n(d(AxC))"),
    concordance = c(51.4, pct_abc, pct_un),
    n = c(NA, n_dir_abc, n_dir_un)
  )
  conc_df$method <- factor(conc_df$method, levels = conc_df$method)

  p1 <- ggplot(conc_df, aes(x = method, y = concordance, fill = method)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = sprintf("%.1f%%", concordance)), vjust = -0.5, size = 4) +
    geom_text(aes(label = ifelse(is.na(n), "", sprintf("n=%d", n))),
              vjust = 1.5, color = "white", size = 3.5) +
    scale_fill_manual(values = c("grey60", "#2166AC", "#B2182B")) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(title = "Loop-ABC Directional Concordance",
         subtitle = "Paired-anchor vs independent overlap",
         x = NULL, y = "Concordance (%)") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank())

  ggsave(file.path(PLOT_DIR, "paired_anchor_concordance.pdf"), p1,
         width = 5, height = 5)
  cat("  Saved paired_anchor_concordance.pdf\n")
}

# --- Plot 2: Loop logFC vs ΔABC scatter ---
p2 <- ggplot(result, aes(x = loop_logFC, y = delta_ABC)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 1.5, color = "#4393C3") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  labs(title = "Loop logFC vs dABC (Paired-Anchor Matches)",
       subtitle = sprintf("n = %d matched triplets", nrow(result)),
       x = "Loop logFC (mutant / control)",
       y = "dABC (KO - WT)") +
  theme_minimal(base_size = 13)

# Add Spearman correlation
cor_test <- cor.test(result$loop_logFC, result$delta_ABC, method = "spearman")
p2 <- p2 + annotate("text", x = Inf, y = Inf,
                     label = sprintf("rho = %.3f\np = %.2e",
                                     cor_test$estimate, cor_test$p.value),
                     hjust = 1.1, vjust = 1.5, size = 4)

ggsave(file.path(PLOT_DIR, "logFC_vs_deltaABC.pdf"), p2,
       width = 6, height = 5)
cat("  Saved logFC_vs_deltaABC.pdf\n")

# --- Plot 3: Loop logFC vs Δ(A×C) scatter ---
p3 <- ggplot(result, aes(x = loop_logFC, y = delta_unnorm)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 1.5, color = "#D6604D") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  labs(title = "Loop logFC vs d(AxC) (Paired-Anchor Matches)",
       subtitle = sprintf("n = %d matched triplets", nrow(result)),
       x = "Loop logFC (mutant / control)",
       y = "d(Activity x Contact) (KO - WT)") +
  theme_minimal(base_size = 13)

cor_test2 <- cor.test(result$loop_logFC, result$delta_unnorm, method = "spearman")
p3 <- p3 + annotate("text", x = Inf, y = Inf,
                     label = sprintf("rho = %.3f\np = %.2e",
                                     cor_test2$estimate, cor_test2$p.value),
                     hjust = 1.1, vjust = 1.5, size = 4)

ggsave(file.path(PLOT_DIR, "logFC_vs_delta_unnorm.pdf"), p3,
       width = 6, height = 5)
cat("  Saved logFC_vs_delta_unnorm.pdf\n")

# --- Plot 4: Combined panel ---
if (nrow(directional_abc) > 0) {
  combined <- (p1 | p2) / (p3 | plot_spacer()) +
    plot_annotation(title = "Paired-Anchor Loop-ABC Analysis",
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))

  ggsave(file.path(PLOT_DIR, "paired_anchor_panel.pdf"), combined,
         width = 11, height = 10)
  cat("  Saved paired_anchor_panel.pdf\n")
}

cat("\n=== Step 9b complete ===\n")
