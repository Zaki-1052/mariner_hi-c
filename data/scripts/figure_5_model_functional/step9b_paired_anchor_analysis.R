# abc/scripts/step9b_paired_anchor_analysis.R
#
# Paired-anchor loop-ABC analysis (v2): tests whether ABC-predicted enhancer-gene
# pairs are physically connected by the SAME Hi-C loop (enhancer at one anchor,
# gene TSS at the other). Runs overlap on ALL 39K loops for visual context,
# then performs statistical analyses on the 2,910 differential loops only.
# Enhancements: FDR stratification, RNA-seq 3-way concordance, K119ub overlay.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

source("data/scripts/_shared/multi_format_output.R") # Original: source("../scripts/utils/multi_format_output.R")

# === CONFIGURATION ===
ALL_LOOPS_FILE  <- "../outputs//250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe" # TODO: not in data/
DIFF_LOOPS_FILE <- "data/upstream/loop_calls/late_characterized_loops.tsv" # Original: ../outputs/250402-late_outputs/merged_loops/non_redundant_loops.tsv
ABC_FILE        <- "data/tsvs/figure_4_abc_analysis/4A_delta_abc_all_pairs.tsv" # Original: results/delta_abc_all_pairs.tsv
TSS_FILE        <- "reference/mm10_tss.bed" # TODO: not in data/
RNASEQ_FILE     <- "data/tsvs/figure_5_model_functional/5B_gene_level_summary.tsv" # Original: results/gene_level_summary.tsv
K119UB_FILE     <- "data/tsvs/figure_4_abc_analysis/4F_k119ub_abc_enhancer_merged.tsv" # Original: results/k119ub_abc_enhancer_merged.tsv
OUT_TSV_DIR     <- "data/tsvs/figure_5_model_functional" # Original: results
PLOT_DIR        <- "data/plots/figure_5_model_functional" # Original: results/paired_anchor_plots
SUPP_PLOT_DIR   <- "data/plots/supplemental" # Original: (new — for paired_anchor_panel)
FIG4_PLOT_DIR   <- "data/plots/figure_4_abc_analysis" # Original: (new — for 4E logFC_vs_deltaABC)
ABC_DELTA_THRESH <- 0.01

# === VALIDATE INPUTS ===
for (f in c(ALL_LOOPS_FILE, DIFF_LOOPS_FILE, ABC_FILE, TSS_FILE, RNASEQ_FILE, K119UB_FILE)) {
  if (!file.exists(f)) stop(sprintf("FATAL: Required input not found: %s", f))
}
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE) # Original: dir.create("results/paired_anchor_plots", ...)
dir.create(OUT_TSV_DIR, recursive = TRUE, showWarnings = FALSE) # Original: (OUT_DIR was "results" — already existed)
dir.create(SUPP_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG4_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Step 9b: Paired-Anchor Loop-ABC Analysis (v2) ===\n\n")

# === STEP 1: LOAD DATA ===
cat("Loading data...\n")

all_loops <- fread(ALL_LOOPS_FILE)
cat(sprintf("  All loops: %d rows\n", nrow(all_loops)))

diff_loops <- fread(DIFF_LOOPS_FILE)
cat(sprintf("  Differential loops: %d rows\n", nrow(diff_loops)))

abc <- fread(ABC_FILE)
cat(sprintf("  ABC pairs: %d rows\n", nrow(abc)))

tss <- fread(TSS_FILE)
setnames(tss, c("chr", "start", "end", "name", "score", "strand", "ensembl_id", "gene_type"))
cat(sprintf("  TSS entries: %d rows\n", nrow(tss)))

rnaseq <- fread(RNASEQ_FILE)
cat(sprintf("  RNA-seq gene summary: %d rows\n", nrow(rnaseq)))

k119ub <- fread(K119UB_FILE)
cat(sprintf("  K119ub enhancer signal: %d rows\n", nrow(k119ub)))

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

# Build diff_loops lookup key for later subsetting
diff_loops[, loop_key := paste(chr1, start1, end1, chr2, start2, end2, sep = "_")]
diff_lookup <- diff_loops[, .(loop_key, loop_id, loop_type, category,
                               diff_significant = as.logical(significant))]
setkey(diff_lookup, loop_key)
cat(sprintf("  Diff loop keys built: %d unique\n", nrow(diff_lookup)))

# === STEP 2: BUILD GRanges FROM ALL LOOPS ===
cat(sprintf("\nBuilding GRanges from ALL %d loops...\n", nrow(all_loops)))

anchor1_gr <- GRanges(seqnames = all_loops$chr1,
                      ranges = IRanges(start = all_loops$start1, end = all_loops$end1))
anchor2_gr <- GRanges(seqnames = all_loops$chr2,
                      ranges = IRanges(start = all_loops$start2, end = all_loops$end2))

enh_gr <- GRanges(seqnames = abc$chr,
                  ranges = IRanges(start = abc$start, end = abc$end))
tss_gr <- GRanges(seqnames = abc$tss_chr,
                  ranges = IRanges(start = abc$tss_start, end = abc$tss_end))

cat(sprintf("  Anchor1/2: %d ranges each\n", length(anchor1_gr)))
cat(sprintf("  Enhancer GRanges: %d ranges\n", length(enh_gr)))
cat(sprintf("  TSS GRanges: %d ranges\n", length(tss_gr)))

# === STEP 3: CASE A — enhancer overlaps anchor1, TSS overlaps anchor2 ===
cat("\nFinding paired-anchor overlaps (Case A: enh->anchor1, TSS->anchor2)...\n")

hits_enh_a1 <- findOverlaps(enh_gr, anchor1_gr)
abc_idx_A <- queryHits(hits_enh_a1)
loop_idx_A <- subjectHits(hits_enh_a1)

pairwise_A <- (as.character(seqnames(tss_gr[abc_idx_A])) == as.character(seqnames(anchor2_gr[loop_idx_A]))) &
              (start(tss_gr[abc_idx_A]) <= end(anchor2_gr[loop_idx_A])) &
              (end(tss_gr[abc_idx_A]) >= start(anchor2_gr[loop_idx_A]))

matched_A <- data.table(abc_row = abc_idx_A[pairwise_A],
                        loop_row = loop_idx_A[pairwise_A],
                        case = "A")
cat(sprintf("  Case A: %d initial overlaps, %d paired matches\n",
            length(abc_idx_A), sum(pairwise_A)))

# === STEP 4: CASE B — enhancer overlaps anchor2, TSS overlaps anchor1 ===
cat("Finding paired-anchor overlaps (Case B: enh->anchor2, TSS->anchor1)...\n")

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
matched_all <- unique(matched_all, by = c("abc_row", "loop_row"))
cat(sprintf("  Total unique (ABC pair, loop) matches: %d\n", nrow(matched_all)))

if (nrow(matched_all) == 0) {
  stop("FATAL: No paired-anchor matches found. Check coordinate systems.")
}

# Build all_result with metadata from both tables
all_result <- data.table(
  loop_chr1      = all_loops$chr1[matched_all$loop_row],
  loop_start1    = all_loops$start1[matched_all$loop_row],
  loop_end1      = all_loops$end1[matched_all$loop_row],
  loop_chr2      = all_loops$chr2[matched_all$loop_row],
  loop_start2    = all_loops$start2[matched_all$loop_row],
  loop_end2      = all_loops$end2[matched_all$loop_row],
  loop_logFC     = all_loops$logFC[matched_all$loop_row],
  loop_FDR       = all_loops$FDR[matched_all$loop_row],
  loop_direction = all_loops$direction[matched_all$loop_row],
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

# Create loop coordinate key for matching to differential loops
all_result[, loop_key := paste(loop_chr1, loop_start1, loop_end1,
                               loop_chr2, loop_start2, loop_end2, sep = "_")]

# Left-join differential loop metadata
all_result <- merge(all_result, diff_lookup, by = "loop_key", all.x = TRUE)
all_result[, is_differential := !is.na(loop_id)]
# For non-differential loops, derive significance from FDR
all_result[is.na(diff_significant), diff_significant := (loop_FDR < 0.05)]

# Split into differential and background subsets
diff_result <- all_result[is_differential == TRUE]
bg_result   <- all_result[is_differential == FALSE]

cat(sprintf("  All matches: %d (from %d unique loops)\n",
            nrow(all_result), length(unique(all_result$loop_key))))
cat(sprintf("  Differential matches: %d (from %d unique differential loops)\n",
            nrow(diff_result), length(unique(diff_result$loop_id))))
cat(sprintf("  Background matches: %d\n", nrow(bg_result)))

# Use diff_result for all statistical analyses (backward-compatible alias)
result <- diff_result

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
            n_unique_loops, nrow(diff_loops), 100 * n_unique_loops / nrow(diff_loops)))
cat(sprintf("  Unique E-G pairs confirmed by loops: %d / %d (%.2f%%)\n",
            n_unique_pairs, nrow(abc), 100 * n_unique_pairs / nrow(abc)))
cat(sprintf("  Unique genes in matched set: %d\n", n_unique_genes))
cat(sprintf("  Unique enhancers in matched set: %d\n", n_unique_enh))

# --- 6b: Directional concordance ---
cat("\n=== DIRECTIONAL CONCORDANCE ===\n")

result[, abc_direction := fifelse(delta_ABC > ABC_DELTA_THRESH, "gained",
                           fifelse(delta_ABC < -ABC_DELTA_THRESH, "lost", "unchanged"))]
result[, unnorm_direction := fifelse(delta_unnorm > 0, "gained",
                               fifelse(delta_unnorm < 0, "lost", "unchanged"))]

# dABC concordance
directional_abc <- result[abc_direction != "unchanged"]
stopifnot(nrow(directional_abc) > 0)

directional_abc[, concordant := (loop_direction == "up_in_mutant" & abc_direction == "gained") |
                                (loop_direction == "down_in_mutant" & abc_direction == "lost")]
n_conc_abc <- sum(directional_abc$concordant)
n_dir_abc  <- nrow(directional_abc)
pct_abc    <- 100 * n_conc_abc / n_dir_abc
binom_abc  <- binom.test(n_conc_abc, n_dir_abc, p = 0.5)

cat(sprintf("  dABC concordance: %d / %d = %.1f%% (p = %.2e, 95%% CI: %.1f-%.1f%%)\n",
            n_conc_abc, n_dir_abc, pct_abc,
            binom_abc$p.value,
            100 * binom_abc$conf.int[1], 100 * binom_abc$conf.int[2]))

# d(AxC) concordance
directional_unnorm <- result[unnorm_direction != "unchanged"]
directional_unnorm[, concordant := (loop_direction == "up_in_mutant" & unnorm_direction == "gained") |
                                   (loop_direction == "down_in_mutant" & unnorm_direction == "lost")]
n_conc_un <- sum(directional_unnorm$concordant)
n_dir_un  <- nrow(directional_unnorm)
pct_un    <- 100 * n_conc_un / n_dir_un
binom_un  <- binom.test(n_conc_un, n_dir_un, p = 0.5)

cat(sprintf("  d(AxC) concordance: %d / %d = %.1f%% (p = %.2e, 95%% CI: %.1f-%.1f%%)\n",
            n_conc_un, n_dir_un, pct_un,
            binom_un$p.value,
            100 * binom_un$conf.int[1], 100 * binom_un$conf.int[2]))

cat("  Reference: Step 9 independent-overlap concordance was 51.4%\n")

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

cat("\n  Concordance by loop type (dABC):\n")
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

# --- 6d: Effect size comparison ---
cat("\n=== EFFECT SIZE COMPARISON ===\n")

matched_delta <- abs(result$delta_ABC)
all_delta     <- abs(abc$delta_ABC)
wt <- wilcox.test(matched_delta, all_delta, alternative = "greater")
cat(sprintf("  Median |dABC| matched: %.4f vs all: %.4f\n",
            median(matched_delta), median(all_delta)))
cat(sprintf("  Wilcoxon p-value (matched > all): %.2e\n", wt$p.value))

matched_unnorm <- abs(result$delta_unnorm)
all_unnorm     <- abs(abc$delta_unnorm)
wt2 <- wilcox.test(matched_unnorm, all_unnorm, alternative = "greater")
cat(sprintf("  Median |d(AxC)| matched: %.6f vs all: %.6f\n",
            median(matched_unnorm), median(all_unnorm)))
cat(sprintf("  Wilcoxon p-value (matched > all): %.2e\n", wt2$p.value))

# --- 6e: FDR stratification ---
cat("\n=== FDR STRATIFICATION ===\n")

result[, fdr_stratum := fifelse(loop_FDR < 0.05, "significant",
                          fifelse(loop_FDR < 0.15, "exploratory", "other"))]

strat_abc <- result[abc_direction != "unchanged"]
strat_abc[, concordant := (loop_direction == "up_in_mutant" & abc_direction == "gained") |
                          (loop_direction == "down_in_mutant" & abc_direction == "lost")]

fdr_strata <- strat_abc[fdr_stratum %in% c("significant", "exploratory"),
  .(n_conc = sum(concordant), n_total = .N,
    pct = 100 * sum(concordant) / .N),
  by = fdr_stratum]

for (i in seq_len(nrow(fdr_strata))) {
  bt <- binom.test(fdr_strata$n_conc[i], fdr_strata$n_total[i], p = 0.5)
  cat(sprintf("  %s (FDR %s): %d/%d = %.1f%% (p = %.2e, 95%% CI: %.1f-%.1f%%)\n",
              fdr_strata$fdr_stratum[i],
              ifelse(fdr_strata$fdr_stratum[i] == "significant", "< 0.05", "0.05-0.15"),
              fdr_strata$n_conc[i], fdr_strata$n_total[i], fdr_strata$pct[i],
              bt$p.value, 100 * bt$conf.int[1], 100 * bt$conf.int[2]))
}

# --- 6f: RNA-seq 3-way concordance ---
cat("\n=== RNA-SEQ 3-WAY CONCORDANCE ===\n")

# Merge RNA-seq DE results onto diff_result
rnaseq_cols <- rnaseq[, .(TargetGene, de_log2FC = log2FC, de_padj = padj)]
result_de <- merge(result, rnaseq_cols, by.x = "target_gene", by.y = "TargetGene", all.x = TRUE)

n_de_matched <- sum(!is.na(result_de$de_log2FC))
cat(sprintf("  Triplets with RNA-seq data: %d / %d (%.1f%%)\n",
            n_de_matched, nrow(result_de), 100 * n_de_matched / nrow(result_de)))

# Filter to triplets where ABC is directional AND gene has significant DE
three_way <- result_de[abc_direction != "unchanged" & !is.na(de_padj) & de_padj < 0.05]
three_way[, de_direction := fifelse(de_log2FC > 0, "up", "down")]

cat(sprintf("  Triplets with directional ABC + significant DE: %d\n", nrow(three_way)))

if (nrow(three_way) > 0) {
  # 2-way: loop-ABC concordance (among this subset)
  three_way[, conc_loop_abc := (loop_direction == "up_in_mutant" & abc_direction == "gained") |
                               (loop_direction == "down_in_mutant" & abc_direction == "lost")]
  # 3-way: loop-ABC-DE all agree
  three_way[, conc_all_three := (loop_direction == "up_in_mutant" & abc_direction == "gained" & de_direction == "up") |
                                (loop_direction == "down_in_mutant" & abc_direction == "lost" & de_direction == "down")]

  n_3way <- nrow(three_way)
  n_conc_2 <- sum(three_way$conc_loop_abc)
  n_conc_3 <- sum(three_way$conc_all_three)
  pct_2way <- 100 * n_conc_2 / n_3way
  pct_3way <- 100 * n_conc_3 / n_3way

  binom_3way <- binom.test(n_conc_3, n_3way, p = 0.125)

  cat(sprintf("  Loop-ABC 2-way concordance (this subset): %d/%d = %.1f%%\n",
              n_conc_2, n_3way, pct_2way))
  cat(sprintf("  Loop-ABC-DE 3-way concordance: %d/%d = %.1f%% (p = %.2e vs 12.5%% chance)\n",
              n_conc_3, n_3way, pct_3way, binom_3way$p.value))

  # Gene-level concordance (deduplicate to strongest dABC per gene)
  gene_level <- three_way[order(-abs(delta_ABC))][!duplicated(target_gene)]
  n_gl <- nrow(gene_level)
  n_gl_conc3 <- sum(gene_level$conc_all_three)
  pct_gl <- 100 * n_gl_conc3 / n_gl
  binom_gl <- binom.test(n_gl_conc3, n_gl, p = 0.125)

  cat(sprintf("  Gene-level 3-way concordance: %d/%d = %.1f%% (p = %.2e vs 12.5%% chance)\n",
              n_gl_conc3, n_gl, pct_gl, binom_gl$p.value))

  # Direction breakdown table
  cat("\n  Direction breakdown (triplet-level):\n")
  combo_table <- three_way[, .N, by = .(loop_direction, abc_direction, de_direction)][order(-N)]
  for (i in seq_len(nrow(combo_table))) {
    cat(sprintf("    Loop=%s, ABC=%s, DE=%s: %d\n",
                combo_table$loop_direction[i], combo_table$abc_direction[i],
                combo_table$de_direction[i], combo_table$N[i]))
  }
} else {
  cat("  No triplets with both directional ABC and significant DE\n")
  three_way <- NULL
}

# --- 6g: K119ub at paired enhancers ---
cat("\n=== K119UB AT PAIRED ENHANCERS ===\n")

k119_cols <- k119ub[, .(chr, start, end, k119_ctrl_mean = ctrl_mean,
                         k119_mut_mean = mut_mean, k119_log2fc = log2fc,
                         k119_signal_class = signal_class)]
result_k119 <- merge(result, k119_cols,
                     by.x = c("enh_chr", "enh_start", "enh_end"),
                     by.y = c("chr", "start", "end"),
                     all.x = TRUE)

n_k119_matched <- sum(!is.na(result_k119$k119_log2fc))
n_k119_quant   <- sum(result_k119$k119_signal_class == "quantifiable", na.rm = TRUE)
cat(sprintf("  Triplets with K119ub data: %d / %d (%.1f%%)\n",
            n_k119_matched, nrow(result_k119), 100 * n_k119_matched / nrow(result_k119)))
cat(sprintf("  Triplets with quantifiable K119ub: %d\n", n_k119_quant))

# K119ub log2fc distribution by ABC direction
k119_quant <- result_k119[k119_signal_class == "quantifiable" & abc_direction != "unchanged"]
if (nrow(k119_quant) > 0) {
  k119_by_dir <- k119_quant[, .(
    median_k119 = median(k119_log2fc),
    mean_k119 = mean(k119_log2fc),
    n = .N
  ), by = abc_direction]

  cat("  K119ub log2fc by ABC direction (quantifiable enhancers):\n")
  for (i in seq_len(nrow(k119_by_dir))) {
    cat(sprintf("    %s: median=%.3f, mean=%.3f, n=%d\n",
                k119_by_dir$abc_direction[i],
                k119_by_dir$median_k119[i],
                k119_by_dir$mean_k119[i],
                k119_by_dir$n[i]))
  }

  # Spearman correlation: K119ub log2fc vs dABC
  cor_k119 <- cor.test(k119_quant$k119_log2fc, k119_quant$delta_ABC, method = "spearman")
  cat(sprintf("  Spearman K119ub log2fc vs dABC: rho = %.3f, p = %.2e (n = %d)\n",
              cor_k119$estimate, cor_k119$p.value, nrow(k119_quant)))

  # Wilcoxon: K119ub log2fc in gained vs lost pairs
  if (all(c("gained", "lost") %in% k119_quant$abc_direction)) {
    wt_k119 <- wilcox.test(k119_quant[abc_direction == "gained"]$k119_log2fc,
                           k119_quant[abc_direction == "lost"]$k119_log2fc)
    cat(sprintf("  Wilcoxon K119ub gained vs lost: p = %.2e\n", wt_k119$p.value))
  }
} else {
  cat("  No quantifiable K119ub data for directional ABC pairs\n")
  k119_quant <- NULL
}

# --- 6h: K119ub stratified by loop type ---
cat("\n=== K119UB BY LOOP TYPE ===\n")

# Use result_k119 which carries loop_type from the diff_lookup merge
k119_by_type <- result_k119[k119_signal_class == "quantifiable"]
k119_by_type_counts <- k119_by_type[, .N, by = loop_type][order(-N)]

# Keep only loop types with n >= 10
reliable_types <- k119_by_type_counts[N >= 10]$loop_type
k119_by_type_filt <- k119_by_type[loop_type %in% reliable_types]

cat(sprintf("  Loop types with n >= 10 quantifiable K119ub: %d / %d\n",
            length(reliable_types), nrow(k119_by_type_counts)))

if (length(reliable_types) > 0 && nrow(k119_by_type_filt) > 0) {
  # Per-type summary: median K119ub log2fc
  k119_type_summary <- k119_by_type_filt[, .(
    median_k119 = median(k119_log2fc),
    mean_k119 = mean(k119_log2fc),
    n = .N
  ), by = loop_type][order(median_k119)]

  cat("  K119ub log2fc by loop type (ordered by median):\n")
  for (i in seq_len(nrow(k119_type_summary))) {
    cat(sprintf("    %s: median=%.3f, mean=%.3f, n=%d\n",
                k119_type_summary$loop_type[i],
                k119_type_summary$median_k119[i],
                k119_type_summary$mean_k119[i],
                k119_type_summary$n[i]))
  }

  # Kruskal-Wallis test across loop types
  if (length(reliable_types) >= 2) {
    kw_test <- kruskal.test(k119_log2fc ~ loop_type, data = k119_by_type_filt)
    cat(sprintf("  Kruskal-Wallis across types: chi2 = %.2f, df = %d, p = %.2e\n",
                kw_test$statistic, kw_test$parameter, kw_test$p.value))
  }

  # Wilcoxon gained vs lost within each type
  k119_type_dir <- k119_by_type_filt[abc_direction %in% c("gained", "lost")]
  if (nrow(k119_type_dir) > 0) {
    cat("  Gained vs Lost K119ub by loop type:\n")
    for (lt in reliable_types) {
      gained_vals <- k119_type_dir[loop_type == lt & abc_direction == "gained"]$k119_log2fc
      lost_vals   <- k119_type_dir[loop_type == lt & abc_direction == "lost"]$k119_log2fc
      if (length(gained_vals) >= 3 && length(lost_vals) >= 3) {
        wt_lt <- wilcox.test(gained_vals, lost_vals)
        cat(sprintf("    %s: gained median=%.3f (n=%d), lost median=%.3f (n=%d), p=%.2e\n",
                    lt, median(gained_vals), length(gained_vals),
                    median(lost_vals), length(lost_vals), wt_lt$p.value))
      }
    }
  }
} else {
  cat("  No loop types with sufficient K119ub data (n >= 10)\n")
  k119_type_summary <- NULL
}

# --- 6i: Distance-dependent concordance ---
cat("\n=== DISTANCE-DEPENDENT CONCORDANCE ===\n")

result[, dist_bin := fcase(
  distance < 100000,  "<100kb",
  distance < 250000,  "100-250kb",
  distance < 500000,  "250-500kb",
  distance < 1000000, "500kb-1Mb",
  default = ">1Mb"
)]
result[, dist_bin := factor(dist_bin,
  levels = c("<100kb", "100-250kb", "250-500kb", "500kb-1Mb", ">1Mb"))]

dist_directional <- result[abc_direction != "unchanged"]
dist_directional[, concordant := (loop_direction == "up_in_mutant" & abc_direction == "gained") |
                                 (loop_direction == "down_in_mutant" & abc_direction == "lost")]

dist_conc <- dist_directional[, .(
  n_conc = sum(concordant),
  n_total = .N,
  pct = 100 * sum(concordant) / .N
), by = dist_bin][order(dist_bin)]

cat("  Concordance by distance bin:\n")
for (i in seq_len(nrow(dist_conc))) {
  bt_dist <- binom.test(dist_conc$n_conc[i], dist_conc$n_total[i], p = 0.5)
  cat(sprintf("    %s: %d/%d = %.1f%% (p = %.2e)\n",
              dist_conc$dist_bin[i], dist_conc$n_conc[i],
              dist_conc$n_total[i], dist_conc$pct[i], bt_dist$p.value))
}

# Spearman: distance vs |delta_ABC| to check magnitude effect
cor_dist <- cor.test(result$distance, abs(result$delta_ABC), method = "spearman")
cat(sprintf("  Spearman distance vs |dABC|: rho = %.3f, p = %.2e\n",
            cor_dist$estimate, cor_dist$p.value))

# --- 6j: Gene Ontology enrichment ---
cat("\n=== GENE ONTOLOGY ENRICHMENT ===\n")

# Extract unique gene symbols split by loop direction
up_genes   <- unique(result[loop_direction == "up_in_mutant"]$target_gene)
down_genes <- unique(result[loop_direction == "down_in_mutant"]$target_gene)
cat(sprintf("  Genes in up-loops: %d, down-loops: %d\n", length(up_genes), length(down_genes)))

# Convert symbols to Entrez IDs
convert_to_entrez <- function(symbols) {
  mapping <- AnnotationDbi::select(org.Mm.eg.db, keys = symbols,
                                   keytype = "SYMBOL", columns = "ENTREZID")
  mapping <- mapping[!is.na(mapping$ENTREZID), ]
  unique(mapping$ENTREZID)
}

up_entrez   <- convert_to_entrez(up_genes)
down_entrez <- convert_to_entrez(down_genes)
cat(sprintf("  Entrez IDs mapped: up=%d, down=%d\n", length(up_entrez), length(down_entrez)))

gene_clusters <- list(Up_in_Mutant = up_entrez, Down_in_Mutant = down_entrez)

# GO Biological Process enrichment
go_result <- NULL
if (length(up_entrez) >= 5 && length(down_entrez) >= 5) {
  go_result <- tryCatch(
    compareCluster(gene_clusters, fun = "enrichGO",
                   OrgDb = org.Mm.eg.db, ont = "BP",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                   readable = TRUE),
    error = function(e) { cat(sprintf("  GO enrichment error: %s\n", e$message)); NULL }
  )

  if (!is.null(go_result) && nrow(go_result@compareClusterResult) > 0) {
    cat(sprintf("  GO BP significant terms: %d\n", nrow(go_result@compareClusterResult)))
    fwrite(as.data.table(go_result@compareClusterResult),
           file.path(OUT_TSV_DIR, "5A_abc_go_enrichment.tsv"), sep = "\t") # Original: file.path(OUT_DIR, "paired_anchor_go_enrichment.tsv")
    cat("  Wrote paired_anchor_go_enrichment.tsv\n")

    # Top 5 per cluster
    go_top <- as.data.table(go_result@compareClusterResult)
    for (cl in unique(go_top$Cluster)) {
      cat(sprintf("  Top GO BP terms (%s):\n", cl))
      top5 <- head(go_top[Cluster == cl][order(p.adjust)], 5)
      for (j in seq_len(nrow(top5))) {
        cat(sprintf("    %s (q=%.2e, %s)\n",
                    top5$Description[j], top5$p.adjust[j], top5$GeneRatio[j]))
      }
    }
  } else {
    cat("  No significant GO BP terms found\n")
    go_result <- NULL
  }
} else {
  cat("  Too few genes for GO enrichment (need >= 5 per direction)\n")
}

# KEGG pathway enrichment
kegg_result <- NULL
if (length(up_entrez) >= 5 && length(down_entrez) >= 5) {
  kegg_result <- tryCatch(
    compareCluster(gene_clusters, fun = "enrichKEGG",
                   organism = "mmu",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05),
    error = function(e) { cat(sprintf("  KEGG enrichment error: %s\n", e$message)); NULL }
  )

  if (!is.null(kegg_result) && nrow(kegg_result@compareClusterResult) > 0) {
    cat(sprintf("  KEGG significant pathways: %d\n", nrow(kegg_result@compareClusterResult)))
    fwrite(as.data.table(kegg_result@compareClusterResult),
           file.path(OUT_TSV_DIR, "5A_abc_kegg_enrichment.tsv"), sep = "\t") # Original: file.path(OUT_DIR, "paired_anchor_kegg_enrichment.tsv")
    cat("  Wrote paired_anchor_kegg_enrichment.tsv\n")

    kegg_top <- as.data.table(kegg_result@compareClusterResult)
    for (cl in unique(kegg_top$Cluster)) {
      cat(sprintf("  Top KEGG pathways (%s):\n", cl))
      top5 <- head(kegg_top[Cluster == cl][order(p.adjust)], 5)
      for (j in seq_len(nrow(top5))) {
        cat(sprintf("    %s (q=%.2e, %s)\n",
                    top5$Description[j], top5$p.adjust[j], top5$GeneRatio[j]))
      }
    }
  } else {
    cat("  No significant KEGG pathways found\n")
    kegg_result <- NULL
  }
} else {
  cat("  Too few genes for KEGG enrichment\n")
}

# === STEP 7: SAVE RESULTS ===
cat("\nSaving results...\n")

# All matches (full background)
fwrite(all_result, file.path(OUT_TSV_DIR, "5C_paired_anchor_all_matches.tsv"), sep = "\t") # Original: file.path(OUT_DIR, "paired_anchor_all_matches.tsv")
cat(sprintf("  Wrote paired_anchor_all_matches.tsv (%d rows)\n", nrow(all_result)))

# Differential matches only (backward-compatible output)
# Select columns matching original output format
diff_out_cols <- c("loop_id", "loop_chr1", "loop_start1", "loop_end1",
                   "loop_chr2", "loop_start2", "loop_end2",
                   "loop_logFC", "loop_FDR", "loop_direction", "loop_type",
                   "diff_significant",
                   "enh_chr", "enh_start", "enh_end",
                   "target_gene", "tss_chr", "tss_start", "tss_end",
                   "abc_score_wt", "abc_score_ko", "delta_ABC", "delta_unnorm",
                   "distance", "match_case")
diff_export <- result[, ..diff_out_cols]
setnames(diff_export, "diff_significant", "loop_significant")
fwrite(diff_export, file.path(OUT_TSV_DIR, "5C_paired_anchor_matches.tsv"), sep = "\t") # Original: file.path(OUT_DIR, "paired_anchor_matches.tsv")
cat(sprintf("  Wrote paired_anchor_matches.tsv (%d rows)\n", nrow(diff_export)))

# Summary text
summary_lines <- c(
  "=== Paired-Anchor Loop-ABC Analysis Summary (v2) ===",
  sprintf("Date: %s", Sys.time()),
  "",
  "--- Input Data ---",
  sprintf("All loops (visual background): %d", nrow(all_loops)),
  sprintf("Differential loops (statistical analysis): %d", nrow(diff_loops)),
  sprintf("ABC E-G pairs: %d", nrow(abc)),
  "",
  "--- Match Statistics ---",
  sprintf("All-loop matches: %d (from %d unique loops)",
          nrow(all_result), length(unique(all_result$loop_key))),
  sprintf("Differential matches: %d (from %d unique diff loops)",
          nrow(result), n_unique_loops),
  sprintf("Unique E-G pairs confirmed by diff loops: %d / %d (%.4f%%)",
          n_unique_pairs, nrow(abc), 100 * n_unique_pairs / nrow(abc)),
  sprintf("Unique genes: %d", n_unique_genes),
  sprintf("Unique enhancers: %d", n_unique_enh),
  "",
  "--- Directional Concordance (differential loops) ---",
  sprintf("dABC concordance: %d / %d = %.1f%% (p = %.2e)",
          n_conc_abc, n_dir_abc, pct_abc, binom_abc$p.value),
  sprintf("d(AxC) concordance: %d / %d = %.1f%% (p = %.2e)",
          n_conc_un, n_dir_un, pct_un, binom_un$p.value),
  "Step 9 independent-overlap concordance: 51.4%",
  ""
)

# FDR stratification
summary_lines <- c(summary_lines, "--- FDR Stratification ---")
for (i in seq_len(nrow(fdr_strata))) {
  summary_lines <- c(summary_lines,
    sprintf("%s: %d/%d = %.1f%%",
            fdr_strata$fdr_stratum[i], fdr_strata$n_conc[i],
            fdr_strata$n_total[i], fdr_strata$pct[i]))
}
summary_lines <- c(summary_lines, "")

# RNA-seq concordance
if (!is.null(three_way) && nrow(three_way) > 0) {
  summary_lines <- c(summary_lines,
    "--- RNA-seq 3-Way Concordance ---",
    sprintf("Triplets with directional ABC + significant DE: %d", n_3way),
    sprintf("Loop-ABC 2-way concordance (this subset): %d/%d = %.1f%%",
            n_conc_2, n_3way, pct_2way),
    sprintf("Loop-ABC-DE 3-way concordance: %d/%d = %.1f%% (p = %.2e vs 12.5%% chance)",
            n_conc_3, n_3way, pct_3way, binom_3way$p.value),
    sprintf("Gene-level 3-way concordance: %d/%d = %.1f%%",
            n_gl_conc3, n_gl, pct_gl),
    ""
  )
}

# K119ub
if (!is.null(k119_quant) && nrow(k119_quant) > 0) {
  summary_lines <- c(summary_lines,
    "--- K119ub at Paired Enhancers ---",
    sprintf("Triplets with quantifiable K119ub: %d", n_k119_quant),
    sprintf("Spearman K119ub log2fc vs dABC: rho = %.3f, p = %.2e",
            cor_k119$estimate, cor_k119$p.value)
  )
  for (i in seq_len(nrow(k119_by_dir))) {
    summary_lines <- c(summary_lines,
      sprintf("  %s: median K119ub log2fc = %.3f (n=%d)",
              k119_by_dir$abc_direction[i],
              k119_by_dir$median_k119[i],
              k119_by_dir$n[i]))
  }
}

# K119ub by loop type
if (!is.null(k119_type_summary) && nrow(k119_type_summary) > 0) {
  summary_lines <- c(summary_lines, "",
    "--- K119ub by Loop Type ---")
  for (i in seq_len(nrow(k119_type_summary))) {
    summary_lines <- c(summary_lines,
      sprintf("  %s: median K119ub log2fc = %.3f (n=%d)",
              k119_type_summary$loop_type[i],
              k119_type_summary$median_k119[i],
              k119_type_summary$n[i]))
  }
  if (exists("kw_test")) {
    summary_lines <- c(summary_lines,
      sprintf("  Kruskal-Wallis p = %.2e", kw_test$p.value))
  }
}

# Distance-dependent concordance
summary_lines <- c(summary_lines, "",
  "--- Distance-Dependent Concordance ---")
for (i in seq_len(nrow(dist_conc))) {
  summary_lines <- c(summary_lines,
    sprintf("  %s: %d/%d = %.1f%%",
            dist_conc$dist_bin[i], dist_conc$n_conc[i],
            dist_conc$n_total[i], dist_conc$pct[i]))
}
summary_lines <- c(summary_lines,
  sprintf("  Spearman distance vs |dABC|: rho = %.3f, p = %.2e",
          cor_dist$estimate, cor_dist$p.value))

# GO enrichment top terms
if (!is.null(go_result) && nrow(go_result@compareClusterResult) > 0) {
  summary_lines <- c(summary_lines, "",
    "--- GO Biological Process Enrichment ---")
  go_top_summary <- as.data.table(go_result@compareClusterResult)
  for (cl in unique(go_top_summary$Cluster)) {
    summary_lines <- c(summary_lines, sprintf("  %s:", cl))
    top5 <- head(go_top_summary[Cluster == cl][order(p.adjust)], 5)
    for (j in seq_len(nrow(top5))) {
      summary_lines <- c(summary_lines,
        sprintf("    %s (q=%.2e)", top5$Description[j], top5$p.adjust[j]))
    }
  }
}

summary_lines <- c(summary_lines, "",
  "--- Effect Size ---",
  sprintf("Median |dABC| matched: %.4f vs all: %.4f (Wilcoxon p = %.2e)",
          median(matched_delta), median(all_delta), wt$p.value),
  sprintf("Median |d(AxC)| matched: %.6f vs all: %.6f (Wilcoxon p = %.2e)",
          median(matched_unnorm), median(all_unnorm), wt2$p.value)
)

writeLines(summary_lines, file.path(OUT_TSV_DIR, "5C_paired_anchor_summary.txt")) # Original: file.path(OUT_DIR, "paired_anchor_summary.txt")
cat(sprintf("  Wrote paired_anchor_summary.txt\n"))

# === STEP 8: PLOTS ===
cat("\nGenerating plots...\n")

# Shared theme
theme_paired <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# --- Plot 1: Concordance comparison bar chart ---
conc_df <- data.frame(
  method = c("Independent\noverlap (Step 9)", "Paired-anchor\n(dABC)", "Paired-anchor\n(d(AxC))"),
  concordance = c(51.4, pct_abc, pct_un),
  n = c(NA, n_dir_abc, n_dir_un)
)
conc_df$method <- factor(conc_df$method, levels = conc_df$method)

p_conc <- ggplot(conc_df, aes(x = method, y = concordance, fill = method)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("%.1f%%", concordance)), vjust = -0.5, size = 3.5) +
  geom_text(aes(label = ifelse(is.na(n), "", sprintf("n=%d", n))),
            vjust = 1.5, color = "white", size = 3) +
  scale_fill_manual(values = c("grey60", "#2166AC", "#B2182B")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Loop-ABC Concordance",
       subtitle = "Paired-anchor vs independent overlap",
       x = NULL, y = "Concordance (%)") +
  theme_paired +
  theme(panel.grid.major.x = element_blank())

save_multiformat_ggplot(p_conc, file.path(PLOT_DIR, "5A_paired_anchor_concordance"), width = 5, height = 5) # Original: file.path(PLOT_DIR, "paired_anchor_concordance")

# --- Plot 2: FDR-stratified concordance ---
fdr_label_map <- c("significant" = "FDR < 0.05", "exploratory" = "FDR 0.05-0.15")
fdr_color_map <- c("FDR < 0.05" = "#2166AC", "FDR 0.05-0.15" = "#92C5DE")
fdr_plot_df <- data.frame(
  stratum = fdr_label_map[fdr_strata$fdr_stratum],
  concordance = fdr_strata$pct,
  n = fdr_strata$n_total,
  stringsAsFactors = FALSE
)
fdr_plot_df$stratum <- factor(fdr_plot_df$stratum,
                              levels = c("FDR < 0.05", "FDR 0.05-0.15"))

p_fdr <- ggplot(fdr_plot_df, aes(x = stratum, y = concordance, fill = stratum)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("%.1f%%\nn=%d", concordance, n)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = fdr_color_map, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Concordance by FDR Stratum",
       subtitle = "dABC direction vs loop direction",
       x = NULL, y = "Concordance (%)") +
  theme_paired +
  theme(panel.grid.major.x = element_blank())

save_multiformat_ggplot(p_fdr, file.path(PLOT_DIR, "5A_fdr_stratified_concordance"), width = 4, height = 5) # Original: file.path(PLOT_DIR, "fdr_stratified_concordance")

# --- Plot 3: Loop logFC vs dABC scatter (with background) ---
# Subsample background for rendering if too large
bg_plot <- if (nrow(bg_result) > 20000) bg_result[sample(.N, 20000)] else bg_result

p_scatter1 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  # Background: all non-differential loops
  geom_point(data = bg_plot, aes(x = loop_logFC, y = delta_ABC),
             color = "grey80", alpha = 0.3, size = 0.8) +
  # Foreground: differential loops colored by direction
  geom_point(data = result[loop_direction == "up_in_mutant"],
             aes(x = loop_logFC, y = delta_ABC),
             color = "#D6604D", alpha = 0.6, size = 1.8) +
  geom_point(data = result[loop_direction == "down_in_mutant"],
             aes(x = loop_logFC, y = delta_ABC),
             color = "#4393C3", alpha = 0.6, size = 1.8) +
  geom_smooth(data = result, aes(x = loop_logFC, y = delta_ABC),
              method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  labs(title = "Loop logFC vs dABC",
       subtitle = sprintf("Grey: %d non-diff matches | Red/Blue: %d diff matches",
                          nrow(bg_result), nrow(result)),
       x = "Loop logFC (mutant / control)",
       y = "dABC (KO - WT)") +
  theme_paired

cor_abc <- cor.test(result$loop_logFC, result$delta_ABC, method = "spearman")
p_scatter1 <- p_scatter1 + annotate("text", x = Inf, y = Inf,
                   label = sprintf("rho = %.3f\np = %.2e",
                                   cor_abc$estimate, cor_abc$p.value),
                   hjust = 1.1, vjust = 1.5, size = 3.5)

save_multiformat_ggplot(p_scatter1, file.path(FIG4_PLOT_DIR, "4E_logFC_vs_deltaABC"), width = 6, height = 5) # Original: file.path(PLOT_DIR, "logFC_vs_deltaABC")

# --- Plot 4: Loop logFC vs d(AxC) scatter (with background) ---
p_scatter2 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(data = bg_plot, aes(x = loop_logFC, y = delta_unnorm),
             color = "grey80", alpha = 0.3, size = 0.8) +
  geom_point(data = result[loop_direction == "up_in_mutant"],
             aes(x = loop_logFC, y = delta_unnorm),
             color = "#D6604D", alpha = 0.6, size = 1.8) +
  geom_point(data = result[loop_direction == "down_in_mutant"],
             aes(x = loop_logFC, y = delta_unnorm),
             color = "#4393C3", alpha = 0.6, size = 1.8) +
  geom_smooth(data = result, aes(x = loop_logFC, y = delta_unnorm),
              method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  labs(title = "Loop logFC vs d(AxC)",
       subtitle = sprintf("Grey: %d non-diff matches | Red/Blue: %d diff matches",
                          nrow(bg_result), nrow(result)),
       x = "Loop logFC (mutant / control)",
       y = "d(Activity x Contact) (KO - WT)") +
  theme_paired

cor_unnorm <- cor.test(result$loop_logFC, result$delta_unnorm, method = "spearman")
p_scatter2 <- p_scatter2 + annotate("text", x = Inf, y = Inf,
                   label = sprintf("rho = %.3f\np = %.2e",
                                   cor_unnorm$estimate, cor_unnorm$p.value),
                   hjust = 1.1, vjust = 1.5, size = 3.5)

save_multiformat_ggplot(p_scatter2, file.path(PLOT_DIR, "5A_logFC_vs_delta_unnorm"), width = 6, height = 5) # Original: file.path(PLOT_DIR, "logFC_vs_delta_unnorm")

# --- Plot 5: RNA-seq 3-way concordance ---
if (!is.null(three_way) && nrow(three_way) > 0) {
  rnaseq_conc_df <- data.frame(
    comparison = c("Loop-ABC\n(2-way)", "Loop-ABC-DE\n(3-way)"),
    concordance = c(pct_2way, pct_3way),
    n = c(n_3way, n_3way),
    chance = c(50, 12.5),
    stringsAsFactors = FALSE
  )
  rnaseq_conc_df$comparison <- factor(rnaseq_conc_df$comparison, levels = rnaseq_conc_df$comparison)

  p_rnaseq <- ggplot(rnaseq_conc_df, aes(x = comparison, y = concordance, fill = comparison)) +
    geom_col(width = 0.5, show.legend = FALSE) +
    geom_hline(data = rnaseq_conc_df, aes(yintercept = chance),
               linetype = "dashed", color = "grey40") +
    geom_text(aes(label = sprintf("%.1f%%\nn=%d", concordance, n)), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("#2166AC", "#762A83")) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(title = "RNA-seq Integration",
         subtitle = "Among pairs with sig. DE (padj < 0.05)",
         x = NULL, y = "Concordance (%)") +
    theme_paired +
    theme(panel.grid.major.x = element_blank())

  save_multiformat_ggplot(p_rnaseq, file.path(PLOT_DIR, "5A_rnaseq_concordance"), width = 4, height = 5) # Original: file.path(PLOT_DIR, "rnaseq_concordance")
} else {
  p_rnaseq <- NULL
  cat("  Skipped rnaseq_concordance.pdf (no data)\n")
}

# --- Plot 6: K119ub at paired enhancers ---
if (!is.null(k119_quant) && nrow(k119_quant) > 0) {
  p_k119ub <- ggplot(k119_quant, aes(x = k119_log2fc, y = delta_ABC, color = abc_direction)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8,
                inherit.aes = FALSE, aes(x = k119_log2fc, y = delta_ABC)) +
    scale_color_manual(values = c("gained" = "#D6604D", "lost" = "#4393C3"),
                       name = "ABC direction") +
    labs(title = "K119ub Signal vs dABC",
         subtitle = sprintf("rho = %.3f, p = %.2e (n = %d quantifiable)",
                            cor_k119$estimate, cor_k119$p.value, nrow(k119_quant)),
         x = "H2AK119ub log2FC (mut/ctrl) at enhancer",
         y = "dABC (KO - WT)") +
    theme_paired +
    theme(legend.position = c(0.85, 0.15),
          legend.background = element_rect(fill = "white", color = "grey80"))

  save_multiformat_ggplot(p_k119ub, file.path(PLOT_DIR, "5A_k119ub_at_paired_enhancers"), width = 6, height = 5) # Original: file.path(PLOT_DIR, "k119ub_at_paired_enhancers")
} else {
  p_k119ub <- NULL
  cat("  Skipped k119ub_at_paired_enhancers.pdf (no data)\n")
}

# --- Plot 7: K119ub by loop type ---
if (!is.null(k119_type_summary) && nrow(k119_type_summary) > 0) {
  # Order loop types by median K119ub for readability
  type_order <- k119_type_summary[order(median_k119)]$loop_type
  k119_plot_data <- k119_by_type_filt[abc_direction %in% c("gained", "lost")]
  k119_plot_data[, loop_type := factor(loop_type, levels = type_order)]

  p_k119_type <- ggplot(k119_plot_data,
                        aes(x = loop_type, y = k119_log2fc, fill = abc_direction)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8, position = position_dodge(width = 0.75)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("gained" = "#D6604D", "lost" = "#4393C3"),
                      name = "ABC direction") +
    labs(title = "K119ub Signal Change by Loop Type",
         subtitle = sprintf("Quantifiable enhancers, types with n >= 10 (Kruskal-Wallis %s)",
                            if (exists("kw_test")) sprintf("p = %.2e", kw_test$p.value) else "N/A"),
         x = "Loop Type", y = "H2AK119ub log2FC (mut/ctrl)") +
    theme_paired +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  save_multiformat_ggplot(p_k119_type, file.path(PLOT_DIR, "5A_k119ub_by_loop_type"), width = 10, height = 6) # Original: file.path(PLOT_DIR, "k119ub_by_loop_type")
} else {
  cat("  Skipped k119ub_by_loop_type (insufficient data)\n")
}

# --- Plot 8: Distance-dependent concordance ---
dist_conc[, label := sprintf("%.1f%%\nn=%d", pct, n_total)]

p_dist <- ggplot(dist_conc, aes(x = dist_bin, y = pct, fill = dist_bin)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = label), vjust = -0.3, size = 3.2) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = "Concordance by Enhancer-Gene Distance",
       subtitle = "Loop direction vs dABC direction",
       x = "Enhancer-Gene Distance", y = "Concordance (%)") +
  theme_paired +
  theme(panel.grid.major.x = element_blank())

save_multiformat_ggplot(p_dist, file.path(PLOT_DIR, "5A_distance_concordance"), width = 6, height = 5) # Original: file.path(PLOT_DIR, "distance_concordance")

# --- Plot 9: GO BP dotplot ---
if (!is.null(go_result) && nrow(go_result@compareClusterResult) > 0) {
  p_go <- dotplot(go_result, showCategory = 15) +
    labs(title = "GO Biological Process Enrichment",
         subtitle = "Paired-anchor loop-connected genes by direction") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(p_go, file.path(PLOT_DIR, "5A_go_bp_dotplot"), width = 10, height = 9) # Original: file.path(PLOT_DIR, "go_bp_dotplot")
} else {
  cat("  Skipped go_bp_dotplot (no significant terms)\n")
}

# --- Plot 10: KEGG dotplot ---
if (!is.null(kegg_result) && nrow(kegg_result@compareClusterResult) > 0) {
  p_kegg <- dotplot(kegg_result, showCategory = 15) +
    labs(title = "KEGG Pathway Enrichment",
         subtitle = "Paired-anchor loop-connected genes by direction") +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  save_multiformat_ggplot(p_kegg, file.path(PLOT_DIR, "5A_kegg_dotplot"), width = 10, height = 9) # Original: file.path(PLOT_DIR, "kegg_dotplot")
} else {
  cat("  Skipped kegg_dotplot (no significant pathways)\n")
}

# --- Plot 11: Combined panel ---
# Build panel from available plots
panel_plots <- list(p_conc, p_fdr, p_scatter1, p_scatter2)
if (!is.null(p_rnaseq)) panel_plots <- c(panel_plots, list(p_rnaseq))
if (!is.null(p_k119ub)) panel_plots <- c(panel_plots, list(p_k119ub))

n_plots <- length(panel_plots)
if (n_plots == 6) {
  combined <- (panel_plots[[1]] | panel_plots[[2]] | panel_plots[[5]]) /
              (panel_plots[[3]] | panel_plots[[4]] | panel_plots[[6]]) +
    plot_annotation(title = "Paired-Anchor Loop-ABC Analysis (v2)",
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  panel_w <- 18
  panel_h <- 11
} else if (n_plots == 5) {
  combined <- (panel_plots[[1]] | panel_plots[[2]] | panel_plots[[5]]) /
              (panel_plots[[3]] | panel_plots[[4]] | plot_spacer()) +
    plot_annotation(title = "Paired-Anchor Loop-ABC Analysis (v2)",
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  panel_w <- 18
  panel_h <- 11
} else {
  combined <- (panel_plots[[1]] | panel_plots[[2]]) /
              (panel_plots[[3]] | panel_plots[[4]]) +
    plot_annotation(title = "Paired-Anchor Loop-ABC Analysis (v2)",
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
  panel_w <- 14
  panel_h <- 11
}

save_multiformat_ggplot(combined, file.path(SUPP_PLOT_DIR, "paired_anchor_panel"), width = panel_w, height = panel_h) # Original: file.path(PLOT_DIR, "paired_anchor_panel")

cat("\n=== Step 9b complete ===\n")
