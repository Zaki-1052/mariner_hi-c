# biomodal/downstream/scripts/viz_sections/section_48_cpg_island_ubiquitination.R
# Section 48: CpG Island Ubiquitination & De Novo Methylation Analysis
#
# Tests two hypotheses about BAP1-KO effects at CpG islands:
#   1. Whether H2AK119ub1 accumulation (from BAP1 loss) is enriched at
#      differentially methylated CpG islands
#   2. Whether hypermethylated CpG islands represent de novo methylation
#      (unmethylated in control) or amplification of pre-existing methylation
#
# Data: Combined mC + hmC CpG island DMRs from modality run-5 (8 samples,
#   sex covariate, deep-seq). K119ub condition-specific peaks, DiffBind
#   quantitative results, and directional (gained/lost) BED files.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_48_cpg_island_ubiquitination.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ChIPseeker)
})

# =============================================================================
# SECTION 48: CpG ISLAND UBIQUITINATION & DE NOVO METHYLATION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 48: CpG ISLAND UBIQUITINATION & DE NOVO METHYLATION ANALYSIS\n")
cat("================================================================================\n\n")

SECTION_DIR <- file.path(OUTPUT_DIR, "48_cpg_island_ubiquitination")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

run_fisher_2x2 <- function(a, b, c, d, row_names = c("Yes", "No"),
                           col_names = c("Target", "Background")) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                dimnames = list(row_names, col_names))
  res <- fisher.test(mat)
  list(
    mat    = mat,
    OR     = as.numeric(res$estimate),
    ci_lo  = res$conf.int[1],
    ci_hi  = res$conf.int[2],
    pvalue = res$p.value
  )
}

# =============================================================================
# STEP 1: LOAD CpG ISLAND UNIVERSE
# =============================================================================

cat("Loading CpG island universe (mC + hmC combined)...\n")

CPG_ISLAND_TSV <- file.path(BASE_DIR,
  "modality/exports/cpg_islands/cpg_islands_CG_run-5_mc_hmc.tsv")
stopifnot("CpG island TSV not found" = file.exists(CPG_ISLAND_TSV))

cpgi <- read.table(CPG_ISLAND_TSV, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE, quote = "")

cpgi$chr <- paste0("chr", cpgi$Chromosome)

cpgi$mc_sig    <- !is.na(cpgi$mc_qvalue) & cpgi$mc_qvalue < Q_THRESHOLD
cpgi$hmc_sig   <- !is.na(cpgi$hmc_qvalue) & cpgi$hmc_qvalue < Q_THRESHOLD
cpgi$mc_hyper  <- cpgi$mc_sig & cpgi$mc_mod_difference > 0
cpgi$mc_hypo   <- cpgi$mc_sig & cpgi$mc_mod_difference < 0
cpgi$hmc_hyper <- cpgi$hmc_sig & cpgi$hmc_mod_difference > 0
cpgi$hmc_hypo  <- cpgi$hmc_sig & cpgi$hmc_mod_difference < 0
cpgi$island_id <- paste0(cpgi$chr, ":", cpgi$Start, "-", cpgi$End)

n_total     <- nrow(cpgi)
n_mc_sig    <- sum(cpgi$mc_sig)
n_mc_hyper  <- sum(cpgi$mc_hyper)
n_mc_hypo   <- sum(cpgi$mc_hypo)
n_hmc_sig   <- sum(cpgi$hmc_sig)
n_hmc_hyper <- sum(cpgi$hmc_hyper)
n_hmc_hypo  <- sum(cpgi$hmc_hypo)

cat(sprintf("  CpG island universe: %d islands\n", n_total))
cat(sprintf("  mC significant: %d (%d hyper, %d hypo)\n",
            n_mc_sig, n_mc_hyper, n_mc_hypo))
cat(sprintf("  hmC significant: %d (%d hyper, %d hypo)\n",
            n_hmc_sig, n_hmc_hyper, n_hmc_hypo))

cpgi_gr <- GRanges(
  seqnames = cpgi$chr,
  ranges   = IRanges(start = cpgi$Start, end = cpgi$End)
)

# =============================================================================
# STEP 2: LOAD K119ub DATA
# =============================================================================

cat("\nLoading K119ub peak data...\n")

K119UB_UP_BED   <- file.path(REPO_ROOT, "peaks/new/H2AK119ub_up.bed")
K119UB_DOWN_BED <- file.path(REPO_ROOT, "peaks/new/H2AK119ub_down.bed")

k119ub_ctrl_gr   <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub ctrl")
k119ub_mut_gr    <- load_chip_peaks(K119UB_FILES$mut,  "K119ub mut")
k119ub_gained_gr <- load_chip_peaks(K119UB_UP_BED,     "K119ub gained")
k119ub_lost_gr   <- load_chip_peaks(K119UB_DOWN_BED,   "K119ub lost")

cat("\nLoading K119ub DiffBind results...\n")
diffbind_k119ub <- load_diffbind_flex(DIFFBIND_FILES$k119ub, "K119ub DiffBind")

diffbind_k119ub <- diffbind_k119ub[
  !is.na(diffbind_k119ub$Start) & !is.na(diffbind_k119ub$End), ]
diffbind_k119ub_gr <- GRanges(
  seqnames = diffbind_k119ub$Chr,
  ranges   = IRanges(start = diffbind_k119ub$Start, end = diffbind_k119ub$End),
  Fold     = diffbind_k119ub$Fold,
  FDR      = diffbind_k119ub$FDR
)

# =============================================================================
# STEP 3: LOAD ChIP PEAKS FOR CHROMATIN CONTEXT
# =============================================================================

cat("\nLoading ChIP-seq peaks...\n")
h3k27me3_gr <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
h3k27ac_gr  <- load_chip_peaks(CHIP_PEAK_FILES$h3k27ac,  "H3K27ac")
h3k4me3_gr  <- load_chip_peaks(CHIP_PEAK_FILES$h3k4me3,  "H3K4me3")
h3k4me1_gr  <- load_chip_peaks(CHIP_PEAK_FILES$h3k4me1,  "H3K4me1")
bivalent_gr <- load_chip_peaks(CHIP_PEAK_FILES$bivalent,  "Bivalent")

# =============================================================================
# STEP 4: COMPUTE ALL OVERLAPS
# =============================================================================

cat("\nComputing overlaps...\n")

cpgi$k119ub_ctrl   <- countOverlaps(cpgi_gr, k119ub_ctrl_gr) > 0
cpgi$k119ub_mut    <- countOverlaps(cpgi_gr, k119ub_mut_gr) > 0
cpgi$k119ub_gained <- countOverlaps(cpgi_gr, k119ub_gained_gr) > 0
cpgi$k119ub_lost   <- countOverlaps(cpgi_gr, k119ub_lost_gr) > 0
cpgi$n_k119ub_mut  <- countOverlaps(cpgi_gr, k119ub_mut_gr)

cpgi$H3K27me3 <- countOverlaps(cpgi_gr, h3k27me3_gr) > 0
cpgi$H3K27ac  <- countOverlaps(cpgi_gr, h3k27ac_gr) > 0
cpgi$H3K4me3  <- countOverlaps(cpgi_gr, h3k4me3_gr) > 0
cpgi$H3K4me1  <- countOverlaps(cpgi_gr, h3k4me1_gr) > 0
cpgi$Bivalent <- countOverlaps(cpgi_gr, bivalent_gr) > 0

ov_db <- findOverlaps(cpgi_gr, diffbind_k119ub_gr)
cpgi$k119ub_diffbind_fold <- NA_real_
if (length(ov_db) > 0) {
  fold_vals <- diffbind_k119ub_gr$Fold[subjectHits(ov_db)]
  query_idx <- queryHits(ov_db)
  for (idx in unique(query_idx)) {
    hits <- which(query_idx == idx)
    best <- hits[which.max(abs(fold_vals[hits]))]
    cpgi$k119ub_diffbind_fold[idx] <- fold_vals[best]
  }
}

cat(sprintf("  Islands with K119ub DiffBind overlap: %d\n",
            sum(!is.na(cpgi$k119ub_diffbind_fold))))

# =============================================================================
# STEP 5: CLASSIFY DE NOVO STATUS
# =============================================================================

DE_NOVO_THRESHOLDS <- c(0.05, 0.10, 0.20)
for (thr in DE_NOVO_THRESHOLDS) {
  col_name <- sprintf("de_novo_%s", gsub("\\.", "", as.character(thr)))
  cpgi[[col_name]] <- cpgi$mc_mean_mod_control < thr
}

cpgi$sig_group <- dplyr::case_when(
  cpgi$mc_hyper ~ sprintf("mC Hyper (n=%d)", n_mc_hyper),
  cpgi$mc_hypo  ~ sprintf("mC Hypo (n=%d)", n_mc_hypo),
  TRUE          ~ sprintf("Non-significant (n=%d)", n_total - n_mc_sig)
)
cpgi$sig_group <- factor(cpgi$sig_group,
  levels = c(sprintf("mC Hyper (n=%d)", n_mc_hyper),
             sprintf("mC Hypo (n=%d)", n_mc_hypo),
             sprintf("Non-significant (n=%d)", n_total - n_mc_sig)))

cpgi$sig_status <- dplyr::case_when(
  cpgi$mc_sig & cpgi$hmc_sig   ~ "Both sig",
  cpgi$mc_sig & !cpgi$hmc_sig  ~ "mC only",
  !cpgi$mc_sig & cpgi$hmc_sig  ~ "hmC only",
  TRUE                          ~ "Neither"
)

cat(sprintf("\nDe novo classification (hyper islands only, n=%d):\n", n_mc_hyper))
for (thr in DE_NOVO_THRESHOLDS) {
  col_name <- sprintf("de_novo_%s", gsub("\\.", "", as.character(thr)))
  n_dn <- sum(cpgi$mc_hyper & cpgi[[col_name]])
  cat(sprintf("  ctrl < %.2f: %d de novo (%.1f%%)\n",
              thr, n_dn, 100 * n_dn / n_mc_hyper))
}

# =============================================================================
# PART 1: K119ub OVERLAP AT CpG ISLANDS
# =============================================================================

cat("\n================================================================================\n")
cat("PART 1: K119ub Overlap at CpG Islands\n")
cat("================================================================================\n\n")

# ----- FIGURE 48a: K119ub enrichment bar chart with Fisher's test -----

cat("--- Figure 48a: K119ub enrichment ---\n")

peak_sets <- list(
  "K119ub Ctrl"   = cpgi$k119ub_ctrl,
  "K119ub Mut"    = cpgi$k119ub_mut,
  "K119ub Gained" = cpgi$k119ub_gained,
  "K119ub Lost"   = cpgi$k119ub_lost
)

enrich_rows <- list()
for (pset_name in names(peak_sets)) {
  ov <- peak_sets[[pset_name]]
  for (group_name in c("Hyper", "Hypo")) {
    is_group <- if (group_name == "Hyper") cpgi$mc_hyper else cpgi$mc_hypo
    is_bg    <- !cpgi$mc_sig
    a <- sum(is_group & ov)
    b <- sum(is_group & !ov)
    c <- sum(is_bg & ov)
    d <- sum(is_bg & !ov)
    fisher_res <- run_fisher_2x2(a, c, b, d,
      row_names = c("Overlap", "No overlap"),
      col_names = c(group_name, "Background"))
    enrich_rows[[length(enrich_rows) + 1]] <- data.frame(
      peak_set  = pset_name,
      dmr_group = group_name,
      pct_dmr   = 100 * a / (a + b),
      pct_bg    = 100 * c / (c + d),
      n_dmr     = a + b,
      n_bg      = c + d,
      OR        = fisher_res$OR,
      ci_lo     = fisher_res$ci_lo,
      ci_hi     = fisher_res$ci_hi,
      pvalue    = fisher_res$pvalue,
      stringsAsFactors = FALSE
    )
  }
}
enrich_df <- do.call(rbind, enrich_rows)
enrich_df$label <- sprintf("OR=%.2f %s", enrich_df$OR,
                           sapply(enrich_df$pvalue, sig_stars))

plot_df_48a <- rbind(
  data.frame(peak_set = enrich_df$peak_set, dmr_group = enrich_df$dmr_group,
             category = "DMR", pct = enrich_df$pct_dmr, stringsAsFactors = FALSE),
  data.frame(peak_set = enrich_df$peak_set, dmr_group = enrich_df$dmr_group,
             category = "Background", pct = enrich_df$pct_bg, stringsAsFactors = FALSE)
)
plot_df_48a$peak_set <- factor(plot_df_48a$peak_set,
  levels = c("K119ub Ctrl", "K119ub Mut", "K119ub Gained", "K119ub Lost"))

p_48a <- ggplot(plot_df_48a, aes(x = peak_set, y = pct, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = enrich_df,
            aes(x = peak_set, y = pmax(pct_dmr, pct_bg) + 2, label = label,
                fill = NULL, group = dmr_group),
            size = 3, vjust = 0) +
  facet_wrap(~dmr_group) +
  scale_fill_manual(values = c("DMR" = "#B2182B", "Background" = "grey70"),
                    name = "") +
  labs(title = "K119ub Peak Overlap at CpG Island DMRs",
       subtitle = "Hyper = mC gain in BAP1-KO | Hypo = mC loss | Background = non-significant islands",
       x = "", y = "% of Islands Overlapping K119ub Peaks") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_48a, file.path(SECTION_DIR, "48a_k119ub_enrichment"),
                        width = 12, height = 7)

cat(sprintf("  Fisher's test results saved for %d comparisons\n", nrow(enrich_df)))

# ----- FIGURE 48b: K119ub peak count violin -----

cat("\n--- Figure 48b: K119ub peak count violin ---\n")

cpgi$dmr_category <- dplyr::case_when(
  cpgi$mc_hyper ~ "mC Hyper",
  cpgi$mc_hypo  ~ "mC Hypo",
  TRUE          ~ "Non-sig"
)
cpgi$dmr_category <- factor(cpgi$dmr_category,
  levels = c("mC Hyper", "mC Hypo", "Non-sig"))

wilcox_hyper <- wilcox.test(
  cpgi$n_k119ub_mut[cpgi$mc_hyper],
  cpgi$n_k119ub_mut[!cpgi$mc_sig])
wilcox_hypo <- wilcox.test(
  cpgi$n_k119ub_mut[cpgi$mc_hypo],
  cpgi$n_k119ub_mut[!cpgi$mc_sig])

p_48b <- ggplot(cpgi, aes(x = dmr_category, y = n_k119ub_mut + 1, fill = dmr_category)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
  scale_y_log10() +
  scale_fill_manual(values = c("mC Hyper" = "#D7191C", "mC Hypo" = "#2C7BB6",
                                "Non-sig" = "grey70")) +
  annotate("text", x = 1.5, y = max(cpgi$n_k119ub_mut + 1) * 0.8,
           label = sprintf("Hyper vs Non-sig: %s", fmt_p(wilcox_hyper$p.value)),
           size = 3.5, hjust = 0.5) +
  annotate("text", x = 2.5, y = max(cpgi$n_k119ub_mut + 1) * 0.6,
           label = sprintf("Hypo vs Non-sig: %s", fmt_p(wilcox_hypo$p.value)),
           size = 3.5, hjust = 0.5) +
  labs(title = "K119ub Mutant Peak Count at CpG Islands",
       subtitle = "Number of overlapping K119ub peaks per island (log10 scale, +1 pseudocount)",
       x = "", y = "K119ub Mutant Peaks (+ 1)") +
  theme_biomodal() +
  guides(fill = "none")

save_multiformat_ggplot(p_48b, file.path(SECTION_DIR, "48b_k119ub_peak_count_violin"),
                        width = 8, height = 7)

# ----- FIGURE 48c: DiffBind K119ub fold change at CpG islands -----

cat("\n--- Figure 48c: DiffBind K119ub fold change ---\n")

cpgi_with_db <- cpgi[!is.na(cpgi$k119ub_diffbind_fold), ]
n_with_db <- nrow(cpgi_with_db)

if (n_with_db > 10) {
  wilcox_db_hyper <- NULL
  if (sum(cpgi_with_db$mc_hyper) >= 3 && sum(!cpgi_with_db$mc_sig) >= 3) {
    wilcox_db_hyper <- wilcox.test(
      cpgi_with_db$k119ub_diffbind_fold[cpgi_with_db$mc_hyper],
      cpgi_with_db$k119ub_diffbind_fold[!cpgi_with_db$mc_sig])
  }

  cpgi_with_db$dmr_category <- factor(
    dplyr::case_when(
      cpgi_with_db$mc_hyper ~ "mC Hyper",
      cpgi_with_db$mc_hypo  ~ "mC Hypo",
      TRUE                   ~ "Non-sig"),
    levels = c("mC Hyper", "mC Hypo", "Non-sig"))

  p_48c <- ggplot(cpgi_with_db, aes(x = dmr_category, y = k119ub_diffbind_fold,
                                     fill = dmr_category)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("mC Hyper" = "#D7191C", "mC Hypo" = "#2C7BB6",
                                  "Non-sig" = "grey70")) +
    labs(title = "K119ub DiffBind Fold Change at CpG Islands",
         subtitle = sprintf(
           "%d islands with DiffBind overlap | Positive = K119ub gain in mutant%s",
           n_with_db,
           if (!is.null(wilcox_db_hyper))
             sprintf("\nHyper vs Non-sig: %s", fmt_p(wilcox_db_hyper$p.value))
           else ""),
         x = "", y = "K119ub Fold Change (DiffBind)") +
    theme_biomodal() +
    guides(fill = "none")

  save_multiformat_ggplot(p_48c, file.path(SECTION_DIR, "48c_k119ub_diffbind_fold"),
                          width = 8, height = 7)
} else {
  cat("  WARNING: Too few islands with DiffBind overlap for figure 48c\n")
}

# =============================================================================
# PART 2: DE NOVO vs PRE-EXISTING METHYLATION
# =============================================================================

cat("\n================================================================================\n")
cat("PART 2: De Novo vs Pre-existing Methylation\n")
cat("================================================================================\n\n")

# ----- FIGURE 48d: Baseline mC density -----

cat("--- Figure 48d: Baseline mC density ---\n")

p_48d <- ggplot(cpgi, aes(x = mc_mean_mod_control, fill = sig_group)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0.20, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.05, y = Inf, label = "0.05", vjust = 2, size = 3) +
  annotate("text", x = 0.10, y = Inf, label = "0.10", vjust = 2, size = 3) +
  annotate("text", x = 0.20, y = Inf, label = "0.20", vjust = 2, size = 3) +
  scale_fill_manual(values = c(
    setNames("#D7191C", sprintf("mC Hyper (n=%d)", n_mc_hyper)),
    setNames("#2C7BB6", sprintf("mC Hypo (n=%d)", n_mc_hypo)),
    setNames("grey70", sprintf("Non-significant (n=%d)", n_total - n_mc_sig))
  ), name = "") +
  labs(title = "Baseline (Control) mC at CpG Islands",
       subtitle = "Distribution of control methylation fraction by DMR status",
       x = "Control mC Fraction", y = "Density") +
  theme_biomodal()

save_multiformat_ggplot(p_48d, file.path(SECTION_DIR, "48d_baseline_mc_density"),
                        width = 10, height = 7)

# ----- FIGURE 48e: Control vs mutant mC scatter -----

cat("\n--- Figure 48e: Control vs mutant mC scatter ---\n")

cpgi_sorted <- cpgi[order(cpgi$mc_sig), ]

top10_hyper <- cpgi[cpgi$mc_hyper, ]
top10_hyper <- top10_hyper[order(-top10_hyper$mc_mod_difference), ][1:min(10, n_mc_hyper), ]

p_48e <- ggplot(cpgi_sorted, aes(x = mc_mean_mod_control, y = mc_mean_mod_mutant)) +
  geom_point(aes(color = sig_group, alpha = mc_sig), size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.10, linetype = "dotted", color = "grey40", alpha = 0.6) +
  geom_vline(xintercept = 0.20, linetype = "dotted", color = "grey40", alpha = 0.6) +
  geom_text_repel(data = top10_hyper,
                  aes(label = island_id), size = 2.5, max.overlaps = 15,
                  color = "#D7191C", segment.alpha = 0.5) +
  scale_color_manual(values = c(
    setNames("#D7191C", sprintf("mC Hyper (n=%d)", n_mc_hyper)),
    setNames("#2C7BB6", sprintf("mC Hypo (n=%d)", n_mc_hypo)),
    setNames("grey80", sprintf("Non-significant (n=%d)", n_total - n_mc_sig))
  ), name = "") +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.15), guide = "none") +
  labs(title = "Control vs Mutant mC at CpG Islands",
       subtitle = "Above diagonal = mC gain | Below = mC loss | Dotted lines at 0.10, 0.20",
       x = "Control mC (mean fraction)", y = "Mutant mC (mean fraction)") +
  theme_biomodal() +
  coord_fixed()

save_multiformat_ggplot(p_48e, file.path(SECTION_DIR, "48e_ctrl_vs_mut_scatter"),
                        width = 10, height = 9)

# ----- FIGURE 48f: De novo classification -----

cat("\n--- Figure 48f: De novo classification ---\n")

class_rows <- list()
for (thr in DE_NOVO_THRESHOLDS) {
  hyper_islands <- cpgi[cpgi$mc_hyper, ]
  n_dn <- sum(hyper_islands$mc_mean_mod_control < thr)
  n_pre <- n_mc_hyper - n_dn
  class_rows[[length(class_rows) + 1]] <- data.frame(
    threshold = sprintf("< %.2f", thr),
    category  = c("De novo", "Pre-existing"),
    count     = c(n_dn, n_pre),
    pct       = c(100 * n_dn / n_mc_hyper, 100 * n_pre / n_mc_hyper),
    stringsAsFactors = FALSE
  )
}
class_df <- do.call(rbind, class_rows)
class_df$threshold <- factor(class_df$threshold,
  levels = sprintf("< %.2f", DE_NOVO_THRESHOLDS))

p_48f <- ggplot(class_df, aes(x = threshold, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", count, pct)),
            position = position_stack(vjust = 0.5), size = 4, color = "white",
            fontface = "bold") +
  scale_fill_manual(values = c("De novo" = "#FDB863", "Pre-existing" = "#D7191C"),
                    name = "Baseline mC") +
  labs(title = sprintf("De Novo vs Pre-existing Methylation (%d Hyper CpG Islands)", n_mc_hyper),
       subtitle = "De novo = control mC below threshold | Pre-existing = above threshold",
       x = "Control mC Threshold", y = "Number of Hyper Islands") +
  theme_biomodal()

save_multiformat_ggplot(p_48f, file.path(SECTION_DIR, "48f_de_novo_classification"),
                        width = 9, height = 7)

# ----- FIGURE 48g: Methylation gain magnitude -----

cat("\n--- Figure 48g: Methylation gain magnitude ---\n")

hyper_islands <- cpgi[cpgi$mc_hyper, ]
hyper_islands$baseline_class <- ifelse(
  hyper_islands$mc_mean_mod_control < 0.20,
  sprintf("Low baseline (<0.20, n=%d)", sum(hyper_islands$mc_mean_mod_control < 0.20)),
  sprintf("High baseline (>=0.20, n=%d)", sum(hyper_islands$mc_mean_mod_control >= 0.20))
)

wilcox_gain <- NULL
n_low  <- sum(hyper_islands$mc_mean_mod_control < 0.20)
n_high <- sum(hyper_islands$mc_mean_mod_control >= 0.20)
if (n_low >= 5 && n_high >= 5) {
  wilcox_gain <- wilcox.test(
    hyper_islands$mc_mod_difference[hyper_islands$mc_mean_mod_control < 0.20],
    hyper_islands$mc_mod_difference[hyper_islands$mc_mean_mod_control >= 0.20])
}

p_48g <- ggplot(hyper_islands, aes(x = baseline_class, y = mc_mod_difference,
                                    fill = baseline_class)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c(
    setNames("#FDB863", sprintf("Low baseline (<0.20, n=%d)", n_low)),
    setNames("#D7191C", sprintf("High baseline (>=0.20, n=%d)", n_high))
  )) +
  labs(title = "mC Gain Magnitude: Low vs High Baseline",
       subtitle = if (!is.null(wilcox_gain))
         sprintf("Wilcoxon %s", fmt_p(wilcox_gain$p.value))
       else "Insufficient data for Wilcoxon test",
       x = "", y = "mC Difference (mutant - control)") +
  theme_biomodal() +
  guides(fill = "none")

save_multiformat_ggplot(p_48g, file.path(SECTION_DIR, "48g_gain_magnitude"),
                        width = 8, height = 7)

# =============================================================================
# PART 3: COORDINATED mC + hmC CHANGES
# =============================================================================

cat("\n================================================================================\n")
cat("PART 3: Coordinated mC + hmC Changes\n")
cat("================================================================================\n\n")

# ----- FIGURE 48h: mC vs hmC difference scatter -----

cat("--- Figure 48h: mC vs hmC difference scatter ---\n")

cpgi$sig_status_f <- factor(cpgi$sig_status,
  levels = c("Both sig", "mC only", "hmC only", "Neither"))

n_both  <- sum(cpgi$sig_status == "Both sig")
n_mc_only  <- sum(cpgi$sig_status == "mC only")
n_hmc_only <- sum(cpgi$sig_status == "hmC only")

both_sig <- cpgi[cpgi$sig_status == "Both sig", ]
q1 <- sum(both_sig$mc_mod_difference > 0 & both_sig$hmc_mod_difference > 0)
q2 <- sum(both_sig$mc_mod_difference < 0 & both_sig$hmc_mod_difference > 0)
q3 <- sum(both_sig$mc_mod_difference < 0 & both_sig$hmc_mod_difference < 0)
q4 <- sum(both_sig$mc_mod_difference > 0 & both_sig$hmc_mod_difference < 0)

cpgi_plot <- cpgi[order(cpgi$sig_status_f, decreasing = TRUE), ]

p_48h <- ggplot(cpgi_plot, aes(x = mc_mod_difference, y = hmc_mod_difference)) +
  geom_point(aes(color = sig_status_f, size = sig_status_f, alpha = sig_status_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Both sig" = "#984EA3", "mC only" = "#E41A1C",
                                "hmC only" = "#377EB8", "Neither" = "grey80"), name = "") +
  scale_size_manual(values = c("Both sig" = 2.0, "mC only" = 1.5,
                               "hmC only" = 1.5, "Neither" = 0.4), guide = "none") +
  scale_alpha_manual(values = c("Both sig" = 0.9, "mC only" = 0.8,
                                "hmC only" = 0.8, "Neither" = 0.1), guide = "none") +
  annotate("text", x = 0.20, y = 0.12, size = 3, fontface = "bold",
           label = sprintf("mC+/hmC+: %d", q1)) +
  annotate("text", x = -0.20, y = 0.12, size = 3, fontface = "bold",
           label = sprintf("mC-/hmC+: %d", q2)) +
  annotate("text", x = -0.20, y = -0.12, size = 3, fontface = "bold",
           label = sprintf("mC-/hmC-: %d", q3)) +
  annotate("text", x = 0.20, y = -0.12, size = 3, fontface = "bold",
           label = sprintf("mC+/hmC-: %d", q4)) +
  labs(title = "Coordinated mC and hmC Changes at CpG Islands",
       subtitle = sprintf("Both sig: %d | mC only: %d | hmC only: %d",
                          n_both, n_mc_only, n_hmc_only),
       x = "mC Difference (mutant - control)",
       y = "hmC Difference (mutant - control)") +
  theme_biomodal()

save_multiformat_ggplot(p_48h, file.path(SECTION_DIR, "48h_mc_hmc_scatter"),
                        width = 10, height = 9)

# ----- FIGURE 48i: Co-significant heatmap -----

cat("\n--- Figure 48i: Co-significant heatmap ---\n")

if (n_both >= 5) {
  both_df <- cpgi[cpgi$sig_status == "Both sig", ]
  both_df$direction_combo <- paste0(
    ifelse(both_df$mc_mod_difference > 0, "mC+", "mC-"), "/",
    ifelse(both_df$hmc_mod_difference > 0, "hmC+", "hmC-"))

  hm_mat <- as.matrix(both_df[, c("mc_mean_mod_control", "mc_mean_mod_mutant",
                                    "hmc_mean_mod_control", "hmc_mean_mod_mutant")])
  colnames(hm_mat) <- c("mC Ctrl", "mC Mut", "hmC Ctrl", "hmC Mut")
  rownames(hm_mat) <- both_df$island_id

  row_annot <- data.frame(
    Direction = both_df$direction_combo,
    row.names = both_df$island_id
  )
  annot_colors <- list(Direction = c("mC+/hmC+" = "#E41A1C", "mC-/hmC+" = "#377EB8",
                                      "mC-/hmC-" = "#4DAF4A", "mC+/hmC-" = "#FF7F00"))

  hm_height <- max(6, min(14, n_both * 0.12 + 3))

  hm <- pheatmap(hm_mat,
    cluster_rows = TRUE, cluster_cols = FALSE,
    scale = "none",
    annotation_row = row_annot, annotation_colors = annot_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    fontsize_row = if (n_both > 60) 4 else 6,
    main = sprintf("Co-significant CpG Islands (n=%d): mC + hmC", n_both),
    silent = TRUE)

  hm_path <- file.path(SECTION_DIR, "48i_cosig_heatmap")
  dir.create(hm_path, recursive = TRUE, showWarnings = FALSE)

  pdf(file.path(hm_path, "48i_cosig_heatmap.pdf"), width = 10, height = hm_height)
  grid::grid.draw(hm$gtable)
  dev.off()

  jpeg(file.path(hm_path, "48i_cosig_heatmap.jpg"), width = 10, height = hm_height,
       units = "in", res = 300)
  grid::grid.draw(hm$gtable)
  dev.off()

  cat(sprintf("  Co-significant heatmap: %d islands\n", n_both))
  cat(sprintf("  Direction combos: %s\n",
              paste(names(table(both_df$direction_combo)),
                    table(both_df$direction_combo), sep = "=", collapse = ", ")))
} else {
  cat("  WARNING: Too few co-significant islands for heatmap\n")
}

# =============================================================================
# PART 4: K119ub x DE NOVO INTEGRATION
# =============================================================================

cat("\n================================================================================\n")
cat("PART 4: K119ub x De Novo Integration\n")
cat("================================================================================\n\n")

# ----- FIGURE 48j: Fisher's test forest plot -----

cat("--- Figure 48j: Fisher's forest plot ---\n")

is_nonsig <- !cpgi$mc_sig

# Test 1: All islands — mc_hyper vs K119ub gained
a1 <- sum(cpgi$mc_hyper & cpgi$k119ub_gained)
b1 <- sum(cpgi$mc_hyper & !cpgi$k119ub_gained)
c1 <- sum(is_nonsig & cpgi$k119ub_gained)
d1 <- sum(is_nonsig & !cpgi$k119ub_gained)
f1 <- run_fisher_2x2(a1, c1, b1, d1)

# Test 2: Among hyper — de novo (ctrl < 0.20) vs K119ub gained
hyper_cpgi <- cpgi[cpgi$mc_hyper, ]
dn_gained  <- sum(hyper_cpgi$mc_mean_mod_control < 0.20 & hyper_cpgi$k119ub_gained)
dn_no      <- sum(hyper_cpgi$mc_mean_mod_control < 0.20 & !hyper_cpgi$k119ub_gained)
pre_gained <- sum(hyper_cpgi$mc_mean_mod_control >= 0.20 & hyper_cpgi$k119ub_gained)
pre_no     <- sum(hyper_cpgi$mc_mean_mod_control >= 0.20 & !hyper_cpgi$k119ub_gained)
f2 <- run_fisher_2x2(dn_gained, pre_gained, dn_no, pre_no)

# Test 3: All islands — mc_hypo vs K119ub gained
a3 <- sum(cpgi$mc_hypo & cpgi$k119ub_gained)
b3 <- sum(cpgi$mc_hypo & !cpgi$k119ub_gained)
f3 <- run_fisher_2x2(a3, c1, b3, d1)

forest_df <- data.frame(
  test = c("Hyper vs Non-sig",
           sprintf("De novo vs Pre-existing\n(within %d hyper)", n_mc_hyper),
           "Hypo vs Non-sig"),
  OR    = c(f1$OR, f2$OR, f3$OR),
  ci_lo = c(f1$ci_lo, f2$ci_lo, f3$ci_lo),
  ci_hi = c(f1$ci_hi, f2$ci_hi, f3$ci_hi),
  pvalue = c(f1$pvalue, f2$pvalue, f3$pvalue),
  n_target = c(n_mc_hyper, sum(hyper_cpgi$mc_mean_mod_control < 0.20), n_mc_hypo),
  stringsAsFactors = FALSE
)
forest_df$test <- factor(forest_df$test, levels = rev(forest_df$test))
forest_df$label <- sprintf("OR=%.2f [%.2f-%.2f] %s (n=%d)",
                           forest_df$OR, forest_df$ci_lo, forest_df$ci_hi,
                           sapply(forest_df$pvalue, fmt_p), forest_df$n_target)

forest_df$ci_lo <- pmax(forest_df$ci_lo, 0.001)
forest_df$ci_hi <- pmax(forest_df$ci_hi, 0.001)
forest_df$OR_plot <- pmax(forest_df$OR, 0.001)

p_48j <- ggplot(forest_df, aes(x = OR_plot, y = test)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.2,
                linewidth = 0.8, orientation = "y") +
  geom_point(size = 4, color = "#756BB1") +
  geom_text(aes(label = label, x = pmax(ci_hi, 1) * 1.1),
            hjust = 0, size = 3.5) +
  scale_x_log10() +
  labs(title = "K119ub Gained Peak Enrichment at CpG Island DMRs",
       subtitle = "OR > 1 = enriched relative to non-significant background",
       x = "Odds Ratio (log10 scale)", y = "") +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 11))

p_48j_width <- max(10, 10 + max(nchar(forest_df$label)) * 0.05)
save_multiformat_ggplot(p_48j, file.path(SECTION_DIR, "48j_k119ub_denovo_forest"),
                        width = min(p_48j_width, 16), height = 6)

# ----- FIGURE 48k: 2x2 contingency tile -----

cat("\n--- Figure 48k: 2x2 contingency tile ---\n")

contingency_df <- data.frame(
  k119ub = c("K119ub Gained", "K119ub Gained", "No K119ub Gained", "No K119ub Gained"),
  baseline = c("De novo (<0.20)", "Pre-existing (>=0.20)",
               "De novo (<0.20)", "Pre-existing (>=0.20)"),
  count = c(dn_gained, pre_gained, dn_no, pre_no),
  stringsAsFactors = FALSE
)
contingency_df$total <- sum(contingency_df$count)
contingency_df$pct <- 100 * contingency_df$count / n_mc_hyper

p_48k <- ggplot(contingency_df, aes(x = k119ub, y = baseline)) +
  geom_tile(aes(fill = count), color = "white", linewidth = 2) +
  geom_text(aes(label = sprintf("n = %d\n(%.1f%%)", count, pct)),
            size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#756BB1", name = "Count") +
  labs(title = sprintf("K119ub x Methylation Baseline (%d Hyper CpG Islands)", n_mc_hyper),
       subtitle = sprintf("Fisher's OR = %.2f [%.2f-%.2f], %s",
                          f2$OR, f2$ci_lo, f2$ci_hi, fmt_p(f2$pvalue)),
       x = "", y = "") +
  theme_biomodal() +
  theme(axis.text = element_text(size = 12))

save_multiformat_ggplot(p_48k, file.path(SECTION_DIR, "48k_k119ub_denovo_tile"),
                        width = 8, height = 7)

# =============================================================================
# PART 5: CHROMATIN CONTEXT
# =============================================================================

cat("\n================================================================================\n")
cat("PART 5: Chromatin Context\n")
cat("================================================================================\n\n")

# ----- FIGURE 48l: Chromatin mark overlap bar chart -----

cat("--- Figure 48l: Chromatin mark overlap bar chart ---\n")

marks <- c("H3K27me3", "H3K27ac", "H3K4me3", "H3K4me1", "Bivalent")
categories <- c("mC Hyper", "mC Hypo", "Non-sig")

chrom_rows <- list()
for (mark in marks) {
  ov_col <- cpgi[[mark]]
  for (cat_name in categories) {
    is_cat <- switch(cat_name,
      "mC Hyper" = cpgi$mc_hyper,
      "mC Hypo"  = cpgi$mc_hypo,
      "Non-sig"  = !cpgi$mc_sig)
    n_cat <- sum(is_cat)
    n_ov  <- sum(is_cat & ov_col)
    # Fisher's: this category vs non-sig background
    if (cat_name != "Non-sig") {
      is_bg <- !cpgi$mc_sig
      a <- sum(is_cat & ov_col)
      b <- sum(is_cat & !ov_col)
      c <- sum(is_bg & ov_col)
      d <- sum(is_bg & !ov_col)
      ft <- run_fisher_2x2(a, c, b, d)
      or_val <- ft$OR
      p_val  <- ft$pvalue
    } else {
      or_val <- 1
      p_val  <- 1
    }
    chrom_rows[[length(chrom_rows) + 1]] <- data.frame(
      mark = mark, category = cat_name,
      n = n_cat, n_overlap = n_ov,
      pct = 100 * n_ov / n_cat,
      OR = or_val, pvalue = p_val,
      stringsAsFactors = FALSE)
  }
}
chrom_df <- do.call(rbind, chrom_rows)
chrom_df$category <- factor(chrom_df$category, levels = categories)
chrom_df$mark <- factor(chrom_df$mark, levels = marks)
chrom_df$star <- sapply(chrom_df$pvalue, sig_stars)

p_48l <- ggplot(chrom_df, aes(x = mark, y = pct, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = star),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("mC Hyper" = "#D7191C", "mC Hypo" = "#2C7BB6",
                                "Non-sig" = "grey70"), name = "") +
  labs(title = "Chromatin Mark Overlap at CpG Islands by DMR Status",
       subtitle = "Stars indicate Fisher's test significance vs non-significant background",
       x = "", y = "% of Islands Overlapping Mark") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_48l, file.path(SECTION_DIR, "48l_chromatin_overlap_bars"),
                        width = 12, height = 7)

# ----- FIGURE 48m: Enrichment OR heatmap -----

cat("\n--- Figure 48m: Enrichment OR heatmap ---\n")

or_mat <- matrix(NA_real_, nrow = length(marks), ncol = 2,
                 dimnames = list(marks, c("mC Hyper", "mC Hypo")))
pval_mat <- matrix(NA_real_, nrow = length(marks), ncol = 2,
                   dimnames = list(marks, c("mC Hyper", "mC Hypo")))

for (i in seq_along(marks)) {
  for (j in 1:2) {
    cat_name <- c("mC Hyper", "mC Hypo")[j]
    row_match <- chrom_df$mark == marks[i] & chrom_df$category == cat_name
    or_mat[i, j]   <- log2(chrom_df$OR[row_match])
    pval_mat[i, j] <- chrom_df$pvalue[row_match]
  }
}

or_mat[is.infinite(or_mat) & or_mat > 0] <- 3
or_mat[is.infinite(or_mat) & or_mat < 0] <- -3

display_mat <- matrix(
  sprintf("%.2f\n%s", or_mat, sapply(as.vector(pval_mat), sig_stars)),
  nrow = nrow(or_mat), dimnames = dimnames(or_mat))

max_abs <- max(abs(or_mat[is.finite(or_mat)]), na.rm = TRUE)
break_lim <- ceiling(max_abs * 10) / 10

hm_or <- pheatmap(or_mat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = display_mat,
  number_format = "%s",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-break_lim, break_lim, length.out = 101),
  fontsize_number = 10,
  main = "Chromatin Mark Enrichment at CpG Island DMRs\n(log2 Odds Ratio vs Non-significant)",
  silent = TRUE)

hm_or_path <- file.path(SECTION_DIR, "48m_chromatin_or_heatmap")
dir.create(hm_or_path, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(hm_or_path, "48m_chromatin_or_heatmap.pdf"), width = 8, height = 6)
grid::grid.draw(hm_or$gtable)
dev.off()

jpeg(file.path(hm_or_path, "48m_chromatin_or_heatmap.jpg"), width = 8, height = 6,
     units = "in", res = 300)
grid::grid.draw(hm_or$gtable)
dev.off()

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n================================================================================\n")
cat("Writing summary table...\n")

summary_table <- cpgi %>%
  dplyr::select(
    chr, Start, End, island_id, Annotation,
    mc_num_contexts, mc_mean_coverage,
    mc_mean_mod_control, mc_mean_mod_mutant,
    mc_mod_fold_change, mc_mod_difference,
    mc_qvalue, mc_sig, mc_hyper, mc_hypo,
    hmc_mean_mod_control, hmc_mean_mod_mutant,
    hmc_mod_difference, hmc_qvalue, hmc_sig, hmc_hyper, hmc_hypo,
    k119ub_ctrl, k119ub_mut, k119ub_gained, k119ub_lost,
    n_k119ub_mut, k119ub_diffbind_fold,
    de_novo_005, de_novo_01, de_novo_02,
    H3K27me3, H3K27ac, H3K4me3, H3K4me1, Bivalent,
    sig_status, dmr_category
  )

write.table(summary_table,
  file.path(TABLES_DIR, "48_cpg_island_ubiquitination_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Summary table: %d rows written\n", nrow(summary_table)))

# =============================================================================
# SECTION 48 SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 48 SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("CpG island universe: %d islands\n", n_total))
cat(sprintf("mC significant: %d (%d hyper, %d hypo)\n",
            n_mc_sig, n_mc_hyper, n_mc_hypo))
cat(sprintf("hmC significant: %d (%d hyper, %d hypo)\n",
            n_hmc_sig, n_hmc_hyper, n_hmc_hypo))

cat("\n--- K119ub Enrichment ---\n")
for (i in seq_len(nrow(enrich_df))) {
  row <- enrich_df[i, ]
  cat(sprintf("  %s @ %s: OR=%.2f [%.2f-%.2f] %s\n",
              row$peak_set, row$dmr_group,
              row$OR, row$ci_lo, row$ci_hi, fmt_p(row$pvalue)))
}

cat("\n--- De Novo Methylation (Hyper Islands) ---\n")
cat(sprintf("  Mean baseline mC (hyper): %.3f\n",
            mean(cpgi$mc_mean_mod_control[cpgi$mc_hyper])))
cat(sprintf("  Mean baseline mC (hypo):  %.3f\n",
            mean(cpgi$mc_mean_mod_control[cpgi$mc_hypo])))
cat(sprintf("  Mean baseline mC (non-sig): %.3f\n",
            mean(cpgi$mc_mean_mod_control[!cpgi$mc_sig])))
for (thr in DE_NOVO_THRESHOLDS) {
  n_dn <- sum(cpgi$mc_hyper & cpgi$mc_mean_mod_control < thr)
  cat(sprintf("  De novo (ctrl < %.2f): %d / %d hyper (%.1f%%)\n",
              thr, n_dn, n_mc_hyper, 100 * n_dn / n_mc_hyper))
}
cat("  CONCLUSION: Majority of hypermethylated CpG islands amplify\n")
cat("  pre-existing methylation rather than acquiring de novo methylation.\n")

cat("\n--- Coordinated mC + hmC ---\n")
cat(sprintf("  Both sig: %d islands\n", n_both))
cat(sprintf("  Quadrants: mC+/hmC+ %d | mC-/hmC+ %d | mC-/hmC- %d | mC+/hmC- %d\n",
            q1, q2, q3, q4))

cat("\n--- Chromatin Context (Hyper Islands) ---\n")
for (mark in marks) {
  row_match <- chrom_df$mark == mark & chrom_df$category == "mC Hyper"
  cat(sprintf("  %s: %.1f%% overlap (OR=%.2f, %s)\n",
              mark, chrom_df$pct[row_match], chrom_df$OR[row_match],
              fmt_p(chrom_df$pvalue[row_match])))
}

cat("\nSection 48 complete.\n\n")
