# biomodal/downstream/scripts/viz_sections/section_47_ctcf_anchor_methylation_overlay.R
# Section 47: CTCF Anchor Methylation Overlay
# Overlays biomodal 5mC/5hmC differential methylation onto CTCF loop anchors
# to test whether lost CTCF loop anchors are hypermethylated — paralleling
# the Flavahan/Bernstein IDH glioma work (Nature 2015).
#
# Distinct from Section 27 (gene-level methylation at all anchors):
#   This section works at the genomic region coordinate level (CpG islands,
#   shores, shelves, promoters) restricted to CTCF-anchored loops.
#
# Sub-analyses:
#   47a: CpG island methylation at CTCF loop anchors (core test)
#   47b: Multi-region comparison (CpG islands, shores, shelves, promoters)
#   47c: Coordinated mC-up / hmC-down at CTCF anchor CpG islands
#   47d: Distance-stratified analysis (controlling loop length confound)
#   47e: Methylation effect size vs loop logFC correlation
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_47_ctcf_anchor_methylation_overlay.R

source("scripts/viz_sections/_shared_config.R")

cat("\n")
cat("================================================================================\n")
cat("SECTION 47: CTCF ANCHOR METHYLATION OVERLAY\n")
cat("================================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

SECTION_DIR <- file.path(OUTPUT_DIR, "47_ctcf_anchor_methylation_overlay")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

DISTANCE_BINS <- c(0, 200e3, 500e3, 1e6, Inf)
DISTANCE_LABELS <- c("<200kb", "200-500kb", "500kb-1Mb", ">1Mb")

ANCHOR_COLORS <- c(
  "Lost CTCF anchor" = "#d73027",
  "Gained CTCF anchor" = "#4575b4",
  "Background" = "grey70"
)

STANDARD_CHROMS <- paste0("chr", c(1:19, "X"))

# =============================================================================
# LOCAL HELPERS
# =============================================================================

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  if (p < 0.001) return(sprintf("p = %.2e", p))
  if (p < 0.01) return(sprintf("p = %.3f", p))
  sprintf("p = %.2f", p)
}

extract_ctcf_anchors_local <- function(loops_df, ctcf_col_prefix = "CTCF_overlap") {
  a1_col <- paste0("anchor1_", ctcf_col_prefix)
  a2_col <- paste0("anchor2_", ctcf_col_prefix)

  a1 <- loops_df %>%
    filter(.data[[a1_col]] == TRUE) %>%
    transmute(chr = chr1, start = start1, end = end1,
              direction = direction, loop_distance = loop_distance,
              loop_id = loop_id)

  a2 <- loops_df %>%
    filter(.data[[a2_col]] == TRUE) %>%
    transmute(chr = chr2, start = start2, end = end2,
              direction = direction, loop_distance = loop_distance,
              loop_id = loop_id)

  anchors <- bind_rows(a1, a2)

  anchors <- anchors %>%
    group_by(chr, start, end, direction) %>%
    summarise(
      loop_distance = max(loop_distance),
      n_loops = n(),
      .groups = "drop"
    )

  anchors
}

region_dmr_to_granges <- function(dmr_df) {
  gr <- GRanges(
    seqnames = dmr_df$chr,
    ranges = IRanges(start = dmr_df$start, end = dmr_df$end)
  )
  mcols(gr)$mod_difference <- dmr_df$mod_difference
  mcols(gr)$dmr_qvalue <- dmr_df$dmr_qvalue
  mcols(gr)$significant <- dmr_df$significant
  mcols(gr)$direction <- dmr_df$direction
  mcols(gr)$dmr_idx <- seq_len(nrow(dmr_df))
  gr
}

run_fisher_2x2 <- function(n_test_pos, n_test_total, n_ref_pos, n_ref_total,
                            test_label, ref_label) {
  mat <- matrix(c(n_test_pos, n_ref_pos,
                  n_test_total - n_test_pos, n_ref_total - n_ref_pos),
                nrow = 2,
                dimnames = list(c("positive", "negative"),
                               c(test_label, ref_label)))
  ft <- fisher.test(mat)
  tibble(
    test_label = test_label,
    ref_label = ref_label,
    test_n = n_test_total,
    test_pos = n_test_pos,
    test_pct = 100 * n_test_pos / n_test_total,
    ref_n = n_ref_total,
    ref_pos = n_ref_pos,
    ref_pct = 100 * n_ref_pos / n_ref_total,
    odds_ratio = as.numeric(ft$estimate),
    ci_lower = ft$conf.int[1],
    ci_upper = ft$conf.int[2],
    p_value = ft$p.value
  )
}

# =============================================================================
# INPUT VALIDATION
# =============================================================================

cat("--- Validating inputs ---\n")
stopifnot("Loop file not found" = file.exists(LOOP_FILES$late))
stopifnot("CpG island mC DMRs not found" = file.exists(DATA_PATHS$cpg_islands_mc))
stopifnot("CpG island hmC DMRs not found" = file.exists(DATA_PATHS$cpg_islands_hmc))
stopifnot("CpG shores mC DMRs not found" = file.exists(DATA_PATHS$cpg_shores_mc))
stopifnot("CpG shores hmC DMRs not found" = file.exists(DATA_PATHS$cpg_shores_hmc))
stopifnot("CpG shelves mC DMRs not found" = file.exists(DATA_PATHS$cpg_shelves_mc))
stopifnot("CpG shelves hmC DMRs not found" = file.exists(DATA_PATHS$cpg_shelves_hmc))
stopifnot("Promoter mC DMRs not found" = file.exists(DATA_PATHS$promoters_mc))
stopifnot("Promoter hmC DMRs not found" = file.exists(DATA_PATHS$promoters_hmc))
cat("  All inputs validated.\n")

# =============================================================================
# DATA LOADING
# =============================================================================

cat("\n--- Loading data ---\n")

loops <- read.table(LOOP_FILES$late, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Loops: %d total (%d lost, %d gained)\n",
            nrow(loops),
            sum(loops$direction == "down_in_mutant"),
            sum(loops$direction == "up_in_mutant")))

cpgi_mc <- load_dmr_bed(DATA_PATHS$cpg_islands_mc)
cpgi_hmc <- load_dmr_bed(DATA_PATHS$cpg_islands_hmc)
shores_mc <- load_dmr_bed(DATA_PATHS$cpg_shores_mc)
shores_hmc <- load_dmr_bed(DATA_PATHS$cpg_shores_hmc)
shelves_mc <- load_dmr_bed(DATA_PATHS$cpg_shelves_mc)
shelves_hmc <- load_dmr_bed(DATA_PATHS$cpg_shelves_hmc)
promoters_mc <- load_dmr_bed(DATA_PATHS$promoters_mc)
promoters_hmc <- load_dmr_bed(DATA_PATHS$promoters_hmc)

cat(sprintf("  CpG islands: %d mC DMRs, %d hmC DMRs\n", nrow(cpgi_mc), nrow(cpgi_hmc)))
cat(sprintf("  CpG shores:  %d mC DMRs, %d hmC DMRs\n", nrow(shores_mc), nrow(shores_hmc)))
cat(sprintf("  CpG shelves: %d mC DMRs, %d hmC DMRs\n", nrow(shelves_mc), nrow(shelves_hmc)))
cat(sprintf("  Promoters:   %d mC DMRs, %d hmC DMRs\n", nrow(promoters_mc), nrow(promoters_hmc)))

# Extract CTCF ChIP+ anchors
cat("\n--- Extracting CTCF anchors ---\n")

anchors_lost <- extract_ctcf_anchors_local(
  loops %>% filter(direction == "down_in_mutant"), "CTCF_overlap")
anchors_gained <- extract_ctcf_anchors_local(
  loops %>% filter(direction == "up_in_mutant"), "CTCF_overlap")

cat(sprintf("  Unique CTCF anchors (ChIP): %d lost, %d gained\n",
            nrow(anchors_lost), nrow(anchors_gained)))

lost_gr <- GRanges(seqnames = anchors_lost$chr,
                   ranges = IRanges(start = anchors_lost$start, end = anchors_lost$end))
mcols(lost_gr)$loop_distance <- anchors_lost$loop_distance
mcols(lost_gr)$n_loops <- anchors_lost$n_loops

gained_gr <- GRanges(seqnames = anchors_gained$chr,
                     ranges = IRanges(start = anchors_gained$start, end = anchors_gained$end))
mcols(gained_gr)$loop_distance <- anchors_gained$loop_distance
mcols(gained_gr)$n_loops <- anchors_gained$n_loops

# Motif-based anchors (for sensitivity in 47a)
anchors_motif_lost <- extract_ctcf_anchors_local(
  loops %>% filter(direction == "down_in_mutant"), "CTCF_motif_overlap")
anchors_motif_gained <- extract_ctcf_anchors_local(
  loops %>% filter(direction == "up_in_mutant"), "CTCF_motif_overlap")
cat(sprintf("  Unique CTCF anchors (motif): %d lost, %d gained\n",
            nrow(anchors_motif_lost), nrow(anchors_motif_gained)))


# =============================================================================
# 47a: CpG ISLAND METHYLATION AT CTCF LOOP ANCHORS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("47a: CpG ISLAND METHYLATION AT CTCF LOOP ANCHORS\n")
cat("================================================================================\n\n")

cpgi_mc_gr <- region_dmr_to_granges(cpgi_mc)

hits_lost <- findOverlaps(cpgi_mc_gr, lost_gr, ignore.strand = TRUE)
hits_gained <- findOverlaps(cpgi_mc_gr, gained_gr, ignore.strand = TRUE)

cpgi_at_lost_idx <- unique(queryHits(hits_lost))
cpgi_at_gained_idx <- unique(queryHits(hits_gained))
cpgi_at_both_idx <- intersect(cpgi_at_lost_idx, cpgi_at_gained_idx)

cpgi_mc$anchor_group <- "Background"
cpgi_mc$anchor_group[cpgi_at_lost_idx] <- "Lost CTCF anchor"
cpgi_mc$anchor_group[cpgi_at_gained_idx] <- "Gained CTCF anchor"
cpgi_mc$anchor_group[cpgi_at_both_idx] <- "Both"

cat(sprintf("  CpG islands at lost CTCF anchors:   %d\n", length(cpgi_at_lost_idx)))
cat(sprintf("  CpG islands at gained CTCF anchors: %d\n", length(cpgi_at_gained_idx)))
cat(sprintf("  CpG islands at both:                %d\n", length(cpgi_at_both_idx)))
cat(sprintf("  CpG islands at neither (background): %d\n",
            sum(cpgi_mc$anchor_group == "Background")))

# Same for hmC
cpgi_hmc_gr <- region_dmr_to_granges(cpgi_hmc)
cpgi_hmc$anchor_group <- "Background"
cpgi_hmc$anchor_group[unique(queryHits(findOverlaps(cpgi_hmc_gr, lost_gr)))] <- "Lost CTCF anchor"
cpgi_hmc$anchor_group[unique(queryHits(findOverlaps(cpgi_hmc_gr, gained_gr)))] <- "Gained CTCF anchor"
both_hmc <- intersect(
  unique(queryHits(findOverlaps(cpgi_hmc_gr, lost_gr))),
  unique(queryHits(findOverlaps(cpgi_hmc_gr, gained_gr)))
)
cpgi_hmc$anchor_group[both_hmc] <- "Both"

# --- Fisher's exact tests (excluding "Both" category) ---
cat("\n  Fisher's exact tests (mC hypermethylation):\n")

fisher_47a <- list()

# mC hyper: lost vs gained
lost_mc <- cpgi_mc %>% filter(anchor_group == "Lost CTCF anchor")
gained_mc <- cpgi_mc %>% filter(anchor_group == "Gained CTCF anchor")
bg_mc <- cpgi_mc %>% filter(anchor_group == "Background")

fisher_47a[[1]] <- run_fisher_2x2(
  sum(lost_mc$significant & lost_mc$direction == "Hypermethylated"), nrow(lost_mc),
  sum(gained_mc$significant & gained_mc$direction == "Hypermethylated"), nrow(gained_mc),
  "Lost_vs_Gained", "mC_hyper")
fisher_47a[[1]]$modality <- "mC"
fisher_47a[[1]]$comparison <- "lost_vs_gained"

cat(sprintf("    Lost vs Gained (mC hyper): OR=%.2f [%.2f-%.2f], %s\n",
            fisher_47a[[1]]$odds_ratio, fisher_47a[[1]]$ci_lower,
            fisher_47a[[1]]$ci_upper, fmt_p(fisher_47a[[1]]$p_value)))

# mC hyper: lost vs background
fisher_47a[[2]] <- run_fisher_2x2(
  sum(lost_mc$significant & lost_mc$direction == "Hypermethylated"), nrow(lost_mc),
  sum(bg_mc$significant & bg_mc$direction == "Hypermethylated"), nrow(bg_mc),
  "Lost_vs_BG", "mC_hyper")
fisher_47a[[2]]$modality <- "mC"
fisher_47a[[2]]$comparison <- "lost_vs_background"

cat(sprintf("    Lost vs Background (mC hyper): OR=%.2f [%.2f-%.2f], %s\n",
            fisher_47a[[2]]$odds_ratio, fisher_47a[[2]]$ci_lower,
            fisher_47a[[2]]$ci_upper, fmt_p(fisher_47a[[2]]$p_value)))

# hmC hypo: lost vs gained
lost_hmc <- cpgi_hmc %>% filter(anchor_group == "Lost CTCF anchor")
gained_hmc <- cpgi_hmc %>% filter(anchor_group == "Gained CTCF anchor")

fisher_47a[[3]] <- run_fisher_2x2(
  sum(lost_hmc$significant & lost_hmc$direction == "Hypomethylated"), nrow(lost_hmc),
  sum(gained_hmc$significant & gained_hmc$direction == "Hypomethylated"), nrow(gained_hmc),
  "Lost_vs_Gained", "hmC_hypo")
fisher_47a[[3]]$modality <- "hmC"
fisher_47a[[3]]$comparison <- "lost_vs_gained"

cat(sprintf("    Lost vs Gained (hmC hypo): OR=%.2f [%.2f-%.2f], %s\n",
            fisher_47a[[3]]$odds_ratio, fisher_47a[[3]]$ci_lower,
            fisher_47a[[3]]$ci_upper, fmt_p(fisher_47a[[3]]$p_value)))

# Motif-based sensitivity
motif_lost_gr <- GRanges(seqnames = anchors_motif_lost$chr,
                         ranges = IRanges(start = anchors_motif_lost$start,
                                          end = anchors_motif_lost$end))
motif_gained_gr <- GRanges(seqnames = anchors_motif_gained$chr,
                           ranges = IRanges(start = anchors_motif_gained$start,
                                            end = anchors_motif_gained$end))

cpgi_motif_group <- rep("Background", nrow(cpgi_mc))
cpgi_motif_group[unique(queryHits(findOverlaps(cpgi_mc_gr, motif_lost_gr)))] <- "Lost"
cpgi_motif_group[unique(queryHits(findOverlaps(cpgi_mc_gr, motif_gained_gr)))] <- "Gained"

motif_lost_n <- sum(cpgi_motif_group == "Lost")
motif_gained_n <- sum(cpgi_motif_group == "Gained")
motif_lost_hyper <- sum(cpgi_mc$significant[cpgi_motif_group == "Lost"] &
                        cpgi_mc$direction[cpgi_motif_group == "Lost"] == "Hypermethylated")
motif_gained_hyper <- sum(cpgi_mc$significant[cpgi_motif_group == "Gained"] &
                          cpgi_mc$direction[cpgi_motif_group == "Gained"] == "Hypermethylated")

fisher_47a[[4]] <- run_fisher_2x2(
  motif_lost_hyper, motif_lost_n,
  motif_gained_hyper, motif_gained_n,
  "Motif_Lost_vs_Gained", "mC_hyper")
fisher_47a[[4]]$modality <- "mC_motif"
fisher_47a[[4]]$comparison <- "lost_vs_gained_motif"

cat(sprintf("    Motif sensitivity (mC hyper): OR=%.2f [%.2f-%.2f], %s\n",
            fisher_47a[[4]]$odds_ratio, fisher_47a[[4]]$ci_lower,
            fisher_47a[[4]]$ci_upper, fmt_p(fisher_47a[[4]]$p_value)))

fisher_47a_df <- bind_rows(fisher_47a)
fisher_47a_df$fdr <- p.adjust(fisher_47a_df$p_value, method = "BH")

write.table(fisher_47a_df, file.path(TABLES_DIR, "47a_fisher_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Wilcoxon tests ---
cat("\n  Wilcoxon rank-sum tests (mod_difference):\n")

wt_mc <- wilcox.test(lost_mc$mod_difference, gained_mc$mod_difference)
cat(sprintf("    mC mod_diff lost vs gained: %s\n", fmt_p(wt_mc$p.value)))
cat(sprintf("      Lost median: %.4f, Gained median: %.4f\n",
            median(lost_mc$mod_difference), median(gained_mc$mod_difference)))

wt_hmc <- wilcox.test(lost_hmc$mod_difference, gained_hmc$mod_difference)
cat(sprintf("    hmC mod_diff lost vs gained: %s\n", fmt_p(wt_hmc$p.value)))
cat(sprintf("      Lost median: %.4f, Gained median: %.4f\n",
            median(lost_hmc$mod_difference), median(gained_hmc$mod_difference)))

# --- Per-CpG-island overlay table ---
overlay_47a <- cpgi_mc %>%
  dplyr::select(chr, start, end, mod_difference, dmr_qvalue, significant, direction, anchor_group) %>%
  dplyr::rename(mc_mod_difference = mod_difference, mc_qvalue = dmr_qvalue,
                mc_significant = significant, mc_direction = direction) %>%
  left_join(
    cpgi_hmc %>%
      dplyr::select(chr, start, end, mod_difference, dmr_qvalue, significant, direction) %>%
      dplyr::rename(hmc_mod_difference = mod_difference, hmc_qvalue = dmr_qvalue,
                    hmc_significant = significant, hmc_direction = direction),
    by = c("chr", "start", "end")
  )

write.table(overlay_47a, file.path(TABLES_DIR, "47a_cpgi_ctcf_anchor_methylation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plots ---
cat("\n  Generating plots...\n")

# Plot 47a.1: Grouped bar chart of % hypermethylated
bar_data <- cpgi_mc %>%
  filter(anchor_group != "Both") %>%
  mutate(anchor_group = factor(anchor_group,
                               levels = c("Lost CTCF anchor", "Gained CTCF anchor", "Background"))) %>%
  group_by(anchor_group) %>%
  summarise(
    n = n(),
    n_hyper = sum(significant & direction == "Hypermethylated"),
    pct_hyper = 100 * n_hyper / n,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%d/%d", n_hyper, n))

p47a_bar <- ggplot(bar_data, aes(x = anchor_group, y = pct_hyper, fill = anchor_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = ANCHOR_COLORS, guide = "none") +
  labs(x = NULL,
       y = "% Significantly Hypermethylated (mC, q<0.05)",
       title = "CpG Island Hypermethylation at CTCF Loop Anchors",
       subtitle = sprintf("Fisher's OR (lost vs gained) = %.2f, %s",
                          fisher_47a[[1]]$odds_ratio,
                          fmt_p(fisher_47a[[1]]$p_value))) +
  theme_biomodal() +
  coord_cartesian(ylim = c(0, max(bar_data$pct_hyper) * 1.2))

save_multiformat_ggplot(p47a_bar,
                        file.path(SECTION_DIR, "47a_cpgi_mc_direction_barchart"),
                        width = 8, height = 6)

# Plot 47a.2: mC mod_difference violin
violin_mc_data <- cpgi_mc %>%
  filter(anchor_group != "Both") %>%
  mutate(anchor_group = factor(anchor_group,
                               levels = c("Lost CTCF anchor", "Gained CTCF anchor", "Background")))

p47a_mc <- ggplot(violin_mc_data, aes(x = anchor_group, y = mod_difference, fill = anchor_group)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = ANCHOR_COLORS, guide = "none") +
  labs(x = NULL,
       y = "mC Modification Difference (mutant - control)",
       title = "5mC Differential Methylation at CpG Islands",
       subtitle = sprintf("Wilcoxon lost vs gained: %s", fmt_p(wt_mc$p.value))) +
  theme_biomodal()

save_multiformat_ggplot(p47a_mc,
                        file.path(SECTION_DIR, "47a_cpgi_mc_violin"),
                        width = 7, height = 6)

# Plot 47a.3: hmC mod_difference violin
violin_hmc_data <- cpgi_hmc %>%
  filter(anchor_group != "Both") %>%
  mutate(anchor_group = factor(anchor_group,
                               levels = c("Lost CTCF anchor", "Gained CTCF anchor", "Background")))

p47a_hmc <- ggplot(violin_hmc_data, aes(x = anchor_group, y = mod_difference, fill = anchor_group)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = ANCHOR_COLORS, guide = "none") +
  labs(x = NULL,
       y = "hmC Modification Difference (mutant - control)",
       title = "5hmC Differential Methylation at CpG Islands",
       subtitle = sprintf("Wilcoxon lost vs gained: %s", fmt_p(wt_hmc$p.value))) +
  theme_biomodal()

save_multiformat_ggplot(p47a_hmc,
                        file.path(SECTION_DIR, "47a_cpgi_hmc_violin"),
                        width = 7, height = 6)


# =============================================================================
# 47b: MULTI-REGION COMPARISON
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("47b: MULTI-REGION COMPARISON\n")
cat("================================================================================\n\n")

region_list <- list(
  "CpG islands" = list(mc = cpgi_mc, hmc = cpgi_hmc),
  "CpG shores"  = list(mc = shores_mc, hmc = shores_hmc),
  "CpG shelves" = list(mc = shelves_mc, hmc = shelves_hmc),
  "Promoters"   = list(mc = promoters_mc, hmc = promoters_hmc)
)

fisher_47b <- list()
wilcox_47b <- list()
idx <- 1

for (region_name in names(region_list)) {
  for (modality in c("mc", "hmc")) {
    dmr_df <- region_list[[region_name]][[modality]]
    dmr_gr <- region_dmr_to_granges(dmr_df)

    at_lost <- unique(queryHits(findOverlaps(dmr_gr, lost_gr)))
    at_gained <- unique(queryHits(findOverlaps(dmr_gr, gained_gr)))

    lost_subset <- dmr_df[at_lost, ]
    gained_subset <- dmr_df[at_gained, ]

    if (modality == "mc") {
      pos_direction <- "Hypermethylated"
      test_label <- "mC_hyper"
    } else {
      pos_direction <- "Hypomethylated"
      test_label <- "hmC_hypo"
    }

    n_lost_pos <- sum(lost_subset$significant & lost_subset$direction == pos_direction)
    n_gained_pos <- sum(gained_subset$significant & gained_subset$direction == pos_direction)

    if (nrow(lost_subset) >= 5 & nrow(gained_subset) >= 5) {
      ft <- run_fisher_2x2(n_lost_pos, nrow(lost_subset),
                            n_gained_pos, nrow(gained_subset),
                            "Lost", "Gained")
      ft$region <- region_name
      ft$modality <- modality
      ft$test_type <- test_label
      fisher_47b[[idx]] <- ft

      wt <- wilcox.test(lost_subset$mod_difference, gained_subset$mod_difference)
      wilcox_47b[[idx]] <- tibble(
        region = region_name, modality = modality,
        lost_n = nrow(lost_subset), gained_n = nrow(gained_subset),
        lost_median = median(lost_subset$mod_difference),
        gained_median = median(gained_subset$mod_difference),
        W = as.numeric(wt$statistic), p_value = wt$p.value
      )

      cat(sprintf("  %s %s: lost %d/%d (%.1f%%) vs gained %d/%d (%.1f%%), OR=%.2f, %s\n",
                  region_name, modality,
                  n_lost_pos, nrow(lost_subset),
                  100 * n_lost_pos / nrow(lost_subset),
                  n_gained_pos, nrow(gained_subset),
                  100 * n_gained_pos / nrow(gained_subset),
                  ft$odds_ratio, fmt_p(ft$p_value)))
    } else {
      cat(sprintf("  %s %s: skipped (lost n=%d, gained n=%d)\n",
                  region_name, modality, nrow(lost_subset), nrow(gained_subset)))
    }
    idx <- idx + 1
  }
}

fisher_47b_df <- bind_rows(fisher_47b)
fisher_47b_df$fdr <- p.adjust(fisher_47b_df$p_value, method = "BH")
write.table(fisher_47b_df, file.path(TABLES_DIR, "47b_multiregion_fisher.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wilcox_47b_df <- bind_rows(wilcox_47b)
wilcox_47b_df$fdr <- p.adjust(wilcox_47b_df$p_value, method = "BH")
write.table(wilcox_47b_df, file.path(TABLES_DIR, "47b_multiregion_wilcoxon.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Plot 47b.1: Forest plot of ORs
cat("\n  Generating plots...\n")

forest_data <- fisher_47b_df %>%
  mutate(
    log2_or = log2(odds_ratio),
    log2_ci_lower = log2(ci_lower),
    log2_ci_upper = log2(ci_upper),
    sig_label = ifelse(fdr < 0.05, "*", ""),
    modality_label = ifelse(modality == "mc", "5mC (hyper)", "5hmC (hypo)"),
    label = paste(region, modality_label, sep = " — ")
  ) %>%
  mutate(label = factor(label, levels = rev(label)))

p47b_forest <- ggplot(forest_data, aes(x = log2_or, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = log2_ci_lower, xmax = log2_ci_upper),
                 height = 0.25, color = "grey30") +
  geom_point(aes(color = fdr < 0.05, size = -log10(fdr)), shape = 16) +
  scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60"),
                     name = "FDR < 0.05") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 6)) +
  labs(x = "log2(Odds Ratio) — Lost vs Gained CTCF Anchors",
       y = NULL,
       title = "Methylation Enrichment at Lost CTCF Anchors",
       subtitle = "Fisher's exact test across region types") +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p47b_forest,
                        file.path(SECTION_DIR, "47b_multiregion_OR_forest"),
                        width = 10, height = 7)

# Plot 47b.2: Heatmap of % hyper/hypo
heatmap_data <- fisher_47b_df %>%
  dplyr::select(region, modality, test_pct, ref_pct) %>%
  pivot_longer(cols = c(test_pct, ref_pct),
               names_to = "anchor", values_to = "pct") %>%
  mutate(anchor = ifelse(anchor == "test_pct", "Lost", "Gained"),
         modality_label = ifelse(modality == "mc", "5mC hyper %", "5hmC hypo %"))

# Add background rates
for (region_name in names(region_list)) {
  for (modality in c("mc", "hmc")) {
    dmr_df <- region_list[[region_name]][[modality]]
    dmr_gr <- region_dmr_to_granges(dmr_df)
    all_anchor_idx <- union(
      unique(queryHits(findOverlaps(dmr_gr, lost_gr))),
      unique(queryHits(findOverlaps(dmr_gr, gained_gr)))
    )
    bg_df <- dmr_df[-all_anchor_idx, ]
    pos_dir <- ifelse(modality == "mc", "Hypermethylated", "Hypomethylated")
    bg_pct <- 100 * sum(bg_df$significant & bg_df$direction == pos_dir) / nrow(bg_df)
    heatmap_data <- bind_rows(heatmap_data, tibble(
      region = region_name,
      modality = modality,
      anchor = "Background",
      pct = bg_pct,
      modality_label = ifelse(modality == "mc", "5mC hyper %", "5hmC hypo %")
    ))
  }
}

heatmap_data <- heatmap_data %>%
  mutate(anchor = factor(anchor, levels = c("Lost", "Gained", "Background")))

p47b_heatmap <- ggplot(heatmap_data, aes(x = anchor, y = region, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), size = 3.5) +
  facet_wrap(~modality_label) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = median(heatmap_data$pct, na.rm = TRUE),
                       name = "% Sig.") +
  labs(x = NULL, y = NULL,
       title = "Significant DMR Rates at CTCF Anchors by Region Type") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_multiformat_ggplot(p47b_heatmap,
                        file.path(SECTION_DIR, "47b_multiregion_pct_heatmap"),
                        width = 9, height = 6)


# =============================================================================
# 47c: COORDINATED mC-UP / hmC-DOWN AT CTCF ANCHOR CpG ISLANDS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("47c: COORDINATED mC-UP / hmC-DOWN AT CTCF ANCHOR CpG ISLANDS\n")
cat("================================================================================\n\n")

cpgi_coord <- inner_join(
  cpgi_mc %>%
    dplyr::select(chr, start, end, mod_difference, dmr_qvalue, significant, direction, anchor_group) %>%
    dplyr::rename(mc_diff = mod_difference, mc_q = dmr_qvalue, mc_sig = significant, mc_dir = direction),
  cpgi_hmc %>%
    dplyr::select(chr, start, end, mod_difference, dmr_qvalue, significant, direction) %>%
    dplyr::rename(hmc_diff = mod_difference, hmc_q = dmr_qvalue, hmc_sig = significant, hmc_dir = direction),
  by = c("chr", "start", "end")
)

cpgi_coord <- cpgi_coord %>%
  mutate(
    coordinated = mc_diff > 0 & hmc_diff < 0,
    coordinated_sig = coordinated & mc_sig
  )

cat(sprintf("  Joined CpG islands: %d\n", nrow(cpgi_coord)))
cat(sprintf("  Coordinated (mC up + hmC down): %d (%.1f%%)\n",
            sum(cpgi_coord$coordinated),
            100 * sum(cpgi_coord$coordinated) / nrow(cpgi_coord)))

# Coordinated rates by anchor group
coord_rates <- cpgi_coord %>%
  filter(anchor_group != "Both") %>%
  group_by(anchor_group) %>%
  summarise(
    n = n(),
    n_coordinated = sum(coordinated),
    pct_coordinated = 100 * n_coordinated / n,
    n_coordinated_sig = sum(coordinated_sig),
    pct_coordinated_sig = 100 * n_coordinated_sig / n,
    .groups = "drop"
  )

cat("\n  Coordinated rates by anchor group:\n")
for (i in 1:nrow(coord_rates)) {
  cat(sprintf("    %s: %d/%d (%.1f%%) coordinated, %d sig-coordinated\n",
              coord_rates$anchor_group[i],
              coord_rates$n_coordinated[i], coord_rates$n[i],
              coord_rates$pct_coordinated[i],
              coord_rates$n_coordinated_sig[i]))
}

# Fisher's: coordinated at lost vs gained
lost_coord <- cpgi_coord %>% filter(anchor_group == "Lost CTCF anchor")
gained_coord <- cpgi_coord %>% filter(anchor_group == "Gained CTCF anchor")
bg_coord <- cpgi_coord %>% filter(anchor_group == "Background")

fisher_47c <- list()

fisher_47c[[1]] <- run_fisher_2x2(
  sum(lost_coord$coordinated), nrow(lost_coord),
  sum(gained_coord$coordinated), nrow(gained_coord),
  "Lost", "Gained")
fisher_47c[[1]]$test <- "coordinated_lost_vs_gained"

fisher_47c[[2]] <- run_fisher_2x2(
  sum(lost_coord$coordinated), nrow(lost_coord),
  sum(bg_coord$coordinated), nrow(bg_coord),
  "Lost", "Background")
fisher_47c[[2]]$test <- "coordinated_lost_vs_background"

fisher_47c[[3]] <- run_fisher_2x2(
  sum(lost_coord$coordinated_sig), nrow(lost_coord),
  sum(gained_coord$coordinated_sig), nrow(gained_coord),
  "Lost", "Gained")
fisher_47c[[3]]$test <- "coordinated_sig_lost_vs_gained"

fisher_47c_df <- bind_rows(fisher_47c)
fisher_47c_df$fdr <- p.adjust(fisher_47c_df$p_value, method = "BH")

cat(sprintf("\n  Fisher's (coordinated, lost vs gained): OR=%.2f, %s\n",
            fisher_47c[[1]]$odds_ratio, fmt_p(fisher_47c[[1]]$p_value)))
cat(sprintf("  Fisher's (coordinated, lost vs BG): OR=%.2f, %s\n",
            fisher_47c[[2]]$odds_ratio, fmt_p(fisher_47c[[2]]$p_value)))

write.table(cpgi_coord, file.path(TABLES_DIR, "47c_cpgi_coordinated_pattern.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fisher_47c_df, file.path(TABLES_DIR, "47c_coordinated_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plots ---
cat("\n  Generating plots...\n")

# Plot 47c.1: Coordinated % bar chart
coord_bar_data <- coord_rates %>%
  mutate(anchor_group = factor(anchor_group,
                               levels = c("Lost CTCF anchor", "Gained CTCF anchor", "Background")),
         label = sprintf("%d/%d", n_coordinated, n))

p47c_bar <- ggplot(coord_bar_data, aes(x = anchor_group, y = pct_coordinated, fill = anchor_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = ANCHOR_COLORS, guide = "none") +
  labs(x = NULL,
       y = "% CpG Islands with Coordinated mC↑ / hmC↓",
       title = "Coordinated Methylation Pattern at CTCF Anchor CpG Islands",
       subtitle = sprintf("Fisher's OR (lost vs gained) = %.2f, %s",
                          fisher_47c[[1]]$odds_ratio,
                          fmt_p(fisher_47c[[1]]$p_value))) +
  theme_biomodal() +
  coord_cartesian(ylim = c(0, max(coord_bar_data$pct_coordinated) * 1.2))

save_multiformat_ggplot(p47c_bar,
                        file.path(SECTION_DIR, "47c_coordinated_barchart"),
                        width = 7, height = 6)

# Plot 47c.2: mC vs hmC scatter
scatter_data <- cpgi_coord %>%
  filter(anchor_group != "Both") %>%
  mutate(anchor_group = factor(anchor_group,
                               levels = c("Lost CTCF anchor", "Gained CTCF anchor", "Background")))

p47c_scatter <- ggplot(scatter_data, aes(x = mc_diff, y = hmc_diff, color = anchor_group)) +
  geom_point(data = scatter_data %>% filter(anchor_group == "Background"),
             alpha = 0.15, size = 0.8) +
  geom_point(data = scatter_data %>% filter(anchor_group != "Background"),
             alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = ANCHOR_COLORS, name = "Anchor Group") +
  labs(x = "mC Modification Difference (mutant - control)",
       y = "hmC Modification Difference (mutant - control)",
       title = "mC vs hmC at CpG Islands by CTCF Anchor Group",
       subtitle = "Upper-left quadrant = coordinated mC↑ / hmC↓ (TET blockade)") +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p47c_scatter,
                        file.path(SECTION_DIR, "47c_scatter_mc_vs_hmc"),
                        width = 9, height = 8)

# Plot 47c.3: Combined effect size violin
cpgi_coord_plot <- cpgi_coord %>%
  filter(anchor_group != "Both") %>%
  mutate(
    combined_effect = abs(mc_diff) + abs(hmc_diff),
    anchor_group = factor(anchor_group,
                          levels = c("Lost CTCF anchor", "Gained CTCF anchor", "Background"))
  )

wt_combined <- wilcox.test(
  cpgi_coord_plot$combined_effect[cpgi_coord_plot$anchor_group == "Lost CTCF anchor"],
  cpgi_coord_plot$combined_effect[cpgi_coord_plot$anchor_group == "Gained CTCF anchor"]
)

p47c_combined <- ggplot(cpgi_coord_plot, aes(x = anchor_group, y = combined_effect, fill = anchor_group)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = ANCHOR_COLORS, guide = "none") +
  labs(x = NULL,
       y = "|mC diff| + |hmC diff|",
       title = "Combined Methylation Effect Size at CpG Islands",
       subtitle = sprintf("Wilcoxon lost vs gained: %s", fmt_p(wt_combined$p.value))) +
  theme_biomodal()

save_multiformat_ggplot(p47c_combined,
                        file.path(SECTION_DIR, "47c_combined_effect_violin"),
                        width = 7, height = 6)


# =============================================================================
# 47d: DISTANCE-STRATIFIED ANALYSIS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("47d: DISTANCE-STRATIFIED ANALYSIS\n")
cat("================================================================================\n\n")

# Assign distance bins to anchors
anchors_lost$distance_bin <- cut(anchors_lost$loop_distance,
                                 breaks = DISTANCE_BINS,
                                 labels = DISTANCE_LABELS,
                                 include.lowest = TRUE)
anchors_gained$distance_bin <- cut(anchors_gained$loop_distance,
                                   breaks = DISTANCE_BINS,
                                   labels = DISTANCE_LABELS,
                                   include.lowest = TRUE)

# For each CpG island at an anchor, inherit the anchor's distance bin
# Use findOverlaps to build the mapping
cpgi_mc_gr_47d <- region_dmr_to_granges(cpgi_mc)

hits_lost_47d <- findOverlaps(cpgi_mc_gr_47d, lost_gr)
hits_gained_47d <- findOverlaps(cpgi_mc_gr_47d, gained_gr)

# Build per-CpG-island table with distance bin from the overlapping anchor
cpgi_lost_dist <- tibble(
  cpgi_idx = queryHits(hits_lost_47d),
  anchor_idx = subjectHits(hits_lost_47d),
  direction = "lost",
  distance_bin = anchors_lost$distance_bin[subjectHits(hits_lost_47d)]
) %>%
  group_by(cpgi_idx) %>%
  slice_min(order_by = anchors_lost$loop_distance[anchor_idx], n = 1, with_ties = FALSE) %>%
  ungroup()

cpgi_gained_dist <- tibble(
  cpgi_idx = queryHits(hits_gained_47d),
  anchor_idx = subjectHits(hits_gained_47d),
  direction = "gained",
  distance_bin = anchors_gained$distance_bin[subjectHits(hits_gained_47d)]
) %>%
  group_by(cpgi_idx) %>%
  slice_min(order_by = anchors_gained$loop_distance[anchor_idx], n = 1, with_ties = FALSE) %>%
  ungroup()

# Merge with DMR data
cpgi_lost_dist$mc_hyper <- cpgi_mc$significant[cpgi_lost_dist$cpgi_idx] &
  cpgi_mc$direction[cpgi_lost_dist$cpgi_idx] == "Hypermethylated"
cpgi_gained_dist$mc_hyper <- cpgi_mc$significant[cpgi_gained_dist$cpgi_idx] &
  cpgi_mc$direction[cpgi_gained_dist$cpgi_idx] == "Hypermethylated"
cpgi_lost_dist$mc_diff <- cpgi_mc$mod_difference[cpgi_lost_dist$cpgi_idx]
cpgi_gained_dist$mc_diff <- cpgi_mc$mod_difference[cpgi_gained_dist$cpgi_idx]

# Per-bin Fisher's tests
strat_results <- list()
strat_idx <- 1
cmh_arrays <- list()

cat("  Per-distance-bin Fisher's tests:\n")
for (bin in DISTANCE_LABELS) {
  lost_bin <- cpgi_lost_dist %>% filter(distance_bin == bin)
  gained_bin <- cpgi_gained_dist %>% filter(distance_bin == bin)

  if (nrow(lost_bin) >= 5 & nrow(gained_bin) >= 5) {
    n_lost_hyper <- sum(lost_bin$mc_hyper)
    n_gained_hyper <- sum(gained_bin$mc_hyper)

    ft <- run_fisher_2x2(n_lost_hyper, nrow(lost_bin),
                          n_gained_hyper, nrow(gained_bin),
                          paste0("lost_", bin), paste0("gained_", bin))
    ft$distance_bin <- bin
    strat_results[[strat_idx]] <- ft

    # Accumulate 2x2 tables for CMH test
    cmh_arrays[[strat_idx]] <- matrix(
      c(n_lost_hyper, nrow(lost_bin) - n_lost_hyper,
        n_gained_hyper, nrow(gained_bin) - n_gained_hyper),
      nrow = 2,
      dimnames = list(c("hyper", "not_hyper"), c("lost", "gained"))
    )

    cat(sprintf("    %s: lost %d/%d (%.1f%%) vs gained %d/%d (%.1f%%), OR=%.2f, %s\n",
                bin, n_lost_hyper, nrow(lost_bin),
                100 * n_lost_hyper / nrow(lost_bin),
                n_gained_hyper, nrow(gained_bin),
                100 * n_gained_hyper / nrow(gained_bin),
                ft$odds_ratio, fmt_p(ft$p_value)))

    strat_idx <- strat_idx + 1
  } else {
    cat(sprintf("    %s: skipped (lost n=%d, gained n=%d)\n",
                bin, nrow(lost_bin), nrow(gained_bin)))
  }
}

strat_df <- bind_rows(strat_results)
strat_df$fdr <- p.adjust(strat_df$p_value, method = "BH")

# Cochran-Mantel-Haenszel test
cmh_result <- NULL
if (length(cmh_arrays) >= 2) {
  cmh_array <- array(
    unlist(cmh_arrays),
    dim = c(2, 2, length(cmh_arrays)),
    dimnames = list(c("hyper", "not_hyper"), c("lost", "gained"),
                    names(cmh_arrays))
  )

  cmh_test <- tryCatch(
    mantelhaen.test(cmh_array, correct = TRUE),
    error = function(e) { cat("    CMH test failed:", e$message, "\n"); NULL }
  )

  if (!is.null(cmh_test)) {
    cmh_result <- tibble(
      test = "Cochran-Mantel-Haenszel",
      common_OR = as.numeric(cmh_test$estimate),
      ci_lower = cmh_test$conf.int[1],
      ci_upper = cmh_test$conf.int[2],
      chi_sq = as.numeric(cmh_test$statistic),
      p_value = cmh_test$p.value
    )
    cat(sprintf("\n  CMH test (all strata): common OR=%.2f [%.2f-%.2f], %s\n",
                cmh_result$common_OR, cmh_result$ci_lower, cmh_result$ci_upper,
                fmt_p(cmh_result$p_value)))
  }
}

# Save results
strat_output <- strat_df
if (!is.null(cmh_result)) {
  strat_output$cmh_OR <- cmh_result$common_OR
  strat_output$cmh_p <- cmh_result$p_value
}
write.table(strat_output, file.path(TABLES_DIR, "47d_distance_stratified_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plots ---
cat("\n  Generating plots...\n")

# Plot 47d.1: Forest plot per distance bin
if (nrow(strat_df) > 0) {
  forest_47d <- strat_df %>%
    mutate(
      log2_or = log2(odds_ratio),
      log2_ci_lower = log2(ci_lower),
      log2_ci_upper = log2(ci_upper),
      distance_bin = factor(distance_bin, levels = rev(DISTANCE_LABELS))
    )

  # Add CMH row if available
  if (!is.null(cmh_result)) {
    cmh_row <- tibble(
      distance_bin = factor("CMH Overall", levels = c("CMH Overall", rev(DISTANCE_LABELS))),
      log2_or = log2(cmh_result$common_OR),
      log2_ci_lower = log2(cmh_result$ci_lower),
      log2_ci_upper = log2(cmh_result$ci_upper),
      fdr = cmh_result$p_value
    )
    forest_47d$distance_bin <- factor(as.character(forest_47d$distance_bin),
                                      levels = c("CMH Overall", rev(DISTANCE_LABELS)))
    forest_47d <- bind_rows(forest_47d, cmh_row)
  }

  p47d_forest <- ggplot(forest_47d, aes(x = log2_or, y = distance_bin)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = log2_ci_lower, xmax = log2_ci_upper),
                   height = 0.25, color = "grey30") +
    geom_point(aes(color = fdr < 0.05), size = 4, shape = 16) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60"),
                       name = "FDR < 0.05") +
    labs(x = "log2(Odds Ratio) — mC Hyper at Lost vs Gained CTCF Anchors",
         y = "Loop Distance Bin",
         title = "Distance-Stratified CpG Island Hypermethylation",
         subtitle = "Controlling for loop length confound (CMH test at bottom)") +
    theme_biomodal()

  save_multiformat_ggplot(p47d_forest,
                          file.path(SECTION_DIR, "47d_distance_stratified_OR"),
                          width = 9, height = 6)
}

# Plot 47d.2: Faceted violin by distance bin
cpgi_dist_combined <- bind_rows(
  cpgi_lost_dist %>% mutate(group = "Lost"),
  cpgi_gained_dist %>% mutate(group = "Gained")
) %>%
  filter(!is.na(distance_bin)) %>%
  mutate(group = factor(group, levels = c("Lost", "Gained")))

if (nrow(cpgi_dist_combined) > 0) {
  p47d_violin <- ggplot(cpgi_dist_combined, aes(x = group, y = mc_diff, fill = group)) +
    geom_violin(alpha = 0.3, color = NA) +
    geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~distance_bin, nrow = 1) +
    scale_fill_manual(values = c("Lost" = "#d73027", "Gained" = "#4575b4"),
                      guide = "none") +
    labs(x = NULL,
         y = "mC Modification Difference",
         title = "mC at CpG Islands by Loop Distance Bin",
         subtitle = "CTCF anchor CpG islands, lost vs gained") +
    theme_biomodal() +
    theme(strip.text = element_text(size = 10))

  save_multiformat_ggplot(p47d_violin,
                          file.path(SECTION_DIR, "47d_mc_violin_by_distance"),
                          width = 12, height = 6)
}


# =============================================================================
# 47e: METHYLATION EFFECT SIZE vs LOOP logFC CORRELATION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("47e: METHYLATION EFFECT SIZE vs LOOP logFC CORRELATION\n")
cat("================================================================================\n\n")

# For each CTCF-CTCF loop, find CpG islands at either anchor
ctcf_ctcf_loops <- loops %>%
  filter(anchor1_CTCF_overlap == TRUE & anchor2_CTCF_overlap == TRUE)

cat(sprintf("  CTCF-CTCF loops: %d (%d lost, %d gained)\n",
            nrow(ctcf_ctcf_loops),
            sum(ctcf_ctcf_loops$direction == "down_in_mutant"),
            sum(ctcf_ctcf_loops$direction == "up_in_mutant")))

# Build anchor GRanges per loop
loop_meth <- lapply(seq_len(nrow(ctcf_ctcf_loops)), function(i) {
  row <- ctcf_ctcf_loops[i, ]

  a1_gr <- GRanges(seqnames = row$chr1, ranges = IRanges(start = row$start1, end = row$end1))
  a2_gr <- GRanges(seqnames = row$chr2, ranges = IRanges(start = row$start2, end = row$end2))
  both_gr <- c(a1_gr, a2_gr)

  mc_hits <- findOverlaps(cpgi_mc_gr, both_gr)
  hmc_hits <- findOverlaps(cpgi_hmc_gr, both_gr)

  mc_idx <- unique(queryHits(mc_hits))
  hmc_idx <- unique(queryHits(hmc_hits))

  tibble(
    loop_id = row$loop_id,
    direction = row$direction,
    logFC = row$logFC,
    loop_distance = row$loop_distance,
    n_cpgi = length(mc_idx),
    mean_mc_diff = if (length(mc_idx) > 0) mean(cpgi_mc$mod_difference[mc_idx]) else NA_real_,
    mean_hmc_diff = if (length(hmc_idx) > 0) mean(cpgi_hmc$mod_difference[hmc_idx]) else NA_real_,
    any_mc_hyper = any(cpgi_mc$significant[mc_idx] & cpgi_mc$direction[mc_idx] == "Hypermethylated"),
    any_hmc_hypo = any(cpgi_hmc$significant[hmc_idx] & cpgi_hmc$direction[hmc_idx] == "Hypomethylated")
  )
})

loop_meth_df <- bind_rows(loop_meth)
loop_meth_with_cpgi <- loop_meth_df %>% filter(!is.na(mean_mc_diff) & n_cpgi > 0)

cat(sprintf("  Loops with CpG islands at anchors: %d / %d (%.1f%%)\n",
            nrow(loop_meth_with_cpgi), nrow(ctcf_ctcf_loops),
            100 * nrow(loop_meth_with_cpgi) / nrow(ctcf_ctcf_loops)))

write.table(loop_meth_df, file.path(TABLES_DIR, "47e_loop_methylation_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Spearman correlations
corr_results <- list()

if (nrow(loop_meth_with_cpgi) >= 10) {
  # Overall mC
  mc_cor <- cor.test(loop_meth_with_cpgi$mean_mc_diff, loop_meth_with_cpgi$logFC,
                     method = "spearman", exact = FALSE)
  corr_results[[1]] <- tibble(
    test = "mC_vs_logFC_all", rho = mc_cor$estimate,
    p_value = mc_cor$p.value, n = nrow(loop_meth_with_cpgi)
  )
  cat(sprintf("  Spearman (mC diff vs logFC, all): rho=%.3f, %s, n=%d\n",
              mc_cor$estimate, fmt_p(mc_cor$p.value), nrow(loop_meth_with_cpgi)))

  # Lost loops only
  lost_loops <- loop_meth_with_cpgi %>% filter(direction == "down_in_mutant")
  if (nrow(lost_loops) >= 10) {
    mc_cor_lost <- cor.test(lost_loops$mean_mc_diff, lost_loops$logFC,
                            method = "spearman", exact = FALSE)
    corr_results[[length(corr_results) + 1]] <- tibble(
      test = "mC_vs_logFC_lost", rho = mc_cor_lost$estimate,
      p_value = mc_cor_lost$p.value, n = nrow(lost_loops)
    )
    cat(sprintf("  Spearman (mC diff vs logFC, lost only): rho=%.3f, %s, n=%d\n",
                mc_cor_lost$estimate, fmt_p(mc_cor_lost$p.value), nrow(lost_loops)))
  }

  # hmC
  hmc_with <- loop_meth_with_cpgi %>% filter(!is.na(mean_hmc_diff))
  if (nrow(hmc_with) >= 10) {
    hmc_cor <- cor.test(hmc_with$mean_hmc_diff, hmc_with$logFC,
                        method = "spearman", exact = FALSE)
    corr_results[[length(corr_results) + 1]] <- tibble(
      test = "hmC_vs_logFC_all", rho = hmc_cor$estimate,
      p_value = hmc_cor$p.value, n = nrow(hmc_with)
    )
    cat(sprintf("  Spearman (hmC diff vs logFC, all): rho=%.3f, %s, n=%d\n",
                hmc_cor$estimate, fmt_p(hmc_cor$p.value), nrow(hmc_with)))
  }

  # Partial correlation controlling for distance
  if (nrow(loop_meth_with_cpgi) >= 20) {
    partial_fit <- lm(logFC ~ mean_mc_diff + log10(loop_distance),
                      data = loop_meth_with_cpgi)
    partial_summary <- summary(partial_fit)
    mc_coef <- partial_summary$coefficients["mean_mc_diff", ]
    corr_results[[length(corr_results) + 1]] <- tibble(
      test = "partial_mC_vs_logFC_adj_distance",
      rho = mc_coef["Estimate"],
      p_value = mc_coef["Pr(>|t|)"],
      n = nrow(loop_meth_with_cpgi)
    )
    cat(sprintf("  Partial (mC diff | distance): beta=%.3f, %s\n",
                mc_coef["Estimate"], fmt_p(mc_coef["Pr(>|t|)"])))
  }
}

corr_results_df <- bind_rows(corr_results)
write.table(corr_results_df, file.path(TABLES_DIR, "47e_correlation_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plots ---
cat("\n  Generating plots...\n")

if (nrow(loop_meth_with_cpgi) >= 10) {
  plot_df <- loop_meth_with_cpgi %>%
    mutate(dir_label = ifelse(direction == "down_in_mutant", "Lost", "Gained"),
           dir_label = factor(dir_label, levels = c("Lost", "Gained")))

  # Plot 47e.1: logFC vs mC scatter
  rho_label <- sprintf("rho = %.3f, %s",
                        corr_results_df$rho[corr_results_df$test == "mC_vs_logFC_all"],
                        fmt_p(corr_results_df$p_value[corr_results_df$test == "mC_vs_logFC_all"]))

  p47e_mc <- ggplot(plot_df, aes(x = mean_mc_diff, y = logFC, color = dir_label)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Lost" = "#d73027", "Gained" = "#4575b4"),
                       name = "Loop Direction") +
    labs(x = "Mean mC Difference at Anchor CpG Islands",
         y = "Loop logFC (mutant / control)",
         title = "Loop Strength vs CpG Island Methylation at CTCF Anchors",
         subtitle = rho_label) +
    theme_biomodal() +
    theme(legend.position = "top")

  save_multiformat_ggplot(p47e_mc,
                          file.path(SECTION_DIR, "47e_logfc_vs_mc_scatter"),
                          width = 8, height = 7)

  # Plot 47e.2: logFC vs hmC scatter
  hmc_plot_df <- plot_df %>% filter(!is.na(mean_hmc_diff))

  if (nrow(hmc_plot_df) >= 10) {
    hmc_rho_row <- corr_results_df %>% filter(test == "hmC_vs_logFC_all")
    hmc_rho_label <- if (nrow(hmc_rho_row) > 0) {
      sprintf("rho = %.3f, %s", hmc_rho_row$rho, fmt_p(hmc_rho_row$p_value))
    } else {
      "hmC correlation not computed"
    }

    p47e_hmc <- ggplot(hmc_plot_df, aes(x = mean_hmc_diff, y = logFC, color = dir_label)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 1, alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = c("Lost" = "#d73027", "Gained" = "#4575b4"),
                         name = "Loop Direction") +
      labs(x = "Mean hmC Difference at Anchor CpG Islands",
           y = "Loop logFC (mutant / control)",
           title = "Loop Strength vs CpG Island Hydroxymethylation at CTCF Anchors",
           subtitle = hmc_rho_label) +
      theme_biomodal() +
      theme(legend.position = "top")

    save_multiformat_ggplot(p47e_hmc,
                            file.path(SECTION_DIR, "47e_logfc_vs_hmc_scatter"),
                            width = 8, height = 7)
  }

  # Plot 47e.3: Partial regression residual plot
  if (nrow(loop_meth_with_cpgi) >= 20) {
    resid_mc <- residuals(lm(mean_mc_diff ~ log10(loop_distance),
                             data = loop_meth_with_cpgi))
    resid_logfc <- residuals(lm(logFC ~ log10(loop_distance),
                                data = loop_meth_with_cpgi))

    partial_df <- tibble(
      resid_mc = resid_mc,
      resid_logfc = resid_logfc,
      dir_label = ifelse(loop_meth_with_cpgi$direction == "down_in_mutant", "Lost", "Gained")
    ) %>%
      mutate(dir_label = factor(dir_label, levels = c("Lost", "Gained")))

    partial_cor <- cor.test(resid_mc, resid_logfc, method = "spearman", exact = FALSE)

    p47e_partial <- ggplot(partial_df, aes(x = resid_mc, y = resid_logfc, color = dir_label)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1, alpha = 0.2, color = "grey30") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = c("Lost" = "#d73027", "Gained" = "#4575b4"),
                         name = "Loop Direction") +
      labs(x = "mC Difference (residual, distance-adjusted)",
           y = "logFC (residual, distance-adjusted)",
           title = "Partial Correlation: Methylation vs Loop Strength",
           subtitle = sprintf("Distance-adjusted Spearman rho = %.3f, %s",
                              partial_cor$estimate, fmt_p(partial_cor$p.value))) +
      theme_biomodal() +
      theme(legend.position = "top")

    save_multiformat_ggplot(p47e_partial,
                            file.path(SECTION_DIR, "47e_residual_partial_correlation"),
                            width = 7, height = 6)
  }
}


# =============================================================================
# SECTION 47 SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 47 SUMMARY\n")
cat("================================================================================\n\n")

cat("47a — CpG Island Methylation at CTCF Anchors (Core Test):\n")
cat(sprintf("  CpG islands at lost CTCF anchors: %d\n", length(cpgi_at_lost_idx)))
cat(sprintf("  CpG islands at gained CTCF anchors: %d\n", length(cpgi_at_gained_idx)))
cat(sprintf("  mC hyper Fisher's (lost vs gained): OR=%.2f, %s\n",
            fisher_47a_df$odds_ratio[1], fmt_p(fisher_47a_df$p_value[1])))
cat(sprintf("  hmC hypo Fisher's (lost vs gained): OR=%.2f, %s\n",
            fisher_47a_df$odds_ratio[3], fmt_p(fisher_47a_df$p_value[3])))
cat(sprintf("  mC Wilcoxon: %s (lost median=%.4f, gained median=%.4f)\n",
            fmt_p(wt_mc$p.value), median(lost_mc$mod_difference),
            median(gained_mc$mod_difference)))

cat("\n47b — Multi-Region Comparison:\n")
for (i in 1:nrow(fisher_47b_df)) {
  cat(sprintf("  %s %s: OR=%.2f, FDR=%.2e\n",
              fisher_47b_df$region[i], fisher_47b_df$modality[i],
              fisher_47b_df$odds_ratio[i], fisher_47b_df$fdr[i]))
}

cat("\n47c — Coordinated Pattern:\n")
cat(sprintf("  Fisher's (coordinated, lost vs gained): OR=%.2f, %s\n",
            fisher_47c_df$odds_ratio[1], fmt_p(fisher_47c_df$p_value[1])))

cat("\n47d — Distance-Stratified:\n")
if (!is.null(cmh_result)) {
  cat(sprintf("  CMH common OR=%.2f [%.2f-%.2f], %s\n",
              cmh_result$common_OR, cmh_result$ci_lower, cmh_result$ci_upper,
              fmt_p(cmh_result$p_value)))
} else {
  cat("  CMH test not computed (insufficient strata)\n")
}

cat("\n47e — logFC vs Methylation Correlation:\n")
if (nrow(corr_results_df) > 0) {
  for (i in 1:nrow(corr_results_df)) {
    cat(sprintf("  %s: rho/beta=%.3f, %s (n=%d)\n",
                corr_results_df$test[i], corr_results_df$rho[i],
                fmt_p(corr_results_df$p_value[i]), corr_results_df$n[i]))
  }
}

cat("\n--- Output files ---\n")
cat(sprintf("Figures (13 panels in %s):\n", SECTION_DIR))
cat("  47a_cpgi_mc_direction_barchart/\n")
cat("  47a_cpgi_mc_violin/\n")
cat("  47a_cpgi_hmc_violin/\n")
cat("  47b_multiregion_OR_forest/\n")
cat("  47b_multiregion_pct_heatmap/\n")
cat("  47c_coordinated_barchart/\n")
cat("  47c_scatter_mc_vs_hmc/\n")
cat("  47c_combined_effect_violin/\n")
cat("  47d_distance_stratified_OR/\n")
cat("  47d_mc_violin_by_distance/\n")
cat("  47e_logfc_vs_mc_scatter/\n")
cat("  47e_logfc_vs_hmc_scatter/\n")
cat("  47e_residual_partial_correlation/\n")
cat(sprintf("\nTables (in %s):\n", TABLES_DIR))
cat("  47a_cpgi_ctcf_anchor_methylation.tsv\n")
cat("  47a_fisher_results.tsv\n")
cat("  47b_multiregion_fisher.tsv\n")
cat("  47b_multiregion_wilcoxon.tsv\n")
cat("  47c_cpgi_coordinated_pattern.tsv\n")
cat("  47c_coordinated_enrichment.tsv\n")
cat("  47d_distance_stratified_results.tsv\n")
cat("  47e_loop_methylation_correlation.tsv\n")
cat("  47e_correlation_results.tsv\n")

cat("\n================================================================================\n")
cat("SECTION 47 COMPLETE\n")
cat("================================================================================\n\n")
