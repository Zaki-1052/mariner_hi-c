# biomodal/downstream/scripts/viz_sections/section_27_methylation_hic_loop_anchor_integration.R
# Section 27: Methylation x Hi-C Loop Anchor Integration
# Tests whether genes at differential loop anchors are enriched for the
# coordinated mC-up/hmC-down pattern, and whether methylation direction
# associates with loop direction (lost vs gained).
#
# Central premise: BAP1-KO -> H2AK119ub accumulation -> TET blockade ->
#   mC up / hmC down at genes whose 3D contacts are disrupted.
#   Loop anchors mark regulatory hubs, so anchor-associated genes should
#   show stronger coordinated methylation changes.
#
# Sub-analyses:
#   27a: Coordinated gene enrichment at loop anchors
#   27b: Stratification by loop direction (lost vs gained)
#   27c: K119ub-loop-methylation convergence
#   27d: Logistic regression (loop features -> hypermethylation)
#   27e: Shared anchor methylation profile
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_27_methylation_hic_loop_anchor_integration.R

source("scripts/viz_sections/_shared_config.R")

# ChIPseeker for K119ub gene annotation in 27c
suppressPackageStartupMessages(library(ChIPseeker))

cat("\n")
cat("================================================================================\n")
cat("SECTION 27: METHYLATION x HI-C LOOP ANCHOR INTEGRATION\n")
cat("================================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input file paths
LOOPS_FILE <- file.path(BASE_DIR, "../../outputs/250402-late_outputs/merged_loops/characterized_loops.tsv")
RATIO_FILE <- file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv")
SHARED_ANCHORS_FILE <- file.path(BASE_DIR, "../../output/shared_anchor_analysis/late/tables/shared_anchors.tsv")
SHARED_ANCHOR_GENES_FILE <- file.path(BASE_DIR, "../../output/shared_anchor_analysis/late/tables/shared_anchor_genes.tsv")

# Validate all inputs exist
stopifnot("Loops file not found" = file.exists(LOOPS_FILE))
stopifnot("Delta-ratio file not found" = file.exists(RATIO_FILE))
stopifnot("Shared anchors file not found" = file.exists(SHARED_ANCHORS_FILE))
stopifnot("Shared anchor genes file not found" = file.exists(SHARED_ANCHOR_GENES_FILE))
stopifnot("K119ub ctrl peaks not found" = file.exists(K119UB_FILES$ctrl))
stopifnot("K119ub mut peaks not found" = file.exists(K119UB_FILES$mut))

# GREAT-style regulatory domain parameters
GREAT_UPSTREAM <- 5000       # 5kb upstream of TSS
GREAT_DOWNSTREAM <- 1000     # 1kb downstream of TSS
GREAT_MAX_EXTENSION <- 100000  # 100kb max extension

# Colors
LOOP_DIRECTION_COLORS <- c(Lost = "#d73027", Gained = "#4575b4", Background = "grey70")

# Section output directory
SECTION_DIR <- file.path(OUTPUT_DIR, "27_methylation_loop_integration")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Build GREAT-style regulatory domains from TxDb
#' Basal: 5kb upstream, 1kb downstream (strand-aware)
#' Extended: up to 100kb, stopping at neighboring gene's basal domain
#' @return data.frame with entrez_id, chr, tss, strand, reg_start, reg_end
build_gene_regulatory_domains <- function() {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  genes_gr <- genes(txdb)

  tss_pos <- ifelse(as.character(strand(genes_gr)) == "-",
                    end(genes_gr), start(genes_gr))

  gene_info <- data.frame(
    entrez_id = names(genes_gr),
    chr = as.character(seqnames(genes_gr)),
    tss = tss_pos,
    strand = as.character(strand(genes_gr)),
    stringsAsFactors = FALSE
  )

  gene_info <- gene_info %>%
    dplyr::arrange(chr, tss) %>%
    dplyr::mutate(
      basal_start = ifelse(strand == "+", tss - GREAT_UPSTREAM, tss - GREAT_DOWNSTREAM),
      basal_end = ifelse(strand == "+", tss + GREAT_DOWNSTREAM, tss + GREAT_UPSTREAM),
      max_start = tss - GREAT_MAX_EXTENSION,
      max_end = tss + GREAT_MAX_EXTENSION
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      prev_basal_end = dplyr::lag(basal_end, default = -Inf),
      next_basal_start = dplyr::lead(basal_start, default = Inf),
      reg_start = pmax(max_start, prev_basal_end, 1),
      reg_end = pmin(max_end, next_basal_start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      reg_start = ifelse(reg_end < reg_start, basal_start, reg_start),
      reg_end = ifelse(reg_end < reg_start, basal_end, reg_end),
      reg_end = pmax(reg_end, reg_start)
    ) %>%
    dplyr::select(entrez_id, chr, tss, strand, reg_start, reg_end)

  return(gene_info)
}

#' Associate genomic anchors (GRanges) with genes via GREAT-style regulatory domains
#' @param anchor_gr GRanges of anchors (with anchor_id, loop_direction metadata)
#' @param gene_domains_gr GRanges of gene regulatory domains (with entrez_id metadata)
#' @return data.frame of anchor-gene associations
associate_anchors_to_genes <- function(anchor_gr, gene_domains_gr) {
  overlaps <- findOverlaps(anchor_gr, gene_domains_gr, ignore.strand = TRUE)

  if (length(overlaps) == 0) {
    stop("No anchor-gene overlaps found. Check coordinate systems.")
  }

  anchor_gene_df <- data.frame(
    anchor_idx = queryHits(overlaps),
    entrez_id = gene_domains_gr$entrez_id[subjectHits(overlaps)],
    anchor_chr = as.character(seqnames(anchor_gr)[queryHits(overlaps)]),
    anchor_start = start(anchor_gr)[queryHits(overlaps)],
    anchor_end = end(anchor_gr)[queryHits(overlaps)],
    stringsAsFactors = FALSE
  )

  # Carry anchor metadata
  if ("anchor_id" %in% names(mcols(anchor_gr))) {
    anchor_gene_df$anchor_id <- anchor_gr$anchor_id[queryHits(overlaps)]
  }
  if ("loop_direction" %in% names(mcols(anchor_gr))) {
    anchor_gene_df$loop_direction <- anchor_gr$loop_direction[queryHits(overlaps)]
  }
  if ("loop_distance" %in% names(mcols(anchor_gr))) {
    anchor_gene_df$loop_distance <- anchor_gr$loop_distance[queryHits(overlaps)]
  }
  if ("anchor_type" %in% names(mcols(anchor_gr))) {
    anchor_gene_df$anchor_type <- anchor_gr$anchor_type[queryHits(overlaps)]
  }

  return(anchor_gene_df)
}

#' Derive differential peaks (gained/lost) from condition-specific peak sets
#' Reused from section_14
derive_differential_peaks <- function(ctrl_gr, mut_gr) {
  mut_hits <- countOverlaps(mut_gr, ctrl_gr)
  ctrl_hits <- countOverlaps(ctrl_gr, mut_gr)
  list(
    gained = mut_gr[mut_hits == 0],
    lost = ctrl_gr[ctrl_hits == 0],
    shared = mut_gr[mut_hits > 0]
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n")

# 1. Differential loops
loops <- read.table(LOOPS_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Loops: %d total\n", nrow(loops)))

# 2. Delta-ratio per gene
ratio_df <- read.table(RATIO_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Delta-ratio genes: %d\n", nrow(ratio_df)))

# 3. Shared anchors
shared_anchors <- read.table(SHARED_ANCHORS_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
shared_anchor_genes <- read.table(SHARED_ANCHOR_GENES_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Shared anchors: %d, shared anchor genes: %d\n",
            nrow(shared_anchors), nrow(shared_anchor_genes)))

# 4. K119ub peaks
cat("  Loading K119ub peaks...\n")
k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub ctrl")
k119ub_mut_gr <- load_chip_peaks(K119UB_FILES$mut, "K119ub mut")
k119ub_diff <- derive_differential_peaks(k119ub_ctrl_gr, k119ub_mut_gr)
cat(sprintf("    K119ub: %d gained, %d lost, %d shared\n",
            length(k119ub_diff$gained), length(k119ub_diff$lost), length(k119ub_diff$shared)))

# =============================================================================
# ENTREZ -> SYMBOL MAPPING (centralized, computed once)
# =============================================================================

cat("\n--- Building Entrez -> Symbol mapping ---\n")

# All unique Entrez IDs from loop anchors
all_entrez <- unique(c(as.character(loops$anchor1_nearest_gene),
                       as.character(loops$anchor2_nearest_gene)))
all_entrez <- all_entrez[!is.na(all_entrez) & all_entrez != ""]

entrez_to_symbol <- tryCatch({
  mapping <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = all_entrez,
                                   columns = "SYMBOL",
                                   keytype = "ENTREZID")
  # Deduplicate (some Entrez map to multiple symbols)
  mapping <- mapping %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  setNames(mapping$SYMBOL, mapping$ENTREZID)
}, error = function(e) {
  stop("Failed to map Entrez IDs to symbols: ", e$message)
})

cat(sprintf("  Mapped %d / %d Entrez IDs to symbols (%.1f%%)\n",
            length(entrez_to_symbol), length(all_entrez),
            100 * length(entrez_to_symbol) / length(all_entrez)))

# =============================================================================
# RE-COMPUTE COORDINATED CHANGES (from mc_dmr + hmc_dmr, for independence)
# =============================================================================

cat("\n--- Computing coordinated methylation changes ---\n")

# All genes with mC DMR data (the background universe)
# Deduplicate by gene symbol (keep first occurrence) to avoid many-to-many joins
mc_all <- mc_dmr %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_sig = significant, mc_direction = direction) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

hmc_all <- hmc_dmr %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_sig = significant, hmc_direction = direction) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Merge mC and hmC for all genes (inner join: genes tested in both)
gene_methylation <- dplyr::inner_join(mc_all, hmc_all, by = "gene")

# Define coordinated pattern: mC up AND hmC down (significant in at least mC)
gene_methylation <- gene_methylation %>%
  dplyr::mutate(
    coordinated = (mc_diff > 0 & hmc_diff < 0),
    coordinated_sig = (mc_sig & mc_diff > 0 & hmc_diff < 0),
    is_hyper = (mc_sig & mc_diff > 0)
  )

# Add delta_ratio from ratio_df
gene_methylation <- gene_methylation %>%
  dplyr::left_join(
    ratio_df %>% dplyr::select(gene, delta_ratio, chromatin_state),
    by = "gene"
  )

cat(sprintf("  Total genes in methylation universe: %d\n", nrow(gene_methylation)))
cat(sprintf("  Coordinated (mC up/hmC down, any significance): %d (%.1f%%)\n",
            sum(gene_methylation$coordinated),
            100 * mean(gene_methylation$coordinated)))
cat(sprintf("  Coordinated + mC significant: %d\n", sum(gene_methylation$coordinated_sig)))
cat(sprintf("  Hypermethylated (mC sig, mc_diff > 0): %d\n", sum(gene_methylation$is_hyper)))

# =============================================================================
# BUILD GREAT-STYLE GENE REGULATORY DOMAINS
# =============================================================================

cat("\n--- Building GREAT-style regulatory domains ---\n")

gene_domains_info <- build_gene_regulatory_domains()
cat(sprintf("  Built regulatory domains for %d genes\n", nrow(gene_domains_info)))
cat(sprintf("  Parameters: %dkb upstream, %dkb downstream, %dkb max extension\n",
            GREAT_UPSTREAM / 1000, GREAT_DOWNSTREAM / 1000, GREAT_MAX_EXTENSION / 1000))

# Convert to GRanges
gene_domains_gr <- GRanges(
  seqnames = gene_domains_info$chr,
  ranges = IRanges(start = gene_domains_info$reg_start, end = gene_domains_info$reg_end),
  entrez_id = gene_domains_info$entrez_id
)

# Map Entrez to Symbol for the domain genes
domain_entrez_ids <- unique(gene_domains_info$entrez_id)
domain_mapping <- tryCatch({
  m <- AnnotationDbi::select(org.Mm.eg.db,
                             keys = domain_entrez_ids,
                             columns = "SYMBOL",
                             keytype = "ENTREZID")
  m %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::distinct(ENTREZID, .keep_all = TRUE)
}, error = function(e) {
  stop("Failed domain Entrez->Symbol mapping: ", e$message)
})
domain_entrez_to_symbol <- setNames(domain_mapping$SYMBOL, domain_mapping$ENTREZID)

# =============================================================================
# ASSOCIATE ANCHORS TO GENES
# =============================================================================

cat("\n--- Associating loop anchors to genes ---\n")

# Build anchor GRanges from loops (both anchors)
anchor1_gr <- GRanges(
  seqnames = loops$anchor1_chr,
  ranges = IRanges(start = loops$anchor1_start, end = loops$anchor1_end),
  anchor_id = paste0(loops$loop_id, "_A1"),
  loop_direction = loops$direction,
  loop_distance = loops$loop_distance,
  anchor_type = loops$anchor1_type,
  loop_logFC = loops$logFC
)

anchor2_gr <- GRanges(
  seqnames = loops$anchor2_chr,
  ranges = IRanges(start = loops$anchor2_start, end = loops$anchor2_end),
  anchor_id = paste0(loops$loop_id, "_A2"),
  loop_direction = loops$direction,
  loop_distance = loops$loop_distance,
  anchor_type = loops$anchor2_type,
  loop_logFC = loops$logFC
)

all_anchors_gr <- c(anchor1_gr, anchor2_gr)
cat(sprintf("  Total anchors: %d (%d loops x 2)\n", length(all_anchors_gr), nrow(loops)))

# GREAT-style association
anchor_gene_assoc <- associate_anchors_to_genes(all_anchors_gr, gene_domains_gr)

# Map entrez to symbol
anchor_gene_assoc$symbol <- domain_entrez_to_symbol[anchor_gene_assoc$entrez_id]
anchor_gene_assoc <- anchor_gene_assoc %>% dplyr::filter(!is.na(symbol))

cat(sprintf("  Anchor-gene associations (GREAT): %d\n", nrow(anchor_gene_assoc)))
cat(sprintf("  Unique genes at anchors: %d\n", length(unique(anchor_gene_assoc$symbol))))

# Nearest-gene (secondary) approach
nearest_genes_entrez <- unique(c(as.character(loops$anchor1_nearest_gene),
                                 as.character(loops$anchor2_nearest_gene)))
nearest_genes_entrez <- nearest_genes_entrez[!is.na(nearest_genes_entrez) & nearest_genes_entrez != ""]
nearest_genes_symbol <- entrez_to_symbol[nearest_genes_entrez]
nearest_genes_symbol <- unique(nearest_genes_symbol[!is.na(nearest_genes_symbol)])
cat(sprintf("  Unique genes (nearest-gene): %d\n", length(nearest_genes_symbol)))

# Deduplicated gene sets for Fisher's tests
great_anchor_genes <- unique(anchor_gene_assoc$symbol)
# Restrict to genes in methylation universe
great_anchor_genes_in_universe <- intersect(great_anchor_genes, gene_methylation$gene)
nearest_anchor_genes_in_universe <- intersect(nearest_genes_symbol, gene_methylation$gene)
cat(sprintf("  GREAT anchor genes in DMR universe: %d\n", length(great_anchor_genes_in_universe)))
cat(sprintf("  Nearest anchor genes in DMR universe: %d\n", length(nearest_anchor_genes_in_universe)))

# Save anchor-gene mapping table
anchor_gene_export <- anchor_gene_assoc %>%
  dplyr::select(anchor_id, anchor_chr, anchor_start, anchor_end,
                entrez_id, symbol, loop_direction, loop_distance, anchor_type)
write.table(anchor_gene_export,
            file.path(TABLES_DIR, "anchor_gene_associations_great.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: anchor_gene_associations_great.tsv\n")

# =============================================================================
# 27a: COORDINATED GENE ENRICHMENT AT LOOP ANCHORS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("27a: COORDINATED GENE ENRICHMENT AT LOOP ANCHORS\n")
cat("================================================================================\n\n")

run_coordinated_enrichment <- function(anchor_genes, method_label, gene_meth) {
  # Mark which genes are at anchors
  gene_meth$at_anchor <- gene_meth$gene %in% anchor_genes

  # 2x2 Fisher's: coordinated x anchor
  a <- sum(gene_meth$at_anchor & gene_meth$coordinated)
  b <- sum(!gene_meth$at_anchor & gene_meth$coordinated)
  c <- sum(gene_meth$at_anchor & !gene_meth$coordinated)
  d <- sum(!gene_meth$at_anchor & !gene_meth$coordinated)

  fisher_mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                       dimnames = list(c("Coordinated", "Not_coordinated"),
                                       c("Anchor", "Non-anchor")))
  fisher_res <- fisher.test(fisher_mat)

  cat(sprintf("  %s: Fisher's OR = %.2f, %s\n",
              method_label, fisher_res$estimate, fmt_p(fisher_res$p.value)))
  cat(sprintf("    Contingency: a=%d, b=%d, c=%d, d=%d\n", a, b, c, d))
  cat(sprintf("    %% coordinated at anchors: %.1f%%, background: %.1f%%\n",
              100 * a / (a + c), 100 * b / (b + d)))

  list(
    method = method_label,
    n_anchor = sum(gene_meth$at_anchor),
    n_coordinated_anchor = a,
    pct_coordinated_anchor = 100 * a / (a + c),
    pct_coordinated_background = 100 * b / (b + d),
    odds_ratio = fisher_res$estimate,
    p_value = fisher_res$p.value,
    ci_lower = fisher_res$conf.int[1],
    ci_upper = fisher_res$conf.int[2]
  )
}

# Run for both methods
res_great <- run_coordinated_enrichment(great_anchor_genes_in_universe, "GREAT-style", gene_methylation)
res_nearest <- run_coordinated_enrichment(nearest_anchor_genes_in_universe, "Nearest-gene", gene_methylation)

# Export table
enrichment_table <- dplyr::bind_rows(
  as.data.frame(res_great, stringsAsFactors = FALSE),
  as.data.frame(res_nearest, stringsAsFactors = FALSE)
)
write.table(enrichment_table,
            file.path(TABLES_DIR, "methylation_loop_anchor_coordinated_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: methylation_loop_anchor_coordinated_enrichment.tsv\n")

# Figure 27a: Grouped bar chart
plot_data_27a <- data.frame(
  Method = rep(c("GREAT-style", "Nearest-gene"), each = 2),
  Group = rep(c("Loop anchor genes", "Background"), 2),
  Pct_coordinated = c(res_great$pct_coordinated_anchor, res_great$pct_coordinated_background,
                      res_nearest$pct_coordinated_anchor, res_nearest$pct_coordinated_background),
  stringsAsFactors = FALSE
)
plot_data_27a$Group <- factor(plot_data_27a$Group, levels = c("Loop anchor genes", "Background"))

p_27a <- ggplot(plot_data_27a, aes(x = Method, y = Pct_coordinated, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Loop anchor genes" = "#4575b4", "Background" = "grey70")) +
  labs(
    title = "Coordinated mC-up / hmC-down Pattern at Loop Anchors",
    subtitle = sprintf("GREAT: OR=%.2f, %s | Nearest: OR=%.2f, %s",
                       res_great$odds_ratio, fmt_p(res_great$p_value),
                       res_nearest$odds_ratio, fmt_p(res_nearest$p_value)),
    x = "Association method", y = "% genes with coordinated pattern",
    fill = NULL
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_27a, file.path(SECTION_DIR, "27a_coordinated_enrichment_at_loop_anchors"),
                        width = 8, height = 6)

# =============================================================================
# 27b: STRATIFICATION BY LOOP DIRECTION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("27b: STRATIFICATION BY LOOP DIRECTION\n")
cat("================================================================================\n\n")

# Build direction-specific gene sets (GREAT-style)
lost_anchor_assoc <- anchor_gene_assoc %>%
  dplyr::filter(loop_direction == "down_in_mutant")
gained_anchor_assoc <- anchor_gene_assoc %>%
  dplyr::filter(loop_direction == "up_in_mutant")

lost_genes <- unique(lost_anchor_assoc$symbol)
gained_genes <- unique(gained_anchor_assoc$symbol)

lost_genes_in_universe <- intersect(lost_genes, gene_methylation$gene)
gained_genes_in_universe <- intersect(gained_genes, gene_methylation$gene)

cat(sprintf("  Lost-loop anchor genes in universe: %d\n", length(lost_genes_in_universe)))
cat(sprintf("  Gained-loop anchor genes in universe: %d\n", length(gained_genes_in_universe)))

# 27b test 1: Fisher's — hypermethylated enrichment at lost vs gained
gene_methylation$at_lost <- gene_methylation$gene %in% lost_genes_in_universe
gene_methylation$at_gained <- gene_methylation$gene %in% gained_genes_in_universe
gene_methylation$at_neither <- !gene_methylation$at_lost & !gene_methylation$at_gained

# Hypermethylated rates
hyper_lost <- sum(gene_methylation$at_lost & gene_methylation$is_hyper)
not_hyper_lost <- sum(gene_methylation$at_lost & !gene_methylation$is_hyper)
hyper_gained <- sum(gene_methylation$at_gained & gene_methylation$is_hyper)
not_hyper_gained <- sum(gene_methylation$at_gained & !gene_methylation$is_hyper)
hyper_bg <- sum(gene_methylation$at_neither & gene_methylation$is_hyper)
not_hyper_bg <- sum(gene_methylation$at_neither & !gene_methylation$is_hyper)

pct_hyper_lost <- 100 * hyper_lost / (hyper_lost + not_hyper_lost)
pct_hyper_gained <- 100 * hyper_gained / (hyper_gained + not_hyper_gained)
pct_hyper_bg <- 100 * hyper_bg / (hyper_bg + not_hyper_bg)

# Fisher: lost vs background
fisher_lost_bg <- fisher.test(matrix(c(hyper_lost, hyper_bg, not_hyper_lost, not_hyper_bg), nrow = 2))
# Fisher: gained vs background
fisher_gained_bg <- fisher.test(matrix(c(hyper_gained, hyper_bg, not_hyper_gained, not_hyper_bg), nrow = 2))
# Fisher: lost vs gained
fisher_lost_gained <- fisher.test(matrix(c(hyper_lost, hyper_gained, not_hyper_lost, not_hyper_gained), nrow = 2))

cat(sprintf("  %% hypermethylated: Lost=%.1f%%, Gained=%.1f%%, Background=%.1f%%\n",
            pct_hyper_lost, pct_hyper_gained, pct_hyper_bg))
cat(sprintf("  Lost vs BG: OR=%.2f, %s\n", fisher_lost_bg$estimate, fmt_p(fisher_lost_bg$p.value)))
cat(sprintf("  Gained vs BG: OR=%.2f, %s\n", fisher_gained_bg$estimate, fmt_p(fisher_gained_bg$p.value)))
cat(sprintf("  Lost vs Gained: OR=%.2f, %s\n", fisher_lost_gained$estimate, fmt_p(fisher_lost_gained$p.value)))

# 27b test 2: Wilcoxon on mc_diff by loop direction
mc_diff_lost <- gene_methylation$mc_diff[gene_methylation$at_lost]
mc_diff_gained <- gene_methylation$mc_diff[gene_methylation$at_gained]
mc_diff_bg <- gene_methylation$mc_diff[gene_methylation$at_neither]
wilcox_mc_lost_gained <- wilcox.test(mc_diff_lost, mc_diff_gained)
wilcox_mc_lost_bg <- wilcox.test(mc_diff_lost, mc_diff_bg)

cat(sprintf("  mc_diff Wilcoxon (lost vs gained): %s\n", fmt_p(wilcox_mc_lost_gained$p.value)))
cat(sprintf("  mc_diff Wilcoxon (lost vs bg): %s\n", fmt_p(wilcox_mc_lost_bg$p.value)))

# 27b test 3: Wilcoxon on delta_ratio by loop direction
dr_lost <- gene_methylation$delta_ratio[gene_methylation$at_lost & !is.na(gene_methylation$delta_ratio)]
dr_gained <- gene_methylation$delta_ratio[gene_methylation$at_gained & !is.na(gene_methylation$delta_ratio)]
dr_bg <- gene_methylation$delta_ratio[gene_methylation$at_neither & !is.na(gene_methylation$delta_ratio)]
wilcox_dr_lost_gained <- wilcox.test(dr_lost, dr_gained)
wilcox_dr_lost_bg <- wilcox.test(dr_lost, dr_bg)

cat(sprintf("  delta_ratio median: Lost=%.4f, Gained=%.4f, BG=%.4f\n",
            median(dr_lost), median(dr_gained), median(dr_bg)))
cat(sprintf("  delta_ratio Wilcoxon (lost vs gained): %s\n", fmt_p(wilcox_dr_lost_gained$p.value)))

# Export direction table
direction_table <- data.frame(
  comparison = c("Lost_vs_BG", "Gained_vs_BG", "Lost_vs_Gained"),
  pct_hyper_group1 = c(pct_hyper_lost, pct_hyper_gained, pct_hyper_lost),
  pct_hyper_group2 = c(pct_hyper_bg, pct_hyper_bg, pct_hyper_gained),
  fisher_or = c(fisher_lost_bg$estimate, fisher_gained_bg$estimate, fisher_lost_gained$estimate),
  fisher_p = c(fisher_lost_bg$p.value, fisher_gained_bg$p.value, fisher_lost_gained$p.value),
  wilcox_mc_diff_p = c(wilcox_mc_lost_bg$p.value, NA, wilcox_mc_lost_gained$p.value),
  wilcox_delta_ratio_p = c(wilcox_dr_lost_bg$p.value, NA, wilcox_dr_lost_gained$p.value),
  median_delta_ratio_group1 = c(median(dr_lost), median(dr_gained), median(dr_lost)),
  median_delta_ratio_group2 = c(median(dr_bg), median(dr_bg), median(dr_gained)),
  n_group1 = c(length(lost_genes_in_universe), length(gained_genes_in_universe), length(lost_genes_in_universe)),
  n_group2 = c(sum(gene_methylation$at_neither), sum(gene_methylation$at_neither), length(gained_genes_in_universe)),
  stringsAsFactors = FALSE
)
write.table(direction_table,
            file.path(TABLES_DIR, "methylation_direction_by_loop_direction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: methylation_direction_by_loop_direction.tsv\n")

# Figure 27b1: Grouped bar (% hyper at lost vs gained vs background)
plot_data_27b1 <- data.frame(
  Group = factor(c("Lost loops", "Gained loops", "Background"),
                 levels = c("Lost loops", "Gained loops", "Background")),
  Pct_hyper = c(pct_hyper_lost, pct_hyper_gained, pct_hyper_bg),
  N = c(length(lost_genes_in_universe), length(gained_genes_in_universe),
        sum(gene_methylation$at_neither)),
  stringsAsFactors = FALSE
)

p_27b1 <- ggplot(plot_data_27b1, aes(x = Group, y = Pct_hyper, fill = Group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("n=%d", N)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Lost loops" = "#d73027",
                               "Gained loops" = "#4575b4",
                               "Background" = "grey70")) +
  labs(
    title = "Hypermethylation Rate by Loop Direction",
    subtitle = sprintf("Lost vs Gained: OR=%.2f, %s",
                       fisher_lost_gained$estimate, fmt_p(fisher_lost_gained$p.value)),
    x = NULL, y = "% genes with significant mC increase"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_27b1, file.path(SECTION_DIR, "27b_mc_direction_by_loop_direction"),
                        width = 7, height = 6)

# Figure 27b2: Violin + boxplot of mc_diff by loop direction
violin_data_mc <- data.frame(
  mc_diff = c(mc_diff_lost, mc_diff_gained, mc_diff_bg),
  Group = factor(
    c(rep("Lost loops", length(mc_diff_lost)),
      rep("Gained loops", length(mc_diff_gained)),
      rep("Background", length(mc_diff_bg))),
    levels = c("Lost loops", "Gained loops", "Background")
  ),
  stringsAsFactors = FALSE
)

p_27b2 <- ggplot(violin_data_mc, aes(x = Group, y = mc_diff, fill = Group)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Lost loops" = "#d73027",
                               "Gained loops" = "#4575b4",
                               "Background" = "grey70")) +
  labs(
    title = "mC Difference by Loop Direction",
    subtitle = sprintf("Lost vs Gained: %s | Lost vs BG: %s",
                       fmt_p(wilcox_mc_lost_gained$p.value), fmt_p(wilcox_mc_lost_bg$p.value)),
    x = NULL, y = "mC difference (mutant - control)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_27b2, file.path(SECTION_DIR, "27b_mc_diff_violin_by_loop_direction"),
                        width = 7, height = 6)

# Figure 27b3: Violin + boxplot of delta_ratio by loop direction
violin_data_dr <- data.frame(
  delta_ratio = c(dr_lost, dr_gained, dr_bg),
  Group = factor(
    c(rep("Lost loops", length(dr_lost)),
      rep("Gained loops", length(dr_gained)),
      rep("Background", length(dr_bg))),
    levels = c("Lost loops", "Gained loops", "Background")
  ),
  stringsAsFactors = FALSE
)

p_27b3 <- ggplot(violin_data_dr, aes(x = Group, y = delta_ratio, fill = Group)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Lost loops" = "#d73027",
                               "Gained loops" = "#4575b4",
                               "Background" = "grey70")) +
  labs(
    title = "Delta Demethylation Ratio by Loop Direction",
    subtitle = sprintf("Lost vs Gained: %s | Median: Lost=%.4f, Gained=%.4f, BG=%.4f",
                       fmt_p(wilcox_dr_lost_gained$p.value),
                       median(dr_lost), median(dr_gained), median(dr_bg)),
    x = NULL, y = "Delta ratio (more negative = greater TET impairment)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_27b3, file.path(SECTION_DIR, "27b_delta_ratio_violin_by_loop_direction"),
                        width = 7, height = 6)

# =============================================================================
# 27c: K119ub-LOOP-METHYLATION CONVERGENCE
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("27c: K119ub-LOOP-METHYLATION CONVERGENCE\n")
cat("================================================================================\n\n")

# Annotate each anchor for K119ub gained overlap
anchor_k119ub_gained <- countOverlaps(all_anchors_gr, k119ub_diff$gained) > 0
anchor_k119ub_lost <- countOverlaps(all_anchors_gr, k119ub_diff$lost) > 0

# Build anchor-level data with K119ub status
anchor_gene_assoc$k119ub_gained_at_anchor <- anchor_k119ub_gained[anchor_gene_assoc$anchor_idx]
anchor_gene_assoc$k119ub_lost_at_anchor <- anchor_k119ub_lost[anchor_gene_assoc$anchor_idx]

# Aggregate to gene level: gene has K119ub gained if ANY of its anchors do
gene_k119ub <- anchor_gene_assoc %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    k119ub_gained = any(k119ub_gained_at_anchor),
    k119ub_lost = any(k119ub_lost_at_anchor),
    .groups = "drop"
  )

# Merge with gene methylation
gene_convergence <- gene_methylation %>%
  dplyr::inner_join(gene_k119ub, by = c("gene" = "symbol"))

cat(sprintf("  Genes with anchor + K119ub data: %d\n", nrow(gene_convergence)))

# Fisher's 2x2: K119ub gained x hypermethylated
a_conv <- sum(gene_convergence$k119ub_gained & gene_convergence$is_hyper)
b_conv <- sum(!gene_convergence$k119ub_gained & gene_convergence$is_hyper)
c_conv <- sum(gene_convergence$k119ub_gained & !gene_convergence$is_hyper)
d_conv <- sum(!gene_convergence$k119ub_gained & !gene_convergence$is_hyper)

fisher_convergence <- fisher.test(matrix(c(a_conv, b_conv, c_conv, d_conv), nrow = 2,
                                         dimnames = list(c("Hyper", "Not_hyper"),
                                                         c("K119ub_gained", "No_K119ub_gained"))))

cat(sprintf("  K119ub gained x Hypermethylated: OR=%.2f, %s\n",
            fisher_convergence$estimate, fmt_p(fisher_convergence$p.value)))
cat(sprintf("    a=%d, b=%d, c=%d, d=%d\n", a_conv, b_conv, c_conv, d_conv))

# Compute O/E for 2x2 heatmap (K119ub direction x mC direction)
k119ub_states <- c("K119ub gained", "K119ub not gained")
mc_states <- c("Hypermethylated", "Not hypermethylated")

oe_matrix <- matrix(NA, nrow = 2, ncol = 2,
                    dimnames = list(k119ub_states, mc_states))
count_matrix <- matrix(c(a_conv, c_conv, b_conv, d_conv), nrow = 2,
                       dimnames = list(k119ub_states, mc_states))

total_n <- sum(count_matrix)
for (i in 1:2) {
  for (j in 1:2) {
    expected <- sum(count_matrix[i, ]) * sum(count_matrix[, j]) / total_n
    oe_matrix[i, j] <- count_matrix[i, j] / expected
  }
}

cat("  O/E matrix:\n")
print(round(oe_matrix, 2))

# Figure 27c1: O/E heatmap
oe_df <- expand.grid(
  K119ub = k119ub_states,
  Methylation = mc_states,
  stringsAsFactors = FALSE
)
oe_df$OE <- as.vector(oe_matrix)
oe_df$Count <- as.vector(count_matrix)
oe_df$Label <- sprintf("O/E=%.2f\n(n=%d)", oe_df$OE, oe_df$Count)

p_27c1 <- ggplot(oe_df, aes(x = Methylation, y = K119ub, fill = OE)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Label), size = 4) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 1, name = "O/E ratio") +
  labs(
    title = "K119ub x Methylation Convergence at Loop Anchors",
    subtitle = sprintf("Fisher's OR = %.2f, %s",
                       fisher_convergence$estimate, fmt_p(fisher_convergence$p.value)),
    x = "Methylation status", y = "K119ub status at anchor"
  ) +
  theme_biomodal() +
  theme(panel.grid = element_blank())

save_multiformat_ggplot(p_27c1, file.path(SECTION_DIR, "27c_k119ub_methylation_loop_convergence_heatmap"),
                        width = 8, height = 6)

# Triple convergence: K119ub gained + hypermethylated + loop lost
gene_convergence$at_lost_loop <- gene_convergence$gene %in% lost_genes_in_universe
triple_count <- sum(gene_convergence$k119ub_gained & gene_convergence$is_hyper & gene_convergence$at_lost_loop)
double_k119ub_hyper <- sum(gene_convergence$k119ub_gained & gene_convergence$is_hyper)
double_k119ub_lost <- sum(gene_convergence$k119ub_gained & gene_convergence$at_lost_loop)
double_hyper_lost <- sum(gene_convergence$is_hyper & gene_convergence$at_lost_loop)

cat(sprintf("\n  Triple convergence (K119ub gained + hyper + lost loop): %d genes\n", triple_count))
cat(sprintf("  K119ub gained + hyper: %d\n", double_k119ub_hyper))
cat(sprintf("  K119ub gained + lost loop: %d\n", double_k119ub_lost))
cat(sprintf("  Hyper + lost loop: %d\n", double_hyper_lost))

# Figure 27c2: Triple convergence bar
convergence_bar_data <- data.frame(
  Category = factor(
    c("K119ub gained +\nhypermethylated +\nlost loop",
      "K119ub gained +\nhypermethylated",
      "K119ub gained +\nlost loop",
      "Hypermethylated +\nlost loop"),
    levels = c("Hypermethylated +\nlost loop",
               "K119ub gained +\nlost loop",
               "K119ub gained +\nhypermethylated",
               "K119ub gained +\nhypermethylated +\nlost loop")
  ),
  Count = c(triple_count, double_k119ub_hyper, double_k119ub_lost, double_hyper_lost),
  stringsAsFactors = FALSE
)

p_27c2 <- ggplot(convergence_bar_data, aes(x = Category, y = Count)) +
  geom_col(fill = "#756BB1", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  labs(
    title = "Multi-Layer Convergence at Loop Anchors",
    subtitle = "Genes with overlapping epigenomic and 3D chromatin changes",
    x = NULL, y = "Number of genes"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 9))

save_multiformat_ggplot(p_27c2, file.path(SECTION_DIR, "27c_triple_convergence_summary"),
                        width = 9, height = 6)

# Export convergence table
convergence_export <- data.frame(
  test = "K119ub_gained_x_Hypermethylated_at_loop_anchors",
  a = a_conv, b = b_conv, c = c_conv, d = d_conv,
  odds_ratio = fisher_convergence$estimate,
  p_value = fisher_convergence$p.value,
  ci_lower = fisher_convergence$conf.int[1],
  ci_upper = fisher_convergence$conf.int[2],
  triple_convergence_count = triple_count,
  double_k119ub_hyper = double_k119ub_hyper,
  double_k119ub_lost = double_k119ub_lost,
  double_hyper_lost = double_hyper_lost,
  stringsAsFactors = FALSE
)
write.table(convergence_export,
            file.path(TABLES_DIR, "k119ub_loop_methylation_convergence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: k119ub_loop_methylation_convergence.tsv\n")

# =============================================================================
# 27d: LOGISTIC REGRESSION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("27d: LOGISTIC REGRESSION\n")
cat("================================================================================\n\n")

# Build modeling data: anchor-gene pairs from GREAT-style mapping
# Need: is_hyper, loop_direction, loop_distance, anchor_type, distance_to_tss
model_df <- anchor_gene_assoc %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::inner_join(gene_methylation %>% dplyr::select(gene, is_hyper, mc_diff, delta_ratio),
                    by = c("symbol" = "gene")) %>%
  dplyr::mutate(
    loop_dir_binary = ifelse(loop_direction == "down_in_mutant", "Lost", "Gained"),
    log_loop_distance = log10(pmax(loop_distance, 1))
  )

# Compute distance to TSS for each anchor-gene pair
# Use gene_domains_info which has TSS positions
model_df <- model_df %>%
  dplyr::left_join(
    gene_domains_info %>%
      dplyr::mutate(symbol = domain_entrez_to_symbol[entrez_id]) %>%
      dplyr::filter(!is.na(symbol)) %>%
      dplyr::select(symbol, tss),
    by = c("symbol")
  ) %>%
  dplyr::mutate(
    anchor_midpoint = (anchor_start + anchor_end) / 2,
    distance_to_tss = abs(anchor_midpoint - tss),
    log_dist_tss = log10(distance_to_tss + 1)
  )

# Clean anchor_type for modeling (set NA to "Other")
model_df$anchor_type[is.na(model_df$anchor_type) | model_df$anchor_type == ""] <- "Other"
model_df$anchor_type <- factor(model_df$anchor_type, levels = CHROMATIN_STATE_ORDER)

# Remove rows with missing outcome
model_df <- model_df %>% dplyr::filter(!is.na(is_hyper))

cat(sprintf("  Modeling data: %d anchor-gene pairs\n", nrow(model_df)))
cat(sprintf("  Hypermethylated: %d (%.1f%%)\n",
            sum(model_df$is_hyper), 100 * mean(model_df$is_hyper)))

# Logistic model: hypermethylation ~ loop direction + distance features
logistic_fit <- tryCatch({
  glm(is_hyper ~ loop_dir_binary + log_loop_distance + anchor_type + log_dist_tss,
      family = binomial, data = model_df)
}, error = function(e) {
  cat("  WARNING: Full logistic model failed:", e$message, "\n")
  cat("  Falling back to reduced model (direction + distance only)\n")
  glm(is_hyper ~ loop_dir_binary + log_loop_distance + log_dist_tss,
      family = binomial, data = model_df)
})

logistic_summary <- summary(logistic_fit)
cat("\n  Logistic regression coefficients:\n")
print(round(logistic_summary$coefficients, 4))

# Extract ORs with 95% CI
or_df <- data.frame(
  term = names(coef(logistic_fit)),
  estimate = coef(logistic_fit),
  or = exp(coef(logistic_fit)),
  or_lower = exp(confint.default(logistic_fit)[, 1]),
  or_upper = exp(confint.default(logistic_fit)[, 2]),
  p_value = logistic_summary$coefficients[, 4],
  stringsAsFactors = FALSE,
  row.names = NULL
)

write.table(or_df,
            file.path(TABLES_DIR, "logistic_regression_methylation_loop.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: logistic_regression_methylation_loop.tsv\n")

# Figure 27d1: Forest plot of ORs
or_plot_df <- or_df %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    term_clean = gsub("loop_dir_binary", "Loop: ", term),
    term_clean = gsub("log_loop_distance", "log10(loop distance)", term_clean),
    term_clean = gsub("log_dist_tss", "log10(dist to TSS)", term_clean),
    term_clean = gsub("anchor_type", "Anchor: ", term_clean),
    sig_label = ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(p_value < 0.05, "*", "ns")))
  )

p_27d1 <- ggplot(or_plot_df, aes(x = or, y = reorder(term_clean, or))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = or_lower, xmax = or_upper), width = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "#4575b4") +
  geom_text(aes(label = sprintf("%.2f %s", or, sig_label)),
            hjust = -0.15, size = 3) +
  scale_x_log10() +
  labs(
    title = "Logistic Regression: Predictors of Hypermethylation",
    subtitle = sprintf("Model: is_hyper ~ loop_direction + log(loop_dist) + anchor_type + log(dist_TSS)\nAIC = %.1f",
                       AIC(logistic_fit)),
    x = "Odds Ratio (log scale)", y = NULL
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 9))

save_multiformat_ggplot(p_27d1, file.path(SECTION_DIR, "27d_logistic_regression_forest_plot"),
                        width = 10, height = 7)

# Linear model: delta_ratio ~ same predictors
model_dr <- model_df %>% dplyr::filter(!is.na(delta_ratio))

linear_fit <- tryCatch({
  lm(delta_ratio ~ loop_dir_binary + log_loop_distance + anchor_type + log_dist_tss,
     data = model_dr)
}, error = function(e) {
  cat("  WARNING: Full linear model failed:", e$message, "\n")
  lm(delta_ratio ~ loop_dir_binary + log_loop_distance + log_dist_tss,
     data = model_dr)
})

linear_summary <- summary(linear_fit)
cat("\n  Linear model (delta_ratio) summary:\n")
cat(sprintf("  R-squared: %.4f, Adj R-squared: %.4f\n",
            linear_summary$r.squared, linear_summary$adj.r.squared))
f_stat <- linear_summary$fstatistic
f_p <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
cat(sprintf("  F-statistic: %.2f, %s\n", f_stat[1], fmt_p(f_p)))

linear_coef_df <- data.frame(
  term = rownames(linear_summary$coefficients),
  estimate = linear_summary$coefficients[, 1],
  std_error = linear_summary$coefficients[, 2],
  t_value = linear_summary$coefficients[, 3],
  p_value = linear_summary$coefficients[, 4],
  stringsAsFactors = FALSE,
  row.names = NULL
)

write.table(linear_coef_df,
            file.path(TABLES_DIR, "linear_model_delta_ratio_loop.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: linear_model_delta_ratio_loop.tsv\n")

# Figure 27d2: Linear model coefficient plot
lm_plot_df <- linear_coef_df %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::mutate(
    term_clean = gsub("loop_dir_binary", "Loop: ", term),
    term_clean = gsub("log_loop_distance", "log10(loop distance)", term_clean),
    term_clean = gsub("log_dist_tss", "log10(dist to TSS)", term_clean),
    term_clean = gsub("anchor_type", "Anchor: ", term_clean),
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error,
    sig_label = ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(p_value < 0.05, "*", "ns")))
  )

p_27d2 <- ggplot(lm_plot_df, aes(x = estimate, y = reorder(term_clean, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "#d73027") +
  geom_text(aes(label = sprintf("%.4f %s", estimate, sig_label)),
            hjust = -0.15, size = 3) +
  labs(
    title = "Linear Model: Predictors of Delta Demethylation Ratio",
    subtitle = sprintf("R\u00b2 = %.4f, F = %.2f, %s",
                       linear_summary$r.squared, f_stat[1], fmt_p(f_p)),
    x = "Coefficient (effect on delta ratio)", y = NULL
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 9))

save_multiformat_ggplot(p_27d2, file.path(SECTION_DIR, "27d_linear_model_coefficients"),
                        width = 10, height = 7)

# =============================================================================
# 27e: SHARED ANCHOR METHYLATION PROFILE
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("27e: SHARED ANCHOR METHYLATION PROFILE\n")
cat("================================================================================\n\n")

# GREAT-style mapping for shared anchors
shared_anchor_gr <- GRanges(
  seqnames = shared_anchors$chr,
  ranges = IRanges(start = shared_anchors$start, end = shared_anchors$end),
  anchor_id = shared_anchors$anchor_id
)

shared_anchor_gene_assoc <- associate_anchors_to_genes(shared_anchor_gr, gene_domains_gr)
shared_anchor_gene_assoc$symbol <- domain_entrez_to_symbol[shared_anchor_gene_assoc$entrez_id]
shared_anchor_gene_assoc <- shared_anchor_gene_assoc %>% dplyr::filter(!is.na(symbol))
shared_great_genes <- unique(shared_anchor_gene_assoc$symbol)
shared_great_in_universe <- intersect(shared_great_genes, gene_methylation$gene)

cat(sprintf("  Shared anchors: %d\n", nrow(shared_anchors)))
cat(sprintf("  GREAT-style genes at shared anchors: %d\n", length(shared_great_genes)))
cat(sprintf("  In methylation universe: %d\n", length(shared_great_in_universe)))

# Also use pre-computed shared anchor genes (42 genes) for validation
pre_shared_genes <- unique(shared_anchor_genes$gene_symbol)
pre_shared_in_universe <- intersect(pre_shared_genes, gene_methylation$gene)
cat(sprintf("  Pre-computed shared anchor genes in universe: %d\n", length(pre_shared_in_universe)))

# Three groups: shared-anchor, non-shared-anchor, background
gene_methylation$at_shared <- gene_methylation$gene %in% shared_great_in_universe
gene_methylation$at_nonshared <- gene_methylation$gene %in% great_anchor_genes_in_universe &
                                  !gene_methylation$at_shared

# Coordinated rates
coord_shared <- sum(gene_methylation$at_shared & gene_methylation$coordinated)
n_shared <- sum(gene_methylation$at_shared)
pct_shared <- 100 * coord_shared / n_shared

coord_nonshared <- sum(gene_methylation$at_nonshared & gene_methylation$coordinated)
n_nonshared <- sum(gene_methylation$at_nonshared)
pct_nonshared <- 100 * coord_nonshared / n_nonshared

coord_bg_27e <- sum(!gene_methylation$at_shared & !gene_methylation$at_nonshared & gene_methylation$coordinated)
n_bg_27e <- sum(!gene_methylation$at_shared & !gene_methylation$at_nonshared)
pct_bg_27e <- 100 * coord_bg_27e / n_bg_27e

cat(sprintf("  %% coordinated: Shared=%.1f%%, Non-shared anchor=%.1f%%, Background=%.1f%%\n",
            pct_shared, pct_nonshared, pct_bg_27e))

# Fisher: shared vs non-shared
fisher_shared <- fisher.test(
  matrix(c(coord_shared, coord_nonshared,
           n_shared - coord_shared, n_nonshared - coord_nonshared), nrow = 2,
         dimnames = list(c("Coordinated", "Not"), c("Shared", "Non-shared")))
)
# Fisher: shared vs background
fisher_shared_bg <- fisher.test(
  matrix(c(coord_shared, coord_bg_27e,
           n_shared - coord_shared, n_bg_27e - coord_bg_27e), nrow = 2,
         dimnames = list(c("Coordinated", "Not"), c("Shared", "Background")))
)

cat(sprintf("  Shared vs Non-shared: OR=%.2f, %s\n",
            fisher_shared$estimate, fmt_p(fisher_shared$p.value)))
cat(sprintf("  Shared vs Background: OR=%.2f, %s\n",
            fisher_shared_bg$estimate, fmt_p(fisher_shared_bg$p.value)))

# Wilcoxon on delta_ratio across 3 groups
dr_shared <- gene_methylation$delta_ratio[gene_methylation$at_shared & !is.na(gene_methylation$delta_ratio)]
dr_nonshared <- gene_methylation$delta_ratio[gene_methylation$at_nonshared & !is.na(gene_methylation$delta_ratio)]
dr_bg_27e <- gene_methylation$delta_ratio[!gene_methylation$at_shared & !gene_methylation$at_nonshared &
                                             !is.na(gene_methylation$delta_ratio)]

wilcox_shared_nonshared <- wilcox.test(dr_shared, dr_nonshared)
wilcox_shared_bg_27e <- wilcox.test(dr_shared, dr_bg_27e)

cat(sprintf("  Delta-ratio median: Shared=%.4f, Non-shared=%.4f, BG=%.4f\n",
            median(dr_shared), median(dr_nonshared), median(dr_bg_27e)))
cat(sprintf("  Shared vs Non-shared Wilcoxon: %s\n", fmt_p(wilcox_shared_nonshared$p.value)))
cat(sprintf("  Shared vs Background Wilcoxon: %s\n", fmt_p(wilcox_shared_bg_27e$p.value)))

# Export table
shared_table <- data.frame(
  group = c("Shared_anchor", "Non-shared_anchor", "Background"),
  n_genes = c(n_shared, n_nonshared, n_bg_27e),
  n_coordinated = c(coord_shared, coord_nonshared, coord_bg_27e),
  pct_coordinated = c(pct_shared, pct_nonshared, pct_bg_27e),
  median_delta_ratio = c(median(dr_shared), median(dr_nonshared), median(dr_bg_27e)),
  fisher_vs_bg_or = c(fisher_shared_bg$estimate, NA, NA),
  fisher_vs_bg_p = c(fisher_shared_bg$p.value, NA, NA),
  fisher_shared_vs_nonshared_or = c(fisher_shared$estimate, NA, NA),
  fisher_shared_vs_nonshared_p = c(fisher_shared$p.value, NA, NA),
  wilcox_vs_bg_p = c(wilcox_shared_bg_27e$p.value, NA, NA),
  stringsAsFactors = FALSE
)
write.table(shared_table,
            file.path(TABLES_DIR, "shared_anchor_methylation_profile.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: shared_anchor_methylation_profile.tsv\n")

# Figure 27e1: Three-group bar chart
plot_data_27e1 <- data.frame(
  Group = factor(c("Shared anchors", "Non-shared anchors", "Background"),
                 levels = c("Shared anchors", "Non-shared anchors", "Background")),
  Pct_coordinated = c(pct_shared, pct_nonshared, pct_bg_27e),
  N = c(n_shared, n_nonshared, n_bg_27e),
  stringsAsFactors = FALSE
)

p_27e1 <- ggplot(plot_data_27e1, aes(x = Group, y = Pct_coordinated, fill = Group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("n=%d", N)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Shared anchors" = "#8c510a",
                               "Non-shared anchors" = "#4575b4",
                               "Background" = "grey70")) +
  labs(
    title = "Coordinated mC-up / hmC-down Rate at Shared vs Non-shared Anchors",
    subtitle = sprintf("Shared vs BG: OR=%.2f, %s | Shared vs Non-shared: OR=%.2f, %s",
                       fisher_shared_bg$estimate, fmt_p(fisher_shared_bg$p.value),
                       fisher_shared$estimate, fmt_p(fisher_shared$p.value)),
    x = NULL, y = "% genes with coordinated pattern"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_27e1, file.path(SECTION_DIR, "27e_shared_anchor_coordinated_rate"),
                        width = 8, height = 6)

# Figure 27e2: Violin + boxplot of delta_ratio across 3 groups
violin_data_27e <- data.frame(
  delta_ratio = c(dr_shared, dr_nonshared, dr_bg_27e),
  Group = factor(
    c(rep("Shared anchors", length(dr_shared)),
      rep("Non-shared anchors", length(dr_nonshared)),
      rep("Background", length(dr_bg_27e))),
    levels = c("Shared anchors", "Non-shared anchors", "Background")
  ),
  stringsAsFactors = FALSE
)

p_27e2 <- ggplot(violin_data_27e, aes(x = Group, y = delta_ratio, fill = Group)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Shared anchors" = "#8c510a",
                               "Non-shared anchors" = "#4575b4",
                               "Background" = "grey70")) +
  labs(
    title = "Delta Demethylation Ratio: Shared vs Non-shared Anchors",
    subtitle = sprintf("Shared vs Non-shared: %s | Shared vs BG: %s",
                       fmt_p(wilcox_shared_nonshared$p.value),
                       fmt_p(wilcox_shared_bg_27e$p.value)),
    x = NULL, y = "Delta ratio (more negative = greater TET impairment)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_27e2, file.path(SECTION_DIR, "27e_shared_anchor_delta_ratio_violin"),
                        width = 8, height = 6)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 27 SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Input: %d differential loops, %d genes in methylation universe\n",
            nrow(loops), nrow(gene_methylation)))
cat(sprintf("GREAT-style anchor genes: %d (in universe: %d)\n",
            length(great_anchor_genes), length(great_anchor_genes_in_universe)))
cat(sprintf("Nearest-gene anchor genes: %d (in universe: %d)\n",
            length(nearest_genes_symbol), length(nearest_anchor_genes_in_universe)))
cat("\n")

cat("27a - Coordinated enrichment:\n")
cat(sprintf("  GREAT: OR=%.2f, %s\n", res_great$odds_ratio, fmt_p(res_great$p_value)))
cat(sprintf("  Nearest: OR=%.2f, %s\n", res_nearest$odds_ratio, fmt_p(res_nearest$p_value)))

cat("\n27b - Direction stratification:\n")
cat(sprintf("  Lost vs Gained hyper rate: OR=%.2f, %s\n",
            fisher_lost_gained$estimate, fmt_p(fisher_lost_gained$p.value)))
cat(sprintf("  Delta-ratio: Lost median=%.4f, Gained median=%.4f\n",
            median(dr_lost), median(dr_gained)))

cat("\n27c - K119ub convergence:\n")
cat(sprintf("  K119ub gained x hyper: OR=%.2f, %s\n",
            fisher_convergence$estimate, fmt_p(fisher_convergence$p.value)))
cat(sprintf("  Triple convergence: %d genes\n", triple_count))

cat("\n27d - Logistic regression:\n")
cat(sprintf("  AIC: %.1f\n", AIC(logistic_fit)))
cat(sprintf("  Linear model R\u00b2: %.4f\n", linear_summary$r.squared))

cat("\n27e - Shared anchors:\n")
cat(sprintf("  Shared vs BG coordinated: OR=%.2f, %s\n",
            fisher_shared_bg$estimate, fmt_p(fisher_shared_bg$p.value)))

cat("\n--- Output files ---\n")
cat("Figures (10 panels in ", SECTION_DIR, "):\n")
cat("  27a_coordinated_enrichment_at_loop_anchors/\n")
cat("  27b_mc_direction_by_loop_direction/\n")
cat("  27b_mc_diff_violin_by_loop_direction/\n")
cat("  27b_delta_ratio_violin_by_loop_direction/\n")
cat("  27c_k119ub_methylation_loop_convergence_heatmap/\n")
cat("  27c_triple_convergence_summary/\n")
cat("  27d_logistic_regression_forest_plot/\n")
cat("  27d_linear_model_coefficients/\n")
cat("  27e_shared_anchor_coordinated_rate/\n")
cat("  27e_shared_anchor_delta_ratio_violin/\n")
cat("\nTables (in ", TABLES_DIR, "):\n")
cat("  anchor_gene_associations_great.tsv\n")
cat("  methylation_loop_anchor_coordinated_enrichment.tsv\n")
cat("  methylation_direction_by_loop_direction.tsv\n")
cat("  k119ub_loop_methylation_convergence.tsv\n")
cat("  logistic_regression_methylation_loop.tsv\n")
cat("  linear_model_delta_ratio_loop.tsv\n")
cat("  shared_anchor_methylation_profile.tsv\n")

cat("\n================================================================================\n")
cat("SECTION 27 COMPLETE\n")
cat("================================================================================\n\n")
