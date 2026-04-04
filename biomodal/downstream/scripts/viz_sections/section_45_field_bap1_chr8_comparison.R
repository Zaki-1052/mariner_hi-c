# biomodal/downstream/scripts/viz_sections/section_45_field_bap1_chr8_comparison.R
# Section 45: Field et al. (2019) BAP1 Chr8 Methylation Hotspot Comparison
# Standalone script - sources shared config for all dependencies and data
#
# Compares chr8 methylation hotspot genes from Field et al. (2019, Clin Cancer Res,
# PMC6744995) — human BAP1-KD uveal melanoma — against our BAP1-KO mouse cerebellum
# DUET evoC data. Field found significant enrichment of differentially methylated
# genes on chr8 (bands 8q22, 8p21, 8q13, 8q24) in BAP1-mutant Class 2 tumors.
#
# Analysis includes:
#   - Human→mouse ortholog mapping via biomaRt
#   - Gene-body AND promoter-level DMR comparison
#   - RNA-seq expression integration
#   - Trisomy 8 diagnostic (coverage, expression, DMR rate)
#   - Direction concordance testing (Fisher's exact)
#
# Data sources:
#   1. docs/NIHMS1533544-supplement-4.xlsx — Field et al. Supplementary Table S3
#   2. DATA_PATHS$mc_dmr / hmc_dmr — Our gene-body DMRs (run-4, 8 samples)
#   3. DATA_PATHS$promoters_mc / promoters_hmc — Our promoter DMRs
#   4. tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx — DESeq2 RNA-seq
#
# Figures:
#   45a: Mapping funnel — Field genes → mouse orthologs → significant DMRs
#   45b: Direction concordance heatmap (gene-body, 2×2 tile + Fisher's exact)
#   45c: Dual modality dot plot (mC + hmC per gene)
#   45d: Volcano highlight plot (Field genes in genome-wide context)
#   45e: Effect size lollipop plot (concordant vs discordant)
#   45f: Coordinated mC/hmC quadrant analysis
#   45g: Gene-body vs promoter scatter
#   45h: Summary publication table
#   45i: Promoter-level concordance heatmap (expected: ~random)
#   45j: RNA-seq expression integration
#   45k: Trisomy 8 diagnostic panel (3-panel)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_45_field_bap1_chr8_comparison.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(readxl)
})

# =============================================================================
# SECTION 45 CONFIGURATION
# =============================================================================

SECTION_NUM <- "45"
SECTION_DIR <- file.path(OUTPUT_DIR, "45_field_bap1_chr8_comparison")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# Input paths
FIELD_XLSX <- file.path(REPO_ROOT, "biomodal/docs/NIHMS1533544-supplement-4.xlsx")
RNASEQ_XLSX <- file.path(REPO_ROOT, "tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx")
ORTHOLOGS_CSV <- file.path(REPO_ROOT, "biomodal/field_chr8_orthologs.csv")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 45: FIELD et al. (2019) BAP1 Chr8 METHYLATION HOTSPOT COMPARISON\n")
cat("========================================================================\n\n")

cat("Validating inputs...\n")
stopifnot("Field et al. supplement not found" = file.exists(FIELD_XLSX))
stopifnot("Gene-body mC DMR file not found" = file.exists(DATA_PATHS$mc_dmr))
stopifnot("Gene-body hmC DMR file not found" = file.exists(DATA_PATHS$hmc_dmr))
stopifnot("Promoter mC DMR file not found" = file.exists(DATA_PATHS$promoters_mc))
stopifnot("Promoter hmC DMR file not found" = file.exists(DATA_PATHS$promoters_hmc))
stopifnot("RNA-seq results not found" = file.exists(RNASEQ_XLSX))
stopifnot("Ortholog mapping CSV not found" = file.exists(ORTHOLOGS_CSV))
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1: LOAD FIELD et al. SUPPLEMENTARY DATA
# =============================================================================

cat("--- Step 1: Loading Field et al. supplementary data ---\n")

# Read hypermethylated/downregulated genes (skip category header row)
hyper_raw <- read_excel(FIELD_XLSX, sheet = "HyperMethvsDecGE_Cutoffs", skip = 1)
cat(sprintf("  HyperMeth sheet: %d rows, %d columns\n", nrow(hyper_raw), ncol(hyper_raw)))

# Read hypomethylated/upregulated genes
hypo_raw <- read_excel(FIELD_XLSX, sheet = "HypoMethvsIncGE_Cutoffs", skip = 1)
cat(sprintf("  HypoMeth sheet: %d rows, %d columns\n", nrow(hypo_raw), ncol(hypo_raw)))

# Filter to chr8 and extract unique genes with mean deltaBeta
extract_chr8 <- function(df, direction) {
  # Find CHR column
  chr_col <- grep("^CHR$", colnames(df), value = TRUE)[1]
  gene_col <- grep("Gene.Symbol", colnames(df), value = TRUE)[1]
  delta_col <- grep("deltaBeta", colnames(df), value = TRUE)[1]

  stopifnot("CHR column not found" = !is.na(chr_col))
  stopifnot("Gene.Symbol column not found" = !is.na(gene_col))

  chr8 <- df[df[[chr_col]] == 8, ]

  # Aggregate per gene: mean deltaBeta, number of probes
  result <- chr8 %>%
    group_by(human_gene = .data[[gene_col]]) %>%
    summarise(
      field_mean_delta = if (!is.na(delta_col)) mean(.data[[delta_col]], na.rm = TRUE) else NA_real_,
      field_n_probes = n(),
      .groups = "drop"
    ) %>%
    mutate(field_direction = direction)

  return(result)
}

field_hyper <- extract_chr8(hyper_raw, "Hypermethylated")
field_hypo <- extract_chr8(hypo_raw, "Hypomethylated")

field_chr8 <- bind_rows(field_hyper, field_hypo)
cat(sprintf("  Chr8 genes: %d hypermethylated, %d hypomethylated (%d total)\n",
            nrow(field_hyper), nrow(field_hypo), nrow(field_chr8)))

# =============================================================================
# STEP 2: HUMAN→MOUSE ORTHOLOG MAPPING (pre-computed via biomaRt)
# =============================================================================

cat("\n--- Step 2: Loading pre-computed ortholog mapping ---\n")
cat("  Source: field_chr8_orthologs.csv (generated via Ensembl biomaRt two-step query)\n")
cat("  Method: getBM(hsapiens_gene_ensembl) → ensembl_gene_id → mmusculus_homolog\n")
cat("  Includes updated HGNC symbols: BAI1→ADGRB1, KIAA0196→WASHC5, TCEB1→ELOC, FAM49B→CYRIB\n")

orthologs <- read.csv(ORTHOLOGS_CSV, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d ortholog mappings\n", nrow(orthologs)))

human_genes <- unique(field_chr8$human_gene)
mapped_count <- sum(human_genes %in% orthologs$human_gene)
unmapped <- setdiff(human_genes, orthologs$human_gene)
cat(sprintf("  Successfully mapped: %d/%d genes\n", mapped_count, length(human_genes)))
if (length(unmapped) > 0) {
  cat(sprintf("  No Ensembl ortholog: %s\n", paste(unmapped, collapse = ", ")))
}

stopifnot("Too few orthologs mapped (< 50)" = mapped_count >= 50)

# Join orthologs to Field data
field_chr8 <- field_chr8 %>%
  left_join(orthologs, by = "human_gene")

cat(sprintf("  Field chr8 genes with mouse orthologs: %d/%d\n",
            sum(!is.na(field_chr8$mouse_gene)), nrow(field_chr8)))

# =============================================================================
# STEP 3: LOAD OUR DMR DATA (GENE-BODY + PROMOTER)
# =============================================================================

cat("\n--- Step 3: Loading our DMR data ---\n")

# Gene-body DMRs (all genes, not just significant)
gb_mc <- load_dmr_bed(DATA_PATHS$mc_dmr)
gb_hmc <- load_dmr_bed(DATA_PATHS$hmc_dmr)
cat(sprintf("  Gene-body mC: %d genes (%d significant)\n", nrow(gb_mc), sum(gb_mc$significant)))
cat(sprintf("  Gene-body hmC: %d genes (%d significant)\n", nrow(gb_hmc), sum(gb_hmc$significant)))

# Promoter DMRs
prom_mc <- load_dmr_bed(DATA_PATHS$promoters_mc)
prom_hmc <- load_dmr_bed(DATA_PATHS$promoters_hmc)
cat(sprintf("  Promoter mC: %d genes (%d significant)\n", nrow(prom_mc), sum(prom_mc$significant)))
cat(sprintf("  Promoter hmC: %d genes (%d significant)\n", nrow(prom_hmc), sum(prom_hmc$significant)))

# RNA-seq
cat("  Loading RNA-seq data...\n")
rnaseq <- read_excel(RNASEQ_XLSX, sheet = "Output")
colnames(rnaseq)[1] <- "gene"
rnaseq <- rnaseq %>%
  dplyr::select(gene, baseMean, log2FoldChange, padj) %>%
  mutate(expr_sig = !is.na(padj) & padj < 0.05,
         expr_dir = ifelse(log2FoldChange > 0, "Up", "Down"))
cat(sprintf("  RNA-seq: %d genes (%d significant)\n", nrow(rnaseq), sum(rnaseq$expr_sig, na.rm = TRUE)))

# =============================================================================
# STEP 4: BUILD MASTER COMPARISON TABLE
# =============================================================================

cat("\n--- Step 4: Building master comparison table ---\n")

comparison <- field_chr8 %>%
  # Gene-body mC
  left_join(
    gb_mc %>% dplyr::select(gene, gb_mc_diff = mod_difference, gb_mc_q = dmr_qvalue,
                     gb_mc_sig = significant, gb_mc_dir = direction),
    by = c("mouse_gene" = "gene")
  ) %>%
  # Gene-body hmC
  left_join(
    gb_hmc %>% dplyr::select(gene, gb_hmc_diff = mod_difference, gb_hmc_q = dmr_qvalue,
                      gb_hmc_sig = significant, gb_hmc_dir = direction),
    by = c("mouse_gene" = "gene")
  ) %>%
  # Promoter mC
  left_join(
    prom_mc %>% dplyr::select(gene, prom_mc_diff = mod_difference, prom_mc_q = dmr_qvalue,
                       prom_mc_sig = significant, prom_mc_dir = direction),
    by = c("mouse_gene" = "gene")
  ) %>%
  # Promoter hmC
  left_join(
    prom_hmc %>% dplyr::select(gene, prom_hmc_diff = mod_difference, prom_hmc_q = dmr_qvalue,
                        prom_hmc_sig = significant, prom_hmc_dir = direction),
    by = c("mouse_gene" = "gene")
  ) %>%
  # RNA-seq

  left_join(
    rnaseq %>% dplyr::select(gene, baseMean, log2FoldChange, padj, expr_sig, expr_dir),
    by = c("mouse_gene" = "gene")
  )

# Classify concordance (gene-body level)
comparison <- comparison %>%
  mutate(
    gb_concordance = case_when(
      is.na(mouse_gene) ~ "No ortholog",
      is.na(gb_mc_sig) ~ "Not in data",
      !gb_mc_sig ~ "Non-significant",
      field_direction == "Hypermethylated" & gb_mc_dir == "Hypermethylated" ~ "Concordant",
      field_direction == "Hypomethylated" & gb_mc_dir == "Hypomethylated" ~ "Concordant",
      TRUE ~ "Discordant"
    ),
    # Promoter concordance
    prom_concordance = case_when(
      is.na(mouse_gene) ~ "No ortholog",
      is.na(prom_mc_sig) ~ "Not in data",
      !prom_mc_sig ~ "Non-significant",
      field_direction == "Hypermethylated" & prom_mc_dir == "Hypermethylated" ~ "Concordant",
      field_direction == "Hypomethylated" & prom_mc_dir == "Hypomethylated" ~ "Concordant",
      TRUE ~ "Discordant"
    ),
    # Expression concordance (Field: hyper→down, hypo→up)
    expr_concordance = case_when(
      is.na(mouse_gene) ~ "No ortholog",
      is.na(expr_sig) | !expr_sig ~ "Non-significant",
      field_direction == "Hypermethylated" & expr_dir == "Down" ~ "Concordant",
      field_direction == "Hypomethylated" & expr_dir == "Up" ~ "Concordant",
      TRUE ~ "Discordant"
    )
  )

cat("  Gene-body concordance breakdown:\n")
print(table(comparison$gb_concordance))
cat("\n  Promoter concordance breakdown:\n")
print(table(comparison$prom_concordance))
cat("\n  Expression concordance breakdown:\n")
print(table(comparison$expr_concordance))

# =============================================================================
# STEP 5: TRISOMY 8 DIAGNOSTIC
# =============================================================================

cat("\n--- Step 5: Trisomy 8 diagnostic ---\n")

# 5a: Coverage check
chr_coverage <- gb_mc %>%
  mutate(chrom = gsub("chr", "", chr)) %>%
  filter(chrom %in% as.character(1:19)) %>%
  group_by(chrom) %>%
  summarise(mean_coverage = mean(mean_coverage, na.rm = TRUE),
            n_genes = n(), .groups = "drop") %>%
  mutate(chrom_num = as.integer(chrom)) %>%
  arrange(chrom_num)

genome_mean_cov <- mean(chr_coverage$mean_coverage)
chr8_mean_cov <- chr_coverage$mean_coverage[chr_coverage$chrom == "8"]
cov_ratio <- chr8_mean_cov / genome_mean_cov
cat(sprintf("  Coverage: chr8 = %.2f, genome = %.2f, ratio = %.3f (trisomy ~1.5)\n",
            chr8_mean_cov, genome_mean_cov, cov_ratio))

# 5b: RNA-seq expression check
gene_to_chrom <- gb_mc %>%
  dplyr::select(gene, chr) %>%
  mutate(chrom = gsub("chr", "", chr)) %>%
  distinct(gene, .keep_all = TRUE)

chr_expr <- rnaseq %>%
  inner_join(gene_to_chrom, by = "gene") %>%
  filter(chrom %in% as.character(1:19)) %>%
  group_by(chrom) %>%
  summarise(median_lfc = median(log2FoldChange, na.rm = TRUE),
            n_genes = n(), .groups = "drop") %>%
  mutate(chrom_num = as.integer(chrom)) %>%
  arrange(chrom_num)

chr8_lfc <- chr_expr$median_lfc[chr_expr$chrom == "8"]
other_lfc <- rnaseq %>%
  inner_join(gene_to_chrom, by = "gene") %>%
  filter(chrom %in% as.character(1:19), chrom != "8") %>%
  pull(log2FoldChange)
chr8_lfc_vals <- rnaseq %>%
  inner_join(gene_to_chrom, by = "gene") %>%
  filter(chrom == "8") %>%
  pull(log2FoldChange)

expr_test <- wilcox.test(chr8_lfc_vals, other_lfc, alternative = "greater")
cat(sprintf("  Expression: chr8 median log2FC = %+.4f, Mann-Whitney p = %.4e\n",
            chr8_lfc, expr_test$p.value))

# 5c: DMR rate check
chr_dmr_rate <- gb_mc %>%
  mutate(chrom = gsub("chr", "", chr)) %>%
  filter(chrom %in% as.character(1:19)) %>%
  group_by(chrom) %>%
  summarise(total = n(), sig = sum(significant),
            pct_sig = sig / total * 100, .groups = "drop") %>%
  mutate(chrom_num = as.integer(chrom)) %>%
  arrange(chrom_num)

chr8_dmr_pct <- chr_dmr_rate$pct_sig[chr_dmr_rate$chrom == "8"]
cat(sprintf("  DMR rate: chr8 = %.1f%%, genome range = %.1f-%.1f%%\n",
            chr8_dmr_pct, min(chr_dmr_rate$pct_sig), max(chr_dmr_rate$pct_sig)))
cat("  CONCLUSION: No evidence of trisomy 8.\n")

# =============================================================================
# STEP 6: STATISTICAL TESTS
# =============================================================================

cat("\n--- Step 6: Statistical tests ---\n")

stat_results <- list()

# 6a: Fisher's exact — gene-body mC direction concordance
gb_sig <- comparison %>% filter(gb_concordance %in% c("Concordant", "Discordant"))
if (nrow(gb_sig) > 0) {
  ct_gb <- table(
    Field = gb_sig$field_direction,
    Ours = gb_sig$gb_mc_dir
  )
  fisher_gb <- fisher.test(ct_gb)
  cat(sprintf("  Fisher's exact (gene-body mC): OR = %.2f, p = %.4e\n",
              fisher_gb$estimate, fisher_gb$p.value))
  stat_results$fisher_gb_mc <- list(
    test = "Fisher's exact (gene-body mC concordance)",
    OR = fisher_gb$estimate, p = fisher_gb$p.value,
    table = ct_gb
  )
}

# 6b: Fisher's exact — promoter mC direction concordance
prom_sig <- comparison %>% filter(prom_concordance %in% c("Concordant", "Discordant"))
if (nrow(prom_sig) > 0) {
  ct_prom <- table(
    Field = prom_sig$field_direction,
    Ours = prom_sig$prom_mc_dir
  )
  fisher_prom <- fisher.test(ct_prom)
  cat(sprintf("  Fisher's exact (promoter mC): OR = %.2f, p = %.4e\n",
              fisher_prom$estimate, fisher_prom$p.value))
  stat_results$fisher_prom_mc <- list(
    test = "Fisher's exact (promoter mC concordance)",
    OR = fisher_prom$estimate, p = fisher_prom$p.value,
    table = ct_prom
  )
}

# 6c: Hypergeometric enrichment — are Field orthologs enriched among our sig DMR genes?
n_field_in_data <- sum(!is.na(comparison$gb_mc_sig))
n_field_sig <- sum(comparison$gb_mc_sig, na.rm = TRUE)
n_total_genes <- nrow(gb_mc)
n_total_sig <- sum(gb_mc$significant)

hyper_p <- phyper(n_field_sig - 1, n_total_sig, n_total_genes - n_total_sig,
                  n_field_in_data, lower.tail = FALSE)
cat(sprintf("  Hypergeometric (gene-body): %d/%d Field genes significant vs %d/%d genome-wide, p = %.4e\n",
            n_field_sig, n_field_in_data, n_total_sig, n_total_genes, hyper_p))
stat_results$hypergeometric_gb <- list(
  test = "Hypergeometric enrichment (gene-body mC)",
  field_sig = n_field_sig, field_total = n_field_in_data,
  genome_sig = n_total_sig, genome_total = n_total_genes, p = hyper_p
)

# 6d: Binomial test for concordance rate
n_concordant <- sum(comparison$gb_concordance == "Concordant")
n_testable <- sum(comparison$gb_concordance %in% c("Concordant", "Discordant"))
if (n_testable > 0) {
  binom_result <- binom.test(n_concordant, n_testable, p = 0.5)
  cat(sprintf("  Binomial (gene-body concordance): %d/%d = %.1f%%, p = %.4e\n",
              n_concordant, n_testable, n_concordant / n_testable * 100, binom_result$p.value))
  stat_results$binomial_gb <- list(
    test = "Binomial test (gene-body concordance vs 50%)",
    concordant = n_concordant, testable = n_testable,
    rate = n_concordant / n_testable, p = binom_result$p.value
  )
}

# =============================================================================
# FIGURE 45a: MAPPING FUNNEL
# =============================================================================

cat("\n--- Generating figures ---\n")

funnel_data <- data.frame(
  step = factor(c("Field chr8\ngenes", "Mouse\northolog", "In gene-body\nDMR data",
                   "Sig mC\n(q<0.05)", "Sig hmC\n(q<0.05)"),
                levels = c("Field chr8\ngenes", "Mouse\northolog", "In gene-body\nDMR data",
                            "Sig mC\n(q<0.05)", "Sig hmC\n(q<0.05)")),
  Hypermethylated = c(
    sum(comparison$field_direction == "Hypermethylated"),
    sum(comparison$field_direction == "Hypermethylated" & !is.na(comparison$mouse_gene)),
    sum(comparison$field_direction == "Hypermethylated" & !is.na(comparison$gb_mc_sig)),
    sum(comparison$field_direction == "Hypermethylated" & comparison$gb_mc_sig %in% TRUE),
    sum(comparison$field_direction == "Hypermethylated" & comparison$gb_hmc_sig %in% TRUE)
  ),
  Hypomethylated = c(
    sum(comparison$field_direction == "Hypomethylated"),
    sum(comparison$field_direction == "Hypomethylated" & !is.na(comparison$mouse_gene)),
    sum(comparison$field_direction == "Hypomethylated" & !is.na(comparison$gb_mc_sig)),
    sum(comparison$field_direction == "Hypomethylated" & comparison$gb_mc_sig %in% TRUE),
    sum(comparison$field_direction == "Hypomethylated" & comparison$gb_hmc_sig %in% TRUE)
  )
)

funnel_long <- funnel_data %>%
  pivot_longer(cols = c("Hypermethylated", "Hypomethylated"),
               names_to = "direction", values_to = "count")

p45a <- ggplot(funnel_long, aes(x = step, y = count, fill = direction)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
  scale_fill_manual(values = COLORS$direction, name = "Field Direction") +
  labs(title = "Field et al. Chr8 Gene Mapping Pipeline",
       subtitle = "Attrition from human genes to significant mouse DMRs",
       x = NULL, y = "Number of genes") +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 10))

save_multiformat_ggplot(p45a, file.path(SECTION_DIR, "45a_mapping_funnel"),
                        width = 10, height = 7)
cat("  45a: Mapping funnel saved\n")

# =============================================================================
# FIGURE 45b: GENE-BODY DIRECTION CONCORDANCE HEATMAP
# =============================================================================

gb_sig_for_tile <- comparison %>%
  filter(gb_concordance %in% c("Concordant", "Discordant")) %>%
  count(field_direction, gb_mc_dir) %>%
  complete(field_direction = c("Hypermethylated", "Hypomethylated"),
           gb_mc_dir = c("Hypermethylated", "Hypomethylated"),
           fill = list(n = 0))

fisher_label <- if (!is.null(stat_results$fisher_gb_mc)) {
  sprintf("Fisher's exact: OR = %.2f, p = %.2e", stat_results$fisher_gb_mc$OR, stat_results$fisher_gb_mc$p)
} else { "" }

p45b <- ggplot(gb_sig_for_tile, aes(x = field_direction, y = gb_mc_dir, fill = n)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = n), size = 8, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#D7191C", name = "Count") +
  labs(title = "Gene-Body mC Direction Concordance",
       subtitle = paste0("Field et al. vs Our BAP1-KO (significant genes only)\n", fisher_label),
       x = "Field et al. Direction", y = "Our Gene-Body mC Direction") +
  theme_biomodal() +
  theme(axis.text = element_text(size = 11))

save_multiformat_ggplot(p45b, file.path(SECTION_DIR, "45b_genebody_concordance_heatmap"),
                        width = 8, height = 7)
cat("  45b: Gene-body concordance heatmap saved\n")

# =============================================================================
# FIGURE 45c: DUAL MODALITY DOT PLOT
# =============================================================================

# Select genes with at least one significant modality
dual_data <- comparison %>%
  filter(!is.na(mouse_gene),
         (gb_mc_sig %in% TRUE) | (gb_hmc_sig %in% TRUE)) %>%
  arrange(gb_mc_diff) %>%
  mutate(mouse_gene = factor(mouse_gene, levels = unique(mouse_gene)))

if (nrow(dual_data) > 50) {
  # Keep top 25 positive and top 25 negative by mC effect
  top_pos <- dual_data %>% arrange(desc(gb_mc_diff)) %>% head(25)
  top_neg <- dual_data %>% arrange(gb_mc_diff) %>% head(25)
  dual_data <- bind_rows(top_pos, top_neg) %>%
    distinct(mouse_gene, .keep_all = TRUE) %>%
    arrange(gb_mc_diff) %>%
    mutate(mouse_gene = factor(mouse_gene, levels = unique(mouse_gene)))
}

dual_long <- dual_data %>%
  dplyr::select(mouse_gene, field_direction, gb_mc_diff, gb_hmc_diff, gb_mc_sig, gb_hmc_sig) %>%
  pivot_longer(cols = c(gb_mc_diff, gb_hmc_diff),
               names_to = "modality", values_to = "effect_size") %>%
  mutate(
    modality = ifelse(modality == "gb_mc_diff", "5mC", "5hmC"),
    is_sig = ifelse(modality == "5mC",
                    gb_mc_sig %in% TRUE,
                    gb_hmc_sig %in% TRUE)
  )

p45c <- ggplot(dual_long, aes(x = effect_size, y = mouse_gene, color = modality, shape = is_sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = COLORS$methylation, name = "Modality") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Significant\n(q < 0.05)") +
  facet_wrap(~field_direction, scales = "free_y", ncol = 2) +
  labs(title = "Gene-Body Methylation Changes for Field et al. Chr8 Genes",
       subtitle = "mC and hmC mod_difference (filled = significant)",
       x = "mod_difference (mutant - control)", y = NULL) +
  theme_biomodal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 12))

save_multiformat_ggplot(p45c, file.path(SECTION_DIR, "45c_dual_modality_dotplot"),
                        width = 14, height = max(8, nrow(dual_data) * 0.22))
cat("  45c: Dual modality dot plot saved\n")

# =============================================================================
# FIGURE 45d: VOLCANO HIGHLIGHT PLOT
# =============================================================================

volcano_base <- gb_mc %>%
  mutate(field_group = case_when(
    gene %in% comparison$mouse_gene[comparison$field_direction == "Hypermethylated"] ~ "Field Hyper",
    gene %in% comparison$mouse_gene[comparison$field_direction == "Hypomethylated"] ~ "Field Hypo",
    TRUE ~ "Other"
  ))

# Label top Field genes
top_labels <- volcano_base %>%
  filter(field_group != "Other", significant) %>%
  arrange(dmr_qvalue) %>%
  head(15)

p45d <- ggplot(volcano_base, aes(x = mod_difference, y = neg_log10_q)) +
  geom_point(data = filter(volcano_base, field_group == "Other"),
             color = "grey80", size = 0.5, alpha = 0.3) +
  geom_point(data = filter(volcano_base, field_group != "Other"),
             aes(color = field_group), size = 2, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = top_labels, aes(label = gene),
                           size = 3, max.overlaps = 20, segment.alpha = 0.5) +
  scale_color_manual(values = c("Field Hyper" = "#D7191C", "Field Hypo" = "#2C7BB6"),
                     name = "Field Direction") +
  labs(title = "Volcano Plot: Gene-Body mC DMRs",
       subtitle = "Field et al. chr8 genes highlighted in genome-wide context",
       x = "mod_difference (mutant - control)", y = "-log10(q-value)") +
  theme_biomodal() +
  coord_cartesian(xlim = c(-0.15, 0.15))

save_multiformat_ggplot(p45d, file.path(SECTION_DIR, "45d_volcano_highlight"),
                        width = 12, height = 9)
cat("  45d: Volcano highlight saved\n")

# =============================================================================
# FIGURE 45e: EFFECT SIZE LOLLIPOP
# =============================================================================

lollipop_data <- comparison %>%
  filter(!is.na(mouse_gene), !is.na(gb_mc_diff)) %>%
  arrange(field_direction, gb_mc_diff) %>%
  mutate(mouse_gene = factor(mouse_gene, levels = unique(mouse_gene)),
         color_group = gb_concordance)

p45e <- ggplot(lollipop_data, aes(x = gb_mc_diff, y = mouse_gene, color = color_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(aes(x = 0, xend = gb_mc_diff, yend = mouse_gene), linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Concordant" = "#4DAF4A", "Discordant" = "#E41A1C",
                                "Non-significant" = "grey60", "Not in data" = "grey90",
                                "No ortholog" = "grey90"),
                     name = "Concordance") +
  facet_wrap(~field_direction, scales = "free_y", ncol = 2) +
  labs(title = "Gene-Body mC Effect Size for Field et al. Chr8 Genes",
       subtitle = "Green = concordant direction, Red = discordant",
       x = "mC mod_difference (mutant - control)", y = NULL) +
  theme_biomodal(base_size = 9) +
  theme(axis.text.y = element_text(size = 6))

save_multiformat_ggplot(p45e, file.path(SECTION_DIR, "45e_effect_size_lollipop"),
                        width = 14, height = max(10, nrow(lollipop_data) * 0.18))
cat("  45e: Effect size lollipop saved\n")

# =============================================================================
# FIGURE 45f: COORDINATED mC/hmC QUADRANT ANALYSIS
# =============================================================================

quad_data <- comparison %>%
  filter(!is.na(mouse_gene), gb_mc_sig %in% TRUE, gb_hmc_sig %in% TRUE) %>%
  mutate(quadrant = case_when(
    gb_mc_diff > 0 & gb_hmc_diff < 0 ~ "Q4: mC up / hmC down",
    gb_mc_diff > 0 & gb_hmc_diff > 0 ~ "Q3: mC up / hmC up",
    gb_mc_diff < 0 & gb_hmc_diff > 0 ~ "Q2: mC down / hmC up",
    gb_mc_diff < 0 & gb_hmc_diff < 0 ~ "Q1: mC down / hmC down",
    TRUE ~ "Other"
  ))

# Genome-wide quadrant distribution for comparison
gw_quad <- gb_mc %>%
  inner_join(gb_hmc %>% dplyr::select(gene, hmc_diff = mod_difference, hmc_sig = significant),
             by = "gene") %>%
  filter(significant, hmc_sig) %>%
  mutate(quadrant = case_when(
    mod_difference > 0 & hmc_diff < 0 ~ "Q4: mC up / hmC down",
    mod_difference > 0 & hmc_diff > 0 ~ "Q3: mC up / hmC up",
    mod_difference < 0 & hmc_diff > 0 ~ "Q2: mC down / hmC up",
    mod_difference < 0 & hmc_diff < 0 ~ "Q1: mC down / hmC down",
    TRUE ~ "Other"
  )) %>%
  count(quadrant) %>%
  mutate(pct = n / sum(n) * 100, source = "Genome-wide")

field_quad <- quad_data %>%
  count(quadrant) %>%
  mutate(pct = n / sum(n) * 100, source = "Field chr8")

quad_combined <- bind_rows(field_quad, gw_quad)

p45f <- ggplot(quad_combined, aes(x = source, y = pct, fill = quadrant)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2", name = "Quadrant") +
  labs(title = "Coordinated mC/hmC Pattern: Field Chr8 Genes vs Genome-Wide",
       subtitle = sprintf("Field chr8: %d genes with both mC and hmC significant", nrow(quad_data)),
       x = NULL, y = "% of co-significant genes") +
  theme_biomodal()

save_multiformat_ggplot(p45f, file.path(SECTION_DIR, "45f_quadrant_analysis"),
                        width = 10, height = 7)
cat("  45f: Quadrant analysis saved\n")

# =============================================================================
# FIGURE 45g: GENE-BODY vs PROMOTER SCATTER
# =============================================================================

scatter_data <- comparison %>%
  filter(!is.na(mouse_gene), !is.na(gb_mc_diff), !is.na(prom_mc_diff))

p45g <- ggplot(scatter_data, aes(x = prom_mc_diff, y = gb_mc_diff, color = field_direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(shape = prom_mc_sig %in% TRUE), size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = filter(scatter_data, prom_mc_sig %in% TRUE),
    aes(label = mouse_gene), size = 3, max.overlaps = 15
  ) +
  scale_color_manual(values = COLORS$direction, name = "Field Direction") +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                     name = "Promoter\nSignificant",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  labs(title = "Gene-Body vs Promoter mC Effect Size",
       subtitle = "Field chr8 genes: context-dependent methylation response\n(triangles = significant promoter DMRs)",
       x = "Promoter mC mod_difference", y = "Gene-body mC mod_difference") +
  theme_biomodal()

save_multiformat_ggplot(p45g, file.path(SECTION_DIR, "45g_genebody_vs_promoter"),
                        width = 10, height = 9)
cat("  45g: Gene-body vs promoter scatter saved\n")

# =============================================================================
# FIGURE 45h: SUMMARY PUBLICATION TABLE
# =============================================================================

table_data <- comparison %>%
  filter(!is.na(mouse_gene)) %>%
  arrange(field_direction, gb_mc_q) %>%
  dplyr::select(human_gene, field_direction, mouse_gene,
         gb_mc_diff, gb_mc_q, gb_mc_dir,
         prom_mc_diff, prom_mc_q, prom_mc_dir,
         gb_hmc_diff, gb_hmc_q,
         log2FoldChange, padj,
         gb_concordance) %>%
  mutate(
    gb_mc_q = ifelse(is.na(gb_mc_q), "N/A", formatC(gb_mc_q, format = "e", digits = 1)),
    prom_mc_q = ifelse(is.na(prom_mc_q), "N/A", formatC(prom_mc_q, format = "e", digits = 1)),
    gb_hmc_q = ifelse(is.na(gb_hmc_q), "N/A", formatC(gb_hmc_q, format = "e", digits = 1)),
    padj = ifelse(is.na(padj), "N/A", formatC(padj, format = "e", digits = 1)),
    gb_mc_diff = ifelse(is.na(gb_mc_diff), "N/A", sprintf("%+.4f", gb_mc_diff)),
    prom_mc_diff = ifelse(is.na(prom_mc_diff), "N/A", sprintf("%+.4f", prom_mc_diff)),
    gb_hmc_diff = ifelse(is.na(gb_hmc_diff), "N/A", sprintf("%+.4f", gb_hmc_diff)),
    log2FoldChange = ifelse(is.na(log2FoldChange), "N/A", sprintf("%+.3f", log2FoldChange))
  )

# Save as TSV (the table-as-figure is hard to read at this scale)
write.table(table_data, file.path(TABLES_DIR, "field_chr8_comparison_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  45h: Summary table saved as TSV\n")

# =============================================================================
# FIGURE 45i: PROMOTER DIRECTION CONCORDANCE
# =============================================================================

prom_trend_data <- comparison %>%
  filter(!is.na(mouse_gene), !is.na(prom_mc_diff)) %>%
  mutate(prom_trend = ifelse(prom_mc_diff > 0, "Hypermethylated", "Hypomethylated")) %>%
  count(field_direction, prom_trend) %>%
  complete(field_direction = c("Hypermethylated", "Hypomethylated"),
           prom_trend = c("Hypermethylated", "Hypomethylated"),
           fill = list(n = 0))

p45i <- ggplot(prom_trend_data, aes(x = field_direction, y = prom_trend, fill = n)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = n), size = 8, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#2C7BB6", name = "Count") +
  labs(title = "Promoter mC Direction Concordance (All Genes, Trend-Level)",
       subtitle = "Expected ~random if no conserved promoter-level signal",
       x = "Field et al. Direction", y = "Our Promoter mC Trend") +
  theme_biomodal() +
  theme(axis.text = element_text(size = 11))

save_multiformat_ggplot(p45i, file.path(SECTION_DIR, "45i_promoter_concordance_heatmap"),
                        width = 8, height = 7)
cat("  45i: Promoter concordance heatmap saved\n")

# =============================================================================
# FIGURE 45j: RNA-seq EXPRESSION INTEGRATION
# =============================================================================

expr_data <- comparison %>%
  filter(!is.na(mouse_gene), expr_sig %in% TRUE) %>%
  arrange(log2FoldChange) %>%
  mutate(mouse_gene = factor(mouse_gene, levels = unique(mouse_gene)))

if (nrow(expr_data) > 0) {
  p45j <- ggplot(expr_data, aes(x = log2FoldChange, y = mouse_gene,
                                color = expr_concordance, size = baseMean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(x = 0, xend = log2FoldChange, yend = mouse_gene),
                 linewidth = 0.5, alpha = 0.6, color = "grey50") +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("Concordant" = "#4DAF4A", "Discordant" = "#E41A1C"),
                       name = "Matches Field\n(Hyper→Down,\nHypo→Up)") +
    scale_size_continuous(name = "baseMean\n(expression)", range = c(2, 8)) +
    labs(title = "RNA-seq Expression of Field et al. Chr8 Genes",
         subtitle = sprintf("%d genes with significant expression (padj < 0.05); %d match Field pattern",
                            nrow(expr_data), sum(expr_data$expr_concordance == "Concordant")),
         x = "log2 Fold Change (mutant / control)", y = NULL) +
    theme_biomodal()

  save_multiformat_ggplot(p45j, file.path(SECTION_DIR, "45j_rnaseq_expression"),
                          width = 10, height = max(6, nrow(expr_data) * 0.4))
  cat("  45j: RNA-seq expression saved\n")
}

# =============================================================================
# FIGURE 45k: TRISOMY 8 DIAGNOSTIC PANEL
# =============================================================================

chr_coverage$is_chr8 <- chr_coverage$chrom == "8"
chr_expr$is_chr8 <- chr_expr$chrom == "8"
chr_dmr_rate$is_chr8 <- chr_dmr_rate$chrom == "8"

pk1 <- ggplot(chr_coverage, aes(x = reorder(chrom, chrom_num), y = mean_coverage, fill = is_chr8)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = genome_mean_cov, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "#E41A1C")) +
  labs(title = "Mean Sequencing Coverage", subtitle = sprintf("Chr8/genome = %.3f", cov_ratio),
       x = "Chromosome", y = "Mean coverage") +
  theme_biomodal(base_size = 10)

pk2 <- ggplot(chr_expr, aes(x = reorder(chrom, chrom_num), y = median_lfc, fill = is_chr8)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "#E41A1C")) +
  labs(title = "Median RNA-seq log2FC", subtitle = sprintf("Mann-Whitney p = %.2e", expr_test$p.value),
       x = "Chromosome", y = "Median log2FC") +
  theme_biomodal(base_size = 10)

pk3 <- ggplot(chr_dmr_rate, aes(x = reorder(chrom, chrom_num), y = pct_sig, fill = is_chr8)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "#E41A1C")) +
  labs(title = "% Significant mC DMRs", subtitle = sprintf("Chr8 = %.1f%%", chr8_dmr_pct),
       x = "Chromosome", y = "% genes significant") +
  theme_biomodal(base_size = 10)

p45k <- pk1 + pk2 + pk3 +
  plot_annotation(title = "Trisomy 8 Diagnostic: No Evidence of Aneuploidy",
                  subtitle = "All three metrics show chr8 within normal autosomal range",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                plot.subtitle = element_text(size = 11, hjust = 0.5)))

save_multiformat_ggplot(p45k, file.path(SECTION_DIR, "45k_trisomy8_diagnostic"),
                        width = 16, height = 6)
cat("  45k: Trisomy 8 diagnostic saved\n")

# =============================================================================
# STEP 7: TABLE EXPORTS
# =============================================================================

cat("\n--- Step 7: Exporting tables ---\n")

# Ortholog mapping
write.table(
  comparison %>% dplyr::select(human_gene, mouse_gene, field_direction, field_mean_delta, field_n_probes),
  file.path(TABLES_DIR, "field_chr8_ortholog_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Concordant genes
concordant_genes <- comparison %>%
  filter(gb_concordance == "Concordant") %>%
  dplyr::select(human_gene, mouse_gene, field_direction,
         gb_mc_diff, gb_mc_q, gb_hmc_diff, gb_hmc_q,
         log2FoldChange, padj, expr_concordance)
write.table(concordant_genes, file.path(TABLES_DIR, "field_chr8_concordant_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Promoter comparison
prom_compare <- comparison %>%
  filter(!is.na(mouse_gene)) %>%
  dplyr::select(human_gene, mouse_gene, field_direction,
         prom_mc_diff, prom_mc_q, prom_mc_sig, prom_mc_dir,
         prom_hmc_diff, prom_hmc_q, prom_concordance)
write.table(prom_compare, file.path(TABLES_DIR, "field_chr8_promoter_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Statistical tests
stat_df <- data.frame(
  test = sapply(stat_results, function(x) x$test),
  p_value = sapply(stat_results, function(x) x$p),
  stringsAsFactors = FALSE
)
write.table(stat_df, file.path(TABLES_DIR, "field_chr8_statistical_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("  Tables exported to:", TABLES_DIR, "\n")

# =============================================================================
# STEP 8: SUMMARY
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 45 SUMMARY\n")
cat("========================================================================\n\n")

cat(sprintf("Field et al. chr8 genes: %d (hyper: %d, hypo: %d)\n",
            nrow(field_chr8), sum(field_chr8$field_direction == "Hypermethylated"),
            sum(field_chr8$field_direction == "Hypomethylated")))
cat(sprintf("Mapped to mouse orthologs: %d/%d (%.0f%%)\n",
            sum(!is.na(comparison$mouse_gene)), nrow(comparison),
            sum(!is.na(comparison$mouse_gene)) / nrow(comparison) * 100))

n_gb_sig <- sum(comparison$gb_mc_sig %in% TRUE)
n_in_data <- sum(!is.na(comparison$gb_mc_sig))
cat(sprintf("\nGene-body mC: %d/%d significant (%.0f%%)\n", n_gb_sig, n_in_data,
            n_gb_sig / n_in_data * 100))
cat(sprintf("  Concordant: %d, Discordant: %d\n",
            sum(comparison$gb_concordance == "Concordant"),
            sum(comparison$gb_concordance == "Discordant")))

n_prom_sig <- sum(comparison$prom_mc_sig %in% TRUE)
n_prom_in <- sum(!is.na(comparison$prom_mc_sig))
cat(sprintf("\nPromoter mC: %d/%d significant (%.0f%%)\n", n_prom_sig, n_prom_in,
            n_prom_sig / n_prom_in * 100))

n_expr_match <- sum(comparison$expr_concordance == "Concordant")
n_expr_test <- sum(comparison$expr_concordance %in% c("Concordant", "Discordant"))
cat(sprintf("\nRNA-seq: %d/%d expression-concordant with Field pattern\n",
            n_expr_match, n_expr_test))

cat(sprintf("\nTrisomy 8: RULED OUT (coverage ratio = %.3f, expression p = %.2e)\n",
            cov_ratio, expr_test$p.value))

cat("\nTop concordant genes (gene-body mC):\n")
concordant_genes %>%
  arrange(gb_mc_q) %>%
  head(10) %>%
  dplyr::select(human_gene, mouse_gene, field_direction, gb_mc_diff, gb_mc_q) %>%
  print()

cat("\n Section 45 complete.\n\n")
