# scripts/figures/figure1_methylation_phenotype.R
# =============================================================================
# FIGURE 1: The Coordinated Gene-Body Methylation Phenotype  (paper R2 + R3)
#
# Question: What is the DNA-methylation phenotype in BAP1-KO cerebellum?
# Answer:   A coordinated, reciprocal gene-body 5mC-up / 5hmC-down shift with a
#           flat total modified-cytosine pool -- the signature of a blocked
#           TET-mediated active-demethylation pathway -- concentrated in active
#           chromatin and reproducible across replicates and sexes.
#
# Panels (per approved master plan):
#   A  Sec 64  Paired dot-line global bulk shift (5mC / 5hmC / modC), lines
#              connect matched ctrl/mut samples; Cohen's d annotated.
#   B  Sec 03  Horizontal % significant-DMR bar across 6 region classes
#              (gene bodies 51.4% >> promoters 8.3%) from pre-loaded region_dmrs.
#   C  Sec 04  Dual 5mC | 5hmC volcano (mod_difference % vs -log10 q) with
#              KEY_GENES labelled and the FDR threshold dashed.
#   D  Sec 05  mC-vs-hmC quadrant scatter (coordinated_changes.tsv); the
#              coordinated mC-up/hmC-down quadrant (Q4 = 6,589, 78.7%) highlighted.
#   E  Sec 07  Mirror-image effect-size densities of significant 5mC vs 5hmC,
#              vertical dashed lines at the (live-computed) medians.
#   F  Sec 42  Side-by-side 5mC | 5hmC PCA of per-sample gene-body fractions;
#              PC1 = condition, PC2 = sex (R3 sex-difference claim).
#   G  Sec 46  Gviz genome-browser at the Syt1 locus (chr10:108,747,000-
#              108,812,000): merged 5mC/5hmC ctrl+mut and H2AK119ub ctrl+mut
#              signal DataTracks + gene-model track. Rendered via base graphics.
#
# Data sources (all run-specific access via _shared_config.R objects/paths):
#   - Pre-loaded objects: mc_dmr, hmc_dmr (deduped by gene), region_dmrs,
#     SAMPLES, KEY_GENES, COLORS, Q_THRESHOLD
#   - Tables: 64_global_methylation_summary.tsv, 64_statistics.tsv,
#             coordinated_changes.tsv
#   - Per-sample fractions: EXTRACT_PATHS$mc_regional_frac / $hmc_regional_frac
#   - BigWigs: METHYLATION_BIGWIGS (cg_mc/hmc ctrl/mut), HISTONE_BIGWIGS (k119ub)
#
# Units note (VERIFIED on disk): mod_difference / mc_diff / hmc_diff are
# FRACTIONS on [-1, 1] (e.g. Syt1 = 0.1731 == "+17.3%"). They are multiplied by
# 100 for percent-axis display in panels C/D/E. The Section-64 summary table is
# ALREADY in percent and is used as-is in panel A.
#
# Output: plots/figures/figure1_methylation_phenotype/  (PDF + SVG + JPG)
#         plus the standalone Gviz sub-panel
#         plots/figures/figure1G_syt1_browser/          (PDF + SVG + JPG)
# =============================================================================

# --- Universal sourcing contract: _shared_config.R first, then _figure_config.R
# Both resolve relative to getwd() == downstream/.
source("scripts/viz_sections/_shared_config.R")
source("scripts/figures/_figure_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

cat("\n")
cat("================================================================================\n")
cat("FIGURE 1: Coordinated Gene-Body Methylation Phenotype (R2 + R3)\n")
cat("================================================================================\n\n")

# Fail loudly if the load-bearing pre-loaded objects are absent.
stopifnot(
  exists("mc_dmr"), exists("hmc_dmr"), exists("region_dmrs"),
  exists("SAMPLES"), exists("KEY_GENES"), exists("COLORS"),
  !is.null(mc_dmr), !is.null(hmc_dmr)
)

# Canonical per-modality colours (5mC red / 5hmC blue) and direction colours.
MC_COL   <- PUB_COLORS$mC      # "#E41A1C"
HMC_COL  <- PUB_COLORS$hmC     # "#377EB8"
HYPER_COL <- PUB_COLORS$Hyper  # "#D7191C"
HYPO_COL  <- PUB_COLORS$Hypo   # "#2C7BB6"

# =============================================================================
# PANEL A -- Section 64: paired global bulk-methylation shift
# =============================================================================
# Per-sample autosomal CG percentages (already in %) + paired t-test Cohen's d.
cat("Panel A: global methylation shift (Section 64)...\n")

global_summary <- read_table_tsv("64_global_methylation_summary.tsv")
global_stats   <- read_table_tsv("64_statistics.tsv")

required_global_cols <- c("sample_id", "condition", "sex", "batch",
                          "5mC_pct", "5hmC_pct", "modC_pct")
if (!all(required_global_cols %in% colnames(global_summary))) {
  stop("Panel A: 64_global_methylation_summary.tsv missing expected columns: ",
       paste(setdiff(required_global_cols, colnames(global_summary)), collapse = ", "))
}
if (!all(c("modality", "delta_mean", "cohens_d") %in% colnames(global_stats))) {
  stop("Panel A: 64_statistics.tsv missing modality/delta_mean/cohens_d columns.")
}

# Long form: one row per (sample, modality). Pair key = sex+batch links the
# matched control/mutant samples so the connecting line is biologically honest.
global_long <- global_summary %>%
  dplyr::transmute(
    sample_id, condition, sex, batch,
    pair = paste(sex, batch, sep = "_"),
    `5mC`  = `5mC_pct`,
    `5hmC` = `5hmC_pct`,
    modC   = `modC_pct`
  ) %>%
  tidyr::pivot_longer(c(`5mC`, `5hmC`, modC),
                      names_to = "modality", values_to = "pct")

modality_levels <- c("5mC", "5hmC", "modC")
global_long$modality  <- factor(global_long$modality, levels = modality_levels)
global_long$condition <- factor(global_long$condition, levels = c("Control", "Mutant"))

# Condition means per modality (drawn as wide cross-bars behind the points).
global_means <- global_long %>%
  dplyr::group_by(modality, condition) %>%
  dplyr::summarise(mean_pct = mean(pct), .groups = "drop")

# Cohen's d + delta annotation strip, one label per modality facet.
stats_lab <- global_stats %>%
  dplyr::filter(modality %in% modality_levels) %>%
  dplyr::transmute(
    modality = factor(modality, levels = modality_levels),
    label = sprintf("Delta = %+.2f%%\nCohen's d = %.2f", delta_mean, cohens_d)
  )

panel_a <- ggplot(global_long, aes(x = condition, y = pct)) +
  geom_line(aes(group = pair), colour = "grey60", linewidth = 0.3) +
  geom_crossbar(
    data = global_means,
    aes(x = condition, y = mean_pct, ymin = mean_pct, ymax = mean_pct,
        colour = condition),
    width = 0.55, fatten = 0, linewidth = 0.6, show.legend = FALSE
  ) +
  geom_point(aes(colour = condition, shape = sex), size = 1.6, stroke = 0.6) +
  geom_text(
    data = stats_lab,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.2, size = 2.1, lineheight = 0.9
  ) +
  facet_wrap(~ modality, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = COLORS$condition, name = "Condition") +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Sex") +
  labs(
    title = "Global CpG methylation shift",
    x = NULL, y = "Autosomal CpG (%)"
  ) +
  theme_pub() +
  theme(legend.position = "right")

# =============================================================================
# PANEL B -- Section 03: % significant DMRs by region class
# =============================================================================
# Replicates region_dmrs semantics exactly: sum(significant) / nrow per region.
# NOTE: "Gene bodies" is deduped-by-gene (per-gene), the other five regions are
# raw per-feature counts (this is the established Section 03 / region_dmrs
# behaviour; the relative dominance of gene bodies is the message).
cat("Panel B: %-significant DMRs by region (Section 03 / region_dmrs)...\n")

region_order <- c("Gene bodies", "CpG shores", "CpG shelves",
                  "Promoters", "CpG islands", "TSS regions")
region_rows <- lapply(region_order, function(rn) {
  df <- region_dmrs[[rn]]
  if (is.null(df)) {
    stop("Panel B: region_dmrs[['", rn, "']] is NULL -- expected all 6 regions.")
  }
  data.frame(
    region  = rn,
    n_sig   = sum(df$significant, na.rm = TRUE),
    n_total = nrow(df),
    stringsAsFactors = FALSE
  )
})
region_df <- do.call(rbind, region_rows)
region_df$pct <- 100 * region_df$n_sig / region_df$n_total
region_df$region <- factor(region_df$region, levels = rev(region_order))
# Region class fill: gene bodies highlighted as the dominant class.
region_df$is_genebody <- region_df$region == "Gene bodies"

panel_b <- ggplot(region_df, aes(x = pct, y = region, fill = is_genebody)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.3) +
  scale_fill_manual(values = c("TRUE" = MC_COL, "FALSE" = "grey60")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     limits = c(0, max(region_df$pct) * 1.18)) +
  labs(
    title = "Significant 5mC DMRs by region",
    x = "Significant DMRs (%)", y = NULL
  ) +
  theme_pub()

# =============================================================================
# PANEL C -- Section 04: dual 5mC | 5hmC volcano
# =============================================================================
# Pre-loaded deduped mc_dmr / hmc_dmr. mod_difference is a FRACTION -> x100.
cat("Panel C: dual volcano (Section 04)...\n")

build_volcano_df <- function(dmr, meth_label) {
  d <- dmr
  d$diff_pct <- d$mod_difference * 100
  d$status <- dplyr::case_when(
    !d$significant ~ "Not Significant",
    d$mod_difference > 0 ~ "Hypermethylated",
    TRUE ~ "Hypomethylated"
  )
  d$status <- factor(d$status,
                     levels = c("Hypermethylated", "Hypomethylated", "Not Significant"))
  d$methylation <- meth_label
  d
}

mc_volc  <- build_volcano_df(mc_dmr,  "5mC")
hmc_volc <- build_volcano_df(hmc_dmr, "5hmC")
volc_df  <- dplyr::bind_rows(mc_volc, hmc_volc)
volc_df$methylation <- factor(volc_df$methylation, levels = c("5mC", "5hmC"))

n_sig_mc  <- sum(mc_dmr$significant,  na.rm = TRUE)
n_sig_hmc <- sum(hmc_dmr$significant, na.rm = TRUE)

# KEY_GENES labels (one row per gene per facet, only if present & significant).
key_label_df <- volc_df %>%
  dplyr::filter(gene %in% KEY_GENES, significant) %>%
  dplyr::group_by(methylation, gene) %>%
  dplyr::slice_max(neg_log10_q, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# Per-facet significance counts for subtitles.
facet_counts <- data.frame(
  methylation = factor(c("5mC", "5hmC"), levels = c("5mC", "5hmC")),
  label = c(sprintf("%s significant (q < %.2f)", format(n_sig_mc, big.mark = ","), Q_THRESHOLD),
            sprintf("%s significant (q < %.2f)", format(n_sig_hmc, big.mark = ","), Q_THRESHOLD))
)

volcano_colors <- c("Hypermethylated" = HYPER_COL,
                    "Hypomethylated"  = HYPO_COL,
                    "Not Significant" = "grey80")

panel_c <- ggplot(volc_df, aes(x = diff_pct, y = neg_log10_q, colour = status)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed",
             colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "grey40", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = key_label_df, aes(label = gene),
    size = 2, colour = "black", fontface = "italic",
    max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.2, box.padding = 0.3, show.legend = FALSE
  ) +
  geom_text(
    data = facet_counts, aes(x = 0, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.3, size = 2.1
  ) +
  facet_wrap(~ methylation, nrow = 1) +
  scale_colour_manual(values = volcano_colors, name = "Direction",
                      drop = FALSE) +
  labs(
    title = "Per-gene differential methylation",
    x = "Methylation change (mutant - control, %)",
    y = expression(-log[10]~"(q-value)")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  theme_pub()

# =============================================================================
# PANEL D -- Section 05: mC-vs-hmC quadrant scatter
# =============================================================================
# coordinated_changes.tsv: co-significant genes (both modalities). mc_diff /
# hmc_diff are FRACTIONS -> x100. Coordinated = mC-up & hmC-down (Q4).
cat("Panel D: coordinated quadrant scatter (Section 05)...\n")

coord <- read_table_tsv("coordinated_changes.tsv")
required_coord_cols <- c("gene", "mc_diff", "hmc_diff", "coordinated_pattern")
if (!all(required_coord_cols %in% colnames(coord))) {
  stop("Panel D: coordinated_changes.tsv missing columns: ",
       paste(setdiff(required_coord_cols, colnames(coord)), collapse = ", "))
}

coord$mc_pct  <- coord$mc_diff  * 100
coord$hmc_pct <- coord$hmc_diff * 100
# Quadrant assignment (sign-based); Q4 = mC up / hmC down = the TET-block class.
coord$quadrant <- dplyr::case_when(
  coord$mc_diff > 0 & coord$hmc_diff < 0 ~ "mC up / hmC down",
  coord$mc_diff > 0 & coord$hmc_diff > 0 ~ "mC up / hmC up",
  coord$mc_diff < 0 & coord$hmc_diff > 0 ~ "mC down / hmC up",
  TRUE                                   ~ "mC down / hmC down"
)
coord$quadrant <- factor(coord$quadrant, levels = c(
  "mC up / hmC down", "mC down / hmC up", "mC up / hmC up", "mC down / hmC down"))

n_coord_true  <- sum(coord$coordinated_pattern, na.rm = TRUE)
n_coord_total <- nrow(coord)
pct_coord     <- 100 * n_coord_true / n_coord_total

quadrant_colors <- c(
  "mC up / hmC down"  = HYPER_COL,   # the coordinated TET-block quadrant
  "mC down / hmC up"  = HYPO_COL,
  "mC up / hmC up"    = "grey55",
  "mC down / hmC down"= "grey75"
)

coord_key <- coord %>% dplyr::filter(gene %in% KEY_GENES)

coord_annot <- sprintf("Coordinated mC up / hmC down\n%s of %s (%.1f%%)",
                       format(n_coord_true, big.mark = ","),
                       format(n_coord_total, big.mark = ","), pct_coord)

panel_d <- ggplot(coord, aes(x = mc_pct, y = hmc_pct, colour = quadrant)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  geom_point(size = 0.4, alpha = 0.45) +
  ggrepel::geom_text_repel(
    data = coord_key, aes(label = gene),
    size = 2, colour = "black", fontface = "italic",
    max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.2, box.padding = 0.3, show.legend = FALSE
  ) +
  annotate("text", x = Inf, y = -Inf, label = coord_annot,
           hjust = 1.05, vjust = -0.4, size = 2.1, lineheight = 0.9) +
  scale_colour_manual(values = quadrant_colors, name = "Quadrant") +
  labs(
    title = "Coordinated 5mC / 5hmC changes",
    x = "5mC change (%)", y = "5hmC change (%)"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  theme_pub()

# =============================================================================
# PANEL E -- Section 07: mirror-image effect-size densities (significant only)
# =============================================================================
# Computed live from the deduped significant DMRs; medians annotated with the
# live-computed values (never hardcoded) so the labels always match the data.
cat("Panel E: effect-size densities (Section 07)...\n")

mc_sig  <- mc_dmr  %>% dplyr::filter(significant)
hmc_sig <- hmc_dmr %>% dplyr::filter(significant)

eff_df <- dplyr::bind_rows(
  data.frame(diff_pct = mc_sig$mod_difference  * 100, methylation = "5mC"),
  data.frame(diff_pct = hmc_sig$mod_difference * 100, methylation = "5hmC")
)
eff_df$methylation <- factor(eff_df$methylation, levels = c("5mC", "5hmC"))

median_mc  <- median(mc_sig$mod_difference  * 100, na.rm = TRUE)
median_hmc <- median(hmc_sig$mod_difference * 100, na.rm = TRUE)

median_lines <- data.frame(
  methylation = factor(c("5mC", "5hmC"), levels = c("5mC", "5hmC")),
  med = c(median_mc, median_hmc)
)
eff_annot <- sprintf("Median 5mC %+.2f%%\nMedian 5hmC %+.2f%%", median_mc, median_hmc)

panel_e <- ggplot(eff_df, aes(x = diff_pct, fill = methylation, colour = methylation)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  geom_vline(data = median_lines, aes(xintercept = med, colour = methylation),
             linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.3) +
  annotate("text", x = Inf, y = Inf, label = eff_annot,
           hjust = 1.05, vjust = 1.3, size = 2.1, lineheight = 0.9) +
  scale_fill_manual(values = c("5mC" = MC_COL, "5hmC" = HMC_COL),
                    name = "Modification") +
  scale_colour_manual(values = c("5mC" = MC_COL, "5hmC" = HMC_COL),
                      name = "Modification") +
  labs(
    title = "Effect-size distribution (significant DMRs)",
    x = "Methylation change (%)", y = "Density"
  ) +
  theme_pub()

# =============================================================================
# PANEL F -- Section 42: PCA of per-sample gene-body methylation (5mC | 5hmC)
# =============================================================================
# Reads the per-sample regional-fraction matrices with check.names = FALSE to
# preserve the ORIGINAL sample-id column names, then derives condition/sex/batch
# by joining to the canonical SAMPLES table. (The section_42 grepl('_F', ...)
# sex test breaks once read.table mangles '-' to '.'; joining to SAMPLES avoids
# that latent bug so the R3 sex axis is correct.)
cat("Panel F: per-sample methylation PCA (Section 42)...\n")

mc_frac_file  <- EXTRACT_PATHS$mc_regional_frac
hmc_frac_file <- EXTRACT_PATHS$hmc_regional_frac
if (!file.exists(mc_frac_file) || !file.exists(hmc_frac_file)) {
  stop("Panel F: per-sample regional-fraction file(s) missing:\n  ",
       mc_frac_file, "\n  ", hmc_frac_file)
}

mc_frac  <- utils::read.table(mc_frac_file,  header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE)
hmc_frac <- utils::read.table(hmc_frac_file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE)

# Map each per-sample column to its sample_id by stripping the modality suffix,
# then to condition/sex/batch via SAMPLES (fail loudly on any mismatch).
build_sample_meta <- function(frac_cols, suffix) {
  ids <- sub(suffix, "", frac_cols)
  meta <- SAMPLES[match(ids, SAMPLES$sample_id), ]
  if (anyNA(meta$sample_id)) {
    stop("Panel F: regional-fraction columns did not map to SAMPLES$sample_id: ",
         paste(frac_cols[is.na(meta$sample_id)], collapse = ", "))
  }
  data.frame(
    col_name  = frac_cols,
    sample_id = meta$sample_id,
    condition = factor(meta$condition, levels = c("Control", "Mutant")),
    sex       = factor(meta$sex, levels = c("Female", "Male")),
    batch     = factor(paste("Batch", meta$batch), levels = c("Batch 1", "Batch 2")),
    short     = meta$short_name,
    stringsAsFactors = FALSE
  )
}

mc_cols  <- colnames(mc_frac)[grepl("_num_mc_region_frac$",  colnames(mc_frac))]
hmc_cols <- colnames(hmc_frac)[grepl("_num_hmc_region_frac$", colnames(hmc_frac))]
if (length(mc_cols) != 8 || length(hmc_cols) != 8) {
  stop("Panel F: expected 8 per-sample columns each; got ",
       length(mc_cols), " (mC) / ", length(hmc_cols), " (hmC).")
}

mc_meta  <- build_sample_meta(mc_cols,  "_num_mc_region_frac")
hmc_meta <- build_sample_meta(hmc_cols, "_num_hmc_region_frac")

# PCA over all gene bodies: genes x samples -> prcomp on samples (t()).
run_methylation_pca <- function(frac, sample_cols, meta, meth_label) {
  mat <- as.matrix(frac[, sample_cols, drop = FALSE])
  # NA -> per-gene mean (matches section_42 handling).
  for (i in seq_len(nrow(mat))) {
    na_idx <- is.na(mat[i, ])
    if (any(na_idx) && !all(na_idx)) mat[i, na_idx] <- mean(mat[i, !na_idx])
  }
  gene_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  mat <- mat[!is.na(gene_vars) & gene_vars > 1e-10, , drop = FALSE]

  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_explained <- summary(pca)$importance[2, ] * 100

  pca_df <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    condition = meta$condition, sex = meta$sex,
    batch = meta$batch, label = meta$short
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, shape = sex)) +
    geom_point(size = 2.2, stroke = 0.8) +
    ggrepel::geom_text_repel(aes(label = label), size = 1.9,
                             show.legend = FALSE, max.overlaps = Inf,
                             min.segment.length = 0, segment.size = 0.2) +
    scale_colour_manual(values = COLORS$condition, name = "Condition") +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Sex") +
    labs(
      title = sprintf("%s PCA", meth_label),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_pub()
}

pca_mc  <- run_methylation_pca(mc_frac,  mc_cols,  mc_meta,  "5mC")
pca_hmc <- run_methylation_pca(hmc_frac, hmc_cols, hmc_meta, "5hmC")
# Collect shared legend on the right of the pair. The 5mC|5hmC pair is itself a
# nested patchwork; wrap_elements() makes it ONE atomic panel so the top-level
# tag_levels = "A" assigns a single 'F' to it rather than tagging both PCA
# subplots (which would break the A-G panel-letter sequence).
panel_f_pair <- (pca_mc | pca_hmc) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")
panel_f <- patchwork::wrap_elements(panel_f_pair)

# =============================================================================
# PANEL G -- Section 46: Gviz genome-browser at the Syt1 locus
# =============================================================================
# Base-graphics render via save_multiformat_base. BigWigs are chr-prefixed
# (verified); the window is the plan's chr10:108,747,000-108,812,000.
cat("Panel G: Syt1 genome-browser track (Section 46)...\n")

suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
})

SYT1_CHR   <- "chr10"
SYT1_START <- 108747000
SYT1_END   <- 108812000

browser_bw <- list(
  mc_ctrl  = METHYLATION_BIGWIGS$cg_mc_ctrl,
  mc_mut   = METHYLATION_BIGWIGS$cg_mc_mut,
  hmc_ctrl = METHYLATION_BIGWIGS$cg_hmc_ctrl,
  hmc_mut  = METHYLATION_BIGWIGS$cg_hmc_mut,
  k119_ctrl = HISTONE_BIGWIGS$k119ub_ctrl,
  k119_mut  = HISTONE_BIGWIGS$k119ub_mut
)
missing_bw <- names(browser_bw)[!file.exists(unlist(browser_bw))]
if (length(missing_bw) > 0) {
  stop("Panel G: required BigWig(s) not found: ",
       paste(missing_bw, collapse = ", "),
       "\nCheck METHYLATION_BIGWIGS / HISTONE_BIGWIGS in _shared_config.R.")
}

make_bw_track <- function(path, name, col, ymax = NULL) {
  Gviz::DataTrack(
    range = path, genome = "mm10", chromosome = SYT1_CHR,
    name = name, type = "histogram", col.histogram = col,
    fill.histogram = col, col.axis = "black", col.title = "black",
    background.title = "white", fontsize = 9, cex.axis = 0.6, cex.title = 0.65,
    ylim = if (is.null(ymax)) NULL else c(0, ymax)
  )
}

g_track  <- Gviz::GenomeAxisTrack(fontsize = 9, cex = 0.6)
# Gene model from the mm10 TxDb (loaded by _shared_config.R).
gene_track <- Gviz::GeneRegionTrack(
  TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
  genome = "mm10", chromosome = SYT1_CHR,
  start = SYT1_START, end = SYT1_END,
  name = "Syt1", showId = FALSE, geneSymbol = FALSE,
  fill = "grey40", col = "grey40",
  background.title = "white", col.title = "black", col.axis = "black",
  fontsize = 9, cex.title = 0.65
)

mc_ctrl_t  <- make_bw_track(browser_bw$mc_ctrl,  "5mC ctrl",   MC_COL)
mc_mut_t   <- make_bw_track(browser_bw$mc_mut,   "5mC mut",    MC_COL)
hmc_ctrl_t <- make_bw_track(browser_bw$hmc_ctrl, "5hmC ctrl",  HMC_COL)
hmc_mut_t  <- make_bw_track(browser_bw$hmc_mut,  "5hmC mut",   HMC_COL)
k119_ctrl_t <- make_bw_track(browser_bw$k119_ctrl, "H2AK119ub ctrl", "#756BB1")
k119_mut_t  <- make_bw_track(browser_bw$k119_mut,  "H2AK119ub mut",  "#756BB1")

browser_tracks <- list(
  g_track,
  mc_ctrl_t, mc_mut_t,
  hmc_ctrl_t, hmc_mut_t,
  k119_ctrl_t, k119_mut_t,
  gene_track
)

# Quoted plot expression (base graphics) for the multi-format saver.
syt1_plot_expr <- quote({
  Gviz::plotTracks(
    browser_tracks,
    from = SYT1_START, to = SYT1_END, chromosome = SYT1_CHR,
    main = "Syt1 (chr10:108,747,000-108,812,000)  |  +17.3% 5mC, -15.0% 5hmC",
    cex.main = 0.8,
    sizes = c(0.8, 1, 1, 1, 1, 1, 1, 0.8)
  )
})

save_multiformat_base(
  syt1_plot_expr,
  base_path = file.path(FIGURE_DIR, "figure1G_syt1_browser"),
  width = 7.1, height = 5.5, dpi = 300, verbose = TRUE, use_subfolders = TRUE
)

# Wrap the rendered browser image as a patchwork element so it can compose into
# the main figure beneath the ggplot panels.
panel_g <- patchwork::wrap_elements(full = grid::grid.grabExpr(
  eval(syt1_plot_expr)
))

# =============================================================================
# COMPOSE -- master plan layout
# =============================================================================
# design = "AABB\nCCCC\nDDEE\nFFFF\nGGGG", heights = c(1, 1.3, 1, 1, 1.4)
cat("Composing Figure 1 (7-panel patchwork)...\n")

design <- "AABB
CCCC
DDEE
FFFF
GGGG"

figure1 <- panel_a + panel_b + panel_c + panel_d + panel_e + panel_f + panel_g +
  patchwork::plot_layout(design = design, heights = c(1, 1.3, 1, 1, 1.4))
# Single plot_annotation call carries BOTH the panel tag sequence and the title:
# chaining two plot_annotation() calls can drop the first's tag_levels, so the
# title and tag_levels = "A" are combined here rather than via add_panel_labels().
figure1 <- figure1 +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Figure 1. Coordinated gene-body 5mC up / 5hmC down phenotype in BAP1-KO cerebellum",
    theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))
  )

# Dimensions: 180 x 290 mm (per plan).
save_figure(figure1, "figure1_methylation_phenotype",
            width_mm = 180, height_mm = 290)

cat("\nFigure 1 complete.\n")
cat(sprintf("  Panel A: delta 5mC %+.2f%% / 5hmC %+.2f%% (Cohen's d %.2f / %.2f)\n",
            global_stats$delta_mean[global_stats$modality == "5mC"],
            global_stats$delta_mean[global_stats$modality == "5hmC"],
            global_stats$cohens_d[global_stats$modality == "5mC"],
            global_stats$cohens_d[global_stats$modality == "5hmC"]))
cat(sprintf("  Panel B: gene bodies %.1f%% vs promoters %.1f%% significant\n",
            region_df$pct[region_df$region == "Gene bodies"],
            region_df$pct[region_df$region == "Promoters"]))
cat(sprintf("  Panel C: %s 5mC / %s 5hmC significant\n",
            format(n_sig_mc, big.mark = ","), format(n_sig_hmc, big.mark = ",")))
cat(sprintf("  Panel D: %s coordinated (%.1f%%)\n",
            format(n_coord_true, big.mark = ","), pct_coord))
cat(sprintf("  Panel E: median 5mC %+.2f%% / 5hmC %+.2f%%\n", median_mc, median_hmc))
cat("  Panels F/G rendered (PCA + Syt1 browser).\n")
cat("================================================================================\n")
