# Permutation Testing Reference

Ground-truth specification for implementing permutation-based enrichment tests across the biomodal DUET evoC differential methylation visualization pipeline. Complements existing Fisher's exact tests (sections 12-33) with genomic permutation tests to address reviewer concerns about spatial non-independence of genomic features.

**Target scripts:** `section_34.R` through `section_37.R` (new files in `scripts/viz_sections/`)
**Dependencies:** `regioneR`, `regioneReloaded` (Bioconductor)
**Local vignette reference:** `downstream/docs/regioneReloaded.Rmd`
**Estimated total runtime:** 2-4 hours (5,000 permutations, 8 cores on HPC)

---

## 1. Motivation

### The Problem: Genomic Non-Independence

Fisher's exact test assumes each observation is independent. In genomic enrichment analyses, this assumption is violated:

- **Adjacent DMRs cluster.** Ten adjacent hypermethylated DMRs all overlapping K119ub-gained peaks are not 10 independent observations — they are likely one locus. Fisher's treats them as 10.
- **Polycomb domains span megabases.** 50 genes in the same Polycomb domain all sharing H3K27me3 and hypermethylation are not 50 independent confirmations.
- **Compartments are megabase-scale.** All genes in the same A compartment share the designation by construction — gene-level Fisher's massively overcounts degrees of freedom.
- **Co-regulated gene neighborhoods exist.** Genes in the same TAD or regulatory neighborhood tend to share chromatin mark directions.

The consequence: Fisher p-values are anti-conservative (too small). O/E ratios as point estimates are unaffected, but our confidence in whether O/E=1.5 is real or noise-driven is inflated.

### The Solution: Permutation Testing

Permutation tests generate an empirical null distribution by randomizing the spatial relationship between features while preserving their genomic structure (interval sizes, chromosome assignment, count). The observed statistic's position in this null distribution gives an honest p-value.

### Two Methodologies Required

| Approach | Use Case | Tool |
|----------|----------|------|
| **Genomic interval permutation** | Testing spatial co-occurrence of region sets (DMR intervals x ChIP peaks, ATAC peaks x chromatin states) | `regioneReloaded::crosswisePermTest()` |
| **Chromosome-stratified label shuffle** | Testing label co-occurrence at gene level (2x2 O/E tables: methylation direction x mark direction) | Custom R implementation |

regioneReloaded cannot be used for gene-level 2x2 tables because the unit there is a gene label (e.g., "mC Up" vs "mC Down"), not a genomic coordinate. The appropriate permutation for gene-level tests shuffles labels within chromosomes to preserve chromosome-level confounding structure.

---

## 2. Package Installation & Genome Setup

### Installation

```r
# regioneReloaded pulls regioneR as a dependency
BiocManager::install("regioneReloaded")

# Verify
library(regioneR)
library(regioneReloaded)
packageVersion("regioneR")         # >= 1.34.0
packageVersion("regioneReloaded")  # >= 1.4.0
```

Neither package is currently in the conda `mariner_env` or `modality` environments. Install on HPC before running sections 34-36.

### Genome Object Construction

```r
# Standard chromosomes only (exclude random/Un contigs)
library(BSgenome.Mmusculus.UCSC.mm10)
genome_full <- getGenomeAndMask("mm10")$genome
standard_chrs <- paste0("chr", c(1:19, "X"))
genome <- genome_full[seqnames(genome_full) %in% standard_chrs]
```

The TxDb (`TxDb.Mmusculus.UCSC.mm10.knownGene`) is already loaded by `_shared_config.R`.

### Shared Permutation Parameters

Every new section should define these at the top:

```r
PERM_NTIMES   <- 5000    # Minimum for publication (regioneReloaded recommends >= 5000)
PERM_CORES    <- 8       # mc.cores for parallel permutation
PERM_PER_CHR  <- TRUE    # Circular permutation within chromosomes
PERM_RANFUN   <- "randomizeRegions"  # Preserves interval sizes exactly
PERM_EVFUN    <- "numOverlaps"       # Default evaluation function
PERM_SEED     <- 42      # Reproducibility
```

**Why `randomizeRegions`?** It performs circular permutation within each chromosome, preserving interval count, size distribution, and chromosome assignment. This is more conservative than `resampleGenome` (which ignores chromosome structure) and more appropriate than `resampleRegions` (which requires a pre-defined universe).

**Why `per.chromosome = TRUE`?** Chromosome-level effects (X-linked genes, gene density gradients, chromosome-specific methylation patterns) would inflate false positives if regions could be shuffled across chromosomes.

---

## 3. Section 34 — DMR x Chromatin Mark Interval Permutation

**File:** `scripts/viz_sections/section_34_permutation_dmr_chromatin_marks.R`
**Validates Fisher tests from:** sections 12a, 14a, 14b, 19f

### Biological Question

Do differentially methylated regions (DMR intervals) spatially co-occur with chromatin mark peaks more than expected by chance, accounting for the clustered nature of DMRs in the genome?

### Region Set Construction

**Alist (2 DMR direction sets):**

```r
# From _shared_config.R: mc_dmr pre-loaded, dmr_to_granges() available
hyper_gr <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant, mod_difference > 0))
hypo_gr  <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant, mod_difference < 0))

Alist <- list(
  "mC Hyper DMRs" = hyper_gr,
  "mC Hypo DMRs"  = hypo_gr
)
```

**Blist (8 chromatin mark peak sets):**

```r
# From _shared_config.R: ATAC_FILES, K119UB_FILES, H3K27AC_FILES, load_chip_peaks()

# Differential ATAC peaks (pre-computed as BED files)
atac_up   <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")

# Condition-specific K119ub peaks
k119ub_ctrl <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Ctrl")
k119ub_mut  <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mut")

# Derived gained/lost K119ub peaks (pattern from section_14, lines 34-47)
# derive_differential_peaks(): gained = mut-only, lost = ctrl-only
diff_k119ub <- derive_differential_peaks(k119ub_ctrl, k119ub_mut)
k119ub_gained <- diff_k119ub$gained
k119ub_lost   <- diff_k119ub$lost

# Condition-specific H3K27ac peaks
h3k27ac_ctrl <- load_chip_peaks(H3K27AC_FILES$ctrl, "H3K27ac Ctrl")
h3k27ac_mut  <- load_chip_peaks(H3K27AC_FILES$mut, "H3K27ac Mut")

Blist <- list(
  "ATAC Up"       = atac_up,
  "ATAC Down"     = atac_down,
  "K119ub Ctrl"   = k119ub_ctrl,
  "K119ub Mut"    = k119ub_mut,
  "K119ub Gained" = k119ub_gained,
  "K119ub Lost"   = k119ub_lost,
  "H3K27ac Ctrl"  = h3k27ac_ctrl,
  "H3K27ac Mut"   = h3k27ac_mut
)
```

### regioneReloaded Call

```r
set.seed(PERM_SEED)
cw_34 <- crosswisePermTest(
  Alist          = Alist,
  Blist          = Blist,
  genome         = genome,
  ranFUN         = PERM_RANFUN,
  evFUN          = PERM_EVFUN,
  ntimes         = PERM_NTIMES,
  mc.cores       = PERM_CORES,
  per.chromosome = PERM_PER_CHR
)
cw_34 <- makeCrosswiseMatrix(cw_34, pvcut = 1)
```

### Figures

| Figure | Description | Function |
|--------|-------------|----------|
| **34a** | Crosswise z-score heatmap (2 x 8 matrix). Blue = depleted, red = enriched. | `plotCrosswiseMatrix(cw_34, matrix_type = "association")` |
| **34b** | Fisher-vs-permutation comparison table. For each original Fisher test: Fisher OR, Fisher p, permutation z-score, permutation p, concordance. | Custom ggplot table/forest plot |
| **34c** | Local z-score curve for strongest association (likely Hyper DMR x K119ub Gained) at +/- 50kb window, showing focal vs broad enrichment. | `multiLocalZscore()` then `plotSingleLZ()` |

### Mapping to Existing Fisher Tests

| Crosswise cell | Original section | Original figure | Description |
|----------------|------------------|-----------------|-------------|
| mC Hyper x ATAC Up | 12 | 12a | Hyper DMR overlap with ATAC-up peaks |
| mC Hyper x ATAC Down | 12 | 12a | Hyper DMR overlap with ATAC-down peaks |
| mC Hypo x ATAC Up | 12 | 12a | Hypo DMR overlap with ATAC-up peaks |
| mC Hypo x ATAC Down | 12 | 12a | Hypo DMR overlap with ATAC-down peaks |
| mC Hyper x K119ub Ctrl | 14 | 14a | Hyper DMR overlap with ctrl K119ub |
| mC Hyper x K119ub Mut | 14 | 14a | Hyper DMR overlap with mut K119ub |
| mC Hypo x K119ub Ctrl | 14 | 14a | Hypo DMR overlap with ctrl K119ub |
| mC Hypo x K119ub Mut | 14 | 14a | Hypo DMR overlap with mut K119ub |
| mC Hyper x K119ub Gained | 14 | 14b | Hyper DMR overlap with gained K119ub |
| mC Hyper x K119ub Lost | 14 | 14b | Hyper DMR overlap with lost K119ub |
| mC Hypo x K119ub Gained | 14 | 14b | Hypo DMR overlap with gained K119ub |
| mC Hypo x K119ub Lost | 14 | 14b | Hypo DMR overlap with lost K119ub |
| mC Hyper x H3K27ac Ctrl | 19 | 19f | Hyper DMR overlap with ctrl H3K27ac |
| mC Hyper x H3K27ac Mut | 19 | 19f | Hyper DMR overlap with mut H3K27ac |
| mC Hypo x H3K27ac Ctrl | 19 | 19f | Hypo DMR overlap with ctrl H3K27ac |
| mC Hypo x H3K27ac Mut | 19 | 19f | Hypo DMR overlap with mut H3K27ac |

### Runtime Estimate
~20-30 minutes at 5,000 permutations, 8 cores.

---

## 4. Section 35 — ATAC x Chromatin Features + Loop Anchor Permutation

**File:** `scripts/viz_sections/section_35_permutation_atac_loops.R`
**Validates Fisher tests from:** sections 13b, 13c, 13d, 27a, 27b, 27e, 31a

### Sub-analysis 35A: ATAC Peaks x ChIP Marks (validates 13b)

**Question:** Are differential ATAC peaks enriched at specific histone mark peaks?

```r
Alist_35a <- list(
  "ATAC Up"   = load_chip_peaks(ATAC_FILES$up, "ATAC Up"),
  "ATAC Down" = load_chip_peaks(ATAC_FILES$down, "ATAC Down")
)

Blist_35a <- list(
  "CTCF"      = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
  "H3K27ac"   = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
  "H3K27me3"  = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  "H3K4me1"   = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
  "H3K4me3"   = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
  "Bivalent"  = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)
```

One `crosswisePermTest()` call: 2 x 6 = 12 pairwise tests.

**Mapping:** Each cell validates one of the 6 Fisher tests in section 13b (ATAC direction x ChIP mark overlap, BH-corrected).

### Sub-analysis 35B: ATAC Peaks x 7 Chromatin States (validates 13c)

**Question:** Are ATAC-up vs ATAC-down peaks differentially distributed across chromatin states?

Construct 7 GRanges objects, one per chromatin state, from the classified genome:

```r
# Get all gene body regions with chromatin state assignments
# Use the classify_chromatin_state() function from _shared_config.R
# on all mc_dmr genes (not just significant ones) to define state regions

# For each state in CHROMATIN_STATE_ORDER, collect the DMR GRanges with that state
all_dmr_gr <- dmr_to_granges(mc_dmr)
chip_peaks <- lapply(CHIP_PEAK_FILES, load_chip_peaks)
overlaps <- compute_chip_overlaps(all_dmr_gr, chip_peaks)
# Need distance_to_tss from ChIPseeker for classification
# (recompute within section, same pattern as section_10)

Blist_35b <- list()
for (state in CHROMATIN_STATE_ORDER) {
  state_mask <- chromatin_states == state
  Blist_35b[[state]] <- all_dmr_gr[state_mask]
}
```

One `crosswisePermTest()` call: 2 x 7 = 14 pairwise tests.

**Mapping:** Validates the 2x7 O/E heatmap in figure 13c.

### Sub-analysis 35C: Loop Anchors x Chromatin Features (validates 13d, 27a/b/e, 31a)

**Question:** Do gained/lost loop anchors show enriched overlap with specific chromatin features?

Construct loop anchor GRanges from `LOOP_FILES$late`:

```r
# Load loop data (pattern from section_27)
loops <- read.table(LOOP_FILES$late, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract anchors as GRanges (each loop has two anchors)
anchor1_gr <- GRanges(seqnames = loops$chr1, ranges = IRanges(start = loops$x1, end = loops$x2))
anchor2_gr <- GRanges(seqnames = loops$chr2, ranges = IRanges(start = loops$y1, end = loops$y2))

# Split by loop direction
gained_mask <- loops$logFC > 0 & loops$FDR < 0.05
lost_mask   <- loops$logFC < 0 & loops$FDR < 0.05

# Combine both anchors per direction, reduce to unique regions
gained_anchors <- reduce(c(anchor1_gr[gained_mask], anchor2_gr[gained_mask]))
lost_anchors   <- reduce(c(anchor1_gr[lost_mask], anchor2_gr[lost_mask]))

Alist_35c <- list(
  "Gained Loop Anchors" = gained_anchors,
  "Lost Loop Anchors"   = lost_anchors
)

Blist_35c <- list(
  "ATAC Up"       = atac_up,
  "ATAC Down"     = atac_down,
  "MeCP2 Up"      = load_chip_peaks(MECP2_FILES$up, "MeCP2 Up"),
  "MeCP2 Down"    = load_chip_peaks(MECP2_FILES$down, "MeCP2 Down"),
  "mC Hyper DMRs" = hyper_gr,
  "mC Hypo DMRs"  = hypo_gr
)
```

One `crosswisePermTest()` call: 2 x 6 = 12 pairwise tests.

**Note on column names:** The exact column names in the loop file depend on the mariner pipeline output. Verify by reading the header of `LOOP_FILES$late` — look for chr1/x1/x2/chr2/y1/y2 or seqnames1/start1/end1/seqnames2/start2/end2 patterns.

### Figures

| Figure | Description |
|--------|-------------|
| **35a** | Crosswise heatmap: ATAC direction x 6 ChIP marks (2 x 6) |
| **35b** | Crosswise heatmap: ATAC direction x 7 chromatin states (2 x 7) |
| **35c** | Crosswise heatmap: Loop anchor direction x chromatin features (2 x 6) |
| **35d** | Comparison table: all Group B+C Fisher tests vs permutation results |
| **35e** | Local z-score curve for strongest loop anchor association |

### Mapping to Existing Fisher Tests

| Sub-analysis | Original section | Original figure | Tests |
|-------------|------------------|-----------------|-------|
| 35A | 13 | 13b | 6 Fisher tests (ATAC up/down x 6 marks), BH-corrected |
| 35B | 13 | 13c | 2x7 O/E heatmap (ATAC direction x chromatin state) |
| 35C | 13 | 13d | 2 Fisher tests (loop direction x ATAC at anchors) |
| 35C | 27 | 27a | 2 Fisher tests (coordinated x loop anchor, GREAT + nearest) |
| 35C | 27 | 27b | 3 Fisher tests (hyper x lost/gained anchor vs background) |
| 35C | 27 | 27e | 2 Fisher tests (coordinated x shared anchor) |
| 35C | 31 | 31a | 2 Fisher tests (loop direction x MeCP2 overlap) |

### Runtime Estimate
~30-45 minutes total (3 crosswisePermTest calls).

---

## 5. Section 36 — Domain-Level Permutation (Compartments + Polycomb)

**File:** `scripts/viz_sections/section_36_permutation_domains.R`
**Validates Fisher tests from:** sections 29, 30

### Sub-analysis 36A: DMR Direction x A/B Compartments (validates section 29)

**Question:** Are hypermethylated DMRs enriched in A-compartment regions, accounting for the megabase-scale structure of compartment domains?

Construct compartment GRanges from HOMER data (pattern from section_29, lines 30-116):

```r
# Load HOMER diffcompartments.txt
COMPARTMENT_FILE <- file.path(BASE_DIR, "../../tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt")
comp_raw <- read.table(COMPARTMENT_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Identify PC1 columns, compute ctrl mean for A/B classification
# (Follow exact pattern from section_29_ab_compartment_methylation_mapping.R)

# Construct GRanges for compartment bins
comp_gr <- GRanges(
  seqnames = comp_raw$chr,
  ranges = IRanges(start = comp_raw$start, end = comp_raw$end),
  mean_ctrl_pc1 = comp_raw$mean_ctrl_pc1
)

# Split into A and B
a_comp <- comp_gr[comp_gr$mean_ctrl_pc1 > 0]
b_comp <- comp_gr[comp_gr$mean_ctrl_pc1 < 0]

# Shift bins (significant compartment switches, FDR < 0.05, |Diff| > 0.30)
b_to_a_shift <- comp_gr[comp_raw$FDR < 0.05 & comp_raw$Difference > 0.30]
a_to_b_shift <- comp_gr[comp_raw$FDR < 0.05 & comp_raw$Difference < -0.30]

Alist_36a <- list(
  "mC Hyper DMRs" = hyper_gr,
  "mC Hypo DMRs"  = hypo_gr,
  "hmC Hypo DMRs" = dmr_to_granges(hmc_dmr %>% dplyr::filter(significant, mod_difference < 0)),
  "hmC Hyper DMRs" = dmr_to_granges(hmc_dmr %>% dplyr::filter(significant, mod_difference > 0))
)

Blist_36a <- list(
  "A Compartment"  = a_comp,
  "B Compartment"  = b_comp,
  "B->A Shift"     = b_to_a_shift,
  "A->B Shift"     = a_to_b_shift
)
```

One `crosswisePermTest()` call: 4 x 4 = 16 pairwise tests.

**Note on HOMER column names:** Verify column names in the diffcompartments.txt file — they depend on HOMER version. Look for `chr`, `start`, `end`, and PC1 sample columns. The FDR and Difference columns come from `getDiffExpression.pl`.

### Sub-analysis 36B: DMR Direction x Polycomb Domains (validates section 30)

**Question:** Are hypermethylated DMRs enriched at Polycomb-marked regions?

```r
Blist_36b <- list(
  "H3K27me3 Peaks" = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  "Bivalent Peaks"  = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)
# Alist_36b = same Alist_36a (4 DMR direction sets)
```

One `crosswisePermTest()` call: 4 x 2 = 8 pairwise tests.

### Figures

| Figure | Description |
|--------|-------------|
| **36a** | Crosswise heatmap: DMR direction x compartment (4 x 4) |
| **36b** | Crosswise heatmap: DMR direction x Polycomb domains (4 x 2) |
| **36c** | Local z-score decay for mC Hyper x A compartment (expected broad, not focal) |
| **36d** | Comparison table: all section 29/30 Fisher tests vs permutation results |

### Mapping to Existing Fisher Tests

| Sub-analysis | Original section | Description | N tests |
|-------------|------------------|-------------|---------|
| 36A | 29 Step 2 | DMR direction x A/B compartment | 4 Fisher tests |
| 36A | 29 Step 3 | Shift direction x DMR direction (interval-level part) | 2 Fisher tests + O/E |
| 36B | 30 | Polycomb definitions x DMR direction | 12+ Fisher tests (BH-corrected) |

### Runtime Estimate
~20-30 minutes (2 crosswisePermTest calls, compartment bins are large so overlaps are fast).

---

## 6. Section 37 — Gene-Level Label-Shuffle Permutation

**File:** `scripts/viz_sections/section_37_permutation_gene_level.R`
**Validates Fisher/O/E from:** sections 12e, 15a, 15b, 15c, 19g, 20d, 27c, 29 Step 3, 31b, 33c

### Why Not regioneReloaded?

These analyses cross-tabulate gene-level binary labels: "mC Up" vs "mC Down" against "ATAC Up" vs "ATAC Down" at the gene level. The unit is a gene symbol with two categorical assignments, not a genomic interval to overlap. regioneReloaded's `crosswisePermTest` tests whether region sets A and B spatially co-occur. That is not the question here — the question is whether the *directional labels* are correlated across genes.

### Methodology: Chromosome-Stratified Label Shuffle

```r
#' Chromosome-stratified permutation test for 2x2 gene-level contingency tables
#'
#' Shuffles col_label assignments within each chromosome, preserving:
#'   - Per-chromosome marginal counts
#'   - Total count of each label
#'   - Genomic neighborhood structure (genes on same chr stay together)
#'
#' @param gene_df data.frame with columns: gene, chr, row_label, col_label
#' @param ntimes Number of permutations (default 10000)
#' @return list with observed_or, null_distribution, empirical_p, z_score, ci_95
permute_gene_2x2 <- function(gene_df, ntimes = 10000) {
  stopifnot(all(c("gene", "chr", "row_label", "col_label") %in% colnames(gene_df)))

  # Observed log2(Fisher OR)
  obs_table <- table(gene_df$row_label, gene_df$col_label)
  obs_fisher <- fisher.test(obs_table)
  obs_log2or <- log2(obs_fisher$estimate)

  # Permutation null: shuffle col_label within chromosomes
  null_log2or <- numeric(ntimes)
  for (i in seq_len(ntimes)) {
    shuffled <- gene_df %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(col_label = sample(col_label)) %>%
      dplyr::ungroup()
    perm_table <- table(shuffled$row_label, shuffled$col_label)
    perm_fisher <- fisher.test(perm_table)
    null_log2or[i] <- log2(perm_fisher$estimate)
  }

  # Empirical p-value (two-sided)
  empirical_p <- (sum(abs(null_log2or) >= abs(obs_log2or)) + 1) / (ntimes + 1)

  # Z-score
  z_score <- (obs_log2or - mean(null_log2or)) / sd(null_log2or)

  # 95% CI from null distribution
  ci_95 <- quantile(null_log2or, c(0.025, 0.975))

  list(
    observed_or = obs_fisher$estimate,
    observed_log2or = obs_log2or,
    fisher_p = obs_fisher$p.value,
    null_distribution = null_log2or,
    empirical_p = empirical_p,
    z_score = z_score,
    ci_95 = ci_95
  )
}
```

**Why shuffle within chromosomes?** Genes on the same chromosome share regulatory context (same TAD neighborhood, similar replication timing, chromosome-level methylation patterns). Free shuffling across the whole genome would overestimate degrees of freedom. Chromosome-stratified shuffling is the gene-level analogue of regioneReloaded's `per.chromosome = TRUE`.

**Why 10,000 permutations (not 5,000)?** Label shuffling is computationally cheap (~0.01s per iteration vs ~1s for genomic interval permutation). 10,000 gives better resolution in the tail for empirical p-values.

### Data Preparation Pattern

Each test requires a `gene_df` with 4 columns. The general pattern for methylation x mark direction:

```r
# Example: hmC direction x MeCP2 direction (validates section 15a)
# Start from hmc_dmr (pre-loaded), join with MeCP2 gene-level data

hmc_sig <- hmc_dmr %>% dplyr::filter(significant)

# MeCP2 gene-level data (re-computed from MECP2_FILES$annotated, same as section 15)
mecp2_gene <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t") %>%
  # ... aggregate nearest-to-TSS peak per gene ...
  dplyr::mutate(mecp2_direction = ifelse(Fold > 0, "MeCP2 Up", "MeCP2 Down"))

gene_df_15a <- hmc_sig %>%
  dplyr::inner_join(mecp2_gene, by = "gene") %>%
  dplyr::filter(mecp2_direction %in% c("MeCP2 Up", "MeCP2 Down")) %>%
  dplyr::transmute(
    gene,
    chr,
    row_label = ifelse(mod_difference < 0, "hmC Down", "hmC Up"),
    col_label = mecp2_direction
  )

result_15a <- permute_gene_2x2(gene_df_15a, ntimes = 10000)
```

### Complete Test Inventory

| Test ID | Row labels | Col labels | Validates | Data source for col_label |
|---------|-----------|------------|-----------|---------------------------|
| 37-01 | mC Up / mC Down | ATAC Up / ATAC Down | 12e | net_atac from ATAC_FILES up/down BEDs, sign of (n_up - n_down) per gene |
| 37-02 | hmC Down / hmC Up | MeCP2 Up / MeCP2 Down | 15a | MECP2_FILES$annotated, nearest-to-TSS peak Fold |
| 37-03 | hmC Down / hmC Up | ATAC Up / ATAC Down | 15b | Same as 37-01 column |
| 37-04 | hmC Down / hmC Up | K119ub Gained / K119ub Lost | 15c | K119UB_FILES ctrl/mut, net overlap (n_mut - n_ctrl) per gene |
| 37-05 | mC Up / mC Down | H3K27ac Gained / H3K27ac Lost | 19g (mC) | H3K27AC_FILES ctrl/mut, net overlap per gene |
| 37-06 | hmC Down / hmC Up | H3K27ac Gained / H3K27ac Lost | 19g (hmC) | Same as 37-05 column |
| 37-07a | mC Up / mC Down | ATAC DiffBind Up / Down | 33c (mC x ATAC) | DIFFBIND_FILES$atac, gene-level via ChIPseeker |
| 37-07b | mC Up / mC Down | K27ac DiffBind Up / Down | 33c (mC x K27ac) | DIFFBIND_FILES$k27ac, gene-level via ChIPseeker |
| 37-07c | mC Up / mC Down | K27me3 DiffBind Up / Down | 33c (mC x K27me3) | DIFFBIND_FILES$k27me3, gene-level via ChIPseeker |
| 37-07d | mC Up / mC Down | K119ub DiffBind Up / Down | 33c (mC x K119ub) | DIFFBIND_FILES$k119ub, gene-level via ChIPseeker |
| 37-08 | mC Up / mC Down | Expr Up / Expr Down | 20d | RNA-seq DEG results (external file, loaded in section 20) |
| 37-09 | K119ub Gained / K119ub Lost | Hyper / Not Hyper | 27c | Loop anchor gene subset, K119ub at anchor |
| 37-10 | B->A Shift / A->B Shift / Stable | mC Hyper / mC Hypo | 29 Step 3 | HOMER compartment shift assignment per gene |
| 37-11 | Loop Gained / Loop Lost | MeCP2 Sig Up / Other | 31b | MeCP2 status at loop-anchor genes |

**Note on 33c:** Section 33 runs 8 gene-level 2x2s (4 marks x 2 methylation perspectives: mC and hmC). Tests 37-07a through 37-07d cover the mC perspective. The hmC perspective can be added as 37-07e through 37-07h if needed (same structure, swap row labels to hmC Down/Up).

**Note on test 37-08:** The RNA-seq DEG results file path is not in `_shared_config.R`. Check section 20 for the exact path — it loads a DESeq2/edgeR results table. This test can be skipped if the file is not available locally.

### Figures

| Figure | Description |
|--------|-------------|
| **37a** | Forest plot of all permutation z-scores, color-coded by original section. 95% CI bars from null distribution. Dashed vertical line at z=0. Stars for empirical p < 0.05. |
| **37b** | Null distribution histograms for 4 strongest effects (2x2 panel grid). Grey histogram = null log2(OR) distribution, red vertical line = observed log2(OR), shaded tail = empirical p-value area. |
| **37c** | Comparison table: test ID, original section/figure, Fisher OR, Fisher p, permutation z, empirical p, concordance label. |
| **37d** | Scatter: observed log2(OR) (x-axis) vs permutation z-score (y-axis). Points should cluster along a positive diagonal if Fisher and permutation agree. Color by concordance status. |

### Runtime Estimate
~15-25 minutes (10,000 permutations x ~15 tests, each iteration is a fast table shuffle).

---

## 7. Complete Mapping Table

Every Tier 1 + Tier 2 Fisher test mapped to its permutation counterpart.

### Tier 1: Interval-Level Tests

| Section | Figure | Test description | Unit | Perm section | Perm cell/sub-analysis |
|---------|--------|------------------|------|-------------|----------------------|
| 12 | 12a | DMR hyper/hypo x ATAC up/down | DMR interval | 34 | crosswise rows 1-2 x cols 1-2 |
| 13 | 13b | ATAC up/down x 6 ChIP marks (BH) | ATAC peak | 35 | 35A crosswise (2 x 6) |
| 13 | 13c | ATAC up/down x 7 chromatin states (O/E) | ATAC peak | 35 | 35B crosswise (2 x 7) |
| 13 | 13d | Loop direction x ATAC at anchors | Loop | 35 | 35C crosswise (ATAC cols) |
| 14 | 14a | DMR hyper/hypo x K119ub ctrl/mut | DMR interval | 34 | crosswise rows 1-2 x cols 3-4 |
| 14 | 14b | DMR hyper/hypo x K119ub gained/lost | DMR interval | 34 | crosswise rows 1-2 x cols 5-6 |
| 19 | 19f | DMR hyper/hypo x H3K27ac ctrl/mut | DMR interval | 34 | crosswise rows 1-2 x cols 7-8 |
| 27 | 27a | Coordinated x loop anchor (GREAT/nearest) | Gene | 35 | 35C (DMR overlap with anchors) |
| 27 | 27b | Hypermethylated x lost/gained/bg anchor | Gene | 35 | 35C (Hyper DMR x anchor direction) |
| 27 | 27e | Coordinated x shared anchor | Gene | 35 | 35C (shared anchor rows) |
| 29 | Step 2 | DMR direction x A/B compartment | Gene | 36 | 36A crosswise (4 x 4) |
| 30 | 30b/d | Polycomb x DMR direction (BH) | Gene | 36 | 36B crosswise (4 x 2) |
| 31 | 31a | Loop direction x MeCP2 overlap | Loop | 35 | 35C (MeCP2 cols) |

### Tier 2: Gene-Level O/E Tests

| Section | Figure | Test description | Unit | Perm section | Perm test ID |
|---------|--------|------------------|------|-------------|-------------|
| 12 | 12e | mC direction x ATAC direction (O/E) | Gene | 37 | 37-01 |
| 15 | 15a | hmC direction x MeCP2 direction (O/E) | Gene | 37 | 37-02 |
| 15 | 15b | hmC direction x ATAC direction (O/E) | Gene | 37 | 37-03 |
| 15 | 15c | hmC direction x K119ub direction (O/E) | Gene | 37 | 37-04 |
| 19 | 19g | mC direction x H3K27ac direction (O/E) | Gene | 37 | 37-05 |
| 19 | 19g | hmC direction x H3K27ac direction (O/E) | Gene | 37 | 37-06 |
| 20 | 20d | mC direction x expression direction (O/E) | Gene | 37 | 37-08 |
| 27 | 27c | K119ub x methylation at anchors (O/E) | Gene | 37 | 37-09 |
| 29 | Step 3 | Compartment shift x mC direction (O/E) | Gene | 37 | 37-10 |
| 31 | 31b | Loop direction x MeCP2 direction (O/E) | Gene | 37 | 37-11 |
| 33 | 33c | 4-mark x mC direction (O/E, 4 tests) | Gene | 37 | 37-07a to 37-07d |

---

## 8. Interpretation Guide

### Z-Score Thresholds

| |z-score| | Approximate p | Interpretation |
|-----------|---------------|----------------|
| > 1.96 | < 0.05 | Significant enrichment/depletion |
| > 2.58 | < 0.01 | Strong enrichment/depletion |
| > 3.29 | < 0.001 | Very strong enrichment/depletion |

regioneReloaded's **normalized z-scores** (`nZS = ZS / sqrt(n)`) are directly comparable across tests with different region set sizes. Use normalized z-scores when comparing across cells in the crosswise matrix.

### Concordance Categories

After running permutation tests, classify each original Fisher result:

| Category | Criteria | Meaning |
|----------|----------|---------|
| **Confirmed** | Fisher p < 0.05 AND permutation p < 0.05, same direction | Non-independence did not meaningfully inflate the Fisher result |
| **Weakened** | Fisher p < 0.05 but permutation p >= 0.05 | Spatial clustering inflated the Fisher result; the association may not be robust |
| **Strengthened** | Fisher p >= 0.05 but permutation p < 0.05 | Rare; possible when Fisher had borderline significance |
| **Concordant NS** | Both p >= 0.05 | No association detected by either method |

### Expected Outcomes

Most Fisher tests in this pipeline have large effect sizes (OR > 2) with many hundreds of observations. Permutation tests will likely **confirm** the majority. The tests most at risk of **weakening** are:
- Moderate O/E enrichments (1.3-1.8) with marginal Fisher p-values
- Tests with small N (e.g., 31b loop direction x MeCP2 with ~300 genes)
- Domain-level tests (sections 29, 30) where the effective sample size after accounting for domain structure is much smaller than the nominal gene count

### Reporting

For each test, report both:
1. The original Fisher OR and p-value (for comparability with existing literature)
2. The permutation z-score and empirical p-value (for reviewer defense)

If a result weakens: report both p-values honestly. Discuss the non-independence explanation. Do not retroactively remove the Fisher result or the original figure.

---

## 9. Runtime & SLURM

### Runtime Estimates

| Section | Method | N permutations | Estimated time (8 cores) |
|---------|--------|----------------|--------------------------|
| 34 | regioneReloaded crosswisePermTest | 5,000 | 20-30 min |
| 35 | regioneReloaded crosswisePermTest (x3) | 5,000 each | 30-45 min |
| 36 | regioneReloaded crosswisePermTest (x2) | 5,000 each | 20-30 min |
| 37 | Custom label shuffle | 10,000 per test x ~15 tests | 15-25 min |
| **Total** | | | **~2-3 hours** |

### SLURM Template

```bash
#!/bin/bash
#SBATCH --job-name=perm_tests
#SBATCH --partition=compute
#SBATCH --account=csd940
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/permutation_%j.out

module load cpu/0.17.3b
module load anaconda3
conda activate mariner_env

cd /expanse/lustre/projects/csd940/zalibhai/biomodal/downstream/

for section in scripts/viz_sections/section_3{4,5,6,7}_*.R; do
  echo "$(date): Running ${section}..."
  Rscript "${section}"
  echo "$(date): Finished ${section}"
  echo ""
done
```

### RDS Caching

Each section should save its raw permutation results as RDS:

```r
# At end of section_34:
saveRDS(cw_34, file.path(TABLES_DIR, "permutation_34_dmr_marks.rds"))

# At start of section_34, check for cache:
cache_path <- file.path(TABLES_DIR, "permutation_34_dmr_marks.rds")
if (file.exists(cache_path)) {
  cat("Loading cached permutation results...\n")
  cw_34 <- readRDS(cache_path)
} else {
  # Run crosswisePermTest...
}
```

This allows re-running the visualization without repeating the expensive permutation step.

### Local Execution

Section 37 (gene-level label shuffle) is fast enough to run locally. Sections 34-36 can also run locally but will take longer without HPC parallelism. Reduce `PERM_NTIMES` to 100 for testing, then scale to 5,000 for final results.

---

## 10. Implementation Checklist

1. [ ] Install `regioneR` and `regioneReloaded` on HPC (`mariner_env` conda environment)
2. [ ] Verify `BSgenome.Mmusculus.UCSC.mm10` is available (needed for genome object)
3. [ ] Implement `section_34_permutation_dmr_chromatin_marks.R`
   - Test with `ntimes = 100` first, verify output structure
   - Scale to `ntimes = 5000` for final run
4. [ ] Implement `section_35_permutation_atac_loops.R`
   - Verify loop file column names from `LOOP_FILES$late` header
   - 3 sub-analyses (35A, 35B, 35C)
5. [ ] Implement `section_36_permutation_domains.R`
   - Verify HOMER compartment file column names
   - 2 sub-analyses (36A, 36B)
6. [ ] Implement `section_37_permutation_gene_level.R`
   - Define `permute_gene_2x2()` helper
   - Run all ~15 tests
   - Verify RNA-seq DEG file availability for test 37-08
7. [ ] Run all 4 sections on HPC via SLURM
8. [ ] Review all comparison tables for concordance classification
9. [ ] Verify cached RDS files saved to `plots/visualizations/tables/`
10. [ ] Update `CLAUDE.md` with new section descriptions (34-37)

---

## 11. Output Directory Structure

```
plots/visualizations/
  34_permutation_dmr_chromatin/
    34a_crosswise_dmr_x_marks/           # PDF + SVG + JPG (via save_multiformat_ggplot)
    34b_fisher_vs_permutation_table/
    34c_local_zscore_strongest/
  35_permutation_atac_loops/
    35a_crosswise_atac_x_chip/
    35b_crosswise_atac_x_chromatin_state/
    35c_crosswise_loops_x_features/
    35d_fisher_vs_permutation_table/
    35e_local_zscore_loop/
  36_permutation_domains/
    36a_crosswise_dmr_x_compartment/
    36b_crosswise_dmr_x_polycomb/
    36c_local_zscore_compartment/
    36d_fisher_vs_permutation_table/
  37_permutation_gene_level/
    37a_zscore_forest_plot/
    37b_null_distribution_top4/
    37c_fisher_vs_permutation_table/
    37d_or_vs_zscore_scatter/

plots/visualizations/tables/
  permutation_34_dmr_marks.rds
  permutation_34_comparison.tsv
  permutation_35_atac_loops.rds
  permutation_35_comparison.tsv
  permutation_36_domains.rds
  permutation_36_comparison.tsv
  permutation_37_gene_level.rds
  permutation_37_comparison.tsv
```

Each comparison TSV has columns: `test_id`, `original_section`, `original_figure`, `description`, `fisher_or`, `fisher_p`, `perm_z`, `perm_p`, `concordance`.
