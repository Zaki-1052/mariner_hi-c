# Hi-C Differential Loop Analysis Results

## Overview

This directory contains results from differential chromatin loop analysis comparing **BAP1-KO mutant vs wildtype** in mouse cerebellum (mm10 genome). The analysis uses biological replicates (n=3 per condition) across two developmental timepoints (Early/P12 and Late/Adult).

**Pipeline:** mariner + edgeR quasi-likelihood GLM with TMM normalization
**Resolutions analyzed:** 5kb, 10kb, 25kb (merged with 10kb tolerance)
**Primary statistical thresholds:** FDR < 0.05, |logFC| > 0.3

---

## Central Biological Finding

**BAP1-KO induces coordinated "loop switching" at Polycomb-marked chromatin domains:**

1. **Long-range loops are lost** (median 1.15 Mb) at specific genomic anchors
2. **Short-range loops are gained** (median 340 kb) at the **same** anchors
3. This switching is **2.5x stronger at H3K27me3-marked (Polycomb) sites** than non-Polycomb sites
4. The pattern is consistent with BAP1's role in H2AK119 deubiquitination affecting Polycomb-mediated long-range chromatin interactions

---

## Directory Structure

```
output/
├── shared_anchor_analysis/       # Core loop switching analysis (Section 1)
│   ├── early/
│   └── late/
│       ├── tables/               # Shared anchors, loops, distances
│       ├── plots/                # Multi-format (PDF/SVG/JPG)
│       ├── apa_subsets/          # RDS files for APA heatmaps
│       └── polycomb_specific/    # Polycomb-enriched subset
│
├── loops_mark_filtered/          # Per-histone-mark CDF/density (Section 4)
│   ├── early/{H3K27ac,H3K27me3,H3K4me1,H3K4me3}/
│   └── late/{H3K27ac,H3K27me3,H3K4me1,H3K4me3,CTCF,Bivalent}/
│
├── chip_distance/                # ChIP-seq signal vs loop distance (Section 3b)
│   ├── early/
│   └── late/
│
├── ctcf_stripe_crossref/         # CTCF-stripe independence test (Section 5c)
│   ├── early/
│   └── late/
│
├── diff_chip_polycomb_enrichment/ # Differential K27me3/H2AK119ub (Section 3f)
│
├── deg_loop_violin/              # RNA-seq integration (Section 2)
│   ├── early/
│   └── late/
│
├── loops_k27me3_filtered/        # K27me3-anchored loop subsets
├── loops_k27me3_global/          # Global K27me3 distance analysis
├── loops_ep_filtered/            # Enhancer-Promoter loops
├── loops_visualization_extended/ # GO/KEGG enrichment, distance plots
└── k27me3_loops_bedpe/           # BEDPE files for Juicebox
```

---

## Analysis Results by Section

### 1. Shared Anchor / Loop Switching Analysis

**Script:** `scripts/shared_anchor_analysis.R`, `scripts/polycomb_shared_anchor_analysis.R`
**Output:** `output/shared_anchor_analysis/{early,late}/`

#### Key Statistics (Late Timepoint)

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Shared anchor regions | **212** | Genomic sites with both lost AND gained loops |
| Lost loops at shared anchors | 283 | 24% of all lost loops |
| Gained loops at shared anchors | 321 | 19% of all gained loops |
| Chromatin state Chi-square p | **3.28e-31** | Anchor chromatin states are non-random |
| Paired Wilcoxon p (lost > gained) | **1.17e-20** | Lost loops are longer than gained |
| **Median lost distance** | **1.15 Mb** | Lost loops are long-range |
| **Median gained distance** | **340 kb** | Gained loops are short-range |
| Distance ratio (lost/gained) | **3.4x** | |
| Cohen's d effect size | 0.31 | Small-medium effect |

#### Polycomb-Specific Results

| Metric | Polycomb Anchors | Non-Polycomb Anchors | Ratio |
|--------|------------------|----------------------|-------|
| N loops | 312 (51.7%) | 292 (48.3%) | - |
| Median lost distance | 1.73 Mb | 685 kb | 2.5x |
| Median gained distance | 370 kb | 360 kb | ~1x |
| **Distance ratio (lost/gained)** | **4.68x** | **1.90x** | **2.46x** |
| Wilcoxon p | 9.63e-36 | 2.15e-07 | - |
| KS test D statistic | 0.713 | 0.268 | - |

**Interaction test:** Two-way ANOVA p < 0.0001 confirms Polycomb-specific amplification of the switching effect.

**Polycomb enrichment:** OR = 2.15 (95% CI: 1.79-2.59), Fisher p = 1.68e-16
Shared anchors are 2x more likely to be Polycomb-marked (51.7% vs 33.2% background).

#### Both-Anchors-Polycomb Subset

When **both** loop anchors are Polycomb-marked, the effect is even stronger:
- Distance ratio: **8.61x** (lost median 3.1 Mb, gained median 360 kb)
- Wilcoxon p = 1.44e-20
- KS D = 0.816

#### Output Files

| File | Description |
|------|-------------|
| `tables/shared_anchors.tsv` | 212 anchor coordinates with chromatin states |
| `tables/shared_anchor_loops.tsv` | 604 loops at shared anchors |
| `tables/paired_distance_stats.tsv` | Per-anchor lost/gained distance comparison |
| `polycomb_specific/tables/polycomb_shared_loops.tsv` | 312 Polycomb-anchored loops |
| `polycomb_specific/tables/interaction_test.tsv` | Statistical interaction results |
| `apa_subsets/shared_lost_longrange.rds` | 141 loops for APA |
| `apa_subsets/shared_gained_shortrange.rds` | 157 loops for APA |

---

### 2. RNA-seq Integration

**Script:** `scripts/deg_loop_anchor_violin.R`
**Output:** `output/deg_loop_violin/{early,late}/`

Maps differentially expressed genes (DEGs) to loop anchors using GREAT-style regulatory domains.

#### Results (Late Timepoint)

- DEGs at shared anchors: 42 genes
- Median log2FC of associated DEGs: 0.050
- One-sample Wilcoxon p: 0.629 (not significant)

**Interpretation:** While there are DEGs at shared anchors, the expression changes are modest and bidirectional. This suggests:
1. Loop switching may affect gene regulation through mechanisms beyond simple on/off changes
2. Post-transcriptional effects or developmental timing may mask transcriptional changes
3. Some loops may be structural rather than regulatory

#### Output Files

| File | Description |
|------|-------------|
| `deg_anchor_genes.tsv` | DEGs mapped to loop anchors |
| `deg_loop_distance_stratified_genes.tsv` | Stratified by loop distance |
| `plots/deg_loop_anchor_violin.pdf` | Main violin plot |
| `plots/deg_loop_violin_longrange_lost_vs_shortrange_gained.pdf` | Key comparison |

---

### 3. Polycomb Loop Analysis

#### 3a-b. K27me3-Filtered Distance Distributions

**Scripts:** `scripts/loop_distance_k27me3_filtered.R`, `scripts/loop_distance_k27me3_global.R`
**Output:** `output/loops_k27me3_filtered/`, `output/loops_k27me3_global/`

CDF and density plots for loops filtered by H3K27me3 overlap at anchors.

**Findings:**
- K27me3-anchored loops show stronger distance shift than non-K27me3 loops
- Bivalent (K4me3+K27me3) loops show intermediate behavior

#### 3c. ChIP-seq Signal vs Loop Distance

**Script:** `scripts/chip_distance_analysis.R`
**Output:** `output/chip_distance/{early,late}/`

**Key Finding:**
- H3K27ac signal **decreases** with increasing loop distance
- H3K27me3 signal **increases** with increasing loop distance
- This supports the model that long-range loops are Polycomb-mediated

#### 3f. Differential H3K27me3/H2AK119ub at Loop Categories

**Script:** `scripts/diff_chip_polycomb_enrichment.R`
**Output:** `output/diff_chip_polycomb_enrichment/`

Uses diffbind differential peaks (summit=400bp) to test whether K27me3/H2AK119ub changes correlate with loop changes.

##### Significant Enrichments (All Loops, FDR < 0.05)

| Comparison | Odds Ratio | FDR | Interpretation |
|------------|------------|-----|----------------|
| **K27me3_down at short_range_gained** | 8.78 | 6.86e-102 | Gained loops lose K27me3 |
| **H2AK119ub_up at long_range_lost vs short_range_gained** | 4.87 | 1.86e-90 | Lost loops gain ubiquitin |
| H2AK119ub_up at long_range_lost | 2.18 | 3.42e-41 | |
| H2AK119ub_down at short_range_gained | 3.20 | 4.46e-41 | |
| K27me3_up at long_range_lost | 3.18 | 3.58e-28 | Lost loops gain K27me3 |

##### Polycomb Shared Anchors

| Comparison | Odds Ratio | FDR |
|------------|------------|-----|
| H2AK119ub_down at long_range_lost vs short_range_gained | 3.92 | 2.42e-03 |

**Biological Interpretation:**
- Lost long-range loops occur at sites gaining H3K27me3 and H2AK119ub
- Gained short-range loops occur at sites losing H3K27me3 and H2AK119ub
- Consistent with BAP1's role as H2AK119 deubiquitinase - loss of BAP1 leads to ubiquitin accumulation which may destabilize long-range Polycomb interactions

**Caveat:** 400bp summit width underestimates broader H3K27me3/H2AK119ub domains; bigWig aggregate analysis recommended for complete picture.

---

### 4. Per-Mark CDF/Density Subsetting

**Script:** `scripts/loop_distance_mark_filtered.R`
**Output:** `output/loops_mark_filtered/{early,late}/{mark}/`

Distance distributions for loops filtered by ChIP-seq mark overlap at anchors.

#### Marks Analyzed

| Mark | One-Anchor Filter | Both-Anchors Filter | Key Finding |
|------|-------------------|---------------------|-------------|
| H3K27ac | Active enhancer | Super-enhancer | 92% of super-enhancer loops lost |
| H3K27me3 | Polycomb | Polycomb-Polycomb | Strong distance shift |
| H3K4me1 | Poised enhancer | - | Moderate shift |
| H3K4me3 | Active promoter | - | Weak shift |
| CTCF | Architectural | Both-CTCF | See Section 5 |
| Bivalent | K4me3+K27me3 | - | Intermediate |

#### Output Files Per Mark

```
output/loops_mark_filtered/late/H3K27ac/
├── 01_cdf_one_anchor.{pdf,svg,jpg}
├── 01_cdf_both_anchors.{pdf,svg,jpg}
├── 03_density_one_anchor.{pdf,svg,jpg}
├── 03_density_both_anchors.{pdf,svg,jpg}
└── filter_summary.txt
```

---

### 5. CTCF / Loop Extrusion Analysis

**Script:** `scripts/ctcf_stripe_crossref.R`
**Output:** `output/ctcf_stripe_crossref/{early,late}/`

Tests whether CTCF-anchored loop loss correlates with stripe defects (cohesin extrusion dysfunction).

#### Results (Late Timepoint)

| Metric | Value |
|--------|-------|
| Total lost loops | 1,187 |
| Lost loops with CTCF at anchor | 890 |
| CTCF lost loops at lost stripes | 17/890 (1.9%) |
| Non-CTCF lost loops at lost stripes | 5/297 (1.7%) |
| **Fisher's exact p-value** | **1.00** |
| Odds ratio | 1.14 (95% CI: 0.40-3.98) |
| Permutation test p-value | 0.18 |

**Interpretation:** No significant enrichment detected. **CTCF loop loss and stripe defects appear to be independent phenomena**, not coordinated by common cohesin dysfunction.

This suggests BAP1-KO affects chromatin loops through Polycomb mechanisms rather than loop extrusion machinery.

---

### 6. TAD Boundary Integration

**Script:** `tads/scripts/boundary_loop_crossref.R`
**Output:** `tads/output/boundary_loop_analysis/`

Cross-references differential loops with differential TAD boundaries from TADCompare.

**Status:** Basic analysis complete; permutation testing needs redo with regioneR/regioneReloaded for proper null distribution.

---

## Statistical Framework

### Differential Analysis

1. **Loop Calling:** HiCCUPS at 5kb, 10kb, 25kb resolutions
2. **Merging:** mariner `mergePairs()` with 10kb tolerance across 6 replicates
3. **Count Extraction:** 5x5 Hi-C matrices via strawr, summed to single value
4. **Normalization:** TMM (Trimmed Mean of M-values)
5. **Model:** Quasi-likelihood GLM with robust dispersion estimation
6. **Testing:** `glmQLFTest` with Benjamini-Hochberg FDR correction

### Quality Control Thresholds

- Within-group replicate correlation: > 0.90
- Library size ratio: < 2-fold
- BCV (biological coefficient of variation): 0.2-0.5

### Enrichment Testing

- Fisher's exact test for 2x2 contingency tables
- Permutation tests (n=1000) for spatial overlap
- Two-way ANOVA for interaction effects
- Kolmogorov-Smirnov for distribution comparisons

---

## Key Data Files

### Loop Data

| Timepoint | File | Loops |
|-----------|------|-------|
| Late | `25042-late_outputs/merged_loops/characterized_loops.tsv` | ~39,000 |
| Early | `250831-early_outputs/merged_loops/characterized_loops.tsv` | Similar |

Columns include: coordinates, logFC, FDR, direction, distance, 5 ChIP-seq overlaps, 7-category anchor type, gene annotations.

### Extended Annotation

| Timepoint | File |
|-----------|------|
| Late | `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv` |
| Early | `peaks/loop_annotation_extended/early/extended_characterized_loops.tsv` |

7 chromatin state categories:
1. Active_Promoter (H3K4me3+, not K27me3, ≤2kb from TSS)
2. Repressed_Promoter (H3K27me3+, not K27ac, ≤2kb from TSS)
3. Bivalent_Promoter (K4me3+K27me3 overlap)
4. Polycomb (H3K27me3+, >2kb from TSS)
5. Active_Enhancer (H3K27ac+, >2kb from TSS)
6. Poised_Enhancer (H3K4me1+, not K27ac/K27me3, >2kb)
7. Other

---

## Outstanding Analyses

### Partially Complete

- [ ] **3d.** APA heatmaps for Polycomb-anchored loops (subset RDS files ready)
- [ ] **3g.** Aggregate ChIP-seq signal heatmaps at loop anchors (bigWig-based)
- [ ] **6d.** Proper permutation analysis for TAD-loop enrichment (regioneR/regioneReloaded)

### Not Started

- [ ] **7.** ABC Model / Enhancer-Gene Linkage
- [ ] **8.** H2AK119ub integration with loop switching

---

## Reproducibility

All analyses were run with:
- R 4.3+, Bioconductor 3.18+
- mariner, edgeR, InteractionSet, GenomicRanges, strawr
- Python 3.8+ (pandas, numpy)
- HDF5Array for out-of-memory matrices

Scripts are parameterized by timepoint (early/late) and designed for SLURM HPC execution.

---

## Summary Interpretation

The data strongly support a model where **BAP1 loss disrupts Polycomb-mediated long-range chromatin interactions**:

1. **Shared anchors** (212 sites) show coordinated loss of long-range and gain of short-range loops
2. This effect is **2.5x stronger at H3K27me3-marked anchors** (Polycomb enrichment OR=2.15)
3. Lost loops occur at sites **gaining H2AK119ub** (consistent with BAP1's deubiquitinase role)
4. Gained short-range loops occur at sites **losing H3K27me3**
5. The phenomenon is **independent of CTCF/cohesin loop extrusion**
6. Gene expression changes are modest, suggesting structural or indirect regulatory effects

The distance shift from 1.15 Mb to 340 kb (3.4x) represents a major reorganization of chromatin architecture, with Polycomb domains driving the transition from distal to proximal interactions.

---

*Generated: 2026-02-01*
*Repository: mariner_hi-c*
*Contact: See CLAUDE.md for project documentation*
