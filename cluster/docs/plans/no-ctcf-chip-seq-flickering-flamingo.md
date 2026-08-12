# Plan: Section 47 — CTCF Anchor Methylation Overlay

## Context

Jesse Dixon connected the lost CTCF-CTCF loop length asymmetry in this BAP1-KO project to the Flavahan/Bernstein IDH glioma work (Nature 2015), where DNA hypermethylation at CpG-rich CTCF sites causes CTCF binding loss and insulation failure. The existing `scripts/cpg_ctcf_analysis.R` (run 2026-04-12) confirmed that lost CTCF loop anchors are enriched for CpG islands (OR=1.53, p=2.9e-6), but tested only annotation/GC content — not actual methylation levels. The biomodal DUET data (run-5, 8 samples) provides region-level differential methylation at CpG islands, shores, shelves, and promoters, which has never been overlaid onto CTCF loop anchors. This section fills that gap.

Section 27 already does methylation × loop anchor integration, but at the **gene level** (GREAT domains, nearest gene) for **all anchors**. Section 47 is distinct: it works at the **genomic region coordinate level** (direct GRanges overlap) restricted to **CTCF-anchored loops**.

## Implementation

Create: `biomodal/downstream/scripts/viz_sections/section_47_ctcf_anchor_methylation_overlay.R`

### Structure (follows section pattern exactly)

```
source("scripts/viz_sections/_shared_config.R")
# Banner, CONFIGURATION, INPUT VALIDATION, DATA LOADING
# Sub-analyses 47a–47e
# SUMMARY block
```

### Data Loading

| Source | Variable | Notes |
|--------|----------|-------|
| `LOOP_FILES$late` | `loops` | 2910 loops, read.table |
| `DATA_PATHS$cpg_islands_mc` | `cpgi_mc` | via `load_dmr_bed()`, 8910 rows |
| `DATA_PATHS$cpg_islands_hmc` | `cpgi_hmc` | via `load_dmr_bed()`, 8910 rows |
| `DATA_PATHS$cpg_shores_mc/hmc` | `shores_mc/hmc` | 32,581 rows |
| `DATA_PATHS$cpg_shelves_mc/hmc` | `shelves_mc/hmc` | 29,094 rows |
| `DATA_PATHS$promoters_mc/hmc` | `promoters_mc/hmc` | 20,417 rows, has gene name |

All coordinate systems will be `chr`-prefixed after `load_dmr_bed()`. No additional prefix handling needed.

### Local Helpers (defined in-script)

- `extract_ctcf_anchors_local(loops_df, ctcf_col_prefix)` — replicate logic from `cpg_ctcf_analysis.R` lines 132-163: extract CTCF+ anchors, deduplicate by `(chr, start, end, direction)`, carry `loop_distance` and `n_loops`
- `fmt_p(p)` — standard p-value formatter
- `region_dmr_to_granges(dmr_df, extra_mcols)` — extends `dmr_to_granges()` to carry mod_difference, dmr_qvalue, significant, direction as mcols

### Sub-Analyses

#### 47a: CpG Island Methylation at CTCF Loop Anchors (Core Test)

1. Extract unique CTCF ChIP+ anchors: lost (~1085) and gained (~2192)
2. Convert anchors + CpG island DMRs to GRanges
3. `findOverlaps(cpgi_gr, lost_anchor_gr)` and `findOverlaps(cpgi_gr, gained_anchor_gr)`
4. Classify each CpG island: at_lost / at_gained / at_both / at_neither (background)
5. **Fisher's exact**: % mC-hypermethylated CpG islands at lost anchors vs gained anchors
6. **Fisher's exact**: % hmC-hypomethylated CpG islands at lost vs gained
7. **Wilcoxon**: mc_mod_difference distribution at lost vs gained vs background
8. **Sensitivity**: repeat Fisher's with motif-based CTCF definition
9. Plots: grouped bar (% hyper by group), violin mC, violin hmC
10. Tables: per-CpG-island overlay, Fisher's results

#### 47b: Multi-Region Comparison

1. For each of 4 region types (CpG islands, shores, shelves, promoters) × 2 modalities (mC, hmC):
   - findOverlaps with lost/gained anchor GRanges
   - Compute % hypermethylated at lost vs gained
   - Fisher's exact test + Wilcoxon
2. BH-correct across all 8 tests
3. Plots: forest plot of ORs, heatmap of % hyper
4. Tables: multi-region Fisher's + Wilcoxon results

#### 47c: Coordinated mC-up / hmC-down at CTCF Anchor CpG Islands

1. Inner join `cpgi_mc` and `cpgi_hmc` on `(chr, start, end)` — same coordinate set
2. Define `coordinated = mc_mod_difference > 0 & hmc_mod_difference < 0`
3. Classify by anchor group (lost/gained/both/neither)
4. **Fisher's**: coordinated enrichment at lost anchors vs gained, vs background
5. Plots: % coordinated bar chart, mC-vs-hmC scatter colored by anchor group, combined effect violin
6. Tables: per-CpG-island coordinated flags, Fisher's enrichment results

#### 47d: Distance-Stratified Analysis

1. Assign each anchor a distance bin: `<200kb`, `200-500kb`, `500kb-1Mb`, `>1Mb`
2. Within each bin: Fisher's on % mC-hyper at lost vs gained CpG islands
3. **Cochran-Mantel-Haenszel test** (`mantelhaen.test()`) across all strata — the definitive distance-controlled test
4. Plots: per-bin forest plot + CMH overall, faceted violin by distance bin
5. Tables: distance-stratified Fisher's + CMH results

#### 47e: Methylation Effect Size vs Loop logFC Correlation

1. For each CTCF-CTCF loop: find CpG islands at either anchor
2. Compute mean `mc_mod_difference` and `hmc_mod_difference` per loop
3. **Spearman correlation**: mean_mc_diff vs logFC (expect: loops with more hypermethylated anchors have more negative logFC)
4. **Partial correlation**: controlling for log(loop_distance)
5. Plots: scatter with loess (mC and hmC), partial regression residual plot
6. Tables: per-loop methylation summary, correlation statistics

### Output Summary

- **9 tables** in `plots/visualizations/tables/` (prefixed `47a_`–`47e_`)
- **13 figure panels** in `plots/visualizations/47_ctcf_anchor_methylation_overlay/` (PDF+SVG+JPG each)

### Key Design Decisions

- **Aggregation**: One row per CpG island (not per anchor). If a CpG island overlaps both a lost and gained anchor, flag as "both" and exclude from primary Fisher's tests (sensitivity).
- **Background**: CpG islands NOT overlapping any CTCF loop anchor (~8500 of 8910 total)
- **No extra libraries needed**: GenomicRanges/IRanges/stats already loaded by `_shared_config.R`
- **CTCF definition**: ChIP-seq primary (column `anchor*_CTCF_overlap`), motif sensitivity in 47a only

### Critical Files

- **Pattern to follow**: `biomodal/downstream/scripts/viz_sections/section_27_methylation_hic_loop_anchor_integration.R`
- **Anchor extraction logic**: `scripts/cpg_ctcf_analysis.R` lines 132-163
- **Shared config**: `biomodal/downstream/scripts/viz_sections/_shared_config.R`
- **Loop data**: `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv`
- **Output utility**: `scripts/utils/multi_format_output.R` (already sourced by shared config)

### Verification

1. Run from `biomodal/downstream/`: `Rscript scripts/viz_sections/section_47_ctcf_anchor_methylation_overlay.R`
2. Check output tables exist in `plots/visualizations/tables/47*.tsv`
3. Check figure directories in `plots/visualizations/47_ctcf_anchor_methylation_overlay/`
4. Sanity check: the 47a Fisher's OR direction should be concordant with `cpg_ctcf_analysis.R` CpG island enrichment (OR>1 for CpG islands at lost anchors), but the methylation test adds whether those CpG islands are actually hypermethylated
5. The 47d CMH test is the definitive result — if significant after distance stratification, the methylation effect is independent of the loop length confound
