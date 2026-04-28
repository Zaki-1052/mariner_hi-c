# Plan: Section 48 — CpG Island Ubiquitination & De Novo Methylation Analysis

## Context

BAP1 is a deubiquitinase removing H2AK119ub1 (PRC1 mark). In BAP1-KO, H2AK119ub1 accumulates, potentially driving chromatin compaction that blocks TET access (mC up / hmC down) or recruits DNMTs for de novo CpG island methylation. This script answers two questions:

1. **Is H2AK119ub1 enriched at differentially methylated CpG islands?**
2. **Are hyper-methylated CpG islands undergoing de novo methylation, or amplifying pre-existing methylation?**

### Critical Data Reality (verified from actual data)

| Metric | Count |
|--------|-------|
| CpG islands tested | 8,910 |
| Significant mC DMRs | 442 (122 hyper + 320 hypo) |
| Significant hmC DMRs | 172 (26 hyper + 146 hypo) |
| Hyper islands with ctrl mC < 0.10 | 12 (10% of 122) |
| Hyper islands with ctrl mC >= 0.20 | 104 (85% of 122) |
| Mean baseline mC of hyper islands | 0.57 |

**Key biological finding:** The vast majority of hypermethylated CpG islands already had substantial mC in control — this is **amplification of existing methylation**, not classical de novo methylation. Only ~12 islands qualify as true de novo. The script must honestly report this rather than forcing a de novo narrative.

## Implementation

### New file: `downstream/scripts/viz_sections/section_48_cpg_island_ubiquitination.R`

Follows the established section pattern: source `_shared_config.R`, load extra packages, section banner, data loading, figures 48a-48m, summary table, console summary.

**Output directory:** `plots/visualizations/48_cpg_island_ubiquitination/`

### Data Sources

| Data | Source | Notes |
|------|--------|-------|
| CpG island universe (mC + hmC) | `modality/exports/cpg_islands/cpg_islands_CG_run-5_mc_hmc.tsv` | 8,910 rows, 22 cols; **bare chr numbers** — must add `chr` prefix |
| K119ub ctrl/mut peaks | `K119UB_FILES$ctrl` / `$mut` from shared config | ~20K peaks each, `chr` prefix |
| K119ub gained/lost peaks | `peaks/new/H2AK119ub_up.bed` / `_down.bed` | 6,164 gained / 1,250 lost, `chr` prefix |
| K119ub DiffBind | `DIFFBIND_FILES$k119ub` from shared config | Quantitative Fold/FDR per peak |
| ChIP peaks (H3K27me3, H3K27ac, H3K4me3) | `CHIP_PEAK_FILES` from shared config | For chromatin context |
| K119ub gene-level signal | `downstream/data/k119ub_gene_signal.tsv` | Per-gene log2fc; for promoter-proximal island subset only |

### Coordinate Handling

CpG island TSV uses bare chromosomes (`1`, `2`, ..., `X`). All ChIP/K119ub BEDs use `chr` prefix. The script adds `chr` prefix immediately on load and calls `seqlevelsStyle(gr) <- "UCSC"` on all GRanges.

### Figures (13 panels across 5 parts)

#### Part 1: K119ub Overlap at CpG Islands

**48a — K119ub enrichment bar chart with Fisher's test**
- For each K119ub peak set (ctrl, mut, gained, lost): fraction of sig mC DMR islands overlapping vs background (non-DMR)
- Split sig DMRs into hyper vs hypo subgroups
- Fisher's exact test for each comparison; annotate OR + p-value
- Uses `countOverlaps()` between `cpgi_gr` and each K119ub GRanges

**48b — K119ub peak count violin by DMR category**
- Count of overlapping K119ub mut peaks per island
- Groups: mC Hyper DMR / mC Hypo DMR / Non-significant
- Wilcoxon test: hyper vs non-sig, hypo vs non-sig
- Optional panel: for promoter-proximal islands (<=2kb from TSS), join to `k119ub_gene_signal.tsv` by gene symbol and show `gb_log2fc`

**48c — DiffBind K119ub fold change at CpG islands**
- `findOverlaps()` between cpgi_gr and DiffBind K119ub peaks
- For each island, take the overlapping DiffBind peak with max |Fold|
- Violin/box split by sig category; Wilcoxon test
- Prediction: hyper islands should show positive K119ub fold (gain in mut)

#### Part 2: De Novo vs Pre-existing Methylation

**48d — Baseline mC density plot**
- Overlapping density of `mc_mean_mod_control` for hyper DMR / hypo DMR / non-significant
- Vertical reference lines at 0.05, 0.10, 0.20 thresholds
- Colors: hyper `#D7191C`, hypo `#2C7BB6`, non-sig `grey70`

**48e — Control vs mutant mC scatter**
- All 8,910 islands; highlight significant DMRs by direction
- Diagonal line (y=x), vertical lines at 0.10 and 0.20
- Points above diagonal = gain, below = loss
- Label top-10 highest-gain hyper islands with `ggrepel`

**48f — De novo classification stacked bars**
- For 122 hyper islands only: classify at thresholds 0.05, 0.10, 0.20
- Stacked bar showing de novo vs pre-existing counts at each threshold
- Include parallel analysis for hypo islands (high baseline loss)

**48g — Methylation gain magnitude comparison**
- For hyper islands: box/violin of `mc_mod_difference` split by de novo (ctrl < 0.20) vs pre-existing (ctrl >= 0.20)
- Wilcoxon test if n >= 5 per group

#### Part 3: Coordinated mC + hmC Changes

**48h — mC vs hmC difference scatter**
- x = `mc_mod_difference`, y = `hmc_mod_difference` for all 8,910 islands
- Color by 4-category sig status: Both sig / mC only / hmC only / Neither
- Quadrant count annotations; dashed lines at x=0, y=0

**48i — Co-significant island heatmap**
- Filter to islands sig in both mC AND hmC
- Matrix: mc_ctrl, mc_mut, hmc_ctrl, hmc_mut
- Row annotation: direction combo (mC gain+hmC loss, concordant, etc.)
- `pheatmap` with row clustering, no column clustering

#### Part 4: K119ub x De Novo Integration (Critical)

**48j — Fisher's test forest plot**
- Three OR estimates with CIs:
  1. All islands: mc_hyper ~ K119ub_gained overlap
  2. Among hyper: de_novo (ctrl < 0.20) ~ K119ub_gained
  3. Among hyper: pre_existing (ctrl >= 0.20) ~ K119ub_gained
- Forest plot with log-scale x-axis
- Note small sample sizes explicitly

**48k — 2x2 contingency tile plot**
- For 122 hyper islands: K119ub gained (yes/no) x methylation baseline (de novo / pre-existing at 0.20)
- Counts, percentages, Fisher's OR + p-value in subtitle

#### Part 5: Chromatin Context

**48l — Chromatin mark overlap bar chart**
- H3K27me3, H3K27ac, H3K4me3 overlap rates at hyper / hypo / non-sig islands
- Fisher's test for each mark x category vs background
- Prediction: hyper enriched for H3K27me3 (Polycomb), depleted for H3K27ac

**48m — Enrichment OR heatmap**
- 3 marks x 3 categories matrix of log2(OR)
- `pheatmap`, diverging blue-white-red palette centered at 0

### Summary Table

**`TABLES_DIR/48_cpg_island_ubiquitination_summary.tsv`** — 8,910 rows with: coordinates, mC/hmC metrics, K119ub overlaps (ctrl/mut/gained/lost), DiffBind fold, de novo classification (3 thresholds), chromatin overlaps, sig status.

### Reused Functions (from `_shared_config.R`)

| Function | Used For |
|----------|----------|
| `load_dmr_bed()` | Loading hmC CpG island BED if needed |
| `load_chip_peaks()` | Loading K119ub and ChIP peak BEDs → GRanges |
| `load_diffbind_flex()` | Loading K119ub DiffBind results |
| `dmr_to_granges()` | Converting DMR data frames to GRanges |
| `theme_biomodal()` | All ggplot2 figures |
| `save_multiformat_ggplot()` | All ggplot outputs (PDF + PNG + SVG) |
| `COLORS$direction` | Hyper/hypo coloring |

### New Helper Functions (defined in section script)

| Function | Purpose |
|----------|---------|
| `load_cpg_island_universe(tsv_path)` | Load TSV, add chr prefix, compute significance flags |
| `make_cpgi_granges(df)` | Data frame → GRanges with all metadata |
| `run_fisher_2x2(a, b, c, d, ...)` | Standard 2x2 Fisher's exact test wrapper |
| `derive_differential_peaks(ctrl_gr, mut_gr)` | Reuse pattern from section_14 |
| `classify_de_novo(df, threshold)` | Classify hyper islands at given baseline threshold |

### Execution Order

1. Source config + load packages (ChIPseeker, ggforce)
2. Load CpG island universe TSV → `cpgi` data frame
3. Build `cpgi_gr` GRanges (with chr prefix)
4. Load K119ub peaks (ctrl, mut, gained, lost) + DiffBind
5. Load ChIP peaks (H3K27me3, H3K27ac, H3K4me3)
6. Compute all overlaps → add columns to `cpgi`
7. Classify de novo at 3 thresholds → add columns
8. Generate figures 48a-48m
9. Write summary table
10. Print console summary

## Verification

```bash
cd /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream
Rscript scripts/viz_sections/section_48_cpg_island_ubiquitination.R
```

Check:
- 13 figures in `plots/visualizations/48_cpg_island_ubiquitination/` (PDF + PNG + SVG each)
- Summary table in `plots/visualizations/tables/48_cpg_island_ubiquitination_summary.tsv` (8,910 rows)
- Console output shows Fisher's test results and summary statistics
- No warnings about chromosome mismatches (coordinate system handling)
