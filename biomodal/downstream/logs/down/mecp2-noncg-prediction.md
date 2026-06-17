# Plan: MeCP2 TAD Methylation Analysis

## Context

The existing TAD methylation pipeline (`compute_tad_methylation.R` → `section_54`) quantifies how 4 methylation contexts (CG 5mC, CG 5hmC, CHG 5mC, CHH 5mC) are organized within 12,141 TAD domains. Key finding: 5hmC shows a 2x inter/intra-TAD variance ratio while 5mC is ~1x, and CHG/CHH are at the noise floor.

PI wants MeCP2 ChIP-seq added as a proxy for non-CG methylation that the DUET kit can't detect. Two additions:
1. **Naive**: Add raw MeCP2 signal as another context in the existing framework
2. **CG-corrected**: Regression-based residual analysis to isolate MeCP2 binding NOT explained by CG methylation, then test whether that residual is TAD-organized

## Data Available (all local)

| Data | Path | Notes |
|------|------|-------|
| MeCP2 ctrl BigWig | `/Users/zakiralibhai/sdsc/bigwigs/MeCP2Ctrl.bw` | 62 MB, merged |
| MeCP2 mut BigWig | `/Users/zakiralibhai/sdsc/bigwigs/MeCP2Mut.bw` | 60 MB, merged |
| CG mC ctrl BigWig | `/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_mc_ctrl.bw` | 212 MB, merged |
| CG mC mut BigWig | `/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_mc_mut.bw` | 213 MB, merged |
| CG hmC ctrl BigWig | `/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_hmc_ctrl.bw` | 187 MB, merged |
| CG hmC mut BigWig | `/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_hmc_mut.bw` | 189 MB, merged |
| MeCP2 peaks (up) | `peaks/mecp2/MeCP2_up.bed` | 7,686 peaks |
| MeCP2 peaks (down) | `peaks/mecp2/MeCP2_down.bed` | 1,200 peaks |
| Existing TAD signal | `downstream/data/tad_methylation_signal_late.tsv` | 12,141 rows × 61 cols |

No HPC compute needed — all BigWigs are local.

---

## Change 1: Naive MeCP2 in TADs

### New file: `downstream/scripts/tad_methylation/add_mecp2_tad_signal.R`

Local R script that augments the existing TSV with MeCP2 columns. Runs on Mac (~5-15 min).

**What it computes** (same 4 metrics as every methylation context, plus peak overlaps):

| Column | Description |
|--------|-------------|
| `mecp2_ctrl_mean` | Mean MeCP2 ChIP signal per TAD (ctrl BigWig) |
| `mecp2_mut_mean` | Mean MeCP2 ChIP signal per TAD (mut BigWig) |
| `mecp2_log2fc` | log2(mut/ctrl) with 5th-percentile pseudocount |
| `mecp2_pseudocount` | The pseudocount used |
| `mecp2_ctrl_cv` | Intra-TAD coefficient of variation (10kb bins, ctrl) |
| `mecp2_mut_cv` | Intra-TAD CV (mut) |
| `mecp2_ctrl_within_var` | Within-domain variance of 10kb bin means (ctrl) |
| `mecp2_mut_within_var` | Within-domain variance (mut) |
| `mecp2_ctrl_insulation` | Boundary insulation score, 50kb flanks (ctrl) |
| `mecp2_mut_insulation` | Boundary insulation score (mut) |
| `mecp2_peak_count` | Total MeCP2 DiffBind peaks overlapping the TAD |
| `mecp2_up_count` | MeCP2-Up peaks (higher in mut) overlapping the TAD |
| `mecp2_down_count` | MeCP2-Down peaks (lower in mut) overlapping the TAD |

**Implementation approach:**
- Reads existing TSV, reconstructs `domain_gr` GRanges from chr/start/end (no TADCompare re-read)
- Replicates the 4 signal extraction functions from `compute_tad_methylation.R` (they're pure functions; avoids cross-script sourcing dependency)
- Uses `GenomicRanges::countOverlaps()` for peak counts
- CLI args: `--input`, `--output`, `--ctrl-bw`, `--mut-bw`, `--peaks-up`, `--peaks-down`
- Overwrites TSV in-place: 61 → 74 columns

### Modify: `section_54_tad_methylation_organization.R`

3-line change to include MeCP2 in all existing figure panels:

```r
CONTEXTS <- c("cg_mc", "cg_hmc", "chg_mc", "chh_mc", "mecp2")

CONTEXT_LABELS — add: mecp2 = "MeCP2"
CONTEXT_COLORS — add: "MeCP2" = "#984EA3"  (purple)
```

Every panel (54a-54j) iterates over `CONTEXTS` generically via `paste0(ctx, "_", cond, "_mean")` etc., so MeCP2 appears automatically. The `scales = "free_y"` faceting already handles the different signal scales (RPM vs fractional methylation).

### Modify: `_shared_config.R`

Add BigWig paths to the existing `MECP2_FILES` list and add a new `METHYLATION_BIGWIGS` list:

```r
MECP2_FILES$ctrl_bw <- "/Users/zakiralibhai/sdsc/bigwigs/MeCP2Ctrl.bw"
MECP2_FILES$mut_bw  <- "/Users/zakiralibhai/sdsc/bigwigs/MeCP2Mut.bw"

METHYLATION_BIGWIGS <- list(
  cg_mc_ctrl  = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_mc_ctrl.bw",
  cg_mc_mut   = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_mc_mut.bw",
  cg_hmc_ctrl = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_hmc_ctrl.bw",
  cg_hmc_mut  = "/Users/zakiralibhai/sdsc/bigwigs/methylation/merged/CG_hmc_mut.bw"
)
```

---

## Change 2: CG-Corrected MeCP2 Proxy (Section 55)

### New file: `downstream/scripts/viz_sections/section_55_mecp2_cg_corrected_proxy.R`

Self-contained viz section. Sources `_shared_config.R`. Reads the MeCP2-augmented TSV from Change 1, PLUS reads BigWigs directly for bin-level regression.

### Statistical Methodology

#### Primary model (Model A): `MeCP2 ~ CG_mC`

Per-TAD OLS regression. Fit on ctrl data only (unperturbed state as baseline):

```r
model_a <- lm(mecp2_ctrl_mean ~ cg_mc_ctrl_mean, data = tad)
```

- **Residuals**: `resid_ctrl = mecp2_ctrl - predicted(cg_mc_ctrl)`
- **Apply to mutant**: Use ctrl-derived coefficients with mut CG values to predict mut MeCP2. `resid_mut = mecp2_mut - predict(model_a, newdata=cg_mc_mut_mean)`. This asks: "In the mutant, how much MeCP2 exceeds what the normal CG-MeCP2 relationship would predict?"
- Positive residuals = candidate non-CG MeCP2 binding sites

#### Secondary model (Model B): `MeCP2 ~ CG_mC + CG_hmC`

Removes ALL CG-context MeCP2 affinity (both 5mC and 5hmC):

```r
model_b <- lm(mecp2_ctrl_mean ~ cg_mc_ctrl_mean + cg_hmc_ctrl_mean, data = tad)
```

Report both sets of residuals. If they tell the same story → robust result. If they diverge → hmC contributes meaningfully and the divergence itself is informative.

#### Model diagnostics

- Residuals vs fitted plot (check linearity)
- Q-Q plot (check normality of residuals)
- If non-linear pattern detected, report it and add a LOESS/GAM sensitivity check (`mgcv::gam(mecp2 ~ s(cg_mc))`)

#### Bin-level variance ratio (rigorous intra-TAD computation)

For the inter/intra-TAD variance ratio of residuals, we need bin-level residuals (not just per-TAD scalars). The script reads BigWigs directly:

1. Tile each TAD into 10kb bins
2. Extract per-bin MeCP2 signal AND per-bin CG mC signal from merged BigWigs
3. Apply TAD-level regression coefficients to predict per-bin MeCP2
4. Compute per-bin residual = actual - predicted
5. Within-TAD variance = `var(bin_residuals)` per domain
6. Inter-TAD variance = `var(per-TAD mean residuals)` across all domains
7. Variance ratio = inter / mean(intra)

This is more rigorous than using raw `mecp2_ctrl_within_var` as the denominator, because the CG correction removes different amounts from different bins (CG varies within TADs, CV ~0.28).

#### CG-stratified residual analysis (the "extra" comparison)

Using the existing `cg_mc_log2fc` column from the TAD signal TSV, stratify TADs by CG methylation change:

1. **Quartile stratification**: Split TADs into quartiles by `cg_mc_log2fc` (Q1 = most hypomethylated, Q4 = most hypermethylated)
2. **Test per stratum**: Are MeCP2 residuals systematically positive in TADs where CG methylation is NOT changing (middle quartiles)?
3. **Key biological question**: If MeCP2 excess exists even at TADs without CG hypermethylation, that's strong evidence for non-CG-driven binding

Additionally, overlap with CG DMR data:
- Load gene-body CG mC DMRs from `DATA_PATHS$mc_dmr` (already in shared config)
- Classify TADs as "contains hypermethylated DMR gene" vs "no DMR" vs "contains hypomethylated DMR gene"
- Compare MeCP2 residuals across these groups
- Wilcoxon tests for each contrast

### TSV Outputs (all saved to `plots/visualizations/tables/`)

| File | Contents |
|------|----------|
| `55_mecp2_regression_coefficients.tsv` | Model A & B: intercept, slope(s), R², adj-R², F-stat, p-value, residual SE, n |
| `55_mecp2_tad_residuals.tsv` | Per-TAD: chr, start, end, mecp2 raw values, cg_mc values, predicted, residual (model A), residual (model B), cg_mc_log2fc quartile, dmr_overlap_flag |
| `55_mecp2_binlevel_variance.tsv` | Per-TAD: within_var_raw, within_var_residual_a, within_var_residual_b (from bin-level computation) |
| `55_mecp2_variance_ratio_summary.tsv` | Variance ratios with bootstrap CIs: raw MeCP2, residual A, residual B, for both ctrl and mut |
| `55_mecp2_cg_stratified_summary.tsv` | Per CG-quartile: n_tads, mean_residual, median_residual, wilcox_p_vs_zero, mean_peak_count |
| `55_mecp2_dmr_stratified_summary.tsv` | Per DMR class (hyper/hypo/none): n_tads, mean_residual, wilcox_p between groups |

### Figure Panels

| Panel | Content |
|-------|---------|
| **55a** | MeCP2 vs CG 5mC scatter + OLS regression line (model A). Annotate R², slope, p. Color by `is_differential`. Shows how much MeCP2 is explained by CG mC. |
| **55b** | Model diagnostics: residuals-vs-fitted + Q-Q plot (2-panel). Quick visual check for non-linearity. |
| **55c** | Residual distributions: density plots of ctrl and mut residuals (model A). Overlay zero line. Wilcoxon test for mut shift from zero. |
| **55d** | **KEY FIGURE**: Variance ratio comparison bar plot with bootstrap CIs. 3 groups side-by-side: CG 5mC, raw MeCP2, residual MeCP2 (model A). Dashed line at ratio=1. Answers: "Is the non-CG component of MeCP2 TAD-organized?" |
| **55e** | Boundary correlation: MeCP2 residual vs `mean_tad_score` (Spearman). Tests whether non-CG MeCP2 preferentially sits at stronger TAD boundaries. |
| **55f** | **CG-stratified residual**: Boxplot of MeCP2 residuals across CG log2FC quartiles (Q1-Q4). One-sample Wilcoxon per quartile testing residual ≠ 0. Key test: are residuals positive even in Q2-Q3 (no CG change)? |
| **55g** | DMR-stratified residual: Violin/boxplot of residuals in TADs with hyper-DMR genes vs hypo-DMR genes vs no DMR. Pairwise Wilcoxon. |
| **55h** | Model B comparison: Same as 55d but overlaying model A vs model B variance ratios. Shows whether including hmC changes the TAD organization conclusion. |
| **55i** | Composite panel (55a + 55d + 55f) for presentation. |

Output directory: `plots/visualizations/55_mecp2_cg_corrected_proxy/`

---

## Files Summary

### Create (2 files)
1. `downstream/scripts/tad_methylation/add_mecp2_tad_signal.R` — local compute, adds 13 MeCP2 columns to TSV
2. `downstream/scripts/viz_sections/section_55_mecp2_cg_corrected_proxy.R` — regression analysis + 9 figure panels + 6 TSV outputs

### Modify (2 files)
3. `downstream/scripts/viz_sections/_shared_config.R` — add BigWig paths to `MECP2_FILES`, add `METHYLATION_BIGWIGS` list
4. `downstream/scripts/viz_sections/section_54_tad_methylation_organization.R` — extend `CONTEXTS`, `CONTEXT_LABELS`, `CONTEXT_COLORS` to include `mecp2`

### Data artifacts produced
5. `downstream/data/tad_methylation_signal_late.tsv` — overwritten: 61 → 74 columns
6. `downstream/plots/visualizations/tables/55_*.tsv` — 6 summary/data TSVs from section 55

---

## Execution Order

```
1. Modify _shared_config.R (add BigWig paths)
2. Create + run add_mecp2_tad_signal.R
   → Produces augmented TSV (74 columns)
3. Modify section_54 (3-line CONTEXTS change)
   → User runs: Rscript scripts/viz_sections/section_54_tad_methylation_organization.R
   → MeCP2 now appears in all 10 existing panels
4. Create section_55
   → User runs: Rscript scripts/viz_sections/section_55_mecp2_cg_corrected_proxy.R
   → 9 figure panels + 6 TSVs
```

The `run_all_sections.sh` glob (`section_*.R`) auto-discovers section_55.

## Verification

1. After step 2: confirm TSV has 74 columns, MeCP2 columns have non-zero medians, peak counts are non-negative integers
2. After step 3: visually check `54a_metatad_signal` and `54d_variance_ratio` now include a purple MeCP2 panel
3. After step 4: check `55_mecp2_regression_coefficients.tsv` for sensible R²; check `55d_variance_ratio` figure; confirm all 6 TSV outputs exist and have expected columns
