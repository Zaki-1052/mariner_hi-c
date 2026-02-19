# Step 10: K119ub-ABC Enhancer Correlation Analysis

## Position in the Pipeline

Step 10 sits between the cross-reference analyses (Steps 9/9b) and the enhancer subset stratified analysis (Step 11). It asks the first direct question about H2AK119ub at ABC-modeled enhancers:

**Does H2AK119ub accumulation at enhancers in BAP1-KO correlate with loss of enhancer-gene linkage (ABC score)?**

BAP1 is a deubiquitinase that removes H2AK119ub. In BAP1-KO, the prediction is straightforward: enhancers that lose ABC activity should show *increased* K119ub (negative correlation between delta-ABC and K119ub log2FC). This step tests that prediction genome-wide across 31,947 quantifiable enhancers.

### Relationship to Other Steps

- **Step 7** established that 58.8% of genes with both significant delta-ABC and differential expression are concordant (binomial p = 6.84e-08), and that unnormalized delta(AxC) outperforms normalized delta-ABC (65.3% vs 58.8%).
- **Step 9** showed that independent loop-ABC overlap concordance is at chance (51.4%).
- **Step 10** (this analysis) adds the K119ub dimension: do ABC changes correlate with the chromatin mark that BAP1 directly regulates?
- **Step 9b** later shows that when geometric constraints are applied (paired-anchor), K119ub correlation jumps to rho = -0.401 (vs -0.023 to -0.090 here).
- **Step 11** further dissects this by classifying enhancers into epigenomic subsets and asking whether K119ub_Only enhancers show contact-level phenotypes.

---

## Input Data

### Preprocessing (HPC)

K119ub signal was extracted on Expanse HPC using `preprocess_k119ub_enhancer_signal.R`:

- **Source:** 8 H2AK119ub ChIP-seq bigWig files (4 ctrl + 4 mut replicates) from `heatmaps/`
- **Regions:** 33,004 unique enhancer coordinates from the 180,423 E-G pairs in `delta_abc_all_pairs.tsv`
- **Method:** Mean bigWig signal per enhancer region via `rtracklayer::import.bw()` + `IRanges::viewMeans()`
- **Normalization:** Median-ratio normalization across all 8 samples (geometric mean of per-sample non-zero medians as reference)
- **Log2FC:** `log2((mut_mean + pseudocount) / (ctrl_mean + pseudocount))` with adaptive pseudocount (5th percentile of non-zero normalized signal)
- **Signal classes:** "quantifiable" (both conditions above pseudocount), "one_condition" (only one above), "no_signal" (both below)

### Analysis Inputs (Local)

| Source | Rows | Description |
|--------|------|-------------|
| `delta_abc_all_pairs.tsv` | 180,423 | All distal E-G pairs with delta-ABC and delta-unnorm |
| `k119ub_enhancer_signal.tsv` | 33,004 | Per-replicate K119ub at unique enhancers (from HPC) |
| `delta_abc_annotated.tsv` | 180,423 | E-G pairs with H3K27ac overlap flag |

---

## Analysis Design

### Part A: Enhancer-Level Aggregation

Each enhancer maps to multiple target genes (mean ~5.5 per enhancer). The script collapses 180,423 E-G pairs to 33,004 unique enhancers by:

- **delta_activity:** Taken from first row per enhancer (identical for all genes from the same enhancer, since activity = ATAC signal at the enhancer itself)
- **mean_delta_abc / mean_delta_unnorm:** Averaged across all target genes per enhancer
- **H3K27ac flag:** From `delta_abc_annotated.tsv` (any True = True for that enhancer)

After joining with K119ub signal: **31,947 enhancers** with quantifiable K119ub in both conditions.

### Part B: Core Spearman Correlations

Five correlations tested, all on n = 31,947 enhancers:

| Correlation | rho | 95% CI | p-value | Interpretation |
|---|---|---|---|---|
| delta_activity vs K119ub log2FC | **-0.039** | [-0.050, -0.028] | 2.0e-12 | Activity loss → K119ub gain |
| mean delta(AxC) vs K119ub log2FC | **-0.090** | [-0.101, -0.079] | 1.6e-58 | **Strongest** — unnormalized best captures the effect |
| mean delta-ABC vs K119ub log2FC | **-0.023** | [-0.034, -0.012] | 3.4e-05 | Weaker after normalization |
| delta_activity vs K119ub mut_mean | **+0.001** | [-0.010, +0.012] | 0.89 | Null control — absolute signal irrelevant |
| K119ub ctrl vs mut (QC) | **+0.837** | [+0.834, +0.841] | ~0 | Replicate agreement — data is solid |

**All three ABC-vs-K119ub correlations are negative**, confirming the biological prediction: enhancers that lose ABC activity in BAP1-KO gain H2AK119ub.

The **unnormalized delta(AxC)** has the strongest correlation (rho = -0.090), consistent with the Step 7 finding that unnormalized scores better predict expression direction (65.3% vs 58.8% concordance). This is expected: the ABC normalization divides by the sum of all Activity x Contact for each gene, compressing real activity changes when BAP1-KO causes widespread chromatin remodeling across many enhancers for the same gene.

The delta_activity vs K119ub_mut_mean correlation is null (rho = +0.001, p = 0.89), confirming that the relationship is with the *change* in K119ub, not with absolute K119ub levels in KO.

### Part C: ABC Category Comparison

Enhancers were classified by mean delta-ABC:

| Category | Median K119ub log2FC | Mean K119ub log2FC | N |
|---|---|---|---|
| **Lost** (delta-ABC < -0.01) | **+0.005** | **+0.011** | 4,861 |
| Unchanged (\|delta-ABC\| <= 0.01) | **-0.037** | **-0.035** | 21,177 |
| Gained (delta-ABC > 0.01) | **-0.015** | **-0.014** | 5,909 |

Kruskal-Wallis: chi-sq = 81.03, p = 2.54e-18.

The **Lost group has the most positive K119ub log2FC** (median +0.005 vs -0.037 for Unchanged). This is the BAP1 mechanism in action: without deubiquitinase activity, K119ub accumulates at enhancers, ATAC/activity drops, and ABC scores decrease.

Pairwise Wilcoxon tests:

| Comparison | p-value | Rank-biserial r |
|---|---|---|
| Lost vs Unchanged | 1.42e-17 | -0.078 |
| Gained vs Unchanged | 5.47e-06 | -0.039 |
| Lost vs Gained | 4.42e-04 | -0.039 |

All comparisons are significant. The Lost vs Unchanged comparison has the strongest effect (r = -0.078), though the effect size is modest.

### Part D: Stratified Analyses

**D1: By enhancer-gene distance.**

| Distance bin | rho | p-value | N |
|---|---|---|---|
| < 50 kb | -0.055 | 2.1e-11 | 14,975 |
| 50-200 kb | -0.063 | 5.1e-09 | 8,582 |
| 200-500 kb | -0.039 | 1.0e-02 | 4,358 |
| 500 kb-1 Mb | -0.050 | 1.8e-02 | 2,216 |
| 1-5 Mb | -0.044 | 6.0e-02 | 1,816 |
| > 5 Mb | — | — | 0 |

The correlation is consistent across all distance bins (<50 kb through 500 kb-1 Mb), ruling out a distance-driven artifact. The 1-5 Mb bin falls below significance (p = 0.06), likely a power issue with n = 1,816 rather than a true biological difference.

**D2: By H3K27ac status.**

| Stratification | Metric | rho | p-value | N |
|---|---|---|---|---|
| H3K27ac+ | delta_activity vs K119ub | **-0.050** | 1.9e-08 | 12,837 |
| H3K27ac+ | delta_unnorm vs K119ub | **-0.082** | 8.6e-21 | 12,837 |
| H3K27ac- | delta_activity vs K119ub | **-0.017** | 1.9e-02 | 19,110 |
| H3K27ac- | delta_unnorm vs K119ub | **-0.069** | 1.5e-21 | 19,110 |

The activity-K119ub correlation is **~3x stronger at H3K27ac+ enhancers** (rho = -0.050) than H3K27ac- (rho = -0.017). This is biologically coherent: BAP1 loss hits active enhancers hardest because those are where deubiquitination was functionally relevant. Active enhancers (H3K27ac+) rely on BAP1 to keep K119ub low; when BAP1 is lost, these enhancers accumulate K119ub and lose activity.

For the unnormalized metric, both H3K27ac+ and H3K27ac- show strong correlations (-0.082 and -0.069), suggesting the contact-level effect captured by delta(AxC) is less dependent on active enhancer status.

---

## Results and Interpretation

### Biological Prediction: Confirmed

The core prediction — negative correlation between ABC changes and K119ub changes — is confirmed across all three ABC metrics. The direction is correct and statistically unambiguous (p-values range from 3.4e-05 to 1.6e-58). Enhancers that lose regulatory function in BAP1-KO (lower ABC score) accumulate H2AK119ub, consistent with BAP1 deubiquitinase activity being required to maintain enhancer accessibility.

### Honest Assessment of Effect Sizes

The rho values (-0.023 to -0.090) are **statistically robust but modest in magnitude**. The strongest correlation (unnormalized delta(AxC), rho = -0.090) explains less than 1% of the variance (R^2 ~ 0.8%). The scatter plots (Panels 01-03) are visually noisy — the red regression line has a perceptible negative slope, but the cloud of points is dominated by scatter.

This is expected for a genome-wide enhancer analysis with 32K data points. Most enhancers are not meaningfully affected by BAP1 loss — the signal is carried by the ~5K "Lost" enhancers, and the category comparison (KW p = 2.5e-18) shows the group-level effect is unambiguous even though the continuous correlation is weak.

**Context from downstream analyses:** Step 9b's paired-anchor analysis shows rho = -0.401 (p = 5.5e-13) when restricting to the 490 enhancer-gene-loop triplets where the geometric constraint confirms a physical loop-mediated connection. The 4.5x increase in correlation strength (from -0.090 to -0.401) demonstrates that the weak genome-wide signal in Step 10 is diluted by the vast majority of enhancers that are not connected to their target genes by specific detectable loops. The biological effect is real but concentrated at a subset of high-confidence regulatory loci.

### QC Check: Passed

K119ub ctrl vs mut correlation: rho = +0.837 (CI: [+0.834, +0.841]). This strong replicate agreement confirms the bigWig data is solid and the median-ratio normalization preserved the biological signal.

---

## Summary of Claims and Evidence Strength

| Claim | Evidence | Strength |
|-------|----------|----------|
| Enhancers that lose ABC score gain K119ub | Spearman rho = -0.023 to -0.090, all p < 3.4e-05 | **Statistically strong**, modest effect size |
| Unnormalized delta(AxC) is the best ABC metric for K119ub correlation | rho = -0.090 vs -0.039 (activity) vs -0.023 (normalized) | **Consistent** with Step 7 findings |
| Lost enhancers have most positive K119ub change | KW p = 2.5e-18, Lost median = +0.005 vs Unchanged -0.037 | **Strong** — clean group separation |
| Correlation is distance-independent | Consistent rho across <50kb to 1Mb bins | **Moderate** — not a distance artifact |
| Effect is stronger at H3K27ac+ enhancers | rho 3x larger at H3K27ac+ for activity metric | **Moderate** — biologically coherent |
| K119ub absolute level is irrelevant | rho = +0.001 (p = 0.89) for mut_mean | **Strong** null — it's the change that matters |

---

## Limitations

1. **ATAC-only ABC mode.** Activity is defined purely from ATAC signal, not H3K27ac (BAMs unavailable). The ABC model benchmarks best with DNase + H3K27ac. ATAC-only mode may underestimate enhancer activity, particularly at elements where H3K27ac and ATAC diverge.

2. **Enhancer-level aggregation averages over genes.** Each enhancer's delta-ABC is the mean across all its target genes. An enhancer with 7 targets gets equal weight to one with 1 target. If the K119ub effect is gene-specific (e.g., stronger at genes with fewer enhancers), this averaging could dilute the signal.

3. **K119ub bigWig normalization.** Median-ratio normalization across 8 samples assumes most regions are unchanged between conditions. For a Polycomb mark in BAP1-KO where genome-wide accumulation is expected, this assumption may not fully hold. The QC correlation (rho = +0.837) suggests the normalization is reasonable, but a systematic upward bias in all KO samples cannot be fully excluded.

4. **H3K27ac annotation is static.** The H3K27ac peaks used for stratification are from a single adult timepoint merge (not condition-specific DiffBind). They tell you which enhancers had active marks in the reference, not whether those marks changed between WT and KO. Step 11 addresses this with DiffBind-derived epigenomic classification.

5. **Correlation does not imply causation.** The negative rho between K119ub gain and ABC loss could reflect: (a) K119ub accumulation causing activity/contact loss (the BAP1 mechanism), (b) activity loss causing K119ub accumulation (loss of active marks exposes chromatin to PRC1), or (c) a shared upstream cause. This analysis cannot distinguish these models.

6. **Weak effect sizes.** Rho values of -0.023 to -0.090 are at the lower end of what is biologically interpretable. While the sample size (31,947) gives extraordinary statistical power, the practical significance is limited. The vast majority of individual enhancers show K119ub and ABC changes that are not meaningfully correlated.

---

## Recommended Figures for Presentation

Based on the analysis results and following the assessment in the original interpretation notes:

### Lead with Panel 06 (Violin + Box) — "The Punchline"

This is the most presentation-friendly figure. It immediately communicates that enhancers which **lost** ABC score in BAP1-KO have a K119ub distribution shifted upward relative to Unchanged/Gained categories. The violin shape shows the full distribution, the box overlay gives medians, and the visual is intuitive without needing to explain correlation coefficients. Your mentor can grasp the result in 5 seconds.

### Follow with Panel 02 (Unnormalized Scatter) — "The Evidence"

This is the strongest quantitative result (rho = -0.090, p = 1.6e-58). The negative slope is visible even with 32K noisy points. Show this second to demonstrate it's not just a categorical artifact — it's a continuous genome-wide relationship. The unnormalized ABC being the strongest predictor also ties back to the Step 7 finding that unnormalized scores outperform normalized (65.3% vs 58.8% concordance with RNA-seq).

### Why Not the Others?

- **Panel 04** (boxplot) shows the same result as 06 but less informatively — the violin adds distributional shape
- **Panel 08** (H3K27ac stratification) is important for the paper but is a secondary analysis — hold it for when asked "is this driven by active enhancers?" (answer: yes, rho is 3x stronger at H3K27ac+ enhancers)
- **Panel 07** (distance bins) is a robustness/controls check, not a lead figure
- **Panel 05** (absolute K119ub in KO) shows that absolute levels don't differ much between categories, reinforcing that it's the *change* that matters — useful as a supplementary control

### The Two-Slide Narrative

> "Enhancers that lose regulatory function in BAP1-KO (lower ABC score) accumulate H2AK119ub (Panel 06). This is a continuous genome-wide relationship (Panel 02, rho = -0.09, p < 10^-58), consistent with BAP1 deubiquitinase activity being required to maintain enhancer accessibility."

---

## Output File Descriptions

### Plots (8 panels, each in PDF + PNG)

| # | File | Content |
|---|------|---------|
| 01 | `01_delta_activity_vs_k119ub` | Density scatter: ATAC activity change vs K119ub log2FC (rho = -0.039) |
| 02 | `02_delta_unnorm_vs_k119ub` | Density scatter: unnormalized delta(AxC) vs K119ub log2FC (rho = -0.090, **strongest**) |
| 03 | `03_delta_abc_vs_k119ub` | Density scatter: normalized delta-ABC vs K119ub log2FC (rho = -0.023) |
| 04 | `04_boxplot_k119ub_by_abc_category` | Boxplot: K119ub log2FC by Lost/Unchanged/Gained (KW p = 2.5e-18) |
| 05 | `05_boxplot_k119ub_mut_by_category` | Boxplot: absolute K119ub signal in KO by category (log1p scale) |
| 06 | `06_violin_k119ub_by_category` | Violin + box: K119ub log2FC distribution by category (**recommended lead figure**) |
| 07 | `07_scatter_by_distance` | Faceted scatter: activity vs K119ub by enhancer-gene distance bin (robustness check) |
| 08 | `08_scatter_by_h3k27ac` | Faceted scatter: activity vs K119ub by H3K27ac status (H3K27ac+ 3x stronger) |

Note: Figures are saved as flat PDF + PNG files (not in the multi-format subfolder structure used in later steps). The Cairo PDF backend produced warnings on macOS due to missing X11 libraries, but files were saved successfully using the fallback R PDF device.

### TSV Files

| File | Rows | Description |
|------|------|-------------|
| `k119ub_abc_correlation_summary.tsv` | 5 | Core Spearman correlations with rho, p-value, n, and 95% CIs |
| `k119ub_abc_enhancer_merged.tsv` | 33,004 | Per-enhancer: ABC metrics + K119ub per-replicate signal + log2FC + H3K27ac flag |

### Interpretation Notes

`results/figures/k119ub_correlation/interpretation.md` contains the original analysis interpretation with figure recommendations, written at the time the analysis was first run.

---

## Bottom Line

**What the data shows:** Enhancers that lose ABC enhancer-gene linkage in BAP1-KO accumulate H2AK119ub, as predicted by the BAP1 deubiquitinase mechanism. The relationship is genome-wide, continuous, distance-independent, and strongest at H3K27ac-marked active enhancers. Unnormalized delta(AxC) captures the effect best (rho = -0.090), consistent with Step 7's finding that ABC normalization compresses real changes.

**What the data does not show:** The correlations are modest (rho < 0.1, R^2 < 1%). While statistically unambiguous with 32K data points, the per-enhancer relationship is noisy — most individual enhancers do not show a clean K119ub-to-ABC correspondence. This is a genome-wide trend, not a per-locus predictor.

**Interpretation:** Step 10 establishes the genome-wide correlation that downstream analyses refine. The signal-to-noise ratio improves dramatically when analyses impose geometric constraints: Step 9b's paired-anchor analysis (rho = -0.401 at enhancers confirmed by loop geometry) and Step 11's enhancer subset analysis (phenotypic classification isolates the K119ub_Only class) both build on this foundation. Step 10's contribution is demonstrating that the BAP1 mechanism operates at the level of individual enhancers across the genome, not just at a handful of outlier loci — even though the per-enhancer effect is small.
