

**Language > Specialist**: R > Epigenomics / Chromatin Biology

## Results Interpretation

The analysis ran successfully on 31,947 quantifiable enhancers. Here's the verdict:

### Biological Prediction: CONFIRMED

All three ΔABC vs ΔK119ub correlations are **negative**, as predicted — enhancers that lose ABC activity in BAP1-KO gain H2AK119ub:

| Correlation | rho | p-value | Interpretation |
|---|---|---|---|
| Δactivity vs K119ub log2FC | **-0.039** | 2.0e-12 | Lost activity → gained ub |
| mean Δ(A×C) vs K119ub log2FC | **-0.090** | 1.6e-58 | Strongest — unnormalized ABC best captures effect |
| mean ΔABC vs K119ub log2FC | **-0.023** | 3.4e-05 | Weaker after normalization (expected) |

The **unnormalized ABC** (Δ(A×C)) has the strongest correlation (rho = -0.090), consistent with your Step 7 finding that unnormalized scores are better predictors. This makes sense: the normalization dilutes the signal by dividing by genome-wide totals.

### Category Comparison: CLEAN SEPARATION

| Category | Median K119ub log2FC | n |
|---|---|---|
| **Lost** (ΔABC < -0.01) | **+0.005** | 4,861 |
| Unchanged | -0.037 | 21,177 |
| Gained (ΔABC > 0.01) | -0.015 | 5,909 |

Kruskal-Wallis p = 2.5e-18. The **Lost group has the most positive K119ub log2FC** (median +0.005 vs -0.037 for unchanged), meaning enhancers that lost ABC score gained K119ub in KO. This is exactly the BAP1 mechanism: without deubiquitinase activity, K119ub accumulates → chromatin closes → ATAC/activity drops → ABC score drops.

### Stratification Insights

**By distance:** The correlation is consistent across all bins (<50kb through 500kb-1Mb), suggesting this isn't a distance-driven artifact. Falls below significance only at >1Mb (n=1,816, p=0.06) — likely power-limited.

**By H3K27ac:** The activity-K119ub correlation is **stronger at H3K27ac+ enhancers** (rho = -0.050) than H3K27ac- (rho = -0.017). This fits: BAP1 loss hits active enhancers hardest because those are where deubiquitination was functionally relevant.

### QC Check: PASSED

K119ub ctrl vs mut correlation: rho = +0.837 — strong replicate agreement confirms the bigWig data is solid.

### Honest Assessment of Effect Sizes

The rho values (-0.02 to -0.09) are **statistically robust but modest in magnitude**. This is expected for a genome-wide enhancer analysis with 32K points — most enhancers aren't meaningfully affected by BAP1 loss. The signal is carried by the ~5K "Lost" enhancers, and the category comparison (Kruskal-Wallis p = 2.5e-18) shows the group-level effect is unambiguous. The scatter plots will look noisy but the trend is real and directionally consistent with the mechanism.

### Cairo Warning

The "failed to load cairo DLL" warnings are cosmetic — PDFs still saved, just using the default R PDF device instead of cairo. The PNGs are fine.


---


I'd recommend showing **two figures**, which together tell the complete story:

### 1. Panel 06 (violin + box) — the "punchline"

This is your most presentation-friendly figure. It immediately communicates: enhancers that **lost** ABC score in BAP1-KO have a K119ub distribution shifted upward relative to unchanged/gained. The violin shape shows the full distribution, the box overlay gives medians, and the visual is intuitive without needing to explain correlation coefficients. Lead with this — your mentor can grasp the result in 5 seconds.

### 2. Panel 02 (Δ(A×C) scatter) — the "evidence"

This is the strongest quantitative result (rho = -0.090, p = 1.6e-58). The negative slope is visible even with 32K noisy points. Show this second to demonstrate it's not just a categorical artifact — it's a continuous genome-wide relationship. The unnormalized ABC being the strongest predictor also ties back to your Step 7 finding nicely.

### Why not the others?

- **Panel 04** (boxplot) shows the same thing as 06 but less informatively — the violin adds distributional shape
- **Panel 08** (H3K27ac stratification) is important for the paper but is a secondary analysis — save it for when the mentor asks "is this driven by active enhancers?" (answer: yes, rho is 3x stronger at H3K27ac+ enhancers)
- **Panel 07** (distance bins) is a controls/robustness check, not a lead figure

### The narrative in two slides

> "Enhancers that lose regulatory function in BAP1-KO (lower ABC score) accumulate H2AK119ub (Panel 06). This is a continuous genome-wide relationship (Panel 02, rho = -0.09, p < 10^-58), consistent with BAP1 deubiquitinase activity being required to maintain enhancer accessibility."