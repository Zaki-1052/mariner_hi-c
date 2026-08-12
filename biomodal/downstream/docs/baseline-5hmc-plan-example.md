# Plan: Section 23 — Baseline 5hmC as Predictor of DMR Susceptibility

## Context

**Problem:** TODO #3 asks whether baseline (wildtype) 5hmC level predicts which genes become DMRs. This tests the "substrate availability" interpretation: genes with more 5hmC have more substrate for TET impediment to affect. Section 18 already showed K119ub is a weak gene-specific predictor (Cliff's delta 0.089) — if baseline 5hmC is stronger, it confirms the revised conclusions about TET-mediated demethylation block.

**Outcome:** A new visualization section (`section_23_baseline_hmc_predictor.R`) with 5 figures (23a–23e), 3 logistic regression models, and exported tables.

---

## Implementation

### New file
`downstream/scripts/viz_sections/section_23_baseline_hmc_predictor.R`

### Data sources (all available locally)

| Source | Path (relative to `downstream/`) | Join key | Provides |
|--------|----------------------------------|----------|----------|
| hmC extract | `modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260221_185106/Extract_hmc_regional-frac_20260221_185106.tsv.gz` | `Name` (col 9) | Per-sample 5hmC regional fractions — average ctrl-M (col 4) + ctrl-F (col 6) → baseline WT 5hmC |
| hmC DMR | Pre-loaded `hmc_dmr` from `_shared_config.R` | `gene` | Binary DMR status (q<0.05), `mod_difference`, `dmr_qvalue` |
| K119ub signal | `data/k119ub_gene_signal.tsv` | `symbol` | `gb_ctrl_signal` (WT gene body K119ub continuous signal) |

### Dependencies
- `pROC` (already installed) — for `roc()`, `auc()`, `ci.auc()`
- All other deps come from `_shared_config.R` (tidyverse, ggplot2, patchwork, ggrepel, GenomicRanges, etc.)

---

### Step-by-step

**Step 1 — Load & compute baseline WT 5hmC (Task 3a)**
- Load `Extract_hmc_regional-frac_*.tsv.gz` (same pattern as section_22 lines 105-122)
- Rename cols 4-7 by position → `hmc_ctrl_M`, `hmc_mut_M`, `hmc_ctrl_F`, `hmc_mut_F`
- `wt_hmc = (hmc_ctrl_M + hmc_ctrl_F) / 2`
- Deduplicate by gene, filter complete cases

**Step 2 — Merge with DMR status**
- From `hmc_dmr`: select `gene, dmr_qvalue, mod_difference, significant`
- `inner_join()` on gene → ~20K genes with WT 5hmC + DMR status

**Step 3 — Load K119ub signal**
- Load `data/k119ub_gene_signal.tsv`, select `symbol` → `gene`, `gb_ctrl_signal`
- Filter to finite values
- `inner_join()` with Step 2 → ~18-20K genes for combined model

**Step 4 — Fit 3 logistic regression models (Tasks 3b + 3c)**

| Model | Formula | Purpose |
|-------|---------|---------|
| A | `significant ~ wt_hmc` | Does baseline 5hmC predict DMR? |
| B | `significant ~ gb_ctrl_signal` | Does K119ub predict DMR? |
| C | `significant ~ wt_hmc + gb_ctrl_signal` | Combined — does 5hmC add info beyond K119ub? |

Extract per model: OR (95% CI), AIC, McFadden pseudo-R², ROC AUC (95% CI via pROC).
Also fit standardized versions (z-scored predictors) for comparable coefficient magnitudes in Model C.

**Step 5 — Chromatin state annotation (for 23e)**
- Reuse `classify_chromatin_state()` from `_shared_config.R`
- Load ChIP peaks, compute overlaps, classify each gene (same pattern as sections 10/22)

---

### Figures

**23a — ROC curves (3 models overlaid)**
- x: 1 − Specificity, y: Sensitivity
- 3 colored lines: Model A (blue `#377EB8`), Model B (purple `#756BB1`), Model C (red `#E41A1C`)
- Legend: model name + AUC, annotation box with AIC + McFadden R²
- Diagonal reference line

**23b — Predicted probability curve (Model A)**
- x: Baseline WT 5hmC, y: P(hmC DMR)
- Points colored by actual DMR status (alpha=0.1 for density)
- Sigmoid via `geom_smooth(method="glm", family=binomial)`
- Rug plot on x-axis, annotation with OR + 95% CI + p-value

**23c — Model comparison (bar/dot chart)**
- 3 rows (A, B, C), two panels via patchwork:
  - Left: AIC bars (lower = better)
  - Right: AUC point + CI whiskers
- Color-coded by model

**23d — Dose-response scatter (Task 3d)**
- x: Baseline WT 5hmC, y: hmC `mod_difference`
- Points colored by significance (red vs grey70)
- Linear trend line, Spearman rho annotation
- Horizontal dashed line at y=0
- Top genes labeled with `ggrepel`
- Expected: negative correlation (high baseline → more loss)

**23e — Dose-response by chromatin state (faceted)**
- Same as 23d, `facet_wrap(~chromatin_state)`
- Per-facet Spearman rho annotation
- Shows whether substrate availability varies by chromatin context

### Tables exported to `TABLES_DIR`
1. `baseline_hmc_predictor_all_genes.tsv` — full merged data with predictions
2. `baseline_hmc_model_comparison.tsv` — 3-row summary (AIC, AUC, CI, R², N)
3. `baseline_hmc_model_coefficients.tsv` — detailed coefficient table (raw + standardized)

---

### Key patterns to follow (from section_22)
- `source("scripts/viz_sections/_shared_config.R")` at top
- `stopifnot()` for all input validation — **no fallbacks**
- `cat()` progress logging with section banner
- `save_multiformat_ggplot()` for all figures
- Same `EXTRACT_DIR` path constant as section_22
- Deduplicate genes with `dplyr::distinct(gene, .keep_all = TRUE)`

### Verification
1. Run: `Rscript scripts/viz_sections/section_23_baseline_hmc_predictor.R` from `downstream/`
2. Check outputs in `plots/visualizations/23_baseline_hmc_predictor/`
3. Verify tables in `plots/visualizations/tables/`
4. Confirm Model A AUC > Model B AUC (expected if substrate availability > K119ub)
5. Confirm negative Spearman rho in 23d (high baseline → more hmC loss)

---

All 5 figures (23a-23e) in PDF/SVG/JPG and all 3 tables generated successfully.

**Key results summary:**

| Model | AIC | AUC | McFadden R² |
|---|---|---|---|
| **A: Baseline 5hmC** | 21,918 | **0.762** | 0.1401 |
| B: K119ub signal | 24,885 | 0.592 | 0.0236 |
| **C: Combined** | **20,784** | **0.800** | **0.1846** |

**Findings confirm the substrate availability hypothesis:**
- Baseline 5hmC is **6× more predictive** than K119ub signal (R² 0.14 vs 0.02)
- Combined model best (AUC=0.800), but most information comes from 5hmC
- Standardized Model C coefficients: 5hmC beta=1.182 vs K119ub beta=-0.850
- **Strong negative dose-response** (Spearman rho = -0.586, p < 2.2e-16): higher baseline 5hmC → more hmC loss in BAP1-KO
- Dose-response varies by chromatin state: strongest in Poised_Enhancer (rho=-0.499), Active_Promoter (rho=-0.466), weakest in Repressed_Promoter (rho=0.052)

---



**Language > Specialist**: R/Statistics > Computational Epigeneticist

Based on the results, the **substrate availability hypothesis is strongly supported**:

## Verdict

**Baseline WT 5hmC is the dominant predictor of which genes become hmC DMRs in BAP1-KO**, far exceeding K119ub chromatin signal (AUC 0.762 vs 0.592, R² 6× larger).

## Biological interpretation

1. **TET impairment acts proportionally to substrate abundance.** The strong negative dose-response (rho = -0.586) means genes with more 5hmC to begin with lose more 5hmC in the knockout. This is exactly what you'd expect if BAP1 loss reduces TET recruitment/activity — genes with more TET-dependent 5hmC have more to lose.

2. **K119ub is upstream but not the bottleneck.** K119ub signal adds predictive value in the combined model (AUC 0.762 → 0.800), but it's secondary. This makes mechanistic sense: H2AK119ub mediates TET recruitment, but whether a gene actually *shows* hmC loss depends more on how much 5hmC (TET product) it had in the first place.

3. **Chromatin context modulates the effect.** The dose-response is strongest at Active_Promoters (rho=-0.466), Poised_Enhancers (rho=-0.499), and Other/unmarked regions (rho=-0.529), but essentially absent at Repressed_Promoters (rho=0.052). This suggests Polycomb-silenced promoters are already in a low-5hmC state with minimal TET activity, so there's little substrate to lose — the floor effect.

4. **Revised mechanistic model:** BAP1 loss → increased H2AK119ub → impaired TET recruitment → reduced 5hmC maintenance, where the **magnitude of 5hmC loss scales with how much active demethylation the gene relied on** in the wildtype state. Genes with high baseline TET conversion are the most vulnerable.

This strengthens your revised conclusions from section 18 — K119ub is a weak gene-specific predictor (Cliff's delta 0.089) because it's the *mechanism*, not the *determinant*. The determinant is how much TET-dependent demethylation the gene was undergoing.