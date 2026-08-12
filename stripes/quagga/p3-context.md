# Stripe Analysis Pipeline - Results & Interpretation

**Status**: Pipeline complete, results analyzed
**Date**: December 2024
**Verdict**: Detection-based differences are likely noise; recommend focusing on medium/high confidence stripes only

---

## Executive Summary

The differential stripe analysis pipeline successfully ran on both timepoints, but the results reveal concerning patterns:

| Metric | Early | Late |
|--------|-------|------|
| Total stripes | 286 | 200 |
| Lost (control_only) | 126 | 83 |
| Gained (mutant_only) | 86 | 73 |
| Unchanged (shared) | 74 | 44 |
| **Directional consistency** | **45.8%** | **53.2%** |
| FDR < 0.05 | 0 | 0 |
| FDR < 0.10 | 0 | 1 |

**Key finding**: ~50% directional consistency is essentially random, suggesting most detection-based differences are noise rather than true BAP1-KO effects.

---

## How logFC Was Calculated

### Step 1: Stripe-to-Coordinate Conversion

Each stripe was converted to a quantifiable interaction:

```
Original stripe:
  anchor: chr1:50,000,000-50,100,000 (stripe anchor region)
  span:   chr1:50,100,000-50,600,000 (where stripe extends)
  direction: 3prime (downstream)

Sampling point calculation:
  anchor_center = (50,000,000 + 50,100,000) / 2 = 50,050,000
  sampling_point = anchor_center + 100kb = 50,150,000

Result: GInteraction between anchor and sampling point
```

### Step 2: Hi-C Matrix Extraction

For each stripe, a 5x5 pixel matrix was extracted at 5kb resolution:

```
                    sampling_point
                         ↓
              ┌─────────────────────────┐
              │  .   .   .   .   .     │  -10kb
              │  .   .   .   .   .     │  -5kb
anchor_center │  .   .   X   .   .     │  center
              │  .   .   .   .   .     │  +5kb
              │  .   .   .   .   .     │  +10kb
              └─────────────────────────┘

Extracted from 6 .hic files with KR normalization
```

### Step 3: Aggregation to Counts

```r
count[stripe_i, sample_j] = sum(5x5_matrix, na.rm = TRUE)
```

This produces a count matrix: `N_stripes x 6_samples`

### Step 4: edgeR Differential Analysis

1. **TMM normalization**: Accounts for library size differences
2. **Dispersion estimation**: Uses 3 replicates per condition
3. **Quasi-likelihood GLM**: `log(counts) ~ intercept + mutant_effect`
4. **logFC = mutant_effect coefficient**: `log2(mutant / control)`

**Interpretation**:
- logFC > 0: Higher signal in BAP1-KO mutant
- logFC < 0: Higher signal in wildtype control

---

## The Classification Problem & Fix

### Original Problem

The original classification logic required FDR < 0.05:

```r
# ORIGINAL (produced 0 lost/gained, all ambiguous)
direction = case_when(
  source == "control_only" & FDR < 0.05 & logFC < 0 ~ "lost",
  source == "mutant_only" & FDR < 0.05 & logFC > 0 ~ "gained",
  ...
  TRUE ~ "ambiguous"
)
```

Since no stripes reached FDR < 0.05, all condition-specific stripes became "ambiguous".

### The Fix: Tiered Classification

Detection is now PRIMARY evidence; FDR + direction provide confidence tiers:

```r
# FIXED (all stripes classified with confidence tier)
direction = case_when(
  source == "control_only" ~ "lost",     # By detection
  source == "mutant_only" ~ "gained",    # By detection
  source == "shared" & FDR < 0.05 & logFC > 0.3 ~ "strengthened",
  source == "shared" & FDR < 0.05 & logFC < -0.3 ~ "weakened",
  source == "shared" ~ "unchanged"
)

direction_confidence = case_when(
  # High: FDR < 0.1 AND correct direction
  source == "control_only" & FDR < 0.1 & logFC < 0 ~ "high",
  source == "mutant_only" & FDR < 0.1 & logFC > 0 ~ "high",

  # Medium: Correct direction (logFC > 0.2)
  source == "control_only" & logFC < -0.2 ~ "medium",
  source == "mutant_only" & logFC > 0.2 ~ "medium",

  # Low: Detection only
  TRUE ~ "low"
)
```

---

## Results Breakdown

### Confidence Tier Distribution

| Timepoint | High | Medium | Low | Total |
|-----------|------|--------|-----|-------|
| Early - Lost | 0 | 9 | 117 | 126 |
| Early - Gained | 0 | 10 | 76 | 86 |
| Late - Lost | 0 | 10 | 73 | 83 |
| Late - Gained | 1 | 6 | 66 | 73 |

**~90% of stripes are low confidence** (detection only, no quantitative support).

### Directional Consistency

"Directional consistency" measures whether the logFC direction matches the detection-based classification:

- **Lost stripes** should have logFC < 0 (lower in mutant)
- **Gained stripes** should have logFC > 0 (higher in mutant)

| Timepoint | Consistent | Inconsistent | % Consistent |
|-----------|------------|--------------|--------------|
| Early | 97 | 115 | **45.8%** |
| Late | 83 | 73 | **53.2%** |

**This is essentially random (50%)**. If stripes were truly lost/gained, we'd expect >90% consistency.

### Fold Change Statistics

| Metric | Early | Late |
|--------|-------|------|
| Median logFC | 0.003 | -0.009 |
| Range | -0.568 to +1.000 | -0.624 to +0.437 |
| FDR < 0.05 | 0 | 0 |
| FDR < 0.10 | 0 | 1 |

The median logFC near zero indicates no systematic shift between conditions.

---

## Interpretation

### What the Data Shows

1. **Detection differences exist**: Many stripes detected in only one condition (126+86=212 early, 83+73=156 late)

2. **Quantitative support is absent**:
   - No stripes reach FDR < 0.05
   - Median logFC ≈ 0
   - Directional consistency ≈ 50%

3. **This is a red flag**: The detection-based classification doesn't hold up to quantitative scrutiny

### Possible Explanations

1. **Quagga detection is noisy at 200kb min_length**
   - Stripes near detection threshold appear stochastically
   - Pooling 3 replicates to merged files creates artificial sensitivity
   - Individual replicates show only 11-20 stripes each

2. **100kb sampling distance may be suboptimal**
   - May not capture the actual stripe signal
   - Stripe intensity might be concentrated closer to the anchor

3. **Biological effect is genuinely weak**
   - BAP1-KO may not strongly affect chromatin stripes
   - Effect may be masked by tissue/cell heterogeneity

4. **Sequencing depth is borderline**
   - ~300M reads per replicate (Quagga paper recommends 500M+)
   - Lower sensitivity for stripe detection

---

## Recommendations

### For Publication/Interpretation

1. **Focus on medium/high confidence stripes only**:
   - Early: 19 stripes (9 lost + 10 gained with directional support)
   - Late: 17 stripes (10 lost + 1 high + 6 medium gained)

2. **Be cautious about "lost/gained" claims**:
   - Detection-based classification lacks quantitative support
   - ~50% directional consistency suggests noise

3. **Report the null finding honestly**:
   - "No statistically significant differential stripes at FDR < 0.05"
   - "Detection-based differences do not show quantitative support"

### For Validation

1. **Manual inspection in Juicebox**:
   - Load the BEDPE files with .hic files
   - Visually compare top candidates between ctrl and mut
   - Check if stripes are genuinely absent vs. just below detection

2. **Cross-reference with CTCF Cut&Run** (when available):
   - CTCF-stripes should have CTCF peaks at anchors
   - Validates that detected stripes are real

3. **Consider alternative quantification**:
   - Try different sampling distances (50kb, 150kb)
   - Integrate signal along entire stripe body
   - Use APA-style aggregate analysis

### For Future Studies

1. **Increase sequencing depth** to 500M+ reads per replicate
2. **Consider lower min_length** if EP-stripes are of interest
3. **Use more stringent detection thresholds** for merged files

---

## Output Files

### Per Timepoint

```
outputs/{timepoint}/res_5kb/
├── 01_unified_stripes.tsv          # Phase 1: All detected stripes
├── 02_stripe_counts.tsv            # Phase 2: Count matrix
├── edgeR_results/
│   ├── 03_all_results.tsv          # Phase 3: edgeR statistics
│   ├── 03_dge_object.rds           # DGEList object
│   └── plots/                      # MDS, BCV, volcano plots
├── 04_final_differential_stripes.tsv   # Phase 4: Final results
├── 04_final_differential_stripes.bedpe # For Juicebox
├── 04_stripes_lost.tsv             # Lost stripes only
├── 04_stripes_gained.tsv           # Gained stripes only
└── 04_summary.txt                  # Summary statistics
```

### Key Columns in Final Output

| Column | Description |
|--------|-------------|
| `stripe_id` | Unique identifier |
| `source` | control_only, mutant_only, shared |
| `direction` | lost, gained, unchanged, strengthened, weakened |
| `direction_confidence` | high, medium, low |
| `direction_consistent` | TRUE/FALSE/NA |
| `logFC` | log2 fold change (mutant/control) |
| `FDR` | Benjamini-Hochberg adjusted p-value |
| `confidence` | Phase 1 detection confidence |

---

## Technical Details

### Parameters Used

| Parameter | Value | Description |
|-----------|-------|-------------|
| Resolution | 5kb | Hi-C bin size |
| Sampling distance | 100kb | Distance into span for quantification |
| Buffer | 2 bins | 5x5 matrix extraction |
| Normalization (Hi-C) | KR | Knight-Ruiz balancing |
| Normalization (edgeR) | TMM | Trimmed Mean of M-values |
| FDR threshold | 0.05 | For shared stripe classification |
| High confidence FDR | 0.10 | For condition-specific confidence |
| Medium confidence logFC | 0.2 | For directional support |

### Scripts Modified for Tiered Classification

| Script | Changes |
|--------|---------|
| `scripts/phase4_integration.R` | Tiered classification logic |
| `scripts/phase4_integration.sb` | New SLURM script |
| `config/stripe_config.yaml` | Tiered threshold parameters |
| `CLAUDE.md` | Updated specification |

---

## Conclusion

The stripe analysis pipeline is technically sound, but the biological results are concerning:

1. **The pipeline works**: All phases complete successfully
2. **The classification fix works**: No more ambiguous stripes
3. **The underlying signal is weak**: ~50% directional consistency, 0 significant stripes

**Bottom line**: Most detection-based stripe differences between BAP1-KO and wildtype appear to be stochastic noise rather than true biological effects. Focus validation efforts on the ~17-19 medium/high confidence stripes per timepoint.

---

*Document generated: December 2024*
*Project: BAP1-KO Differential Chromatin Stripe Analysis*
*Organism: Mouse (mm10)*
