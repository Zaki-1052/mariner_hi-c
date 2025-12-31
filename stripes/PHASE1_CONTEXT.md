# Stripe Analysis Pipeline - Phase 1 Context Document

**Status**: Phase 1 implementation complete and ready to run on Expanse
**Date**: December 2024
**Next Step**: Test Phase 1 on both timepoints, then implement Phase 2

---

## What Has Been Implemented

### Files Created

1. **`stripes/config/stripe_config.yaml`** - Configuration file with all paths and parameters
2. **`stripes/scripts/phase1_detection.R`** - Phase 1 detection script (modular, uses config)
3. **`stripes/scripts/phase1_detection.sb`** - SLURM batch script for running Phase 1

### Directory Structure

```
/expanse/lustre/projects/csd940/zalibhai/stripes/
├── config/
│   └── stripe_config.yaml          # Configuration with paths
├── scripts/
│   ├── phase1_detection.R          # Phase 1 R script
│   └── phase1_detection.sb         # SLURM batch script
├── stripe-data/                    # Input BEDPE files (already exists)
│   ├── 5kb/
│   │   ├── 250402/ (late)
│   │   └── 250831/ (early)
│   └── 10kb/
│       ├── 250402/
│       └── 250831/
├── StripeCaller/                   # Hi-C files (already exists)
│   └── data/hic/
│       ├── 250402/
│       └── 250831/
└── outputs/                        # Will be created by scripts
    ├── early/res_5kb/
    └── late/res_5kb/
```

---

## Path Structure in Config

All paths are built dynamically from config values using helper functions:

**Stripe BEDPE files**:
```
{stripe_data.base_path}/{resolution}kb/{timepoint}/{sample}_stripes.bedpe
→ /expanse/.../stripe-data/5kb/250831/ctrl_merged_stripes.bedpe
```

**Hi-C files** (for Phase 2):
```
{hic_files.base_path}/{timepoint}/{sample}.hic
→ /expanse/.../StripeCaller/data/hic/250831/ctrl_M1.hic
```

**Outputs**:
```
{outputs.base_dir}/{timepoint}/res_{resolution}kb/
→ /expanse/.../stripes/outputs/early/res_5kb/
```

---

## How to Run Phase 1

### On Expanse HPC

```bash
cd /expanse/lustre/projects/csd940/zalibhai/stripes

# Option 1: Run both timepoints in parallel (SLURM array)
sbatch scripts/phase1_detection.sb

# Option 2: Run individual timepoints
Rscript scripts/phase1_detection.R early
Rscript scripts/phase1_detection.R late
```

### Expected Runtime
- ~30-60 minutes per timepoint
- Minimal memory usage (~8GB)

### Expected Output

For each timepoint, Phase 1 creates:

```
outputs/{timepoint}/res_5kb/
├── 01_unified_stripes.rds  # R data structure
└── 01_unified_stripes.tsv  # Human-readable table
```

**Columns in output**:
- `stripe_id`: Unique identifier (stripe_0001, stripe_0002, ...)
- `chr, anchor_x1, anchor_x2`: Anchor coordinates
- `span_y1, span_y2`: Span coordinates
- `direction_type`: "3prime" or "5prime"
- `source`: "control_only", "mutant_only", "shared"
- `pval_ctrl, pval_mut`: Quagga detection p-values
- `n_ctrl_reps, n_mut_reps`: Replicate support counts (0-3)
- `in_10kb`: Boolean - validated at 10kb resolution
- `confidence`: "high", "medium", "low"
- `anchor_center`: Center of anchor (for Phase 2 quantification)

---

## Expected Stripe Counts

Based on data inspection:

| Timepoint | Ctrl Merged | Mut Merged | Expected Union |
|-----------|-------------|------------|----------------|
| **Early (250831)** | 201 | 161 | ~250-300 |
| **Late (250402)** | 128 | 118 | ~150-200 |

Union includes:
- Shared stripes (same anchor + direction in both)
- Control-only stripes
- Mutant-only stripes

---

## Key Implementation Details

### 1. Direction-Aware Anchor Matching

**Critical**: Stripes with same anchor but different direction are treated as **SEPARATE stripes** (not shared).

```r
# 3' stripe: span downstream (y1 >= x2)
# 5' stripe: span upstream (y2 <= x1)

# Matching logic:
same_direction = gr_ctrl$direction == gr_mut$direction
shared_stripes = matches[same_direction, ]  # Only same direction = shared
```

### 2. Confidence Scoring

From PRD specifications:

```r
confidence = case_when(
  n_reps >= 2 & in_10kb ~ "high",
  n_reps >= 1 | in_10kb | pval < 1e-10 ~ "medium",
  TRUE ~ "low"
)
```

Where `n_reps` is:
- For control_only: `n_ctrl_reps`
- For mutant_only: `n_mut_reps`
- For shared: `max(n_ctrl_reps, n_mut_reps)`

### 3. Modular Path Construction

All paths built from config using helper functions:

```r
get_stripe_bedpe_path(resolution_kb, timepoint, sample, config)
get_output_dir(timepoint, resolution_kb, config)
```

No hardcoded paths in the script.

---

## Validation Checklist for Phase 1 Output

After running, validate that:

- [ ] Output files exist for both timepoints
- [ ] Union set contains expected number of stripes (~150-300 per timepoint)
- [ ] Source classification looks reasonable:
  - Should have mix of control_only, mutant_only, and shared
  - Shared should be largest category (typically 40-60%)
- [ ] Direction inference worked:
  - Mix of 3prime and 5prime stripes
  - No "ambiguous" in direction_type column
- [ ] Confidence distribution reasonable:
  - Mix of high/medium/low
  - High confidence stripes should have n_reps >= 2 AND in_10kb = TRUE
- [ ] No stripes incorrectly classified as "shared" when directions differ

---

## Next Steps After Phase 1 Validation

Once Phase 1 outputs look good:

1. **Implement Phase 2: Quantification**
   - Convert stripes to GInteractions (anchor + 100kb sampling point)
   - Use mariner `pullHicMatrices()` to extract 5×5 matrices
   - Sum matrices to count matrix

2. **Implement Phase 3: edgeR**
   - TMM normalization
   - GLM with robust dispersion
   - Skip filtering (small feature set)

3. **Implement Phase 4: Integration**
   - Merge Phase 1 annotations + Phase 3 statistics
   - Classify direction: lost/gained/strengthened/weakened

---

## Files to Reference for Next Phases

**Templates from loop pipeline**:
- `scripts/extract_counts.R` - For Phase 2 mariner extraction pattern
- `scripts/aggregate.R` - For Phase 2 matrix summation
- `scripts/edgeR.R` - For Phase 3 differential analysis
- `scripts/downstream_analysis.R` - For Phase 4 integration

**Key adaptation for Phase 2**: Stripe→GInteractions conversion
```r
# Calculate sampling point based on direction
sampling_point = ifelse(
  direction == "3prime",
  anchor_center + 100000,  # Downstream
  anchor_center - 100000   # Upstream
)

# Create GInteractions
anchor1 = GRanges(chr, anchor_x1:anchor_x2)  # Stripe anchor
anchor2 = GRanges(chr, sampling_point)        # Sampling point
gi = GInteractions(anchor1, anchor2)
```

---

## Configuration Reference

Key parameters in `stripe_config.yaml`:

```yaml
detection:
  anchor_tolerance_bp: 50000  # 50kb matching tolerance

quantification:
  sampling_distance_bp: 100000  # 100kb into span
  buffer_bins: 2                # 5×5 matrix extraction

edger:
  skip_filtering: true          # Keep all stripes
  normalization_method: "TMM"
  fdr_primary: 0.05
```

---

## Important Notes

1. **Run from project base**: Scripts expect to be run from `/expanse/.../stripes/` where `config/` exists
2. **Phase 1 is independent**: Does not need .hic files (only BEDPE)
3. **Output location**: All outputs go to `/expanse/.../stripes/outputs/` not local
4. **Array job**: SLURM script runs both timepoints in parallel

---

## Questions/Issues

If Phase 1 fails or produces unexpected results:

1. Check SLURM log: `logs/phase1_<jobid>_<taskid>.out`
2. Verify BEDPE files exist at expected paths
3. Check config file is at `config/stripe_config.yaml`
4. Ensure running from correct directory

---

**Ready to run**: Phase 1 is complete and ready for testing on Expanse.
