# Stripe Analysis Pipeline - Phases 2-4 Context Document

**Status**: Implementation complete
**Date**: December 2024
**Implements**: Differential stripe quantification and analysis per CLAUDE.md specification

---

## Files Created

| Phase | File | Purpose |
|-------|------|---------|
| 2 | `scripts/phase2_quantification.R` | Stripe-to-GInteractions + mariner extraction + aggregation |
| 2 | `scripts/phase2_quantification.sb` | SLURM batch script (array job) |
| 3 | `scripts/phase3_edgeR.R` | Differential analysis with edgeR QLF |
| 3 | `scripts/phase3_edgeR.sb` | SLURM batch script |
| 4 | `scripts/phase4_integration.R` | Merge annotations + classify direction |
| All | `scripts/run_stripe_phases234.sb` | Master script (runs Phases 2-4 sequentially) |

---

## How to Run

### On Expanse HPC

```bash
cd /expanse/lustre/projects/csd940/zalibhai/stripes

# Option 1: Run master script (recommended)
sbatch scripts/run_stripe_phases234.sb

# Option 2: Run phases individually
sbatch scripts/phase2_quantification.sb  # ~30-60 min
sbatch --dependency=afterok:<jobid> scripts/phase3_edgeR.sb  # ~5 min
# Phase 4 runs automatically or manually:
Rscript scripts/phase4_integration.R early
Rscript scripts/phase4_integration.R late
```

### Expected Runtime
- Phase 2: ~30-60 minutes per timepoint (Hi-C extraction)
- Phase 3: ~5-10 minutes per timepoint (edgeR analysis)
- Phase 4: ~5 minutes per timepoint (integration)
- **Total: ~1-2 hours per timepoint**

---

## Output Structure

```
outputs/{timepoint}/res_5kb/
├── 01_unified_stripes.tsv          # Phase 1 output
├── 01_unified_stripes.rds
├── 02_stripe_counts.tsv            # Phase 2: Count matrix
├── 02_stripe_counts.rds
├── 02_extraction_metadata.rds
├── 02_extracted/                    # HDF5 matrices
├── temp_hdf5/                       # Temporary extraction files
├── edgeR_results/                   # Phase 3 outputs
│   ├── 03_dge_object.rds
│   ├── 03_all_results.tsv
│   ├── 03_significant_FDR05.tsv
│   ├── 03_summary.txt
│   └── plots/
│       ├── mds_plot.pdf
│       ├── bcv_plot.pdf
│       ├── volcano_plot.pdf
│       └── ma_plot.pdf
├── 04_final_differential_stripes.tsv   # Phase 4 final output
├── 04_final_differential_stripes.rds
├── 04_final_differential_stripes.bedpe # Juicebox visualization
├── 04_summary.txt
├── 04_stripes_lost.tsv              # Per-direction subsets
├── 04_stripes_gained.tsv
├── 04_stripes_strengthened.tsv
└── 04_stripes_weakened.tsv
```

---

## Key Implementation Details

### Phase 2: Stripe-to-GInteractions Conversion

```r
# Calculate sampling point based on direction
sampling_point <- ifelse(
  stripes$direction_type == "3prime",
  stripes$anchor_center + 100000,  # Downstream
  stripes$anchor_center - 100000   # Upstream
)

# Create GInteractions
anchor1 <- GRanges(chr, anchor_x1:anchor_x2)  # Stripe anchor
anchor2 <- GRanges(chr, sampling_point)        # 100kb into span
gi <- GInteractions(anchor1, anchor2)

# Bin and buffer
binned <- assignToBins(gi, binSize = 5000, pos1 = "center", pos2 = "center")
buffered <- pixelsToMatrices(binned, buffer = 2)  # 5x5 matrix

# Extract with mariner
pixels <- pullHicMatrices(buffered, hicFiles, binSize = 5000, norm = "KR", ...)
```

### Phase 3: edgeR (Skip Filtering)

```r
# Key difference from loop pipeline: SKIP FILTERING
if (config$edger$skip_filtering) {
  # Retain ALL stripes (small feature set ~150-300)
  dge <- calcNormFactors(dge, method = "TMM")
} else {
  keep <- filterByExpr(dge, ...)
  dge <- dge[keep, ]
}

# Standard QL-GLM workflow
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)  # Mutant effect
```

### Phase 4: Direction Classification

Per CLAUDE.md specification:

```r
direction <- case_when(
  source == "control_only" & FDR < 0.05 & logFC < 0 ~ "lost",
  source == "mutant_only" & FDR < 0.05 & logFC > 0 ~ "gained",
  source == "shared" & FDR < 0.05 & logFC > 0.3 ~ "strengthened",
  source == "shared" & FDR < 0.05 & logFC < -0.3 ~ "weakened",
  source == "shared" ~ "unchanged",
  TRUE ~ "ambiguous"
)
```

---

## Validation Checklist

### After Phase 2
- [ ] Count matrix dimensions: ~286 x 6 (early) or ~200 x 6 (late)
- [ ] Within-group correlation > 0.90
- [ ] NA percentage < 5%
- [ ] `02_stripe_counts.rds` exists

### After Phase 3
- [ ] MDS plot shows condition separation
- [ ] BCV in reasonable range (0.2-0.5)
- [ ] Significant stripes at FDR < 0.05 (expect 10-50)
- [ ] `edgeR_results/03_all_results.tsv` exists

### After Phase 4
- [ ] All stripes classified (no NAs in direction)
- [ ] Mix of lost/gained/strengthened/weakened
- [ ] BEDPE loads in Juicebox
- [ ] `04_final_differential_stripes.tsv` exists

---

## Configuration Parameters Used

From `config/stripe_config.yaml`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `quantification.sampling_distance_bp` | 100000 | 100kb into span |
| `quantification.buffer_bins` | 2 | 5x5 matrix extraction |
| `quantification.normalization` | "KR" | Knight-Ruiz preferred |
| `edger.skip_filtering` | true | Retain all stripes |
| `edger.normalization_method` | "TMM" | Trimmed Mean of M |
| `edger.fdr_primary` | 0.05 | Primary significance threshold |
| `classification.logFC_threshold` | 0.3 | For strengthened/weakened |

---

## Troubleshooting

### Phase 2 Fails
1. Check Hi-C files exist at expected paths
2. Verify `readHicNormTypes()` returns KR or VC
3. Check memory (increase to 128G if needed)
4. Review `logs/phase2_*.out` for errors

### Phase 3 Low Power
1. Check within-group correlations from Phase 2
2. Review BCV plot - high BCV reduces power
3. Consider increasing replicate count if possible

### Juicebox Not Loading BEDPE
1. Verify BEDPE has header line
2. Check coordinates are 0-based
3. Ensure chr names match .hic file (chr1 vs 1)

---

## Reference Templates

These scripts adapted patterns from the loop pipeline:

| Template | Source | Key Adaptation |
|----------|--------|----------------|
| Matrix extraction | `scripts/extract_counts.R` | Stripe geometry, sampling point |
| Aggregation | `scripts/aggregate.R` | Same 5x5 sum strategy |
| edgeR analysis | `scripts/edgeR.R` | Skip filtering |
| Config pattern | `stripes/scripts/phase1_detection.R` | YAML loading, path helpers |

---

**Pipeline Complete**: Phases 2-4 are ready for execution on Expanse.
