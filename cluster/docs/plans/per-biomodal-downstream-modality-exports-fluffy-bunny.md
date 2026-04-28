# Plan: SLURM Script for UCSC CpG Island Quantification

## Context

Challana asked Zakir to run quants over `sorted_CpGmm10.bed` (16,009 UCSC CpG islands for mm10) — this is the full UCSC set (~2x larger than the gencode vM25 annotation used in run-5's 8,910 islands). The goal is feature extraction only (count, mean, regional-frac) against the CG zarr store from run-5 (8 samples, deep-seq).

## Critical Format Issue

The new BED file uses `chr`-prefixed chromosomes (`chr1`, `chr2`, ...) while the zarr store indexes on bare numbers (`1`, `2`, ...). The script must strip the `chr` prefix and add the header + Annotation column that modality expects.

## Implementation

Create one file: `downstream/modality/run_cpg_ucsc_quants.sb`

### SLURM Parameters
- Partition: `shared`
- CPUs: 8 (quants alone don't need 16)
- Memory: 64GB (conservative for 16k regions × 8 samples)
- Time: 2:00:00 (expected runtime ~15 min, 8x buffer)
- Account: `csd940`
- Log: `logs/cpg_ucsc_%j.out`

### Script Steps

**Step 0 — Preprocess BED:**
- `awk` to strip `chr` prefix, add header (`Chromosome\tStart\tEnd\tAnnotation`), add `CpG_island` annotation column
- Write to `sorted_CpGmm10.annotation.bed` in the HPC modality working dir
- Idempotent (skip if already exists)
- Validate line count (16,010 = 16,009 data + 1 header) and chromosome format

**Step 1 — Count (total_c):**
```bash
modality -v get count \
    --zarr-path "$ZARR" \
    --sample-sheet "$METADATA" \
    --bedfile "$PROCESSED_BED" \
    --fields total_c \
    --output-dir "$REGION_RESULTS_DIR"
```

**Step 2 — Mean (total_c):**
```bash
modality -v get mean \
    --zarr-path "$ZARR" \
    --sample-sheet "$METADATA" \
    --bedfile "$PROCESSED_BED" \
    --fields total_c \
    --order-by-group condition \
    --output-dir "$REGION_RESULTS_DIR"
```

**Step 3 — Regional fraction (mc + hmc, min-coverage 10):**
```bash
modality -v get regional-frac \
    --zarr-path "$ZARR" \
    --sample-sheet "$METADATA" \
    --bedfile "$PROCESSED_BED" \
    --fields mc hmc \
    --min-coverage 10 \
    --order-by-group condition \
    --output-dir "$REGION_RESULTS_DIR"
```

### Paths

| Variable | Value |
|----------|-------|
| WORK_DIR | `/expanse/lustre/projects/csd940/zalibhai/biomodal/modality` |
| REPO_DIR | `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/biomodal/downstream/modality` |
| ZARR | `/expanse/lustre/projects/csd940/zalibhai/biomodal/evoC-run/output/duet-1.5.0_evoC_Bap1_run_6bp/sample_outputs/zarr_store/CG/evoC_Bap1_run.genome.GRCm38_primary_assembly.dedup.CG.zarrz` |
| METADATA | `${WORK_DIR}/metadata.tsv` |
| RAW_BED | `${REPO_DIR}/sorted_CpGmm10.bed` |
| PROCESSED_BED | `${WORK_DIR}/sorted_CpGmm10.annotation.bed` |
| Output dir | `${WORK_DIR}/outputs/cpg_ucsc_CG/Results/sorted_CpGmm10.annotation/` |

### Style Notes
- No `set -euo pipefail` (per user preference)
- Progress logging with timestamps
- Input validation (file existence checks)
- No DMR calling — quants only
- Conda env: `modality`
- `cd` to WORK_DIR (matches run_modality.sb pattern)

## Verification

After submitting on HPC:
1. Check `logs/cpg_ucsc_*.out` for completion message
2. Confirm 3 `Extract_*` directories in `outputs/cpg_ucsc_CG/Results/sorted_CpGmm10.annotation/`
3. Verify `.tsv.gz` files exist and have 16,009 data rows
4. Spot-check that chromosome column uses bare numbers (not `chr`)

## Files to Create
- `downstream/modality/run_cpg_ucsc_quants.sb` (new)
