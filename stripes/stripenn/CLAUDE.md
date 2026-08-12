# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@README.md

## Project Overview

Stripenn-based differential chromatin stripe analysis pipeline for BAP1-KO vs wildtype mouse cerebellum (mm10). Companion to the Quagga stripe pipeline (`stripes/quagga/`) and the main mariner loop pipeline (repo root). Uses Stripenn's Canny edge detection + pixel saturation instead of Quagga's Gabor filtering, with a built-in `score` command that replaces mariner-based quantification.

**HPC:** SDSC Expanse (SLURM). **Genome:** mm10. **Conditions:** BAP1-KO mutant vs wildtype control, n=3 replicates per condition. **Timepoints:** 250831 (early/P12), 250402 (late/adult). **Resolutions:** 5kb and 10kb.

## Running the Pipeline

### Full pipeline (Stages 0-7, from HPC login node)

```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn/scripts
bash run_full_stripenn.sh              # all stages
bash run_full_stripenn.sh --skip-stage0  # skip .hic->.mcool conversion
```

### Individual stages

```bash
# Stage 0: .hic -> .mcool (uses hictk conda env, not mariner_env)
bash submit_00_convert.sh

# Stages 1-7 use mariner_env
bash submit_01_call.sh        # stripe detection on merged samples
bash submit_02_union.sh       # build union set across conditions
bash submit_03_score.sh       # quantify replicates against union
bash submit_04_edgeR.sh       # differential analysis
bash submit_05_integration.sh # classify direction + confidence
bash submit_06_compare.sh     # cross-resolution comparison
sbatch scripts/stripenn_visualizations.sb  # visualization + annotation
```

### Local-only scripts (no HPC needed)

```bash
python scripts/generate_bedpe.py  # simple 15-column BEDPEs for JuiceBox
```

### Dependency chaining

Submit wrappers accept `--dependency=afterok:<jobid>` as first arg. `run_full_stripenn.sh` wires this automatically. If a job fails, all downstream jobs get `DependencyNeverSatisfied` — cancel stale jobs with `scancel -u $USER --state=PENDING`, fix, and resubmit from the failed stage.

## Two-Root Path Convention

| Root | Path | Contents |
|------|------|----------|
| CODE_DIR | `/expanse/.../mariner_hi-c/stripes/stripenn` | Scripts, config, small TSVs (GitHub-synced) |
| DATA_DIR | `/expanse/.../stripes/stripenn` | mcool files, outputs, logs (NOT synced) |

All scripts use absolute paths referencing these two roots. The CODE_DIR is this repo directory; DATA_DIR holds large data.

## Two Conda Environments

| Env | Stages | Why |
|-----|--------|-----|
| `hictk` | Stage 0 only | hictk for .hic v9 -> .mcool conversion |
| `mariner_env` | Stages 1-7 | Python 3.13 + stripenn (patched) + R packages |

## Critical: Stripenn Site-Packages Patches

Two patches to `getStripe.py` in `mariner_env`'s site-packages fix NumPy 2.0 unsigned integer overflow on Python 3.13. **If stripenn is ever reinstalled, patches must be reapplied:**

```bash
conda activate mariner_env
SITE_PKG=$(python -c "import stripenn; import os; print(os.path.dirname(stripenn.__file__))")
sed -i 's/        self.resol = resol/        self.resol = int(resol)/' "${SITE_PKG}/getStripe.py"
sed -i 's/chrsize = self.chromnames2sizes\[chr\]/chrsize = int(self.chromnames2sizes[chr])/' "${SITE_PKG}/getStripe.py"
```

Without these patches, `stripenn compute` crashes with `ValueError: End coordinate less than start` due to unsigned integer underflow in coordinate arithmetic.

## Architecture

8 stages, each with a SLURM `.sb` script, an optional `.R` script, and a `submit_*.sh` wrapper.

```
Stage 0  00_convert_hic_to_cool.sb   .hic -> .mcool (hictk env)       16 jobs
Stage 1  01_call_stripes.sb          stripenn compute on merged        8 jobs
Stage 2  02_build_union.R/.sb        anchor matching + union set       4 jobs
Stage 3  03_score_replicates.sb      stripenn score per replicate      24 jobs
Stage 4  04_edgeR.R/.sb              TMM + robust QL-GLM              4 jobs
Stage 5  05_integration.R/.sb        direction + confidence tiers      4 jobs
Stage 6  06_compare_resolutions.R/.sb  5kb vs 10kb concordance        2 jobs
Stage 7  stripenn_visualizations.R/.sb  plots + ChIP annotation        1 job
```

Each `submit_*.sh` iterates over the parameter space (timepoints x resolutions x samples) and emits job IDs for dependency chaining.

### Key data flow

1. `.hic` -> `.mcool` (hictk convert, Stage 0)
2. `stripenn compute` on merged ctrl/mut mcools -> `result_filtered.tsv` per condition (Stage 1)
3. Merge ctrl + mut calls with 50kb anchor tolerance -> `union_stripes.tsv` + `.bedpe` (Stage 2)
4. `stripenn score` each replicate against union BEDPE -> `O_Sum_added` per stripe per sample (Stage 3)
5. `O_Sum_added` count matrix -> edgeR QL-GLM -> `all_results.tsv` with logFC/FDR (Stage 4)
6. Classify: control_only->"lost", mutant_only->"gained", shared->FDR direction (Stage 5)
7. Cross-resolution overlap and concordance (Stage 6)
8. Volcano plots, ChIP-seq annotation, GO/KEGG enrichment, annotated BEDPEs (Stage 7)

## Configuration

`config/stripenn_config.yaml` — single source of truth for paths, stripenn parameters, edgeR thresholds, ChIP-seq peak files, and SLURM resources. All R and SLURM scripts load this file.

Key parameters:
- `stripenn_compute.seed: 42` (not stripenn's default 123456789)
- `edger.skip_filtering: true` (feature sets are small: ~1,800-7,400)
- `detection.anchor_tolerance_bp: 50000`
- `filtering.exclude_chromosomes: [chrX, chrY, chrM]` (applied post-hoc in R, not by stripenn)

## Stripenn CLI Gotchas

- Output is `result_filtered.tsv` (NOT `.txt` despite upstream README)
- `stripenn score` has NO `-h`/`--header` flag — it auto-detects headers
- `--seed 42` must be explicit; stripenn default is 123456789
- `yes Y |` prefix required for `stripenn compute` to handle interactive overwrite prompt
- Stage 1 appears stuck after chromosome progress bars — this is normal. The `pvalue()` function runs silently for 1-3+ hours per maxpixel percentile. Total 5kb runtime can be 8-20+ hours.
- `export PYTHONWARNINGS="ignore::FutureWarning"` suppresses cosmetic pandas deprecation warnings

## 250402 KR Normalization Issue

The 250402 merged .hic files lack native KR normalization. ICE weights were computed via `cooler balance --name KR` on the two affected mcool files (`ctrl_merged.mcool`, `mut_merged.mcool`). All other 14 mcool files have native KR from hictk conversion. See README.md Section 5.2 for full history.

## Differences from Quagga Pipeline

| Aspect | Quagga | Stripenn |
|--------|--------|----------|
| Input format | .hic directly | .mcool (converted via hictk) |
| Quantification | mariner `pullHicMatrices()` 5x5 at 100kb | `stripenn score` -> `O_Sum_added` |
| Bioconductor deps | mariner, InteractionSet, HDF5Array, strawr | None for detection/scoring |
| Replicate handling | Call per-replicate for confidence | Score union against replicates |
| Parallelism | Array jobs for timepoints | Individual sbatch per (tp, res, sample) with `--dependency=afterok` |

R scripts in Stages 2/4/5/6/7 were ported from Quagga equivalents (see README.md Section 10 for mapping).

## Current State

All Stages 0-7 complete. Stage 7 ran locally 2026-04-21 producing all visualization, annotation, and enrichment outputs. Results documents fully updated at `outputs/{tp}/{tp}_results.md`.

Key results:
- 250402 (late): 7,371 union stripes at 5kb, 2,320 significant (31.5%), more gained than lost
- 250831 (early): 4,008 union stripes at 5kb, 96 significant (2.4%), weak signal
- All effect sizes minimal (max |logFC| = 0.389), BCV very low (0.011-0.021)
- Cross-res logFC correlation: r=0.850 (late), r=0.808 (early)

## Key Files

- `README.md` — full implementation history, bugs, fixes, verification, directory trees
- `config/stripenn_config.yaml` — all parameters and paths
- `outputs/{tp}/{tp}_results.md` — per-timepoint results with top stripes and JuiceBox coordinates
- `docs/stripenn-analysis-prompt.md` — template for analyzing BEDPE files
- `repo/` — upstream stripenn source (reference only, not used at runtime)
- `src/` — empty (placeholder)
