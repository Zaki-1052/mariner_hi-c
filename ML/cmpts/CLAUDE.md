# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

CALDER2 + SNIPER subcompartment analysis for BAP1-KO mouse cerebellum (mm10). Decomposes coarse HOMER A/B compartment calls into subcompartments (A.1, A.2, B.1, B.2) to distinguish true compartment flips from quantitative within-compartment weakening. CALDER2 is the primary tool (reference-guided, works natively on mm10); SNIPER is secondary (neural network, trained on hg19, requires retraining for mm10).

Two parallel tracks:
- **Track A (CALDER2):** A1 calling → A2 differential → A3 epigenomic validation → A4 HOMER integration
- **Track B (SNIPER):** B1 crop map generation → B2+ training/application (in progress)

## Running the Pipeline

### Track A — Full pipeline (HPC, from login node)

```bash
cd /expanse/.../mariner_hi-c/ML/cmpts
bash scripts/run_full_calder2.sh                        # both timepoints, A1-A4
bash scripts/run_full_calder2.sh --tp early             # early (250831) only
bash scripts/run_full_calder2.sh --tp early --start A2  # skip A1
```

Stages chain via `--dependency=afterok`. A3 and A4 run in parallel after A2.

### Track A — Individual stages

```bash
# A1: CALDER2 subcompartment calling (per sample)
bash scripts/A1_submit_calder2.sh
# Or: sbatch scripts/A1_run_calder2.sb <timepoint> <sample>

# A2: Differential analysis (per timepoint)
sbatch scripts/A2_run.sb <timepoint>

# A3: Epigenomic validation (both TPs in one job)
sbatch scripts/A3_run.sb

# A4: HOMER compartment integration (both TPs in one job)
sbatch scripts/A4_run.sb
```

### Track B — SNIPER (secondary)

```bash
# B1: Generate mm10 crop map from ctrl_merged.hic
sbatch scripts/B1_generate_cropmap.sb <timepoint>
```

### Prerequisites

- A0 env setup: `bash scripts/A0_setup_calder2_env.sh` (one-time, creates `calder2_env` conda env)
- SNIPER env: `sniper_env` with PyTorch + scipy (created separately)
- CALDER2 R package deps: `Rscript scripts/utils/install_calder2_deps.R`

## Two-Root Path Convention

| Root | Path | Contents |
|------|------|----------|
| CODE_DIR | `/expanse/.../mariner_hi-c/ML/cmpts` | Scripts, config, small TSVs, figures (GitHub-synced) |
| DATA_DIR | `/expanse/.../sniper` | Large intermediates: CALDER2 .Rds files, SNIPER dumps (NOT synced) |

A1 writes large intermediates to DATA_DIR, then copies small result files (<1MB) into `outputs/` under CODE_DIR for local rsync. All scripts take `<data_root>` and `<code_root>` as CLI args.

## Conda Environments

| Env | Used by | Key packages |
|-----|---------|--------------|
| `calder2_env` | Track A (A1-A4) | R + CALDER2, data.table, rtracklayer, pheatmap |
| `sniper_env` | Track B | Python 3.x, PyTorch, scipy, numpy, pandas |

## Output Structure

```
outputs/
├── calder2/
│   ├── early/         # 250831 results
│   │   ├── 250831/    # raw CALDER2 output per sample
│   │   ├── 250831_subcompartment_labels_100kb.tsv
│   │   ├── 250831_transition_matrix.tsv
│   │   └── 250831_enrichment_matrix.tsv
│   └── late/          # 250402 results (same structure)
├── integration/
│   ├── early/         # HOMER × CALDER2 cross-tabulations
│   ├── late/
│   └── combined_weakening_decomposition.tsv
└── sniper/            # (Track B outputs, when complete)
```

## Key Constants

- **Resolution:** CALDER2 runs at 50kb, reports at 100kb. SNIPER runs at 100kb.
- **Subcompartment labels:** `A.1`, `A.2`, `B.1`, `B.2` (ordered). No B.3 in mouse.
- **Chromosomes:** chr1-chr19 autosomes only (X, Y, M excluded).
- **Timepoints:** `250402` = late/adult, `250831` = early/P12.
- **Samples:** `ctrl_merged`, `mut_merged` (merged biological replicates).

## Configuration

`config/sniper_config.yaml` centralizes all HPC paths, timepoint mappings, BigWig filenames, and SLURM resource settings. Scripts currently hard-code paths in their constants sections rather than reading this YAML — the config exists as a reference/future migration target.

## Conventions

- **A2 bins via plurality vote:** CALDER2 produces variable-width segments. A2 bins them to uniform 100kb by assigning the most-covered subcompartment label per bin.
- **Multi-format figures:** R scripts use `scripts/utils/multi_format_output.R` (`save_multiformat_ggplot()`) to save PNG + PDF + SVG + JPG in per-figure subdirectories.
- **SNIPER is hg19-only by design.** The `repos/SNIPER/` code is the upstream repo with added `utilities/mm10_config.py`, `utilities/data_processing_mm10.py`, and `utilities/interchromosome_matrix_mm10.py` for mm10 adaptation. B1 generates a mm10 crop map; further stages require retraining the neural network.
- **CALDER2 reference correlation:** All chromosomes must pass >0.59 correlation with the mm10 reference. Early timepoint systematically correlates better (~0.77 vs ~0.69).
- **Script CLI pattern:** All R scripts take positional args `<timepoint> <data_root> <code_root>` (A1, A2) or `<data_root> <code_root>` (A3, A4 process both TPs). SLURM `.sb` wrappers parse the timepoint and pass the roots.
