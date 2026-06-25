# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Biomodal DUET evoC differential methylation analysis (5mC and 5hmC) in BAP1-KO mutant vs wildtype control mouse cerebellum (mm10/GRCm38). Uses biomodal's dual-epigenetic sequencing at 6bp resolution. The downstream analysis has grown into an 85-section R visualization pipeline integrating methylation with Hi-C loops, CUT&RUN marks, ATAC-seq, RNA-seq, and chromatin state annotations.

**Samples:** 8 total (4 control + 4 mutant, 2 batches, 1 per sex per condition)
**Primary run:** run-5 (deep-seq, 8 samples, sex covariate included)
**Key biological finding:** 92.3% of co-significant genes show coordinated mC↑/hmC↓, indicating a TET-mediated demethylation block upon BAP1 loss

## Module Layout

| Directory | Purpose | Languages |
|-----------|---------|-----------|
| `upstream/` | Biomodal DUET v1.5.0 Nextflow pipeline — FASTQ → Zarr stores | Nextflow, Bash |
| `downstream/` | **Primary working area.** Modality XPLR DMR calling + 85-section R visualization pipeline | R, Bash |
| `downstream/modality/` | Modality XPLR CLI configs, SLURM scripts, and run outputs (run-1 through run-5) | Bash |
| `methylseq/` | nf-core/methylseq pipeline (vendored fork for EM-seq/WGBS comparison) | Nextflow |
| `non_cg_analysis/` | Non-CG methylation comparisons with Luo et al. 2017 neuronal atlas | Python |
| `wgbs/` | WGBS validation experiment (samplesheet + Nextflow config) | Nextflow, Bash |
| `docs/` | Literature reviews, QC reports, analysis notes | Markdown |

Each subdirectory with substantial code has its own `CLAUDE.md` — **read `downstream/CLAUDE.md` before working in `downstream/`**, it documents the visualization architecture in detail.

## Running the Pipeline

### Upstream (DUET — HPC only)

```bash
# conda activate env_nf
sbatch upstream/evoc_run.sb
```

~48hr runtime. Produces Zarr stores in `upstream/duet-1.5.0_evoC_Bap1_run_6bp/sample_outputs/zarr_store/`.

### Downstream DMR Calling (HPC only)

```bash
# conda activate modality
sbatch downstream/modality/scripts/run_modality.sb CG    # Primary CpG
sbatch downstream/modality/scripts/run_modality.sb CHG   # Optional
sbatch downstream/modality/scripts/run_modality.sb CHH   # Optional
```

Copies `config_${CONTEXT}.txt` → `config.txt`, then runs `core-workflow-v1.3.sh` (4 steps: BioQC → Feature Extraction → DMR Calling → DMR Visualization).

### Visualization Pipeline (local R)

```bash
# Must run from downstream/ directory
cd downstream/

# All 85 sections sequentially
bash scripts/viz_sections/run_all_sections.sh

# Individual section
Rscript scripts/viz_sections/section_04_volcano_plots.R

# Batch runners for grouped sections
bash scripts/viz_sections/run_sections_22_26.sh   # demethylation ratio + TET
bash scripts/viz_sections/run_sections_61_63.sh   # stoichiometry + heatmap
bash scripts/viz_sections/run_sections_74_78.sh   # neuronal + aging

# Permutation tests (HPC — heavy compute)
sbatch scripts/viz_sections/run_permutation.sb
```

### Updating for a New Modality Run

Follow `downstream/docs/UPDATE.md`. Four files have hardcoded run paths:
1. `downstream/scripts/viz_sections/_shared_config.R` — `DATA_PATHS`, `EXTRACT_PATHS`, `CHG_DATA_PATHS`
2. `downstream/modality/scripts/process_all_dmr.sh` — `INPUT_DIR`
3. `downstream/modality/scripts/visualize_dmr_results.R` — `load_dmr_data()`
4. `downstream/scripts/viz_sections/compare_shallow_vs_deep.R` — `SHALLOW_PATHS`

## Architecture

### Central Config (`downstream/scripts/viz_sections/_shared_config.R`)

~620-line file sourced by every section script. Contains:
- **Path blocks** — Hardcoded paths to run-5 DMR BED files, per-sample extractions, and non-CG contexts. All paths embed modality timestamps (e.g., `DMR_20260402_191818`).
- **External data paths** — ChIP peaks, MeCP2 DiffBind, ATAC/K119ub/K27ac BEDs, DiffBind results, BigWigs, Hi-C loops. These do NOT change between runs.
- **Helper functions** — `load_dmr_bed()`, `load_diffbind_flex()`, `dedup_by_gene()`, `classify_chromatin_state()`, `patch_chooseHclustMet()`, `theme_biomodal()`
- **Pre-loaded data** — `mc_dmr`, `hmc_dmr` (deduplicated by gene), `bioqc`, `upstream`, `region_dmrs`
- **Color/label standards** — `COLORS`, `CHROMATIN_STATE_ORDER`, `CHROMATIN_STATE_COLORS`

### Section Script Convention

Each `section_NN_name.R` is self-contained: sources `_shared_config.R`, reads additional data if needed, produces plots in `downstream/plots/visualizations/{section_num}_{name}/` in 4 formats (PNG + PDF + SVG + JPG) via `multi_format_output.R`, and saves tables to `downstream/plots/visualizations/tables/`.

Section scripts must NOT contain hardcoded run paths — verify with:
```bash
grep -rn 'run-[0-9]' downstream/scripts/viz_sections/section_*.R  # should return 0 results
```

### DMR BED Format (modality output)

13-14 tab-delimited columns: `chr`, `start`, `end`, `num_contexts`, `mean_coverage`, `mean_mod_group1`, `mean_mod_group2`, `mod_fold_change`, `mod_difference`, `test_statistic`, `dmr_pvalue`, `dmr_qvalue`, `annotation`, [`gene`]

### Chromatin State Classification (7 categories, shared with Hi-C pipeline)

1. **Active_Promoter** — H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
2. **Repressed_Promoter** — H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
3. **Bivalent_Promoter** — K4me3+K27me3 overlap (pre-computed BED)
4. **Polycomb** — H3K27me3+ AND >2kb from TSS
5. **Active_Enhancer** — H3K27ac+ AND >2kb from TSS
6. **Poised_Enhancer** — H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
7. **Unmarked** — No marks

## Conda Environments

| Env | Used by | Key detail |
|-----|---------|------------|
| `env_nf` | Upstream DUET pipeline | Nextflow + biomodal CLI |
| `modality` | Downstream DMR calling | modality XPLR CLI |

Visualization R scripts use the system R 4.5.2 (`/usr/local/bin/Rscript`) directly — no conda env needed.

## Analysis Run History

| Run | Samples | Sex Covariate | Notes |
|-----|---------|---------------|-------|
| run-1 | 4 (shallow-seq) | Yes | No significant DMRs — sex confounded at n=2/group |
| run-2 | 4 (shallow-seq) | No | First significant results |
| run-3 | 4 (deep-seq) | No | Improved power |
| run-4 | 8 (deep-seq) | No | Added batch 2 replicates |
| run-5 | 8 (deep-seq) | Yes | **Current primary** — sex covariate properly powered |

## Gotchas

- **DMR paths embed timestamps.** Files like `DMR_mc_control__mutant_20260402_191818.bed` must match both the parent directory name and the filename. See `docs/UPDATE.md` for the update procedure.
- **Working directory matters.** Visualization scripts assume `getwd()` is `downstream/`. Running from the wrong directory breaks all relative path resolution.
- **MeCP2 is CUT&RUN, not ChIP.** Use "signal" or "mark" in labels — never "ChIP".
- **Non-CG contexts are noise.** CHG/CHH show <1% methylation, no significant DMRs. They exist for completeness.
- **Sections 34-36 (permutation tests)** require `regioneReloaded` and use `patch_chooseHclustMet()` to fix a crash with ≤2 rows.
- **Section 77 (MeCP2 aging trajectory)** is blocked on external DiffBind data from collaborator Jai — the script includes `stopifnot()` guards.
- **BigWig canonical source** is `/Users/zakiralibhai/sdsc/bigwigs/` — do NOT use `peaks/bigwigs/macs2.narrow.aug18.dedup/` (has 0-byte mutant files).

## HPC Paths (Expanse)

```
Base: /expanse/lustre/projects/csd940/zalibhai/biomodal/
├── evoC-run/input/          # Raw FASTQ files
├── evoC-run/output/         # DUET pipeline outputs
└── modality/                # Downstream DMR analysis
    ├── mm10/                # Reference region BED files
    ├── outputs/run-{N}/     # Per-run results
    └── config_*.txt         # Context-specific configs
```

## Documentation Hierarchy

1. **This file** — Top-level pipeline overview and entry point
2. `downstream/CLAUDE.md` — Detailed visualization architecture, section groupings, config structure
3. `downstream/docs/UPDATE.md` — Step-by-step run update procedure with checklist
4. `downstream/docs/FIGURES.md` — Figure descriptions and presentation mapping
5. `upstream/biomodal_docs.md` — DUET pipeline v1.5.0 reference
6. `downstream/modality/biomodal-workflow.md` — Modality XPLR v1.3 documentation
7. `downstream/docs/biomodal-interpretation-guide.md` — DMR interpretation guide

@downstream/CLAUDE.md
@downstream/Methylation_Paper.md
@downstream/TODOS.md
@downstream/PLAN.md
