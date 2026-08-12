# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Downstream analysis and visualization of biomodal DUET evoC differential methylation (5mC and 5hmC) between BAP1-KO mutant and wildtype control mice (mm10/GRCm38). Uses the modality XPLR CLI for GLM-based DMR calling, followed by a modular R visualization pipeline (85 section scripts).

**Samples:** 8 total (4 control + 4 mutant, 2 batches × 1 per sex per condition)
**Primary context:** CG (CpG) — CHG/CHH show minimal signal (<1% methylation)
**Current primary run:** run-5 (deep-seq, 8 samples, sex covariate included)
**Key finding:** 92.3% of co-significant genes exhibit coordinated mC↑/hmC↓, indicating TET-mediated demethylation block

## Running the Pipeline

### DMR Analysis (HPC with SLURM)

```bash
cd downstream/modality
sbatch scripts/run_modality.sb CG    # Primary CpG analysis
sbatch scripts/run_modality.sb CHG   # Optional
sbatch scripts/run_modality.sb CHH   # Optional
```

Requires conda environment `modality`. The SLURM script copies `config_${CONTEXT}.txt` → `config.txt` then runs `core-workflow-v1.3.sh`.

### Visualization Pipeline

```bash
# Run all sections sequentially (must run from downstream/ directory)
cd downstream/
bash scripts/viz_sections/run_all_sections.sh

# Run individual section
Rscript scripts/viz_sections/section_04_volcano_plots.R

# Batch runners for grouped sections
bash scripts/viz_sections/run_sections_22_26.sh   # demethylation ratio + TET
bash scripts/viz_sections/run_sections_61_63.sh   # stoichiometry + heatmap
bash scripts/viz_sections/run_sections_74_78.sh   # neuronal + aging

# Permutation tests (HPC, heavy compute)
sbatch scripts/viz_sections/run_permutation.sb
```

All section scripts source `_shared_config.R` which loads data, defines helpers, and sets paths. **Working directory must be `downstream/`.**

### Updating for a New Modality Run

Follow `docs/UPDATE.md` — covers the 4 files with hardcoded run paths, timestamp discovery, and verification grep.

## Architecture

### Two-Stage Pipeline

```
Zarr stores (from upstream DUET pipeline)
    ↓ [modality XPLR CLI: core-workflow-v1.3.sh]
    ├── Step 1: Biological QC (correlations, PCA)
    ├── Step 2: Feature Extraction (count, mean, regional-frac)
    ├── Step 3: DMR Calling (GLM per region × context)
    └── Step 4: DMR Visualization (volcano plots, filtered BED)
    ↓
Post-processing + Visualization (R scripts)
    ├── _shared_config.R (centralized config, data loading, helpers)
    └── section_01 through section_78 (independent visualization modules)
```

### Modular Visualization System (`scripts/viz_sections/`)

85 section scripts, each independent and sourcing `_shared_config.R` for shared state. Major groups:

| Sections | Theme |
|----------|-------|
| 01–09 | Core DMR analysis: QC, correlation, statistics, volcano, coordinated changes, top genes, enrichment |
| 10–10f | Chromatin state classification (7-category system shared with Hi-C pipeline) |
| 11–19 | CUT&RUN/ChIP integration: MeCP2, ATAC-seq, H2AK119ub, H3K27ac |
| 20–21 | Coordinated RNA-seq, discordant mC/hmC |
| 22–26 | Demethylation ratio, baseline hmC, DNMT3A prediction, TET-KO comparison |
| 27–31 | Hi-C loop anchor integration, A/B compartment mapping, Polycomb targets, MeCP2 loops |
| 32–37 | CHG exploratory, multi-mark DiffBind, permutation tests (sections 34–36 use regioneReloaded) |
| 38–41 | H3K36me2/me3 analysis: gene body, boundary, combined, DNMT3A vs DNMT3B |
| 42–45 | Max-significance gene lists, CG exploratory, allele-specific, Field BAP1 chr8 |
| 46–50 | Genome browser loci, CTCF anchor overlay, CpG island ubiquitination, HOMER motif enrichment |
| 51–60 | Non-CG methylation: MeCP2 CG/non-CG, CpG-resolution gene body, TAD organization, Ecker validation |
| 61–63 | Stoichiometry mechanism, multifeature regression, MeCP2 master heatmap |
| 64–72 | Global methylation levels, subcompartment methylation, MeCP2/K119ub unmethylated, neuronal characterization |
| 73–78 | Neuronal chromatin remodeling, gene-set overlaps, MeCP2 aging trajectory, stoichiometry broad |
| `compare_shallow_vs_deep.R` | Standalone: shallow-seq vs deep-seq comparison |

### Shared Config (`_shared_config.R`)

Central configuration loaded by every section script (620 lines). Contains:

**Path blocks (change when updating runs — see `docs/UPDATE.md`):**
- `DATA_PATHS` — Hardcoded paths to run-5 DMR BED files (13 entries: gene body + 5 regions × mC/hmC + BioQC JSON)
- `EXTRACT_PATHS` — Per-sample regional fractions (gene body mC/hmC)
- `CHG_DATA_PATHS`, `CHH_DATA_PATHS` — Non-CG context DMR paths
- `CHG_EXTRACT_PATHS`, `CHH_EXTRACT_PATHS` — Non-CG per-sample extractions

**External data (do NOT change between runs):**
- `CHIP_PEAK_FILES` — ChIP-seq BED paths for chromatin state classification
- `MECP2_FILES` — MeCP2 DiffBind results and BigWigs
- `ATAC_FILES`, `K119UB_FILES`, `H3K27AC_FILES` — Condition-specific peak BEDs
- `DIFFBIND_FILES` — Quantitative differential binding (ATAC, K27ac, K27me3, K119ub)
- `H3K36ME2_FILES`, `H3K36ME3_FILES` — H3K36me2/me3 DiffBind results
- `METHYLATION_BIGWIGS`, `HISTONE_BIGWIGS`, `ECKER_BIGWIGS` — Signal tracks
- `LOOP_FILES` — Hi-C loop annotations

**Helpers:**
- `load_dmr_bed()` — Parses modality DMR BED format (13-14 columns), adds significance/direction
- `load_diffbind_flex()` — Loads DiffBind results with flexible column schema (seqnames or Summit_Chr)
- `dedup_by_gene()` — Keeps lowest q-value row per gene
- `load_chip_peaks()` — BED → GRanges
- `dmr_to_granges()` — DMR data frame → GRanges with metadata
- `compute_chip_overlaps()` — Computes overlap booleans for 6 histone marks
- `classify_chromatin_state()` — 7-category priority system matching Hi-C pipeline
- `patch_chooseHclustMet()` — Bug fix for regioneReloaded with ≤2 rows (used by sections 34–36)
- `theme_biomodal()` — Consistent ggplot2 theme
- `COLORS` — Standardized color palettes for condition, direction, methylation type, marks
- `CHROMATIN_STATE_ORDER`, `CHROMATIN_STATE_COLORS` — Consistent state labels and colors

**Pre-loaded data:** `mc_dmr`, `hmc_dmr` (deduplicated by gene), `bioqc` (JSON), `upstream` (CSV), `region_dmrs` (6 regions)

### DMR BED Format (modality output)

13-14 tab-delimited columns: `chr`, `start`, `end`, `num_contexts`, `mean_coverage`, `mean_mod_group1`, `mean_mod_group2`, `mod_fold_change`, `mod_difference`, `test_statistic`, `dmr_pvalue`, `dmr_qvalue`, `annotation`, [`gene`]

### Chromatin State Classification (7 Categories)

Used in section_10 and shared with the Hi-C loop analysis pipeline:
1. **Active_Promoter** — H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS
2. **Repressed_Promoter** — H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS
3. **Bivalent_Promoter** — K4me3+K27me3 overlap (pre-computed)
4. **Polycomb** — H3K27me3+ AND >2kb from TSS
5. **Active_Enhancer** — H3K27ac+ AND >2kb from TSS
6. **Poised_Enhancer** — H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb
7. **Unmarked** — No marks

## Configuration

### Modality Config (`config_CG.txt`)

Key parameters:
- `Zarr` — Path to upstream Zarr store
- `Group_Column=condition` — Metadata column for contrast
- `Condition_Order=control mutant` — Reference vs test
- `Covariates=sex` — Include sex covariate (run-5); removing it yields more DMRs but confounds sex/genotype
- `Depth_Filter=10` — Minimum coverage threshold
- `Overdispersion=False` — GLM overdispersion correction toggle

### Analysis Run History

| Run | Samples | Sex Covariate | Notes |
|-----|---------|---------------|-------|
| run-1 | 4 (shallow-seq) | Yes | No significant DMRs — sex confounded at n=2/group |
| run-2 | 4 (shallow-seq) | No | First significant results; archived at `plots/visualizations_run2_shallow-seq/` |
| run-3 | 4 (deep-seq) | No | Improved power; archived at `plots/visualizations_run_3_deep-4/` |
| run-4 | 8 (deep-seq) | No | Added batch 2 replicates; archived at `plots/visualizations_run4_deep-8/` |
| run-5 | 8 (deep-seq) | Yes | **Current primary** — sex covariate properly powered with 8 samples |

## Output Locations

- **DMR results:** `modality/outputs/run-5/outputs_CG/Results/{region}/DMR_*/`
- **QC reports:** `modality/outputs/run-5/outputs_CG/Results/BioQC_*/`
- **Visualization plots:** `plots/visualizations/{section_num}_{name}/`
- **Export tables:** `plots/visualizations/tables/`
- **Archived prior runs:** `plots/visualizations_run2_shallow-seq/`, `plots/visualizations_run_3_deep-4/`, `plots/visualizations_run4_deep-8/`

## Dependencies

**R packages (loaded by `_shared_config.R`):** tidyverse, ggplot2, patchwork, ggrepel, RColorBrewer, scales, pheatmap, jsonlite, clusterProfiler, enrichplot, org.Mm.eg.db, ggVennDiagram, GenomicRanges, rtracklayer, TxDb.Mmusculus.UCSC.mm10.knownGene

**Additional R packages (loaded by individual sections as needed):** regioneReloaded (sections 34–36), Gviz (section 46)

**External:** modality XPLR CLI (conda env `modality`), bedToBigBed (UCSC utility at `modality/scripts/bedToBigBed`)

**Shared utility:** Sources `multi_format_output.R` from `../../scripts/utils/` (generates PDF + PNG + SVG)

## Key Caveats

- DMR BED filenames contain timestamps (e.g., `DMR_mc_control__mutant_20260402_191818.bed`) — paths in `_shared_config.R` are hardcoded to run-5 timestamps. See `docs/UPDATE.md` for the update procedure.
- `compare_shallow_vs_deep.R` hardcodes both old and new run timestamps for side-by-side comparison — update `SHALLOW_PATHS` when rotating runs.
- Section scripts must NOT contain hardcoded run paths — all data access goes through `_shared_config.R`. Verify with: `grep -rn 'run-[0-9]' scripts/viz_sections/section_*.R` (should return zero results).
- Sections 34–36 (permutation tests) require `regioneReloaded` and use `patch_chooseHclustMet()` from `_shared_config.R` to fix a crash with ≤2 rows.
- Non-CpG contexts (CHG/CHH) show <1% methylation and no significant DMRs.
- MeCP2 data is CUT&RUN, not ChIP-seq — use "signal" or "mark" in labels, never "ChIP".
