# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Downstream analysis and visualization of biomodal DUET evoC differential methylation (5mC and 5hmC) between BAP1-KO mutant and wildtype control mice (mm10/GRCm38). Uses the modality XPLR CLI for GLM-based DMR calling, followed by a modular R visualization pipeline.

**Samples:** 4 total (2 control + 2 mutant, 1 per sex per condition)
**Primary context:** CG (CpG) — CHG/CHH show minimal signal
**Key finding:** 92.3% of co-significant genes exhibit coordinated mC↑/hmC↓, indicating TET-mediated demethylation block

## Running the Pipeline

### DMR Analysis (HPC with SLURM)

```bash
cd downstream/modality
sbatch run_modality.sb CG    # Primary CpG analysis
sbatch run_modality.sb CHG   # Optional
sbatch run_modality.sb CHH   # Optional
```

Requires conda environment `modality`. The SLURM script copies `config_${CONTEXT}.txt` → `config.txt` then runs `core-workflow-v1.3.sh`.

### Visualization Pipeline

```bash
# Run all 19 sections sequentially (must run from downstream/ directory)
cd downstream/
bash scripts/viz_sections/run_all_sections.sh

# Run individual section
Rscript scripts/viz_sections/section_04_volcano_plots.R
```

All section scripts source `_shared_config.R` which loads data, defines helpers, and sets paths. Working directory must be `downstream/`.

### Post-Processing (mC/hmC Comparison)

```bash
cd modality/DMR_Genes_CG/
bash compare_mc_hmc.sh           # Compare mC vs hmC directions
bash significant_to_bigbed.sh    # Convert BED → BigBed for genome browser
```

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
    └── section_01 through section_10 (independent visualization modules)
```

### Modular Visualization System (`scripts/viz_sections/`)

Each section script is independent and sources `_shared_config.R` for shared state:

| Section | Purpose |
|---------|---------|
| `section_01_qc_overview.R` | Sample quality, coverage, upstream metrics |
| `section_02_correlation.R` | mC/hmC sample correlation heatmaps |
| `section_03_dmr_statistics.R` | DMR counts across 6 genomic regions |
| `section_04_volcano_plots.R` | Significance vs effect size plots |
| `section_05_coordinated_changes.R` | mC↑/hmC↓ coordinated pattern analysis |
| `section_06_top_genes.R` | Ranked gene tables (Syt1, Trpm3, etc.) |
| `section_07_effect_size.R` | Methylation change distributions |
| `section_08_enrichment.R` | GO/KEGG functional enrichment |
| `section_09_summary.R` | Key findings and summary tables |
| `section_10_chromatin_state.R` | ChIP-seq integration, 7-state classification |
| `section_11` through `section_19` | MeCP2, ATAC-seq, H2AK119ub, H3K27ac integrations |
| `compare_shallow_vs_deep.R` | Standalone: shallow-seq (run-2) vs deep-seq (run-3) comparison |

### Shared Config (`_shared_config.R`)

Central configuration loaded by every section script. Contains:
- `DATA_PATHS` — Hardcoded paths to run-3 (deep-seq) DMR BED files (timestamped filenames)
- `CHIP_PEAK_FILES` — ChIP-seq BED paths for chromatin state classification
- `load_dmr_bed()` — Parses modality DMR BED format (13-14 columns)
- `theme_biomodal()` — Consistent ggplot2 theme
- `classify_chromatin_state()` — 7-category priority system matching Hi-C pipeline
- `COLORS` — Standardized color palettes for condition, direction, methylation type
- Pre-loaded data: `mc_dmr`, `hmc_dmr`, `bioqc` (JSON), `region_dmrs` (6 regions)

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
7. **Other** — No marks

## Configuration

### Modality Config (`config_CG.txt`)

Key parameters:
- `Zarr` — Path to upstream Zarr store
- `Group_Column=condition` — Metadata column for contrast
- `Condition_Order=control mutant` — Reference vs test
- `Covariates=sex` — Optional; removing sex covariate yields significant DMRs (run-2)
- `Depth_Filter=10` — Minimum coverage threshold
- `Overdispersion=False` — GLM overdispersion correction toggle

### Three Analysis Runs

- **run-1** (with sex covariate): No significant DMRs — sex confounded with genotype at n=2/group
- **run-2** (sex removed, shallow-seq): Significant DMRs but cannot fully separate sex vs genotype effects. Visualization outputs archived at `plots/visualizations_run2_shallow-seq/`
- **run-3** (sex removed, deep-seq): Current primary results — deeper sequencing for improved statistical power. Visualization pipeline and `_shared_config.R` now point to run-3

## Output Locations

- **DMR results:** `modality/outputs/run-3/outputs_CG/Results/{region}/DMR_*/`
- **QC reports:** `modality/outputs/run-3/outputs_CG/Results/BioQC_*/`
- **Archived run-2 plots:** `plots/visualizations_run2_shallow-seq/`
- **Visualization plots:** `plots/visualizations/{01-10}_*/`
- **Export tables:** `plots/visualizations/tables/`
- **BigBed files:** `modality/DMR_Genes_CG/bigbed/`

## Dependencies

**R packages:** tidyverse, ggplot2, patchwork, ggrepel, RColorBrewer, scales, pheatmap, jsonlite, clusterProfiler, enrichplot, org.Mm.eg.db, ggVennDiagram, GenomicRanges, rtracklayer, TxDb.Mmusculus.UCSC.mm10.knownGene

**External:** modality XPLR CLI (conda env `modality`), bedToBigBed (UCSC utility)

**Shared utility:** Sources `multi_format_output.R` from `../../scripts/utils/` (generates PDF + PNG + SVG)

## Key Caveats

- DMR BED filenames contain timestamps (e.g., `DMR_mc_control__mutant_20260221_190322.bed`) — paths in `_shared_config.R` are hardcoded to run-3 timestamps; `compare_shallow_vs_deep.R` hardcodes run-2 timestamps for comparison
- With n=2 per condition, sex and genotype effects are confounded — results from run-2 (no sex covariate) are used but should be interpreted with this caveat
- Non-CpG contexts (CHG/CHH) show <1% methylation and no significant DMRs
