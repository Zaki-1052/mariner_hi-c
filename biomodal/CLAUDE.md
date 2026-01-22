# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Biomodal DUET evoC Differential Methylation Analysis Pipeline** - Analyzes differential methylation (5mC and 5hmC) between BAP1-KO mutant and wildtype control mice using biomodal's dual-epigenetic sequencing technology with 6bp resolution on the mm10/GRCm38 mouse genome.

**Key Technologies:**
- Upstream: biomodal DUET pipeline v1.5.0 (Nextflow-based)
- Downstream: modality XPLR CLI v1.0+ (GLM-based DMR calling)
- Data format: Zarr stores (HDF5-backed methylation data)
- HPC: SLURM job scheduler on Expanse cluster

## Running the Pipeline

### Upstream (DUET Pipeline) - Raw Sequencing Processing

```bash
# Submit to SLURM
sbatch upstream/evoc_run.sb

# Or run directly (not recommended - requires 48hr runtime)
biomodal run duet \
  --input-path /path/to/fastq \
  --output-path /path/to/output \
  --meta-file metadata.csv \
  --run-name evoC_Bap1_run \
  --mode 6bp \
  --additional-profile deep_seq \
  --chg-chh-contexts \
  --reference-genome-name GRCm38
```

**Required conda environment:** `env_nf`

### Downstream (Modality XPLR) - DMR Analysis

```bash
# Submit to SLURM for specific methylation context
sbatch downstream/modality/run_modality.sb CG    # CpG context (primary)
sbatch downstream/modality/run_modality.sb CHG   # CHG context
sbatch downstream/modality/run_modality.sb CHH   # CHH context

# Or run manually
cd downstream/modality
cp config_CG.txt config.txt
./core-workflow-v1.3.sh
```

**Required conda environment:** `modality`

## Pipeline Architecture

### Data Flow

```
FASTQ files (4 samples × 2 lanes)
    ↓ [biomodal DUET v1.5.0]
Zarr stores (CG/CHG/CHH contexts)
    ↓ [modality XPLR]
DMR results (BED files, volcano plots, HTML reports)
```

### Upstream Outputs (duet-1.5.0_evoC_Bap1_run_6bp/)

| Directory | Contents |
|-----------|----------|
| `sample_outputs/bams/` | Deduplicated BAM files |
| `sample_outputs/zarr_store/` | Methylation data by context (CG, CHG, CHH) |
| `sample_outputs/variant_call_files/` | Germline VCF files |
| `reports/` | Summary metrics (CSV, XLSX) |
| `diagnostics/` | FastQC reports |

### Downstream Outputs (modality/outputs/)

```
outputs/run-{N}/outputs_{CONTEXT}/
├── Results/
│   ├── BioQC_*/                    # QC metrics + JSON
│   └── {region}/                   # Per-region analysis
│       ├── Extract_*/              # Feature extraction (count, mean, regional-frac)
│       ├── DMR_*/                  # DMR calling results (BED files)
│       └── DMR_Report_*/           # Visualization outputs
└── Reports/
    ├── biological_qc_report*.html  # QC visualizations
    └── {region}/                   # HTML reports per region
```

## Configuration

### Config File Parameters (config.txt)

| Parameter | Description | Example |
|-----------|-------------|---------|
| `Zarr` | Path to zarr store | `/path/to/CG/sample.zarrz` |
| `Metadata` | Sample metadata TSV | `metadata.tsv` |
| `Group_Column` | Condition column name | `condition` |
| `Condition_Order` | Reference vs test | `control mutant` |
| `Covariates` | Optional covariates (space-separated) | `sex` or blank |
| `Regions_Directory` | BED files for regions | `mm10/` |
| `Depth_Filter` | Minimum coverage threshold | `10` |
| `Overdispersion` | Apply overdispersion correction | `True` or `False` |

### Metadata Format (metadata.tsv)

```
sample_id	condition	sex
evoC-Bap1-ctrl-M	control	male
evoC-Bap1-mut-M	mutant	male
evoC-Bap1-ctrl-F	control	female
evoC-Bap1-mut-F	mutant	female
```

### Reference Regions (mm10/)

- `gencode.vM25.mouse.genes.annotation.bed.gz` - Gene bodies
- `gencode.vM25.mouse.promoters.annotation.bed.gz` - Promoter regions
- `gencode.vM25.mouse.cpg_islands.annotation.bed.gz` - CpG islands
- `gencode.vM25.mouse.cpg_shores.annotation.bed.gz` - CpG shores
- `gencode.vM25.mouse.cpg_shelves.annotation.bed.gz` - CpG shelves
- `gencode.vM25.mouse.tss_region.annotation.bed.gz` - TSS regions

## Key Samples

| Sample ID | Condition | Sex | Genotype |
|-----------|-----------|-----|----------|
| evoC-Bap1-ctrl-F | Control | Female | Bap1-ff |
| evoC-Bap1-ctrl-M | Control | Male | Bap1-ff |
| evoC-Bap1-mut-F | Mutant | Female | Bap1-ff-cre |
| evoC-Bap1-mut-M | Mutant | Male | Bap1-ff-cre |

## DMR Calling Workflow (core-workflow-v1.3.sh)

The pipeline executes 4 sequential steps:

1. **Biological QC** - Sample correlations, PCA, coverage statistics
2. **Feature Extraction** - Per-region methylation statistics (count, mean, regional-frac)
3. **DMR Calling** - GLM-based differential methylation (Benjamini-Hochberg FDR)
4. **DMR Visualization** - Volcano plots, p-value histograms, filtered BED files

Each step uses the `modality` CLI:
```bash
modality biological-qc ...
modality get count/mean/regional-frac ...
modality dmr call ...
modality dmr plot ...
```

## Critical Analysis Considerations

### Sex Confounding Issue

With n=2 per condition (4 samples total), sex effects are confounded with BAP1 effects:
- **Run-1 (with sex covariate)**: No significant DMRs
- **Run-2 (without sex covariate)**: Significant DMRs but cannot separate sex vs genotype effects

### Methylation Contexts

| Context | Signal Level | DMR Analysis |
|---------|--------------|--------------|
| **CG (CpG)** | Primary (~82% methylated) | Main analysis |
| **CHG** | Low (<1% methylated) | No significant results |
| **CHH** | Very low (<1% methylated) | Not analyzed |

### Expected Correlation Ranges

| Context | Within-Group | Between-Group |
|---------|--------------|---------------|
| CG mC | 0.76-0.79 | 0.76-0.77 |
| CG hmC | 0.48-0.51 | 0.47-0.49 |
| CHG | 0.36-0.52 | Highly variable (noise) |

## Documentation Hierarchy

1. `CLAUDE.md` - This file (pipeline overview)
2. `upstream/biomodal_docs.md` - DUET pipeline v1.5.0 guide
3. `downstream/biomodal-interpretation-guide.md` - DMR interpretation
4. `downstream/modality/biomodal-workflow.md` - Modality XPLR v1.3 documentation
5. `downstream/modality/outputs/analysis_summary.md` - Results summary

## Common Tasks

### Re-run DMR Analysis with Different Covariates

1. Edit `downstream/modality/config_CG.txt`
2. Modify `Covariates=` line (blank for none, or `sex` to include)
3. Run: `sbatch downstream/modality/run_modality.sb CG`

### Examine Top DMRs

```bash
# Gene body mC DMRs (sorted by q-value)
zcat outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_*/num_mc_dmr_results.bed.gz | \
  sort -k11,11g | head -20

# Filter to q<0.05
zcat outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_Report_*/num_mc_dmr_results_max_q_0.05.bed.gz
```

### View QC Reports

HTML reports are located in:
- `outputs/run-N/outputs_CG/Reports/biological_qc_report*.html` - Sample correlations, PCA
- `outputs/run-N/outputs_CG/Reports/{region}/*.html` - Per-region volcano plots

## Resource Requirements

| Pipeline Stage | CPUs | Memory | Runtime |
|----------------|------|--------|---------|
| DUET upstream | 16 | 128GB | ~48hr |
| Modality downstream | 16 | 128GB | ~12hr |

## HPC Paths (Expanse)

```
Base: /expanse/lustre/projects/csd940/zalibhai/biomodal/
├── evoC-run/input/          # Raw FASTQ files
├── evoC-run/output/         # DUET pipeline outputs
└── modality/                # Downstream analysis
    ├── mm10/                # Reference regions
    ├── outputs/             # DMR results
    └── config_*.txt         # Context-specific configs
```
