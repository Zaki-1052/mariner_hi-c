# WGBS nf-core/methylseq Setup — Session Summary

## What was done

Set up nf-core/methylseq (v4.2.0) with Bismark for whole-genome bisulfite sequencing of BAP1-KO mouse cerebellum (mm10/GRCm38). The goal is non-CG methylation extraction (CpG + CHG + CHH) from WGBS data.

## Files created

All files are in `biomodal/wgbs/`:

| File | Purpose |
|------|---------|
| `nextflow.config` | Expanse SLURM + Singularity configuration with resource limits and Bismark-specific process overrides |
| `run_methylseq.sb` | SLURM submission script — runs the full pipeline with all relevant flags |
| `download_fastqs.sb` | Template SLURM script for downloading FASTQs from sequencing facility FTP — **needs FTP URL and file list filled in** |
| `samplesheet.csv` | Empty template — **needs FASTQ paths filled in** |

## Key parameter choices

- **`--cytosine_report`** — the critical flag for non-CG. Runs Bismark's `coverage2cytosine --CX` to produce genome-wide CpG/CHG/CHH cytosine-level reports
- **`--comprehensive`** — merges strand-specific calls into context-dependent output files
- **`--save_reference`** — persists the built Bismark index for reuse
- **`--save_align_intermeds`** — keeps deduplicated BAMs
- **`--save_trimmed`** / **`--unmapped`** — saves intermediate reads
- **`--run_qualimap`** / **`--run_preseq`** — enables extra QC

## Why no separate trimming

nf-core/methylseq runs Trim Galore internally (FastQC → Trim Galore → Bismark). Pre-trimming with fastp would double-trim and is not needed. Feed raw FASTQs directly.

## What still needs to happen

1. **Check Nextflow version** on Expanse — pipeline v4.2.0 requires Nextflow >= 25.04.0:
   ```bash
   conda activate env_nf
   nextflow -version
   # If too old: conda create -n env_nf -c bioconda nextflow=25.04.3
   ```

2. **Get FASTQ download link** from sequencing facility — update `FTP_BASE` and file list in `download_fastqs.sb`

3. **Download FASTQs** on Expanse:
   ```bash
   sbatch download_fastqs.sb
   ```

4. **Fill in `samplesheet.csv`** with absolute paths to the downloaded FASTQs. Format:
   ```csv
   sample,fastq_1,fastq_2,genome
   ctrl_F,/expanse/.../fastqs/ctrl_F_R1.fastq.gz,/expanse/.../fastqs/ctrl_F_R2.fastq.gz,
   ```
   If a sample was sequenced across multiple lanes, use the same sample name on multiple rows — the pipeline will concatenate them automatically.

5. **Deploy to Expanse**:
   ```bash
   # Push to GitHub from Mac
   # On Expanse:
   cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c && git pull
   mkdir -p /expanse/lustre/projects/csd940/zalibhai/methylseq/logs
   cp biomodal/wgbs/nextflow.config /expanse/lustre/projects/csd940/zalibhai/methylseq/
   cp biomodal/wgbs/run_methylseq.sb /expanse/lustre/projects/csd940/zalibhai/methylseq/
   cp biomodal/wgbs/download_fastqs.sb /expanse/lustre/projects/csd940/zalibhai/methylseq/
   cp biomodal/wgbs/samplesheet.csv /expanse/lustre/projects/csd940/zalibhai/methylseq/
   ```

6. **Run the pipeline**:
   ```bash
   cd /expanse/lustre/projects/csd940/zalibhai/methylseq
   sbatch run_methylseq.sb
   ```

## Expanse directory layout (target)

```
/expanse/lustre/projects/csd940/zalibhai/methylseq/
├── nextflow.config
├── run_methylseq.sb
├── download_fastqs.sb
├── samplesheet.csv
├── logs/
├── fastqs/          # downloaded raw FASTQs
└── results/         # pipeline output (created by nextflow)
```

## Context

- The cloned nf-core/methylseq repo is at `biomodal/methylseq/` (v4.2.0, not modified)
- Existing biomodal evoC analysis used a different pipeline (biomodal DUET + modality XPLR) for 5mC/5hmC DMR calling on 4 samples
- This new WGBS pipeline will use Bismark for conventional bisulfite-seq methylation, with the specific goal of extracting non-CG (CHG/CHH) methylation levels
- The `non_cg_analysis/` directory in `biomodal/` contains downstream analysis scripts for non-CG data
