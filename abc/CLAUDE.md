# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**ABC (Activity-By-Contact) Enhancer-Gene Linkage Analysis** for BAP1-KO vs wildtype mouse cerebellum (mm10). Uses the Broad Institute's ABC model to link enhancers to target genes, then computes differential ABC scores (ΔABC) to identify enhancer-gene pairs where linkage strength changes upon BAP1 knockout.

**Core hypothesis:** BAP1-KO eliminates H2AK119ub-mediated suppression of enhancer activity, leading to measurable changes in enhancer-gene linkage that correlate with differential gene expression and H2AK119ub signal changes.

**ABC formula:** `ABC score = (Enhancer Activity × Contact Frequency) / Σ(Activity × Contact)` over all candidate elements within 5 Mb of each gene.

## Running the Pipeline

### Prerequisites

- **HPC (Expanse):** `conda activate abc-env` — Snakemake, bedtools, MACS2, samtools
- **Local (macOS):** Python 3.8+ with pandas/numpy/scipy/matplotlib, R with ggplot2/patchwork/ggpointdensity

### Step-by-step execution

Steps 0–5 and 9 run on HPC (Expanse via SLURM). Steps 6–8, 10 run locally.

```bash
# HPC — Reference generation (one-time)
sbatch scripts/step0a_gene_annotations.sb
bash scripts/step0e_consensus_to_narrowpeak.sh

# HPC — Run ABC pipeline via Snakemake
sbatch scripts/step4_run_abc.sb

# HPC — QC validation
sbatch scripts/step5_qc.sb

# Local — Delta computation and integration
python scripts/step6_delta_abc.py
python scripts/step7_rnaseq_integration.py

# HPC — Cross-reference with loops and H3K27ac
python scripts/step9_cross_reference.py

# HPC — Extract K119ub signal at enhancers
Rscript scripts/preprocess_k119ub_enhancer_signal.R \
  --abc-pairs results/delta_abc_all_pairs.tsv \
  --bw-dir /expanse/lustre/projects/csd940/zalibhai/loop-class/heatmaps/ \
  --output results/k119ub_enhancer_signal.tsv

# Local — K119ub correlation analysis (8-panel visualization)
cd abc && Rscript scripts/step10_k119ub_abc_correlation.R
```

No single entry point script exists — steps are run individually. Steps 6 and 7 have no CLI args (paths hardcoded at top of each script).

## Architecture

### Critical Design Decision: Consensus Enhancer Universe

Both WT and KO ABC runs use the **same** `consensus_all.bed` (75,371 ATAC peaks from Batch 1 union). This is injected via symlink into both biosamples' peak directories. Using identical candidate regions ensures ΔABC measures true E-P linkage changes, not artifacts from different enhancer definitions.

### Data Flow

```
consensus_all.bed ──┐
ATAC tagAlign ──────┤  Step 4: ABC Snakemake
Hi-C .hic files ────┤  (ATAC-only mode, no H3K27ac BAMs)
mm10 references ────┘
        │
        ▼
AllPutative.tsv.gz (WT + KO)
        │
        ▼  Step 6: Outer join, remove self-promoters, compute ΔABC
delta_abc_all_pairs.tsv (180K E-G pairs)
        │
   ┌────┴────────────────┐
   ▼                     ▼
Step 7: + RNA-seq     Step 9: + loops + H3K27ac
   │                     │
   ▼                     ▼
gene_level_summary    delta_abc_annotated
(13,588 genes)        (180K pairs + flags)
                         │
                         ▼  Step 10: + K119ub bigWigs
                      k119ub_abc_correlation
                      (8-panel figures + stats)
```

### Key Files

| File | Rows | Description |
|------|------|-------------|
| `results/delta_abc_all_pairs.tsv` | 180,423 | All distal E-G pairs with ΔABC, Δ(A×C), activity, contact |
| `results/gene_level_summary.tsv` | 13,588 | One row per gene: strongest enhancer, aggregate stats, DE results |
| `results/delta_abc_annotated.tsv` | 180,423 | Same as above + H3K27ac overlap flags |
| `results/delta_abc_with_rnaseq.tsv` | 113,453 | E-G pairs merged with log2FC and padj |
| `results/k119ub_enhancer_signal.tsv` | ~15K | Per-replicate K119ub at enhancers (from HPC preprocessing) |
| `results/k119ub_abc_enhancer_merged.tsv` | ~9K | Enhancer aggregates with K119ub + ABC deltas |

### Script Dependencies

- **step6** reads from HPC output path (`ABC-Enhancer-Gene-Prediction/results/*/Predictions/`). Transfer `delta_abc_all_pairs.tsv` to local after running.
- **step7** reads `delta_abc_all_pairs.tsv` (local) + RNA-seq Excel from `tads/` sibling directory.
- **step9** requires bedtools (HPC) and reads loops from parent `mariner_hi-c/` directory.
- **step10** reads from `results/` (relative to `abc/` working directory). Must `cd abc` first.
- **preprocess_k119ub** runs on HPC, needs rtracklayer + GenomicRanges. Output transferred to local.

## ABC Pipeline Configuration

Key mm10-specific changes from the default hg38 config:

```yaml
genome_size: mm              # MACS2 genome size (not hs)
use_qnorm: False             # K562 DHS quantile ref inappropriate for mouse ATAC
threshold: 0.02              # ABC score cutoff for thresholded predictions
```

Running in **ATAC-only mode** — activity is ATAC signal only because H3K27ac BAMs are unavailable (only bigWigs exist). H3K27ac peaks are used downstream for annotation in Step 9, not for ABC scoring.

## Key Thresholds

| Parameter | Value | Used in |
|-----------|-------|---------|
| ABC score threshold | 0.02 | Step 6 (at least one condition ≥ 0.02) |
| ΔABC significance | 0.01 | Steps 7, 9 (|ΔABC| > 0.01 = gained/lost) |
| RNA-seq DE | padj < 0.05, |log2FC| > 0.5 | Step 7 concordance |

## Key Findings

- **Unnormalized Δ(A×C) outperforms normalized ΔABC** for concordance with DE: 65.3% vs 58.8%. ABC normalization compresses real activity changes when BAP1-KO causes widespread chromatin remodeling.
- **940 dysregulated genes** with both significant ΔABC and DE; 58.8% concordant (binomial p = 6.84e-08).
- **H3K27ac+ enhancers** show asymmetric loss (1.17× lost/gained ratio) vs symmetric H3K27ac- enhancers.
- **Loop-ABC directional concordance at chance (51.4%)** — expected due to methodology differences (loop strength ≠ ABC score), not model failure.

## Reference Files

Located in `reference/`:
- `mm10_genes.bed` — 28,230 genes from GENCODE vM25 (protein_coding + lincRNA + processed_transcript)
- `mm10_tss.bed` — 500bp TSS regions (same gene set)
- `mm10_merged_blacklist.bed` — Lab + ENCODE blacklists merged (3,581 regions)
- `mm10_ubiquitous_genes.tsv` — Top 5% by baseMean (730 genes, used by ABC)

## Documentation

- `docs/abc-prd-v2.md` — Full product requirements document
- `docs/abc-execution-plan.md` — Step-by-step implementation plan (~1600 lines)
- `docs/abc-analysis-context.md` — QC findings and results from Steps 5–9
- `docs/methods.md` — Scientific methodology
