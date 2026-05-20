# ABC Enhancer-Gene Linkage Analysis

Differential Activity-By-Contact (ABC) analysis of enhancer-gene linkage in BAP1-KO vs wildtype mouse cerebellum (mm10). Uses the [Broad Institute ABC model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) to link distal enhancers to target genes, then computes differential ABC scores to identify enhancer-gene pairs where linkage strength changes upon BAP1 knockout.

## Biological Context

BAP1 is the catalytic subunit of the PR-DUB complex, a deubiquitinase that removes H2AK119ub from chromatin. In BAP1-KO cerebellum, loss of deubiquitinase activity leads to aberrant accumulation of H2AK119ub at regulatory elements, altering enhancer activity and downstream gene expression. This pipeline quantifies those changes through the lens of 3D enhancer-gene linkage.

**ABC score** = (Enhancer Activity x Contact Frequency) / Sum(Activity x Contact) over all candidate elements within 5 Mb of each gene, where activity is the geometric mean of ATAC-seq and H3K27ac signal: `sqrt(ATAC.RPM x H3K27ac.RPM)`.

## Pipeline Overview

```
Step 0:  Reference generation (gene annotations, consensus peaks, blacklists)
Step 4:  Run ABC model via Snakemake (WT + KO, shared enhancer universe)
Step 5:  QC validation of ABC predictions
Step 6:  Compute ΔABC between conditions (outer join, self-promoter removal)
Step 7:  Integrate with RNA-seq differential expression
Step 9:  Cross-reference with Hi-C loops and H3K27ac peaks
Step 9b: Paired-anchor analysis (loop-ABC physical overlap)
Step 10: K119ub-ABC correlation at enhancers (DiffBind)
Step 11: Enhancer phenotypic subset analysis + HOMER motif enrichment
Step 12: Activity vs contact decomposition scatter plots
Step 13: Discordant gene analysis + GO/KEGG enrichment + K119ub concordance
Step 14: Loop type decomposition at enhancer subsets
Step 15: ABC component (activity vs contact) vs expression
Step 16: K119ub-filtered scatter (FDR < 0.05 sites)
```

## Requirements

| Environment | Software |
|---|---|
| HPC (SDSC Expanse) | `conda activate abc-env` — Snakemake, bedtools, MACS2, samtools, HOMER |
| HPC (R preprocessing) | `conda activate mariner_env` — rtracklayer, GenomicRanges |
| Local (Python) | Python 3.8+, pandas, numpy, scipy, matplotlib, openpyxl |
| Local (R) | R 4.x, ggplot2, patchwork, ggpointdensity, GenomicRanges, scales, clusterProfiler, org.Mm.eg.db, data.table |

## Running the Pipeline

Steps 0-5 and 9 run on HPC via SLURM. Steps 6-8 and 10+ run locally.

### Reference generation (one-time)

```bash
# Generate mm10 gene annotations from GENCODE vM25
sbatch scripts/step0a_gene_annotations.sb

# Convert consensus ATAC peaks to narrowPeak format
bash scripts/step0e_consensus_to_narrowpeak.sh

# Generate ubiquitous gene list (top 5% baseMean)
python scripts/step0d_ubiquitous_genes.py

# Fix ATAC tagAlign compression for ABC input (bgzip + tabix)
sbatch scripts/step0c_fix_tagalign.sb
```

### Core ABC pipeline

```bash
# Run ABC model (Snakemake, ~40 hr wall time)
sbatch scripts/step4_run_abc.sb

# QC: column headers, prediction counts, score distributions
sbatch scripts/step5_qc.sb
```

### Delta computation and integration (local)

```bash
# Compute ΔABC and Δ(A×C) between KO and WT
python scripts/step6_delta_abc.py

# Merge with RNA-seq DE, gene-level aggregation, concordance analysis
python scripts/step7_rnaseq_integration.py
```

### Cross-referencing (HPC)

```bash
# Annotate with Hi-C loops (bedtools closest) and H3K27ac peaks
python scripts/step9_cross_reference.py

# Extract per-replicate K119ub signal at enhancers
sbatch scripts/preprocess_k119ub_enhancer_signal.sb
```

### Downstream analysis (local, from `abc/` directory)

```bash
cd abc

# Paired-anchor loop-ABC overlap analysis
Rscript scripts/step9b_paired_anchor_analysis.R

# K119ub-ABC correlation (DiffBind-based, 12 panels)
Rscript scripts/step10_k119ub_abc_correlation.R

# Enhancer phenotypic subset stratification
Rscript scripts/step11_enhancer_subset_analysis.R

# HOMER motif enrichment (HPC)
sbatch scripts/step11_homer_motif_enrichment.sb

# HOMER motif visualization (local)
Rscript scripts/step11b_homer_motif_visualization.R

# Activity vs contact decomposition
Rscript scripts/step12_activity_contact_scatter.R
Rscript scripts/step12b_promoter_distal_scatter.R

# Discordant gene characterization
Rscript scripts/step13_discordant_gene_analysis.R
Rscript scripts/step13b_go_enrichment.R
Rscript scripts/step13c_k119ub_concordance.R

# Loop type at enhancer subsets
Rscript scripts/step14_loop_type_at_enhancer_subsets.R

# ABC component decomposition vs expression (Jesse Dixon suggestion)
Rscript scripts/step15_abc_component_vs_expression.R

# K119ub filtered scatter (FDR < 0.05 DiffBind sites)
Rscript scripts/step16_k119ub_filtered_scatter.R
```

## Key Design Decisions

### Consensus enhancer universe

Both WT and KO ABC runs use the **same** candidate element set: `consensus_all.bed` (75,371 ATAC peaks from Batch 1 union across both conditions). This is symlinked into both biosamples' peak directories so that ΔABC reflects true changes in enhancer-gene linkage strength, not artifacts from different enhancer definitions between conditions.

### Activity metric

Enhancer activity uses the geometric mean of ATAC-seq and H3K27ac: `sqrt(ATAC.RPM x H3K27ac.RPM)`. ATAC provides accessibility; H3K27ac distinguishes active enhancers from merely open chromatin. Each condition has 1 merged ATAC tagAlign + 5 H3K27ac BAM replicates.

### Quantile normalization disabled

The ABC default quantile normalization reference is derived from K562 DHS (human). This is inappropriate for mm10 ATAC-seq data, so `use_qnorm: False` in the Snakemake config.

### Unnormalized Δ(A×C) vs normalized ΔABC

ABC normalization divides by the genome-wide sum of (Activity x Contact), which compresses real activity changes when BAP1-KO causes widespread chromatin remodeling. The unnormalized Δ(A×C) shows stronger concordance with DE (65.3% vs 58.8%). Both metrics are computed and carried through all analyses.

## Thresholds

| Parameter | Value | Context |
|---|---|---|
| ABC score threshold | 0.02 | At least one condition must reach 0.02 (Step 6 filter) |
| ΔABC significance | 0.01 | \|ΔABC\| > 0.01 classifies gained/lost (Steps 7, 9) |
| RNA-seq DE | padj < 0.05, \|log2FC\| > 0.5 | Concordance with ΔABC direction |
| DiffBind K119ub FDR | 0.05 | Step 16 filtered scatter |

## Enhancer Phenotypic Classes

Enhancers are classified into four categories based on DiffBind results across three ChIP marks (K119ub, H3K27ac, ATAC):

| Class | Count | Definition |
|---|---|---|
| Activity_Lost | 7,503 | Significant change in H3K27ac and/or ATAC (activity markers) |
| K119ub_Only | 2,479 | Significant K119ub change but no activity change — key biological class |
| Activity_Gain | 2,851 | Gained activity markers |
| Stable | 42,864 | No significant change in any mark |

The **K119ub_Only** class is the most biologically interesting: these enhancers gain ubiquitination without detectable activity changes, testing whether K119ub alone is sufficient to alter 3D contact or downstream gene expression.

## Output Files

### Core results (`results/`)

| File | Rows | Description |
|---|---|---|
| `delta_abc_all_pairs.tsv` | 179,709 | All distal E-G pairs with ΔABC, Δ(A×C), per-condition activity and contact |
| `delta_abc_annotated.tsv` | 179,709 | Same + H3K27ac overlap flags from Step 9 |
| `delta_abc_significant.tsv` | 86,629 | Pairs where \|ΔABC\| > 0.01 |
| `delta_abc_with_rnaseq.tsv` | 113,943 | E-G pairs merged with RNA-seq log2FC and padj |
| `gene_level_summary.tsv` | 13,749 | One row per gene: strongest enhancer, aggregate stats, DE results |
| `k119ub_enhancer_signal.tsv` | 37,431 | Per-replicate K119ub signal at enhancers (HPC preprocessing) |
| `k119ub_abc_enhancer_merged.tsv` | 35,777 | Enhancer-level K119ub + ABC deltas joined |
| `k119ub_abc_correlation_summary.tsv` | 5 | Correlation statistics across stratifications |
| `discordant_gene_characteristics.tsv` | 957 | Per-gene stats for concordance classification |
| `loops_with_gene_assignments.tsv` | 2,910 | Differential loops with nearest-gene annotations |
| `paired_anchor_matches.tsv` | 513 | Loop-ABC pairs physically connected at both anchors |
| `paired_anchor_all_matches.tsv` | 17,020 | Paired-anchor matches using all 39K loops |

### Figures (`results/figures/`)

| Directory | Description |
|---|---|
| `k119ub_correlation/` | 12-panel K119ub vs ΔABC correlation (Step 10) |
| `k119ub_filtered_scatter/` | FDR-filtered K119ub scatter (Step 16) |
| `activity_contact_scatter/` | Activity vs contact decomposition (Steps 12, 12b) |
| `abc_component_expression/` | Component correlations with gene expression (Step 15) |
| `discordant_analysis/` | Discordant gene characterization + GO/KEGG (Steps 13, 13b, 13c) |
| `concordance_pie/` | Concordance pie charts (Step 7) |

### Enhancer subset results (`results/enhancer_subset_analysis/`)

Stratified analyses from Step 11: per-class loop metrics, ABC score distributions, contact profiles, and HOMER motif enrichment results.

### Loop type results (`results/loop_type_at_subsets/`)

Step 14 output: 7-category chromatin state annotation of loop partner anchors at each enhancer class.

## Reference Files

All in `reference/`:

| File | Description |
|---|---|
| `mm10_genes.bed` | 28,230 genes from GENCODE vM25 (protein_coding, lincRNA, processed_transcript) |
| `mm10_tss.bed` | 500 bp TSS regions (same gene set) |
| `mm10_merged_blacklist.bed` | Lab + ENCODE blacklists merged (3,581 regions) |
| `mm10_ubiquitous_genes.tsv` | Top 5% by baseMean (730 genes, used by ABC for normalization) |

## Project Structure

```
abc/
├── ABC-Enhancer-Gene-Prediction/  # Broad Institute ABC model (Snakemake pipeline)
├── scripts/
│   ├── step0a_gene_annotations.sb       # GENCODE vM25 → mm10 gene/TSS BED
│   ├── step0c_fix_tagalign.sb           # bgzip + tabix ATAC tagAligns
│   ├── step0d_ubiquitous_genes.py       # Top 5% expressed genes
│   ├── step0e_consensus_to_narrowpeak.sh # Consensus ATAC peaks → narrowPeak
│   ├── step4_run_abc.sb                 # Snakemake ABC pipeline
│   ├── step5_qc.sb                      # QC column headers + distributions
│   ├── step6_delta_abc.py               # ΔABC computation (outer join)
│   ├── step7_rnaseq_integration.py      # RNA-seq merge + concordance
│   ├── step9_cross_reference.py         # Loop + H3K27ac annotation
│   ├── step9b_paired_anchor_analysis.R  # Loop-ABC physical overlap
│   ├── preprocess_k119ub_enhancer_signal.R   # HPC: K119ub bigWig extraction
│   ├── preprocess_k119ub_enhancer_signal.sb  # SLURM wrapper
│   ├── step10_k119ub_abc_correlation.R  # K119ub-ABC DiffBind correlation
│   ├── step11_enhancer_subset_analysis.R     # Enhancer class stratification
│   ├── step11_homer_motif_enrichment.sb      # HOMER motif enrichment (HPC)
│   ├── step11b_homer_motif_visualization.R   # HOMER result visualization
│   ├── step12_activity_contact_scatter.R     # Activity vs contact scatter
│   ├── step12b_promoter_distal_scatter.R     # Promoter vs distal stratification
│   ├── step13_discordant_gene_analysis.R     # Discordant gene characterization
│   ├── step13b_go_enrichment.R               # GO/KEGG: concordant vs discordant
│   ├── step13c_k119ub_concordance.R          # K119ub by concordance status
│   ├── step14_loop_type_at_enhancer_subsets.R # Loop types at enhancer classes
│   ├── step15_abc_component_vs_expression.R  # Component decomposition vs DE
│   └── step16_k119ub_filtered_scatter.R      # FDR-filtered K119ub scatter
├── enhancer_subsets/          # 4 enhancer class TSVs (from DiffBind classification)
├── reference/                 # mm10 gene annotations, blacklists, ubiquitous genes
├── results/                   # All output TSVs and figures
│   ├── figures/               # Visualization outputs (PDF + SVG + JPG)
│   ├── enhancer_subset_analysis/  # Step 11 stratified results + HOMER
│   ├── loop_type_at_subsets/      # Step 14 loop decomposition
│   └── paired_anchor_plots/       # Step 9b visualizations
├── docs/                      # Methods, execution plan, analysis context
└── logs/                      # SLURM job logs
```

## Documentation

| Document | Description |
|---|---|
| `docs/abc-prd-v2.md` | Full product requirements document |
| `docs/abc-execution-plan.md` | Step-by-step implementation plan (~1,600 lines) |
| `docs/abc-analysis-context.md` | QC findings and results from Steps 5-9 |
| `docs/methods.md` | ABC model scientific methodology (from Broad) |
| `docs/step10-k119ub-abc-correlation.md` | Step 10 analysis documentation |
| `docs/step11-enhancer-subset-analysis.md` | Enhancer subset analysis documentation |
| `docs/step9b-paired-anchor-analysis.md` | Paired-anchor loop analysis documentation |
