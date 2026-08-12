# Steps 5-7 & 9: ABC Pipeline QC, Delta-ABC, RNA-seq Integration, and Cross-Reference

## Position in the Pipeline

These steps form the foundation of the ABC enhancer-gene linkage analysis. Everything downstream (Steps 9b, 10, 11) depends on the outputs generated here:

- **Step 5** (QC): Validates the ABC pipeline ran correctly on both conditions
- **Step 6** (Delta-ABC): Computes differential ABC scores between BAP1-KO and WT
- **Step 7** (RNA-seq Integration): Merges delta-ABC with gene expression data and tests concordance
- **Step 9** (Cross-Reference): Overlaps ABC enhancers with Hi-C differential loops and H3K27ac peaks

### What Came Before

Steps 0-4 ran on HPC (Expanse) and are not documented here:

- **Step 0**: Reference generation — mm10 gene annotations, TSS bed, blacklists, ubiquitous genes, consensus ATAC peaks
- **Step 0e**: Convert consensus ATAC peaks to narrowPeak format for ABC input
- **Step 4**: Run the Broad Institute ABC model via Snakemake for both WT and KO conditions
- Both conditions used the **same** `consensus_all.bed` (75,371 ATAC peaks from Batch 1 union), ensuring delta-ABC measures true E-G linkage changes, not artifacts from different enhancer definitions

---

## Step 5: QC Validation

### Pipeline Configuration

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Organism | Mouse mm10, GENCODE vM25 | - |
| Tool | ABC v1.1.2 (Snakemake) | - |
| Mode | ATAC-only | H3K27ac BAMs unavailable; only bigWigs exist |
| Hi-C | Condition-specific .hic at 5kb | Cell-type specific (best mode per ABC docs) |
| Consensus peaks | 75,371 | Union of Batch 1 ctrl/mut ATAC, >=2/3 reps per genotype |
| genome_size | mm | MACS2 flag for mouse (not hg38 default) |
| use_qnorm | False | Reference is human K562 DHS — cross-species and cross-assay |
| threshold | 0.02 | Set explicitly because disabling qnorm invalidates auto-calibration |

### Reference File Statistics

| File | Count | Notes |
|------|-------|-------|
| mm10_genes.bed | 28,230 | 21,833 protein_coding, 5,620 lincRNA, 777 processed_transcript |
| mm10_tss.bed | 28,230 TSS | All 500bp, matches gene count |
| mm10_merged_blacklist.bed | 3,581 | Lab (254) + ENCODE (3,435), merged |
| mm10_ubiquitous_genes.tsv | 730 | Top 5% by baseMean from RNA-seq |
| consensus_all.narrowPeak | 75,371 | 20 standard chromosomes (chr1-19, chrX) |

Spot-checks: Actb present in ubiquitous list (rank 575). Gapdh absent (rank 1,054, just below 95th percentile cutoff — not a bug). Rpl13a absent (rank 11,925, unexpectedly low — likely ribosomal gene mapping issue). The ubiquitous list is used only for flagging, not filtering.

### Per-Condition Predictions (ABC >= 0.02)

| Metric | WT | KO |
|--------|----|----|
| Total E-G pairs | 122,212 | 119,965 |
| Self-promoter entries | 29,152 | 29,152 |
| Distal E-G pairs | 93,060 | 90,813 |
| Unique genes | 27,492 | 27,492 |
| Mean distal enhancers/gene | 3.38 | 3.30 |
| Median ABC score (distal) | 0.0343 | 0.0344 |
| Median distance (distal) | 127,108 bp | 115,795 bp |
| Chromosomes | 20 | 20 |

**100% gene overlap** between conditions (27,492 genes in both) — confirms consensus peak injection worked correctly. **2,247 fewer E-G pairs in KO** (1.8% reduction) — modest net loss of enhancer-gene connections crossing threshold.

### Hi-C Sparsity

| Distance bin | Zero hic_contact (WT) |
|--------------|-----------------------|
| < 50 kb | 0.1% |
| 50-200 kb | 1.1% |
| 200-500 kb | 7.0% |
| 500 kb-1 Mb | 22.0% |
| 1-2 Mb | 43.4% |
| 2-5 Mb | 65.1% |
| **Overall** | **47.3% (WT) / 49.1% (KO)** |

The headline number (~48% zeros) sounds alarming but is dominated by long-range pairs. Among thresholded distal pairs (ABC >= 0.02), only 3.1-3.4% have zero Hi-C contact. **98.7% of pairs that matter have real Hi-C data** in at least one condition. Hi-C sparsity is not a concern for the differential analysis.

### QC Readiness Check

All 10 checks passed: reference files exist, column headers identical between WT/KO, 20 chromosomes in each, biotype filter clean. Ready for Step 6.

---

## Step 6: Delta-ABC Computation

### Join Statistics

- AllPutative pairs in both conditions: 7,942,829 (92.9% overlap)
- WT only: 354,726 / KO only: 256,770
- After threshold (ABC >= 0.02 in either condition): **180,423 distal E-G pairs**
- Unique genes: 20,265

### Delta-ABC Distribution

| Statistic | Value |
|-----------|-------|
| Mean | -0.000459 |
| Median | -0.000042 |
| Std | 0.031158 |
| Q5 / Q95 | -0.037 / +0.035 |
| Gained (delta-ABC > 0.01) | 44,060 (24.4%) |
| Lost (delta-ABC < -0.01) | 47,499 (26.3%) |
| Unchanged (\|delta\| <= 0.01) | 88,864 (49.3%) |

Slightly negative mean indicates a modest net loss of ABC scores in KO. 7.8% excess of lost over gained links.

### Distance Dependence

| Distance | N | Mean delta | Lost/Gained ratio |
|----------|---|-----------|-------------------|
| < 50 kb | 44,234 | -0.0019 | 1.02 |
| 50-200 kb | 47,044 | -0.0025 | 1.33 |
| 200-500 kb | 32,060 | -0.0030 | 1.41 |
| 500 kb-1 Mb | 18,721 | -0.0035 | 1.40 |
| 1-5 Mb | 31,720 | -0.0032 | 1.27 |

Loss of E-G connections preferentially affects **mid-to-long range interactions** (50 kb-1 Mb), with a lost/gained ratio up to 1.41. Short-range (<50 kb) is near-symmetric (ratio 1.02). This is biologically consistent: BAP1-KO disrupts long-range enhancer-gene connections mediated by chromatin loops while leaving proximal interactions (which are less dependent on 3D structure) relatively intact.

---

## Step 7: RNA-seq Integration

### Gene Matching

| Metric | Value |
|--------|-------|
| ABC genes | 20,265 |
| RNA-seq genes | 16,572 |
| Matched | 13,588 (67.1% of ABC) |
| E-G pairs with RNA-seq | 113,453 |

Unmatched ABC genes are primarily lincRNAs and processed transcripts not in the DESeq2 output.

### Gene-Level Concordance — The Core RNA-seq Result

Concordance tests whether the direction of enhancer-gene linkage change (delta-ABC) matches the direction of gene expression change (log2FC). Tested at the gene level (one observation per gene, using the strongest delta-ABC enhancer) to avoid pseudo-replication.

**Normalized delta-ABC:**

| Metric | Value |
|--------|-------|
| Genes with \|delta-ABC\| > 0.01 | 12,665 |
| Genes with DE (padj < 0.05, \|LFC\| > 0.5) | 978 |
| Genes with both | 940 |
| **Concordant** | **553 (58.8%)** |
| **Binomial p-value** | **6.84e-08** |

**Unnormalized delta(AxC):**

| Metric | Value |
|--------|-------|
| Genes with both changed | 911 |
| **Concordant** | **595 (65.3%)** |
| **Binomial p-value** | **1.64e-20** |

The unnormalized score substantially outperforms the normalized score (65.3% vs 58.8%). This is the first indication that ABC normalization is counterproductive for differential analysis.

### 2x2 Directional Breakdown (Normalized)

|  | Upregulated | Downregulated |
|---|---|---|
| Gained enhancer | 154 | 252 |
| Lost enhancer | 135 | 399 |

Lost+downregulated (399) is the dominant quadrant (42%). Total downregulated: 651. Total upregulated: 289. Ratio 2.25:1. The asymmetry toward downregulation is consistent with BAP1-KO predominantly disrupting enhancer function rather than activating it.

### Scatter Plots — Strongest Enhancer per Gene

The `delta_abc_vs_log2fc.png` figure contains two panels showing delta-ABC (or delta-unnorm) vs RNA-seq log2FC for 13,588 genes, with DE genes (padj < 0.05) highlighted.

**Left panel (Normalized delta-ABC):** Pearson r = 0.035 (p = 0.13, **non-significant**). Spearman rho = 0.092. The point cloud is compressed into a tight vertical stripe near delta-ABC = 0. The normalization is destroying the signal. This makes biological sense: when BAP1-KO causes widespread activity changes across many enhancers for the same gene, the denominator (sum of AxC for all elements within 5 Mb) shifts proportionally, so the normalized ratio barely moves even when the absolute regulatory input changes substantially.

**Right panel (Unnormalized delta(AxC)):** Pearson r = 0.269 (p ~ 10^-32). Spearman rho = 0.300 (p ~ 10^-40). The diagonal trend is visible — genes gaining enhancer activity tend to be upregulated, and vice versa. Spearman > Pearson indicates a monotonic but nonlinear relationship, expected since ABC scores and log2FC operate on different scales. An r ~ 0.3 for a single-enhancer-to-expression correlation is solid for a chromatin model — ~9% of variance in DE explained from enhancer activity changes alone.

### Scatter Plots — Sum of All Enhancers per Gene

The `sum_delta_abc_vs_log2fc.png` tests whether summing delta-ABC across all enhancers per gene (rather than taking just the strongest) improves concordance.

**Normalized sum:** Pearson r = 0.031 (dead). Spearman rho = 0.114 (barely improved from 0.092). The normalization forces sum-ABC ~ 1 per gene by design, so sum-delta-ABC ~ 0 by construction. The signal is mathematically suppressed.

**Unnormalized sum:** Pearson r = 0.451 (p ~ 10^-92). **Spearman rho = 0.582 (p ~ 10^-172)**. This is the strongest correlation in the entire analysis. The jump from rho = 0.300 (strongest enhancer) to rho = 0.582 (sum) means the aggregate captures ~34% of rank variance, up from ~9%.

**Biological interpretation:** BAP1-KO doesn't just flip one enhancer per gene on/off — it remodels the entire regulatory landscape around affected genes. Summing delta(AxC) across all enhancers captures the cumulative regulatory input change, and that aggregate shift predicts expression direction far better than any single enhancer.

### Honest Assessment

The unnormalized-sum result (rho = 0.582) is genuinely strong for a chromatin model. But several caveats apply:

1. **Shared ATAC data.** The ABC activity component is ATAC signal. If ATAC and RNA-seq correlate for reasons unrelated to enhancer function (e.g., both reflect chromatin accessibility at promoters), some concordance is expected by construction. The contact component of ABC is from Hi-C (independent), but for short-range pairs, the powerlaw contribution is large.

2. **Many "strongest enhancers" are promoters.** Several top dysregulated genes have the strongest-enhancer classified as "promoter" with small distances (Cbr1 at 1,789bp, Kctd15 at 1,152bp). These are consensus ATAC peaks overlapping the promoter itself. The truly interesting hits for the BAP1 story are long-range connections like Col15a1 (4.8 Mb, intergenic, log2FC = 1.82) and Robo1 (2 Mb, 15 enhancers, log2FC = -1.74).

3. **Normalization is not wrong, just inappropriate here.** Per-gene normalization is designed for single-condition enhancer-gene mapping (identifying which enhancers regulate which genes). For differential analysis where a perturbation causes global chromatin remodeling, it destroys the signal. This is a methodological finding worth noting in the paper.

---

## Step 9: Cross-Reference with Loops and H3K27ac

### 9A: Differential Loop Overlap

- Total differential loops: 2,910
- ABC enhancers overlapping anchor1: 1,394/33,004 (4.2%)
- ABC enhancers overlapping anchor2: 1,330/33,004 (4.0%)
- ABC target genes in loop annotations: 2,433/20,265 (12.0%)
- **Directional concordance (loop vs delta-ABC): 1,209/2,352 (51.4%)** — essentially at chance

**Why 51.4% is not a failure:** Loop annotations use nearest-gene assignment (which gene is closest to the anchor), while ABC uses activity x contact scoring (which enhancer-gene pair has the strongest regulatory linkage). These are fundamentally different ways of assigning genes to genomic regions. The chance-level concordance reflects this methodology mismatch, not a failure of either dataset. Step 9b later resolves this by requiring the enhancer at one anchor and the gene TSS at the other, achieving 89.7% concordance.

### 9B: H3K27ac Peak Overlap

| Category | Enhancers | E-G pairs | Mean \|delta-ABC\| |
|----------|-----------|-----------|-------------------|
| H3K27ac+ | 12,993 (39.4%) | 98,440 | 0.01805 |
| H3K27ac- | 20,011 (60.6%) | 81,983 | 0.01598 |

**H3K27ac-stratified gain/loss asymmetry:**

| | Gained | Lost | Lost/Gained | Mean delta-ABC |
|---|---|---|---|---|
| **H3K27ac+** | 24,346 | 28,409 | **1.17** | **-0.00102** |
| **H3K27ac-** | 19,714 | 19,089 | **0.97** | +0.00021 |

The gain/loss asymmetry is concentrated at H3K27ac-marked enhancers (1.17 lost/gained ratio). H3K27ac- enhancers are essentially symmetric. This is consistent with BAP1-KO preferentially disrupting active enhancers (where deubiquitination was functionally relevant) rather than globally affecting all enhancers.

**Important caveat:** H3K27ac peaks are from a single merged adult timepoint dataset, not condition-specific DiffBind. The annotation identifies enhancers with active marks in the reference, not whether those marks changed between WT and KO. Step 11 addresses this with DiffBind-derived classification.

---

## Summary of Claims and Evidence Strength

| Claim | Evidence | Strength |
|-------|----------|----------|
| ABC pipeline ran correctly for both conditions | 100% gene overlap, 10/10 QC checks passed | **Strong** |
| Net loss of E-G connections in BAP1-KO | 7.8% excess lost vs gained links | **Moderate** — real but small |
| Loss preferentially affects mid-to-long range | Lost/gained ratio 1.02 (<50kb) to 1.41 (200-500kb) | **Strong** — consistent gradient |
| Unnormalized delta(AxC) outperforms normalized delta-ABC | 65.3% vs 58.8% concordance; rho 0.300 vs 0.092 | **Strong** — consistent across all analyses |
| Sum across enhancers outperforms strongest enhancer | rho 0.582 vs 0.300 (unnormalized) | **Strong** — 34% vs 9% variance explained |
| Per-gene normalization is inappropriate for differential ABC | Pearson r = 0.035 (non-significant) for normalized; 0.269 for unnormalized | **Strong** — mathematically expected |
| 940 dysregulated genes are concordant above chance | Binomial p = 6.84e-08 | **Strong** statistically |
| Active (H3K27ac+) enhancers are preferentially affected | 1.17 lost/gained ratio at H3K27ac+ vs 0.97 at H3K27ac- | **Moderate** — consistent but modest |
| Independent loop-ABC overlap is at chance | 51.4% concordance | **Strong** null — reflects methodology mismatch |

---

## Recommended Figures for Presentation

### Lead with the Sum Scatter (Right Panel of `sum_delta_abc_vs_log2fc.png`) — "The Headline Result"

The unnormalized-sum scatter (Spearman rho = 0.582) is the strongest single correlation in the entire ABC pipeline. The diagonal structure is visually obvious — blue DE points track from lower-left (lost enhancer input, downregulated) to upper-right (gained enhancer input, upregulated). This is the figure that establishes the ABC model *works* for predicting differential gene expression in BAP1-KO, and it does so using the cumulative enhancer landscape, not just the single strongest link.

### Follow with the 4-Panel Comparison — "The Methodological Insight"

Show both `delta_abc_vs_log2fc.png` (strongest enhancer) and `sum_delta_abc_vs_log2fc.png` (sum) as a 2x2 comparison: normalized vs unnormalized (columns) x strongest vs sum (rows). The visual progression from upper-left (dead signal, r = 0.035) to lower-right (strong signal, rho = 0.582) tells the methodological story: (a) normalization kills the signal, and (b) summing across enhancers captures the cumulative regulatory change. This is actually a useful contribution to the field — most ABC papers use normalized scores for single-condition mapping, not differential comparisons.

### Why Not the Others?

- The QC log results are essential background but not presentation figures — mention "10/10 QC checks passed" verbally
- The Hi-C sparsity data is a robustness control — mention "98.7% of thresholded pairs have real Hi-C data" to preempt reviewer questions
- The Step 9 cross-reference (51.4% concordance) is important context for Step 9b but should be presented as the "before" in the 9b concordance bar chart (Panel 1), not as a standalone figure

### The Narrative

> "We computed differential ABC scores for 180,423 enhancer-gene pairs between BAP1-KO and WT. Per-gene normalization mathematically suppresses the signal (left panels, r ~ 0.03), but the unnormalized sum of Activity x Contact across all enhancers per gene predicts expression direction with Spearman rho = 0.582 (lower-right panel). BAP1-KO remodels the entire regulatory landscape around affected genes, and this cumulative change predicts which genes are up- or downregulated."

---

## Limitations

1. **ATAC-only ABC mode.** Activity is ATAC signal only (H3K27ac BAMs unavailable). ABC benchmarks best with DNase + H3K27ac. ATAC-only may underestimate enhancer activity, particularly where H3K27ac and ATAC diverge. The ABC documentation notes bulk ATAC-seq performs worse than DNase-seq in CRISPR benchmarks.

2. **QNorm disabled.** The quantile normalization reference is human K562 DHS — inappropriate for mouse ATAC. Disabling qnorm means ABC scores are not calibrated to the CRISPR benchmark, and the 0.02 threshold is set manually rather than by auto-calibration. Relative comparisons (delta-ABC) are unaffected, but absolute ABC scores should not be compared to other studies that use qnorm.

3. **~8% of genes have 2+ self-promoter entries** due to overlapping consensus ATAC peaks at promoter regions. This slightly inflates their ABC denominators, deflating distal enhancer scores. Effect cancels in delta-ABC (same peaks in both conditions).

4. **Hi-C sparsity at long range.** ~65% of 2-5 Mb pairs have zero raw Hi-C contact. For these, ABC falls back to powerlaw pseudocount (identical between conditions), so long-range delta-ABC reflects only delta-activity, not delta-contact. This primarily affects the AllPutative universe; 98.7% of thresholded pairs are fine.

5. **RNA-seq is one timepoint.** Adult cerebellum steady-state. If BAP1-KO effects on gene expression are developmental, this misses them. Step 11's concordance failure for K119ub_Only enhancers may partly reflect this.

6. **Gene matching is 67.1%.** One-third of ABC genes lack RNA-seq data (mostly non-coding). This limits the concordance analysis to the protein-coding-biased subset.

---

## Output File Descriptions

### Figures

| File | Content |
|------|---------|
| `figures/delta_abc_vs_log2fc.pdf/png` | 2-panel scatter: strongest enhancer delta-ABC (normalized, left) and delta-unnorm (right) vs RNA-seq log2FC. DE genes (padj < 0.05) highlighted. Pearson/Spearman stats annotated. |
| `figures/sum_delta_abc_vs_log2fc.pdf/png` | 2-panel scatter: sum of all enhancers' delta-ABC (normalized, left) and delta-unnorm (right) vs log2FC. Contains the headline rho = 0.582 result. |
| `figures/interpretation.md` | Original analysis interpretation with figure recommendations |
| `figures/k119ub_correlation/` | Step 10 K119ub-ABC correlation (8 panels) — documented separately in `step10-k119ub-abc-correlation.md` |

### TSV Files

| File | Rows | Description |
|------|------|-------------|
| `delta_abc_all_pairs.tsv` | 180,423 | All thresholded distal E-G pairs with delta-ABC, delta-unnorm, activity/contact components, Hi-C flags |
| `delta_abc_significant.tsv` | 91,559 | Subset with \|delta-ABC\| > 0.01 |
| `delta_abc_with_rnaseq.tsv` | 113,453 | E-G pairs merged with RNA-seq (log2FC, padj, baseMean) |
| `gene_level_summary.tsv` | 13,588 | One row per gene: strongest delta-ABC enhancer, aggregate stats (sum, mean, counts), DE results, dysregulated flag |
| `delta_abc_annotated.tsv` | 180,423 | E-G pairs annotated with `has_H3K27ac` flag from Step 9 |
| `loops_with_gene_assignments.tsv` | - | Loops with nearest-gene annotations from Step 9 |

### Key Column Inventory for `delta_abc_all_pairs.tsv`

- `chr`, `start`, `end` — enhancer coordinates
- `TargetGene` — gene symbol
- `distance` — enhancer-to-TSS distance
- `class` — genic/intergenic/promoter
- `ABC.Score_WT`, `ABC.Score_KO` — normalized scores per condition
- `ABC.Score.Numerator_WT`, `ABC.Score.Numerator_KO` — unnormalized AxC
- `activity_base_WT`, `activity_base_KO` — ATAC activity component
- `hic_contact_pl_scaled_adj_WT`, `hic_contact_pl_scaled_adj_KO` — adjusted Hi-C contact
- `delta_ABC`, `delta_unnorm` — KO minus WT deltas
- `hic_data_wt`, `hic_data_ko`, `has_hic_either` — Hi-C coverage flags

---

## Bottom Line

**What the data shows:** The ABC model successfully identifies enhancer-gene linkage changes in BAP1-KO. The unnormalized sum of Activity x Contact across all enhancers per gene predicts expression direction with rho = 0.582 — a strong result for a chromatin model. Loss of enhancer-gene connections preferentially affects active (H3K27ac+) enhancers and mid-to-long range interactions, consistent with BAP1 maintaining active enhancer function.

**What the data does not show:** Per-gene normalized ABC scores are ineffective for differential analysis (r = 0.035, non-significant). The ABC normalization, designed for single-condition enhancer-gene mapping, is mathematically counterproductive when a perturbation causes global chromatin remodeling. Independent loop-ABC overlap is at chance (51.4%), but this reflects methodology mismatch (nearest-gene vs AxC scoring), not model failure.

**Interpretation:** The ABC pipeline produces a reliable enhancer-gene linkage map for differential analysis, provided the unnormalized score is used. The 180,423 delta-ABC pairs and 13,588 gene-level summaries form the input for all downstream analyses (Steps 9b, 10, 11), where the findings are progressively refined by integrating Hi-C loop geometry, H2AK119ub signal, and epigenomic classification.
