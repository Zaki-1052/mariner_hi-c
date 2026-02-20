# Research Compilation: Regulation of Chromatin Conformation by the Histone Deubiquitinase BAP1 in the Brain

**Author:** Zakir Alibhai
**Period:** ~20 weeks (October 2025 - February 2026)
**Affiliation:** Ferguson Lab, UCSD
**Organism:** Mouse (mm10 genome)
**Conditions:** BAP1-KO (CRISPR deletion) vs Wildtype Control
**Tissue:** Cerebellum
**Timepoints:** Early (P12) and Late (P60/Adult)
**Replicates:** n=3 biological replicates per condition per timepoint

---

## Table of Contents

1. [Project Overview & Biological Context](#1-project-overview--biological-context)
2. [Hi-C Data Processing Pipeline (HPC)](#2-hi-c-data-processing-pipeline-hpc)
3. [Differential Chromatin Loop Analysis (Mariner/edgeR)](#3-differential-chromatin-loop-analysis-marineredger)
4. [Loop Annotation & ChIP-seq Integration](#4-loop-annotation--chip-seq-integration)
5. [Loop Distance Shift & "Loop Rewriting" Phenomenon](#5-loop-distance-shift--loop-rewriting-phenomenon)
6. [Shared Anchor Analysis (Loop Switching)](#6-shared-anchor-analysis-loop-switching)
7. [H2AK119ub Integration (BAP1 Substrate)](#7-h2ak119ub-integration-bap1-substrate)
8. [Differential ChIP-seq Polycomb Enrichment](#8-differential-chip-seq-polycomb-enrichment)
9. [DEG-Loop Anchor Integration](#9-deg-loop-anchor-integration)
10. [TAD Boundary Analysis (TADCompare)](#10-tad-boundary-analysis-tadcompare)
11. [Chromatin Compartment (PC1) Analysis](#11-chromatin-compartment-pc1-analysis)
12. [Chromatin Stripe Analysis (Quagga)](#12-chromatin-stripe-analysis-quagga)
13. [CTCF-Stripe Cross-Reference](#13-ctcf-stripe-cross-reference)
14. [ABC Enhancer-Gene Linkage Model](#14-abc-enhancer-gene-linkage-model)
15. [Machine Learning: LoopBin Unsupervised Clustering](#15-machine-learning-loopbin-unsupervised-clustering)
16. [Biomodal Methylation Analysis (5mC/5hmC)](#16-biomodal-methylation-analysis-5mc5hmc)
17. [ATAC-seq Consensus Peak Analysis](#17-atac-seq-consensus-peak-analysis)
18. [HOMER Motif Analysis](#18-homer-motif-analysis)
19. [Aggregate Peak Analysis (APA)](#19-aggregate-peak-analysis-apa)
20. [Hi-C Processing Infrastructure (HPC)](#20-hi-c-processing-infrastructure-hpc)
21. [Technical Innovations & Methodological Contributions](#21-technical-innovations--methodological-contributions)
22. [Complete Results Summary Table](#22-complete-results-summary-table)
23. [Key Biological Findings & Unified Model](#23-key-biological-findings--unified-model)
24. [Publication-Ready Outputs](#24-publication-ready-outputs)
25. [Tools & Technologies Used](#25-tools--technologies-used)

---

## 1. Project Overview & Biological Context

### The Biological Question

**BAP1** (BRCA1-Associated Protein 1) is a histone deubiquitinase that removes H2AK119ub (monoubiquitination of histone H2A at lysine 119), a Polycomb Repressive Complex 1 (PRC1)-deposited mark. BAP1 is essential for proper chromatin regulation, and its loss leads to accumulation of H2AK119ub at regulatory elements. This project investigates how BAP1 loss reshapes three-dimensional chromatin architecture in the developing and adult mouse cerebellum.

### Central Hypothesis

BAP1 loss disrupts chromatin 3D organization at multiple scales -- loops, TADs, stripes, and compartments -- through accumulated H2AK119ub at regulatory elements, leading to a "loop rewriting" phenomenon where long-range developmental loops are lost and replaced by shorter-range ectopic contacts.

### Experimental Design

| Parameter | Details |
|-----------|---------|
| **Organism** | Mus musculus (mm10 genome) |
| **Tissue** | Cerebellum |
| **Genotypes** | Wildtype (Control) vs BAP1-KO (Mutant, CRISPR) |
| **Early timepoint** | P12 (postnatal day 12, designated 250831) |
| **Late timepoint** | P60/Adult (designated 250402) |
| **Biological replicates** | n=3 per condition per timepoint |
| **Samples** | ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3 |
| **Hi-C protocol** | Arima Hi-C |
| **Sequencing depth** | ~600M-1B reads per replicate |

### Epigenomic Data Available

| Mark | Type | Source |
|------|------|--------|
| H3K27ac | Active enhancer/promoter | CUT&RUN |
| H3K27me3 | Polycomb repression | CUT&RUN |
| H3K4me1 | Enhancer (poised/active) | CUT&RUN |
| H3K4me3 | Active promoter | CUT&RUN |
| H2AK119ub | PRC1 ubiquitination (BAP1 substrate) | CUT&RUN |
| CTCF | Structural insulator | CUT&RUN + DNA motifs |
| MeCP2 | Methyl-CpG binding | CUT&RUN |
| ATAC-seq | Chromatin accessibility | ATAC-seq (n=3-5/condition) |
| RNA-seq | Gene expression | Bulk RNA-seq (both timepoints) |
| Methylation | 5mC/5hmC | Biomodal DUET evoC |

---

## 2. Hi-C Data Processing Pipeline (HPC)

### Infrastructure & Computational Environment

All primary Hi-C processing was performed on **SDSC Expanse** (San Diego Supercomputer Center) via SLURM job scheduling.

### Pipeline Overview

**Two complementary pipelines were maintained:**

#### A. Juicer Pipeline (Shallow Sequencing QC)
- Used for initial shallow sequencing QC (<600M reads)
- Produces .hic files with VC, KR normalization
- QC metrics: alignment rate (80-90%), PCR duplicates (<20%), MAPQ filtering
- Enabled HiCCUPS loop calling, Arrowhead TAD calling, APA

#### B. Nextflow-HiC Pipeline (Deep Sequencing, Primary)
- Used for deep sequencing production runs (600M-1B reads)
- Based on HiC-Pro pipeline: bowtie2 alignment + cooler contact maps
- Outputs: cool/mcool files, ValidPairs, pairix files, compartment eigenvectors
- Includes fastp trimming, MultiQC reports, insulation analysis
- Runtime: ~25-40 hours per sample depending on read depth

### Post-Processing Steps

1. **fastp trimming** of adapter sequences before alignment
2. **Nextflow-HiC** pipeline execution with bowtie2 alignment
3. **format_4DN.sb** -- Format pairs for Juicer compatibility
4. **juicer_pre.sb** -- Generate Juicer v2.0 .hic files from pairs
5. **sortmerge_contacts.sh** -- Merge replicates for combined analysis
6. **HiCCUPS** -- Loop calling on individual and merged .hic files
7. **Cooler utilities** -- Balancing (ICE), zoomify, dump, info

### HOMER Tag Directory Generation

For differential TAD and compartment analysis:
1. Merged BAM files across batches (`merge_and_splitbam.sb`)
2. Split into R1/R2 alignment files (`unmergebam.sb`)
3. Generated tag directories for each sample (`maketagdir.sb`)
4. Used for downstream HOMER TAD/PC1/motif analysis

### Key Data Outputs

| Output | Location | Purpose |
|--------|----------|---------|
| .hic files (Juicer v2.0) | juicerpre/ | Loop calling, APA, straw extraction |
| .cool/.mcool files | Contact_maps/ | GENOVA, cooler analysis |
| ValidPairs | hicpro/Valid_pairs/ | Alternate contact map generation |
| BAM files | hicpro/Mapping/ | HOMER tag directories |
| HiCCUPS loops | hiccups_results/ | Input to mariner pipeline |
| MultiQC reports | Multiqc/ | QC metrics |

---

## 3. Differential Chromatin Loop Analysis (Mariner/edgeR)

### Pipeline Architecture

Developed a complete **11-step replicate-aware differential loop analysis pipeline** using the mariner R/Bioconductor package and edgeR statistical framework. This is the core analytical contribution of the project.

### The 11-Step Workflow

```
Input: HiCCUPS loop calls (BEDPE) + Hi-C contact matrices (.hic)
  |
[1] prep_loops.R        -> Merge 6 replicates -> consensus loop positions
[2] extract_counts.R    -> Extract 5x5 Hi-C matrices from .hic files
[3] aggregate.R         -> Sum matrices -> count matrix (loops x samples)
[4] qc-val.R            -> Quality control validation (9 QC diagnostics)
[5] edgeR.R             -> Quasi-likelihood GLM differential analysis
[6] compare_resolutions.R -> Cross-resolution comparison & overlap
[7] generate_final_results.py -> Stringent filters (|logFC|>0.3, FDR<0.03)
[8] convert_final_bedpe.sh -> BEDPE for Juicebox visualization
[9] downstream_analysis.R -> Merge resolutions, annotate with ChIP-seq
[10] visualizations.R    -> Publication plots, enrichment, volcano
[11] apa_analysis.R      -> Aggregate Peak Analysis heatmaps
  |
Output: Differential loops, annotations, statistics, visualizations
```

### Multi-Resolution Analysis

All analyses run at three resolutions to capture different loop types:

| Resolution | Best For | Total Tested | Significant (FDR<0.05) | Final (Stringent) |
|------------|----------|-------------|------------------------|-------------------|
| **5 kb** | Short-range E-P | 17,982 | 1,766 (9.8%) | 1,120 (6.2%) |
| **10 kb** | Intermediate | 22,632 | 3,981 (17.6%) | 1,532 (6.8%) |
| **25 kb** | Long-range, inter-TAD | 20,398 | 4,774 (23.4%) | 1,189 (5.8%) |

### Statistical Framework

- **Method:** Quasi-likelihood negative binomial GLM (edgeR)
- **Normalization:** TMM (Trimmed Mean of M-values)
- **Design:** `~group` (intercept + mutant effect)
- **Testing:** `glmQLFTest` for mutant vs control
- **Dispersion:** Data-driven (no fixed BCV), `robust=TRUE` for outlier protection
- **Multiple testing:** Benjamini-Hochberg FDR

### Key Quality Metrics (Late Timepoint, 5kb)

| Metric | Value |
|--------|-------|
| Control library size | 2,864,059 |
| Mutant library size | 3,820,723 |
| Common BCV | 0.038 |
| Median tagwise BCV | 0.019 |
| Within-group correlation | >0.95 |

### Merged Results (Late Timepoint)

After cross-resolution merging with overlap removal:

| Metric | Value |
|--------|-------|
| **Total non-redundant differential loops** | **2,910** |
| Up in BAP1-KO (gained/strengthened) | 1,723 (59.2%) |
| Down in BAP1-KO (lost/weakened) | 1,187 (40.8%) |
| From 5kb resolution | 509 |
| From 10kb resolution | 1,335 |
| From 25kb resolution | 1,066 |

### Early vs Late Timepoint Comparison

| Metric | Early (P12) | Late (P60) | Fold Change |
|--------|-------------|------------|-------------|
| Total differential loops | 165 | 2,910 | **17.6x** |
| % significant (FDR<0.05) | 0.4-1.4% | 9.8-23.4% | ~21x |
| Direction bias | 57% DOWN | 59% UP | Reversal |
| Dominant loop type | Repressed_Promoter | Mixed/diverse | Progression |

---

## 4. Loop Annotation & ChIP-seq Integration

### 8-Category Chromatin State Classification

Developed a **priority-based 8-category classification system** for loop anchors using 5 histone marks + CTCF motifs:

| Priority | Category | Definition |
|----------|----------|-----------|
| 1 | **Active_Promoter** | H3K4me3+ AND NOT K27me3 AND <=2kb TSS |
| 2 | **Repressed_Promoter** | H3K27me3+ AND NOT K27ac AND <=2kb TSS |
| 3 | **Bivalent_Promoter** | K4me3+K27me3 overlap (developmental poised) |
| 4 | **Polycomb** | H3K27me3+ AND >2kb TSS (distal repressive) |
| 5 | **Active_Enhancer** | H3K27ac+ AND >2kb TSS |
| 6 | **Poised_Enhancer** | H3K4me1+ AND NOT K27ac AND NOT K27me3 AND >2kb |
| 7 | **CTCF_Site** | CTCF DNA motif+ (114,081 genome-wide motifs) |
| 8 | **Other** | No marks / structural elements |

### Anchor Type Distribution (Late, 2,910 Loops)

| Category | Anchor 1 | Anchor 2 |
|----------|---------|---------|
| Active_Promoter | 8.7% | 8.8% |
| Repressed_Promoter | 9.2% | 10.0% |
| Bivalent_Promoter | 0.9% | 0.8% |
| Polycomb | 16.0% | 15.8% |
| Active_Enhancer | 10.8% | 12.0% |
| Poised_Enhancer | 21.2% | 19.1% |
| Other | 33.2% | 33.4% |

### Top Loop Type Combinations (Late)

| Loop Type | Count | % |
|-----------|-------|---|
| Other-Other | 491 | 16.9% |
| Poised_Enhancer-Other | 386 | 13.3% |
| Polycomb-Other | 237 | 8.1% |
| Polycomb-Polycomb | 175 | 6.0% |
| Repressed_Promoter-Polycomb | 172 | 5.9% |
| Active_Enhancer-Poised_Enhancer | 162 | 5.6% |

### ChIP-seq Peak Files Used

| Mark | Early Peaks | Late Peaks |
|------|------------|------------|
| H3K27ac | 18,178 | 15,105 |
| H3K27me3 | 12,473 | 15,809 |
| H3K4me1 | 93,859 | 113,781 |
| H3K4me3 | 10,370 | 6,581 |
| Bivalent (K4me3+K27me3) | 931 | 318 |
| CTCF (DNA motifs) | 114,081 | 114,081 |
| CTCF (ChIP-seq) | 32,487 | 32,487 |

---

## 5. Loop Distance Shift & "Loop Rewriting" Phenomenon

### Core Discovery

**BAP1-KO preferentially loses long-range chromatin loops and gains shorter-range loops** -- a phenomenon termed "loop rewriting."

### Distance Statistics (Late Timepoint)

| Metric | Lost Loops | Gained Loops | Ratio |
|--------|-----------|-------------|-------|
| Count | 1,187 | 1,723 | 0.69 |
| **Median distance** | **625 kb** | **320 kb** | **1.95x** |
| <100 kb | 8.1% | 12.0% | 0.68 |
| 100-500 kb | 35.3% | 58.7% | 0.60 |
| 500 kb - 1 Mb | 21.8% | 18.8% | 1.16 |
| **>1 Mb** | **34.8%** | **10.5%** | **3.31x** |

### Statistical Significance

| Test | Statistic | P-value |
|------|-----------|---------|
| Kolmogorov-Smirnov | D=0.279 | <2.2e-16 |
| Wilcoxon Rank-Sum | W=1,356,071 | <2.2e-16 |
| Chi-Square (distance categories) | chi2=294.94, df=3 | <2.2e-16 |
| Spearman (distance vs logFC) | rho=-0.244 | <2.2e-16 |

### Biological Interpretation

Long-range loops (>1 Mb) are **3.3x enriched for loss** in BAP1-KO. This is consistent with a model where:
1. BAP1 loss causes H2AK119ub accumulation at regulatory elements
2. Excess ubiquitination destabilizes long-range Polycomb-mediated contacts
3. Chromatin architecture reorganizes into shorter-range, potentially ectopic interactions
4. Creates a "regulatory traffic jam" where developmental loops collapse

---

## 6. Shared Anchor Analysis (Loop Switching)

### Hypothesis

Specific genomic loci act as "anchor hubs" that participate in BOTH lost and gained loops -- losing a long-range contact and gaining a short-range one at the same site.

### Results (Late Timepoint)

| Metric | Value |
|--------|-------|
| Shared anchor regions identified | 212 |
| Lost loops at shared anchors | 283 (23.8% of all lost) |
| Gained loops at shared anchors | 321 (18.6% of all gained) |
| Total loops at shared anchors | 604 (20.8% of differential) |
| **Direction support (lost > gained distance)** | **83.0%** |
| Median lost loop distance | 1,150 kb |
| Median gained loop distance | 340 kb |
| **Median distance difference** | **694 kb (3.4x)** |
| Paired Wilcoxon signed-rank test | **p=1.17e-20** |

### ChIP-seq Enrichment at Shared Anchors

| Mark | Odds Ratio | P-value | Direction |
|------|-----------|---------|-----------|
| **H3K27me3** | **2.04** | **1.75e-24** | ENRICHED |
| **Bivalent** | **1.74** | **0.013** | ENRICHED |
| H3K4me1 | 1.27 | 4.6e-04 | Enriched |
| H3K27ac | 0.68 | 3.9e-06 | DEPLETED |
| H3K4me3 | 1.05 | 0.64 | No difference |

### Polycomb-Specific Shared Anchor Analysis

- 312 Polycomb-classified shared loops analyzed separately
- Polycomb anchors show stronger distance shift than non-Polycomb anchors
- Confirms Polycomb-mediated mechanism for loop rewriting

### Assessment

Loop switching is **partially supported**: explains ~20% of differential loops with strong statistical significance and expected H3K27me3 enrichment. However, 22% of shared anchors show opposite pattern, and no significant gene expression change detected at shared anchor genes (p=0.63).

---

## 7. H2AK119ub Integration (BAP1 Substrate)

### Rationale

BAP1 is the deubiquitinase for H2AK119ub. In BAP1-KO, H2AK119ub accumulates at Polycomb target sites. This analysis tests whether H2AK119ub changes predict loop changes.

### K119ub Peak Data

| Category | Peaks |
|----------|-------|
| K119ub_up (gained in mutant) | 6,164 |
| K119ub_down (lost in mutant) | 1,250 |
| K119ub_ctrl (control) | 20,592 |
| K119ub_mut (mutant) | 20,399 |

### Key Results

#### Enrichment Analysis (11/12 tests significant)
- K119ub_up at long_range_lost vs short_range_gained: **OR=4.87, FDR=3.73e-90**
- Interpretation: Loops with elevated K119ub are enriched for being lost at long range

#### Continuous Signal Correlation
| Anchor Type | Spearman rho | P-value | N |
|-------------|-------------|---------|---|
| Active anchors | **-0.314** | 7.46e-19 | 764 |
| Other anchors | -0.013 | 0.676 | 973 |
| Polycomb/Repressive | +0.177 | 2.82e-09 | 1,118 |

#### Logistic Regression: P(lost) ~ K119ub_fc + log(distance)
- K119ub_fc coefficient: 2.370 (p=1.71e-91)
- **Odds Ratio: 10.70** -- higher K119ub change strongly predicts loop loss

#### Shared Anchor K119ub Enrichment
| K119ub Category | Shared vs Non-Shared OR | FDR |
|-----------------|------------------------|-----|
| K119ub_down | **OR=2.54** | **1.95e-04** |
| K119ub_ctrl | **OR=2.41** | **2.49e-17** |
| K119ub_mut | **OR=1.99** | **2.20e-11** |

---

## 8. Differential ChIP-seq Polycomb Enrichment

### Analysis Design

Tested whether differential H3K27me3 and H2AK119ub peaks specifically overlap loop categories (long-range lost, short-range gained, unchanged).

### Significant Enrichments (All Loops, FDR < 0.05)

| Test | OR | FDR |
|------|-----|-----|
| K27me3_down at short_range_gained vs unchanged | **8.78** | **6.86e-102** |
| H2AK119ub_up at long_range_lost vs short_range_gained | **4.87** | **1.86e-90** |
| H2AK119ub_down at short_range_gained vs unchanged | **3.20** | **4.46e-41** |
| K27me3_up at long_range_lost vs short_range_gained | **5.05** | **2.85e-28** |
| K27me3_up at long_range_lost vs unchanged | **3.18** | **3.58e-28** |

### Interpretation

- Long-range lost loops: Enriched for K27me3_UP and K119ub_UP (Polycomb accumulation)
- Short-range gained loops: Enriched for K27me3_DOWN and K119ub_DOWN (Polycomb depletion)
- Supports mechanistic model: BAP1 loss redistributes Polycomb marks, destabilizing long-range contacts while enabling short-range ones

---

## 9. DEG-Loop Anchor Integration

### Design

Connected differential gene expression (RNA-seq) with loop anchor changes using GREAT-style regulatory domain assignment (5kb upstream, 1kb downstream, 100kb max extension).

### Results (Late Timepoint)

| Metric | Value |
|--------|-------|
| Total DEGs near differential anchors | 457 |
| DEGs near lost anchors | 230 |
| DEGs near gained anchors | ~200+ |

### Chromatin State-Stratified Analysis

Violin plots showing DEG log2FC distribution stratified by:
- Loop direction (lost vs gained)
- Distance category
- Anchor chromatin state (Active_Promoter, Polycomb, etc.)

### Permutation Testing

1,000 random anchor assignments as baseline to establish significance of observed DEG associations.

---

## 10. TAD Boundary Analysis (TADCompare)

### Pipeline

6-step TADCompare (Cresswell & Dozmorov 2020) pipeline for differential TAD boundary detection at 10kb resolution:

1. Extract sparse contact matrices from .hic files (straw)
2. Run TADCompare on merged contact matrices
3. ConsensusTADs robustness assessment across individual replicates
4. Shift distance calculation (bedtools closest)
5. Blacklist filtering (ENCODE mm10 + lab-specific)
6. Comprehensive visualization suite (40+ plots)

### Results

| Metric | Late | Early |
|--------|------|-------|
| Total boundaries detected | 23,496 | 24,801 |
| **Differential boundaries** | **3,822 (16.3%)** | **5,055 (20.4%)** |
| Shifted | 1,000 (4.3%) | 1,175 (4.7%) |
| Strength Change | 1,416 (6.0%) | 1,114 (4.5%) |
| Complex | 922 (3.9%) | 1,121 (4.5%) |
| Merge | 838 (3.6%) | 780 (3.1%) |
| Split | 640 (2.7%) | 865 (3.5%) |
| Median shift distance | 0 kb | N/A |
| Mean shift distance | 6.2 kb | 4.6 kb |
| After blacklist filtering | 21,755 | 23,104 |

### ChIP-seq Classification of Boundaries

| State | Early | Late |
|-------|-------|------|
| Active_Promoter | 15.4% | 11.2% |
| CTCF_Site | 24.2% | 12.7% |
| Poised_Enhancer | 22.7% | 19.7% |
| Polycomb | 5.2% | 6.2% |

### HOMER TAD Inclusion Ratio Analysis

| Threshold | TADs Meeting Criteria | Direction |
|-----------|----------------------|-----------|
| Relaxed (FDR<0.15, |Diff|>0.15) | 923 (18.2%) | 809 increased, 114 decreased |
| Standard (FDR<0.05, |Diff|>0.30) | 101 (2.0%) | 84 increased, 17 decreased |

### Boundary-Loop Cross-Reference

- Lost loops: 2.15x odds of being within 10kb of a boundary (p=0.039)
- But only 46.7% directional concordance -- boundary and loop changes partially independent

### Key Finding

**TAD structure is relatively stable** with BAP1 loss (only 16-20% differential boundaries), suggesting BAP1's primary effects are at the loop and compartment level, not the TAD boundary level.

---

## 11. Chromatin Compartment (PC1) Analysis

### HOMER Differential Compartment Analysis

| Threshold | Regions | Interpretation |
|-----------|---------|---------------|
| Relaxed (FDR<0.15, |Diff|>0.15) | 60,388 (59.4%) | Massive shifts |
| Standard (FDR<0.05, |Diff|>0.30) | 44,703 (44.0%) | 6,462 B->A, 4,184 A->B |

### Key Findings

- **44% of the genome shows significant compartment shifts** -- far more dramatic than TAD changes
- 1,830 unique genes in significantly shifted regions
- Strongest shifts on chrX (Ctag2, Fmr1os, Mir509 loci: 2.7-2.8 difference units)
- Net direction: More B->A (active) shifts than A->B (repressive)

### Biological Interpretation

Unlike TADs (relatively stable), **chromatin compartments undergo dramatic reorganization** with BAP1 loss. Widespread shifts toward more active chromatin are consistent with loss of Polycomb-mediated repression (BAP1-KO fails to remove H2AK119ub at proper developmental timepoints).

---

## 12. Chromatin Stripe Analysis (Quagga)

### Pipeline

4-phase differential stripe analysis using Quagga (Feng et al., Genome Research 2025):

1. **Detection & Union Set** -- Load merged 5kb stripes (min_length=200kb), classify as control_only/mutant_only/shared
2. **Hi-C Quantification** -- Convert to GInteractions, extract 5x5 matrices via mariner
3. **edgeR Differential** -- Quasi-likelihood GLM (same framework as loops)
4. **Integration & Classification** -- Tiered confidence (high/medium/low)

### Results

| Metric | Early | Late |
|--------|-------|------|
| Total stripes | 286 | 200 |
| Lost (control_only) | 126 (44%) | 83 (42%) |
| Gained (mutant_only) | 86 (30%) | 73 (36%) |
| Shared | 74 (26%) | 44 (22%) |
| FDR < 0.05 | 0 | 0 |
| FDR < 0.10 | 0 | 1 |
| BCV | 0.067 (6.7%) | 0.062 (6.2%) |

### The Single High-Confidence Stripe

**stripe_0014 (GAINED, Late):**
- Location: chr10:90,660,000-91,150,000
- Gene: *Apaf1* (Apoptotic Protease Activating Factor 1)
- Anchor: Active_Enhancer (H3K27ac+, H3K4me1+)
- FDR: 0.0753 (only stripe below 0.10)
- logFC: 0.34

### Key Finding

**BAP1-KO does NOT strongly alter stripe architecture.** The null finding is informative -- stripe-level organization is largely preserved, suggesting BAP1's effects manifest at finer scales (individual loops) or broader scales (compartments), not at the stripe/cohesin-extrusion level. The remarkably low BCV (~6.7%) indicates excellent technical reproducibility and sufficient power to detect real changes.

---

## 13. CTCF-Stripe Cross-Reference

### Analysis

Tested whether lost CTCF-anchored loops co-occur with stripe defects, which would indicate common cohesin/extrusion dysfunction.

### Results

- Loop anchors overlapped with stripe anchors (25kb tolerance)
- Fisher's exact test for enrichment of CTCF loops at stripe sites
- Venn diagrams showing overlap between CTCF loops and stripes
- Minimal coordination between loop and stripe changes

---

## 14. ABC Enhancer-Gene Linkage Model

### Overview

Implemented the Broad Institute's **Activity-By-Contact (ABC) model** to predict enhancer-gene regulatory linkages and quantify how BAP1 loss alters these connections.

**ABC Formula:** `ABC = (Enhancer Activity x Contact Frequency) / Sum(Activity x Contact) within 5 Mb`

### Critical Design Innovation

Both WT and KO conditions use the **identical consensus enhancer universe** (75,371 ATAC peaks) to ensure DABC measures true linkage changes, not enhancer definition artifacts.

### Scale of Analysis

| Metric | WT | KO |
|--------|----|----|
| Total E-G pairs (>=0.02) | 122,212 | 119,965 |
| Distal E-G pairs | 93,060 | 90,813 |
| Unique genes | 27,492 | 27,492 |
| Median distance (distal) | 127 kb | 116 kb |

### Delta-ABC Distribution

| Category | Count | % |
|----------|-------|---|
| Gained (D > 0.01) | 44,060 | 24.4% |
| Lost (D < -0.01) | 47,499 | 26.3% |
| Unchanged | 88,864 | 49.3% |

### Headline Results

#### RNA-seq Integration (Strongest Correlation)
- **Unnormalized D(AxC) sum across enhancers: Spearman rho = 0.582** (p ~ 10^-172)
- ~34% of rank variance explained
- **Methodological contribution:** Unnormalized outperforms normalized DABC (65.3% vs 58.8% concordance)

#### Paired-Anchor Loop-ABC Concordance (Geometric Constraint)
- Independent loop-ABC overlap: 51.4% (at chance)
- **With geometric constraint (enhancer at one anchor, TSS at other): 89.7%** (p=1.67e-48)
- **3-way concordance (loop + ABC + RNA-seq): 88.2%** (p=1.69e-45)

#### K119ub at Paired Enhancers
- **Spearman rho = -0.401** (p=5.48e-13) -- 4.5x stronger than genome-wide
- Lost enhancers: median K119ub log2FC = +0.321
- Gained enhancers: median K119ub log2FC = -0.074

### Enhancer Subset Stratified Analysis (Step 11)

| Class | N | Definition | Loop logFC | RNA Concordance |
|-------|---|-----------|-----------|----------------|
| Activity_Lost | 7,503 | K27ac down | -0.088 | **59.8%** |
| **K119ub_Only** | **2,479** | K119ub up, marks unchanged | **-0.054** | **49.4% (chance)** |
| Activity_Gain | 2,851 | K27ac up | +0.066 | **61.9%** |
| Stable | 42,864 | No changes | -0.013 | 48.8% |

### Critical Finding

**K119ub_Only enhancers show real but sub-functional loop weakening** (3.7% reduction, p<2.2e-16) that does NOT translate to detectable gene expression changes (concordance at chance). This establishes that K119ub accumulation alone is insufficient for transcriptional consequence in steady-state adult tissue -- suggesting the effect requires additional cooperative disruption or manifests during developmental transitions.

---

## 15. Machine Learning: LoopBin Unsupervised Clustering

### Method

**LoopBin** -- Variational Autoencoder with Deep Embedding (VADE) model for unsupervised loop classification:
- Input: ~384 features per loop (256 from 16x16 Hi-C matrices + ~128 from Cut&Tag epigenetic signals)
- Architecture: 500-500-2000-10 (latent) encoder/decoder
- Gaussian Mixture Model for cluster assignment
- 7-cluster final model after merging small clusters

### Results

| Cluster | Loops | % Total |
|---------|-------|---------|
| 0 | 617 | 25.5% |
| 6 | 452 | 18.7% |
| 9 | 505 | 20.9% |
| 3 | 252 | 10.4% |
| 8 | 247 | 10.2% |
| 7 | 139 | 5.8% |
| 4 | 109 | 4.5% |

### Key Finding

**Loops don't disappear when BAP1 is degraded -- they shift between clusters.** BAP1 loss alters loop contact pattern properties rather than eliminating loops entirely. This provides a complementary view to the statistical (edgeR) approach:
- edgeR asks: "which loops change in strength?"
- LoopBin asks: "how do loops change their character?"

### Cluster Transitions

- Control Cluster 1 -> Degron Cluster 2: 21 loops
- Control Cluster 1 -> Degron Cluster 6: 31 loops
- Cluster 5 shows strongest DEG correlation (log2FC 2.34 in condition-specific loops)

---

## 16. Biomodal Methylation Analysis (5mC/5hmC)

### Technology

**Biomodal DUET evoC** platform -- dual-epigenetic sequencing providing simultaneous 5mC and 5hmC quantification at 6bp resolution.

### Pipeline

1. **Upstream:** DUET Nextflow pipeline v1.5.0 (FASTQ -> Zarr stores)
2. **Downstream:** DMR calling via modality XPLR CLI (GLM with B-H FDR)

### Results

- CpG methylation: ~82% methylated (genome-wide)
- CHG/CHH contexts: <1% (expected for mammalian genome)
- **Critical limitation:** Sex confounding with n=2 per condition -- with sex covariate, no significant DMRs; without, confounded

### 17-Section Visualization Pipeline

Developed a comprehensive **17-section R visualization pipeline** (`downstream/scripts/viz_sections/section_*.R`) covering:
- DMR statistics, correlations, effect sizes
- Coordinated changes, enrichment, chromatin state
- Cross-data integration: K119ub, ATAC, MeCP2, loop-methylation correlation
- Each section is self-contained, auto-included via `run_all_sections.sh`

---

## 17. ATAC-seq Consensus Peak Analysis

### Consensus Generation

10 raw narrowPeak files processed:
- Batch 1 (250310): n=3 control + n=3 mutant (primary)
- Batch 2 (250731): n=2 control + n=2 mutant (uncertain identity, possibly technical replicates)

### Method

`bedtools multiIntersectBed` with >=2 replicate threshold, merged overlapping regions.

### Results

| Category | Peaks |
|----------|-------|
| Control consensus (Batch 1) | ~66,000 |
| Mutant consensus (Batch 1) | ~68,000 |
| Union (all accessible) | ~71,000 |
| **ATAC_up (more accessible in mutant)** | **7,620** |
| **ATAC_down (more accessible in control)** | **3,744** |

### Use in ABC Model

The union consensus (75,371 peaks) served as the candidate enhancer universe for the ABC model, ensuring identical enhancer definitions across conditions.

---

## 18. HOMER Motif Analysis

### CTCF Motif Extraction

- Extracted canonical CTCF(Zf) motifs from HOMER's pre-computed mm10 database
- **114,081 CTCF motifs** genome-wide
- Used for loop anchor classification (CTCF_Site category) across both timepoints
- Design decision: DNA motifs preferred over ChIP-seq for cross-timepoint consistency

### Motif Enrichment at Enhancer Subsets

HOMER motif enrichment run on 7 enhancer subset BED files (4 phenotypic classes + 3 K119ub tertiles) -- visualization pending completion.

---

## 19. Aggregate Peak Analysis (APA)

### Method

APA averages Hi-C contact matrices across sets of loops to visualize aggregate enrichment:
- Extract matrices around loops (50kb windows)
- Average across loops per condition
- Calculate P2LL (Peak to Lower-Left) enrichment score
- Statistical comparison via Wilcoxon test

### Scale

| Metric | Value |
|--------|-------|
| Resolutions processed | 3 (5kb, 10kb, 25kb) |
| Loop sets analyzed | 6 (resolution-specific + merged, per direction) |
| Total successful analyses | 12/12 |
| Total APA PDF outputs | 155 |
| Runtime | ~2 hours per resolution |

### Results

APA heatmaps confirm differential contact frequency changes:
- **Upregulated loops:** Higher enrichment in mutant
- **Downregulated loops:** Higher enrichment in control
- Visual confirmation of the statistical findings from edgeR

---

## 20. Hi-C Processing Infrastructure (HPC)

### Complete HPC Pipeline Documentation

Documented the entire Hi-C processing workflow from raw sequencing data to analysis-ready files, including:

1. **Data retrieval** from IGM FTP server and import to SDSC Expanse
2. **Juicer pipeline** for shallow sequencing QC
3. **Nextflow-HiC pipeline** for deep sequencing production
4. **fastp trimming** for adapter removal
5. **juicerpre** for Juicer v2.0 .hic file generation
6. **HiCCUPS** for loop calling
7. **Cooler utilities** (balance, dump, zoomify, info)
8. **HOMER tag directory** generation for TAD/PC1 analysis
9. **GENOVA** configuration for loci visualization and aggregate analysis
10. **multiHiCompare** setup for replicate-aware Hi-C comparison
11. **dchic** setup for differential compartment analysis
12. **hictk** for file format conversion and merging

### Computing Resources

- **SDSC Expanse:** Primary HPC (SLURM, project csd940)
- **AWS EC2:** GENOVA and multiHiCompare (m5.4xlarge, upgraded to m5.8xlarge for high-memory jobs)
- **Local:** Downstream analysis, visualization, documentation

---

## 21. Technical Innovations & Methodological Contributions

### 1. Multi-Resolution Differential Loop Analysis

First implementation of a complete multi-resolution (5kb/10kb/25kb) replicate-aware differential loop pipeline using mariner + edgeR, with cross-resolution validation and non-redundant merging.

### 2. 8-Category Chromatin State Classification

Priority-based classification using 5 histone marks + CTCF motifs, with CTCF motif-based classification for cross-timepoint consistency (rather than condition-specific ChIP-seq).

### 3. Unnormalized Delta(AxC) for Differential ABC

Discovery that unnormalized D(AxC) sum outperforms normalized DABC for differential analysis (rho 0.582 vs 0.035), because per-gene normalization compresses real activity changes in widespread remodeling.

### 4. Geometric Constraint for Loop-ABC Concordance

Requiring enhancer at one anchor and TSS at the other increases loop-ABC concordance from 51.4% (chance) to 89.7% (highly significant) -- demonstrating that geometric co-localization is critical for meaningful cross-method integration.

### 5. K119ub Dose-Response at Enhancer Subsets

Stratified analysis showing K119ub_Only enhancers produce real but sub-functional contact changes (3.7% reduction, OR=10.70 in logistic regression) that don't reach transcriptional consequence -- establishing a threshold model for Polycomb-mediated chromatin disruption.

### 6. Tiered Stripe Confidence Classification

Detection-primary confidence tiers (high/medium/low) for stripe analysis where traditional FDR-based classification fails due to small feature sets.

### 7. Loop Switching Quantification

Formal quantification of the "shared anchor" phenomenon with paired distance statistics, ChIP-seq enrichment, and Polycomb-specific subsetting.

---

## 22. Complete Results Summary Table

| Analysis | Key Metric | Value |
|----------|-----------|-------|
| **Differential Loops (Late)** | Non-redundant differential | 2,910 |
| | Up in BAP1-KO | 1,723 (59.2%) |
| | Down in BAP1-KO | 1,187 (40.8%) |
| **Differential Loops (Early)** | Non-redundant differential | 165 |
| **Loop Rewriting** | Median distance lost/gained | 625 kb / 320 kb (1.95x) |
| | >1Mb lost enrichment | 3.31x |
| **Shared Anchors** | Hub regions identified | 212 |
| | Direction support | 83.0% |
| | H3K27me3 enrichment OR | 2.04 (p=1.75e-24) |
| **K119ub Integration** | Logistic regression OR | 10.70 (p=1.71e-91) |
| | Active anchor correlation | rho=-0.314 |
| **TAD Boundaries** | Differential (Late) | 3,822 (16.3%) |
| | Differential (Early) | 5,055 (20.4%) |
| **Compartments** | Significant shifts (standard) | 44,703 (44.0%) |
| | Unique genes affected | 1,830 |
| **Stripes** | Significant (FDR<0.10) | 1 (out of ~200-286) |
| **ABC Model** | E-G pairs analyzed | 180,423 |
| | RNA-seq concordance (unnorm sum) | rho=0.582 |
| | Paired-anchor 3-way concordance | 88.2% |
| | K119ub at paired enhancers | rho=-0.401 |
| **LoopBin ML** | Loop clusters | 7 |
| | Cluster transitions | 66 condition-specific |
| **ATAC-seq** | Consensus peaks | ~71,000 |
| | Differential (up/down) | 7,620 / 3,744 |
| **Total Visualization Outputs** | PDFs generated | 500+ |
| **Total Data Files** | TSV/RDS/BEDPE | 200+ |

---

## 23. Key Biological Findings & Unified Model

### The Multi-Scale Chromatin Disruption Model

BAP1 loss in the cerebellum causes a **progressive, multi-scale disruption** of 3D chromatin architecture:

#### Scale 1: Compartments (Most Affected)
- 44% of the genome shows significant compartment shifts
- Net B->A (more active) direction consistent with loss of Polycomb repression
- Widespread but graded effect

#### Scale 2: Loops (Primary Effect)
- 2,910 differential loops at late timepoint (17.6x increase from early)
- "Loop rewriting": long-range loops lost, short-range gained
- Lost loops enriched for H3K27me3 and H2AK119ub accumulation
- 212 shared anchor hubs where loops switch partners

#### Scale 3: TAD Boundaries (Moderate Effect)
- 16-20% of boundaries differential
- Boundaries relatively stable but not immune
- Lost loops near boundaries (2.15x enrichment)

#### Scale 4: Stripes (Minimal Effect)
- No significant widespread changes
- Stripe/cohesin extrusion architecture preserved

### Temporal Dynamics

| Feature | Early (P12) | Late (P60) | Interpretation |
|---------|-------------|------------|----------------|
| Loops | 165 differential | 2,910 differential | 17.6x amplification |
| Direction | 57% DOWN | 59% UP | Phase reversal |
| Loop types | Repressed_Promoter | Mixed/diverse | Broader involvement |
| TAD boundaries | 20.4% differential | 16.3% differential | Slight stabilization |
| Compartments | N/A (same analysis) | 44% shifted | Major reorganization |

### Mechanistic Chain

1. **BAP1 loss** -> H2AK119ub accumulates (cannot be removed)
2. **K119ub accumulation** -> Destabilizes long-range Polycomb contacts (OR=10.70)
3. **Long-range loop loss** -> Reorganization into shorter-range contacts ("loop rewriting")
4. **At enhancers**: K119ub_Only causes 3.7% contact weakening (sub-functional)
5. **Activity loss** (K27ac down): Causes functional E-G disruption (59.8% concordance)
6. **Gene expression**: 940 dysregulated genes with concordant DABC and DE
7. **Temporal progression**: Early subtle -> Late massive (regulatory cascade)

### Key Loci of Interest

- **Syt1/Nav3** (chr10): Most impacted region, essentially collapsed by P60
- **Apaf1** (chr10): Only significant stripe change (apoptosis regulator)
- **Hox clusters**: Developmental gene regulatory hubs affected
- **Sox2/6/9, Foxg1, Shh**: Gained loop targets (developmental transcription factors)
- **BDNF, Mef2c, Pax6**: Lost loop targets (neuronal function genes)

---

## 24. Publication-Ready Outputs

### Visualization Count

| Category | Count |
|----------|-------|
| APA heatmaps | 155 PDFs |
| Volcano plots | 19 PDFs |
| Loop classification plots | 4 PDFs |
| GO/KEGG enrichment | 4 PDFs |
| H2AK119ub integration plots | 20 PDFs |
| Polycomb enrichment plots | 9 PDFs |
| Shared anchor plots | 8 PDFs |
| DEG-loop violin plots | 8+ PDFs |
| TAD analysis plots | 40+ PDFs |
| Stripe analysis plots | 20+ PDFs |
| ABC analysis plots | 30+ PDFs |
| LoopBin ML visualizations | 15+ PDFs |
| Biomodal methylation plots | 17 sections |
| Distance shift plots | 6+ PDFs |
| **Total** | **500+ publication-ready PDFs** |

### Key Data Tables

| File | Rows | Content |
|------|------|---------|
| characterized_loops.tsv | 2,910 (or 39,344 all) | Annotated differential loops |
| delta_abc_all_pairs.tsv | 180,423 | All E-G pairs with DABC |
| gene_level_summary.tsv | 13,588 | Per-gene regulatory summary |
| tadcompare_final_filtered.tsv | ~21,000-23,000 | Annotated TAD boundaries |
| unified_stripes.rds | 200-286 | Classified differential stripes |
| merged_log_data.npy | 2,910 x 384 | LoopBin feature matrix |

---

## 25. Tools & Technologies Used

### Languages
- **R** (primary): Bioconductor packages, ggplot2, statistical modeling
- **Python**: pandas, numpy, ABC model, feature extraction
- **Bash**: Pipeline orchestration, HPC job submission, file processing
- **Java**: Juicer tools backend

### Key R/Bioconductor Packages
- **mariner**: Hi-C loop analysis framework (pullHicMatrices, mergePairs, binPairs)
- **edgeR**: Differential analysis (quasi-likelihood GLM)
- **InteractionSet**: Genomic interaction data structures (GInteractions)
- **GenomicRanges**: Genomic interval operations
- **HDF5Array**: Out-of-memory array storage
- **strawr**: .hic file reading via Straw API
- **TADCompare**: Differential TAD boundary detection
- **ChIPseeker**: Peak annotation and genomic feature distribution
- **clusterProfiler**: GO/KEGG functional enrichment
- **EnhancedVolcano**: Publication volcano plots
- **GENOVA**: Hi-C visualization and aggregate analysis

### Bioinformatics Tools
- **HiCCUPS**: Loop calling (Juicer suite)
- **Juicer**: Hi-C processing, APA, Arrowhead
- **Nextflow-HiC**: Production Hi-C pipeline (bowtie2 + cooler)
- **HOMER**: TAD analysis, compartments, motif enrichment
- **Quagga**: Stripe detection (Gabor filtering + Poisson statistics)
- **ABC Model**: Broad Institute enhancer-gene prediction
- **deepTools**: ChIP-seq heatmap visualization
- **bedtools**: Genomic interval operations
- **cooler**: Hi-C contact map manipulation
- **hictk**: Format conversion and merging
- **fastp**: Read trimming
- **Biomodal DUET evoC**: Dual 5mC/5hmC methylation
- **modality XPLR**: DMR calling

### ML/Deep Learning
- **PyTorch**: VADE model (Variational Autoencoder + Deep Embedding)
- **scikit-learn**: GMM, silhouette scoring, t-SNE

### Infrastructure
- **SDSC Expanse**: Primary HPC (SLURM)
- **AWS EC2**: GENOVA and multiHiCompare (m5.4x/8xlarge)
- **Git/GitHub**: Version control and documentation

---

*Document compiled: February 2026*
*Repository: mariner_hi-c*
*Total codebase: ~45 analysis scripts, 500+ output files, 200+ data tables*
