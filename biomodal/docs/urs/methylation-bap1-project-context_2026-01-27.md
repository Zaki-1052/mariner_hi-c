# BAP1-KO Mouse Cerebellum Project: Complete Experimental Context

## Purpose of This Document

This document provides comprehensive context for understanding an ongoing research project studying BAP1 knockout in mouse cerebellum. It synthesizes experimental design, key findings from methylation and chromatin loop analyses, and the emerging biological narrative. Use this alongside the separate literature review document to understand how these findings fit into existing knowledge.

---

## 1. Project Overview

### What is BAP1?

BAP1 (BRCA1-Associated Protein 1) is a deubiquitinase that serves as the catalytic subunit of the **PR-DUB (Polycomb Repressive-Deubiquitinase) complex**. Its primary function is removing monoubiquitin from histone H2A lysine 119 (H2AK119ub). BAP1 is well-characterized as a tumor suppressor in cancer contexts but has been minimally studied in the brain.

### The Experimental Model

**Organism**: Mouse (*Mus musculus*)  
**Genome**: GRCm38/mm10  
**Tissue**: Cerebellum  
**Genotypes**:
- Control (wildtype): Bap1-ff (floxed but no Cre)
- Mutant: Bap1-ff-cre (conditional BAP1 knockout)

### Lab Context

The Ferguson Lab studies how genes control brain formation and function, with particular interest in ubiquitin signaling and degradative pathways in neurodevelopmental and neurodegenerative disorders. Methods include genomics (CUT&RUN, RNA-seq, ATAC-seq, Hi-C, Hi-ChIP) and proteomics (TMT-MS) in genetically modified mice.

A potential collaboration with the Rao Lab (La Jolla Institute) is relevant—they discovered TET enzyme activities and study TET-mediated DNA demethylation, which intersects with findings from this project.

---

## 2. Two Parallel Analyses

This project involves two major genomic analyses on the same BAP1-KO model:

| Analysis | Technology | Sample Size | Primary Finding |
|----------|------------|-------------|-----------------|
| **DNA Methylation** | Biomodal DUET evoC (6bp resolution) | n=2 per condition (4 total) | Coordinated 5mC↑ / 5hmC↓ at gene bodies |
| **Chromatin Loops** | Hi-C with HiCCUPS loop calling | n=3 per condition (6 total) | Differential looping, Polycomb anchor enrichment |

---

## 3. Methylation Analysis: DUET evoC Results

### Technology

Biomodal DUET evoC provides simultaneous, base-resolution quantification of both 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC) at 6bp resolution—enabling detection of both marks at the same genomic loci.

### Sample Information

| Sample ID | Condition | Sex | Genotype | Mice Pooled |
|-----------|-----------|-----|----------|-------------|
| evoC-Bap1-ctrl-F | Control | Female | Bap1-ff | 3398-3379-3138 |
| evoC-Bap1-ctrl-M | Control | Male | Bap1-ff | 3139-3147-3387 |
| evoC-Bap1-mut-F | Mutant | Female | Bap1-ff-cre | 3380-3388-3381 |
| evoC-Bap1-mut-M | Mutant | Male | Bap1-ff-cre | 3391-3392-3386 |

### Critical Caveat: Sex Confounding

With n=2 per condition:
- **Run-1 (with sex covariate)**: No significant results
- **Run-2 (without sex covariate)**: Significant results

This indicates sex effects are confounded with BAP1 effects. More replicates are being generated.

### Key Statistics

**5mC (5-methylcytosine) Changes:**
- Total genes tested: 20,979
- Significant (q < 0.05): **8,836 (42.1%)**
- Hypermethylated: 6,635 (75.1%)
- Hypomethylated: 2,201 (24.9%)
- Mean change (significant): **+2.27%**

**5hmC (5-hydroxymethylcytosine) Changes:**
- Total genes tested: 20,979
- Significant (q < 0.05): **9,930 (47.3%)**
- Increased: 1,307 (13.2%)
- Decreased: 8,623 (86.8%)
- Mean change (significant): **-2.08%**

### The Coordinated Pattern

**Critical finding**: 5mC and 5hmC change in opposite directions at the same loci.

- Genes significant in both mC AND hmC: 6,750
- Genes with **mC↑ / hmC↓ pattern**: **5,708 (84.6%)**
- Genes with opposite pattern (mC↓ / hmC↑): 1,042 (15.4%)

This 84.6% concordance is remarkably high and suggests a specific mechanistic cause rather than random dysregulation.

### Most Affected Genes

| Gene | mC Change | hmC Change | Combined Effect | Function |
|------|-----------|------------|-----------------|----------|
| **Syt1** | +17.9% | -15.3% | 33.2% | Synaptic vesicle exocytosis, Ca²⁺ sensor |
| Tmem238 | +20.3% | -11.1% | 31.4% | Transmembrane protein |
| Sap30 | +22.4% | -6.3% | 28.7% | Sin3A-associated protein, chromatin |
| Gclm | +16.3% | -11.4% | 27.7% | Glutamate-cysteine ligase, antioxidant |
| Mcu | +14.7% | -9.8% | 24.5% | Mitochondrial calcium uniporter |

**Syt1 (Synaptotagmin-1)** is the most significantly affected gene and encodes the primary calcium sensor for fast synchronous neurotransmitter release.

### Genomic Distribution

**Primary affected regions**: Gene bodies (NOT promoters or CpG islands)

| Region | mC Significant | hmC Significant |
|--------|----------------|-----------------|
| **Gene bodies** | **8,836 (42.1%)** | **9,930 (47.3%)** |
| CpG shores | 772 (2.8%) | 55 (0.2%) |
| CpG shelves | 397 (1.5%) | 65 (0.3%) |
| CpG islands | 23 (1.9%) | 2 (0.2%) |
| Promoters | 28 (0.2%) | 4 (0.03%) |

### Chromatin State at DMRs

DMRs were classified using ChIP-seq peaks from Late cerebellum:

| Chromatin State | Count | Percentage |
|-----------------|-------|------------|
| **Active_Promoter** | 4,413 | **49.9%** |
| Repressed_Promoter | 1,211 | 13.7% |
| Other | 3,049 | 34.5% |
| Bivalent_Promoter | 101 | 1.1% |

**Key observation**: 94.0% of Active Promoter DMRs are hypermethylated. Hypomethylated DMRs are predominantly at Repressed Promoters (93.8%).

### Methylation Context

| Context | Baseline Methylation | Significant DMRs |
|---------|---------------------|------------------|
| **CG (CpG)** | ~82% modified | **8,836 mC / 9,930 hmC** |
| CHG | <1% modified | 0 |
| CHH | <1% modified | Not analyzed |

Non-CpG methylation shows no significant changes, confirming specificity to CpG methylation.

### Functional Enrichment (GO/KEGG)

Top enriched pathways among affected genes:
1. RNA splicing
2. Golgi vesicle transport
3. Vesicle organization
4. Macroautophagy
5. Autophagy
6. Endocytosis

Vesicular transport and RNA processing functions dominate.

---

## 4. Hi-C Chromatin Loop Analysis

### Technology and Pipeline

**Method**: Hi-C with HiCCUPS loop calling, analyzed using the mariner R/Bioconductor package with edgeR statistical framework.

**Key features**:
- Multi-resolution analysis: 5kb, 10kb, and 25kb
- Biological replicates: n=3 per condition (6 samples total)
- Quasi-likelihood GLM for robust differential testing

### Sample Structure

**Control (wildtype)**: ctrl_M1, ctrl_M2, ctrl_M3  
**Mutant (BAP1-KO)**: mut_M1, mut_M2, mut_M3

Each sample has HiCCUPS loop calls and Hi-C contact matrices (.hic files).

### Pipeline Workflow

1. **prep_loops.R** → Merge 6 replicates, create consensus loop positions
2. **extract_counts.R** → Extract 5×5 Hi-C matrices at each loop
3. **aggregate.R** → Sum matrices → count matrix (loops × samples)
4. **qc-val.R** → Quality control validation
5. **edgeR.R** → Quasi-likelihood GLM differential analysis
6. **compare_resolutions.R** → Cross-resolution comparison
7. **downstream_analysis.R** → Merge resolutions, annotate with genes & ChIP-seq
8. **visualizations.R** → Publication plots, enrichment analysis
9. **apa_analysis.R** → Aggregate Peak Analysis heatmaps

### Statistical Method

- Design: ~group (Control vs Mutant)
- Test: Quasi-likelihood GLM with robust dispersion
- Multiple testing: Benjamini-Hochberg FDR
- Primary threshold: FDR < 0.05
- Stringent threshold: FDR < 0.03 AND |logFC| > 0.3

### Power Improvement from Replicates

- **Without replicates (n=1)**: ~2 significant loops
- **With replicates (n=3)**: 500-5,000+ significant loops
- **Improvement**: 250-2,500× increase in statistical power

---

## 5. Loop Anchor Chromatin State Classification

### 8-Category System

Loop anchors are classified using ChIP-seq data (H3K27ac, H3K27me3, H3K4me1, H3K4me3) in priority order:

| Priority | Category | Definition |
|----------|----------|------------|
| 1 | **Active_Promoter** | H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS |
| 2 | **Repressed_Promoter** | H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS |
| 3 | **Bivalent_Promoter** | H3K4me3+H3K27me3 overlap (pre-computed) |
| 4 | **Polycomb** | H3K27me3+ AND >2kb from TSS |
| 5 | **Active_Enhancer** | H3K27ac+ AND >2kb from TSS |
| 6 | **Poised_Enhancer** | H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb |
| 7 | **CTCF_Site** | CTCF motif+ (early) or CTCF ChIP+ (late) |
| 8 | **Other** | No histone marks AND no CTCF evidence |

### Loop Type Combinations

Loops are classified by combining anchor types, yielding 36 possible combinations:
- Active_Promoter-Active_Promoter
- Active_Promoter-Polycomb
- Polycomb-CTCF_Site
- etc.

### Timepoint-Aware CTCF Classification

- **Early timepoint**: CTCF motif only (no ChIP-seq available)
- **Late timepoint**: CTCF ChIP-seq only (actual binding data)

This avoids mixing developmental stages by using adult CTCF ChIP for early data.

---

## 6. The Emerging Biological Story

### Grad Mentor's Interpretation

> "I think our story's clearing up at the 'shared' anchors - long polycomb loops are collapsing inwards to a closer TAD checkpoint, making a higher density TAD"

This suggests:
- Differential loops involve **shared anchors** (same genomic position in multiple loops)
- **Long-range Polycomb loops** are being lost
- Loops are **collapsing inward** to nearer CTCF/TAD boundary sites
- Result: **higher density TADs** (more loops within shorter distances)

### Mechanistic Chain (Proposed, Requires Validation)

Based on known BAP1 biology and observed results:

```
BAP1 loss
    ↓
H2AK119ub accumulation (BAP1 normally removes this mark)
    ↓
Chromatin compaction / altered Polycomb targeting
    ↓
Two downstream effects:

METHYLATION ARM:
    → TET enzyme access restricted
    → 5mC cannot be oxidized to 5hmC
    → 5mC accumulates, 5hmC depletes
    → Coordinated mC↑ / hmC↓ pattern at gene bodies

CHROMATIN ARCHITECTURE ARM:
    → Long-range Polycomb loops destabilized
    → Loops collapse to nearer TAD boundaries
    → Higher density TAD structure
```

### What's Expected vs. Potentially Novel

**Expected based on literature**:
- High 5hmC levels in brain/cerebellum
- Gene body methylation as primary site of 5hmC
- BAP1 loss causing chromatin changes via H2AK119ub
- Neuronal genes being vulnerable (high baseline 5hmC)

**Potentially novel**:
- The specific coordinated mC↑/hmC↓ pattern (84.6% concordance)
- Connection between BAP1 and TET-mediated demethylation
- Syt1 as top affected gene (first link to synaptic transmission)
- Active Promoter specificity for methylation changes
- The Polycomb loop collapse/TAD densification story (if validated)

---

## 7. Key Numbers to Remember

### Methylation Analysis

| Metric | Value |
|--------|-------|
| Sample size | n=2 per condition (sex confounded) |
| Genes with coordinated mC↑/hmC↓ | **5,708 (84.6%)** |
| Mean 5mC change | **+2.27%** |
| Mean 5hmC change | **-2.08%** |
| Top gene (Syt1) mC change | +17.9% |
| Top gene (Syt1) hmC change | -15.3% |
| DMRs at Active Promoter | **49.9%** |

### Hi-C Analysis

| Metric | Value |
|--------|-------|
| Sample size | n=3 per condition |
| Resolutions analyzed | 5kb, 10kb, 25kb |
| Power improvement vs n=1 | 250-2,500× |

### ChIP-seq Peak Counts (Late Cerebellum)

| Mark | Count |
|------|-------|
| H3K27ac | ~15,000 |
| H3K27me3 | ~15,800 |
| H3K4me1 | ~113,000 |
| H3K4me3 | ~6,500 (single) / ~9,600 (consensus) |
| Bivalent | ~320 (single) / ~690 (consensus) |
| CTCF | 32,487 peaks |

---

## 8. Known Limitations and Caveats

1. **Methylation sample size**: n=2 per condition with sex confounding. More replicates in progress. Current FDR thresholds should be interpreted cautiously.

2. **BAP1 in brain is understudied**: Most BAP1 literature is from cancer. Mechanistic parallels may not translate directly to neurons.

3. **Correlation ≠ Causation**: The coordinated methylation pattern is consistent with TET blockade but doesn't prove it. Additional experiments needed (TET expression, H2AK119ub ChIP-seq in these samples).

4. **Timepoint mixing for ChIP-seq annotations**: Early timepoint uses different ChIP-seq sources than late. CTCF is handled via timepoint-aware classification (motif vs ChIP).

5. **Single tissue**: Results are from cerebellum. Effects may be tissue-specific given cerebellum's uniquely high 5hmC levels.

---

## 9. File Locations and Data Structure

### Methylation Data

```
biomodal/
├── upstream/duet-1.5.0_evoC_Bap1_run_6bp/    # DUET pipeline outputs
│   ├── reports/                               # Summary metrics
│   └── sample_outputs/zarr_store/             # Methylation data
└── downstream/modality/outputs/run-2/         # DMR analysis
    ├── outputs_CG/Results/                    # CG context results
    └── outputs_CHG/Results/                   # CHG context results
```

### Hi-C Data

```
mariner_hi-c/
├── outputs/
│   ├── res_5kb/, res_10kb/, res_25kb/        # Per-resolution intermediates
│   ├── edgeR_results_res_*/                   # Differential analysis
│   ├── merged_loops/                          # Cross-resolution merged
│   ├── visualizations/                        # Publication figures
│   └── apa_results/                           # Aggregate Peak Analysis
└── peaks/
    ├── beds/                                  # ChIP-seq peak files
    └── loop_annotation_extended/              # Chromatin state annotations
```

---

## 10. Questions This Project Aims to Address

1. How does BAP1 loss affect DNA methylation homeostasis in neurons?
2. Is there a connection between H2AK119ub accumulation and TET enzyme access?
3. What happens to chromatin 3D architecture when BAP1/PR-DUB is lost?
4. Are Polycomb-mediated long-range loops specifically affected?
5. What is the relationship between methylation changes and loop changes?
6. Which genes/pathways are most affected and what are the functional consequences?

---

*Document prepared: January 2026*  
*For use as context in AI-assisted analysis discussions*
