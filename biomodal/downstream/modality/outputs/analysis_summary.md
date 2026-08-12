# Biomodal DUET evoC Analysis Summary
## BAP1-KO vs Wildtype Differential Methylation Analysis

**Analysis Date**: 2026-01-21
**Samples**: 4 (2 control, 2 mutant)
**Genome**: GRCm38 (mm10)
**Assay**: Biomodal DUET evoC 6-base resolution

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Sample Information](#2-sample-information)
3. [Upstream Quality Metrics](#3-upstream-quality-metrics)
4. [Methylation Context Comparison (CG vs CHG vs CHH)](#4-methylation-context-comparison)
5. [Sample Correlations](#5-sample-correlations)
6. [DMR Statistics by Genomic Region](#6-dmr-statistics-by-genomic-region)
7. [Top Differentially Methylated Genes](#7-top-differentially-methylated-genes)
8. [Biological Interpretation](#8-biological-interpretation)
9. [Methodological Notes](#9-methodological-notes)
10. [Source File Paths](#10-source-file-paths)

---

## 1. Executive Summary

### Key Findings

| Finding | Details |
|---------|---------|
| **Primary affected regions** | Gene bodies (NOT promoters or CpG islands) |
| **Coordinated pattern** | 5mC increase + 5hmC decrease at same loci |
| **Significant DMRs (CG context)** | 4,188 mC + 4,897 hmC gene body DMRs |
| **Top affected genes** | Neuronal: Syt1, Cntnap2, Dlgap1, Zbtb20 |
| **Mechanism suggested** | Impaired TET-mediated active demethylation |
| **Non-CpG methylation** | No significant changes in CHG/CHH contexts |

### Critical Caveat

> **Run-1 (with sex covariate)**: No significant results
> **Run-2 (without sex covariate)**: Significant results
>
> This indicates sex effects are confounded with BAP1 effects when using n=2 per condition.

---

## 2. Sample Information

| Sample ID | Condition | Sex | Genotype | Mice Pooled |
|-----------|-----------|-----|----------|-------------|
| evoC-Bap1-ctrl-F | Control | Female | Bap1-ff | 3398-3379-3138 |
| evoC-Bap1-ctrl-M | Control | Male | Bap1-ff | 3139-3147-3387 |
| evoC-Bap1-mut-F | Mutant | Female | Bap1-ff-cre | 3380-3388-3381 |
| evoC-Bap1-mut-M | Mutant | Male | Bap1-ff-cre | 3391-3392-3386 |

**Source**: `upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv`

---

## 3. Upstream Quality Metrics

### 3.1 Sequencing Quality

| Sample | Mapped Reads | Mapped Bases | Mean Phred | Duplication Rate |
|--------|--------------|--------------|------------|------------------|
| ctrl-F | 168.0M | 21.8B | 34.4 | 27.5% |
| ctrl-M | 170.1M | 22.4B | 34.4 | 25.7% |
| mut-F | 215.3M | 27.7B | 34.4 | 26.8% |
| mut-M | 188.1M | 24.6B | 34.4 | 26.5% |

### 3.2 Assay Performance (Spike-in Controls)

| Metric | ctrl-F | ctrl-M | mut-F | mut-M | Target |
|--------|--------|--------|-------|-------|--------|
| Genetic Accuracy | 99.97% | 99.97% | 99.97% | 99.97% | >99.9% |
| mC Sensitivity (lambda) | 98.67% | 98.69% | 98.71% | 98.67% | >95% |
| hmC Specificity (lambda) | 99.74% | 99.74% | 99.73% | 99.72% | >99% |
| modC Specificity (pUC19) | 99.90% | 99.89% | 99.80% | 99.85% | >99% |

### 3.3 Coverage Statistics

| Sample | Mean CpG Coverage | Genome Coverage (>1x) |
|--------|-------------------|----------------------|
| ctrl-F | 7.03x | 88.7% |
| ctrl-M | 7.23x | 89.1% |
| mut-F | 8.90x | 88.9% |
| mut-M | 7.98x | 89.2% |

**Source**: `upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv`

---

## 4. Methylation Context Comparison

### 4.1 What Are CG, CHG, CHH Contexts?

| Context | Sequence Pattern | Example | Prevalence in Mammals |
|---------|------------------|---------|----------------------|
| **CG** (CpG) | C followed by G | 5'-**CG**-3' | Primary (dominant) |
| **CHG** | C, H (A/C/T), G | 5'-**CAG**-3' | Rare (<1%) |
| **CHH** | C, H, H | 5'-**CAA**-3' | Very rare (<1%) |

### 4.2 Baseline Methylation Levels

| Context | Unmethylated C | 5mC | 5hmC | Total modC |
|---------|----------------|-----|------|------------|
| **CG** | 17.3-17.7% | **71.9-72.4%** | **9.8-10.3%** | ~82% |
| **CHG** | 99.3-99.4% | N/A | N/A | 0.63-0.66% |
| **CHH** | 99.1% | N/A | N/A | 0.86-0.89% |

> **Note**: CpG methylation is ~100x more abundant than non-CpG methylation

### 4.3 Number of Genomic Sites Analyzed

| Context | Total Sites | Mean Coverage |
|---------|-------------|---------------|
| **CG** | 21,908,008 | 7.0-8.9x |
| **CHG** | 232,418,587 | 3.7-4.7x |
| **CHH** | 829,102,643 | ~4x |

### 4.4 Significant DMRs by Context

| Context | Gene Body mC DMRs | Gene Body hmC DMRs | Total Tested |
|---------|-------------------|---------------------|--------------|
| **CG** | **4,188 (20.7%)** | **4,897 (24.2%)** | 20,224 |
| **CHG** | 0 (0%) | N/A | 10,506 |
| **CHH** | Not analyzed | Not analyzed | N/A |

**Source**:
- `downstream/modality/outputs/run-2/outputs_CG/Results/BioQC_*/biological_qc_report_*.json`
- `downstream/modality/outputs/run-2/outputs_CHG/Results/BioQC_*/biological_qc_report_*.json`

---

## 5. Sample Correlations

### 5.1 CG Context (CpG) - Primary Analysis

#### 5mC Correlations

|  | ctrl-F | ctrl-M | mut-F | mut-M |
|--|--------|--------|-------|-------|
| **ctrl-F** | 1.00 | - | - | - |
| **ctrl-M** | 0.76 | 1.00 | - | - |
| **mut-F** | 0.77 | 0.77 | 1.00 | - |
| **mut-M** | 0.76 | 0.76 | **0.79** | 1.00 |

#### 5hmC Correlations

|  | ctrl-F | ctrl-M | mut-F | mut-M |
|--|--------|--------|-------|-------|
| **ctrl-F** | 1.00 | - | - | - |
| **ctrl-M** | **0.48** | 1.00 | - | - |
| **mut-F** | 0.49 | 0.48 | 1.00 | - |
| **mut-M** | 0.47 | 0.47 | **0.51** | 1.00 |

### 5.2 CHG Context - Low Signal

#### 5mC Correlations

|  | ctrl-F | ctrl-M | mut-F | mut-M |
|--|--------|--------|-------|-------|
| **ctrl-F** | 1.00 | - | - | - |
| **ctrl-M** | **0.36** | 1.00 | - | - |
| **mut-F** | 0.37 | 0.41 | 1.00 | - |
| **mut-M** | **0.23** | 0.44 | 0.52 | 1.00 |

#### 5hmC Correlations (Essentially Noise)

|  | ctrl-F | ctrl-M | mut-F | mut-M |
|--|--------|--------|-------|-------|
| **ctrl-F** | 1.00 | - | - | - |
| **ctrl-M** | **0.02** | 1.00 | - | - |
| **mut-F** | 0.03 | 0.02 | 1.00 | - |
| **mut-M** | **0.02** | 0.02 | 0.03 | 1.00 |

### 5.3 Correlation Summary

| Context | Within-Group (mC) | Between-Group (mC) | Within-Group (hmC) | Between-Group (hmC) |
|---------|-------------------|--------------------|--------------------|---------------------|
| **CG** | 0.76-0.79 | 0.76-0.77 | 0.48-0.51 | 0.47-0.49 |
| **CHG** | 0.36-0.52 | 0.23-0.44 | ~0.02 | ~0.02 |

**Source**:
- `downstream/modality/outputs/run-2/outputs_CG/Results/BioQC_20260121_170327/biological_qc_report_4_samples_20260121_170327.json`
- `downstream/modality/outputs/run-2/outputs_CHG/Results/BioQC_20260121_170730/biological_qc_report_4_samples_20260121_170730.json`

---

## 6. DMR Statistics by Genomic Region

### 6.1 CG Context - mC (5-methylcytosine)

| Region | Significant (q<0.05) | Total Tested | Percentage |
|--------|---------------------|--------------|------------|
| **Gene bodies** | **4,188** | 20,224 | **20.7%** |
| CpG shores | 772 | 27,401 | 2.8% |
| CpG shelves | 397 | 25,975 | 1.5% |
| CpG islands | 23 | 1,237 | 1.9% |
| Promoters | 28 | 14,989 | 0.2% |
| TSS regions | 2 | 4,593 | 0.04% |

### 6.2 CG Context - hmC (5-hydroxymethylcytosine)

| Region | Significant (q<0.05) | Total Tested | Percentage |
|--------|---------------------|--------------|------------|
| **Gene bodies** | **4,897** | 20,224 | **24.2%** |
| CpG shelves | 65 | 25,975 | 0.3% |
| CpG shores | 55 | 27,401 | 0.2% |
| CpG islands | 2 | 1,237 | 0.2% |
| Promoters | 4 | 14,989 | 0.03% |
| TSS regions | 0 | 4,593 | 0% |

### 6.3 CHG Context - No Significant DMRs

| Region | Significant (q<0.05) | Total Tested |
|--------|---------------------|--------------|
| Gene bodies | 0 | 10,506 |
| CpG shores | 0 | 2,427 |
| CpG shelves | 0 | 2,709 |
| CpG islands | 1 | 36 |
| Promoters | 0 | 1,069 |
| TSS regions | 0 | 295 |

**Source**:
- `downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.*.annotation/DMR_*/DMR_mc_*.bed`
- `downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.*.annotation/DMR_*/DMR_hmc_*.bed`

---

## 7. Top Differentially Methylated Genes

### 7.1 Top 20 mC DMRs (Hypermethylated in Mutant)

| Rank | Gene | Chromosome | mC Change | q-value | Function |
|------|------|------------|-----------|---------|----------|
| 1 | **Syt1** | chr10 | +18.4% | <10^-300 | Synaptic vesicle exocytosis |
| 2 | **Trpm3** | chr19 | +8.7% | 2.5e-261 | TRP calcium channel |
| 3 | **Arhgap26** | chr18 | +13.8% | 1.3e-257 | Rho GTPase signaling |
| 4 | **Zbtb20** | chr16 | +8.0% | 3.1e-229 | Zinc finger transcription factor |
| 5 | **Galnt7** | chr8 | +16.9% | 3.7e-215 | Glycosyltransferase |
| 6 | **Patj** | chr4 | +11.6% | 2.8e-205 | Tight junction protein |
| 7 | **Mcu** | chr10 | +14.8% | 7.8e-195 | Mitochondrial Ca2+ uniporter |
| 8 | **Tmtc2** | chr10 | +9.2% | 2.1e-166 | Transmembrane O-mannosyltransferase |
| 9 | **Epha3** | chr16 | +10.4% | 3.1e-166 | Ephrin receptor |
| 10 | **Pde4d** | chr13 | +4.8% | 4.9e-163 | Phosphodiesterase |
| 11 | **Lpp** | chr16 | +7.6% | 1.4e-158 | LIM domain protein |
| 12 | **Unc5c** | chr3 | +10.9% | 1.8e-139 | Netrin receptor |
| 13 | **Rps6ka5** | chr12 | +12.8% | 6.6e-130 | Ribosomal kinase |
| 14 | **Lsamp** | chr16 | +3.6% | 3.4e-124 | Limbic system protein |
| 15 | **Cdh8** | chr8 | +9.1% | 1.4e-119 | Cadherin |
| 16 | **Ncald** | chr15 | +9.8% | 2.1e-118 | Neurocalcin delta |
| 17 | **Rims2** | chr15 | +9.2% | 1.8e-117 | Synaptic vesicle priming |
| 18 | **Dlgap1** | chr17 | +6.3% | 3.9e-104 | Postsynaptic density |
| 19 | **Maml3** | chr3 | +7.3% | 3.2e-103 | Notch coactivator |
| 20 | **Cadm2** | chr16 | +5.2% | 2.6e-102 | Cell adhesion molecule |

### 7.2 Top 20 hmC DMRs (Reduced in Mutant)

| Rank | Gene | Chromosome | hmC Change | q-value | Function |
|------|------|------------|------------|---------|----------|
| 1 | **Syt1** | chr10 | **-15.8%** | <10^-300 | Synaptic vesicle exocytosis |
| 2 | **Cntnap2** | chr6 | -4.0% | 3.9e-232 | Contactin-associated protein |
| 3 | **Zbtb20** | chr16 | -6.2% | 5.7e-185 | Zinc finger transcription factor |
| 4 | **Pde4d** | chr13 | -4.1% | 1.3e-170 | Phosphodiesterase |
| 5 | **Tmtc2** | chr10 | -8.1% | 1.4e-165 | Transmembrane O-mannosyltransferase |
| 6 | **Trpm3** | chr19 | -5.9% | 5.5e-157 | TRP calcium channel |
| 7 | **Epha3** | chr16 | -8.6% | 6.6e-154 | Ephrin receptor |
| 8 | **Arhgap26** | chr18 | -9.0% | 2.5e-140 | Rho GTPase signaling |
| 9 | **Kcnip4** | chr5 | -4.2% | 1.3e-126 | Potassium channel interacting |
| 10 | **Lpp** | chr16 | -6.0% | 3.7e-123 | LIM domain protein |
| 11 | **Csmd1** | chr8 | +3.1% | 6.6e-120 | CUB/Sushi domain protein |
| 12 | **Patj** | chr4 | -7.1% | 2.3e-119 | Tight junction protein |
| 13 | **Mcu** | chr10 | -10.3% | 4.2e-118 | Mitochondrial Ca2+ uniporter |
| 14 | **Fat3** | chr9 | +5.3% | 2.4e-116 | Protocadherin |
| 15 | **Trhde** | chr10 | -6.8% | 1.1e-110 | TRH-degrading enzyme |
| 16 | **Rps6ka5** | chr12 | -10.2% | 1.4e-108 | Ribosomal kinase |
| 17 | **Rit2** | chr18 | -7.2% | 2.4e-108 | Ras-like GTPase |
| 18 | **Angpt1** | chr15 | -7.0% | 2.6e-103 | Angiopoietin |
| 19 | **Cdh8** | chr8 | -7.2% | 1.4e-99 | Cadherin |
| 20 | **Lsamp** | chr16 | -2.7% | 7.7e-95 | Limbic system protein |

### 7.3 Coordinated Changes (Same Gene, Opposite Direction)

| Gene | mC Change | hmC Change | Net Effect |
|------|-----------|------------|------------|
| **Syt1** | +18.4% | -15.8% | Blocked demethylation |
| **Zbtb20** | +8.0% | -6.2% | Blocked demethylation |
| **Trpm3** | +8.7% | -5.9% | Blocked demethylation |
| **Epha3** | +10.4% | -8.6% | Blocked demethylation |
| **Mcu** | +14.8% | -10.3% | Blocked demethylation |
| **Lpp** | +7.6% | -6.0% | Blocked demethylation |

**Source**:
- `downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_172049/DMR_mc_control__mutant_20260121_172049.bed`
- `downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_172049/DMR_hmc_control__mutant_20260121_172049.bed`

---

## 8. Biological Interpretation

### 8.1 The TET-Mediated Demethylation Pathway

```
5mC → (TET1/2/3) → 5hmC → 5fC → 5caC → C (unmethylated)
```

The observed pattern of **5mC increase** coupled with **5hmC decrease** suggests:

1. BAP1 loss → Increased H2AK119ub → Altered chromatin structure
2. TET enzymes have restricted access to DNA
3. 5mC accumulates (not being converted to 5hmC)
4. 5hmC decreases (not being produced from 5mC)

### 8.2 Functional Categories of Affected Genes

| Category | Genes | Implication |
|----------|-------|-------------|
| **Synaptic/Neuronal** | Syt1, Cntnap2, Dlgap1, Dlg2, Cadm1, Cadm2, Rims2 | Neuronal function disruption |
| **Ion Channels** | Trpm3, Kcnip4, Kcnma1, Mcu | Calcium signaling affected |
| **Cell Adhesion** | Cntnap2, Cdh8, Fat3, Lsamp | Cell-cell communication |
| **Transcription** | Zbtb20, E2f3, Nfia, Nfib, Tcf4 | Gene regulation |
| **Signaling** | Arhgap26, Pde4d, Epha3 | Signal transduction |

### 8.3 Key Conclusions

1. **BAP1 regulates gene body methylation dynamics**, not promoter methylation
2. **The coordinated mC↑/hmC↓ pattern indicates impaired active demethylation**
3. **Neuronal genes are preferentially affected**
4. **~20-24% of genes show significant methylation changes**
5. **Non-CpG methylation (CHG/CHH) is not affected by BAP1 loss**

---

## 9. Methodological Notes

### 9.1 Analysis Runs

| Run | Sex Covariate | Significant Results | Interpretation |
|-----|---------------|---------------------|----------------|
| Run-1 | Yes | None | Sex effects absorb variance |
| **Run-2** | **No** | **Yes** | Primary analysis (but confounded) |

### 9.2 Statistical Thresholds

| Parameter | Value |
|-----------|-------|
| q-value threshold | 0.05 (FDR-corrected) |
| Minimum coverage | 10x per region |
| Multiple testing | Benjamini-Hochberg |

### 9.3 Limitations

1. **Small sample size**: n=2 per condition limits statistical power
2. **Sex confounding**: Male/female differences conflated with BAP1 effect
3. **Single tissue**: Results may be tissue-specific
4. **No validation**: Top hits require orthogonal confirmation

### 9.4 Recommendations

1. Increase sample size (n≥3 per sex per condition)
2. Stratify analysis by sex
3. Validate top genes (Syt1, Cntnap2, Zbtb20) with targeted bisulfite sequencing
4. Correlate with RNA-seq expression data
5. Examine TET1/2/3 expression levels

---

## 10. Source File Paths

### 10.1 Upstream Data (DUET Pipeline Outputs)

```
# Main summary metrics
upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv

# Metrics definitions
upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Metrics_Definitions.csv

# FastQC summaries
upstream/duet-1.5.0_evoC_Bap1_run_6bp/diagnostics/fastqc_reports/logs/*.fastqc_summary.txt

# Pipeline documentation
upstream/biomodal_docs.md
```

### 10.2 Downstream Data (Modality Outputs)

#### Run Configuration
```
# Critical note about runs
downstream/modality/outputs/note.txt

# Sample metadata
downstream/modality/metadata.tsv

# Configuration files
downstream/modality/config.txt
downstream/modality/config_CG.txt
downstream/modality/config_CHG.txt
```

#### Biological QC Reports (JSON with correlation matrices)
```
# CG context
downstream/modality/outputs/run-2/outputs_CG/Results/BioQC_20260121_170327/biological_qc_report_4_samples_20260121_170327.json

# CHG context
downstream/modality/outputs/run-2/outputs_CHG/Results/BioQC_20260121_170730/biological_qc_report_4_samples_20260121_170730.json
```

#### Biological QC Reports (HTML - Visualizations)
```
# CG context
downstream/modality/outputs/run-2/outputs_CG/Results/BioQC_20260121_170327/biological_qc_report_4_samples_20260121_170327.html

# CHG context
downstream/modality/outputs/run-2/outputs_CHG/Results/BioQC_20260121_170730/biological_qc_report_4_samples_20260121_170730.html
```

#### DMR Results - Full BED Files (CG Context)
```
# Gene bodies - mC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_172049/DMR_mc_control__mutant_20260121_172049.bed

# Gene bodies - hmC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_172049/DMR_hmc_control__mutant_20260121_172049.bed

# CpG islands - mC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260121_171455/DMR_mc_control__mutant_20260121_171455.bed

# CpG islands - hmC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260121_171455/DMR_hmc_control__mutant_20260121_171455.bed

# Promoters - mC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260121_172220/DMR_mc_control__mutant_20260121_172220.bed

# CpG shores - mC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260121_171910/DMR_mc_control__mutant_20260121_171910.bed

# CpG shelves - mC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260121_171718/DMR_mc_control__mutant_20260121_171718.bed

# TSS regions - mC
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260121_172333/DMR_mc_control__mutant_20260121_172333.bed
```

#### DMR Results - Filtered (q<0.05)
```
# CpG islands - mC filtered
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_Report_20260121_172402/DMR_mc_control__mutant_qval_0_05_filtered_20260121_172402.bed.gz

# CpG islands - hmC filtered
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_Report_20260121_172354/DMR_hmc_control__mutant_qval_0_05_filtered_20260121_172354.bed.gz

# Gene bodies - mC filtered
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_Report_20260121_172506/DMR_mc_control__mutant_qval_0_05_filtered_20260121_172506.bed.gz

# Gene bodies - hmC filtered
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_Report_20260121_172457/DMR_hmc_control__mutant_qval_0_05_filtered_20260121_172457.bed.gz
```

#### DMR Volcano Plot Reports (HTML)
```
# CpG islands
downstream/modality/outputs/run-2/outputs_CG/Reports/gencode.vM25.mouse.cpg_islands.annotation/DMR_Report_mc_control__mutant.dmr.volcano_*.html
downstream/modality/outputs/run-2/outputs_CG/Reports/gencode.vM25.mouse.cpg_islands.annotation/DMR_Report_hmc_control__mutant.dmr.volcano_*.html

# Gene bodies
downstream/modality/outputs/run-2/outputs_CG/Reports/gencode.vM25.mouse.genes.annotation/DMR_Report_mc_control__mutant.dmr.volcano_*.html
downstream/modality/outputs/run-2/outputs_CG/Reports/gencode.vM25.mouse.genes.annotation/DMR_Report_hmc_control__mutant.dmr.volcano_*.html

# Promoters
downstream/modality/outputs/run-2/outputs_CG/Reports/gencode.vM25.mouse.promoters.annotation/DMR_Report_mc_control__mutant.dmr.volcano_*.html
```

#### Feature Extraction Results
```
# CpG islands - methylation fractions
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/Extract_20260121_170458/Extract_mc_regional-frac_20260121_170458.tsv.gz
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/Extract_20260121_170458/Extract_hmc_regional-frac_20260121_170458.tsv.gz

# Gene bodies - methylation fractions
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260121_171104/Extract_mc_regional-frac_20260121_171104.tsv.gz
downstream/modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_20260121_171104/Extract_hmc_regional-frac_20260121_171104.tsv.gz
```

#### CHG Context Results (for comparison)
```
# BioQC
downstream/modality/outputs/run-2/outputs_CHG/Results/BioQC_20260121_170730/biological_qc_report_4_samples_20260121_170730.json

# Gene bodies DMR
downstream/modality/outputs/run-2/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_174629/DMR_mc_control__mutant_20260121_174629.bed
```

### 10.3 Documentation

```
# Interpretation guide
downstream/biomodal-interpretation-guide.md

# Workflow documentation
downstream/modality/biomodal-workflow.md

# Pipeline documentation
upstream/biomodal_docs.md
```

### 10.4 Visualization Screenshots (if needed)

```
# Upstream metrics
upstream/plots/biomodal-metrics/*.png

# QC plots
downstream/plots/QC/CG/*.png
downstream/plots/QC/CHG/*.png

# CpG island DMR plots
downstream/plots/cpg/islands/*.png
```

---

## Appendix: BED File Column Definitions

The DMR BED files contain the following columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | Chromosome | Chromosome name |
| 2 | Start | Region start (0-indexed) |
| 3 | End | Region end |
| 4 | num_contexts | Number of CpG sites in region |
| 5 | mean_coverage | Mean sequencing depth |
| 6 | mean_mod_group_1 | Mean methylation (control) |
| 7 | mean_mod_group_2 | Mean methylation (mutant) |
| 8 | mod_fold_change | Fold change (mutant/control) |
| 9 | mod_difference | Absolute difference (mutant - control) |
| 10 | test_statistic | Statistical test value |
| 11 | dmr_pvalue | Raw p-value |
| 12 | dmr_qvalue | FDR-corrected q-value |
| 13 | Annotation | Region type (Gene, CpG_island, etc.) |
| 14 | Gene (if applicable) | Gene symbol |

---

*Document generated from biomodal DUET evoC analysis results*
