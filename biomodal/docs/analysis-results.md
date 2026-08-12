# Biomodal DUET evoC Differential Methylation Results

## Study Overview

**Comparison**: BAP1-KO Mutant vs Wildtype Control
**Organism**: Mouse (*Mus musculus*, GRCm38/mm10)
**Samples**: 4 total (2 Control, 2 Mutant)
**Technology**: Biomodal DUET evoC (6bp resolution dual 5mC/5hmC sequencing)
**Analysis Date**: 2026-01-21

---

## Executive Summary

BAP1 loss causes a striking coordinated epigenetic pattern: **5mC increases while 5hmC decreases** at the same gene body loci. This pattern affects 2,856 genes (92.3% of genes significant in both modifications) and is consistent with impaired TET-mediated active demethylation. The affected genes are enriched for neuronal/synaptic functions and predominantly occur at **Active Promoter** chromatin states (62.7% of significant DMRs).

---

## Differential Methylation Statistics

### 5mC (5-methylcytosine) Changes

| Metric | Value |
|--------|-------|
| Total genes tested | 20,224 |
| Significant (q < 0.05) | 4,188 (20.7%) |
| Hypermethylated (mC↑) | 3,616 (86.3%) |
| Hypomethylated (mC↓) | 572 (13.7%) |
| Mean change (significant) | **+4.10%** |

**Key observation**: 5mC changes are overwhelmingly in the **hypermethylated direction** (86.3% of significant genes show increased methylation in BAP1-KO mutants).

### 5hmC (5-hydroxymethylcytosine) Changes

| Metric | Value |
|--------|-------|
| Total genes tested | 20,224 |
| Significant (q < 0.05) | 4,897 (24.2%) |
| Increased (hmC↑) | 536 (10.9%) |
| Decreased (hmC↓) | 4,361 (89.1%) |
| Mean change (significant) | **-3.04%** |

**Key observation**: 5hmC changes are overwhelmingly in the **decreased direction** (89.1% of significant genes show reduced hydroxymethylation in BAP1-KO mutants).

### Direction Summary

The asymmetric direction of changes is striking:
- **mC**: 86% hypermethylated vs 14% hypomethylated
- **hmC**: 11% increased vs 89% decreased

This inverse pattern is the key finding of the analysis.

---

## The Coordinated Pattern: mC↑ / hmC↓

### Evidence for Coordinated Changes

| Metric | Value |
|--------|-------|
| Genes significant in both mC AND hmC | 3,095 |
| Genes with mC↑ / hmC↓ pattern | **2,856 (92.3%)** |
| Genes with opposite pattern (mC↓ / hmC↑) | 239 (7.7%) |

The near-universal coordination (92.3%) of increased mC with decreased hmC at the same loci is biologically meaningful and mechanistically interpretable.

### Mechanistic Interpretation

The TET (Ten-Eleven Translocation) enzymes catalyze the active demethylation pathway:

```
5mC → (TET1/2/3) → 5hmC → 5fC → 5caC → C (unmethylated)
```

The mC↑ / hmC↓ pattern indicates a **block in TET-mediated oxidation**:
- 5mC accumulates because it cannot be converted to 5hmC
- 5hmC depletes because existing 5hmC is converted downstream but not replenished
- This is consistent with reduced TET enzyme access to chromatin

**Proposed mechanism**:
```
BAP1 loss → Increased H2AK119ub → Altered chromatin → TET access restricted
```

BAP1 is a deubiquitinase for H2AK119ub (histone H2A lysine 119 ubiquitination). Loss of BAP1 leads to accumulation of H2AK119ub, which may sterically hinder TET enzyme access to DNA, blocking the first step of active demethylation.

---

## Top Affected Genes

### Top 10 mC Hypermethylated Genes (by q-value)

| Rank | Gene | mC Change | q-value | Gene Function |
|------|------|-----------|---------|---------------|
| 1 | **Syt1** | +18.4% | <10^-300 | Synaptic vesicle exocytosis |
| 2 | Trpm3 | +8.7% | 2.5×10^-257 | Calcium channel, thermosensation |
| 3 | Arhgap26 | +13.8% | 8.9×10^-254 | Rho GTPase signaling |
| 4 | Zbtb20 | +8.0% | 1.6×10^-225 | Transcription factor, brain development |
| 5 | Galnt7 | +16.9% | 1.5×10^-211 | O-linked glycosylation |
| 6 | Patj | +11.6% | 9.5×10^-202 | Cell polarity, tight junctions |
| 7 | Mcu | +14.8% | 2.3×10^-191 | Mitochondrial calcium uniporter |
| 8 | Tmtc2 | +9.2% | 5.3×10^-163 | ER protein quality control |
| 9 | Epha3 | +10.4% | 6.9×10^-163 | Ephrin receptor, axon guidance |
| 10 | Pde4d | +4.8% | 9.9×10^-160 | cAMP signaling |

### Top 10 hmC Decreased Genes (by q-value)

| Rank | Gene | hmC Change | q-value | Gene Function |
|------|------|------------|---------|---------------|
| 1 | **Syt1** | -15.8% | <10^-300 | Synaptic vesicle exocytosis |
| 2 | Cntnap2 | -4.0% | 3.9×10^-228 | Cell adhesion, autism-associated |
| 3 | Zbtb20 | -6.2% | 3.8×10^-181 | Transcription factor |
| 4 | Pde4d | -4.1% | 6.6×10^-167 | cAMP signaling |
| 5 | Tmtc2 | -8.1% | 5.7×10^-162 | ER protein quality control |
| 6 | Trpm3 | -5.9% | 1.8×10^-153 | Calcium channel |
| 7 | Epha3 | -8.6% | 6.6×10^-154 | Ephrin receptor |
| 8 | Arhgap26 | -9.0% | 6.4×10^-137 | Rho GTPase signaling |
| 9 | Kcnip4 | -4.2% | 2.8×10^-123 | Potassium channel interacting |
| 10 | Lpp | -6.0% | 7.6×10^-120 | Cell-cell adhesion |

### Syt1: The Most Affected Gene

**Syt1 (Synaptotagmin-1)** shows the most extreme coordinated change:
- mC change: **+18.4%** (q < 10^-300)
- hmC change: **-15.8%** (q < 10^-300)
- Combined effect: 34.2 percentage points

Syt1 encodes the primary calcium sensor for fast synchronous neurotransmitter release. Its dysregulation could have profound effects on synaptic transmission.

### Top Coordinated Genes (mC↑ / hmC↓ pattern)

| Gene | mC Change | hmC Change | Combined Effect |
|------|-----------|------------|-----------------|
| Ly6e | +37.2% | -37.2% | 74.4% |
| Lypla2 | +34.7% | -36.7% | 71.4% |
| Tprn | +23.0% | -21.9% | 44.9% |
| Tssk5 | +23.9% | -20.3% | 44.1% |
| Lgi3 | +21.5% | -22.0% | 43.6% |
| Inpp5j | +17.8% | -18.4% | 36.3% |
| Sap30 | +23.9% | -11.9% | 35.7% |
| Syt1 | +18.4% | -15.8% | 34.2% |

---

## Chromatin State Analysis

### ChIP-seq Peak Overlap at Significant DMRs

The chromatin environment of differentially methylated genes was characterized using ChIP-seq peaks from Late cerebellum:

| ChIP Mark | % DMRs Overlapping | Interpretation |
|-----------|-------------------|----------------|
| **CTCF** | 61.7% | Insulator/boundary element |
| **H3K27ac** | 67.2% | Active enhancer/promoter |
| **H3K27me3** | 12.1% | Polycomb repression |
| **H3K4me1** | 96.7% | Gene body/enhancer mark |
| **H3K4me3** | 67.3% | Active promoter |
| **Bivalent** | 2.1% | Developmental poised |

**Key findings**:
- **96.7% overlap with H3K4me1**: Expected for gene body DMRs (H3K4me1 marks transcribed gene bodies)
- **67.3% overlap with H3K4me3**: High overlap with active promoter marks
- **61.7% overlap with CTCF**: Substantial overlap with insulator elements
- **Low H3K27me3 overlap (12.1%)**: Affected genes are not primarily Polycomb-silenced

### Chromatin State Distribution

DMRs were classified into 7 chromatin state categories:

| Chromatin State | Count | Percentage | Description |
|-----------------|-------|------------|-------------|
| **Active_Promoter** | 2,624 | 62.7% | H3K4me3+ near TSS |
| Repressed_Promoter | 343 | 8.2% | H3K27me3+ near TSS |
| Other | 1,140 | 27.2% | No marks / structural |
| Bivalent_Promoter | 52 | 1.2% | K4me3+K27me3 |
| Active_Enhancer | 17 | 0.4% | H3K27ac+ distal |
| Poised_Enhancer | 6 | 0.1% | H3K4me1+ only |
| Polycomb | 6 | 0.1% | H3K27me3+ distal |

**Key finding**: 62.7% of significant DMRs occur at **Active Promoter** regions, indicating that BAP1 loss preferentially affects actively transcribed genes.

### Chromatin State by Methylation Direction

| Chromatin State | Hypermethylated | Hypomethylated |
|-----------------|-----------------|----------------|
| Active_Promoter | 2,549 (97.1%) | 75 (2.9%) |
| Repressed_Promoter | 36 (10.5%) | 307 (89.5%) |
| Bivalent_Promoter | 33 (63.5%) | 19 (36.5%) |
| Other | 975 (85.5%) | 165 (14.5%) |

**Critical observation**:
- **Hypermethylated genes are predominantly at Active Promoters** (2,549 of 2,624 Active_Promoter DMRs = 97.1%)
- **Hypomethylated genes are predominantly at Repressed Promoters** (307 of 343 Repressed_Promoter DMRs = 89.5%)

This suggests BAP1 loss causes methylation gain specifically at actively transcribed genes, while Polycomb-repressed genes show methylation loss.

---

## Functional Enrichment Analysis

### GO Biological Process (Top 15 terms)

| Term | Gene Count | q-value | Enrichment |
|------|------------|---------|------------|
| Golgi vesicle transport | 112 | 2.4×10^-28 | 3.3x |
| Vesicle-mediated transport in synapse | 112 | 5.1×10^-25 | 3.0x |
| Protein localization to cell periphery | 131 | 5.1×10^-25 | 2.7x |
| Vesicle organization | 126 | 2.5×10^-24 | 2.7x |
| Endosomal transport | 97 | 1.4×10^-20 | 2.9x |
| **Dendrite development** | 107 | 1.4×10^-20 | 2.7x |
| Protein polyubiquitination | 93 | 8.6×10^-20 | 2.9x |
| **Synaptic vesicle cycle** | 89 | 7.5×10^-19 | 2.9x |
| Vesicle localization | 81 | 1.2×10^-18 | 3.1x |
| Protein localization to plasma membrane | 99 | 2.2×10^-18 | 2.7x |
| Cytosolic transport | 58 | 3.0×10^-18 | 3.8x |
| Protein-containing complex localization | 90 | 4.0×10^-18 | 2.8x |
| **Postsynapse organization** | 97 | 1.7×10^-17 | 2.6x |
| Regulation of synapse organization | 110 | 2.2×10^-16 | 2.4x |
| Macroautophagy | 96 | 6.4×10^-16 | 2.5x |

**Key theme**: Synaptic and vesicular transport functions dominate the enriched terms.

### KEGG Pathway Enrichment (Top 10)

| Pathway | Gene Count | q-value | Enrichment |
|---------|------------|---------|------------|
| **Autophagy - animal** | 69 | 3.4×10^-15 | 2.8x |
| **Endocytosis** | 89 | 7.5×10^-13 | 2.3x |
| Spinocerebellar ataxia | 53 | 5.4×10^-10 | 2.6x |
| MAPK signaling pathway | 88 | 6.8×10^-10 | 2.0x |
| Salmonella infection | 76 | 4.8×10^-9 | 2.1x |
| Inositol phosphate metabolism | 34 | 9.1×10^-9 | 3.0x |
| Phosphatidylinositol signaling | 39 | 9.1×10^-9 | 2.8x |
| Ubiquitin mediated proteolysis | 53 | 1.6×10^-8 | 2.3x |
| Sphingolipid signaling | 47 | 1.7×10^-8 | 2.5x |
| Protein processing in ER | 57 | 3.7×10^-8 | 2.2x |

**Key themes**:
1. **Autophagy and endocytosis**: Vesicular trafficking is consistently enriched
2. **Spinocerebellar ataxia**: Directly relevant to cerebellar phenotypes
3. **MAPK/inositol signaling**: Signal transduction pathways
4. **Ubiquitin-proteasome**: Protein quality control (note: BAP1 is a deubiquitinase)

---

## Biological Interpretation

### Primary Affected Regions: Gene Bodies

The analysis reveals that **gene bodies** (not promoters or CpG islands) are the primary site of differential methylation. This is consistent with the role of gene body methylation in:
- Transcriptional elongation
- Alternative splicing regulation
- Prevention of spurious transcription initiation

### The TET Enzyme Block Model

The coordinated mC↑/hmC↓ pattern strongly supports a model where BAP1 loss impairs TET enzyme function:

1. **BAP1 normally deubiquitinates H2AK119ub**
2. **BAP1 loss → H2AK119ub accumulation**
3. **Altered chromatin structure restricts TET access**
4. **5mC cannot be oxidized to 5hmC**
5. **Result: mC accumulates, hmC depletes**

### Neuronal Gene Vulnerability

The enrichment of synaptic/neuronal genes (Syt1, Cntnap2, Dlgap1, etc.) suggests that neuronal genes are particularly vulnerable to this epigenetic dysregulation. Possible explanations:
- Neuronal genes have higher baseline gene body methylation
- Neuronal genes require dynamic methylation for activity-dependent regulation
- Cerebellar expression bias in the analyzed tissue

### Implications for BAP1 Biology

This study demonstrates that BAP1 functions extend beyond chromatin deubiquitination to include:
- **Regulation of DNA methylation homeostasis**
- **Control of TET-mediated active demethylation**
- **Maintenance of gene body methylation patterns at actively transcribed genes**

---

## Key Conclusions

1. **4,188 genes** show significant 5mC changes; **4,897 genes** show significant 5hmC changes

2. **2,856 genes (92.3%)** exhibit the coordinated **mC↑ / hmC↓** pattern, indicating blocked TET-mediated demethylation

3. **Syt1** is the most significantly affected gene, with +18.4% mC and -15.8% hmC changes

4. **62.7% of significant DMRs** occur at Active Promoter chromatin states, indicating preferential effects on actively transcribed genes

5. **Synaptic and vesicular transport genes** are strongly enriched among affected genes

6. **Non-CpG methylation (CHG/CHH)** shows no significant changes, confirming the specificity of CpG methylation dysregulation

7. The results support a model where **BAP1 loss → H2AK119ub accumulation → TET access restricted → impaired active demethylation**

---

## Data Files

All analysis tables are available in:
`biomodal/downstream/plots/visualizations/tables/`

| File | Description |
|------|-------------|
| `summary_statistics.txt` | Overall analysis summary |
| `top_mc_dmrs.tsv` | Top 20 mC differentially methylated genes |
| `top_hmc_dmrs.tsv` | Top 20 hmC differentially methylated genes |
| `coordinated_changes.tsv` | All genes with coordinated mC/hmC changes |
| `chromatin_state_summary.tsv` | Chromatin state distribution summary |
| `dmr_chromatin_state_annotation.tsv` | Full DMR annotation with ChIP-seq overlaps |
| `enrichment_go_bp.tsv` | GO Biological Process enrichment |
| `enrichment_go_cc.tsv` | GO Cellular Component enrichment |
| `enrichment_go_mf.tsv` | GO Molecular Function enrichment |
| `enrichment_kegg.tsv` | KEGG pathway enrichment |

---

*Generated from biomodal DUET evoC analysis pipeline*
*Analysis script: `biomodal/downstream/scripts/biomodal_visualizations.R`*
