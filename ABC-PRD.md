# PRD: ABC Model Enhancer-Gene Linkage Analysis

## Purpose

Use the Activity-By-Contact (ABC) model to link enhancers to genes, then demonstrate that enhancer-gene linkages are dysregulated in BAP1-KO conditions and correlate with differentially expressed genes.

The ABC model formula: **ABC score = (Enhancer Activity) × (Contact Frequency)**

Where:
- Activity = H3K27ac signal at enhancer
- Contact = Hi-C contact frequency between enhancer and promoter

---

## Tool Selection

**Broad Institute ABC Implementation**
- Repository: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
- Version: v1.1.2 (March 2024)
- Documentation: https://abc-enhancer-gene-prediction.readthedocs.io/
- Reference paper: Fulco et al. 2019 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6886585/)

### Installation

Conda/mamba environment creation (~15 minutes with mamba). Dependencies include Python 3.9+, samtools, bedtools, MACS2, juicertools for Hi-C processing. Snakemake workflow supports parallel execution via `-j` flag.

### Tool Input Requirements

| Input | Format | Required |
|-------|--------|----------|
| Chromatin accessibility | DNase-seq BAM or ATAC-seq tagAlign | Yes |
| H3K27ac ChIP-seq | bigWig | Yes (or ATAC-only mode) |
| Hi-C contact data | `.hic` format or BEDPE | Optional (powerlaw distance model substitutes) |
| Gene annotations | BED/TSV with TSS coordinates | Yes |

### Tool Outputs

- `EnhancerPredictions_threshold_*.tsv` containing:
  - Enhancer-gene pairs
  - ABC scores
  - Activity components
  - Contact frequencies
  - Distances
- Recommended ABC score threshold: **0.015** (achieves ~70% recall)
- Expected output: ~3 distal enhancers per expressed gene

---

## Existing Project Resources

### H3K27ac Peaks (Enhancer Activity Proxy)
- `peaks/beds/H3K27acCerebellumLate2.bed`
- `peaks/beds/H3K27acCerebellumEarly2.bed`

### ATAC-seq Data
- Differential peaks available as up/down BED files
- **Consensus ATAC bed files:** Challana to provide (intersect across replicates, merge across genotype)
- Individual peak files on instance: `/data2/rs_256/Func_annotation_v2/subtracted_bedfiles`
- Will need conversion to tagAlign format for ABC input

### Loop Contact Data
- `characterized_loops.tsv`
  - Contains: `logFC`, `anchor1_H3K27ac_overlap`, `anchor2_H3K27ac_overlap`, `anchor1_type`, `loop_type`
  - Gene annotations included: `anchor1_nearest_gene`, `anchor2_nearest_gene`, `anchor1_distance_to_tss`

### RNA-seq Differential Expression
- `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
- `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx`

### H2AK119ub ChIP-seq (bigWig)
Location: `heatmaps/`

**Ctrl (WT) — 4 usable replicates:**
- `H2AK119ubCtrl1.bw`
- `H2AK119ubCtrl2.bw`
- `H2AK119ubCtrl3.bw`
- `H2AK119ubCtrl4.bw`

**Mut (KO) — 4 usable replicates:**
- `H2AK119ubMut1.bw`
- `H2AK119ubMut2.bw`
- `H2AK119ubMut3.bw`
- `H2AK119ubMut4.bw`

**Excluded (bad agreement with prior data):**
- `H2AK119ubCtrl5-badagreementwithpriordata.bw`
- `H2AK119ubCtrl6-badagreementwithpriordata.bw`
- `H2AK119ubMut5-badagreementwithpriordata.bw`
- `H2AK119ubMut6-badagreementwithpriordata.bw`

---

## Tasks

### 7a. Implement ABC Model Scoring
Compute ABC score = (enhancer Activity) × (Contact frequency) for each enhancer-gene pair.
- Activity = H3K27ac signal at enhancer
- Contact = Hi-C contact frequency between enhancer and promoter (from loop logFC or raw counts)

### 7b. Compute Delta Enhancer-Promoter Contact Per Gene
Using the differential loop data, calculate change in contact between each promoter and its linked enhancers.

### 7c. Correlate Delta Contacts with Differential Gene Expression
Plot delta E-P contact vs log2FC from RNA-seq.

### 7d. Identify Dysregulated Enhancer-Gene Linkages
Find E-P pairs where both contact and expression change significantly.

### 7e. Super-Enhancer to Gene Linkages
Identify loops where both anchors are H3K27ac+ (super-enhancer signature), link to target genes, assess expression changes.

### 7f. Correlate Delta E-P Contacts with H2AK119ub Levels (2nd Goal)
Test hypothesis: "ub is buffer to stop k27ac contact" — once ubiquitination threshold is reached, contacts form.

Steps:
1. Quantify H2AK119ub signal at enhancers and promoters using `bigWigAverageOverBed` or deepTools `computeMatrix`
2. Calculate differential ub levels (Mut vs Ctrl) at E-P anchor regions
3. Correlate Δub with ΔE-P contact (from task 7b)
4. Test whether loss of ub (in BAP1-KO) corresponds to gain of E-P contact

---

## Workflow

The ABC model has no built-in differential analysis. Predictions must be run separately per condition, then integrated.

### Critical: Consensus Enhancer Universe (NOT Cell-Type Specific)

**IMPORTANT:** Many papers use ABC to define cell-type specific enhancers by calling separate enhancer sets per condition. **This is NOT what we want.** Defining separate enhancers between ctrl and mut would confound differential analysis.

**Correct approach:**
1. Define a **consensus enhancer universe** shared across conditions
2. Generate ATAC-seq consensus: **intersect across replicates, merge across genotype**
   - Double-check that intersecting with a particular sample isn't removing too many peaks
3. Filter to elements within 5-100 kb (or up to 200 kb) of gene TSSs
4. Quantify enhancer **activity using H3K27ac AND ATAC-seq signal** per condition
5. Run ABC with this fixed enhancer set for both WT and KO, contact from Hi-C

This ensures we're comparing the same genomic elements across conditions.

### Recommended Workflow (CCBB Integration Approach)

```
1. Define consensus enhancer universe
   - Intersect ATAC-seq peaks across replicates (within each genotype)
   - Merge across genotypes (WT ∪ KO)
   - QC: Check that intersecting doesn't remove too many peaks from any sample
   - Filter to regions within 5-200 kb of gene TSSs
   - Cap to top 100-150k enhancers by H3K27ac activity (if runtime/compute issues)
   - Note: Challana will provide ATAC bed files

2. Run ABC separately per condition
   - Input: Consensus ATAC set + H3K27ac + Hi-C for WT
   - Input: Consensus ATAC set + H3K27ac + Hi-C for KO
   - Activity = H3K27ac signal AND ATAC signal
   - Contact = Hi-C
   - Use same consensus enhancer set for both conditions

3. Filter ABC predictions
   - Keep links with ABC ≥ 0.02 per condition

4. Compute delta ABC
   - ΔABC = ABC_KO - ABC_WT for each enhancer-gene pair

5. Integrate with RNA-seq
   - Join ΔABC with gene log2FC from RNA-seq
   - Summarize directional concordance (ΔABC↑ + log2FC↑, etc.)

6. Aggregate to gene level
   - Strongest link per gene (max |ΔABC|)
   - Total ABC per gene (sum of all links)
   - Number of perturbed links per gene

7. (Later) Correlate with H2AK119ub - interpret results first before proceeding
```

### Alternative Workflow (Original - Modified)

```
1. DiffBind → Identify differential H3K27ac regions (for downstream interpretation)
2. ABC Model → Generate enhancer-gene predictions per condition (using consensus enhancer set)
3. Custom integration →
   - Calculate ΔABC scores between conditions
   - Intersect with DiffBind differential peaks
   - Correlate enhancer logFC with target gene logFC from RNA-seq
4. Filter → Prioritize concordant enhancer-gene pairs
```

### Differential Analysis Tools (if needed)

| Tool | Purpose | Notes |
|------|---------|-------|
| DiffBind | Differential H3K27ac peaks | Accepts MACS2 peaks + BAMs, uses DESeq2/edgeR internally |
| MAnorm2 | Alternative to DiffBind | Better when groups have different within-group variability |
| dcHiC/HiC-DC+ | Differential Hi-C loops | Outputs differential loops with fold-change and significance |

### Super-Enhancer Identification

| Tool | Method |
|------|--------|
| ROSE | Stitches H3K27ac peaks within 12.5 kb (excluding TSS ±2.5 kb), ranks by signal |
| HOMER findPeaks | Same ROSE methodology, use `-style super` flag |

---

## Technical Constraints

- ABC predictions depend heavily on accurate TSS annotations; errors >5 kb degrade prediction quality
- Model assumes independent, additive enhancer contributions
- Only considers elements within 5 Mb on the same chromosome
- H3K27ac HiChIP was attempted but not working; Hi-C + H3K27ac peak overlap is the alternative approach

---

## Project-Specific Notes

From meeting notes, two distinct goals:
1. **1st goal:** Link change in E-P contacts to differentially expressed genes
2. **2nd goal:** Tie change in delta contacts to ubiquitinated histone (H2AK119ub)

H2AK119ub ChIP-seq data IS available as bigWigs in `heatmaps/` (4 Ctrl + 4 Mut replicates).

Meeting hypothesis: "ub is buffer to stop k27ac contact" — once ubiquitination threshold is reached, contacts form.

---

## Companion Resources

- Pre-computed ABC predictions for 131 cell types: https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/
- BENGI benchmark datasets for validation: ChIA-PET, Hi-C, eQTL gold standards
- ABC-Max pipeline for variant prioritization: https://github.com/EngreitzLab/ABC-Max-pipeline
