# PRD: ABC Model Enhancer-Gene Linkage Analysis (v2)

## Purpose

Use the Activity-By-Contact (ABC) model to link enhancers to genes, then demonstrate that enhancer-gene linkages are dysregulated in BAP1-KO conditions and correlate with differentially expressed genes.

The ABC model formula:

**ABC score = (Enhancer Activity × Contact Frequency) / Σ(Activity × Contact) over all candidate elements within 5 Mb**

Where:
- Activity = H3K27ac signal (and/or ATAC-seq signal) at enhancer
- Contact = Hi-C contact frequency between enhancer and promoter
- Denominator normalizes per gene so scores represent relative contribution

**Reference paper:** Fulco et al. 2019 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC6886585/>

---

## Organism & Genome

- **Organism:** Mouse (Mus musculus)
- **Genome build:** mm10
- **Chromosome sizes file:** `/expanse/lustre/projects/csd940/zalibhai/taiji-new/reference/mm10/mm10.chrom.sizes`

---

## Tool Selection

**Broad Institute ABC Implementation**
- Repository: <https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction>
- Version: v1.1.2 (March 2024)
- Documentation: <https://abc-enhancer-gene-prediction.readthedocs.io/>
- Installation location: `/expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction/`

### Installation

Conda/mamba environment creation (~15 min with mamba). Dependencies include Python 3.9+, samtools, bedtools, MACS2, juicertools for Hi-C processing. Snakemake workflow supports parallel execution via `-j` flag.

### Tool Input Requirements

| Input | Format | Required | Status |
|-------|--------|----------|--------|
| Chromatin accessibility | ATAC-seq tagAlign | Yes | ✅ Ready |
| H3K27ac ChIP-seq | BAM (preferred) or bigWig | Yes (or ATAC-only mode) | ⚠️ bigWig only — see notes |
| Hi-C contact data | `.hic` format | Optional (powerlaw substitutes) | ✅ Ready |
| Gene annotations | BED6 + Ensembl ID | Yes | ❌ Needs generation |
| Chromosome sizes | TSV (chr, size) | Yes | ✅ Ready |

### Tool Outputs

- `EnhancerPredictions_threshold_*.tsv` containing:
  - Enhancer-gene pairs
  - ABC scores
  - Activity components
  - Contact frequencies
  - Distances
- Recommended ABC score threshold: **0.015** (achieves ~70% recall); project will use **≥ 0.02**
- Expected output: ~3 distal enhancers per expressed gene

---

## HPC Environment

- **Cluster:** SDSC Expanse
- **Scheduler:** SLURM
- **Project directory:** `/expanse/lustre/projects/csd940/zalibhai/abc/`
- **Account:** `csd940`
- **Partition:** `shared`

### Standard SBATCH Parameters

All job scripts must use these parameters (adjust resources as needed, do not omit any field):

```bash
#!/bin/bash
#SBATCH --job-name=abc_<step_name>
#SBATCH --output=logs/abc_<step_name>_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=40:00:00
#SBATCH --account=csd940
```

The user will manually adjust resource requests if they exceed allocation. Always include:

```bash
mkdir -p /expanse/lustre/projects/csd940/zalibhai/abc/logs
source ~/.bashrc
conda activate abc-env
```

---

## Input Data Inventory

### 1. ATAC-seq (Accessibility — Candidate Enhancers)

**TagAlign files (ABC input format, merged across replicates per condition):**
- `/expanse/lustre/projects/csd940/zalibhai/abc/input/ctrl_ATAC.tagAlign.gz`
- `/expanse/lustre/projects/csd940/zalibhai/abc/input/mut_ATAC.tagAlign.gz`

Generated from merged BAMs via:
```bash
bedtools bamtobed -i ctrl_merged_ATAC.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print}' | gzip > ctrl_ATAC.tagAlign.gz
bedtools bamtobed -i mut_merged_ATAC.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print}' | gzip > mut_ATAC.tagAlign.gz
```

Original BAMs were merged from **3 replicates per condition** (Batch 1).

**Raw replicate narrowPeak files (10 total: 3+3 Batch 1, 2+2 Batch 2):**
Location: `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/raw/`

Batch 1 (250310) — primary, n=3 per genotype:
- `250310_ATAC_N702_N501_BAP1_MATH1_Ctr1_P100_S1_L006_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250310_ATAC_N702_N503_BAP1_MATH1_Ctr2_P100_S3_L006_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250310_ATAC_N703_N501_BAP1_MATH1_Ctr3_P100_S5_L006_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250310_ATAC_N702_N502_BAP1_MATH1_mut1_P100_S2_L006_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250310_ATAC_N702_N504_BAP1_MATH1_Mut2_P100_S4_L006_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250310_ATAC_N703_N502_BAP1_MATH1_Mut3_P100_S6_L006_aligned_reads_peaks.narrowPeak_subtracted.bed`

Batch 2 (250731) — secondary, n=2 per genotype:
- `250731_ATAC__N701_N501_BAP1_MATH1_Ctr1_3037_S41_L007_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250731_ATAC__N701_N503_BAP1_MATH1_Ctr2_3093_S43_L007_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250731_ATAC__N701_N502_BAP1_MATH1_mut1_3036_S42_L007_aligned_reads_peaks.narrowPeak_subtracted.bed`
- `250731_ATAC__N701_N504_BAP1_MATH1_Mut2_3092_S44_L007_aligned_reads_peaks.narrowPeak_subtracted.bed`

**Consensus peak files (pre-generated):**
Location: `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/atac_seq/`
- `consensus_control.bed` — Batch 1, ≥2/3 replicates
- `consensus_mutant.bed` — Batch 1, ≥2/3 replicates
- `consensus_all.bed` — **Union of control + mutant consensus (PRIMARY — use for ABC)**
- `consensus_control_all.bed` — All samples, ≥2/5 replicates
- `consensus_mutant_all.bed` — All samples, ≥2/5 replicates
- `consensus_all_combined.bed` — Union of all-sample consensus
- `ATAC_up.bed` — Differential peaks (up in mutant)
- `ATAC_down.bed` — Differential peaks (down in mutant)

The consensus generation script is at `peaks/atac_seq/generate_consensus.sh` (attached separately). It uses `bedtools multiIntersectBed` with threshold ≥2 replicates, then merges across genotypes.

**For ABC, use `consensus_all.bed` as the consensus enhancer universe.** This is the union of Batch 1 control and mutant consensus peaks — exactly what the PRD requires (intersect within genotype, merge across genotypes).

### 2. H3K27ac ChIP-seq (Enhancer Activity)

**bigWig files (merged across replicates per condition):**
- `/expanse/lustre/projects/csd940/zalibhai/abc/input/H3K27ac_ctrl_merged.bw`
- `/expanse/lustre/projects/csd940/zalibhai/abc/input/H3K27ac_mut_merged.bw`

**⚠️ IMPORTANT:** BAM files are NOT available for H3K27ac. Only bigWigs exist. The standard ABC Snakemake pipeline expects BAM input for H3K27ac (it computes signal internally from BAMs). This requires one of:

1. **ATAC-only mode:** Run ABC using only ATAC-seq for both activity and accessibility. H3K27ac is technically optional in the ABC pipeline. Activity is then defined purely from ATAC signal. This is the simplest path.
2. **Custom quantification:** Use `bigWigAverageOverBed` or deepTools to quantify H3K27ac signal at candidate enhancers outside the pipeline, then feed pre-computed activity values into a modified ABC scoring step.
3. **Obtain BAMs:** Check if H3K27ac BAMs can be sourced from the original processing.

**Decision needed:** Determine whether to run ATAC-only mode or pursue a workaround for bigWig input. The ABC pipeline does support running without H3K27ac.

**H3K27ac peak BED files (for downstream cross-referencing, NOT direct ABC input):**
- `peaks/beds/H3K27acCerebellumLate2.bed` — **Late timepoint (adult) = PRIMARY**
- `peaks/beds/H3K27acCerebellumEarly2.bed` — Early timepoint (young), included for reference only

### 3. Hi-C Contact Data

**Per-condition `.hic` files at 5kb resolution:**
- WT: `/expanse/lustre/projects/csd940/zalibhai/nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/resorted_ctrl.hic`
- KO: `/expanse/lustre/projects/csd940/zalibhai/nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/resorted_mut.hic`

Available resolutions: 5kb, 10kb, 25kb. **ABC requires 5kb.**

These are the raw Hi-C contact matrices needed by ABC to query contact frequency between any enhancer-promoter pair — distinct from the pre-called differential loops.

**Pre-called differential loops (for downstream cross-referencing, NOT direct ABC input):**
- `characterized_loops.tsv`
  - Contains: `logFC`, `anchor1_H3K27ac_overlap`, `anchor2_H3K27ac_overlap`, `anchor1_type`, `loop_type`
  - Gene annotations: `anchor1_nearest_gene`, `anchor2_nearest_gene`, `anchor1_distance_to_tss`
  - Generated by a custom pipeline (not relevant to ABC setup)
  - Will be used in downstream integration after ABC predictions are generated

### 4. RNA-seq Differential Expression

**Primary (adult timepoint):**
- `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
- Standard DESeq2 output

| Column | Content | Example |
|--------|---------|---------|
| `ensembl_gene_id` | Gene symbols (NOT Ensembl IDs despite column name) | `Xkr4` |
| `log2FoldChange` | log2 fold change (KO vs WT) | `0.5737` |
| `padj` | Adjusted p-value (BH) | `0.00523` |

Also contains: `baseMean`, `lfcSE`, `stat`, `pvalue`, per-sample normalized counts (`ctrl2-10`, `mut2-10`).

**Secondary (young timepoint — NOT used in this analysis):**
- `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx`

### 5. Gene Annotations

**Status: ❌ Needs generation/reformatting for ABC**

Existing GENCODE vM25 (mm10) annotations:
- `gencode.vM25.mouse.genes.annotation.bed.gz` — Gene bodies (e.g., `1  4807787  4848409  Gene  Lypla1`)
- `gencode.vM25.mouse.tss_region.annotation.bed.gz` — TSS ±200bp (e.g., `1  4807587  4807987  TSS  Lypla1`)

**Two issues to fix before use with ABC:**
1. **No `chr` prefix** — These use `1, 2, 3...` but Hi-C/ATAC data uses `chr1, chr2, chr3...`. Requires `sed 's/^/chr/'` or similar.
2. **ABC expects specific format** — ABC gene annotation needs: BED6 (chr, start, end, gene, score, strand) + Ensembl ID as 7th column. Current files lack strand info.

**Recommended approach:** Download GENCODE vM25 GTF and extract a properly formatted BED6 + Ensembl ID file. The ABC repo may have a helper script, or this is straightforward from the GTF:
```
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
```

The gene annotation must also be used to generate the `genome_tss` file (500bp TSS regions with Ensembl ID), which ABC uses as the promoter reference for computing E-P contact.

### 6. H2AK119ub ChIP-seq (Deferred — 2nd Goal)

**Available as bigWig files in `heatmaps/`:**

Ctrl (WT) — 4 usable replicates:
- `H2AK119ubCtrl1.bw`, `H2AK119ubCtrl2.bw`, `H2AK119ubCtrl3.bw`, `H2AK119ubCtrl4.bw`

Mut (KO) — 4 usable replicates:
- `H2AK119ubMut1.bw`, `H2AK119ubMut2.bw`, `H2AK119ubMut3.bw`, `H2AK119ubMut4.bw`

Excluded (bad agreement with prior data):
- `H2AK119ubCtrl5-badagreementwithpriordata.bw`, `H2AK119ubCtrl6-badagreementwithpriordata.bw`
- `H2AK119ubMut5-badagreementwithpriordata.bw`, `H2AK119ubMut6-badagreementwithpriordata.bw`

**This data is deferred until after ABC results from Goal 1 are interpreted.**

---

## Analysis Scope

### Timepoint
- **Primary:** Adult (late) timepoint only
- **Mention:** Early timepoint exists but is not analyzed here

### Two Distinct Goals
1. **1st goal (this analysis):** Link change in E-P contacts to differentially expressed genes via ABC model
2. **2nd goal (deferred):** Correlate delta E-P contacts with H2AK119ub levels — test hypothesis that "ub is buffer to stop k27ac contact"

---

## Critical Design Decision: Consensus Enhancer Universe

**IMPORTANT:** Many papers use ABC to define cell-type specific enhancers by calling separate enhancer sets per condition. **This is NOT what we want.** Defining separate enhancers between ctrl and mut would confound differential analysis.

**Correct approach:**
1. Use a **consensus enhancer universe** shared across both conditions
2. The consensus is already generated: **`consensus_all.bed`** = union of Batch 1 control and mutant consensus ATAC peaks (intersected at ≥2/3 replicates within genotype, then merged across genotypes)
3. Quantify enhancer activity (H3K27ac and/or ATAC signal) at this fixed set per condition
4. Run ABC with this fixed enhancer set for both WT and KO
5. Compare ABC scores across conditions for the same enhancer-gene pairs

---

## Workflow

### Step 0: Prepare Gene Annotations for mm10

ABC ships with hg38 annotations. For mm10, generate:
1. Download GENCODE vM25 GTF
2. Extract gene annotation BED: `chr  start  end  gene_symbol  score  strand  ensembl_id`
3. Extract TSS file: 500bp regions centered on TSS, same column format
4. Ensure `chr` prefix on all chromosomes
5. Place in `/expanse/lustre/projects/csd940/zalibhai/abc/reference/`

### Step 1: Prepare ABC Configuration

Create `config-biosamples.tsv` with two rows (one per condition):

| biosample | DHS | ATAC | H3K27ac | default_accessibility_feature | HiC_file | HiC_type | HiC_resolution |
|-----------|-----|------|---------|-------------------------------|----------|----------|----------------|
| WT_cerebellum | | /path/to/ctrl_ATAC.tagAlign.gz | (see H3K27ac note) | ATAC | /path/to/resorted_ctrl.hic | hic | 5000 |
| KO_cerebellum | | /path/to/mut_ATAC.tagAlign.gz | (see H3K27ac note) | ATAC | /path/to/resorted_mut.hic | hic | 5000 |

**H3K27ac column:** If BAMs are not available and ATAC-only mode is used, leave this column blank. If a workaround for bigWig input is implemented, populate accordingly.

Update `config/config.yaml`:
- `biosamplesTable`: point to the new TSV
- `chrom_sizes`: `/expanse/lustre/projects/csd940/zalibhai/taiji-new/reference/mm10/mm10.chrom.sizes`
- `genes`: path to generated mm10 gene BED
- `genome_tss`: path to generated mm10 TSS BED
- `regions_blocklist`: mm10 blocklist (ENCODE, needs download if not present)
- Adjust `params_macs` genome size for mouse (`mm`)

### Step 2: ATAC-seq Index Requirement

ABC requires a tabix `.tbi` index for tagAlign files. Verify these exist:
```bash
ls -la /expanse/lustre/projects/csd940/zalibhai/abc/input/ctrl_ATAC.tagAlign.gz.tbi
ls -la /expanse/lustre/projects/csd940/zalibhai/abc/input/mut_ATAC.tagAlign.gz.tbi
```
If missing, generate with: `tabix -p bed ctrl_ATAC.tagAlign.gz`
(File must be bgzipped, not just gzipped — if generated with `gzip`, will need `zcat | bgzip` first.)

### Step 3: Run ABC Per Condition

Run the Snakemake pipeline separately for WT and KO:
```bash
conda activate abc-env
snakemake -j16  # uses 16 cores from SLURM allocation
```

Or submit as SLURM job with the standard parameters above.

Both conditions must use the same consensus enhancer set. If ABC's internal peak calling selects different enhancers per run, the consensus ATAC peaks from `consensus_all.bed` should be provided as pre-defined candidate regions (check ABC docs for how to supply pre-called peaks instead of letting MACS2 call them internally).

**Runtime/compute note:** If the candidate enhancer set is too large (>150k), cap to top 100-150k by H3K27ac activity to manage compute. This is controlled in ABC's configuration.

### Step 4: Filter ABC Predictions

Per condition, filter to links with **ABC score ≥ 0.02**.

Output files:
- `results/WT_cerebellum/Predictions/EnhancerPredictions_threshold_0.02.tsv`
- `results/KO_cerebellum/Predictions/EnhancerPredictions_threshold_0.02.tsv`

### Step 5: Compute Delta ABC

For each enhancer-gene pair present in either condition:
```
ΔABC = ABC_KO − ABC_WT
```

This requires joining WT and KO prediction tables on (enhancer coordinates, gene). Pairs present in only one condition get the other condition's score set to 0.

### Step 6: Integrate with RNA-seq

1. Read `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
2. Join ΔABC with gene `log2FoldChange` on gene symbol (`ensembl_gene_id` column, which actually contains gene symbols)
3. Summarize directional concordance:
   - ΔABC ↑ + log2FC ↑ (gained enhancer contact → upregulated gene)
   - ΔABC ↓ + log2FC ↓ (lost enhancer contact → downregulated gene)
   - Discordant pairs

### Step 7: Aggregate to Gene Level

For each gene, compute:
- **Strongest link:** max |ΔABC| across all linked enhancers
- **Total ABC per gene:** Σ ABC scores (per condition)
- **Number of perturbed links:** count of enhancers with |ΔABC| > threshold

### Step 8: Downstream Cross-Referencing

After ABC predictions are generated, cross-reference with:
- `characterized_loops.tsv` — do ABC-predicted E-P pairs overlap with called differential loops?
- H3K27ac peak BEDs — which enhancers overlap differential H3K27ac peaks?
- Super-enhancer identification (Task 7e) — loops where both anchors are H3K27ac+

---

## Tasks Mapping

| Task | Description | Dependencies |
|------|-------------|-------------|
| **7a** | Run ABC model per condition → enhancer-gene links with scores | Steps 0-4 |
| **7b** | Compute ΔABC per enhancer-gene pair | Step 5 |
| **7c** | Correlate ΔABC with RNA-seq log2FC (scatter plot) | Step 6 |
| **7d** | Identify dysregulated E-G pairs (both |ΔABC| and |log2FC| significant) | Steps 6-7 |
| **7e** | Super-enhancer to gene linkages (both anchors H3K27ac+) | Step 8 + H3K27ac peaks |
| **7f** | Correlate with H2AK119ub (DEFERRED) | After 7a-7e interpreted |

---

## Open Questions / Blockers

### 1. H3K27ac BAM Availability (Medium Priority)

Only bigWig files exist for H3K27ac. Options in order of preference:
1. **Run ATAC-only mode** — simplest, well-supported by ABC. Activity = ATAC signal only. H3K27ac peaks can still be used for downstream annotation.
2. **Custom activity quantification** — use `bigWigAverageOverBed` on H3K27ac bigWigs to quantify signal at consensus enhancers, then integrate into ABC scoring outside the Snakemake pipeline.
3. **Obtain BAMs** — check if original H3K27ac BAMs exist somewhere.

### 2. TagAlign Index Files

The tagAlign files were generated with `gzip` (not `bgzip`). ABC requires tabix-indexed tagAlign, which requires bgzip compression. If `.tbi` index files don't exist alongside the tagAlign files, will need:
```bash
zcat ctrl_ATAC.tagAlign.gz | bgzip > ctrl_ATAC.tagAlign.bgz
tabix -p bed ctrl_ATAC.tagAlign.bgz
```

### 3. Consensus Peaks as ABC Input

Need to verify how to supply pre-defined candidate enhancer regions (from `consensus_all.bed`) to the ABC pipeline instead of letting it call peaks internally. The standard pipeline runs MACS2 on the ATAC data — we want to bypass this and use our pre-generated consensus set to ensure both conditions use identical candidate elements.

### 4. mm10 Blocklist

ABC requires a regions blocklist. Need to obtain the ENCODE mm10 blocklist:
```
https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
```

### 5. Ubiquitous Genes List

ABC flags ubiquitously expressed genes (they typically lack distal enhancers). Need an mm10-appropriate list. The ABC repo may include one for human — will need a mouse equivalent, or this parameter can potentially be left empty.

---

## File Path Summary

```
/expanse/lustre/projects/csd940/zalibhai/
├── abc/
│   ├── ABC-Enhancer-Gene-Prediction/    # Cloned repo
│   ├── input/
│   │   ├── ctrl_ATAC.tagAlign.gz        # WT ATAC (merged, 3 reps)
│   │   ├── mut_ATAC.tagAlign.gz         # KO ATAC (merged, 3 reps)
│   │   ├── H3K27ac_ctrl_merged.bw       # WT H3K27ac (bigWig only)
│   │   └── H3K27ac_mut_merged.bw        # KO H3K27ac (bigWig only)
│   ├── reference/                        # To be populated
│   │   ├── mm10_genes.bed               # Generate from GENCODE vM25
│   │   ├── mm10_tss.bed                 # Generate from GENCODE vM25
│   │   └── mm10-blacklist.v2.bed        # Download from ENCODE
│   └── logs/
├── mariner_hi-c/
│   ├── peaks/
│   │   ├── atac_seq/
│   │   │   ├── consensus_all.bed        # ← CONSENSUS ENHANCER UNIVERSE
│   │   │   ├── consensus_control.bed
│   │   │   ├── consensus_mutant.bed
│   │   │   ├── ATAC_up.bed
│   │   │   ├── ATAC_down.bed
│   │   │   └── raw/                     # 10 narrowPeak files
│   │   └── beds/
│   │       ├── H3K27acCerebellumLate2.bed   # Adult H3K27ac peaks
│   │       └── H3K27acCerebellumEarly2.bed  # Young H3K27ac peaks
│   ├── characterized_loops.tsv          # Differential loops (downstream)
│   └── tads/
│       ├── adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx
│       └── young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx
├── nf-hic/250402_Bap1_deepseq/juicerpre/merged/hic/
│   ├── resorted_ctrl.hic                # WT Hi-C (5kb resolution)
│   └── resorted_mut.hic                 # KO Hi-C (5kb resolution)
├── taiji-new/reference/mm10/
│   └── mm10.chrom.sizes
└── heatmaps/                            # H2AK119ub bigWigs (deferred)
    ├── H2AK119ubCtrl[1-4].bw
    └── H2AK119ubMut[1-4].bw
```

---

## Technical Constraints

- ABC predictions depend heavily on accurate TSS annotations; errors >5 kb degrade prediction quality
- Model assumes independent, additive enhancer contributions
- Only considers elements within 5 Mb on the same chromosome — no interchromosomal predictions
- H3K27ac HiChIP was attempted but is not working; Hi-C + H3K27ac peak overlap is the alternative approach
- ABC's Snakemake pipeline is designed for hg38 — config changes needed for mm10 (genome size param in MACS2, reference files)
- Elements are defined as ~500bp regions centered on ATAC peaks

---

## Companion Resources

- Pre-computed ABC predictions for 131 cell types: <https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/>
- BENGI benchmark datasets: ChIA-PET, Hi-C, eQTL gold standards
- ABC-Max pipeline for variant prioritization: <https://github.com/EngreitzLab/ABC-Max-pipeline>
- ENCODE mm10 blocklist: <https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz>
- GENCODE vM25 GTF: <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz>
