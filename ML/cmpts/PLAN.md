# Plan: CALDER2 + SNIPER Subcompartment Calling Pipeline for BAP1-KO Cerebellum

**Companion documents:** `ML/cmpts/SNIPER/README.md` (SNIPER tool), `ML/cmpts/CALDER2/README.md` (CALDER2 tool), `loops/docs/dixon-meeting-summary.md` (biological context), `tads/CLAUDE.md` (existing compartment pipeline).

## Context

The Dixon meeting (2026-04-10) identified subcompartment calling as a key analysis — Jesse recommended SNIPER over DCIC (which gave iffy results) to look for B1↔B2 or A1↔A2 switches rather than full A↔B flips. Our HOMER-based A/B compartment analysis at 25kb found 8,189 significant differential bins (2,704 A→B, 5,485 B→A), but Jesse's observation was that these may be quantitative weakening within subcompartments rather than true flips. Subcompartment resolution would reveal whether "B→A" transitions are actually B2→B1 (within B) or B1→A2 (genuine).

**Critical finding during planning:** SNIPER's pre-trained models are hardcoded for human hg19 (chr1-22) and cannot be applied to mouse mm10 (chr1-19). The neural network input dimensions (13,393/13,594 bins) are baked into the .h5 weights. Five hard blockers: hardcoded chromosome ranges, mismatched model dimensions, hg19-specific crop maps, hg19-specific coordinate maps, and hg19-default BED output. Retraining from scratch on mm10 is required, which needs ground-truth subcompartment labels.

**Solution:** CALDER2 as primary tool (R package, works natively on mm10, no training needed), SNIPER retrained on mm10 as secondary validation (using CALDER2 labels as ground truth). Both timepoints (250402=late/adult, 250831=early/P12). Integration with existing HOMER A/B compartment calls.

**Scope:** Both timepoints, merged samples only (ctrl_merged, mut_merged) for primary analysis. 100kb resolution (subcompartment standard). Autosomes only (chr1-19).

```
                          PIPELINE DEPENDENCY GRAPH

  Track A: CALDER2 (Primary) ─────────────────────────────────────────────────
                                                                              |
  A0: Env setup (interactive, once)                                           |
       │                                                                      |
       v                                                                      |
  A1: Run CALDER2 ──── 4 jobs ──── {250402,250831} × {ctrl,mut}_merged       |
       │                                                                      |
       v                                                                      |
  A2: Differential subcompartments ── 2 jobs (1 per tp) ──────────────────────|
       │                                                                      |
       ├──────────────┬──────────────┐                                        |
       v              v              v                                        |
  A3: Epigenomic   A4: HOMER     B2: Generate                                |
   validation      integration   SNIPER labels                               |
                                     │                                        |
  Track B: SNIPER (Secondary) ───────┼────────────────────────────────────────|
                                     │                                        |
  B0: Python 3.7 + TF 1.15 env      │                                        |
  B1: Adapt SNIPER for mm10 ────────>│                                        |
       │                             v                                        |
       └────────────────────> B3: Train SNIPER (2 jobs, 1 per tp)             |
                                     │                                        |
                                     v                                        |
                              B4: Apply SNIPER ── 4 jobs ── merged × tp       |
                                     │                                        |
                                     v                                        |
                              B5: Concordance (SNIPER vs CALDER2)             |
                                                                              |
  Track C: Integration (after A + B) ─────────────────────────────────────────|
   C1: Combined comparison (subcompartment tracks, transition matrices)       |
   C2: Epigenomic enrichment heatmap (Fig 2c equivalent)                      |
   C3: Loop/stripe × subcompartment integration                               |
   C4: HOMER A/B × subcompartment decomposition                              |
```

---

## Track A: CALDER2 (Primary Subcompartment Calling)

### A0 — Environment Setup — PENDING

**Script:** `scripts/A0_setup_calder2_env.sh` (interactive, run once on Expanse login node)

CALDER2 (`library(CALDER)`) is an R package with C++ compiled code (RcppArmadillo + OpenMP) and Bioconductor dependencies. It accepts .hic files directly via `strawr`, ships an mm10 reference BED (`inst/extdata/mm10_all_sub_compartments.bed`, 15,364 bins), and is fully mm10-aware via `genome='mm10'`.

**Dependencies (from DESCRIPTION):**
- R ≥ 3.5.2
- Bioconductor: `GenomicRanges`, `rhdf5`
- CRAN: `strawr` (≥0.0.9), `data.table`, `ape`, `dendextend`, `fitdistrplus`, `igraph`, `Matrix`, `rARPACK`, `factoextra`, `fields`, `ggplot2`, `optparse`, `R.utils`, `doParallel`
- LinkingTo: `Rcpp`, `RcppArmadillo`
- System: C++11 compiler with OpenMP, BLAS/LAPACK

**Install approach:**
```bash
conda create -n calder2_env r-base=4.2 r-rcpp r-rcpparmadillo -c conda-forge -y
conda activate calder2_env
Rscript -e '
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install(c("GenomicRanges", "rhdf5"))
  install.packages(c("strawr", "data.table", "ape", "dendextend", "fitdistrplus",
                      "igraph", "Matrix", "rARPACK", "factoextra", "fields",
                      "ggplot2", "optparse", "R.utils", "doParallel"))
  install.packages("/path/to/ML/cmpts/CALDER2", repos=NULL, type="source")
'
```

Alternatively: `conda install -c bioconda r-calder2` (if available).

**Verification:** `Rscript -e "library(CALDER); cat('CALDER loaded OK\n')"` exits 0.

### A1 — Run CALDER2 on Each Sample — PENDING

**Script:** `scripts/A1_run_calder2.R` + `scripts/A1_run_calder2.sb` + `scripts/A1_submit_calder2.sh`

4 jobs: {250402, 250831} × {ctrl_merged, mut_merged}.

**CALDER2 function call:**
```r
library(CALDER)

CALDER(
  contact_file_hic    = hic_path,            # .hic file, strawr extracts KR-normalized contacts
  chrs                = as.character(1:19),   # autosomes only, no "chr" prefix (CALDER strips it)
  bin_size            = 50000,               # 50kb — CALDER2 auto-extends to {50kb, 100kb}
  save_dir            = out_dir,
  save_intermediate_data = TRUE,             # needed for sub_domains and debugging
  genome              = "mm10",              # activates mm10 reference for A/B polarity
  n_cores             = 8,                   # parallel per-chromosome processing
  sub_domains         = FALSE,               # skip nested boundaries (not needed for subcompartments)
  single_binsize_only = FALSE                # use multi-resolution optimization (50kb + 100kb)
)
```

**Resolution choice:** `bin_size=50000` (50kb). CALDER2 automatically extends to `{50000, 100000}` when `genome='mm10'` and `single_binsize_only=FALSE`. The output `all_sub_compartments.tsv` contains bins at the optimal resolution per region. For uniform 100kb output (needed for SNIPER integration), we'll aggregate in A2.

**Key outputs per sample:**
```
DATA_DIR/calder2/{tp}/{sample}/
  sub_compartments/
    all_sub_compartments.bed     # BED9, hierarchical labels (e.g., A.1.1, B.2.2.1.2...)
    all_sub_compartments.tsv     # 6-col: chr, pos_start, pos_end, comp_name, comp_rank, continuous_rank
    cor_with_ref.txt             # per-chr correlation with mm10 reference
  intermediate_data/             # per-chr per-resolution .Rds files
```

**CALDER2 label hierarchy:**
- Depth 1: A, B (2 classes — standard compartments)
- Depth 2: A.1, A.2, B.1, B.2 (4 classes — primary subcompartments)
- Depth 3: A.1.1, A.1.2, A.2.1, A.2.2, B.1.1, B.1.2, B.2.1, B.2.2 (8 classes)
- Deeper: up to ~10 levels for fine-grained domains

For SNIPER compatibility (5 classes: A1/A2/B1/B2/B3), we use depth 2 (4 classes). The mm10 reference has no B.3 at depth 2 (only A.1=4,498, A.2=3,692, B.1=3,134, B.2=4,040 bins). B3 (a specialized human subcompartment on chr19) is not expected in mouse. We'll use 4-class labeling throughout.

**SLURM:** `--cpus-per-task=8 --mem=64G --time=12:00:00 --account=csd940 --partition=shared`

**Verification:**
- All 4 `all_sub_compartments.tsv` files exist with rows for all 19 autosomes
- Per-sample label distribution at depth 2: A.1 ~25-35%, A.2 ~20-25%, B.1 ~15-25%, B.2 ~20-30% (rough expectation from mm10 reference: A.1=29%, A.2=24%, B.1=20%, B.2=26%)
- `cor_with_ref.txt` shows positive correlations (>0.3) for most chromosomes

### A2 — Differential Subcompartment Analysis — PENDING

**Script:** `scripts/A2_differential_subcompartments.R` + `scripts/A2_run.sb`

Dependencies: All 4 A1 jobs complete.

**Logic:**
1. Load all 4 `all_sub_compartments.tsv` files
2. Truncate hierarchical labels to depth 2 (e.g., `A.1.1.2.1` → `A.1`) for 4-class analysis
3. Bin to uniform 100kb resolution: for variable-width CALDER2 output, assign each 100kb bin the label of the largest overlapping CALDER2 segment (plurality vote)
4. For each timepoint, join ctrl and mut labels on (chr, 100kb_bin_start)
5. Classify transitions: `ctrl_label → mut_label` for each bin
6. Compute transition matrix (4×4: A.1, A.2, B.1, B.2)
7. Chi-squared test on transition counts vs expected (no-change null)

**Outputs:**
```
ML/cmpts/sniper/outputs/calder2/
  {tp}_subcompartment_labels_100kb.tsv     # chr, start, end, ctrl_label, mut_label, continuous_rank_ctrl, continuous_rank_mut
  {tp}_transition_matrix.tsv               # 4x4 counts
  {tp}_transition_summary.tsv              # counts per transition type
  {tp}_transition_sankey/                  # multi-format figure
  {tp}_subcompartment_genome_pct/          # multi-format figure
```

**SLURM:** `--cpus-per-task=4 --mem=16G --time=02:00:00`

**Verification:**
- Transition matrices sum to total number of 100kb autosomal bins (~25,600 for mm10)
- Chi-squared p-value < 0.05 (late timepoint expected to show significant transitions)
- Key question answered: what fraction of bins change subcompartment? If <5%, the A/B system is largely stable; if >10%, there are substantial sub-compartment rearrangements

### A3 — Epigenomic Validation — PENDING

**Script:** `scripts/A3_epigenomic_validation.R` + `scripts/A3_run.sb`

Replicates SNIPER paper Figure 2c: fold-enrichment heatmap of epigenomic marks per subcompartment.

**Available marks for validation (10 marks × ctrl/mut = 20 BigWigs):**

| Mark | BigWig | Subcompartment expectation |
|------|--------|--------------------------|
| H3K27ac | `{H3K27acCtrl,H3K27acMut}.bw` | A.1 >> A.2 > B.1 > B.2 |
| H3K4me3 | `{H3K4me3Ctrl,H3K4me3Mut}.bw` | A.1 >> A.2 > B.1 > B.2 |
| H3K27me3 | `{H3K27me3Ctrl,H3K27me3Mut}.bw` | B.1 >> B.2 > A.2 > A.1 |
| H2AK119ub | `{H2AK119ubCtrl,H2AK119ubMut}.bw` | B.1 enriched (Polycomb) |
| ATAC-seq | `{ATACctrl,ATACmut}.bw` | A.1 >> A.2 > B.1 > B.2 |
| RNA-seq | `{RNActrl,RNAmut}.bw` | A.1 > A.2 >> B.1 ≥ B.2 |
| DNA methylation | `{DNAmethylationCtrl,DNAmethylationMut}.bw` | Complex |
| H3K27me1 | `{H3K27me1Ctrl,H3K27me1Mut}.bw` | Not well characterized |

Missing: H3K36me3 BigWig (peaks exist but no signal track), H3K79me2, H4K20me1, LMNB1. B.2 vs B.3 validation is limited without these.

BigWig path on HPC: `/expanse/lustre/projects/csd940/zalibhai/bigwigs/`

**Logic:**
1. Load `{tp}_subcompartment_labels_100kb.tsv` for ctrl samples
2. For each BigWig, use `rtracklayer::import.bw(bw, which=GRanges(bins))` to extract mean signal per 100kb bin
3. Group by depth-2 label, compute median signal per mark per subcompartment
4. Compute fold-enrichment: `median_in_subcomp / median_across_all_bins`
5. Plot as heatmap (rows=marks, cols=subcompartments A.1→B.2, color=RdBu fold-enrichment)

**SLURM:** `--cpus-per-task=8 --mem=32G --time=04:00:00`

**Verification:**
- H3K27ac gradient monotonically decreasing A.1 → B.2 in ctrl
- H3K27me3 enriched in B.1 (>1.5× fold) — this is the key Polycomb compartment
- If gradients are flat or inverted, something is wrong with compartment polarity

### A4 — Integration with HOMER A/B Compartments — PENDING

**Script:** `scripts/A4_integrate_homer_compartments.R` + `scripts/A4_run.sb`

**Input:**
- HOMER: `tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv` (104,071 × 25kb bins)
- CALDER2: `{tp}_subcompartment_labels_100kb.tsv` (from A2)

**Logic:**
1. Aggregate HOMER 25kb bins to 100kb: mean PC1 per 100kb window (4 bins per window)
2. Join HOMER 100kb bins with CALDER2 labels via `GenomicRanges::findOverlaps()`
3. For the 8,189 significant HOMER differential bins (FDR<0.05, |Diff|>0.30):
   - Extract ctrl_label and mut_label from CALDER2
   - Cross-tabulate: HOMER direction × CALDER2 transition type
4. Key decomposition:
   - 2,704 HOMER A→B bins: how many are A.1→B.1 vs A.2→B.1 vs A.1→B.2?
   - 5,485 HOMER B→A bins: how many are B.1→A.2 vs B.2→A.2 vs B.2→A.1?
   - "Weakening" bins (same A or same B, but significant PC1 change): A.1→A.2 or B.1→B.2?

**Outputs:**
```
ML/cmpts/sniper/outputs/integration/
  homer_calder2_crosstab.tsv
  homer_weakening_decomposition.tsv
  homer_calder2_sankey/            # multi-format figure
  homer_calder2_dotplot/           # multi-format figure
```

**SLURM:** `--cpus-per-task=4 --mem=16G --time=02:00:00`

---

## Track B: SNIPER (Secondary Validation)

### B0 — Environment Setup — PENDING

**Script:** `scripts/B0_setup_sniper_env.sh` (interactive, run once)

SNIPER requires legacy Python + TensorFlow. Using TF 1.15 instead of 1.12 (backward-compatible, supports Python 3.7, last TF 1.x release).

```bash
conda create -n sniper_env python=3.7 -c conda-forge -y
conda activate sniper_env
pip install tensorflow==1.15.5 keras==2.2.4 numpy==1.15.4 scipy==1.1.0 \
            pandas==0.24.2 h5py==2.8.0 PyYAML==4.2b2
```

**Note:** TF 1.15 GPU requires CUDA 10.0 (Expanse has CUDA 11+). Use CPU-only TF 1.15. Training on CPU is feasible — the mm10 inter-chromosomal matrix is ~12.8k × 11.8k entries, and the autoencoder has ~15M parameters. Training takes 12-20h on CPU.

**juicer_tools:** SNIPER calls `java -jar {juicer_tools} dump observed KR {hic} {chr1} {chr2} BP 100000 {output}`. Need juicer_tools.jar on HPC. Check existing paths: `lab/pipeline-scripts/juicer.sb` uses Singularity container (`juicer_2.0.1.sif`), but SNIPER needs the JAR directly. Verify/locate: the tads pipeline extracts via `straw` (R), but SNIPER needs the Java `dump` command.

**Verification:** `python -c "import tensorflow as tf; print(tf.__version__)"` → `1.15.5`

### B1 — Adapt SNIPER for mm10 — PENDING

**Script:** `scripts/B1_adapt_sniper_mm10.py` (one-time utility, run interactively or as SLURM job)

This creates mm10-adapted copies of SNIPER's core modules. We modify copies, not the original SNIPER source, to preserve the upstream repo.

**Files to create in `ML/cmpts/SNIPER/`:**

| File | Purpose |
|------|---------|
| `utilities/mm10_config.py` | mm10 chromosome constants: `ODD_CHROMS = [1,3,5,7,9,11,13,15,17,19]`, `EVEN_CHROMS = [2,4,6,8,10,12,14,16,18]`, `SIZES_FILE = 'data/mm10.chrom.sizes'` |
| `utilities/data_processing_mm10.py` | Patched `hicToMat()` and `predictionsToBed()` using mm10 chromosome ranges and `mm10.chrom.sizes` |
| `utilities/interchromosome_matrix_mm10.py` | Patched `construct()` using mm10 chromosome ranges |
| `sniper_train_mm10.py` | mm10 training entry point (imports mm10 variants) |
| `sniper_apply_mm10.py` | mm10 application entry point |

**Specific code changes (copies, not patches to originals):**

1. **`data_processing_mm10.py`** — `hicToMat()`:
   - Line 38: `range(1,23,2)` → `range(1,20,2)` (odd: 1,3,5,...,19)
   - Line 39: `range(2,23,2)` → `range(2,20,2)` (even: 2,4,6,...,18)
   - Lines 63-64: same change in autoremove loop
   - `predictionsToBed()` line 136: default `sizes_file='data/mm10.chrom.sizes'`
   - Line 137: `np.int` → `int` (numpy deprecation fix for numpy ≥1.20)

2. **`interchromosome_matrix_mm10.py`** — `construct()`:
   - Line 10: default `sizes_file='data/mm10.chrom.sizes'`
   - Line 15: `range(1,23,2)` → `range(1,20,2)`
   - Line 21: `range(2,23,2)` → `range(2,20,2)`

3. **Crop map generation** — `B1_adapt_sniper_mm10.py` also generates:
   - `crop_map/mm10_cropMap.mat`: `rowMap` (N_odd × 3) and `colMap` (N_even × 3) mapping each retained bin to (original_index, chrom_number, bin_position)
   - `crop_map/mm10_cropIndices.mat`: `odd_indices` and `even_indices` boolean arrays

   The crop map is generated from the ctrl_merged inter-chromosomal matrix by:
   a. Extract inter-chromosomal contacts via juicer_tools dump (KR normalized, 100kb)
   b. Build the full 12,808 × 11,831 matrix
   c. Compute per-row and per-column sums
   d. Mark bins with coverage < 1st percentile OR in ENCODE blacklist (`tads/mm10-blacklist.v2.bed`) as excluded
   e. Save retained indices as crop map

**Matrix dimensions (mm10 at 100kb):**
- Odd chromosomes (1,3,5,7,9,11,13,15,17,19): 1955+1601+1519+1455+1246+1221+1205+1041+950+615 = **12,808 bins**
- Even chromosomes (2,4,6,8,10,12,14,16,18): 1822+1566+1498+1295+1307+1202+1250+983+908 = **11,831 bins**
- Note: asymmetric (10 odd, 9 even), unlike hg19 (11 odd, 11 even)

**SLURM for crop map generation:** `--cpus-per-task=8 --mem=64G --time=04:00:00` (juicer_tools dump of 10×9=90 chromosome pairs takes ~2h)

**Verification:**
- `mm10_cropMap.mat` exists with `rowMap` shape (~11k, 3) and `colMap` shape (~10k, 3)
- `mm10_cropIndices.mat` exists with `odd_indices` and `even_indices` arrays
- Test extraction: `python -c "from utilities.data_processing_mm10 import hicToMat; print('OK')"` from SNIPER directory

### B2 — Generate Ground Truth Labels from CALDER2 — PENDING

**Script:** `scripts/B2_generate_labels_from_calder2.R` + `scripts/B2_run.sb`

Dependencies: A1 (CALDER2 output) + B1 (crop map).

**Logic:**
1. Load CALDER2 `all_sub_compartments.tsv` for ctrl_merged (one per timepoint)
2. Truncate labels to depth 2: A.1→0, A.2→1, B.1→2, B.2→3 (4-class, no B3 in mm10)
3. Bin to 100kb and assign labels (plurality vote for variable-width CALDER2 segments)
4. For each 100kb bin in the odd-chromosome order (chr1, chr3, ..., chr19), look up label
5. Apply `mm10_cropIndices` mask to keep only non-sparse bins
6. Save as `mm10_labels_{tp}.mat` with struct fields `rows` (odd labels) and `cols` (even labels)

**Mapping:** CALDER2 depth-2 → SNIPER integer labels:
| CALDER2 | SNIPER int | SNIPER name |
|---------|-----------|-------------|
| A.1 | 0 | A1 |
| A.2 | 1 | A2 |
| B.1 | 2 | B1 |
| B.2 | 3 | B2 |

SNIPER's classifier has 5 output neurons (A1/A2/B1/B2/B3). Since mm10 lacks B3, we train with 4 classes. Modify the classifier output from 5 → 4 neurons in the mm10 training script, or keep 5 and accept that B3 will never be predicted (the bootstrap resampling handles class imbalance, but an empty class would cause errors). **Decision: use 4-class classifier.**

**SLURM:** `--cpus-per-task=4 --mem=16G --time=01:00:00`

**Verification:**
- `.mat` files have `rows` and `cols` arrays matching cropped bin counts
- Label distribution matches CALDER2 genome fractions (±5%)

### B3 — Train SNIPER on Control Merged — PENDING

**Script:** `scripts/B3_train_sniper.sb` + `scripts/B3_submit_train.sh`

2 jobs: one per timepoint.

**Training protocol (from SNIPER paper):**
1. Extract inter-chromosomal contacts from ctrl_merged.hic via juicer_tools dump (KR, 100kb)
2. Build contact probability matrix: P_ij = exp(-1 / (C_ij + 1e-10))
3. Trim with mm10 crop indices
4. Train denoising autoencoder (DAE): input=sparse probabilities, target=same (self-supervised on this single sample since we don't have a high-coverage reference like GM12878)
5. Extract latent variables from encoder
6. Train MLP classifier on latent variables using CALDER2-derived labels

**Architecture (from `pipeline/models.py`, dimensions adapt to mm10):**
- DAE: Input(N_even_bins) → Dense(1024,ReLU) → Dropout(0.25) → Dense(512,ReLU) → Dense(256,ReLU) → Dropout(0.25) → Dense(128,sigmoid) [latent] → Dense(256,ReLU) → Dense(512,ReLU) → Dense(1024,ReLU) → Dense(N_even_bins,sigmoid)
- Classifier: Input(128) → Dense(64,ReLU) → Dropout(0.25) → Dense(16,ReLU) → Dropout(0.25) → Dense(4,softmax) [modified from 5 to 4 for mm10]

**Training split:** `inputM[:7000] / inputM[7000:]` (hardcoded in `pipeline/training.py`). mm10 has ~12,808 odd bins pre-crop → after crop ~9-11k → safely above 7000.

**SLURM:** `--cpus-per-task=16 --mem=96G --time=24:00:00 --account=csd940 --partition=shared`

**Outputs per timepoint:**
```
DATA_DIR/sniper_mm10/models_{tp}/
  odd_chrm_autoencoder.h5, odd_chrm_encoder.h5, odd_chrm_classifier.h5
  even_chrm_autoencoder.h5, even_chrm_encoder.h5, even_chrm_classifier.h5
```

**Verification:** All 6 .h5 files exist (each >1MB). Training log shows decreasing loss.

### B4 — Apply SNIPER to All Samples — PENDING

**Script:** `scripts/B4_apply_sniper.sb` + `scripts/B4_submit_apply.sh`

4 jobs: {ctrl_merged, mut_merged} × {250402, 250831}. Uses timepoint-specific trained models.

**SLURM:** `--cpus-per-task=8 --mem=64G --time=08:00:00`

**Outputs:**
```
DATA_DIR/sniper_mm10/predictions/{tp}/{sample}/
  predictions.mat      # softmax probabilities (N_bins × 4)
  predictions.bed      # BED9 with subcompartment labels
```

### B5 — SNIPER-CALDER2 Concordance — PENDING

**Script:** `scripts/B5_concordance_analysis.R` + `scripts/B5_run.sb`

**Logic:**
1. Load SNIPER `predictions.bed` and CALDER2 `{tp}_subcompartment_labels_100kb.tsv` for ctrl_merged
2. Map to common 100kb bins via `GenomicRanges::findOverlaps()`
3. Compute confusion matrix (SNIPER × CALDER2), Cohen's kappa, per-class precision/recall
4. For differential analysis: compare ctrl→mut transitions between both tools

**Expected concordance:** Since SNIPER is trained on CALDER2 labels, concordance on ctrl should be high (kappa > 0.7). On mut (unseen), concordance tests whether SNIPER's inter-chromosomal approach generalizes beyond training data.

**SLURM:** `--cpus-per-task=4 --mem=16G --time=02:00:00`

---

## Track C: Visualization & Integration

### C1 — Combined Subcompartment Comparison — PENDING

**Script:** `scripts/C1_combined_comparison.R`

- Genome-wide subcompartment karyotype track (CALDER2 ctrl vs mut, both timepoints)
- Transition matrix heatmap: ctrl → mut, counts colored by log10(n)
- Developmental comparison: early vs late transition rates
- Subcompartment stability analysis: what fraction of bins change between timepoints?

### C2 — Epigenomic Enrichment (Publication Figure) — PENDING

**Script:** `scripts/C2_epigenomic_enrichment.R`

Full SNIPER paper Figure 2c equivalent using CALDER2 calls + all available marks. Uses `ComplexHeatmap` for publication quality.

### C3 — Loop/Stripe × Subcompartment Integration — PENDING

**Script:** `scripts/C3_loops_stripes_integration.R`

Key biological question: do gained Polycomb loops (cluster pipeline clust5) sit in B.1 (facultative heterochromatin)?

**Logic:**
1. Load `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt` (38,948 loops × 6 clusters)
2. For each loop anchor midpoint, look up CALDER2 subcompartment label at 100kb
3. Cross-tabulate: cluster × subcompartment (chi-squared)
4. Test: clust5 (gained, 97% up_in_mutant) enriched in B.1 vs other clusters?
5. Repeat for stripes: `stripes/stripenn/outputs/{tp}/cross_res_merged.tsv`

### C4 — HOMER A/B Decomposition — PENDING

**Script:** `scripts/C4_homer_decomposition.R`

The key integrative analysis: decompose HOMER's coarse A/B transitions into subcompartment-resolution switches.

---

## Configuration

**`config/sniper_config.yaml`:**
```yaml
paths:
  code_root: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts/sniper"
  data_root: "/expanse/lustre/projects/csd940/zalibhai/sniper"
  hic_root: "/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic"
  sniper_source: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts/SNIPER"
  bigwig_dir: "/expanse/lustre/projects/csd940/zalibhai/bigwigs"

timepoints:
  early: "250831"
  late: "250402"

samples:
  merged: ["ctrl_merged", "mut_merged"]
  groups: ["ctrl", "mut"]

reference:
  genome: "mm10"
  chromosomes: ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"]
  chrom_sizes: "data/mm10.chrom.sizes"
  blacklist: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/mm10-blacklist.v2.bed"

calder2:
  bin_size: 50000
  genome: "mm10"
  conda_env: "calder2_env"
  n_cores: 8

sniper:
  resolution: 100000
  conda_env: "sniper_env"
  n_classes: 4
  juicer_tools: "/path/to/juicer_tools.jar"

existing:
  homer_compartments: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv"
  loop_clusters: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt"
  stripe_results: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn/outputs"

bigwigs:
  H3K27ac_ctrl: "H3K27acCtrl.bw"
  H3K27ac_mut: "H3K27acMut.bw"
  H3K27me3_ctrl: "H3K27me3Ctrl.bw"
  H3K27me3_mut: "H3K27me3Mut.bw"
  H3K4me3_ctrl: "H3K4me3Ctrl.bw"
  H3K4me3_mut: "H3K4me3Mut.bw"
  H2AK119ub_ctrl: "H2AK119ubCtrl.bw"
  H2AK119ub_mut: "H2AK119ubMut.bw"
  ATAC_ctrl: "ATACctrl.bw"
  ATAC_mut: "ATACmut.bw"
  RNA_ctrl: "RNActrl.bw"
  RNA_mut: "RNAmut.bw"
  DNAmethylation_ctrl: "DNAmethylationCtrl.bw"
  DNAmethylation_mut: "DNAmethylationMut.bw"
  H3K27me1_ctrl: "H3K27me1Ctrl.bw"
  H3K27me1_mut: "H3K27me1Mut.bw"

slurm:
  account: "csd940"
  partition: "shared"
  A1_calder2: { cpus: 8, mem: "64G", time: "12:00:00" }
  A2_diff: { cpus: 4, mem: "16G", time: "02:00:00" }
  A3_epigenomic: { cpus: 8, mem: "32G", time: "04:00:00" }
  A4_homer: { cpus: 4, mem: "16G", time: "02:00:00" }
  B1_cropmap: { cpus: 8, mem: "64G", time: "04:00:00" }
  B3_train: { cpus: 16, mem: "96G", time: "24:00:00" }
  B4_apply: { cpus: 8, mem: "64G", time: "08:00:00" }

filtering:
  exclude_chromosomes: ["X", "Y", "M"]
```

---

## Directory Structure

```
ML/cmpts/sniper/                          # CODE_DIR (GitHub-synced)
├── config/
│   └── sniper_config.yaml
├── scripts/
│   ├── A0_setup_calder2_env.sh           # Interactive: conda + R package install
│   ├── A1_run_calder2.R                  # CALDER2 wrapper R script
│   ├── A1_run_calder2.sb                 # SLURM worker (8c/64G/12h)
│   ├── A1_submit_calder2.sh              # Submit 4 jobs, print job IDs to stdout
│   ├── A2_differential_subcompartments.R # Ctrl vs mut transitions per tp
│   ├── A2_run.sb                         # SLURM wrapper (4c/16G/2h)
│   ├── A3_epigenomic_validation.R        # BigWig signal × subcompartment heatmap
│   ├── A3_run.sb                         # SLURM wrapper (8c/32G/4h)
│   ├── A4_integrate_homer_compartments.R # HOMER 25kb A/B × CALDER2 100kb
│   ├── A4_run.sb                         # SLURM wrapper (4c/16G/2h)
│   ├── B0_setup_sniper_env.sh            # Interactive: Python 3.7 + TF 1.15
│   ├── B1_adapt_sniper_mm10.py           # Generate mm10 variants + crop map
│   ├── B1_generate_cropmap.sb            # SLURM worker (8c/64G/4h)
│   ├── B2_generate_labels_from_calder2.R # CALDER2 → SNIPER .mat labels
│   ├── B2_run.sb                         # SLURM wrapper (4c/16G/1h)
│   ├── B3_train_sniper.sb                # SLURM worker (16c/96G/24h)
│   ├── B3_submit_train.sh               # Submit 2 jobs (1 per tp)
│   ├── B4_apply_sniper.sb               # SLURM worker (8c/64G/8h)
│   ├── B4_submit_apply.sh               # Submit 4 jobs
│   ├── B5_concordance_analysis.R         # SNIPER vs CALDER2 comparison
│   ├── B5_run.sb                         # SLURM wrapper (4c/16G/2h)
│   ├── C1_combined_comparison.R          # Genome-wide tracks + transitions
│   ├── C2_epigenomic_enrichment.R        # Publication-quality Fig 2c equiv
│   ├── C3_loops_stripes_integration.R    # Loop cluster × subcompartment
│   ├── C4_homer_decomposition.R          # A/B → subcompartment breakdown
│   ├── C_run.sb                          # SLURM wrapper for all C scripts
│   ├── run_full_calder2.sh               # Track A master driver
│   ├── run_full_sniper.sh                # Track B master driver
│   ├── run_full_integration.sh           # Track C master driver
│   └── utils/
│       └── multi_format_output.R         # save_multiformat_ggplot() utility
├── outputs/                              # Small results (GitHub-synced)
│   ├── calder2/                          # A2 transition TSVs, summaries
│   ├── sniper/                           # B5 concordance results
│   └── integration/                      # C1-C4 integration tables
└── docs/
    └── PLAN.md                           # This document

ML/cmpts/SNIPER/                          # SNIPER source (upstream clone)
├── data/
│   ├── mm10.chrom.sizes                  # Already exists
│   ├── hg19.chrom.sizes                  # Original
│   └── labels.mat                        # GM12878 hg19 labels (not used for mm10)
├── crop_map/
│   ├── cropMap.mat                       # hg19 original
│   ├── cropIndices.mat                   # hg19 original
│   ├── mm10_cropMap.mat                  # NEW: generated by B1
│   └── mm10_cropIndices.mat              # NEW: generated by B1
├── utilities/
│   ├── data_processing.py               # Original (unchanged)
│   ├── data_processing_mm10.py           # NEW: mm10-adapted copy
│   ├── interchromosome_matrix.py         # Original (unchanged)
│   ├── interchromosome_matrix_mm10.py    # NEW: mm10-adapted copy
│   ├── mm10_config.py                    # NEW: chromosome constants
│   └── input.py                          # Original (may need mm10 arg additions)
├── sniper_train_mm10.py                  # NEW: mm10 training entry point
├── sniper_apply_mm10.py                  # NEW: mm10 application entry point
└── pipeline/
    ├── models.py                         # Original (dynamic input dims, no change needed)
    ├── training.py                       # Original (7000 split, check mm10 bin count)
    └── application.py                    # Original (unchanged)

/expanse/.../sniper/                      # DATA_DIR (HPC-only, never synced)
├── calder2/
│   ├── 250402/
│   │   ├── ctrl_merged/                  # CALDER2 output
│   │   │   ├── sub_compartments/
│   │   │   │   ├── all_sub_compartments.bed
│   │   │   │   ├── all_sub_compartments.tsv
│   │   │   │   └── cor_with_ref.txt
│   │   │   └── intermediate_data/
│   │   └── mut_merged/                   # Same structure
│   └── 250831/                           # Same structure
├── sniper_mm10/
│   ├── dump_{tp}_{sample}/               # juicer_tools intermediate .txt
│   ├── models_{tp}/                      # 6 × .h5 trained models
│   ├── predictions/{tp}/{sample}/        # predictions.mat + predictions.bed
│   └── mm10_labels_{tp}.mat              # CALDER2-derived training labels
├── integration/                          # C-stage outputs (figures, large TSVs)
└── logs/                                 # SLURM logs
```

---

## Master Driver Scripts

### `run_full_calder2.sh` (Track A)

```bash
#!/bin/bash
# ML/cmpts/sniper/scripts/run_full_calder2.sh
# Master driver for Track A: CALDER2 subcompartment calling.
# Chains stages A1-A4 via SLURM --dependency=afterok.
# Run from login node: bash run_full_calder2.sh
#
# Prerequisites: A0_setup_calder2_env.sh completed successfully.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts/sniper"
cd "${CODE_DIR}/scripts"

echo "==========================================="
echo "CALDER2 Subcompartment Pipeline (Track A)"
echo "Started: $(date)"
echo "==========================================="

# A1: CALDER2 — 4 jobs
JIDS_A1=$(bash A1_submit_calder2.sh 2>/dev/null | paste -sd:)
if [ -z "${JIDS_A1}" ]; then echo "ERROR: A1 submission failed." >&2; exit 1; fi
echo "A1 jobs: ${JIDS_A1}"

# A2: Differential — 1 job (depends on all A1)
A2_JID=$(sbatch --parsable --dependency=afterok:${JIDS_A1} A2_run.sb)
echo "A2 job: ${A2_JID}"

# A3 + A4: parallel, both depend on A2
A3_JID=$(sbatch --parsable --dependency=afterok:${A2_JID} A3_run.sb)
A4_JID=$(sbatch --parsable --dependency=afterok:${A2_JID} A4_run.sb)
echo "A3 job: ${A3_JID} | A4 job: ${A4_JID}"

echo "==========================================="
echo "Track A submitted. Total: 4+1+1+1 = 7 jobs"
echo "Monitor: squeue -u $USER"
echo "==========================================="
```

### `run_full_sniper.sh` (Track B)

Same pattern. B1 (cropmap) → B2 → B3 → B4 → B5. B2 depends on A2 completion (manual check or explicit dependency).

---

## Files to Create

| File | Track | Type | Description |
|------|-------|------|-------------|
| `config/sniper_config.yaml` | — | Config | All paths, params, SLURM hints |
| `scripts/A0_setup_calder2_env.sh` | A | Setup | conda + R package install |
| `scripts/A1_run_calder2.R` | A | R | CALDER2 wrapper |
| `scripts/A1_run_calder2.sb` | A | SLURM | Worker: 8c/64G/12h |
| `scripts/A1_submit_calder2.sh` | A | Wrapper | Submit 4 jobs |
| `scripts/A2_differential_subcompartments.R` | A | R | Ctrl vs mut transitions |
| `scripts/A2_run.sb` | A | SLURM | Wrapper |
| `scripts/A3_epigenomic_validation.R` | A | R | BigWig × subcompartment |
| `scripts/A3_run.sb` | A | SLURM | Wrapper |
| `scripts/A4_integrate_homer_compartments.R` | A | R | HOMER integration |
| `scripts/A4_run.sb` | A | SLURM | Wrapper |
| `scripts/B0_setup_sniper_env.sh` | B | Setup | Python 3.7 + TF 1.15 |
| `scripts/B1_adapt_sniper_mm10.py` | B | Python | mm10 adaptation + crop map |
| `scripts/B1_generate_cropmap.sb` | B | SLURM | Worker: 8c/64G/4h |
| `scripts/B2_generate_labels_from_calder2.R` | B | R | CALDER2 → .mat labels |
| `scripts/B2_run.sb` | B | SLURM | Wrapper |
| `scripts/B3_train_sniper.sb` | B | SLURM | Worker: 16c/96G/24h |
| `scripts/B3_submit_train.sh` | B | Wrapper | Submit 2 jobs |
| `scripts/B4_apply_sniper.sb` | B | SLURM | Worker: 8c/64G/8h |
| `scripts/B4_submit_apply.sh` | B | Wrapper | Submit 4 jobs |
| `scripts/B5_concordance_analysis.R` | B | R | SNIPER vs CALDER2 |
| `scripts/B5_run.sb` | B | SLURM | Wrapper |
| `scripts/C1_combined_comparison.R` | C | R | Genome-wide tracks |
| `scripts/C2_epigenomic_enrichment.R` | C | R | Fig 2c equivalent |
| `scripts/C3_loops_stripes_integration.R` | C | R | Loop × subcompartment |
| `scripts/C4_homer_decomposition.R` | C | R | A/B breakdown |
| `scripts/C_run.sb` | C | SLURM | Wrapper for all C |
| `scripts/run_full_calder2.sh` | A | Driver | Track A master |
| `scripts/run_full_sniper.sh` | B | Driver | Track B master |
| `scripts/run_full_integration.sh` | C | Driver | Track C master |
| `scripts/utils/multi_format_output.R` | — | Utility | save_multiformat_ggplot() |

## Files to Modify (in SNIPER source, as copies)

| Original | mm10 Copy | Changes |
|----------|-----------|---------|
| `utilities/data_processing.py` | `utilities/data_processing_mm10.py` | `range(1,23,2)` → `range(1,20,2)` (×4 locations); `sizes_file` default → `mm10.chrom.sizes`; `np.int` → `int` |
| `utilities/interchromosome_matrix.py` | `utilities/interchromosome_matrix_mm10.py` | `range(1,23,2)` → `range(1,20,2)` (×2); `sizes_file` default → `mm10.chrom.sizes` |
| `pipeline/models.py` | (inline in `sniper_train_mm10.py`) | Classifier output: 5 → 4 neurons |

---

## Verification Plan

| Stage | Criterion | Expected |
|-------|-----------|----------|
| A0 | `library(CALDER)` loads | Exit 0 |
| A1 | 4× `all_sub_compartments.tsv` exist, all 19 chroms | ~25,600 bins at 100kb |
| A1 | Label distribution at depth 2 | A.1 ~25-35%, A.2 ~20-25%, B.1 ~15-25%, B.2 ~20-30% |
| A1 | Reference correlation | >0.3 per chromosome |
| A2 | Transition matrix non-uniform | Chi-squared p < 0.05 (late) |
| A3 | H3K27ac gradient A.1 → B.2 | Monotonic decrease |
| A3 | H3K27me3 enriched in B.1 | >1.5× fold enrichment |
| A4 | Cross-tab significant | Chi-squared p << 0.05 |
| B0 | TF 1.15 importable | `tensorflow.__version__` = 1.15.5 |
| B1 | mm10 crop map files exist | rowMap ~(11k,3), colMap ~(10k,3) |
| B2 | Label .mat files valid | `rows`/`cols` fields match crop dims |
| B3 | 6 .h5 model files per tp | Each >1MB, decreasing training loss |
| B4 | `predictions.bed` for all 4 samples | Correct bin count |
| B5 | SNIPER-CALDER2 concordance on ctrl | Cohen's kappa > 0.6 |
| C1-C4 | All figures in 4 formats | Non-zero file sizes |

---

## Critical Path and Timeline

```
Day 1: A0 (setup) → A1 submitted (runs overnight, ~8h)
Day 2: A1 complete → A2 (30min) → A3+A4 (2h) → Track A DONE
        Parallel: B0 (setup) + B1 (manual, ~2h interactive + 2h SLURM cropmap)
Day 2 evening: B2 (15min) → B3 submitted (runs overnight, ~16h)
Day 3: B3 complete → B4 (4h) → B5 (20min) → Track B DONE
Day 3 afternoon: C1-C4 (2h total) → PIPELINE COMPLETE
```

Total wall time: ~3 days with no failures. Track A alone completes in 1.5 days.

---

## Open Questions for Execution

1. **juicer_tools.jar location on Expanse** — SNIPER calls `java -jar {juicer_tools} dump`. The existing `lab/pipeline-scripts/juicer.sb` uses a Singularity container. Need to confirm standalone JAR location or extract from container.

2. **H3K36me3 BigWig** — Available as peaks (SEACR) but not as signal track in `sdsc/bigwigs/`. May need to generate from BAM, or use peak overlap fraction instead of continuous signal for A3 validation.

3. **CALDER2 conda install** — Verify whether `conda install -c bioconda r-calder2` works on Expanse, or if local source install (`install.packages("path/to/CALDER2", repos=NULL, type="source")`) is needed. The cloned repo at `ML/cmpts/CALDER2/` can be installed directly.

4. **TF 1.15 + Python 3.7 on Expanse** — Verify conda-forge still has Python 3.7 packages. If not, fall back to Python 3.6 + TF 1.12.2 (CPU).
