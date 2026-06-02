# Plan: CALDER2 + SNIPER Subcompartment Calling Pipeline for BAP1-KO Cerebellum

**Companion documents:** `ML/cmpts/repos/SNIPER/README.md` (SNIPER tool), `ML/cmpts/repos/CALDER2/README.md` (CALDER2 tool), `loops/docs/dixon-meeting-summary.md` (biological context), `tads/CLAUDE.md` (existing compartment pipeline).

**Directory convention:** Upstream tool repos live under `ML/cmpts/repos/` (SNIPER, CALDER2). Pipeline scaffolding (config, scripts, outputs, docs) lives directly under `ML/cmpts/`. This avoids cluttering the cloned repos with our pipeline code.

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
  B0: Python 3.6 + TF 1.12 env      │                                        |
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

### A0 — Environment Setup — DONE (2026-05-26)

**Scripts:** `scripts/A0_setup_calder2_env.sh` (conda env) + `scripts/utils/install_calder2_deps.R` (R package installation). Run interactively on Expanse login node.

**Status:** `calder2_env` conda env created and all packages installed. CALDER v2.0 loads successfully with mm10 reference BED (15,364 bins). Key deps verified: strawr 0.0.92, GenomicRanges 1.50.2, rhdf5 2.42.1, data.table 1.18.4, ggplot2 4.0.3, doParallel 1.0.17, igraph 2.0.3.

**Addendum (A3 prerequisite):** `rtracklayer` must also be installed for BigWig signal extraction in A3. Run interactively: `Rscript -e 'BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)'`

### A1 — Run CALDER2 on Each Sample — DONE (2026-05-27)

**Script:** `scripts/A1_run_calder2.R` + `scripts/A1_run_calder2.sb` + `scripts/A1_submit_calder2.sh`

4 jobs: {250402, 250831} × {ctrl_merged, mut_merged}.

**CALDER2 source patch (2026-05-27):** `repos/CALDER2/R/compartment_data_generation_fun.R:42-45` — changed `|`/`&` to `||`/`&&` and `class(x)==` to `inherits(x,)` in the normalization fallback chain. The 250402 merged .hic files lack KR normalization (only have VC_SQRT, VC). CALDER2's fallback logic (KR → VC_SQRT → VC) was broken: `try()` returns a try-error object, `nrow()` on it returns NULL, `NULL < 100` = `logical(0)`, and `TRUE | logical(0)` = `logical(0)` in R, crashing the `if` statement. The patch uses short-circuit `||`/`&&` so the `nrow()` branch is never evaluated when `try()` catches the error. Verified: chr1 completes with VC_SQRT fallback, correlation 0.68 against mm10 reference. **CALDER must be reinstalled from source after syncing the patch.**

**CALDER2 function call:**
```r
library(CALDER)

CALDER(
  contact_file_hic    = hic_path,            # .hic file, strawr extracts KR-normalized contacts
  chrs                = as.character(1:19),   # autosomes only, no "chr" prefix (CALDER strips it)
  bin_size            = 50000,               # 50kb — CALDER2 auto-extends to {50kb, 100kb}
  save_dir            = work_dir,
  save_intermediate_data = TRUE,             # needed for sub_domains and debugging
  genome              = "mm10",              # activates mm10 reference for A/B polarity
  n_cores             = 8,                   # parallel per-chromosome processing
  sub_domains         = FALSE,               # skip nested boundaries (not needed for subcompartments)
  single_binsize_only = FALSE                # use multi-resolution optimization (50kb + 100kb)
)
```

**Resolution choice:** `bin_size=50000` (50kb). CALDER2 automatically extends to `{50000, 100000}` when `genome='mm10'` and `single_binsize_only=FALSE`. The output `all_sub_compartments.tsv` contains bins at the optimal resolution per region. For uniform 100kb output (needed for SNIPER integration), we'll aggregate in A2.

**Key outputs per sample (in repo, rsync-able):**
```
ML/cmpts/outputs/calder2/{tp}/{sample}/
  sub_compartments/
    all_sub_compartments.bed     # BED9, hierarchical labels (e.g., A.1.1, B.2.2.1.2...)
    all_sub_compartments.tsv     # 6-col: chr, pos_start, pos_end, comp_name, comp_rank, continous_rank
    cor_with_ref.txt             # per-chr correlation with mm10 reference
    cor_with_ref.ALL.txt         # all bin sizes
    cor_with_ref.pdf             # correlation plot
```

**Large intermediates (HPC only, not synced):**
```
DATA_DIR/calder2/{tp}/{sample}/intermediate_data/   # per-chr per-resolution .Rds files (~1-3GB)
```

**SLURM logs (in repo):** `ML/cmpts/logs/A1_calder2_{tp}_{sample}_{jobid}.out`

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
- SLURM logs in `ML/cmpts/logs/` show exit 0

### A2 — Differential Subcompartment Analysis — DONE (2026-05-27)

**Script:** `scripts/A2_differential_subcompartments.R` + `scripts/A2_run.sb` + `scripts/A2_submit_differential.sh`

2 jobs: {250402, 250831}. Dependencies: All 4 A1 jobs complete.

**Logic (per timepoint):**
1. Load 2 `all_sub_compartments.tsv` files (ctrl + mut for this timepoint)
2. Truncate hierarchical labels to depth 2 (e.g., `A.1.1.2.1` → `A.1`) for 4-class analysis
3. Bin to uniform 100kb resolution: for variable-width CALDER2 output, assign each 100kb bin the label of the largest overlapping CALDER2 segment (plurality vote)
4. Join ctrl and mut labels on (chr, 100kb_bin_start)
5. Classify transitions: `ctrl_label → mut_label` for each bin
6. Compute transition matrix (4×4: A.1, A.2, B.1, B.2)
7. Chi-squared independence test + Cramer's V effect size

**Outputs (per timepoint):**
```
ML/cmpts/outputs/calder2/
  {tp}_subcompartment_labels_100kb.tsv     # chr, start, end, ctrl_label, mut_label, continous_rank_ctrl, continous_rank_mut, label_changed
  {tp}_transition_matrix.tsv               # 4x4 counts
  {tp}_transition_summary.tsv              # counts per transition type
  {tp}_transition_heatmap/                 # multi-format figure (pdf/svg/png/jpg)
  {tp}_transition_sankey/                  # multi-format figure
  {tp}_subcompartment_genome_pct/          # multi-format figure
```

**SLURM:** `--cpus-per-task=4 --mem=16G --time=02:00:00` (1 job per timepoint, 2 jobs total)

**Additional dependency:** `ggalluvial` R package (installed on first run if missing)

**Verification:**
- Transition matrices sum to ~24,639 (exact total of 100kb autosomal bins for mm10 chr1-19)
- Chi-squared p-value < 0.05 (late timepoint expected to show significant transitions)
- Key question answered: what fraction of bins change subcompartment? If <5%, the A/B system is largely stable; if >10%, there are substantial sub-compartment rearrangements

**Results (2026-05-27):**
- 250402 (late/adult): 15.3% bins changed (3,659/23,853), X²=39,957, p≈0, V=0.747
- 250831 (early/P12): 18.3% bins changed (4,371/23,823), X²=39,122, p≈0, V=0.740
- Both timepoints show substantial rearrangements (>10%). Early timepoint has more transitions than late (unexpected — may reflect developmental plasticity at P12).
- ~800 bins per timepoint uncallable (centromeric/telomeric gaps).

### A3 — Epigenomic Validation — DONE (2026-05-27)

**Script:** `scripts/A3_epigenomic_validation.R` + `scripts/A3_run.sb`

Single SLURM job processes both timepoints (BigWigs are not timepoint-specific). Produces ctrl validation, mut validation, and differential heatmaps per timepoint. Prerequisites: `rtracklayer` installed in `calder2_env`, A2 outputs exist, `merge_h3k36me3_bigwigs.sb` complete.

Replicates SNIPER paper Figure 2c: fold-enrichment heatmap of epigenomic marks per subcompartment.

**Available marks for validation (9 marks × ctrl/mut = 18 BigWigs):**

| Mark | BigWig | Subcompartment expectation |
|------|--------|--------------------------|
| H3K27ac | `{H3K27acCtrl,H3K27acMut}.bw` | A.1 >> A.2 > B.1 > B.2 |
| H3K4me3 | `{H3K4me3Ctrl,H3K4me3Mut}.bw` | A.1 >> A.2 > B.1 > B.2 |
| H3K36me3 | `{H3K36me3Ctrl,H3K36me3Mut}.bw` | A.1 > A.2 >> B.1 > B.2 (active gene bodies) |
| H3K27me3 | `{H3K27me3Ctrl,H3K27me3Mut}.bw` | B.1 >> B.2 > A.2 > A.1 |
| H2AK119ub | `{H2AK119ubCtrl,H2AK119ubMut}.bw` | B.1 enriched (Polycomb) |
| ATAC-seq | `{ATACctrl,ATACmut}.bw` | A.1 >> A.2 > B.1 > B.2 |
| RNA-seq | `{RNActrl,RNAmut}.bw` | A.1 > A.2 >> B.1 ≥ B.2 |
| DNA methylation | `{DNAmethylationCtrl,DNAmethylationMut}.bw` | Complex |
| H3K27me1 | `{H3K27me1Ctrl,H3K27me1Mut}.bw` | Not well characterized |

Missing: H3K79me2, H4K20me1, LMNB1. B.2 vs B.3 validation is limited without these.

BigWig path on HPC: `/expanse/lustre/projects/csd940/zalibhai/bigwigs/`

**Logic:**
1. Load `{tp}_subcompartment_labels_100kb.tsv` for ctrl samples
2. For each BigWig, use `rtracklayer::import.bw(bw, which=GRanges(bins))` to extract mean signal per 100kb bin
3. Group by depth-2 label, compute median signal per mark per subcompartment
4. Compute fold-enrichment: `median_in_subcomp / median_across_all_bins`
5. Plot as heatmap (rows=marks, cols=subcompartments A.1→B.2, color=RdBu fold-enrichment)

**SLURM:** `--cpus-per-task=8 --mem=32G --time=04:00:00` (single job, both TPs)

**Outputs (per timepoint):**
```
ML/cmpts/outputs/calder2/
  {tp}_bin_signals.tsv                    # ~24k rows × 23 cols (reusable by C-stage)
  {tp}_enrichment_matrix.tsv              # long format: mark, subcompartment, fold, log2
  {tp}_differential_matrix.tsv            # ctrl_label × log2(mut/ctrl) per mark
  {tp}_enrichment_heatmap_ctrl/           # 4 formats: pdf/svg/png/jpg
  {tp}_enrichment_heatmap_mut/
  {tp}_enrichment_heatmap_diff/
  {tp}_enrichment_combined/               # 3-panel combined figure
```

**Verification:**
- H3K27ac gradient monotonically decreasing A.1 → B.2 in ctrl
- H3K27me3 enriched in B.1 (>1.5× fold) — this is the key Polycomb compartment
- H2AK119ub enriched in B.1 (Polycomb signature)
- No NaN/Inf in enrichment results
- If gradients are flat or inverted, something is wrong with compartment polarity

**Results (2026-05-27):**
- Both timepoints validate: H3K27ac monotonically decreasing A.1→B.2 (late: 3.89→1.58→0.89→0.36; early: 2.74→1.24→0.68→0.29). All active marks (H3K4me3, ATAC, RNA) show strong A.1 enrichment.
- 250402 (late): H3K27me3 B.1=2.25× (strong Polycomb compartment); H2AK119ub B.1=1.25×. All checks passed.
- 250831 (early): H3K27me3 B.1=1.41× (marginally below 1.5 threshold — Polycomb may be less concentrated at P12); H2AK119ub B.1=0.86× (below 1.0, not enriched at early timepoint). RNA B.2=0.00 (no transcription in constitutive heterochromatin).
- Runtime: 2.8 min total (1.5 + 1.4 min per TP). 23,853/23,823 callable bins per TP.
- SLURM log: `logs/a3_epigenomic_49836042.out`

### A4 — Integration with HOMER A/B Compartments — DONE (2026-05-28)

**Script:** `scripts/A4_integrate_homer_compartments.R` + `scripts/A4_run.sb`

Single SLURM job processes both timepoints. Uses data.table arithmetic join (no GenomicRanges needed — 25kb bins tile 100kb bins exactly). Classifies each significant HOMER bin as a true compartment flip, within-compartment shift, subcompartment-stable, or uncallable.

**Input:**
- HOMER late: `tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv` (104,071 × 25kb bins)
- HOMER early: `tads/tad-pc-analysis/output/compartment_analysis_early/compartment_all_annotated.tsv` (101,684 × 25kb bins)
- CALDER2: `outputs/calder2/{tp}_subcompartment_labels_100kb.tsv` (from A2)

**Coordinate join:** `floor(Start / 100000) * 100000 + 1` maps HOMER 0-based Start to CALDER2 1-based bin_start. Merge on `(Chr, calder_bin_start) = (chr, bin_start)`.

**Logic:**
1. Load HOMER 25kb bins, filter to autosomes (chr1-19)
2. Left-join each 25kb bin to its containing CALDER2 100kb bin
3. For significant bins (adj_pvalue<0.05, |Difference|>0.30), classify transition:
   - `True_A_to_B`: CALDER2 ctrl=A.x → mut=B.y
   - `True_B_to_A`: CALDER2 ctrl=B.x → mut=A.y
   - `Within_A_shift`: A.1↔A.2 (quantitative weakening/strengthening)
   - `Within_B_shift`: B.1↔B.2
   - `Stable`: same CALDER2 label despite significant HOMER change
   - `Uncallable`: CALDER2 label is NA (centromeric/telomeric gaps)
4. Cross-tabulate: HOMER direction × CALDER2 ctrl_label × mut_label
5. Aggregate HOMER to 100kb (mean PC1, n_sig_bins) for C4 reuse

**Outputs (per timepoint):**
```
ML/cmpts/outputs/integration/
  {tp}_homer_calder2_joined.tsv          # all autosomal 25kb bins with CALDER2 labels
  {tp}_homer_100kb_aggregated.tsv        # HOMER aggregated to 100kb (C4 handoff)
  {tp}_homer_calder2_crosstab.tsv        # sig bins: direction × ctrl × mut counts
  {tp}_homer_weakening_decomposition.tsv # summary by transition category
  {tp}_homer_calder2_sankey/             # multi-format 3-axis alluvial figure
  {tp}_homer_calder2_dotplot/            # multi-format bubble plot
  combined_weakening_decomposition.tsv   # both TPs stacked
```

**SLURM:** `--cpus-per-task=4 --mem=16G --time=02:00:00` (single job, both TPs, ~5 min expected)

**Results (2026-05-28):**
- Runtime: 0.2 min total. SLURM log: `logs/a4_homer_integration_49838312.out`
- 250402 (late): 7,575 significant autosomal bins (8,189 total minus 614 on chrX/Y). 0% uncallable. 28.4% true compartment flips (2,152/7,575). B→A bins: 34.7% true B→A, 32.6% stable, 18.4% within-B, 14.2% within-A. A→B bins: 58.6% stable, 14.0% true A→B, 15.2% within-A, 11.3% within-B.
- 250831 (early): 5,036 significant autosomal bins. 0% uncallable. 29.9% true flips (1,508/5,036).
- **Polarity discrepancy (early timepoint, initial run):** The early HOMER `direction` column was systematically inverted relative to CALDER2 orientation — HOMER's PC1 eigenvector had the wrong sign for the 250831 analysis. CALDER2 orients against an mm10 reference BED; HOMER does not validate eigenvector polarity. This manifested as HOMER "B→A" bins showing 36.1% `True_A_to_B` by CALDER2 (and only 0.1% `True_B_to_A`).
- **Fix (2026-05-28):** Added gene-density polarity validation to `compartment_volcano_plot.R` (both copies) and `A4_integrate_homer_compartments.R`. The check compares mean ctrl PC1 for genic vs intergenic bins (A compartment must be gene-rich). If flipped, `ctrl_avg_PC1`, `mut_avg_PC1`, and `Difference` are negated before direction assignment. A4 re-run needed to produce corrected decomposition.
- Key biological finding: Jesse's hypothesis confirmed at late timepoint — ~71% of significant HOMER bins are NOT true compartment flips. They are quantitative weakening within subcompartments (within-A/B shifts) or subcompartment-stable (PC1 change without label change).

**Verification:**
- All checks passed for both timepoints
- 0% uncallable (excellent coordinate join coverage)
- All 13 output files written (4 TSV + 2 figure dirs per TP + 1 combined TSV)
- Polarity fix must be validated by re-running A4 on HPC (early TP should now show HOMER/CALDER2 direction agreement)

---

## Track B: SNIPER (Secondary Validation)

### B0 — Environment Setup — DONE (2026-05-28)

**Script:** `scripts/B0_setup_sniper_env.sh` (interactive, run once)

SNIPER requires Python 3.6 + TensorFlow-GPU 1.12.0 (pinned exact versions from upstream `requirements.txt`). Python 3.7 is explicitly unsupported by TF 1.12 via pip.

**Status:** `sniper_env` conda env created. TF-GPU 1.12.0 verified on Tesla V100-SXM2-32GB (gpu-debug partition, 31.7GB VRAM). No system CUDA modules needed — CUDA 9.0 runtime and cuDNN 7.x installed via conda into the env.

**Setup steps (run interactively on login node, except GPU test on gpu-debug):**
```bash
conda create -n sniper_env python=3.6.7 -c conda-forge -y
conda activate sniper_env
grep -v PyYAML repos/SNIPER/requirements.txt | pip install -r /dev/stdin
pip install PyYAML==5.4.1              # 3.13 fails to compile on Expanse gcc 7.2.0; 5.4.1 has a pre-built wheel
conda install cudatoolkit=9.0 -c nvidia
conda install "cudnn>=7,<8" -c nvidia  # TF 1.12 needs cuDNN 7.x; version resolved by conda against cudatoolkit 9.0
```

**GPU partition for SLURM jobs:** `--partition=gpu-shared --gpus=1` (or `gpu-debug` for interactive testing, 30 min max). No `module load cuda` needed — conda env provides CUDA/cuDNN runtime libraries.

**juicer_tools:** SNIPER calls `java -jar {juicer_tools} dump observed KR {hic} {chr1} {chr2} BP 100000 {output}`. The JAR is inside the Singularity container at `/cm/shared/apps/containers/singularity/juicer/juicer_2.0.1.sif` — path `/opt/juicer/CPU/common/juicer_tools.jar` (also at `/opt/scripts/common/juicer_tools.jar`, used by `abc/scripts/addnorm.sb`). SNIPER's `data_processing.py:45` shells out via `os.system()`, so the mm10 adapter will need to wrap calls through `singularity exec --bind /scratch,/expanse {container} java -jar {jar} dump ...`.

**Verification:** `python -c "import tensorflow as tf; print(tf.__version__); print(tf.test.is_gpu_available())"` → `1.12.0`, `True`

### B1 — Adapt SNIPER for mm10 — DONE (2026-06-01)

**Script:** `scripts/B1_adapt_sniper_mm10.py` + `scripts/B1_generate_cropmap.sb` (SLURM job)

This creates mm10-adapted copies of SNIPER's core modules. We modify copies, not the original SNIPER source, to preserve the upstream repo.

**Status (2026-05-28):** All 7 files written and import-tested locally. The 5 mm10-adapted modules (`utilities/mm10_config.py`, `utilities/data_processing_mm10.py`, `utilities/interchromosome_matrix_mm10.py`, `sniper_train_mm10.py`, `sniper_apply_mm10.py`) pass all local verification tests. The circular import between `data_processing_mm10` ↔ `interchromosome_matrix_mm10` is resolved via lazy import in `constructAndSave()`. Crop map generation (`B1_adapt_sniper_mm10.py` + `B1_generate_cropmap.sb`) awaits HPC execution. Both timepoints can run concurrently — output filenames include the timepoint suffix:
```bash
sbatch scripts/B1_generate_cropmap.sb 250402
sbatch scripts/B1_generate_cropmap.sb 250831
```

**Files to create in `ML/cmpts/repos/SNIPER/`:**

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

3. **Crop map generation** — `B1_adapt_sniper_mm10.py` also generates (timepoint-suffixed via `--timepoint`):
   - `crop_map/mm10_cropMap_{tp}.mat`: `rowMap` (N_odd × 3) and `colMap` (N_even × 3) mapping each retained bin to (original_index, chrom_number, bin_position)
   - `crop_map/mm10_cropIndices_{tp}.mat`: `odd_indices` and `even_indices` boolean arrays

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
- `mm10_cropMap_{tp}.mat` exists with `rowMap` shape (~11k, 3) and `colMap` shape (~10k, 3)
- `mm10_cropIndices_{tp}.mat` exists with `odd_indices` and `even_indices` arrays
- Test extraction: `python -c "from utilities.data_processing_mm10 import hicToMat; print('OK')"` from SNIPER directory

**Results (2026-06-01):**
- Both timepoints complete. Dimensions nearly identical: rowMap ~(10,686, 3), colMap ~(9,462–9,463, 3). ~80–83% bin retention after blacklist + sparse filtering.
- Each timepoint had 1 empty inter-chromosomal dump (different pairs); handled with zeros, not a blocker.

### B2 — Generate Ground Truth Labels from CALDER2 — DONE (2026-06-01)

**Script:** `scripts/B2_generate_labels_from_calder2.py` + `scripts/B2_run.sb`

Dependencies: A1 (CALDER2 output) + B1 (crop map).

**Implementation note (2026-06-01):** Changed from R (`.R`) to Python (`.py`) running in `sniper_env`. The .mat I/O requirement (reading B1 crop indices via `scipy.io.loadmat`, writing labels for B3's `scipy.io.loadmat`) made Python the natural choice. Avoids adding `R.matlab` dependency to `calder2_env` and eliminates cross-library .mat format compatibility risk. The plurality vote logic is reimplemented in pandas/numpy, mirroring A2's R `bin_plurality_vote()` exactly.

**Logic:**
1. Load CALDER2 `all_sub_compartments.tsv` for ctrl_merged (one per timepoint)
2. Truncate labels to depth 2 via regex `^([AB]\.[12])`: A.1→0, A.2→1, B.1→2, B.2→3 (4-class, no B3 in mm10)
3. Bin to 100kb via plurality vote (expand segments to bins, compute bp overlap, pick max-overlap label per bin)
4. For each 100kb bin in odd-chromosome order (chr1, chr3, ..., chr19), look up label; separately for even (chr2, chr4, ..., chr18)
5. Apply `mm10_cropIndices_{tp}` mask to keep only non-sparse bins; impute any gaps with chromosome-mode label
6. Save as `mm10_labels_{tp}.mat` with fields `rows` (odd labels, shape 1×K_odd) and `cols` (even labels, shape 1×K_even)

**Coordinate conversion:** CALDER2 1-based → SNIPER 0-based: `bin_0based = (pos_start - 1) // 100000`. Cross-validated against A2's `ceiling(pos_start / BIN_SIZE)` (1-based bin N = 0-based bin N-1).

**Mapping:** CALDER2 depth-2 → SNIPER integer labels:
| CALDER2 | SNIPER int | SNIPER name |
|---------|-----------|-------------|
| A.1 | 0 | A1 |
| A.2 | 1 | A2 |
| B.1 | 2 | B1 |
| B.2 | 3 | B2 |

**Decision: use 4-class classifier** (mm10 has no B.3). `mm10_config.N_CLASSES=4` and `Classifier_mm10` output layer uses 4 softmax neurons.

**SLURM:** `--cpus-per-task=4 --mem=16G --time=01:00:00`

**Verification (7 checks in script):**
1. Input existence: CALDER2 TSV and crop indices `.mat` exist and non-empty
2. Chromosome completeness: All 19 autosomes present in CALDER2 TSV
3. Label distribution: Per-class percentages; warn if any <5% or >60%
4. Dimension match: `len(rows) == len(odd_indices)`, `len(cols) == len(even_indices)`
5. No invalid labels: All values in {0,1,2,3}
6. Coverage: All retained bins labeled (uncallable bins imputed)
7. Training split: K_odd > 7000 and K_even > 7000

**Output:** `{CODE_DIR}/outputs/sniper/mm10_labels_{tp}.mat`

**Results (2026-06-01):**
- Runtime: ~8s per timepoint. SLURM logs: `logs/b2_labels_49961411.out` (late), `logs/b2_labels_49961412.out` (early).
- 250402 (late): rows=10,686, cols=9,463. 50+23 uncallable bins imputed. Labels: A.1=31.7%/33.5%, A.2=14.6%/14.2%, B.1=13.6%/13.9%, B.2=40.1%/38.4% (odd/even). B.2-heavy vs mm10 reference (~26%).
- 250831 (early): rows=10,686, cols=9,462. 45+25 uncallable bins imputed. Labels: A.1=33.1%/35.5%, A.2=20.0%/19.3%, B.1=16.4%/14.7%, B.2=30.5%/30.6% (odd/even). More balanced, closer to mm10 reference.
- All 7 validation checks passed for both timepoints. Dimension match exact. All labels in {0..3}.
- **Note:** Required `.astype(np.int64)` on bin arithmetic for Python 3.6 / old pandas compatibility (floor division produced float64).

### B3 — Train SNIPER on Control Merged — DONE (2026-06-01)

**Script:** `scripts/B3_train_sniper.sb` + `scripts/B3_submit_train.sh`

2 jobs: one per timepoint. Dependencies: B1 (crop maps) + B2 (labels).

**Training protocol (self-supervised, from SNIPER paper):**
1. Extract inter-chromosomal contacts from ctrl_merged.hic via juicer_tools dump (KR, 100kb) — 90 chromosome pairs (10 odd × 9 even), ~2h
2. Build contact probability matrix: P_ij = exp(-1 / (C_ij + 1e-10))
3. Trim with mm10 crop indices
4. Train denoising autoencoder (DAE): input=target=same ctrl_merged.hic (self-supervised — no high-coverage reference like GM12878 available)
5. Extract latent variables from encoder (128-dim)
6. Train MLP classifier on latent variables using CALDER2-derived labels (B2 output)

**Architecture (from `pipeline/models.py`, dimensions adapt to mm10):**
- DAE: Input(N_even_bins) → Dense(1024,ReLU) → Dropout(0.25) → Dense(512,ReLU) → Dense(256,ReLU) → Dropout(0.25) → Dense(128,sigmoid) [latent] → Dense(256,ReLU) → Dense(512,ReLU) → Dense(1024,ReLU) → Dense(N_even_bins,sigmoid)
- Classifier: Input(128) → Dense(64,ReLU) → Dropout(0.25) → Dense(16,ReLU) → Dropout(0.25) → Dense(4,softmax) [modified from 5 to 4 for mm10]

**Training split:** `inputM[:7000] / inputM[7000:]` (TRAIN_VAL_SPLIT=7000 in `mm10_config.py`). mm10 has ~12,808 odd bins pre-crop → after crop ~10,686 → safely above 7000.

**SLURM:** `--cpus-per-task=8 --mem=64G --gpus=1 --time=08:00:00 --account=csd940 --partition=gpu-shared`

**Invocation:**
```bash
# Both timepoints:
cd /expanse/.../mariner_hi-c/ML/cmpts
bash scripts/B3_submit_train.sh

# Chain after B2:
bash scripts/B3_submit_train.sh --dependency=afterok:${B2_JID_1}:${B2_JID_2}

# Single timepoint:
sbatch scripts/B3_train_sniper.sb 250402
```

**dump_dir dual purpose:** `sniper_train_mm10.py` saves both juicer dump intermediates AND .h5 models to `params['dump_dir']`. The SLURM wrapper sets dump_dir to `DATA_DIR/sniper_mm10/dump_train_{tp}` (large intermediates on Lustre), then copies the 6 small .h5 files to `CODE_DIR/outputs/sniper/models_{tp}/` post-training. The `-sm` flag also saves `input_matrix.mat` in dump_dir for B4 `-usemat` reuse. The `-ar` flag removes the 90 intermediate .txt files after matrix construction.

**DAE epoch fix (2026-06-01):** The upstream SNIPER code ships with `epochs=10` for the DAE, but the published paper (Xiong & Ma 2019, Methods p.10 and Fig 1a) specifies 25 epochs for both the DAE and classifier. Fixed `sniper_train_mm10.py` to use `epochs=25` for the DAE (classifier was already 25). This is especially important for our self-supervised protocol (input=target) where the DAE needs more training to learn meaningful latent structure without a separate high-coverage target.

**Import path fix (2026-06-01):** Added `sys.path` resolution (3-line `__file__`-relative pattern, same as B1) to `sniper_train_mm10.py` and `sniper_apply_mm10.py`. Both scripts previously depended on the caller to set PYTHONPATH; they are now self-contained. The SLURM wrapper also sets `PYTHONPATH` as belt-and-suspenders.

**Outputs per timepoint:**
```
CODE_DIR/outputs/sniper/models_{tp}/
  odd_chrm_autoencoder.h5, odd_chrm_encoder.h5, odd_chrm_classifier.h5
  even_chrm_autoencoder.h5, even_chrm_encoder.h5, even_chrm_classifier.h5
DATA_DIR/sniper_mm10/dump_train_{tp}/
  input_matrix.mat               # saved inter-chromosomal matrix (for B4 -usemat)
  target_matrix.mat              # identical to input (self-supervised)
```

**Verification:** All 6 .h5 files exist in models_{tp}/ (each >1MB). Training log shows decreasing DAE loss over 25 epochs and classifier trains successfully.

**Crop map path fix (2026-06-01):** `sniper_train_mm10.py` and `sniper_apply_mm10.py` originally hardcoded `crop_map/mm10_cropMap.mat` (no timepoint suffix), but B1 generates `mm10_cropMap_{tp}.mat`. Fixed: both scripts now require `-tp <timepoint>` flag and construct the correct timepoint-suffixed filenames. B3/B4 wrappers must pass `-tp ${TP}`.

**Results (2026-06-01):**
- 250402 (late): DAE converged by epoch ~20 (train_loss=0.458, val_loss=0.503). Val_loss minimum at epoch 21 (0.5034), slight uptick by epoch 25 (0.5094) — expected plateau behavior. All 6 models copied (~1MB each). Training time ~2s/epoch on V100.
- 250831 (early): Both timepoints complete. All 12 model files (6 per TP) in `outputs/sniper/models_{tp}/`.
- SLURM logs: `logs/b3_train_250402_*.out`, `logs/b3_train_250831_*.out`

### B4 — Apply SNIPER to All Samples — DONE (2026-06-01)

**Script:** `scripts/B4_apply_sniper.sb` + `scripts/B4_submit_apply.sh`

4 jobs: {ctrl_merged, mut_merged} × {250402, 250831}. Uses timepoint-specific trained models. Takes two positional args: `<timepoint> <sample>`.

**Bug fix (pre-requisite):** `sniper_apply_mm10.py:apply_on_mat_mm10()` was missing `predictionsToBed()` — the `.mat` code path saved `predictions.mat` but not `predictions.bed`. Fixed by adding 4 lines mirroring the `.hic` code path. All imports and `params['cropMap']` were already available.

**Input selection:** ctrl_merged reuses B3's saved inter-chromosomal matrix (`dump_train_{tp}/input_matrix.mat`, saved by B3's `-sm` flag), skipping ~2h of juicer dumps. The `.mat` extension auto-triggers `apply_on_mat_mm10()` via `get_application_params()`. mut_merged uses the `.hic` path (no pre-extracted matrix exists) and saves its matrix with `-sm` for potential reuse.

**SLURM:** `--cpus-per-task=8 --mem=64G --gpus=1 --time=04:00:00 --partition=gpu-shared`

Resource adjustment from original plan: 64G (was 32G) because mut_merged builds the full 12,808×11,831 matrix in memory, same as B3. 4h (was 2h) because juicer dumps alone take ~2h for 90 chromosome pairs.

**Invocation:**
```bash
# Both timepoints, both samples:
cd /expanse/.../mariner_hi-c/ML/cmpts
bash scripts/B4_submit_apply.sh

# Chained after B3:
B3_JIDS=$(bash scripts/B3_submit_train.sh | paste -sd:)
bash scripts/B4_submit_apply.sh --dependency=afterok:${B3_JIDS}

# Single job:
sbatch scripts/B4_apply_sniper.sb 250402 mut_merged
```

**Outputs:**
```
CODE_DIR/outputs/sniper/predictions/{tp}/{sample}/
  predictions.mat      # softmax probabilities (N_bins × 4)
  predictions.bed      # BED9 with subcompartment labels (~20,148 lines)
DATA_DIR/sniper_mm10/dump_apply_{tp}_{sample}/
  input_matrix.mat     # saved inter-chromosomal matrix (mut_merged only, via -sm)
```

**Verification:**
- All 4 `predictions.mat` files exist and are >1MB
- All 4 `predictions.bed` files have ~20,148 lines (10,686 odd + 9,462 even bins)
- BED contains all 4 labels (A1, A2, B1, B2) — no degenerate single-class predictions
- SLURM logs show exit 0
- ctrl_merged BED should closely match CALDER2 labels (model trained on this data)

**Results (2026-06-01):**
- All 4 jobs completed successfully. SLURM logs: `logs/b4_apply_{tp}_{sample}_*.out`
- ctrl_merged jobs used mat path (~16s each); mut_merged used .hic path with juicer dumps (~3.5 min each — much faster than B3's ~2h because inference skips training epochs)
- 250402 (late): ctrl=20,149 lines, mut=20,149 lines. Label distribution ctrl: B2=7,048, A1=5,540, A2=3,854, B1=3,707. Label distribution mut: B2=6,579, A1=5,184, A2=4,381, B1=4,005. Shift toward A2/B1 in mutant consistent with Polycomb derepression.
- 250831 (early): ctrl=20,148 lines, mut=20,148 lines. Label distribution ctrl: A1=6,553, A2=4,791, B2=4,522, B1=4,282. Label distribution mut: A1=6,580, A2=5,431, B2=4,210, B1=3,927. More balanced; A1 dominant in both conditions.
- 1-bin difference (20,149 vs 20,148) between timepoints reflects crop map differences (B1 generated per-timepoint crop indices from ctrl_merged.hic).
- Keras "No training configuration found" warning is benign — models loaded for inference only.
- Juicer "Development mode is enabled" warnings (90 per mut_merged job) are benign — standard juicer_tools stderr output.

### B5 — SNIPER-CALDER2 Concordance — READY (2026-06-01)

**Script:** `scripts/B5_concordance_analysis.R` + `scripts/B5_run.sb`

Both timepoints processed in a single job. Uses `calder2_env` (R + data.table). No GenomicRanges — coordinate join via data.table merge after converting SNIPER 0-based BED start to CALDER2 1-based bin_start (+1L). Label mapping: SNIPER `A1`/`A2`/`B1`/`B2` → CALDER2 `A.1`/`A.2`/`B.1`/`B.2`. Cohen's kappa computed manually (no `irr` dependency).

**Logic:**
1. Load SNIPER `predictions.bed` (BED9, 0-based) and CALDER2 `{tp}_subcompartment_labels_100kb.tsv` (1-based) for both ctrl and mut
2. Left-join on SNIPER bins (sparse ~20,149) → CALDER2 labels (dense ~24,639), excluding CALDER2 NA bins
3. Per-condition (ctrl, mut): confusion matrix (4×4), Cohen's kappa, per-class precision/recall/F1, accuracy
4. Differential transition concordance: intersect bins with all 4 labels (SNIPER ctrl+mut, CALDER2 ctrl+mut), classify into 5 categories: Both_stable, Both_change_agree, Both_change_disagree, SNIPER_only, CALDER2_only

**Expected concordance:** Ctrl kappa > 0.5 (trained on these labels). Mut kappa likely lower but > 0.3 if SNIPER generalizes.

**Outputs per timepoint (7 TSVs + 5 figure dirs):**
```
ML/cmpts/outputs/sniper/
  {tp}_concordance_ctrl.tsv              # per-bin: chr, bin_start, sniper_label, calder_label, concordant
  {tp}_concordance_mut.tsv               # same for mut
  {tp}_confusion_matrix_ctrl.tsv         # 4×4 matrix
  {tp}_confusion_matrix_mut.tsv          # 4×4 matrix
  {tp}_per_class_metrics.tsv             # long format: class, condition, precision, recall, f1
  {tp}_transition_concordance.tsv        # per-bin with all 4 labels + agreement category
  {tp}_transition_agreement_summary.tsv  # counts per category
  {tp}_confusion_heatmap_ctrl/           # multi-format figure
  {tp}_confusion_heatmap_mut/            # multi-format figure
  {tp}_discordant_alluvial/              # multi-format figure
  {tp}_per_class_concordance/            # multi-format figure
  {tp}_transition_agreement/             # multi-format figure
  combined_concordance_summary.tsv       # both TPs stacked: accuracy, kappa, per-class F1
```

**SLURM:** `--cpus-per-task=4 --mem=16G --time=02:00:00`

**Verification (7 checks):**
1. Join coverage ≥ 95% of SNIPER bins
2. Ctrl accuracy > 70%
3. Ctrl kappa > 0.5
4. No empty classes in confusion matrix
5. Mut kappa > 0.3
6. Both_stable > 80% of transition bins
7. All output files exist

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
  code_root: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts"
  data_root: "/expanse/lustre/projects/csd940/zalibhai/sniper"
  hic_root: "/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic"
  sniper_source: "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts/repos/SNIPER"
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
  juicer_container: "/cm/shared/apps/containers/singularity/juicer/juicer_2.0.1.sif"
  juicer_tools_jar: "/opt/juicer/CPU/common/juicer_tools.jar"

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
  B3_train: { cpus: 8, mem: "32G", gpus: 1, time: "04:00:00", partition: "gpu-shared" }
  B4_apply: { cpus: 8, mem: "32G", gpus: 1, time: "02:00:00", partition: "gpu-shared" }

filtering:
  exclude_chromosomes: ["X", "Y", "M"]
```

---

## Directory Structure

```
ML/cmpts/                                 # CODE_DIR (GitHub-synced)
├── PLAN.md                               # This document
├── config/
│   └── sniper_config.yaml
├── logs/                                 # SLURM logs (rsync to local, inspect in VS Code)
│   └── A1_calder2_{tp}_{sample}_{jobid}.out
├── scripts/
│   ├── A0_setup_calder2_env.sh           # Interactive: conda + R package install ✅
│   ├── A1_run_calder2.R                  # CALDER2 wrapper R script
│   ├── A1_run_calder2.sb                 # SLURM worker (8c/64G/12h)
│   ├── A1_submit_calder2.sh              # Submit 4 jobs, print job IDs to stdout
│   ├── A2_differential_subcompartments.R # Ctrl vs mut transitions (1 tp)
│   ├── A2_run.sb                         # SLURM worker (4c/16G/2h)
│   ├── A2_submit_differential.sh         # Submit 2 jobs, print job IDs to stdout
│   ├── A3_epigenomic_validation.R        # BigWig signal × subcompartment heatmap
│   ├── A3_run.sb                         # SLURM wrapper (8c/32G/4h)
│   ├── A4_integrate_homer_compartments.R # HOMER 25kb A/B × CALDER2 100kb
│   ├── A4_run.sb                         # SLURM wrapper (4c/16G/2h)
│   ├── B0_setup_sniper_env.sh            # Interactive: Python 3.6 + TF-GPU 1.12
│   ├── B1_adapt_sniper_mm10.py           # Generate mm10 variants + crop map
│   ├── B1_generate_cropmap.sb            # SLURM worker (8c/64G/4h)
│   ├── B2_generate_labels_from_calder2.py # CALDER2 → SNIPER .mat labels
│   ├── B2_run.sb                         # SLURM wrapper (4c/16G/1h)
│   ├── B3_train_sniper.sb                # SLURM worker (8c/64G/1GPU/8h, gpu-shared)
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
├── outputs/                              # All results <100MB (GitHub-synced, rsync-able)
│   ├── calder2/                          # A1 subcompartment calls + A2 transition TSVs
│   │   ├── {tp}/{sample}/sub_compartments/  # TSV, BED, cor_with_ref (copied from DATA_DIR)
│   │   └── ...
│   ├── sniper/                           # All Track B small outputs (<100MB)
│   │   ├── mm10_labels_{tp}.mat          # B2: CALDER2-derived training labels
│   │   ├── models_{tp}/                  # B3: 6 × .h5 trained models (~5MB each)
│   │   ├── predictions/{tp}/{sample}/    # B4: predictions.mat + predictions.bed
│   │   └── ...                           # B5: concordance results
│   └── integration/                      # C1-C4 integration tables
├── docs/                                 # Context docs
│
└── repos/                                # Upstream tool clones (do NOT add pipeline files here)
    ├── CALDER2/                           # CALDER2 v2.0 R package source
    └── SNIPER/                            # SNIPER upstream (mm10 files added by B1)
        ├── data/mm10.chrom.sizes          # Already exists
        ├── crop_map/
        │   ├── mm10_cropMap_{tp}.mat      # NEW: generated by B1 (per timepoint)
        │   └── mm10_cropIndices_{tp}.mat  # NEW: generated by B1 (per timepoint)
        ├── utilities/
        │   ├── data_processing_mm10.py    # NEW: mm10-adapted copy
        │   ├── interchromosome_matrix_mm10.py  # NEW: mm10-adapted copy
        │   └── mm10_config.py             # NEW: chromosome constants
        ├── sniper_train_mm10.py           # NEW: mm10 training entry point
        └── sniper_apply_mm10.py           # NEW: mm10 application entry point

/expanse/.../sniper/                      # DATA_DIR (HPC-only, large intermediates >100MB)
├── calder2/
│   ├── 250402/
│   │   ├── ctrl_merged/                  # CALDER2 work directory
│   │   │   ├── sub_compartments/         # small results (copied to repo by R script)
│   │   │   └── intermediate_data/        # large .Rds files (~1-3GB, stays here)
│   │   └── mut_merged/                   # Same structure
│   └── 250831/                           # Same structure
├── sniper_mm10/
│   └── dump_{tp}_{sample}/               # juicer_tools intermediate .txt (large, stays here)
├── integration/                          # C-stage outputs (figures, large TSVs)
└── logs/                                 # SLURM logs
```

---

## Master Driver Scripts

### `run_full_calder2.sh` (Track A)

```bash
#!/bin/bash
# ML/cmpts/scripts/run_full_calder2.sh
# Master driver for Track A: CALDER2 subcompartment calling.
# Chains stages A1-A4 via SLURM --dependency=afterok.
# Run from login node: bash run_full_calder2.sh
#
# Prerequisites: A0_setup_calder2_env.sh completed successfully.

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/ML/cmpts"
cd "${CODE_DIR}/scripts"

echo "==========================================="
echo "CALDER2 Subcompartment Pipeline (Track A)"
echo "Started: $(date)"
echo "==========================================="

# A1: CALDER2 — 4 jobs
JIDS_A1=$(bash A1_submit_calder2.sh 2>/dev/null | paste -sd:)
if [ -z "${JIDS_A1}" ]; then echo "ERROR: A1 submission failed." >&2; exit 1; fi
echo "A1 jobs: ${JIDS_A1}"

# A2: Differential — 2 jobs (1 per tp, depends on all A1)
JIDS_A2=$(bash A2_submit_differential.sh --dependency=afterok:${JIDS_A1} 2>/dev/null | paste -sd:)
if [ -z "${JIDS_A2}" ]; then echo "ERROR: A2 submission failed." >&2; exit 1; fi
echo "A2 jobs: ${JIDS_A2}"

# A3 + A4: parallel, both depend on all A2 jobs
A3_JID=$(sbatch --parsable --dependency=afterok:${JIDS_A2} A3_run.sb)
A4_JID=$(sbatch --parsable --dependency=afterok:${JIDS_A2} A4_run.sb)
echo "A3 job: ${A3_JID} | A4 job: ${A4_JID}"

echo "==========================================="
echo "Track A submitted. Total: 4+2+1+1 = 8 jobs"
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
| `scripts/A2_differential_subcompartments.R` | A | R | Ctrl vs mut transitions (1 tp) |
| `scripts/A2_run.sb` | A | SLURM | Worker: 4c/16G/2h |
| `scripts/A2_submit_differential.sh` | A | Wrapper | Submit 2 jobs |
| `scripts/A3_epigenomic_validation.R` | A | R | BigWig × subcompartment |
| `scripts/A3_run.sb` | A | SLURM | Wrapper |
| `scripts/A4_integrate_homer_compartments.R` | A | R | HOMER integration |
| `scripts/A4_run.sb` | A | SLURM | Wrapper |
| `scripts/B0_setup_sniper_env.sh` | B | Setup | Python 3.6 + TF-GPU 1.12 |
| `scripts/B1_adapt_sniper_mm10.py` | B | Python | mm10 adaptation + crop map |
| `scripts/B1_generate_cropmap.sb` | B | SLURM | Worker: 8c/64G/4h |
| `scripts/B2_generate_labels_from_calder2.py` | B | Python | CALDER2 → .mat labels |
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
| A1 | 4× `all_sub_compartments.tsv` exist, all 19 chroms | Variable-width segments (~237k rows); 24,639 bins after A2 100kb binning |
| A1 | Label distribution at depth 2 | A.1 ~25-35%, A.2 ~20-25%, B.1 ~15-25%, B.2 ~20-30% |
| A1 | Reference correlation | >0.3 per chromosome |
| A2 | Transition matrix non-uniform | Chi-squared p < 0.05 (late) |
| A3 | H3K27ac gradient A.1 → B.2 | Monotonic decrease |
| A3 | H3K27me3 enriched in B.1 | >1.5× fold enrichment |
| A4 | Cross-tab significant | Chi-squared p << 0.05 |
| B0 | TF-GPU 1.12 importable, GPU detected | `tensorflow.__version__` = 1.12.0, `is_gpu_available()` = True |
| B1 | mm10 crop map files exist (per tp) | `mm10_cropMap_{tp}.mat` rowMap ~(11k,3), colMap ~(10k,3) |
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

1. **juicer_tools.jar location on Expanse** — RESOLVED (2026-05-28). JAR is inside Singularity container `/cm/shared/apps/containers/singularity/juicer/juicer_2.0.1.sif` at `/opt/juicer/CPU/common/juicer_tools.jar`. SNIPER's mm10 adapter will call via `singularity exec`. Pattern established in `abc/scripts/addnorm.sb`.

2. **H3K36me3 BigWig** — RESOLVED (2026-05-26). Per-replicate BigWigs (6 ctrl + 6 mut, `_norm.bw` + `_rnorm.bw`) synced from EC2 (`/media/rs_256/normalization/`) to `sdsc/bigwigs/h3k36me3/`. On HPC at `/expanse/lustre/projects/csd940/zalibhai/bigwigs/h3k36me3/`.

3. **CALDER2 conda install** — Verify whether `conda install -c bioconda r-calder2` works on Expanse, or if local source install (`install.packages("path/to/CALDER2", repos=NULL, type="source")`) is needed. The cloned repo at `ML/cmpts/repos/CALDER2/` can be installed directly.

4. **TF 1.12 + Python 3.6 + CUDA/cuDNN on Expanse** — RESOLVED (2026-05-28). System CUDA modules (9.1, 10.2, 11.x) don't match TF 1.12's requirement for CUDA 9.0 exactly (`libcublas.so.9.0`). Solution: `conda install cudatoolkit=9.0 cudnn>=7,<8 -c nvidia` provides CUDA 9.0 runtime + cuDNN 7.x inside the conda env. No `module load cuda` needed. Verified: TF 1.12.0 GPU on Tesla V100-SXM2-32GB via `gpu-debug` partition.
