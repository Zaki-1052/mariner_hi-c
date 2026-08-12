# Plan: Adapt Popay Hi-C Loop Clustering Pipeline for BAP1-KO Cerebellum

**Companion document:** `CONTEXT-CLUSTER.md` contains full biological context, meeting notes, paper summaries, tool references, available data files, and module APIs. Read that file first for the "why"; this file covers the "how".

## Context

The Dixon meeting (2026-04-10) identified the Popay anchor-vs-span ChromHMM analysis as the highest-priority outstanding item for the BAP1-KO paper. Popay's code is in `cluster/`. The goal: reproduce Popay Figure 2f — ChromHMM chromatin state enrichment at loop anchors vs. loop spans across k-means clusters — adapted to our BAP1-KO vs wildtype mouse cerebellum data. This answers the central question: is Polycomb enrichment at anchors (sensitivity model) or across the loop body (extrusion impediment model)?

**Scope:** Late/adult timepoint (250402). All 39,344 nonredundant merged loops (multi-resolution). Replicates averaged to 2 columns (ctrl_merged, mut_merged). ChromHMM segmentation generated from our 5 ChIP-seq marks.

```
                          PIPELINE DEPENDENCY GRAPH

  Phase 0: Bug fixes ─────────────────────────────────────────────────────────
   (plotting.py, cooltools_called.py, chromHMM_heatmap.py,                   |
    deeptools_plotting.py, cluster_tools.py, deepTools_pipeline.py,          |
    bedpe_analysis.py)                                                       |
                                                                             |
  Phase 1: Data Prep ──────────────────────┐                                 |
   01_build_loop_count_file.py             |                                 |
   02_build_mm10_gene_annotation.R         |                                 |
                       ┌───────────────────┤                                 |
                       v                   v                                 |
  Phase 2: ChromHMM           Phase 3: K-means clustering                   |
   03_chromhmm_segmentation.sh  04_clustering.py                            |
   (BinarizeBed -> LearnModel    (elbow -> Cluster 3.0 -> sort)             |
    -> manual state naming)             |                                    |
           |                            |                                    |
           └─────────────┬──────────────┘                                    |
                         v                                                   |
                 Phase 4: Downstream analyses <──────────────────────────────┘
                  05_grouped_analyses.py
                  |-- 4.1 Loop size
                  |-- 4.2 Loop classification (CTCF-only, no RAD21)
                  |-- 4.3 ChIP signal in anchors (raw RPKM, no normalization)
                  |-- 4.4 ChromHMM anchor vs span enrichment  <-- KEY (Fig 2f)
                  |-- 4.5 ChromHMM proportions stacked bar
                  |-- 4.6 Gene annotation
                  |-- 4.7 DiffBind relationship (K27ac, K27me3, K119ub)
                  `-- 4.8 Cluster x differential status crosstab (NEW)
                         |
                         v
                 Phase 5: deepTools metagene
                  (H3K27me3/H3K27ac/H2AK119ub at anchors +/-5kb per cluster)
                         |
                         v
                 Phase 6: Cooltools pileup (DEFERRED -- needs mcool sync)
```

---

## Phase 0: Bug Fixes in Popay Pipeline Modules — DONE (2026-04-26)

Four Python modules crash on import due to hardcoded paths; one has a removed matplotlib style.

**Status:** All 8 fixes applied across 7 files. All modules import successfully under `cluster` conda env.

### 0.1 — Fix hardcoded `custom_params.json` path (4 files)

All four load `/Users/tessapopay/example_data/custom_params.json` at module top level. The file exists locally at `cluster/custom_params.json`. Replace with path relative to module location using `os.path.join(os.path.dirname(os.path.abspath(__file__)), 'custom_params.json')`.

Three of the four files need `import os` added (it's missing from their imports):

| File | Line(s) | Needs `import os` |
|------|---------|-------------------|
| `cluster/plotting.py` | 24-25 | Yes |
| `cluster/cooltools_called.py` | 18-19 | Yes |
| `cluster/chromHMM_heatmap.py` | 9-10 | No (already has it at line 6) |
| `cluster/deeptools_plotting.py` | 8-9 | Yes |

### 0.2 — Fix `plt.style.use('seaborn-poster')` crash

**File:** `cluster/cooltools_called.py:16`
**Bug:** `'seaborn-poster'` removed in matplotlib 3.6+.
**Fix:** Replace with `plt.style.use('seaborn-v0_8-poster')` (the renamed equivalent in matplotlib 3.6+). This preserves the larger font/figure sizing Popay intended, which `custom_params` alone does not provide (custom_params sets tick/label sizes but not the base figure/font scale).

### 0.3 — Fix `elbow()` blank PNG

**File:** `cluster/cluster_tools.py:18-19`
**Bug:** `plt.show()` (line 18) clears figure before `plt.savefig()` (line 19).
**Fix:** Swap order — `savefig` first, then `show`. Add `plt.close()` after.

### 0.4 — Fix `joint()` logX zero-check on wrong column

**File:** `cluster/plotting.py:501-502`
**Bug:** The `if logX:` block checks `data[data[ycol] == 0]` but should check `data[data[xcol] == 0]`.
**Fix:** Change `ycol` to `xcol` on line 502.

### 0.5 — Remove hardcoded biology overrides

| File | Lines | What to remove |
|------|-------|---------------|
| `cluster/plotting.py` | 132-139 | 8 `if dataset == 'TFAP2A'`/`'MYC'`/etc. ylim overrides |
| `cluster/deeptools_plotting.py` | 108 | `if bam == 'HA-NIPBL': line.set_ylim([0,2])` |

### 0.6 — Fix hg38 genome size in deepTools_pipeline.py

**File:** `cluster/deepTools_pipeline.py:15`
**Bug:** `--effectiveGenomeSize 2913022398` is hg38.
**Fix:** Add `genome='mm10'` parameter to `bam_coverage()`, with lookup dict `{'mm10': 2494787188, 'hg38': 2913022398}`.

### 0.7 — Fix cooltools default genome

**File:** `cluster/cooltools_called.py:27`, `mcool_pileup()` signature
**Fix:** Change `genome='hg38'` default to `genome='mm10'`.

### 0.8 — Fix bedpe_analysis.py default gene path and FPKM guard

**File:** `cluster/bedpe_analysis.py`
- Line 67: Change `default_gene_path = '~/example_data/gencode.v25.annotation_pp.bed'` to `os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'mm10_knownGene_pp.bed')`
- Lines 142-159: Add `if FPKM_df is not None:` guard around the FPKM `.to_csv()` and plotting calls. Currently crashes if FPKM_df is None because it tries to group/plot on a column that doesn't exist.

### Phase 0 Verification

```bash
cd cluster
python3 -c "import plotting; import cluster_tools; import chromHMM_heatmap; import deeptools_plotting; import bedpe_analysis; print('All imports OK')"
```

**Verified 2026-04-26:** All 7 modules (plotting, cluster_tools, chromHMM_heatmap, deeptools_plotting, bedpe_analysis, cooltools_called, deepTools_pipeline) import cleanly. No remnants of hardcoded Popay paths, deprecated styles, or hg38 defaults remain (confirmed via grep).

---

## Phase 1: Data Preparation — DONE (2026-04-26)

**Status:** All three artifacts produced and verified end-to-end. Loop count file (`cluster/data/late_merged_loop_counts.txt`, 39,344 × 8), differential-stats sidecar (`cluster/data/late_merged_loop_metadata.tsv`, 39,344 × 16), and mm10 promoter BED (`cluster/data/mm10_knownGene_pp.bed`, 24,515 × 7) all match expected schemas. Output directory tree in `cluster/bap1_late/` is created.

**Corrections applied during execution (parent plan was wrong about these — Phase 3+ must follow the corrections, not the original plan):**

- **Column naming.** Count file uses `ctrl_merge` / `mut_merge` (NOT `ctrl_merged` / `mut_merged`). Popay's downstream notebook code does `treatment_group.str.replace('_merge','')`; the `_merged` suffix from the parent plan would produce broken `ctrld` / `mutd` display labels. Phase 3 must use `normalize_col = 'ctrl_merge'` and `palette = {'ctrl_merge': 'darkgrey', 'mut_merge': 'forestgreen'}`.
- **Filter threshold (Phase 3.1) must be re-derived empirically.** The parent plan's `filter_threshold = 0.008` is calibrated for Popay's KR-balanced contact frequencies (~0.05 scale). Our mariner-aggregated averaged counts span a much larger scale: ctrl_merge median 451, min 5.50, 1%-tile 49, 99%-tile 4,527 overall — and the per-resolution distributions differ ~6× (5kb median 171, 10kb 403, 25kb 984). A single absolute-value filter would disproportionately remove 5kb loops. Phase 3 should choose a percentile-based or per-resolution threshold from the actual `ctrl_merge` distribution.
- **R environment.** `mariner_env` does not exist on this Mac (`conda env list` shows only base, cleave-pipeline, cluster, hic, hictk). All required Bioconductor packages — `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`, `GenomicRanges`, `GenomicFeatures`, `GenomeInfoDb` — are in **system R 4.5.2** at `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library`. Invoke as `/usr/local/bin/Rscript`. Note: `keepSeqlevels()` is exported from `GenomeInfoDb`, not `GenomicRanges` — load `GenomeInfoDb` explicitly.
- **Python invocation.** `conda run -n cluster python3` does not activate the env reliably (pandas import fails). Use the absolute interpreter path `/opt/homebrew/anaconda3/envs/cluster/bin/python3` for all Phase 2+ Python scripts. Verified: pandas 1.5.3, bioframe 0.6.1, pyBigWig, numpy 1.24.4, sklearn all importable.
- **`pd.to_csv` line terminator.** pandas 1.5.3 has deprecated `line_terminator` in favor of `lineterminator` (one word). Use `lineterminator='\n'` for forward compatibility.
- **`merged_all_loops_nonredundant.bedpe` schema.** 32 columns (header confirmed). `kept_from_resolution` is **column 30** (kb units: 5/10/25), not column 29 as easy off-by-one might suggest. Distribution: 5kb=7,901 / 10kb=14,553 / 25kb=16,890. Coordinates in this merged file are byte-identical to the per-resolution edgeR results (verified: row 2 = `loop_2` in 5kb edgeR), enabling clean coord-keyed joins via the loop_id-bridge pattern.

### 1.1 — Build loop count file in Popay format

**Script:** `cluster/scripts/01_build_loop_count_file.py`

**Use the multi-resolution merged loop list** (not 5kb alone):
- `outputs/250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe` — 39,344 nonredundant loops (5kb: 7,901; 10kb: 14,553; 25kb: 16,890). Each loop has a `kept_from_resolution` field indicating which resolution's data to use.

**Per-resolution count matrices (for pulling raw counts):**
- `outputs/250402-late_outputs/res_5kb/06_counts_matrix.tsv` (17,983 loops x 6 replicates)
- `outputs/250402-late_outputs/res_10kb/06_counts_matrix.tsv` (22,639 loops x 6 replicates)
- `outputs/250402-late_outputs/res_25kb/06_counts_matrix.tsv` (20,421 loops x 6 replicates)
- Corresponding `edgeR_results_res_{RES}kb/primary_analysis/all_results_primary.tsv` (maps loop_id -> coordinates)

**Logic:**
1. Load merged BEDPE — provides nonredundant coordinates + `resolution`, `kept_from_resolution`, `logFC`, `FDR`, `direction`
2. For each resolution, load the count matrix (index_col=0) and edgeR results to build a loop_id -> (coordinates, counts) mapping
3. For each loop in the merged set, pull counts from its `kept_from_resolution` count matrix via coordinate matching (chr1/start1/end1/chr2/start2/end2)
4. Average replicates: `ctrl_merged = mean(ctrl_M1, ctrl_M2, ctrl_M3)`, `mut_merged = mean(mut_M1, mut_M2, mut_M3)`
5. Output Popay-format: `chr1  x1  x2  chr2  y1  y2  ctrl_merged  mut_merged`
6. Output metadata sidecar with edgeR stats for Phase 4.8

**Outputs:**
- `cluster/data/late_merged_loop_counts.txt` (~39,344 rows, 8 columns, tab-separated with header)
- `cluster/data/late_merged_loop_metadata.tsv` (same rows, differential stats + resolution)

### 1.2 — Generate mm10 gene annotation BED

**Script:** `cluster/scripts/02_build_mm10_gene_annotation.R`

mm10 equivalent of Popay's `gencode.v25.annotation_pp.bed` (7 columns: chr, start, end, gene_id, score, strand, gene_name; TSS +/-750bp).

**Logic:**
1. Load `TxDb.Mmusculus.UCSC.mm10.knownGene` and `org.Mm.eg.db`
2. Extract gene coordinates, compute promoters (upstream=750, downstream=750)
3. Map Entrez IDs to gene symbols
4. Output 7-column BED matching Popay format

**Output:** `cluster/data/mm10_knownGene_pp.bed` (~24K gene entries)

### 1.3 — Create output directory structure

```bash
mkdir -p cluster/{scripts,data}
mkdir -p cluster/bap1_late/{cluster3,figures/{annotation,ChIP_intersect,deeptools,deeptools_input},chromHMM/{binarized,learned_model,anchor_input,span_input,peak_beds},cooltools}
```

### Phase 1 Verification

```bash
bash cluster/scripts/run_phase1.sh
```

**Verified 2026-04-26:**

- `cluster/data/late_merged_loop_counts.txt`: 39,344 rows × 8 cols. Columns exactly `chr1 x1 x2 chr2 y1 y2 ctrl_merge mut_merge`. No NaN. All `ctrl_merge` and `mut_merge` strictly positive (mins 5.50 / 6.82). Row 1 spot-check matches `loop_2` in `res_5kb/06_counts_matrix.tsv` exactly: 146.09 / 150.30 = arithmetic mean of `[ctrl_M1, ctrl_M2, ctrl_M3]` and `[mut_M1, mut_M2, mut_M3]` respectively.
- `cluster/data/late_merged_loop_metadata.tsv`: 39,344 × 16. Inner-joins to count file on `(chr1, x1, x2, chr2, y1, y2)` with zero loss. Direction split: unchanged=31,363, up_in_mutant=4,253, down_in_mutant=3,728 — matches CONTEXT-CLUSTER §11. Resolution split: 5kb=7,901 / 10kb=14,553 / 25kb=16,890.
- `cluster/data/mm10_knownGene_pp.bed`: 24,515 × 7. All rows exactly 1,500 bp wide (`end - start == 1500`). Score column constant 0. Strand split: 12,318 (+) / 12,197 (-). Standard chroms only (chr1–19, chrX, chrY); chrM has no genes in TxDb. 19 entries lack a symbol mapping → `gene_name` left blank, Entrez ID retained in `gene_id` (do not drop these — `bedpe_analysis.bedtools_annotation()` tolerates blanks).
- Notebook compatibility smoke test: `from cluster_tools import comparison_type; comparison_type(['ctrl_merge','mut_merge'])` returns `('multiple', ['ctrl', 'mut'])` (set order may vary). `[c.replace('_merge','') for c in data_cols]` yields `['ctrl', 'mut']` — clean labels for `treatment_group` axis.
- Output directory tree present under `cluster/bap1_late/{cluster3, figures/{annotation,ChIP_intersect,deeptools,deeptools_input}, chromHMM/{binarized,learned_model,anchor_input,span_input,peak_beds}, cooltools}`.

**Files created in Phase 1:**

| File | Purpose |
|------|---------|
| `cluster/scripts/01_build_loop_count_file.py` | Build Popay-format count file + metadata sidecar via loop_id-bridge join through edgeR results |
| `cluster/scripts/02_build_mm10_gene_annotation.R` | mm10 promoter BED, 7-col gencode-style, TSS ± 750 bp |
| `cluster/scripts/run_phase1.sh` | Driver: mkdir + run both scripts + summary |

---

## Phase 2: ChromHMM Segmentation — DONE (2026-04-27)

**Status:** 12-state model learned from 5 ChIP-seq peak BEDs on 21 standard chromosomes. Segmentation (`cerebellum_late_12_segments.bed`, 335,569 segments), emission matrix (`emissions_12.txt`, 12 × 5), and rename file (`12state_rename_cerebellum.txt`, 12 biological names) all produced and verified. Smoke-tested rename file via OverlapEnrichment on all 39,344 loop anchors — biologically coherent enrichments (Active_Promoter 5.31×, Quiescent 0.87×). Model shows no degeneracy: all 12 states have unique emission profiles (every pair differs by ≥1 fully-on vs fully-off mark), clean convergence in 32 iterations, no state >50% of segments.

**Corrections applied during execution (Phase 4+ must follow these, not the original plan):**

- **Filtered chromsizes.** Used `cluster/bap1_late/chromHMM/mm10_standard.txt` (21 chroms: chr1-19 + chrX + chrY) instead of the full `cluster/ChromHMM/CHROMSIZES/mm10.txt` (which includes chrM, chr*_random, chrUn_* contigs). Peak BEDs were also filtered to these 21 chroms. This avoids segmenting ~45 contigs that no loops or genes ever touch.
- **Peak BED chromosome contamination.** 4 of 5 source BED files (all except H3K4me3) contained peaks on `chr*_random`, `chrUn_*`, or `chrM`. Filtered during staging: H3K27ac 15,105→15,100; H3K27me3 15,809→15,801; H3K4me1 113,781→113,722; CTCF 32,487→32,474; H3K4me3 unchanged at 6,581.
- **State ID prefix is `E`, not `U`.** ChromHMM v1.27 emits state IDs as `E1`–`E12` in segments.bed and emissions_12.txt. Popay's example uses `U1`–`U12` (from an older version or relabeling). The rename file uses `E` prefix verbatim.
- **Runtime: 93 seconds, not 30-90 min.** LearnModel converged in 32 iterations (total 71 sec) with 4 threads. The mm10 genome with 5 binary marks produces a small enough state space that convergence is fast.
- **Segments.bed has 335,569 lines, not ~14M.** The original plan estimate (~14M = genome/200bp) was the per-bin count; the segments BED merges adjacent bins of the same state into contiguous regions.
- **No blacklist pre-filter applied.** ChromHMM treats blacklist regions correctly as no-signal bins. Blacklist filtering deferred to Phase 5 (deepTools `--blackListFileName`).
- **Rename file uses underscored names** (e.g. `Active_Promoter`, not `Active promoter`) for cross-script string-matching with the project's existing 7-category annotation taxonomy and Phase 4 palette dicts.
- **ChromHMM wrapper `-mx` deprecation warning.** Java 25 emits a benign warning about the deprecated `-mx8G` flag in `cluster/ChromHMM/chromhmm`. Functionally harmless; could update to `-Xmx8G` but not in scope.
- **Peak BEDs are control-condition data.** The 5 ChIP-seq peak files fed into ChromHMM are from wildtype/control samples, published in the PLOS Genetics paper. The model therefore represents the **control chromatin landscape**, not a condition-agnostic or pooled reference. Downstream enrichment (Phase 4.4) asks "what control-condition chromatin states do the differential loop anchors/spans sit in?"

### 2.1 — Prepare peak BEDs for BinarizeBed

**Script:** `cluster/scripts/03_chromhmm_segmentation.sh`

5 peak BEDs (verified on disk):
| Mark | Source file | Peaks |
|------|-----------|-------|
| H3K27ac | `peaks/beds/H3K27acCerebellumLate2.bed` | 15,105 |
| H3K27me3 | `peaks/beds/H3K27me3CerebellumLate1.bed` | 15,809 |
| H3K4me1 | `peaks/beds/H3K4me1CerebellumLate1.bed` | 113,781 |
| H3K4me3 | `peaks/beds/H3K4me3CerebellumLate2.bed` | 6,581 |
| CTCF | `peaks/CTCF.bed` | 32,487 |

Stage 3-column BEDs into `cluster/bap1_late/chromHMM/peak_beds/` and create cellmarkfiletable.

### 2.2 — Run ChromHMM BinarizeBed

ChromHMM.jar verified working. mm10.txt chromsizes present at `cluster/ChromHMM/CHROMSIZES/mm10.txt`.

```bash
java -mx4G -jar cluster/ChromHMM/ChromHMM.jar BinarizeBed \
    -peaks \
    cluster/ChromHMM/CHROMSIZES/mm10.txt \
    cluster/bap1_late/chromHMM/peak_beds/ \
    cluster/bap1_late/chromHMM/cellmarkfiletable.txt \
    cluster/bap1_late/chromHMM/binarized/
```

### 2.3 — Run ChromHMM LearnModel (12 states)

```bash
java -mx8G -jar cluster/ChromHMM/ChromHMM.jar LearnModel \
    -p 4 \
    -l cluster/ChromHMM/CHROMSIZES/mm10.txt \
    cluster/bap1_late/chromHMM/binarized/ \
    cluster/bap1_late/chromHMM/learned_model/ \
    12 \
    mm10
```

Runtime: ~30-90 min. Key outputs: `cerebellum_late_12_segments.bed`, `emissions_12.txt`.

### 2.4 — Interpret emission matrix and create rename file

**Manual step.** Inspect `emissions_12.txt` to assign biological names in 2-column TSV format matching `cluster/clustering_example_data/12state_rename.txt` (state_id e.g. `U1`, biological_name e.g. `Active_Promoter`).

**Output:** `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt`

### Phase 2 Verification
```bash
wc -l cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed
# Expect ~335K lines (region-merged segments, not per-bin)
```

**Verified 2026-04-27:**

- `cluster/bap1_late/chromHMM/peak_beds/{H3K27ac,H3K27me3,H3K4me1,H3K4me3,CTCF}.bed`: 5 staged 3-col BEDs, standard chroms only. Total 183,678 peaks.
- `cluster/bap1_late/chromHMM/mm10_standard.txt`: 21 chroms, sum 2,725,519,400 bp.
- `cluster/bap1_late/chromHMM/cellmarkfiletable.txt`: 5 rows, cell = `cerebellum_late`.
- `cluster/bap1_late/chromHMM/binarized/`: 21 binary files (`cerebellum_late_chr{1..19,X,Y}_binary.txt`).
- `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed`: 335,569 segments across all 21 chroms. 13,627,597 total 200bp bins covered.
- `cluster/bap1_late/chromHMM/learned_model/emissions_12.txt`: 12 states × 5 marks. All emission probabilities are cleanly "digital" (near 0 or 1), no blurred intermediate states. Convergence: 32 iterations, log-likelihood -1,617,268.978.
- State distribution (segment counts): E1=136,333 / E3=96,673 / E10=21,252 / E11=20,365 / E9=16,476 / E4=12,851 / E2=6,896 / E12=6,763 / E8=6,214 / E5=5,511 / E7=3,383 / E6=2,852. All 12 states meaningfully populated. E1 Quiescent = 40.6% of segments / 92.0% of genome (expected for 5-mark model).
- `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt`: 12 lines, tab-separated, state_id `E{N}` → underscored biological name. Ordered thematically (active → repressive → structural → quiescent):

| State | Genome % | Marks (≥0.95 prob) | Name | TSS ×fold | CpG ×fold |
|-------|----------|--------------------|------|-----------|-----------|
| E8 | 0.35% | K27ac + K4me3 | Active_Promoter | 45.07 | 83.76 |
| E5 | 0.20% | K27ac + K4me3 + K4me1 | Active_Promoter_Flank | 19.28 | 21.15 |
| E6 | 0.08% | K4me3 + K4me1 | Poised_Promoter | 19.47 | 20.66 |
| E7 | 0.11% | K4me3 (+ slight K27me3) | Weak_Promoter | 51.97 | 88.89 |
| E9 | 0.36% | K27ac alone | Strong_Enhancer | 7.18 | 12.70 |
| E4 | 0.77% | K27ac + K4me1 | Active_Enhancer | 3.00 | 2.33 |
| E3 | 3.67% | K4me1 alone | Poised_Enhancer | 2.32 | 1.60 |
| E12 | 0.31% | K4me1 + K27me3 | Bivalent_Enhancer | 23.33 | 49.75 |
| E11 | 1.81% | K27me3 alone | Polycomb | 2.93 | 3.79 |
| E2 | 0.09% | K4me1 + CTCF | CTCF_Boundary | 9.44 | 8.98 |
| E10 | 0.30% | CTCF alone | Insulator | 4.66 | 3.94 |
| E1 | 91.97% | none | Quiescent | 0.49 | 0.21 |

- Smoke-test OverlapEnrichment (all 39,344 anchor1 regions): Active_Promoter 5.31×, Active_Promoter_Flank 4.90×, Poised_Promoter 4.28×, Weak_Promoter 4.32×, Strong_Enhancer 3.67×, Active_Enhancer 3.45×, Poised_Enhancer 2.06×, Bivalent_Enhancer 2.02×, Polycomb 1.56×, CTCF_Boundary 4.11×, Insulator 3.41×, Quiescent 0.87× (depleted). Biologically coherent — active/structural states enriched at loop anchors, quiescent depleted.
- **Model degeneracy assessment: k=12 is well-justified.** Every pair of states differs by ≥1 mark being fully on vs fully off (no near-duplicate emissions). No state >50% of segments. Small genome-% states (E2 0.09%, E6 0.08%, E7 0.11%) have very high TSS/CpG enrichments (up to 89×), proving they mark real genomic features, not noise. A k=10 model would merge E5↔E8 (losing TSS-flank vs core distinction) and E6↔E7 (losing K4me1-flanked vs isolated promoter distinction) — both biologically meaningful.

**Files created in Phase 2:**

| File | Purpose |
|------|---------|
| `cluster/scripts/03_chromhmm_segmentation.sh` | Worker: stage BEDs, filter chromsizes, BinarizeBed + LearnModel(k=12) |
| `cluster/scripts/run_phase2.sh` | Driver: call worker with status check, summary |
| `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` | State ID → biological name mapping (12 lines) |
| `cluster/bap1_late/chromHMM/mm10_standard.txt` | 21-row filtered chromsizes |
| `cluster/bap1_late/chromHMM/cellmarkfiletable.txt` | 5-row mark → BED table |
| `cluster/phase2.txt` | Full run log |

---

## Phase 3: K-means Clustering — DONE (2026-04-27)

**Status:** Two runs preserved on disk for full reproducibility.
- **v1** (no ratio bounds, k=6) — diagnostic; produced a degenerate 2-loop cluster driven by 2 statistical outliers. Outputs at `cluster/bap1_late/cluster3/{elbow_plot,k-6}_v1_no-ratio-bound/`. Log: `cluster/phase3.txt`.
- **v2** (with symmetric ratio bounds `0.333 ≤ mut/ctrl ≤ 3.0`, k=6) — **canonical result**, used by Phase 4+. Outputs at `cluster/bap1_late/cluster3/{elbow_plot,k-6}/`. Log: `cluster/phase3_v2.txt`.

Six well-populated clusters in v2 (12,298 / 10,970 / 8,738 / 3,916 / 667 / 2,359), no degeneracy, biologically clean — captures 99.3% of edgeR's down_in_mutant calls and 99.5% of up_in_mutant. clust5 = strong gain (97% up_in_mutant), clust6 = strong loss (78% down_in_mutant), clust3/4 = moderate loss/gain, clust1/2 = unchanged bulk.

**Corrections applied during execution (Phase 4+ must follow these, not the original plan):**

- **Filter percentile, not Popay's `0.008` absolute.** The parent plan's `filter_threshold = 0.008` was calibrated for Popay's KR-balanced contact frequencies (~0.05 scale). Our mariner-aggregated counts are ~5,000× larger and per-resolution medians differ ~6× (5kb=171, 10kb=403, 25kb=984). Use **per-resolution 1%-tile** threshold computed empirically. Empirical thresholds: 5kb=38.0, 10kb=61.5, 25kb=56.5. Drops 394 loops (1%), keeps 38,950.
- **Cluster 3.0 binary lives at `/usr/local/bin/cluster`**, NOT `~/apps/cluster-1.59/src/cluster` as the parent plan §3 / §15 claimed. Verified via `which cluster`. PATH conflict-free with the conda env named `cluster` (envs don't add a binary by that name).
- **Symmetric mut/ctrl ratio bound is necessary** to prevent k-means degeneracy. v1 (no bounds) put one centroid at ratio≈7.3 to capture 2 extreme outlier loops (`chr8:69440000` at ratio=4.6, `chr7:104225000` at ratio=10.08), wasting one slot on 2 loops while leaving the loss side under-resolved. Both outliers are flagged `unchanged` by edgeR (FDR>0.25) — replicate-noise driven, not real biology. Adding `--max-ratio 3.0 --min-ratio 0.333` drops these 2 loops and lets k-means produce a balanced 6-cluster solution.
- **`sort_clusters` orders by raw signal mean, NOT mut/ctrl ratio.** It calls `clusters_df[data_cols].to_numpy().mean()` on the un-normalized values. So `clust1` has the highest mean ctrl_merge+mut_merge across loops (median ctrl=555 in v2), not the highest mut/ctrl ratio. The directional information lives in the `direction` column of the metadata sidecar, not in the cluster index.
- **Multi-format figure output.** Mirroring `scripts/utils/multi_format_output.R` (Illustrator + slides + publication workflow), each figure is saved in its own subfolder as `{name}.{png,pdf,svg,jpg}`. Implemented as a Python context-manager `multi_format_savefig()` that monkey-patches `plt.savefig` so Popay's plotting modules emit all four formats with no edits to her code. Utility lives at `cluster/scripts/utils/multi_format_output.py`.
- **`Y_range=None` (auto-fit) for line/strip/box, not Popay's hardcoded `[0,1]` / `[0,0.4]`.** Her ranges were tuned for NIPBL depletion (loop-loss-only); our BAP1-KO data spans 0.39–2.31 (post-bound) with both gains and losses, so hardcoded ranges would clip the gain side.
- **Python 3.8 in the conda env requires `Optional[float]` instead of `float | None`** (PEP 604 syntax requires 3.10+). Use `from typing import Optional`.
- **Cluster 3.0 default jobname behavior.** With `-f input_matrix.txt`, output `.kgg` lands at `<input_dir>/input_matrix_K_G{k}.kgg` (same dir as input). The script reads it from there.

### 3.1 — Elbow plot

```python
loop_count_file = 'cluster/data/late_merged_loop_counts.txt'
normalize_col = 'ctrl_merge'             # NOT 'ctrl_merged' — see Phase 1 corrections
filter_pct    = 0.01                      # per-resolution percentile, NOT absolute 0.008
min_ratio     = 0.333                     # drop mut/ctrl < 0.333 (symmetric)
max_ratio     = 3.0                       # drop mut/ctrl > 3.0 (removes 2 outliers)
```

Load TSV → join metadata (for `resolution_kb`) → filter per-resolution percentile → apply ratio bounds → normalize all data columns by `ctrl_merge` → `cluster_tools.elbow()` for k=1..14.

**Output:** `cluster/bap1_late/cluster3/elbow_plot/elbow_plot.{png,pdf,svg,jpg}` (subfolder due to multi-format).

**v2 elbow shape:** k=1 SSE=550, k=4 SSE≈100, k=6 SSE≈50, k=8 SSE≈30, asymptote ≈ 10. Bend visible around k=4-5; k=6 sits at inflection point and matches Popay precedent. Defensible.

### 3.2 — Run Cluster 3.0 k-means

Logic (from HiC_cluster3 cell-5):
1. Load + filter + normalize (as above)
2. Build `id = chr1-x1-x2-chr2-y1-y2`
3. Save normalized matrix to `input_matrix.txt`
4. Shell out: `/usr/local/bin/cluster -f input_matrix.txt -r 100 -g 7 -k {k}`
5. Read `input_matrix_K_G{k}.kgg`, merge with raw-counts df on `id`
6. `cluster_tools.sort_clusters()` → assigns `clust1..clustK` ordered by descending raw signal mean (NOT ratio — see corrections)
7. Save to `combined-clusters.txt`

**Output:** `cluster/bap1_late/cluster3/k-{k}/data/combined-clusters.txt`
Format: `GROUP  chr1  x1  x2  chr2  y1  y2  ctrl_merge  mut_merge` (9 cols, header row, ~38,948 rows in v2)

### 3.3 — Cluster visualization

Logic (from HiC_cluster3 cell-9): max-normalize per row, melt to long format, generate heatmap/lineplot/boxplot/stripplot. Our data is 2 conditions, so `comparison_type(['ctrl_merge','mut_merge']) → ('multiple', ['ctrl','mut'])`.

Palette covers both pre- and post-`_merge`-strip keys:
```python
palette = {
    'ctrl_merge': 'darkgrey', 'mut_merge': 'forestgreen',  # for heat() — uses raw data_cols
    'ctrl':       'darkgrey', 'mut':       'forestgreen',  # for box/strip — _merge stripped
}
```

Each plot wrapped in `multi_format_savefig()` and `figure_subfolder()` so output lands as:
- `cluster/bap1_late/cluster3/k-{k}/figures/heatmap/heatmap.{png,pdf,svg,jpg}`
- `cluster/bap1_late/cluster3/k-{k}/figures/lineplot/lineplot.{png,pdf,svg,jpg}`
- `cluster/bap1_late/cluster3/k-{k}/figures/stripplot/stripplot.{png,pdf,svg,jpg}`
- `cluster/bap1_late/cluster3/k-{k}/figures/boxplot/boxplot.{png,pdf,svg,jpg}`

### Phase 3 Verification

```bash
LOG=cluster/phase3_v2.txt bash cluster/scripts/run_phase3.sh 6 0.01 0.333 3.0
```

**Verified 2026-04-27 (v2 — canonical):**

- `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt`: 38,949 lines (38,948 rows + header). All 9 expected columns. Coordinates inner-join cleanly to Phase 1 metadata sidecar.
- **Cluster sizes (no degenerate cluster):** clust1=12,298 (31.6%) / clust2=10,970 (28.2%) / clust3=8,738 (22.4%) / clust4=3,916 (10.1%) / clust5=667 (1.7%) / clust6=2,359 (6.1%). Smallest cluster (clust5) has 667 loops — sufficient for stable ChromHMM enrichment in Phase 4.4.
- **Per-cluster mut/ctrl ratio (median):** clust1=0.93 / clust2=1.01 / clust3=0.86 / clust4=1.12 / clust5=1.34 / clust6=0.76. Range bounded by `[0.333, 3.0]` filter — max observed 2.31 in clust5.
- **Per-cluster median logFC (edgeR):** clust1=-0.03 / clust2=+0.09 / clust3=-0.14 / clust4=+0.25 / clust5=+0.50 / clust6=-0.31.
- **cluster × direction crosstab (row %):**
  | cluster | down_in_mutant | unchanged | up_in_mutant |
  |---------|---------------|-----------|--------------|
  | clust1  | 0.0           | **100.0** | 0.0          |
  | clust2  | 0.0           | 92.2      | 7.8          |
  | clust3  | **21.2**      | 78.8      | 0.0          |
  | clust4  | 0.0           | 30.2      | **69.8**     |
  | clust5  | 0.0           | 2.7       | **97.3**     |
  | clust6  | **78.5**      | 21.5      | 0.0          |
- **Differential capture:** 3,703 of 3,728 down_in_mutant loops (99.3%) and 4,232 of 4,253 up_in_mutant (99.5%) are present in v2 (the 2 outliers and the 1%-tile-filter drops account for the rest).
- **Resolution distribution per cluster:** all 6 clusters draw from all 3 resolutions in roughly proportional ratios (no cluster monopolizes one resolution). E.g., clust5: 5kb=112 / 10kb=379 / 25kb=176; clust1: 5kb=2,274 / 10kb=4,258 / 25kb=5,766.
- **Multi-format outputs:** all 4 figure types × 4 formats present in their respective subfolders. Heatmap.svg ≈ 12 MB, stripplot.svg ≈ 10 MB (one vector element per data point × 38k loops). PDF/PNG/JPG sizes 60–110 KB.
- **Elbow plot v2:** SSE drops from 550 (k=1) to 100 (k=4) to 50 (k=6) to 30 (k=8); inflection at k=4-5 with continued gradual improvement. k=6 defensible.
- **Multi-format wrapper smoke test:** `multi_format_savefig()` correctly emits 4 sibling files per figure, restores `plt.savefig` on context exit. Verified before integrating.

**v1 vs v2 comparison (preserved on disk):**

| Metric                  | v1 (no bounds)              | v2 (bounds applied)              |
|-------------------------|-----------------------------|-----------------------------------|
| Loops clustered         | 38,950                      | 38,948 (−2 outliers)              |
| clust1 size             | **2** (degenerate)          | 12,298 (well-populated)           |
| All 6 clusters > 100    | No (clust1 = 2)             | Yes (smallest = 667)              |
| Smallest cluster        | 2                           | 667                               |
| Median ratio range      | 0.80–7.34                   | 0.76–1.34                         |
| Differential capture    | 99.3% / 99.5%               | 99.3% / 99.5%                     |
| Output dir              | `*_v1_no-ratio-bound/`      | `elbow_plot/`, `k-6/` (canonical) |
| Log file                | `cluster/phase3.txt`        | `cluster/phase3_v2.txt`           |

**Files created in Phase 3:**

| File | Purpose |
|------|---------|
| `cluster/scripts/utils/multi_format_output.py` | Context-manager savefig wrapper + figure_subfolder helper (Python equivalent of `scripts/utils/multi_format_output.R`); reused by Phase 4-5 |
| `cluster/scripts/04_clustering.py` | Elbow + Cluster 3.0 k-means + sort + 4 visualizations; CLI flags `--elbow-only --k --filter-pct --min-ratio --max-ratio` |
| `cluster/scripts/run_phase3.sh` | Driver: positional args `K FILTER_PCT MIN_RATIO MAX_RATIO`, env var `LOG` for log path, `tee` to `$LOG` |
| `cluster/phase3.txt` | v1 run log (no ratio bounds — diagnostic) |
| `cluster/phase3_v2.txt` | v2 run log (with bounds — canonical) |
| `cluster/bap1_late/cluster3/elbow_plot_v1_no-ratio-bound/` | v1 elbow plot (preserved) |
| `cluster/bap1_late/cluster3/elbow_plot/` | v2 elbow plot (canonical) |
| `cluster/bap1_late/cluster3/k-6_v1_no-ratio-bound/` | v1 cluster outputs (preserved) |
| `cluster/bap1_late/cluster3/k-6/` | v2 cluster outputs (canonical, used by Phase 4+) |

---
