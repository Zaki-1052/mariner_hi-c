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

## Phase 2: ChromHMM Segmentation

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
# Expect ~14M lines (mm10 / 200bp)
```

---

## Phase 3: K-means Clustering

**Script:** `cluster/scripts/04_clustering.py` (adapted from HiC_cluster3.ipynb)

Cluster 3.0 is installed on the user's Mac and available on PATH as `cluster`. Use it faithfully as Popay does — subprocess call with `-r 100 -g 7 -k {k}`.

### 3.1 — Elbow plot

```python
loop_count_file = 'cluster/data/late_merged_loop_counts.txt'
normalize_col = 'ctrl_merged'
filter_threshold = 0.008
```

Load TSV -> filter rows where ctrl_merged > threshold -> normalize all data columns by ctrl_merged -> `cluster_tools.elbow()` for k=1-14.

**Output:** `cluster/bap1_late/cluster3/elbow_plot.png`

**Decision point:** Inspect elbow to choose k. Start with k=6 (matching Popay).

### 3.2 — Run Cluster 3.0 k-means

Logic (from HiC_cluster3 cell-5):
1. Normalize to ctrl_merged, filter as above
2. Create `id` column: `chr1-x1-x2-chr2-y1-y2`
3. Save normalized matrix to `input_matrix.txt`
4. Shell out to Cluster 3.0: `cluster -f input_matrix.txt -r 100 -g 7 -k {k}`
5. Read `.kgg` output, merge with original coordinates
6. Sort clusters by descending mean signal via `cluster_tools.sort_clusters()`
7. Save to `combined-clusters.txt`

**Output:** `cluster/bap1_late/cluster3/k-{k}/data/combined-clusters.txt`
Format: `GROUP  chr1  x1  x2  chr2  y1  y2  ctrl_merged  mut_merged`

### 3.3 — Cluster visualization

Logic (from HiC_cluster3 cell-9): max-normalize per row, melt to long format, generate heatmap/lineplot/boxplot/stripplot.

Our data has 2 conditions (not 3 like Popay's DMSO/4hr/24hr). `comparison_type()` returns `'multiple'` since column names have no `_` separator creating pairwise structure. Palette:
```python
palette = {'ctrl_merged': 'darkgrey', 'mut_merged': 'forestgreen'}
```

**Outputs:** `cluster/bap1_late/cluster3/k-{k}/figures/{heatmap,lineplot,boxplot,stripplot}.{png,pdf}`

### Phase 3 Verification
```bash
wc -l cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt
# Expect ~39,344 + 1 header
cut -f1 cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt | sort | uniq -c
# Check 6 clusters with reasonable distribution (no cluster < 1%)
```

---

## Phase 4: Downstream Analyses

**Script:** `cluster/scripts/05_grouped_analyses.py` (adapted from grouped_loops_figures.ipynb)

All analyses use `group_dict` built from `combined-clusters.txt` split by GROUP column.

### 4.1 — Loop size

**Reuse:** `bedpe_analysis.loop_size(out_dir, bedpe_dict=group_dict)` — no changes needed.
**Output:** `cluster/bap1_late/figures/loop_size.{stats.txt,png,pdf}`

### 4.2 — Loop classification (CTCF-only, no RAD21)

Popay's classification requires CTCF+RAD21 at both anchors for "structural". We have no RAD21 data, so use CTCF alone as the structural marker.

Modify the classifier_dict to omit RAD21 and adjust the classification rules:
- **structural:** CTCF at both anchors, no enhancer/promoter at either (`CTCF==2 & EorP<2`)
- **CRE:** enhancer or promoter at both anchors, no structural marker (`CTCF<2 & EorP==2`)
- **mixed:** one anchor CTCF, other anchor CRE (reflected by remaining combinations)
- **unclassified:** everything else

The bioframe overlap logic stays identical — just remove RAD21 from `classifier_dict` and replace `(loop_df['CTCF'] == 2) & (loop_df['RAD21'] == 2)` with `(loop_df['CTCF'] == 2)`.

**Input files:**
- CTCF: `peaks/CTCF.bed` (32,487 peaks)
- Enhancer: `peaks/beds/H3K27acCerebellumLate2.bed` (15,105 peaks)
- Promoter: `cluster/data/mm10_knownGene_pp.bed` (Phase 1 output)

**Output:** `cluster/bap1_late/figures/loop_classification.{png,pdf}`

### 4.3 — ChIP signal in anchors (raw RPKM, no RAD21 normalization)

Popay normalizes all ChIP RPKM to RAD21 signal at RAD21 peaks overlapping the anchors. Since we lack RAD21:

**Approach:** Extract raw RPKM signal directly at loop anchors using pyBigWig (mean signal per anchor, then mean of both anchors per loop). Skip the bed_intersect filtering step entirely — compute signal at ALL anchor regions, not just those overlapping a reference peak. This is biologically valid: we're asking "what is the histone mark level at loop anchors per cluster?" rather than "what is the mark level at cohesin-bound loop anchors?". The structural information is already captured by the CTCF-based classification in 4.2.

Per-mark boxplots across clusters (Kruskal-Wallis + pairwise Wilcoxon), same statistical framework as Popay.

**Input BigWigs (verified on disk at `peaks/bigwigs/macs2.narrow.aug18.dedup/`):**
- ctrl H3K27ac: `index_13_ctrl_1_H3K27ac_S25_L001_aligned_reads.sorted.bw`
- ctrl H3K27me3: `index_25_ctrl_1_H3K27me3_S37_L001_aligned_reads.sorted.bw`
- mut H3K27ac: `index_19_mut_1_H3K27ac_S46_L002_R2_001.fastq.gz_aligned_reads.sorted.bw`
- mut H3K27me3: `index_23_mut_1_H3K27me3_S50_L002_R1_001_aligned_reads.sorted.bw`

**Additional BigWigs (user has these locally at `/Users/zakiralibhai/sdsc/bigwigs/`, ask to sync):**
- H2AK119ubCtrl.bw, H2AK119ubMut.bw — enables K119ub anchor signal per cluster
- H3K27me1Ctrl.bw, H3K27me1Mut.bw — H3K27me1 (not in macs2 bigwigs)

**Output:** `cluster/bap1_late/figures/ChIP_intersect/anchor_ChIP_box.{png,pdf}`, `anchor_ChIP.stats.txt`

### 4.4 — ChromHMM anchor vs span enrichment (KEY — Fig 2f)

**Inputs:**
- Segmentation: `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` (Phase 2)
- Rename: `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` (Phase 2)
- Loop clusters: `combined-clusters.txt` (Phase 3)

**Logic (from notebook cell-23):**

1. For each cluster, write two BED files:
   - **Anchor BED:** Both anchors concatenated — `[chr1,x1,x2]` + `[chr2,y1,y2]` renamed to 3-col
   - **Span BED:** Full loop extent — `[chr1, x1, y2]` (matching notebook cell-23: `df[['chr1','x1','y2']]`)

2. Write `coordlistfile.txt` listing BED filenames

3. Run ChromHMM OverlapEnrichment for anchors and spans:
```bash
java -mx4G -jar cluster/ChromHMM/ChromHMM.jar OverlapEnrichment \
    -noimage -uniformscale \
    -m cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt \
    -f cluster/bap1_late/chromHMM/coordlistfile.txt \
    -colfields 0,1,2 \
    cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed \
    cluster/bap1_late/chromHMM/anchor_input/ \
    cluster/bap1_late/chromHMM/anchor
```
   (Repeat with `span_input/` and `span`)

4. Generate heatmaps via `chromHMM_heatmap.heatmap_plot(path=..., normalize=False)`

**Outputs:**
- `cluster/bap1_late/chromHMM/{anchor,span}.txt` — fold enrichment tables
- `cluster/bap1_late/chromHMM/{anchor,span}.{png,pdf}` — heatmaps

### 4.5 — ChromHMM proportions stacked bar

Logic (from notebook cell-25): bioframe overlap of each loop anchor with segmentation, find most prominent state per anchor by overlap length, plot proportions per cluster.

**Adaptations:**
- Update `palette` dict from Popay's 12-state names to our cerebellum names (from rename file)
- Dynamically read cluster list from `group_dict.keys()` (not hardcoded `clust1-6`)
- Map state exclusions (Popay excludes `U12`/`U9`) using our rename file

**Output:** `cluster/bap1_late/figures/chromHMM_anchor.{png,pdf}`

### 4.6 — Gene annotation

**Reuse:** `bedpe_analysis.bedtools_annotation()` with `cluster/data/mm10_knownGene_pp.bed`.

Set `FPKM=None` — Popay's pipeline optionally takes an RNA-seq FPKM DataFrame to rank/filter annotated genes by expression. We don't have RNA-seq FPKM data in the format her code expects. With `FPKM=None`, the function still produces per-cluster gene lists (via bedtools closest on the gene annotation BED) — it just skips the FPKM-based expression boxplots and the `.to_csv()` call that writes expression-ranked gene tables. The Phase 0 fix (0.8) adds a `if FPKM_df is not None:` guard around those downstream calls so the function doesn't crash when FPKM is absent.

**Output:** `cluster/bap1_late/figures/annotation/{clust1..clustN}_annotation.txt`

### 4.7 — DiffBind relationship

Logic (from notebook cells 17-19): overlaps differential ChIP peaks with loop anchors per cluster.

**Input files (verified on disk):**
- `peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt`

**Column mapping (our files vs Popay's):**
- Popay: `diffbind[['Chr','Start','End','Fold','sig']]`
- Ours: `diffbind[['Peak_Chr','Peak_Start','Peak_End','Fold','sig']]` (columns 4-6 for coordinates, col 12 for Fold, col 14 for FDR)
- Update bioframe overlap `cols2` references and drop-column names accordingly

**Outputs:**
- `cluster/bap1_late/figures/ChIP_intersect/differential_binding.{png,pdf}` (proportions per mark)
- `cluster/bap1_late/figures/ChIP_intersect/ChIP_FC_*.{png,pdf}` (fold-change boxplots)

### 4.8 — Cluster x differential status crosstab (NEW)

Not in Popay's pipeline. Uses metadata sidecar from Phase 1.

**Logic:**
1. Join `combined-clusters.txt` with `late_merged_loop_metadata.tsv` on coordinates
2. Contingency table: cluster x direction (up_in_mutant / down_in_mutant / not_significant)
3. Chi-squared test
4. Stacked bar chart via `plotting.stacked()`

**Output:** `cluster/bap1_late/figures/cluster_differential_status.{png,pdf,stats.txt}`

---

## Phase 5: deepTools Metagene Analysis

H3K27me3/H3K27ac/H2AK119ub signal across loop anchors +/-5kb per cluster.

### 5.1 — Prepare per-cluster anchor BED files

For each cluster, extract both anchors (deduplicated) into 3-column BED:
`cluster/bap1_late/figures/deeptools_input/{clust1..clustN}_anchors.bed`

### 5.2 — Run deepTools bed_pileup

Uses `deepTools_pipeline.bed_pileup()` which shells out to `computeMatrix reference-point` + `deeptools_plotting.heatmap_plot()`.

**BigWig dict — need user to confirm/sync these files:**
- H3K27me3: ctrl + mut (from `peaks/bigwigs/macs2.narrow.aug18.dedup/`)
- H3K27ac: ctrl + mut (same directory)
- H2AK119ub: ctrl + mut (user has at `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` — needs sync to `cluster/data/bigwigs/` or similar)

**Blacklist:** `tads/mm10-blacklist.v2.bed` (verified on disk)

**Output:** `cluster/bap1_late/figures/deeptools/histone_anchors.{png,pdf}`

---

## Phase 6: Cooltools Pileup (DEFERRED — needs mcool sync)

mcool files are on HPC at `/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/data/cool/250402/`. Samples: ctrl_merged.mcool, mut_merged.mcool (plus 6 individual replicates).

User can sync the merged mcools (~2-5GB each) to Mac for local cooltools pileup. Write stub script `cluster/scripts/06_cooltools_pileup.py` with the correct `cooltools_called.mcool_pileup()` call using `genome='mm10'` and the cluster BEDPE dict, with a configurable `mcool_dict` pointing to local paths.

---

## Directory Structure

```
cluster/
├── ChromHMM/                          # Existing — ChromHMM.jar + CHROMSIZES/
├── clustering_example_data/           # Existing — Popay example files
├── custom_params.json                 # Existing — matplotlib style params
├── *.py                               # Existing — Popay modules (modified in Phase 0)
│
├── scripts/                           # NEW — all runnable scripts
│   ├── 00_fix_popay_modules.sh        #   Phase 0: apply bug fixes
│   ├── 01_build_loop_count_file.py    #   Phase 1: merged loops → Popay format
│   ├── 02_build_mm10_gene_annotation.R#   Phase 1: mm10 gene BED
│   ├── 03_chromhmm_segmentation.sh    #   Phase 2: BinarizeBed + LearnModel
│   ├── 04_clustering.py               #   Phase 3: elbow + Cluster 3.0 + viz
│   ├── 05_grouped_analyses.py         #   Phase 4-5: all downstream analyses
│   ├── 06_cooltools_pileup.py         #   Phase 6: cooltools pileup (stub)
│   ├── run_phase0.sh                  #   Runner: apply fixes + verify imports
│   ├── run_phase1.sh                  #   Runner: data prep
│   ├── run_phase2.sh                  #   Runner: ChromHMM segmentation
│   ├── run_phase3.sh                  #   Runner: clustering
│   ├── run_phase4.sh                  #   Runner: downstream analyses
│   ├── run_phase5.sh                  #   Runner: deepTools metagene
│   └── sync_from_hpc.sh              #   rsync helper for HPC data
│
├── data/                              # NEW — input data files
│   ├── late_merged_loop_counts.txt    #   Phase 1 output (Popay format)
│   ├── late_merged_loop_metadata.tsv  #   Phase 1 output (edgeR stats sidecar)
│   ├── mm10_knownGene_pp.bed          #   Phase 1 output (gene annotation)
│   └── bigwigs/                       #   Synced BigWigs for deepTools
│       ├── H2AK119ubCtrl.bw           #   (synced from HPC via sync_from_hpc.sh)
│       └── H2AK119ubMut.bw
│
└── bap1_late/                         # NEW — all analysis outputs
    ├── cluster3/
    │   ├── elbow_plot.png
    │   └── k-{k}/
    │       ├── data/
    │       │   ├── input_matrix.txt
    │       │   └── combined-clusters.txt
    │       └── figures/
    │           ├── heatmap.{png,pdf}
    │           ├── lineplot.{png,pdf}
    │           ├── boxplot.{png,pdf}
    │           └── stripplot.{png,pdf}
    ├── figures/
    │   ├── annotation/                # Gene lists per cluster
    │   ├── ChIP_intersect/            # DiffBind + ChIP signal plots
    │   ├── deeptools/                 # Metagene heatmaps
    │   ├── deeptools_input/           # Per-cluster anchor BEDs
    │   ├── loop_size.{png,pdf}
    │   ├── loop_classification.{png,pdf}
    │   ├── chromHMM_anchor.{png,pdf}
    │   └── cluster_differential_status.{png,pdf}
    ├── chromHMM/
    │   ├── peak_beds/                 # 3-col BEDs staged for BinarizeBed
    │   ├── cellmarkfiletable.txt
    │   ├── binarized/                 # BinarizeBed output
    │   ├── learned_model/             # LearnModel output
    │   │   ├── cerebellum_late_12_segments.bed
    │   │   └── emissions_12.txt
    │   ├── 12state_rename_cerebellum.txt  # Manual state naming
    │   ├── anchor_input/              # Per-cluster anchor BEDs for OverlapEnrichment
    │   ├── span_input/                # Per-cluster span BEDs for OverlapEnrichment
    │   ├── {anchor,span}.txt          # Fold enrichment tables
    │   └── {anchor,span}.{png,pdf}    # Heatmaps (Fig 2f equivalent)
    └── cooltools/                     # Phase 6 outputs (after mcool sync)
```

---

## Runner Scripts

Each phase gets a `run_phase{N}.sh` bash script in `cluster/scripts/` that runs the steps for that phase sequentially with error checking. These are meant to be run manually one phase at a time (not chained).

### `run_phase0.sh` — Apply bug fixes + verify

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Bug fixes are applied via edits (not a script), but this verifies they work:
python3 -c "
import plotting
import cluster_tools
import chromHMM_heatmap
import deeptools_plotting
import bedpe_analysis
print('All imports OK')
"
```

### `run_phase1.sh` — Data preparation

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data bap1_late/{cluster3,figures/{annotation,ChIP_intersect,deeptools,deeptools_input},chromHMM/{binarized,learned_model,anchor_input,span_input,peak_beds},cooltools}
python3 scripts/01_build_loop_count_file.py
Rscript scripts/02_build_mm10_gene_annotation.R
echo "Phase 1 complete. Verify:"
wc -l data/late_merged_loop_counts.txt data/mm10_knownGene_pp.bed
```

### `run_phase2.sh` — ChromHMM segmentation

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/03_chromhmm_segmentation.sh
echo "Phase 2 complete. Verify:"
wc -l bap1_late/chromHMM/learned_model/*_12_segments.bed
echo "Now inspect emissions_12.txt and create 12state_rename_cerebellum.txt"
```

### `run_phase3.sh` — Clustering

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/04_clustering.py
echo "Phase 3 complete. Verify:"
wc -l bap1_late/cluster3/k-*/data/combined-clusters.txt
```

### `run_phase4.sh` — Downstream analyses

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/05_grouped_analyses.py
echo "Phase 4 complete. Check bap1_late/figures/ and bap1_late/chromHMM/"
```

### `run_phase5.sh` — deepTools metagene

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Requires BigWigs in data/bigwigs/ — run sync_from_hpc.sh first if needed
python3 scripts/05_grouped_analyses.py --deeptools-only
echo "Phase 5 complete. Check bap1_late/figures/deeptools/"
```

### `sync_from_hpc.sh` — rsync HPC data

```bash
#!/bin/bash
set -euo pipefail
# Sync BigWigs and mcools from SDSC Expanse (ssh alias: expanse)
# Run from repo root: bash cluster/scripts/sync_from_hpc.sh

REMOTE="expanse"
LOCAL_BW="cluster/data/bigwigs"
LOCAL_MCOOL="cluster/data/mcools"
mkdir -p "$LOCAL_BW" "$LOCAL_MCOOL"

echo "=== Syncing BigWigs ==="
# H2AK119ub (not in peaks/bigwigs/)
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H2AK119ubCtrl.bw" "$LOCAL_BW/"
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H2AK119ubMut.bw" "$LOCAL_BW/"
# H3K27me1
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H3K27me1Ctrl.bw" "$LOCAL_BW/"
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H3K27me1Mut.bw" "$LOCAL_BW/"

echo "=== Syncing mcools (Phase 6 only — large files, ~2-5GB each) ==="
echo "Uncomment below when ready for cooltools pileup:"
# rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/data/cool/250402/ctrl_merged.mcool" "$LOCAL_MCOOL/"
# rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/data/cool/250402/mut_merged.mcool" "$LOCAL_MCOOL/"

echo "Done. Files synced to:"
ls -lh "$LOCAL_BW"/ 2>/dev/null || echo "  (no bigwigs yet)"
ls -lh "$LOCAL_MCOOL"/ 2>/dev/null || echo "  (no mcools yet)"
```

**Note on BigWig paths in `sync_from_hpc.sh`:** The remote paths above assume `/expanse/lustre/projects/csd940/zalibhai/bigwigs/` matches what's on the user's Mac at `/Users/zakiralibhai/sdsc/bigwigs/`. If the HPC paths differ, update the script before running. The H3K27ac/H3K27me3 BigWigs are already on disk at `peaks/bigwigs/macs2.narrow.aug18.dedup/` and don't need syncing.

---

## Files to Create

| File | Phase | Description |
|------|-------|-------------|
| **Analysis scripts** | | |
| `cluster/scripts/01_build_loop_count_file.py` | 1.1 | Merged loops → Popay format count file |
| `cluster/scripts/02_build_mm10_gene_annotation.R` | 1.2 | mm10 promoter BED from TxDb |
| `cluster/scripts/03_chromhmm_segmentation.sh` | 2 | BinarizeBed + LearnModel pipeline |
| `cluster/scripts/04_clustering.py` | 3 | Elbow plot + Cluster 3.0 k-means + visualization |
| `cluster/scripts/05_grouped_analyses.py` | 4-5 | All downstream analyses (ChromHMM, ChIP, DiffBind, deepTools) |
| `cluster/scripts/06_cooltools_pileup.py` | 6 | Cooltools pileup stub (after mcool sync) |
| **Runner scripts** | | |
| `cluster/scripts/run_phase0.sh` | 0 | Verify bug fixes |
| `cluster/scripts/run_phase1.sh` | 1 | Data prep + directory setup |
| `cluster/scripts/run_phase2.sh` | 2 | ChromHMM segmentation |
| `cluster/scripts/run_phase3.sh` | 3 | Clustering |
| `cluster/scripts/run_phase4.sh` | 4 | Downstream analyses |
| `cluster/scripts/run_phase5.sh` | 5 | deepTools metagene |
| `cluster/scripts/sync_from_hpc.sh` | — | rsync BigWigs + mcools from Expanse |

## Files to Modify

| File | Changes |
|------|---------|
| `cluster/plotting.py` | Add `import os`; fix custom_params path (L24-25); remove gene ylim overrides (L132-139); fix joint() logX bug (L502: `ycol`→`xcol`) |
| `cluster/cooltools_called.py` | Add `import os`; replace `plt.style.use('seaborn-poster')` with `'seaborn-v0_8-poster'` (L16); fix custom_params path (L18-19); change default genome to mm10 (L27) |
| `cluster/chromHMM_heatmap.py` | Fix custom_params path (L9-10) |
| `cluster/deeptools_plotting.py` | Add `import os`; fix custom_params path (L8-9); remove HA-NIPBL ylim override (L108) |
| `cluster/cluster_tools.py` | Fix elbow() savefig/show order (L18-19) |
| `cluster/deepTools_pipeline.py` | Add genome param to bam_coverage() with mm10 size (L15) |
| `cluster/bedpe_analysis.py` | Fix default gene path (L67); add FPKM_df None guard (L142-159) |

## Verification Plan

1. **Phase 0:** `run_phase0.sh` — all module imports succeed without error
2. **Phase 1:** `late_merged_loop_counts.txt` has ~39,344 rows, 8 columns; `mm10_knownGene_pp.bed` has ~24K entries
3. **Phase 2:** Segmentation BED covers all autosomes + chrX; emissions matrix shows interpretable mark combinations
4. **Phase 3:** Elbow plot shows clear bend; combined-clusters.txt has ~39K rows, no cluster < 1%
5. **Phase 4.4 (key result):** ChromHMM anchor vs span heatmaps show differential Polycomb enrichment across clusters
6. **Phase 4.8:** Chi-squared p < 0.05 for cluster x differential status
7. **Phase 5:** deepTools heatmaps show H3K27me3/H3K27ac differences across clusters at anchors