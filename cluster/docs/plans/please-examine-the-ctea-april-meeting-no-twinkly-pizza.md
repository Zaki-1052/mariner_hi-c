# Plan: Adapt Popay Hi-C Loop Clustering Pipeline for BAP1-KO Cerebellum

## Context

The Dixon meeting (2026-04-10) identified the Popay anchor-vs-span ChromHMM analysis as the highest-priority outstanding item for the BAP1-KO paper. Tessa Popay shared her code (2026-04-21), which is now in `cluster/`. The goal is to reproduce Popay's Figure 2f — ChromHMM chromatin state enrichment at loop anchors vs. loop spans across k-means clusters — adapted to our BAP1-KO vs wildtype mouse cerebellum data. This addresses the central biological question: is Polycomb enrichment at the anchors themselves (sensitivity model) or across the loop body (extrusion impediment model)?

**Scope:** Late/adult timepoint (250402) first. All ~18K loops from 5kb count matrix. Replicates merged to 2 columns (ctrl_merged, mut_merged). ChromHMM segmentation generated from our ChIP-seq peaks.

---

## Phase 0: Bug Fixes in Popay Pipeline Modules

Four Python modules crash on import due to hardcoded paths and a removed matplotlib style. These must be fixed before any downstream code runs.

### 0.1 — Fix hardcoded `custom_params.json` path (4 files)

All four files load `/Users/tessapopay/example_data/custom_params.json` at module import time. Replace with path relative to module location. The local copy already exists at `cluster/custom_params.json`.

| File | Line | Fix |
|------|------|-----|
| `cluster/plotting.py` | 24-25 | Use `os.path.join(os.path.dirname(os.path.abspath(__file__)), 'custom_params.json')` |
| `cluster/cooltools_called.py` | 18-19 | Same pattern |
| `cluster/chromHMM_heatmap.py` | 9-10 | Same pattern |
| `cluster/deeptools_plotting.py` | 8-9 | Same pattern |

### 0.2 — Fix `plt.style.use('seaborn-poster')` crash

**File:** `cluster/cooltools_called.py:16`
**Bug:** `'seaborn-poster'` was removed in matplotlib 3.6. Our env has 3.7.5.
**Fix:** Delete line 16. The custom_params via `sns.set_theme()` already controls styling.

### 0.3 — Fix `elbow()` blank PNG

**File:** `cluster/cluster_tools.py:18-19`
**Bug:** `plt.show()` clears the figure before `plt.savefig()`.
**Fix:** Swap the two lines — `plt.savefig()` first, then `plt.show()`.

### 0.4 — Fix `joint()` logX zero-check on wrong column

**File:** `cluster/plotting.py` (in `joint()` function body, the `if logX:` block)
**Bug:** The zero-value filter checks `data[ycol]` instead of `data[xcol]`.
**Fix:** Change the column reference from `ycol` to `xcol` in the `if logX:` block.

### 0.5 — Remove hardcoded biology overrides

| File | Location | What to remove |
|------|----------|---------------|
| `cluster/plotting.py` | `line()` function | 8 `if dataset == 'TFAP2A'`-style ylim overrides (Popay-specific genes) |
| `cluster/deeptools_plotting.py` | `heatmap_plot()` | `if bam == 'HA-NIPBL': line.set_ylim([0,2])` |

### 0.6 — Fix hg38 genome size in deepTools_pipeline.py

**File:** `cluster/deepTools_pipeline.py:15`
**Bug:** `--effectiveGenomeSize 2913022398` hardcoded (hg38 value).
**Fix:** Add `genome` parameter to `bam_coverage()`, default `'mm10'`, with dict lookup: `{'mm10': 2494787188, 'hg38': 2913022398}`.

### 0.7 — Fix cooltools default genome

**File:** `cluster/cooltools_called.py`, `mcool_pileup()` signature
**Fix:** Change default `genome='hg38'` to `genome='mm10'`.

### Verification

```bash
conda activate cluster
cd /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster
python -c "import plotting; import cluster_tools; import chromHMM_heatmap; import deeptools_plotting; print('All imports OK')"
# cooltools_called requires cooler/cooltools — test separately:
python -c "import cooltools_called; print('cooltools OK')"
```

---

## Phase 1: Data Preparation

### 1.1 — Build loop count file in Popay format

**Script:** `cluster/scripts/01_build_loop_count_file.py`

**Inputs:**
- `outputs/250402-late_outputs/res_5kb/06_counts_matrix.tsv` (17,983 loops x 6 replicate columns)
- `outputs/250402-late_outputs/edgeR_results_res_5kb/primary_analysis/all_results_primary.tsv` (same loops with chr1/start1/end1/chr2/start2/end2 coordinates)

**Logic:**
1. Load count matrix (loop_id index, ctrl_M1/M2/M3, mut_M1/M2/M3)
2. Load edgeR results (loop_id, chr1, start1, end1, chr2, start2, end2, logFC, FDR, direction, ...)
3. Inner join on loop_id
4. Average replicates: `ctrl_merged = mean(ctrl_M1, ctrl_M2, ctrl_M3)`, `mut_merged = mean(mut_M1, mut_M2, mut_M3)`
5. Output columns: `chr1  x1  x2  chr2  y1  y2  ctrl_merged  mut_merged` (rename start1->x1, end1->x2, start2->y1, end2->y2)
6. Also output a metadata sidecar preserving the edgeR stats (logFC, FDR, direction, category) for each loop_id, keyed by coordinate string — used in Phase 4 to annotate clusters with differential status

**Output:** `cluster/data/late_5kb_loop_counts.txt` (17,983 rows, 8 columns, tab-separated with header)

### 1.2 — Generate mm10 gene annotation BED

**Script:** `cluster/scripts/02_build_mm10_gene_annotation.R`

**Goal:** Produce mm10 equivalent of Popay's `gencode.v25.annotation_pp.bed` (hg38 Gencode promoter-proximal BED, 7 columns: chr start end gene_id score strand gene_name, TSS ± 750bp).

**Logic:**
1. Load `TxDb.Mmusculus.UCSC.mm10.knownGene` and `org.Mm.eg.db` (both already installed in mariner_env)
2. Extract gene coordinates, compute promoters (upstream=750, downstream=750)
3. Map Entrez IDs to gene symbols via `org.Mm.eg.db`
4. Output 7-column BED matching Popay's format

**Output:** `cluster/data/mm10_knownGene_pp.bed`

Also update `bedpe_analysis.py` default gene path to point here (or pass explicitly in notebooks).

### 1.3 — Create output directory structure

```
cluster/
  data/                          # Phase 1 outputs
    late_5kb_loop_counts.txt
    late_5kb_loop_metadata.tsv
    mm10_knownGene_pp.bed
  bap1_late/                     # Working directory for late timepoint
    cluster3/                    # Phase 3 outputs (elbow, k-means, figures)
    chromHMM/                    # Phase 2 + 4 outputs
      binarized/
      learned_model/
      anchor_input/
      span_input/
    figures/                     # Phase 4-5 outputs
      annotation/
      ChIP_intersect/
      deeptools/
    cooltools/                   # Phase 6 outputs (HPC)
  scripts/                      # New scripts
```

---

## Phase 2: ChromHMM Segmentation Generation

**Goal:** Produce a full genome-wide mm10 cerebellum chromatin state segmentation in 4-column BED format (chr start end state_label) for use with ChromHMM OverlapEnrichment.

### 2.1 — Prepare peak BEDs for BinarizeBed

**Problem:** Our peak BEDs are split across two directories (`peaks/beds/` and `peaks/`), and BinarizeBed needs all BEDs in one `inputbamdir`. Also, our BEDs have 6+ columns but BinarizeBed with `-peaks` only reads chr/start/end.

**Script:** `cluster/scripts/03_chromhmm_segmentation.sh`

**Step 2.1a:** Create a staging directory with 3-column BED copies (or symlinks):
```
cluster/bap1_late/chromHMM/peak_beds/
  H3K27ac.bed      <- cut -f1-3 peaks/beds/H3K27acCerebellumLate2.bed
  H3K27me3.bed     <- cut -f1-3 peaks/beds/H3K27me3CerebellumLate1.bed
  H3K4me1.bed      <- cut -f1-3 peaks/beds/H3K4me1CerebellumLate1.bed
  H3K4me3.bed      <- cut -f1-3 peaks/beds/H3K4me3CerebellumLate2.bed
  CTCF.bed         <- cut -f1-3 peaks/CTCF.bed
```

**Step 2.1b:** Create cellmarkfiletable (3 columns: cell, mark, filename):
```
cluster/bap1_late/chromHMM/cellmarkfiletable.txt:
cerebellum_late	H3K27ac	H3K27ac.bed
cerebellum_late	H3K27me3	H3K27me3.bed
cerebellum_late	H3K4me1	H3K4me1.bed
cerebellum_late	H3K4me3	H3K4me3.bed
cerebellum_late	CTCF	CTCF.bed
```

### 2.2 — Run ChromHMM BinarizeBed

```bash
cluster/ChromHMM/chromhmm BinarizeBed \
    -peaks \
    cluster/ChromHMM/CHROMSIZES/mm10.txt \
    cluster/bap1_late/chromHMM/peak_beds/ \
    cluster/bap1_late/chromHMM/cellmarkfiletable.txt \
    cluster/bap1_late/chromHMM/binarized/
```

The `-peaks` flag treats BEDs as peak calls directly — any 200bp bin overlapping a peak gets a 1. This produces one binary file per chromosome named `cerebellum_late_chr*_binary.txt`.

### 2.3 — Run ChromHMM LearnModel

```bash
cluster/ChromHMM/chromhmm LearnModel \
    -p 4 \
    -l cluster/ChromHMM/CHROMSIZES/mm10.txt \
    cluster/bap1_late/chromHMM/binarized/ \
    cluster/bap1_late/chromHMM/learned_model/ \
    12 \
    mm10
```

**Runtime:** ~30-90 min locally (5 marks, mm10, 200bp bins). `-p 4` uses 4 threads.

**Key outputs:**
- `cerebellum_late_12_segments.bed` — the segmentation (chr, start, end, E1-E12)
- `emissions_12.txt` — emission probability matrix (12 states x 5 marks)
- `transitions_12.txt` — state transition matrix
- `model_12.txt` — model parameters
- Enrichment heatmaps against built-in mm10 annotations (RefSeqTSS, CpGIsland, etc.)

### 2.4 — Interpret emission matrix and create rename file

**Manual step** after LearnModel completes. Inspect `emissions_12.txt` to assign biological names:

| Emission pattern | State name |
|-----------------|------------|
| High H3K4me3 + H3K27ac | Active_Promoter |
| High H3K4me1 + H3K27ac (low H3K4me3) | Active_Enhancer |
| High H3K4me1 only | Poised_Enhancer |
| High H3K27me3 only | Polycomb |
| High H3K4me3 + H3K27me3 | Bivalent |
| High CTCF | Insulator |
| Low all marks | Low_Signal |
| (other combinations) | Named by dominant marks |

Create `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` in same format as `cluster/clustering_example_data/12state_rename.txt` (2-column TSV: state_id, biological_name).

### Verification

```bash
wc -l cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed
# Should be ~14M lines (mm10 genome / 200bp bins)
head -3 cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed
# Should show: chr1  0  XXXXX  E1 (or similar)
```

---

## Phase 3: K-means Clustering

### 3.1 — Elbow plot

**Script:** `cluster/scripts/04_clustering.py` (standalone script adapted from HiC_cluster3.ipynb cells 1-3)

```python
loop_count_file = 'cluster/data/late_5kb_loop_counts.txt'
normalize_col = 'ctrl_merged'
filter_threshold = 0.008  # Popay default — will pass nearly all loops since normalized values center ~1.0
```

Logic (from HiC_cluster3 cell-3):
1. Load loop count TSV
2. Filter: keep rows where `ctrl_merged > filter_threshold`
3. Normalize: divide all data columns by `ctrl_merged`
4. Run `cluster_tools.elbow()` for k=1-14
5. Save elbow plot

**Output:** `cluster/bap1_late/cluster3/elbow_plot.png`

**Decision point:** Inspect elbow plot to choose k. Expect k=4-8 based on 2-condition data. Start with k=6 (matching Popay).

### 3.2 — Run Cluster 3.0 k-means

Logic (from HiC_cluster3 cell-5):
1. Normalize to `ctrl_merged`, filter as above
2. Create `id` column: `chr1-x1-x2-chr2-y1-y2`
3. Save normalized matrix to `input_matrix.txt`
4. Run Cluster 3.0: `~/apps/cluster-1.59/src/cluster -f input_matrix.txt -r 100 -g 7 -k {k}`
5. Read `.kgg` output, merge with original coordinates
6. Sort clusters by descending mean signal via `cluster_tools.sort_clusters()`
7. Save to `combined-clusters.txt`

**Output:** `cluster/bap1_late/cluster3/k-{k}/data/combined-clusters.txt`
Format: `GROUP  chr1  x1  x2  chr2  y1  y2  ctrl_merged  mut_merged`

### 3.3 — Cluster visualization

Logic (from HiC_cluster3 cell-9):
1. Max-normalize per row
2. Melt to long format
3. Generate heatmap, lineplot, boxplot, stripplot

**Key adaptation:** Update palette dict from Popay's `{'DMSO_merge': 'darkgrey', 'dTAG': 'forestgreen'}` to `{'ctrl_merged': 'darkgrey', 'mut_merged': 'forestgreen'}`.

**Outputs:**
- `cluster/bap1_late/cluster3/k-{k}/figures/heatmap.{png,pdf}`
- `cluster/bap1_late/cluster3/k-{k}/figures/lineplot.{png,pdf}`
- `cluster/bap1_late/cluster3/k-{k}/figures/boxplot.{png,pdf}`
- `cluster/bap1_late/cluster3/k-{k}/figures/stripplot.{png,pdf}`

### Verification

```bash
wc -l cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt
# Should be ~17,983 + 1 header = ~17,984
cut -f1 cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt | sort | uniq -c
# Should show 6 clusters with reasonable distribution
```

---

## Phase 4: Downstream Analyses

**Script:** `cluster/scripts/05_grouped_analyses.py` (standalone script adapted from grouped_loops_figures.ipynb)

All analyses use `group_dict` built from `combined-clusters.txt` split by GROUP column.

### 4.1 — Loop size analysis

**Reuse:** `bedpe_analysis.loop_size(out_dir, bedpe_dict=group_dict)`
**No changes needed** — function computes loop size as `y2 - x1` per loop.
**Output:** `cluster/bap1_late/figures/loop_size.{stats.txt,png,pdf}`

### 4.2 — Loop classification (structural/CRE/mixed)

**Adaptation needed:** No RAD21 data available. Use CTCF alone as structural marker.

**Input files:**
- CTCF: `peaks/CTCF.bed` (32,487 peaks)
- Enhancer: `peaks/beds/H3K27acCerebellumLate2.bed` (active enhancers by K27ac)
- Promoter: `cluster/data/mm10_knownGene_pp.bed` (Phase 1 output)

**Modified classification rules:**
- structural: CTCF at both anchors, no enhancer/promoter at either
- CRE: enhancer or promoter at both anchors (no CTCF requirement)
- mixed: one anchor CTCF, other anchor CRE
- unclassified: everything else

**Output:** `cluster/bap1_late/figures/loop_classification.{png,pdf}`

### 4.3 — ChIP signal in anchors

**Input BigWigs** (from `peaks/bigwigs/macs2.narrow.aug18.dedup/`):
- ctrl H3K27ac: `index_13_ctrl_1_H3K27ac_S25_L001_aligned_reads.sorted.bw`
- ctrl H3K27me3: `index_25_ctrl_1_H3K27me3_S37_L001_aligned_reads.sorted.bw`
- mut H3K27ac: `index_19_mut_1_H3K27ac_S46_L002_R2_001.fastq.gz_aligned_reads.sorted.bw`
- mut H3K27me3: `index_23_mut_1_H3K27me3_S50_L002_R1_001_aligned_reads.sorted.bw`

**Adaptation:** Popay normalizes all ChIP to RAD21. We have no RAD21, so normalize H3K27ac to H3K27me3 (or skip normalization and plot raw RPKM).

**Alternative approach:** Use existing `peaks/k119ub_anchor_signal.tsv` (ctrl_signal, mut_signal per anchor) to add K119ub signal per cluster — merge by anchor coordinates with group_dict.

**Output:** `cluster/bap1_late/figures/ChIP_intersect/anchor_ChIP_box.{png,pdf}`, `anchor_ChIP.stats.txt`

### 4.4 — ChromHMM anchor vs. span enrichment (KEY ANALYSIS)

**The Fig 2f equivalent.** This is the highest-value output.

**Inputs:**
- Segmentation: `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` (Phase 2)
- Rename: `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` (Phase 2)
- Loop clusters: `combined-clusters.txt` (Phase 3)

**Logic (from grouped_loops_figures.ipynb cell-23):**

1. For each cluster, write anchor BED (both anchors concatenated) and span BED (chr1, x2, y1 — the interval between anchors):
   ```
   cluster/bap1_late/chromHMM/anchor_input/{clust1..clustN}.bed
   cluster/bap1_late/chromHMM/span_input/{clust1..clustN}.bed
   ```

2. Write `coordlistfile.txt` listing BED filenames (one per line)

3. Run ChromHMM OverlapEnrichment for anchors:
   ```bash
   cluster/ChromHMM/chromhmm OverlapEnrichment \
       -noimage \
       -m cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt \
       -uniformscale \
       -f cluster/bap1_late/chromHMM/coordlistfile.txt \
       -colfields 0,1,2 \
       cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed \
       cluster/bap1_late/chromHMM/anchor_input/ \
       cluster/bap1_late/chromHMM/anchor
   ```

4. Repeat for spans (same command, swap `anchor_input/` → `span_input/`, `anchor` → `span`)

5. Generate heatmaps via `chromHMM_heatmap.heatmap_plot()`:
   ```python
   heatmap_plot(path='cluster/bap1_late/chromHMM/anchor.txt', normalize=False)
   heatmap_plot(path='cluster/bap1_late/chromHMM/span.txt', normalize=False)
   ```

**Outputs:**
- `cluster/bap1_late/chromHMM/anchor.txt` — fold enrichment table (states x clusters)
- `cluster/bap1_late/chromHMM/span.txt` — fold enrichment table
- `cluster/bap1_late/chromHMM/anchor.{png,pdf}` — heatmap
- `cluster/bap1_late/chromHMM/span.{png,pdf}` — heatmap

### 4.5 — ChromHMM proportions stacked bar

**Logic (from grouped_loops_figures.ipynb cell-25):**
Uses bioframe to overlap each loop anchor with the segmentation, finds the most prominent state per anchor, then plots proportions per cluster as a stacked bar chart.

**Adaptation:** Update the hardcoded `palette` dict from Popay's 12-state names to our cerebellum state names (from `12state_rename_cerebellum.txt`). Update the cluster list from `['clust1',...,'clust6']` to match our k.

**Output:** `cluster/bap1_late/figures/chromHMM_anchor.{png,pdf}`

### 4.6 — Gene annotation

**Logic (from grouped_loops_figures.ipynb cell-11):**
Uses `bedpe_analysis.bedtools_annotation()` to intersect loop anchors with gene promoters.

**Adaptation:** Pass `cluster/data/mm10_knownGene_pp.bed` as the gene annotation file. Set `FPKM=None` (no RNA-seq FPKM file in Popay format available locally — could generate from our DESeq2 results but not critical for first pass).

**Output:** `cluster/bap1_late/figures/annotation/{clust1..clustN}_annotation.txt`

### 4.7 — DiffBind relationship

**Logic (from grouped_loops_figures.ipynb cells 17-19):**
Overlaps differential ChIP peaks with loop anchors per cluster.

**Input files:**
- `peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt`

**Adaptation:** Column names differ from Popay's format. Our diffbind has `Peak_Chr, Peak_Start, Peak_End, Fold, FDR` vs Popay's `Chr, Start, End, Fold, FDR`. Update the bioframe overlap column references.

**Outputs:**
- `cluster/bap1_late/figures/ChIP_intersect/differential_binding.{png,pdf}` (per mark)
- `cluster/bap1_late/figures/ChIP_intersect/ChIP_FC_*.{png,pdf}` (fold-change boxplots)

### 4.8 — Cluster × differential status crosstab

**New analysis** (not in Popay's pipeline but essential for our story).

Using the metadata sidecar from Phase 1 (loop_id → logFC, FDR, direction), join with cluster assignments and produce:
1. Contingency table: cluster × direction (up_in_mutant / down_in_mutant / unchanged)
2. Chi-squared test for non-random distribution
3. Stacked bar chart of differential status per cluster

This directly answers: which clusters are enriched for gained vs. lost loops?

**Output:** `cluster/bap1_late/figures/cluster_differential_status.{png,pdf,stats.txt}`

---

## Phase 5: deepTools Metagene Analysis

**Goal:** H3K27me3/H3K27ac signal across loop anchors ± 5kb per cluster (the Dixon suggestion).

### 5.1 — Prepare per-cluster anchor BED files

For each cluster, extract both anchors (deduplicated) into a 3-column BED:
```
cluster/bap1_late/figures/deeptools_input/{clust1..clustN}_anchors.bed
```

### 5.2 — Run deepTools bed_pileup (local, H3K27me3 + H3K27ac)

```python
from deepTools_pipeline import bed_pileup

bigWig_dict = {
    'ctrl_H3K27me3': 'peaks/bigwigs/macs2.narrow.aug18.dedup/index_25_ctrl_1_H3K27me3_S37_L001_aligned_reads.sorted.bw',
    'mut_H3K27me3': 'peaks/bigwigs/macs2.narrow.aug18.dedup/index_23_mut_1_H3K27me3_S50_L002_R1_001_aligned_reads.sorted.bw',
    'ctrl_H3K27ac': 'peaks/bigwigs/macs2.narrow.aug18.dedup/index_13_ctrl_1_H3K27ac_S25_L001_aligned_reads.sorted.bw',
    'mut_H3K27ac': 'peaks/bigwigs/macs2.narrow.aug18.dedup/index_19_mut_1_H3K27ac_S46_L002_R2_001.fastq.gz_aligned_reads.sorted.bw',
}

bed_dict = {group: f'cluster/bap1_late/figures/deeptools_input/{group}_anchors.bed'
            for group in cluster_order}

bed_pileup(
    bed_dict=bed_dict,
    bigWig_dict=bigWig_dict,
    out_dir='cluster/bap1_late/figures/deeptools',
    blacklisted_regions='tads/mm10-blacklist.v2.bed',
    up_down=5000,
    out_name='histone_anchors'
)
```

**Output:** `cluster/bap1_late/figures/deeptools/histone_anchors.{png,pdf}`

### 5.3 — K119ub metagene (DEFERRED — requires HPC)

No K119ub BigWig files exist locally. This step requires either:
- HPC access to K119ub BAMs/BigWigs at `/expanse/lustre/projects/csd940/ctea/`
- Or generating BigWigs from BAMs on HPC and transferring

**HPC script to write later:** `cluster/scripts/06_k119ub_metagene_hpc.sh`

---

## Phase 6: Cooltools Pileup / APA per Cluster (DEFERRED — requires HPC)

mcool files are only on HPC. This produces off-diagonal APA heatmaps per cluster showing aggregate Hi-C contact enrichment.

**HPC script to write later:** `cluster/scripts/07_cooltools_pileup_hpc.sh`

Uses `cooltools_called.mcool_pileup()` with `genome='mm10'`, our mcool paths, and the cluster BEDPE dict.

---

## Files to Create

| File | Phase | Type |
|------|-------|------|
| `cluster/scripts/01_build_loop_count_file.py` | 1 | Python script |
| `cluster/scripts/02_build_mm10_gene_annotation.R` | 1 | R script |
| `cluster/scripts/03_chromhmm_segmentation.sh` | 2 | Bash script |
| `cluster/scripts/04_clustering.py` | 3 | Python script |
| `cluster/scripts/05_grouped_analyses.py` | 4-5 | Python script |

## Files to Modify (Phase 0 bug fixes)

| File | Changes |
|------|---------|
| `cluster/plotting.py` | Fix custom_params path, remove gene-specific ylim overrides, fix joint() logX bug |
| `cluster/cooltools_called.py` | Remove seaborn-poster style, fix custom_params path, change default genome to mm10 |
| `cluster/chromHMM_heatmap.py` | Fix custom_params path |
| `cluster/deeptools_plotting.py` | Fix custom_params path, remove HA-NIPBL ylim override |
| `cluster/cluster_tools.py` | Fix elbow() show-before-savefig |
| `cluster/deepTools_pipeline.py` | Add genome parameter to bam_coverage() with mm10 genome size |

## Verification Plan

1. **Phase 0:** All module imports succeed without error
2. **Phase 1:** `late_5kb_loop_counts.txt` has 17,983 rows, 8 columns; `mm10_knownGene_pp.bed` has ~24K gene entries
3. **Phase 2:** Segmentation BED covers all autosomes + chrX; `emissions_12.txt` shows interpretable mark combinations; all 12 states have non-trivial genome coverage
4. **Phase 3:** Elbow plot shows clear bend; `combined-clusters.txt` has ~17K rows with reasonable cluster size distribution (no cluster < 1% of total)
5. **Phase 4.4:** ChromHMM heatmaps show differential enrichment patterns between clusters (the key biological result — expect Polycomb-enriched clusters to be enriched for gained/strengthened loops)
6. **Phase 4.8:** Chi-squared p-value < 0.05 for cluster × differential status (clusters are not randomly distributed across loop directions)
7. **Phase 5:** deepTools heatmaps show visible H3K27me3 and H3K27ac signal differences across clusters at loop anchors
