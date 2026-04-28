# Plan: Lab-Meeting Summary Figures for BAP1-KO Loop Clustering

## Context

Phases 0-4 of the Popay-adapted clustering pipeline are complete, producing 21 figures across 8 sub-analyses. Each figure shows ONE dimension of the 6 k-means clusters (loop size, classification, ChromHMM enrichment, etc.), requiring the audience to mentally cross-reference 5+ separate plots to understand what each cluster biologically represents. The existing ChromHMM proportions plot also has a poor aspect ratio (very wide, squished text).

For the lab meeting, we need integrated composite figures where the biological story is immediately obvious: gained loops (clust5) sit in expanding Polycomb domains (anchor AND span enriched), while lost loops (clust6) have Polycomb only at anchors (sensitivity/anchor-disruption model). Two mechanisms from the same BAP1-KO perturbation.

**User preferences (confirmed via questions):**
1. X-axis: Biological gradient — strong loss -> strong gain
2. Labels: Biological names + cluster ID + n= counts
3. Layout: Dashboard composite + separate mechanism figure + feature heatmap
4. Extras: K119ub signal, gene examples, all included

## Script: `cluster/scripts/08_summary_figures.py`

New standalone script (07 is taken by cooltools_pileup.py). Reads pre-computed Phase 4 outputs plus 2 fresh BigWig queries (~24s for K119ub). Produces 3 figures.

**Driver:** `cluster/scripts/run_summary.sh` (simple tee-to-log wrapper, same pattern as run_phase4.sh)

## Data Sources (all pre-computed except K119ub)

| Data | Source file | How loaded |
|------|------------|------------|
| Cluster assignments + counts | `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt` | `pd.read_csv(sep='\t')` |
| logFC, direction | `cluster/data/late_merged_loop_metadata.tsv` | join on 6 coord cols |
| ChromHMM anchor enrichment | `cluster/bap1_late/chromHMM/anchor.txt` | `pd.read_csv(sep='\t', index_col=0)`, drop Base row |
| ChromHMM span enrichment | `cluster/bap1_late/chromHMM/span.txt` | same |
| Anchor proportions | `cluster/bap1_late/figures/chromHMM_anchor/chromHMM_anchor.proportions.tsv` | `pd.read_csv(sep='\t', index_col=0)` |
| Direction %s | `cluster/bap1_late/figures/cluster_differential_status.stats.txt` | parse "Row %:" block |
| Classification %s | `cluster/bap1_late/figures/loop_classification/loop_classification.stats.txt` | parse "Percent stacked:" block |
| Gene annotations | `cluster/bap1_late/figures/annotation/clust{1..6}_annotation.txt` | 20-col TSV, extract chr1_gene_name + chr2_gene_name |
| K119ub signal (FRESH) | `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` | pyBigWig queries on deduplicated anchors, ~24s |

## Cluster Ordering and Labels

**Biological gradient (x-axis left to right):**

| Position | Cluster | Biological label |
|----------|---------|-----------------|
| 1 (left) | clust6 | Strong loss\n(anchor-disrupted)\nclust6, n=2,359 |
| 2 | clust3 | Moderate loss\nclust3, n=8,738 |
| 3 | clust1 | Unchanged (high)\nclust1, n=12,298 |
| 4 | clust2 | Unchanged\nclust2, n=10,970 |
| 5 | clust4 | Moderate gain\nclust4, n=3,916 |
| 6 (right) | clust5 | Strong gain\n(Polycomb domain)\nclust5, n=667 |

## Figure 1: Cluster Summary Dashboard

**Output:** `cluster/bap1_late/figures/summary_figures/dashboard/cluster_dashboard.{png,pdf,svg,jpg}`

**Layout:** `figsize=(5.5, 14)`, single column, 6 panels via `GridSpec(6, 1, height_ratios=[0.8, 1.2, 1.4, 1.2, 1.2, 1.2], hspace=0.12)`. All panels share x-axis; only bottom panel shows x-tick labels. Each panel gets a letter label (A-F) in upper-left.

| Panel | What | Data source | Visual |
|-------|------|-------------|--------|
| **A** | Median logFC per cluster | metadata join -> groupby median | Diverging bar: red (#d62728) for negative, blue (#1f77b4) for positive. Horizontal line at 0. Values annotated on bars. |
| **B** | % direction (down/unchanged/up) | `cluster_differential_status.stats.txt` Row % block | Stacked horizontal bar. Colors: red/gold/green matching existing palette. |
| **C** | **KEY: Polycomb anchor vs span** | `anchor.txt` row 11_Polycomb + `span.txt` row 11_Polycomb | Grouped vertical bar per cluster: anchor (solid dark purple #7b2d8e) vs span (hatched light purple #c9a0dc). Dashed line at y=1.0 (background). Fold values annotated. |
| **D** | Loop size distribution | Compute `(y1 - x2) / 1000` from combined-clusters.txt | Boxplot (seaborn), `showfliers=False`. Median annotated. |
| **E** | Loop classification | `loop_classification.stats.txt` Percent stacked block | Stacked bar: structural (grey #7f7f7f) / CRE (green #2ca02c) / mixed (orange #ff7f0e) / unclassified (light grey #c7c7c7). |
| **F** | H2AK119ub signal | Fresh BigWig query (~24s) | Grouped bar: ctrl (grey #aaaaaa) vs mut (forestgreen #2ca02c). Bottom panel — has full x-axis labels. |

**Key implementation details:**
- Create bottom axis first (`ax_f`), then create all other axes with `sharex=ax_f`
- All upper panels: `ax.tick_params(labelbottom=False)`
- Bottom panel sets multi-line biological labels with `rotation=45, ha='right', fontsize=7`
- Use `ax.bar()` directly with explicit `x` positions (integers 0-5), NOT through Popay's plotting modules (which create their own figures)

## Figure 2: Anchor vs Span Mechanism

**Output:** `cluster/bap1_late/figures/summary_figures/mechanism/anchor_vs_span_mechanism.{png,pdf,svg,jpg}`

**Layout:** `figsize=(8, 5.5)`, `GridSpec(1, 2, wspace=0.4)`. Left panel = clust5 (gained), right = clust6 (lost).

Each panel shows all 12 ChromHMM states as **horizontal bars** (states on y-axis, fold enrichment on x-axis):
- Anchor enrichment: solid bar in state-specific color
- Span enrichment: hatched bar (`///`) in same color, alpha=0.6
- Vertical dashed line at x=1.0 (genome background)
- States ordered bottom-to-top by emission number (Quiescent at bottom, Bivalent_Enhancer at top — Polycomb states visually prominent at top)

**Annotations:**
- Left panel title: "Gained loops (clust5, n=667)\nPolycomb at anchor AND span"
- Right panel title: "Lost loops (clust6, n=2,359)\nPolycomb at anchor ONLY"
- Text box in each panel with top-3 gene examples from annotation files
- Mechanism label: "Domain compaction" (left) / "Anchor disruption" (right)

## Figure 3: Feature Summary Heatmap

**Output:** `cluster/bap1_late/figures/summary_figures/heatmap/feature_heatmap.{png,pdf,svg,jpg}`

**Layout:** `figsize=(8, 4)`, single `seaborn.heatmap`.

**Rows** (6): Clusters in biological gradient order (clust6 at top, clust5 at bottom). Row labels = biological names.

**Columns** (12):

| Column | Source | Raw value |
|--------|--------|-----------|
| Median logFC | metadata join | float |
| % Gained | diff status Row % | float |
| % Lost | diff status Row % | float |
| Median size (kb) | compute from coords | int |
| % Structural | classification Pct | float |
| % CRE | classification Pct | float |
| Polycomb anchor (fold) | anchor.txt 11_Polycomb | float |
| Polycomb span (fold) | span.txt 11_Polycomb | float |
| Bivalent anchor (fold) | anchor.txt 12_Bivalent_Enhancer | float |
| % Polycomb (proportion) | chromHMM_anchor.proportions.tsv | float |
| K119ub ctrl (median) | BigWig query | float |
| K119ub mut (median) | BigWig query | float |

**Normalization:** Column-wise z-score via `scipy.stats.zscore(axis=0)`. Colormap `RdBu_r`, `center=0`, `vmin=-2, vmax=2`. Each cell annotated with raw value (2 decimal places). Column labels rotated 35 degrees.

**Gene annotation panel:** Right margin text listing top-3 genes per cluster row.

## Script Structure

```
08_summary_figures.py
├── Imports + path constants (same pattern as 05_grouped_analyses.py)
├── BIO_ORDER = ['clust6','clust3','clust1','clust2','clust4','clust5']
├── BIOL_LABELS dict (filled with n= after data load)
├── Color constants (LOGFC_COLORS, DIRECTION_COLORS, etc.)
│
├── load_cluster_data()          # combined-clusters + metadata join
├── load_enrichment_tables()     # anchor.txt + span.txt, drop Base row
├── load_proportions()           # chromHMM_anchor.proportions.tsv
├── parse_stats_table()          # generic stats .txt block parser
├── compute_k119ub()             # pyBigWig queries on 2 BigWigs
├── extract_top_genes(n=5)       # annotation files -> top gene names
│
├── make_dashboard(...)          # Figure 1: 6-panel composite
├── make_mechanism_figure(...)   # Figure 2: clust5 vs clust6 anchor/span
├── make_feature_heatmap(...)    # Figure 3: z-scored summary heatmap
│
└── main()                       # load all -> compute K119ub -> render 3 figs
```

## Stats File Parsing

The two stats files have different formats:

**`cluster_differential_status.stats.txt`** — pandas `to_string()` output with named index:
```
Row %:
direction  down_in_mutant  unchanged  up_in_mutant
GROUP                                             
clust1               0.00     100.00          0.00
...
```
Parse: find "Row %:" line, skip 2 header lines (column names + index name), read data lines as whitespace-separated.

**`loop_classification.stats.txt`** — pandas `to_string()` with unnamed index:
```
Percent stacked:
              clust1  clust2  clust3  clust4  clust5  clust6
structural     32.20   40.39   31.90   52.27   54.72   29.33
...
```
Parse: find "Percent stacked:" line, line[0] = column headers, lines[1:] = data rows.

A single `parse_stats_block(path, block_label)` helper handles both by detecting whether the second line is an index-name-only line.

## K119ub Computation

Follow the exact `_bw_anchor_mean()` pattern from `05_grouped_analyses.py`:
1. Load combined-clusters.txt, split by GROUP
2. For each cluster, concatenate chr1 anchors + chr2 anchors, deduplicate
3. Open BigWig, query `bw.stats(chrom, start, end, type="mean")` per anchor
4. Wrap in try/except (RuntimeError, ValueError) -> np.nan
5. Average chr1+chr2 anchor signal per loop
6. Return per-cluster median as DataFrame

BigWig paths: `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw`

## Gene Extraction

Annotation files: `cluster/bap1_late/figures/annotation/clust{N}_annotation.txt`
- 20 columns, tab-separated, with header
- Gene name columns: `chr1_gene_name` (col 12) and `chr2_gene_name` (col 19)
- Extract unique non-empty gene names, count frequency (a gene at both anchors counts twice), take top-N per cluster
- Return `Dict[str, List[str]]`

## Key Implementation Notes

1. **No Popay plotting modules** — use matplotlib/seaborn directly for full layout control. `plotting.stacked()` creates its own figure; incompatible with GridSpec composite.
2. **`multi_format_savefig`** wraps each `fig.savefig()` call, NOT the entire function.
3. **`bbox_inches='tight'`** on every savefig — critical for multi-line x-labels.
4. **ChromHMM anchor.txt has a "Base" row** (line 14) that must be dropped before plotting.
5. **Column names in anchor.txt/span.txt** are `clust1.bed` through `clust6.bed` — strip `.bed` suffix.
6. **seaborn 0.13.2** is available in the cluster env. Use `sns.boxplot()` for panel D.
7. **Python 3.8** in cluster env — use `typing.Dict`, `typing.List`, `typing.Tuple` (not PEP 604).
8. **custom_params.json** — load from `cluster/custom_params.json` for consistent styling.
9. **Annotation file col 1** (clust5) may have rows with blank gene names — filter with `str != ''` and `notna()`.

## Verification

1. Run: `/opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/08_summary_figures.py`
2. Check 3 figure directories exist with 4 formats each:
   ```
   ls cluster/bap1_late/figures/summary_figures/{dashboard,mechanism,heatmap}/
   ```
3. Visual checks:
   - Dashboard: x-axis labels show biological names, gradient flows loss->gain left->right
   - Dashboard panel C: clust5 has BOTH anchor AND span bars above 1.0; clust6 has only anchor above 1.0
   - Mechanism: clust5 Polycomb bars clearly elevated in both solid+hatched; clust6 only solid
   - Heatmap: z-score colors show clear gradient from loss (left/top) to gain (right/bottom)
4. K119ub values: clust5 should show elevated mut vs ctrl signal (Polycomb domain = ub accumulation)
5. Gene names: non-empty for all 6 clusters, no NaN artifacts
