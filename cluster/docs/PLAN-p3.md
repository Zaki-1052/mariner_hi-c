---

## Phase 7: Summary Figures for Lab Meeting — DONE (2026-04-27)

**Status:** Three composite figures produced for lab-meeting / publication context, integrating Phase 4 outputs into integrated visualizations. Single Python orchestrator at `cluster/scripts/08_summary_figures.py` produces (1) Cluster Summary Dashboard (6-panel vertical), (2) Anchor vs Span Mechanism (clust5 vs clust6 focused), (3) Feature Summary Heatmap (z-scored, all clusters × all features). All inputs are pre-computed by Phases 1–4 except H2AK119ub anchor signal, which is queried fresh from BigWig files (~30s for ~63K unique anchors). Total runtime ~1 min on Mac. All 3 figures emit PNG + PDF + SVG + JPG via the Figure.savefig patch, plus a `feature_values.tsv` sidecar for the heatmap. Driver: `cluster/scripts/run_summary.sh`.

**Corrections applied during execution (Phase 8+ must follow these, not the original plan):**

- **Established biological cluster ordering:** `BIO_ORDER = ['clust6', 'clust3', 'clust1', 'clust2', 'clust4', 'clust5']` (loss → unchanged → gain). This ordering is reused by Phase 8 visualization scripts. The numeric clust1..clust6 ordering from Phase 3 v2 (descending mean signal) is preserved on disk; biological reorder happens only at the figure-rendering layer.
- **`BIO_NAMES` dict is project-canonical:** `clust6 = 'Strong loss\n(anchor-disrupted)'`, `clust3 = 'Moderate loss'`, `clust1 = 'Unchanged\n(high signal)'`, `clust2 = 'Unchanged'`, `clust4 = 'Moderate gain'`, `clust5 = 'Strong gain\n(Polycomb domain)'`. Used by `08_summary_figures.py` Panel-F labels and by `09_oriented_anchor_metagene.py` / `visualize_orientation_asymmetry.py`. Future scripts that label clusters in figures should reuse these strings.
- **Polycomb state row in `anchor.txt` / `span.txt` uses the OverlapEnrichment-prefix `'11_'`** (i.e. row label = `'11_Polycomb'` not just `'Polycomb'`). The state-id-prefix-stripping done in Phase 4.5 strips the leading `11_` for proportions display, but the OverlapEnrichment `.txt` files retain it — Panel C indexes via `anchor_enrich.loc['11_Polycomb', BIO_ORDER]`. Same for `'12_Bivalent_Enhancer'`. If the rename file changes, update the literal index strings here.
- **Phase 4 stats text files (`*.stats.txt`) parsed via fragile `StringIO + re` pattern.** `parse_direction_pct()` finds the `'Row %:\n'` marker and splits subsequent lines; `parse_classification_pct()` finds `'Percent stacked:\n'`. If Phase 4 changes the text emitted to those files (e.g. removes the marker, swaps column order), Phase 7 will silently produce wrong numbers. Defensive: call `df.index.name = None` on the parsed crosstab to avoid the `KeyError('xcol')` in `plotting.stacked` (same fix as Phase 4.8).
- **Annotation files use `chr1_gene_name` and `chr2_gene_name` columns** (not `gene_name` / `nearest_gene` etc). `extract_top_genes()` reads both and `.value_counts()` to get the most-frequent genes per cluster. Tolerates blanks and works through `on_bad_lines='warn'` for older pandas.
- **Heatmap uses pandas `apply(scipy.stats.zscore)` per column** then `.clip(-2.5, 2.5)`. This z-scores within each feature across the 6 clusters — NOT across all features. Raw values are still annotated in cells (with format depending on column type: `%`, `kb`, generic). RdBu_r colormap centered at 0.
- **`feature_values.tsv` sidecar** records the raw (non-z-scored) feature matrix for Fig 3, in case the heatmap needs to be regenerated externally or features need to be added to a downstream analysis.
- **K119ub query is the only "live" computation** — all other panels use precomputed Phase 4 artifacts. ~63K unique anchor regions × 2 BigWigs takes ~30s; the script falls back gracefully (Panel F shows "K119ub BigWigs not found" placeholder text) if `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` are missing.
- **`run_summary.sh` is intentionally minimal** — no `PYTHONUNBUFFERED=1` or `python3 -u` like Phase 5 needs, since this script is fast and stdout-only. Logs to `cluster/phase8_summary.txt` if `LOG` env var unset (note: even though the phase number is 7 in this plan, the driver's default log path predates the renumbering and uses `phase8_summary.txt`; not worth changing since it's just a path string).

### 7.1 — Loaders (`load_*` functions in `08_summary_figures.py`)

| Loader | Source | Output |
|--------|--------|--------|
| `load_cluster_data()` | `combined-clusters.txt` + `late_merged_loop_metadata.tsv` | Per-loop df with `GROUP`, `logFC`, `direction`, `loop_size_kb` |
| `load_enrichment_tables()` | `bap1_late/chromHMM/{anchor,span}.txt` | Two 12×6 fold-enrichment matrices (state × cluster) |
| `load_proportions()` | `chromHMM_anchor.proportions.tsv` | State proportion-per-cluster table |
| `parse_direction_pct()` | `cluster_differential_status.stats.txt` | 6×3 crosstab (cluster × direction) |
| `parse_classification_pct()` | `loop_classification.stats.txt` | 4×6 crosstab (class × cluster) |
| `compute_k119ub()` | BigWigs at `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` | 6×2 median signal matrix |
| `extract_top_genes(n=5)` | `bap1_late/figures/annotation/clust*_annotation.txt` | Dict[cluster → top-n gene-name list] |

### 7.2 — Figure 1: Cluster Summary Dashboard (6-panel vertical)

`make_dashboard()` produces 6 stacked panels with shared x-axis, biological cluster ordering on x-axis labels (Panel F only):

| Panel | What | Source |
|-------|------|--------|
| A: Median logFC | Bar chart, signed, centered at 0 | edgeR via metadata |
| B: Direction (%) | Stacked bar (lost / unchanged / gained) | Phase 4.8 stats |
| C: **Polycomb anchor vs span (KEY)** | Grouped bar, anchor solid / span hatched, with "KEY RESULT" annotation | Phase 4.4 enrichment |
| D: Loop size (kb) | Box plot, median annotated | Phase 1 loop_size_kb |
| E: Loop type (%) | Stacked bar (structural / CRE / mixed / unclassified) | Phase 4.2 stats |
| F: H2AK119ub anchor signal | Grouped bar, ctrl grey / mut green | Live BigWig query |

### 7.3 — Figure 2: Anchor vs Span Mechanism (clust5 vs clust6, 2-panel)

`make_mechanism_figure()` produces a side-by-side comparison of clust5 vs clust6 with:
- horizontal bars per ChromHMM state (anchor solid + span hatched)
- "Domain compaction" / "Anchor disruption" mechanism labels
- Top 3 genes per cluster annotated inline
- Color-coded by ChromHMM state (`STATE_COLORS` dict)

### 7.4 — Figure 3: Feature Summary Heatmap

`make_feature_heatmap()` z-scores 12 features per cluster, renders as RdBu_r heatmap with:
- raw values annotated in cells (formatted by feature type)
- top-3 genes per cluster as italic annotations on the right edge
- Sidecar `feature_values.tsv` with the un-zscored matrix

### Phase 7 Verification

**Verified 2026-04-27 (11:25 PDT, ~1 min runtime):**

- `cluster/bap1_late/figures/summary_figures/{dashboard,mechanism,heatmap}/` — 3 subfolders all present.
- `dashboard/cluster_dashboard.{png,pdf,svg,jpg}`: 316KB / 52KB / 137KB / 369KB respectively. 6-panel composite (A through F) renders cleanly with shared x-axis.
- `mechanism/anchor_vs_span_mechanism.{png,pdf,svg,jpg}`: 219KB / 68KB / 121KB / 278KB. 2-panel side-by-side with mechanism labels and top-genes annotations.
- `heatmap/feature_heatmap.{png,pdf,svg,jpg}`: 286KB / 62KB / 125KB / 226KB. 6×12 z-scored heatmap with raw-value annotations and per-cluster top-3 gene lists.
- `heatmap/feature_values.tsv`: 6 rows × 12 columns. Sample values for clust5 (Strong gain): logFC=+0.50, %Gained=97.3%, size=330kb, Structural=54.7%, CRE=2.85%, Polycomb-anchor=6.59×, Polycomb-span=3.03×, Bivalent-anchor=7.91×, %Polycomb-anchor=86.96%, K119ub-ctrl=1.37, K119ub-mut=1.25. clust6 (Strong loss): logFC=-0.31, %Lost=78.5%, size=550kb, Structural=29.3%, CRE=30.9%, Polycomb-anchor=2.09×, Polycomb-span=0.94×, %Polycomb-anchor=25.5%. Confirms Phase 4.4 KEY result and adds gene context.
- K119ub fresh-query worked: ~63K unique anchors × 2 BigWigs in ~30s. Per-cluster medians: clust6 ctrl=1.16/mut=1.58 (mut up at lost-loop anchors), clust5 ctrl=1.37/mut=1.25 (both high, slight ctrl > mut), clust1 ctrl=1.06/mut=0.96 (both modest). Pattern is mixed — H2AK119ub_mut > ctrl at clust3 and clust6 (loss-biased), but not at clust4/clust5 (gain-biased) — consistent with the Phase 4.3 Kruskal-Wallis omnibus showing increased mut-side variance, not a uniform mut-up shift.
- All figures emit 4 formats per the Figure.savefig patch. Phase 4's correction (`Figure.savefig` not `plt.savefig`) carries through here — confirmed by file listings.

**Files created in Phase 7:**

| File | Purpose |
|------|---------|
| `cluster/scripts/08_summary_figures.py` | Phase 7 orchestrator: 3 composite figures (dashboard, mechanism, heatmap); reuses Phase 4 artifacts; live K119ub query |
| `cluster/scripts/run_summary.sh` | Driver: invokes `08_summary_figures.py` via cluster-env python; tees stdout; lists outputs at end |

**Outputs created in Phase 7 (under `cluster/bap1_late/figures/summary_figures/`):**

| Path | What |
|------|------|
| `dashboard/cluster_dashboard.{png,pdf,svg,jpg}` | 6-panel cluster summary (logFC / direction / Polycomb anchor-vs-span / size / classification / K119ub) |
| `mechanism/anchor_vs_span_mechanism.{png,pdf,svg,jpg}` | clust5 vs clust6 horizontal-bar comparison with mechanism labels |
| `heatmap/feature_heatmap.{png,pdf,svg,jpg}` | z-scored 6×12 heatmap with raw-value annotations + top-3 genes |
| `heatmap/feature_values.tsv` | Raw (non-z-scored) feature matrix for Fig 3 — 6 rows × 12 cols |

**Re-running:**
```bash
LOG=cluster/phase8_summary.txt bash cluster/scripts/run_summary.sh
# Note: log file path is `phase8_summary.txt` even though the phase is 7 in this plan;
# the path was set before phase renumbering and isn't worth changing.
```

---

## Phase 8: Oriented Anchor Metagene + Asymmetry Quantification — DONE (2026-04-27)

**Status:** Three-step pipeline tests whether K27me3 (and the other 7 BigWig marks) is asymmetrically enriched on the EXTERIOR (away from loop body) vs INTERIOR (toward loop body) side of loop anchors. (1) `09_oriented_anchor_metagene.py` builds per-cluster oriented BED6 files (anchor1=strand+ / anchor2=strand-) and runs deepTools `bed_pileup` with the same 8 BigWigs and ±5kb window as Phase 5. deepTools natively flips signal arrays for strand- regions, so plot left = -5kb (exterior) / right = +5kb (interior) for both anchor types. (2) `quantify_orientation_asymmetry.py` reads the resulting `oriented_anchors_values` matrix, computes per-cluster per-mark exterior vs interior means, and runs Wilcoxon signed-rank tests. Outputs `asymmetry_quantification.tsv` (48 rows = 6 clusters × 8 marks). (3) `visualize_orientation_asymmetry.py` produces a focused dual-panel figure (per-cluster K27me3 metagene profiles at top; ext/int ratio bar chart with significance asterisks at bottom). Total runtime: ~1.7h for the metagene step (matches Phase 5), <1 min for quantify, <1 min for visualize.

**KEY biological finding:** K27me3 is **INTERIOR-enriched** at clust4 (mod gain) and clust5 (strong gain) anchors — Ext/Int ratio 0.93 / 0.91 with p<0.005 in both ctrl and mut. Polycomb spreads INTO the loop body, supporting the Polycomb-domain compaction model: clust5 anchors sit AT THE EDGE of expanding Polycomb domains. Conversely, **clust6 (strong loss) shows NO K27me3 asymmetry** (Ext/Int = 1.02 ctrl, 1.02 mut, p ≈ 0.52 both) — Polycomb gain at lost-loop anchors is symmetric, NOT a side-specific extrusion barrier. This refines the Phase 4.4 mechanism: clust6's anchor-Polycomb enrichment (2.09×) is from anchors flanked by Polycomb on BOTH sides (consistent with anchor-disruption / "sensitivity model"), while clust5's (6.59×) reflects the anchor sitting at the boundary of a one-sided Polycomb domain (consistent with extrusion impediment / domain-compaction).

**Corrections applied during execution:**

- **Strand encoding convention:** anchor1 (left, x1<y1) gets strand `+`; anchor2 (right) gets strand `-`. After deepTools native strand-flip, plot left ALWAYS = exterior (away from loop body), plot right ALWAYS = interior (toward loop body). This is the simplest correct encoding — alternatives (e.g. centering on anchor midpoint with custom flipping) are more error-prone.
- **Dedup on `(chrom, start, end, strand)`** — different from Phase 5 which dedups on `(chrom, start, end)`. Hub anchors that participate in different loops as both anchor1 AND anchor2 keep BOTH orientation entries (they appear once with `+` and once with `-`). For clust1: 12,298 loops × 2 = 24,596 raw anchors → 20,739 oriented BED rows (930 hub anchors retain dual orientation; 80.5% of raw anchors deduplicated to single-orientation entries). Per-cluster oriented anchor counts: clust1=20,739 (930 hub) / clust2=18,752 (589) / clust3=14,424 (665) / clust4=6,887 (104) / clust5=1,198 (11) / clust6=3,964 (147).
- **Post-computeMatrix retention rate matches Phase 5 (~98%).** computeMatrix dropped some rows due to blacklist overlap or BigWig signal gaps. Rows surviving for asymmetry quantification: clust1=20,373 / clust2=18,423 / clust3=14,137 / clust4=6,749 / clust5=1,171 / clust6=3,860 (1.5–2.4% dropped per cluster, similar to Phase 5).
- **8 BigWigs identical to Phase 5** (4 marks × ctrl/mut; H3K27ac, H3K27me3, H2AK119ub, H3K27me1). Same `vmax_groups` pairing for ctrl/mut visual consistency. Same color_dict (Blues/Reds/Greens/Purples). Same blacklist (mm10 v2). Same ±5kb window with default `--binSize 10`. Same `pileup_type='referencePoint'`. Net effect: oriented metagene heatmap structure mirrors Phase 5's, but rows are flipped for strand-`-` anchors so spatial asymmetry becomes interpretable.
- **Same intermediate-file caveat as Phase 5.** `oriented_anchors_values` (1.6GB tab text) and `oriented_anchors` (100MB gz binary) are intermediate `computeMatrix` outputs needed only by `heatmap_plot` and the quantify script. Add to `.gitignore` along with the Phase 5 intermediates:
  ```
  cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors
  cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values
  ```
- **`xticklabels=['-5kb (exterior)', 'anchor', '+5kb (interior)']`** passed to `bed_pileup()` — deepTools 3.5.5 honors this and labels the heatmap x-axis correctly. Without this override, the default labels are `['-5.0Kb', '0', '+5.0Kb']` which is correct numerically but loses the biological interpretation.
- **Asymmetry index sign convention:** `asymmetry_index = (exterior_mean − interior_mean) / (exterior_mean + interior_mean)`. POSITIVE = exterior-enriched, NEGATIVE = interior-enriched. K27me3 interior enrichment for clust5 → asymmetry_index = -0.045 (ctrl) / -0.048 (mut). Reported in `asymmetry_quantification.tsv`.
- **Wilcoxon signed-rank** uses paired `(exterior_per_anchor, interior_per_anchor)` differences with `alternative='two-sided'` and `nonzero` filter (drops anchors where ext = int exactly). Bonferroni not applied per-mark because each test answers a distinct biological question; per-mark FDR or family-wise correction can be added if needed for publication.
- **Plot styling:** `visualize_orientation_asymmetry.py` uses BIO_ORDER from Phase 7 (loss → unchanged → gain). Shaded background (blue=exterior, peach=interior) makes the asymmetry direction visually obvious. Significance asterisks: `*` p<0.05, `**` p<0.01, `***` p<0.001. Y-axis range `[0.885, 1.055]` is hand-tuned for our K27me3 ratio range; a more general script would auto-fit.
- **Driver `run_oriented_metagene.sh` mirrors `run_phase5.sh`:** prepends cluster env bin to PATH (for `computeMatrix` subprocess), sets `PYTHONUNBUFFERED=1`, invokes `python3 -u`, tees to `cluster/oriented_metagene.txt`. The quantify and visualize scripts are run separately by the user (no driver) — they're fast (<1 min each) and produce stdout-only diagnostic tables that the user can capture as needed.

### 8.1 — `09_oriented_anchor_metagene.py` (the metagene step, ~1.7h on Mac)

`build_oriented_anchor_beds()` (line 80):
1. Read `combined-clusters.txt`
2. For each cluster, build anchor1 BED6 with strand=`+`, anchor2 BED6 with strand=`-`
3. Concat → dedup on `(chrom, start, end, strand)` → sort by `(chrom, start, end, strand)`
4. Write `cluster{1..6}_oriented_anchors.bed`

`bed_pileup()` (deepTools_pipeline.py:25) called with all 8 BigWigs, all 6 oriented BEDs, `up_down=5000`, `pileup_type='referencePoint'`, blacklist, paired ctrl/mut vmax groups, custom xticklabels.

Outputs: `oriented_anchors_values.{png,pdf,svg,jpg}` (the heatmap), `oriented_anchors_values` (1.6GB text matrix), `oriented_anchors` (100MB gz binary).

### 8.2 — `quantify_orientation_asymmetry.py` (the analysis step, <1 min)

`parse_header()` reads the first 3 lines of `oriented_anchors_values` to extract cluster names, sizes, BigWig labels, bin size, and bin count. Then loads the data matrix, splits each BigWig's bins into exterior half (-5kb to 0) and interior half (0 to +5kb), computes per-anchor mean signal in each half, runs Wilcoxon, produces `asymmetry_quantification.tsv`.

Output table columns: `mark, cluster, direction, n_anchors, exterior_mean, interior_mean, ext_int_ratio, asymmetry_index, wilcoxon_pval`.

Stdout: human-readable per-mark significance tables + a `CLUST5 vs CLUST6 COMPARISON` block summarizing the key biological contrast.

### 8.3 — `visualize_orientation_asymmetry.py` (the figure step, <1 min)

Produces `orientation_asymmetry_figure.{png,pdf,svg,jpg}` — a dual-panel composite:
- **Top:** 1×6 grid of K27me3 metagene profiles, one panel per cluster, in BIO_ORDER. Each panel shows ctrl (blue) and mut (red) mean ± SEM lines with shaded backgrounds (blue exterior, peach interior). Per-panel inset: n, ctrl Ext/Int ratio + significance, mut Ext/Int ratio + significance.
- **Bottom:** Bar chart of K27me3 ctrl + mut Ext/Int ratios across all 6 clusters with significance asterisks. Reference line at 1.0 (symmetric); above = exterior-enriched, below = interior-enriched.

### Phase 8 Verification

**Verified 2026-04-27:**

- `cluster/bap1_late/figures/deeptools_input/clust{1..6}_oriented_anchors.bed`: 6 BED6 files, 35KB–606KB. All sorted, all deduplicated on `(chrom, start, end, strand)`. Strand counts roughly balanced (anchor1+ ≈ anchor2-; minor skew from hub-anchor dedup).
- `cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values.{png,pdf,svg,jpg}`: PNG=2.3MB, PDF=1.1MB, SVG=665KB, JPG=399KB. Heatmap structure matches Phase 5's (6 cluster rows × 8 BigWig columns) with reference point at center; clear exterior-vs-interior asymmetry visible by eye for clust4/clust5 K27me3 panels.
- `oriented_anchors_values` 1.6GB tab-text matrix and `oriented_anchors` 100MB gz binary present (intermediate; should be gitignored).
- `cluster/bap1_late/figures/deeptools/oriented_anchors/asymmetry_quantification.tsv`: 48 rows (8 marks × 6 clusters) + header. All 8 marks × 6 clusters covered. Wilcoxon p-values populated for all rows.
- **K27me3 quantitative results (matches the orientation_asymmetry_figure):**
  | Cluster | Direction | n | ctrl Ext/Int | ctrl p | mut Ext/Int | mut p |
  |---------|-----------|---|--------------|--------|-------------|-------|
  | clust1 | unchanged | 20,373 | 1.015 | 0.541 | 1.002 | 0.575 |
  | clust2 | ~unchanged | 18,423 | 0.952 | 5.9e-6 *** | 0.956 | 0.016 * |
  | clust3 | mod loss | 14,137 | 0.968 | 0.433 | 0.983 | 0.554 |
  | clust4 | mod gain | 6,749 | 0.928 | 5.1e-5 *** | 0.922 | 8.6e-6 *** |
  | **clust5** | **strong gain** | **1,171** | **0.914** | **0.004 \*\*** | **0.909** | **0.004 \*\*** |
  | **clust6** | **strong loss** | **3,860** | **1.020** | **0.522 ns** | **1.023** | **0.571 ns** |

  Confirms KEY finding: K27me3 interior enrichment is the gain signature (clust5, clust4, weakly clust2); clust6 has NO asymmetry (symmetric Polycomb at lost-loop anchors).

- **`orientation_asymmetry_figure.{png,pdf,svg,jpg}`:** PNG=269KB, PDF=192KB, SVG=755KB, JPG=200KB. Top row of 6 panels shows per-cluster K27me3 profiles with ctrl+mut overlay. Bottom row bar chart shows the Ext/Int ratios with asterisks marking p<0.05.

**Files created in Phase 8:**

| File | Purpose |
|------|---------|
| `cluster/scripts/09_oriented_anchor_metagene.py` | Phase 8 orchestrator step 1: per-cluster oriented BED6 prep + bed_pileup with 8 BigWigs at ±5kb / strand-aware orientation |
| `cluster/scripts/quantify_orientation_asymmetry.py` | Phase 8 step 2: parse `oriented_anchors_values`, split bins into exterior/interior halves, Wilcoxon test → `asymmetry_quantification.tsv` |
| `cluster/scripts/visualize_orientation_asymmetry.py` | Phase 8 step 3: dual-panel K27me3 figure (per-cluster profiles + ext/int bar chart) |
| `cluster/scripts/run_oriented_metagene.sh` | Driver for step 1 only: PATH + PYTHONUNBUFFERED setup; tees to `cluster/oriented_metagene.txt`. Steps 2 and 3 are run manually by user. |

**Outputs created in Phase 8 (under `cluster/bap1_late/`):**

| Path | What |
|------|------|
| `figures/deeptools_input/clust{1..6}_oriented_anchors.bed` | 6 deduplicated BED6 files (1,198–20,739 rows; strand-aware) |
| `figures/deeptools/oriented_anchors/oriented_anchors` | computeMatrix gz binary matrix (100MB; intermediate, gitignore) |
| `figures/deeptools/oriented_anchors/oriented_anchors_values` | computeMatrix tab-text matrix (1.6GB; intermediate, gitignore) |
| `figures/deeptools/oriented_anchors/oriented_anchors_values.{png,pdf,svg,jpg}` | Combined heatmap (8 BigWigs × 6 clusters, strand-aware) |
| `figures/deeptools/oriented_anchors/asymmetry_quantification.tsv` | 48 rows × 9 cols: per-(mark,cluster) Wilcoxon results |
| `figures/deeptools/oriented_anchors/orientation_asymmetry_figure.{png,pdf,svg,jpg}` | K27me3 dual-panel figure (per-cluster profiles + Ext/Int bars) |

**Re-running:**
```bash
# Step 1 (~1.7h): metagene
LOG=cluster/oriented_metagene.txt bash cluster/scripts/run_oriented_metagene.sh

# Step 2 (<1 min): quantify
/opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/quantify_orientation_asymmetry.py | tee cluster/orientation_asymmetry_stats.txt

# Step 3 (<1 min): visualize (after Step 2)
/opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/visualize_orientation_asymmetry.py
```

---

## Directory Structure

```
cluster/
├── ChromHMM/                          # Existing — ChromHMM.jar + CHROMSIZES/
├── clustering_example_data/           # Existing — Popay example files
├── custom_params.json                 # Existing — matplotlib style params
├── *.py                               # Existing — Popay modules (modified in Phase 0)
│
├── scripts/                           # all runnable scripts
│   ├── 01_build_loop_count_file.py    #   Phase 1: merged loops → Popay format
│   ├── 02_build_mm10_gene_annotation.R#   Phase 1: mm10 gene BED
│   ├── 03_chromhmm_segmentation.sh    #   Phase 2: BinarizeBed + LearnModel
│   ├── 04_clustering.py               #   Phase 3: elbow + Cluster 3.0 + viz
│   ├── 05_grouped_analyses.py         #   Phase 4: 8 downstream sub-analyses (DONE)
│   ├── 06_deeptools_metagene.py       #   Phase 5: anchor metagene (DONE)
│   ├── 07_cooltools_pileup.py         #   Phase 6: cooltools pileup orchestrator (DONE on HPC)
│   ├── 07_cooltools_pileup.sb         #   Phase 6: SLURM wrapper for Expanse
│   ├── cooler_balance.sb              #   Phase 6 prereq: ICE balance mcools (NEW)
│   ├── 08_summary_figures.py          #   Phase 7: lab-meeting composite figures (DONE)
│   ├── 09_oriented_anchor_metagene.py #   Phase 8 step 1: strand-aware anchor metagene (DONE)
│   ├── quantify_orientation_asymmetry.py  # Phase 8 step 2: ext/int Wilcoxon → TSV (DONE)
│   ├── visualize_orientation_asymmetry.py # Phase 8 step 3: K27me3 dual-panel figure (DONE)
│   ├── run_phase1.sh                  #   Runner: data prep
│   ├── run_phase2.sh                  #   Runner: ChromHMM segmentation
│   ├── run_phase3.sh                  #   Runner: clustering
│   ├── run_phase4.sh                  #   Runner: downstream analyses
│   ├── run_phase5.sh                  #   Runner: deepTools metagene
│   ├── run_summary.sh                 #   Runner: Phase 7 summary figures
│   ├── run_oriented_metagene.sh       #   Runner: Phase 8 step 1 only
│   └── utils/
│       └── multi_format_output.py     #   Figure.savefig patch for 4-format output
│
├── logs/                              # NEW — HPC SLURM stdout (synced from Expanse)
│   ├── cooler_balance_<jobid>.out     #   Phase 6 prereq logs
│   └── phase6_pileup_<jobid>.out      #   Phase 6 pileup attempts (one per submission)
│
├── data/                              # input data files
│   ├── late_merged_loop_counts.txt    #   Phase 1 output (Popay format)
│   ├── late_merged_loop_metadata.tsv  #   Phase 1 output (edgeR stats sidecar)
│   └── mm10_knownGene_pp.bed          #   Phase 1 output (gene annotation)
│   # Phase 5/7/8 BigWigs source: /Users/zakiralibhai/sdsc/bigwigs/ (not in repo)
│   # Phase 6 mcools source:    /expanse/.../stripes/stripenn/data/cool/250402/ (HPC)
│
└── bap1_late/                         # all analysis outputs
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
    │   ├── annotation/                # Phase 4.6: per-cluster gene lists
    │   ├── ChIP_intersect/            # Phase 4.3, 4.7: DiffBind + ChIP signal
    │   ├── loop_size/                 # Phase 4.1: per-cluster loop-size dist.
    │   ├── loop_classification/       # Phase 4.2: structural/CRE/mixed/uncl.
    │   ├── chromHMM_anchor/           # Phase 4.5: ChromHMM proportions
    │   ├── cluster_differential_status/  # Phase 4.8: cluster x edgeR direction
    │   ├── deeptools_input/           # Phase 5: anchor BEDs + Phase 8: oriented anchor BEDs
    │   ├── deeptools/
    │   │   ├── histone_anchors/       # Phase 5: combined heatmap (DONE)
    │   │   │   └── histone_anchors_values.{png,pdf,svg,jpg}
    │   │   └── oriented_anchors/      # Phase 8: strand-aware metagene (NEW, DONE)
    │   │       ├── oriented_anchors_values.{png,pdf,svg,jpg}
    │   │       ├── asymmetry_quantification.tsv  # 48 rows = 8 marks × 6 clusters
    │   │       └── orientation_asymmetry_figure.{png,pdf,svg,jpg}  # K27me3 dual-panel
    │   └── summary_figures/           # Phase 7: lab-meeting figures (NEW, DONE)
    │       ├── dashboard/             #   cluster_dashboard.{png,pdf,svg,jpg}
    │       ├── mechanism/             #   anchor_vs_span_mechanism.{png,pdf,svg,jpg}
    │       └── heatmap/               #   feature_heatmap.{png,pdf,svg,jpg} + feature_values.tsv
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
    └── cooltools/                     # Phase 6 outputs (DONE)
        └── obs_exp_contacts/          #   obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}
```

---

## Runner Scripts

Each phase gets a `run_phase{N}.sh` bash script in `cluster/scripts/` that runs the steps for that phase sequentially with error checking. These are meant to be run manually one phase at a time (not chained).

### `run_phase0.sh` — Apply bug fixes + verify

```bash
#!/bin/bash
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
cd "$(dirname "$0")/.."
bash scripts/03_chromhmm_segmentation.sh
echo "Phase 2 complete. Verify:"
wc -l bap1_late/chromHMM/learned_model/*_12_segments.bed
echo "Now inspect emissions_12.txt and create 12state_rename_cerebellum.txt"
```

### `run_phase3.sh` — Clustering

```bash
#!/bin/bash
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
set -e
cd "$(dirname "$0")/../.."
CLUSTER_ENV_BIN=/opt/homebrew/anaconda3/envs/cluster/bin
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=cluster/scripts/06_deeptools_metagene.py
LOG=${LOG:-cluster/phase5.txt}
{
  "$PYTHON" -u "$SCRIPT"
} 2>&1 | tee "$LOG"
```

### `run_summary.sh` — Phase 7 summary figures

```bash
#!/bin/bash
set -e
cd "$(dirname "$0")/../.."
PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=cluster/scripts/08_summary_figures.py
LOG=${LOG:-cluster/phase8_summary.txt}
{
  "$PYTHON" "$SCRIPT"
  echo "--- Output inventory ---"
  for dir in dashboard mechanism heatmap; do
    ls -lh "cluster/bap1_late/figures/summary_figures/${dir}/" 2>/dev/null
  done
} 2>&1 | tee "$LOG"
```

### `run_oriented_metagene.sh` — Phase 8 step 1 only (metagene)

```bash
#!/usr/bin/env bash
# Driver for 09_oriented_anchor_metagene.py only.
# Run quantify + visualize separately after this finishes (~1.7h).
set -e
cd "$(dirname "$0")/../.."
CLUSTER_ENV_BIN=/opt/homebrew/anaconda3/envs/cluster/bin
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=cluster/scripts/09_oriented_anchor_metagene.py
LOG=${LOG:-cluster/oriented_metagene.txt}
{
  "$PYTHON" -u "$SCRIPT"
} 2>&1 | tee "$LOG"
```

### `cooler_balance.sb` — Phase 6 prerequisite (HPC SLURM)

ICE-balance mcools at a given resolution (writes the `weight` column required by `cooltools.expected_cis`). Runs ctrl + mut concurrently with 2 nproc each.

```bash
#!/bin/bash
#SBATCH --account=csd940 --partition=shared --cpus-per-task=4 --mem=64G --time=02:00:00
#SBATCH --output=cluster/logs/cooler_balance_%j.out
MCOOL_CTRL="$1"; MCOOL_MUT="$2"; RES="${3:-10000}"
source ~/.bashrc; conda activate cluster
cooler balance --nproc 2 "${MCOOL_CTRL}::/resolutions/${RES}" &
cooler balance --nproc 2 "${MCOOL_MUT}::/resolutions/${RES}"  &
wait
```

Usage:
```bash
sbatch cluster/scripts/cooler_balance.sb \
  /expanse/.../ctrl_merged.mcool \
  /expanse/.../mut_merged.mcool \
  10000
```

### `sync_from_hpc.sh` — rsync HPC data

```bash
#!/bin/bash
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
| `cluster/scripts/04_clustering.py` | 3 | Elbow plot + Cluster 3.0 k-means + visualization (incl. `--min-ratio` / `--max-ratio` for outlier removal) |
| `cluster/scripts/utils/multi_format_output.py` | 3+ | Context-manager savefig wrapper + figure_subfolder helper (PNG + PDF + SVG + JPG, subfolder per figure); shared by Phase 4-5 |
| `cluster/scripts/05_grouped_analyses.py` | 4 | All Phase 4 downstream analyses (ChromHMM, ChIP, DiffBind, gene annotation) |
| `cluster/scripts/06_deeptools_metagene.py` | 5 | Per-cluster anchor BED prep + bed_pileup (8 BigWigs x 6 clusters, single combined heatmap) |
| `cluster/scripts/07_cooltools_pileup.py` | 6 | Loads combined-clusters.txt → cooler-convention bedpe_dict; calls mcool_pileup with viewframe built inline from cooler chromsizes (no UCSC fetch) |
| `cluster/scripts/07_cooltools_pileup.sb` | 6 | SLURM wrapper for Expanse: csd940 / shared / 16 cpus / 64G / 4h, two-root paths |
| `cluster/scripts/cooler_balance.sb` | 6 (prereq) | NEW: SLURM wrapper to ICE-balance mcools at a given resolution (concurrent ctrl+mut). Required when source mcool only has KR weights (no `weight` column from cooler balance). |
| `cluster/scripts/08_summary_figures.py` | 7 | Lab-meeting summary figures: cluster dashboard (6-panel), anchor-vs-span mechanism (clust5/clust6), feature summary heatmap. Reuses Phase 4 artifacts; live K119ub query. |
| `cluster/scripts/09_oriented_anchor_metagene.py` | 8 | Strand-aware anchor metagene: anchor1=+/anchor2=- → deepTools native flip → exterior/interior asymmetry-ready matrix (8 BigWigs × 6 clusters) |
| `cluster/scripts/quantify_orientation_asymmetry.py` | 8 | Reads oriented_anchors_values, splits ±5kb window into ext/int halves, runs Wilcoxon signed-rank → asymmetry_quantification.tsv (48 rows) |
| `cluster/scripts/visualize_orientation_asymmetry.py` | 8 | K27me3-focused dual-panel figure (per-cluster metagene profiles + ext/int bar chart with significance) |
| **Runner scripts** | | |
| `cluster/scripts/run_phase1.sh` | 1 | Data prep + directory setup |
| `cluster/scripts/run_phase2.sh` | 2 | ChromHMM segmentation |
| `cluster/scripts/run_phase3.sh` | 3 | Clustering — accepts `K FILTER_PCT MIN_RATIO MAX_RATIO` positional + `LOG` env var |
| `cluster/scripts/run_phase4.sh` | 4 | Downstream analyses — accepts `ANALYSES` positional (e.g. `4.4` or `all`) + `LOG` env var |
| `cluster/scripts/run_phase5.sh` | 5 | deepTools metagene — prepends cluster env bin to PATH, sets PYTHONUNBUFFERED=1, invokes `python3 -u` |
| `cluster/scripts/run_summary.sh` | 7 | Summary figures — minimal driver, tees stdout to `cluster/phase8_summary.txt` (legacy filename, predates phase renumbering) |
| `cluster/scripts/run_oriented_metagene.sh` | 8 | Oriented metagene step 1 only — same PATH/PYTHONUNBUFFERED setup as Phase 5; tees to `cluster/oriented_metagene.txt`. Steps 2–3 run manually. |

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

1. **Phase 0:** ✅ All Popay modules import without error (`plotting`, `cluster_tools`, `chromHMM_heatmap`, `deeptools_plotting`, `bedpe_analysis`, `cooltools_called`, `deepTools_pipeline`)
2. **Phase 1:** ✅ `late_merged_loop_counts.txt` has 39,344 rows, 8 columns; `mm10_knownGene_pp.bed` has 24,515 entries
3. **Phase 2:** ✅ Segmentation BED covers all 21 standard chroms (chr1-19 + chrX + chrY), 335,569 segments; 12-state emissions matrix biologically interpretable (Active_Promoter through Quiescent)
4. **Phase 3:** ✅ Elbow plot shows bend at k=4-5; v2 combined-clusters.txt has 38,948 rows, all 6 clusters > 600 loops, captures 99%+ of edgeR differential calls; multi-format outputs (PNG + PDF + SVG + JPG) per figure subfolder
5. **Phase 4.4 (KEY result):** ✅ ChromHMM anchor vs span heatmaps show differential Polycomb enrichment across clusters — clust5 anchor 6.59x + span 3.03x (extrusion impediment); clust6 anchor 2.09x + span 0.94x (sensitivity model)
6. **Phase 4.8:** ✅ Chi-squared p ≈ 0 for cluster x differential status (clust5=97% up, clust6=78% down, clust1=100% unchanged)
7. **Phase 5:** ✅ deepTools metagene shows clust5 anchors enriched for H3K27me3 + H2AK119ub + H3K27me1 with low H3K27ac (Polycomb signature); clust1-3 K27ac-dominant; clust4-5 paired vmax shows visible mut increases in repressive marks. Heatmap PDF/PNG/SVG/JPG all present at `cluster/bap1_late/figures/deeptools/histone_anchors/`
8. **Phase 6:** ✅ Cooltools pileup completed on Expanse (job `48577675`, 18:43 PDT 2026-04-27). Required ICE-balancing the mcools first (new `cooler_balance.sb` script, job `48574993`, ~52 min) and replacing the §6.3 make_arms_robust patch with an inline cooler-chromsizes viewframe. Output: `cluster/bap1_late/cooltools/obs_exp_contacts/obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}`. Visual confirmation: clust5 mut > ctrl central pixel (gain), clust6 mut < ctrl (loss), clust1 ctrl ≈ mut (unchanged).
9. **Phase 7:** ✅ Three composite summary figures produced — cluster dashboard (6-panel), anchor-vs-span mechanism (clust5/clust6), feature summary heatmap. Output: `cluster/bap1_late/figures/summary_figures/{dashboard,mechanism,heatmap}/`. `feature_values.tsv` records the underlying 6×12 feature matrix for reuse.
10. **Phase 8:** ✅ Oriented anchor metagene + asymmetry quantification completed. KEY result: K27me3 INTERIOR-enriched at clust5 (Ext/Int=0.91, p=0.004 ctrl/mut) and clust4 (0.92–0.93, p<5e-5) — Polycomb-domain compaction at gained loops. clust6 (lost loops): NO K27me3 asymmetry (Ext/Int=1.02, p=0.52) — symmetric Polycomb at anchor-disrupted loops, NOT a side-specific extrusion barrier. Output: `cluster/bap1_late/figures/deeptools/oriented_anchors/{oriented_anchors_values,asymmetry_quantification.tsv,orientation_asymmetry_figure}.*`.