# Phase 3 Plan: K-means Clustering of BAP1-KO Hi-C Loops via Cluster 3.0

**Companion docs (already produced & verified):**
- `cluster/CONTEXT-CLUSTER.md` — full biological context, Dixon meeting, Popay collaboration
- `cluster/PLAN-CLUSTER.md` — master plan (Phases 0-6); Phases 0-2 marked DONE
- `cluster/phase1.txt`, `cluster/phase2.txt` — execution logs for completed phases
- `cluster/HiC_cluster3.ipynb` — Popay's reference notebook (cells 1-5 cover Phase 3 logic)

---

## 1. Context

**Why we are doing this.** The Dixon meeting (2026-04-10) and CTEA April directives both flagged the central unresolved question for this paper: is the H3K27me3 / H2AK119ub enrichment that drives BAP1-KO loop dysregulation localized at loop **anchors** (sensitivity model — repressive chromatin disrupts CTCF binding directly) or distributed across the loop **body / span** (extrusion-impediment model — Polycomb spreading blocks cohesin extrusion before the long loop completes)? Tessa Popay shared the exact pipeline she used to answer this for NIPBL depletion (Popay et al., Nat Genet 2026, Fig 2f) and we are adapting it to our system.

The downstream Fig 2f equivalent (Phase 4.4) requires loops grouped into k-means clusters by their normalized contact-frequency profile across conditions — which is what this phase produces. **Phase 3 is the bridge between data prep (Phases 1-2) and the mechanistic ChromHMM anchor-vs-span analysis (Phase 4.4).**

**What Phases 1-2 produced (verified, reused as-is):**
- `cluster/data/late_merged_loop_counts.txt` — 39,344 loops × 8 cols (`chr1 x1 x2 chr2 y1 y2 ctrl_merge mut_merge`), Popay format
- `cluster/data/late_merged_loop_metadata.tsv` — 39,344 × 16 cols, includes `resolution_kb` for per-resolution operations
- `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` — 335,569 segments, 12 biologically-named states (Active_Promoter, Polycomb, etc.)

**What Phase 3 produces:**
- Elbow plot (k=1..14) for k selection
- Per-loop k-means cluster assignments via Cluster 3.0 binary (`-r 100 -g 7 -k {k}`)
- `combined-clusters.txt` = `(GROUP, chr1, x1, x2, chr2, y1, y2, ctrl_merge, mut_merge)` with `GROUP ∈ {clust1..clustK}`, sorted descending by mean signal
- 4 visualizations per k: heatmap (subplot per cluster), lineplot, stripplot, boxplot

---

## 2. Key Design Decisions (with rationale)

### 2.1 Filter threshold — **per-resolution 1%-tile of `ctrl_merge`**

The parent plan's `filter_threshold = 0.008` is mis-calibrated (it was Popay's filter on KR-balanced counts in 0.05–0.30 range). Our mariner-aggregated averaged counts span ~5.5 to 10,422, with **per-resolution medians differing ~6×** (5kb=171, 10kb=403, 25kb=984). Phase 1 already flagged this and required Phase 3 to choose empirically.

**Empirical thresholds (computed from actual data):**

| Resolution | 1%-tile | 2.5%-tile | 5%-tile |
|------------|---------|-----------|---------|
| 5 kb       | 38.0    | 47.8      | 59.1    |
| 10 kb      | 61.5    | 84.4      | 107.2   |
| 25 kb      | 56.5    | 98.7      | 159.0   |
| **Total dropped** | **394 (1.0%)** | **985 (2.5%)** | **1,968 (5.0%)** |

**Recommendation: 1%-tile per resolution.** Rationale:
- mariner's `filterByExpr(min.count=5)` already filtered weak loops upstream; we are applying a second-line cut only against the noisiest tail
- 1%-tile drops just 394 loops (essentially the noise-floor) yet **directly addresses the mut/ctrl ratio outlier problem** — the max ratio of 10.08× and 99.9%-tile of 1.68× are driven by loops with very low `ctrl_merge` denominators (the same loops the 1%-tile cut removes)
- Per-resolution preserves the 5kb/10kb/25kb balance — a global threshold would disproportionately filter 5kb (whose absolute counts are naturally smaller)
- Conservative: keeps 99% of the dataset for clustering

Expose `--filter-pct` CLI flag so this can be re-run at 5%-tile or any other percentile if the elbow plot shows pathology (e.g., one cluster absorbing extreme outliers).

### 2.2 k selection — **start with k=6, generate elbow first**

Popay used k=6 for her RAD21-dependence analysis (matching the published Fig 1c-d clusters: min-dep / low-dep / med-dep / mixed-dep / high-dep / max-dep). Our biology is different (developmental BAP1 loss vs. acute NIPBL degron) and our data is effectively 1D in normalized space (since `ctrl_norm = ctrl_merge / ctrl_merge ≡ 1.0`, so k-means is clustering on `mut/ctrl` ratio alone), so the optimal k may differ.

**Workflow:**
1. Generate elbow plot first (`--elbow-only`)
2. User inspects `cluster/bap1_late/cluster3/elbow_plot/elbow_plot.{png,pdf,svg,jpg}`, picks k
3. Re-run with `--k {chosen_k}` (default=6)

The script writes outputs to `cluster/bap1_late/cluster3/k-{k}/`, so multiple k values can be explored side-by-side.

### 2.3 Cluster 3.0 binary path — **`/usr/local/bin/cluster`** (NOT `~/apps/cluster-1.59/src/cluster`)

`PLAN-CLUSTER.md §3` and §15 both claim the binary is at `~/apps/cluster-1.59/src/cluster`, but `which cluster` resolves to `/usr/local/bin/cluster` and the apps path does not exist. Use the absolute path `/usr/local/bin/cluster` for clarity and to avoid PATH/conda-env-name confusion (the conda environment is also called `cluster`).

Verified CLI flags: `-f filename -r runs -g distance_measure -k k_clusters`.

### 2.4 Cluster naming — **`clust1..clustK`** sorted descending by mean signal

Empirically verified by reading `cluster_tools.sort_clusters` (cluster_tools.py:22-30). It returns `(replace_dict, cluster_order)` where:
- `replace_dict` maps int Cluster-3.0 IDs (1..k) to strings `clust1..clustK`
- Order is **descending mean signal across all `data_cols`** — so `clust1` has the highest mean (gained / mut-up loops, since ctrl_norm≡1 makes mean ∝ mut/ctrl) and `clustK` has the lowest (lost / mut-down loops)
- `cluster_order` is a Python list of these strings in sort order

Example from Popay (`cluster/clustering_example_data/combined-clusters.txt`): clust1=1171, clust2=3164, clust3=3931, clust4=1757, clust5=4124, clust6=2662 — exactly the RAD21-dependency clusters from the paper.

### 2.5 Visualization palette — **two-key dict for both stripped and unstripped labels**

Notebook cell-9 strips `_merge` from `treatment_group` before plotting. Phase 1 verified: `comparison_type(['ctrl_merge','mut_merge']) → ('multiple', ['ctrl','mut'])` and `[c.replace('_merge','') for c in data_cols] → ['ctrl','mut']`. To survive both pre- and post-replacement lookups in `heat()` (uses raw `data_cols`) and `box`/`strip` (use stripped `treatment_group`):

```python
palette = {
    'ctrl_merge': 'darkgrey', 'mut_merge': 'forestgreen',  # for heat()
    'ctrl':       'darkgrey', 'mut':       'forestgreen',  # for box/strip after _merge stripped
}
```

(NOT `ctrl_merged` / `mut_merged` — that was an error in the parent plan. Phase 1 corrections noted this.)

### 2.6 Y-axis ranges — **`Y_range=None` (auto-fit)** instead of Popay's hard-coded `[0,1]` / `[0,0.4]`

Notebook cell-9 hard-codes `Y_range=[0,1]` for line/strip and `[0,0.4]` for box. These were tuned for NIPBL depletion (loop-loss-only biology, ratios <1). **Our BAP1-KO data has both gained AND lost loops** with mut/ctrl ratios spanning 0.39–10.08 (99.9%-tile=1.68). Hard-coded `[0,1]` would clip all gained-loop signal; `[0,0.4]` would clip the entire ctrl distribution (always 1.0 after normalization). Use `Y_range=None` to auto-fit; can be re-tightened after first inspection.

### 2.7 Output directory — **`cluster/bap1_late/cluster3/`**

Notebook auto-derives output path from the input file's directory (`os.path.dirname(loop_count_file) + '/' + file_name + '/cluster3'` → `cluster/data/late_merged_loop_counts/cluster3/`). Our convention puts all Phase 2-6 outputs under `cluster/bap1_late/`, so the script overrides this with an explicit `OUT_DIR` constant pointing at `cluster/bap1_late/cluster3/`. The `cluster/bap1_late/cluster3/` directory was already created by Phase 1's `run_phase1.sh` mkdir.

### 2.8 Multi-format figure output — **PDF + SVG + JPG + PNG, each figure in its own subfolder**

Mirror the project utility `scripts/utils/multi_format_output.R` (the `save_multiformat_*` family with `use_subfolders=TRUE`). Popay's `plotting.py` and `cluster_tools.py` only emit `.png` / `.pdf` and dump everything into a flat `out_dir`. Phase 3 (and all later phases) need:
- **PDF** — publication standard (vector, already produced by Popay)
- **SVG** — Illustrator-editable vector (NEW)
- **JPG** — Google Slides / presentation use (NEW)
- **PNG** — kept as-is (Popay-native; cheap to retain)
- **Subfolder per figure** — `{fig_dir}/{name}/{name}.{pdf,svg,jpg,png}` instead of `{fig_dir}/{name}.{png,pdf}`

**Strategy: a context-manager monkey-patch on `plt.savefig`**, not modifying Popay's modules:
- During the patched scope, any `plt.savefig('foo.png', ...)` call ALSO writes `foo.svg`, `foo.jpg`, and `foo.pdf` from the same active figure (so we capture the same render, not a re-rendering)
- The `.pdf` Popay subsequently writes overwrites our `.pdf` (identical vector output, harmless)
- For functions that only write PNG (like `cluster_tools.elbow`), the wrapper still produces all four formats — exactly what we want
- Subfolder is pre-created and passed AS `out_dir` to the Popay function, so Popay's existing `out_dir + '/' + out_name + '.png'` string concat lands files inside the subfolder without modification

This keeps Popay's modules untouched (Phase 0 contract preserved) while delivering the multi-format / subfolder layout that matches the project's R utility and downstream workflow.

---

## 3. Implementation

### 3.1 New file: `cluster/scripts/utils/multi_format_output.py`

Python equivalent of `scripts/utils/multi_format_output.R`. Provides one context manager and one helper, both reusable by Phase 4-5 scripts:

```python
"""Multi-format plot output utility.

Python equivalent of scripts/utils/multi_format_output.R. Augments matplotlib's
plt.savefig() to emit PDF + SVG + JPG (alongside the original PNG) when used
inside the context manager, and pre-creates a subfolder per figure so output
matches the R utility's use_subfolders=TRUE layout.

Usage:
    from multi_format_output import multi_format_savefig, figure_subfolder

    fig_dir = figure_subfolder('cluster/bap1_late/figures', 'heatmap')
    with multi_format_savefig():
        plotting.heat(out_dir=str(fig_dir), out_name='heatmap', ...)
    # produces: cluster/bap1_late/figures/heatmap/heatmap.{png,pdf,svg,jpg}
"""
import os
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt


@contextmanager
def multi_format_savefig():
    """Inside this scope, plt.savefig('foo.png') ALSO writes foo.{svg,jpg,pdf}.

    Triggers only on .png (the leading format Popay's plotting.py writes), then
    re-savefigs the SAME active figure as SVG (vector, dpi-independent), JPG
    (raster @ same dpi), and PDF (vector). Popay's subsequent .pdf savefig
    overwrites ours with identical vector output (harmless).

    For functions that emit only .png (e.g. cluster_tools.elbow), the wrapper
    still yields all four formats — matching the R utility's PDF+SVG+JPG set.
    Original savefig is restored on exit.
    """
    original = plt.savefig

    def patched(fname, *args, **kwargs):
        base, ext = os.path.splitext(str(fname))
        result = original(fname, *args, **kwargs)
        if ext.lower() == '.png':
            vector_kwargs = {k: v for k, v in kwargs.items() if k != 'dpi'}
            original(base + '.svg', *args, **vector_kwargs)
            original(base + '.pdf', *args, **vector_kwargs)
            jpg_kwargs = dict(kwargs)
            jpg_kwargs.setdefault('dpi', 300)
            original(base + '.jpg', *args, **jpg_kwargs)
        return result

    plt.savefig = patched
    try:
        yield
    finally:
        plt.savefig = original


def figure_subfolder(parent_dir, name) -> Path:
    """Mirror multi_format_output.R use_subfolders=TRUE: create and return
    {parent_dir}/{name}/. Caller passes this as the plotting function's
    out_dir so files land at {parent_dir}/{name}/{name}.{ext}.
    """
    sub = Path(parent_dir) / name
    sub.mkdir(parents=True, exist_ok=True)
    return sub
```

### 3.2 New file: `cluster/scripts/04_clustering.py`

A single self-contained script with two modes (`--elbow-only` for k selection, full run otherwise). Modular function structure:

```python
#!/usr/bin/env python3
"""K-means clustering of nonredundant Hi-C loops via Cluster 3.0.

Adapted from Popay HiC_cluster3.ipynb (cells 1-9) for BAP1-KO mouse cerebellum.
Per-resolution percentile filter, ctrl_merge normalization, sort_clusters by
descending mean signal, four visualization types matching Popay Fig 1d.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR))         # for cluster_tools, plotting
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))# for multi_format_output

from cluster_tools import elbow, sort_clusters, comparison_type
from plotting import heat, line, box, strip, randomize_hex
from multi_format_output import multi_format_savefig, figure_subfolder

# -------- config (paths & constants) --------
LOOP_COUNT_FILE = REPO_ROOT / 'cluster/data/late_merged_loop_counts.txt'
METADATA_FILE   = REPO_ROOT / 'cluster/data/late_merged_loop_metadata.tsv'
OUT_DIR         = REPO_ROOT / 'cluster/bap1_late/cluster3'
CLUSTER_BIN     = '/usr/local/bin/cluster'      # confirmed via `which cluster`
NORMALIZE_COL   = 'ctrl_merge'
DATA_COLS       = ['ctrl_merge', 'mut_merge']   # explicit; matches Phase 1 output schema
COORD_COLS      = ['chr1','x1','x2','chr2','y1','y2']

PALETTE = {
    'ctrl_merge': 'darkgrey', 'mut_merge': 'forestgreen',
    'ctrl':       'darkgrey', 'mut':       'forestgreen',
}

# -------- core --------
def load_and_filter(filter_pct: float) -> pd.DataFrame:
    """Load count file + metadata; apply per-resolution percentile filter on ctrl_merge."""
    df   = pd.read_csv(LOOP_COUNT_FILE, sep='\t')
    meta = pd.read_csv(METADATA_FILE, sep='\t', usecols=COORD_COLS + ['resolution_kb'])
    df   = df.merge(meta, on=COORD_COLS, how='left')
    assert df['resolution_kb'].notna().all(), 'metadata join lost rows'

    thresholds = df.groupby('resolution_kb')[NORMALIZE_COL].quantile(filter_pct).to_dict()
    keep = df.apply(lambda r: r[NORMALIZE_COL] >= thresholds[r['resolution_kb']], axis=1)

    print(f'\nPer-resolution filter: keep {NORMALIZE_COL} >= {filter_pct*100:.1f}%-tile')
    for res in sorted(thresholds):
        n_total = (df['resolution_kb'] == res).sum()
        n_drop  = ((df['resolution_kb'] == res) & ~keep).sum()
        print(f'  res={res:>2}kb  thr={thresholds[res]:8.2f}  '
              f'drop={n_drop:>4}  keep={n_total-n_drop:>5} ({(n_total-n_drop)/n_total*100:.1f}%)')
    print(f'  TOTAL  keep {keep.sum()} of {len(df)} ({keep.sum()/len(df)*100:.1f}%)\n')

    return df[keep].drop(columns=['resolution_kb']).reset_index(drop=True)


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Divide DATA_COLS by NORMALIZE_COL; ctrl→1.0, mut→mut/ctrl. NaN→0."""
    out = df.copy()
    out[DATA_COLS] = out[DATA_COLS].div(out[NORMALIZE_COL], axis=0).fillna(0)
    return out


def make_id(df: pd.DataFrame) -> pd.Series:
    """Composite row id: chr1-x1-x2-chr2-y1-y2 (matches notebook cell-5)."""
    return (df['chr1'] + '-' + df['x1'].astype(str) + '-' + df['x2'].astype(str) + '-'
            + df['chr2'] + '-' + df['y1'].astype(str) + '-' + df['y2'].astype(str))


def run_elbow(filter_pct: float) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_norm = normalize(load_and_filter(filter_pct))
    # subfolder so all 4 formats land in cluster/bap1_late/cluster3/elbow_plot/
    elbow_dir = figure_subfolder(OUT_DIR, 'elbow_plot')
    with multi_format_savefig():
        elbow(out_dir=str(elbow_dir), count_matrix=df_norm, data_cols=DATA_COLS)
    print(f'Elbow plot: {elbow_dir}/elbow_plot.{{png,pdf,svg,jpg}}')


def run_clustering(k: int, filter_pct: float) -> None:
    df       = load_and_filter(filter_pct)
    df_norm  = normalize(df)
    df_norm['id'] = make_id(df_norm)

    k_dir    = OUT_DIR / f'k-{k}'
    data_dir = k_dir  / 'data'
    fig_dir  = k_dir  / 'figures'
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = data_dir / 'input_matrix.txt'
    df_norm[['id'] + DATA_COLS].to_csv(matrix_path, sep='\t', index=False, lineterminator='\n')

    cmd = [CLUSTER_BIN, '-f', str(matrix_path), '-r', '100', '-g', '7', '-k', str(k)]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    kgg = pd.read_csv(data_dir / f'input_matrix_K_G{k}.kgg', sep='\t')
    df['id'] = make_id(df)                                 # rebuild id on un-normalized df
    clusters_df = kgg.merge(df, on='id').drop(columns=['id'])

    replace_dict, cluster_order = sort_clusters(clusters_df, DATA_COLS)
    for orig, new in replace_dict.items():
        clusters_df['GROUP'] = clusters_df['GROUP'].replace(orig, new)

    out_path = data_dir / 'combined-clusters.txt'
    clusters_df.to_csv(out_path, sep='\t', index=False, lineterminator='\n')
    print(f'Wrote {out_path}')

    print('\nCluster sizes (descending mean signal — clust1=highest mut/ctrl):')
    for c, n in clusters_df['GROUP'].value_counts().reindex(cluster_order).items():
        print(f'  {c:>8}: {n:>5} ({n/len(clusters_df)*100:.1f}%)')

    visualize(clusters_df, cluster_order, fig_dir)


def visualize(clusters_df: pd.DataFrame, cluster_order: list, fig_dir: Path) -> None:
    """Heat / line / strip / box per notebook cell-9, comparison='multiple' branch."""
    clusters_df = clusters_df.sort_values('GROUP').copy()

    # max-normalize per row for heatmap
    clusters_norm = clusters_df.copy()
    clusters_norm['_max']     = clusters_norm[DATA_COLS].max(axis=1)
    clusters_norm[DATA_COLS]  = clusters_norm[DATA_COLS].div(clusters_norm['_max'], axis=0).fillna(0)
    clusters_norm = clusters_norm.drop(columns=['_max'])

    info_cols = [c for c in clusters_norm.columns if c not in DATA_COLS]
    stacked_df   = clusters_df.melt(id_vars=info_cols, value_vars=DATA_COLS,
                                    var_name='treatment_group', value_name='balanced counts')
    stacked_norm = clusters_norm.melt(id_vars=info_cols, value_vars=DATA_COLS,
                                      var_name='treatment_group',
                                      value_name='max-normalized balanced counts')
    stacked_df['treatment_group']   = stacked_df['treatment_group'].str.replace('_merge', '')
    stacked_norm['treatment_group'] = stacked_norm['treatment_group'].str.replace('_merge', '')

    line_palette = {c: randomize_hex() for c in cluster_order}

    # Each plot in its own subfolder; multi_format_savefig adds svg/jpg/pdf siblings.
    with multi_format_savefig():
        heatmap_dir = figure_subfolder(fig_dir, 'heatmap')
        heat(count_matrix=clusters_norm, data_cols=DATA_COLS, out_dir=str(heatmap_dir),
             measure='max-normalized balanced counts', data_type='chromatin loops',
             title='', order=cluster_order, subplot_col='GROUP', subplots=True,
             out_name='heatmap', proportional=False, add_n=True)

        lineplot_dir = figure_subfolder(fig_dir, 'lineplot')
        line(melted_df=stacked_norm, ycol='max-normalized balanced counts',
             xcol='treatment_group', hue_col='GROUP', out_dir=str(lineplot_dir),
             ycol_measure='max-normalized balanced counts', xcol_measure=None,
             title='', out_name='lineplot', palette=line_palette,
             subplots=False, subplot_col=None, sharey=False, Y_range=None,
             add_n='hue_col', order=None)

        stripplot_dir = figure_subfolder(fig_dir, 'stripplot')
        strip(melted_df=stacked_df, ycol='balanced counts', xcol='GROUP',
              out_dir=str(stripplot_dir), measure='balanced counts', title='',
              out_name='stripplot', sharey=False, Y_range=None,
              hue_col='treatment_group', subplots=False, subplot_col=None,
              palette=PALETTE, add_n='xcol', order=None)

        boxplot_dir = figure_subfolder(fig_dir, 'boxplot')
        box(melted_df=stacked_df, ycol='balanced counts', xcol='GROUP',
            out_dir=str(boxplot_dir), measure='balanced counts', title='',
            out_name='boxplot', sharey=False, Y_range=None,
            hue_col='treatment_group', subplots=False, subplot_col=None,
            palette=PALETTE, add_n='xcol', order=None)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--elbow-only', action='store_true', help='Generate elbow plot only')
    p.add_argument('--k', type=int, default=6, help='k for k-means (default 6)')
    p.add_argument('--filter-pct', type=float, default=0.01,
                   help='Per-resolution percentile floor for ctrl_merge (default 0.01)')
    args = p.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.elbow_only:
        run_elbow(args.filter_pct)
    else:
        run_clustering(args.k, args.filter_pct)


if __name__ == '__main__':
    main()
```

### 3.3 New file: `cluster/scripts/run_phase3.sh`

Driver script — runs elbow first, prints inspection prompt, then runs full clustering. Per `feedback_user_preferences` memory, **does NOT use `set -euo pipefail`** (use plain `set -e` only).

```bash
#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."   # repo root

PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=cluster/scripts/04_clustering.py
LOG=cluster/phase3.txt

K=${1:-6}
FILTER_PCT=${2:-0.01}

{
  echo "============================================================"
  echo "Phase 3: K-means clustering"
  echo "Repo root: $(pwd)"
  echo "Started:   $(date)"
  echo "k=${K}  filter_pct=${FILTER_PCT}"
  echo "============================================================"

  echo
  echo "[1/2] Generating elbow plot (k=1..14)..."
  "$PYTHON" "$SCRIPT" --elbow-only --filter-pct "$FILTER_PCT"

  echo
  echo "[2/2] Running clustering with k=${K}..."
  "$PYTHON" "$SCRIPT" --k "$K" --filter-pct "$FILTER_PCT"

  echo
  echo "============================================================"
  echo "Phase 3 outputs"
  echo "============================================================"
  ls -lh cluster/bap1_late/cluster3/elbow_plot/ 2>/dev/null || true
  ls -lh "cluster/bap1_late/cluster3/k-${K}/data/combined-clusters.txt" 2>/dev/null || true
  for sub in heatmap lineplot stripplot boxplot; do
    ls -lh "cluster/bap1_late/cluster3/k-${K}/figures/${sub}/" 2>/dev/null || true
  done

  echo "Phase 3 finished: $(date)"
} 2>&1 | tee "$LOG"
```

Usage:
- Default (k=6, 1%-tile filter): `bash cluster/scripts/run_phase3.sh`
- Override k: `bash cluster/scripts/run_phase3.sh 5`
- Override both: `bash cluster/scripts/run_phase3.sh 6 0.05`

---

## 4. Critical files

| Path | Status | Role in Phase 3 |
|------|--------|------------------|
| `cluster/data/late_merged_loop_counts.txt` | exists (Phase 1) | input — 39,344 loops, 8 cols |
| `cluster/data/late_merged_loop_metadata.tsv` | exists (Phase 1) | input — provides `resolution_kb` for per-res filter |
| `cluster/cluster_tools.py` | exists, Phase 0 fixes verified | imports `elbow`, `sort_clusters`, `comparison_type` |
| `cluster/plotting.py` | exists, Phase 0 fixes verified | imports `heat`, `line`, `box`, `strip`, `randomize_hex` |
| `cluster/custom_params.json` | exists | matplotlib styling (loaded by `plotting.py`) |
| `/usr/local/bin/cluster` | exists, executable | Cluster 3.0 binary (`-f -r -g -k`) |
| `/opt/homebrew/anaconda3/envs/cluster/bin/python3` | exists | pandas 1.5.3, numpy 1.24.4, seaborn 0.13.2, matplotlib 3.7.5, sklearn 1.3.2 |
| `cluster/bap1_late/cluster3/` | exists (empty, Phase 1 mkdir) | output root |
| `scripts/utils/multi_format_output.R` | exists (reference; not invoked from Python) | format/subfolder conventions to mirror |
| `cluster/scripts/utils/multi_format_output.py` | **NEW (this phase)** | Python equivalent — context-manager savefig wrapper + subfolder helper |
| `cluster/scripts/04_clustering.py` | **NEW (this phase)** | clustering script |
| `cluster/scripts/run_phase3.sh` | **NEW (this phase)** | driver |
| `cluster/phase3.txt` | **NEW (this phase)** | execution log |

---

## 5. Verification (run after `bash cluster/scripts/run_phase3.sh`)

### 5.1 Filter and clustering produced expected counts
```bash
# combined-clusters.txt should have ~38,950 rows (header + 38,950) at 1%-tile filter
wc -l cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt
# expect: 38951

# Cluster size distribution: all clusters should have >1% of rows (>389 loops)
awk 'NR>1 {print $1}' cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt \
    | sort | uniq -c | sort -k2,2
# expect: 6 lines, each count > 389; ordered clust1..clust6
```

### 5.2 Elbow plot is interpretable
- Open `cluster/bap1_late/cluster3/elbow_plot/elbow_plot.png` (or .svg/.pdf for higher fidelity)
- Look for visible bend; confirm k=6 is at or beyond the bend (Popay's heuristic)

### 5.3 Visualizations rendered (multi-format, subfolder layout)
```bash
ls cluster/bap1_late/cluster3/elbow_plot/
# expect: elbow_plot.{png,pdf,svg,jpg}

for sub in heatmap lineplot stripplot boxplot; do
  echo "=== $sub ==="
  ls cluster/bap1_late/cluster3/k-6/figures/${sub}/
done
# expect: each subfolder contains {sub}.{png,pdf,svg,jpg}  (4 files per figure)
```

**Total figure files for k=6 run:** 4 plots × 4 formats + elbow × 4 formats = **20 image files**.
SVG sizes for heatmap/stripplot may be 10-50 MB (one vector element per data point × ~38k loops). Acceptable for an Illustrator-edit workflow.

### 5.4 Cluster ordering is monotonic in mean signal
```bash
/opt/homebrew/anaconda3/envs/cluster/bin/python3 -c "
import pandas as pd
df = pd.read_csv('cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt', sep='\t')
print(df.groupby('GROUP')[['ctrl_merge','mut_merge']].mean().sort_index())
print('mean per cluster (sum of ctrl+mut):')
print(df.groupby('GROUP')[['ctrl_merge','mut_merge']].mean().sum(axis=1).sort_values(ascending=False))
"
# expect: clust1 has highest mean, monotonically decreasing through clustK
# expect: mut/ctrl ratio decreases monotonically clust1 -> clustK
```

### 5.5 Sanity check vs. metadata sidecar (cross-check biology)
Loops with up_in_mutant direction (gained) should concentrate in early (high-mean) clusters; down_in_mutant (lost) in late (low-mean) clusters. This is a soft check — significant loops are only ~20% of the dataset.
```bash
/opt/homebrew/anaconda3/envs/cluster/bin/python3 -c "
import pandas as pd
clusters = pd.read_csv('cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt', sep='\t')
meta = pd.read_csv('cluster/data/late_merged_loop_metadata.tsv', sep='\t')
m = clusters.merge(meta[['chr1','x1','x2','chr2','y1','y2','direction']],
                   on=['chr1','x1','x2','chr2','y1','y2'])
print(pd.crosstab(m['GROUP'], m['direction'], normalize='index').round(3))
"
# expect: clust1 enriched for up_in_mutant; clustK enriched for down_in_mutant
```

---

## 6. Risks & known unknowns

1. **Effective 1D clustering.** With only 2 conditions and ctrl-normalization (ctrl_norm≡1), k-means is clustering on a 1D axis (mut/ctrl ratio). This is mathematically valid and produces interpretable groups, but downstream interpretations should remember each cluster ≈ a fold-change band, not an independently-discovered chromatin state. The biological meaning emerges in Phase 4 when we overlay ChromHMM, ChIP-seq, and gene annotation per cluster.

2. **Outlier ratios driving cluster centroids.** Max mut/ctrl ratio is 10.08; the 1%-tile filter is the primary mitigation. If a cluster ends up with <50 loops (mostly outliers), rerun with `--filter-pct 0.05` or higher.

3. **Cluster 3.0 non-determinism.** k-means with `-r 100` picks the best of 100 random initializations but is not seeded. Reruns may produce slightly different cluster assignments at the boundaries. The 100-run consensus minimizes this; cluster sizes should be stable to ~5%.

4. **Phase 4 dependency on `combined-clusters.txt` schema.** Phase 4.4-4.8 all read this file via `group_dict = clusters.groupby('GROUP')`. Confirmed format matches Popay's example exactly: `GROUP \t chr1 \t x1 \t x2 \t chr2 \t y1 \t y2 \t ctrl_merge \t mut_merge`.

5. **Y_range deferred.** Set to None initially (auto-fit). After visual inspection, if box/strip plots have too much whitespace from the long tail, can re-tune to e.g. `[0, 1.5]` per condition. Not blocking.

6. **Out of scope for this phase:**
   - Phase 4 downstream analyses (loop size, classification, ChromHMM enrichment, etc.)
   - Choosing the final k — this plan generates elbow + a default k=6 run; user reviews and re-runs at chosen k if different
   - The user must inspect the elbow plot to confirm k before treating clusters as final
