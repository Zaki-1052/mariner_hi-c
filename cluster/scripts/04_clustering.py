#!/usr/bin/env python3
# cluster/scripts/04_clustering.py
"""K-means clustering of nonredundant Hi-C loops via Cluster 3.0.

Adapted from Popay HiC_cluster3.ipynb (cells 1-9) for BAP1-KO mouse cerebellum.
Per-resolution percentile filter on ctrl_merge, ctrl_merge normalization,
sort_clusters by descending mean signal, four visualizations matching Popay
Fig 1d. Each figure written in its own subfolder as PNG + PDF + SVG + JPG via
the multi_format_savefig context manager.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR))           # cluster_tools, plotting
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))  # multi_format_output

from cluster_tools import elbow, sort_clusters
from plotting import heat, line, box, strip, randomize_hex
from multi_format_output import multi_format_savefig, figure_subfolder

# -------- config --------
LOOP_COUNT_FILE = REPO_ROOT / 'cluster/data/late_merged_loop_counts.txt'
METADATA_FILE   = REPO_ROOT / 'cluster/data/late_merged_loop_metadata.tsv'
OUT_DIR         = REPO_ROOT / 'cluster/bap1_late/cluster3'
CLUSTER_BIN     = '/usr/local/bin/cluster'
NORMALIZE_COL   = 'ctrl_merge'
DATA_COLS       = ['ctrl_merge', 'mut_merge']
COORD_COLS      = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']

PALETTE = {
    'ctrl_merge': 'darkgrey', 'mut_merge': 'forestgreen',
    'ctrl':       'darkgrey', 'mut':       'forestgreen',
}


def load_and_filter(filter_pct: float) -> pd.DataFrame:
    # Per-resolution percentile filter on ctrl_merge — addresses the ~6x scale
    # difference between 5/10/25 kb count distributions.
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
              f'drop={n_drop:>4}  keep={n_total - n_drop:>5} '
              f'({(n_total - n_drop) / n_total * 100:.1f}%)')
    print(f'  TOTAL  keep {keep.sum()} of {len(df)} '
          f'({keep.sum() / len(df) * 100:.1f}%)\n')

    return df[keep].drop(columns=['resolution_kb']).reset_index(drop=True)


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    # Divide DATA_COLS by NORMALIZE_COL: ctrl -> 1.0, mut -> mut/ctrl
    out = df.copy()
    out[DATA_COLS] = out[DATA_COLS].div(out[NORMALIZE_COL], axis=0).fillna(0)
    return out


def make_id(df: pd.DataFrame) -> pd.Series:
    # Composite row id chr1-x1-x2-chr2-y1-y2 for joining .kgg back to coords
    return (df['chr1'] + '-' + df['x1'].astype(str) + '-' + df['x2'].astype(str)
            + '-' + df['chr2'] + '-' + df['y1'].astype(str) + '-' + df['y2'].astype(str))


def run_elbow(filter_pct: float) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df_norm = normalize(load_and_filter(filter_pct))
    elbow_dir = figure_subfolder(OUT_DIR, 'elbow_plot')
    with multi_format_savefig():
        elbow(out_dir=str(elbow_dir), count_matrix=df_norm, data_cols=DATA_COLS)
    print(f'Elbow plot: {elbow_dir}/elbow_plot.{{png,pdf,svg,jpg}}')


def run_clustering(k: int, filter_pct: float) -> None:
    df       = load_and_filter(filter_pct)
    df_norm  = normalize(df)
    df_norm['id'] = make_id(df_norm)

    k_dir    = OUT_DIR / f'k-{k}'
    data_dir = k_dir / 'data'
    fig_dir  = k_dir / 'figures'
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = data_dir / 'input_matrix.txt'
    df_norm[['id'] + DATA_COLS].to_csv(matrix_path, sep='\t', index=False,
                                       lineterminator='\n')

    cmd = [CLUSTER_BIN, '-f', str(matrix_path), '-r', '100', '-g', '7', '-k', str(k)]
    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)

    kgg_path = data_dir / f'input_matrix_K_G{k}.kgg'
    kgg = pd.read_csv(kgg_path, sep='\t')

    df['id'] = make_id(df)
    clusters_df = kgg.merge(df, on='id').drop(columns=['id'])

    replace_dict, cluster_order = sort_clusters(clusters_df, DATA_COLS)
    for orig, new in replace_dict.items():
        clusters_df['GROUP'] = clusters_df['GROUP'].replace(orig, new)

    out_path = data_dir / 'combined-clusters.txt'
    clusters_df.to_csv(out_path, sep='\t', index=False, lineterminator='\n')
    print(f'Wrote {out_path}')

    print('\nCluster sizes (descending mean signal — clust1=highest mut/ctrl):')
    counts = clusters_df['GROUP'].value_counts().reindex(cluster_order)
    for c, n in counts.items():
        print(f'  {c:>8}: {n:>5} ({n / len(clusters_df) * 100:.1f}%)')

    visualize(clusters_df, cluster_order, fig_dir)


def visualize(clusters_df: pd.DataFrame, cluster_order: list, fig_dir: Path) -> None:
    # Heat / line / strip / box per HiC_cluster3.ipynb cell-9 'multiple' branch.
    # Each plot wrapped in its own subfolder via figure_subfolder + multi_format_savefig.
    clusters_df = clusters_df.sort_values('GROUP').copy()

    # Max-normalize per row for heatmap
    clusters_norm = clusters_df.copy()
    clusters_norm['_max']    = clusters_norm[DATA_COLS].max(axis=1)
    clusters_norm[DATA_COLS] = clusters_norm[DATA_COLS].div(
        clusters_norm['_max'], axis=0).fillna(0)
    clusters_norm = clusters_norm.drop(columns=['_max'])

    info_cols = [c for c in clusters_norm.columns if c not in DATA_COLS]
    stacked_df = clusters_df.melt(
        id_vars=info_cols, value_vars=DATA_COLS,
        var_name='treatment_group', value_name='balanced counts')
    stacked_norm = clusters_norm.melt(
        id_vars=info_cols, value_vars=DATA_COLS,
        var_name='treatment_group',
        value_name='max-normalized balanced counts')
    stacked_df['treatment_group']   = stacked_df['treatment_group'].str.replace('_merge', '')
    stacked_norm['treatment_group'] = stacked_norm['treatment_group'].str.replace('_merge', '')

    line_palette = {c: randomize_hex() for c in cluster_order}

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
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--elbow-only', action='store_true',
                   help='Generate elbow plot only (k=1..14)')
    p.add_argument('--k', type=int, default=6,
                   help='k for k-means clustering (default 6)')
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
