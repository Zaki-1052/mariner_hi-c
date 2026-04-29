# cluster/scripts/11_histone_anchors_metagene.py
"""Regenerate Phase 5 histone anchor metagene as clean profile figure.

Reads the existing computeMatrix output (histone_anchors_values, 1.6GB)
and renders a 4-row x 6-column grid of per-cluster ctrl/mut profile lines
for all 4 histone marks (H3K27ac, H3K27me3, H2AK119ub, H3K27me1).

Runtime: ~1-2 min (dominated by loading the matrix).
"""
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))
from multi_format_output import multi_format_savefig  # noqa: E402

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as f:
    plt.rcParams.update(json.load(f))

_out = os.environ.get('CLUSTER_OUT_DIR', 'outputs/bap1_late')
SDSC_DIR    = Path(os.environ.get('CLUSTER_SDSC_DIR', '/Users/zakiralibhai/sdsc'))
VALUES_FILE = SDSC_DIR / 'histone_anchors_values'
OUT_DIR     = CLUSTER_DIR / _out / 'figures/deeptools/histone_anchors'

BIO_ORDER = ['clust6', 'clust3', 'clust1', 'clust2', 'clust4', 'clust5']
BIO_LABELS = {
    'clust6': 'clust6\nloss (78%)',
    'clust3': 'clust3\nmod loss',
    'clust1': 'clust1\nunchanged',
    'clust2': 'clust2\n~unchanged',
    'clust4': 'clust4\nmod gain',
    'clust5': 'clust5\ngain (97%)',
}

CTRL_COLOR = '#2166ac'
MUT_COLOR  = '#b2182b'

MARK_ROWS = [
    ('H3K27ac',   [('H3K27ac_ctrl', 0),   ('H3K27ac_mut', 1)]),
    ('H3K27me3',  [('H3K27me3_ctrl', 2),  ('H3K27me3_mut', 3)]),
    ('H2AK119ub', [('H2AK119ub_ctrl', 4), ('H2AK119ub_mut', 5)]),
    ('H3K27me1',  [('H3K27me1_ctrl', 6),  ('H3K27me1_mut', 7)]),
]

BIGWIG_INDICES = {bw: idx for _, pairs in MARK_ROWS for bw, idx in pairs}

MARK_SHORT = {
    'H3K27ac': 'K27ac', 'H3K27me3': 'K27me3',
    'H2AK119ub': 'K119ub', 'H3K27me1': 'K27me1',
}


def parse_header(filepath):
    # type: (Path) -> Dict
    with open(filepath) as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()

    tokens1 = [t.replace('#', '') for t in line1.split('\t') if ':' in t]
    cluster_names = [t.split(':')[0].replace('_anchors.bed', '') for t in tokens1]
    cluster_sizes = [int(t.split(':')[1]) for t in tokens1]

    params = {}
    for p in line2.replace('#', '').split('\t'):
        if ':' in p:
            k, v = p.split(':', 1)
            params[k.strip()] = int(v.strip())
    bin_size = params['bin size']
    n_bins = (params['upstream'] + params['downstream']) // bin_size

    group_set = set(tokens1)
    bigwig_labels = [t for t in line3.split('\t') if t not in group_set]
    bigwig_order = list(dict.fromkeys(bigwig_labels))

    return {
        'cluster_names': cluster_names,
        'cluster_sizes': cluster_sizes,
        'bin_size': bin_size,
        'n_bins': n_bins,
        'bigwig_order': bigwig_order,
    }


def compute_profiles(filepath, header, bigwig_indices):
    # type: (Path, Dict, Dict[str, int]) -> Dict[Tuple[str, str], Dict]
    n_bins = header['n_bins']
    total = sum(header['cluster_sizes'])
    n_cols = len(bigwig_indices) * n_bins
    print('Loading {} x {} matrix...'.format(total, n_cols))
    data = pd.read_csv(filepath, sep='\t', header=None, skiprows=3)
    arr = data.select_dtypes(include=[np.number]).values.astype(np.float32)
    print('Loaded: {}'.format(arr.shape))

    labels = np.concatenate([[name] * size for name, size in
                             zip(header['cluster_names'], header['cluster_sizes'])])[:arr.shape[0]]

    result = {}  # type: Dict[Tuple[str, str], Dict]
    for bw_name, bw_idx in bigwig_indices.items():
        cols = slice(bw_idx * n_bins, (bw_idx + 1) * n_bins)
        for cluster in header['cluster_names']:
            mask = labels == cluster
            rows = arr[mask][:, cols]
            n = mask.sum()
            mean = np.nanmean(rows, axis=0)
            sem = np.nanstd(rows, axis=0) / np.sqrt(n)
            result[(bw_name, cluster)] = {'mean': mean, 'sem': sem, 'n': n}
    return result


def render_figure(profiles, header, out_dir):
    # type: (Dict, Dict, Path) -> None
    n_bins = header['n_bins']
    x = np.arange(n_bins) * header['bin_size'] - 5000 + header['bin_size'] // 2

    row_ymaxes = {}  # type: Dict[str, float]
    for mark_name, bw_pairs in MARK_ROWS:
        ymax = 0.0
        for bw_name, _ in bw_pairs:
            for cluster in BIO_ORDER:
                peak = float(profiles[(bw_name, cluster)]['mean'].max())
                if peak > ymax:
                    ymax = peak
        row_ymaxes[mark_name] = ymax * 1.12

    fig = plt.figure(figsize=(9.5, 8.5))
    gs = GridSpec(4, 6, figure=fig,
                  hspace=0.40, wspace=0.12,
                  left=0.07, right=0.97, top=0.94, bottom=0.06)

    for row_idx, (mark_name, bw_pairs) in enumerate(MARK_ROWS):
        y_max = row_ymaxes[mark_name]

        for col_idx, cluster in enumerate(BIO_ORDER):
            ax = fig.add_subplot(gs[row_idx, col_idx])

            ax.axvspan(-5000, 0, alpha=0.06, color='#4393c3', zorder=0)
            ax.axvspan(0, 5000, alpha=0.06, color='#f4a582', zorder=0)
            ax.axvline(0, color='#888888', linewidth=0.3, linestyle='--', zorder=1)

            for bw_name, _ in bw_pairs:
                is_ctrl = bw_name.endswith('_ctrl')
                color = CTRL_COLOR if is_ctrl else MUT_COLOR
                p = profiles[(bw_name, cluster)]
                ax.plot(x, p['mean'], color=color, linewidth=0.6, zorder=2)
                ax.fill_between(x, p['mean'] - p['sem'], p['mean'] + p['sem'],
                                color=color, alpha=0.12, zorder=2)

            ax.set_xlim(-5000, 5000)
            ax.set_ylim(0, y_max)
            ax.tick_params(axis='both', length=2, pad=1)

            if row_idx == 0:
                ax.set_title(BIO_LABELS[cluster], fontsize=5.5, pad=2, linespacing=1.1)
                n_anchors = profiles[(bw_pairs[0][0], cluster)]['n']
                ax.text(0.97, 0.97, 'n={:,}'.format(n_anchors),
                        transform=ax.transAxes, fontsize=3.8, ha='right', va='top',
                        family='monospace',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  alpha=0.7, linewidth=0.3))

            if col_idx == 0:
                ax.set_ylabel('{} (RPKM)'.format(mark_name), fontsize=6)
            else:
                ax.yaxis.set_ticklabels([])

            if row_idx == 3:
                ax.set_xticks([-5000, 0, 5000])
                ax.set_xticklabels(['-5kb', 'anchor', '+5kb'], fontsize=4.5)
            else:
                ax.tick_params(labelbottom=False)

            if row_idx == 0 and col_idx == 5:
                ax.plot([], [], color=CTRL_COLOR, linewidth=0.8, label='ctrl')
                ax.plot([], [], color=MUT_COLOR, linewidth=0.8, label='mut')
                ax.legend(fontsize=4.5, loc='upper left', frameon=False,
                          handlelength=1.2, handletextpad=0.3)

    fig.suptitle('Histone mark profiles at loop anchors (per cluster)',
                 fontsize=9, fontweight='bold', y=0.98)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / 'histone_anchors_metagene.png'
    with multi_format_savefig():
        fig.savefig(str(out_path), dpi=300)
    plt.close(fig)
    print('Saved: {}.{{png,pdf,svg,jpg}}'.format(out_path.stem))


def main():
    # type: () -> None
    if not VALUES_FILE.exists():
        raise FileNotFoundError('Values file not found: {}'.format(VALUES_FILE))

    print('[1/3] Parsing header...')
    header = parse_header(VALUES_FILE)
    print('  Clusters: {}'.format(
        ', '.join('{} (n={})'.format(c, s)
                  for c, s in zip(header['cluster_names'], header['cluster_sizes']))))
    print('  BigWigs:  {}'.format(header['bigwig_order']))
    print('  Bins:     {} x {}bp = {}bp window'.format(
        header['n_bins'], header['bin_size'],
        header['n_bins'] * header['bin_size']))

    print('\n[2/3] Computing profiles...')
    profiles = compute_profiles(VALUES_FILE, header, BIGWIG_INDICES)

    print('\n[3/3] Rendering figure...')
    render_figure(profiles, header, OUT_DIR)

    print('\nDone. Output in: {}'.format(OUT_DIR))


if __name__ == '__main__':
    main()
