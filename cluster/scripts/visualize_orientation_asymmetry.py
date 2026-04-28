#!/usr/bin/env python3
# cluster/scripts/visualize_orientation_asymmetry.py
"""Visualize K27me3 exterior/interior asymmetry from oriented anchor metagene.

Top: per-cluster K27me3 metagene profiles (1x6 grid, biological order).
Bottom: Ext/Int ratio bar chart with significance.
"""
import json
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))
from multi_format_output import multi_format_savefig  # noqa: E402

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as f:
    plt.rcParams.update(json.load(f))

VALUES_FILE = CLUSTER_DIR / 'bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values'
ASYM_TSV    = CLUSTER_DIR / 'bap1_late/figures/deeptools/oriented_anchors/asymmetry_quantification.tsv'
OUT_DIR     = CLUSTER_DIR / 'bap1_late/figures/deeptools/oriented_anchors'

BIO_ORDER  = ['clust6', 'clust3', 'clust1', 'clust2', 'clust4', 'clust5']
BIO_LABELS = {
    'clust6': 'clust6\nloss (78%)',
    'clust3': 'clust3\nmod loss',
    'clust1': 'clust1\nunchanged',
    'clust2': 'clust2\n~unchanged',
    'clust4': 'clust4\nmod gain',
    'clust5': 'clust5\ngain (97%)',
}
DIRECTION_SHORT = {
    'clust6': 'loss', 'clust3': 'mod loss', 'clust1': 'unchanged',
    'clust2': '~unchanged', 'clust4': 'mod gain', 'clust5': 'gain',
}

CTRL_COLOR = '#2166ac'
MUT_COLOR  = '#b2182b'


def parse_header(filepath):
    # type: (Path) -> Dict
    with open(filepath) as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
    tokens1 = [t.replace('#', '') for t in line1.split('\t') if ':' in t]
    cluster_names = [t.split(':')[0].replace('_oriented_anchors.bed', '') for t in tokens1]
    cluster_sizes = [int(t.split(':')[1]) for t in tokens1]
    params = {}
    for p in line2.replace('#', '').split('\t'):
        if ':' in p:
            k, v = p.split(':', 1)
            params[k.strip()] = int(v.strip())
    bin_size = params['bin size']
    n_bins = (params['upstream'] + params['downstream']) // bin_size
    return {
        'cluster_names': cluster_names,
        'cluster_sizes': cluster_sizes,
        'bin_size': bin_size,
        'n_bins': n_bins,
    }


def compute_profiles(filepath, header, bigwig_indices):
    # type: (Path, Dict, Dict[str, int]) -> Dict[Tuple[str, str], Dict]
    """Compute per-cluster per-bin mean and SEM for selected BigWig columns."""
    n_bins = header['n_bins']
    print(f'Loading {sum(header["cluster_sizes"])} x {8 * n_bins} matrix...')
    data = pd.read_csv(filepath, sep='\t', header=None, skiprows=3)
    arr = data.select_dtypes(include=[np.number]).values.astype(np.float64)
    print(f'Loaded: {arr.shape}')

    labels = np.concatenate([[name] * size for name, size in
                             zip(header['cluster_names'], header['cluster_sizes'])])[:arr.shape[0]]

    result = {}
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


def main():
    # type: () -> None
    header = parse_header(VALUES_FILE)
    n_bins = header['n_bins']
    x = np.arange(n_bins) * header['bin_size'] - 5000 + header['bin_size'] // 2

    bigwig_indices = {'H3K27me3_ctrl': 2, 'H3K27me3_mut': 3}
    profiles = compute_profiles(VALUES_FILE, header, bigwig_indices)

    asym_df = pd.read_csv(ASYM_TSV, sep='\t')

    fig = plt.figure(figsize=(7.5, 4.2))
    gs = GridSpec(2, 6, figure=fig, height_ratios=[1.3, 1.0],
                  hspace=0.55, wspace=0.15,
                  left=0.08, right=0.95, top=0.92, bottom=0.10)

    y_max = max(profiles[(bw, c)]['mean'].max()
                for bw in bigwig_indices for c in header['cluster_names']) * 1.12

    for i, cluster in enumerate(BIO_ORDER):
        ax = fig.add_subplot(gs[0, i])

        ax.axvspan(-5000, 0, alpha=0.06, color='#4393c3', zorder=0)
        ax.axvspan(0, 5000, alpha=0.06, color='#f4a582', zorder=0)
        ax.axvline(0, color='#888888', linewidth=0.3, linestyle='--', zorder=1)

        for bw_name, color in [('H3K27me3_ctrl', CTRL_COLOR), ('H3K27me3_mut', MUT_COLOR)]:
            p = profiles[(bw_name, cluster)]
            ax.plot(x, p['mean'], color=color, linewidth=0.6, zorder=2)
            ax.fill_between(x, p['mean'] - p['sem'], p['mean'] + p['sem'],
                            color=color, alpha=0.12, zorder=2)

        ax.set_xlim(-5000, 5000)
        ax.set_ylim(0, y_max)
        ax.set_title(BIO_LABELS[cluster], fontsize=5.5, pad=2, linespacing=1.1)
        ax.set_xticks([-5000, 0, 5000])
        ax.set_xticklabels(['-5kb\next', '0', '+5kb\nint'], fontsize=4.5)
        ax.tick_params(axis='both', length=2, pad=1)

        row_c = asym_df[(asym_df['mark'] == 'H3K27me3_ctrl') & (asym_df['cluster'] == cluster)]
        row_m = asym_df[(asym_df['mark'] == 'H3K27me3_mut') & (asym_df['cluster'] == cluster)]
        ratio_c = row_c['ext_int_ratio'].values[0]
        pval_c = row_c['wilcoxon_pval'].values[0]
        ratio_m = row_m['ext_int_ratio'].values[0]
        pval_m = row_m['wilcoxon_pval'].values[0]

        def sig(p):
            if pd.isna(p): return ''
            if p < 0.001: return '***'
            if p < 0.01:  return '**'
            if p < 0.05:  return '*'
            return 'ns'

        n_anchors = row_c['n_anchors'].values[0]
        ax.text(0.97, 0.97,
                f'n={n_anchors:,}\n'
                f'ctrl E/I={ratio_c:.2f} {sig(pval_c)}\n'
                f'mut  E/I={ratio_m:.2f} {sig(pval_m)}',
                transform=ax.transAxes, fontsize=3.8, ha='right', va='top',
                family='monospace', linespacing=1.3,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, linewidth=0.3))

        if i == 0:
            ax.set_ylabel('H3K27me3 (RPKM)', fontsize=6)
        else:
            ax.yaxis.set_ticklabels([])

        if i == 5:
            ax.plot([], [], color=CTRL_COLOR, linewidth=0.8, label='ctrl')
            ax.plot([], [], color=MUT_COLOR, linewidth=0.8, label='mut')
            ax.legend(fontsize=4.5, loc='upper left', frameon=False,
                      handlelength=1.2, handletextpad=0.3)

    ax_bar = fig.add_subplot(gs[1, :])

    bar_width = 0.28
    x_pos = np.arange(len(BIO_ORDER))

    for j, (mark, color, label) in enumerate([
        ('H3K27me3_ctrl', CTRL_COLOR, 'K27me3 ctrl'),
        ('H3K27me3_mut', MUT_COLOR, 'K27me3 mut'),
    ]):
        vals = []
        pvals = []
        for cluster in BIO_ORDER:
            row = asym_df[(asym_df['mark'] == mark) & (asym_df['cluster'] == cluster)]
            vals.append(row['ext_int_ratio'].values[0])
            pvals.append(row['wilcoxon_pval'].values[0])

        offset = (j - 0.5) * bar_width
        bars = ax_bar.bar(x_pos + offset, vals, bar_width, color=color,
                          edgecolor='black', linewidth=0.3, alpha=0.85, label=label)

        for k, (v, p) in enumerate(zip(vals, pvals)):
            if p < 0.05:
                s = '***' if p < 0.001 else '**' if p < 0.01 else '*'
                y_off = 0.004 if v >= 1.0 else -0.008
                va = 'bottom' if v >= 1.0 else 'top'
                ax_bar.text(x_pos[k] + offset, v + y_off, s,
                            ha='center', va=va, fontsize=5, fontweight='bold')

    ax_bar.axhline(1.0, color='#888888', linewidth=0.5, linestyle='--', zorder=0)
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels([BIO_LABELS[c].replace('\n', ' — ') for c in BIO_ORDER], fontsize=5)
    ax_bar.set_ylabel('Ext / Int ratio', fontsize=6)
    ax_bar.set_ylim(0.885, 1.055)
    ax_bar.legend(fontsize=5, loc='lower left', frameon=False,
                  handlelength=1.2, handletextpad=0.3)
    ax_bar.tick_params(axis='both', length=2, pad=1)

    ax_bar.text(-0.6, 1.035, 'exterior-\nenriched', fontsize=4, color='#666666',
                ha='center', va='center', style='italic')
    ax_bar.text(-0.6, 0.94, 'interior-\nenriched', fontsize=4, color='#666666',
                ha='center', va='center', style='italic')

    fig.suptitle('H3K27me3 orientation asymmetry at loop anchors',
                 fontsize=8, fontweight='bold', y=0.98)

    out_path = OUT_DIR / 'orientation_asymmetry_figure.png'
    with multi_format_savefig():
        fig.savefig(str(out_path), dpi=300)
    plt.close(fig)
    print(f'\nSaved: {out_path.stem}.{{png,pdf,svg,jpg}}')


if __name__ == '__main__':
    main()
