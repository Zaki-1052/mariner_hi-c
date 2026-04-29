#!/usr/bin/env python3
# cluster/scripts/10_clust6_subgroup_asymmetry.py
"""Split clust6 (strong loss) by loop length and rerun asymmetry analysis.

Clust6 is heterogeneous: 29% structural + 31% CRE. This script splits at
800kb (1,438 short / 921 long), builds oriented anchor BED6 files for each
sub-group, runs computeMatrix, quantifies exterior/interior K27me3 asymmetry,
and produces a focused 2-panel comparison figure.

Runtime: ~2-4 min (only ~2,800 oriented anchors vs 63k for full pipeline).
"""
import json
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR / 'modules'))
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))

from deepTools_pipeline import bed_pileup          # noqa: E402
from multi_format_output import multi_format_savefig, figure_subfolder  # noqa: E402

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as f:
    plt.rcParams.update(json.load(f))

import os as _os
_out = _os.environ.get('CLUSTER_OUT_DIR', 'outputs/bap1_late')
_k = _os.environ.get('CLUSTER_K', '6')
CLUSTER_FILE = Path(_os.environ.get('CLUSTER_COMBINED',
    str(CLUSTER_DIR / '{}/cluster3/k-{}/data/combined-clusters.txt'.format(_out, _k))))
BIGWIG_BASE  = Path(_os.environ.get('CLUSTER_BIGWIG_DIR', '/Users/zakiralibhai/sdsc/bigwigs'))
BLACKLIST    = REPO_ROOT / 'tads/mm10-blacklist.v2.bed'

ANCHOR_BED_DIR = CLUSTER_DIR / _out / 'figures/deeptools_input'
DEEPTOOLS_DIR  = CLUSTER_DIR / _out / 'figures/deeptools'

SIZE_THRESHOLD = 800_000
OUT_NAME = 'clust6_subgroups'

SUB_ORDER = ['clust6_short', 'clust6_long']
SUB_LABELS = {
    'clust6_short': 'clust6 short\n(<800kb)',
    'clust6_long':  'clust6 long\n(>=800kb)',
}
SUB_DIRECTION = {
    'clust6_short': 'short (<800kb)',
    'clust6_long':  'long (>=800kb)',
}

BIGWIG_DICT = {
    'H3K27ac_ctrl':   BIGWIG_BASE / 'H3K27acCtrl.bw',
    'H3K27ac_mut':    BIGWIG_BASE / 'H3K27acMut.bw',
    'H3K27me3_ctrl':  BIGWIG_BASE / 'H3K27me3Ctrl.bw',
    'H3K27me3_mut':   BIGWIG_BASE / 'H3K27me3Mut.bw',
    'H2AK119ub_ctrl': BIGWIG_BASE / 'H2AK119ubCtrl.bw',
    'H2AK119ub_mut':  BIGWIG_BASE / 'H2AK119ubMut.bw',
    'H3K27me1_ctrl':  BIGWIG_BASE / 'H3K27me1Ctrl.bw',
    'H3K27me1_mut':   BIGWIG_BASE / 'H3K27me1Mut.bw',
}

VMAX_GROUPS = [
    ['H3K27ac_ctrl',   'H3K27ac_mut'],
    ['H3K27me3_ctrl',  'H3K27me3_mut'],
    ['H2AK119ub_ctrl', 'H2AK119ub_mut'],
    ['H3K27me1_ctrl',  'H3K27me1_mut'],
]

COLOR_DICT = {
    'H3K27ac_ctrl':   'Blues',   'H3K27ac_mut':    'Blues',
    'H3K27me3_ctrl':  'Reds',    'H3K27me3_mut':   'Reds',
    'H2AK119ub_ctrl': 'Greens',  'H2AK119ub_mut':  'Greens',
    'H3K27me1_ctrl':  'Purples', 'H3K27me1_mut':   'Purples',
}

XTICKLABELS = ['-5kb (exterior)', 'anchor', '+5kb (interior)']
CTRL_COLOR = '#2166ac'
MUT_COLOR  = '#b2182b'


def build_subgroup_beds(cluster_file, out_dir):
    # type: (Path, Path) -> Dict[str, str]
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(cluster_file, sep='\t')
    clust6 = df[df['GROUP'] == 'clust6'].copy()
    clust6['loop_size'] = clust6['y2'] - clust6['x1']

    short_df = clust6[clust6['loop_size'] < SIZE_THRESHOLD]
    long_df  = clust6[clust6['loop_size'] >= SIZE_THRESHOLD]

    print(f'  clust6 total: {len(clust6)} loops')
    print(f'  clust6_short (<{SIZE_THRESHOLD // 1000}kb): {len(short_df)} loops '
          f'(median {short_df["loop_size"].median() / 1000:.0f}kb)')
    print(f'  clust6_long (>={SIZE_THRESHOLD // 1000}kb): {len(long_df)} loops '
          f'(median {long_df["loop_size"].median() / 1000:.0f}kb)')

    bed_dict = {}
    for name, sub in [('clust6_short', short_df), ('clust6_long', long_df)]:
        a1 = sub[['chr1', 'x1', 'x2']].rename(
            columns={'chr1': 'chrom', 'x1': 'start', 'x2': 'end'})
        a1['name']   = '.'
        a1['score']  = 0
        a1['strand'] = '+'

        a2 = sub[['chr2', 'y1', 'y2']].rename(
            columns={'chr2': 'chrom', 'y1': 'start', 'y2': 'end'})
        a2['name']   = '.'
        a2['score']  = 0
        a2['strand'] = '-'

        oriented = (
            pd.concat([a1, a2])
              .drop_duplicates(subset=['chrom', 'start', 'end', 'strand'])
              .sort_values(['chrom', 'start', 'end', 'strand'])
              .reset_index(drop=True)
        )
        bed_path = out_dir / '{}_oriented_anchors.bed'.format(name)
        oriented[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            bed_path, sep='\t', header=False, index=False, lineterminator='\n')

        n_unique = len(
            pd.concat([a1[['chrom', 'start', 'end']], a2[['chrom', 'start', 'end']]])
              .drop_duplicates())
        n_hub = len(oriented) - n_unique
        print(f'  {name}: {len(sub):>5} loops -> {len(oriented):>5} oriented anchors '
              f'({n_hub} hub) -> {bed_path.name}')
        bed_dict[name] = str(bed_path)
    return bed_dict


def parse_header(filepath):
    # type: (Path) -> Dict
    with open(filepath) as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()

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
    data = pd.read_csv(filepath, sep='\t', header=None, skiprows=3)
    arr = data.select_dtypes(include=[np.number]).values.astype(np.float64)

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


def quantify_asymmetry(values_file, header, out_dir):
    # type: (Path, Dict, Path) -> pd.DataFrame
    n_bins = header['n_bins']
    n_half = n_bins // 2
    bigwigs = header['bigwig_order']
    clusters = header['cluster_names']
    sizes = header['cluster_sizes']

    print(f'  Sub-groups: {clusters}')
    print(f'  Sizes:      {sizes}')
    print(f'  BigWigs:    {len(bigwigs)}')
    print(f'  Bins:       {n_bins} per BigWig, half={n_half}')

    data = pd.read_csv(values_file, sep='\t', header=None, skiprows=3)
    arr = data.select_dtypes(include=[np.number]).values.astype(np.float64)
    print(f'  Matrix:     {arr.shape[0]} rows x {arr.shape[1]} cols')

    cluster_labels = np.concatenate([[name] * size for name, size in zip(clusters, sizes)])
    cluster_labels = cluster_labels[:arr.shape[0]]

    results = []  # type: List[Dict]
    for bw_idx, bw_name in enumerate(bigwigs):
        col_start = bw_idx * n_bins
        ext_cols = slice(col_start, col_start + n_half)
        int_cols = slice(col_start + n_half, col_start + n_bins)

        for cluster in clusters:
            mask = cluster_labels == cluster
            rows = arr[mask]

            ext_per_row = np.nanmean(rows[:, ext_cols], axis=1)
            int_per_row = np.nanmean(rows[:, int_cols], axis=1)

            valid = ~(np.isnan(ext_per_row) | np.isnan(int_per_row))
            ext_v = ext_per_row[valid]
            int_v = int_per_row[valid]

            ext_mean = float(np.mean(ext_v))
            int_mean = float(np.mean(int_v))
            denom = ext_mean + int_mean
            asym = (ext_mean - int_mean) / denom if denom > 0 else 0.0

            diff = ext_v - int_v
            nonzero = diff[diff != 0]
            if len(nonzero) >= 10:
                _, pval = stats.wilcoxon(nonzero, alternative='two-sided')
            else:
                pval = float('nan')

            results.append({
                'mark': bw_name,
                'cluster': cluster,
                'direction': SUB_DIRECTION.get(cluster, ''),
                'n_anchors': int(valid.sum()),
                'exterior_mean': round(ext_mean, 5),
                'interior_mean': round(int_mean, 5),
                'ext_int_ratio': round(ext_mean / int_mean, 4) if int_mean > 0 else float('nan'),
                'asymmetry_index': round(asym, 5),
                'wilcoxon_pval': pval,
            })

    df = pd.DataFrame(results)
    out_tsv = out_dir / 'asymmetry_quantification.tsv'
    df.to_csv(out_tsv, sep='\t', index=False, lineterminator='\n')
    print(f'  Wrote: {out_tsv}')
    return df


def sig(p):
    # type: (float) -> str
    if pd.isna(p): return ''
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'


def print_summary(asym_df):
    # type: (pd.DataFrame) -> None
    def p_str(p):
        if pd.isna(p):  return '       N/A'
        if p < 0.001:   return '{:>10.2e}'.format(p)
        return '{:>10.4f}'.format(p)

    for mark in ['H3K27me3_ctrl', 'H3K27me3_mut']:
        sub = asym_df[asym_df['mark'] == mark]
        print('\n{}'.format('=' * 85))
        print(' {}'.format(mark))
        print('{}'.format('=' * 85))
        print('{:<15} {:<15} {:>6} {:>10} {:>10} {:>8} {:>8} {:>10} {:>3}'.format(
            'Cluster', 'Dir', 'n', 'Exterior', 'Interior', 'Ext/Int', 'Asym', 'p-value', ''))
        print('-' * 85)
        for _, r in sub.iterrows():
            print('{:<15} {:<15} {:>6} {:>10.4f} {:>10.4f} {:>8.4f} {:>+8.5f} {} {}'.format(
                r['cluster'], r['direction'], r['n_anchors'],
                r['exterior_mean'], r['interior_mean'],
                r['ext_int_ratio'], r['asymmetry_index'],
                p_str(r['wilcoxon_pval']), sig(r['wilcoxon_pval'])))

    for mark in ['H3K27ac_ctrl', 'H3K27ac_mut', 'H2AK119ub_ctrl', 'H2AK119ub_mut',
                 'H3K27me1_ctrl', 'H3K27me1_mut']:
        sub = asym_df[asym_df['mark'] == mark]
        print('\n{}'.format('=' * 85))
        print(' {}'.format(mark))
        print('{}'.format('=' * 85))
        print('{:<15} {:<15} {:>6} {:>10} {:>10} {:>8} {:>8} {:>10} {:>3}'.format(
            'Cluster', 'Dir', 'n', 'Exterior', 'Interior', 'Ext/Int', 'Asym', 'p-value', ''))
        print('-' * 85)
        for _, r in sub.iterrows():
            print('{:<15} {:<15} {:>6} {:>10.4f} {:>10.4f} {:>8.4f} {:>+8.5f} {} {}'.format(
                r['cluster'], r['direction'], r['n_anchors'],
                r['exterior_mean'], r['interior_mean'],
                r['ext_int_ratio'], r['asymmetry_index'],
                p_str(r['wilcoxon_pval']), sig(r['wilcoxon_pval'])))


def visualize(values_file, header, asym_df, out_dir):
    # type: (Path, Dict, pd.DataFrame, Path) -> None
    n_bins = header['n_bins']
    x = np.arange(n_bins) * header['bin_size'] - 5000 + header['bin_size'] // 2

    bigwig_indices = {'H3K27me3_ctrl': 2, 'H3K27me3_mut': 3}
    profiles = compute_profiles(values_file, header, bigwig_indices)

    fig = plt.figure(figsize=(3.5, 4.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.3, 1.0],
                  hspace=0.55, wspace=0.15,
                  left=0.13, right=0.95, top=0.90, bottom=0.12)

    y_max = max(profiles[(bw, c)]['mean'].max()
                for bw in bigwig_indices for c in header['cluster_names']) * 1.12

    for i, cluster in enumerate(SUB_ORDER):
        if cluster not in header['cluster_names']:
            continue
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
        ax.set_title(SUB_LABELS[cluster], fontsize=5.5, pad=2, linespacing=1.1)
        ax.set_xticks([-5000, 0, 5000])
        ax.set_xticklabels(['-5kb\next', '0', '+5kb\nint'], fontsize=4.5)
        ax.tick_params(axis='both', length=2, pad=1)

        row_c = asym_df[(asym_df['mark'] == 'H3K27me3_ctrl') & (asym_df['cluster'] == cluster)]
        row_m = asym_df[(asym_df['mark'] == 'H3K27me3_mut') & (asym_df['cluster'] == cluster)]
        ratio_c = row_c['ext_int_ratio'].values[0]
        pval_c = row_c['wilcoxon_pval'].values[0]
        ratio_m = row_m['ext_int_ratio'].values[0]
        pval_m = row_m['wilcoxon_pval'].values[0]
        n_anchors = row_c['n_anchors'].values[0]

        ax.text(0.97, 0.97,
                'n={:,}\nctrl E/I={:.2f} {}\nmut  E/I={:.2f} {}'.format(
                    n_anchors, ratio_c, sig(pval_c), ratio_m, sig(pval_m)),
                transform=ax.transAxes, fontsize=3.8, ha='right', va='top',
                family='monospace', linespacing=1.3,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, linewidth=0.3))

        if i == 0:
            ax.set_ylabel('H3K27me3 (RPKM)', fontsize=6)
        else:
            ax.yaxis.set_ticklabels([])

        if i == 1:
            ax.plot([], [], color=CTRL_COLOR, linewidth=0.8, label='ctrl')
            ax.plot([], [], color=MUT_COLOR, linewidth=0.8, label='mut')
            ax.legend(fontsize=4.5, loc='upper left', frameon=False,
                      handlelength=1.2, handletextpad=0.3)

    ax_bar = fig.add_subplot(gs[1, :])
    bar_width = 0.28
    x_pos = np.arange(len(SUB_ORDER))

    for j, (mark, color, label) in enumerate([
        ('H3K27me3_ctrl', CTRL_COLOR, 'K27me3 ctrl'),
        ('H3K27me3_mut', MUT_COLOR, 'K27me3 mut'),
    ]):
        vals = []
        pvals = []
        for cluster in SUB_ORDER:
            row = asym_df[(asym_df['mark'] == mark) & (asym_df['cluster'] == cluster)]
            vals.append(row['ext_int_ratio'].values[0])
            pvals.append(row['wilcoxon_pval'].values[0])

        offset = (j - 0.5) * bar_width
        ax_bar.bar(x_pos + offset, vals, bar_width, color=color,
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
    ax_bar.set_xticklabels([SUB_LABELS[c].replace('\n', ' ') for c in SUB_ORDER], fontsize=5)
    ax_bar.set_ylabel('Ext / Int ratio', fontsize=6)

    all_ratios = asym_df[asym_df['mark'].str.startswith('H3K27me3')]['ext_int_ratio'].values
    y_min = max(0.88, float(np.min(all_ratios)) * 0.97)
    y_max_bar = min(1.15, float(np.max(all_ratios)) * 1.04)
    ax_bar.set_ylim(y_min, y_max_bar)

    ax_bar.legend(fontsize=5, loc='lower left', frameon=False,
                  handlelength=1.2, handletextpad=0.3)
    ax_bar.tick_params(axis='both', length=2, pad=1)

    ax_bar.text(-0.45, y_max_bar - 0.005, 'exterior-\nenriched', fontsize=4,
                color='#666666', ha='center', va='top', style='italic')
    ax_bar.text(-0.45, y_min + 0.005, 'interior-\nenriched', fontsize=4,
                color='#666666', ha='center', va='bottom', style='italic')

    fig.suptitle('H3K27me3 asymmetry: clust6 short vs long loops',
                 fontsize=8, fontweight='bold', y=0.98)

    out_path = out_dir / 'clust6_subgroup_asymmetry.png'
    with multi_format_savefig():
        fig.savefig(str(out_path), dpi=300)
    plt.close(fig)
    print(f'  Figure: {out_path.stem}.{{png,pdf,svg,jpg}}')


def main():
    # type: () -> None
    for path, label in [(CLUSTER_FILE, 'Cluster file'),
                        (BIGWIG_BASE,  'BigWig dir'),
                        (BLACKLIST,    'Blacklist BED')]:
        if not path.exists():
            raise FileNotFoundError('{} not found: {}'.format(label, path))
    if shutil.which('computeMatrix') is None:
        raise RuntimeError(
            'computeMatrix not on PATH. Prepend '
            '/opt/homebrew/anaconda3/envs/cluster/bin to PATH '
            '(see run_clust6_subgroups.sh).')

    bigwig_dict = {}
    for label, path in BIGWIG_DICT.items():
        if not path.exists():
            raise FileNotFoundError('BigWig missing: {}'.format(path))
        bigwig_dict[label] = str(path)

    print('\n[1/4] Building clust6 sub-group oriented BED6 files...')
    bed_dict = build_subgroup_beds(CLUSTER_FILE, ANCHOR_BED_DIR)

    print('\n[2/4] Running computeMatrix ({} BigWigs x {} sub-groups)...'.format(
        len(bigwig_dict), len(bed_dict)))
    out_dir = figure_subfolder(DEEPTOOLS_DIR, 'clust6_subgroups')
    print('  out_dir: {}'.format(out_dir))

    with multi_format_savefig():
        bed_pileup(
            bed_dict=bed_dict,
            bigWig_dict=bigwig_dict,
            out_dir=str(out_dir),
            out_name=OUT_NAME,
            blacklisted_regions=str(BLACKLIST),
            up_down=5000,
            color_dict=COLOR_DICT,
            vmax_groups=VMAX_GROUPS,
            line_measure='mean',
            pileup_type='referencePoint',
            xticklabels=XTICKLABELS,
        )

    values_file = out_dir / '{}_values'.format(OUT_NAME)
    if not values_file.exists():
        raise FileNotFoundError('computeMatrix output missing: {}'.format(values_file))

    print('\n[3/4] Quantifying exterior/interior asymmetry...')
    header = parse_header(values_file)
    asym_df = quantify_asymmetry(values_file, header, out_dir)

    print('\n[4/4] Generating comparison figure...')
    visualize(values_file, header, asym_df, out_dir)

    print_summary(asym_df)

    print('\n\nClust6 sub-group asymmetry analysis complete.')
    print('Outputs in: {}'.format(out_dir))


if __name__ == '__main__':
    main()
