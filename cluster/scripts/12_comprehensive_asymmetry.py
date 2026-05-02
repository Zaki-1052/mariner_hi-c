#!/usr/bin/env python3
# cluster/scripts/12_comprehensive_asymmetry.py
"""Comprehensive interior/exterior asymmetry for H2AK119ub, H3K27ac, PC1, insulation.

Two-scale design:
  Run A (histone): H2AK119ub_{ctrl,mut} + H3K27ac_{ctrl,mut} at +-5kb, binSize=10bp
  Run B (compartment): PC1_{ctrl,mut} + insulation_{ctrl,mut} at +-50kb, binSize=1000bp

Clusters: clust5, clust6, clust6_short (<800kb), clust6_long (>=800kb).
All 6 clusters included in computeMatrix for completeness; focused figures on
the 4-group contrast.

Steps:
  [0/6] Validate inputs
  [1/6] Compute H3K27ac phasing tracks (pyBigWig per-25kb-bin mean -> TSV)
  [2/6] Compute PC1 BigWigs (cooltools eigs-cis --bigwig)
  [3/6] Compute insulation BigWigs (cooltools insulation --bigwig)
  [4/6] Run computeMatrix x2 (histone +-5kb + compartment +-50kb)
  [5/6] Quantify ext/int asymmetry (Wilcoxon)
  [6/6] Visualize (3 figures, 4 formats each)

Designed for SDSC Expanse via 12_comprehensive_asymmetry.sb.
"""
import argparse
import json
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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

from multi_format_output import multi_format_savefig, figure_subfolder  # noqa: E402

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as _f:
    plt.rcParams.update(json.load(_f))

import os as _os
_out = _os.environ.get('CLUSTER_OUT_DIR', 'outputs/bap1_late')
_k = _os.environ.get('CLUSTER_K', '6')
CLUSTER_FILE = Path(_os.environ.get('CLUSTER_COMBINED',
    str(CLUSTER_DIR / '{}/cluster3/k-{}/data/combined-clusters.txt'.format(_out, _k))))
BLACKLIST    = REPO_ROOT / 'tads/mm10-blacklist.v2.bed'

ANCHOR_BED_DIR = CLUSTER_DIR / _out / 'figures/deeptools_input'
OUT_BASE       = CLUSTER_DIR / _out / 'figures/deeptools/comprehensive_asymmetry'

SIZE_THRESHOLD       = 800_000
MCOOL_RES_PC1        = 25000
MCOOL_RES_INS        = 10000
INS_WINDOW_BP        = 200_000
HISTONE_UP_DOWN      = 5000
COMPARTMENT_UP_DOWN  = 50_000
COMPARTMENT_BIN_SIZE = 1000

FOCUS_ORDER = ['clust5', 'clust6', 'clust6_short', 'clust6_long']
ALL_CLUSTERS = ['clust1', 'clust2', 'clust3', 'clust4', 'clust5', 'clust6',
                'clust6_short', 'clust6_long']

CLUSTER_LABELS = {
    'clust5':       'clust5\ngain (97%)',
    'clust6':       'clust6\nloss (78%)',
    'clust6_short': 'clust6 short\n(<800kb)',
    'clust6_long':  'clust6 long\n(>=800kb)',
    'clust1':       'clust1\nunchanged',
    'clust2':       'clust2\n~unchanged',
    'clust3':       'clust3\nmod loss',
    'clust4':       'clust4\nmod gain',
}

CLUSTER_DIRECTION = {
    'clust5': 'strong gain', 'clust6': 'strong loss',
    'clust6_short': 'loss <800kb', 'clust6_long': 'loss >=800kb',
    'clust1': 'unchanged', 'clust2': '~unchanged',
    'clust3': 'mod loss', 'clust4': 'mod gain',
}

CTRL_COLOR = '#2166ac'
MUT_COLOR  = '#b2182b'


# ---------------------------------------------------------------------------
# Oriented BED construction (reused from scripts 09/10)
# ---------------------------------------------------------------------------

def build_oriented_beds(cluster_file, out_dir):
    # type: (Path, Path) -> Dict[str, str]
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(cluster_file, sep='\t')

    bed_dict = {}  # type: Dict[str, str]
    for cluster in ['clust1', 'clust2', 'clust3', 'clust4', 'clust5', 'clust6']:
        sub = df[df['GROUP'] == cluster]
        if len(sub) == 0:
            continue
        bed_path = out_dir / '{}_oriented_anchors.bed'.format(cluster)
        if bed_path.exists():
            bed_dict[cluster] = str(bed_path)
            continue
        oriented = _make_oriented(sub)
        oriented.to_csv(bed_path, sep='\t', header=False, index=False, lineterminator='\n')
        print('  {} -> {}'.format(cluster, bed_path.name))
        bed_dict[cluster] = str(bed_path)

    clust6 = df[df['GROUP'] == 'clust6'].copy()
    clust6['loop_size'] = clust6['y2'] - clust6['x1']
    for name, sub in [('clust6_short', clust6[clust6['loop_size'] < SIZE_THRESHOLD]),
                      ('clust6_long',  clust6[clust6['loop_size'] >= SIZE_THRESHOLD])]:
        bed_path = out_dir / '{}_oriented_anchors.bed'.format(name)
        if bed_path.exists():
            bed_dict[name] = str(bed_path)
            continue
        oriented = _make_oriented(sub)
        oriented.to_csv(bed_path, sep='\t', header=False, index=False, lineterminator='\n')
        print('  {} -> {} ({} loops -> {} anchors)'.format(
            name, bed_path.name, len(sub), len(oriented)))
        bed_dict[name] = str(bed_path)

    return bed_dict


def _make_oriented(sub):
    # type: (pd.DataFrame) -> pd.DataFrame
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
    return (
        pd.concat([a1, a2])
          .drop_duplicates(subset=['chrom', 'start', 'end', 'strand'])
          .sort_values(['chrom', 'start', 'end', 'strand'])
          .reset_index(drop=True)
    )[['chrom', 'start', 'end', 'name', 'score', 'strand']]


# ---------------------------------------------------------------------------
# PC1 and insulation computation
# ---------------------------------------------------------------------------

def compute_phasing_track(mcool_path, h3k27ac_bw, resolution, out_tsv):
    # type: (str, str, int, Path) -> Path
    import cooler
    import pyBigWig

    print('  Building phasing track from {} at {}kb...'.format(
        Path(h3k27ac_bw).name, resolution // 1000))
    clr = cooler.Cooler('{}::/resolutions/{}'.format(mcool_path, resolution))
    bins = clr.bins()[['chrom', 'start', 'end']][:]

    bw = pyBigWig.open(str(h3k27ac_bw))
    vals = np.full(len(bins), np.nan)
    for chrom in clr.chromnames:
        mask = bins['chrom'] == chrom
        chrom_bins = bins[mask]
        if len(chrom_bins) == 0:
            continue
        try:
            raw = bw.stats(chrom, int(chrom_bins['start'].iloc[0]),
                           int(chrom_bins['end'].iloc[-1]),
                           type='mean', nBins=len(chrom_bins))
        except RuntimeError:
            continue
        arr = np.array([v if v is not None else np.nan for v in raw])
        vals[mask.values] = arr
    bw.close()

    bins['h3k27ac'] = vals
    valid = bins.dropna(subset=['h3k27ac'])
    valid.to_csv(str(out_tsv), sep='\t', header=True, index=False, lineterminator='\n')
    print('    {} bins with signal (of {} total)'.format(len(valid), len(bins)))
    return out_tsv


def compute_pc1_bigwig(mcool_path, phasing_tsv, resolution, out_prefix):
    # type: (str, Path, int, Path) -> Path
    out_bw = Path('{}.cis.bw'.format(out_prefix))
    print('  Running cooltools eigs-cis -> {}'.format(out_bw.name))
    cmd = [
        'cooltools', 'eigs-cis',
        '{}::/resolutions/{}'.format(mcool_path, resolution),
        '--phasing-track', '{}::h3k27ac'.format(phasing_tsv),
        '--n-eigs', '1',
        '--clr-weight-name', 'weight',
        '--out-prefix', str(out_prefix),
        '--bigwig',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print('  STDERR: {}'.format(result.stderr[-500:] if result.stderr else '(none)'))
        raise RuntimeError('cooltools eigs-cis failed (exit={})'.format(result.returncode))
    if not out_bw.exists():
        raise FileNotFoundError('Expected BigWig not produced: {}'.format(out_bw))
    print('    OK: {}'.format(out_bw))
    return out_bw


def compute_insulation_bigwig(mcool_path, resolution, window_bp, out_prefix, nproc):
    # type: (str, int, int, Path, int) -> Path
    out_bw = Path('{}.{}.bw'.format(out_prefix, window_bp))
    print('  Running cooltools insulation -> {}'.format(out_bw.name))
    cmd = [
        'cooltools', 'insulation',
        '{}::/resolutions/{}'.format(mcool_path, resolution),
        str(window_bp),
        '-p', str(nproc),
        '-o', str(out_prefix),
        '--bigwig',
        '--clr-weight-name', 'weight',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print('  STDERR: {}'.format(result.stderr[-500:] if result.stderr else '(none)'))
        raise RuntimeError('cooltools insulation failed (exit={})'.format(result.returncode))
    if not out_bw.exists():
        raise FileNotFoundError('Expected BigWig not produced: {}'.format(out_bw))
    print('    OK: {}'.format(out_bw))
    return out_bw


# ---------------------------------------------------------------------------
# computeMatrix wrapper
# ---------------------------------------------------------------------------

def run_computematrix(bw_dict, bed_dict, out_dir, name, up_down,
                      blacklist, nproc, bin_size=None, missing_as_zero=True):
    # type: (Dict[str, str], Dict[str, str], Path, str, int, str, int, Optional[int], bool) -> Path
    out_dir.mkdir(parents=True, exist_ok=True)
    out_gz     = out_dir / name
    out_values = out_dir / '{}_values'.format(name)

    bw_files  = ' '.join(str(v) for v in bw_dict.values())
    bw_labels = ' '.join(bw_dict.keys())
    bed_files = ' '.join(str(v) for v in bed_dict.values())

    cmd = ('computeMatrix reference-point --referencePoint center'
           ' -S {bw_files}'
           ' -R {bed_files}'
           ' -o {out_gz}'
           ' -bl {bl}'
           ' -p {nproc}'
           ' --outFileNameMatrix {out_values}'
           ' -b {up} -a {up}'
           ' --sortRegions keep'
           ' --samplesLabel {labels}').format(
        bw_files=bw_files, bed_files=bed_files,
        out_gz=out_gz, bl=blacklist, nproc=nproc,
        out_values=out_values, up=up_down, labels=bw_labels)

    if missing_as_zero:
        cmd += ' --missingDataAsZero'
    if bin_size is not None:
        cmd += ' --binSize {}'.format(bin_size)

    print('  CMD: {}'.format(cmd[:200]))
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print('  STDERR: {}'.format(result.stderr[-500:] if result.stderr else '(none)'))
        raise RuntimeError('computeMatrix failed (exit={})'.format(result.returncode))
    if not out_values.exists():
        raise FileNotFoundError('computeMatrix output missing: {}'.format(out_values))
    print('    OK: {}'.format(out_values))
    return out_values


# ---------------------------------------------------------------------------
# Parsing and quantification (reused from scripts 09/10)
# ---------------------------------------------------------------------------

def parse_header(filepath):
    # type: (Path,) -> Dict
    with open(filepath) as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()

    tokens1 = [t.replace('#', '') for t in line1.split('\t') if ':' in t]
    cluster_names = [t.split(':')[0].replace('_oriented_anchors.bed', '') for t in tokens1]
    cluster_sizes = [int(t.split(':')[1]) for t in tokens1]

    params = {}  # type: Dict[str, int]
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


def quantify_asymmetry(values_file, header, direction_map, out_tsv):
    # type: (Path, Dict, Dict[str, str], Path) -> pd.DataFrame
    n_bins = header['n_bins']
    n_half = n_bins // 2
    bigwigs = header['bigwig_order']
    clusters = header['cluster_names']
    sizes = header['cluster_sizes']

    print('  Sub-groups: {}'.format(clusters))
    print('  Sizes:      {}'.format(sizes))
    print('  BigWigs:    {}'.format(len(bigwigs)))
    print('  Bins:       {} per BigWig, half={}'.format(n_bins, n_half))

    data = pd.read_csv(values_file, sep='\t', header=None, skiprows=3)
    arr = data.select_dtypes(include=[np.number]).values.astype(np.float64)
    print('  Matrix:     {} rows x {} cols'.format(arr.shape[0], arr.shape[1]))

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

            ext_mean = float(np.mean(ext_v)) if len(ext_v) > 0 else float('nan')
            int_mean = float(np.mean(int_v)) if len(int_v) > 0 else float('nan')
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
                'direction': direction_map.get(cluster, ''),
                'n_anchors': int(valid.sum()),
                'exterior_mean': round(ext_mean, 5),
                'interior_mean': round(int_mean, 5),
                'ext_int_ratio': round(ext_mean / int_mean, 4) if int_mean != 0 else float('nan'),
                'asymmetry_index': round(asym, 5),
                'wilcoxon_pval': pval,
            })

    df = pd.DataFrame(results)
    df.to_csv(str(out_tsv), sep='\t', index=False, lineterminator='\n')
    print('  Wrote: {}'.format(out_tsv))
    return df


def compute_profiles(filepath, header, bigwig_indices):
    # type: (Path, Dict, Dict[str, int]) -> Dict[Tuple[str, str], Dict]
    n_bins = header['n_bins']
    data = pd.read_csv(filepath, sep='\t', header=None, skiprows=3)
    arr = data.select_dtypes(include=[np.number]).values.astype(np.float64)

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
            sem = np.nanstd(rows, axis=0) / np.sqrt(n) if n > 1 else np.zeros_like(mean)
            result[(bw_name, cluster)] = {'mean': mean, 'sem': sem, 'n': n}
    return result


# ---------------------------------------------------------------------------
# Visualization helpers
# ---------------------------------------------------------------------------

def _sig(p):
    # type: (float,) -> str
    if pd.isna(p):  return ''
    if p < 0.001:   return '***'
    if p < 0.01:    return '**'
    if p < 0.05:    return '*'
    return 'ns'


def _profile_panel(ax, x, profiles, cluster, bw_ctrl, bw_mut, asym_df):
    # type: (...) -> None
    ax.axvspan(x[0], 0, alpha=0.06, color='#4393c3', zorder=0)
    ax.axvspan(0, x[-1], alpha=0.06, color='#f4a582', zorder=0)
    ax.axvline(0, color='#888888', linewidth=0.3, linestyle='--', zorder=1)

    for bw_name, color in [(bw_ctrl, CTRL_COLOR), (bw_mut, MUT_COLOR)]:
        key = (bw_name, cluster)
        if key not in profiles:
            continue
        p = profiles[key]
        ax.plot(x, p['mean'], color=color, linewidth=0.6, zorder=2)
        ax.fill_between(x, p['mean'] - p['sem'], p['mean'] + p['sem'],
                        color=color, alpha=0.12, zorder=2)

    row_c = asym_df[(asym_df['mark'] == bw_ctrl) & (asym_df['cluster'] == cluster)]
    row_m = asym_df[(asym_df['mark'] == bw_mut) & (asym_df['cluster'] == cluster)]
    if len(row_c) > 0 and len(row_m) > 0:
        n = row_c['n_anchors'].values[0]
        rc = row_c['ext_int_ratio'].values[0]
        pc = row_c['wilcoxon_pval'].values[0]
        rm = row_m['ext_int_ratio'].values[0]
        pm = row_m['wilcoxon_pval'].values[0]
        ax.text(0.97, 0.97,
                'n={:,}\nctrl E/I={:.2f} {}\nmut  E/I={:.2f} {}'.format(
                    n, rc, _sig(pc), rm, _sig(pm)),
                transform=ax.transAxes, fontsize=3.5, ha='right', va='top',
                family='monospace', linespacing=1.3,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          alpha=0.7, linewidth=0.3))


def _bar_panel(ax, asym_df, marks_colors, clusters, y_label):
    # type: (...) -> None
    bar_width = 0.28
    x_pos = np.arange(len(clusters))

    for j, (mark, color, label) in enumerate(marks_colors):
        vals = []
        pvals = []
        for c in clusters:
            row = asym_df[(asym_df['mark'] == mark) & (asym_df['cluster'] == c)]
            if len(row) > 0:
                vals.append(row['ext_int_ratio'].values[0])
                pvals.append(row['wilcoxon_pval'].values[0])
            else:
                vals.append(float('nan'))
                pvals.append(float('nan'))

        offset = (j - 0.5) * bar_width
        ax.bar(x_pos + offset, vals, bar_width, color=color,
               edgecolor='black', linewidth=0.3, alpha=0.85, label=label)

        for k, (v, p) in enumerate(zip(vals, pvals)):
            if not pd.isna(p) and p < 0.05:
                s = '***' if p < 0.001 else '**' if p < 0.01 else '*'
                y_off = 0.004 if v >= 1.0 else -0.008
                va = 'bottom' if v >= 1.0 else 'top'
                ax.text(x_pos[k] + offset, v + y_off, s,
                        ha='center', va=va, fontsize=5, fontweight='bold')

    ax.axhline(1.0, color='#888888', linewidth=0.5, linestyle='--', zorder=0)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([CLUSTER_LABELS.get(c, c).replace('\n', ' ') for c in clusters],
                       fontsize=5)
    ax.set_ylabel(y_label, fontsize=6)
    ax.legend(fontsize=4.5, loc='best', frameon=False,
              handlelength=1.2, handletextpad=0.3)
    ax.tick_params(axis='both', length=2, pad=1)


# ---------------------------------------------------------------------------
# Figure functions
# ---------------------------------------------------------------------------

def visualize_histone(values_file, header, asym_df, out_dir):
    # type: (Path, Dict, pd.DataFrame, Path) -> None
    n_bins = header['n_bins']
    x = np.arange(n_bins) * header['bin_size'] - HISTONE_UP_DOWN + header['bin_size'] // 2

    bw_order = header['bigwig_order']
    bw_idx = {name: i for i, name in enumerate(bw_order)}
    profiles = compute_profiles(values_file, header, bw_idx)

    clusters = [c for c in FOCUS_ORDER if c in header['cluster_names']]
    n_cols = len(clusters)

    fig = plt.figure(figsize=(2.0 * n_cols, 5.5))
    gs = GridSpec(3, n_cols, figure=fig, height_ratios=[1.3, 1.3, 1.0],
                  hspace=0.55, wspace=0.15,
                  left=0.12, right=0.95, top=0.90, bottom=0.10)

    mark_rows = [
        ('H2AK119ub_ctrl', 'H2AK119ub_mut', 'H2AK119ub (RPKM)'),
        ('H3K27ac_ctrl',   'H3K27ac_mut',   'H3K27ac (RPKM)'),
    ]

    for row_i, (bw_ctrl, bw_mut, ylabel) in enumerate(mark_rows):
        y_max = 0
        for c in clusters:
            for bw in [bw_ctrl, bw_mut]:
                key = (bw, c)
                if key in profiles:
                    y_max = max(y_max, np.nanmax(profiles[key]['mean']))
        y_max *= 1.15

        for col_i, cluster in enumerate(clusters):
            ax = fig.add_subplot(gs[row_i, col_i])
            _profile_panel(ax, x, profiles, cluster, bw_ctrl, bw_mut, asym_df)
            ax.set_xlim(-HISTONE_UP_DOWN, HISTONE_UP_DOWN)
            ax.set_ylim(0, y_max)
            if row_i == 0:
                ax.set_title(CLUSTER_LABELS.get(cluster, cluster),
                             fontsize=5.5, pad=2, linespacing=1.1)
            ax.set_xticks([-HISTONE_UP_DOWN, 0, HISTONE_UP_DOWN])
            ax.set_xticklabels(['-5kb\next', '0', '+5kb\nint'], fontsize=4.5)
            ax.tick_params(axis='both', length=2, pad=1)
            if col_i == 0:
                ax.set_ylabel(ylabel, fontsize=6)
            else:
                ax.yaxis.set_ticklabels([])
            if col_i == n_cols - 1 and row_i == 0:
                ax.plot([], [], color=CTRL_COLOR, linewidth=0.8, label='ctrl')
                ax.plot([], [], color=MUT_COLOR, linewidth=0.8, label='mut')
                ax.legend(fontsize=4.5, loc='upper left', frameon=False,
                          handlelength=1.2, handletextpad=0.3)

    ax_bar = fig.add_subplot(gs[2, :])
    _bar_panel(ax_bar, asym_df, [
        ('H2AK119ub_ctrl', CTRL_COLOR, 'K119ub ctrl'),
        ('H2AK119ub_mut',  MUT_COLOR,  'K119ub mut'),
    ], clusters, 'Ext / Int ratio')

    all_ratios = asym_df[asym_df['mark'].str.contains('H2AK119ub')]['ext_int_ratio'].dropna()
    if len(all_ratios) > 0:
        y_lo = max(0.88, float(all_ratios.min()) * 0.97)
        y_hi = min(1.15, float(all_ratios.max()) * 1.04)
        ax_bar.set_ylim(y_lo, y_hi)

    fig.suptitle('H2AK119ub and H3K27ac orientation asymmetry at loop anchors',
                 fontsize=8, fontweight='bold', y=0.98)

    out_path = out_dir / 'histone_asymmetry.png'
    with multi_format_savefig():
        fig.savefig(str(out_path), dpi=300)
    plt.close(fig)
    print('  Figure: {}'.format(out_path.stem))


def visualize_compartment(values_file, header, asym_df, out_dir):
    # type: (Path, Dict, pd.DataFrame, Path) -> None
    n_bins = header['n_bins']
    x = np.arange(n_bins) * header['bin_size'] - COMPARTMENT_UP_DOWN + header['bin_size'] // 2

    bw_order = header['bigwig_order']
    bw_idx = {name: i for i, name in enumerate(bw_order)}
    profiles = compute_profiles(values_file, header, bw_idx)

    clusters = [c for c in FOCUS_ORDER if c in header['cluster_names']]
    n_cols = len(clusters)

    fig = plt.figure(figsize=(2.0 * n_cols, 5.5))
    gs = GridSpec(3, n_cols, figure=fig, height_ratios=[1.3, 1.3, 1.0],
                  hspace=0.55, wspace=0.15,
                  left=0.12, right=0.95, top=0.90, bottom=0.10)

    mark_rows = [
        ('PC1_ctrl', 'PC1_mut', 'PC1 eigenvector'),
        ('insulation_ctrl', 'insulation_mut', 'log2 insulation\n(200kb window)'),
    ]

    for row_i, (bw_ctrl, bw_mut, ylabel) in enumerate(mark_rows):
        y_vals = []
        for c in clusters:
            for bw in [bw_ctrl, bw_mut]:
                key = (bw, c)
                if key in profiles:
                    m = profiles[key]['mean']
                    y_vals.extend([np.nanmin(m), np.nanmax(m)])
        if y_vals:
            y_lo = min(y_vals) * 1.1 if min(y_vals) < 0 else min(y_vals) * 0.9
            y_hi = max(y_vals) * 1.1
        else:
            y_lo, y_hi = -1, 1

        for col_i, cluster in enumerate(clusters):
            ax = fig.add_subplot(gs[row_i, col_i])
            _profile_panel(ax, x, profiles, cluster, bw_ctrl, bw_mut, asym_df)
            ax.set_xlim(-COMPARTMENT_UP_DOWN, COMPARTMENT_UP_DOWN)
            ax.set_ylim(y_lo, y_hi)
            if row_i == 0:
                ax.axhline(0, color='#888888', linewidth=0.3, linestyle=':', zorder=0)
                ax.set_title(CLUSTER_LABELS.get(cluster, cluster),
                             fontsize=5.5, pad=2, linespacing=1.1)
            ax.set_xticks([-COMPARTMENT_UP_DOWN, 0, COMPARTMENT_UP_DOWN])
            ax.set_xticklabels(['-50kb\next', '0', '+50kb\nint'], fontsize=4.5)
            ax.tick_params(axis='both', length=2, pad=1)
            if col_i == 0:
                ax.set_ylabel(ylabel, fontsize=6)
            else:
                ax.yaxis.set_ticklabels([])
            if col_i == n_cols - 1 and row_i == 0:
                ax.plot([], [], color=CTRL_COLOR, linewidth=0.8, label='ctrl')
                ax.plot([], [], color=MUT_COLOR, linewidth=0.8, label='mut')
                ax.legend(fontsize=4.5, loc='upper left', frameon=False,
                          handlelength=1.2, handletextpad=0.3)

    ax_bar = fig.add_subplot(gs[2, :])
    _bar_panel(ax_bar, asym_df, [
        ('PC1_ctrl', CTRL_COLOR, 'PC1 ctrl'),
        ('PC1_mut',  MUT_COLOR,  'PC1 mut'),
    ], clusters, 'Ext / Int ratio')

    fig.suptitle('PC1 and insulation score orientation asymmetry at loop anchors',
                 fontsize=8, fontweight='bold', y=0.98)

    out_path = out_dir / 'compartment_asymmetry.png'
    with multi_format_savefig():
        fig.savefig(str(out_path), dpi=300)
    plt.close(fig)
    print('  Figure: {}'.format(out_path.stem))


def visualize_summary(histone_asym, compartment_asym, out_dir):
    # type: (pd.DataFrame, pd.DataFrame, Path) -> None
    combined = pd.concat([histone_asym, compartment_asym], ignore_index=True)
    focus = ['clust5', 'clust6']

    mark_groups = [
        ('H2AK119ub_ctrl', 'H2AK119ub_mut', 'H2AK119ub'),
        ('H3K27ac_ctrl',   'H3K27ac_mut',   'H3K27ac'),
        ('PC1_ctrl',       'PC1_mut',        'PC1'),
        ('insulation_ctrl','insulation_mut',  'Insulation'),
    ]

    fig = plt.figure(figsize=(7.0, 3.5))
    gs = GridSpec(1, len(mark_groups), figure=fig, wspace=0.35,
                  left=0.08, right=0.95, top=0.85, bottom=0.18)

    for col_i, (ctrl, mut, title) in enumerate(mark_groups):
        ax = fig.add_subplot(gs[0, col_i])
        bar_width = 0.22
        x_pos = np.arange(len(focus))

        for j, (mark, color, label) in enumerate([
            (ctrl, CTRL_COLOR, 'ctrl'), (mut, MUT_COLOR, 'mut'),
        ]):
            vals = []
            pvals = []
            for c in focus:
                row = combined[(combined['mark'] == mark) & (combined['cluster'] == c)]
                if len(row) > 0:
                    vals.append(row['ext_int_ratio'].values[0])
                    pvals.append(row['wilcoxon_pval'].values[0])
                else:
                    vals.append(float('nan'))
                    pvals.append(float('nan'))

            offset = (j - 0.5) * bar_width
            ax.bar(x_pos + offset, vals, bar_width, color=color,
                   edgecolor='black', linewidth=0.3, alpha=0.85, label=label)

            for k, (v, p) in enumerate(zip(vals, pvals)):
                if not pd.isna(p) and p < 0.05:
                    s = '***' if p < 0.001 else '**' if p < 0.01 else '*'
                    y_off = 0.004 if v >= 1.0 else -0.008
                    va = 'bottom' if v >= 1.0 else 'top'
                    ax.text(x_pos[k] + offset, v + y_off, s,
                            ha='center', va=va, fontsize=6, fontweight='bold')

        ax.axhline(1.0, color='#888888', linewidth=0.5, linestyle='--', zorder=0)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(['clust5\n(gain)', 'clust6\n(loss)'], fontsize=6)
        ax.set_title(title, fontsize=7, fontweight='bold')
        ax.tick_params(axis='both', length=2, pad=1)
        if col_i == 0:
            ax.set_ylabel('Ext / Int ratio', fontsize=6)
            ax.legend(fontsize=5, loc='best', frameon=False,
                      handlelength=1.0, handletextpad=0.3)

    fig.suptitle('Asymmetry summary: gained vs lost loops',
                 fontsize=9, fontweight='bold', y=0.97)

    out_path = out_dir / 'asymmetry_summary.png'
    with multi_format_savefig():
        fig.savefig(str(out_path), dpi=300)
    plt.close(fig)
    print('  Figure: {}'.format(out_path.stem))


# ---------------------------------------------------------------------------
# Summary print
# ---------------------------------------------------------------------------

def print_summary(histone_df, compartment_df):
    # type: (pd.DataFrame, pd.DataFrame) -> None
    combined = pd.concat([histone_df, compartment_df], ignore_index=True)

    def p_str(p):
        if pd.isna(p):  return '       N/A'
        if p < 0.001:   return '{:>10.2e}'.format(p)
        return '{:>10.4f}'.format(p)

    for mark in combined['mark'].unique():
        sub = combined[combined['mark'] == mark]
        print('\n{}'.format('=' * 90))
        print(' {}'.format(mark))
        print('{}'.format('=' * 90))
        print('{:<15} {:<15} {:>6} {:>10} {:>10} {:>8} {:>8} {:>10} {:>3}'.format(
            'Cluster', 'Dir', 'n', 'Exterior', 'Interior', 'Ext/Int', 'Asym', 'p-value', ''))
        print('-' * 90)
        for _, r in sub.iterrows():
            print('{:<15} {:<15} {:>6} {:>10.4f} {:>10.4f} {:>8.4f} {:>+8.5f} {} {}'.format(
                r['cluster'], r['direction'], r['n_anchors'],
                r['exterior_mean'], r['interior_mean'],
                r['ext_int_ratio'], r['asymmetry_index'],
                p_str(r['wilcoxon_pval']), _sig(r['wilcoxon_pval'])))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # type: () -> None
    ap = argparse.ArgumentParser(
        description='Comprehensive interior/exterior asymmetry analysis')
    ap.add_argument('--mcool-ctrl', required=True, help='Path to ctrl mcool')
    ap.add_argument('--mcool-mut',  required=True, help='Path to mut mcool')
    ap.add_argument('--bigwig-dir', default='/Users/zakiralibhai/sdsc/bigwigs',
                    help='Directory with H2AK119ub/H3K27ac BigWigs')
    ap.add_argument('--nproc', type=int, default=4,
                    help='Number of processes for cooltools/computeMatrix')
    ap.add_argument('--skip-compute', action='store_true',
                    help='Skip PC1/insulation BigWig generation (use existing)')
    args = ap.parse_args()

    # -- [0/6] Validate inputs ------------------------------------------------
    print('\n[0/6] Validating inputs...')
    for path, label in [(CLUSTER_FILE, 'Cluster file'),
                        (BLACKLIST,    'Blacklist BED')]:
        if not path.exists():
            raise FileNotFoundError('{} not found: {}'.format(label, path))
    for label, path_str in [('ctrl mcool', args.mcool_ctrl),
                            ('mut mcool',  args.mcool_mut)]:
        if not Path(path_str).exists():
            raise FileNotFoundError('{} not found: {}'.format(label, path_str))

    bw_dir = Path(args.bigwig_dir)
    histone_bws = OrderedDict([
        ('H2AK119ub_ctrl', bw_dir / 'H2AK119ubCtrl.bw'),
        ('H2AK119ub_mut',  bw_dir / 'H2AK119ubMut.bw'),
        ('H3K27ac_ctrl',   bw_dir / 'H3K27acCtrl.bw'),
        ('H3K27ac_mut',    bw_dir / 'H3K27acMut.bw'),
    ])
    for label, path in histone_bws.items():
        if not path.exists():
            raise FileNotFoundError('BigWig not found: {} -> {}'.format(label, path))
    print('  Histone BigWigs: OK ({})'.format(len(histone_bws)))

    if shutil.which('computeMatrix') is None:
        raise RuntimeError(
            'computeMatrix not on PATH. Activate the cluster conda env.')
    if not args.skip_compute and shutil.which('bedGraphToBigWig') is None:
        raise RuntimeError(
            'bedGraphToBigWig not on PATH. Install: '
            'conda install -c bioconda ucsc-bedgraphtobigwig')
    print('  Tools on PATH: OK')

    # -- Build oriented BED files if needed -----------------------------------
    print('\n  Checking oriented anchor BEDs...')
    bed_dict = build_oriented_beds(CLUSTER_FILE, ANCHOR_BED_DIR)
    focus_beds = OrderedDict()
    for c in FOCUS_ORDER:
        if c in bed_dict:
            focus_beds[c] = bed_dict[c]
        else:
            bed_path = ANCHOR_BED_DIR / '{}_oriented_anchors.bed'.format(c)
            if bed_path.exists():
                focus_beds[c] = str(bed_path)
            else:
                raise FileNotFoundError(
                    'Oriented BED not found for {}: {}'.format(c, bed_path))
    print('  BEDs: {}'.format(list(focus_beds.keys())))

    out_dir = OUT_BASE
    out_dir.mkdir(parents=True, exist_ok=True)
    bw_out = out_dir / 'bigwigs'
    bw_out.mkdir(parents=True, exist_ok=True)

    # -- [1/6] Phasing tracks ------------------------------------------------
    if not args.skip_compute:
        print('\n[1/6] Computing H3K27ac phasing tracks...')
        ctrl_phasing = compute_phasing_track(
            args.mcool_ctrl, str(histone_bws['H3K27ac_ctrl']),
            MCOOL_RES_PC1, bw_out / 'ctrl_phasing_25kb.tsv')
        mut_phasing = compute_phasing_track(
            args.mcool_mut, str(histone_bws['H3K27ac_mut']),
            MCOOL_RES_PC1, bw_out / 'mut_phasing_25kb.tsv')
    else:
        ctrl_phasing = bw_out / 'ctrl_phasing_25kb.tsv'
        mut_phasing  = bw_out / 'mut_phasing_25kb.tsv'
        print('\n[1/6] Skipped (--skip-compute)')

    # -- [2/6] PC1 BigWigs ---------------------------------------------------
    pc1_ctrl_bw = Path('{}.cis.bw'.format(bw_out / 'PC1_ctrl'))
    pc1_mut_bw  = Path('{}.cis.bw'.format(bw_out / 'PC1_mut'))

    if not args.skip_compute:
        print('\n[2/6] Computing PC1 BigWigs (cooltools eigs-cis)...')
        pc1_ctrl_bw = compute_pc1_bigwig(
            args.mcool_ctrl, ctrl_phasing, MCOOL_RES_PC1, bw_out / 'PC1_ctrl')
        pc1_mut_bw = compute_pc1_bigwig(
            args.mcool_mut, mut_phasing, MCOOL_RES_PC1, bw_out / 'PC1_mut')

        import cooler
        ctrl_vecs = pd.read_csv(str(bw_out / 'PC1_ctrl.cis.vecs.tsv'), sep='\t')
        mut_vecs  = pd.read_csv(str(bw_out / 'PC1_mut.cis.vecs.tsv'), sep='\t')
        merged = ctrl_vecs.merge(mut_vecs, on=['chrom', 'start', 'end'],
                                 suffixes=('_ctrl', '_mut'))
        e1_cols = [c for c in merged.columns if c.startswith('E1')]
        if len(e1_cols) >= 2:
            valid = merged.dropna(subset=e1_cols[:2])
            corr = valid[e1_cols[0]].corr(valid[e1_cols[1]])
            print('  PC1 ctrl-mut correlation: r={:.4f} (n={})'.format(corr, len(valid)))
    else:
        print('\n[2/6] Skipped (--skip-compute)')

    # -- [3/6] Insulation BigWigs ---------------------------------------------
    ins_ctrl_bw = Path('{}.{}.bw'.format(bw_out / 'insulation_ctrl_10kb', INS_WINDOW_BP))
    ins_mut_bw  = Path('{}.{}.bw'.format(bw_out / 'insulation_mut_10kb', INS_WINDOW_BP))

    if not args.skip_compute:
        print('\n[3/6] Computing insulation BigWigs (cooltools insulation)...')
        ins_ctrl_bw = compute_insulation_bigwig(
            args.mcool_ctrl, MCOOL_RES_INS, INS_WINDOW_BP,
            bw_out / 'insulation_ctrl_10kb', args.nproc)
        ins_mut_bw = compute_insulation_bigwig(
            args.mcool_mut, MCOOL_RES_INS, INS_WINDOW_BP,
            bw_out / 'insulation_mut_10kb', args.nproc)
    else:
        print('\n[3/6] Skipped (--skip-compute)')

    for bw, label in [(pc1_ctrl_bw, 'PC1 ctrl'), (pc1_mut_bw, 'PC1 mut'),
                      (ins_ctrl_bw, 'insulation ctrl'), (ins_mut_bw, 'insulation mut')]:
        if not bw.exists():
            raise FileNotFoundError('{} BigWig not found: {}. '
                                    'Run without --skip-compute first.'.format(label, bw))

    # -- [4/6] computeMatrix x2 ----------------------------------------------
    print('\n[4/6] Running computeMatrix (histone +-5kb + compartment +-50kb)...')

    histone_bw_paths = OrderedDict([
        ('H2AK119ub_ctrl', str(histone_bws['H2AK119ub_ctrl'])),
        ('H2AK119ub_mut',  str(histone_bws['H2AK119ub_mut'])),
        ('H3K27ac_ctrl',   str(histone_bws['H3K27ac_ctrl'])),
        ('H3K27ac_mut',    str(histone_bws['H3K27ac_mut'])),
    ])

    compartment_bw_paths = OrderedDict([
        ('PC1_ctrl',       str(pc1_ctrl_bw)),
        ('PC1_mut',        str(pc1_mut_bw)),
        ('insulation_ctrl', str(ins_ctrl_bw)),
        ('insulation_mut',  str(ins_mut_bw)),
    ])

    focus_bed_paths = OrderedDict([(k, str(v)) for k, v in focus_beds.items()])

    print('\n  Histone run (+-5kb, {} BigWigs x {} groups)...'.format(
        len(histone_bw_paths), len(focus_bed_paths)))
    histone_values = run_computematrix(
        histone_bw_paths, focus_bed_paths, out_dir, 'histone_anchors',
        HISTONE_UP_DOWN, str(BLACKLIST), args.nproc,
        missing_as_zero=True)

    print('\n  Compartment run (+-50kb/1kb bins, {} BigWigs x {} groups)...'.format(
        len(compartment_bw_paths), len(focus_bed_paths)))
    compartment_values = run_computematrix(
        compartment_bw_paths, focus_bed_paths, out_dir, 'compartment_anchors',
        COMPARTMENT_UP_DOWN, str(BLACKLIST), args.nproc,
        bin_size=COMPARTMENT_BIN_SIZE, missing_as_zero=False)

    # -- [5/6] Quantify asymmetry ---------------------------------------------
    print('\n[5/6] Quantifying exterior/interior asymmetry...')

    print('\n  Histone asymmetry:')
    histone_header = parse_header(histone_values)
    histone_asym = quantify_asymmetry(
        histone_values, histone_header, CLUSTER_DIRECTION,
        out_dir / 'asymmetry_histone.tsv')

    print('\n  Compartment asymmetry:')
    compartment_header = parse_header(compartment_values)
    compartment_asym = quantify_asymmetry(
        compartment_values, compartment_header, CLUSTER_DIRECTION,
        out_dir / 'asymmetry_compartment.tsv')

    # -- [6/6] Visualize ------------------------------------------------------
    print('\n[6/6] Generating figures...')

    print('\n  Histone figure:')
    visualize_histone(histone_values, histone_header, histone_asym, out_dir)

    print('\n  Compartment figure:')
    visualize_compartment(compartment_values, compartment_header, compartment_asym, out_dir)

    print('\n  Summary figure:')
    visualize_summary(histone_asym, compartment_asym, out_dir)

    # -- Summary tables -------------------------------------------------------
    print_summary(histone_asym, compartment_asym)

    print('\n\nComprehensive asymmetry analysis complete.')
    print('Outputs in: {}'.format(out_dir))


if __name__ == '__main__':
    main()
