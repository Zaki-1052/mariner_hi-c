#!/usr/bin/env python3
# cluster/scripts/quantify_orientation_asymmetry.py
"""Quantify exterior vs interior signal asymmetry from oriented anchor metagene.

Reads the oriented_anchors_values matrix, splits each BigWig's bins into
exterior half (-5kb to 0) and interior half (0 to +5kb), computes per-anchor
mean signal in each half, then reports asymmetry index and Wilcoxon test.
"""
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
VALUES_FILE = CLUSTER_DIR / 'bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values'
OUT_DIR     = CLUSTER_DIR / 'bap1_late/figures/deeptools/oriented_anchors'

CLUSTER_DIRECTION = {
    'clust1': 'unchanged',
    'clust2': '~unchanged',
    'clust3': 'mod loss',
    'clust4': 'mod gain',
    'clust5': 'strong gain',
    'clust6': 'strong loss',
}


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


def main():
    # type: () -> None
    if not VALUES_FILE.exists():
        print(f'ERROR: {VALUES_FILE} not found', file=sys.stderr)
        sys.exit(1)

    header = parse_header(VALUES_FILE)
    n_bins = header['n_bins']
    n_half = n_bins // 2
    bigwigs = header['bigwig_order']
    clusters = header['cluster_names']
    sizes = header['cluster_sizes']

    print(f'Clusters: {clusters}')
    print(f'Sizes:    {sizes} (total {sum(sizes)})')
    print(f'BigWigs:  {bigwigs}')
    print(f'Bins:     {n_bins} per BigWig, half={n_half}')
    print(f'\nLoading {sum(sizes)} x {len(bigwigs)*n_bins} matrix...')

    data = pd.read_csv(VALUES_FILE, sep='\t', header=None, skiprows=3)
    numeric = data.select_dtypes(include=[np.number])
    arr = numeric.values.astype(np.float64)
    print(f'Loaded: {arr.shape[0]} rows x {arr.shape[1]} cols')

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
                'direction': CLUSTER_DIRECTION.get(cluster, ''),
                'n_anchors': int(valid.sum()),
                'exterior_mean': round(ext_mean, 5),
                'interior_mean': round(int_mean, 5),
                'ext_int_ratio': round(ext_mean / int_mean, 4) if int_mean > 0 else float('nan'),
                'asymmetry_index': round(asym, 5),
                'wilcoxon_pval': pval,
            })

    df = pd.DataFrame(results)
    out_tsv = OUT_DIR / 'asymmetry_quantification.tsv'
    df.to_csv(out_tsv, sep='\t', index=False, lineterminator='\n')
    print(f'\nFull table: {out_tsv}')

    def sig_str(p):
        if pd.isna(p): return '   '
        if p < 0.001:  return '***'
        if p < 0.01:   return '** '
        if p < 0.05:   return '*  '
        return 'ns '

    def p_str(p):
        if pd.isna(p):  return '       N/A'
        if p < 0.001:   return f'{p:>10.2e}'
        return f'{p:>10.4f}'

    for mark in ['H3K27me3_ctrl', 'H3K27me3_mut']:
        sub = df[df['mark'] == mark]
        print(f'\n{"="*85}')
        print(f' {mark}')
        print(f'{"="*85}')
        print(f'{"Cluster":<10} {"Dir":<13} {"n":>6} {"Exterior":>10} {"Interior":>10} '
              f'{"Ext/Int":>8} {"Asym":>8} {"p-value":>10} {"":>3}')
        print('-' * 85)
        for _, r in sub.iterrows():
            print(f'{r["cluster"]:<10} {r["direction"]:<13} {r["n_anchors"]:>6} '
                  f'{r["exterior_mean"]:>10.4f} {r["interior_mean"]:>10.4f} '
                  f'{r["ext_int_ratio"]:>8.4f} {r["asymmetry_index"]:>+8.5f} '
                  f'{p_str(r["wilcoxon_pval"])} {sig_str(r["wilcoxon_pval"])}')

    for mark in ['H3K27ac_ctrl', 'H3K27ac_mut', 'H2AK119ub_ctrl', 'H2AK119ub_mut',
                 'H3K27me1_ctrl', 'H3K27me1_mut']:
        sub = df[df['mark'] == mark]
        print(f'\n{"="*85}')
        print(f' {mark}')
        print(f'{"="*85}')
        print(f'{"Cluster":<10} {"Dir":<13} {"n":>6} {"Exterior":>10} {"Interior":>10} '
              f'{"Ext/Int":>8} {"Asym":>8} {"p-value":>10} {"":>3}')
        print('-' * 85)
        for _, r in sub.iterrows():
            print(f'{r["cluster"]:<10} {r["direction"]:<13} {r["n_anchors"]:>6} '
                  f'{r["exterior_mean"]:>10.4f} {r["interior_mean"]:>10.4f} '
                  f'{r["ext_int_ratio"]:>8.4f} {r["asymmetry_index"]:>+8.5f} '
                  f'{p_str(r["wilcoxon_pval"])} {sig_str(r["wilcoxon_pval"])}')

    print(f'\n{"="*85}')
    print(' CLUST5 vs CLUST6 COMPARISON (key biological contrast)')
    print(f'{"="*85}')
    for cluster in ['clust5', 'clust6']:
        sub = df[df['cluster'] == cluster]
        print(f'\n  --- {cluster} ({CLUSTER_DIRECTION[cluster]}) ---')
        print(f'  {"Mark":<20} {"Exterior":>10} {"Interior":>10} {"Ext/Int":>8} {"Asym":>8}')
        for _, r in sub.iterrows():
            print(f'  {r["mark"]:<20} {r["exterior_mean"]:>10.4f} {r["interior_mean"]:>10.4f} '
                  f'{r["ext_int_ratio"]:>8.4f} {r["asymmetry_index"]:>+8.5f}')


if __name__ == '__main__':
    main()
