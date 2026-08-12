# cluster/scripts/14_mechanism_9mark_clust6_split.py
"""Anchor vs Span mechanism figure for clust6 short/long subgroups (9-mark model).

Splits clust6 (lost loops) at the 800kb size threshold into short (<800kb) and
long (>=800kb) subgroups, generates per-subgroup ChromHMM OverlapEnrichment,
and produces a 2-panel mechanism figure comparing the two subgroups.

Companion to 13_mechanism_9mark.py which shows the full clust5 vs clust6.

Usage:
  /opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/14_mechanism_9mark_clust6_split.py
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT = CLUSTER_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))
from multi_format_output import multi_format_savefig, figure_subfolder
from pipeline_config import get_size_threshold

import os as _os
_out = _os.environ.get('CLUSTER_OUT_DIR', 'outputs/bap1_late')
_base = REPO_ROOT / 'cluster' / _out

CLUSTER_FILE = _base / 'cluster3/k-6/data/combined-clusters.txt'
CHROMHMM_9M = _base / 'chromHMM_9mark_intersect'
SEGMENT_BED = CHROMHMM_9M / 'learned_model_18/cerebellum_late_18_segments.bed'
RENAME_FILE = CHROMHMM_9M / '18state_rename_cerebellum.txt'
ANCHOR_INPUT = CHROMHMM_9M / 'anchor_input_split'
SPAN_INPUT = CHROMHMM_9M / 'span_input_split'
ANCHOR_SPLIT = CHROMHMM_9M / 'anchor_split_18.txt'
SPAN_SPLIT = CHROMHMM_9M / 'span_split_18.txt'
ANNOT_DIR = _base / 'figures/annotation'
OUT_DIR = _base / 'figures/summary_figures'

SIZE_THRESHOLD = get_size_threshold()

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as f:
    _custom = json.load(f)
plt.rcParams.update(_custom)
sns.set_theme(style='ticks', rc=_custom)

SPLIT_GROUPS = ['clust6_short', 'clust6_long']

STATE_COLORS = {
    'K119ub_Only':               '#663399',
    'Polycomb_K119ub':           '#4a0e4e',
    'Repressed_Enhancer_K119ub': '#6b3a6b',
    'K119ub_Poised_Enhancer':    '#8b5e83',
    'Poised_Enhancer':           '#c2e085',
    'Active_Enhancer_K9ac':      '#3cb371',
    'Active_Enhancer':           '#ffcd00',
    'Active_Enhancer_K119ub':    '#556b2f',
    'Strong_Enhancer':           '#faca00',
    'Weak_Enhancer':             '#dda0dd',
    'Active_Promoter':           '#ff0000',
    'K9ac_Promoter':             '#e3877c',
    'Quiescent':                 '#e0e0e0',
    'CTCF_Open':                 '#a9a9a9',
    'ATAC_Enhancer':             '#a8e6cf',
    'Insulator':                 '#00bfff',
    'Polycomb':                  '#808080',
    'Heterochromatin':           '#1a1a2e',
}

DISPLAY_ORDER = [
    'Active_Promoter',
    'K9ac_Promoter',
    'Active_Enhancer',
    'Active_Enhancer_K119ub',
    'Active_Enhancer_K9ac',
    'Strong_Enhancer',
    'Weak_Enhancer',
    'Poised_Enhancer',
    'K119ub_Poised_Enhancer',
    'Repressed_Enhancer_K119ub',
    'Polycomb_K119ub',
    'Polycomb',
    'K119ub_Only',
    'Heterochromatin',
    'Insulator',
    'CTCF_Open',
    'ATAC_Enhancer',
    'Quiescent',
]


def build_split_beds():
    # type: () -> Tuple[Dict[str, int], List[str], List[str]]
    """Split clust6 by loop size and write anchor/span BEDs.

    Returns (counts, anchor_filenames, span_filenames).
    """
    df = pd.read_csv(CLUSTER_FILE, sep='\t')
    clust6 = df[df['GROUP'] == 'clust6'].copy()
    clust6['loop_size'] = clust6['y2'] - clust6['x1']

    assert (clust6['chr1'] == clust6['chr2']).all(), 'Trans loops in clust6'

    short = clust6[clust6['loop_size'] < SIZE_THRESHOLD]
    long = clust6[clust6['loop_size'] >= SIZE_THRESHOLD]
    groups = {'clust6_short': short, 'clust6_long': long}
    counts = {name: len(sub) for name, sub in groups.items()}

    anchor_fns = []  # type: List[str]
    span_fns = []  # type: List[str]

    for kind, out_dir in [('anchor', ANCHOR_INPUT), ('span', SPAN_INPUT)]:
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, sub in groups.items():
            bed_path = out_dir / '{}.bed'.format(name)
            if kind == 'anchor':
                out_df = pd.concat([
                    sub[['chr1', 'x1', 'x2']],
                    sub[['chr2', 'y1', 'y2']].rename(
                        columns={'chr2': 'chr1', 'y1': 'x1', 'y2': 'x2'}),
                ], ignore_index=True)
            else:
                out_df = sub[['chr1', 'x1', 'y2']].copy()
            out_df.to_csv(bed_path, sep='\t', header=False, index=False)
            fn = '{}.bed'.format(name)
            if kind == 'anchor':
                anchor_fns.append(fn)
            else:
                span_fns.append(fn)

    print('  clust6_short: {} loops'.format(counts['clust6_short']))
    print('  clust6_long:  {} loops'.format(counts['clust6_long']))
    return counts, anchor_fns, span_fns


def run_overlap_enrichment(coord_kind, filenames, input_dir, output_prefix):
    # type: (str, List[str], Path, Path) -> Path
    out_txt = Path(str(output_prefix) + '.txt')
    if out_txt.exists() and out_txt.stat().st_size > 0:
        print('  Skipping {} OverlapEnrichment (output exists)'.format(coord_kind))
        return out_txt

    coordlist = CHROMHMM_9M / 'coordlistfile_{}_split.txt'.format(coord_kind)
    coordlist.write_text('\n'.join(filenames) + '\n')

    cmd = [
        'chromhmm', 'OverlapEnrichment',
        '-noimage', '-uniformscale',
        '-m', str(RENAME_FILE),
        '-f', str(coordlist),
        '-colfields', '0,1,2',
        str(SEGMENT_BED),
        str(input_dir) + '/',
        str(output_prefix),
    ]
    print('  Running: {}'.format(' '.join(cmd)))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print('  STDOUT:', proc.stdout)
        print('  STDERR:', proc.stderr)
        raise RuntimeError(
            'OverlapEnrichment failed for {} (exit {})'.format(coord_kind, proc.returncode))
    if not out_txt.exists() or out_txt.stat().st_size == 0:
        raise RuntimeError('OverlapEnrichment produced no output: {}'.format(out_txt))
    return out_txt


def load_enrichment(path):
    # type: (Path) -> pd.DataFrame
    df = pd.read_csv(path, sep='\t', index_col=0)
    df = df[~df.index.str.startswith('Base')]
    df.columns = [c.replace('.bed', '') for c in df.columns]
    names = []
    for s in df.index:
        parts = s.split('_', 1)
        names.append(parts[1] if len(parts) == 2 else s)
    df.index = names
    return df


def extract_top_genes(cluster_file, n=5):
    # type: (Path, int) -> Dict[str, List[str]]
    """Extract top genes per clust6 subgroup by joining annotation with the split."""
    annot_path = ANNOT_DIR / 'clust6_annotation.txt'
    if not annot_path.exists():
        return {'clust6_short': [], 'clust6_long': []}
    try:
        annot = pd.read_csv(annot_path, sep='\t', on_bad_lines='warn')
    except TypeError:
        annot = pd.read_csv(annot_path, sep='\t', error_bad_lines=False, warn_bad_lines=True)
    if annot.empty:
        return {'clust6_short': [], 'clust6_long': []}

    clusters = pd.read_csv(cluster_file, sep='\t')
    clust6 = clusters[clusters['GROUP'] == 'clust6'].copy()
    clust6['loop_size'] = clust6['y2'] - clust6['x1']
    clust6['subgroup'] = np.where(
        clust6['loop_size'] < SIZE_THRESHOLD, 'clust6_short', 'clust6_long')

    coord_cols = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']
    merged = annot.merge(clust6[coord_cols + ['subgroup']], on=coord_cols, how='inner')

    result = {}  # type: Dict[str, List[str]]
    for sg in SPLIT_GROUPS:
        sub = merged[merged['subgroup'] == sg]
        genes = []  # type: List[str]
        for col in ['chr1_gene_name', 'chr2_gene_name']:
            if col in sub.columns:
                vals = sub[col].dropna().astype(str)
                vals = vals[vals.str.strip() != '']
                genes.extend(vals.tolist())
        if genes:
            counts = pd.Series(genes).value_counts()
            result[sg] = counts.head(n).index.tolist()
        else:
            result[sg] = []
    return result


def make_mechanism_figure_9mark_split(anchor_enrich, span_enrich, top_genes,
                                      counts):
    # type: (pd.DataFrame, pd.DataFrame, Dict[str, List[str]], Dict[str, int]) -> None
    print('\n=== Mechanism Figure (9-mark, clust6 short vs long) ===')

    ordered = [s for s in DISPLAY_ORDER if s in anchor_enrich.index]

    a_df = anchor_enrich.loc[ordered]
    s_df = span_enrich.loc[ordered]

    fig, (ax_short, ax_long) = plt.subplots(1, 2, figsize=(11, 7.5), sharey=True)
    fig.subplots_adjust(wspace=0.06, left=0.24, right=0.96, top=0.86, bottom=0.07)

    y = np.arange(len(ordered))
    h = 0.35

    n_short = '{:,}'.format(counts['clust6_short'])
    n_long = '{:,}'.format(counts['clust6_long'])
    thresh_kb = SIZE_THRESHOLD // 1000

    panels = [
        (ax_short, 'clust6_short',
         'Lost loops (clust6 short, n={})\n<{}kb, K119ub at active enhancers'.format(
             n_short, thresh_kb),
         'Active anchor disruption'),
        (ax_long, 'clust6_long',
         'Lost loops (clust6 long, n={})\n≥{}kb, PRC1+PRC2 at anchors'.format(
             n_long, thresh_kb),
         'Polycomb anchor silencing'),
    ]

    for ax, clust, title, mech in panels:
        a_vals = a_df[clust].values.astype(float)
        s_vals = s_df[clust].values.astype(float)

        for yi, sn in enumerate(ordered):
            color = STATE_COLORS.get(sn, '#cccccc')
            ax.barh(yi + h/2, a_vals[yi], height=h, color=color,
                    edgecolor='#333333', linewidth=0.4)
            ax.barh(yi - h/2, s_vals[yi], height=h, color=color,
                    edgecolor='#333333', linewidth=0.4, alpha=0.5, hatch='///')

        ax.axvline(1.0, color='grey', linestyle='--', linewidth=0.8, alpha=0.7, zorder=0)
        ax.set_title(title, fontsize=8.5, fontweight='bold', pad=8)
        ax.set_xlabel('Fold enrichment over genome', fontsize=8)
        ax.tick_params(axis='x', labelsize=7)

        ax.text(0.97, 0.04, mech, transform=ax.transAxes, fontsize=8,
                fontstyle='italic', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                          edgecolor='#cccccc', alpha=0.9))

        genes = top_genes.get(clust, [])[:3]
        if genes:
            gene_str = 'Top genes: ' + ', '.join(genes)
            ax.text(0.97, 0.14, gene_str, transform=ax.transAxes, fontsize=6,
                    ha='right', va='bottom', fontstyle='italic',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                              edgecolor='#cccccc', alpha=0.9))

    display_labels = [s.replace('_', ' ') for s in ordered]
    ax_short.set_yticks(y)
    ax_short.set_yticklabels(display_labels, fontsize=6.5)
    ax_long.tick_params(axis='y', length=0)

    legend_handles = [
        Patch(facecolor='#999999', edgecolor='black', label='Anchor'),
        Patch(facecolor='#999999', edgecolor='black', alpha=0.5, hatch='///', label='Span'),
    ]
    ax_long.legend(handles=legend_handles, loc='center right', fontsize=7, framealpha=0.9)

    fig.suptitle(
        'ChromHMM Enrichment: Anchor vs Loop Span (9-mark, 18-state)\n'
        'Clust6 lost loops split by loop size at {}kb'.format(thresh_kb),
        fontsize=10, fontweight='bold', y=0.97)

    for ax, letter in [(ax_short, 'A'), (ax_long, 'B')]:
        ax.text(-0.02 if letter == 'B' else -0.18, 1.08, letter,
                transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
        sns.despine(ax=ax)

    sub = figure_subfolder(OUT_DIR, 'mechanism_9mark_clust6_split')
    with multi_format_savefig():
        fig.savefig(str(sub / 'anchor_vs_span_mechanism_9mark_clust6_split.png'),
                    dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('  Saved: {}/'.format(sub))


def main():
    print('9-mark mechanism figure -- clust6 short/long split')
    print('  Size threshold: {} bp ({}kb)'.format(SIZE_THRESHOLD, SIZE_THRESHOLD // 1000))
    print('  Cluster file:   {}'.format(CLUSTER_FILE))
    print('  Segment BED:    {}'.format(SEGMENT_BED))
    print('  Rename file:    {}'.format(RENAME_FILE))

    for f in [CLUSTER_FILE, SEGMENT_BED, RENAME_FILE]:
        if not f.exists():
            raise FileNotFoundError('Missing: {}'.format(f))

    if not shutil.which('chromhmm'):
        raise RuntimeError('chromhmm not found on PATH')

    counts, anchor_fns, span_fns = build_split_beds()

    anchor_txt = run_overlap_enrichment(
        'anchor', anchor_fns, ANCHOR_INPUT,
        CHROMHMM_9M / 'anchor_split_18')
    span_txt = run_overlap_enrichment(
        'span', span_fns, SPAN_INPUT,
        CHROMHMM_9M / 'span_split_18')

    anchor_enrich = load_enrichment(anchor_txt)
    span_enrich = load_enrichment(span_txt)
    print('  States: {}'.format(anchor_enrich.index.tolist()))

    top_genes = extract_top_genes(CLUSTER_FILE, n=5)
    make_mechanism_figure_9mark_split(anchor_enrich, span_enrich, top_genes, counts)

    print('\nDone.')


if __name__ == '__main__':
    main()
