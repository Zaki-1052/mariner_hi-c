# cluster/scripts/13_mechanism_9mark.py
"""Anchor vs Span mechanism figure using the 9-mark 18-state ChromHMM model.

Parallel to the 12-state mechanism figure in 08_summary_figures.py but with
18 states that resolve K119ub-specific biology (Polycomb_K119ub,
Repressed_Enhancer_K119ub, Active_Enhancer_K119ub).

Usage:
  /opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/13_mechanism_9mark.py
"""

import json
import sys
from pathlib import Path
from typing import Dict, List

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

import os as _os
_out = _os.environ.get('CLUSTER_OUT_DIR', 'outputs/bap1_late')
_base = REPO_ROOT / 'cluster' / _out

ANCHOR_18 = _base / 'chromHMM_9mark_intersect/anchor_18.txt'
SPAN_18 = _base / 'chromHMM_9mark_intersect/span_18.txt'
ANNOT_DIR = _base / 'figures/annotation'
OUT_DIR = _base / 'figures/summary_figures'

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as f:
    _custom = json.load(f)
plt.rcParams.update(_custom)
sns.set_theme(style='ticks', rc=_custom)

BIO_ORDER = ['clust6', 'clust3', 'clust1', 'clust2', 'clust4', 'clust5']

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


def extract_top_genes(n=5):
    # type: (int) -> Dict[str, List[str]]
    result = {}  # type: Dict[str, List[str]]
    for clust in BIO_ORDER:
        path = ANNOT_DIR / '{}_annotation.txt'.format(clust)
        if not path.exists():
            result[clust] = []
            continue
        try:
            df = pd.read_csv(path, sep='\t', on_bad_lines='warn')
        except TypeError:
            df = pd.read_csv(path, sep='\t', error_bad_lines=False, warn_bad_lines=True)
        genes = []  # type: List[str]
        for col in ['chr1_gene_name', 'chr2_gene_name']:
            if col in df.columns:
                vals = df[col].dropna().astype(str)
                vals = vals[vals.str.strip() != '']
                genes.extend(vals.tolist())
        if genes:
            counts = pd.Series(genes).value_counts()
            result[clust] = counts.head(n).index.tolist()
        else:
            result[clust] = []
    return result


def make_mechanism_figure_9mark(anchor_enrich, span_enrich, top_genes):
    # type: (pd.DataFrame, pd.DataFrame, Dict[str, List[str]]) -> None
    print('\n=== Mechanism Figure (9-mark, 18-state) ===')

    ordered = [s for s in DISPLAY_ORDER if s in anchor_enrich.index]

    a_df = anchor_enrich.loc[ordered]
    s_df = span_enrich.loc[ordered]

    fig, (ax_gain, ax_loss) = plt.subplots(1, 2, figsize=(11, 7.5), sharey=True)
    fig.subplots_adjust(wspace=0.06, left=0.24, right=0.96, top=0.86, bottom=0.07)

    y = np.arange(len(ordered))
    h = 0.35

    panels = [
        (ax_gain, 'clust5',
         'Gained loops (clust5, n=667)\nPRC1+PRC2 at anchor AND span',
         'Polycomb domain compaction'),
        (ax_loss, 'clust6',
         'Lost loops (clust6, n=2,359)\nK119ub at anchor ONLY',
         'Anchor disruption'),
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
    ax_gain.set_yticks(y)
    ax_gain.set_yticklabels(display_labels, fontsize=6.5)
    ax_loss.tick_params(axis='y', length=0)

    legend_handles = [
        Patch(facecolor='#999999', edgecolor='black', label='Anchor'),
        Patch(facecolor='#999999', edgecolor='black', alpha=0.5, hatch='///', label='Span'),
    ]
    ax_loss.legend(handles=legend_handles, loc='center right', fontsize=7, framealpha=0.9)

    fig.suptitle('ChromHMM Enrichment: Anchor vs Loop Span (9-mark, 18-state)\n'
                 'H2AK119ub resolves PRC1-specific states at BAP1-KO loop anchors',
                 fontsize=10, fontweight='bold', y=0.97)

    for ax, letter in [(ax_gain, 'A'), (ax_loss, 'B')]:
        ax.text(-0.02 if letter == 'B' else -0.18, 1.08, letter,
                transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
        sns.despine(ax=ax)

    sub = figure_subfolder(OUT_DIR, 'mechanism_9mark')
    with multi_format_savefig():
        fig.savefig(str(sub / 'anchor_vs_span_mechanism_9mark.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('  Saved: {}/'.format(sub))


def main():
    print('9-mark mechanism figure')
    print('  Anchor enrichment: {}'.format(ANCHOR_18))
    print('  Span enrichment:   {}'.format(SPAN_18))

    for f in [ANCHOR_18, SPAN_18]:
        if not f.exists():
            raise FileNotFoundError('Missing: {}'.format(f))

    anchor_enrich = load_enrichment(ANCHOR_18)
    span_enrich = load_enrichment(SPAN_18)
    print('  States: {}'.format(anchor_enrich.index.tolist()))

    top_genes = extract_top_genes(n=5)
    make_mechanism_figure_9mark(anchor_enrich, span_enrich, top_genes)

    print('\nDone.')


if __name__ == '__main__':
    main()
