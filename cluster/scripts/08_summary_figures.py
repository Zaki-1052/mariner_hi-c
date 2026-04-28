# cluster/scripts/08_summary_figures.py
"""Lab-meeting summary figures for BAP1-KO loop clustering.

Produces 3 composite figures integrating Phase 4 results:
  1. Cluster Summary Dashboard (6-panel vertical composite)
  2. Anchor vs Span Mechanism (clust5 vs clust6 focused)
  3. Feature Summary Heatmap (z-scored, all clusters x all features)

All data is pre-computed by Phase 4 except H2AK119ub signal
which is queried fresh from BigWig files (~30s).

Usage:
  /opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/08_summary_figures.py
"""

import json
import sys
import warnings
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns
from scipy.stats import zscore as scipy_zscore

SCRIPT_DIR = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT = CLUSTER_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))
from multi_format_output import multi_format_savefig, figure_subfolder

# ---------------------------------------------------------------------------
# Path constants
# ---------------------------------------------------------------------------
CLUSTER_FILE = REPO_ROOT / 'cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt'
METADATA_FILE = REPO_ROOT / 'cluster/data/late_merged_loop_metadata.tsv'
ANCHOR_ENRICH = REPO_ROOT / 'cluster/bap1_late/chromHMM/anchor.txt'
SPAN_ENRICH = REPO_ROOT / 'cluster/bap1_late/chromHMM/span.txt'
PROPORTIONS_FILE = REPO_ROOT / 'cluster/bap1_late/figures/chromHMM_anchor/chromHMM_anchor.proportions.tsv'
DIFF_STATS_FILE = REPO_ROOT / 'cluster/bap1_late/figures/cluster_differential_status.stats.txt'
CLASS_STATS_FILE = REPO_ROOT / 'cluster/bap1_late/figures/loop_classification/loop_classification.stats.txt'
ANNOT_DIR = REPO_ROOT / 'cluster/bap1_late/figures/annotation'
BIGWIG_DIR = Path('/Users/zakiralibhai/sdsc/bigwigs')
OUT_DIR = REPO_ROOT / 'cluster/bap1_late/figures/summary_figures'
COORD_COLS = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']

# ---------------------------------------------------------------------------
# Biological ordering and labels
# ---------------------------------------------------------------------------
BIO_ORDER = ['clust6', 'clust3', 'clust1', 'clust2', 'clust4', 'clust5']
BIO_NAMES = {
    'clust6': 'Strong loss\n(anchor-disrupted)',
    'clust3': 'Moderate loss',
    'clust1': 'Unchanged\n(high signal)',
    'clust2': 'Unchanged',
    'clust4': 'Moderate gain',
    'clust5': 'Strong gain\n(Polycomb domain)',
}
BIO_LABELS = {}  # type: Dict[str, str]

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
LOGFC_POS = '#1f77b4'
LOGFC_NEG = '#d62728'
DIR_COLORS = {
    'down_in_mutant': '#d62728',
    'unchanged':      '#bcbd22',
    'up_in_mutant':   '#2ca02c',
}
CLASS_COLORS = {
    'structural':    '#7f7f7f',
    'CRE':           '#2ca02c',
    'mixed':         '#ff7f0e',
    'unclassified':  '#c7c7c7',
}
POLY_ANCHOR = '#7b2d8e'
POLY_SPAN   = '#c9a0dc'
K119_CTRL   = '#aaaaaa'
K119_MUT    = '#2ca02c'
STATE_COLORS = {
    'Active_Promoter':       '#ff0000',
    'Active_Promoter_Flank': '#ff4500',
    'Poised_Promoter':       '#c83fff',
    'Weak_Promoter':         '#ff69b4',
    'Strong_Enhancer':       '#faca00',
    'Active_Enhancer':       '#ffcd00',
    'Poised_Enhancer':       '#c2e085',
    'Bivalent_Enhancer':     '#bdb76b',
    'Polycomb':              '#808080',
    'CTCF_Boundary':         '#0abab5',
    'Insulator':             '#00bfff',
    'Quiescent':             '#e0e0e0',
}

with open(CLUSTER_DIR / 'modules' / 'custom_params.json') as f:
    _custom = json.load(f)
plt.rcParams.update(_custom)
sns.set_theme(style='ticks', rc=_custom)


# ===================================================================
# Data loaders
# ===================================================================

def load_cluster_data():
    # type: () -> pd.DataFrame
    clusters = pd.read_csv(CLUSTER_FILE, sep='\t')
    meta = pd.read_csv(METADATA_FILE, sep='\t',
                        usecols=COORD_COLS + ['logFC', 'direction'])
    df = clusters.merge(meta, on=COORD_COLS, how='inner')
    df['loop_size_kb'] = (df['y1'] - df['x2']) / 1000.0
    for g in BIO_ORDER:
        n = int((df['GROUP'] == g).sum())
        BIO_LABELS[g] = '{}\n{}, n={:,}'.format(BIO_NAMES[g], g, n)
    print('  Loaded {} loops across {} clusters'.format(len(df), df['GROUP'].nunique()))
    return df


def load_enrichment_tables():
    # type: () -> Tuple[pd.DataFrame, pd.DataFrame]
    dfs = []
    for path in [ANCHOR_ENRICH, SPAN_ENRICH]:
        df = pd.read_csv(path, sep='\t', index_col=0)
        df = df[~df.index.str.startswith('Base')]
        df.columns = [c.replace('.bed', '') for c in df.columns]
        dfs.append(df)
    print('  Loaded enrichment tables: {} states x {} clusters'.format(
        dfs[0].shape[0], dfs[0].shape[1]))
    return dfs[0], dfs[1]


def load_proportions():
    # type: () -> pd.DataFrame
    df = pd.read_csv(PROPORTIONS_FILE, sep='\t', index_col=0)
    return df


def parse_direction_pct():
    # type: () -> pd.DataFrame
    text = DIFF_STATS_FILE.read_text()
    idx = text.find('Row %:\n')
    block = text[idx + len('Row %:\n'):]
    lines = [l for l in block.strip().split('\n') if l.strip()]
    # lines[0] = header with columns.name prefix, lines[1] = index name, lines[2:] = data
    clean = lines[0] + '\n' + '\n'.join(lines[2:])
    df = pd.read_csv(StringIO(clean), sep=r'\s+', index_col=0)
    df.index.name = None
    return df


def parse_classification_pct():
    # type: () -> pd.DataFrame
    text = CLASS_STATS_FILE.read_text()
    idx = text.find('Percent stacked:\n')
    block = text[idx + len('Percent stacked:\n'):]
    lines = [l for l in block.strip().split('\n') if l.strip()]
    clean = '\n'.join(lines)
    df = pd.read_csv(StringIO(clean), sep=r'\s+', index_col=0)
    df.index.name = None
    return df


def compute_k119ub(cluster_df):
    # type: (pd.DataFrame) -> Optional[pd.DataFrame]
    bw_paths = {
        'ctrl': BIGWIG_DIR / 'H2AK119ubCtrl.bw',
        'mut':  BIGWIG_DIR / 'H2AK119ubMut.bw',
    }
    for name, path in bw_paths.items():
        if not path.exists():
            print('  WARNING: {} not found, skipping K119ub'.format(path))
            return None

    a1 = cluster_df[['chr1', 'x1', 'x2']].rename(
        columns={'chr1': 'chrom', 'x1': 'start', 'x2': 'end'})
    a2 = cluster_df[['chr2', 'y1', 'y2']].rename(
        columns={'chr2': 'chrom', 'y1': 'start', 'y2': 'end'})
    all_anchors = pd.concat([a1, a2]).drop_duplicates().reset_index(drop=True)
    print('  {} unique anchors to query'.format(len(all_anchors)))

    def _query_bw(bw_path, anchors):
        # type: (Path, pd.DataFrame) -> Dict[str, float]
        bw = pyBigWig.open(str(bw_path))
        sig = {}
        for _, row in anchors.iterrows():
            key = '{}:{}-{}'.format(row['chrom'], int(row['start']), int(row['end']))
            try:
                s = bw.stats(row['chrom'], int(row['start']), int(row['end']), type='mean')
                sig[key] = s[0] if s and s[0] is not None else np.nan
            except (RuntimeError, ValueError):
                sig[key] = np.nan
        bw.close()
        return sig

    cluster_df = cluster_df.copy()
    cluster_df['_k1'] = cluster_df['chr1'] + ':' + cluster_df['x1'].astype(str) + '-' + cluster_df['x2'].astype(str)
    cluster_df['_k2'] = cluster_df['chr2'] + ':' + cluster_df['y1'].astype(str) + '-' + cluster_df['y2'].astype(str)

    results = {}
    for cond, bw_path in bw_paths.items():
        print('  Querying {}...'.format(bw_path.name))
        sig_map = _query_bw(bw_path, all_anchors)
        s1 = cluster_df['_k1'].map(sig_map)
        s2 = cluster_df['_k2'].map(sig_map)
        cluster_df['_sig'] = (s1 + s2) / 2.0
        results[cond] = cluster_df.groupby('GROUP')['_sig'].median()

    cluster_df.drop(columns=['_k1', '_k2', '_sig'], inplace=True, errors='ignore')
    result = pd.DataFrame(results)
    print('  K119ub medians:\n{}'.format(result.reindex(BIO_ORDER).to_string()))
    return result


def extract_top_genes(n=5):
    # type: (int) -> Dict[str, List[str]]
    result = {}
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


# ===================================================================
# Figure 1: Cluster Summary Dashboard
# ===================================================================

def make_dashboard(cluster_df, direction_pct, class_pct, anchor_enrich,
                   span_enrich, k119ub_df):
    # type: (...) -> None
    print('\n=== [Fig 1] Cluster Summary Dashboard ===')
    x = np.arange(len(BIO_ORDER))
    bar_w = 0.6

    fig = plt.figure(figsize=(6, 16))
    gs = gridspec.GridSpec(6, 1, height_ratios=[0.8, 1.2, 1.5, 1.2, 1.2, 1.2],
                           hspace=0.15, top=0.96, bottom=0.10, left=0.14, right=0.95)
    axes = [fig.add_subplot(gs[i]) for i in range(6)]

    # Panel A: Median logFC
    ax = axes[0]
    logfc = cluster_df.groupby('GROUP')['logFC'].median().reindex(BIO_ORDER)
    colors_a = [LOGFC_NEG if v < 0 else LOGFC_POS for v in logfc]
    ax.bar(x, logfc.values, width=bar_w, color=colors_a, edgecolor='black', linewidth=0.5)
    ax.axhline(0, color='black', linewidth=0.5)
    for xi, v in zip(x, logfc.values):
        offset = 0.02 if v >= 0 else -0.02
        ax.text(xi, v + offset, '{:.2f}'.format(v), ha='center',
                va='bottom' if v >= 0 else 'top', fontsize=6, fontweight='bold')
    ax.set_ylabel('Median logFC', fontsize=8)
    ax.set_title('BAP1-KO Loop Cluster Summary Dashboard',
                 fontsize=10, fontweight='bold', pad=8)

    # Panel B: Direction proportions
    ax = axes[1]
    dir_order = ['down_in_mutant', 'unchanged', 'up_in_mutant']
    dir_labels = {'down_in_mutant': 'Lost', 'unchanged': 'Unchanged', 'up_in_mutant': 'Gained'}
    bottom_b = np.zeros(len(BIO_ORDER))
    for d in dir_order:
        vals = direction_pct.reindex(BIO_ORDER)[d].values
        ax.bar(x, vals, width=bar_w, bottom=bottom_b, color=DIR_COLORS[d],
               edgecolor='white', linewidth=0.3, label=dir_labels[d])
        bottom_b = bottom_b + vals
    ax.set_ylabel('Direction (%)', fontsize=8)
    ax.set_ylim(0, 108)
    ax.legend(fontsize=6, loc='upper right', framealpha=0.9, ncol=3)

    # Panel C: Polycomb anchor vs span — KEY
    ax = axes[2]
    anchor_poly = anchor_enrich.loc['11_Polycomb', BIO_ORDER].values.astype(float)
    span_poly = span_enrich.loc['11_Polycomb', BIO_ORDER].values.astype(float)
    w = 0.28
    bars_anc = ax.bar(x - w/2, anchor_poly, width=w, color=POLY_ANCHOR,
                       edgecolor='black', linewidth=0.5, label='Anchor', zorder=3)
    bars_spn = ax.bar(x + w/2, span_poly, width=w, color=POLY_SPAN,
                       edgecolor='black', linewidth=0.5, hatch='///',
                       label='Span', zorder=3)
    ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.8, alpha=0.7, zorder=1)
    ax.set_ylabel('Polycomb fold\nenrichment', fontsize=8)
    for xi, va, vs in zip(x, anchor_poly, span_poly):
        ax.text(xi - w/2, va + 0.15, '{:.1f}x'.format(va), ha='center', va='bottom',
                fontsize=5.5, fontweight='bold', color=POLY_ANCHOR)
        ax.text(xi + w/2, vs + 0.15, '{:.1f}x'.format(vs), ha='center', va='bottom',
                fontsize=5.5, color='#666666')
    ax.legend(fontsize=6, loc='upper left', framealpha=0.9)
    ax.text(0.98, 0.95, 'KEY RESULT', transform=ax.transAxes, fontsize=7,
            fontweight='bold', ha='right', va='top', color=POLY_ANCHOR,
            bbox=dict(boxstyle='round,pad=0.2', facecolor='#f3e5f5', alpha=0.8))

    # Panel D: Loop size
    ax = axes[3]
    size_data = [cluster_df.loc[cluster_df['GROUP'] == g, 'loop_size_kb'].values
                 for g in BIO_ORDER]
    bp = ax.boxplot(size_data, positions=x, widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(linewidth=0.7),
                    capprops=dict(linewidth=0.7))
    for patch in bp['boxes']:
        patch.set_facecolor('#6baed6')
        patch.set_edgecolor('black')
        patch.set_linewidth(0.5)
    meds = [np.median(d) for d in size_data]
    for xi, m in zip(x, meds):
        ax.text(xi, m + 30, '{:.0f}'.format(m), ha='center', va='bottom',
                fontsize=5.5, fontweight='bold')
    ax.set_ylabel('Loop size (kb)', fontsize=8)

    # Panel E: Loop classification
    ax = axes[4]
    class_order = ['structural', 'CRE', 'mixed', 'unclassified']
    bottom_e = np.zeros(len(BIO_ORDER))
    for ct in class_order:
        vals = class_pct.loc[ct, BIO_ORDER].values.astype(float)
        ax.bar(x, vals, width=bar_w, bottom=bottom_e, color=CLASS_COLORS[ct],
               edgecolor='white', linewidth=0.3, label=ct)
        bottom_e = bottom_e + vals
    ax.set_ylabel('Loop type (%)', fontsize=8)
    ax.set_ylim(0, 108)
    ax.legend(fontsize=6, loc='upper right', framealpha=0.9, ncol=2)

    # Panel F: K119ub signal
    ax = axes[5]
    if k119ub_df is not None:
        ctrl_v = k119ub_df.reindex(BIO_ORDER)['ctrl'].values
        mut_v = k119ub_df.reindex(BIO_ORDER)['mut'].values
        ax.bar(x - w/2, ctrl_v, width=w, color=K119_CTRL, edgecolor='black',
               linewidth=0.5, label='Ctrl')
        ax.bar(x + w/2, mut_v, width=w, color=K119_MUT, edgecolor='black',
               linewidth=0.5, label='Mut (BAP1-KO)')
        ax.legend(fontsize=6, loc='upper left', framealpha=0.9)
    else:
        ax.text(0.5, 0.5, 'K119ub BigWigs not found', transform=ax.transAxes,
                ha='center', va='center', fontsize=8, color='grey')
    ax.set_ylabel('H2AK119ub\nanchor signal', fontsize=8)

    # Shared x-axis: labels only on bottom
    for i in range(5):
        axes[i].set_xticks(x)
        axes[i].set_xticklabels([])
    axes[5].set_xticks(x)
    axes[5].set_xticklabels([BIO_LABELS[g] for g in BIO_ORDER],
                             rotation=45, ha='right', fontsize=6)

    for i, a in enumerate(axes):
        a.text(-0.06, 1.06, chr(65 + i), transform=a.transAxes,
               fontsize=12, fontweight='bold', va='top')
        a.tick_params(axis='both', labelsize=7)
        sns.despine(ax=a)

    sub = figure_subfolder(OUT_DIR, 'dashboard')
    with multi_format_savefig():
        fig.savefig(str(sub / 'cluster_dashboard.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('  Saved: {}/'.format(sub))


# ===================================================================
# Figure 2: Anchor vs Span Mechanism
# ===================================================================

def make_mechanism_figure(anchor_enrich, span_enrich, top_genes):
    # type: (...) -> None
    print('\n=== [Fig 2] Anchor vs Span Mechanism ===')

    states = anchor_enrich.index.tolist()
    state_names = []
    for s in states:
        parts = s.split('_', 1)
        state_names.append(parts[1] if len(parts) == 2 else s)

    fig, (ax_gain, ax_loss) = plt.subplots(1, 2, figsize=(10, 6), sharey=True)
    fig.subplots_adjust(wspace=0.06, left=0.20, right=0.96, top=0.86, bottom=0.07)

    y = np.arange(len(states))
    h = 0.35

    panels = [
        (ax_gain, 'clust5',
         'Gained loops (clust5, n=667)\nPolycomb at anchor AND span',
         'Domain compaction'),
        (ax_loss, 'clust6',
         'Lost loops (clust6, n=2,359)\nPolycomb at anchor ONLY',
         'Anchor disruption'),
    ]

    for ax, clust, title, mech in panels:
        a_vals = anchor_enrich[clust].values.astype(float)
        s_vals = span_enrich[clust].values.astype(float)

        for yi, sn in enumerate(state_names):
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

    ax_gain.set_yticks(y)
    ax_gain.set_yticklabels(state_names, fontsize=7)
    ax_loss.tick_params(axis='y', length=0)

    legend_handles = [
        Patch(facecolor='#999999', edgecolor='black', label='Anchor'),
        Patch(facecolor='#999999', edgecolor='black', alpha=0.5, hatch='///', label='Span'),
    ]
    ax_loss.legend(handles=legend_handles, loc='center right', fontsize=7, framealpha=0.9)

    fig.suptitle('ChromHMM Enrichment: Anchor vs Loop Span\n'
                 'Two mechanisms of BAP1-KO loop dysregulation',
                 fontsize=10, fontweight='bold', y=0.97)

    for ax, letter in [(ax_gain, 'A'), (ax_loss, 'B')]:
        ax.text(-0.02 if letter == 'B' else -0.15, 1.08, letter,
                transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
        sns.despine(ax=ax)

    sub = figure_subfolder(OUT_DIR, 'mechanism')
    with multi_format_savefig():
        fig.savefig(str(sub / 'anchor_vs_span_mechanism.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('  Saved: {}/'.format(sub))


# ===================================================================
# Figure 3: Feature Summary Heatmap
# ===================================================================

def make_feature_heatmap(cluster_df, direction_pct, class_pct, anchor_enrich,
                         span_enrich, proportions, k119ub_df, top_genes):
    # type: (...) -> None
    print('\n=== [Fig 3] Feature Summary Heatmap ===')

    feat = {}  # type: Dict[str, np.ndarray]

    feat['Median\nlogFC'] = (cluster_df.groupby('GROUP')['logFC'].median()
                             .reindex(BIO_ORDER).values)
    feat['%\nGained'] = direction_pct.reindex(BIO_ORDER)['up_in_mutant'].values
    feat['%\nLost'] = direction_pct.reindex(BIO_ORDER)['down_in_mutant'].values

    feat['Median\nsize (kb)'] = (cluster_df.groupby('GROUP')['loop_size_kb'].median()
                                  .reindex(BIO_ORDER).values)

    feat['%\nStructural'] = class_pct.loc['structural', BIO_ORDER].values.astype(float)
    feat['%\nCRE'] = class_pct.loc['CRE', BIO_ORDER].values.astype(float)

    feat['Polycomb\nanchor'] = anchor_enrich.loc['11_Polycomb', BIO_ORDER].values.astype(float)
    feat['Polycomb\nspan'] = span_enrich.loc['11_Polycomb', BIO_ORDER].values.astype(float)
    feat['Bivalent\nanchor'] = anchor_enrich.loc['12_Bivalent_Enhancer', BIO_ORDER].values.astype(float)

    if 'Polycomb' in proportions.index:
        feat['% Polycomb\n(anchors)'] = proportions.loc['Polycomb', BIO_ORDER].values.astype(float)

    if k119ub_df is not None:
        feat['K119ub\nctrl'] = k119ub_df.reindex(BIO_ORDER)['ctrl'].values
        feat['K119ub\nmut'] = k119ub_df.reindex(BIO_ORDER)['mut'].values

    raw_df = pd.DataFrame(feat, index=BIO_ORDER)

    z_vals = raw_df.apply(lambda col: scipy_zscore(col.astype(float)), axis=0)
    z_vals = z_vals.clip(-2.5, 2.5)

    annot_strs = raw_df.copy()
    for col in annot_strs.columns:
        if 'size' in col.lower():
            annot_strs[col] = annot_strs[col].apply(lambda v: '{:.0f}'.format(v))
        elif '%' in col:
            annot_strs[col] = annot_strs[col].apply(lambda v: '{:.1f}'.format(v))
        else:
            annot_strs[col] = annot_strs[col].apply(lambda v: '{:.2f}'.format(v))

    row_labels = ['{}\n({})'.format(BIO_NAMES[g].split('\n')[0], g) for g in BIO_ORDER]

    ncols = len(raw_df.columns)
    fig_w = max(8, ncols * 0.75 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, 4.5))

    sns.heatmap(z_vals.values.astype(float), ax=ax, cmap='RdBu_r', center=0,
                vmin=-2.5, vmax=2.5,
                annot=annot_strs.values, fmt='', annot_kws={'fontsize': 6},
                xticklabels=raw_df.columns.tolist(),
                yticklabels=row_labels,
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'Z-score', 'shrink': 0.7})

    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha='right', fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
    ax.set_title('Cluster Feature Summary\n(z-score normalized, raw values annotated)',
                 fontsize=10, fontweight='bold', pad=10)

    for i, g in enumerate(BIO_ORDER):
        genes = top_genes.get(g, [])[:3]
        if genes:
            ax.text(ncols + 0.3, i + 0.5, ', '.join(genes), fontsize=5,
                    va='center', ha='left', fontstyle='italic', color='#555555')

    sub = figure_subfolder(OUT_DIR, 'heatmap')
    with multi_format_savefig():
        fig.savefig(str(sub / 'feature_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('  Saved: {}/'.format(sub))

    raw_df.to_csv(str(sub / 'feature_values.tsv'), sep='\t')
    print('  Raw values: {}/feature_values.tsv'.format(sub))


# ===================================================================
# Main
# ===================================================================

def main():
    print('=' * 60)
    print('Summary Figures for Lab Meeting')
    print('=' * 60)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print('\n--- Loading pre-computed data ---')
    cluster_df = load_cluster_data()
    anchor_enrich, span_enrich = load_enrichment_tables()
    proportions = load_proportions()
    direction_pct = parse_direction_pct()
    class_pct = parse_classification_pct()
    top_genes = extract_top_genes(n=5)

    print('\n--- Computing K119ub signal (~30s) ---')
    k119ub_df = compute_k119ub(cluster_df)

    print('\n--- Rendering figures ---')
    make_dashboard(cluster_df, direction_pct, class_pct, anchor_enrich,
                   span_enrich, k119ub_df)
    make_mechanism_figure(anchor_enrich, span_enrich, top_genes)
    make_feature_heatmap(cluster_df, direction_pct, class_pct, anchor_enrich,
                         span_enrich, proportions, k119ub_df, top_genes)

    print('\n' + '=' * 60)
    print('All summary figures complete.')
    print('Output: {}/'.format(OUT_DIR))
    print('=' * 60)


if __name__ == '__main__':
    main()
