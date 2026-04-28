#!/usr/bin/env python3
# cluster/scripts/05_grouped_analyses.py
"""Phase 4 grouped downstream analyses for BAP1-KO Hi-C loop clusters.

Adapts Popay grouped_loops_figures.ipynb for cerebellum data. Eight sub-analyses
run on the canonical 6-cluster solution from Phase 3 (combined-clusters.txt).
Sub-analyses are gated by --analyses (comma-separated; default 'all').

Sub-analyses:
  4.1  loop_size           — per-cluster loop-size distribution (Popay cell-3)
  4.2  loop_classification — CTCF-only structural / CRE / mixed / unclassified
  4.3  anchor_chip         — raw RPKM at anchors per cluster (4 marks x ctrl/mut)
  4.4  chromhmm_span       — Fig 2f-equivalent anchor-vs-span enrichment (KEY)
  4.5  chromhmm_proportions — most-prominent state per anchor, grouped bar
  4.6  gene_annotation     — bedtools-based per-cluster gene lists
  4.7  diffbind            — DiffBind peak overlap with cluster anchors (3 marks
                              x 2 thresholds: project 0.3 + Popay 0.0)
  4.8  cluster_diff        — cluster x edgeR direction crosstab (NEW)
"""
import argparse
import shutil
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use('Agg')

import bioframe as bf
import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import chi2_contingency

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR / 'modules'))  # bedpe_analysis, plotting, etc.
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))  # multi_format_output

import bedpe_analysis        # noqa: E402
import plotting              # noqa: E402
import statistics_functions  # noqa: E402
from cluster_tools import sort_by_strings  # noqa: E402
from chromHMM_heatmap import heatmap_plot  # noqa: E402
from multi_format_output import multi_format_savefig, figure_subfolder  # noqa: E402

# -------- inputs --------
CLUSTER_FILE  = REPO_ROOT / 'cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt'
METADATA_FILE = REPO_ROOT / 'cluster/data/late_merged_loop_metadata.tsv'
PROMOTER_BED  = REPO_ROOT / 'cluster/data/mm10_knownGene_pp.bed'
SEGMENT_BED   = REPO_ROOT / 'cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed'
RENAME_FILE   = REPO_ROOT / 'cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt'

CTCF_BED      = REPO_ROOT / 'peaks/CTCF.bed'
ENHANCER_BED  = REPO_ROOT / 'peaks/beds/H3K27acCerebellumLate2.bed'

DIFFBIND_FILES = {
    'K27ac':  REPO_ROOT / 'peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt',
    'K27me3': REPO_ROOT / 'peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt',
    'K119ub': REPO_ROOT / 'peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt',
}

# All BigWigs sourced from sdsc — peaks/bigwigs/macs2.narrow.aug18.dedup mut files
# are partially 0-byte (corrupted), per data-format audit. sdsc has all 16 marks.
BIGWIG_BASE = Path('/Users/zakiralibhai/sdsc/bigwigs')
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

# -------- outputs --------
OUT_BASE     = REPO_ROOT / 'cluster/bap1_late'
FIG_BASE     = OUT_BASE / 'figures'
CHROMHMM_DIR = OUT_BASE / 'chromHMM'

# -------- constants --------
CLUSTER_ORDER = ['clust1', 'clust2', 'clust3', 'clust4', 'clust5', 'clust6']
COORD_COLS    = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']

# 12-state ChromHMM cerebellum model (E-prefix). Quiescent excluded from
# proportions plots — 92% of genome in this state, would dominate.
STATE_COLORS = {
    'Active_Promoter':       '#d62728',
    'Active_Promoter_Flank': '#ff9896',
    'Poised_Promoter':       '#ff7f0e',
    'Weak_Promoter':         '#ffbb78',
    'Strong_Enhancer':       '#9467bd',
    'Active_Enhancer':       '#2ca02c',
    'Poised_Enhancer':       '#98df8a',
    'Bivalent_Enhancer':     '#bcbd22',
    'Polycomb':              '#17becf',
    'CTCF_Boundary':         '#7f7f7f',
    'Insulator':             '#c7c7c7',
}

DIRECTION_COLORS = {
    'down_in_mutant': '#d62728',
    'unchanged':      '#bcbd22',
    'up_in_mutant':   '#2ca02c',
}

CLASSIFICATION_COLORS = {
    'structural':   '#7f7f7f',
    'CRE':          '#2ca02c',
    'mixed':        '#ff7f0e',
    'unclassified': '#c7c7c7',
}

DIFFBIND_COLORS = {
    'no peak':   '#c7c7c7',
    'non-sig.':  '#7f7f7f',
    'decreased': '#1f77b4',
    'increased': '#d62728',
}

VALID_ANALYSES = {'4.1', '4.2', '4.3', '4.4', '4.5', '4.6', '4.7', '4.8'}

# Default suspended layout for 4.3 anchor-ChIP boxplot (8 panels, 4 marks x 2 conditions)
ANCHOR_CHIP_PANEL_ORDER = list(BIGWIG_DICT.keys())


# ============================================================
# Shared loaders
# ============================================================
def load_clusters() -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame], List[str]]:
    """Load combined-clusters.txt → sorted bedpe_df + group_dict (coord-only) + actual_order."""
    if not CLUSTER_FILE.exists():
        raise FileNotFoundError(f'Cluster file not found: {CLUSTER_FILE}')
    df = pd.read_csv(CLUSTER_FILE, sep='\t')
    bedpe_df = sort_by_strings(df=df, sort_col='GROUP', order=CLUSTER_ORDER)
    actual_order = [g for g in CLUSTER_ORDER if g in bedpe_df['GROUP'].unique()]
    group_dict = {}  # type: Dict[str, pd.DataFrame]
    for g in actual_order:
        sub = bedpe_df[bedpe_df['GROUP'] == g][COORD_COLS].drop_duplicates(ignore_index=True)
        group_dict[g] = sub
    print(f'Loaded {len(bedpe_df)} loops in {len(actual_order)} clusters: '
          + ', '.join(f'{g}={len(group_dict[g])}' for g in actual_order))
    return bedpe_df, group_dict, actual_order


def load_state_rename() -> pd.DataFrame:
    """Load 12state_rename → sort by state-id integer (regex strips E-prefix)."""
    df = pd.read_csv(RENAME_FILE, sep='\t', header=None, names=['short', 'state'])
    df['sort_col'] = df['short'].str.replace(r'^[A-Z]', '', regex=True).astype(int)
    df.sort_values('sort_col', inplace=True, ignore_index=True)
    return df


def parse_analyses(arg: str) -> List[str]:
    if arg == 'all':
        return sorted(VALID_ANALYSES)
    requested = {x.strip() for x in arg.split(',') if x.strip()}
    unknown = requested - VALID_ANALYSES
    if unknown:
        raise SystemExit(
            f'Unknown analyses: {", ".join(sorted(unknown))}. '
            f'Valid: {", ".join(sorted(VALID_ANALYSES))}, all'
        )
    return sorted(requested)


# ============================================================
# 4.1 Loop size
# ============================================================
def run_loop_size(group_dict: dict) -> None:
    print('\n=== [4.1] Loop size per cluster ===')
    sub = figure_subfolder(FIG_BASE, 'loop_size')
    with multi_format_savefig():
        bedpe_analysis.loop_size(out_dir=str(sub), bedpe_dict=group_dict, logY=False)
    print(f'  Outputs: {sub}/loop_size{{,_strip}}.{{png,pdf,svg,jpg}} + loop_size.stats.txt')


# ============================================================
# 4.2 Loop classification (CTCF-only, no RAD21)
# ============================================================
def run_loop_classification(group_dict: dict) -> None:
    """CTCF + enhancer + promoter overlap → structural / CRE / mixed / unclassified.

    Adapts Popay cell-5: drop RAD21 from classifier (we have no RAD21 ChIP);
    use CTCF-alone for structural rule.
    """
    print('\n=== [4.2] Loop classification (CTCF-only) ===')
    sub = figure_subfolder(FIG_BASE, 'loop_classification')

    classifier_dict = {
        'CTCF':     pd.read_csv(CTCF_BED,     sep='\t', header=None,
                                names=['chr', 'start', 'end'], usecols=[0, 1, 2]),
        'enhancer': pd.read_csv(ENHANCER_BED, sep='\t', header=None,
                                names=['chr', 'start', 'end'], usecols=[0, 1, 2]),
        'promoter': pd.read_csv(PROMOTER_BED, sep='\t', header=None,
                                names=['chr', 'start', 'end'], usecols=[0, 1, 2]),
    }

    classified_df = pd.DataFrame()
    for loop_group, loop_df in group_dict.items():
        loop_df = loop_df.copy()
        for anchor in [['chr1', 'x1', 'x2'], ['chr2', 'y1', 'y2']]:
            for classifier, class_df in classifier_dict.items():
                ov = bf.overlap(df1=loop_df[anchor], df2=class_df,
                                cols1=anchor, cols2=['chr', 'start', 'end'])
                col_name = anchor[0] + '_' + classifier
                ov[col_name] = (ov['start_'].notnull()).astype(int)
                ov = (ov[anchor + [col_name]]
                      .sort_values(by=col_name, ascending=False)
                      .drop_duplicates(subset=anchor))
                loop_df = loop_df.merge(ov, on=anchor, how='left')

        loop_df['chr1_EorP'] = loop_df[['chr1_promoter', 'chr1_enhancer']].max(axis=1)
        loop_df['chr2_EorP'] = loop_df[['chr2_promoter', 'chr2_enhancer']].max(axis=1)
        for c in ['EorP', 'CTCF']:
            loop_df[c] = loop_df[['chr1_' + c, 'chr2_' + c]].sum(axis=1)
        loop_df = loop_df[['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'EorP', 'CTCF']]

        # CTCF-only structural rule (Popay's needed CTCF==2 AND RAD21==2)
        loop_df['classification'] = 'unclassified'
        loop_df.loc[(loop_df['CTCF'] == 2) & (loop_df['EorP'] < 2), 'classification'] = 'structural'
        loop_df.loc[(loop_df['CTCF'] < 2)  & (loop_df['EorP'] == 2), 'classification'] = 'CRE'
        loop_df.loc[(loop_df['CTCF'] == 2) & (loop_df['EorP'] == 2), 'classification'] = 'mixed'

        for c in ['structural', 'CRE', 'mixed', 'unclassified']:
            classified_df.loc[c, loop_group] = (loop_df['classification'] == c).sum()

    print('\nRaw counts:')
    print(classified_df.astype(int).to_string())
    chi2, p, dof, expected = chi2_contingency(classified_df.values)
    print(f'\nchi2 = {chi2:.2f}  p = {p:.3e}  dof = {dof}')

    pct = 100 * classified_df / classified_df.sum()
    print('\nPercent stacked:')
    print(pct.round(1).to_string())

    # Save stats + percent table
    stats_path = sub / 'loop_classification.stats.txt'
    with stats_path.open('w') as f:
        f.write(f'chi2 = {chi2:.4f}\np = {p:.6e}\ndof = {dof}\n\n')
        f.write('Raw counts:\n')
        f.write(classified_df.astype(int).to_string())
        f.write('\n\nPercent stacked:\n')
        f.write(pct.round(2).to_string())
        f.write('\n')

    # Stacked bar (cluster x-axis, classification stacks)
    palette = {k: CLASSIFICATION_COLORS[k] for k in pct.index}
    with multi_format_savefig():
        plotting.stacked(
            count_table=pct.transpose(),
            measure='loop proportion (%)',
            title='', out_dir=str(sub),
            out_name='loop_classification',
            palette=palette,
        )
    print(f'  Outputs: {sub}/loop_classification.{{png,pdf,svg,jpg}} + loop_classification.stats.txt')


# ============================================================
# 4.3 ChIP signal at anchors (raw RPKM, no RAD21 normalization)
# ============================================================
def _bw_anchor_mean(bw, chrom: str, start: int, end: int) -> float:
    """Return mean BigWig signal in [chrom:start-end], NaN if region is missing."""
    try:
        vals = bw.values(chrom, int(start), int(end))
    except (RuntimeError, ValueError):
        return np.nan
    if vals is None or len(vals) == 0:
        return np.nan
    arr = np.asarray(vals, dtype=float)
    if np.all(np.isnan(arr)):
        return np.nan
    return float(np.nanmean(arr))


def run_anchor_chip(group_dict: dict) -> None:
    """Raw mean BigWig signal per anchor per loop, faceted by mark x condition."""
    print('\n=== [4.3] ChIP signal at anchors (raw RPKM) ===')
    sub = figure_subfolder(FIG_BASE / 'ChIP_intersect', 'anchor_ChIP_box')
    stats_path = FIG_BASE / 'ChIP_intersect' / 'anchor_ChIP.stats.txt'

    # Verify all BigWigs exist + non-empty
    missing = [name for name, p in BIGWIG_DICT.items()
               if not p.exists() or p.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(f'BigWigs missing/empty: {missing}')

    # Build all_groups: every loop with its cluster label
    all_groups = pd.concat(
        [df.assign(group=g) for g, df in group_dict.items()],
        ignore_index=True,
    )
    print(f'  Querying {len(all_groups)} loops x {len(BIGWIG_DICT)} BigWigs')

    # For each BigWig: dedupe each anchor side, query pyBigWig, merge back, average
    # the two anchors per loop. Dedup matters — hub anchors participate in many
    # loops and would otherwise be over-weighted in per-cluster means.
    for name, bw_path in BIGWIG_DICT.items():
        print(f'    {name}: {bw_path.name}')
        bw = pyBigWig.open(str(bw_path))
        for anchor_cols in [['chr1', 'x1', 'x2'], ['chr2', 'y1', 'y2']]:
            anchors = (all_groups[anchor_cols]
                       .drop_duplicates(subset=anchor_cols)
                       .reset_index(drop=True))
            sig_col = f'{name}_{anchor_cols[0]}'
            anchors[sig_col] = [
                _bw_anchor_mean(bw, row[anchor_cols[0]],
                                row[anchor_cols[1]], row[anchor_cols[2]])
                for _, row in anchors.iterrows()
            ]
            all_groups = all_groups.merge(anchors, on=anchor_cols, how='left')
        all_groups[name] = all_groups[[f'{name}_chr1', f'{name}_chr2']].mean(axis=1)
        all_groups = all_groups.drop(columns=[f'{name}_chr1', f'{name}_chr2'])
        bw.close()

    # Long-form for plotting
    all_groups = all_groups[COORD_COLS + ['group'] + list(BIGWIG_DICT)]
    melted = all_groups.melt(
        id_vars=COORD_COLS + ['group'],
        value_vars=list(BIGWIG_DICT),
        var_name='ChIP', value_name='signal',
    ).dropna(subset=['signal']).reset_index(drop=True)

    # Stats: per mark+condition, 6-cluster Kruskal-Wallis + pairwise Wilcoxon
    all_stats = pd.DataFrame()
    for chip in ANCHOR_CHIP_PANEL_ORDER:
        names, lists = [], []
        for g in CLUSTER_ORDER:
            sel = melted.loc[(melted['ChIP'] == chip) & (melted['group'] == g), 'signal'].tolist()
            if not sel:
                continue
            names.append(f'{chip}_{g}')
            lists.append(sel)
        if len(lists) < 2:
            continue
        s = statistics_functions.kruskal_wilcoxon(data_names=names, data_list=lists)
        all_stats = pd.concat([all_stats, s], ignore_index=True)
    all_stats.to_csv(stats_path, sep='\t', index=False)

    # Filter empty subplot panels (defensive — plot_defaults_b can crash on
    # zero-tick axes if a panel is fully empty).
    nonempty_chips = [c for c in ANCHOR_CHIP_PANEL_ORDER
                      if (melted['ChIP'] == c).any()]

    # Faceted boxplot: 4 marks x 2 conditions = 8 panels, ncols=2 → 4x2 layout
    with multi_format_savefig():
        plotting.box(
            melted_df=melted, xcol='group', ycol='signal',
            out_dir=str(sub), measure='mean anchor ChIP signal (RPKM)',
            title='', out_name='anchor_ChIP_box',
            subplots=True, subplot_col='ChIP', order=nonempty_chips,
            ncols=2, sharey=False, Y_range=None,
        )
    print(f'  Outputs: {sub}/anchor_ChIP_box.{{png,pdf,svg,jpg}} + {stats_path.name}')


# ============================================================
# 4.4 ChromHMM anchor vs span enrichment (KEY — Fig 2f equivalent)
# ============================================================
def _write_per_cluster_bed(group_dict: dict, kind: str, out_dir: Path) -> List[str]:
    """Write per-cluster {anchor,span} BEDs. kind in {'anchor','span'}.

    anchor BED = pool of both anchors (chr1,x1,x2) + (chr2,y1,y2).
    span BED   = full loop extent (chr1, x1, y2) — REQUIRES cis loops.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    filenames = []
    for group, df in group_dict.items():
        # cis-only assertion — protects span coords from chr1 != chr2 loops
        if (df['chr1'] != df['chr2']).any():
            n_trans = (df['chr1'] != df['chr2']).sum()
            raise ValueError(
                f'{group} has {n_trans} trans loops (chr1 != chr2). '
                f'Span BED logic only valid for cis loops.'
            )
        if kind == 'anchor':
            out_df = pd.concat([
                df[['chr1', 'x1', 'x2']],
                df[['chr2', 'y1', 'y2']].rename(columns={'chr2': 'chr1', 'y1': 'x1', 'y2': 'x2'}),
            ], ignore_index=True)
        elif kind == 'span':
            out_df = df[['chr1', 'x1', 'y2']].copy()
        else:
            raise ValueError(f'Unknown kind: {kind}')
        bed_path = out_dir / f'{group}.bed'
        out_df.to_csv(bed_path, sep='\t', header=False, index=False)
        filenames.append(f'{group}.bed')
    return filenames


def _run_overlap_enrichment(coord_kind: str, filenames: list) -> Path:
    """Run ChromHMM OverlapEnrichment for a coord_kind ∈ {'anchor','span'}.

    Returns the output enrichment .txt path.
    """
    coordlist = CHROMHMM_DIR / f'coordlistfile_{coord_kind}.txt'
    coordlist.write_text('\n'.join(filenames) + '\n')

    input_dir = CHROMHMM_DIR / f'{coord_kind}_input'
    output_prefix = CHROMHMM_DIR / coord_kind
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
    print(f'  Running: {" ".join(cmd)}')
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print('  STDOUT:', proc.stdout)
        print('  STDERR:', proc.stderr)
        raise RuntimeError(f'OverlapEnrichment failed for {coord_kind} (exit {proc.returncode})')
    out_txt = CHROMHMM_DIR / f'{coord_kind}.txt'
    if not out_txt.exists() or out_txt.stat().st_size == 0:
        raise RuntimeError(f'OverlapEnrichment produced no output: {out_txt}')
    return out_txt


def run_chromhmm_anchor_vs_span(group_dict: dict) -> None:
    """Popay Fig 2f-equivalent: ChromHMM state enrichment at anchors vs spans.

    Tests the central biological question: is Polycomb at anchors (sensitivity
    model) or across loop body (extrusion impediment model)?
    """
    print('\n=== [4.4] ChromHMM anchor vs span enrichment (KEY) ===')
    CHROMHMM_DIR.mkdir(parents=True, exist_ok=True)

    # Verify segment + rename inputs
    for f in [SEGMENT_BED, RENAME_FILE]:
        if not f.exists():
            raise FileNotFoundError(f'Phase 2 artifact missing: {f}')

    # 1. Write per-cluster BEDs
    anchor_files = _write_per_cluster_bed(group_dict, 'anchor', CHROMHMM_DIR / 'anchor_input')
    span_files   = _write_per_cluster_bed(group_dict, 'span',   CHROMHMM_DIR / 'span_input')

    # 2. Run OverlapEnrichment for anchor + span
    anchor_txt = _run_overlap_enrichment('anchor', anchor_files)
    span_txt   = _run_overlap_enrichment('span',   span_files)

    # 3. Heatmaps via Popay's chromHMM_heatmap.heatmap_plot
    #    Output lands at {input_dir}/{stem}.{png,pdf} = cluster/bap1_late/chromHMM/{anchor,span}.{png,pdf}
    #    Our patch adds .svg + .jpg for the .png save (the .pdf save passes through).
    with multi_format_savefig():
        heatmap_plot(path=str(anchor_txt), normalize=False)
        heatmap_plot(path=str(span_txt),   normalize=False)

    print(f'  Enrichment tables: {anchor_txt}, {span_txt}')
    print(f'  Heatmaps: {CHROMHMM_DIR}/{{anchor,span}}.{{png,pdf,svg,jpg}}')

    # Print a quick Polycomb-row summary to log for biological QC
    for label, path in [('anchor', anchor_txt), ('span', span_txt)]:
        df = pd.read_csv(path, sep='\t', header=0, index_col=0)
        if 'Genome %' in df.columns:
            df = df.drop(columns=['Genome %'])
        if 'Polycomb' in df.index:
            print(f'\n  Polycomb fold-enrichment ({label}):')
            print('    ' + df.loc['Polycomb'].round(2).to_string().replace('\n', '\n    '))


# ============================================================
# 4.5 ChromHMM proportions
# ============================================================
def run_chromhmm_proportions(bedpe_df: pd.DataFrame, group_order: list) -> None:
    """Per-loop most-prominent state by overlap length → proportion grouped bar.

    Adapts Popay cell-25 with the E-prefix sort fix (str.replace 'U'->'E' regex)
    and a dynamic palette built from our underscored state names.
    """
    print('\n=== [4.5] ChromHMM proportions per cluster ===')
    sub = figure_subfolder(FIG_BASE, 'chromHMM_anchor')

    segment_bed = pd.read_csv(SEGMENT_BED, sep='\t', header=None,
                              names=['chr', 'start', 'end', 'state'])
    rename_df = load_state_rename()  # cols: short, state, sort_col

    bedpe_cp = bedpe_df.copy()
    bedpe_cp['loop'] = (bedpe_cp['chr1'] + '-' + bedpe_cp['x1'].astype(str)
                        + '-' + bedpe_cp['y2'].astype(str))

    # Pool both anchors as 1 long table, keyed by loop + GROUP
    anchor1 = bedpe_cp[['chr1', 'x1', 'x2', 'loop', 'GROUP']].rename(
        columns={'chr1': 'chr', 'x1': 'start', 'x2': 'end'})
    anchor2 = bedpe_cp[['chr2', 'y1', 'y2', 'loop', 'GROUP']].rename(
        columns={'chr2': 'chr', 'y1': 'start', 'y2': 'end'})
    cluster_df = pd.concat([anchor1, anchor2], ignore_index=True)

    overlap_df = bf.overlap(
        df1=segment_bed, df2=cluster_df,
        cols1=['chr', 'start', 'end'], cols2=['chr', 'start', 'end'],
        how='inner', return_overlap=True,
    )
    overlap_df['overlap_length'] = overlap_df['overlap_end'] - overlap_df['overlap_start']

    # Most-prominent state per anchor-instance (max overlap)
    max_overlap = overlap_df.sort_values(by='overlap_length', ascending=False).drop_duplicates(
        subset=['loop_'], ignore_index=True)

    # Map E-id → biological name; exclude Quiescent (E1, 92% of genome)
    excluded_short = {rename_df.loc[rename_df['state'] == 'Quiescent', 'short'].values[0]}

    prop = pd.DataFrame()
    for group in group_order:
        subset = max_overlap[max_overlap['GROUP_'] == group]
        for short in rename_df['short'].unique():
            if short in excluded_short:
                continue
            biol = rename_df.loc[rename_df['short'] == short, 'state'].values[0]
            prop.loc[biol, group] = (subset['state'] == short).sum()
        prop[group] = 100 * prop[group] / prop[group].sum()

    print('\nPercent state per cluster (Quiescent excluded):')
    print(prop.round(2).to_string())

    # Grouped bar — use plotting.bar (Popay cell-25's active code path)
    melted = prop.transpose().reset_index().rename(columns={'index': 'cluster'})
    melted = melted.melt(id_vars='cluster', var_name='state', value_name='proportion')

    palette = {state: STATE_COLORS.get(state, '#000000') for state in prop.index}
    with multi_format_savefig():
        plotting.bar(
            melted_df=melted, xcol='cluster', ycol='proportion',
            out_dir=str(sub), out_name='chromHMM_anchor',
            measure='proportion of loop anchors (%)', title='',
            hue_col='state', palette=palette,
        )

    # Also save the proportion table
    prop.to_csv(sub / 'chromHMM_anchor.proportions.tsv', sep='\t', index=True)
    print(f'  Outputs: {sub}/chromHMM_anchor.{{png,pdf,svg,jpg}} + chromHMM_anchor.proportions.tsv')


# ============================================================
# 4.6 Gene annotation
# ============================================================
def run_gene_annotation(group_dict: dict) -> None:
    """bedtools-based per-cluster promoter overlap. FPKM-free."""
    print('\n=== [4.6] Gene annotation per cluster ===')
    if shutil.which('bedtools') is None:
        raise RuntimeError(
            'bedtools not on PATH. Install via brew or conda: '
            '`conda install -n cluster -c bioconda bedtools`'
        )

    sub = FIG_BASE / 'annotation'
    sub.mkdir(parents=True, exist_ok=True)
    temp_dir = sub / '_temp'
    temp_dir.mkdir(parents=True, exist_ok=True)

    annotation = bedpe_analysis.bedtools_annotation(
        out_dir=str(sub),
        bedpe_dict=group_dict,
        FPKM_df=None,
        temp_dir=str(temp_dir),
    )

    n_total = 0
    for group, df in annotation.items():
        out_path = sub / f'{group}_annotation.txt'
        df.to_csv(out_path, sep='\t', index=False)
        # Count unique non-blank gene names per cluster
        for col in df.columns:
            if col.endswith('_gene_name'):
                pass  # gene_name col exists per anchor
        all_gene_cols = [c for c in df.columns if c.endswith('_gene_name')]
        if all_gene_cols:
            unique_genes = pd.unique(df[all_gene_cols].values.ravel())
            unique_genes = [g for g in unique_genes if pd.notna(g) and g != '']
            print(f'  {group}: {len(df)} loops, {len(unique_genes)} unique genes')
            n_total += len(unique_genes)
    print(f'  Outputs: {sub}/{{clust1..clust6}}_annotation.txt')


# ============================================================
# 4.7 DiffBind relationship
# ============================================================
def _load_diffbind(path: Path, fdr_cutoff: float, fc_cutoff: float) -> pd.DataFrame:
    """Load DiffBind file, rename our schema → Popay schema, tag sig."""
    df = pd.read_csv(path, sep='\t')
    # Rename Peak_* → Chr/Start/End so cell-17/19 logic works unchanged
    df = df.rename(columns={'Peak_Chr': 'Chr', 'Peak_Start': 'Start', 'Peak_End': 'End'})
    df['sig'] = 'non-sig.'
    df.loc[(df['FDR'] < fdr_cutoff) & (df['Fold'] < -fc_cutoff), 'sig'] = 'decreased'
    df.loc[(df['FDR'] < fdr_cutoff) & (df['Fold'] >  fc_cutoff), 'sig'] = 'increased'
    return df


def _run_diffbind_for(mark: str, db: pd.DataFrame, group_dict: dict,
                      fc_cutoff: float, sub: Path) -> None:
    """One mark × one cutoff: proportions stacked + FC boxplot."""
    fc_tag = f'fc{fc_cutoff:.1f}'.replace('.', 'p')
    prop_dir = figure_subfolder(sub, f'differential_binding_{mark}_{fc_tag}')
    box_dir  = figure_subfolder(sub, f'ChIP_FC_{mark}_{fc_tag}')
    stats_path = sub / f'diffbind_stats_{mark}_{fc_tag}.txt'

    # ---- 4.7a proportions (cell-17 logic, with anchor2 typo fix) ----
    sig_summary = pd.DataFrame()
    for group, df in group_dict.items():
        out_df = df.copy()
        # overlap each anchor with diffbind peaks (suffix per chr)
        out_df = bf.overlap(out_df, db[['Chr', 'Start', 'End', 'Fold', 'sig']],
                            cols1=['chr1', 'x1', 'x2'], cols2=['Chr', 'Start', 'End'],
                            suffixes=('', '_chr1'))
        out_df = bf.overlap(out_df, db[['Chr', 'Start', 'End', 'Fold', 'sig']],
                            cols1=['chr2', 'y1', 'y2'], cols2=['Chr', 'Start', 'End'],
                            suffixes=('', '_chr2'))
        out_df = out_df.drop(columns=['Chr_chr1', 'Start_chr1', 'End_chr1',
                                      'Chr_chr2', 'Start_chr2', 'End_chr2'])
        # CRITICAL FIX: Popay cell-17 has typo value_vars=['sig_chr1','sig_chr1']
        # — both chr1, anchor2 silently ignored. Fixed here.
        melt_df = out_df.melt(id_vars=COORD_COLS,
                              value_vars=['sig_chr1', 'sig_chr2'],
                              var_name='anchor', value_name='sig')
        # Priority: increased > decreased > non-sig. > NaN (no peak) — pick the
        # most informative tag per loop across both anchors.
        sorted_melt = sort_by_strings(
            df=melt_df, sort_col='sig',
            order=['no peak', 'non-sig.', 'decreased', 'increased'],
        )
        sorted_melt = sorted_melt.dropna(subset=['sig']).reset_index(drop=True)
        # For each loop, keep the LAST (highest priority) tag — sort_by_strings
        # already encoded order, so highest int-encoded value wins. Re-decode
        # done by sort_by_strings; we just keep_last per loop.
        deduped = sorted_melt.drop_duplicates(subset=COORD_COLS, keep='last',
                                              ignore_index=True)
        # Loops where neither anchor overlapped any peak don't appear in deduped
        n_loops = df.shape[0]
        n_with_peak = deduped.shape[0]
        sig_summary.loc['no peak', group]    = n_loops - n_with_peak
        sig_summary.loc['non-sig.', group]   = (deduped['sig'] == 'non-sig.').sum()
        sig_summary.loc['decreased', group]  = (deduped['sig'] == 'decreased').sum()
        sig_summary.loc['increased', group]  = (deduped['sig'] == 'increased').sum()

    chi2, p, dof, _ = chi2_contingency(sig_summary.values)
    pct = 100 * sig_summary / sig_summary.sum()

    palette_prop = {k: DIFFBIND_COLORS[k] for k in pct.index}
    with multi_format_savefig():
        plotting.stacked(
            count_table=pct.transpose(),
            measure='proportion of loops (%)',
            title='', out_dir=str(prop_dir),
            out_name=f'differential_binding_{mark}_{fc_tag}',
            palette=palette_prop,
        )

    # ---- 4.7b mean fold-change boxplot (cell-19 logic) ----
    fc_summary = pd.DataFrame()
    data_names, data_list = [], []
    for group, df in group_dict.items():
        out_df = df.copy()
        chr1_df = bf.overlap(out_df, db[['Chr', 'Start', 'End', 'Fold', 'sig']],
                             cols1=['chr1', 'x1', 'x2'], cols2=['Chr', 'Start', 'End'],
                             suffixes=('', '_chr1'))
        chr1_df = (chr1_df[COORD_COLS + ['Fold_chr1']].dropna(subset=['Fold_chr1'])
                   .groupby(by=COORD_COLS, as_index=False).max())

        chr2_df = bf.overlap(out_df, db[['Chr', 'Start', 'End', 'Fold', 'sig']],
                             cols1=['chr2', 'y1', 'y2'], cols2=['Chr', 'Start', 'End'],
                             suffixes=('', '_chr2'))
        chr2_df = (chr2_df[COORD_COLS + ['Fold_chr2']].dropna(subset=['Fold_chr2'])
                   .groupby(by=COORD_COLS, as_index=False).max())

        out_df = out_df.merge(chr1_df, on=COORD_COLS, how='left')
        out_df = out_df.merge(chr2_df, on=COORD_COLS, how='left')
        out_df['mean_fold'] = out_df[['Fold_chr1', 'Fold_chr2']].mean(axis=1)
        out_df = out_df.dropna(subset=['mean_fold'])
        out_df['group'] = group
        fc_summary = pd.concat([fc_summary, out_df], ignore_index=True)
        if not out_df.empty:
            data_names.append(group)
            data_list.append(out_df['mean_fold'].tolist())

    kw_stats = statistics_functions.kruskal_wilcoxon(data_names=data_names,
                                                      data_list=data_list)

    with multi_format_savefig():
        plotting.box(
            melted_df=fc_summary, xcol='group', ycol='mean_fold',
            out_dir=str(box_dir),
            measure=f'mean log2 fold-change ({mark})',
            title='', out_name=f'ChIP_FC_{mark}_{fc_tag}',
            Y_range=None,
        )

    # Save stats
    with stats_path.open('w') as f:
        f.write(f'mark = {mark}\nFDR < 0.05, |Fold| > {fc_cutoff}\n\n')
        f.write(f'Proportions chi2 = {chi2:.4f}  p = {p:.6e}  dof = {dof}\n\n')
        f.write('Counts:\n' + sig_summary.astype(int).to_string() + '\n\n')
        f.write('Percent stacked:\n' + pct.round(2).to_string() + '\n\n')
        f.write('Fold-change stats (Kruskal + Wilcoxon):\n')
        f.write(kw_stats.to_string(index=False) + '\n')

    print(f'    {mark} fc>{fc_cutoff}: chi2={chi2:.2f}, p={p:.2e}; '
          f'fc-data n={len(fc_summary)} loops')


def run_diffbind(group_dict: dict, fc_cutoffs=(0.3, 0.0), fdr_cutoff=0.05) -> None:
    print('\n=== [4.7] DiffBind relationship (3 marks x 2 thresholds) ===')
    sub = FIG_BASE / 'ChIP_intersect'
    sub.mkdir(parents=True, exist_ok=True)

    for mark, path in DIFFBIND_FILES.items():
        if not path.exists():
            raise FileNotFoundError(f'DiffBind file missing: {path}')
        for fc in fc_cutoffs:
            db = _load_diffbind(path, fdr_cutoff, fc)
            print(f'  Mark={mark}  FDR<{fdr_cutoff}  |Fold|>{fc}: '
                  f'{(db["sig"] != "non-sig.").sum()} sig peaks of {len(db)}')
            _run_diffbind_for(mark, db, group_dict, fc, sub)


# ============================================================
# 4.8 Cluster x differential status crosstab (NEW)
# ============================================================
def run_cluster_diff_crosstab(bedpe_df: pd.DataFrame, group_order: list) -> None:
    print('\n=== [4.8] Cluster x edgeR direction crosstab (NEW) ===')
    sub = figure_subfolder(FIG_BASE, 'cluster_differential_status')

    meta = pd.read_csv(METADATA_FILE, sep='\t', usecols=COORD_COLS + ['direction'])
    joined = bedpe_df.merge(meta, on=COORD_COLS, how='inner')
    if len(joined) != len(bedpe_df):
        warnings.warn(f'Metadata join lost {len(bedpe_df) - len(joined)} rows')

    crosstab = pd.crosstab(joined['GROUP'], joined['direction'])
    # Reorder rows by cluster order, columns to ensure direction order
    crosstab = crosstab.reindex(index=[g for g in group_order if g in crosstab.index],
                                columns=['down_in_mutant', 'unchanged', 'up_in_mutant'],
                                fill_value=0)

    chi2, p, dof, _ = chi2_contingency(crosstab.values)
    pct = 100 * crosstab.div(crosstab.sum(axis=1), axis=0)

    print('\nCounts:')
    print(crosstab.to_string())
    print('\nRow %:')
    print(pct.round(1).to_string())
    print(f'\nchi2 = {chi2:.2f}  p = {p:.3e}  dof = {dof}')

    # Save stats
    stats_path = sub.parent / 'cluster_differential_status.stats.txt'
    with stats_path.open('w') as f:
        f.write(f'chi2 = {chi2:.4f}\np = {p:.6e}\ndof = {dof}\n\n')
        f.write('Counts:\n' + crosstab.to_string() + '\n\n')
        f.write('Row %:\n' + pct.round(2).to_string() + '\n')

    # Stacked bar: cluster on x-axis (rows of crosstab → x), direction stacks.
    # plotting.stacked melts via reset_index(drop=False).rename({'index':'xcol'})
    # — only works when the index is unnamed. pd.crosstab sets index.name='GROUP',
    # so we clear it. Same trick avoids KeyError('xcol') in the upstream melt.
    pct.index.name = None
    pct.columns.name = None
    palette = {k: DIRECTION_COLORS[k] for k in pct.columns}
    with multi_format_savefig():
        plotting.stacked(
            count_table=pct,  # index=cluster (x-axis), columns=direction (stacks)
            measure='proportion of loops (%)',
            title='', out_dir=str(sub),
            out_name='cluster_differential_status',
            palette=palette,
        )
    print(f'  Outputs: {sub}/cluster_differential_status.{{png,pdf,svg,jpg}} + '
          f'{stats_path.name}')


# ============================================================
# Main dispatch
# ============================================================
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--analyses', default='all',
                   help='Comma-separated subset to run, or "all" (default). '
                        f'Valid: {", ".join(sorted(VALID_ANALYSES))}.')
    args = p.parse_args()

    selected = parse_analyses(args.analyses)
    print(f'Selected analyses: {", ".join(selected)}')

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    FIG_BASE.mkdir(parents=True, exist_ok=True)
    (FIG_BASE / 'ChIP_intersect').mkdir(parents=True, exist_ok=True)
    CHROMHMM_DIR.mkdir(parents=True, exist_ok=True)

    bedpe_df, group_dict, group_order = load_clusters()

    # Order: KEY result first, then cheap, then expensive.
    if '4.4' in selected:
        run_chromhmm_anchor_vs_span(group_dict)
    if '4.1' in selected:
        run_loop_size(group_dict)
    if '4.8' in selected:
        run_cluster_diff_crosstab(bedpe_df, group_order)
    if '4.5' in selected:
        run_chromhmm_proportions(bedpe_df, group_order)
    if '4.6' in selected:
        run_gene_annotation(group_dict)
    if '4.2' in selected:
        run_loop_classification(group_dict)
    if '4.7' in selected:
        run_diffbind(group_dict)
    if '4.3' in selected:
        run_anchor_chip(group_dict)

    print('\n=== All requested Phase 4 analyses complete ===')


if __name__ == '__main__':
    main()
