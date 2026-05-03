#!/usr/bin/env python3
# cluster/scripts/09_oriented_anchor_metagene.py
"""Oriented anchor metagene -- test K27me3 exterior/interior asymmetry.

Builds per-cluster oriented BED6 files (anchor1 -> strand+, anchor2 -> strand-)
then runs deepTools bed_pileup with the same 8 BigWigs as Phase 5.
deepTools natively flips the signal array for strand '-' regions so that the
left side of the plot (-b / 'upstream') always = EXTERIOR (away from loop body)
and right side (+a / 'downstream') always = INTERIOR (toward loop body).

Tests: clust6 (lost loops, Polycomb at anchors) should show asymmetric K27me3
enrichment on the exterior side only (extrusion impediment model).
Contrast: clust5 (gained Polycomb loops) expects symmetric K27me3 on both sides.

Reads:
  cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt
  /Users/zakiralibhai/sdsc/bigwigs/                           (8 BigWigs)
  tads/mm10-blacklist.v2.bed

Writes:
  cluster/bap1_late/figures/deeptools_input/clust{1..6}_oriented_anchors.bed
  cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values{,.png,.pdf,.svg,.jpg}
"""
import shutil
import sys
from pathlib import Path
from typing import Dict

import matplotlib
matplotlib.use('Agg')

import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR / 'modules'))  # deepTools_pipeline, deeptools_plotting
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))

from deepTools_pipeline import bed_pileup        # noqa: E402
from multi_format_output import multi_format_savefig, figure_subfolder  # noqa: E402
from pipeline_config import (                     # noqa: E402
    get_cluster_order, build_bigwig_dict, build_vmax_groups, build_color_dict,
)

import os as _os
_out = _os.environ.get('CLUSTER_OUT_DIR', 'outputs/bap1_late')
_k = _os.environ.get('CLUSTER_K', '6')
CLUSTER_FILE  = Path(_os.environ.get('CLUSTER_COMBINED',
    str(CLUSTER_DIR / '{}/cluster3/k-{}/data/combined-clusters.txt'.format(_out, _k))))
BIGWIG_BASE   = Path(_os.environ.get('CLUSTER_BIGWIG_DIR', '/Users/zakiralibhai/sdsc/bigwigs'))
BLACKLIST     = REPO_ROOT / 'tads/mm10-blacklist.v2.bed'

ANCHOR_BED_DIR = CLUSTER_DIR / _out / 'figures/deeptools_input'
DEEPTOOLS_DIR  = CLUSTER_DIR / _out / 'figures/deeptools'

CLUSTER_ORDER = get_cluster_order()
BIGWIG_DICT   = build_bigwig_dict(BIGWIG_BASE)
VMAX_GROUPS   = build_vmax_groups(BIGWIG_DICT)
COLOR_DICT    = build_color_dict(BIGWIG_DICT)

XTICKLABELS = ['-5kb (exterior)', 'anchor', '+5kb (interior)']


def build_oriented_anchor_beds(cluster_file, out_dir):
    # type: (Path, Path) -> Dict[str, str]
    """Build per-cluster oriented BED6 files.

    Anchor1 (left, x1 < y1) gets strand '+': deepTools downstream = rightward = loop interior.
    Anchor2 (right) gets strand '-': deepTools flips so downstream = leftward = loop interior.
    Result: plot left = exterior, plot right = interior for both anchor types.

    Dedup on (chrom, start, end, strand) so hub anchors appearing as both
    anchor1 and anchor2 within a cluster keep both orientation entries."""
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(cluster_file, sep='\t')
    bed_dict = {}
    for cluster in CLUSTER_ORDER:
        sub = df[df['GROUP'] == cluster]
        if len(sub) == 0:
            print(f'  WARNING: cluster {cluster} has 0 loops -- skipping')
            continue
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
        bed_path = out_dir / f'{cluster}_oriented_anchors.bed'
        oriented[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            bed_path, sep='\t', header=False, index=False, lineterminator='\n')

        n_unique_positions = len(
            pd.concat([a1[['chrom', 'start', 'end']], a2[['chrom', 'start', 'end']]])
              .drop_duplicates())
        n_hub = len(oriented) - n_unique_positions
        print(f'  {cluster}: {len(sub):>6} loops -> {len(oriented):>6} oriented anchors '
              f'({n_hub} hub anchors with dual orientation) -> {bed_path.name}')
        bed_dict[cluster] = str(bed_path)
    return bed_dict


def main():
    # type: () -> None
    for path, label in [(CLUSTER_FILE, 'Cluster file'),
                        (BIGWIG_BASE,  'BigWig dir'),
                        (BLACKLIST,    'Blacklist BED')]:
        if not path.exists():
            raise FileNotFoundError(f'{label} not found: {path}')
    if shutil.which('computeMatrix') is None:
        raise RuntimeError(
            'computeMatrix not on PATH. Prepend '
            '/opt/homebrew/anaconda3/envs/cluster/bin to PATH '
            '(see run_oriented_metagene.sh).')

    bigwig_dict = {}
    for label, path in BIGWIG_DICT.items():
        if not path.exists():
            raise FileNotFoundError(f'BigWig missing: {path}')
        bigwig_dict[label] = str(path)

    print(f'\n[1/2] Building per-cluster oriented BED6 in {ANCHOR_BED_DIR}...')
    bed_dict = build_oriented_anchor_beds(CLUSTER_FILE, ANCHOR_BED_DIR)

    print(f'\n[2/2] Running deepTools bed_pileup '
          f'({len(bigwig_dict)} BigWigs x {len(bed_dict)} clusters)...')
    print(f'  Orientation: anchor1=strand+, anchor2=strand-')
    print(f'  After deepTools flip: left=-5kb(exterior), right=+5kb(interior)')
    out_dir = figure_subfolder(DEEPTOOLS_DIR, 'oriented_anchors')
    print(f'  out_dir:   {out_dir}')
    print(f'  blacklist: {BLACKLIST}')

    with multi_format_savefig():
        bed_pileup(
            bed_dict=bed_dict,
            bigWig_dict=bigwig_dict,
            out_dir=str(out_dir),
            out_name='oriented_anchors',
            blacklisted_regions=str(BLACKLIST),
            up_down=5000,
            color_dict=COLOR_DICT,
            vmax_groups=VMAX_GROUPS,
            line_measure='mean',
            pileup_type='referencePoint',
            xticklabels=XTICKLABELS,
        )

    print('\nOriented metagene complete. Outputs:')
    print(f'  {out_dir}/oriented_anchors_values.{{png,pdf,svg,jpg}}')


if __name__ == '__main__':
    main()
