#!/usr/bin/env python3
# cluster/scripts/06_deeptools_metagene.py
"""Phase 5: deepTools metagene at loop anchors per cluster.

Builds per-cluster anchor BEDs (both halves of each loop, deduplicated to 3-col)
then runs deepTools `bed_pileup` to produce one combined heatmap with rows = 6
clusters and cols = 8 BigWigs (4 marks x ctrl/mut). vmax_groups pairs ctrl/mut
for each mark to enable fair visual comparison; color_dict assigns mark-specific
colormaps. Output via multi_format_savefig context for .png/.pdf/.svg/.jpg
sibling emission.

Reads:
  cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt   (Phase 3 v2)
  /Users/zakiralibhai/sdsc/bigwigs/                            (8 BigWigs)
  tads/mm10-blacklist.v2.bed                                   (computeMatrix -bl)

Writes:
  cluster/bap1_late/figures/deeptools_input/clust{1..6}_anchors.bed
  cluster/bap1_late/figures/deeptools/histone_anchors/histone_anchors_values{,.png,.pdf,.svg,.jpg}
"""
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR))            # deepTools_pipeline, deeptools_plotting
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))   # multi_format_output

from deepTools_pipeline import bed_pileup        # noqa: E402
from multi_format_output import multi_format_savefig, figure_subfolder  # noqa: E402

# -------- inputs --------
CLUSTER_FILE = CLUSTER_DIR / 'bap1_late/cluster3/k-6/data/combined-clusters.txt'
BIGWIG_BASE  = Path('/Users/zakiralibhai/sdsc/bigwigs')
BLACKLIST    = REPO_ROOT / 'tads/mm10-blacklist.v2.bed'

# -------- outputs --------
ANCHOR_BED_DIR = CLUSTER_DIR / 'bap1_late/figures/deeptools_input'
DEEPTOOLS_DIR  = CLUSTER_DIR / 'bap1_late/figures/deeptools'

# -------- constants --------
CLUSTER_ORDER = ['clust1', 'clust2', 'clust3', 'clust4', 'clust5', 'clust6']

# 4 marks x 2 conditions = 8 BigWigs (BAP1-relevant: PR-DUB substrate K119ub,
# PRC2 product K27me3, K27ac active marker, K27me1 PRC2 intermediate)
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

# Pair ctrl + mut so the 98th-percentile vmax is shared across conditions per mark.
# Without paired vmax, mut signal would render against its own scale and visual
# differences ctrl vs mut would be artifacts of the colormap rather than biology.
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


def build_anchor_beds(cluster_file: Path, out_dir: Path) -> dict:
    """For each cluster, write a deduplicated 3-col BED of both anchor halves.
    Returns dict: cluster_name -> str(absolute path to BED).
    Hub anchors that participate in multiple loops within a cluster appear once."""
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(cluster_file, sep='\t')
    bed_dict = {}
    for cluster in CLUSTER_ORDER:
        sub = df[df['GROUP'] == cluster]
        if len(sub) == 0:
            print(f'  WARNING: cluster {cluster} has 0 loops -- skipping')
            continue
        a1 = sub[['chr1', 'x1', 'x2']].rename(columns={'chr1': 'chrom', 'x1': 'start', 'x2': 'end'})
        a2 = sub[['chr2', 'y1', 'y2']].rename(columns={'chr2': 'chrom', 'y1': 'start', 'y2': 'end'})
        anchors = (
            pd.concat([a1, a2])
              .drop_duplicates(subset=['chrom', 'start', 'end'])
              .sort_values(['chrom', 'start', 'end'])
              .reset_index(drop=True)
        )
        bed_path = out_dir / f'{cluster}_anchors.bed'
        anchors.to_csv(bed_path, sep='\t', header=False, index=False)
        bed_dict[cluster] = str(bed_path)
        print(f'  {cluster}: {len(sub):>6} loops -> {len(anchors):>6} unique anchors -> {bed_path.name}')
    return bed_dict


def main() -> None:
    if not CLUSTER_FILE.exists():
        raise FileNotFoundError(f'Cluster file not found: {CLUSTER_FILE}')
    if not BIGWIG_BASE.exists():
        raise FileNotFoundError(f'BigWig dir not found: {BIGWIG_BASE}')
    if not BLACKLIST.exists():
        raise FileNotFoundError(f'Blacklist BED not found: {BLACKLIST}')
    if shutil.which('computeMatrix') is None:
        raise RuntimeError(
            'computeMatrix not on PATH. Either activate the cluster env, or '
            'prepend /opt/homebrew/anaconda3/envs/cluster/bin to PATH '
            '(see run_phase5.sh).'
        )

    bigwig_dict = {}
    for label, path in BIGWIG_DICT.items():
        if not path.exists():
            raise FileNotFoundError(f'BigWig missing: {path}')
        bigwig_dict[label] = str(path)

    print(f'\n[1/2] Building per-cluster anchor BEDs in {ANCHOR_BED_DIR}...')
    bed_dict = build_anchor_beds(CLUSTER_FILE, ANCHOR_BED_DIR)

    print(f'\n[2/2] Running deepTools bed_pileup '
          f'({len(bigwig_dict)} BigWigs x {len(bed_dict)} clusters)...')
    out_dir = figure_subfolder(DEEPTOOLS_DIR, 'histone_anchors')
    print(f'  out_dir:     {out_dir}')
    print(f'  blacklist:   {BLACKLIST}')
    print(f'  up_down:     +/-5kb (referencePoint=center)')
    print(f'  vmax_groups: pair ctrl/mut per mark')

    with multi_format_savefig():
        bed_pileup(
            bed_dict=bed_dict,
            bigWig_dict=bigwig_dict,
            out_dir=str(out_dir),
            out_name='histone_anchors',
            blacklisted_regions=str(BLACKLIST),
            up_down=5000,
            color_dict=COLOR_DICT,
            vmax_groups=VMAX_GROUPS,
            line_measure='mean',
            pileup_type='referencePoint',
        )

    print('\nPhase 5 complete. Outputs:')
    print(f'  {out_dir}/histone_anchors_values.{{png,pdf,svg,jpg}}')


if __name__ == '__main__':
    main()
