#!/usr/bin/env python3
# cluster/scripts/07_cooltools_pileup.py
"""Phase 6: cooltools off-diagonal pileup over per-cluster BEDPEs.

Per-cluster pileup of log2(obs/exp) Hi-C contact frequency, computed at 10kb
resolution over a +/-500kb window around each loop midpoint, separately for
ctrl_merged and mut_merged mcools. Cluster ordering matches Phase 3 v2 output.

Reads:
  cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt   (Phase 3 v2)
  ${MCOOL_CTRL}, ${MCOOL_MUT}                                 (CLI args)

Writes:
  cluster/bap1_late/cooltools/obs_exp_contacts/obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}

NOTE: cooltools_called.make_arms() requires outbound HTTPS to UCSC (bioframe
fetches mm10 chromsizes + centromeres). If the compute node has no network
access, pre-fetch on a login node and patch make_arms to load from cache.
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')           # disable GUI; headless on HPC compute nodes

import bioframe
import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT   = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR))            # cooltools_called
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))   # multi_format_output

import cooltools_called                         # noqa: E402
from cooltools_called import mcool_pileup        # noqa: E402
from multi_format_output import multi_format_savefig, figure_subfolder  # noqa: E402


def make_arms_robust(genome: str = 'mm10') -> pd.DataFrame:
    """Robust replacement for cooltools_called.make_arms.

    bioframe>=0.5 deprecated UCSC online fetch and ships a 'local' provider
    that needs bundled centromere data. mm10 centromere data is not bundled
    in some installations -> fetch_centromeres raises ValueError.

    Fallback: whole-chromosome arms (one "arm" per chromosome from 0 to
    chromsize). cooltools.expected_cis with single-arm-per-chrom view_df
    computes whole-chromosome expected, which is fine for our +/-500kb
    pileup use case (centromere-aware splits matter only for trans or
    centromere-spanning cis contacts).
    """
    chromsizes = bioframe.fetch_chromsizes(genome)
    try:
        cens = bioframe.fetch_centromeres(genome)
        arms = bioframe.make_chromarms(chromsizes, cens)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        print(f'  make_arms: fetch_centromeres failed '
              f'({type(exc).__name__}: {exc}); '
              f'falling back to whole-chromosome arms.', flush=True)
        arms = pd.DataFrame({
            'chrom': list(chromsizes.index),
            'start': 0,
            'end':   [int(v) for v in chromsizes.values],
            'name':  list(chromsizes.index),
        })
    return arms


# Monkey-patch cooltools_called.make_arms BEFORE mcool_pileup is invoked so that
# the call inside mcool_pileup uses the robust version. Module-level function
# lookup in Python is dynamic, so replacing the attribute in the module dict
# redirects subsequent calls.
cooltools_called.make_arms = make_arms_robust

CLUSTER_FILE  = CLUSTER_DIR / 'bap1_late/cluster3/k-6/data/combined-clusters.txt'
OUT_BASE      = CLUSTER_DIR / 'bap1_late/cooltools'
CLUSTER_ORDER = ['clust1', 'clust2', 'clust3', 'clust4', 'clust5', 'clust6']


def load_bedpe_dict(path: Path) -> dict:
    """Read combined-clusters.txt -> dict[cluster_name -> DataFrame].
    cooltools.pileup requires the cooler/cooltools 6-column convention
    (chrom1/start1/end1/chrom2/start2/end2), NOT the project's chr1/x1/x2."""
    df = pd.read_csv(path, sep='\t')
    bedpe_dict = {}
    for c in CLUSTER_ORDER:
        sub = df[df['GROUP'] == c].copy()
        if len(sub) == 0:
            print(f'  WARNING: cluster {c} has 0 loops -- skipping')
            continue
        sub = sub.rename(columns={
            'chr1': 'chrom1', 'x1': 'start1', 'x2': 'end1',
            'chr2': 'chrom2', 'y1': 'start2', 'y2': 'end2',
        })
        bedpe_dict[c] = sub[
            ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']
        ].reset_index(drop=True)
    return bedpe_dict


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--mcool-ctrl', required=True, help='Path to ctrl mcool')
    ap.add_argument('--mcool-mut',  required=True, help='Path to mut mcool')
    ap.add_argument('--resolution', type=int, default=10000,
                    help='Pileup binsize in bp (default 10000 = 10kb, Popay default)')
    ap.add_argument('--flank',      type=int, default=500000,
                    help='Window flank in bp (default 500000 = +/-500kb, Popay default)')
    ap.add_argument('--out-dir',    default=str(OUT_BASE), help='Output base directory')
    args = ap.parse_args()

    if not CLUSTER_FILE.exists():
        raise FileNotFoundError(f'Cluster file not found: {CLUSTER_FILE}')
    for label, path in [('ctrl', args.mcool_ctrl), ('mut', args.mcool_mut)]:
        if not Path(path).exists():
            raise FileNotFoundError(f'{label} mcool not found: {path}')

    print(f'[1/2] Loading clusters from {CLUSTER_FILE}...', flush=True)
    bedpe_dict = load_bedpe_dict(CLUSTER_FILE)
    for c, df in bedpe_dict.items():
        print(f'  {c}: {len(df)} loops', flush=True)

    mcool_dict = {'ctrl': args.mcool_ctrl, 'mut': args.mcool_mut}
    out_base = Path(args.out_dir)
    out_base.mkdir(parents=True, exist_ok=True)
    sub = figure_subfolder(out_base, 'obs_exp_contacts')

    print(f'\n[2/2] cooltools pileup at {args.resolution}bp / flank +/-{args.flank}bp', flush=True)
    print(f'  out_dir:    {sub}')
    print(f'  mcool_ctrl: {args.mcool_ctrl}')
    print(f'  mcool_mut:  {args.mcool_mut}')
    print('  NOTE: bioframe.fetch_chromsizes("mm10") + fetch_centromeres requires '
          'outbound HTTPS to UCSC. If this fails on a compute node, pre-fetch '
          'arms on the login node and patch make_arms to load from cache.', flush=True)

    with multi_format_savefig():
        mcool_pileup(
            mcool_dict=mcool_dict,
            bedpe_dict=bedpe_dict,
            out_dir=str(sub),
            out_name='obs_exp_contacts',
            flank=args.flank,
            resolution=args.resolution,
            over_expected=True,
            split_diagonal=False,    # latent bug if True; see Plan-p1 corrections
            v_range=[-1, 2],         # Popay notebook value (log2obs/exp range)
            genome='mm10',           # Phase 0 fix; redundant with default
        )

    print('\nPhase 6 complete. Outputs:', flush=True)
    print(f'  {sub}/obs_exp_contacts_ctrl.{{png,pdf,svg,jpg}}')
    print(f'  {sub}/obs_exp_contacts_mut.{{png,pdf,svg,jpg}}')


if __name__ == '__main__':
    main()
