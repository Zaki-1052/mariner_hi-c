# ML/cmpts/scripts/B1_adapt_sniper_mm10.py
# Generate mm10 crop map for SNIPER from ctrl_merged.hic.
# Dumps inter-chromosomal contacts via juicer_tools, builds the full matrix,
# filters sparse and blacklisted bins, and saves cropMap + cropIndices .mat files.
import os
import sys
import argparse
import numpy as np
import pandas as pd
from subprocess import call
from scipy.io import savemat

_script_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.abspath(os.path.join(_script_dir, '..', '..', '..'))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from ML.cmpts.repos.SNIPER.utilities.mm10_config import (
    ODD_CHROMS, EVEN_CHROMS,
    N_ODD_BINS, N_EVEN_BINS,
    ODD_BIN_COUNTS, EVEN_BIN_COUNTS,
    ODD_OFFSETS, EVEN_OFFSETS,
    RESOLUTION,
    JUICER_CONTAINER, JUICER_JAR
)

JUICER_CMD = 'singularity exec --bind /scratch,/expanse {0} java -jar {1}'.format(
    JUICER_CONTAINER, JUICER_JAR
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate mm10 crop map for SNIPER from ctrl_merged.hic'
    )
    parser.add_argument('hic_path', help='Path to ctrl_merged.hic file')
    parser.add_argument('blacklist_path', help='Path to mm10-blacklist.v2.bed')
    parser.add_argument('output_dir', help='Directory for mm10_cropMap.mat and mm10_cropIndices.mat')
    parser.add_argument('--tmp-dir', default='.', help='Directory for juicer_tools .txt intermediates')
    parser.add_argument('--timepoint', default=None, help='Timepoint label (e.g. 250402) appended to output filenames')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing .txt files')
    parser.add_argument('--sizes-file', default=None,
                        help='Path to mm10.chrom.sizes (default: data/mm10.chrom.sizes relative to SNIPER repo)')
    return parser.parse_args()


def dump_interchromosomal_contacts(hic_path, tmp_dir, overwrite):
    """Dump 90 inter-chromosomal contact files via juicer_tools (10 odd x 9 even)."""
    os.makedirs(tmp_dir, exist_ok=True)
    total = len(ODD_CHROMS) * len(EVEN_CHROMS)
    done = 0

    for chrm1 in ODD_CHROMS:
        for chrm2 in EVEN_CHROMS:
            output_path = os.path.join(tmp_dir, 'cropmap_chrm{0}_chrm{1}.txt'.format(chrm1, chrm2))

            if os.path.isfile(output_path) and not overwrite:
                if os.path.getsize(output_path) > 0:
                    done += 1
                    continue

            cmd = '{0} dump observed KR {1} chr{2} chr{3} BP {4} {5}'.format(
                JUICER_CMD, hic_path, chrm1, chrm2, RESOLUTION, output_path
            )
            call([cmd], shell=True)

            if not os.path.isfile(output_path) or os.path.getsize(output_path) == 0:
                print('WARNING: empty dump for chr{0} x chr{1}'.format(chrm1, chrm2))

            done += 1
            if done % 10 == 0:
                print('  Dumped {0}/{1} chromosome pairs'.format(done, total))

    print('  Dump complete: {0}/{1} pairs'.format(done, total))


def build_matrix(tmp_dir, sizes_file):
    """Build the full inter-chromosomal matrix from dump files."""
    from ML.cmpts.repos.SNIPER.utilities.interchromosome_matrix_mm10 import construct
    print('Building inter-chromosomal matrix...')
    M = construct(hic_dir=tmp_dir, prefix='cropmap', sizes_file=sizes_file)
    print('  Matrix shape: {0}'.format(M.shape))
    assert M.shape == (N_ODD_BINS, N_EVEN_BINS), (
        'Expected ({0}, {1}), got {2}'.format(N_ODD_BINS, N_EVEN_BINS, M.shape)
    )
    return M


def compute_coverage_filter(M):
    """Identify sparse bins (bottom 1st percentile of coverage)."""
    row_sums = np.nansum(M, axis=1)
    col_sums = np.nansum(M, axis=0)

    row_thresh = np.percentile(row_sums, 1)
    col_thresh = np.percentile(col_sums, 1)

    sparse_rows = row_sums <= row_thresh
    sparse_cols = col_sums <= col_thresh

    print('  Coverage filter:')
    print('    Row threshold (1st pctile): {0:.4f}'.format(row_thresh))
    print('    Col threshold (1st pctile): {0:.4f}'.format(col_thresh))
    print('    Sparse rows: {0}'.format(sparse_rows.sum()))
    print('    Sparse cols: {0}'.format(sparse_cols.sum()))

    return sparse_rows, sparse_cols


def compute_blacklist_filter(blacklist_path):
    """Identify bins overlapping ENCODE blacklist regions."""
    bl_df = pd.read_csv(blacklist_path, sep='\t', header=None,
                        usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    bl_df = bl_df[bl_df['chrom'].str.match(r'^chr\d+$')]

    bl_row_mask = np.zeros(N_ODD_BINS, dtype=bool)
    bl_col_mask = np.zeros(N_EVEN_BINS, dtype=bool)

    for c_idx, c in enumerate(ODD_CHROMS):
        chrom = 'chr{0}'.format(c)
        chrom_bl = bl_df[bl_df['chrom'] == chrom]
        if chrom_bl.empty:
            continue
        offset = ODD_OFFSETS[c_idx]
        n_bins = ODD_BIN_COUNTS[c_idx]
        bin_starts = np.arange(n_bins) * RESOLUTION
        bin_ends = bin_starts + RESOLUTION
        bl_starts = chrom_bl['start'].values
        bl_ends = chrom_bl['end'].values
        overlaps = (bin_starts[:, np.newaxis] < bl_ends) & (bl_starts < bin_ends[:, np.newaxis])
        bl_row_mask[offset:offset + n_bins] = overlaps.any(axis=1)

    for c_idx, c in enumerate(EVEN_CHROMS):
        chrom = 'chr{0}'.format(c)
        chrom_bl = bl_df[bl_df['chrom'] == chrom]
        if chrom_bl.empty:
            continue
        offset = EVEN_OFFSETS[c_idx]
        n_bins = EVEN_BIN_COUNTS[c_idx]
        bin_starts = np.arange(n_bins) * RESOLUTION
        bin_ends = bin_starts + RESOLUTION
        bl_starts = chrom_bl['start'].values
        bl_ends = chrom_bl['end'].values
        overlaps = (bin_starts[:, np.newaxis] < bl_ends) & (bl_starts < bin_ends[:, np.newaxis])
        bl_col_mask[offset:offset + n_bins] = overlaps.any(axis=1)

    print('  Blacklist filter:')
    print('    Blacklisted rows: {0}'.format(bl_row_mask.sum()))
    print('    Blacklisted cols: {0}'.format(bl_col_mask.sum()))

    return bl_row_mask, bl_col_mask


def build_crop_arrays(retain_rows, retain_cols):
    """Build rowMap, colMap, odd_indices, and even_indices arrays."""
    odd_global_indices = np.where(retain_rows)[0]
    even_global_indices = np.where(retain_cols)[0]

    odd_offsets_arr = np.array(ODD_OFFSETS + [N_ODD_BINS])
    odd_chrom_idxs = np.searchsorted(odd_offsets_arr, odd_global_indices, side='right') - 1
    rowMap = np.column_stack([
        odd_global_indices,
        np.array(ODD_CHROMS)[odd_chrom_idxs],
        odd_global_indices - np.array(ODD_OFFSETS)[odd_chrom_idxs]
    ]).astype(np.int64)

    even_offsets_arr = np.array(EVEN_OFFSETS + [N_EVEN_BINS])
    even_chrom_idxs = np.searchsorted(even_offsets_arr, even_global_indices, side='right') - 1
    colMap = np.column_stack([
        even_global_indices,
        np.array(EVEN_CHROMS)[even_chrom_idxs],
        even_global_indices - np.array(EVEN_OFFSETS)[even_chrom_idxs]
    ]).astype(np.int64)

    return rowMap, colMap, odd_global_indices, even_global_indices


def save_crop_maps(output_dir, rowMap, colMap, odd_indices, even_indices, timepoint=None):
    """Save cropMap and cropIndices as .mat files."""
    os.makedirs(output_dir, exist_ok=True)

    suffix = '_{0}'.format(timepoint) if timepoint else ''
    cropmap_path = os.path.join(output_dir, 'mm10_cropMap{0}.mat'.format(suffix))
    cropidx_path = os.path.join(output_dir, 'mm10_cropIndices{0}.mat'.format(suffix))

    savemat(cropmap_path, {'rowMap': rowMap, 'colMap': colMap})
    savemat(cropidx_path, {
        'odd_indices': odd_indices.reshape(1, -1),
        'even_indices': even_indices.reshape(1, -1)
    })

    print('  Saved: {0}'.format(cropmap_path))
    print('  Saved: {0}'.format(cropidx_path))

    return cropmap_path, cropidx_path


def main():
    args = parse_args()

    sniper_root = os.path.abspath(os.path.join(_script_dir, '..', 'repos', 'SNIPER'))
    if args.sizes_file is not None:
        sizes_file = args.sizes_file
    else:
        sizes_file = os.path.join(sniper_root, 'data', 'mm10.chrom.sizes')

    if not os.path.isfile(args.hic_path):
        print('ERROR: .hic file not found: {0}'.format(args.hic_path))
        sys.exit(1)
    if not os.path.isfile(args.blacklist_path):
        print('ERROR: blacklist not found: {0}'.format(args.blacklist_path))
        sys.exit(1)
    if not os.path.isfile(sizes_file):
        print('ERROR: chrom.sizes not found: {0}'.format(sizes_file))
        sys.exit(1)

    print('=========================================')
    print('B1: mm10 Crop Map Generation')
    print('=========================================')
    print('  HiC:       {0}'.format(args.hic_path))
    print('  Blacklist: {0}'.format(args.blacklist_path))
    print('  Output:    {0}'.format(args.output_dir))
    print('  Tmp dir:   {0}'.format(args.tmp_dir))
    print('  Sizes:     {0}'.format(sizes_file))
    print('  Overwrite: {0}'.format(args.overwrite))
    print('')

    # ── Step 1: Dump inter-chromosomal contacts ──
    print('Step 1: Dumping inter-chromosomal contacts (90 pairs)...')
    dump_interchromosomal_contacts(args.hic_path, args.tmp_dir, args.overwrite)

    dump_files = [f for f in os.listdir(args.tmp_dir)
                  if f.startswith('cropmap_chrm') and f.endswith('.txt')]
    non_empty = [f for f in dump_files
                 if os.path.getsize(os.path.join(args.tmp_dir, f)) > 0]
    print('  Dump inventory: {0} files, {1} non-empty'.format(
        len(dump_files), len(non_empty)))
    if len(non_empty) < len(ODD_CHROMS) * len(EVEN_CHROMS):
        missing = set('cropmap_chrm{0}_chrm{1}.txt'.format(o, e)
                      for o in ODD_CHROMS for e in EVEN_CHROMS) - set(dump_files)
        empty = set(dump_files) - set(non_empty)
        if missing:
            print('  Missing files: {0}'.format(sorted(missing)))
        if empty:
            print('  Empty files:   {0}'.format(sorted(empty)))

    # ── Step 2: Build full matrix ──
    print('Step 2: Building inter-chromosomal matrix...')
    M = build_matrix(args.tmp_dir, sizes_file)

    # ── Step 3: Coverage filter ──
    print('Step 3: Computing coverage filter...')
    sparse_rows, sparse_cols = compute_coverage_filter(M)

    # ── Step 4: Blacklist filter ──
    print('Step 4: Computing blacklist filter...')
    bl_row_mask, bl_col_mask = compute_blacklist_filter(args.blacklist_path)

    # ── Step 5: Combine and build crop arrays ──
    print('Step 5: Building crop map arrays...')
    exclude_rows = sparse_rows | bl_row_mask
    exclude_cols = sparse_cols | bl_col_mask
    retain_rows = ~exclude_rows
    retain_cols = ~exclude_cols

    rowMap, colMap, odd_indices, even_indices = build_crop_arrays(retain_rows, retain_cols)

    # ── Step 6: Save ──
    print('Step 6: Saving crop map files...')
    save_crop_maps(args.output_dir, rowMap, colMap, odd_indices, even_indices, args.timepoint)

    # ── Summary ──
    print('')
    print('=========================================')
    print('mm10 crop map generated successfully.')
    print('  Odd bins:  {0} retained / {1} total ({2:.1f}%)'.format(
        len(odd_indices), N_ODD_BINS, 100.0 * len(odd_indices) / N_ODD_BINS))
    print('  Even bins: {0} retained / {1} total ({2:.1f}%)'.format(
        len(even_indices), N_EVEN_BINS, 100.0 * len(even_indices) / N_EVEN_BINS))
    print('  Excluded rows: {0} sparse + {1} blacklist = {2} total'.format(
        sparse_rows.sum(), bl_row_mask.sum(), exclude_rows.sum()))
    print('  Excluded cols: {0} sparse + {1} blacklist = {2} total'.format(
        sparse_cols.sum(), bl_col_mask.sum(), exclude_cols.sum()))
    print('  rowMap shape: {0}'.format(rowMap.shape))
    print('  colMap shape: {0}'.format(colMap.shape))
    print('=========================================')


if __name__ == '__main__':
    main()
