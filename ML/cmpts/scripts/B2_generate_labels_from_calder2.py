# ML/cmpts/scripts/B2_generate_labels_from_calder2.py
# Stage B2: Convert CALDER2 subcompartment calls to SNIPER training labels.
#
# Reads variable-width CALDER2 segments (ctrl_merged), bins to 100kb via
# plurality vote, applies B1 crop mask, and saves as .mat for B3 training.
#
# Usage:
#   python B2_generate_labels_from_calder2.py <timepoint> <data_dir> <code_dir>
#     <timepoint> : 250402 | 250831
#     <data_dir>  : HPC data root (e.g. /expanse/.../sniper)
#     <code_dir>  : repo root (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

import os
import sys

import numpy as np
import pandas as pd
from scipy.io import loadmat, savemat

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'repos', 'SNIPER'))
from utilities.mm10_config import (
    ODD_CHROMS, EVEN_CHROMS,
    ODD_BIN_COUNTS, EVEN_BIN_COUNTS,
    ODD_OFFSETS, EVEN_OFFSETS,
    N_ODD_BINS, N_EVEN_BINS,
    RESOLUTION, N_CLASSES, TRAIN_VAL_SPLIT
)

LABEL_MAP = {'A.1': 0, 'A.2': 1, 'B.1': 2, 'B.2': 3}
LABEL_NAMES = {v: k for k, v in LABEL_MAP.items()}

TP_DIRS = {'250402': 'late', '250831': 'early'}

ALL_CHROMS = ['chr{}'.format(c) for c in sorted(ODD_CHROMS + EVEN_CHROMS)]
ODD_CHROM_STRS = ['chr{}'.format(c) for c in ODD_CHROMS]
EVEN_CHROM_STRS = ['chr{}'.format(c) for c in EVEN_CHROMS]

EXPECTED_COLS = ['chr', 'pos_start', 'pos_end', 'comp_name', 'comp_rank', 'continous_rank']


def load_calder2_tsv(tsv_path):
    if not os.path.isfile(tsv_path):
        print('ERROR: CALDER2 TSV not found: {}'.format(tsv_path), file=sys.stderr)
        sys.exit(1)
    if os.path.getsize(tsv_path) < 100:
        print('ERROR: CALDER2 TSV suspiciously small: {}'.format(tsv_path), file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(tsv_path, sep='\t')

    missing = set(EXPECTED_COLS) - set(df.columns)
    if missing:
        print('ERROR: Missing columns in TSV: {}'.format(missing), file=sys.stderr)
        sys.exit(1)

    chroms_found = set(df['chr'].unique())
    chroms_expected = set(ALL_CHROMS)
    missing_chroms = chroms_expected - chroms_found
    if missing_chroms:
        print('ERROR: Missing chromosomes in CALDER2 TSV: {}'.format(
            sorted(missing_chroms)), file=sys.stderr)
        sys.exit(1)

    df = df[df['chr'].isin(chroms_expected)].copy()
    print('  Loaded {} segments across {} chromosomes'.format(len(df), df['chr'].nunique()))
    return df


def truncate_to_depth2(comp_name_series):
    truncated = comp_name_series.str.extract(r'^([AB]\.[12])', expand=False)
    n_invalid = truncated.isna().sum()
    if n_invalid > 0:
        bad = comp_name_series[truncated.isna()].unique()[:5]
        print('WARNING: {} segments have non-standard labels: {}'.format(
            n_invalid, list(bad)), file=sys.stderr)
    return truncated


def plurality_vote_100kb(df, chrom_strs, offsets, bin_counts):
    df = df[df['chr'].isin(chrom_strs)].copy()
    if df.empty:
        print('ERROR: No segments found for chromosomes {}'.format(chrom_strs[:3]),
              file=sys.stderr)
        sys.exit(1)

    df['bin_first'] = ((df['pos_start'] - 1) // RESOLUTION).astype(np.int64)
    df['bin_last'] = ((df['pos_end'] - 1) // RESOLUTION).astype(np.int64)

    n_bins_each = (df['bin_last'] - df['bin_first'] + 1).to_numpy().astype(np.int64)
    repeated = df.loc[df.index.repeat(n_bins_each)].reset_index(drop=True)

    bin_indices = np.concatenate([
        np.arange(f, l + 1, dtype=np.int64)
        for f, l in zip(df['bin_first'].to_numpy(), df['bin_last'].to_numpy())
    ])
    repeated['bin_0based'] = bin_indices

    seg_start_0 = repeated['pos_start'].to_numpy() - 1
    seg_end_0 = repeated['pos_end'].to_numpy()
    bin_start_0 = repeated['bin_0based'].to_numpy() * RESOLUTION
    bin_end_0 = bin_start_0 + RESOLUTION

    repeated['overlap'] = (
        np.minimum(seg_end_0, bin_end_0) - np.maximum(seg_start_0, bin_start_0)
    )
    repeated = repeated[repeated['overlap'] > 0].copy()

    by_label = (repeated
                .groupby(['chr', 'bin_0based', 'label_d2'], sort=False)['overlap']
                .sum()
                .reset_index())

    winner_idx = by_label.groupby(['chr', 'bin_0based'], sort=False)['overlap'].idxmax()
    winners = by_label.loc[winner_idx, ['chr', 'bin_0based', 'label_d2']].copy()

    chrom_to_offset = dict(zip(chrom_strs, offsets))
    winners['global_idx'] = (
        winners['chr'].map(chrom_to_offset).to_numpy() + winners['bin_0based'].to_numpy()
    ).astype(np.int64)
    winners['label_int'] = winners['label_d2'].map(LABEL_MAP).astype(np.int32)

    n_total = sum(bin_counts)
    oob = winners[winners['global_idx'] >= n_total]
    if len(oob) > 0:
        print('WARNING: {} bins exceed expected total {} (max idx {})'.format(
            len(oob), n_total, winners['global_idx'].max()), file=sys.stderr)

    return dict(zip(winners['global_idx'].to_numpy(), winners['label_int'].to_numpy()))


def build_label_array(vote_dict, crop_indices, n_total_bins, offsets, bin_counts, group_name):
    full_labels = np.full(n_total_bins, -1, dtype=np.int32)

    for global_idx, label_int in vote_dict.items():
        if 0 <= global_idx < n_total_bins:
            full_labels[global_idx] = label_int

    n_voted = np.sum(full_labels >= 0)
    print('  {} bins labeled out of {} total (before crop)'.format(n_voted, n_total_bins))

    retained = full_labels[crop_indices]
    n_uncallable = np.sum(retained == -1)

    if n_uncallable > 0:
        print('  {} uncallable bins in crop mask — imputing with chromosome mode'.format(
            n_uncallable))
        if n_uncallable > 50:
            print('ERROR: Too many uncallable bins ({}) in {} crop mask'.format(
                n_uncallable, group_name), file=sys.stderr)
            sys.exit(1)

        chrom_modes = {}
        for i, offset in enumerate(offsets):
            end = offset + bin_counts[i]
            chrom_bins = full_labels[offset:end]
            valid = chrom_bins[chrom_bins >= 0]
            if len(valid) > 0:
                counts = np.bincount(valid, minlength=N_CLASSES)
                chrom_modes[i] = np.argmax(counts)

        for j in np.where(retained == -1)[0]:
            global_idx = crop_indices[j]
            chrom_idx = np.searchsorted(offsets, global_idx, side='right') - 1
            if chrom_idx in chrom_modes:
                retained[j] = chrom_modes[chrom_idx]

    still_missing = np.sum(retained == -1)
    if still_missing > 0:
        print('ERROR: {} bins still unlabeled after imputation'.format(still_missing),
              file=sys.stderr)
        sys.exit(1)

    return retained


def validate(rows, cols, odd_indices, even_indices):
    checks_passed = 0
    checks_total = 7

    print('\n--- Validation ---')

    # 1. Input existence (already validated upstream, so this is a formality)
    print('[1/7] Input existence: PASS (validated during load)')
    checks_passed += 1

    # 2. Chromosome completeness (validated during load)
    print('[2/7] Chromosome completeness: PASS (validated during load)')
    checks_passed += 1

    # 3. Label distribution
    print('[3/7] Label distribution:')
    all_ok = True
    for name, arr, group in [('rows (odd)', rows, 'odd'), ('cols (even)', cols, 'even')]:
        counts = np.bincount(arr, minlength=N_CLASSES)
        pcts = 100.0 * counts / len(arr)
        parts = []
        for i in range(N_CLASSES):
            parts.append('  {}={:.1f}%'.format(LABEL_NAMES[i], pcts[i]))
        print('  {}: {}'.format(name, ', '.join(parts)))
        if np.any(pcts < 5) or np.any(pcts > 60):
            print('  WARNING: Unusual label distribution in {}'.format(group))
            all_ok = False
    status = 'PASS' if all_ok else 'WARN'
    print('  Status: {}'.format(status))
    checks_passed += 1

    # 4. Dimension match
    ok = (len(rows) == len(odd_indices)) and (len(cols) == len(even_indices))
    print('[4/7] Dimension match: rows={} vs odd_indices={}, cols={} vs even_indices={} — {}'.format(
        len(rows), len(odd_indices), len(cols), len(even_indices),
        'PASS' if ok else 'FAIL'))
    if ok:
        checks_passed += 1

    # 5. No invalid labels
    valid_range = set(range(N_CLASSES))
    rows_ok = set(np.unique(rows)).issubset(valid_range)
    cols_ok = set(np.unique(cols)).issubset(valid_range)
    ok = rows_ok and cols_ok
    print('[5/7] All labels in {{0..{}}}: {}'.format(N_CLASSES - 1,
        'PASS' if ok else 'FAIL'))
    if ok:
        checks_passed += 1

    # 6. Coverage (informational — already checked in build_label_array)
    print('[6/7] Coverage: rows={} bins, cols={} bins (all labeled)'.format(
        len(rows), len(cols)))
    checks_passed += 1

    # 7. Training split
    ok = (len(rows) > TRAIN_VAL_SPLIT) and (len(cols) > TRAIN_VAL_SPLIT)
    print('[7/7] Training split: rows={} > {}, cols={} > {} — {}'.format(
        len(rows), TRAIN_VAL_SPLIT, len(cols), TRAIN_VAL_SPLIT,
        'PASS' if ok else 'FAIL'))
    if ok:
        checks_passed += 1

    print('\nResult: {}/{} checks passed'.format(checks_passed, checks_total))

    if checks_passed < checks_total:
        print('ERROR: Validation failed', file=sys.stderr)
        sys.exit(1)


def main():
    if len(sys.argv) != 4:
        print('Usage: python B2_generate_labels_from_calder2.py <timepoint> <data_dir> <code_dir>')
        sys.exit(1)

    tp = sys.argv[1]
    data_dir = sys.argv[2]
    code_dir = sys.argv[3]

    if tp not in TP_DIRS:
        print('ERROR: Unknown timepoint {}. Expected 250402 or 250831'.format(tp),
              file=sys.stderr)
        sys.exit(1)

    tp_dir = TP_DIRS[tp]

    print('===========================================')
    print('B2: Generate SNIPER Labels from CALDER2')
    print('===========================================')
    print('Timepoint:  {} ({})'.format(tp, tp_dir))
    print('Data dir:   {}'.format(data_dir))
    print('Code dir:   {}'.format(code_dir))
    print('')

    tsv_path = os.path.join(
        code_dir, 'outputs', 'calder2', tp_dir, tp,
        'ctrl_merged', 'sub_compartments', 'all_sub_compartments.tsv'
    )
    crop_idx_path = os.path.join(
        code_dir, 'repos', 'SNIPER', 'crop_map',
        'mm10_cropIndices_{}.mat'.format(tp)
    )
    output_dir = os.path.join(code_dir, 'outputs', 'sniper')
    output_path = os.path.join(output_dir, 'mm10_labels_{}.mat'.format(tp))

    print('CALDER2 TSV: {}'.format(tsv_path))
    print('Crop idx:    {}'.format(crop_idx_path))
    print('Output:      {}'.format(output_path))
    print('')

    # --- Load inputs ---

    print('Loading CALDER2 ctrl_merged subcompartments...')
    df = load_calder2_tsv(tsv_path)

    print('Truncating labels to depth 2...')
    df['label_d2'] = truncate_to_depth2(df['comp_name'])
    df = df.dropna(subset=['label_d2'])

    print('Loading crop indices...')
    if not os.path.isfile(crop_idx_path):
        print('ERROR: Crop indices not found: {}'.format(crop_idx_path), file=sys.stderr)
        sys.exit(1)

    crop_mat = loadmat(crop_idx_path)
    odd_indices = crop_mat['odd_indices'].flatten().astype(np.int64)
    even_indices = crop_mat['even_indices'].flatten().astype(np.int64)
    print('  odd_indices: {} bins, even_indices: {} bins'.format(
        len(odd_indices), len(even_indices)))

    # --- Plurality vote ---

    print('\nBinning odd chromosomes ({}) to 100kb...'.format(
        ', '.join('chr{}'.format(c) for c in ODD_CHROMS)))
    odd_votes = plurality_vote_100kb(
        df, ODD_CHROM_STRS, ODD_OFFSETS, ODD_BIN_COUNTS)
    print('  {} bins labeled'.format(len(odd_votes)))

    print('Binning even chromosomes ({}) to 100kb...'.format(
        ', '.join('chr{}'.format(c) for c in EVEN_CHROMS)))
    even_votes = plurality_vote_100kb(
        df, EVEN_CHROM_STRS, EVEN_OFFSETS, EVEN_BIN_COUNTS)
    print('  {} bins labeled'.format(len(even_votes)))

    # --- Build label arrays ---

    print('\nBuilding odd-chromosome label array...')
    rows = build_label_array(
        odd_votes, odd_indices, N_ODD_BINS,
        ODD_OFFSETS, ODD_BIN_COUNTS, 'odd')

    print('Building even-chromosome label array...')
    cols = build_label_array(
        even_votes, even_indices, N_EVEN_BINS,
        EVEN_OFFSETS, EVEN_BIN_COUNTS, 'even')

    # --- Validate ---

    validate(rows, cols, odd_indices, even_indices)

    # --- Save ---

    os.makedirs(output_dir, exist_ok=True)

    savemat(output_path, {
        'rows': rows.reshape(1, -1),
        'cols': cols.reshape(1, -1)
    })

    print('\nSaved: {}'.format(output_path))
    print('  rows: shape={}, dtype={}'.format(rows.shape, rows.dtype))
    print('  cols: shape={}, dtype={}'.format(cols.shape, cols.dtype))
    print('\nB2 complete.')


if __name__ == '__main__':
    main()
