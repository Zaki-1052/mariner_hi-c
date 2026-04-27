#!/usr/bin/env python3
# cluster/scripts/01_build_loop_count_file.py
# Build Popay-format loop count file + edgeR stats sidecar from the
# nonredundant multi-resolution merged BEDPE plus per-resolution count
# matrices. Adapts BAP1-KO cerebellum data (250402 late) to the input
# format expected by Popay's HiC_cluster3.ipynb / clustering pipeline.

import os
import sys

import pandas as pd


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

MERGED_BEDPE = os.path.join(
    REPO_ROOT,
    'outputs/250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe',
)
COUNTS_TEMPLATE = os.path.join(
    REPO_ROOT,
    'outputs/250402-late_outputs/res_{kb}kb/06_counts_matrix.tsv',
)
EDGER_TEMPLATE = os.path.join(
    REPO_ROOT,
    'outputs/250402-late_outputs/edgeR_results_res_{kb}kb/primary_analysis/all_results_primary.tsv',
)

OUT_COUNTS = os.path.join(REPO_ROOT, 'cluster/data/late_merged_loop_counts.txt')
OUT_META = os.path.join(REPO_ROOT, 'cluster/data/late_merged_loop_metadata.tsv')

EXPECTED_TOTAL = 39344
EXPECTED_PER_RES = {5: 7901, 10: 14553, 25: 16890}

CTRL_REPS = ['ctrl_M1', 'ctrl_M2', 'ctrl_M3']
MUT_REPS = ['mut_M1', 'mut_M2', 'mut_M3']
COORD_COLS = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
COORD_INT_COLS = ['start1', 'end1', 'start2', 'end2']
META_COLS = [
    'logFC', 'FDR', 'PValue', 'direction', 'significant', 'category',
    'resolution', 'resolution_kb', 'kept_from_resolution', 'is_multi_resolution',
]


def fail(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def load_per_resolution_counts(kb):
    counts_path = COUNTS_TEMPLATE.format(kb=kb)
    edger_path = EDGER_TEMPLATE.format(kb=kb)

    counts = pd.read_csv(counts_path, sep='\t', index_col=0)
    counts.index.name = 'loop_id'

    expected_reps = CTRL_REPS + MUT_REPS
    missing = [c for c in expected_reps if c not in counts.columns]
    if missing:
        fail(f"{counts_path} missing replicate columns: {missing}")

    counts['ctrl_merge'] = counts[CTRL_REPS].mean(axis=1)
    counts['mut_merge'] = counts[MUT_REPS].mean(axis=1)
    counts = counts[['ctrl_merge', 'mut_merge']].reset_index()

    edger = pd.read_csv(edger_path, sep='\t', usecols=['loop_id'] + COORD_COLS)

    joined = edger.merge(counts, on='loop_id', how='inner', validate='one_to_one')
    for col in COORD_INT_COLS:
        joined[col] = joined[col].astype('int64')
    return joined


def main():
    print(f"[1/5] Loading merged BEDPE: {MERGED_BEDPE}")
    merged = pd.read_csv(MERGED_BEDPE, sep='\t', header=0)
    print(f"      rows: {len(merged)}")

    if len(merged) != EXPECTED_TOTAL:
        fail(f"Expected {EXPECTED_TOTAL} merged loops, got {len(merged)}")

    for col in COORD_COLS + ['kept_from_resolution']:
        if col not in merged.columns:
            fail(f"Merged BEDPE missing column: {col}")

    for col in COORD_INT_COLS:
        merged[col] = merged[col].astype('int64')
    merged['kept_from_resolution'] = merged['kept_from_resolution'].astype('int64')

    res_counts = merged['kept_from_resolution'].value_counts().to_dict()
    print(f"      kept_from_resolution counts: {res_counts}")
    for kb, expected in EXPECTED_PER_RES.items():
        observed = res_counts.get(kb, 0)
        if observed != expected:
            fail(f"kept_from_resolution={kb}: expected {expected}, got {observed}")

    print(f"[2/5] Loading per-resolution counts (5/10/25 kb)...")
    per_res = {}
    for kb in [5, 10, 25]:
        per_res[kb] = load_per_resolution_counts(kb)
        print(f"      res {kb}kb: {len(per_res[kb])} loops with merged-replicate counts")

    print(f"[3/5] Per-resolution coordinate join...")
    pieces = []
    for kb in [5, 10, 25]:
        sub = merged[merged['kept_from_resolution'] == kb].copy()
        joined = sub.merge(
            per_res[kb][COORD_COLS + ['ctrl_merge', 'mut_merge']],
            on=COORD_COLS,
            how='left',
            validate='one_to_one',
        )
        unmatched = joined[joined['ctrl_merge'].isna()]
        if len(unmatched) > 0:
            print(f"FATAL: {len(unmatched)} merged-BEDPE rows at "
                  f"kept_from_resolution={kb}kb did not match per-resolution "
                  f"counts. First 10 unmatched coordinates:", file=sys.stderr)
            print(unmatched[COORD_COLS].head(10).to_string(), file=sys.stderr)
            sys.exit(1)
        pieces.append(joined)

    final = pd.concat(pieces, ignore_index=True)
    if len(final) != EXPECTED_TOTAL:
        fail(f"Concat row count: {len(final)} (expected {EXPECTED_TOTAL})")

    print(f"[4/5] Writing count file: {OUT_COUNTS}")
    counts_out = (
        final[COORD_COLS + ['ctrl_merge', 'mut_merge']]
        .rename(columns={'start1': 'x1', 'end1': 'x2', 'start2': 'y1', 'end2': 'y2'})
        [['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'ctrl_merge', 'mut_merge']]
    )

    n_nan = int(counts_out.isna().sum().sum())
    if n_nan != 0:
        fail(f"NaN cells in count output: {n_nan}")
    if not (counts_out['ctrl_merge'] > 0).all():
        fail("non-positive ctrl_merge values present")
    if not (counts_out['mut_merge'] > 0).all():
        fail("non-positive mut_merge values present")

    counts_out.to_csv(OUT_COUNTS, sep='\t', index=False, lineterminator='\n')
    print(f"      wrote {len(counts_out)} rows × {counts_out.shape[1]} cols")

    print(f"[5/5] Writing metadata sidecar: {OUT_META}")
    missing_meta = [c for c in META_COLS if c not in final.columns]
    if missing_meta:
        fail(f"Merged BEDPE missing metadata columns: {missing_meta}")

    meta = (
        final[COORD_COLS + META_COLS]
        .rename(columns={'start1': 'x1', 'end1': 'x2', 'start2': 'y1', 'end2': 'y2'})
        [['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2'] + META_COLS]
    )
    meta.to_csv(OUT_META, sep='\t', index=False, lineterminator='\n')
    print(f"      wrote {len(meta)} rows × {meta.shape[1]} cols")

    print()
    print("=== Summary ===")
    print(f"Count file: {OUT_COUNTS}")
    print(f"  shape: {counts_out.shape}")
    print("  ctrl_merge distribution:")
    print(counts_out['ctrl_merge']
          .describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99]).to_string())
    print("  mut_merge distribution:")
    print(counts_out['mut_merge']
          .describe(percentiles=[0.01, 0.05, 0.50, 0.95, 0.99]).to_string())
    print()
    print(f"Metadata sidecar: {OUT_META}")
    print(f"  shape: {meta.shape}")
    print(f"  direction counts: {meta['direction'].value_counts().to_dict()}")
    print(f"  significant counts: {meta['significant'].value_counts().to_dict()}")
    print()
    print("Phase 1.1 complete.")


if __name__ == '__main__':
    main()
