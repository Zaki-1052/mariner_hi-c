#!/usr/bin/env python3
# scripts/create_nonredundant_bedpe.py
"""
Create non-redundant BEDPE from merged_all_results.tsv

This script converts the properly de-duplicated merged_all_results.tsv
(created by downstream_analysis.R with overlap removal) to BEDPE format.

The original convert_final_bedpe.py just concatenates per-resolution files
without removing overlaps. This script uses the correct non-redundant source.

Usage:
    python scripts/create_nonredundant_bedpe.py [input_tsv] [output_bedpe]

Defaults:
    input:  25042-late_outputs/merged_loops/merged_all_results.tsv
    output: 25042-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe
"""

import pandas as pd
import sys
from pathlib import Path


def convert_merged_results_to_bedpe(input_tsv: str, output_bedpe: str) -> None:
    """Convert merged_all_results.tsv to BEDPE format."""

    print("=" * 60)
    print("Non-redundant BEDPE Conversion")
    print("=" * 60)

    # Load non-redundant results
    print(f"\nLoading: {input_tsv}")
    df = pd.read_csv(input_tsv, sep='\t')
    print(f"  Loaded {len(df)} non-redundant loops")

    # Check required columns
    required_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Define BEDPE output columns (coordinates first, then useful metadata)
    coord_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']

    # Select informative columns for BEDPE (avoid internal indices)
    exclude_cols = {'loop_id', 'anchor1_nearest_gene_idx', 'anchor2_nearest_gene_idx'}
    meta_cols = [c for c in df.columns if c not in coord_cols and c not in exclude_cols]

    # Reorder: coordinates first
    output_cols = coord_cols + meta_cols
    bedpe = df[output_cols].copy()

    # Ensure coordinates are integers
    for col in ['start1', 'end1', 'start2', 'end2']:
        bedpe[col] = bedpe[col].astype(int)

    # Create output directory if needed
    Path(output_bedpe).parent.mkdir(parents=True, exist_ok=True)

    # Write BEDPE with header
    bedpe.to_csv(output_bedpe, sep='\t', index=False, header=True, na_rep='.')

    # Summary statistics
    print(f"\nOutput: {output_bedpe}")
    print(f"  Total loops: {len(bedpe)}")
    print(f"  File size: {Path(output_bedpe).stat().st_size / 1024:.1f} KB")

    if 'logFC' in bedpe.columns:
        up = (bedpe['logFC'] > 0).sum()
        down = (bedpe['logFC'] < 0).sum()
        print(f"\n  logFC range: [{bedpe['logFC'].min():.3f}, {bedpe['logFC'].max():.3f}]")
        print(f"  Up in mutant: {up} ({100*up/len(bedpe):.1f}%)")
        print(f"  Down in mutant: {down} ({100*down/len(bedpe):.1f}%)")

    if 'FDR' in bedpe.columns:
        sig_005 = (bedpe['FDR'] < 0.05).sum()
        sig_003 = (bedpe['FDR'] < 0.03).sum()
        print(f"\n  FDR < 0.05: {sig_005} ({100*sig_005/len(bedpe):.1f}%)")
        print(f"  FDR < 0.03: {sig_003} ({100*sig_003/len(bedpe):.1f}%)")

    if 'resolution_kb' in bedpe.columns:
        print(f"\n  Resolution breakdown:")
        for res in sorted(bedpe['resolution_kb'].unique()):
            count = (bedpe['resolution_kb'] == res).sum()
            print(f"    {int(res)}kb: {count} ({100*count/len(bedpe):.1f}%)")

    if 'is_multi_resolution' in bedpe.columns:
        multi = bedpe['is_multi_resolution'].sum()
        print(f"\n  Multi-resolution loops: {multi} ({100*multi/len(bedpe):.1f}%)")

    print("\n" + "=" * 60)
    print("Conversion complete!")
    print("=" * 60 + "\n")


def main():
    # Default paths
    default_input = "25042-late_outputs/merged_loops/merged_all_results.tsv"
    default_output = "25042-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe"

    # Parse arguments
    if len(sys.argv) >= 3:
        input_tsv = sys.argv[1]
        output_bedpe = sys.argv[2]
    elif len(sys.argv) == 2:
        input_tsv = sys.argv[1]
        output_bedpe = default_output
    else:
        input_tsv = default_input
        output_bedpe = default_output

    # Validate input exists
    if not Path(input_tsv).exists():
        print(f"ERROR: Input file not found: {input_tsv}")
        sys.exit(1)

    convert_merged_results_to_bedpe(input_tsv, output_bedpe)


if __name__ == '__main__':
    main()
