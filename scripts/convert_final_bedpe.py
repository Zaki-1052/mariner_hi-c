#!/usr/bin/env python3
# scripts/convert_final_bedpe.py
"""
Convert final_results.tsv files (|logFC| > 0.3, FDR < 0.03) to BEDPE format.

Creates resolution-specific BEDPE files plus a merged union across all resolutions.

Usage:
    python scripts/convert_final_bedpe.py

Outputs:
    - outputs/bedpe_final/5kb_final.bedpe
    - outputs/bedpe_final/10kb_final.bedpe
    - outputs/bedpe_final/25kb_final.bedpe
    - outputs/bedpe_final/merged_all_resolutions.bedpe
"""

import pandas as pd
import sys
from pathlib import Path


def load_final_results(resolution_kb):
    """Load final_results.tsv for a specific resolution."""
    filepath = f"outputs/edgeR_results_res_{resolution_kb}kb/primary_analysis/final_results.tsv"

    if not Path(filepath).exists():
        print(f"  WARNING: File not found: {filepath}")
        return None

    df = pd.read_csv(filepath, sep='\t')
    print(f"  Loaded {len(df)} loops from {resolution_kb}kb resolution")

    # Add resolution column for tracking
    df['resolution_kb'] = resolution_kb

    return df


def convert_to_bedpe(df, resolution_kb=None):
    """
    Convert dataframe to BEDPE format.

    Format: chr1, start1, end1, chr2, start2, end2, [additional columns...]
    """
    # Required coordinate columns
    coord_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']

    # Additional columns to include (exclude internal IDs)
    additional_cols = [col for col in df.columns
                      if col not in coord_cols and col != 'loop_id']

    # Reorder: coordinates first, then all other columns
    output_cols = coord_cols + additional_cols
    bedpe = df[output_cols].copy()

    # Ensure coordinates are integers
    for col in ['start1', 'end1', 'start2', 'end2']:
        bedpe[col] = bedpe[col].astype(int)

    if resolution_kb:
        print(f"  Output columns ({resolution_kb}kb): {', '.join(output_cols)}")

    return bedpe


def write_bedpe(bedpe, output_path, description):
    """Write BEDPE dataframe to file with summary."""
    # Create output directory
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Write with header for Juicebox tooltip parsing
    bedpe.to_csv(
        output_path,
        sep='\t',
        index=False,
        header=True,
        na_rep='.'
    )

    print(f"\n✓ {description}")
    print(f"  File: {output_path}")
    print(f"  Loops: {len(bedpe)}")
    print(f"  Size: {Path(output_path).stat().st_size / 1024:.1f} KB")

    # Print statistics
    if 'logFC' in bedpe.columns:
        print(f"  logFC range: [{bedpe['logFC'].min():.3f}, {bedpe['logFC'].max():.3f}]")
        up_count = (bedpe['logFC'] > 0).sum()
        down_count = (bedpe['logFC'] < 0).sum()
        print(f"  Up-regulated: {up_count} ({up_count/len(bedpe)*100:.1f}%)")
        print(f"  Down-regulated: {down_count} ({down_count/len(bedpe)*100:.1f}%)")

    if 'FDR' in bedpe.columns:
        print(f"  FDR range: [{bedpe['FDR'].min():.6f}, {bedpe['FDR'].max():.6f}]")


def merge_resolutions(res_data):
    """
    Merge loops from multiple resolutions into a union.

    Concatenates all loops across resolutions. Same genomic loop may appear
    multiple times with different bin coordinates (5kb vs 10kb vs 25kb binning).
    """
    print("\n" + "="*50)
    print("Merging resolutions into union")
    print("="*50)

    all_loops = []
    for res_kb, df in res_data.items():
        if df is not None:
            print(f"  Adding {len(df)} loops from {res_kb}kb")
            all_loops.append(df)

    if not all_loops:
        print("  ERROR: No data to merge!")
        return None

    # Concatenate all resolutions
    merged = pd.concat(all_loops, ignore_index=True)

    print(f"\n  Total loops in union: {len(merged)}")
    print(f"  Resolutions included: {', '.join([str(k) for k in res_data.keys() if res_data[k] is not None])}")

    # Resolution breakdown
    if 'resolution_kb' in merged.columns:
        print("\n  Breakdown by resolution:")
        for res_kb in sorted(merged['resolution_kb'].unique()):
            count = (merged['resolution_kb'] == res_kb).sum()
            print(f"    {res_kb}kb: {count} loops ({count/len(merged)*100:.1f}%)")

    return merged


def main():
    """Main execution function."""
    print("="*50)
    print("Final Results BEDPE Conversion")
    print("="*50)
    print("\nInput: final_results.tsv (|logFC| > 0.3, FDR < 0.03)")
    print("Output: Resolution-specific + merged BEDPE files\n")

    output_dir = "outputs/bedpe_final"
    resolutions = [5, 10, 25]  # kb

    # Load data from all resolutions
    print("Loading final_results.tsv files...")
    print("-" * 50)
    res_data = {}
    for res_kb in resolutions:
        print(f"\n{res_kb}kb resolution:")
        res_data[res_kb] = load_final_results(res_kb)

    # Check if any data loaded
    if all(df is None for df in res_data.values()):
        print("\nERROR: No final_results.tsv files found!")
        print("Make sure you've run the edgeR analysis first.")
        sys.exit(1)

    print("\n" + "="*50)
    print("Converting to BEDPE format")
    print("="*50)

    # Convert each resolution to BEDPE
    for res_kb in resolutions:
        if res_data[res_kb] is not None:
            print(f"\n{res_kb}kb resolution:")
            bedpe = convert_to_bedpe(res_data[res_kb], res_kb)
            output_path = f"{output_dir}/{res_kb}kb_final.bedpe"
            write_bedpe(bedpe, output_path, f"{res_kb}kb final results")

    # Merge all resolutions
    merged = merge_resolutions(res_data)

    if merged is not None:
        bedpe_merged = convert_to_bedpe(merged)
        output_path = f"{output_dir}/merged_all_resolutions.bedpe"
        write_bedpe(bedpe_merged, output_path, "Merged union (all resolutions)")

    print("\n" + "="*50)
    print("Conversion Complete!")
    print("="*50)
    print(f"\nOutput directory: {output_dir}/")
    print("\nGenerated files:")
    bedpe_files = sorted(Path(output_dir).glob("*.bedpe"))
    for f in bedpe_files:
        print(f"  - {f.name}")

    print("\nTo visualize in Juicebox:")
    print("  1. Open .hic file in Juicebox")
    print("  2. Annotations > Load Basic Annotations")
    print("  3. Select BEDPE file")
    print("  4. Hover over loops to see statistics\n")


if __name__ == '__main__':
    main()
