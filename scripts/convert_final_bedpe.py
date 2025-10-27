#!/usr/bin/env python3
# scripts/convert_final_bedpe.py
"""
Convert edgeR results to BEDPE format for Juicebox visualization.

Processes two file types:
1. final_results.tsv (|logFC| > 0.3, FDR < 0.03) - differential loops only
2. all_results_primary.tsv - all tested loops

Creates resolution-specific BEDPE files plus merged unions across all resolutions.

Usage:
    python scripts/convert_final_bedpe.py

Outputs:
    outputs/bedpe_final/
    - 5kb_final.bedpe, 10kb_final.bedpe, 25kb_final.bedpe
    - merged_final_resolutions.bedpe
    - 5kb_all_loops.bedpe, 10kb_all_loops.bedpe, 25kb_all_loops.bedpe
    - merged_all_loops.bedpe
"""

import pandas as pd
import sys
from pathlib import Path


def load_results_file(resolution_kb, file_type):
    """
    Load TSV results file for a specific resolution.

    Args:
        resolution_kb: Resolution in kb (5, 10, or 25)
        file_type: Either 'final_results' or 'all_results_primary'
    """
    filepath = f"outputs/edgeR_results_res_{resolution_kb}kb/primary_analysis/{file_type}.tsv"

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


def process_file_type(file_type, file_description, output_suffix, resolutions, output_dir):
    """
    Process a specific file type (final_results or all_results_primary).

    Args:
        file_type: TSV filename (e.g., 'final_results' or 'all_results_primary')
        file_description: Human-readable description
        output_suffix: Suffix for output files (e.g., 'final' or 'all_loops')
        resolutions: List of resolutions in kb
        output_dir: Output directory path
    """
    print("\n" + "="*70)
    print(f"Processing: {file_description}")
    print("="*70)

    # Load data from all resolutions
    print(f"\nLoading {file_type}.tsv files...")
    print("-" * 70)
    res_data = {}
    for res_kb in resolutions:
        print(f"\n{res_kb}kb resolution:")
        res_data[res_kb] = load_results_file(res_kb, file_type)

    # Check if any data loaded
    if all(df is None for df in res_data.values()):
        print(f"\nWARNING: No {file_type}.tsv files found!")
        print("Skipping this file type...")
        return

    print("\n" + "="*70)
    print(f"Converting {file_description} to BEDPE format")
    print("="*70)

    # Convert each resolution to BEDPE
    for res_kb in resolutions:
        if res_data[res_kb] is not None:
            print(f"\n{res_kb}kb resolution:")
            bedpe = convert_to_bedpe(res_data[res_kb], res_kb)
            output_path = f"{output_dir}/{res_kb}kb_{output_suffix}.bedpe"
            write_bedpe(bedpe, output_path, f"{res_kb}kb {file_description}")

    # Merge all resolutions
    merged = merge_resolutions(res_data)

    if merged is not None:
        bedpe_merged = convert_to_bedpe(merged)
        output_path = f"{output_dir}/merged_{output_suffix}.bedpe"
        write_bedpe(bedpe_merged, output_path, f"Merged {file_description} (all resolutions)")


def main():
    """Main execution function."""
    print("="*70)
    print("edgeR Results BEDPE Conversion")
    print("="*70)
    print("\nConverts edgeR TSV results to BEDPE format for Juicebox visualization")
    print("\nFile types processed:")
    print("  1. final_results.tsv - Differential loops (|logFC| > 0.3, FDR < 0.03)")
    print("  2. all_results_primary.tsv - All tested loops")
    print("\nOutput: Resolution-specific + merged union BEDPE files\n")

    output_dir = "outputs/bedpe_final"
    resolutions = [5, 10, 25]  # kb

    # Process final_results.tsv (differential loops only)
    process_file_type(
        file_type='final_results',
        file_description='differential loops (|logFC| > 0.3, FDR < 0.03)',
        output_suffix='final',
        resolutions=resolutions,
        output_dir=output_dir
    )

    # Process all_results_primary.tsv (all tested loops)
    process_file_type(
        file_type='all_results_primary',
        file_description='all tested loops',
        output_suffix='all_loops',
        resolutions=resolutions,
        output_dir=output_dir
    )

    print("\n" + "="*70)
    print("Conversion Complete!")
    print("="*70)
    print(f"\nOutput directory: {output_dir}/")
    print("\nGenerated files:")
    bedpe_files = sorted(Path(output_dir).glob("*.bedpe"))
    if bedpe_files:
        for f in bedpe_files:
            size_kb = f.stat().st_size / 1024
            print(f"  - {f.name:<35} ({size_kb:>8.1f} KB)")
    else:
        print("  (No BEDPE files generated)")

    print("\nTo visualize in Juicebox:")
    print("  1. Open .hic file in Juicebox")
    print("  2. Annotations > Load Basic Annotations")
    print("  3. Select BEDPE file")
    print("  4. Hover over loops to see statistics\n")


if __name__ == '__main__':
    main()
