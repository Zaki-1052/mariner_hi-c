#!/usr/bin/env python3
# scripts/generate_final_results.py
"""
Generate final_results.tsv from existing all_results_primary.tsv files.

Applies stringent thresholds without re-running edgeR analysis:
- |logFC| > 0.3
- FDR < 0.03

Usage:
    python scripts/generate_final_results.py

Processes all resolutions: 5kb, 10kb, 25kb
"""

import pandas as pd
import sys
from pathlib import Path


def filter_final_results(all_results_df, logfc_threshold=0.3, fdr_threshold=0.03):
    """
    Apply stringent filters to generate final results.

    Args:
        all_results_df: DataFrame with all edgeR results
        logfc_threshold: Minimum absolute logFC (default: 0.3)
        fdr_threshold: Maximum FDR (default: 0.03)

    Returns:
        Filtered DataFrame
    """
    # Apply filters
    final_df = all_results_df[
        (abs(all_results_df['logFC']) > logfc_threshold) &
        (all_results_df['FDR'] < fdr_threshold)
    ].copy()

    return final_df


def process_resolution(resolution_kb, logfc_threshold=0.3, fdr_threshold=0.03):
    """
    Process a single resolution: load, filter, save.

    Args:
        resolution_kb: Resolution in kb (5, 10, or 25)
        logfc_threshold: Minimum absolute logFC
        fdr_threshold: Maximum FDR

    Returns:
        True if successful, False if file not found
    """
    # Define paths
    results_dir = f"outputs/edgeR_results_res_{resolution_kb}kb/primary_analysis"
    input_file = Path(results_dir) / "all_results_primary.tsv"
    output_file = Path(results_dir) / "final_results.tsv"

    # Check if input exists
    if not input_file.exists():
        print(f"  WARNING: Input file not found: {input_file}")
        return False

    # Load all results
    print(f"  Loading: {input_file.name}")
    all_results = pd.read_csv(input_file, sep='\t')
    print(f"    Total loops: {len(all_results)}")

    # Apply filters
    print(f"  Applying filters: |logFC| > {logfc_threshold}, FDR < {fdr_threshold}")
    final_results = filter_final_results(all_results, logfc_threshold, fdr_threshold)

    n_filtered = len(final_results)
    pct_retained = (n_filtered / len(all_results)) * 100 if len(all_results) > 0 else 0

    print(f"    Loops passing filters: {n_filtered} ({pct_retained:.1f}%)")

    if n_filtered == 0:
        print(f"    WARNING: No loops passed filters!")
        return False

    # Statistics
    up_count = (final_results['logFC'] > 0).sum()
    down_count = (final_results['logFC'] < 0).sum()
    logfc_min = final_results['logFC'].min()
    logfc_max = final_results['logFC'].max()
    fdr_min = final_results['FDR'].min()
    fdr_max = final_results['FDR'].max()

    print(f"    Direction: {up_count} up, {down_count} down")
    print(f"    logFC range: [{logfc_min:.3f}, {logfc_max:.3f}]")
    print(f"    FDR range: [{fdr_min:.6f}, {fdr_max:.6f}]")

    # Write output
    print(f"  Writing: {output_file.name}")
    final_results.to_csv(output_file, sep='\t', index=False)

    file_size_kb = output_file.stat().st_size / 1024
    print(f"    File size: {file_size_kb:.1f} KB")
    print(f"  ✓ Successfully generated final_results.tsv for {resolution_kb}kb")

    return True


def main():
    """Main execution function."""
    print("="*70)
    print("Generate Final Results from edgeR Output")
    print("="*70)
    print("\nFilters: |logFC| > 0.3 AND FDR < 0.03")
    print("Input:   all_results_primary.tsv (per resolution)")
    print("Output:  final_results.tsv (per resolution)")
    print()

    # Configuration
    resolutions = [5, 10, 25]  # kb
    logfc_threshold = 0.3
    fdr_threshold = 0.03

    # Process each resolution
    success_count = 0
    for res_kb in resolutions:
        print("="*70)
        print(f"Resolution: {res_kb}kb")
        print("="*70)

        success = process_resolution(res_kb, logfc_threshold, fdr_threshold)
        if success:
            success_count += 1

        print()

    # Summary
    print("="*70)
    print("Summary")
    print("="*70)
    print(f"Resolutions processed successfully: {success_count}/{len(resolutions)}")

    if success_count == 0:
        print("\nERROR: No final_results.tsv files were generated!")
        print("Make sure you've run the edgeR analysis first.")
        sys.exit(1)
    elif success_count < len(resolutions):
        print(f"\nWARNING: Only {success_count} out of {len(resolutions)} resolutions processed.")
        print("Check that edgeR analysis completed for all resolutions.")
    else:
        print("\n✓ All resolutions processed successfully!")

    print("\nGenerated files:")
    for res_kb in resolutions:
        output_file = Path(f"outputs/edgeR_results_res_{res_kb}kb/primary_analysis/final_results.tsv")
        if output_file.exists():
            size_kb = output_file.stat().st_size / 1024
            print(f"  - {output_file} ({size_kb:.1f} KB)")

    print("\nNext steps:")
    print("  1. Convert to BEDPE: bash scripts/convert_final_bedpe.sh")
    print("  2. Visualize in Juicebox\n")


if __name__ == '__main__':
    main()
