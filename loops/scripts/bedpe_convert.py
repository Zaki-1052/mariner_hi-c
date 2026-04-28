#!/usr/bin/env python3
# scripts/convert_to_bedpe.py
"""
Convert edgeR loop analysis results (TSV) to BEDPE format for Juicebox visualization.

Outputs coordinates (chr1, start1, end1, chr2, start2, end2) followed by all
statistical columns for Juicebox tooltip display.

Usage:
    python scripts/convert_to_bedpe.py <input.tsv> <output.bedpe> [options]

Examples:
    # All differential loops (up + down, FDR<0.05)
    python scripts/convert_to_bedpe.py \\
        outputs/edgeR_results_res_5kb/primary_analysis/all_results_primary.tsv \\
        outputs/bedpe/5kb_all_differential.bedpe \\
        --fdr-cutoff 0.05 --exclude-unchanged

    # Only up-regulated loops
    python scripts/convert_to_bedpe.py \\
        outputs/edgeR_results_res_5kb/primary_analysis/significant_loops_fdr05.tsv \\
        outputs/bedpe/5kb_up.bedpe \\
        --direction up_in_mutant

    # Only down-regulated loops
    python scripts/convert_to_bedpe.py \\
        outputs/edgeR_results_res_5kb/primary_analysis/significant_loops_fdr05.tsv \\
        outputs/bedpe/5kb_down.bedpe \\
        --direction down_in_mutant
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Convert edgeR loop results to BEDPE format for Juicebox',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'input_tsv',
        type=str,
        help='Input TSV file from edgeR analysis'
    )

    parser.add_argument(
        'output_bedpe',
        type=str,
        help='Output BEDPE file path'
    )

    parser.add_argument(
        '--fdr-cutoff',
        type=float,
        default=None,
        help='Filter loops by FDR threshold (e.g., 0.05)'
    )

    parser.add_argument(
        '--min-logfc',
        type=float,
        default=None,
        help='Filter loops by minimum absolute logFC (e.g., 0.5)'
    )

    parser.add_argument(
        '--direction',
        choices=['up', 'down', 'up_in_mutant', 'down_in_mutant', 'unchanged'],
        default=None,
        help='Filter loops by direction of change'
    )

    parser.add_argument(
        '--exclude-unchanged',
        action='store_true',
        help='Exclude loops with direction=unchanged (keep only differential)'
    )

    parser.add_argument(
        '--significant-only',
        action='store_true',
        help='Only include statistically significant loops (significant=TRUE)'
    )

    parser.add_argument(
        '--keep-columns',
        type=str,
        default=None,
        help='Comma-separated list of columns to keep after coordinates (default: keep all)'
    )

    return parser.parse_args()


def validate_input_file(filepath):
    """Validate input TSV file exists and has required columns."""
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    df = pd.read_csv(filepath, sep='\t', nrows=1)
    required_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    return True


def load_loops(filepath):
    """Load loop data from TSV file."""
    print(f"Loading loops from: {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    print(f"  Loaded {len(df)} loops")
    return df


def filter_loops(df, args):
    """Apply filtering criteria to loop dataframe."""
    initial_count = len(df)

    # Filter by significance
    if args.significant_only and 'significant' in df.columns:
        df = df[df['significant'] == True]
        print(f"  After significance filter: {len(df)} loops")

    # Filter by FDR
    if args.fdr_cutoff is not None and 'FDR' in df.columns:
        df = df[df['FDR'] < args.fdr_cutoff]
        print(f"  After FDR < {args.fdr_cutoff} filter: {len(df)} loops")

    # Filter by logFC
    if args.min_logfc is not None and 'logFC' in df.columns:
        df = df[abs(df['logFC']) >= args.min_logfc]
        print(f"  After |logFC| >= {args.min_logfc} filter: {len(df)} loops")

    # Exclude unchanged loops
    if args.exclude_unchanged and 'direction' in df.columns:
        df = df[df['direction'] != 'unchanged']
        print(f"  After excluding unchanged: {len(df)} loops")

    # Filter by direction
    if args.direction is not None and 'direction' in df.columns:
        df = df[df['direction'] == args.direction]
        print(f"  After direction = {args.direction} filter: {len(df)} loops")

    filtered_count = initial_count - len(df)
    if filtered_count > 0:
        print(f"  Filtered out {filtered_count} loops ({filtered_count/initial_count*100:.1f}%)")

    return df


def convert_to_bedpe(df, args):
    """
    Convert dataframe to BEDPE format for Juicebox.

    Format: chr1, start1, end1, chr2, start2, end2, [additional columns...]
    First 6 columns are coordinates, remaining columns show in Juicebox tooltip.
    """
    print(f"\nConverting {len(df)} loops to BEDPE format...")

    # Start with coordinate columns
    coord_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']

    # Determine which additional columns to include
    if args.keep_columns:
        # User specified columns
        additional_cols = [col.strip() for col in args.keep_columns.split(',')]
        additional_cols = [col for col in additional_cols if col in df.columns]
    else:
        # Keep all columns except coordinates and loop_id
        additional_cols = [col for col in df.columns
                          if col not in coord_cols and col != 'loop_id']

    # Create output dataframe
    output_cols = coord_cols + additional_cols
    bedpe = df[output_cols].copy()

    # Ensure coordinates are integers
    for col in ['start1', 'end1', 'start2', 'end2']:
        bedpe[col] = bedpe[col].astype(int)

    print(f"  Output columns: {', '.join(output_cols)}")
    print(f"  Tooltip will show: {', '.join(additional_cols)}")

    return bedpe


def write_bedpe(bedpe, output_path):
    """Write BEDPE dataframe to file."""
    print(f"\nWriting BEDPE to: {output_path}")

    # Create output directory if needed
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    # Write with header (so Juicebox can parse column names for tooltip)
    bedpe.to_csv(
        output_path,
        sep='\t',
        index=False,
        header=True,
        na_rep='.'
    )

    print(f"  Wrote {len(bedpe)} loops")
    print(f"  File size: {Path(output_path).stat().st_size / 1024:.1f} KB")

    # Print summary statistics if available
    if 'logFC' in bedpe.columns:
        print(f"\nlogFC statistics:")
        print(f"  Min: {bedpe['logFC'].min():.3f}")
        print(f"  Max: {bedpe['logFC'].max():.3f}")
        print(f"  Mean: {bedpe['logFC'].mean():.3f}")

    if 'FDR' in bedpe.columns:
        sig_count = (bedpe['FDR'] < 0.05).sum()
        print(f"\nFDR statistics:")
        print(f"  FDR < 0.05: {sig_count} loops ({sig_count/len(bedpe)*100:.1f}%)")


def main():
    """Main execution function."""
    args = parse_arguments()

    try:
        # Validate input file
        validate_input_file(args.input_tsv)

        # Load loops
        df = load_loops(args.input_tsv)

        # Apply filters
        df = filter_loops(df, args)

        if len(df) == 0:
            print("\nWARNING: No loops remaining after filtering!")
            print("Consider relaxing filter criteria.")
            sys.exit(1)

        # Convert to BEDPE
        bedpe = convert_to_bedpe(df, args)

        # Write output
        write_bedpe(bedpe, args.output_bedpe)

        print("\n✓ Conversion complete!")
        print(f"\nTo load in Juicebox:")
        print(f"  1. Open your .hic file in Juicebox")
        print(f"  2. Go to Annotations > Load Basic Annotations")
        print(f"  3. Select: {args.output_bedpe}")
        print(f"  4. Hover over loops to see tooltip with all columns")

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
