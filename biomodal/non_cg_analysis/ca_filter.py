# non_cg_analysis/ca_filter.py
"""
CA-context filtering for EM-seq methylKit files.

Reads CHH and CHG methylKit files, looks up each position in the reference
genome to determine trinucleotide context, writes CA-only sites to new files,
and produces summary statistics.

Requires pysam and a FASTA reference with .fai index.

Usage:
    python ca_filter.py \
        --ref /path/to/mm10+controls.fa \
        --input-dir /path/to/methylKit/files \
        --output-dir ./ca_filtered/
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pysam

AUTOSOMES = {f"chr{i}" for i in range(1, 20)}
SEX_CHROMS = {"chrX", "chrY"}
SPIKE_INS = {"phage_lambda", "plasmid_puc19c", "phage_T4", "phage_Xp12"}

METHYLKIT_BASENAMES = [
    "250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHH.methylKit.gz",
    "250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHG.methylKit.gz",
    "250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHH.methylKit.gz",
    "250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHG.methylKit.gz",
]


def categorize_chrom(chrom: str) -> str:
    if chrom in AUTOSOMES:
        return "autosomes"
    if chrom in SEX_CHROMS:
        return "sex"
    if chrom in SPIKE_INS:
        return chrom
    return "other"


def is_ca_context(fasta: pysam.FastaFile, chrom: str, pos_1based: int, strand: str) -> bool:
    """Check if a cytosine site is in CA context.

    methylKit positions are 1-based. pysam.fetch uses 0-based half-open intervals.

    For forward strand (F): C at pos, next base at pos+1
        -> fetch(chrom, pos, pos+1) to get the base immediately after C
        -> CA if next_base == 'A'

    For reverse strand (R): G on forward strand at pos (complement of C)
        -> The C on reverse strand 'sees' the base to its 3' side on the reverse strand,
           which is the base to the 5' side on the forward strand (pos - 1 on forward, 1-based)
        -> fetch(chrom, pos-2, pos-1) to get the forward-strand base at position pos-1
        -> CA on reverse strand if that forward base == 'T' (complement of A)
    """
    try:
        if strand == "F":
            # Next base after C on forward strand (0-based: pos to pos+1)
            next_base = fasta.fetch(chrom, pos_1based, pos_1based + 1).upper()
            return next_base == "A"
        elif strand == "R":
            # Base before G on forward strand = base at pos-1 (1-based)
            # 0-based fetch: pos_1based - 2 to pos_1based - 1
            if pos_1based < 2:
                return False
            prev_base = fasta.fetch(chrom, pos_1based - 2, pos_1based - 1).upper()
            return prev_base == "T"
    except (ValueError, KeyError):
        return False
    return False


def process_file(fasta: pysam.FastaFile, input_path: Path, output_path: Path) -> dict:
    """Filter a single methylKit file for CA-context sites.

    Returns summary stats dict.
    """
    stats = {
        "total_sites": 0,
        "ca_sites": 0,
        "categories": defaultdict(lambda: {
            "total_sites": 0, "ca_sites": 0,
            "total_cov": 0, "ca_cov": 0,
            "total_meth_reads": 0.0, "ca_meth_reads": 0.0,
            "per_chr": defaultdict(lambda: {
                "total_sites": 0, "ca_sites": 0,
                "total_cov": 0, "ca_cov": 0,
                "total_meth_reads": 0.0, "ca_meth_reads": 0.0,
            }),
        }),
    }

    with gzip.open(input_path, "rt") as fin, gzip.open(output_path, "wt") as fout:
        header = fin.readline()
        fout.write(header)

        for line_num, line in enumerate(fin, start=1):
            parts = line.rstrip("\n").split("\t")
            chrom = parts[1]
            pos = int(parts[2])
            strand = parts[3]
            coverage = int(parts[4])
            freq_c = float(parts[5])
            meth_reads = coverage * freq_c / 100.0

            cat = categorize_chrom(chrom)
            stats["total_sites"] += 1
            stats["categories"][cat]["total_sites"] += 1
            stats["categories"][cat]["total_cov"] += coverage
            stats["categories"][cat]["total_meth_reads"] += meth_reads
            stats["categories"][cat]["per_chr"][chrom]["total_sites"] += 1
            stats["categories"][cat]["per_chr"][chrom]["total_cov"] += coverage
            stats["categories"][cat]["per_chr"][chrom]["total_meth_reads"] += meth_reads

            ca = is_ca_context(fasta, chrom, pos, strand)
            if ca:
                fout.write(line)
                stats["ca_sites"] += 1
                stats["categories"][cat]["ca_sites"] += 1
                stats["categories"][cat]["ca_cov"] += coverage
                stats["categories"][cat]["ca_meth_reads"] += meth_reads
                stats["categories"][cat]["per_chr"][chrom]["ca_sites"] += 1
                stats["categories"][cat]["per_chr"][chrom]["ca_cov"] += coverage
                stats["categories"][cat]["per_chr"][chrom]["ca_meth_reads"] += meth_reads

            if line_num % 10_000_000 == 0:
                pct_ca = stats["ca_sites"] / stats["total_sites"] * 100 if stats["total_sites"] > 0 else 0
                print(f"    {line_num:>12,} lines processed, {stats['ca_sites']:>10,} CA sites ({pct_ca:.1f}%)",
                      flush=True)

    return stats


def write_summary(all_stats: dict, output_dir: Path):
    """Write ca_summary_stats.tsv with per-file, per-category breakdown."""
    outpath = output_dir / "ca_summary_stats.tsv"
    with open(outpath, "w") as fh:
        cols = [
            "file", "category", "chromosome",
            "total_sites", "ca_sites", "ca_fraction",
            "total_cov", "ca_cov",
            "total_rate_pct", "ca_rate_pct",
        ]
        fh.write("\t".join(cols) + "\n")

        for filename, stats in sorted(all_stats.items()):
            for cat_name, cat_data in sorted(stats["categories"].items()):
                # Aggregate row for this category
                total_rate = (cat_data["total_meth_reads"] / cat_data["total_cov"] * 100
                              if cat_data["total_cov"] > 0 else 0)
                ca_rate = (cat_data["ca_meth_reads"] / cat_data["ca_cov"] * 100
                           if cat_data["ca_cov"] > 0 else 0)
                ca_frac = (cat_data["ca_sites"] / cat_data["total_sites"]
                           if cat_data["total_sites"] > 0 else 0)

                fh.write("\t".join([
                    filename, cat_name, "ALL",
                    str(cat_data["total_sites"]), str(cat_data["ca_sites"]), f"{ca_frac:.6f}",
                    str(cat_data["total_cov"]), str(cat_data["ca_cov"]),
                    f"{total_rate:.6f}", f"{ca_rate:.6f}",
                ]) + "\n")

                # Per-chromosome rows
                for chrom, chr_data in sorted(cat_data["per_chr"].items()):
                    chr_total_rate = (chr_data["total_meth_reads"] / chr_data["total_cov"] * 100
                                     if chr_data["total_cov"] > 0 else 0)
                    chr_ca_rate = (chr_data["ca_meth_reads"] / chr_data["ca_cov"] * 100
                                  if chr_data["ca_cov"] > 0 else 0)
                    chr_ca_frac = (chr_data["ca_sites"] / chr_data["total_sites"]
                                   if chr_data["total_sites"] > 0 else 0)

                    fh.write("\t".join([
                        filename, cat_name, chrom,
                        str(chr_data["total_sites"]), str(chr_data["ca_sites"]), f"{chr_ca_frac:.6f}",
                        str(chr_data["total_cov"]), str(chr_data["ca_cov"]),
                        f"{chr_total_rate:.6f}", f"{chr_ca_rate:.6f}",
                    ]) + "\n")

    print(f"  Summary written to {outpath}")


def main():
    parser = argparse.ArgumentParser(description="Filter methylKit files for CA-context sites")
    parser.add_argument("--ref", type=Path, required=True,
                        help="Reference FASTA (must have .fai index)")
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Directory containing methylKit.gz files")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for CA-filtered files and summary")
    args = parser.parse_args()

    if not args.ref.exists():
        sys.exit(f"ERROR: Reference not found: {args.ref}")
    fai_path = Path(str(args.ref) + ".fai")
    if not fai_path.exists():
        sys.exit(f"ERROR: Reference index not found: {fai_path} (run samtools faidx)")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    fasta = pysam.FastaFile(str(args.ref))
    all_stats = {}

    for basename in METHYLKIT_BASENAMES:
        input_path = args.input_dir / basename
        if not input_path.exists():
            print(f"[WARN] Skipping missing file: {input_path}")
            continue

        # Output filename: replace CHH/CHG with CA_CHH/CA_CHG
        out_basename = basename.replace("_CHH.", "_CA_CHH.").replace("_CHG.", "_CA_CHG.")
        output_path = args.output_dir / out_basename

        print(f"\nProcessing {basename} ...")
        stats = process_file(fasta, input_path, output_path)
        all_stats[basename] = stats

        ca_pct = stats["ca_sites"] / stats["total_sites"] * 100 if stats["total_sites"] > 0 else 0
        print(f"  Total sites: {stats['total_sites']:,}, CA sites: {stats['ca_sites']:,} ({ca_pct:.1f}%)")

    fasta.close()

    print("\nWriting summary ...")
    write_summary(all_stats, args.output_dir)
    print("Done.")


if __name__ == "__main__":
    main()
