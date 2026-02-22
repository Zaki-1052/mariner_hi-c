# non_cg_analysis/non_cg_methylation.py
"""
Compare non-CG methylation (CHH, CHG, CA-context) between EM-seq and evoC.

Genome-wide aggregate rates AND per-gene comparisons.
EM-seq data from MethylDackel methylKit files; evoC from duet summary + DMR BEDs.

Usage:
    python non_cg_methylation.py [--ca-dir PATH]  # optional CA-filtered data from cluster
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths (relative to project root)
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent

METHYLKIT_DIR = PROJECT_ROOT / "em-seq_output" / "methylDackelExtracts"
METHYLKIT_FILES = {
    "ctrl1_CHH": METHYLKIT_DIR / "250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHH.methylKit.gz",
    "ctrl1_CHG": METHYLKIT_DIR / "250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHG.methylKit.gz",
    "ctrl2_CHH": METHYLKIT_DIR / "250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHH.methylKit.gz",
    "ctrl2_CHG": METHYLKIT_DIR / "250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHG.methylKit.gz",
}

EVOC_SUMMARY_CSV = (
    PROJECT_ROOT
    / "biomodal/upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv"
)

EVOC_DMR_BEDS = {
    "CHH": (
        PROJECT_ROOT
        / "biomodal/downstream/modality/outputs/run-2/outputs_CHH/Results"
        / "gencode.vM25.mouse.genes.annotation/DMR_20260121_184551"
        / "DMR_mc_control__mutant_20260121_184551.bed"
    ),
    "CHG": (
        PROJECT_ROOT
        / "biomodal/downstream/modality/outputs/run-2/outputs_CHG/Results"
        / "gencode.vM25.mouse.genes.annotation/DMR_20260121_174629"
        / "DMR_mc_control__mutant_20260121_174629.bed"
    ),
}

GENE_BED = (
    PROJECT_ROOT
    / "biomodal/downstream/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz"
)

OUTPUT_DIR = Path(__file__).resolve().parent

AUTOSOMES = {f"chr{i}" for i in range(1, 20)}
SEX_CHROMS = {"chrX", "chrY"}
SPIKE_INS = {"phage_lambda", "plasmid_puc19c", "phage_T4", "phage_Xp12"}


# ---------------------------------------------------------------------------
# Step 2: Parse EM-seq methylKit files → genome-wide non-CG rates
# ---------------------------------------------------------------------------

def categorize_chrom(chrom: str) -> str:
    """Assign chromosome to a reporting category."""
    if chrom in AUTOSOMES:
        return "autosomes"
    if chrom in SEX_CHROMS:
        return "sex"
    if chrom in SPIKE_INS:
        return chrom
    return "other"


def parse_methylkit(filepath: Path) -> dict:
    """Stream a gzipped methylKit file and accumulate per-category stats.

    Returns dict[category] -> {meth_reads, total_cov, n_sites, per_chr: {chr -> {meth_reads, total_cov, n_sites}}}
    """
    stats = defaultdict(lambda: {"meth_reads": 0.0, "total_cov": 0, "n_sites": 0,
                                  "per_chr": defaultdict(lambda: {"meth_reads": 0.0, "total_cov": 0, "n_sites": 0})})

    with gzip.open(filepath, "rt") as fh:
        header = fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[1]
            coverage = int(parts[4])
            freq_c = float(parts[5])

            meth_reads = coverage * freq_c / 100.0
            cat = categorize_chrom(chrom)

            stats[cat]["meth_reads"] += meth_reads
            stats[cat]["total_cov"] += coverage
            stats[cat]["n_sites"] += 1
            stats[cat]["per_chr"][chrom]["meth_reads"] += meth_reads
            stats[cat]["per_chr"][chrom]["total_cov"] += coverage
            stats[cat]["per_chr"][chrom]["n_sites"] += 1

    return dict(stats)


def compute_rate(d: dict) -> float:
    """Weighted methylation rate (%) from accumulated stats."""
    if d["total_cov"] == 0:
        return 0.0
    return d["meth_reads"] / d["total_cov"] * 100.0


def run_emseq_genomewide() -> pd.DataFrame:
    """Parse all EM-seq methylKit files and return a summary DataFrame."""
    rows = []
    per_chr_data = {}  # for per-chromosome breakdown

    for label, fpath in METHYLKIT_FILES.items():
        sample, context = label.rsplit("_", 1)
        print(f"  Parsing {fpath.name} ...", flush=True)
        stats = parse_methylkit(fpath)

        auto = stats.get("autosomes", {"meth_reads": 0, "total_cov": 0, "n_sites": 0, "per_chr": {}})
        sex = stats.get("sex", {"meth_reads": 0, "total_cov": 0, "n_sites": 0, "per_chr": {}})

        rows.append({
            "method": "EM-seq",
            "sample": sample,
            "context": context,
            "autosomes_rate_pct": compute_rate(auto),
            "autosomes_sites": auto["n_sites"],
            "autosomes_cov": auto["total_cov"],
            "sex_rate_pct": compute_rate(sex),
            "lambda_rate_pct": compute_rate(stats.get("phage_lambda", {"meth_reads": 0, "total_cov": 0, "n_sites": 0})),
            "lambda_sites": stats.get("phage_lambda", {"n_sites": 0})["n_sites"],
            "puc19_rate_pct": compute_rate(stats.get("plasmid_puc19c", {"meth_reads": 0, "total_cov": 0, "n_sites": 0})),
            "puc19_sites": stats.get("plasmid_puc19c", {"n_sites": 0})["n_sites"],
        })

        # Store per-chromosome data for autosome breakdown
        key = (sample, context)
        per_chr_data[key] = auto.get("per_chr", {})

    return pd.DataFrame(rows), per_chr_data


# ---------------------------------------------------------------------------
# Step 3: evoC baseline from summary CSV and DMR BEDs
# ---------------------------------------------------------------------------

def extract_evoc_genomewide() -> pd.DataFrame:
    """Extract evoC genome-wide non-CG rates from the duet summary CSV."""
    rows = []
    with open(EVOC_SUMMARY_CSV) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sid = row["sample_id"]
            # Only control samples for fair comparison
            if "ctrl" not in sid:
                continue
            for ctx in ["chh", "chg"]:
                modc = float(row[f"modality_summary_{ctx}_autosomes_modc"])
                rows.append({
                    "method": "evoC",
                    "sample": sid,
                    "context": ctx.upper(),
                    "autosomes_rate_pct": modc,
                    "autosomes_sites": np.nan,
                    "autosomes_cov": np.nan,
                    "sex_rate_pct": np.nan,
                    "lambda_mc_sens": float(row["dqs_lambda_mc_sensitivity"]),
                    "puc19_mc_spec": float(row["dqs_puc19_mc_specificity"]),
                    "lambda_rate_pct": np.nan,
                    "lambda_sites": np.nan,
                    "puc19_rate_pct": np.nan,
                    "puc19_sites": np.nan,
                })
    return pd.DataFrame(rows)


def load_evoc_per_gene(context: str) -> pd.DataFrame:
    """Load evoC per-gene methylation from DMR BED (mean_mod_group_1 = control mean)."""
    bed_path = EVOC_DMR_BEDS[context]
    df = pd.read_csv(bed_path, sep="\t")
    # mean_mod_group_1 is the control group mean (already a fraction like 0.00086)
    # Convert to percentage for consistency
    df["evoc_rate_pct"] = df["mean_mod_group_1"] * 100.0
    df["evoc_coverage"] = df["mean_coverage"]
    df["evoc_num_contexts"] = df["num_contexts"]
    return df[["Chromosome", "Start", "End", "Name", "evoc_rate_pct", "evoc_coverage", "evoc_num_contexts"]]


# ---------------------------------------------------------------------------
# Step 5d: Per-gene EM-seq methylation
# ---------------------------------------------------------------------------

def load_gene_bed() -> pd.DataFrame:
    """Load gencode gene annotation BED. Columns: Chromosome, Start, End, Annotation, Name.
    Note: Chromosome is numeric (1, 2, ...) without 'chr' prefix.
    """
    df = pd.read_csv(GENE_BED, sep="\t", compression="gzip")
    return df


def compute_emseq_per_gene(gene_bed: pd.DataFrame, context: str, min_cov: int = 1) -> pd.DataFrame:
    """Aggregate EM-seq methylation per gene body for a given context.

    Averages ctrl1 and ctrl2 per-gene weighted methylation.
    Uses interval-based lookup: for each methylKit site, find overlapping gene(s).

    Strategy: Build a per-chromosome sorted gene array, binary search for overlaps.
    """
    # Gene BED uses '1', '2', ... but methylKit uses 'chr1', 'chr2', ...
    genes_by_chr = defaultdict(list)
    for _, row in gene_bed.iterrows():
        chrom_num = str(row["Chromosome"])
        chrom_key = f"chr{chrom_num}"
        genes_by_chr[chrom_key].append({
            "start": int(row["Start"]),
            "end": int(row["End"]),
            "name": row["Name"],
        })

    # Sort genes by start position for binary search
    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort(key=lambda g: g["start"])

    # Accumulate per-gene stats from both control samples
    gene_stats = defaultdict(lambda: {"meth_reads": 0.0, "total_cov": 0, "n_sites": 0})

    for sample in ["ctrl1", "ctrl2"]:
        label = f"{sample}_{context}"
        fpath = METHYLKIT_FILES[label]
        print(f"    Per-gene: streaming {fpath.name} ...", flush=True)

        with gzip.open(fpath, "rt") as fh:
            fh.readline()  # skip header
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                chrom = parts[1]
                pos = int(parts[2])  # 1-based
                coverage = int(parts[4])
                freq_c = float(parts[5])

                if coverage < min_cov:
                    continue
                if chrom not in genes_by_chr:
                    continue

                # Find overlapping genes (pos is 1-based; gene coords are 0-based BED)
                for gene in genes_by_chr[chrom]:
                    if pos < gene["start"]:
                        break  # sorted, no more overlaps
                    if gene["start"] <= pos <= gene["end"]:
                        gene_stats[gene["name"]]["meth_reads"] += coverage * freq_c / 100.0
                        gene_stats[gene["name"]]["total_cov"] += coverage
                        gene_stats[gene["name"]]["n_sites"] += 1

    # Convert to DataFrame
    rows = []
    for gene_name, st in gene_stats.items():
        if st["total_cov"] > 0:
            rows.append({
                "Name": gene_name,
                "emseq_rate_pct": st["meth_reads"] / st["total_cov"] * 100.0,
                "emseq_sites": st["n_sites"],
                "emseq_total_cov": st["total_cov"],
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5c: Coverage-stratified analysis
# ---------------------------------------------------------------------------

def coverage_stratified_rates(context: str, thresholds: list[int] = None) -> pd.DataFrame:
    """Compute autosome methylation at different coverage thresholds."""
    if thresholds is None:
        thresholds = [1, 2, 3, 5]

    rows = []
    for sample in ["ctrl1", "ctrl2"]:
        label = f"{sample}_{context}"
        fpath = METHYLKIT_FILES[label]
        print(f"    Cov-stratified: streaming {fpath.name} ...", flush=True)

        # Accumulate per-threshold
        accum = {t: {"meth_reads": 0.0, "total_cov": 0, "n_sites": 0} for t in thresholds}

        with gzip.open(fpath, "rt") as fh:
            fh.readline()
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                chrom = parts[1]
                if chrom not in AUTOSOMES:
                    continue
                coverage = int(parts[4])
                freq_c = float(parts[5])
                meth_reads = coverage * freq_c / 100.0

                for t in thresholds:
                    if coverage >= t:
                        accum[t]["meth_reads"] += meth_reads
                        accum[t]["total_cov"] += coverage
                        accum[t]["n_sites"] += 1

        for t in thresholds:
            rate = compute_rate(accum[t])
            rows.append({
                "sample": sample,
                "context": context,
                "min_cov": t,
                "rate_pct": rate,
                "n_sites": accum[t]["n_sites"],
                "total_cov": accum[t]["total_cov"],
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CA-filtered data (optional, from cluster output)
# ---------------------------------------------------------------------------

def load_ca_summary(ca_dir: Path) -> pd.DataFrame | None:
    """Load CA-filtered summary stats produced by ca_filter.py on the cluster."""
    summary_path = ca_dir / "ca_summary_stats.tsv"
    if not summary_path.exists():
        return None
    return pd.read_csv(summary_path, sep="\t")


def load_existing_summary(output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load previously saved non_cg_summary.tsv, split into EM-seq and evoC DataFrames."""
    path = output_dir / "non_cg_summary.tsv"
    if not path.exists():
        sys.exit(f"ERROR: Existing summary not found: {path}\nRun the full pipeline first before using --ca-only.")
    df = pd.read_csv(path, sep="\t")
    return df[df["method"] == "EM-seq"].copy(), df[df["method"] == "evoC"].copy()


# Filename patterns → (sample, context)
_FILENAME_MAP = {
    "250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHH.methylKit.gz": ("ctrl1", "CHH"),
    "250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHG.methylKit.gz": ("ctrl1", "CHG"),
    "250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHH.methylKit.gz": ("ctrl2", "CHH"),
    "250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHG.methylKit.gz": ("ctrl2", "CHG"),
}


def extract_ca_autosome_rates(ca_df: pd.DataFrame) -> pd.DataFrame:
    """Extract CA-context autosome rates from ca_summary_stats.tsv.

    Returns DataFrame with columns: sample, context, ca_rate_pct, total_rate_pct,
    ca_sites, total_sites, ca_cov, total_cov, ca_fraction.
    """
    # Filter to autosomes aggregate rows only
    auto_all = ca_df[(ca_df["category"] == "autosomes") & (ca_df["chromosome"] == "ALL")]

    rows = []
    for _, row in auto_all.iterrows():
        mapping = _FILENAME_MAP.get(row["file"])
        if mapping is None:
            continue
        sample, context = mapping
        rows.append({
            "sample": sample,
            "context": context,
            "ca_rate_pct": row["ca_rate_pct"],
            "total_rate_pct": row["total_rate_pct"],
            "ca_sites": int(row["ca_sites"]),
            "total_sites": int(row["total_sites"]),
            "ca_cov": int(row["ca_cov"]),
            "total_cov": int(row["total_cov"]),
            "ca_fraction": row["ca_fraction"],
        })

    return pd.DataFrame(rows)


def extract_ca_lambda_rates(ca_df: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Extract CA-context lambda FP rates from ca_summary_stats.tsv.

    Returns dict mapping (sample, context) -> CA-specific lambda false-positive rate (%).
    """
    lam_all = ca_df[(ca_df["category"] == "phage_lambda") & (ca_df["chromosome"] == "ALL")]

    rates = {}
    for _, row in lam_all.iterrows():
        mapping = _FILENAME_MAP.get(row["file"])
        if mapping is None:
            continue
        rates[mapping] = row["ca_rate_pct"]

    return rates


# ---------------------------------------------------------------------------
# Step 7: Visualization
# ---------------------------------------------------------------------------

def plot_genomewide_comparison(emseq_df: pd.DataFrame, evoc_df: pd.DataFrame,
                                ca_rates_df: pd.DataFrame | None, output_dir: Path):
    """Bar chart: CHH, CHG (and CA-filtered) methylation by method/sample."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for idx, context in enumerate(["CHH", "CHG"]):
        ax = axes[idx]

        labels = []
        rates = []
        colors = []

        # EM-seq bars
        em_data = emseq_df[emseq_df["context"] == context]
        for _, r in em_data.iterrows():
            labels.append(f"EM-seq {r['sample']}")
            rates.append(r["autosomes_rate_pct"])
            colors.append("#4C72B0")

        # CA-filtered bars (if available)
        if ca_rates_df is not None and not ca_rates_df.empty:
            ca_ctx = ca_rates_df[ca_rates_df["context"] == context]
            for _, r in ca_ctx.iterrows():
                labels.append(f"CA-only {r['sample']}")
                rates.append(r["ca_rate_pct"])
                colors.append("#55A868")

        # evoC bars
        ev_data = evoc_df[evoc_df["context"] == context]
        for _, r in ev_data.iterrows():
            labels.append(f"evoC {r['sample'].replace('evoC-Bap1-', '')}")
            rates.append(r["autosomes_rate_pct"])
            colors.append("#DD8452")

        rates = np.array(rates)
        x = np.arange(len(labels))
        bars = ax.bar(x, rates, color=colors, edgecolor="black", linewidth=0.5)

        for bar, rate in zip(bars, rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                    f"{rate:.3f}%", ha="center", va="bottom", fontsize=8)

        ax.set_title(f"{context} Methylation (Autosomes)", fontsize=13, fontweight="bold")
        ax.set_ylabel("Methylation (%)" if idx == 0 else "")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
        ax.set_ylim(0, max(rates) * 1.5 if max(rates) > 0 else 2)

    title = "Non-CG Methylation: EM-seq vs evoC (Control Samples)"
    if ca_rates_df is not None and not ca_rates_df.empty:
        title += "\n(includes CA-context sub-rates)"
    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    outpath = output_dir / "non_cg_comparison.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_per_chromosome(per_chr_data: dict, output_dir: Path):
    """Per-autosome non-CG rates for EM-seq control samples."""
    chroms = [f"chr{i}" for i in range(1, 20)]
    contexts = ["CHH", "CHG"]
    samples = ["ctrl1", "ctrl2"]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    for ctx_idx, context in enumerate(contexts):
        ax = axes[ctx_idx]
        x = np.arange(len(chroms))
        width = 0.35

        for s_idx, sample in enumerate(samples):
            key = (sample, context)
            chr_dict = per_chr_data.get(key, {})
            rates = []
            for c in chroms:
                d = chr_dict.get(c, {"meth_reads": 0, "total_cov": 0})
                rates.append(compute_rate(d) if d["total_cov"] > 0 else 0)
            offset = (s_idx - 0.5) * width
            ax.bar(x + offset, rates, width, label=sample, edgecolor="black", linewidth=0.3)

        ax.set_ylabel(f"{context} Methylation (%)")
        ax.set_title(f"Per-Autosome {context} Methylation (EM-seq)", fontweight="bold")
        ax.legend()

    axes[-1].set_xticks(np.arange(len(chroms)))
    axes[-1].set_xticklabels(chroms, rotation=45, ha="right")
    fig.tight_layout()
    outpath = output_dir / "per_chromosome_rates.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_per_gene_scatter(merged_df: pd.DataFrame, context: str, output_dir: Path,
                           min_emseq_sites: int = 3):
    """Scatter: EM-seq per-gene % vs evoC per-gene % for a given context."""
    df = merged_df[merged_df["emseq_sites"] >= min_emseq_sites].copy()
    if df.empty:
        print(f"  [WARN] No genes with >= {min_emseq_sites} EM-seq {context} sites; skipping scatter plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df["evoc_rate_pct"], df["emseq_rate_pct"], alpha=0.3, s=10, c="#4C72B0")

    max_val = max(df["evoc_rate_pct"].max(), df["emseq_rate_pct"].max(), 1)
    ax.plot([0, max_val * 1.1], [0, max_val * 1.1], "k--", alpha=0.4, linewidth=0.8)

    ax.set_xlabel(f"evoC {context} Methylation (%, control mean)")
    ax.set_ylabel(f"EM-seq {context} Methylation (%, ctrl1+ctrl2 weighted)")
    ax.set_title(
        f"Per-Gene {context} Methylation: EM-seq vs evoC\n"
        f"(n={len(df)} genes with >= {min_emseq_sites} EM-seq sites)",
        fontweight="bold",
    )
    fig.tight_layout()
    outpath = output_dir / f"per_gene_scatter_{context}.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_ca_comparison(emseq_df: pd.DataFrame, ca_rates_df: pd.DataFrame, output_dir: Path):
    """Grouped bar chart: original vs CA-only methylation rate per sample/context."""
    contexts = ["CHH", "CHG"]
    samples = ["ctrl1", "ctrl2"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    for idx, context in enumerate(contexts):
        ax = axes[idx]
        x = np.arange(len(samples))
        width = 0.3

        # Original rates from emseq_df
        orig_rates = []
        for sample in samples:
            row = emseq_df[(emseq_df["sample"] == sample) & (emseq_df["context"] == context)]
            orig_rates.append(row["autosomes_rate_pct"].values[0] if len(row) > 0 else 0)

        # CA-only rates
        ca_rates = []
        for sample in samples:
            row = ca_rates_df[(ca_rates_df["sample"] == sample) & (ca_rates_df["context"] == context)]
            ca_rates.append(row["ca_rate_pct"].values[0] if len(row) > 0 else 0)

        bars1 = ax.bar(x - width / 2, orig_rates, width, label=f"All {context}", color="#4C72B0",
                       edgecolor="black", linewidth=0.5)
        bars2 = ax.bar(x + width / 2, ca_rates, width, label=f"CA-only {context}", color="#55A868",
                       edgecolor="black", linewidth=0.5)

        for bar, rate in zip(bars1, orig_rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"{rate:.3f}%", ha="center", va="bottom", fontsize=9)
        for bar, rate in zip(bars2, ca_rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"{rate:.3f}%", ha="center", va="bottom", fontsize=9)

        ax.set_title(f"{context} Context (Autosomes)", fontsize=13, fontweight="bold")
        ax.set_ylabel("Methylation (%)" if idx == 0 else "")
        ax.set_xticks(x)
        ax.set_xticklabels(samples, fontsize=11)
        ax.legend(fontsize=9)
        all_rates = orig_rates + ca_rates
        ax.set_ylim(0, max(all_rates) * 1.5 if max(all_rates) > 0 else 2)

    fig.suptitle("EM-seq Non-CG Methylation: All Sites vs CA-Context Only",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    outpath = output_dir / "ca_vs_total_comparison.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def plot_corrected_comparison(emseq_df: pd.DataFrame, ca_rates_df: pd.DataFrame,
                              ca_lambda_rates: dict[tuple[str, str], float],
                              evoc_df: pd.DataFrame, output_dir: Path):
    """Bar chart comparing FP-corrected EM-seq rates (all-context and CA-only) vs evoC."""
    contexts = ["CHH", "CHG"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

    for idx, context in enumerate(contexts):
        ax = axes[idx]
        labels = []
        rates = []
        colors = []
        hatches = []

        samples = ["ctrl1", "ctrl2"]
        for sample in samples:
            # Raw all-context rate
            em_row = emseq_df[(emseq_df["sample"] == sample) & (emseq_df["context"] == context)]
            raw_all = em_row["autosomes_rate_pct"].values[0] if len(em_row) > 0 else 0
            lam_all = em_row["lambda_rate_pct"].values[0] if len(em_row) > 0 else 0

            # Raw CA-only rate
            ca_row = ca_rates_df[(ca_rates_df["sample"] == sample) & (ca_rates_df["context"] == context)]
            raw_ca = ca_row["ca_rate_pct"].values[0] if len(ca_row) > 0 else 0
            lam_ca = ca_lambda_rates.get((sample, context), 0)

            corrected_all = max(raw_all - lam_all, 0)
            corrected_ca = max(raw_ca - lam_ca, 0)

            # Corrected all-context
            labels.append(f"{sample}\nall {context}")
            rates.append(corrected_all)
            colors.append("#4C72B0")
            hatches.append("")

            # Corrected CA-only
            labels.append(f"{sample}\nCA-only")
            rates.append(corrected_ca)
            colors.append("#55A868")
            hatches.append("")

        # evoC rates
        ev_data = evoc_df[evoc_df["context"] == context]
        for _, r in ev_data.iterrows():
            short_name = r["sample"].replace("evoC-Bap1-", "")
            labels.append(f"{short_name}\nevoC")
            rates.append(r["autosomes_rate_pct"])
            colors.append("#DD8452")
            hatches.append("")

        rates_arr = np.array(rates)
        x = np.arange(len(labels))
        bars = ax.bar(x, rates_arr, color=colors, edgecolor="black", linewidth=0.5)

        for bar, rate in zip(bars, rates_arr):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"{rate:.3f}%", ha="center", va="bottom", fontsize=8)

        ax.set_title(f"{context} Methylation (Autosomes)", fontsize=13, fontweight="bold")
        ax.set_ylabel("Methylation (%)" if idx == 0 else "")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, max(rates_arr) * 1.6 if max(rates_arr) > 0 else 2)

    fig.suptitle("FP-Corrected Non-CG Methylation: EM-seq vs evoC\n"
                 "(EM-seq rates corrected by subtracting lambda FP rate)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    outpath = output_dir / "corrected_comparison.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def print_ca_summary_table(emseq_df: pd.DataFrame, ca_rates_df: pd.DataFrame,
                           ca_lambda_rates: dict[tuple[str, str], float] | None = None):
    """Print original vs CA-filtered rates side-by-side, with optional FP correction."""
    has_correction = ca_lambda_rates is not None and len(ca_lambda_rates) > 0
    width = 130 if has_correction else 100

    print("\n" + "=" * width)
    print("CA-CONTEXT METHYLATION COMPARISON (EM-seq, Autosomes)")
    print("=" * width)

    header = (f"{'Sample':<8} {'Context':<6} {'All Rate(%)':<13} {'CA Rate(%)':<13} "
              f"{'CA/All Ratio':<13} {'CA Sites':<14} {'CA Fraction':<13}")
    if has_correction:
        header += f"{'Lambda FP(%)':<13} {'Corrected(%)':<13}"
    print(header)
    print("-" * width)

    for _, ca_row in ca_rates_df.iterrows():
        sample = ca_row["sample"]
        context = ca_row["context"]
        em_row = emseq_df[(emseq_df["sample"] == sample) & (emseq_df["context"] == context)]
        orig_rate = em_row["autosomes_rate_pct"].values[0] if len(em_row) > 0 else 0
        ca_rate = ca_row["ca_rate_pct"]
        ratio = ca_rate / orig_rate if orig_rate > 0 else 0
        line = (f"{sample:<8} {context:<6} {orig_rate:<13.4f} {ca_rate:<13.4f} "
                f"{ratio:<13.2f} {ca_row['ca_sites']:>13,} {ca_row['ca_fraction']:<13.4f}")
        if has_correction:
            lam_fp = ca_lambda_rates.get((sample, context), 0)
            corrected = max(ca_rate - lam_fp, 0)
            line += f"{lam_fp:<13.4f} {corrected:<13.4f}"
        print(line)

    print("=" * width)


# ---------------------------------------------------------------------------
# Step 7: Summary output
# ---------------------------------------------------------------------------

def print_summary_table(emseq_df: pd.DataFrame, evoc_df: pd.DataFrame):
    """Print a formatted comparison table to console."""
    print("\n" + "=" * 90)
    print("GENOME-WIDE NON-CG METHYLATION COMPARISON")
    print("=" * 90)

    header = f"{'Method':<10} {'Sample':<22} {'Context':<6} {'Rate(%)':<10} {'Sites':<14} {'Coverage':<14}"
    print(header)
    print("-" * 90)

    for _, r in emseq_df.iterrows():
        sites = f"{r['autosomes_sites']:,.0f}" if not np.isnan(r.get("autosomes_sites", np.nan)) else "N/A"
        cov = f"{r['autosomes_cov']:,.0f}" if not np.isnan(r.get("autosomes_cov", np.nan)) else "N/A"
        print(f"{'EM-seq':<10} {r['sample']:<22} {r['context']:<6} {r['autosomes_rate_pct']:<10.4f} {sites:<14} {cov:<14}")

    print("-" * 90)
    for _, r in evoc_df.iterrows():
        print(f"{'evoC':<10} {r['sample']:<22} {r['context']:<6} {r['autosomes_rate_pct']:<10.4f} {'N/A':<14} {'N/A':<14}")

    print("\n" + "=" * 90)
    print("SPIKE-IN FALSE-POSITIVE RATES (EM-seq)")
    print("=" * 90)
    header2 = f"{'Sample':<10} {'Context':<6} {'Lambda(%)':<12} {'Lambda Sites':<14} {'pUC19(%)':<12} {'pUC19 Sites':<14}"
    print(header2)
    print("-" * 90)
    for _, r in emseq_df.iterrows():
        lam_sites = f"{r['lambda_sites']:,.0f}" if not np.isnan(r.get("lambda_sites", np.nan)) else "N/A"
        puc_sites = f"{r['puc19_sites']:,.0f}" if not np.isnan(r.get("puc19_sites", np.nan)) else "N/A"
        print(f"{r['sample']:<10} {r['context']:<6} {r['lambda_rate_pct']:<12.4f} {lam_sites:<14} {r['puc19_rate_pct']:<12.4f} {puc_sites:<14}")

    print("\n" + "=" * 90)
    print("evoC QC METRICS (Control Samples)")
    print("=" * 90)
    for _, r in evoc_df.iterrows():
        if "lambda_mc_sens" in r and not np.isnan(r.get("lambda_mc_sens", np.nan)):
            print(f"  {r['sample']}: lambda_mc_sensitivity={r['lambda_mc_sens']:.4f}, puc19_mc_specificity={r['puc19_mc_spec']:.4f}")


def save_summary_tsv(emseq_df: pd.DataFrame, evoc_df: pd.DataFrame, output_dir: Path):
    """Save machine-readable summary."""
    combined = pd.concat([emseq_df, evoc_df], ignore_index=True)
    outpath = output_dir / "non_cg_summary.tsv"
    combined.to_csv(outpath, sep="\t", index=False)
    print(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Non-CG methylation comparison: EM-seq vs evoC")
    parser.add_argument("--ca-dir", type=Path, default=None,
                        help="Directory with CA-filtered output from ca_filter.py (optional)")
    parser.add_argument("--ca-only", action="store_true",
                        help="Only run CA analysis: load existing non_cg_summary.tsv, integrate CA data, regenerate plots")
    parser.add_argument("--skip-per-gene", action="store_true",
                        help="Skip per-gene analysis (slow with 300M-line files)")
    parser.add_argument("--skip-cov-stratified", action="store_true",
                        help="Skip coverage-stratified analysis (requires extra pass)")
    parser.add_argument("--min-gene-sites", type=int, default=3,
                        help="Minimum EM-seq sites per gene for scatter plot (default: 3)")
    args = parser.parse_args()

    # --ca-only: fast path that skips all heavy methylKit parsing
    if args.ca_only:
        if args.ca_dir is None:
            sys.exit("ERROR: --ca-only requires --ca-dir to be set.")
        print("[CA-only mode] Loading existing data and integrating CA rates ...\n")

        print("[1/4] Loading existing non_cg_summary.tsv ...")
        emseq_df, evoc_df = load_existing_summary(OUTPUT_DIR)
        print(f"  EM-seq rows: {len(emseq_df)}, evoC rows: {len(evoc_df)}")

        print("\n[2/4] Loading CA-filtered summary ...")
        ca_raw = load_ca_summary(args.ca_dir)
        if ca_raw is None:
            sys.exit(f"ERROR: No ca_summary_stats.tsv found in {args.ca_dir}")
        ca_rates_df = extract_ca_autosome_rates(ca_raw)
        print(f"  CA autosome rate rows: {len(ca_rates_df)}")

        ca_lambda_rates = extract_ca_lambda_rates(ca_raw)
        print(f"  CA lambda FP rates: {ca_lambda_rates}")

        print("\n[3/4] CA comparison table ...")
        print_ca_summary_table(emseq_df, ca_rates_df, ca_lambda_rates)

        print("\n[4/4] Generating plots ...")
        plot_genomewide_comparison(emseq_df, evoc_df, ca_rates_df, OUTPUT_DIR)
        plot_ca_comparison(emseq_df, ca_rates_df, OUTPUT_DIR)
        plot_corrected_comparison(emseq_df, ca_rates_df, ca_lambda_rates, evoc_df, OUTPUT_DIR)

        print("\nDone (CA-only mode).")
        return

    # Verify input files exist
    print("[1/7] Verifying input files ...")
    for label, fpath in METHYLKIT_FILES.items():
        if not fpath.exists():
            sys.exit(f"ERROR: Missing methylKit file: {fpath}")
    if not EVOC_SUMMARY_CSV.exists():
        sys.exit(f"ERROR: Missing evoC summary: {EVOC_SUMMARY_CSV}")
    for ctx, bed in EVOC_DMR_BEDS.items():
        if not bed.exists():
            sys.exit(f"ERROR: Missing evoC DMR BED ({ctx}): {bed}")
    if not GENE_BED.exists():
        sys.exit(f"ERROR: Missing gene BED: {GENE_BED}")
    print("  All input files found.")

    # Step 2: EM-seq genome-wide rates
    print("\n[2/7] Computing EM-seq genome-wide non-CG rates ...")
    emseq_df, per_chr_data = run_emseq_genomewide()

    # Step 3: evoC baseline
    print("\n[3/7] Extracting evoC baseline ...")
    evoc_df = extract_evoc_genomewide()

    # Step 5a: Genome-wide comparison table
    print("\n[4/7] Genome-wide comparison ...")
    print_summary_table(emseq_df, evoc_df)

    # Step 5b: Per-chromosome breakdown
    print("\n[5/7] Per-chromosome breakdown ...")
    plot_per_chromosome(per_chr_data, OUTPUT_DIR)

    # Step 5c: Coverage-stratified analysis
    if not args.skip_cov_stratified:
        print("\n[5c/7] Coverage-stratified analysis ...")
        for context in ["CHH", "CHG"]:
            cov_df = coverage_stratified_rates(context)
            print(f"\n  {context} coverage-stratified rates (autosomes):")
            print(f"  {'Sample':<8} {'MinCov':<8} {'Rate(%)':<10} {'Sites':<14} {'TotalCov':<14}")
            for _, r in cov_df.iterrows():
                print(f"  {r['sample']:<8} {r['min_cov']:<8} {r['rate_pct']:<10.4f} {r['n_sites']:>13,.0f} {r['total_cov']:>13,.0f}")
    else:
        print("\n[5c/7] Skipping coverage-stratified analysis (--skip-cov-stratified)")

    # Step 5d: Per-gene comparison
    if not args.skip_per_gene:
        print("\n[6/7] Per-gene comparison ...")
        gene_bed = load_gene_bed()

        for context in ["CHH", "CHG"]:
            print(f"\n  --- {context} per-gene ---")
            emseq_genes = compute_emseq_per_gene(gene_bed, context)
            evoc_genes = load_evoc_per_gene(context)

            if emseq_genes.empty:
                print(f"  [WARN] No EM-seq per-gene data for {context}; skipping.")
                continue

            merged = pd.merge(emseq_genes, evoc_genes, on="Name", how="inner")
            print(f"  Merged genes: {len(merged)} (EM-seq: {len(emseq_genes)}, evoC: {len(evoc_genes)})")

            if not merged.empty:
                plot_per_gene_scatter(merged, context, OUTPUT_DIR, min_emseq_sites=args.min_gene_sites)
    else:
        print("\n[6/7] Skipping per-gene analysis (--skip-per-gene)")

    # CA-filtered data (if available)
    ca_rates_df = None
    ca_lambda_rates = None
    if args.ca_dir:
        print(f"\n  Loading CA-filtered data from {args.ca_dir} ...")
        ca_raw = load_ca_summary(args.ca_dir)
        if ca_raw is not None:
            ca_rates_df = extract_ca_autosome_rates(ca_raw)
            ca_lambda_rates = extract_ca_lambda_rates(ca_raw)
            print_ca_summary_table(emseq_df, ca_rates_df, ca_lambda_rates)
        else:
            print("  [WARN] No ca_summary_stats.tsv found in --ca-dir")

    # Step 7: Plots and output files
    print("\n[7/7] Generating plots and saving outputs ...")
    plot_genomewide_comparison(emseq_df, evoc_df, ca_rates_df, OUTPUT_DIR)
    if ca_rates_df is not None:
        plot_ca_comparison(emseq_df, ca_rates_df, OUTPUT_DIR)
    if ca_rates_df is not None and ca_lambda_rates is not None:
        plot_corrected_comparison(emseq_df, ca_rates_df, ca_lambda_rates, evoc_df, OUTPUT_DIR)
    save_summary_tsv(emseq_df, evoc_df, OUTPUT_DIR)

    print("\nDone.")


if __name__ == "__main__":
    main()
