#!/usr/bin/env python3
# abc/scripts/step7_rnaseq_integration.py
"""
Integrate ΔABC with RNA-seq log2FC, aggregate to gene level,
and identify dysregulated genes.

Combines original Steps 7 and 8:
  - Merge ΔABC pairs with RNA-seq DE results
  - Gene-level summary (strongest ΔABC per gene + aggregate stats)
  - Concordance analysis (gene-level to avoid pseudo-replication)
  - Dysregulated gene identification and top-gene reporting
  - Scatter plots (normalized and unnormalized)

Usage: python step7_rnaseq_integration.py
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# === CONFIGURATION ===
DELTA_ABC = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/abc/results/delta_abc_all_pairs.tsv"
RNASEQ = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
OUTPUT_DIR = "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/abc/results"
GENE_COL = "TargetGene"

# Thresholds
ABC_DELTA_THRESH = 0.01
RNA_PADJ = 0.05
RNA_LFC = 0.5

os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)


def load_data():
    """Load delta-ABC and RNA-seq, report overlap."""
    delta = pd.read_csv(DELTA_ABC, sep="\t")
    rna = pd.read_excel(RNASEQ)
    rna = rna.rename(columns={"ensembl_gene_id": "gene_symbol"})

    print(f"ΔABC pairs: {len(delta):,}")
    print(f"ΔABC unique genes: {delta[GENE_COL].nunique():,}")
    print(f"RNA-seq genes: {len(rna):,}")
    return delta, rna


def merge_with_rnaseq(delta, rna):
    """Inner join ΔABC pairs with RNA-seq on gene symbol.

    Renames log2FoldChange -> log2FC consistently across both outputs.
    """
    rna_cols = ["gene_symbol", "log2FoldChange", "padj", "baseMean"]
    merged = pd.merge(
        delta, rna[rna_cols],
        left_on=GENE_COL, right_on="gene_symbol", how="inner",
    )
    merged.drop(columns=["gene_symbol"], inplace=True)
    merged.rename(columns={"log2FoldChange": "log2FC"}, inplace=True)

    n_abc = delta[GENE_COL].nunique()
    n_matched = merged[GENE_COL].nunique()
    print(f"\nGenes in ABC: {n_abc:,}")
    print(f"Genes in RNA-seq: {rna['gene_symbol'].nunique():,}")
    print(f"Genes matched: {n_matched:,} ({100*n_matched/n_abc:.1f}% of ABC genes)")
    print(f"E-G pairs with RNA-seq: {len(merged):,}")
    return merged


def build_gene_summary(merged):
    """Aggregate to gene level: strongest ΔABC enhancer + per-gene stats.

    Carries through strongest-enhancer metadata (distance, class, Hi-C flag)
    and adds aggregate metrics (mean delta, enhancer counts, sum scores).
    """
    # Row with max |delta_ABC| per gene
    idx_max = merged.groupby(GENE_COL)["delta_ABC"].apply(
        lambda x: x.abs().idxmax()
    )
    strongest = merged.loc[idx_max].copy()
    strongest = strongest.rename(columns={
        "delta_ABC": "max_delta_abc",
        "delta_unnorm": "max_delta_unnorm",
        "chr": "top_enh_chr",
        "start": "top_enh_start",
        "end": "top_enh_end",
        "distance": "top_enh_distance",
        "class": "top_enh_class",
    })

    # Per-gene aggregate stats from all enhancers
    agg = merged.groupby(GENE_COL).agg(
        n_enhancers=("delta_ABC", "size"),
        mean_delta_abc=("delta_ABC", "mean"),
        sum_delta_abc=("delta_ABC", "sum"),
        sum_delta_unnorm=("delta_unnorm", "sum"),
        n_gained=("delta_ABC", lambda x: (x > ABC_DELTA_THRESH).sum()),
        n_lost=("delta_ABC", lambda x: (x < -ABC_DELTA_THRESH).sum()),
        sum_abc_wt=("ABC.Score_WT", "sum"),
        sum_abc_ko=("ABC.Score_KO", "sum"),
    ).reset_index()

    gene_summary = pd.merge(strongest, agg, on=GENE_COL, how="left")
    print(f"Gene-level summary: {len(gene_summary):,} genes")
    return gene_summary


def concordance_analysis(gene_summary):
    """Gene-level concordance between ΔABC direction and log2FC direction."""
    gs = gene_summary.copy()
    gs["abc_changed"] = gs["max_delta_abc"].abs() > ABC_DELTA_THRESH
    gs["rna_changed"] = (gs["padj"] < RNA_PADJ) & (gs["log2FC"].abs() > RNA_LFC)

    both = gs[gs["abc_changed"] & gs["rna_changed"]].copy()

    print(f"\n=== Gene-Level Concordance ===")
    print(f"  Genes with |ΔABC| > {ABC_DELTA_THRESH}: {gs['abc_changed'].sum():,}")
    print(f"  Genes with DE (padj<{RNA_PADJ}, |LFC|>{RNA_LFC}): {gs['rna_changed'].sum():,}")
    print(f"  Genes with both: {len(both):,}")

    if len(both) == 0:
        print("  WARNING: No genes with both changed. Try relaxing thresholds.")
        return

    both["concordant"] = np.sign(both["max_delta_abc"]) == np.sign(both["log2FC"])
    n_conc = both["concordant"].sum()
    n_disc = (~both["concordant"]).sum()
    n_tot = len(both)
    binom_p = stats.binomtest(n_conc, n_tot, 0.5).pvalue

    print(f"  Concordant (same direction): {n_conc} ({100*n_conc/n_tot:.1f}%)")
    print(f"  Discordant: {n_disc} ({100*n_disc/n_tot:.1f}%)")
    print(f"  Binomial test p-value (vs 50%): {binom_p:.2e}")

    gained_up = ((both["max_delta_abc"] > 0) & (both["log2FC"] > 0)).sum()
    lost_down = ((both["max_delta_abc"] < 0) & (both["log2FC"] < 0)).sum()
    gained_down = ((both["max_delta_abc"] > 0) & (both["log2FC"] < 0)).sum()
    lost_up = ((both["max_delta_abc"] < 0) & (both["log2FC"] > 0)).sum()
    print(f"\n  Gained enhancer + upregulated:   {gained_up}")
    print(f"  Lost enhancer + downregulated:   {lost_down}")
    print(f"  Gained enhancer + downregulated: {gained_down}")
    print(f"  Lost enhancer + upregulated:     {lost_up}")


def identify_dysregulated(gene_summary):
    """Flag dysregulated genes and report top hits."""
    gs = gene_summary
    gs["dysregulated"] = (
        (gs["max_delta_abc"].abs() > ABC_DELTA_THRESH)
        & (gs["padj"] < RNA_PADJ)
        & (gs["log2FC"].abs() > RNA_LFC)
    )

    n_dys = gs["dysregulated"].sum()
    dys = gs[gs["dysregulated"]]
    n_conc = (np.sign(dys["max_delta_abc"]) == np.sign(dys["log2FC"])).sum()

    print(f"\n=== Dysregulated Genes ===")
    print(f"  Total: {n_dys} (|ΔABC|>{ABC_DELTA_THRESH} AND padj<{RNA_PADJ} AND |LFC|>{RNA_LFC})")
    print(f"  Concordant: {n_conc}, Discordant: {n_dys - n_conc}")

    if n_dys > 0:
        top_cols = [GENE_COL, "max_delta_abc", "log2FC", "padj",
                    "n_enhancers", "n_gained", "n_lost",
                    "top_enh_distance", "top_enh_class"]
        # Filter to columns that exist
        top_cols = [c for c in top_cols if c in dys.columns]
        top = dys.reindex(dys["max_delta_abc"].abs().sort_values(ascending=False).index).head(20)[top_cols]
        print(f"\n  Top 20 by |ΔABC|:")
        print(top.to_string(index=False))

    return gs


def make_scatter_plots(gene_summary):
    """Two-panel scatter: normalized ΔABC and unnormalized Δ(A×C) vs log2FC."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sig = gene_summary["padj"] < RNA_PADJ
    ns = ~sig

    def pct_lim(series, pad=0.1):
        lo, hi = series.quantile(0.005), series.quantile(0.995)
        rng = hi - lo
        return lo - pad * rng, hi + pad * rng

    # --- Panel A: Normalized ΔABC ---
    ax = axes[0]
    ax.scatter(gene_summary.loc[ns, "max_delta_abc"],
               gene_summary.loc[ns, "log2FC"],
               alpha=0.3, s=10, c="gray", label=f"padj ≥ {RNA_PADJ}", rasterized=True)
    ax.scatter(gene_summary.loc[sig, "max_delta_abc"],
               gene_summary.loc[sig, "log2FC"],
               alpha=0.5, s=15, c="crimson", label=f"padj < {RNA_PADJ}", rasterized=True)
    ax.axhline(0, color="black", lw=0.5, ls="--")
    ax.axvline(0, color="black", lw=0.5, ls="--")
    ax.set_xlabel("ΔABC (strongest enhancer per gene)")
    ax.set_ylabel("log2FC (KO vs WT)")
    ax.set_title("Normalized ABC Score")
    ax.legend(fontsize=8)

    de = gene_summary[sig & gene_summary["max_delta_abc"].notna()
                       & gene_summary["log2FC"].notna()]
    if len(de) > 2:
        r, p = stats.pearsonr(de["max_delta_abc"], de["log2FC"])
        rho, sp = stats.spearmanr(de["max_delta_abc"], de["log2FC"])
        ax.text(0.05, 0.95,
                f"Pearson r={r:.3f}, p={p:.2e}\nSpearman ρ={rho:.3f}, p={sp:.2e}\n(DE genes, n={len(de)})",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    # --- Panel B: Unnormalized Δ(A×C) ---
    ax = axes[1]
    xlim = pct_lim(gene_summary["max_delta_unnorm"].dropna())
    ax.scatter(gene_summary.loc[ns, "max_delta_unnorm"],
               gene_summary.loc[ns, "log2FC"],
               alpha=0.3, s=10, c="gray", label=f"padj ≥ {RNA_PADJ}", rasterized=True)
    ax.scatter(gene_summary.loc[sig, "max_delta_unnorm"],
               gene_summary.loc[sig, "log2FC"],
               alpha=0.5, s=15, c="steelblue", label=f"padj < {RNA_PADJ}", rasterized=True)
    ax.axhline(0, color="black", lw=0.5, ls="--")
    ax.axvline(0, color="black", lw=0.5, ls="--")
    ax.set_xlim(xlim)
    ax.set_xlabel("Δ(Activity × Contact) (strongest enhancer per gene)")
    ax.set_ylabel("log2FC (KO vs WT)")
    ax.set_title("Unnormalized Score")
    ax.legend(fontsize=8)

    if len(de) > 2:
        r2, p2 = stats.pearsonr(de["max_delta_unnorm"], de["log2FC"])
        rho2, sp2 = stats.spearmanr(de["max_delta_unnorm"], de["log2FC"])
        ax.text(0.05, 0.95,
                f"Pearson r={r2:.3f}, p={p2:.2e}\nSpearman ρ={rho2:.3f}, p={sp2:.2e}\n(DE genes, n={len(de)})",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        fig.savefig(f"{OUTPUT_DIR}/figures/delta_abc_vs_log2fc.{ext}",
                    dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {OUTPUT_DIR}/figures/delta_abc_vs_log2fc.pdf")


def make_sum_scatter_plots(gene_summary):
    """Two-panel scatter: sum of all enhancers' ΔABC and Δ(A×C) vs log2FC."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sig = gene_summary["padj"] < RNA_PADJ
    ns = ~sig

    def pct_lim(series, pad=0.1):
        lo, hi = series.quantile(0.005), series.quantile(0.995)
        rng = hi - lo
        return lo - pad * rng, hi + pad * rng

    # --- Panel A: Normalized sum ΔABC ---
    ax = axes[0]
    ax.scatter(gene_summary.loc[ns, "sum_delta_abc"],
               gene_summary.loc[ns, "log2FC"],
               alpha=0.3, s=10, c="gray", label=f"padj \u2265 {RNA_PADJ}", rasterized=True)
    ax.scatter(gene_summary.loc[sig, "sum_delta_abc"],
               gene_summary.loc[sig, "log2FC"],
               alpha=0.5, s=15, c="crimson", label=f"padj < {RNA_PADJ}", rasterized=True)
    ax.axhline(0, color="black", lw=0.5, ls="--")
    ax.axvline(0, color="black", lw=0.5, ls="--")
    xlim_a = pct_lim(gene_summary["sum_delta_abc"].dropna())
    ax.set_xlim(xlim_a)
    ax.set_xlabel("\u03a3 \u0394ABC (all enhancers per gene)")
    ax.set_ylabel("log2FC (KO vs WT)")
    ax.set_title("Normalized ABC Score (Sum)")
    ax.legend(fontsize=8)

    de = gene_summary[sig & gene_summary["sum_delta_abc"].notna()
                       & gene_summary["log2FC"].notna()]
    if len(de) > 2:
        r, p = stats.pearsonr(de["sum_delta_abc"], de["log2FC"])
        rho, sp = stats.spearmanr(de["sum_delta_abc"], de["log2FC"])
        ax.text(0.05, 0.95,
                f"Pearson r={r:.3f}, p={p:.2e}\nSpearman \u03c1={rho:.3f}, p={sp:.2e}\n(DE genes, n={len(de)})",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    # --- Panel B: Unnormalized sum Δ(A×C) ---
    ax = axes[1]
    xlim_b = pct_lim(gene_summary["sum_delta_unnorm"].dropna())
    ax.scatter(gene_summary.loc[ns, "sum_delta_unnorm"],
               gene_summary.loc[ns, "log2FC"],
               alpha=0.3, s=10, c="gray", label=f"padj \u2265 {RNA_PADJ}", rasterized=True)
    ax.scatter(gene_summary.loc[sig, "sum_delta_unnorm"],
               gene_summary.loc[sig, "log2FC"],
               alpha=0.5, s=15, c="steelblue", label=f"padj < {RNA_PADJ}", rasterized=True)
    ax.axhline(0, color="black", lw=0.5, ls="--")
    ax.axvline(0, color="black", lw=0.5, ls="--")
    ax.set_xlim(xlim_b)
    ax.set_xlabel("\u03a3 \u0394(Activity \u00d7 Contact) (all enhancers per gene)")
    ax.set_ylabel("log2FC (KO vs WT)")
    ax.set_title("Unnormalized Score (Sum)")
    ax.legend(fontsize=8)

    de2 = gene_summary[sig & gene_summary["sum_delta_unnorm"].notna()
                        & gene_summary["log2FC"].notna()]
    if len(de2) > 2:
        r2, p2 = stats.pearsonr(de2["sum_delta_unnorm"], de2["log2FC"])
        rho2, sp2 = stats.spearmanr(de2["sum_delta_unnorm"], de2["log2FC"])
        ax.text(0.05, 0.95,
                f"Pearson r={r2:.3f}, p={p2:.2e}\nSpearman \u03c1={rho2:.3f}, p={sp2:.2e}\n(DE genes, n={len(de2)})",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()
    for ext in ["pdf", "png"]:
        fig.savefig(f"{OUTPUT_DIR}/figures/sum_delta_abc_vs_log2fc.{ext}",
                    dpi=300 if ext == "pdf" else 150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR}/figures/sum_delta_abc_vs_log2fc.pdf")


def print_correlation_comparison(gene_summary):
    """Print strongest-link vs sum correlations side by side."""
    sig = gene_summary["padj"] < RNA_PADJ
    de = gene_summary[sig].dropna(subset=["log2FC", "max_delta_abc",
                                           "max_delta_unnorm", "sum_delta_abc",
                                           "sum_delta_unnorm"])
    if len(de) < 3:
        print("WARNING: Too few DE genes for correlation comparison.")
        return

    metrics = [
        ("Normalized (strongest)",   "max_delta_abc"),
        ("Normalized (sum)",         "sum_delta_abc"),
        ("Unnormalized (strongest)", "max_delta_unnorm"),
        ("Unnormalized (sum)",       "sum_delta_unnorm"),
    ]

    print(f"\n{'='*65}")
    print(f"  Correlation Comparison: Strongest Enhancer vs Sum (n={len(de)} DE genes)")
    print(f"{'='*65}")
    print(f"  {'Metric':<28s} {'Pearson r':>10s} {'Spearman ρ':>12s}")
    print(f"  {'-'*28} {'-'*10} {'-'*12}")
    for label, col in metrics:
        r, _ = stats.pearsonr(de[col], de["log2FC"])
        rho, _ = stats.spearmanr(de[col], de["log2FC"])
        print(f"  {label:<28s} {r:>10.3f} {rho:>12.3f}")
    print(f"{'='*65}")


def main():
    delta, rna = load_data()
    merged = merge_with_rnaseq(delta, rna)
    gene_summary = build_gene_summary(merged)
    concordance_analysis(gene_summary)
    gene_summary = identify_dysregulated(gene_summary)
    make_scatter_plots(gene_summary)
    make_sum_scatter_plots(gene_summary)
    print_correlation_comparison(gene_summary)

    # Save
    merged.to_csv(f"{OUTPUT_DIR}/delta_abc_with_rnaseq.tsv",
                   sep="\t", index=False, float_format="%.6f")
    gene_summary.to_csv(f"{OUTPUT_DIR}/gene_level_summary.tsv",
                         sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved: {OUTPUT_DIR}/delta_abc_with_rnaseq.tsv ({len(merged):,} pairs)")
    print(f"Saved: {OUTPUT_DIR}/gene_level_summary.tsv ({len(gene_summary):,} genes)")


if __name__ == "__main__":
    main()
