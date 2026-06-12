# non_cg_analysis/atlas_comparison.py
"""
Compare evoC non-CG methylation with Liu et al. 2023 mouse brain atlas.

Source: Supplementary Table 2 (MOESM5_ESM.csv) of:
  Liu, Zeng, Zhou et al. (2023) Nature 624, 366-377.
  doi:10.1038/s41586-023-06805-y

Usage:
    python atlas_comparison.py
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ATLAS_CSV = SCRIPT_DIR / "41586_2023_6805_MOESM5_ESM.csv"

PROJECT_ROOT = SCRIPT_DIR.parent.parent
EVOC_SUMMARY_CSV = (
    PROJECT_ROOT
    / "biomodal/upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv"
)

CEREBELLAR_KEYWORDS = ["cb ", "granule", "purk", "cerebel", "bergmann", "cbn "]


def load_atlas(path):
    stats = defaultdict(lambda: {"mCH": [], "mCG": [], "mCCC": []})
    n = 0
    with open(path) as f:
        for row in csv.DictReader(f):
            n += 1
            sub = row["Subclass"].strip()
            if not sub:
                continue
            stats[sub]["mCH"].append(float(row["mCHFrac"]))
            stats[sub]["mCG"].append(float(row["mCGFrac"]))
            stats[sub]["mCCC"].append(float(row["mCCCFrac"]))
    print(f"Atlas: {n} cells, {len(stats)} subclasses")
    return dict(stats)


def load_evoc(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append({
                "sample": row["sample_id"],
                "condition": "ctrl" if "ctrl" in row["sample_id"] else "mut",
                "CHH": float(row["modality_summary_chh_autosomes_modc"]),
                "CHG": float(row["modality_summary_chg_autosomes_modc"]),
                "CG": float(row["modality_summary_cg_autosomes_modc"]),
                "FP": (1.0 - float(row["dqs_puc19_modc_specificity"])) * 100.0,
            })
    return rows


def mean(vals):
    return sum(vals) / len(vals) if vals else 0.0


def linreg(x, y):
    n = len(x)
    xm, ym = mean(x), mean(y)
    ss_xx = sum((xi - xm) ** 2 for xi in x)
    ss_xy = sum((xi - xm) * (yi - ym) for xi, yi in zip(x, y))
    ss_yy = sum((yi - ym) ** 2 for yi in y)
    if ss_xx == 0 or ss_yy == 0:
        return 0.0, ym, 0.0, 0.0
    slope = ss_xy / ss_xx
    intercept = ym - slope * xm
    r = ss_xy / (ss_xx * ss_yy) ** 0.5
    return slope, intercept, r, r * r


def is_cerebellar(name):
    low = name.lower()
    return any(kw in low for kw in CEREBELLAR_KEYWORDS)


def print_section(title):
    print()
    print("=" * 90)
    print(title)
    print("=" * 90)


def plot_conversion_confound(evoc, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

    for idx, ctx in enumerate(["CHH", "CHG"]):
        ax = axes[idx]
        fp = [r["FP"] for r in evoc]
        meth = [r[ctx] for r in evoc]
        cond = [r["condition"] for r in evoc]

        slope, intercept, r, r2 = linreg(fp, meth)

        for i, row in enumerate(evoc):
            color = "#2166ac" if cond[i] == "ctrl" else "#b2182b"
            marker = "o" if cond[i] == "ctrl" else "s"
            label = row["sample"].replace("evoC-Bap1-", "")
            ax.scatter(fp[i], meth[i], c=color, marker=marker, s=60, zorder=3)
            ax.annotate(label, (fp[i], meth[i]), fontsize=6.5,
                        xytext=(4, 4), textcoords="offset points")

        x_line = np.linspace(min(fp) - 0.01, max(fp) + 0.01, 100)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, "k--", linewidth=1, alpha=0.7)

        eq = f"y = {slope:.3f}x + {intercept:.4f}"
        stats_text = f"r = {r:.3f}, R² = {r2:.3f}, n = {len(fp)}"
        ax.text(0.05, 0.95, f"{eq}\n{stats_text}",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="gray", alpha=0.9))

        ax.set_xlabel("pUC19 false-positive rate (%)", fontsize=10)
        ax.set_ylabel(f"{ctx} autosomal modC (%)", fontsize=10)
        ax.set_title(f"{ctx} modC vs conversion failure", fontsize=11,
                     fontweight="bold")

        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#2166ac",
                   markersize=8, label="control"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="#b2182b",
                   markersize=8, label="mutant"),
        ]
        ax.legend(handles=legend_elements, fontsize=8, loc="lower right")

    fig.suptitle("evoC non-CG methylation vs conversion efficiency (8 libraries)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    outpath = output_dir / "conversion_confound.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outpath}")


def main():
    if not ATLAS_CSV.exists():
        sys.exit(f"ERROR: {ATLAS_CSV}")

    stats = load_atlas(ATLAS_CSV)

    # --- Cerebellar subclasses ---
    print_section("CEREBELLAR SUBCLASSES (Liu et al. 2023, Supp Table 2)")
    print(f"{'Subclass':<35s} {'n':>6s} {'mCH%':>8s} {'mCG%':>8s} {'mCCC%':>8s}")
    print("-" * 70)
    cb = [(s, d) for s, d in stats.items() if is_cerebellar(s)]
    cb.sort(key=lambda x: mean(x[1]["mCH"]))
    for s, d in cb:
        print(f"{s:<35s} {len(d['mCH']):>6d} {100*mean(d['mCH']):>8.4f} "
              f"{100*mean(d['mCG']):>8.4f} {100*mean(d['mCCC']):>8.4f}")

    # --- Full ranking ---
    print_section("ALL 274 SUBCLASSES RANKED BY GLOBAL mCH")
    ranked = sorted(stats.items(), key=lambda x: mean(x[1]["mCH"]))
    print(f"{'#':>4s} {'Subclass':<42s} {'n':>6s} {'mCH%':>8s} {'mCG%':>8s} {'mCCC%':>8s}")
    print("-" * 80)
    for i, (s, d) in enumerate(ranked, 1):
        tag = ""
        if s == "CB Granule Glut":
            tag = " ***"
        elif s == "CBX Purkinje Gaba":
            tag = " **"
        elif s == "CBN Neurod2 Pvalb Glut":
            tag = " *"
        print(f"{i:>4d} {s:<42s} {len(d['mCH']):>6d} {100*mean(d['mCH']):>8.4f} "
              f"{100*mean(d['mCG']):>8.4f} {100*mean(d['mCCC']):>8.4f}{tag}")

    # --- evoC all 8 samples ---
    if not EVOC_SUMMARY_CSV.exists():
        print(f"WARNING: {EVOC_SUMMARY_CSV} not found")
        return

    evoc = load_evoc(EVOC_SUMMARY_CSV)

    print_section("evoC SAMPLES (sorted by FP rate)")
    print(f"{'Sample':<28s} {'Cond':>5s} {'CHH%':>8s} {'CHG%':>8s} "
          f"{'CG%':>8s} {'FP%':>8s}")
    print("-" * 70)
    evoc_sorted = sorted(evoc, key=lambda r: r["FP"])
    for r in evoc_sorted:
        print(f"{r['sample']:<28s} {r['condition']:>5s} {r['CHH']:>8.4f} "
              f"{r['CHG']:>8.4f} {r['CG']:>8.4f} {r['FP']:>8.4f}")

    ctrl = [r for r in evoc if r["condition"] == "ctrl"]
    mut = [r for r in evoc if r["condition"] == "mut"]
    print()
    print(f"  ctrl mean (n={len(ctrl)}):  CHH={mean([r['CHH'] for r in ctrl]):.4f}%  "
          f"CHG={mean([r['CHG'] for r in ctrl]):.4f}%  "
          f"FP={mean([r['FP'] for r in ctrl]):.4f}%")
    print(f"  mut  mean (n={len(mut)}):   CHH={mean([r['CHH'] for r in mut]):.4f}%  "
          f"CHG={mean([r['CHG'] for r in mut]):.4f}%  "
          f"FP={mean([r['FP'] for r in mut]):.4f}%")

    # --- Regression: apparent mCH vs FP ---
    print_section("LINEAR REGRESSION: apparent modC = slope * FP + intercept")
    all_fp = [r["FP"] for r in evoc]
    for ctx in ["CHH", "CHG"]:
        all_ctx = [r[ctx] for r in evoc]
        slope, intercept, r, r2 = linreg(all_fp, all_ctx)
        print(f"  {ctx}:  slope={slope:.4f}  intercept={intercept:.4f}  r={r:.4f}  R²={r2:.4f}  n={len(evoc)}")

    fp_gap = mean([r["FP"] for r in mut]) - mean([r["FP"] for r in ctrl])
    chh_gap = mean([r["CHH"] for r in mut]) - mean([r["CHH"] for r in ctrl])
    chg_gap = mean([r["CHG"] for r in mut]) - mean([r["CHG"] for r in ctrl])
    chh_slope = linreg(all_fp, [r["CHH"] for r in evoc])[0]
    chg_slope = linreg(all_fp, [r["CHG"] for r in evoc])[0]

    print()
    print(f"  Genotype deltas (mut - ctrl):")
    print(f"    dFP  = {fp_gap:+.4f}%")
    print(f"    dCHH = {chh_gap:+.4f}%    slope * dFP = {chh_slope * fp_gap:+.4f}%")
    print(f"    dCHG = {chg_gap:+.4f}%    slope * dFP = {chg_slope * fp_gap:+.4f}%")

    # --- Floor correction: atlas vs evoC ---
    cb_mCH = 100 * mean(stats["CB Granule Glut"]["mCH"])
    cb_mCG = 100 * mean(stats["CB Granule Glut"]["mCG"])
    cb_mCCC = 100 * mean(stats["CB Granule Glut"]["mCCC"])
    evoc_chh = mean([r["CHH"] for r in ctrl])
    evoc_chg = mean([r["CHG"] for r in ctrl])
    evoc_cg = mean([r["CG"] for r in ctrl])
    evoc_fp = mean([r["FP"] for r in ctrl])

    print_section("FLOOR CORRECTION: atlas mCCC vs evoC pUC19 FP")
    print(f"  {'':>25s} {'Atlas CB Gran':>14s} {'evoC ctrl':>14s}")
    print(f"  {'':>25s} {'(bisulfite)':>14s} {'(enzymatic)':>14s}")
    print(f"  {'-'*55}")
    print(f"  {'Raw mCH / CHH':>25s} {cb_mCH:>13.4f}% {evoc_chh:>13.4f}%")
    print(f"  {'FP floor':>25s} {cb_mCCC:>13.4f}% {evoc_fp:>13.4f}%")
    print(f"  {'Corrected':>25s} {max(cb_mCH - cb_mCCC, 0):>13.4f}% {max(evoc_chh - evoc_fp, 0):>13.4f}%")
    print()
    print(f"  {'Raw mCG / CG':>25s} {cb_mCG:>13.4f}% {evoc_cg:>13.4f}%")

    # --- Plot ---
    print()
    plot_conversion_confound(evoc, SCRIPT_DIR)


if __name__ == "__main__":
    main()
