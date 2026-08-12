#!/usr/bin/env python3
# abc/scripts/step6_delta_abc.py
"""
Compute ΔABC and Δ(Activity × Contact) between KO and WT conditions.

Column names verified from Step 5 QC output (2026-02-17).
Self-promoter entries (forced ABC=1.0) are excluded before delta computation.

Usage: python step6_delta_abc.py
"""

import pandas as pd
import numpy as np
import sys
import os

# === CONFIGURATION (verified from Step 5 QC) ===
ABC_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction"
OUTPUT_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

WT_PRED = f"{ABC_DIR}/results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz"
KO_PRED = f"{ABC_DIR}/results/KO_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz"

# Column names (confirmed from QC Section 2)
ENHANCER_COLS = ["chr", "start", "end"]
GENE_COL = "TargetGene"
ABC_SCORE_COL = "ABC.Score"
UNNORM_COL = "ABC.Score.Numerator"
SELF_PROMOTER_COL = "isSelfPromoter"
DISTANCE_COL = "distance"
CLASS_COL = "class"
ACTIVITY_COL = "activity_base"
HIC_CONTACT_COL = "hic_contact"
HIC_ADJ_COL = "hic_contact_pl_scaled_adj"

JOIN_COLS = ENHANCER_COLS + [GENE_COL]

# Columns shared between conditions (same enhancer-gene pair = same values)
# Coalesce from WT after join to avoid duplication
SHARED_COLS = [DISTANCE_COL, CLASS_COL]

# Per-condition columns to carry through the join
PERCOND_COLS = [ABC_SCORE_COL, UNNORM_COL, ACTIVITY_COL, HIC_ADJ_COL, HIC_CONTACT_COL]

ABC_THRESHOLD = 0.02
DELTA_SIG_THRESHOLD = 0.01


def load_and_filter(path, label):
    """Load AllPutative predictions, verify columns, remove self-promoters."""
    print(f"Loading {label} predictions: {path}")
    df = pd.read_csv(path, sep="\t")
    print(f"  Raw: {len(df):,} E-G pairs, {df[GENE_COL].nunique():,} genes")

    # Verify required columns exist
    required = JOIN_COLS + [ABC_SCORE_COL, UNNORM_COL, SELF_PROMOTER_COL]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: Missing columns in {label}: {missing}")
        print(f"  Available: {df.columns.tolist()}")
        sys.exit(1)

    # Filter out self-promoters (forced ABC=1.0, would create delta≈0 noise)
    n_before = len(df)
    n_sp = (df[SELF_PROMOTER_COL] == True).sum()
    df = df[df[SELF_PROMOTER_COL] != True].copy()
    print(f"  Removed {n_sp:,} self-promoter entries ({100*n_sp/n_before:.1f}%)")
    print(f"  Distal: {len(df):,} E-G pairs, {df[GENE_COL].nunique():,} genes")

    return df


def main():
    wt = load_and_filter(WT_PRED, "WT")
    ko = load_and_filter(KO_PRED, "KO")

    # === SELECT COLUMNS FOR JOIN ===
    # Shared cols in both so coalesce has a fallback for KO-only pairs
    wt_cols = JOIN_COLS + SHARED_COLS + PERCOND_COLS
    ko_cols = JOIN_COLS + SHARED_COLS + PERCOND_COLS
    # Only select columns that actually exist
    wt_slim = wt[[c for c in wt_cols if c in wt.columns]].copy()
    ko_slim = ko[[c for c in ko_cols if c in ko.columns]].copy()

    # === OUTER JOIN ===
    print("\nJoining WT and KO on enhancer-gene pairs...")
    merged = pd.merge(
        wt_slim, ko_slim,
        on=JOIN_COLS,
        how="outer",
        suffixes=("_WT", "_KO"),
        indicator=True,
    )

    merge_counts = merged["_merge"].value_counts()
    n_both = merge_counts.get("both", 0)
    n_wt_only = merge_counts.get("left_only", 0)
    n_ko_only = merge_counts.get("right_only", 0)
    print(f"  Both conditions: {n_both:,}")
    print(f"  WT only:         {n_wt_only:,}")
    print(f"  KO only:         {n_ko_only:,}")
    print(f"  Overlap: {100*n_both/len(merged):.1f}%")

    if 100 * n_both / len(merged) < 90:
        print("WARNING: <90% overlap — check if candidate regions differ between conditions.")

    # === FILL MISSING SCORES WITH 0 ===
    score_cols_wt = [f"{c}_WT" for c in PERCOND_COLS if f"{c}_WT" in merged.columns]
    score_cols_ko = [f"{c}_KO" for c in PERCOND_COLS if f"{c}_KO" in merged.columns]
    for col in score_cols_wt + score_cols_ko:
        merged[col] = merged[col].fillna(0)

    # Coalesce shared cols: prefer WT value, fall back to KO for KO-only pairs
    for col in SHARED_COLS:
        if f"{col}_WT" in merged.columns and f"{col}_KO" in merged.columns:
            # merge created suffixed versions since col was in both slims
            merged[col] = merged[f"{col}_WT"].fillna(merged[f"{col}_KO"])
            merged.drop(columns=[f"{col}_WT", f"{col}_KO"], inplace=True)
        # If col only came from wt_slim (not in ko_slim), it's already named correctly

    # === COMPUTE DELTAS ===
    abc_wt = f"{ABC_SCORE_COL}_WT"
    abc_ko = f"{ABC_SCORE_COL}_KO"
    unnorm_wt = f"{UNNORM_COL}_WT"
    unnorm_ko = f"{UNNORM_COL}_KO"

    merged["delta_ABC"] = merged[abc_ko] - merged[abc_wt]
    merged["delta_unnorm"] = merged[unnorm_ko] - merged[unnorm_wt]

    # === FLAG Hi-C SPARSITY ===
    # QC showed 3.1-3.4% of thresholded distal pairs lack real Hi-C contact.
    # Flag pairs where BOTH conditions have zero raw Hi-C (powerlaw-only).
    hic_wt = f"{HIC_CONTACT_COL}_WT"
    hic_ko = f"{HIC_CONTACT_COL}_KO"
    if hic_wt in merged.columns and hic_ko in merged.columns:
        merged["hic_data_wt"] = merged[hic_wt] > 0
        merged["hic_data_ko"] = merged[hic_ko] > 0
        merged["has_hic_either"] = merged["hic_data_wt"] | merged["hic_data_ko"]

    # === APPLY THRESHOLD ===
    # Keep pairs where EITHER condition has ABC >= threshold
    mask = (merged[abc_wt] >= ABC_THRESHOLD) | (merged[abc_ko] >= ABC_THRESHOLD)
    filtered = merged[mask].copy()
    filtered.drop(columns=["_merge"], inplace=True)

    print(f"\nTotal distal E-G pairs (union, AllPutative): {len(merged):,}")
    print(f"Pairs with ABC >= {ABC_THRESHOLD} in at least one condition: {len(filtered):,}")
    print(f"Unique genes: {filtered[GENE_COL].nunique():,}")

    # === SUMMARY STATISTICS ===
    print("\n=== ΔABC Summary (thresholded distal pairs) ===")
    d = filtered["delta_ABC"]
    print(f"  n:      {len(d):,}")
    print(f"  Mean:   {d.mean():.6f}")
    print(f"  Median: {d.median():.6f}")
    print(f"  Std:    {d.std():.6f}")
    print(f"  Q5/Q95: {d.quantile(0.05):.5f} / {d.quantile(0.95):.5f}")

    n_gained = (d > DELTA_SIG_THRESHOLD).sum()
    n_lost = (d < -DELTA_SIG_THRESHOLD).sum()
    n_unchanged = (d.abs() <= DELTA_SIG_THRESHOLD).sum()
    print(f"  Gained (ΔABC > {DELTA_SIG_THRESHOLD}):    {n_gained:,}")
    print(f"  Lost (ΔABC < -{DELTA_SIG_THRESHOLD}):     {n_lost:,}")
    print(f"  Unchanged (|Δ| <= {DELTA_SIG_THRESHOLD}): {n_unchanged:,}")

    print("\n=== Δ(Activity × Contact) Summary ===")
    du = filtered["delta_unnorm"]
    print(f"  Mean:   {du.mean():.6f}")
    print(f"  Median: {du.median():.6f}")
    print(f"  Std:    {du.std():.6f}")

    # Hi-C data coverage in thresholded pairs
    if "has_hic_either" in filtered.columns:
        n_hic = filtered["has_hic_either"].sum()
        print(f"\n=== Hi-C Coverage (thresholded pairs) ===")
        print(f"  With Hi-C in at least one condition: {n_hic:,}/{len(filtered):,} ({100*n_hic/len(filtered):.1f}%)")
        print(f"  Powerlaw-only (no Hi-C either cond): {len(filtered)-n_hic:,}")

    # Asymmetry check
    if abs(d.mean()) > 0.005:
        print(f"\n⚠️  ΔABC mean is notably non-zero ({d.mean():.6f}).")
        print("   This may indicate global activity shift in KO vs WT.")
        print("   The unnormalized Δ(A×C) may be more informative for directionality.")

    # Distance-stratified summary
    if DISTANCE_COL in filtered.columns:
        print("\n=== ΔABC by Distance Bin ===")
        bins = [0, 50_000, 200_000, 500_000, 1_000_000, 5_000_000]
        labels = ["<50kb", "50-200kb", "200-500kb", "500kb-1Mb", "1-5Mb"]
        filtered["dist_bin"] = pd.cut(filtered[DISTANCE_COL], bins=bins, labels=labels, right=False)
        dist_summary = filtered.groupby("dist_bin", observed=True).agg(
            n=("delta_ABC", "size"),
            mean_delta=("delta_ABC", "mean"),
            n_gained=("delta_ABC", lambda x: (x > DELTA_SIG_THRESHOLD).sum()),
            n_lost=("delta_ABC", lambda x: (x < -DELTA_SIG_THRESHOLD).sum()),
        )
        for idx, row in dist_summary.iterrows():
            print(f"  {idx:>10s}: n={int(row['n']):>6,}  mean_Δ={row['mean_delta']:>+.5f}  "
                  f"gained={int(row['n_gained']):>5,}  lost={int(row['n_lost']):>5,}")
        filtered.drop(columns=["dist_bin"], inplace=True)

    # === SAVE ===
    out_all = f"{OUTPUT_DIR}/delta_abc_all_pairs.tsv"
    out_sig = f"{OUTPUT_DIR}/delta_abc_significant.tsv"

    filtered.to_csv(out_all, sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved all thresholded pairs: {out_all} ({len(filtered):,} pairs)")

    sig = filtered[d.abs() > DELTA_SIG_THRESHOLD]
    sig.to_csv(out_sig, sep="\t", index=False, float_format="%.6f")
    print(f"Saved significant pairs (|ΔABC| > {DELTA_SIG_THRESHOLD}): {out_sig} ({len(sig):,} pairs)")

    # Column inventory for Steps 7-9
    print("\n=== Output Column Inventory ===")
    for i, col in enumerate(filtered.columns, 1):
        print(f"  {i:>2}. {col}")


if __name__ == "__main__":
    main()
