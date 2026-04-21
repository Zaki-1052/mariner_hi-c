# stripes/stripenn/scripts/generate_bedpe.py
# Generate JuiceBox BEDPE files from cross_res_merged.tsv for each timepoint.
# Produces 5 BEDPE files per timepoint: 3 significance tiers + diagonal + rectangle formats.

import argparse
import math
import os
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# Color scheme: direction x confidence → RGB string
# ---------------------------------------------------------------------------
COLOR_MAP = {
    ("lost",    "high"):   "0,0,139",       # dark blue
    ("lost",    "medium"): "65,105,225",     # royal blue
    ("lost",    "low"):    "135,206,235",    # sky blue
    ("gained",  "high"):   "139,0,0",        # dark red
    ("gained",  "medium"): "220,20,60",      # crimson
    ("gained",  "low"):    "255,160,122",    # light salmon
}
COLOR_DEFAULT = "128,128,128"  # gray for unchanged / other

BEDPE_HEADER = [
    "chr1", "x1", "x2", "chr2", "y1", "y2",
    "name", "score", "strand1", "strand2", "color",
    "direction", "logFC", "FDR", "source",
]


def compute_score(fdr_series: pd.Series) -> pd.Series:
    """Return -log10(FDR + 1e-300), capped at 300."""
    return (-fdr_series.clip(lower=0).add(1e-300).apply(math.log10)).clip(upper=300)


def assign_color(direction: str, confidence: str) -> str:
    return COLOR_MAP.get((direction, confidence), COLOR_DEFAULT)


def build_bedpe(df: pd.DataFrame, x1_col="pos1", x2_col="pos2",
                y1_col="pos3", y2_col="pos4") -> pd.DataFrame:
    """Construct the 15-column BEDPE dataframe from a filtered stripe dataframe."""
    out = pd.DataFrame()
    out["chr1"]      = df["chr"]
    out["x1"]        = df[x1_col]
    out["x2"]        = df[x2_col]
    out["chr2"]      = df["chr"]       # intrachromosomal: chr == chr2
    out["y1"]        = df[y1_col]
    out["y2"]        = df[y2_col]
    out["name"]      = df["stripe_id"]
    out["score"]     = compute_score(df["FDR"]).round(4)
    out["strand1"]   = "."
    out["strand2"]   = "."
    out["color"]     = df.apply(
        lambda r: assign_color(r["direction"], r["direction_confidence"]), axis=1
    )
    out["direction"] = df["direction"]
    out["logFC"]     = df["logFC"].round(4)
    out["FDR"]       = df["FDR"].apply(lambda v: f"{v:.4e}")
    out["source"]    = df["source"]
    return out.reset_index(drop=True)


def write_bedpe(df_bedpe: pd.DataFrame, path: str, label: str) -> int:
    """Write BEDPE with header; return row count."""
    df_bedpe.to_csv(path, sep="\t", index=False, header=True)
    n = len(df_bedpe)
    print(f"  {label:40s}: {n:6d} stripes  -> {os.path.basename(path)}")
    return n


def process_timepoint(tsv_path: str, tp: str, out_dir: str) -> None:
    """Read cross_res_merged.tsv and write all 5 BEDPE files for one timepoint."""
    print(f"\n[{tp}] Loading {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"  Total rows: {len(df)}")

    os.makedirs(out_dir, exist_ok=True)

    # --- Boolean masks ---------------------------------------------------
    is_differential = df["direction"].isin(["lost", "gained"])
    is_sig05        = df["significant_FDR05"] == True  # noqa: E712  (numpy bool)
    is_highconf     = df["direction_confidence"] == "high"
    is_concordant   = df["resolution_support"] == "both_concordant"

    # --- Tier 1: High-confidence -----------------------------------------
    mask_t1 = is_differential & is_sig05 & is_highconf
    bedpe_t1 = build_bedpe(df[mask_t1])
    write_bedpe(bedpe_t1, os.path.join(out_dir, f"{tp}_stripes_highconf.bedpe"),
                "Tier1 high-confidence (sig+high)")

    # --- Tier 2: All significant -----------------------------------------
    mask_t2 = is_differential & is_sig05
    bedpe_t2 = build_bedpe(df[mask_t2])
    write_bedpe(bedpe_t2, os.path.join(out_dir, f"{tp}_stripes_allsig.bedpe"),
                "Tier2 all-significant (sig)")

    # --- Tier 3: Cross-resolution concordant -----------------------------
    mask_t3 = is_differential & is_concordant
    bedpe_t3 = build_bedpe(df[mask_t3])
    write_bedpe(bedpe_t3, os.path.join(out_dir, f"{tp}_stripes_concordant.bedpe"),
                "Tier3 cross-res concordant")

    # --- Diagonal (from Tier 1): full stripe extent on both axes ---------
    df_t1 = df[mask_t1].copy()
    df_t1["diag_min"] = df_t1[["pos1", "pos3"]].min(axis=1)
    df_t1["diag_max"] = df_t1[["pos2", "pos4"]].max(axis=1)
    bedpe_diag = build_bedpe(df_t1,
                             x1_col="diag_min", x2_col="diag_max",
                             y1_col="diag_min", y2_col="diag_max")
    write_bedpe(bedpe_diag, os.path.join(out_dir, f"{tp}_stripes_diagonal.bedpe"),
                "Diagonal (Tier1, full extent both axes)")

    # --- Rectangle (from Tier 1): narrow anchor x, full extent y --------
    df_t1["extent_min"] = df_t1[["pos1", "pos3"]].min(axis=1)
    df_t1["extent_max"] = df_t1[["pos2", "pos4"]].max(axis=1)
    bedpe_rect = build_bedpe(df_t1,
                             x1_col="pos1",       x2_col="pos2",
                             y1_col="extent_min", y2_col="extent_max")
    write_bedpe(bedpe_rect, os.path.join(out_dir, f"{tp}_stripes_rectangle.bedpe"),
                "Rectangle (Tier1, anchor x full-extent y)")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate JuiceBox BEDPE files from Stripenn cross_res_merged.tsv"
    )
    p.add_argument(
        "--base-dir",
        default="/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/stripes/stripenn/outputs",
        help="Base directory containing {timepoint}/cross_res_merged.tsv sub-directories",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    timepoints = ["250402", "250831"]

    for tp in timepoints:
        tsv_path = os.path.join(args.base_dir, tp, "cross_res_merged.tsv")
        if not os.path.isfile(tsv_path):
            print(f"ERROR: {tsv_path} not found — skipping {tp}", file=sys.stderr)
            continue
        out_dir = os.path.join(args.base_dir, tp)
        process_timepoint(tsv_path, tp, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
