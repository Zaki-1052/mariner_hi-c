#!/usr/bin/env python3
# abc/scripts/step9_cross_reference.py
"""
Cross-reference ABC predictions with:
  A) Differential Hi-C loops (characterized_loops.tsv)
     - Assigns nearest gene to each loop anchor via bedtools closest
       against mm10_tss.bed (replaces unusable Entrez ID columns)
  B) H3K27ac peaks (H3K27acCerebellumLate2.bed)
     - Annotates enhancers with H3K27ac overlap
     - Directional gain/loss analysis stratified by H3K27ac status

Usage: python step9_cross_reference.py
"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
from scipy import stats

# === CONFIGURATION ===
OUTPUT_DIR = "/expanse/lustre/projects/csd940/zalibhai/abc/results"
LOOPS = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/characterized_loops.tsv"
H3K27AC = "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/beds/H3K27acCerebellumLate2.bed"
TSS_BED = "/expanse/lustre/projects/csd940/zalibhai/abc/reference/mm10_tss.bed"
DELTA_ABC = f"{OUTPUT_DIR}/delta_abc_all_pairs.tsv"

ENHANCER_CHR = "chr"
ENHANCER_START = "start"
ENHANCER_END = "end"
GENE_COL = "TargetGene"
ENH_COORDS = [ENHANCER_CHR, ENHANCER_START, ENHANCER_END]

ABC_DELTA_THRESH = 0.01


def write_sorted_bed(df, cols, path):
    """Write a sorted BED file for bedtools."""
    bed = df[cols].drop_duplicates().sort_values(cols[:2])
    bed.to_csv(path, sep="\t", header=False, index=False)
    return len(bed)


def bedtools_intersect(a_bed, b_bed, flags="-u"):
    """Run bedtools intersect, return stdout lines."""
    result = subprocess.run(
        ["bedtools", "intersect", "-a", a_bed, "-b", b_bed] + flags.split(),
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  bedtools error: {result.stderr.strip()}")
        return []
    return result.stdout.strip().split("\n") if result.stdout.strip() else []


def assign_nearest_genes(anchor_df, chr_col, start_col, end_col, tss_bed_sorted):
    """Use bedtools closest to assign nearest gene symbol to loop anchors."""
    anchor_tmp = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    anchors = anchor_df[[chr_col, start_col, end_col]].drop_duplicates()
    anchors_sorted = anchors.sort_values([chr_col, start_col])
    anchors_sorted.to_csv(anchor_tmp, sep="\t", header=False, index=False)
    anchor_tmp.close()

    result = subprocess.run(
        ["bedtools", "closest",
         "-a", anchor_tmp.name,
         "-b", tss_bed_sorted,
         "-d", "-t", "first"],
        capture_output=True, text=True,
    )
    os.unlink(anchor_tmp.name)

    if result.returncode != 0:
        print(f"  bedtools closest error: {result.stderr.strip()}")
        return pd.DataFrame()

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        rows.append({
            "anchor_chr": parts[0],
            "anchor_start": int(parts[1]),
            "anchor_end": int(parts[2]),
            "nearest_gene": parts[6],
            "tss_distance": int(parts[11]),
        })

    return pd.DataFrame(rows)


def cross_reference_loops(delta, enhancers_bed):
    """9A: Overlap ABC enhancers with differential loop anchors."""
    print("\n=== 9A: Differential Loop Overlap ===")

    loops = pd.read_csv(LOOPS, sep="\t")
    print(f"Differential loops: {len(loops):,}")

    # Prepare sorted TSS BED (skip header, sort)
    tss_sorted = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
    subprocess.run(
        f"tail -n +2 {TSS_BED} | sort -k1,1 -k2,2n",
        shell=True, stdout=tss_sorted,
    )
    tss_sorted.close()

    # Assign nearest gene to each anchor
    print("  Assigning nearest genes to loop anchors...")
    a1_genes = assign_nearest_genes(loops, "chr1", "start1", "end1", tss_sorted.name)
    a2_genes = assign_nearest_genes(loops, "chr2", "start2", "end2", tss_sorted.name)
    os.unlink(tss_sorted.name)

    print(f"  Anchor1 genes assigned: {len(a1_genes):,}")
    print(f"  Anchor2 genes assigned: {len(a2_genes):,}")

    # Merge gene assignments back into loops
    loops = loops.merge(
        a1_genes.rename(columns={"nearest_gene": "anchor1_gene", "tss_distance": "anchor1_tss_dist"}),
        left_on=["chr1", "start1", "end1"],
        right_on=["anchor_chr", "anchor_start", "anchor_end"],
        how="left",
    ).drop(columns=["anchor_chr", "anchor_start", "anchor_end"], errors="ignore")

    loops = loops.merge(
        a2_genes.rename(columns={"nearest_gene": "anchor2_gene", "tss_distance": "anchor2_tss_dist"}),
        left_on=["chr2", "start2", "end2"],
        right_on=["anchor_chr", "anchor_start", "anchor_end"],
        how="left",
    ).drop(columns=["anchor_chr", "anchor_start", "anchor_end"], errors="ignore")

    print("  Sample anchor gene assignments:")
    for _, row in loops.head(5).iterrows():
        print(f"    {row['loop_id']}: {row.get('anchor1_gene','?')} ({row.get('anchor1_tss_dist',0):.0f}bp) "
              f"<-> {row.get('anchor2_gene','?')} ({row.get('anchor2_tss_dist',0):.0f}bp)")

    # --- Enhancer-anchor overlap ---
    n_enh = sum(1 for _ in open(enhancers_bed))

    a1_bed = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    a2_bed = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    write_sorted_bed(loops, ["chr1", "start1", "end1"], a1_bed)
    write_sorted_bed(loops, ["chr2", "start2", "end2"], a2_bed)

    for label, anchor_bed in [("anchor1", a1_bed), ("anchor2", a2_bed)]:
        hits = bedtools_intersect(enhancers_bed, anchor_bed)
        print(f"  ABC enhancers overlapping {label}: {len(hits):,} / {n_enh:,}")

    os.unlink(a1_bed)
    os.unlink(a2_bed)

    # --- Gene-level overlap ---
    abc_genes = set(delta[GENE_COL].unique())
    loop_genes = set(loops["anchor1_gene"].dropna()) | set(loops["anchor2_gene"].dropna())
    overlap_genes = abc_genes & loop_genes
    print(f"  Unique genes in loop annotations: {len(loop_genes):,}")
    print(f"  ABC target genes in loop annotations: {len(overlap_genes):,} / {len(abc_genes):,}")

    # --- Directional concordance ---
    gene_directions = {}
    for _, row in loops.iterrows():
        d = row.get("direction")
        if pd.isna(d):
            continue
        for g in [row.get("anchor1_gene"), row.get("anchor2_gene")]:
            if pd.notna(g) and g in abc_genes:
                gene_directions.setdefault(g, []).append(d)

    gene_dabc = delta.groupby(GENE_COL)["delta_ABC"].apply(
        lambda x: x.loc[x.abs().idxmax()]
    ).to_dict()

    conc, disc, ambig = 0, 0, 0
    for g, dirs in gene_directions.items():
        if g not in gene_dabc:
            continue
        dabc = gene_dabc[g]
        n_up = sum(1 for d in dirs if d == "up_in_mutant")
        n_down = sum(1 for d in dirs if d == "down_in_mutant")
        if n_up > n_down:
            majority = "up_in_mutant"
        elif n_down > n_up:
            majority = "down_in_mutant"
        else:
            ambig += 1
            continue

        if (majority == "up_in_mutant" and dabc > 0) or \
           (majority == "down_in_mutant" and dabc < 0):
            conc += 1
        else:
            disc += 1

    total_dir = conc + disc
    if total_dir > 0:
        print(f"\n  Directional concordance (loop vs ΔABC): {conc}/{total_dir} "
              f"({100*conc/total_dir:.1f}%)")
    if ambig > 0:
        print(f"  Ambiguous (tied loop directions): {ambig}")

    # Save annotated loops
    loops_out = f"{OUTPUT_DIR}/loops_with_gene_assignments.tsv"
    loops.to_csv(loops_out, sep="\t", index=False)
    print(f"  Saved: {loops_out}")


def cross_reference_h3k27ac(delta, enhancers_bed):
    """9B: Annotate ABC enhancers with H3K27ac peak overlap.

    Includes directional gain/loss analysis stratified by H3K27ac status
    to test whether BAP1-KO preferentially disrupts active enhancers.
    """
    print("\n=== 9B: H3K27ac Peak Overlap ===")

    n_enh = sum(1 for _ in open(enhancers_bed))

    hit_lines = bedtools_intersect(enhancers_bed, H3K27AC, "-u")
    n_with = len(hit_lines)
    n_without = n_enh - n_with
    print(f"  ABC enhancers overlapping H3K27ac: {n_with:,} / {n_enh:,}")
    print(f"  ABC enhancers without H3K27ac:     {n_without:,} / {n_enh:,}")

    overlap_set = set()
    for line in hit_lines:
        parts = line.split("\t")
        overlap_set.add((parts[0], int(parts[1]), int(parts[2])))

    enh_key = delta[ENH_COORDS].copy()
    enh_key["_key"] = list(zip(enh_key[ENHANCER_CHR],
                                enh_key[ENHANCER_START].astype(int),
                                enh_key[ENHANCER_END].astype(int)))
    delta["has_H3K27ac"] = enh_key["_key"].isin(overlap_set)

    with_k27 = delta[delta["has_H3K27ac"]]
    without_k27 = delta[~delta["has_H3K27ac"]]
    print(f"\n  E-G pairs at H3K27ac+ enhancers: {len(with_k27):,}")
    print(f"  E-G pairs at H3K27ac- enhancers: {len(without_k27):,}")
    print(f"  Mean |ΔABC| H3K27ac+: {with_k27['delta_ABC'].abs().mean():.5f}")
    print(f"  Mean |ΔABC| H3K27ac-: {without_k27['delta_ABC'].abs().mean():.5f}")

    # --- Directional gain/loss by H3K27ac status ---
    # Tests whether BAP1-KO loss of enhancer connections concentrates
    # at H3K27ac-marked (active) enhancers, as predicted by:
    #   BAP1 loss → increased H2AK119ub → Polycomb repression
    #   → preferential disruption of active enhancers
    print("\n=== 9C: H3K27ac-Stratified Gain/Loss Asymmetry ===")

    for label, sub in [("H3K27ac+", with_k27), ("H3K27ac-", without_k27)]:
        n_gained = (sub["delta_ABC"] > ABC_DELTA_THRESH).sum()
        n_lost = (sub["delta_ABC"] < -ABC_DELTA_THRESH).sum()
        n_unchanged = len(sub) - n_gained - n_lost
        ratio = n_lost / n_gained if n_gained > 0 else float("inf")
        mean_d = sub["delta_ABC"].mean()

        print(f"\n  {label} (n={len(sub):,} E-G pairs):")
        print(f"    Gained (ΔABC > {ABC_DELTA_THRESH}):  {n_gained:,}")
        print(f"    Lost (ΔABC < -{ABC_DELTA_THRESH}):   {n_lost:,}")
        print(f"    Unchanged:                 {n_unchanged:,}")
        print(f"    Lost/Gained ratio:         {ratio:.2f}")
        print(f"    Mean ΔABC:                 {mean_d:+.5f}")

    # Statistical test: is the loss bias significantly stronger for H3K27ac+ vs H3K27ac-?
    # Chi-square test on the 2x2 table of (gained/lost) × (H3K27ac+/H3K27ac-)
    k27_gained = (with_k27["delta_ABC"] > ABC_DELTA_THRESH).sum()
    k27_lost = (with_k27["delta_ABC"] < -ABC_DELTA_THRESH).sum()
    nok27_gained = (without_k27["delta_ABC"] > ABC_DELTA_THRESH).sum()
    nok27_lost = (without_k27["delta_ABC"] < -ABC_DELTA_THRESH).sum()

    contingency = np.array([[k27_gained, k27_lost],
                            [nok27_gained, nok27_lost]])
    chi2, chi2_p, dof, expected = stats.chi2_contingency(contingency)

    print(f"\n  Chi-square test (gain/loss asymmetry differs by H3K27ac status):")
    print(f"    Contingency table:")
    print(f"                 Gained    Lost")
    print(f"      H3K27ac+   {k27_gained:>6,}  {k27_lost:>6,}  (ratio {k27_lost/k27_gained:.2f})")
    print(f"      H3K27ac-   {nok27_gained:>6,}  {nok27_lost:>6,}  (ratio {nok27_lost/nok27_gained:.2f})")
    print(f"    Chi2 = {chi2:.2f}, p = {chi2_p:.2e}")

    # Mann-Whitney test: is the ΔABC distribution shifted more negative for H3K27ac+?
    u_stat, mw_p = stats.mannwhitneyu(
        with_k27["delta_ABC"], without_k27["delta_ABC"],
        alternative="less",
    )
    print(f"\n  Mann-Whitney U test (H3K27ac+ ΔABC < H3K27ac- ΔABC):")
    print(f"    U = {u_stat:,.0f}, p = {mw_p:.2e}")

    return delta


def main():
    delta = pd.read_csv(DELTA_ABC, sep="\t")
    print(f"ΔABC pairs: {len(delta):,}")

    enhancers_bed = tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name
    n_enh = write_sorted_bed(delta, ENH_COORDS, enhancers_bed)
    print(f"Unique ABC enhancers: {n_enh:,}")

    cross_reference_loops(delta, enhancers_bed)
    delta = cross_reference_h3k27ac(delta, enhancers_bed)

    os.unlink(enhancers_bed)

    out = f"{OUTPUT_DIR}/delta_abc_annotated.tsv"
    delta.to_csv(out, sep="\t", index=False, float_format="%.6f")
    print(f"\nSaved: {out} ({len(delta):,} pairs)")


if __name__ == "__main__":
    main()
