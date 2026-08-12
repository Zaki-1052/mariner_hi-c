# Plan: DEG-Centric Enhancer Contact Analysis (TODO Section 3)

## Context

TODO Section 3 asks for a DEG-centric analysis that quantifies **Hi-C contact strength** between DEGs and enhancers, comparing against matched invariant genes. The existing `se_deg_proximity.R` (Section 1a) measures genomic **proximity** (how many SEs/enhancers fall within 2 Mb of a DEG) but not actual contact frequency. The ABC pipeline has already computed per-enhancer-gene delta-ABC scores for 179K pairs (`abc/results/delta_abc_all_pairs.tsv`) and per-gene summaries (`abc/results/gene_level_summary.tsv`). The missing piece is joining these ABC contact strength metrics to the DEG vs. invariant framework from `se_utils.R`.

Meeting note: "maybe will clear up ABC plot" — the genome-wide ABC scatter is noisy; filtering to pre-matched DEG/invariant sets should produce a cleaner signal.

---

## Approach: New Script `loops/scripts/se_abc_integration.R`

A new standalone script (not an extension of `se_deg_proximity.R`) because it has different inputs (ABC TSVs, not BED proximity windows) and a different analytical framing (contact strength, not genomic distance).

---

## Inputs

| File | Purpose |
|---|---|
| `abc/results/gene_level_summary.tsv` (13,749 genes) | Per-gene ABC summary: `max_delta_abc`, `sum_delta_abc`, `n_gained`, `n_lost` |
| `abc/results/delta_abc_all_pairs.tsv` (179K pairs) | Per-enhancer-gene: `delta_ABC`, `hic_contact_pl_scaled_adj_{WT,KO}`, `class` |
| `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` | Late DEGs |
| `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx` | Early DEGs |
| `peaks/Superenhancers_P60.bed`, `peaks/Superenhancers_encode.bed` | SE coordinates |
| `peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt` | K27ac DiffBind (Fold, FDR) |

Reuses: `loops/scripts/utils/se_utils.R` (DEG loading, invariant matching, K27ac classification, color schemes), `loops/scripts/utils/multi_format_output.R` (plot saving).

## Output Directory

`loops/output/superenhancer_analysis/3a_deg_abc_contact/{late,early}/`

---

## Analysis Steps

### Step 1: Load DEGs and Invariant Genes
- Call `load_rna_for_se()` + `select_invariant_genes()` from `se_utils.R` (same logic as `se_deg_proximity.R`)
- Build `gene_tbl` with `gene_class` = DEG_up / DEG_down / invariant
- Expected: ~553 DEG_up, ~1,055 DEG_down, ~4,698 invariant (late)

### Step 2: Gene-Level ABC Join (TODO 3a)
- Read `gene_level_summary.tsv`, rename `TargetGene` → `gene_symbol`
- Drop its `log2FC`/`padj` columns (use RNA-seq values from gene_tbl)
- Inner join with `gene_tbl` on `gene_symbol` (~5,600-5,800 matches)
- This gives per-gene: `sum_delta_abc`, `max_delta_abc`, `n_gained`, `n_lost`, `gene_class`

### Step 3: Gene-Level Plots + Statistics (TODO 3a)
1. **`delta_abc_violin`** — Violin/box of `sum_delta_abc` by gene_class (DEG_up / DEG_down / invariant), Wilcoxon comparisons. The primary "cleaner signal" plot.
2. **`max_delta_abc_violin`** — Same for `max_delta_abc` (single strongest enhancer link)
3. **`n_gained_lost_bar`** — Mean n_gained vs n_lost per gene_class with error bars
4. **`delta_abc_vs_lfc_scatter`** — sum_delta_abc vs log2FoldChange for DEGs only, Spearman rho. The focused version of the noisy genome-wide ABC scatter.

Stats: 3 pairwise Wilcoxon tests (sum_delta_abc, max_delta_abc), Kruskal-Wallis omnibus, Spearman rho

### Step 4: Enhancer-Level ABC Pairs Join
- Read `delta_abc_all_pairs.tsv`
- **Filter out promoter class** (`class != "promoter"`, removes ~88K → ~92K distal pairs remain)
- Filter to genes in `gene_tbl`, attach `gene_class` column
- Compute `log2fc_contact` from `hic_contact_pl_scaled_adj_{WT,KO}` with pseudocount

### Step 5: SE Classification of Enhancers (TODO 3b)
- Build GRanges from ABC enhancer coordinates (chr/start/end)
- `findOverlaps()` against P60 + ENCODE SE GRanges
- Tag each pair: `is_se` (TRUE/FALSE), `se_source` (P60/ENCODE/NA)

### Step 6: K27ac DiffBind Join (TODO 3b)
- Load DiffBind GRanges via `load_k27ac_diffbind()`
- For k27ac class: use `classify_k27ac_change()` from `se_utils.R`
- For k27ac **magnitude** (needed for correlation plot): do `findOverlaps()` inline and retain the `fold` value from the best-hit DiffBind peak (the utility function only returns class, not fold)

### Step 7: Enhancer-Level Plots + Statistics (TODO 3b)
5. **`se_vs_regular_contact_violin`** — log2fc_contact for SE vs regular enhancers, faceted by gene_class
6. **`se_contact_proportion_bar`** — Proportion SE vs regular by gene_class, Fisher test
7. **`k27ac_contact_scatter`** — K27ac Fold (x) vs log2fc_contact (y), colored by gene_class, Spearman rho. Directly answers: "does K27ac change magnitude correlate with contact frequency change?"
8. **`k27ac_class_contact_violin`** — log2fc_contact by k27ac_class (gained/stable/lost)
9. **`se_k27ac_subclass_bar`** — k27ac class distribution among SE contacts by gene_class, Fisher concordance test

Stats: Fisher tests for SE enrichment and K27ac concordance, Wilcoxon for SE vs regular contact, Spearman K27ac fold vs contact

### Step 8: Write Outputs
- `tables/deg_abc_gene_level.tsv` — gene-level joined table
- `tables/deg_abc_pairs_classified.tsv` — enhancer-level with SE + K27ac classification
- `tables/se_contact_summary.tsv` — aggregated counts by gene_class x is_se x k27ac_class
- `statistics/abc_contact_statistics.txt` — all test results

---

## Key Implementation Notes

- **Promoter filter is critical**: `delta_abc_all_pairs.tsv` includes promoter self-contacts (~88K rows). These must be excluded — only `genic` and `intergenic` pairs are actual enhancer-gene contacts.
- **Use `hic_contact_pl_scaled_adj_*`**, not `hic_contact_*`, for log2FC. The `_pl_scaled_adj` column is the KR-normalized, pseudocounted, power-law-scaled value used in ABC scoring.
- **Gene symbols align safely**: both RNA-seq Excel and ABC gene_level_summary use UCSC-style symbols (confirmed by existing pipeline).
- **SE BEDs are broad regions (10s of kb)**, ABC enhancers are 500bp summits — `findOverlaps` correctly identifies which ABC enhancers fall within SE regions.
- **No new packages needed**: all dependencies (tidyverse, GenomicRanges, ggpubr, patchwork) already used in se_deg_proximity.R.

---

## Verification

1. Run: `Rscript loops/scripts/se_abc_integration.R --timepoint late`
2. Check output exists in `loops/output/superenhancer_analysis/3a_deg_abc_contact/late/`
3. Verify gene counts in stats file match expected (~553 DEG_up, ~1,055 DEG_down, ~4,698 invariant)
4. Inspect `delta_abc_violin` — expect DEG_down to show more negative sum_delta_abc than invariant
5. Inspect `k27ac_contact_scatter` — expect positive Spearman rho (gained K27ac ↔ gained contact)
6. Run with `--timepoint early` for early timepoint
7. Spot-check `deg_abc_pairs_classified.tsv` — confirm no `class == "promoter"` rows, confirm `is_se` flags are populated
