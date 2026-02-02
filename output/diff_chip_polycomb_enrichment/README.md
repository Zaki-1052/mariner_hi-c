# Task 3f: Differential H3K27me3/H2AK119ub Overlap with Loop Categories

## Overview

Assessment of whether differential H3K27me3 and H2AK119ub peaks (from DiffBind) overlap with specific loop categories.

**Script:** `scripts/diff_chip_polycomb_enrichment.R`
**Analysis Date:** 2026-02-01

---

## Input Data

### Loop Data
| Dataset | Source | Count |
|---------|--------|-------|
| All loops | `25042-late_outputs/merged_loops/merged_all_results.tsv` | 39,344 |
| Polycomb shared | `output/shared_anchor_analysis/late/polycomb_specific/tables/polycomb_shared_loops.tsv` | 312 |

### Differential Peak Files (DiffBind, 400bp summits)
| Mark | Direction | File | Peaks |
|------|-----------|------|-------|
| H3K27me3 | Down in mutant | `peaks/new/adult_K27me3_down.bed` | 1,563 |
| H3K27me3 | Up in mutant | `peaks/new/adult_K27me3_up.bed` | 1,373 |
| H2AK119ub | Down in mutant | `peaks/new/H2AK119ub_down.bed` | 1,250 |
| H2AK119ub | Up in mutant | `peaks/new/H2AK119ub_up.bed` | 6,164 |

---

## Loop Category Definitions

| Category | Definition | All Loops | Polycomb Shared |
|----------|------------|-----------|-----------------|
| Long-range lost | `down_in_mutant` AND distance > 500kb | 1,558 (4.0%) | 121 (38.8%) |
| Short-range gained | `up_in_mutant` AND distance < 500kb | 3,115 (7.9%) | 111 (35.6%) |
| Unchanged | FDR >= 0.1 | 28,507 (72.5%) | 0 (0.0%) |
| Other differential | Remaining FDR < 0.1 | 6,164 (15.7%) | 80 (25.6%) |

**Note:** Polycomb shared loops are all differential by definition (from task 3c), hence 0 unchanged.

---

## Results: All Loops (n=39,344)

### Overlap Percentages by Loop Category

| Mark | Long-range Lost | Short-range Gained | Unchanged | Other Differential |
|------|-----------------|-------------------|-----------|-------------------|
| K27me3_down | 2.2% (34/1558) | 7.5% (233/3115) | 0.9% (260/28507) | 1.4% (89/6164) |
| K27me3_up | 9.4% (147/1558) | 2.0% (63/3115) | 3.2% (905/28507) | 6.3% (387/6164) |
| H2AK119ub_down | 4.9% (77/1558) | 7.4% (231/3115) | 2.4% (696/28507) | 2.5% (153/6164) |
| H2AK119ub_up | 34.9% (543/1558) | 9.9% (308/3115) | 19.7% (5618/28507) | 29.2% (1802/6164) |

### Fisher's Exact Test Results (vs Unchanged as Reference)

All 12 tests significant at FDR < 0.05:

| Mark | Comparison | Odds Ratio | 95% CI | FDR | Direction |
|------|------------|------------|--------|-----|-----------|
| K27me3_down | Long-range lost vs Unchanged | 2.42 | [1.64, 3.49] | 1.6e-05 | **Enriched** |
| K27me3_down | Short-range gained vs Unchanged | 8.78 | [7.30, 10.57] | 6.9e-102 | **Enriched** |
| K27me3_up | Long-range lost vs Unchanged | 3.18 | [2.63, 3.82] | 3.6e-28 | **Enriched** |
| K27me3_up | Short-range gained vs Unchanged | 0.63 | [0.48, 0.82] | 2.6e-04 | **Depleted** |
| H2AK119ub_down | Long-range lost vs Unchanged | 2.08 | [1.61, 2.65] | 5.4e-08 | **Enriched** |
| H2AK119ub_down | Short-range gained vs Unchanged | 3.20 | [2.73, 3.74] | 4.5e-41 | **Enriched** |
| H2AK119ub_up | Long-range lost vs Unchanged | 2.18 | [1.95, 2.43] | 3.4e-41 | **Enriched** |
| H2AK119ub_up | Short-range gained vs Unchanged | 0.45 | [0.39, 0.50] | 1.7e-45 | **Depleted** |

### Direct Comparison: Long-range Lost vs Short-range Gained

| Mark | Odds Ratio | FDR | Direction |
|------|------------|-----|-----------|
| K27me3_down | 0.28 | 4.3e-15 | Lost < Gained |
| K27me3_up | 5.05 | 2.8e-28 | Lost > Gained |
| H2AK119ub_down | 0.65 | 1.1e-03 | Lost < Gained |
| H2AK119ub_up | 4.87 | 1.9e-90 | Lost > Gained |

---

## Results: Polycomb Shared Loops (n=312)

### Overlap Percentages

| Mark | Long-range Lost (n=121) | Short-range Gained (n=111) | Other Differential (n=80) |
|------|------------------------|---------------------------|--------------------------|
| K27me3_down | 24.0% (29) | 19.8% (22) | 13.8% (11) |
| K27me3_up | 5.0% (6) | 1.8% (2) | 5.0% (4) |
| H2AK119ub_down | 24.0% (29) | 7.2% (8) | 18.8% (15) |
| H2AK119ub_up | 14.0% (17) | 9.0% (10) | 17.5% (14) |

### Fisher's Exact Test (1 significant)

| Mark | Comparison | Odds Ratio | FDR | Direction |
|------|------------|------------|-----|-----------|
| H2AK119ub_down | Long-range lost vs Short-range gained | 3.92 | 2.4e-03 | **Enriched at lost** |

---

## Key Data Patterns

### Pattern 1: K27me3_up and H2AK119ub_up enriched at long-range lost, depleted at short-range gained
- K27me3_up: OR=3.18 at lost, OR=0.63 at gained (5.0x difference)
- H2AK119ub_up: OR=2.18 at lost, OR=0.45 at gained (4.8x difference)

### Pattern 2: K27me3_down and H2AK119ub_down enriched at BOTH categories (but different magnitudes)
- K27me3_down: OR=2.42 at lost, OR=8.78 at gained (strongest enrichment at short-range gained)
- H2AK119ub_down: OR=2.08 at lost, OR=3.20 at gained

### Pattern 3: Polycomb shared loops show H2AK119ub_down preferentially at lost loops
- 24.0% of long-range lost vs 7.2% of short-range gained have H2AK119ub_down overlap
- This is the only significant enrichment in the Polycomb shared subset (OR=3.92, FDR=0.002)

---

## Caveats

1. **400bp summit width** from DiffBind may underestimate broader H3K27me3/H2AK119ub domains
2. **Expected false negatives** due to narrow peak calls
3. **H2AK119ub files** are not timepoint-specific (combined across conditions)
4. **Polycomb shared loops have no unchanged category** - all are differential by definition

---

## Output Files

```
output/diff_chip_polycomb_enrichment/
├── tables/
│   ├── overlap_summary_all_loops.tsv
│   ├── overlap_summary_polycomb_shared.tsv
│   ├── enrichment_tests_all_loops.tsv
│   ├── enrichment_tests_polycomb_shared.tsv
│   ├── all_loops_with_diff_chip_overlap.tsv
│   └── polycomb_loops_with_diff_chip_overlap.tsv
├── plots/
│   ├── 01_overlap_barplot_all_loops/
│   ├── 01_overlap_barplot_polycomb_shared/
│   ├── 02_enrichment_dotplot_all_loops/
│   ├── 02_enrichment_dotplot_polycomb_shared/
│   ├── 03_heatmap_overlap_all_loops/
│   └── 03_heatmap_overlap_polycomb_shared/
├── enrichment_analysis_summary.txt
└── README.md
```

---

## Follow-up (Task 3g)

Aggregate ChIP-seq signal heatmaps at loop anchors using bigWigs to visualize broader domains beyond the 400bp summit limitation.
