# Refactor Prompt: Re-path `data/scripts/` to Source from `data/`

## Context

The scripts in `data/scripts/` are **copies** of the original analysis scripts from across the repository (`scripts/`, `tads/scripts/`, `abc/scripts/`, `peaks/scripts/`). They currently reference input/output paths relative to their original locations (e.g., `output/`, `outputs/`, `peaks/`, `tads/results/`, `abc/results/`).

**Goal:** Update these copied scripts so they:
1. Read input data from `data/tsvs/` and `data/upstream/` instead of original output directories
2. Write output plots to `data/plots/` instead of original plot directories
3. Source `multi_format_output.R` from `data/scripts/_shared/` instead of `scripts/utils/`
4. Become self-contained: runnable from the repo root using only `data/` contents

**This is a mechanical refactor** — no logic changes, no new analyses. Just find-and-replace of file paths.

## Rules

1. **ONLY modify files inside `data/scripts/`**. The originals in `scripts/`, `tads/scripts/`, `abc/scripts/`, `peaks/scripts/` must remain untouched.
2. **Do not change any analysis logic.** Only update `source()`, `read.*()`, `readRDS()`, `write.*()`, `saveRDS()`, `pdf()`, `ggsave()`, `save_multiformat_ggplot()`, and directory-creation paths.
3. **Preserve the original path in a comment** next to each change, e.g.:
   ```r
   # Original: output/h2ak119ub_loop_integration/late/k119ub_enrichment_summary.tsv
   k119ub_file <- "data/tsvs/figure_3_epigenetic_integration/3A_k119ub_enrichment_global.tsv"
   ```
4. **If a script reads an intermediate file that is NOT in `data/`**, leave the path unchanged and add a `# TODO: not in data/` comment. Do not create fallbacks or stubs.

## Step 0: Shared Utility

Nearly every script sources `multi_format_output.R`. In all copied scripts, update:

```r
# FIND (variations):
source("scripts/utils/multi_format_output.R")
source("../scripts/utils/multi_format_output.R")
source(file.path(scripts_dir, "utils/multi_format_output.R"))

# REPLACE WITH:
source("data/scripts/_shared/multi_format_output.R")
```

## Step 1: Per-Script Path Mapping

For each script below, grep for `source()`, `read.table()`, `read.csv()`, `read_tsv()`, `readRDS()`, `fread()`, `write.table()`, `write.csv()`, `write_tsv()`, `saveRDS()`, `pdf()`, `svg()`, `jpeg()`, `ggsave()`, `save_multiformat_ggplot()`, `dir.create()`, and `file.path()` calls. Update paths according to the mappings below.

### Approach

For each script:
1. Read the script
2. Identify all I/O paths (source, read, write)
3. Map each path to its `data/` equivalent using the tables below
4. Make the substitution, preserving original as comment
5. Verify no paths are missed by grepping for common path fragments like `output/`, `outputs/`, `results/`, `peaks/`

---

### `_shared/multi_format_output.R`
No changes needed — this is a utility that other scripts source.

---

### `figure_1_tads_boundaries_compartments/`

#### `tad_volcano_plot.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/tad-pc-analysis/output/tad_analysis/{early,late}/` (reads) | `data/tsvs/figure_1_tads_boundaries_compartments/` |
| `tads/tad-pc-analysis/output/tad_analysis/` (writes: TSV, PDF) | `data/tsvs/figure_1_tads_boundaries_compartments/` (TSV), `data/plots/figure_1_tads_boundaries_compartments/` (plots) |
| Output filenames: `tad_volcano_{relaxed,standard}.{pdf,svg}`, `tad_significant_*.tsv`, `tad_all_annotated.tsv` | Prefix with `1B_{early,late}_` per the naming convention |

#### `tad_visualizations.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/results/{early,late}/final/` (reads boundary BEDs) | `data/tsvs/figure_1_tads_boundaries_compartments/` |
| `tads/results/visualizations/late/syt1_nav3_focus/` (writes) | `data/plots/figure_1_tads_boundaries_compartments/` (1E plots), `data/tsvs/figure_1_tads_boundaries_compartments/` (1E TSVs) |
| `tads/results/visualizations/late/` (GO dotplots, boundary genes) | `data/plots/figure_5_model_functional/` (5A plots), `data/tsvs/figure_5_model_functional/` (5A TSVs) |

#### `compartment_volcano_plot.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/tad-pc-analysis/output/compartment_analysis/` (reads) | `data/tsvs/figure_1_tads_boundaries_compartments/` |
| Output volcanos and TSVs | `data/plots/figure_1_tads_boundaries_compartments/` (1D plots), `data/tsvs/figure_1_tads_boundaries_compartments/` (1D TSVs) |

#### `compartment_genome_percentage.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/tad-pc-analysis/output/compartment_analysis/` (reads) | `data/tsvs/figure_1_tads_boundaries_compartments/` |
| Output pie/bar charts and summary | `data/plots/figure_1_tads_boundaries_compartments/` (1D plots), `data/tsvs/figure_1_tads_boundaries_compartments/` (1D summary) |

#### `tad_chip_classification.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/results/{early,late}/final/` (reads boundary data) | `data/tsvs/figure_1_tads_boundaries_compartments/` |
| `peaks/beds/` (reads ChIP peak files) | Keep as-is or `# TODO: not in data/` |
| `tads/results/visualizations/chip/{early,late}/` (writes) | `data/plots/figure_3_epigenetic_integration/` (3B plots), `data/tsvs/figure_3_epigenetic_integration/` (3B TSVs), `data/tsvs/figure_1_tads_boundaries_compartments/` (1F TSVs) |

#### `boundary_loop_crossref.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/results/late/boundary_loop_analysis/` (reads/writes) | `data/tsvs/figure_1_tads_boundaries_compartments/` (1F TSVs), `data/plots/figure_1_tads_boundaries_compartments/` (1F permutation plot) |
| Loop input files from `outputs/` or `peaks/` | `data/upstream/loop_calls/` or `# TODO: not in data/` |

---

### `figure_2_loop_rewiring/`

#### `edgeR.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `outputs/*/edgeR_results_res_*kb/primary_analysis/` (writes) | `data/tsvs/figure_2_loop_rewiring/` (2A TSVs) |
| `outputs/*/edgeR_results_res_*kb/plots/` (writes) | `data/plots/figure_2_loop_rewiring/` |
| Intermediate RDS files from `outputs/res_*kb/` (reads) | `# TODO: not in data/` (pipeline intermediates) |

#### `visualizations.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `outputs/*/edgeR_results_res_*kb/` (reads edgeR results) | `data/tsvs/figure_2_loop_rewiring/` |
| `outputs/*/visualizations/volcano/` (writes) | `data/plots/figure_2_loop_rewiring/` (2A merged volcanos) |

#### `loop_distance_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/loops_visualization_extended/{early,late}/` (writes) | `data/plots/figure_2_loop_rewiring/` (2B CDF), `data/tsvs/figure_2_loop_rewiring/` (2B distance summary) |
| `output/loops_visualization_extended/late/` (GO outputs) | `data/plots/figure_5_model_functional/` (5A GO comparison), `data/tsvs/figure_5_model_functional/` (5A GO TSVs) |
| Loop rewriting summary | `data/plots/supplemental/` |
| Characterized loops input | `data/upstream/loop_calls/` or `data/tsvs/figure_2_loop_rewiring/2H_late_characterized_loops.tsv` |

#### `apa_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `outputs/*/apa_results/res_*kb/merged_loops/` (writes) | `data/plots/figure_2_loop_rewiring/apa/` (2C heatmaps), `data/tsvs/figure_2_loop_rewiring/` (2C enrichment/stats TSVs) |
| `.hic` files on HPC (reads) | `# TODO: not in data/` (HPC-only) |

#### `loop_distance_k27me3_filtered.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/loops_k27me3_filtered/{early,late}/` (writes) | `data/plots/figure_2_loop_rewiring/` (2E K27me3 CDFs) |

#### `loop_distance_mark_filtered.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/loops_mark_filtered/late/{H3K27ac,H3K27me3,CTCF}/` (writes) | `data/plots/figure_2_loop_rewiring/` (2E, 2I mark-specific CDFs) |

#### `annotate_loops_extended_peaks.R` (from `peaks/scripts/`)
| Original Path Pattern | New Path |
|----------------------|----------|
| `peaks/loop_annotation_extended/{early,late}/` (writes) | `data/tsvs/figure_2_loop_rewiring/` (2F, 2G TSVs), `data/plots/figure_2_loop_rewiring/` (2F, 2G plots) |
| `peaks/beds/` (reads ChIP peaks) | `# TODO: not in data/` |

#### `annotate_loops_extended.R` (from `scripts/`)
| Original Path Pattern | New Path |
|----------------------|----------|
| Similar to above but may use different paths | Check and map accordingly |

#### `deg_loop_anchor_violin.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/deg_loop_violin/late/` (writes) | `data/plots/figure_2_loop_rewiring/` (2H violin), `data/plots/figure_5_model_functional/` (5B violins), `data/tsvs/figure_5_model_functional/` (5B TSVs) |
| RNA-seq input | `data/upstream/rna_seq/adult_rnaseq_results.xlsx` |

---

### `figure_3_epigenetic_integration/`

#### `h2ak119ub_loop_integration.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `peaks/k119ub_anchor_signal.tsv` (reads) | `data/upstream/chip_peaks/k119ub_anchor_signal.tsv` |
| `output/h2ak119ub_loop_integration/late/` (writes all 3A outputs) | `data/plots/figure_3_epigenetic_integration/` (3A plots), `data/tsvs/figure_3_epigenetic_integration/` (3A TSVs) |
| Characterized loops input | `data/upstream/loop_calls/late_characterized_loops.tsv` |

#### `preprocess_k119ub_anchor_signal.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| DiffBind output paths (reads) | `# TODO: not in data/` (HPC intermediate) |
| `peaks/k119ub_anchor_signal.tsv` (writes) | `data/upstream/chip_peaks/k119ub_anchor_signal.tsv` |

#### `diff_chip_polycomb_enrichment.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/diff_chip_polycomb_enrichment/` (writes) | `data/plots/figure_3_epigenetic_integration/` (3C enrichment dotplots), `data/tsvs/figure_3_epigenetic_integration/` (3C enrichment TSVs) |

#### `loop_compartment_crossref.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/loop_compartment_crossref/` (writes) | `data/plots/figure_3_epigenetic_integration/` (3C loop_cmpt plots), `data/tsvs/figure_3_epigenetic_integration/` (3C loop_compartment TSV) |
| Compartment data input | `data/tsvs/figure_1_tads_boundaries_compartments/1D_compartment_all_annotated.tsv` |

#### `timepoint_comparison.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/results/visualizations/comparison/` (writes) | `data/plots/figure_3_epigenetic_integration/` (3D plots), `data/tsvs/figure_3_epigenetic_integration/` (3D stats TSV) |

---

### `figure_4_abc_analysis/`

#### `step12_activity_contact_scatter.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/activity_contact_scatter/` (writes) | `data/plots/figure_4_abc_analysis/` (4A plots) |
| `abc/results/delta_abc_*.tsv` (reads) | `data/tsvs/figure_4_abc_analysis/4A_delta_abc_*.tsv` |

#### `step12b_promoter_distal_scatter.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/activity_contact_scatter/` (writes) | `data/plots/figure_4_abc_analysis/4A_raw_delta_promoter_distal.svg` |

#### `concordance_pie_chart.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/concordance_pie/` (writes) | `data/plots/figure_4_abc_analysis/` (4B pie charts) |

#### `step13_discordant_gene_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/discordant_analysis/` (writes) | `data/plots/figure_4_abc_analysis/` (4B discordant plots) |
| `abc/results/discordant_gene_characteristics.tsv` (writes) | `data/tsvs/figure_4_abc_analysis/4B_discordant_gene_characteristics.tsv` |

#### `step13b_go_enrichment.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/discordant_analysis/` (writes GO/KEGG plots) | `data/plots/figure_4_abc_analysis/` (4B GO plot) |
| `abc/results/figures/discordant_analysis/go_bp_enrichment_results.tsv` | `data/tsvs/figure_4_abc_analysis/4B_discordant_go_bp.tsv` |
| `abc/results/figures/discordant_analysis/kegg_enrichment_results.tsv` | `data/tsvs/figure_4_abc_analysis/4B_discordant_kegg.tsv` |

#### `step13c_k119ub_concordance.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/discordant_analysis/` (writes K119ub plots) | `data/plots/figure_4_abc_analysis/` (4B K119ub plots) |

#### `step10_k119ub_abc_correlation.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/figures/k119ub_correlation/` (writes) | `data/plots/figure_4_abc_analysis/` (4C volcano, 4F contingency/boxplot) |
| `abc/results/k119ub_abc_correlation_summary.tsv` (writes) | `data/tsvs/figure_4_abc_analysis/4F_k119ub_abc_correlation_summary.tsv` |
| `abc/results/k119ub_abc_enhancer_merged.tsv` (writes) | `data/tsvs/figure_4_abc_analysis/4F_k119ub_abc_enhancer_merged.tsv` |
| K119ub enhancer signal (reads) | `data/upstream/chip_peaks/k119ub_enhancer_signal.tsv` |

#### `step11_enhancer_subset_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/enhancer_subset_analysis/` (writes all 4D/4E/4F outputs) | `data/plots/figure_4_abc_analysis/` (4D/4E/4F plots), `data/tsvs/figure_4_abc_analysis/` (4D/4E/4F TSVs) |
| `abc/enhancer_subsets/` (writes class files) | `data/tsvs/figure_4_abc_analysis/4D_enhancer_classes_*.tsv` |

---

### `figure_5_model_functional/`

#### `deg_tad_violin.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `tads/results/visualizations/late/deg_violin/` (writes) | `data/plots/figure_5_model_functional/` (5B TAD violin), `data/tsvs/figure_5_model_functional/` (5B boundary genes) |
| RNA-seq input | `data/upstream/rna_seq/adult_rnaseq_results.xlsx` |

#### `step9b_paired_anchor_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/paired_anchor_plots/` (writes) | `data/plots/figure_4_abc_analysis/` (4E logFC_vs_deltaABC), `data/plots/figure_5_model_functional/` (5A GO/KEGG), `data/plots/supplemental/` (paired_anchor_panel) |
| `abc/results/paired_anchor_*.tsv` (writes) | `data/tsvs/figure_5_model_functional/` (5A/5C TSVs) |
| `abc/results/gene_level_summary.tsv` (writes) | `data/tsvs/figure_5_model_functional/5B_gene_level_summary.tsv` |

#### `network_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| Output network figures and TSVs | `data/plots/figure_5_model_functional/` (5C plots), `data/tsvs/figure_5_model_functional/` (5C TSVs) |
| Reads loop data, boundary genes, ABC summary | Map to `data/tsvs/` or `data/upstream/` equivalents |

#### `structural_heatmap.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| Output heatmaps and top-50 TSVs | `data/plots/figure_5_model_functional/` (5D plots), `data/tsvs/figure_5_model_functional/` (5D TSVs) |

---

### `supplemental/`

#### `shared_anchor_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/shared_anchor_analysis/late/` (writes) | `data/tsvs/supplemental/` (shared_anchor TSVs), `data/plots/supplemental/` (shared_anchor plots) |

#### `polycomb_shared_anchor_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/shared_anchor_analysis/late/polycomb_specific/` (writes) | `data/tsvs/supplemental/`, `data/plots/supplemental/` |

#### `shared_anchor_boundary_analysis.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/shared_anchor_boundary_analysis/late/` (writes) | `data/tsvs/supplemental/shared_boundary_*.tsv`, `data/plots/supplemental/` |

#### `ctcf_stripe_crossref.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/ctcf_stripe_crossref/late/` (writes) | `data/tsvs/supplemental/ctcf_stripe_*.tsv`, `data/plots/supplemental/ctcf_stripe_*.svg` |

#### `apa_shared_anchors.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `output/apa_shared_anchors/late/` (writes) | `data/plots/supplemental/apa_shared_*.svg`, `data/tsvs/supplemental/` |

#### `step11b_homer_motif_visualization.R`
| Original Path Pattern | New Path |
|----------------------|----------|
| `abc/results/enhancer_subset_analysis/homer_motif_viz/` (writes) | `data/plots/supplemental/homer_narrative_composite.svg` |
| `abc/results/enhancer_subset_analysis/homer_motif_summary_stats.tsv` (writes) | `data/tsvs/supplemental/homer_motif_summary_stats.tsv` |

---

## Step 2: Cross-Reference Checklist

For **each** of the 39 scripts, verify:

- [ ] All `source()` calls updated to `data/scripts/_shared/multi_format_output.R`
- [ ] All `read.*()` / `readRDS()` / `fread()` input paths point to `data/tsvs/`, `data/upstream/`, or are marked `# TODO: not in data/`
- [ ] All write paths (TSV, plot) point to `data/tsvs/` or `data/plots/`
- [ ] All `dir.create()` calls updated to `data/` subdirectories
- [ ] No references to `output/`, `outputs/`, `abc/results/`, `tads/results/`, or `peaks/` remain (unless marked TODO)
- [ ] Original path preserved as comment next to each change
- [ ] Script runs without error from repo root (after creating necessary output dirs)

## Step 3: Automated Verification

After all changes, run:

```bash
# Check for any remaining references to original output directories
for script in $(find data/scripts/ -name "*.R"); do
  echo "=== $script ==="
  grep -n 'output/' "$script" | grep -v '# Original:' | grep -v '# TODO:'
  grep -n 'outputs/' "$script" | grep -v '# Original:' | grep -v '# TODO:'
  grep -n 'abc/results/' "$script" | grep -v '# Original:' | grep -v '# TODO:'
  grep -n 'tads/results/' "$script" | grep -v '# Original:' | grep -v '# TODO:'
  grep -n 'peaks/' "$script" | grep -v '# Original:' | grep -v '# TODO:'
done
```

Any matches indicate missed path updates.

## File-Level Cross-Reference

Every file in `data/tsvs/` and `data/plots/` should trace back to exactly one script in `data/scripts/`. Use `data/INDEX.md` as the authoritative mapping — each row has a `Script` column identifying the source.

---

## Notes

- **HPC-only inputs** (`.hic` files, raw bigwigs, DiffBind results) cannot be bundled into `data/`. Scripts that read these should have their HPC paths marked `# TODO: not in data/` and will not be fully runnable locally.
- **Pipeline intermediates** (RDS files from `extract_counts.R`, `aggregate.R`, etc.) are not in `data/`. Scripts that depend on these intermediates will need the full pipeline run first.
- **The `data/upstream/` directory** contains the key shared reference files (RNA-seq results, K119ub signals, loop calls) that many scripts need. These should be the primary read targets for cross-figure inputs.
