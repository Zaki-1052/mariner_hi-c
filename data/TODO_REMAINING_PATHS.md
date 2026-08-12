# TODO: Remaining HPC-Only Paths in `data/scripts/`

After the path refactor, only **12 TODO markers remain across 3 scripts**. All reference HPC pipeline intermediates that require the full analysis pipeline to generate.

---

## `edgeR.R` (8 TODOs)

Pipeline RDS intermediates and output directories that require running the full mariner pipeline (`prep_loops.R` -> `extract_counts.R` -> `aggregate.R`):

| Line | Path | What |
|------|------|------|
| 54 | `outputs/res_{5,10,25}kb/` | Input directory (resolution-specific) |
| 55 | `06_counts_matrix.rds` | HDF5-backed count matrix |
| 56 | `03_binned.rds` | Binned loop coordinates |
| 57 | `02_merged.rds` | Merged consensus loops |
| 58 | `qc_report/qc_report_summary.rds` | QC summary |
| 62 | `outputs/edgeR_results_res_*kb/` | Output base dir (logs/RDS) |
| 66 | `hiccups_comparison/` | Comparison output subdir |
| 70 | `logs/` | Log output subdir |

## `apa_analysis.R` (3 TODOs)

| Line | Path | What |
|------|------|------|
| 71 | `.hic` files | Contact matrices (multi-GB, HPC-only) |
| 794 | `outputs/res_*kb/03_binned.rds` | Pipeline RDS intermediate |
| 796 | `primary_analysis/all_results_primary.tsv` | Per-resolution edgeR output |

## `apa_shared_anchors.R` (1 TODO)

| Line | Path | What |
|------|------|------|
| 47 | `config/paths_config.yaml` | HPC config with .hic file paths |

---

## Resolution

These scripts are included in `data/scripts/` for reference but require the full HPC pipeline to execute. The actual analyses were run from the original script locations (`scripts/edgeR.R`, `scripts/apa_analysis.R`, `scripts/apa_shared_anchors.R`).

All other paths across all 38 scripts now point to either:
- `data/tsvs/`, `data/plots/`, `data/upstream/` (bundled in data/)
- Repo-relative paths (e.g., `peaks/beds/`, `outputs/250402-late_outputs/`) marked with `# Note: repo-relative path, not bundled in data/`
