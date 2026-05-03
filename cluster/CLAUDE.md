# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

@docs/CONTEXT-cluster.md
@docs/PLAN-p1.md
@docs/PLAN-p2.md
@docs/PLAN-p3.md
@docs/RESULTS-cluster.md
@README.md

## What This Is

Adaptation of the Popay et al. (Nat Genet 2026) Hi-C loop clustering pipeline for BAP1-KO mouse cerebellum (mm10, late/adult timepoint 250402). Answers whether Polycomb enrichment at differential loop anchors vs. loop spans supports an anchor-disruption or extrusion-impediment model. The key result is ChromHMM state enrichment heatmaps (anchor vs. span) across k-means clusters — the BAP1-KO equivalent of Popay Figure 2f.

The original Popay code (hg38/hTERT-RPE1/NIPBL depletion) was cloned from `github.com/tpopay/HiC-clustering` and adapted in Phase 0. Eight Python modules in `modules/` are the modified Popay library; numbered scripts in `scripts/` are the BAP1-KO pipeline.

## Environment

**Conda env:** `cluster` (Python 3.8.18). Lockfile: `cluster.yml`. Minimal spec: `environment_mac.yml` / `environment_linux.yml`.

**Python interpreter (Mac):** `/opt/homebrew/anaconda3/envs/cluster/bin/python3`
- `conda run -n cluster python3` does NOT reliably activate the env (pandas import fails). Always use the absolute path.
- Python 3.8: PEP 585 generics (`list[X]`, `dict[X,Y]`) are NOT available. Use `from typing import List, Dict, Tuple, Optional`.
- pandas 1.5.3: use `lineterminator='\n'` (single word), not deprecated `line_terminator`.

**R (Phase 1 only):** system R 4.5.2 invoked as `/usr/local/bin/Rscript`. Bioconductor packages in `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library`. `keepSeqlevels()` is in `GenomeInfoDb`, not `GenomicRanges`.

**External binaries on PATH:**
- `cluster` — Cluster 3.0 at `/usr/local/bin/cluster` (k-means engine)
- `chromhmm` — wrapper in `ChromHMM/` dir, added to PATH via `~/.zshrc`
- `bedtools` — v2.31.1

Key Python packages: matplotlib 3.7.5, seaborn 0.13.2, scikit-learn 1.3.2, scipy 1.10.1, bioframe 0.6.1, pyBigWig, cooler 0.9.3, cooltools 0.6.1, deeptools 3.5.5.

## Running the Pipeline

All commands run from `cluster/` (the working directory the runners cd to via `$(dirname "$0")/..`). Invoking via `bash cluster/scripts/run_phaseN.sh` from the repo root works too — runners always end up cwd'd to `cluster/` regardless. Phase drivers do NOT use `set -euo pipefail` (benign Java/awk warnings would abort); they use explicit `$?` checks. Logs default to `docs/phaseN.txt` (relative to `cluster/`, i.e. `cluster/docs/phaseN.txt`).

```bash
cd cluster   # (optional — runners cd themselves)

# Phase 1: Data prep (~1 min)
bash scripts/run_phase1.sh

# Phase 2: ChromHMM segmentation (30–90 min for LearnModel)
bash scripts/run_phase2.sh
# MANUAL STEP: edit outputs/bap1_late/chromHMM/12state_rename_cerebellum.txt after

# Phase 3: K-means clustering (~2 min)
bash scripts/run_phase3.sh [K] [FILTER_PCT] [MIN_RATIO] [MAX_RATIO]
# Canonical: bash scripts/run_phase3.sh 6 0.01 0.333 3.0

# Phase 4: Downstream analyses (~6 min)
bash scripts/run_phase4.sh [all|4.1,4.4,...]
# Run only the key result: bash scripts/run_phase4.sh 4.4

# Phase 5: deepTools metagene (90–120 min)
bash scripts/run_phase5.sh

# Phase 6: Cooltools pileup (HPC only — needs mcool files on Expanse)
sbatch scripts/07_cooltools_pileup.sb        # MUST run cooler_balance.sb first if mcool only has KR

# Phase 7: Lab-meeting summary figures (~1 min)
bash scripts/run_summary.sh

# Phase 8: Oriented anchor metagene + asymmetry (~1.7 h for step 1)
bash scripts/run_oriented_metagene.sh
/opt/homebrew/anaconda3/envs/cluster/bin/python3 scripts/quantify_orientation_asymmetry.py
/opt/homebrew/anaconda3/envs/cluster/bin/python3 scripts/visualize_orientation_asymmetry.py

# Phase 9: Clust6 subgroup asymmetry (short vs long, ~2-4 min)
bash scripts/run_clust6_subgroups.sh

# Phase 10: Histone anchors metagene profile regeneration (~1-2 min)
/opt/homebrew/anaconda3/envs/cluster/bin/python3 scripts/11_histone_anchors_metagene.py

# Phase 11: Comprehensive asymmetry — H2AK119ub, H3K27ac, PC1, insulation (HPC only)
sbatch scripts/12_comprehensive_asymmetry.sb

# Phase 2b: 9-mark ChromHMM expansion (HPC for segmentation, Mac for downstream)
sbatch scripts/03b_chromhmm_9mark.sb intersect   # or: union, both
sbatch scripts/03b_chromhmm_9mark.sb union
# MANUAL: write 18state_rename_cerebellum.txt in chromHMM_9mark_{intersect,union}/
bash scripts/run_phase4_9mark.sh intersect 18     # Phase 4.4+4.5 for 9-mark model
```

## Pipeline Architecture

```
Phase 0: Bug fixes in inherited Popay modules (DONE)
Phase 1: 01_build_loop_count_file.py  →  data/late_merged_loop_counts.txt (39,344 loops)
         02_build_mm10_gene_annotation.R  →  data/mm10_knownGene_pp.bed
              ┌──────────┴──────────┐
Phase 2: ChromHMM               Phase 3: K-means
  03_chromhmm_segmentation.sh     04_clustering.py
  (BinarizeBed → LearnModel)      (elbow → Cluster3.0 → sort)
              └──────────┬──────────┘
Phase 4: 05_grouped_analyses.py (8 sub-analyses)
  4.1 loop size  │  4.2 loop classification  │  4.3 anchor ChIP signal
  4.4 ChromHMM anchor vs span (KEY RESULT)   │  4.5 ChromHMM proportions
  4.6 gene annotation  │  4.7 DiffBind       │  4.8 cluster × differential
Phase 5: 06_deeptools_metagene.py (anchor ±5kb per cluster × 8 BigWigs)
Phase 6: 07_cooltools_pileup.py (off-diagonal Hi-C pileup — HPC only;
         requires cooler_balance.sb prereq if mcool lacks `weight` column)
Phase 7: 08_summary_figures.py (3 composite lab-meeting figures —
         dashboard / mechanism / heatmap)
Phase 8: 09_oriented_anchor_metagene.py (strand-aware anchor metagene)
       + quantify_orientation_asymmetry.py  (Wilcoxon ext/int → TSV)
       + visualize_orientation_asymmetry.py (K27me3 dual-panel figure)
Phase 9: 10_clust6_subgroup_asymmetry.py (clust6 short/long split + asymmetry)
Phase 10: 11_histone_anchors_metagene.py (clean profile figure from Phase 5 matrix)
Phase 11: 12_comprehensive_asymmetry.py (H2AK119ub, H3K27ac, PC1, insulation
           asymmetry — HPC only, needs mcools for cooltools eigs-cis + insulation)
Phase 2b: 03b_chromhmm_9mark_segmentation.sh (9-mark ChromHMM, 15+18 states)
        + run_phase4_9mark.sh (Phase 4.4/4.5 rerun with 9-mark env vars)
```

Phases 2 and 3 are independent of each other (both depend on Phase 1). Phase 4 requires both 2 and 3. Phases 5–11 all depend on Phase 3's clustering output. Phase 11 additionally requires mcools on HPC.

## Module Architecture

Eight `.py` files in `modules/` are the inherited Popay library (modified in Phase 0). Scripts in `scripts/` add `cluster/modules/` to `sys.path` and call into them:

| Module | Purpose | Key functions |
|--------|---------|---------------|
| `plotting.py` | All statistical figure types | `heat()`, `line()`, `box()`, `strip()`, `stacked()`, `bar()`, `joint()` |
| `cluster_tools.py` | Clustering utilities | `elbow()`, `sort_clusters()`, `sort_by_strings()` |
| `statistics_functions.py` | Non-parametric tests | `kruskal_wilcoxon()` — Kruskal-Wallis + pairwise Wilcoxon with Bonferroni |
| `chromHMM_heatmap.py` | ChromHMM enrichment heatmaps | `heatmap_plot()` — output lands adjacent to input file, not configurable |
| `deeptools_plotting.py` | deepTools matrix visualization | `heatmap_plot()` — line + heatmap from `computeMatrix` output |
| `deepTools_pipeline.py` | deepTools subprocess wrapper | `bed_pileup()` — shells out to `computeMatrix` |
| `bedpe_analysis.py` | BEDPE-level analyses | `loop_size()`, `bedtools_annotation()` |
| `cooltools_called.py` | Cooltools pileup wrapper | `mcool_pileup()` |

All modules load `modules/custom_params.json` at import time (matplotlib rcParams for publication figures: Arial, small labels, TrueType embedding). `plotting.py` uses `mpl_toolkits.axes_grid1.Divider` for fixed-pixel axis sizing with dimensions computed from label text width via PIL `ImageFont.textbbox`. Scripts that need `custom_params.json` directly (e.g. `08_summary_figures.py`, `visualize_orientation_asymmetry.py`) reference `CLUSTER_DIR / 'modules' / 'custom_params.json'`.

`scripts/utils/multi_format_output.py` provides a `multi_format_savefig()` context manager that monkey-patches `Figure.savefig` to auto-emit `.svg`, `.pdf`, `.jpg` alongside any `.png` save.

## Column Naming Convention

The loop count file uses `ctrl_merge` / `mut_merge` (NOT `_merged`). Downstream code does `treatment_group.str.replace('_merge', '')` to get display labels `ctrl` / `mut`. Using `_merged` would produce broken labels `ctrld` / `mutd`.

## Data Locations

**In-repo (`cluster/data/`):**
- `late_merged_loop_counts.txt` — 39,344 × 8 (chr1, x1, x2, chr2, y1, y2, ctrl_merge, mut_merge)
- `late_merged_loop_metadata.tsv` — 39,344 × 16 (coords + edgeR logFC/FDR/direction + resolution info)
- `mm10_knownGene_pp.bed` — 24,515 promoter regions (TSS ± 750bp)

**External (not in repo):**
- BigWigs: `/Users/zakiralibhai/sdsc/bigwigs/{H3K27ac,H3K27me3,H2AK119ub,H3K27me1}{Ctrl,Mut}.bw`
  - Do NOT use `peaks/bigwigs/macs2.narrow.aug18.dedup/` — has 0-byte mutant files
- ChIP peaks: `peaks/CTCF.bed`, `peaks/beds/H3K27acCerebellumLate2.bed`, etc.
- DiffBind: `peaks/diffbind/{K27ac,K27me3,K119ub}_diffbind_results_summit_appended_ap.txt`
- Blacklist: `tads/mm10-blacklist.v2.bed`
- mcools (HPC only): `/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/data/cool/250402/{ctrl,mut}_merged.mcool`

**Outputs (`cluster/outputs/bap1_late/`):**
- `cluster3/k-6/data/combined-clusters.txt` — canonical clustering (38,948 loops × 6 clusters)
- `chromHMM/learned_model/cerebellum_late_12_segments.bed` — 12-state segmentation
- `chromHMM/12state_rename_cerebellum.txt` — manual E1–E12 → biological name mapping
- `chromHMM/{anchor,span}.txt` — OverlapEnrichment fold-enrichment tables
- `figures/` — all Phase 4/5/7/8 figure outputs organized by analysis type
- `cooltools/obs_exp_contacts/obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}` — Phase 6 pileup outputs

**Logs (`cluster/docs/`):** captured stdout/stderr from each phase driver; the canonical log directory after the 2026-04-27 reorg. SLURM logs from HPC are in `cluster/logs/`.

## Canonical Clustering (k=6)

Clusters sorted by descending mean mut/ctrl signal. Biological ordering for summary figures: `[clust6, clust3, clust1, clust2, clust4, clust5]` (loss → unchanged → gain).

| Cluster | n | Direction | Key biology |
|---------|---|-----------|-------------|
| clust1 | 12,298 | 100% unchanged | High-signal unchanged |
| clust2 | 10,970 | 92% unchanged | Moderate unchanged |
| clust3 | 8,738 | 79% unchanged, 21% down | Moderate loss |
| clust4 | 3,916 | 70% up | Moderate gain |
| clust5 | 667 | 97% up | Strong gain — Polycomb domain compaction (anchor 6.59× / span 3.03×) |
| clust6 | 2,359 | 78% down | Strong loss — anchor disruption (anchor 2.09× / span 0.94×) |

## ChromHMM 12-State Model (5 marks, original)

| ID | State | Notes |
|----|-------|-------|
| E1 | Quiescent | 92% of genome; excluded from proportion plots |
| E2 | CTCF_Boundary | |
| E3 | Poised_Enhancer | |
| E4 | Active_Enhancer | |
| E5 | Active_Promoter_Flank | |
| E6 | Poised_Promoter | |
| E7 | Weak_Promoter | |
| E8 | Active_Promoter | |
| E9 | Strong_Enhancer | |
| E10 | Insulator | |
| E11 | Polycomb | |
| E12 | Bivalent_Enhancer | |

## ChromHMM 18-State Model (9 marks, expanded)

9-mark expansion: original 5 + H2AK119ub, ATAC, H3K9ac, H3K9me3. Intersect consensus for H3K9 replicates. H3K9me3 is early-only (limitation). Outputs under `outputs/bap1_late/chromHMM_9mark_intersect/`.

| ID | State | Marks ON |
|----|-------|----------|
| E1 | K119ub_Only | K119ub |
| E2 | Polycomb_K119ub | K119ub + K27me3 |
| E3 | Repressed_Enhancer_K119ub | K4me1 + K119ub + K27me3 |
| E4 | K119ub_Poised_Enhancer | K4me1 + K119ub |
| E5 | Poised_Enhancer | K4me1 |
| E6 | Active_Enhancer_K9ac | K4me1 + K9ac |
| E7 | Active_Enhancer | K4me1 + K9ac + K27ac |
| E8 | Active_Enhancer_K119ub | K4me1 + K9ac + K27ac + K119ub |
| E9 | Strong_Enhancer | K4me1 + K27ac |
| E10 | Weak_Enhancer | K27ac |
| E11 | Active_Promoter | K9ac + K4me3 + K27ac |
| E12 | K9ac_Promoter | K9ac |
| E13 | Quiescent | (none) |
| E14 | CTCF_Open | ATAC |
| E15 | ATAC_Enhancer | ATAC + K4me1 |
| E16 | Insulator | CTCF |
| E17 | Polycomb | K27me3 |
| E18 | Heterochromatin | K9me3 |

`05_grouped_analyses.py` is parameterizable via env vars to use either model:
- `CLUSTER_NSTATES` (default `12`), `CLUSTER_CHROMHMM_SUBDIR` (default `chromHMM`), `CLUSTER_MODEL_SUBDIR` (default `learned_model`), `CLUSTER_RENAME_SUFFIX` (default `cerebellum`), `CLUSTER_ENRICH_SUFFIX` (default `''`).
- 9-mark invocation: `bash scripts/run_phase4_9mark.sh intersect 18`

## Critical Gotchas

**Runtime:**
- Phase 5 `computeMatrix` takes 90–120 min over 63k anchors × 8 BigWigs. Set `PYTHONUNBUFFERED=1` or the log appears hung (block-buffered stdout during subprocess). `run_phase5.sh` handles this.
- Phase 5 intermediates (`histone_anchors` 96MB gz, `histone_anchors_values` 1.6GB) should NOT be committed.
- `run_phase5.sh` prepends cluster env bin to PATH so `computeMatrix` subprocess can find deepTools executables.

**API quirks in Popay modules:**
- `chromHMM_heatmap.heatmap_plot(path)` saves output adjacent to `path` (e.g., `dirname(path)/stem.png`), not to a configurable output directory.
- `plotting.stacked()` requires `df.index.name = None` — `pd.crosstab` sets a named index that breaks the internal `reset_index().rename({'index': 'xcol'})`.
- `bedpe_analysis.bedtools_annotation()` requires non-None `temp_dir` even when `FPKM_df=None`.
- `cooltools_called.mcool_pileup()` expects BEDPE DataFrames with columns `chrom1/start1/end1/chrom2/start2/end2` (not the project's `chr1/x1/x2`).
- `bioframe.fetch_centromeres('mm10')` fails on bioframe 0.6.1. `cooltools_called.mcool_pileup` builds the viewframe inline from `cooler.Cooler(mcool).chromsizes` (avoids the bioframe call entirely). The earlier `make_arms_robust` monkey-patch in `07_cooltools_pileup.py` was superseded — the script no longer monkey-patches anything.
- `cooltools.expected_cis` requires the cooler `weight` column from `cooler balance`. Source mcools with only `KR` weights from juicer-pre will fail with `cooler is not balanced`. Run `cluster/scripts/cooler_balance.sb` first.

**Data integrity:**
- All 38,948 retained loops are cis (`chr1 == chr2`). Span BED construction (`df[['chr1', 'x1', 'y2']]`) is only valid for cis loops — an assertion guards this.
- DiffBind files use `Peak_Chr/Peak_Start/Peak_End`; Popay code expects `Chr/Start/End`. Renamed at load time.
- Filter threshold for Phase 3 must be percentile-based (not absolute) because per-resolution count scales differ ~6× (5kb median 171, 10kb 403, 25kb 984).

## Reference Documents

- `README.md` — project overview, quick start, key result, output structure
- `CLAUDE.md` — this file (AI assistant context, environment, gotchas, module API)

All planning, results, and log files live under `cluster/docs/` (post-2026-04-27 reorg):

- `docs/CONTEXT-CLUSTER.md` — biological context, meeting notes, Popay paper summary, data inventory
- `docs/RESULTS-cluster.md` — statistical outputs by phase, biological interpretation, open questions
- `docs/plan-p1.md` — session-level plan for Phases 0–3 with corrections and verification
- `docs/plan-p2.md` — session-level plan for Phases 4–6 with corrections and verification
- `docs/plan-p3.md` — session-level plan for Phases 7–8 with corrections and verification
- `docs/phaseN.txt` — captured stdout/stderr from each phase driver execution
- `example/` — Popay's original hg38 example files (reference schemas)
