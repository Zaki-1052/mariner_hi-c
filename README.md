# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

BAP1-KO mouse cerebellum (mm10) 3D genome analysis project. Studies how loss of the H2AK119ub deubiquitinase BAP1 disrupts chromatin architecture — loops, stripes, TADs, compartments, enhancer-gene linkage, and DNA methylation — across two developmental timepoints (P12 early/250831, adult late/250402). Three biological replicates per condition (mutant vs wildtype control).

The central biological finding: BAP1 loss causes two simultaneous mechanisms — Polycomb-domain compaction (gained loops) and anchor disruption (lost loops) — with H2AK119ub as the upstream signal.

## Repository Layout

Each subdirectory is a self-contained analysis module with its own `CLAUDE.md` and/or `README.md`. Read those before working in a submodule.

| Directory | What | Languages | Compute |
|-----------|------|-----------|---------|
| `loops/` | **Core pipeline.** Mariner-based differential loop calling (edgeR QLF-GLM, 3 resolutions). Entry point: `sbatch scripts/run_full_pipeline.sb` or per-step `Rscript scripts/{step}.R {resolution}` | R, Python, Bash | HPC (SLURM) |
| `cluster/` | Popay-style k-means loop clustering + ChromHMM anchor-vs-span enrichment. 14 phases. Entry point: `bash scripts/run_phase{N}.sh` | Python 3.8 (cluster conda env), R, Bash | Mac + HPC |
| `stripes/stripenn/` | Stripenn differential stripe analysis (8 stages). Entry point: `bash scripts/run_full_stripenn.sh` | Python, R, Bash | HPC |
| `stripes/quagga/` | Quagga stripe pipeline (alternative caller) | R, Python, Bash | HPC |
| `tads/` | TADCompare differential TAD boundary analysis. Entry point: `sbatch scripts/run_full_pipeline.sb` | R, Bash | HPC |
| `peaks/` | ChIP-seq peak annotation of loop anchors (8-category chromatin state system) | R | Local |
| `abc/` | ABC (Activity-By-Contact) enhancer-gene linkage via Broad Institute Snakemake pipeline + delta-ABC | Python, R, Bash | HPC + Local |
| `biomodal/` | Biomodal DUET evoC differential methylation (5mC/5hmC). Upstream: Nextflow. Downstream: modality XPLR + R visualization | Nextflow, R, Bash | HPC |
| `ML/cmpts/` | SNIPER (hg19-only, retrained) + CALDER2 (primary for mm10) subcompartment analysis | Python, R | HPC |
| `ML/loops/` | LoopBin unsupervised deep learning loop clustering (VADE model) | Python (TensorFlow) | HPC |
| `data/` | Shared upstream data files, TSVs, plotting scripts | R | Local |
| `poster/` | Conference poster figures | R | Local |
| `lab/` | Lab utilities — CUT&RUN pipeline, .hic file auditor, archival tools | Bash, Python | Mixed |
| `notes/` | Obsidian vault, meeting notes (gitignored) | Markdown | — |
| `figs/` | Static figures, AI-generated art | — | — |

## Conda Environments

| Env | Used by | Key packages |
|-----|---------|-------------|
| `cluster` | `cluster/` | Python 3.8, pandas 1.5.3, matplotlib 3.7.5, cooltools 0.6.1, deeptools 3.5.5, bioframe 0.6.1 |
| `mariner_env` | `loops/`, `stripes/stripenn/` (stages 1-7) | R + Bioconductor (mariner, edgeR, InteractionSet), Python 3.13, stripenn (patched) |
| `hictk` | `stripes/stripenn/` (stage 0 only) | hictk for .hic v9 → .mcool conversion |
| `env_nf` | `biomodal/` upstream | Nextflow, biomodal DUET pipeline |
| `modality` | `biomodal/` downstream DMR calling | modality XPLR CLI |
| `abc-env` | `abc/` | Snakemake, bedtools, MACS2, samtools |

System R 4.5.2 (`/usr/local/bin/Rscript`) is used by `cluster/` Phase 1 and various local R scripts. `conda run -n cluster python3` does NOT reliably activate the cluster env — use `/opt/homebrew/anaconda3/envs/cluster/bin/python3` directly.

## HPC (SDSC Expanse)

SLURM account: `csd940`, partition: `shared`. Two-root path convention used by `stripes/stripenn/`:

| Root | Path | Synced |
|------|------|--------|
| CODE | `/expanse/.../mariner_hi-c/stripes/stripenn` | Yes (GitHub) |
| DATA | `/expanse/.../stripes/stripenn` | No (large files) |

.hic files: `/expanse/.../stripes/StripeCaller/data/hic/{250402,250831}/`
mcool files: `/expanse/.../stripes/stripenn/data/cool/{250402,250831}/`

## Key External Data Paths (Mac)

- BigWigs: `/Users/zakiralibhai/sdsc/bigwigs/` — canonical source for H3K27ac, H3K27me3, H2AK119ub, H3K27me1 (ctrl + mut). Do NOT use `peaks/bigwigs/macs2.narrow.aug18.dedup/` — has 0-byte mutant files.
- ChIP peak BEDs: `peaks/beds/` (H3K27ac, H3K27me3, H3K4me1, H3K4me3) and `peaks/CTCF.bed`
- DiffBind results: `peaks/diffbind/` (K27ac, K27me3, K119ub)
- Blacklist: `tads/mm10-blacklist.v2.bed`
- DEG lists: `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`, `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx`

## External Tool Dependencies

| Tool | Path / Install | Used by |
|------|---------------|---------|
| ChromHMM v1.27 | `cluster/ChromHMM/ChromHMM.jar`, wrapper on PATH via `~/.zshrc` | `cluster/` Phase 2 |
| Cluster 3.0 | `/usr/local/bin/cluster` | `cluster/` Phase 3 (k-means) |
| bedtools 2.31.1 | `/opt/homebrew/bin/bedtools` | Multiple modules |
| Java 25 | System | ChromHMM |
| hictk | conda `hictk` env | .hic → .mcool conversion |

## Multi-Format Figure Convention

All figure-generating scripts save in 4 formats (PNG + PDF + SVG + JPG) in per-figure subfolders. Python scripts use `cluster/scripts/utils/multi_format_output.py` (`multi_format_savefig()` context manager). R scripts use `scripts/utils/multi_format_output.R` (`save_multiformat_ggplot()`).

## Conventions and Gotchas

- **Column naming:** Loop count files use `ctrl_merge` / `mut_merge` (NOT `_merged`). Downstream display labels are derived via `.str.replace('_merge', '')` → `ctrl` / `mut`. Using `_merged` produces broken labels.
- **Python 3.8 (cluster env):** No PEP 585 generics — use `from typing import List, Dict, Tuple, Optional` instead of `list[X]`, `dict[X,Y]`.
- **pandas 1.5.3 (cluster env):** Use `lineterminator='\n'` (one word), not deprecated `line_terminator`.
- **No `set -euo pipefail` in pipeline drivers.** Benign Java/awk warnings would abort the pipeline. Drivers use explicit `$?` checks instead.
- **Stripenn patches required.** Two NumPy 2.0 unsigned integer overflow fixes in `getStripe.py` must be reapplied if stripenn is reinstalled. See `stripes/stripenn/CLAUDE.md` §2.
- **ChromHMM state prefix is `E` (not `U`).** v1.27 emits `E1`–`E12`, not Popay's `U1`–`U12`.
- **Bivalent = Repressed Enhancer.** PI says K4me1+K27me3 = `Repressed_Enhancer`, not `Bivalent_Enhancer` (though the 5-mark 12-state model retains the `Bivalent_Enhancer` name for backward compatibility).
- **deepTools computeMatrix** takes 90-120 min for ~63k anchors × 8 BigWigs. Set `PYTHONUNBUFFERED=1` or logs appear hung. `run_phase5.sh` handles this.
- **Intermediate deepTools files are LARGE.** `*_values` (1.6GB) and `*` gz (96MB) under `cluster/bap1_late/figures/deeptools/` should not be committed.
- **mcools need ICE balancing** before `cooltools.expected_cis` works. The 250402 merged mcools had only KR weights from juicer-pre. Run `cluster/scripts/cooler_balance.sb` first.
