# Stripenn Pipeline Implementation Plan (Stages 1-6 + Orchestration)

## Context

The BAP1-KO differential chromatin stripe analysis pipeline needs a parallel Stripenn-based track alongside the existing Quagga pipeline. Stage 0 (`.hic` -> `.mcool` conversion via hictk) is **complete**. The mcool files are at `${DATA_DIR}/data/cool/{250831,250402}/{sample}.mcool` with 5kb and 10kb resolutions.

This plan implements **Stages 1-6** plus the master orchestration driver — **16 new files** total — following patterns from the existing Quagga pipeline scripts and the completed Stage 0 scripts.

---

## Critical Findings from Exploration

1. **Output filename is `result_filtered.tsv`** (not `.txt` as the README says) — confirmed from `stripenn.py:156`
2. **`stripenn score` has NO `-h`/`--header` flag** — the PLAN.md's score command is wrong. `score.py:11-21` auto-detects headers by checking if the first row is numeric vs string
3. **Score output adds 6 columns**: `pvalue_added, Stripiness_added, O_Mean_added, O_Sum_added, O/E_Mean_added, O/E_Total_added` — `O_Sum_added` is the edgeR count input
4. **`--seed 42`** must be passed explicitly (stripenn default is 123456789, config says 42)
5. **Chromosome filtering** must be post-hoc in R — stripenn score.py only filters chrY and scaffolds, not chrX/chrM

---

## Implementation Order & Files

All files go in: `stripes/stripenn/scripts/`

### Stage 1: Call Stripes on Merged Files (2 files)

**`01_call_stripes.sb`** — SLURM script, CLI: `<timepoint> <resolution_bp> <sample>`
- Template: `00_convert_hic_to_cool.sb` (same header structure, arg validation, CODE_DIR/DATA_DIR)
- Uses `mariner_env` (not hictk)
- SLURM: 16 cpus, 96G mem, 8h (from config `stage_01_call`)
- Key command:
  ```bash
  yes Y | stripenn compute \
    --cool "${MCOOL}::/resolutions/${RES}" \
    --out "${OUT_DIR}" \
    --norm KR -k all \
    -m 0.95,0.96,0.97,0.98,0.99 \
    -p 0.1 -n ${SLURM_CPUS_PER_TASK} \
    --seed 42 -l 10 -w 8 -b 3
  ```
- Input: `${DATA_DIR}/data/cool/${TP}/${SAMPLE}.mcool`
- Output: `${DATA_DIR}/outputs/${TP}/res_${RES_KB}kb/calls/${SAMPLE}/result_filtered.tsv`
- Verify: output exists and is non-empty, count stripes

**`submit_01_call.sh`** — Wrapper, submits 8 jobs (2tp x 2res x 2 merged)
- Template: `submit_00_convert.sh` (same dependency pass-through, sbatch --parsable)
- Samples: `ctrl_merged`, `mut_merged` only

### Stage 2: Build Union Set (3 files)

**`02_build_union.R`** — R script, CLI: `Rscript 02_build_union.R <timepoint> <resolution_bp>`
- **Port from**: `stripes/scripts/phase1_detection.R` — reuse these functions:
  - `stripes_to_granges()` (lines 127-139) — adapt column names for stripenn format
  - `match_anchors()` (lines 144-166) — reuse verbatim
  - `build_union_set()` (lines 171-249) — adapt for stripenn's 12-column output, add Stripiness/Mean columns
- **Key differences from Quagga phase1**:
  - Input format: stripenn's `result_filtered.tsv` has header + 12 cols (`chr,pos1,pos2,chr2,pos3,pos4,length,width,Mean,maxpixel,pvalue,Stripiness`)
  - No replicate support counting (replicates scored in Stage 3)
  - No 10kb validation (done in Stage 6)
  - Simpler output: just union coordinates + source + detection stats
- Chromosome filtering: exclude chrX, chrY, chrM from config
- Outputs:
  - `union_stripes.tsv` — full annotations (stripe_id, chr, pos1-pos4, direction, source, pval, Stripiness, etc.)
  - `union_stripes.bedpe` — 6-column with header (`chr\tpos1\tpos2\tchr2\tpos3\tpos4`) for `stripenn score` input
- Config: loaded from `${CODE_DIR}/config/stripenn_config.yaml` (absolute path)

**`02_build_union.sb`** — SLURM wrapper, 4 cpus, 16G, 1h
- Validates Stage 1 outputs exist before running R

**`submit_02_union.sh`** — Submits 4 jobs (2tp x 2res)

### Stage 3: Score Replicates (2 files)

**`03_score_replicates.sb`** — SLURM script, CLI: `<timepoint> <resolution_bp> <sample>`
- SLURM: 12 cpus, 64G mem, 4h (from config `stage_03_score`)
- Key command (**NO `-h` flag** — score auto-detects headers):
  ```bash
  stripenn score \
    --cool "${MCOOL}::/resolutions/${RES}" \
    -c "${UNION_BEDPE}" \
    --norm KR \
    -n ${SLURM_CPUS_PER_TASK} \
    --seed 42 \
    -o "${OUT}"
  ```
- Input: `${DATA_DIR}/data/cool/${TP}/${SAMPLE}.mcool` + `union_stripes.bedpe`
- Output: `${DATA_DIR}/outputs/${TP}/res_${RES_KB}kb/scores/${SAMPLE}.scores.tsv`
- Verify: output row count matches union BEDPE row count

**`submit_03_score.sh`** — Submits 24 jobs (2tp x 2res x 6 replicates)
- Samples: `ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3`

### Stage 4: edgeR Differential Analysis (3 files)

**`04_edgeR.R`** — R script, CLI: `Rscript 04_edgeR.R <timepoint> <resolution_bp>`
- **Port from**: `stripes/scripts/phase3_edgeR.R` — reuse:
  - `save_multiformat()` helper (lines 20-31)
  - Full edgeR workflow (skip filtering, TMM, robust QL-GLM)
  - Diagnostic plots: MDS, BCV, QL dispersion, volcano
- **Data loading** (new for stripenn):
  ```r
  # Read 6 score files, extract O_Sum_added, build count matrix
  for (sample in config$samples$replicates) {
    score_file <- file.path(scores_dir, paste0(sample, ".scores.tsv"))
    scores <- read.delim(score_file)
    count_matrix[, sample] <- round(scores$O_Sum_added)
  }
  ```
- Chromosome filtering: remove rows where chr in exclude_chromosomes
- `O_Sum_added` rounded to integer for edgeR
- Skip `filterByExpr()` (small feature set, config `skip_filtering: true`)
- Outputs in `${DATA_DIR}/outputs/${TP}/res_${RES_KB}kb/04_edgeR/`:
  - `all_results.tsv` + `.rds`
  - `dge_object.rds`, `fit_object.rds`
  - `plots/` — mds, bcv, ql_dispersion, volcano (PDF/SVG/JPEG)
  - `summary.txt`

**`04_edgeR.sb`** — SLURM wrapper, 4 cpus, 16G, 1h

**`submit_04_edgeR.sh`** — Submits 4 jobs (2tp x 2res)

### Stage 5: Integration & Classification (3 files)

**`05_integration.R`** — R script, CLI: `Rscript 05_integration.R <timepoint> <resolution_bp>`
- **Port from**: `stripes/scripts/phase4_integration.R` — reuse classification logic
- Merge `union_stripes.tsv` (source, pval, Stripiness) + `all_results.tsv` (logFC, FDR) by stripe_id
- Direction classification:
  - control_only -> "lost", mutant_only -> "gained"
  - shared + FDR<0.05 + logFC>0.3 -> "strengthened"
  - shared + FDR<0.05 + logFC<-0.3 -> "weakened"
  - shared otherwise -> "unchanged"
- Confidence tiers (high/medium/low) based on FDR + logFC agreement with source
- Directional consistency flag
- Outputs:
  - `05_final_differential.tsv` + `.rds` + `.bedpe` (Juicebox-compatible)
  - Per-category: `05_stripes_lost.tsv`, `05_stripes_gained.tsv`, etc.
  - `05_summary.txt`

**`05_integration.sb`** — SLURM wrapper, 4 cpus, 16G, 1h

**`submit_05_integration.sh`** — Submits 4 jobs (2tp x 2res)

### Stage 6: Cross-Resolution Comparison (3 files)

**`06_compare_resolutions.R`** — R script, CLI: `Rscript 06_compare_resolutions.R <timepoint>`
- **Port from**: `scripts/compare_resolutions.R` — simplified for 2 resolutions (5kb, 10kb)
- Load `05_final_differential.tsv` from both `res_5kb/` and `res_10kb/`
- Anchor-based matching (50kb tolerance) between resolutions — reuse `match_anchors()` from Stage 2
- Classify:
  - **High-confidence**: significant at BOTH resolutions with concordant direction
  - **Exploratory**: significant at only one resolution
- Plots: 2-set Venn, logFC correlation scatter, bar chart
- Output: `${DATA_DIR}/outputs/${TP}/cross_res_merged.tsv` + plots

**`06_compare_resolutions.sb`** — SLURM wrapper, 4 cpus, 16G, 1h

**`submit_06_compare.sh`** — Submits 2 jobs (2tp)

### Orchestration (1 file)

**`run_full_stripenn.sh`** — Master driver (bash, NOT .sb — runs from login node)
- Chains all stages with `--dependency=afterok:...`
- Stage 0 can be included (already done, will skip via idempotent check) or skipped
- Pattern: capture job IDs from each `submit_*.sh`, pass to next stage
- Prints progress to stderr, job IDs to stdout

---

## Key Patterns to Reuse (with file paths)

| Pattern | Source File | Lines | Target |
|---------|------------|-------|--------|
| SLURM header + arg validation | `scripts/00_convert_hic_to_cool.sb` | 1-52 | All `.sb` files |
| Wrapper dependency pass-through | `scripts/submit_00_convert.sh` | 18-22 | All `submit_*.sh` |
| Direction inference | `stripes/scripts/phase1_detection.R` | 104-116 | `02_build_union.R` |
| `stripes_to_granges()` | `stripes/scripts/phase1_detection.R` | 127-139 | `02_build_union.R` |
| `match_anchors()` | `stripes/scripts/phase1_detection.R` | 144-166 | `02_build_union.R`, `06_compare.R` |
| `build_union_set()` | `stripes/scripts/phase1_detection.R` | 171-249 | `02_build_union.R` |
| `save_multiformat()` | `stripes/scripts/phase3_edgeR.R` | 20-31 | `04_edgeR.R` |
| edgeR QL-GLM workflow | `stripes/scripts/phase3_edgeR.R` | 149-287 | `04_edgeR.R` |
| Tiered classification | `stripes/scripts/phase4_integration.R` | 139-185 | `05_integration.R` |
| Cross-res matching | `scripts/compare_resolutions.R` | 237-315 | `06_compare_resolutions.R` |

---

## Verification Plan

| Stage | Check | Expected |
|-------|-------|----------|
| 1 | `wc -l result_filtered.tsv` per call | 100-500 stripes per merged sample |
| 2 | `wc -l union_stripes.bedpe` | >= max(ctrl, mut) count (union) |
| 2 | source distribution | shared + control_only + mutant_only = total |
| 3 | score file row count | matches union BEDPE (minus header) |
| 3 | `O_Sum_added` values | non-negative, mostly > 0 |
| 4 | MDS plot | ctrl/mut cluster separation |
| 4 | BCV | 0.2-0.5 range |
| 5 | source x direction cross-tab | control_only all "lost", mutant_only all "gained" |
| 6 | high-confidence subset | proper subset of all significant |

---

## Implementation Sequence

Build and test each stage sequentially — each depends on the previous stage's output:

1. **Stage 1**: `01_call_stripes.sb` + `submit_01_call.sh` (simplest — just wraps a CLI)
2. **Stage 2**: `02_build_union.R` + `.sb` + `submit` (core matching logic — most complex R script)
3. **Stage 3**: `03_score_replicates.sb` + `submit` (simple CLI wrapper)
4. **Stage 4**: `04_edgeR.R` + `.sb` + `submit` (closest port from existing code)
5. **Stage 5**: `05_integration.R` + `.sb` + `submit` (classification logic port)
6. **Stage 6**: `06_compare_resolutions.R` + `.sb` + `submit` (last analysis stage)
7. **Orchestration**: `run_full_stripenn.sh` (write last, test with manual stage submission first)

For each stage, write all files (`.sb`, `.R` if applicable, `submit_*.sh`), then the user can test on HPC before moving to the next stage.
