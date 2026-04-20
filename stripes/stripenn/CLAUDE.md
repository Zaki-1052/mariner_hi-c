# Stripenn Differential Stripe Analysis Pipeline: Context Document

## Document Purpose

This document captures all implementation progress, bugs encountered, fixes applied, and remaining work for the Stripenn-based differential chromatin stripe analysis pipeline. It is a companion to the Quagga pipeline (`stripes/quagga/CLAUDE.md`) and serves as context for continuing development across AI sessions.

**Last updated:** 2026-04-20 (initial session)

---

## 1. Why Stripenn

The existing Quagga-based pipeline (`stripes/quagga/`) found ~150-300 stripes per timepoint but very few significant differential calls (FDR<0.05: ~0-10). The PI suggested trying Stripenn as an alternative stripe caller. Stripenn uses image-processing (Canny edge detection, pixel saturation thresholds) rather than Gabor filtering, and includes a built-in `score` function that quantifies stripe signal in any sample against a coordinate set — eliminating the need for mariner-based quantification.

**Key advantage over Quagga pipeline:** Stripenn's `score` command directly quantifies stripes, producing `O_Sum_added` (sum of observed contacts) per stripe per sample. This replaces the Quagga pipeline's Phase 2 (mariner extraction + 5x5 matrix aggregation), making the pipeline simpler and purely Python/R with no Bioconductor dependency for quantification.

---

## 2. Environment Setup

### Two conda environments (important)

| Env | Used by | Why |
|-----|---------|-----|
| `hictk` | Stage 0 only (.hic -> .mcool conversion) | Source .hic files are Juicer format v9; `hic2cool` cannot parse v9 |
| `mariner_env` | Stages 1-6 (stripenn, R, edgeR) | Python 3.13 + stripenn + R packages |

### Stripenn installation in mariner_env

Installed with `pip install --no-deps stripenn` because stripenn's `pandas<2.0.0` pin cannot be satisfied on Python 3.13. Pandas 2.x works in practice for stripenn's actual usage. Runtime deps installed separately: `opencv-python`, `scikit-image`, `typer`, `joblib`, `tqdm`, `matplotlib`.

### Critical patches applied to installed stripenn (site-packages)

Two patches were applied to the installed `getStripe.py` to fix Python 3.13 / NumPy 2.0 incompatibilities:

**Patch 1 — Integer overflow fix (CRITICAL):**
```python
# Line 20: self.resol = resol
# Changed to:
self.resol = int(resol)
```
**Root cause:** `Lib.binsize` (from cooler/h5py) returns a numpy unsigned integer. Under NumPy 2.0's NEP 50 type promotion rules, `Python_int - numpy_uint` wraps around instead of going negative. This corrupts genomic coordinates in `main_null_calc()`, producing regions where start > end, which crashes cooler's `fetch()` with `ValueError: End coordinate less than start`.

**Patch 2 — Chromosome size overflow fix:**
```python
# All occurrences of: chrsize = self.chromnames2sizes[chr]
# Changed to:
chrsize = int(self.chromnames2sizes[chr])
```
Same unsigned integer issue for chromosome sizes used in coordinate arithmetic.

**How to verify patches are applied:**
```bash
conda activate mariner_env
SITE_PKG=$(python -c "import stripenn; import os; print(os.path.dirname(stripenn.__file__))")
grep 'self.resol = ' "${SITE_PKG}/getStripe.py"
# Should show: self.resol = int(resol)
```

**If patches are lost** (e.g., after reinstalling stripenn):
```bash
SITE_PKG=$(python -c "import stripenn; import os; print(os.path.dirname(stripenn.__file__))")
sed -i 's/        self.resol = resol/        self.resol = int(resol)/' "${SITE_PKG}/getStripe.py"
sed -i 's/chrsize = self.chromnames2sizes\[chr\]/chrsize = int(self.chromnames2sizes[chr])/' "${SITE_PKG}/getStripe.py"
```

### FutureWarning suppression

Pandas deprecation warnings (`Series.__getitem__ treating keys as positions`) are silenced in `01_call_stripes.sb` and `03_score_replicates.sb` via:
```bash
export PYTHONWARNINGS="ignore::FutureWarning"
```
These warnings are cosmetic and do not affect results.

---

## 3. Path Convention

Every script uses absolute paths. Two root directories:

| Root | Path | Synced to GitHub? |
|------|------|-------------------|
| CODE_DIR | `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn` | Yes (scripts, config, small TSVs) |
| DATA_DIR | `/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn` | No (mcools, outputs, logs) |

Source .hic files: `/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/{250831,250402}/`

---

## 4. Pipeline Architecture

### 4.1 Overview

7 stages, 19 scripts total. Each stage has a SLURM `.sb` script (accepts CLI args), an optional R script, and a `submit_*.sh` wrapper that iterates over the parameter space and supports `--dependency=afterok:...` chaining.

```
Stage 0: hictk convert     (.hic -> .mcool)           16 jobs  [COMPLETE]
Stage 1: stripenn compute   (call stripes on merged)    8 jobs
Stage 2: 02_build_union.R  (merge ctrl/mut calls)      4 jobs
Stage 3: stripenn score     (score 6 replicates)       24 jobs
Stage 4: 04_edgeR.R        (differential analysis)     4 jobs
Stage 5: 05_integration.R  (classify direction)        4 jobs
Stage 6: 06_compare_res.R  (5kb vs 10kb comparison)    2 jobs
```

**Master orchestration:** `run_full_stripenn.sh [--skip-stage0]` chains all stages via SLURM dependencies.

### 4.2 Stage Details

**Stage 0 — hictk convert** (COMPLETE)
- Script: `00_convert_hic_to_cool.sb` + `submit_00_convert.sh`
- Converts .hic to .mcool with `--resolutions 5000 10000`
- Uses `hictk` conda env
- Idempotent: skips if output exists and is non-empty
- Output: `${DATA_DIR}/data/cool/{tp}/{sample}.mcool`

**Stage 1 — stripenn compute** (stripe detection)
- Script: `01_call_stripes.sb` + `submit_01_call.sh`
- Runs on merged files only (ctrl_merged, mut_merged) at both resolutions
- `yes Y |` prefix handles interactive overwrite prompt
- `--seed 42` (not stripenn's default 123456789)
- SLURM: 16 cpus, 96G, 8h
- Output: `${DATA_DIR}/outputs/{tp}/res_{kb}kb/calls/{sample}/result_filtered.tsv`
- **Output is `result_filtered.tsv`** (NOT `.txt` despite what the stripenn README says)

**Stage 2 — Build union set** (anchor matching)
- Scripts: `02_build_union.R` + `02_build_union.sb` + `submit_02_union.sh`
- Loads ctrl + mut `result_filtered.tsv`, matches anchors with 50kb tolerance
- Direction-aware: same anchor + different direction = separate stripes
- Filters chrX, chrY, chrM per config
- Outputs:
  - `union_stripes.tsv` — full annotations (stripe_id, source, pvals, Stripiness, etc.)
  - `union_stripes.bedpe` — 6-column with header for `stripenn score` input

**Stage 3 — Score replicates** (quantification)
- Script: `03_score_replicates.sb` + `submit_03_score.sh`
- Runs `stripenn score` on each of 6 replicates against the union BEDPE
- **`stripenn score` has NO `-h`/`--header` flag** — it auto-detects headers
- SLURM: 12 cpus, 64G, 4h
- Output: `${DATA_DIR}/outputs/{tp}/res_{kb}kb/scores/{sample}.scores.tsv`
- Key output column: `O_Sum_added` (sum of observed contacts = edgeR count input)

**Stage 4 — edgeR differential analysis**
- Scripts: `04_edgeR.R` + `04_edgeR.sb` + `submit_04_edgeR.sh`
- Builds count matrix from `O_Sum_added` (rounded to integer) across 6 replicates
- Skips `filterByExpr()` (small feature set ~150-300 stripes)
- TMM normalization, robust QL-GLM (`estimateDisp` + `glmQLFit` + `glmQLFTest`)
- Diagnostic plots: MDS, BCV, QL dispersion, volcano, MA
- Output: `${DATA_DIR}/outputs/{tp}/res_{kb}kb/04_edgeR/`

**Stage 5 — Integration & classification**
- Scripts: `05_integration.R` + `05_integration.sb` + `submit_05_integration.sh`
- Merges union annotations with edgeR stats by stripe_id
- Direction: control_only->"lost", mutant_only->"gained", shared->"strengthened"/"weakened"/"unchanged"
- Confidence tiers (high/medium/low) based on FDR + logFC agreement
- Directional consistency flag
- Output: `05_final_differential.tsv`, `.rds`, `.bedpe` (Juicebox), per-category TSVs

**Stage 6 — Cross-resolution comparison**
- Scripts: `06_compare_resolutions.R` + `06_compare_resolutions.sb` + `submit_06_compare.sh`
- Per timepoint, matches 5kb and 10kb results by anchor overlap
- High-confidence: significant at both resolutions with concordant direction
- Output: `cross_res_merged.tsv`, logFC correlation plot, direction bar chart

### 4.3 Configuration

All parameters in `config/stripenn_config.yaml`:
- Stripenn compute params: KR norm, all chroms, canny=2.0, minL=10, maxW=8, pvalue=0.1, seed=42
- Detection: 50kb anchor tolerance
- edgeR: skip filtering, TMM, robust, FDR primary=0.05, exploratory=0.10
- Classification: logFC threshold=0.3, tiered confidence
- Filtering: exclude chrX, chrY, chrM
- SLURM resource hints per stage

---

## 5. Bugs Encountered and Resolved

### 5.1 NumPy 2.0 unsigned integer overflow (CRITICAL — fixed)

**Symptom:** `stripenn compute` crashes during "3.2. Constituting background" with `ValueError: End coordinate less than start`, preceded by many `RuntimeWarning: overflow encountered in scalar subtract` warnings.

**Root cause:** `Lib.binsize` returns numpy unsigned int. Under NumPy 2.0 NEP 50, `Python_int - numpy_uint` wraps around (unsigned underflow) instead of promoting to signed. The wrapped value (~4.3 billion) bypasses the `if test_region_start0 <= 1` safety check and creates an invalid genomic region.

**Fix:** Patch `getStripe.py` in site-packages — `self.resol = int(resol)` and `chrsize = int(...)`. See Section 2 for details.

### 5.2 250402 merged .hic files missing KR normalization (fixed)

**Symptom:** hictk convert initially failed with `unable to read "GW_SCALE" weights`. After adding `--normalization-methods KR`, hictk warned `KR normalization vector is missing. SKIPPING!`. The resulting mcool files had pixel data but no KR weights.

**Root cause:** The 250402 merged .hic files were created by a different lab member using a merging tool/version that did not compute KR normalization. These files only have: GW_SCALE, INTER_SCALE, SCALE, VC, VC_SQRT. The 250831 merged files (created separately) do have KR.

**Attempted fix 1:** `juicer_tools addNorm` (via `abc/scripts/addnorm.sb`) — ran successfully (exit code 0, "Finished writing norms"), but hictk still could not read the KR vectors. Likely a Juicer v9 format incompatibility between juicer_tools 2.0.1's writer and hictk v2.1.4's reader.

**Final fix:** `cooler balance --name KR` computed ICE balancing weights directly on the mcool files and named the column "KR" so stripenn's `--norm KR` finds it. Applied via `fix_250402_balance.sb`. ICE and KR are functionally equivalent iterative matrix balancing algorithms.

**Affected files (only these two):**
- `data/cool/250402/ctrl_merged.mcool` — KR from cooler balance
- `data/cool/250402/mut_merged.mcool` — KR from cooler balance

All other 14 mcool files have native KR from the source .hic files.

### 5.3 Incomplete mcool conversions (fixed)

**Symptom:** Some Stage 1 jobs failed with `KeyError: No cooler found at ... Coolers found in ['/resolutions/5000']` — the 10kb resolution was missing.

**Root cause:** The initial hictk convert for 250402 merged files crashed on GW_SCALE after writing 5kb but before writing 10kb. The Stage 0 idempotent check (`if [ -s "${MCOOL}" ]`) saw the file existed and skipped it on re-runs.

**Fix:** Delete incomplete mcool files, reconvert with `--normalization-methods KR` to skip GW_SCALE (via `fix_250402_merged.sb`), then balance.

### 5.4 Pandas FutureWarning (cosmetic — suppressed)

**Symptom:** Verbose `FutureWarning: Series.__getitem__ treating keys as positions is deprecated` in stripenn output logs.

**Fix:** `export PYTHONWARNINGS="ignore::FutureWarning"` in Stage 1 and Stage 3 `.sb` scripts.

---

## 6. Current State (as of 2026-04-20)

### 6.1 Data files

| File | Status |
|------|--------|
| 250831 mcools (8 files) | Complete, both resolutions, native KR |
| 250402 individual replicate mcools (6 files) | Complete, both resolutions, native KR |
| 250402 ctrl_merged.mcool | Complete, both resolutions, KR from cooler balance |
| 250402 mut_merged.mcool | Complete, both resolutions, KR from cooler balance |

### 6.2 Pipeline execution

| Stage | Status | Notes |
|-------|--------|-------|
| 0 (hictk convert) | COMPLETE | All 16 mcool files verified |
| 1 (stripenn compute) | PARTIAL | 250831 5kb runs completed successfully after numpy fix. Other jobs need (re)run |
| 2 (build union) | NOT RUN | |
| 3 (score replicates) | NOT RUN | |
| 4 (edgeR) | NOT RUN | |
| 5 (integration) | NOT RUN | |
| 6 (cross-resolution) | NOT RUN | |

### 6.3 To resume

```bash
# Cancel any stale jobs from broken dependency chains
scancel -u $USER --state=PENDING

# Re-run the full pipeline from Stage 1
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn/scripts
bash run_full_stripenn.sh --skip-stage0
```

If the dependency chain breaks (e.g., one job fails), you can resubmit stages manually:
```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn/scripts

# Check which stages completed
ls ${DATA_DIR}/outputs/250831/res_5kb/

# Resubmit from the failed stage onward (no deps = run immediately)
bash submit_02_union.sh
# Then chain the rest:
JIDS=$(bash submit_02_union.sh | paste -sd:)
JIDS=$(bash submit_03_score.sh --dependency=afterok:${JIDS} | paste -sd:)
# ... etc
```

---

## 7. Verification Checklist

| Stage | Check | Expected |
|-------|-------|----------|
| 0 | `hictk metadata *.mcool` | All 16 files show resolutions [5000, 10000] |
| 0 | KR present | `cooler dump -t bins -H -c chrom,start,end,KR <uri> \| awk '$4!=""' \| head` shows numeric values |
| 1 | `wc -l result_filtered.tsv` | 100-500 stripes per merged sample |
| 2 | `wc -l union_stripes.bedpe` | >= max(ctrl, mut) count |
| 2 | source distribution | shared + control_only + mutant_only = total |
| 3 | score row count | matches union BEDPE minus header |
| 3 | `O_Sum_added` values | non-negative, mostly > 0 |
| 4 | MDS plot | ctrl/mut cluster separation |
| 4 | BCV | 0.2-0.5 range |
| 5 | source x direction | control_only all "lost", mutant_only all "gained" |
| 6 | high-confidence set | proper subset of all significant |

---

## 8. Key Corrections to PLAN.md

The original `PLAN.md` had several inaccuracies discovered during implementation:

1. **Output filename:** `result_filtered.tsv` (not `.txt`) — confirmed from `stripenn.py:156`
2. **`stripenn score` has NO `-h`/`--header` flag** — the PLAN's score command included `-h` which doesn't exist. `score.py` auto-detects headers by checking if the first row is numeric.
3. **`--seed 42`** must be passed explicitly — stripenn's default is 123456789
4. **Chromosome filtering** must be post-hoc in R — stripenn's `score.py` only filters chrY and certain scaffolds, not chrX or chrM

---

## 9. Directory Structure

### CODE root (GitHub-synced): `stripes/stripenn/`

```
stripenn/
├── CLAUDE.md                          # This document
├── PLAN.md                            # Original architectural plan
├── config/
│   └── stripenn_config.yaml           # All paths, params, thresholds
├── repo/                              # Upstream stripenn source (reference only)
│   ├── stripenn/
│   │   ├── cli.py                     # CLI entry point (compute, score, seeimage)
│   │   ├── stripenn.py                # compute() implementation
│   │   ├── score.py                   # getScore() implementation
│   │   └── getStripe.py               # Core stripe detection (patched in site-packages)
│   └── README.rst
└── scripts/
    ├── 00_convert_hic_to_cool.sb      # Stage 0: hictk convert
    ├── submit_00_convert.sh
    ├── 01_call_stripes.sb             # Stage 1: stripenn compute
    ├── submit_01_call.sh
    ├── 02_build_union.R               # Stage 2: anchor matching + union
    ├── 02_build_union.sb
    ├── submit_02_union.sh
    ├── 03_score_replicates.sb         # Stage 3: stripenn score
    ├── submit_03_score.sh
    ├── 04_edgeR.R                     # Stage 4: differential analysis
    ├── 04_edgeR.sb
    ├── submit_04_edgeR.sh
    ├── 05_integration.R               # Stage 5: classify direction
    ├── 05_integration.sb
    ├── submit_05_integration.sh
    ├── 06_compare_resolutions.R       # Stage 6: 5kb vs 10kb
    ├── 06_compare_resolutions.sb
    ├── submit_06_compare.sh
    ├── run_full_stripenn.sh           # Master orchestration driver
    ├── fix_250402_merged.sb           # One-off fix (delete after use)
    └── fix_250402_balance.sb          # One-off fix (delete after use)
```

### DATA root (NOT synced, on HPC): `/expanse/.../stripes/stripenn/`

```
stripenn/
├── data/
│   └── cool/
│       ├── 250831/                    # Early timepoint
│       │   ├── ctrl_M1.mcool         # Individual replicates (native KR)
│       │   ├── ctrl_M2.mcool
│       │   ├── ctrl_M3.mcool
│       │   ├── ctrl_merged.mcool     # Merged (native KR)
│       │   ├── mut_M1.mcool
│       │   ├── mut_M2.mcool
│       │   ├── mut_M3.mcool
│       │   └── mut_merged.mcool      # Merged (native KR)
│       └── 250402/                    # Late timepoint
│           ├── ctrl_M1.mcool         # Individual replicates (native KR)
│           ├── ...
│           ├── ctrl_merged.mcool     # Merged (KR from cooler balance)
│           └── mut_merged.mcool      # Merged (KR from cooler balance)
│
├── outputs/
│   └── {250831,250402}/              # Per timepoint
│       ├── res_5kb/                   # Per resolution
│       │   ├── calls/
│       │   │   ├── ctrl_merged/
│       │   │   │   └── result_filtered.tsv    # Stage 1 output
│       │   │   └── mut_merged/
│       │   │       └── result_filtered.tsv
│       │   ├── union_stripes.tsv      # Stage 2: full annotations
│       │   ├── union_stripes.bedpe    # Stage 2: 6-col for score input
│       │   ├── scores/                # Stage 3: per-replicate scores
│       │   │   ├── ctrl_M1.scores.tsv
│       │   │   ├── ctrl_M2.scores.tsv
│       │   │   ├── ctrl_M3.scores.tsv
│       │   │   ├── mut_M1.scores.tsv
│       │   │   ├── mut_M2.scores.tsv
│       │   │   └── mut_M3.scores.tsv
│       │   ├── 04_edgeR/              # Stage 4
│       │   │   ├── all_results.tsv
│       │   │   ├── significant_FDR05.tsv
│       │   │   ├── count_matrix.tsv
│       │   │   ├── dge_object.rds
│       │   │   ├── fit_object.rds
│       │   │   ├── summary.txt
│       │   │   └── plots/
│       │   │       ├── mds_plot.{pdf,svg,jpg}
│       │   │       ├── bcv_plot.*
│       │   │       ├── ql_dispersion_plot.*
│       │   │       ├── volcano_plot.*
│       │   │       └── ma_plot.*
│       │   ├── 05_final_differential.tsv   # Stage 5
│       │   ├── 05_final_differential.rds
│       │   ├── 05_final_differential.bedpe # Juicebox-compatible
│       │   ├── 05_stripes_lost.tsv
│       │   ├── 05_stripes_gained.tsv
│       │   ├── 05_stripes_strengthened.tsv
│       │   ├── 05_stripes_weakened.tsv
│       │   └── 05_summary.txt
│       ├── res_10kb/                  # Same structure as res_5kb/
│       ├── cross_res_merged.tsv       # Stage 6
│       ├── cross_res_summary.txt
│       └── cross_res_plots/
│           ├── logFC_correlation.*
│           └── direction_by_resolution.*
│
└── logs/                              # SLURM output logs
    ├── 00_convert_{tp}_{sample}_%j.out
    ├── 01_call_{tp}_{kb}kb_{sample}_%j.out
    ├── 02_union_{tp}_{kb}kb_%j.out
    ├── 03_score_{tp}_{kb}kb_{sample}_%j.out
    ├── 04_edgeR_{tp}_{kb}kb_%j.out
    ├── 05_integration_{tp}_{kb}kb_%j.out
    └── 06_crossres_{tp}_%j.out
```

---

## 10. Files Reference

### Pipeline scripts (CODE_DIR/scripts/)

| File | Stage | Type | Jobs |
|------|-------|------|------|
| `00_convert_hic_to_cool.sb` | 0 | SLURM | 16 |
| `submit_00_convert.sh` | 0 | Wrapper | |
| `01_call_stripes.sb` | 1 | SLURM | 8 |
| `submit_01_call.sh` | 1 | Wrapper | |
| `02_build_union.R` | 2 | R | |
| `02_build_union.sb` | 2 | SLURM | 4 |
| `submit_02_union.sh` | 2 | Wrapper | |
| `03_score_replicates.sb` | 3 | SLURM | 24 |
| `submit_03_score.sh` | 3 | Wrapper | |
| `04_edgeR.R` | 4 | R | |
| `04_edgeR.sb` | 4 | SLURM | 4 |
| `submit_04_edgeR.sh` | 4 | Wrapper | |
| `05_integration.R` | 5 | R | |
| `05_integration.sb` | 5 | SLURM | 4 |
| `submit_05_integration.sh` | 5 | Wrapper | |
| `06_compare_resolutions.R` | 6 | R | |
| `06_compare_resolutions.sb` | 6 | SLURM | 2 |
| `submit_06_compare.sh` | 6 | Wrapper | |
| `run_full_stripenn.sh` | all | Driver | |
| `fix_250402_merged.sb` | fix | One-off | DELETE after use |
| `fix_250402_balance.sb` | fix | One-off | DELETE after use |

### Config

- `config/stripenn_config.yaml` — All paths, parameters, thresholds, SLURM resources

### Architecture docs

- `PLAN.md` — Original architectural plan (has some inaccuracies; see Section 8)
- `CLAUDE.md` — This document (current state and context)

### Ported logic from Quagga pipeline

| Stripenn script | Ported from | What was reused |
|----------------|-------------|-----------------|
| `02_build_union.R` | `stripes/quagga/scripts/phase1_detection.R` | `stripes_to_granges()`, `match_anchors()`, `build_union_set()`, direction inference |
| `04_edgeR.R` | `stripes/quagga/scripts/phase3_edgeR.R` | `save_multiformat()`, full edgeR QL-GLM workflow, diagnostic plots |
| `05_integration.R` | `stripes/quagga/scripts/phase4_integration.R` | Tiered direction classification, confidence scoring, BEDPE output |
| `06_compare_resolutions.R` | `scripts/compare_resolutions.R` | Anchor matching across resolutions, concordance analysis |

---

## 10. Key Differences from Quagga Pipeline

| Aspect | Quagga Pipeline | Stripenn Pipeline |
|--------|----------------|-------------------|
| **Caller** | Quagga (Gabor filtering + Poisson) | Stripenn (Canny edge detection + pixel saturation) |
| **Quantification** | mariner `pullHicMatrices()` — 5x5 matrices at 100kb sampling point | `stripenn score` — built-in scoring, outputs O_Sum_added directly |
| **Bioconductor deps** | mariner, InteractionSet, HDF5Array, strawr | None (pure Python for detection/scoring) |
| **Input format** | .hic files directly | .mcool (converted from .hic via hictk) |
| **Replicate handling** | Call stripes on individual replicates for confidence | Score union set against replicates (no per-replicate calling) |
| **10kb validation** | Done in Phase 1 (detection) | Done in Stage 6 (cross-resolution) |
| **Confidence scoring** | Phase 1: n_reps + in_10kb + pval | Stage 5: FDR + logFC agreement with source |
| **Resolutions** | 5kb primary, 10kb validation only | Both 5kb and 10kb run independently, compared in Stage 6 |
| **SLURM pattern** | Array jobs (`--array=0-1`) for timepoints | Individual sbatch per (tp, res, sample) with wrapper scripts |
| **Parallelism** | Sequential phases, parallel within | Each job is its own SLURM submission, chained with `--dependency=afterok` |

### Why both pipelines?

The Quagga pipeline found very few significant differential calls (FDR<0.05: ~0-10 stripes). Stripenn may detect different/additional stripes due to its distinct algorithm, and its built-in `score` function provides a cleaner quantification path. Running both allows comparison of results and higher confidence in stripes called by both methods.

---

## 11. Troubleshooting

### Stage 1: stripenn compute hangs or crashes

**"End coordinate less than start":** NumPy overflow — verify patches are applied (Section 2).

**Interactive prompt hangs (no `yes Y`):** If `yes Y |` isn't piping correctly (rare), delete the output directory manually before rerunning.

**Zero stripes found:** Check that the mcool has data at the requested resolution: `cooler dump -r chr1:0-1000000 <uri> | wc -l`. If zero, the conversion failed silently.

### Stage 3: stripenn score fails with "KR not found"

The mcool doesn't have KR normalization weights. Check:
```bash
cooler dump -t bins -H -c chrom,start,end,KR <mcool_uri> | awk -F'\t' '$4 != ""' | head -3
```
If empty, run `cooler balance --name KR` on the affected file (see Section 5.2 for the 250402 merged file precedent).

### Stage 4: edgeR — no significant stripes

Expected with ~150-300 features and high detection noise. Check:
- BCV plot — if BCV > 0.5, biological variability is high
- MDS plot — if ctrl/mut don't separate, the biological effect is weak
- Use exploratory FDR < 0.10 threshold
- Focus on detection-based categories (lost/gained) rather than FDR gates

### Stage 4: edgeR — all logFC near zero

O_Sum_added values may be too similar across conditions. Check the count matrix:
```bash
head -2 ${DATA_DIR}/outputs/{tp}/res_5kb/04_edgeR/count_matrix.tsv
```
If counts are very low (<10), the stripe signal is weak at the scored regions.

### Broken dependency chain

If any job in the chain fails, all downstream jobs get `DependencyNeverSatisfied`. Fix the failed job, cancel stale jobs, and resubmit from that stage:
```bash
scancel -u $USER --state=PENDING
# Then resubmit from the failed stage onward (Section 6.3)
```

### Checking job status and logs

```bash
squeue -u $USER                           # Running/pending jobs
sacct -u $USER --starttime=today -X       # Today's completed jobs
cat ${DATA_DIR}/logs/01_call_*_5kb_*.out  # Stage 1 logs
```

---

## 12. Dependencies

### Python (in mariner_env)

- `stripenn` (pip install --no-deps, patched — see Section 2)
- `cooler` (mcool I/O)
- `numpy`, `scipy`, `pandas` (>=2.0, despite stripenn's pin)
- `opencv-python`, `scikit-image` (image processing for stripenn)
- `typer`, `joblib`, `tqdm`, `matplotlib`

### R (in mariner_env)

- `edgeR` (differential analysis — Stage 4)
- `GenomicRanges`, `IRanges` (anchor matching — Stage 2)
- `yaml` (config loading)
- `dplyr` (data manipulation)
- `ggplot2`, `svglite` (diagnostic plots)
- `VennDiagram` (cross-resolution comparison — Stage 6)

### System (HPC)

- `hictk` (separate conda env, Stage 0 only)
- SLURM job scheduler
- 96G RAM for Stage 1, 64G for Stages 0/3, 16G for Stages 2/4/5/6

---

## 13. Downstream (after pipeline completes)

The existing `stripes/quagga/scripts/stripe_visualizations.R` can be pointed at `05_final_differential.tsv` (same schema) for volcano plots, length distributions, ChIP-seq anchor annotation, and GO/KEGG enrichment — no rewrite needed.

The `cross_res_merged.tsv` from Stage 6 provides high-confidence stripe calls for publication-ready analyses.

### Planned downstream analyses

- Comparison with Quagga pipeline results (stripes called by both methods = highest confidence)
- ChIP-seq anchor annotation (H3K27ac, H3K27me3, H3K4me1 — timepoint-specific)
- GO/KEGG functional enrichment of genes near differential stripe anchors
- JuiceBox visualization using `.bedpe` files (color-coded by direction and confidence)
- Integration with loop differential results from the main mariner pipeline

---

**Last Updated:** 2026-04-20
**Project:** BAP1-KO Differential Chromatin Stripe Analysis (Stripenn track)
**Organism:** Mouse (mm10)
**HPC:** SDSC Expanse
