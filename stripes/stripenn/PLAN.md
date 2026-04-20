# Stripenn Pipeline Setup — Differential Architectural Stripe Analysis

## Context

The existing Quagga-based pipeline (`stripes/README.md`) found ~150–300 stripes per timepoint but very few significant differential calls (FDR<0.05: ~0–10). PI suggested trying **Stripenn** (Yoon & Vahedi, upstream at `stripes/stripenn/`) — an image-processing-based stripe caller already cloned into the repo.

Goal: build a parallel pipeline under `stripes/stripenn/` that (1) converts `.hic` to `.cool`, (2) calls stripes per condition-merged file with Stripenn, (3) scores the union set against each replicate to get a replicate-aware count-like matrix, and (4) feeds that into edgeR for differential testing — mirroring the structure of the existing Quagga / mariner pipelines so downstream annotation/visualization scripts can be reused.

Inputs already on HPC at `/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/{250831,250402}/` — 8 `.hic` per timepoint (`{ctrl,mut}_{M1,M2,M3,merged}.hic`).

Stripenn only accepts `.cool`/`.mcool`, so an upstream conversion step is required.

---

## User Decisions

- **Env:** install into `mariner_env`.
- **Resolutions:** run at **both 5 kb and 10 kb** (cross-resolution comparison like the mariner loop pipeline).
- **Diff strategy:** full replicate-aware edgeR (score union stripes against 6 replicates).
- **Scope:** both timepoints in a single master SLURM job.
- **Parallelism:** one `sbatch` submission per sample/stripe-call/score job — a wrapper `.sh` iterates and submits concurrent jobs. Each `.sb` accepts CLI args instead of encoding the sample list inside the script. Follows existing `phase1_detection.sb`/`phase2_quantification.sb` conventions (absolute paths, `cd /expanse/.../mariner_hi-c/stripes`, `conda activate mariner_env`, `mkdir -p logs`, exit-code checks), but replaces `--array=0-1` with argparse-style positional args + a wrapper.
- **Path split (important):**
  - **CODE root** (synced to GitHub): `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn/` — scripts, config, small final summary TSVs.
  - **DATA root** (not synced, >100 MB): `/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/` — `data/cool/` mcools, `outputs/` raw stripenn compute/score results, intermediate RDS, plots.
  - Every path in every script is absolute; no relative paths. Final small TSVs / final volcano PDFs can be manually rsynced to the repo later.

Total workload: 2 tp × 2 res × (2 merged calls + 6 replicate scores) = **8 `stripenn compute` jobs + 24 `stripenn score` jobs** on top of 16 `hic2cool` jobs. With SLURM concurrency each job is its own submission and runs in parallel on its own node/partition slot.

## Architecture

**CODE root** — `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn/`:

```
stripenn/                                     (in GitHub repo)
├── config/
│   └── stripenn_config.yaml                  (paths, resolutions, thresholds, sample lists)
├── scripts/
│   ├── 00_convert_hic_to_cool.sb             (CLI: <timepoint> <sample>)
│   ├── submit_00_convert.sh                  (wrapper: sbatch × 16)
│   ├── 01_call_stripes.sb                    (CLI: <timepoint> <resolution> <sample>)
│   ├── submit_01_call.sh                     (wrapper: sbatch × 8)
│   ├── 02_build_union.R                      (CLI args: <timepoint> <resolution>)
│   ├── 02_build_union.sb                     (CLI: <timepoint> <resolution>)
│   ├── submit_02_union.sh                    (wrapper: sbatch × 4)
│   ├── 03_score_replicates.sb                (CLI: <timepoint> <resolution> <sample>)
│   ├── submit_03_score.sh                    (wrapper: sbatch × 24)
│   ├── 04_edgeR.R + 04_edgeR.sb              (CLI: <timepoint> <resolution>)
│   ├── submit_04_edgeR.sh                    (wrapper: sbatch × 4)
│   ├── 05_integration.R + 05_integration.sb  (CLI: <timepoint> <resolution>)
│   ├── submit_05_integration.sh              (wrapper: sbatch × 4)
│   ├── 06_compare_resolutions.R + .sb        (CLI: <timepoint>)
│   └── submit_06_compare.sh                  (wrapper: sbatch × 2)
└── (small final TSVs rsync'd back here manually as needed)
```

**DATA root** — `/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/` (NOT in repo, already created by user):

```
stripenn/
├── data/
│   └── cool/
│       ├── 250831/{ctrl,mut}_{M1..3,merged}.mcool
│       └── 250402/{ctrl,mut}_{M1..3,merged}.mcool
├── outputs/
│   └── {timepoint}/                          (250831|250402)
│       └── res_{5,10}kb/
│           ├── calls/{ctrl,mut}_merged/      (stripenn compute result_*.txt)
│           ├── union_stripes.tsv + .bedpe
│           ├── scores/{sample}.scores.out    (6 replicates)
│           ├── 04_edgeR/                     (RDS, TSV, volcano/MDS/BCV PDFs)
│           └── 05_final_differential.tsv + .bedpe
│   └── {timepoint}/cross_res_merged.tsv      (step 6)
└── logs/                                     (SLURM stdout/err per submission)
```

Both paths are set as shell variables at the top of every `.sb`:
```bash
CODE_DIR=/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn
DATA_DIR=/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn
```

---

## Step 1 — Environment & Installation

**Recommended:** install into `mariner_env` (already has R/edgeR for downstream; just needs Python deps).

```bash
conda activate mariner_env
# Stripenn + deps (README explicitly says install via PyPI, NOT git clone)
pip install stripenn hic2cool
# hic2cool pulls cooler, h5py — stripenn pulls opencv-python, scikit-image, scipy
```

Verify:
```bash
stripenn --help
hic2cool --help
python -c "import cooler, cv2, skimage; print('ok')"
```

The cloned `stripes/stripenn/` source is kept only as reference — **do not** `pip install -e ./stripes/stripenn/` per the README warning.

---

## Step 2 — `.hic` → `.mcool` Conversion

Uses `hic2cool convert` to produce multi-resolution `.mcool` preserving KR weights from the original `.hic`.

`00_convert_hic_to_cool.sb` accepts two positional args — timepoint and sample:
```bash
# inside the .sb
TIMEPOINT=$1         # 250831 | 250402
SAMPLE=$2            # ctrl_M1 | ... | mut_merged
HIC=/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/${TIMEPOINT}/${SAMPLE}.hic
MCOOL=${DATA_DIR}/data/cool/${TIMEPOINT}/${SAMPLE}.mcool
mkdir -p "$(dirname ${MCOOL})"
hic2cool convert "${HIC}" "${MCOOL}" -p ${SLURM_CPUS_PER_TASK}
```

Wrapper `submit_00_convert.sh` (in `scripts/`):
```bash
#!/bin/bash
CODE_DIR=/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn
for TP in 250831 250402; do
  for SAMPLE in ctrl_M1 ctrl_M2 ctrl_M3 ctrl_merged mut_M1 mut_M2 mut_M3 mut_merged; do
    sbatch --job-name=hic2cool_${TP}_${SAMPLE} \
           --output=${CODE_DIR}/../../../stripes/stripenn/logs/hic2cool_${TP}_${SAMPLE}_%j.out \
           ${CODE_DIR}/scripts/00_convert_hic_to_cool.sb ${TP} ${SAMPLE}
  done
done
```

Submits 16 concurrent jobs (SLURM schedules as nodes free). Expected wall time per job: ~10–30 min. Output: `.mcool` with all native resolutions + KR weights.

Target resolutions for stripenn: **5 kb and 10 kb** (both embedded in the mcool; stripenn selects via `::/resolutions/{res}`).

---

## Step 3 — Call Stripes on Merged Files

`01_call_stripes.sb` takes `<timepoint> <resolution_bp> <sample>` and runs one `stripenn compute`:
```bash
TIMEPOINT=$1   # 250831 | 250402
RES=$2         # 5000 | 10000
SAMPLE=$3      # ctrl_merged | mut_merged
RES_KB=$((RES/1000))
MCOOL=${DATA_DIR}/data/cool/${TIMEPOINT}/${SAMPLE}.mcool
OUT=${DATA_DIR}/outputs/${TIMEPOINT}/res_${RES_KB}kb/calls/${SAMPLE}/
mkdir -p "${OUT}"
yes Y | stripenn compute \
  --cool ${MCOOL}::/resolutions/${RES} \
  --out ${OUT} \
  -k all \
  -m 0.95,0.96,0.97,0.98,0.99 \
  -p 0.1 \
  -n ${SLURM_CPUS_PER_TASK} \
  --norm KR
```

(`yes Y` handles stripenn's interactive overwrite prompt when rerunning.)

Wrapper `submit_01_call.sh` submits **8 jobs**:
```bash
for TP in 250831 250402; do
  for RES in 5000 10000; do
    for S in ctrl_merged mut_merged; do
      sbatch 01_call_stripes.sb ${TP} ${RES} ${S}
    done
  done
done
```

Chromosome filtering (`chrX,chrY,chrM` exclusion from existing config) can be done post-hoc or via `-k` listing autosomes.

Output: `outputs/{tp}/calls/{ctrl,mut}_merged/result_filtered.txt` — TSV with 12 columns (`chr, pos1..pos4, length, width, Mean, maxpixel, pvalue, Stripiness`).

---

## Step 4 — Build Union (`02_build_union.R`)

Mirror `phase1_detection.R` logic:
- Load both `result_filtered.txt` per timepoint.
- Compute anchor match with 50 kb tolerance (config: `anchor_tolerance_bp: 50000`).
- Classify `control_only`, `mutant_only`, `shared`.
- Output union as 6-column BEDPE (stripenn `score` expects that format: `chr,pos1,pos2,chr2,pos3,pos4`).

Output: `outputs/{tp}/union_stripes.tsv` + `union_stripes.bedpe`.

---

## Step 5 — Score Replicates

`03_score_replicates.sb` takes `<timepoint> <resolution_bp> <sample>`:
```bash
TIMEPOINT=$1; RES=$2; SAMPLE=$3
RES_KB=$((RES/1000))
MCOOL=${DATA_DIR}/data/cool/${TIMEPOINT}/${SAMPLE}.mcool
UNION=${DATA_DIR}/outputs/${TIMEPOINT}/res_${RES_KB}kb/union_stripes.bedpe
OUT=${DATA_DIR}/outputs/${TIMEPOINT}/res_${RES_KB}kb/scores/${SAMPLE}.scores.out
mkdir -p "$(dirname ${OUT})"
stripenn score \
  --cool ${MCOOL}::/resolutions/${RES} \
  -c ${UNION} \
  -h \
  -n ${SLURM_CPUS_PER_TASK} \
  --norm KR \
  -o ${OUT}
```

Wrapper `submit_03_score.sh` submits **24 jobs** (2 tp × 2 res × 6 replicates):
```bash
for TP in 250831 250402; do
  for RES in 5000 10000; do
    for S in ctrl_M1 ctrl_M2 ctrl_M3 mut_M1 mut_M2 mut_M3; do
      sbatch 03_score_replicates.sb ${TP} ${RES} ${S}
    done
  done
done
```

Output columns include `O_Sum_added`, `O_Mean_added`, `O/E_Mean_added`, `O/E_Total_added` — `O_Sum_added` is the natural input for edgeR (sum of raw-like contacts per stripe).

---

## Step 6 — edgeR Differential (`04_edgeR.R`)

Build count matrix: rows = union stripes, cols = 6 replicates. Use `O_Sum_added` as counts (round to integer). Mirror `scripts/phase3_edgeR.R` statistical framework:

```r
dge <- DGEList(counts = count_matrix, group = c("ctrl","ctrl","ctrl","mut","mut","mut"))
# Skip filterByExpr — feature count is small (~200-400)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~group)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)
results <- topTags(qlf, n = Inf)$table
```

Output: `04_edgeR/all_results.tsv`, MDS/BCV/volcano PDFs, `dge.rds`.

---

## Step 7 — Integration (`05_integration.R`)

Join Step 4 detection annotations (`source = control_only|mutant_only|shared`) with Step 6 stats. Classify using the existing tiered scheme (see `stripes/README.md` lines 212–225):

- `lost`: `source == control_only` (tiered on FDR/logFC)
- `gained`: `source == mutant_only` (tiered)
- `strengthened`: `shared` + FDR<0.05 + logFC>0.3
- `weakened`: `shared` + FDR<0.05 + logFC<-0.3
- `unchanged`: everything else

Output: `05_final_differential.tsv`, `05_final.bedpe` (Juicebox), per-category TSVs.

**Downstream**: the existing `scripts/stripe_visualizations.R` can be pointed at this new TSV (same schema) for volcano plots, length distribution, ChIP-seq anchor annotation, GO/KEGG enrichment — no rewrite needed.

---

## Step 8 — Cross-Resolution Comparison (`06_compare_resolutions.R`)

Per timepoint, merge 5 kb and 10 kb final differential TSVs (following `scripts/compare_resolutions.R` from the mariner loop pipeline):
- Overlap analysis (Venn / Jaccard) of stripes called at both resolutions.
- Stripes significant at both resolutions → **high-confidence** set.
- Resolution-specific calls → **exploratory** set.

Output: `outputs/{tp}/cross_res_merged.tsv` + overlap plots.

---

## Critical Files Referenced / Reused

- `stripes/README.md` — pipeline conventions, tiered classification logic
- `stripes/config/stripe_config.yaml` — sample names, paths, thresholds to copy/adapt
- `stripes/scripts/phase1_detection.R` — anchor matching with 50kb tolerance (port logic to `02_build_union.R`)
- `stripes/scripts/phase3_edgeR.R` — edgeR template (Skip-filter + TMM + robust QL-GLM)
- `stripes/scripts/phase4_integration.R` — classification rules
- `stripes/scripts/stripe_visualizations.R` — reusable as-is for Step 7 output
- `loop-prep/loopbin_prep.sb` (user-supplied) — SLURM + parallel `&`/`wait` pattern
- `stripes/stripenn/src/stripenn/cli.py` lines 9–29 — compute args reference
- `stripes/stripenn/src/stripenn/score.py` — score output columns

---

## Verification Plan

1. **After install**: `stripenn compute --help` and `hic2cool --help` both print usage.
2. **After conversion** (Step 2): `cooler ls data/cool/250831/ctrl_merged.mcool` shows expected resolutions; `cooler info ::/resolutions/10000` shows nnz > 0 and KR weight column present.
3. **Smoke test** (Step 3): run `stripenn compute` on a single chromosome (`-k chr19`) on `ctrl_merged` 250831 — expect completion in <15 min and a populated `result_filtered.txt`.
4. **After Step 3 full run**: both merged calls produce 100–500 stripes each (Stripenn is more permissive than Quagga; count may be higher).
5. **After Step 5**: every replicate score file has the same number of rows as union stripes; `O_Sum_added` > 0 for the vast majority.
6. **After Step 6**: MDS plot separates ctrl vs mut clusters; BCV ~0.2–0.5; volcano shows some differential features.
7. **End-to-end**: run `stripe_visualizations.R` pointed at `05_final_differential.tsv` and confirm all downstream plots render.

---

## Resolved User Decisions

- Env: `mariner_env` (shared with R downstream).
- Resolutions: 5 kb AND 10 kb (cross-res comparison in Step 8).
- Differential strategy: full replicate-aware edgeR.
- Scope: both timepoints in one master SLURM job.

## Orchestration — Staged Submission with `--dependency=afterok`

No single monolithic `.sb`. Instead, each stage's wrapper `.sh` submits its jobs with a dependency on the previous stage. Top-level driver `run_full_stripenn.sh`:

```bash
#!/bin/bash
CODE_DIR=/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn
cd ${CODE_DIR}/scripts

# Stage 0: hic2cool × 16
JIDS_00=$(bash submit_00_convert.sh | awk '{print $NF}' | paste -sd:)

# Stage 1: stripenn compute × 8 (depends on all stage 0 jobs)
JIDS_01=$(bash submit_01_call.sh --dependency=afterok:${JIDS_00} | awk '{print $NF}' | paste -sd:)

# Stage 2: build_union × 4
JIDS_02=$(bash submit_02_union.sh --dependency=afterok:${JIDS_01} | awk '{print $NF}' | paste -sd:)

# Stage 3: score × 24 (needs union_stripes.bedpe)
JIDS_03=$(bash submit_03_score.sh --dependency=afterok:${JIDS_02} | awk '{print $NF}' | paste -sd:)

# Stage 4: edgeR × 4
JIDS_04=$(bash submit_04_edgeR.sh --dependency=afterok:${JIDS_03} | awk '{print $NF}' | paste -sd:)

# Stage 5: integration × 4
JIDS_05=$(bash submit_05_integration.sh --dependency=afterok:${JIDS_04} | awk '{print $NF}' | paste -sd:)

# Stage 6: cross-resolution × 2
bash submit_06_compare.sh --dependency=afterok:${JIDS_05}
```

Each `submit_*.sh` accepts an optional trailing `--dependency=…` arg that it splices into its `sbatch` calls. Within a stage, all jobs run in parallel; across stages, SLURM enforces ordering. Expected wall time: **8–16 h** (vs ~30 h if fully serial) — stage 0 is the longest bottleneck (~30 min/file × 16 files on queue).

Each stage's `.sb` uses the existing repo conventions: absolute paths, `source ~/.bashrc; conda activate mariner_env`, `logs/` under `${DATA_DIR}/logs/`, exit-code validation (matching `phase1_detection.sb:40-67`).
