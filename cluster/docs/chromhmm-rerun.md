# Plan: ChromHMM 9-Mark Expansion (Intersect + Union × 15 + 18 States)

## Context

The existing 5-mark ChromHMM (H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF) produced a 12-state model for the Phase 4.4 KEY result (anchor-vs-span Polycomb enrichment). We're adding 4 marks (H2AK119ub, ATAC, H3K9ac, H3K9me3) and learning both 15-state and 18-state models with two replicate-consensus strategies (intersect and union) for the H3K9 marks. This yields **4 models** total for comparison.

**Scope:** ChromHMM segmentation + Phase 4.4/4.5 downstream only. Clustering (Phase 3, count-based) is mark-independent — reuses existing k=6. All outputs in parallel `chromHMM_9mark_{intersect,union}/` folders; existing 5-mark/12-state outputs untouched.

---

## Step 0: Create Consensus BEDs for H3K9 Replicates

Both H3K9 marks have 2 replicates. Create **both** intersect (conservative) and union (permissive) consensus.

### Intersect consensus (peaks in BOTH replicates)

```bash
# H3K9ac: Late reps
bedtools intersect -a peaks/h3k9/H3K9acCerebellumLate1.bed \
                   -b peaks/h3k9/H3K9acCerebellumLate2.bed -u \
  > peaks/h3k9/H3K9ac_consensus_intersect.bed

# H3K9me3: Early reps
bedtools intersect -a peaks/h3k9/H3K9me3CerebellumEarly1.bed \
                   -b peaks/h3k9/H3K9me3CerebellumEarly2.bed -u \
  > peaks/h3k9/H3K9me3_consensus_intersect.bed
```

Expected: H3K9ac ~12K-16K peaks; H3K9me3 ~8K-11K peaks.

### Union consensus (peaks in EITHER replicate)

```bash
# H3K9ac: concat → sort → merge overlapping
cat peaks/h3k9/H3K9acCerebellumLate1.bed peaks/h3k9/H3K9acCerebellumLate2.bed \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n \
  | bedtools merge -i - > peaks/h3k9/H3K9ac_consensus_union.bed

# H3K9me3: same approach
cat peaks/h3k9/H3K9me3CerebellumEarly1.bed peaks/h3k9/H3K9me3CerebellumEarly2.bed \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n \
  | bedtools merge -i - > peaks/h3k9/H3K9me3_consensus_union.bed
```

Expected: H3K9ac ~23K-25K peaks (Late1+Late2 union); H3K9me3 ~16K-20K peaks (Early1+Early2 union). Union uses `bedtools merge` after concatenation to collapse overlapping intervals into single regions.

---

## Step 1: New Segmentation Script

**Create:** `cluster/scripts/03b_chromhmm_9mark_segmentation.sh`

Mirrors `03_chromhmm_segmentation.sh` (same 5-stage structure, explicit `$?` checks, no `set -euo pipefail`). Parameterized by env vars so the same script handles both consensus approaches.

### 9-Mark cellmarkfiletable (same for both approaches)

```
cerebellum_late    H3K27ac      H3K27ac.bed
cerebellum_late    H3K27me3     H3K27me3.bed
cerebellum_late    H3K4me1      H3K4me1.bed
cerebellum_late    H3K4me3      H3K4me3.bed
cerebellum_late    CTCF         CTCF.bed
cerebellum_late    H2AK119ub    H2AK119ub.bed
cerebellum_late    ATAC         ATAC.bed
cerebellum_late    H3K9ac       H3K9ac.bed
cerebellum_late    H3K9me3      H3K9me3.bed
```

### Source mark array (env-var overridable for H3K9 consensus swap)

```bash
MARKS=(H3K27ac H3K27me3 H3K4me1 H3K4me3 CTCF H2AK119ub ATAC H3K9ac H3K9me3)
SOURCES=(
  "${CLUSTER_PEAK_K27AC:-peaks/beds/H3K27acCerebellumLate2.bed}"
  "${CLUSTER_PEAK_K27ME3:-peaks/beds/H3K27me3CerebellumLate1.bed}"
  "${CLUSTER_PEAK_K4ME1:-peaks/beds/H3K4me1CerebellumLate1.bed}"
  "${CLUSTER_PEAK_K4ME3:-peaks/beds/H3K4me3CerebellumLate2.bed}"
  "${CLUSTER_PEAK_CTCF:-peaks/CTCF.bed}"
  "${CLUSTER_PEAK_K119UB:-peaks/intersect/P51_K119ub_ctrl_intersect.bed}"
  "${CLUSTER_PEAK_ATAC:-peaks/atac_seq/consensus_control.bed}"
  "${CLUSTER_PEAK_K9AC}"      # NO default — must be set by caller
  "${CLUSTER_PEAK_K9ME3}"     # NO default — must be set by caller
)
```

H3K9 marks have **no default** — the sbatch wrapper must explicitly set them to the intersect or union consensus BED for each run. This prevents accidentally using the wrong consensus.

### Output directory structure (per consensus approach)

```
cluster/outputs/bap1_late/chromHMM_9mark_{intersect,union}/
    peak_beds/                          # 9 staged 3-col BEDs
    mm10_standard.txt                   # 21-chrom chromsizes
    cellmarkfiletable.txt               # 9 rows
    binarized/                          # cerebellum_late_chr*_binary.txt
    learned_model_15/                   # LearnModel k=15
        cerebellum_late_15_segments.bed
        emissions_15.txt
        model_15.txt
    learned_model_18/                   # LearnModel k=18
        cerebellum_late_18_segments.bed
        emissions_18.txt
        model_18.txt
    15state_rename_cerebellum.txt       # MANUAL
    18state_rename_cerebellum.txt       # MANUAL
```

### Script structure

1. `CHM_DIR` derived from `CLUSTER_CHROMHMM_SUBDIR` env var (e.g., `chromHMM_9mark_intersect`)
2. BinarizeBed runs once (shared binarized/ dir)
3. LearnModel loops over `NSTATES_LIST=(15 18)`, each into `learned_model_${NSTATES}/`
4. `-p ${SLURM_CPUS_PER_TASK:-4}` for auto-scaling on HPC
5. Rename file guidance printed at end uses `_cerebellum` suffix (matching existing convention)

---

## Step 2: SLURM sbatch Script

**Create:** `cluster/scripts/03b_chromhmm_9mark.sb`

Accepts a positional argument for the consensus method (`intersect`, `union`, or `both`). Submit twice for concurrent execution:

```bash
# Concurrent submission (recommended — halves wall time):
sbatch cluster/scripts/03b_chromhmm_9mark.sb intersect
sbatch cluster/scripts/03b_chromhmm_9mark.sb union

# Or sequential in one job:
sbatch cluster/scripts/03b_chromhmm_9mark.sb both
```

```
#SBATCH --job-name=chromhmm_9mark
#SBATCH --output=cluster/logs/chromhmm_9mark_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1 --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --account=csd940
```

**Script body:**

1. Parse `$1` (default: `both`). Validate ∈ {intersect, union, both}.
2. Banner (job ID, node, start time, consensus method)
3. `source ~/.bashrc && conda activate cluster`
4. Pre-flight: verify `java`, `bedtools`, all 7 non-H3K9 source BEDs exist
5. **Create consensus BEDs** for the requested method(s) only:
   - If `intersect` or `both`: create `H3K9ac_consensus_intersect.bed` + `H3K9me3_consensus_intersect.bed`
   - If `union` or `both`: create `H3K9ac_consensus_union.bed` + `H3K9me3_consensus_union.bed`
   - Verify outputs non-empty
6. **Run segmentation** for each requested method:
   ```bash
   # For intersect:
   export CLUSTER_CHROMHMM_SUBDIR=chromHMM_9mark_intersect
   export CLUSTER_PEAK_K9AC=peaks/h3k9/H3K9ac_consensus_intersect.bed
   export CLUSTER_PEAK_K9ME3=peaks/h3k9/H3K9me3_consensus_intersect.bed
   bash cluster/scripts/03b_chromhmm_9mark_segmentation.sh

   # For union:
   export CLUSTER_CHROMHMM_SUBDIR=chromHMM_9mark_union
   export CLUSTER_PEAK_K9AC=peaks/h3k9/H3K9ac_consensus_union.bed
   export CLUSTER_PEAK_K9ME3=peaks/h3k9/H3K9me3_consensus_union.bed
   bash cluster/scripts/03b_chromhmm_9mark_segmentation.sh
   ```
7. Verification: check segment BEDs + emission matrices for the method(s) run
8. Print rename file guidance

**Runtime estimate per method:** BinarizeBed ~1 min. LearnModel × 2 (15 + 18 states): 5-mark/12-state took 93s on Mac; 9-mark with 15-18 states → ~15-60 min each. Per-method total: ~1-2h. 6h wall time per job provides headroom. Concurrent submission means both methods finish in ~2h real time.

---

## Step 3: Parameterize `05_grouped_analyses.py`

**Modify:** `cluster/scripts/05_grouped_analyses.py` (lines 53-64, 88-90, 98-110, 412, 429, 487)

### 3a. Path parameterization (lines 53-64, 88-90)

Add env vars after the existing block at line 56:

```python
_nstates       = _os.environ.get('CLUSTER_NSTATES', '12')
_chromhmm_sub  = _os.environ.get('CLUSTER_CHROMHMM_SUBDIR', 'chromHMM')
_model_sub     = _os.environ.get('CLUSTER_MODEL_SUBDIR', 'learned_model')
_rename_suffix = _os.environ.get('CLUSTER_RENAME_SUFFIX', 'cerebellum')
_enrich_suffix = _os.environ.get('CLUSTER_ENRICH_SUFFIX', '')
```

Update path derivations:

```python
SEGMENT_BED  = REPO_ROOT / 'cluster/{}/{}/{}/{}_{}_segments.bed'.format(
    _out, _chromhmm_sub, _model_sub, _cell, _nstates)
RENAME_FILE  = REPO_ROOT / 'cluster/{}/{}/{}state_rename_{}.txt'.format(
    _out, _chromhmm_sub, _nstates, _rename_suffix)
...
CHROMHMM_DIR = OUT_BASE / _chromhmm_sub
```

**Why `_rename_suffix` is separate from `_cell`:** The segmentation script (`03_chromhmm_segmentation.sh` line 191) hardcodes `_cerebellum` as the rename file suffix — a human convention not derived from the `$CELL` env var (`cerebellum_late`). The existing on-disk file is `12state_rename_cerebellum.txt`, but the current code line 64 would resolve to `12state_rename_cerebellum_late.txt` — a pre-existing mismatch. Adding `_rename_suffix` defaulting to `cerebellum` fixes this cleanly: it matches the on-disk file for the existing 12-state model, and the new 9-mark segmentation script follows the same `_cerebellum` convention. The segment BED continues using `_cell` because ChromHMM generates that filename from the cellmarkfiletable cell name (`cerebellum_late`).

**Backward compatibility check:**
| Variable | Default | Existing 12-state path | 9-mark 15-state intersect path |
|----------|---------|------------------------|--------------------------------|
| `_nstates` | `12` | `12` | `15` |
| `_chromhmm_sub` | `chromHMM` | `chromHMM/` | `chromHMM_9mark_intersect/` |
| `_model_sub` | `learned_model` | `learned_model/` | `learned_model_15/` |
| `_rename_suffix` | `cerebellum` | `12state_rename_cerebellum.txt` ✓ | `15state_rename_cerebellum.txt` |
| `_enrich_suffix` | `''` | `anchor.txt`, `span.txt` | `anchor_15.txt`, `span_15.txt` |
| SEGMENT_BED | — | `learned_model/cerebellum_late_12_segments.bed` ✓ | `learned_model_15/cerebellum_late_15_segments.bed` |

All defaults resolve to existing paths. No breakage.

### 3b. OverlapEnrichment output naming (line 412, 429)

Suffix the output prefix and txt check:

```python
output_prefix = CHROMHMM_DIR / f'{coord_kind}{_enrich_suffix}'
# ...
out_txt = CHROMHMM_DIR / f'{coord_kind}{_enrich_suffix}.txt'
```

Without this, running 15-state then 18-state would overwrite `anchor.txt`/`span.txt`. With suffix: `anchor_15.txt`, `anchor_18.txt`, etc.

### 3c. Phase 4.5 figure subfolder (line 487)

```python
sub = figure_subfolder(FIG_BASE, f'chromHMM_anchor{_enrich_suffix}')
```

Produces `figures/chromHMM_anchor_15/` for 15-state, `figures/chromHMM_anchor_18/` for 18-state.

### 3d. STATE_COLORS expansion (lines 98-110)

Add provisional entries for new state types. These are placeholders — actual names will be assigned in the manual rename step after inspecting emissions:

```python
# States anticipated from 9-mark model (H3K9me3, H2AK119ub, ATAC additions)
'Heterochromatin':       '#1a1a2e',   # dark navy — H3K9me3-only
'Polycomb_K119ub':       '#4a0e4e',   # dark purple — K27me3 + K119ub
'Polycomb_K9me3':        '#8c564b',   # brown — K27me3 + K9me3
'Open_Chromatin':        '#00ced1',   # dark turquoise — ATAC-only
'ATAC_Enhancer':         '#a8e6cf',   # light green — ATAC + K27ac + K4me1
'Active_K9ac':           '#e377c2',   # pink — K9ac + K27ac
```

States not in the dict render as black (`#000000`) via `.get()` fallback — functional but will need updating after the actual state names are decided.

---

## Step 4: Downstream Driver Script

**Create:** `cluster/scripts/run_phase4_9mark.sh`

```bash
#!/usr/bin/env bash
# Usage:
#   bash cluster/scripts/run_phase4_9mark.sh                 # all 4 models
#   bash cluster/scripts/run_phase4_9mark.sh intersect 15    # specific model
#   bash cluster/scripts/run_phase4_9mark.sh union both      # union × 15+18
```

For each (consensus, nstates) combination, sets:

```bash
export CLUSTER_NSTATES=${NSTATES}
export CLUSTER_CHROMHMM_SUBDIR=chromHMM_9mark_${CONSENSUS}
export CLUSTER_MODEL_SUBDIR=learned_model_${NSTATES}
export CLUSTER_ENRICH_SUFFIX=_${NSTATES}
export CLUSTER_RENAME_SUFFIX=cerebellum
```

Then calls `$PYTHON scripts/05_grouped_analyses.py --analyses 4.4,4.5`.

Loop over all requested combinations. Outputs land in their respective dirs:
- `chromHMM_9mark_intersect/anchor_15.txt`, `span_15.txt`, `anchor_15.{png,...}`
- `chromHMM_9mark_intersect/anchor_18.txt`, `span_18.txt`, `anchor_18.{png,...}`
- `chromHMM_9mark_union/anchor_15.txt`, `span_15.txt`, etc.
- `figures/chromHMM_anchor_15/` per chromHMM subdir

Can run on Mac (fast — OverlapEnrichment takes seconds) or HPC.

---

## Step 5: Manual Rename Step (After LearnModel)

Inspect all 4 emission matrices:
```
chromHMM_9mark_intersect/learned_model_15/emissions_15.txt
chromHMM_9mark_intersect/learned_model_18/emissions_18.txt
chromHMM_9mark_union/learned_model_15/emissions_15.txt
chromHMM_9mark_union/learned_model_18/emissions_18.txt
```

Each is `NSTATES × 9` (cols: H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF, H2AK119ub, ATAC, H3K9ac, H3K9me3).

**New biology the 9-mark model can resolve:**
- `Heterochromatin` — H3K9me3 high, everything else low (constitutive)
- `Polycomb_K119ub` — H3K27me3 + H2AK119ub high (PRC1+PRC2 co-marked, key for BAP1-KO)
- `Polycomb_K9me3` — H3K27me3 + H3K9me3 (facultative → constitutive transition)
- `Open_Chromatin` — ATAC high, no histone marks (accessible but unmarked)
- `ATAC_Enhancer` — ATAC + H3K27ac + H3K4me1 (active accessible enhancer)
- `Active_K9ac` — H3K9ac + K27ac + K4me3 (active promoter with K9ac)

**Model comparison criteria** (across all 4 models):
1. Log-likelihood from `model_{N}.txt` (compare via BIC = -2*LL + k*ln(n))
2. No near-duplicate emission profiles (all pairs differ by ≥1 mark on/off)
3. Quiescent fragmentation check (>2 quiescent-like states = too many states)
4. Intersect vs union: does the union approach resolve additional states, or just add noise?

Create 4 rename files:
```
chromHMM_9mark_intersect/15state_rename_cerebellum.txt
chromHMM_9mark_intersect/18state_rename_cerebellum.txt
chromHMM_9mark_union/15state_rename_cerebellum.txt
chromHMM_9mark_union/18state_rename_cerebellum.txt
```

Format: `E{N}\t{Biological_Name}` (underscored names for cross-script compatibility).

---

## File Summary

### Files to create

| File | Purpose |
|------|---------|
| `peaks/h3k9/H3K9ac_consensus_intersect.bed` | Intersect of 2 late reps (created by sbatch) |
| `peaks/h3k9/H3K9me3_consensus_intersect.bed` | Intersect of 2 early reps (created by sbatch) |
| `peaks/h3k9/H3K9ac_consensus_union.bed` | Union merge of 2 late reps (created by sbatch) |
| `peaks/h3k9/H3K9me3_consensus_union.bed` | Union merge of 2 early reps (created by sbatch) |
| `cluster/scripts/03b_chromhmm_9mark_segmentation.sh` | 9-mark segmentation worker |
| `cluster/scripts/03b_chromhmm_9mark.sb` | SLURM wrapper (runs both intersect + union) |
| `cluster/scripts/run_phase4_9mark.sh` | Driver for 4.4/4.5 rerun |

### Files to modify

| File | Changes |
|------|---------|
| `cluster/scripts/05_grouped_analyses.py` | Add 5 env vars (`CLUSTER_NSTATES`, `CLUSTER_CHROMHMM_SUBDIR`, `CLUSTER_MODEL_SUBDIR`, `CLUSTER_RENAME_SUFFIX`, `CLUSTER_ENRICH_SUFFIX`). Update `SEGMENT_BED`, `RENAME_FILE`, `CHROMHMM_DIR`, enrichment output prefix/txt, 4.5 figure subfolder. Add `STATE_COLORS` entries. All backward-compatible. |

### Files NOT modified

- `cluster/scripts/03_chromhmm_segmentation.sh`
- `cluster/outputs/bap1_late/chromHMM/` (entire existing directory)
- All Phase 3 clustering outputs

---

## Verification Checklist

### After consensus BEDs
- [ ] 4 consensus BEDs exist and are non-empty
- [ ] Intersect counts < union counts for both marks

### After segmentation (all 4 models)
- [ ] 2 × `peak_beds/` dirs each with 9 `.bed` files, 3-col, standard chroms
- [ ] 2 × `binarized/` dirs with 21 binary files; headers have 10 cols (cell + 9 marks)
- [ ] 4 emission files: each `NSTATES × 9` + header
- [ ] 4 segment BEDs covering all 21 chroms
- [ ] Existing `chromHMM/` directory mtimes unchanged

### After Phase 4.4/4.5 rerun (all 4 models)
- [ ] 4 × {anchor, span} enrichment tables + heatmaps
- [ ] 4 × proportions TSVs
- [ ] Polycomb/K119ub states show clust5/clust6 enrichment (biological sanity)
- [ ] Original `chromHMM/anchor.txt` untouched

---

## Caveats

### Timepoint limitation (document in figure legends)

| Mark | Timepoint | Notes |
|------|-----------|-------|
| H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF, H3K9ac | Late (P60) | Matches Hi-C |
| H2AK119ub | Late (P51) | 9 days from P60, acceptable |
| ATAC | Combined (ctrl consensus) | Closest available |
| **H3K9me3** | **Early** | No late peaks exist — limitation |

Only H3K9me3 is a true timepoint mismatch. H2AK119ub at P51 is close enough to P60 to be considered late.

### Ctrl-only peaks
All 9 marks use wildtype/control peak calls. ChromHMM represents the control chromatin landscape.

### Intersect vs union tradeoff
Intersect (peaks in both reps) is more conservative and higher confidence. Union (peaks in either rep) captures weaker/replicate-specific signals and will produce more peaks. Comparing both approaches shows whether the extra union peaks meaningfully change state assignments at loop anchors, or just add noise to sparse states.
