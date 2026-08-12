# Phase 2 Plan: ChromHMM 12-state Segmentation for BAP1-KO Cerebellum

**Companion docs:** `cluster/CONTEXT-CLUSTER.md` (biology + meeting context), `cluster/PLAN-CLUSTER.md` (overall pipeline). This plan refines and operationalizes Phase 2 of `PLAN-CLUSTER.md` (lines 205–262).

---

## Context

Phase 1 is complete: 39,344 nonredundant merged loops are in Popay format (`cluster/data/late_merged_loop_counts.txt`) with a metadata sidecar carrying edgeR stats, plus a 24,515-row mm10 promoter BED. The mariner Popay-collaboration pipeline now needs a chromatin-state segmentation of the mouse cerebellum to drive the headline analysis: **Phase 4.4 ChromHMM anchor-vs-span enrichment (Popay Fig 2f equivalent).**

That figure is the mechanistic answer the Dixon meeting asked for — it tells us whether Polycomb in BAP1-KO sits *at* loop anchors (sensitivity model: H3K27me3 directly disrupts CTCF binding) or *across the loop body* (extrusion-impediment model: K27me3 spreading blocks cohesin extrusion of long loops). All of Phase 4.4–4.5 read directly from the segmentation BED produced here.

Phase 2 takes the 5 cerebellum ChIP peak BEDs (H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF) → ChromHMM BinarizeBed → LearnModel(k=12) → `cerebellum_late_12_segments.bed` + `emissions_12.txt`. The final step is human inspection of the emission matrix to write a state-id-to-biological-name rename file that matches the project's existing 7-category chromatin taxonomy.

---

## Pre-flight verification (already done)

| Check | Result |
|---|---|
| 5 peak BEDs exist with expected line counts | ✅ K27ac 15,105 / K27me3 15,809 / K4me1 113,781 / K4me3 6,581 / CTCF 32,487 |
| All BEDs have ≥3 cols (BinarizeBed only reads cols 1-3) | ✅ K27ac/K27me3/K4me1/K4me3 = 6 cols, CTCF = 10 cols |
| Chromosome content | ⚠️ K27ac, K27me3, K4me1, CTCF contain `chr*_random`, `chrUn_*`, `chrM` records — must be filtered to chr1-19 + chrX + chrY before BinarizeBed (K4me3 is already clean) |
| `cluster/ChromHMM/ChromHMM.jar` present (3.8 MB) | ✅ |
| Wrapper script `cluster/ChromHMM/chromhmm` works (-mx8G) | ✅ |
| Java 25 OpenJDK installed | ✅ (`-mx` deprecation warning is benign — leave as-is) |
| `cluster/ChromHMM/CHROMSIZES/mm10.txt` covers all 21 standard chroms | ✅ (also lists chrM + many random/chrUn we will exclude) |
| `cluster/bap1_late/chromHMM/{peak_beds,binarized,learned_model,...}` exist (created in Phase 1) and are empty | ✅ no partial outputs to clean |
| No prior `03_chromhmm_segmentation.sh` exists | ✅ writing fresh |
| Popay rename example at `cluster/clustering_example_data/12state_rename.txt` | ✅ format confirmed: `<state_id><TAB><name>`, 12 lines, tab-separated, names use spaces (we will use underscores for cross-script compatibility — see below) |

---

## Key scientific decisions

1. **Standard-chromosome filter (chr1-19 + chrX + chrY) on both peak BEDs and chromsizes.** All four ChIP files except K4me3 contain peaks on `chr*_random` and `chrUn_*` contigs, plus K27me3 has 1 chrM peak. None of our 39,344 loops or 24,515 promoters live on these contigs, so segmenting them produces ChromHMM bins that no downstream analysis ever touches. We drop them at staging time.

2. **No mm10-blacklist filtering of peak BEDs.** Popay's pipeline doesn't pre-filter peaks against blacklist; ChromHMM treats blacklist regions correctly as "no signal" when the underlying peak callers already excluded them. Blacklist will be applied later in Phase 5 (deepTools `--blackListFileName`) where signal averaging actually needs it.

3. **k=12 states.** Direct comparability with Popay Fig 2f. Caveat: with only 5 marks (vs Popay's 12 marks), expect 2–4 redundant low-signal/quiescent states. Acceptable — Popay's own rename file already has a "Low signal" state. If post-hoc the model is degenerate (>50% bins in one state, or several states with near-identical emissions), we revisit with k=10 or use `chromhmm CompareModels` — not in scope for Phase 2.

4. **Cell name = `cerebellum_late`** — becomes the prefix for ChromHMM output files (`cerebellum_late_12_segments.bed`, `cerebellum_late_chr1_binary.txt`, ...). Distinguishes from a possible future P12/early run.

5. **State-name convention: underscores, not spaces** (e.g. `Active_Promoter`, not `Active promoter`). The project's existing annotation taxonomy (`annotate_loops_extended.R` per `CLAUDE.md`) uses underscored names: `Active_Promoter`, `Repressed_Promoter`, `Bivalent_Promoter`, `Polycomb`, `Active_Enhancer`, `Poised_Enhancer`, `Other`. Matching this lets Phase 4 plotting code do clean string comparisons against existing palettes/dicts.

6. **`set -e` not used in the runner.** Per project memory and matching `run_phase1.sh`: explicit `$?` checking after each long step rather than `set -euo pipefail` (which would abort on benign warnings from Java/awk).

---

## Files to create

| Path | Purpose |
|---|---|
| `cluster/scripts/03_chromhmm_segmentation.sh` | Worker: stage BEDs, filter chromsizes, write cellmarkfiletable, BinarizeBed, LearnModel, verify |
| `cluster/scripts/run_phase2.sh` | Driver: cd to repo root, call the worker, print final summary (mirrors `run_phase1.sh` style) |
| `cluster/bap1_late/chromHMM/peak_beds/{H3K27ac,H3K27me3,H3K4me1,H3K4me3,CTCF}.bed` | Staged 3-col BEDs filtered to standard chroms (worker output) |
| `cluster/bap1_late/chromHMM/mm10_standard.txt` | 21-row chromsizes (chr1-19 + chrX + chrY only) — used for both BinarizeBed and LearnModel `-l` (worker output) |
| `cluster/bap1_late/chromHMM/cellmarkfiletable.txt` | 5-row tab table mapping marks to BED filenames (worker output) |
| `cluster/bap1_late/chromHMM/binarized/cerebellum_late_chr*_binary.txt` | 21 binary files from BinarizeBed (one per chrom) |
| `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` | Final segmentation (~0.5–5M lines) |
| `cluster/bap1_late/chromHMM/learned_model/emissions_12.txt` | 12×5 emission probability matrix |
| `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` | **Manually authored after LearnModel completes** — maps state IDs to biological names |
| `cluster/bap1_late/chromHMM/phase2.txt` | Run log (created via `tee`/redirect when invoking `run_phase2.sh`) |

No edits to existing files — Phase 0 already fixed the Popay modules; Phase 1 produced the data.

---

## Implementation steps

### Step 1 — `03_chromhmm_segmentation.sh` (worker script)

Mirrors the `run_phase1.sh` style: explicit status checks, no `set -e`, banner-style logging, paths derived from `BASH_SOURCE` to find repo root.

Outline:

```bash
#!/usr/bin/env bash
# cluster/scripts/03_chromhmm_segmentation.sh
# Phase 2 worker: ChromHMM 12-state segmentation from 5 cerebellum ChIP peak BEDs.
# Input  : peaks/{beds/H3K27acCerebellumLate2,beds/H3K27me3CerebellumLate1,
#          beds/H3K4me1CerebellumLate1,beds/H3K4me3CerebellumLate2,CTCF}.bed
# Output : cluster/bap1_late/chromHMM/learned_model/{cerebellum_late_12_segments.bed,
#          emissions_12.txt,transitions_12.txt,model_12.txt,webpage_12.html}
# Runtime: ~30-90 min (LearnModel is the bottleneck; BinarizeBed is ~1-2 min)

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CHROMHMM="${REPO_ROOT}/cluster/ChromHMM/chromhmm"
CHM_DIR="${REPO_ROOT}/cluster/bap1_late/chromHMM"
PEAK_DIR="${CHM_DIR}/peak_beds"
BIN_DIR="${CHM_DIR}/binarized"
MODEL_DIR="${CHM_DIR}/learned_model"
CMT="${CHM_DIR}/cellmarkfiletable.txt"
CHROMSIZES_SRC="${REPO_ROOT}/cluster/ChromHMM/CHROMSIZES/mm10.txt"
CHROMSIZES_DST="${CHM_DIR}/mm10_standard.txt"
CELL="cerebellum_late"
NSTATES=12
ASSEMBLY="mm10"

cd "${REPO_ROOT}"
```

The worker performs five sub-steps:

**[1/5] Stage peak BEDs** — for each mark, awk-filter to `chr1-19|chrX|chrY`, take cols 1-3, sort and dedupe, write to `${PEAK_DIR}/<MARK>.bed`. Echo input→output line counts so any unexpectedly aggressive filter is visible.

```bash
declare -A SOURCES=(
  [H3K27ac]="peaks/beds/H3K27acCerebellumLate2.bed"
  [H3K27me3]="peaks/beds/H3K27me3CerebellumLate1.bed"
  [H3K4me1]="peaks/beds/H3K4me1CerebellumLate1.bed"
  [H3K4me3]="peaks/beds/H3K4me3CerebellumLate2.bed"
  [CTCF]="peaks/CTCF.bed"
)
for MARK in H3K27ac H3K27me3 H3K4me1 H3K4me3 CTCF; do
  SRC="${SOURCES[$MARK]}"
  DST="${PEAK_DIR}/${MARK}.bed"
  IN_LINES=$(wc -l < "${SRC}")
  awk 'BEGIN{OFS="\t"} $1 ~ /^(chr[0-9]+|chrX|chrY)$/ {print $1,$2,$3}' "${SRC}" \
    | LC_ALL=C sort -k1,1 -k2,2n -u > "${DST}"
  OUT_LINES=$(wc -l < "${DST}")
  echo "    ${MARK}: ${IN_LINES} -> ${OUT_LINES}"
done
```

Expected post-filter line counts (only the contaminated files lose rows; clean ones are unchanged modulo dedupe):
- H3K27ac: 15,105 → ~15,090 (drops a few `chr*_random` + `chrUn_JH584304`)
- H3K27me3: 15,809 → ~15,790 (drops chrM peak + a few random/chrUn)
- H3K4me1: 113,781 → ~113,750
- H3K4me3: 6,581 → 6,581 (already clean)
- CTCF: 32,487 → ~32,475

**[2/5] Filter chromsizes** — `awk '$1 ~ /^(chr[0-9]+|chrX|chrY)$/' "${CHROMSIZES_SRC}" > "${CHROMSIZES_DST}"`. Verify exactly 21 lines. Sum of lengths ≈ 2.73 Gb.

**[3/5] Build cellmarkfiletable.txt** — five tab-separated lines via `printf`. Format: `<cell><TAB><mark><TAB><bedfile>`. The bedfile names are bare (no path prefix) because BinarizeBed receives the directory separately.

```bash
{
  printf '%s\tH3K27ac\tH3K27ac.bed\n'  "${CELL}"
  printf '%s\tH3K27me3\tH3K27me3.bed\n' "${CELL}"
  printf '%s\tH3K4me1\tH3K4me1.bed\n'  "${CELL}"
  printf '%s\tH3K4me3\tH3K4me3.bed\n'  "${CELL}"
  printf '%s\tCTCF\tCTCF.bed\n'         "${CELL}"
} > "${CMT}"
```

**[4/5] BinarizeBed** —

```bash
"${CHROMHMM}" BinarizeBed -peaks \
    "${CHROMSIZES_DST}" \
    "${PEAK_DIR}" \
    "${CMT}" \
    "${BIN_DIR}"
status=$?
[ ${status} -eq 0 ] || { echo "ERROR: BinarizeBed exited ${status}" >&2; exit ${status}; }
n_bin=$(ls "${BIN_DIR}" | wc -l)
echo "    binary files produced: ${n_bin} (expect 21)"
```

Expected output: 21 files named `cerebellum_late_chr{1..19,X,Y}_binary.txt`. Each file has a 2-line header (cell name; tab-separated mark labels) followed by N rows of 5 binary values (one per 200bp bin per chrom).

**[5/5] LearnModel(k=12)** — the long step. Use `-p 4` for parallel chromosome processing, `-l` to give proper segmentation lengths.

```bash
echo "    LearnModel started: $(date)  -- expect 30-90 min"
"${CHROMHMM}" LearnModel \
    -p 4 \
    -l "${CHROMSIZES_DST}" \
    "${BIN_DIR}" \
    "${MODEL_DIR}" \
    "${NSTATES}" \
    "${ASSEMBLY}"
status=$?
[ ${status} -eq 0 ] || { echo "ERROR: LearnModel exited ${status}" >&2; exit ${status}; }
echo "    LearnModel finished: $(date)"
```

Final verification block at end of script:

```bash
SEG_BED="${MODEL_DIR}/${CELL}_${NSTATES}_segments.bed"
EMI="${MODEL_DIR}/emissions_${NSTATES}.txt"
[ -s "${SEG_BED}" ] && echo "  OK ${SEG_BED}: $(wc -l < "${SEG_BED}") lines" \
                    || echo "  MISSING ${SEG_BED}"
[ -s "${EMI}" ] && cat "${EMI}" || echo "  MISSING ${EMI}"
echo ""
echo "  Chrom coverage in segments:"
cut -f1 "${SEG_BED}" | sort -u
echo ""
echo "  State distribution (col 4 of segments.bed):"
awk '{print $4}' "${SEG_BED}" | sort | uniq -c | sort -rn
echo ""
echo "NEXT (manual): inspect ${EMI} and write rename file at:"
echo "  ${CHM_DIR}/12state_rename_cerebellum.txt"
```

### Step 2 — `run_phase2.sh` (driver)

Mirrors `run_phase1.sh` exactly: cd to repo root, banner, call worker with status check, end banner. Keep it ~25 lines.

### Step 3 — Execute (NOT in plan mode)

```bash
bash cluster/scripts/run_phase2.sh 2>&1 | tee cluster/bap1_late/chromHMM/phase2.txt
```

### Step 4 — Manual emission interpretation (POST-AUTOMATION)

After LearnModel finishes, `emissions_12.txt` will have shape (12 states × 5 marks) of probabilities in [0,1]. Each row is one state; columns are H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF (in the order ChromHMM emits — verify by reading the header line).

Decision tree for naming each state from its emission row (probabilities relative to other states):

| Emission signature | Suggested name |
|---|---|
| K4me3 high + K27ac high | `Active_Promoter` |
| K4me3 high + K27me3 high (no/low K27ac) | `Bivalent_Promoter` |
| K4me3 high + K27ac low (no K27me3) | `Weak_Promoter` or `Poised_Promoter` |
| K27ac high + K4me1 high (no K4me3) | `Active_Enhancer` |
| K27ac very high, K4me1 lower | `Strong_Enhancer` |
| K4me1 high + K27ac low + K27me3 low | `Poised_Enhancer` |
| K4me1 high + K27me3 high | `Bivalent_Enhancer` (uncommon but possible) |
| K27me3 high alone (other marks low) | `Polycomb` |
| CTCF high alone | `Insulator` |
| CTCF high + active marks | `Active_CTCF` |
| All marks ≤ 0.05 | `Quiescent` (likely 1-3 such states; suffix as `Quiescent_1`, `Quiescent_2` to keep IDs unique, or fold to one name — both work for OverlapEnrichment) |

Format of `12state_rename_cerebellum.txt` (matches Popay's example, tab-separated, no header):

```
E1<TAB>Active_Promoter
E2<TAB>Bivalent_Promoter
...
E12<TAB>Quiescent
```

The state ID prefix (`E1` vs `1` vs `U1`) depends on what ChromHMM v1.27 actually emits in `emissions_12.txt` and `_segments.bed` col 4 — read the file before writing the rename. ChromHMM v1.27 default is `E1`-`E12` (no leading zero). Use whatever prefix appears verbatim.

Sanity for naming: if two states have nearly identical emissions (Euclidean distance < 0.1 across all 5 marks), it's fine to give them the same biological name — OverlapEnrichment in Phase 4.4 will collapse rows with identical labels in the heatmap.

---

## Verification

### Per-step checks (in worker script)
1. After staging: each `${PEAK_DIR}/<MARK>.bed` has only 3 cols, `cut -f1 | sort -u` returns subset of {chr1..chr19, chrX, chrY}.
2. After chromsizes filter: `wc -l ${CHROMSIZES_DST}` returns exactly 21.
3. After BinarizeBed: `ls ${BIN_DIR} | wc -l` returns exactly 21.
4. After LearnModel: `${MODEL_DIR}/cerebellum_late_12_segments.bed` exists and is non-empty; `${MODEL_DIR}/emissions_12.txt` exists.

### Post-hoc model sanity (printed by worker)
- All 21 standard chroms appear in `cut -f1 segments.bed | sort -u`.
- State distribution: `awk '{print $4}' segments.bed | sort | uniq -c` should show 12 distinct states. **No state >50% of bins** (would mean quiescent-collapse — probably means peak BEDs were over-filtered or one mark dominates noise). **No state <0.5% of bins** (over-fitting; harmless but uninterpretable).
- segments.bed line count: expect ~0.5M–5M lines (region-merged across 200bp bins; the original PLAN-CLUSTER estimate of "~14M" was per-bin, not per-segment — segments.bed merges adjacent identical-state bins).

### Cross-check vs Popay
`cluster/clustering_example_data/12state_rename.txt` was learned from 12 ENCODE marks in RPE-1 cells. We don't have H3K36me3, so we will not have "Transcriptional elongation" / "Transcriptional transition" states. We will have a richer Polycomb/CTCF/enhancer breakdown. Different state vocabulary is expected and biologically correct.

### Smoke test for Phase 4.4 readiness
After writing rename file, run a minimal OverlapEnrichment as a smoke test (optional — not part of Phase 2 deliverable, but cheap):

```bash
# Only if user wants to validate before Phase 3:
echo "${REPO_ROOT}/cluster/data/late_merged_loop_counts.txt" > /tmp/coordlist.txt
"${CHROMHMM}" OverlapEnrichment -noimage -uniformscale \
    -m "${CHM_DIR}/12state_rename_cerebellum.txt" \
    -f /tmp/coordlist.txt \
    -colfields 0,1,2 \
    "${MODEL_DIR}/cerebellum_late_12_segments.bed" \
    /tmp/ \
    /tmp/all_loops_smoke
# Inspect /tmp/all_loops_smoke.txt — fold-enrichment values should be finite,
# states with active marks should be enriched at K27ac-rich loop coords.
```

---

## Critical files referenced (verified to exist on disk)

| Role | Path |
|---|---|
| ChromHMM jar | `cluster/ChromHMM/ChromHMM.jar` |
| ChromHMM wrapper | `cluster/ChromHMM/chromhmm` |
| Source chromsizes | `cluster/ChromHMM/CHROMSIZES/mm10.txt` |
| Peak BED — H3K27ac | `peaks/beds/H3K27acCerebellumLate2.bed` |
| Peak BED — H3K27me3 | `peaks/beds/H3K27me3CerebellumLate1.bed` |
| Peak BED — H3K4me1 | `peaks/beds/H3K4me1CerebellumLate1.bed` |
| Peak BED — H3K4me3 | `peaks/beds/H3K4me3CerebellumLate2.bed` |
| Peak BED — CTCF | `peaks/CTCF.bed` |
| Style reference | `cluster/scripts/run_phase1.sh` |
| Rename format reference | `cluster/clustering_example_data/12state_rename.txt` |
| Project taxonomy reference (for naming) | `CLAUDE.md` §"7-category" + `scripts/annotate_loops_extended.R` |
| Phase 1 outputs (consumers of segmentation in Phase 4.4) | `cluster/data/late_merged_loop_counts.txt`, `cluster/data/late_merged_loop_metadata.tsv` |

---

## Followups (not in scope for Phase 2)

- **Phase 3** (k-means clustering) reads only `cluster/data/late_merged_loop_counts.txt` — independent of Phase 2 outputs, can run in parallel if desired.
- **Phase 4.4** (anchor-vs-span ChromHMM enrichment) is the headline consumer; needs `cerebellum_late_12_segments.bed`, `12state_rename_cerebellum.txt`, AND `combined-clusters.txt` from Phase 3.
- **Phase 4.5** (ChromHMM proportions stacked bar) — same inputs.
- The rename file's biological names will get reused as palette keys in `chromHMM_heatmap.heatmap_plot()` and `plotting.stacked()`; Phase 4 plotting code may need a palette dict that lists the names we choose here. Note the chosen names in the eventual Phase 4 plan so palette/legend code matches exactly.
