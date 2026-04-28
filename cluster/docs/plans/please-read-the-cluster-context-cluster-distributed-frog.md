# Phase 4 Plan: Downstream Grouped Analyses for BAP1-KO Hi-C Loop Clusters

## Context

The Popay anchor-vs-span ChromHMM analysis is the **#1 outstanding item** for the BAP1-KO Hi-C paper, identified by both Jesse Dixon and Tessa Popay as the highest-priority unresolved question (CONTEXT-CLUSTER §2, §3). It mechanistically distinguishes two competing models:

- **Anchor-specific Polycomb enrichment** → sensitivity model: CTCF binding sites directly disrupted by repressive chromatin
- **Loop-body/span Polycomb enrichment** → extrusion impediment model: Polycomb spreading blocks cohesin loop extrusion

Phases 0–3 are complete (Popay module bug fixes, loop count file + mm10 gene BED, 12-state ChromHMM segmentation, k=6 k-means clusters with v2 ratio bounds). Phase 4 takes the canonical 6 clusters at `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt` (38,948 loops) and runs 8 grouped downstream analyses adapted from Popay's `grouped_loops_figures.ipynb`. The KEY result is **4.4: ChromHMM anchor-vs-span enrichment heatmaps** (Popay Fig 2f equivalent) — the paper-ready figure that answers the central biological question.

User-confirmed scope (this round):
- **BigWig source:** `/Users/zakiralibhai/sdsc/bigwigs/` exclusively (mariner_hi-c `macs2.narrow.aug18.dedup` paths are partially 0-byte — discard).
- **Phase 4.3 marks:** H3K27ac, H3K27me3, H2AK119ub, H3K27me1 (4 marks × ctrl+mut = 8 BigWigs).
- **Phase 4.7 DiffBind:** emit BOTH project-standard (FDR<0.05, |Fold|>0.3) AND Popay-default (FDR<0.05, |Fold|>0) figures.
- **Implementation:** all 8 sub-analyses (4.1–4.8) in one round, gated by `--analyses` CLI flag.

---

## Architecture

**One Python script** at `cluster/scripts/05_grouped_analyses.py` orchestrating 8 sub-analysis functions, **one driver** at `cluster/scripts/run_phase4.sh`, **one utility patch** to `cluster/scripts/utils/multi_format_output.py`. Mirrors run_phase3.sh exactly (LOG env var → tee, set -e only, no -uo pipefail per project memory).

**Sub-analysis execution order** (designed to fail-fast on KEY result, then cheap-fast, then expensive):
1. **4.4** ChromHMM anchor-vs-span (KEY — Fig 2f equivalent)
2. **4.1** Loop size + 4.8 cluster×differential crosstab (cheap, validates clustering)
3. **4.5** ChromHMM proportions + 4.6 gene annotation
4. **4.2** Loop classification + 4.7 DiffBind relationship
5. **4.3** ChIP signal at anchors (slowest — 8 BigWigs × ~38k anchors)

Shared inputs are loaded once and passed by reference to each sub-analysis function: `bedpe_df_sorted` (DataFrame, 38,948 × 9), `group_dict` (dict of 6 DataFrames), `cluster_order` (list `['clust1', ..., 'clust6']`), `palette` (state name → color dict, derived dynamically from rename file).

---

## Phase 4.0 — Patch `multi_format_output.py` to also intercept `Figure.savefig`

**Why:** `plotting.stacked()` (used in 4.2 and 4.7) calls `fig.savefig(...)` on `cluster/plotting.py:630-631`, NOT `plt.savefig(...)`. The current `multi_format_savefig()` only patches `plt.savefig`, so stacked-bar outputs would bypass 4-format emission.

**Edit:** Extend `multi_format_savefig()` to additionally patch `matplotlib.figure.Figure.savefig` at the class level. Within the patch:
- Trigger on PNG-extension target only (skip when called with .pdf to avoid double-emission, since `stacked()` calls savefig twice — once for .png, once for .pdf — and we want one round of 4-format emission).
- Restore the original `Figure.savefig` in `finally` block.

**File:** `cluster/scripts/utils/multi_format_output.py`

---

## Phase 4.1 — Loop size

**Function:** `run_loop_size(group_dict, out_dir)`

**Logic:** Pre-create `cluster/bap1_late/figures/loop_size/` subfolder (per `figure_subfolder()` pattern), then call `bedpe_analysis.loop_size(out_dir=str(subfolder), bedpe_dict=group_dict, logY=False)` inside `multi_format_savefig()`. Function internally calls `plotting.box()` then `plotting.strip()` on `(max(x1,y1) - min(x1,y1)) / 1000` per cluster, plus writes `loop_size.stats.txt` (Kruskal-Wallis + pairwise Wilcoxon).

**Outputs:** `cluster/bap1_late/figures/loop_size/loop_size.{png,pdf,svg,jpg}`, `loop_size_strip.{png,pdf,svg,jpg}`, `loop_size.stats.txt`.

**Expected biology:** clust6 (loss-dominant) median loop size > clust5 (gain-dominant). Mirrors the abstract's "lost loops median 625 kb vs gained 320 kb."

---

## Phase 4.2 — Loop classification (CTCF-only, no RAD21)

**Function:** `run_loop_classification(group_dict, out_dir)`

**Adaptations from Popay cell-5:**
- **Drop RAD21** from `classifier_dict` (we have no RAD21 ChIP).
- **Reclassify** with CTCF-only structural rule:
  - structural: `(CTCF == 2) & (EorP < 2)`
  - CRE: `(CTCF < 2) & (EorP == 2)`
  - mixed: `(CTCF == 2) & (EorP == 2)`
  - unclassified: everything else
- **Peak file column count mismatch:** `peaks/CTCF.bed` is 10-col narrowPeak, `peaks/beds/H3K27acCerebellumLate2.bed` is 6-col. Use `usecols=[0,1,2]` for all peak loads (uniform 3-col).
- **Inputs:**
  - CTCF: `peaks/CTCF.bed`
  - enhancer: `peaks/beds/H3K27acCerebellumLate2.bed`
  - promoter: `cluster/data/mm10_knownGene_pp.bed`
- Run `chi2_contingency` on the resulting 4 × 6 count table.
- `plotting.stacked(count_table=classified_df.transpose(), ...)` for the 100% stacked bar (transpose so cluster is x-axis, classification stacks vertically).

**Output:** `cluster/bap1_late/figures/loop_classification/loop_classification.{png,pdf,svg,jpg}` + `loop_classification.stats.txt` (chi-squared output).

---

## Phase 4.3 — ChIP signal at anchors (raw RPKM, no RAD21 normalization)

**Function:** `run_anchor_chip(group_dict, out_dir, bigwig_dict)`

**Adaptations from Popay cell-7:**
- **No `bed_intersect` filter** — we lack RAD21 peaks. Compute mean BigWig signal at every anchor region directly.
- **Critical: dedupe anchors before pyBigWig query.** A single 10kb anchor region may participate in multiple loops (hub anchors). Without dedup the same region is queried N times, over-weighting hubs in per-cluster means. Apply `.drop_duplicates(subset=anchor)` on the union of anchor1 + anchor2 BEFORE the pyBigWig loop. Re-merge back to loops afterward.
- **No RAD21 normalization** (source data unavailable). Report raw RPKM.
- **8 BigWigs from `/Users/zakiralibhai/sdsc/bigwigs/`:**
  - H3K27ac: `H3K27acCtrl.bw`, `H3K27acMut.bw`
  - H3K27me3: `H3K27me3Ctrl.bw`, `H3K27me3Mut.bw`
  - H2AK119ub: `H2AK119ubCtrl.bw`, `H2AK119ubMut.bw`
  - H3K27me1: `H3K27me1Ctrl.bw`, `H3K27me1Mut.bw`
- **Per loop:** mean of (anchor1 mean signal, anchor2 mean signal). Average via `.groupby(anchor).mean()` to handle sub-anchor positional duplicates.
- **Stats:** for each mark, run `statistics_functions.kruskal_wilcoxon(data_names=['{mark}_{ctrl/mut}_clust1', ..., '{mark}_{ctrl/mut}_clust6'], data_list=...)` — 12-way comparison per mark, 4 marks total.
- **Plot:** single faceted boxplot via `plotting.box(melted_df, xcol='group', ycol='signal', subplots=True, subplot_col='ChIP_condition', ...)`. The subplot column encodes mark × condition (8 panels: K27ac_ctrl, K27ac_mut, K27me3_ctrl, K27me3_mut, K119ub_ctrl, K119ub_mut, K27me1_ctrl, K27me1_mut). Pass `Y_range=None` for auto-fit (Popay's `[0,1]` was tuned for hg38 RAD21-normalized data).
- **Pitfall guard:** if any panel ends up with zero rows, `plot_defaults_b` may throw on `get_yticks()[1]`. Filter empty subsets before plotting.

**Outputs:**
- `cluster/bap1_late/figures/ChIP_intersect/anchor_ChIP_box/anchor_ChIP_box.{png,pdf,svg,jpg}`
- `cluster/bap1_late/figures/ChIP_intersect/anchor_ChIP.stats.txt`

---

## Phase 4.4 — ChromHMM anchor-vs-span enrichment **(KEY — Fig 2f equivalent)**

**Function:** `run_chromhmm_anchor_vs_span(group_dict, out_dir, segment_path, rename_path)`

**Logic (per Popay cell-23):**
1. For each cluster, write two BED files:
   - **Anchor BED** (`cluster/bap1_late/chromHMM/anchor_input/{clust}.bed`): concatenate `df[['chr1','x1','x2']]` + `df[['chr2','y1','y2']].rename(columns={'chr2':'chr1','y1':'x1','y2':'x2'})`. Both anchors pooled (validated by Plan agent — concatenation is correct, separating would halve sample sizes and break OverlapEnrichment fold-enrichment math).
   - **Span BED** (`cluster/bap1_late/chromHMM/span_input/{clust}.bed`): `df[['chr1','x1','y2']]` (full loop extent). Pre-validate: assert `(df['chr1'] == df['chr2']).all()` to confirm cis-only loops; chrM/trans loops would produce meaningless span coords.
2. Write `coordlistfile.txt` listing the 6 cluster BED filenames in order. Controls heatmap column order.
3. For each of `{anchor, span}`, run ChromHMM OverlapEnrichment as `subprocess.run(..., check=True, capture_output=True)`:
   ```bash
   chromhmm OverlapEnrichment \
       -noimage -uniformscale \
       -m cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt \
       -f cluster/bap1_late/chromHMM/coordlistfile.txt \
       -colfields 0,1,2 \
       cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed \
       cluster/bap1_late/chromHMM/{anchor,span}_input/ \
       cluster/bap1_late/chromHMM/{anchor,span}
   ```
   - Use the `chromhmm` wrapper on PATH (per CONTEXT-CLUSTER §8 / `~/.zshrc`).
   - On non-zero exit, surface stderr and abort the analysis (no silent FileNotFoundError downstream).
   - Verify both `anchor.txt` and `span.txt` exist + are non-empty before proceeding.
4. Inside `multi_format_savefig()`, call `chromHMM_heatmap.heatmap_plot(path=anchor_txt, normalize=False)` then again for `span_txt`. Function auto-derives output filename from input stem and writes to the input's directory (no `out_dir` arg). So heatmaps land at `cluster/bap1_late/chromHMM/{anchor,span}.{png,pdf,svg,jpg}` — adjacent to the .txt enrichment tables.
5. **Note:** `heatmap_plot()` calls `plt.savefig` twice internally (once per format). Our patch only fires 4-format emission on the .png call (.pdf call passes through cleanly), so we get exactly one set of 4 outputs per heatmap. Validated.

**Outputs:**
- `cluster/bap1_late/chromHMM/{anchor,span}.txt` — 12 × 6 fold enrichment matrices
- `cluster/bap1_late/chromHMM/{anchor,span}.{png,pdf,svg,jpg}` — Fig 2f-equivalent heatmaps

**Expected biology:** If extrusion-impediment model holds, clust6 (loss-dominant) should show high Polycomb enrichment in **span** but normal in **anchor**, mirroring Popay's mixed-dependency cluster (Fig 2g).

---

## Phase 4.5 — ChromHMM proportions

**Function:** `run_chromhmm_proportions(bedpe_df_sorted, out_dir, segment_path, rename_path)`

**Adaptations from Popay cell-25:**
- **CRITICAL FIX:** `rename_df['sort_col'].str.replace('U','').astype(int)` will fail with our E-prefixed states. Change `'U'` → `'E'` (or use a regex that strips any leading letter prefix).
- **Use `plotting.bar()` not `plotting.stacked()`** — cell-25's active code path is `plotting.bar(melted_df, xcol='index', ycol='proportion', hue_col='state', palette=palette, ...)`. The stacked() call is commented out.
- **Excluded states:** Popay excludes `U12` (Heterochromatin) and `U9` (Low signal). Our equivalents are E1 (Quiescent, 92% of genome) and arguably none other since our 12 states are all biologically meaningful. **Exclude only E1 (Quiescent)** to avoid the proportions being dominated by background.
- **Dynamic palette construction** from rename file (per Plan agent guidance — Popay's hardcoded space-name palette won't match our underscore-name states). Build a 12-color dict at runtime keyed by our biological state names:
  ```python
  state_colors = {
      'Active_Promoter': 'tab:red',
      'Active_Promoter_Flank': 'salmon',
      'Poised_Promoter': 'tab:orange',
      'Weak_Promoter': 'gold',
      'Strong_Enhancer': 'tab:purple',
      'Active_Enhancer': 'tab:green',
      'Poised_Enhancer': 'lightgreen',
      'Bivalent_Enhancer': 'olive',
      'Polycomb': 'tab:cyan',
      'CTCF_Boundary': 'tab:gray',
      'Insulator': 'lightgray',
      # Quiescent excluded
  }
  ```
- **Method:** `bf.overlap` of cluster anchors with segments → per-loop most-prominent state by overlap length → proportion table per cluster. Cell-25's pattern reused verbatim with the fixes above.

**Output:** `cluster/bap1_late/figures/chromHMM_anchor/chromHMM_anchor.{png,pdf,svg,jpg}` (grouped bar).

---

## Phase 4.6 — Gene annotation

**Function:** `run_gene_annotation(group_dict, out_dir)`

**Logic:**
- Call `bedpe_analysis.bedtools_annotation(out_dir=out_dir/'annotation', bedpe_dict=group_dict, FPKM_df=None, temp_dir=str(out_dir/'annotation/_temp'))`.
- **CRITICAL:** `temp_dir` must be a real path (not None) or function crashes at line 64. Per Phase 0 fix the `FPKM_df=None` branch now skips the FPKM blocks safely.
- Function returns `annotation` dict (deep copy of group_dict with gene columns merged). Iterate to write per-cluster TSV: `for group, df in annotation.items(): df.to_csv(out_dir/'annotation'/f'{group}_annotation.txt', sep='\t', index=False)`.
- Verify `bedtools` is on PATH inside the cluster conda env at the start of the function (`shutil.which('bedtools')` check + clear error if missing).

**Outputs:**
- `cluster/bap1_late/figures/annotation/{clust1..clust6}_annotation.txt`
- `cluster/bap1_late/figures/annotation/_temp/` (transient — keep for debug, easy to clean later)

---

## Phase 4.7 — DiffBind relationship

**Function:** `run_diffbind(group_dict, out_dir, diffbind_dict, fc_cutoffs=[0.3, 0.0])`

**Per user request: emit BOTH thresholds (project-standard 0.3 + Popay default 0.0) as separate figures.**

**Adaptations from Popay cells 17 + 19:**
- **Column rename at load time** to avoid the 8-10 scattered references the Plan agent flagged. After loading each diffbind file, `df.rename(columns={'Peak_Chr':'Chr', 'Peak_Start':'Start', 'Peak_End':'End'})` — then the rest of Popay's code works unchanged.
- **CRITICAL FIX:** cell-17 has `value_vars=['sig_chr1','sig_chr1']` (typo — both are chr1, anchor2 silently ignored). Fix to `['sig_chr1','sig_chr2']`.
- **Three input files** at `peaks/diffbind/`:
  - `K27ac_diffbind_results_summit_appended_ap.txt`
  - `K27me3_diffbind_results_summit_appended_ap.txt`
  - `K119ub_diffbind_results_summit_appended_ap.txt`
- For each (diffbind file, FC cutoff) pair:
  - Tag peaks: `decreased` (FDR<0.05 & Fold<-cutoff), `increased` (FDR<0.05 & Fold>cutoff), `non-sig.` otherwise.
  - **Cell-17 logic** (proportions): per cluster, count anchors overlapping `sig` / `non-sig` / `no peak` peaks → 3 × 6 table → `chi2_contingency` → `plotting.stacked(count_table.transpose())`.
  - **Cell-19 logic** (FC magnitudes): per cluster, mean Fold across both anchors → 6-way Kruskal-Wallis + pairwise Wilcoxon → `plotting.box`.

**Outputs (per mark, per FC cutoff — 6 figure sets total):**
- `cluster/bap1_late/figures/ChIP_intersect/differential_binding_{mark}_fc{cutoff}/differential_binding_{mark}_fc{cutoff}.{png,pdf,svg,jpg}` (proportions)
- `cluster/bap1_late/figures/ChIP_intersect/ChIP_FC_{mark}_fc{cutoff}/ChIP_FC_{mark}_fc{cutoff}.{png,pdf,svg,jpg}` (boxplots)
- `cluster/bap1_late/figures/ChIP_intersect/diffbind_stats_{mark}_fc{cutoff}.txt` (chi-squared + Kruskal-Wallis + Wilcoxon)

3 marks × 2 cutoffs × (proportions + boxplots) = 12 figure subfolders + 6 stats files.

---

## Phase 4.8 — Cluster × differential status crosstab (NEW)

**Function:** `run_cluster_diff_crosstab(bedpe_df_sorted, metadata_path, out_dir)`

Not in Popay's pipeline. Validates the clustering against edgeR's independent differential calls.

**Logic:**
1. Load `cluster/data/late_merged_loop_metadata.tsv` (16 columns including `direction`).
2. Inner-join `bedpe_df_sorted` with metadata on the 6 coordinate columns (`chr1, x1, x2, chr2, y1, y2`). Verify zero loss (Phase 1 confirmed clean join).
3. Build a 6 × 3 contingency table: rows = clust1..clust6, cols = `down_in_mutant`/`unchanged`/`up_in_mutant`.
4. `chi2_contingency(crosstab)` → expect p≈0 given Phase 3's clean separation (clust5: 97% up, clust6: 78.5% down).
5. Plot `plotting.stacked(count_table=crosstab.transpose(), out_name='cluster_differential_status')` — 100% stacked bars per cluster colored by direction.
6. Save the raw + percent crosstab as TSV alongside the figure.

**Outputs:**
- `cluster/bap1_late/figures/cluster_differential_status/cluster_differential_status.{png,pdf,svg,jpg}`
- `cluster/bap1_late/figures/cluster_differential_status.stats.txt` (chi-squared + crosstab)

---

## Driver script: `cluster/scripts/run_phase4.sh`

Mirrors `run_phase3.sh` exactly:
- `set -e` only (no -uo pipefail — project memory).
- `cd "$(dirname "$0")/../.."` to repo root.
- `LOG=${LOG:-cluster/phase4.txt}`, body wrapped in `{ ... } 2>&1 | tee "$LOG"`.
- `PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3` (per Phase 1 corrections — conda activation unreliable).
- Positional arg passes through to `--analyses` (default `all`).
- Verification block at end: `ls -lh` for all 8 expected output groups.

Invocation:
```bash
LOG=cluster/phase4.txt bash cluster/scripts/run_phase4.sh             # all 8
LOG=cluster/phase4_test.txt bash cluster/scripts/run_phase4.sh 4.4    # KEY result only
LOG=cluster/phase4_chip.txt bash cluster/scripts/run_phase4.sh 4.3,4.7  # subset
```

---

## Files to Create

| File | Purpose |
|------|---------|
| `cluster/scripts/05_grouped_analyses.py` | Phase 4 main script — 8 sub-analysis functions + argparse + shared loader |
| `cluster/scripts/run_phase4.sh` | Driver script — LOG env var + tee, mirrors run_phase3.sh |

## Files to Modify

| File | Change |
|------|--------|
| `cluster/scripts/utils/multi_format_output.py` | Add `Figure.savefig` patch alongside existing `plt.savefig` patch (Phase 4.0). PNG-extension trigger only to avoid double emission. Restore in `finally`. |

## Files to Read/Reuse (no changes)

| File | Why |
|------|-----|
| `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt` | Input — 38,948 clustered loops |
| `cluster/data/late_merged_loop_metadata.tsv` | Input — edgeR direction labels for 4.8 |
| `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` | Input — 4.4, 4.5 |
| `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` | Input — 4.4, 4.5 (state naming) |
| `cluster/data/mm10_knownGene_pp.bed` | Input — 4.2 (promoter overlap), 4.6 (gene annotation) |
| `peaks/CTCF.bed`, `peaks/beds/H3K27acCerebellumLate2.bed` | Input — 4.2 (loop classification) |
| `peaks/diffbind/{K27ac,K27me3,K119ub}_diffbind_results_summit_appended_ap.txt` | Input — 4.7 |
| `/Users/zakiralibhai/sdsc/bigwigs/{H3K27ac,H3K27me3,H2AK119ub,H3K27me1}{Ctrl,Mut}.bw` | Input — 4.3 |
| `cluster/bedpe_analysis.py` | Reuse `loop_size`, `bedtools_annotation` |
| `cluster/plotting.py` | Reuse `box`, `strip`, `stacked`, `bar` |
| `cluster/chromHMM_heatmap.py` | Reuse `heatmap_plot` |
| `cluster/statistics_functions.py` | Reuse `kruskal_wilcoxon` |
| `cluster/cluster_tools.py` | Reuse `sort_by_strings` |
| `cluster/scripts/04_clustering.py` | Reference for argparse / logging / palette / multi_format_savefig pattern |
| `cluster/scripts/utils/multi_format_output.py` | Use `multi_format_savefig`, `figure_subfolder` (after Phase 4.0 patch) |

---

## Verification Plan

After running `LOG=cluster/phase4.txt bash cluster/scripts/run_phase4.sh`:

1. **Phase 4.0 (utility patch):** Smoke test — call `plotting.stacked` inside the `with multi_format_savefig():` block on a tiny synthetic 3×3 count table, confirm 4 sibling files (.png, .pdf, .svg, .jpg) emitted, `Figure.savefig` restored after exit (assert `matplotlib.figure.Figure.savefig is _original`).
2. **Phase 4.1:** `loop_size.stats.txt` shows clust6 median > clust5 median (loss-dominant cluster has longer loops). All 4 figure formats present.
3. **Phase 4.2:** Stacked bar shows ≥1 of `(structural, CRE, mixed, unclassified)` differs across clusters; chi-squared p < 0.001.
4. **Phase 4.3:** Faceted boxplot has 8 panels (4 marks × 2 conditions). H3K27me3 panel for clust6 (loss-dominant) is visibly elevated relative to clust5 (gain-dominant), or H2AK119ub similarly elevated. Stats file lists 4 omnibus + ≤60 pairwise tests.
5. **Phase 4.4 (KEY):** Both `anchor.txt` and `span.txt` produce a 12-row × 6-column matrix. Heatmaps differ between anchor and span (specifically: Polycomb row should be more enriched in span than anchor for at least one cluster, supporting extrusion-impediment model). Verify by `awk '$1 == "Polycomb"' anchor.txt span.txt`.
6. **Phase 4.5:** Grouped bar shows 11 states (excluding Quiescent) with proportions summing to 100% per cluster. clust6 should have visibly higher Polycomb proportion than clust5.
7. **Phase 4.6:** 6 per-cluster annotation TSVs present, each with non-zero gene count and `gene_name` column.
8. **Phase 4.7:** 12 figure subfolders (3 marks × 2 cutoffs × 2 plot types) + 6 stats files. Project-standard threshold should produce qualitatively similar but cleaner figures than Popay default.
9. **Phase 4.8:** chi-squared p ≪ 0.001 (expected given Phase 3 v2's 97%/78.5% directional capture). Stacked bar shows clust5 dominantly green (up_in_mutant), clust6 dominantly orange (down_in_mutant).

After all 8 succeed, append a "Phase 4 — DONE" section to `cluster/PLAN-CLUSTER.md` (in `plan-p2.md` style: status summary, corrections-to-plan during execution, file inventory, verification snapshot). Following `cluster/feedback_phase_log_capture.md` memory: tee output is captured at `cluster/phase4.txt` for later reference.
