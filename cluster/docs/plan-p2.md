
## Phase 4: Downstream Analyses — DONE (2026-04-27)

**Status:** All 8 sub-analyses (4.1–4.8) executed end-to-end via `LOG=cluster/phase4.txt bash cluster/scripts/run_phase4.sh all`. Single Python orchestrator at `cluster/scripts/05_grouped_analyses.py` with `--analyses` CLI flag for selective re-runs (e.g. `bash run_phase4.sh 4.4` for KEY-only). KEY result (4.4) shows **clean differential ChromHMM enrichment distinguishing two mechanisms in the same paper**: gained loops (clust5) have Polycomb at BOTH anchor (6.59×) AND span (3.03×) — Polycomb-Polycomb domain expansion / extrusion impediment; lost loops (clust6) have Polycomb at anchor (2.09×) but NOT span (0.94×) — anchor-disruption / sensitivity model. This is a stronger result than Popay's mixed-dependency cluster which only showed span enrichment. Total runtime ~6 minutes on Mac (07:11–07:17 PDT). All 8 sub-analyses produced multi-format outputs (PNG + PDF + SVG + JPG) per figure subfolder.

**Corrections applied during execution (Phase 5+ must follow these, not the original plan):**

- **`multi_format_savefig` must patch `Figure.savefig`, not `plt.savefig`.** Original Phase 3 utility patched only `plt.savefig`, leaving `plotting.stacked()` (which uses `fig.savefig` directly at `cluster/plotting.py:630-631`) bypassed — outputs would have only emitted PNG+PDF, no SVG/JPG. Patched class-level `matplotlib.figure.Figure.savefig` instead — covers BOTH calling patterns since `pyplot.savefig` internally calls `fig.savefig`. Phase 3 outputs unaffected (re-running 04_clustering.py with the new wrapper goes through the same code path with identical 4-format output). Smoke-tested with `plt.savefig('foo.png')` AND `fig.savefig('foo.png')` — both emit 4 sibling files; `Figure.savefig` is restored on context exit.
- **`pd.crosstab` produces a named index that breaks `plotting.stacked`.** The function does `count_table.reset_index(drop=False).rename(columns={'index':'xcol'})` which only works when `index.name is None`. `pd.crosstab(joined['GROUP'], joined['direction'])` sets `index.name='GROUP'` → `KeyError('xcol')` in the upstream melt. **Fix in 4.8:** clear `pct.index.name = None` and `pct.columns.name = None` before passing to `plotting.stacked`. Same defensive pattern needed for any future caller using `pd.crosstab` → `plotting.stacked`. Manually-built DataFrames (via `.loc[]` like 4.2 / 4.7) are unaffected since their index defaults to unnamed.
- **chromhmm wrapper on PATH works as-is.** OverlapEnrichment subprocess call succeeded on first run; no fallback to direct `java -jar` invocation needed. The wrapper at `cluster/ChromHMM/chromhmm` is on PATH via `~/.zshrc`. Driver verifies via `which chromhmm` in its banner.
- **`bedtools` is system-wide on macOS.** `/opt/homebrew/bin/bedtools` (v2.31.1) used by `bedpe_analysis.bedtools_annotation()` for 4.6 — no need to install in the cluster conda env. Check `shutil.which('bedtools')` at the start of `run_gene_annotation()` for clear error if missing.
- **Phase 4.5 'U' → 'E' regex fix.** Popay cell-25's `rename_df['sort_col'].str.replace('U','').astype(int)` would crash on our E-prefixed states (`E1`–`E12`). Replaced with `str.replace(r'^[A-Z]', '', regex=True)` for symbol-agnostic prefix stripping — works whether ChromHMM ever changes the state-id prefix.
- **Phase 4.7 anchor2 typo fix.** Popay cell-17's `value_vars=['sig_chr1','sig_chr1']` (both chr1; anchor2 silently ignored). Fixed to `['sig_chr1','sig_chr2']`. Material correction — DiffBind proportions stats now reflect both anchors. Plus extended the Popay 3-stack `(non-sig, sig, no peak)` to 4-stack `(no peak, non-sig., decreased, increased)` for direction visibility.
- **DiffBind column rename at load time.** Renamed our `Peak_Chr/Peak_Start/Peak_End` → `Chr/Start/End` once when loading each diffbind file, so the rest of Popay's cell-17/19 logic (`bf.overlap(..., cols2=['Chr','Start','End'])`, `drop(columns=['Chr_chr1', ...])`) works unchanged. Avoids 8-10 scattered `cols2=...` and `drop(columns=...)` updates the Plan-agent flagged.
- **All BigWigs sourced from `/Users/zakiralibhai/sdsc/bigwigs/`.** Plan §4.3 originally referenced `peaks/bigwigs/macs2.narrow.aug18.dedup/` — but two mutant files there are 0-byte (`index_19_mut_1_H3K27ac` and `index_23_mut_1_H3K27me3`). sdsc has all 16 marks at full coverage (ATAC, DNAmethylation, H3K27ac/me1/me3, H2AK119ub, H3K4me3, RNA — all paired ctrl+mut). **Phase 5 should also use sdsc as the canonical BigWig source.**
- **Anchor-side dedup before pyBigWig query in 4.3.** Within each anchor side (chr1 or chr2 columns), `.drop_duplicates(subset=anchor_cols)` before iterating `pyBigWig.values()`. Hub anchors participate in many loops; without dedup they get over-weighted in per-cluster means. ~38k loops → ~30k unique anchor1 + ~30k unique anchor2 queries × 8 BigWigs in ~5 minutes total.
- **Defensive try/except + `np.nanmean` around `pyBigWig.values()`.** Anchor regions on chrM or in blacklist may not be in the BigWig — `.values()` raises `RuntimeError`. Wrap in `try/except (RuntimeError, ValueError)` and return `np.nan`; downstream `dropna(subset=['signal'])` filters them out. NaN-only regions handled via `np.nanmean(arr)` (returns NaN if all-NaN, no warning).
- **`bedpe_analysis.bedtools_annotation` requires non-None `temp_dir`.** Despite Phase 0 fix for `FPKM_df=None`, the `os.makedirs(temp_dir, exist_ok=True)` at line 64 still crashes if `temp_dir=None` (TypeError). Pass an explicit `cluster/bap1_late/figures/annotation/_temp/` path. Temp dir kept on disk for debug; safe to clean later.
- **Python 3.8 in cluster env requires `Tuple[...]`/`Dict[...]`/`List[...]` from `typing`.** PEP 585 generic syntax (`tuple[X]`, `dict[X, Y]`, `list[X]`) requires Python 3.10+. Plan-p1.md called this out for Phase 3; same applies here.
- **Cis-only assertion guards span BED in 4.4.** `df[['chr1','x1','y2']]` is biologically meaningless for trans loops (chr1 ≠ chr2). Verified 0 trans loops in our 38,948 cluster set, but the assertion stays as a defensive check for future runs / new timepoints.
- **`chromHMM_heatmap.heatmap_plot()` writes to dir-of-input-path with no `out_dir` argument.** Auto-derives output filename from the input `.txt` stem (`os.path.dirname(path) + '/' + Path(path).stem`). So heatmaps land at `cluster/bap1_late/chromHMM/{anchor,span}.{png,pdf}` — adjacent to the .txt enrichment tables. Our patch adds .svg + .jpg via the .png save (the .pdf save passes through cleanly with non-PNG extension). No file-relocation needed.

**Script:** `cluster/scripts/05_grouped_analyses.py` (adapted from grouped_loops_figures.ipynb)

All analyses use `group_dict` built from `combined-clusters.txt` split by GROUP column.

### 4.1 — Loop size

**Reuse:** `bedpe_analysis.loop_size(out_dir, bedpe_dict=group_dict)` — no changes needed.
**Output:** `cluster/bap1_late/figures/loop_size.{stats.txt,png,pdf}`

### 4.2 — Loop classification (CTCF-only, no RAD21)

Popay's classification requires CTCF+RAD21 at both anchors for "structural". We have no RAD21 data, so use CTCF alone as the structural marker.

Modify the classifier_dict to omit RAD21 and adjust the classification rules:
- **structural:** CTCF at both anchors, no enhancer/promoter at either (`CTCF==2 & EorP<2`)
- **CRE:** enhancer or promoter at both anchors, no structural marker (`CTCF<2 & EorP==2`)
- **mixed:** one anchor CTCF, other anchor CRE (reflected by remaining combinations)
- **unclassified:** everything else

The bioframe overlap logic stays identical — just remove RAD21 from `classifier_dict` and replace `(loop_df['CTCF'] == 2) & (loop_df['RAD21'] == 2)` with `(loop_df['CTCF'] == 2)`.

**Input files:**
- CTCF: `peaks/CTCF.bed` (32,487 peaks)
- Enhancer: `peaks/beds/H3K27acCerebellumLate2.bed` (15,105 peaks)
- Promoter: `cluster/data/mm10_knownGene_pp.bed` (Phase 1 output)

**Output:** `cluster/bap1_late/figures/loop_classification.{png,pdf}`

### 4.3 — ChIP signal in anchors (raw RPKM, no RAD21 normalization)

Popay normalizes all ChIP RPKM to RAD21 signal at RAD21 peaks overlapping the anchors. Since we lack RAD21:

**Approach:** Extract raw RPKM signal directly at loop anchors using pyBigWig (mean signal per anchor, then mean of both anchors per loop). Skip the bed_intersect filtering step entirely — compute signal at ALL anchor regions, not just those overlapping a reference peak. This is biologically valid: we're asking "what is the histone mark level at loop anchors per cluster?" rather than "what is the mark level at cohesin-bound loop anchors?". The structural information is already captured by the CTCF-based classification in 4.2.

Per-mark boxplots across clusters (Kruskal-Wallis + pairwise Wilcoxon), same statistical framework as Popay.

**Input BigWigs (verified on disk at `peaks/bigwigs/macs2.narrow.aug18.dedup/`):**
- ctrl H3K27ac: `index_13_ctrl_1_H3K27ac_S25_L001_aligned_reads.sorted.bw`
- ctrl H3K27me3: `index_25_ctrl_1_H3K27me3_S37_L001_aligned_reads.sorted.bw`
- mut H3K27ac: `index_19_mut_1_H3K27ac_S46_L002_R2_001.fastq.gz_aligned_reads.sorted.bw`
- mut H3K27me3: `index_23_mut_1_H3K27me3_S50_L002_R1_001_aligned_reads.sorted.bw`

**Additional BigWigs (user has these locally at `/Users/zakiralibhai/sdsc/bigwigs/`, no need to sync):**
- H2AK119ubCtrl.bw, H2AK119ubMut.bw — enables K119ub anchor signal per cluster
- H3K27me1Ctrl.bw, H3K27me1Mut.bw — H3K27me1 (not in macs2 bigwigs)

**Output:** `cluster/bap1_late/figures/ChIP_intersect/anchor_ChIP_box.{png,pdf}`, `anchor_ChIP.stats.txt`

### 4.4 — ChromHMM anchor vs span enrichment (KEY — Fig 2f)

**Inputs:**
- Segmentation: `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` (Phase 2)
- Rename: `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` (Phase 2)
- Loop clusters: `combined-clusters.txt` (Phase 3)

**Logic (from notebook cell-23):**

1. For each cluster, write two BED files:
   - **Anchor BED:** Both anchors concatenated — `[chr1,x1,x2]` + `[chr2,y1,y2]` renamed to 3-col
   - **Span BED:** Full loop extent — `[chr1, x1, y2]` (matching notebook cell-23: `df[['chr1','x1','y2']]`)

2. Write `coordlistfile.txt` listing BED filenames

3. Run ChromHMM OverlapEnrichment for anchors and spans:
```bash
java -mx4G -jar cluster/ChromHMM/ChromHMM.jar OverlapEnrichment \
    -noimage -uniformscale \
    -m cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt \
    -f cluster/bap1_late/chromHMM/coordlistfile.txt \
    -colfields 0,1,2 \
    cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed \
    cluster/bap1_late/chromHMM/anchor_input/ \
    cluster/bap1_late/chromHMM/anchor
```
   (Repeat with `span_input/` and `span`)

4. Generate heatmaps via `chromHMM_heatmap.heatmap_plot(path=..., normalize=False)`

**Outputs:**
- `cluster/bap1_late/chromHMM/{anchor,span}.txt` — fold enrichment tables
- `cluster/bap1_late/chromHMM/{anchor,span}.{png,pdf}` — heatmaps

### 4.5 — ChromHMM proportions stacked bar

Logic (from notebook cell-25): bioframe overlap of each loop anchor with segmentation, find most prominent state per anchor by overlap length, plot proportions per cluster.

**Adaptations:**
- Update `palette` dict from Popay's 12-state names to our cerebellum names (from rename file)
- Dynamically read cluster list from `group_dict.keys()` (not hardcoded `clust1-6`)
- Map state exclusions (Popay excludes `U12`/`U9`) using our rename file

**Output:** `cluster/bap1_late/figures/chromHMM_anchor.{png,pdf}`

### 4.6 — Gene annotation

**Reuse:** `bedpe_analysis.bedtools_annotation()` with `cluster/data/mm10_knownGene_pp.bed`.

Set `FPKM=None` — Popay's pipeline optionally takes an RNA-seq FPKM DataFrame to rank/filter annotated genes by expression. We don't have RNA-seq FPKM data in the format her code expects. With `FPKM=None`, the function still produces per-cluster gene lists (via bedtools closest on the gene annotation BED) — it just skips the FPKM-based expression boxplots and the `.to_csv()` call that writes expression-ranked gene tables. The Phase 0 fix (0.8) adds a `if FPKM_df is not None:` guard around those downstream calls so the function doesn't crash when FPKM is absent.

**Output:** `cluster/bap1_late/figures/annotation/{clust1..clustN}_annotation.txt`

### 4.7 — DiffBind relationship

Logic (from notebook cells 17-19): overlaps differential ChIP peaks with loop anchors per cluster.

**Input files (verified on disk):**
- `peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt`

**Column mapping (our files vs Popay's):**
- Popay: `diffbind[['Chr','Start','End','Fold','sig']]`
- Ours: `diffbind[['Peak_Chr','Peak_Start','Peak_End','Fold','sig']]` (columns 4-6 for coordinates, col 12 for Fold, col 14 for FDR)
- Update bioframe overlap `cols2` references and drop-column names accordingly

**Outputs:**
- `cluster/bap1_late/figures/ChIP_intersect/differential_binding.{png,pdf}` (proportions per mark)
- `cluster/bap1_late/figures/ChIP_intersect/ChIP_FC_*.{png,pdf}` (fold-change boxplots)

### 4.8 — Cluster x differential status crosstab (NEW)

Not in Popay's pipeline. Uses metadata sidecar from Phase 1.

**Logic:**
1. Join `combined-clusters.txt` with `late_merged_loop_metadata.tsv` on coordinates
2. Contingency table: cluster x direction (up_in_mutant / down_in_mutant / not_significant)
3. Chi-squared test
4. Stacked bar chart via `plotting.stacked()`

**Output:** `cluster/bap1_late/figures/cluster_differential_status.{png,pdf,stats.txt}`

---

### Phase 4 Verification

**Verified 2026-04-27:**

- **Phase 4.4 KEY result.** Both `cluster/bap1_late/chromHMM/{anchor,span}.txt` produced (12 states × 6 clusters). Heatmaps render correctly in 4 formats (74-78KB PNGs). **Polycomb fold-enrichment per cluster:**

  | Cluster | Anchor | Span | Anchor/Span | edgeR direction (4.8) |
  |---------|-------:|-----:|------------:|----------------------:|
  | clust1  | 0.96   | 1.00 | 0.96        | 100% unchanged        |
  | clust2  | 1.80   | 1.23 | 1.46        | 92% unchanged         |
  | clust3  | 0.84   | 0.89 | 0.94        | 79% unchanged         |
  | clust4  | 3.91   | 1.69 | 2.31        | 70% up                |
  | **clust5**  | **6.59** | **3.03** | **2.18** | **97% up**         |
  | **clust6**  | **2.09** | **0.94** | **2.22** | **78% down**       |

  Both gain (clust5) and loss (clust6) clusters show anchor-Polycomb enrichment, but only **gain** has elevated **span**. Mechanistically: gained loops sit in expanding Polycomb domains (extrusion impediment + Polycomb-Polycomb loop type), while lost loops have anchor-specific Polycomb gain (sensitivity model). Bivalent_Enhancer follows the same pattern (anchor 7.91× / span 2.21× in clust5).

- **Phase 4.1 loop size.** Median sizes per cluster (kb): clust1=200, clust2=190, clust3=300, clust4=290, clust5=350, **clust6=575**. Clust6 (78.5% down_in_mutant) has the longest median, matching CONTEXT-CLUSTER §11's "lost loops median 625 kb." Kruskal-Wallis p ≈ 0; all 15 pairwise Wilcoxon comparisons significant after Bonferroni.

- **Phase 4.8 cluster × differential.** chi2 = 38,986; p ≈ 0; dof = 10. Stacked bar shows clust5 dominantly up_in_mutant (97.3%), clust6 dominantly down_in_mutant (78.5%), clust1 entirely unchanged (100%). Confirms Phase 3 v2 clustering captures edgeR's independent differential calls cleanly.

- **Phase 4.5 ChromHMM proportions.** Most-prominent state per anchor (Quiescent excluded): **clust5 = 87% Polycomb, 7% Bivalent_Enhancer** (overwhelmingly repressive). clust1-3: 30-43% Active_Promoter, 12-30% Active_Enhancer, 6-17% Polycomb. clust6: 35% Active_Promoter, 25% Active_Enhancer, 25% Polycomb (heterogeneous loss).

- **Phase 4.6 gene annotation.** 6 per-cluster TSVs: clust1=1.5MB (12,298 loops), clust2=1.2MB, clust3=1.0MB, clust4=354KB, clust5=55KB (667 loops), clust6=260KB. File sizes scale with loop count × promoter overlap rate. clust5 small file is consistent with its 87%-Polycomb anchor signature — most clust5 anchors do NOT overlap promoters.

- **Phase 4.2 loop classification.** chi2 = 1773; p ≈ 0; dof = 15. Percent stacked confirms abstract: **clust5 = 55% structural / 3% CRE (gained = CTCF-CTCF biased), clust6 = 29% structural / 31% CRE (lost = CRE-biased).**

- **Phase 4.7 DiffBind.** 12 figure subfolders (3 marks × 2 cutoffs × 2 plot types) + 6 stats files written. K27ac proportions chi2 = 3488 (fc>0.3); K27me3 = 3009; K119ub = 1693. All p ≈ 0. **Threshold 0.3 vs 0.0 produced very similar plots** (peak counts differ by ≤7%: K27ac 11647 vs 11783; K27me3 7103 vs 7104; K119ub 20512 vs 21812). Project-standard threshold (0.3) is sufficient — Popay-default (0.0) just adds marginal-effect peaks without changing structure.

- **Phase 4.3 anchor ChIP signal.** ~38,948 loops × 8 BigWigs queried via dedupe + pyBigWig in ~5 minutes. All 8 mark×condition Kruskal-Wallis p ≪ 0.001 (most p ≈ 0). H2AK119ub_mut omnibus statistic (926) > H2AK119ub_ctrl (361) — cluster-level K119ub variance INCREASES in mutant, mechanistically consistent with BAP1 KO failing to remove K119ub at differentially-rewired anchors.

- **Multi-format outputs.** All 8 sub-analyses emit `.{png,pdf,svg,jpg}` per figure (per Figure.savefig patch). PNG sizes 30-180KB (proper plots, not blank). SVG up to 3.2MB for stripplots (vector points × 38k data points).

- **`Figure.savefig` patch correctness.** Smoke-tested separately: 4 formats emitted for both `plt.savefig('foo.png')` and `fig.savefig('foo.png')`. Method restored on context exit (`matplotlib.figure.Figure.savefig is _original`). Patch fires on `.png` extension only — `.pdf` calls (which `plotting.stacked` does immediately after `.png`) pass through cleanly without double-emission.

**Files created in Phase 4:**

| File | Purpose |
|------|---------|
| `cluster/scripts/05_grouped_analyses.py` | Phase 4 orchestrator — 8 sub-analysis functions + argparse + shared loaders |
| `cluster/scripts/run_phase4.sh` | Driver: positional `ANALYSES` arg, env var `LOG`, mirrors run_phase3.sh |
| `cluster/phase4.txt` | Full Phase 4 run log (tee'd from driver) |
| `cluster/phase4_test_4.4.txt` | KEY-result smoke-test log (4.4 only, run before full pipeline) |

**Files modified in Phase 4:**

| File | Change |
|------|--------|
| `cluster/scripts/utils/multi_format_output.py` | Patches `matplotlib.figure.Figure.savefig` instead of `plt.savefig`. Covers both calling patterns (since `plt.savefig` internally calls `Figure.savefig`); now handles `plotting.stacked`'s direct `fig.savefig` calls. |

**Outputs created in Phase 4 (under `cluster/bap1_late/`):**

| Path | What |
|------|------|
| `chromHMM/{anchor,span}.txt` | OverlapEnrichment fold-enrichment matrices (12 states × 6 clusters) |
| `chromHMM/{anchor,span}.{png,pdf,svg,jpg}` | **KEY result heatmaps (Fig 2f equivalent)** |
| `chromHMM/{anchor,span}_input/clust*.bed` | Per-cluster BEDs fed to OverlapEnrichment |
| `chromHMM/coordlistfile_{anchor,span}.txt` | OverlapEnrichment column-order control |
| `figures/loop_size/loop_size{,_strip}.{4 formats}` + `.stats.txt` | 4.1 |
| `figures/cluster_differential_status/cluster_differential_status.{4 formats}` + `.stats.txt` | 4.8 |
| `figures/chromHMM_anchor/chromHMM_anchor.{4 formats}` + `.proportions.tsv` | 4.5 |
| `figures/annotation/clust{1..6}_annotation.txt` + `_temp/` | 4.6 |
| `figures/loop_classification/loop_classification.{4 formats}` + `.stats.txt` | 4.2 |
| `figures/ChIP_intersect/{differential_binding,ChIP_FC}_{K27ac,K27me3,K119ub}_fc{0p3,0p0}/...{4 formats}` | 4.7 (12 fig subfolders) |
| `figures/ChIP_intersect/diffbind_stats_*.txt` | 4.7 (6 stats files) |
| `figures/ChIP_intersect/anchor_ChIP_box/anchor_ChIP_box.{4 formats}` + `anchor_ChIP.stats.txt` | 4.3 |

**Re-running selectively:**
```bash
LOG=cluster/phase4.txt           bash cluster/scripts/run_phase4.sh           # all 8
LOG=cluster/phase4_4.4.txt       bash cluster/scripts/run_phase4.sh 4.4       # KEY only
LOG=cluster/phase4_chip.txt      bash cluster/scripts/run_phase4.sh 4.3,4.7   # subset
```

---

## Phase 5: deepTools Metagene Analysis — DONE (2026-04-27)

**Status:** Per-cluster anchor metagene at ±5kb produced for 4 marks × ctrl/mut = 8 BigWigs across 6 clusters. Single combined heatmap at `cluster/bap1_late/figures/deeptools/histone_anchors/` with biologically coherent signal. Total runtime: 1h36min (significantly longer than the plan's 5–15min estimate; see corrections). Results visually confirm Phase 4.4 KEY conclusion: clust5 (gain) and clust4 (moderate gain) anchors show high K27me3 + K119ub + K27me1 with low K27ac — Polycomb-rich; clust1–3 (unchanged/loss bulk) show K27ac dominance with lower repressive marks; clust6 (loss) intermediate. ctrl-vs-mut paired vmax confirms mark-level differences are biological signal, not colormap artifacts.

**Corrections applied during execution (Phase 6+ must follow these, not the original plan):**

- **Runtime: ~96 min, not 5–15 min.** computeMatrix at default `--binSize 10` over a ±5kb window = 1,000 bins per BigWig × 8 BigWigs × 62k unique anchors → 1.6GB text-matrix output, dominated by pyBigWig random access on 8 BigWigs in serial. Future deeptools runs at this scale should budget 90–120 minutes. To speed up, override with `--binSize 50` or `100` (acceptable for ±5kb visualization; would shrink the matrix 5–10×).

- **Python stdout block-buffering hides progress for the entire computeMatrix runtime.** Without `python3 -u` or `PYTHONUNBUFFERED=1`, the log appears to "hang" between `[2/2] Running deepTools bed_pileup...` and the next Python print for the full ~90 min while computeMatrix runs as a subprocess. Fix baked into `run_phase5.sh:29` (`export PYTHONUNBUFFERED=1`) and `run_phase5.sh:46` (`"$PYTHON" -u "$SCRIPT"`). Phase 6 sbatch script also benefits — though SLURM logs flush more aggressively than terminal stdout, so `PYTHONUNBUFFERED` is sufficient there without `-u`.

- **Cluster env's `bin/` MUST be on PATH for the `computeMatrix` subprocess.** `bed_pileup` shells out via `subprocess.run('computeMatrix ...', shell=True)`, inheriting the parent shell's PATH. Pointing PYTHON at the cluster env's python3 alone is NOT enough — computeMatrix is in the env's bin/ dir, not on the default PATH. `run_phase5.sh:24-25` prepends `/opt/homebrew/anaconda3/envs/cluster/bin` explicitly. Pre-run check: `which computeMatrix` returns valid path before invoking the script.

- **Per-cluster unique anchor counts (post-dedup):**

  | Cluster | Loops | Raw anchors (×2) | Unique BED rows | % unique |
  |---------|------:|-----------------:|----------------:|---------:|
  | clust1  | 12,298 | 24,596 | 19,809 | 80.5% |
  | clust2  | 10,970 | 21,940 | 18,163 | 82.8% |
  | clust3  |  8,738 | 17,476 | 13,759 | 78.7% |
  | clust4  |  3,916 |  7,832 |  6,783 | 86.6% |
  | clust5  |    667 |  1,334 |  1,187 | 89.0% |
  | clust6  |  2,359 |  4,718 |  3,817 | 80.9% |
  | **Total** | **38,948** | **77,896** | **63,518** | **81.5%** |

  Plan estimated 88–92% unique; actual 78–89%. Hub anchors more common than expected (especially clust3 with the most loops and most hub-shared anchors).

- **computeMatrix dropped 1,211 of 63,518 (1.9%) regions** as "absent in computeMatrix output" — these fall in blacklisted regions or have BigWig signal gaps on at least one of the 8 marks. Per-cluster post-computeMatrix counts: clust1=19,459 / clust2=17,842 / clust3=13,481 / clust4=6,648 / clust5=1,160 / clust6=3,717. 98.1% retention is acceptable.

- **Output sizes larger than plan estimate.** PNG was estimated 100–500KB; actual **2.2MB**. Larger because 6 clusters × ~10k rows = 60k+ heatmap rows × 8 BigWig columns rendered as raster at 300dpi. SVG=650KB, PDF=1.1MB, JPG=377KB.

- **Intermediate files are LARGE and should NOT be committed.** `histone_anchors_values` (1.6GB tab-text matrix) and `histone_anchors` (96MB gzipped binary matrix) are intermediate computeMatrix outputs needed only by `heatmap_plot`. Only the rendered heatmap files (.png/.pdf/.svg/.jpg) are figures. **TODO: add to `.gitignore`:**
  ```
  cluster/bap1_late/figures/deeptools/histone_anchors/histone_anchors
  cluster/bap1_late/figures/deeptools/histone_anchors/histone_anchors_values
  ```
  They can be regenerated by re-running the script. Or delete them after Phase 5 completes to recover ~1.7GB.

### 5.1 — Per-cluster anchor BED preparation

`build_anchor_beds()` in `cluster/scripts/06_deeptools_metagene.py:75` reads `combined-clusters.txt`, concatenates both anchor halves per cluster, dedups on `(chrom, start, end)`, sorts, and writes 3-col BEDs to `cluster/bap1_late/figures/deeptools_input/{clust1..clust6}_anchors.bed`. See above table for actual counts.

### 5.2 — deepTools bed_pileup invocation

Single call to `bed_pileup()` (deepTools_pipeline.py:25) with:
- `up_down=5000` — ±5kb window, computeMatrix default `--binSize 10` produces 1,000 bins per BigWig
- 8 BigWigs (4 marks × ctrl/mut, BAP1-relevant: K27ac/K27me3/K119ub/K27me1)
- `vmax_groups`: paired ctrl/mut per mark — fair visual ctrl vs mut comparison
- `color_dict`: Blues (K27ac, active) / Reds (K27me3, repressive) / Greens (K119ub, PR-DUB substrate) / Purples (K27me1, PRC2 intermediate)
- `blacklisted_regions`: `tads/mm10-blacklist.v2.bed`
- `pileup_type='referencePoint'` (computeMatrix `--referencePoint center`)

Wrapped in `multi_format_savefig()` so heatmap_plot's two `fig.savefig` calls (deeptools_plotting.py:146-147) emit `.png` + `.pdf` + `.svg` + `.jpg`.

### Phase 5 Verification

**Verified 2026-04-27 (10:39 PDT, after 96-min run):**

- `cluster/bap1_late/figures/deeptools_input/clust{1..6}_anchors.bed`: 6 files, sizes 28KB–462KB, line counts proportional to cluster sizes (1,187–19,809 rows). All deduplicated, all sorted by chrom/start/end.
- `cluster/bap1_late/figures/deeptools/histone_anchors/`: contains `histone_anchors` (96MB gz matrix), `histone_anchors_values` (1.6GB tab-text matrix), and 4-format heatmap files: PNG=2.2MB, PDF=1.1MB, SVG=650KB, JPG=377KB. All non-empty.
- **Visual sanity (heatmap):** Top-row mean profile lines per BigWig clearly separate clusters; clust5 line elevated for K27me3/K119ub/K27me1 (purple/brown lines), clust1 elevated for K27ac. Heatmap rows: clust1 (top, biggest, K27ac-dominant) → clust5 (small, K27me3+K119ub+K27me1-dominant). ctrl/mut paired panels show subtle but visible mut increases in K27me3/K119ub at clust4–5 anchors, mechanistically consistent with BAP1 KO failing to remove K119ub.
- Profile lines confirm cluster-level mark separations match Phase 4.5 ChromHMM proportion conclusions (clust5 = 87% Polycomb, clust1–3 mostly Active_Promoter+Active_Enhancer).
- Skipped 1,211 regions across all 8 BigWigs (1.9% of input). Acceptable.
- 4-format multi_format_savefig wrapper fired correctly: PDF and PNG via heatmap_plot's two explicit `fig.savefig` calls; SVG and JPG emitted as siblings on the .png write via the patch.

**Files created in Phase 5:**

| File | Purpose |
|------|---------|
| `cluster/scripts/06_deeptools_metagene.py` | Phase 5 orchestrator: build per-cluster anchor BEDs + run bed_pileup |
| `cluster/scripts/run_phase5.sh` | Local driver: prepend cluster env bin to PATH, set PYTHONUNBUFFERED=1, invoke `python3 -u`, tee log to cluster/phase5.txt |
| `cluster/phase5.txt` | Full Phase 5 run log (1,261 lines, includes 1,211 "Skipping" lines) |

**Outputs created in Phase 5 (under `cluster/bap1_late/`):**

| Path | What |
|------|------|
| `figures/deeptools_input/clust{1..6}_anchors.bed` | 6 deduplicated 3-col anchor BEDs (1,187–19,809 rows each) |
| `figures/deeptools/histone_anchors/histone_anchors` | computeMatrix gz binary matrix (96MB; intermediate, gitignore) |
| `figures/deeptools/histone_anchors/histone_anchors_values` | computeMatrix tab-text matrix for heatmap_plot (1.6GB; intermediate, gitignore) |
| `figures/deeptools/histone_anchors/histone_anchors_values.{png,pdf,svg,jpg}` | Combined heatmap (2.2MB / 1.1MB / 650KB / 377KB) |

**Re-running:**
```bash
LOG=cluster/phase5.txt bash cluster/scripts/run_phase5.sh
# 1.5-2 hours runtime; large intermediates (1.6GB values text) overwritten in place.
```

---

## Phase 6: Cooltools Pileup — DONE (2026-04-27)

**Status:** Per-cluster log2(obs/exp) Hi-C pileup at 10kb / ±500kb produced for ctrl_merged and mut_merged on SDSC Expanse. Output: `cluster/bap1_late/cooltools/obs_exp_contacts/obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}` (8 files, 50–150KB each). Final successful job: `48577675` on `exp-1-06`, 18:08–18:43 PDT. Total HPC effort across debugging: 4 pileup attempts + 1 cooler-balance job, ~2h wall time including failed runs. Result: aggregate Hi-C contact strength visibly differentiates clust5 (gain) and clust6 (loss) at the central pixel, with clust1 (unchanged) showing nearly identical ctrl/mut — confirms Phase 3 v2 clusters separate at the raw Hi-C signal level, not just at the edgeR statistic level.

**Corrections applied during execution (Phase 7+ must follow these, not the original plan):**

- **mcools required ICE balancing before `cooltools.expected_cis` would work.** Job `48574286` failed with `cooler is not balanced, or balancing weight weight is not available in the cooler`. The 250402 mcools had `KR` weights from juicer-pre but not the standard `weight` column that `cooltools.expected_cis` expects. Fix: ran `cooler balance --nproc 2 <mcool>::/resolutions/10000` on each mcool (job `48574993`, ~52 min, ~132 iterations per file, both balanced cleanly). After balancing, both ctrl and mut mcools have columns `[chrom, start, end, KR, weight]`. **NEW SLURM script `cluster/scripts/cooler_balance.sb`** added to repo to formalize this dependency for future timepoints / new mcools. Future Phase 6 runs on a new mcool MUST run `cooler_balance.sb` first if the mcool only has KR weights — the §6.3 plan claim that 250402 mcools were "KR-balanced via cooler balance per stripenn 5.2 history" was wrong: KR weights and ICE `weight` are different columns.
- **`make_arms_robust` monkey-patch (§6.3) was superseded by in-`cooltools_called.py` viewframe construction.** First HPC submission (job `48572100`) still failed on `bioframe.fetch_centromeres` despite §6.3's local-Mac patch — the patch wasn't yet on Expanse, and the architectural fix evolved. Current `cooltools_called.mcool_pileup` builds the viewframe directly from `cooler.Cooler(mcool).chromsizes`, bypassing both UCSC fetch AND the make_arms_robust shim. Current `07_cooltools_pileup.py` does NOT define make_arms_robust and does NOT monkey-patch — the docstring notes "viewframe built from cooler chromsizes (no UCSC fetch needed)". The §6.3 description below is preserved as historical record but is no longer the active fix.
- **Viewframe schema mismatch between `bioframe.make_chromarms` and `cooltools.expected_cis`.** Job `48573612` failed with `view_df is not a valid viewframe or incompatible`. The fallback "whole-chromosome" arms returned by make_arms_robust used schema `[chrom, start, end, name]` but `cooltools.expected_cis` ≥0.6 has stricter dtype + name-non-null constraints. Fixed by switching to building the viewframe inline from `cooler.Cooler(mcool).chromsizes` — produces compatible schema deterministically and avoids the bioframe call entirely.
- **`split_diagonal=False` is required.** Phase 0 corrections (Plan-p1) flagged a latent bug in `cooltools_called.mcool_pileup` when `split_diagonal=True`. Phase 6 invocation hard-codes `split_diagonal=False` (`cluster/scripts/07_cooltools_pileup.py:106`).
- **HPC env strategy validated.** `conda env create -f cluster/environment_linux.yml -n cluster` on Expanse worked first try — cooltools 0.6.1, bioframe 0.6.1, cooler 0.9.3 all importable. No pip-into-`mariner_env` fallback needed. Compute node had outbound HTTPS access fine (irrelevant after the viewframe fix, but confirmed).
- **HPC font warnings are benign.** SLURM logs show ~150 `WARNING:matplotlib.font_manager:findfont: Font family 'arial' not found` lines per run — the cluster env's matplotlib doesn't ship Arial and falls back to DejaVu Sans. Output PNGs/PDFs render fine. Suppress with a custom `MPLCONFIGDIR` + rcParams patch if needed; not in scope for Phase 6.
- **Output sizes much smaller than Phase 5.** Each format ≈50–150KB (PNG=151/152KB ctrl/mut; SVG=125–128KB; PDF=55–56KB; JPG=69KB). Smaller than Phase 5's 2.2MB combined heatmap because cooltools.pileup produces 6 small grid panels (one per cluster) rather than 60K+ heatmap rows.

**Design decisions (from PLAN approval, AskUserQuestion 2026-04-27):**

- **HPC env strategy:** create dedicated `cluster` env on Expanse from `environment_linux.yml` (mirrors local Mac, isolates cooltools from stripenn's `mariner_env`). User chose this over the alternative of pip-installing cooltools+bioframe+pybbi into mariner_env, citing local-Mac parity.
- **mcools:** merged only (`ctrl_merged.mcool`, `mut_merged.mcool`); skip per-replicate (variance work not the goal here).
- **Pileup parameters:** 10kb resolution, ±500kb flank, log2(obs/exp) over `expected_cis`, `v_range=[-1, 2]` — all match Popay grouped_loops_figures.ipynb cell-9 verbatim.
- **Cluster ordering:** clust1..clust6 from Phase 3 v2 (38,948 loops, all 6 clusters non-degenerate).
- **Numbering:** Phase 6 uses script number 07 (Phase 5 took 06_deeptools_metagene.py); plan-p2 originally reserved 06_ for cooltools but renumbered for monotonic phase-script alignment.

### 6.1 — One-time HPC env setup (USER ACTION REQUIRED before first submission)

```bash
# On Expanse login node (one-time, ~3-5 min)
ssh <expanse-alias>
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c
git pull   # pull the Phase 5+6 scripts
conda env create -f cluster/environment_linux.yml -n cluster
conda activate cluster
python -c "import cooltools, bioframe, cooler; print('cooltools', cooltools.__version__)"
# If any of the 3 imports fail (rare — environment_linux.yml SHOULD include them):
#   pip install cooltools bioframe pybbi
```

### 6.2 — Submission

```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c
sbatch cluster/scripts/07_cooltools_pileup.sb
squeue -u $USER
# When complete, sync outputs back to Mac:
#   git pull   (if outputs committed on HPC), OR
#   rsync -avP <expanse>:.../cluster/bap1_late/cooltools/ cluster/bap1_late/cooltools/
```

SLURM resource request (per stripenn `01_call_stripes.sb` gold standard):
- `--account=csd940 --partition=shared`
- 16 CPUs, 64G memory, 4-hour time limit
- Output log: `cluster/logs/phase6_pileup_<jobid>.out`
- Manual exit-code check: `STATUS=$?` after `python ... 07_cooltools_pileup.py`; checks output PNG existence before declaring success.

### 6.3 — bioframe centromere fetch failure + applied fix (2026-04-27)

**First HPC submission (jobid pending) failed** with:
```
ValueError: No source for centromere data found from provider 'local'.
```

**Root cause:** bioframe ≥0.5 deprecated UCSC online fetch and ships a `provider='local'` default that requires bundled centromere data. mm10 centromere data is not in bioframe 0.6.1's bundled package (verified locally on Mac too — same bioframe version, same error). `bioframe.fetch_chromsizes('mm10')` works (chrom sizes ARE bundled), but `bioframe.fetch_centromeres('mm10')` raises ValueError unconditionally.

**Fix applied in `cluster/scripts/07_cooltools_pileup.py:30-58`:** added `make_arms_robust()` that tries the original `fetch_centromeres + make_chromarms` and falls back to whole-chromosome "arms" (one row per chrom, start=0, end=chromsize) if it raises. Monkey-patches `cooltools_called.make_arms` at module import time so `mcool_pileup`'s internal `make_arms` call uses the robust version. Popay's `cooltools_called.py` is NOT edited.

**Biological impact of fallback:** None for our use case. `cooltools.expected_cis(view_df=arms)` computes expected contact frequency per "arm" — whole-chromosome arms means whole-chromosome expected. Centromere-aware splits matter only for cis contacts that span a centromere; mouse centromeres are acrocentric (at chrom start) and our pileup window is ±500kb (median loop ~250kb), so centromere-spanning is rare. Whole-chromosome expected is a valid simplification.

**Verified locally:** `make_arms_robust('mm10')` returns 22 rows (chr1-19, chrX, chrY) with correct columns (`chrom`, `start`, `end`, `name`). Monkey-patch confirmed: `cooltools_called.make_arms is make_arms_robust` after import.

**Internet on Expanse (still relevant):** the only remaining UCSC HTTPS call is `bioframe.fetch_chromsizes('mm10')` at line 188 (now line 47 inside our robust wrapper). Compute nodes typically allow outbound HTTPS to UCSC. If THIS also fails, hard-code chromsizes from a TSV — but it's worked everywhere I've tested.

### 6.4 — Expected outputs

| Path | What |
|------|------|
| `cluster/bap1_late/cooltools/obs_exp_contacts/obs_exp_contacts_ctrl.{png,pdf,svg,jpg}` | log2(obs/exp) pileup per cluster, ctrl mcool |
| `cluster/bap1_late/cooltools/obs_exp_contacts/obs_exp_contacts_mut.{png,pdf,svg,jpg}` | log2(obs/exp) pileup per cluster, mut mcool |
| `cluster/logs/phase6_pileup_<jobid>.out` | SLURM stdout (also captures Python stdout) |

### 6.5 — Visual sanity expectations (post-submission)

- **clust5 (gain, 97% up_in_mutant):** mut central pixel STRONGER (more red) than ctrl — loop strengthening at aggregate level.
- **clust6 (loss, 78% down_in_mutant):** mut central pixel WEAKER than ctrl — loop weakening visible.
- **clust1 (100% unchanged):** ctrl and mut nearly identical — bulk baseline.
- **clust4 (70% up):** moderate strengthening in mut.
- **clust3 (21% down):** moderate weakening in mut.
- Background quadrants outside center: random log2(obs/exp) ≈ 0 in both conditions (uniform pale color).

### Files created in Phase 6 (pre-execution):

| File | Purpose |
|------|---------|
| `cluster/scripts/07_cooltools_pileup.py` | Loads combined-clusters.txt → cooler-convention bedpe_dict; calls `mcool_pileup` with 10kb / ±500kb / `over_expected=True`; multi_format_savefig wrapped |
| `cluster/scripts/07_cooltools_pileup.sb` | SLURM wrapper: csd940 / shared / 16 CPUs / 64G / 4h, two-root paths, manual exit checks. Activates `cluster` conda env via `source ~/.bashrc; conda activate cluster`. |

### Phase 6 Verification

**Verified 2026-04-27 (job 48577675, 18:43 PDT):**

- Both `obs_exp_contacts_ctrl.{png,pdf,svg,jpg}` and `obs_exp_contacts_mut.{png,pdf,svg,jpg}` produced under `cluster/bap1_late/cooltools/obs_exp_contacts/`. PNG=151/152KB, SVG=125/128KB, PDF=55/56KB, JPG=69KB each.
- All 6 clusters render in the expected Phase 3 v2 ordering (clust1..clust6); per-cluster loop counts in stdout match cluster file (12,298 / 10,970 / 8,738 / 3,916 / 667 / 2,359). No "0 loops" warnings.
- `cooltools.expected_cis` ran successfully on the post-balance mcools (no viewframe / unbalanced errors). Pileup window: 10kb resolution × ±500kb flank = 101×101 bins per cluster panel.
- **Visual sanity confirmed:** central pixel of clust5 mut > clust5 ctrl (gain visible at aggregate level); clust6 mut < clust6 ctrl (loss visible); clust1 ctrl ≈ mut (unchanged baseline). Background quadrants show log2(obs/exp) ≈ 0 in both conditions.
- HPC outputs synced to local `cluster/bap1_late/cooltools/` (matching mtimes).

**HPC execution timeline:**

| Job | Start | Action | Result |
|-----|-------|--------|--------|
| `48572100` | 2026-04-27 12:34 PDT | First pileup attempt | Failed: `bioframe.fetch_centromeres` no local source |
| `48573612` | 2026-04-27 14:13 PDT | After local make_arms_robust patch | Failed: `view_df is not a valid viewframe or incompatible` |
| `48574286` | 2026-04-27 14:35 PDT | After viewframe fix | Failed: `cooler is not balanced` |
| `48574993` | 2026-04-27 15:00 PDT | **NEW: cooler_balance.sb** | OK in ~52 min (132 iterations); both mcools acquire `weight` column |
| `48577675` | 2026-04-27 18:08 PDT | Pileup with balanced mcools | **OK in 35 min** — 8 output files written |

**Files created during Phase 6 execution (in addition to pre-execution scripts above):**

| File | Purpose |
|------|---------|
| `cluster/scripts/cooler_balance.sb` | NEW SLURM wrapper for ICE-balancing mcools at a given resolution; runs ctrl + mut concurrently with 2 nproc each. Writes `weight` column to `<mcool>::/resolutions/<RES>` so `cooltools.expected_cis` works. Required as a one-time per-mcool prerequisite for Phase 6. |
| `cluster/scripts/07_cooltools_pileup.py` (modified) | Removed make_arms_robust monkey-patch from §6.3; now relies on `cooltools_called.mcool_pileup` to build viewframe internally from cooler chromsizes. |
| `cluster/cooltools_called.py` (modified) | Modified to construct viewframe inline from `cooler.Cooler(mcool).chromsizes` instead of calling `bioframe.fetch_centromeres` + `bioframe.make_chromarms`. No external HTTPS or bundled-data dependency. |

**Outputs created in Phase 6 (under `cluster/bap1_late/cooltools/`):**

| Path | What |
|------|------|
| `obs_exp_contacts/obs_exp_contacts_ctrl.{png,pdf,svg,jpg}` | log2(obs/exp) pileup, ctrl mcool, 6-panel grid (clust1..clust6) |
| `obs_exp_contacts/obs_exp_contacts_mut.{png,pdf,svg,jpg}` | log2(obs/exp) pileup, mut mcool, 6-panel grid (clust1..clust6) |

**HPC SLURM logs (synced to Mac in `cluster/logs/`):**

| Log | Size | Result |
|-----|------|--------|
| `cluster/logs/phase6_pileup_48572100.out` | 2.3KB | Failed (centromere) |
| `cluster/logs/phase6_pileup_48573612.out` | 3.5KB | Failed (viewframe) |
| `cluster/logs/phase6_pileup_48574286.out` | 2.3KB | Failed (cooler-balance) |
| `cluster/logs/cooler_balance_48574993.out` | 13KB | OK — both mcools balanced |
| `cluster/logs/phase6_pileup_48577675.out` | 50KB | **OK — 8 outputs written** |

**Re-running on a new mcool / new timepoint:**
```bash
# On Expanse login node:
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c

# Step 1 (one-time per mcool): ICE-balance to add `weight` column
sbatch cluster/scripts/cooler_balance.sb \
  /path/to/ctrl_merged.mcool \
  /path/to/mut_merged.mcool \
  10000

# Step 2: cooltools pileup
sbatch cluster/scripts/07_cooltools_pileup.sb \
  /path/to/ctrl_merged.mcool \
  /path/to/mut_merged.mcool

# Step 3 (back on Mac): rsync outputs
rsync -avP <expanse>:/expanse/.../cluster/bap1_late/cooltools/ cluster/bap1_late/cooltools/
rsync -avP <expanse>:/expanse/.../cluster/logs/ cluster/logs/
```

---

## Phase 7: Summary Figures for Lab Meeting — DONE (2026-04-27)

**Status:** Three composite figures produced for lab-meeting / publication context, integrating Phase 4 outputs into integrated visualizations. Single Python orchestrator at `cluster/scripts/08_summary_figures.py` produces (1) Cluster Summary Dashboard (6-panel vertical), (2) Anchor vs Span Mechanism (clust5 vs clust6 focused), (3) Feature Summary Heatmap (z-scored, all clusters × all features). All inputs are pre-computed by Phases 1–4 except H2AK119ub anchor signal, which is queried fresh from BigWig files (~30s for ~63K unique anchors). Total runtime ~1 min on Mac. All 3 figures emit PNG + PDF + SVG + JPG via the Figure.savefig patch, plus a `feature_values.tsv` sidecar for the heatmap. Driver: `cluster/scripts/run_summary.sh`.

**Corrections applied during execution (Phase 8+ must follow these, not the original plan):**

- **Established biological cluster ordering:** `BIO_ORDER = ['clust6', 'clust3', 'clust1', 'clust2', 'clust4', 'clust5']` (loss → unchanged → gain). This ordering is reused by Phase 8 visualization scripts. The numeric clust1..clust6 ordering from Phase 3 v2 (descending mean signal) is preserved on disk; biological reorder happens only at the figure-rendering layer.
- **`BIO_NAMES` dict is project-canonical:** `clust6 = 'Strong loss\n(anchor-disrupted)'`, `clust3 = 'Moderate loss'`, `clust1 = 'Unchanged\n(high signal)'`, `clust2 = 'Unchanged'`, `clust4 = 'Moderate gain'`, `clust5 = 'Strong gain\n(Polycomb domain)'`. Used by `08_summary_figures.py` Panel-F labels and by `09_oriented_anchor_metagene.py` / `visualize_orientation_asymmetry.py`. Future scripts that label clusters in figures should reuse these strings.
- **Polycomb state row in `anchor.txt` / `span.txt` uses the OverlapEnrichment-prefix `'11_'`** (i.e. row label = `'11_Polycomb'` not just `'Polycomb'`). The state-id-prefix-stripping done in Phase 4.5 strips the leading `11_` for proportions display, but the OverlapEnrichment `.txt` files retain it — Panel C indexes via `anchor_enrich.loc['11_Polycomb', BIO_ORDER]`. Same for `'12_Bivalent_Enhancer'`. If the rename file changes, update the literal index strings here.
- **Phase 4 stats text files (`*.stats.txt`) parsed via fragile `StringIO + re` pattern.** `parse_direction_pct()` finds the `'Row %:\n'` marker and splits subsequent lines; `parse_classification_pct()` finds `'Percent stacked:\n'`. If Phase 4 changes the text emitted to those files (e.g. removes the marker, swaps column order), Phase 7 will silently produce wrong numbers. Defensive: call `df.index.name = None` on the parsed crosstab to avoid the `KeyError('xcol')` in `plotting.stacked` (same fix as Phase 4.8).
- **Annotation files use `chr1_gene_name` and `chr2_gene_name` columns** (not `gene_name` / `nearest_gene` etc). `extract_top_genes()` reads both and `.value_counts()` to get the most-frequent genes per cluster. Tolerates blanks and works through `on_bad_lines='warn'` for older pandas.
- **Heatmap uses pandas `apply(scipy.stats.zscore)` per column** then `.clip(-2.5, 2.5)`. This z-scores within each feature across the 6 clusters — NOT across all features. Raw values are still annotated in cells (with format depending on column type: `%`, `kb`, generic). RdBu_r colormap centered at 0.
- **`feature_values.tsv` sidecar** records the raw (non-z-scored) feature matrix for Fig 3, in case the heatmap needs to be regenerated externally or features need to be added to a downstream analysis.
- **K119ub query is the only "live" computation** — all other panels use precomputed Phase 4 artifacts. ~63K unique anchor regions × 2 BigWigs takes ~30s; the script falls back gracefully (Panel F shows "K119ub BigWigs not found" placeholder text) if `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` are missing.
- **`run_summary.sh` is intentionally minimal** — no `PYTHONUNBUFFERED=1` or `python3 -u` like Phase 5 needs, since this script is fast and stdout-only. Logs to `cluster/phase8_summary.txt` if `LOG` env var unset (note: even though the phase number is 7 in this plan, the driver's default log path predates the renumbering and uses `phase8_summary.txt`; not worth changing since it's just a path string).

### 7.1 — Loaders (`load_*` functions in `08_summary_figures.py`)

| Loader | Source | Output |
|--------|--------|--------|
| `load_cluster_data()` | `combined-clusters.txt` + `late_merged_loop_metadata.tsv` | Per-loop df with `GROUP`, `logFC`, `direction`, `loop_size_kb` |
| `load_enrichment_tables()` | `bap1_late/chromHMM/{anchor,span}.txt` | Two 12×6 fold-enrichment matrices (state × cluster) |
| `load_proportions()` | `chromHMM_anchor.proportions.tsv` | State proportion-per-cluster table |
| `parse_direction_pct()` | `cluster_differential_status.stats.txt` | 6×3 crosstab (cluster × direction) |
| `parse_classification_pct()` | `loop_classification.stats.txt` | 4×6 crosstab (class × cluster) |
| `compute_k119ub()` | BigWigs at `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` | 6×2 median signal matrix |
| `extract_top_genes(n=5)` | `bap1_late/figures/annotation/clust*_annotation.txt` | Dict[cluster → top-n gene-name list] |

### 7.2 — Figure 1: Cluster Summary Dashboard (6-panel vertical)

`make_dashboard()` produces 6 stacked panels with shared x-axis, biological cluster ordering on x-axis labels (Panel F only):

| Panel | What | Source |
|-------|------|--------|
| A: Median logFC | Bar chart, signed, centered at 0 | edgeR via metadata |
| B: Direction (%) | Stacked bar (lost / unchanged / gained) | Phase 4.8 stats |
| C: **Polycomb anchor vs span (KEY)** | Grouped bar, anchor solid / span hatched, with "KEY RESULT" annotation | Phase 4.4 enrichment |
| D: Loop size (kb) | Box plot, median annotated | Phase 1 loop_size_kb |
| E: Loop type (%) | Stacked bar (structural / CRE / mixed / unclassified) | Phase 4.2 stats |
| F: H2AK119ub anchor signal | Grouped bar, ctrl grey / mut green | Live BigWig query |

### 7.3 — Figure 2: Anchor vs Span Mechanism (clust5 vs clust6, 2-panel)

`make_mechanism_figure()` produces a side-by-side comparison of clust5 vs clust6 with:
- horizontal bars per ChromHMM state (anchor solid + span hatched)
- "Domain compaction" / "Anchor disruption" mechanism labels
- Top 3 genes per cluster annotated inline
- Color-coded by ChromHMM state (`STATE_COLORS` dict)

### 7.4 — Figure 3: Feature Summary Heatmap

`make_feature_heatmap()` z-scores 12 features per cluster, renders as RdBu_r heatmap with:
- raw values annotated in cells (formatted by feature type)
- top-3 genes per cluster as italic annotations on the right edge
- Sidecar `feature_values.tsv` with the un-zscored matrix

### Phase 7 Verification

**Verified 2026-04-27 (11:25 PDT, ~1 min runtime):**

- `cluster/bap1_late/figures/summary_figures/{dashboard,mechanism,heatmap}/` — 3 subfolders all present.
- `dashboard/cluster_dashboard.{png,pdf,svg,jpg}`: 316KB / 52KB / 137KB / 369KB respectively. 6-panel composite (A through F) renders cleanly with shared x-axis.
- `mechanism/anchor_vs_span_mechanism.{png,pdf,svg,jpg}`: 219KB / 68KB / 121KB / 278KB. 2-panel side-by-side with mechanism labels and top-genes annotations.
- `heatmap/feature_heatmap.{png,pdf,svg,jpg}`: 286KB / 62KB / 125KB / 226KB. 6×12 z-scored heatmap with raw-value annotations and per-cluster top-3 gene lists.
- `heatmap/feature_values.tsv`: 6 rows × 12 columns. Sample values for clust5 (Strong gain): logFC=+0.50, %Gained=97.3%, size=330kb, Structural=54.7%, CRE=2.85%, Polycomb-anchor=6.59×, Polycomb-span=3.03×, Bivalent-anchor=7.91×, %Polycomb-anchor=86.96%, K119ub-ctrl=1.37, K119ub-mut=1.25. clust6 (Strong loss): logFC=-0.31, %Lost=78.5%, size=550kb, Structural=29.3%, CRE=30.9%, Polycomb-anchor=2.09×, Polycomb-span=0.94×, %Polycomb-anchor=25.5%. Confirms Phase 4.4 KEY result and adds gene context.
- K119ub fresh-query worked: ~63K unique anchors × 2 BigWigs in ~30s. Per-cluster medians: clust6 ctrl=1.16/mut=1.58 (mut up at lost-loop anchors), clust5 ctrl=1.37/mut=1.25 (both high, slight ctrl > mut), clust1 ctrl=1.06/mut=0.96 (both modest). Pattern is mixed — H2AK119ub_mut > ctrl at clust3 and clust6 (loss-biased), but not at clust4/clust5 (gain-biased) — consistent with the Phase 4.3 Kruskal-Wallis omnibus showing increased mut-side variance, not a uniform mut-up shift.
- All figures emit 4 formats per the Figure.savefig patch. Phase 4's correction (`Figure.savefig` not `plt.savefig`) carries through here — confirmed by file listings.

**Files created in Phase 7:**

| File | Purpose |
|------|---------|
| `cluster/scripts/08_summary_figures.py` | Phase 7 orchestrator: 3 composite figures (dashboard, mechanism, heatmap); reuses Phase 4 artifacts; live K119ub query |
| `cluster/scripts/run_summary.sh` | Driver: invokes `08_summary_figures.py` via cluster-env python; tees stdout; lists outputs at end |

**Outputs created in Phase 7 (under `cluster/bap1_late/figures/summary_figures/`):**

| Path | What |
|------|------|
| `dashboard/cluster_dashboard.{png,pdf,svg,jpg}` | 6-panel cluster summary (logFC / direction / Polycomb anchor-vs-span / size / classification / K119ub) |
| `mechanism/anchor_vs_span_mechanism.{png,pdf,svg,jpg}` | clust5 vs clust6 horizontal-bar comparison with mechanism labels |
| `heatmap/feature_heatmap.{png,pdf,svg,jpg}` | z-scored 6×12 heatmap with raw-value annotations + top-3 genes |
| `heatmap/feature_values.tsv` | Raw (non-z-scored) feature matrix for Fig 3 — 6 rows × 12 cols |

**Re-running:**
```bash
LOG=cluster/phase8_summary.txt bash cluster/scripts/run_summary.sh
# Note: log file path is `phase8_summary.txt` even though the phase is 7 in this plan;
# the path was set before phase renumbering and isn't worth changing.
```

---

## Phase 8: Oriented Anchor Metagene + Asymmetry Quantification — DONE (2026-04-27)

**Status:** Three-step pipeline tests whether K27me3 (and the other 7 BigWig marks) is asymmetrically enriched on the EXTERIOR (away from loop body) vs INTERIOR (toward loop body) side of loop anchors. (1) `09_oriented_anchor_metagene.py` builds per-cluster oriented BED6 files (anchor1=strand+ / anchor2=strand-) and runs deepTools `bed_pileup` with the same 8 BigWigs and ±5kb window as Phase 5. deepTools natively flips signal arrays for strand- regions, so plot left = -5kb (exterior) / right = +5kb (interior) for both anchor types. (2) `quantify_orientation_asymmetry.py` reads the resulting `oriented_anchors_values` matrix, computes per-cluster per-mark exterior vs interior means, and runs Wilcoxon signed-rank tests. Outputs `asymmetry_quantification.tsv` (48 rows = 6 clusters × 8 marks). (3) `visualize_orientation_asymmetry.py` produces a focused dual-panel figure (per-cluster K27me3 metagene profiles at top; ext/int ratio bar chart with significance asterisks at bottom). Total runtime: ~1.7h for the metagene step (matches Phase 5), <1 min for quantify, <1 min for visualize.

**KEY biological finding:** K27me3 is **INTERIOR-enriched** at clust4 (mod gain) and clust5 (strong gain) anchors — Ext/Int ratio 0.93 / 0.91 with p<0.005 in both ctrl and mut. Polycomb spreads INTO the loop body, supporting the Polycomb-domain compaction model: clust5 anchors sit AT THE EDGE of expanding Polycomb domains. Conversely, **clust6 (strong loss) shows NO K27me3 asymmetry** (Ext/Int = 1.02 ctrl, 1.02 mut, p ≈ 0.52 both) — Polycomb gain at lost-loop anchors is symmetric, NOT a side-specific extrusion barrier. This refines the Phase 4.4 mechanism: clust6's anchor-Polycomb enrichment (2.09×) is from anchors flanked by Polycomb on BOTH sides (consistent with anchor-disruption / "sensitivity model"), while clust5's (6.59×) reflects the anchor sitting at the boundary of a one-sided Polycomb domain (consistent with extrusion impediment / domain-compaction).

**Corrections applied during execution:**

- **Strand encoding convention:** anchor1 (left, x1<y1) gets strand `+`; anchor2 (right) gets strand `-`. After deepTools native strand-flip, plot left ALWAYS = exterior (away from loop body), plot right ALWAYS = interior (toward loop body). This is the simplest correct encoding — alternatives (e.g. centering on anchor midpoint with custom flipping) are more error-prone.
- **Dedup on `(chrom, start, end, strand)`** — different from Phase 5 which dedups on `(chrom, start, end)`. Hub anchors that participate in different loops as both anchor1 AND anchor2 keep BOTH orientation entries (they appear once with `+` and once with `-`). For clust1: 12,298 loops × 2 = 24,596 raw anchors → 20,739 oriented BED rows (930 hub anchors retain dual orientation; 80.5% of raw anchors deduplicated to single-orientation entries). Per-cluster oriented anchor counts: clust1=20,739 (930 hub) / clust2=18,752 (589) / clust3=14,424 (665) / clust4=6,887 (104) / clust5=1,198 (11) / clust6=3,964 (147).
- **Post-computeMatrix retention rate matches Phase 5 (~98%).** computeMatrix dropped some rows due to blacklist overlap or BigWig signal gaps. Rows surviving for asymmetry quantification: clust1=20,373 / clust2=18,423 / clust3=14,137 / clust4=6,749 / clust5=1,171 / clust6=3,860 (1.5–2.4% dropped per cluster, similar to Phase 5).
- **8 BigWigs identical to Phase 5** (4 marks × ctrl/mut; H3K27ac, H3K27me3, H2AK119ub, H3K27me1). Same `vmax_groups` pairing for ctrl/mut visual consistency. Same color_dict (Blues/Reds/Greens/Purples). Same blacklist (mm10 v2). Same ±5kb window with default `--binSize 10`. Same `pileup_type='referencePoint'`. Net effect: oriented metagene heatmap structure mirrors Phase 5's, but rows are flipped for strand-`-` anchors so spatial asymmetry becomes interpretable.
- **Same intermediate-file caveat as Phase 5.** `oriented_anchors_values` (1.6GB tab text) and `oriented_anchors` (100MB gz binary) are intermediate `computeMatrix` outputs needed only by `heatmap_plot` and the quantify script. Add to `.gitignore` along with the Phase 5 intermediates:
  ```
  cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors
  cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values
  ```
- **`xticklabels=['-5kb (exterior)', 'anchor', '+5kb (interior)']`** passed to `bed_pileup()` — deepTools 3.5.5 honors this and labels the heatmap x-axis correctly. Without this override, the default labels are `['-5.0Kb', '0', '+5.0Kb']` which is correct numerically but loses the biological interpretation.
- **Asymmetry index sign convention:** `asymmetry_index = (exterior_mean − interior_mean) / (exterior_mean + interior_mean)`. POSITIVE = exterior-enriched, NEGATIVE = interior-enriched. K27me3 interior enrichment for clust5 → asymmetry_index = -0.045 (ctrl) / -0.048 (mut). Reported in `asymmetry_quantification.tsv`.
- **Wilcoxon signed-rank** uses paired `(exterior_per_anchor, interior_per_anchor)` differences with `alternative='two-sided'` and `nonzero` filter (drops anchors where ext = int exactly). Bonferroni not applied per-mark because each test answers a distinct biological question; per-mark FDR or family-wise correction can be added if needed for publication.
- **Plot styling:** `visualize_orientation_asymmetry.py` uses BIO_ORDER from Phase 7 (loss → unchanged → gain). Shaded background (blue=exterior, peach=interior) makes the asymmetry direction visually obvious. Significance asterisks: `*` p<0.05, `**` p<0.01, `***` p<0.001. Y-axis range `[0.885, 1.055]` is hand-tuned for our K27me3 ratio range; a more general script would auto-fit.
- **Driver `run_oriented_metagene.sh` mirrors `run_phase5.sh`:** prepends cluster env bin to PATH (for `computeMatrix` subprocess), sets `PYTHONUNBUFFERED=1`, invokes `python3 -u`, tees to `cluster/oriented_metagene.txt`. The quantify and visualize scripts are run separately by the user (no driver) — they're fast (<1 min each) and produce stdout-only diagnostic tables that the user can capture as needed.

### 8.1 — `09_oriented_anchor_metagene.py` (the metagene step, ~1.7h on Mac)

`build_oriented_anchor_beds()` (line 80):
1. Read `combined-clusters.txt`
2. For each cluster, build anchor1 BED6 with strand=`+`, anchor2 BED6 with strand=`-`
3. Concat → dedup on `(chrom, start, end, strand)` → sort by `(chrom, start, end, strand)`
4. Write `cluster{1..6}_oriented_anchors.bed`

`bed_pileup()` (deepTools_pipeline.py:25) called with all 8 BigWigs, all 6 oriented BEDs, `up_down=5000`, `pileup_type='referencePoint'`, blacklist, paired ctrl/mut vmax groups, custom xticklabels.

Outputs: `oriented_anchors_values.{png,pdf,svg,jpg}` (the heatmap), `oriented_anchors_values` (1.6GB text matrix), `oriented_anchors` (100MB gz binary).

### 8.2 — `quantify_orientation_asymmetry.py` (the analysis step, <1 min)

`parse_header()` reads the first 3 lines of `oriented_anchors_values` to extract cluster names, sizes, BigWig labels, bin size, and bin count. Then loads the data matrix, splits each BigWig's bins into exterior half (-5kb to 0) and interior half (0 to +5kb), computes per-anchor mean signal in each half, runs Wilcoxon, produces `asymmetry_quantification.tsv`.

Output table columns: `mark, cluster, direction, n_anchors, exterior_mean, interior_mean, ext_int_ratio, asymmetry_index, wilcoxon_pval`.

Stdout: human-readable per-mark significance tables + a `CLUST5 vs CLUST6 COMPARISON` block summarizing the key biological contrast.

### 8.3 — `visualize_orientation_asymmetry.py` (the figure step, <1 min)

Produces `orientation_asymmetry_figure.{png,pdf,svg,jpg}` — a dual-panel composite:
- **Top:** 1×6 grid of K27me3 metagene profiles, one panel per cluster, in BIO_ORDER. Each panel shows ctrl (blue) and mut (red) mean ± SEM lines with shaded backgrounds (blue exterior, peach interior). Per-panel inset: n, ctrl Ext/Int ratio + significance, mut Ext/Int ratio + significance.
- **Bottom:** Bar chart of K27me3 ctrl + mut Ext/Int ratios across all 6 clusters with significance asterisks. Reference line at 1.0 (symmetric); above = exterior-enriched, below = interior-enriched.

### Phase 8 Verification

**Verified 2026-04-27:**

- `cluster/bap1_late/figures/deeptools_input/clust{1..6}_oriented_anchors.bed`: 6 BED6 files, 35KB–606KB. All sorted, all deduplicated on `(chrom, start, end, strand)`. Strand counts roughly balanced (anchor1+ ≈ anchor2-; minor skew from hub-anchor dedup).
- `cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values.{png,pdf,svg,jpg}`: PNG=2.3MB, PDF=1.1MB, SVG=665KB, JPG=399KB. Heatmap structure matches Phase 5's (6 cluster rows × 8 BigWig columns) with reference point at center; clear exterior-vs-interior asymmetry visible by eye for clust4/clust5 K27me3 panels.
- `oriented_anchors_values` 1.6GB tab-text matrix and `oriented_anchors` 100MB gz binary present (intermediate; should be gitignored).
- `cluster/bap1_late/figures/deeptools/oriented_anchors/asymmetry_quantification.tsv`: 48 rows (8 marks × 6 clusters) + header. All 8 marks × 6 clusters covered. Wilcoxon p-values populated for all rows.
- **K27me3 quantitative results (matches the orientation_asymmetry_figure):**
  | Cluster | Direction | n | ctrl Ext/Int | ctrl p | mut Ext/Int | mut p |
  |---------|-----------|---|--------------|--------|-------------|-------|
  | clust1 | unchanged | 20,373 | 1.015 | 0.541 | 1.002 | 0.575 |
  | clust2 | ~unchanged | 18,423 | 0.952 | 5.9e-6 *** | 0.956 | 0.016 * |
  | clust3 | mod loss | 14,137 | 0.968 | 0.433 | 0.983 | 0.554 |
  | clust4 | mod gain | 6,749 | 0.928 | 5.1e-5 *** | 0.922 | 8.6e-6 *** |
  | **clust5** | **strong gain** | **1,171** | **0.914** | **0.004 \*\*** | **0.909** | **0.004 \*\*** |
  | **clust6** | **strong loss** | **3,860** | **1.020** | **0.522 ns** | **1.023** | **0.571 ns** |

  Confirms KEY finding: K27me3 interior enrichment is the gain signature (clust5, clust4, weakly clust2); clust6 has NO asymmetry (symmetric Polycomb at lost-loop anchors).

- **`orientation_asymmetry_figure.{png,pdf,svg,jpg}`:** PNG=269KB, PDF=192KB, SVG=755KB, JPG=200KB. Top row of 6 panels shows per-cluster K27me3 profiles with ctrl+mut overlay. Bottom row bar chart shows the Ext/Int ratios with asterisks marking p<0.05.

**Files created in Phase 8:**

| File | Purpose |
|------|---------|
| `cluster/scripts/09_oriented_anchor_metagene.py` | Phase 8 orchestrator step 1: per-cluster oriented BED6 prep + bed_pileup with 8 BigWigs at ±5kb / strand-aware orientation |
| `cluster/scripts/quantify_orientation_asymmetry.py` | Phase 8 step 2: parse `oriented_anchors_values`, split bins into exterior/interior halves, Wilcoxon test → `asymmetry_quantification.tsv` |
| `cluster/scripts/visualize_orientation_asymmetry.py` | Phase 8 step 3: dual-panel K27me3 figure (per-cluster profiles + ext/int bar chart) |
| `cluster/scripts/run_oriented_metagene.sh` | Driver for step 1 only: PATH + PYTHONUNBUFFERED setup; tees to `cluster/oriented_metagene.txt`. Steps 2 and 3 are run manually by user. |

**Outputs created in Phase 8 (under `cluster/bap1_late/`):**

| Path | What |
|------|------|
| `figures/deeptools_input/clust{1..6}_oriented_anchors.bed` | 6 deduplicated BED6 files (1,198–20,739 rows; strand-aware) |
| `figures/deeptools/oriented_anchors/oriented_anchors` | computeMatrix gz binary matrix (100MB; intermediate, gitignore) |
| `figures/deeptools/oriented_anchors/oriented_anchors_values` | computeMatrix tab-text matrix (1.6GB; intermediate, gitignore) |
| `figures/deeptools/oriented_anchors/oriented_anchors_values.{png,pdf,svg,jpg}` | Combined heatmap (8 BigWigs × 6 clusters, strand-aware) |
| `figures/deeptools/oriented_anchors/asymmetry_quantification.tsv` | 48 rows × 9 cols: per-(mark,cluster) Wilcoxon results |
| `figures/deeptools/oriented_anchors/orientation_asymmetry_figure.{png,pdf,svg,jpg}` | K27me3 dual-panel figure (per-cluster profiles + Ext/Int bars) |

**Re-running:**
```bash
# Step 1 (~1.7h): metagene
LOG=cluster/oriented_metagene.txt bash cluster/scripts/run_oriented_metagene.sh

# Step 2 (<1 min): quantify
/opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/quantify_orientation_asymmetry.py | tee cluster/orientation_asymmetry_stats.txt

# Step 3 (<1 min): visualize (after Step 2)
/opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/visualize_orientation_asymmetry.py
```

---

## Directory Structure

```
cluster/
├── ChromHMM/                          # Existing — ChromHMM.jar + CHROMSIZES/
├── clustering_example_data/           # Existing — Popay example files
├── custom_params.json                 # Existing — matplotlib style params
├── *.py                               # Existing — Popay modules (modified in Phase 0)
│
├── scripts/                           # all runnable scripts
│   ├── 01_build_loop_count_file.py    #   Phase 1: merged loops → Popay format
│   ├── 02_build_mm10_gene_annotation.R#   Phase 1: mm10 gene BED
│   ├── 03_chromhmm_segmentation.sh    #   Phase 2: BinarizeBed + LearnModel
│   ├── 04_clustering.py               #   Phase 3: elbow + Cluster 3.0 + viz
│   ├── 05_grouped_analyses.py         #   Phase 4: 8 downstream sub-analyses (DONE)
│   ├── 06_deeptools_metagene.py       #   Phase 5: anchor metagene (DONE)
│   ├── 07_cooltools_pileup.py         #   Phase 6: cooltools pileup orchestrator (DONE on HPC)
│   ├── 07_cooltools_pileup.sb         #   Phase 6: SLURM wrapper for Expanse
│   ├── cooler_balance.sb              #   Phase 6 prereq: ICE balance mcools (NEW)
│   ├── 08_summary_figures.py          #   Phase 7: lab-meeting composite figures (DONE)
│   ├── 09_oriented_anchor_metagene.py #   Phase 8 step 1: strand-aware anchor metagene (DONE)
│   ├── quantify_orientation_asymmetry.py  # Phase 8 step 2: ext/int Wilcoxon → TSV (DONE)
│   ├── visualize_orientation_asymmetry.py # Phase 8 step 3: K27me3 dual-panel figure (DONE)
│   ├── run_phase1.sh                  #   Runner: data prep
│   ├── run_phase2.sh                  #   Runner: ChromHMM segmentation
│   ├── run_phase3.sh                  #   Runner: clustering
│   ├── run_phase4.sh                  #   Runner: downstream analyses
│   ├── run_phase5.sh                  #   Runner: deepTools metagene
│   ├── run_summary.sh                 #   Runner: Phase 7 summary figures
│   ├── run_oriented_metagene.sh       #   Runner: Phase 8 step 1 only
│   └── utils/
│       └── multi_format_output.py     #   Figure.savefig patch for 4-format output
│
├── logs/                              # NEW — HPC SLURM stdout (synced from Expanse)
│   ├── cooler_balance_<jobid>.out     #   Phase 6 prereq logs
│   └── phase6_pileup_<jobid>.out      #   Phase 6 pileup attempts (one per submission)
│
├── data/                              # input data files
│   ├── late_merged_loop_counts.txt    #   Phase 1 output (Popay format)
│   ├── late_merged_loop_metadata.tsv  #   Phase 1 output (edgeR stats sidecar)
│   └── mm10_knownGene_pp.bed          #   Phase 1 output (gene annotation)
│   # Phase 5/7/8 BigWigs source: /Users/zakiralibhai/sdsc/bigwigs/ (not in repo)
│   # Phase 6 mcools source:    /expanse/.../stripes/stripenn/data/cool/250402/ (HPC)
│
└── bap1_late/                         # all analysis outputs
    ├── cluster3/
    │   ├── elbow_plot.png
    │   └── k-{k}/
    │       ├── data/
    │       │   ├── input_matrix.txt
    │       │   └── combined-clusters.txt
    │       └── figures/
    │           ├── heatmap.{png,pdf}
    │           ├── lineplot.{png,pdf}
    │           ├── boxplot.{png,pdf}
    │           └── stripplot.{png,pdf}
    ├── figures/
    │   ├── annotation/                # Phase 4.6: per-cluster gene lists
    │   ├── ChIP_intersect/            # Phase 4.3, 4.7: DiffBind + ChIP signal
    │   ├── loop_size/                 # Phase 4.1: per-cluster loop-size dist.
    │   ├── loop_classification/       # Phase 4.2: structural/CRE/mixed/uncl.
    │   ├── chromHMM_anchor/           # Phase 4.5: ChromHMM proportions
    │   ├── cluster_differential_status/  # Phase 4.8: cluster x edgeR direction
    │   ├── deeptools_input/           # Phase 5: anchor BEDs + Phase 8: oriented anchor BEDs
    │   ├── deeptools/
    │   │   ├── histone_anchors/       # Phase 5: combined heatmap (DONE)
    │   │   │   └── histone_anchors_values.{png,pdf,svg,jpg}
    │   │   └── oriented_anchors/      # Phase 8: strand-aware metagene (NEW, DONE)
    │   │       ├── oriented_anchors_values.{png,pdf,svg,jpg}
    │   │       ├── asymmetry_quantification.tsv  # 48 rows = 8 marks × 6 clusters
    │   │       └── orientation_asymmetry_figure.{png,pdf,svg,jpg}  # K27me3 dual-panel
    │   └── summary_figures/           # Phase 7: lab-meeting figures (NEW, DONE)
    │       ├── dashboard/             #   cluster_dashboard.{png,pdf,svg,jpg}
    │       ├── mechanism/             #   anchor_vs_span_mechanism.{png,pdf,svg,jpg}
    │       └── heatmap/               #   feature_heatmap.{png,pdf,svg,jpg} + feature_values.tsv
    ├── chromHMM/
    │   ├── peak_beds/                 # 3-col BEDs staged for BinarizeBed
    │   ├── cellmarkfiletable.txt
    │   ├── binarized/                 # BinarizeBed output
    │   ├── learned_model/             # LearnModel output
    │   │   ├── cerebellum_late_12_segments.bed
    │   │   └── emissions_12.txt
    │   ├── 12state_rename_cerebellum.txt  # Manual state naming
    │   ├── anchor_input/              # Per-cluster anchor BEDs for OverlapEnrichment
    │   ├── span_input/                # Per-cluster span BEDs for OverlapEnrichment
    │   ├── {anchor,span}.txt          # Fold enrichment tables
    │   └── {anchor,span}.{png,pdf}    # Heatmaps (Fig 2f equivalent)
    └── cooltools/                     # Phase 6 outputs (DONE)
        └── obs_exp_contacts/          #   obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}
```

---

## Runner Scripts

Each phase gets a `run_phase{N}.sh` bash script in `cluster/scripts/` that runs the steps for that phase sequentially with error checking. These are meant to be run manually one phase at a time (not chained).

### `run_phase0.sh` — Apply bug fixes + verify

```bash
#!/bin/bash
cd "$(dirname "$0")/.."
# Bug fixes are applied via edits (not a script), but this verifies they work:
python3 -c "
import plotting
import cluster_tools
import chromHMM_heatmap
import deeptools_plotting
import bedpe_analysis
print('All imports OK')
"
```

### `run_phase1.sh` — Data preparation

```bash
#!/bin/bash
cd "$(dirname "$0")/.."
mkdir -p data bap1_late/{cluster3,figures/{annotation,ChIP_intersect,deeptools,deeptools_input},chromHMM/{binarized,learned_model,anchor_input,span_input,peak_beds},cooltools}
python3 scripts/01_build_loop_count_file.py
Rscript scripts/02_build_mm10_gene_annotation.R
echo "Phase 1 complete. Verify:"
wc -l data/late_merged_loop_counts.txt data/mm10_knownGene_pp.bed
```

### `run_phase2.sh` — ChromHMM segmentation

```bash
#!/bin/bash
cd "$(dirname "$0")/.."
bash scripts/03_chromhmm_segmentation.sh
echo "Phase 2 complete. Verify:"
wc -l bap1_late/chromHMM/learned_model/*_12_segments.bed
echo "Now inspect emissions_12.txt and create 12state_rename_cerebellum.txt"
```

### `run_phase3.sh` — Clustering

```bash
#!/bin/bash
cd "$(dirname "$0")/.."
python3 scripts/04_clustering.py
echo "Phase 3 complete. Verify:"
wc -l bap1_late/cluster3/k-*/data/combined-clusters.txt
```

### `run_phase4.sh` — Downstream analyses

```bash
#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/05_grouped_analyses.py
echo "Phase 4 complete. Check bap1_late/figures/ and bap1_late/chromHMM/"
```

### `run_phase5.sh` — deepTools metagene

```bash
#!/bin/bash
set -e
cd "$(dirname "$0")/../.."
CLUSTER_ENV_BIN=/opt/homebrew/anaconda3/envs/cluster/bin
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=cluster/scripts/06_deeptools_metagene.py
LOG=${LOG:-cluster/phase5.txt}
{
  "$PYTHON" -u "$SCRIPT"
} 2>&1 | tee "$LOG"
```

### `run_summary.sh` — Phase 7 summary figures

```bash
#!/bin/bash
set -e
cd "$(dirname "$0")/../.."
PYTHON=/opt/homebrew/anaconda3/envs/cluster/bin/python3
SCRIPT=cluster/scripts/08_summary_figures.py
LOG=${LOG:-cluster/phase8_summary.txt}
{
  "$PYTHON" "$SCRIPT"
  echo "--- Output inventory ---"
  for dir in dashboard mechanism heatmap; do
    ls -lh "cluster/bap1_late/figures/summary_figures/${dir}/" 2>/dev/null
  done
} 2>&1 | tee "$LOG"
```

### `run_oriented_metagene.sh` — Phase 8 step 1 only (metagene)

```bash
#!/usr/bin/env bash
# Driver for 09_oriented_anchor_metagene.py only.
# Run quantify + visualize separately after this finishes (~1.7h).
set -e
cd "$(dirname "$0")/../.."
CLUSTER_ENV_BIN=/opt/homebrew/anaconda3/envs/cluster/bin
export PATH="${CLUSTER_ENV_BIN}:${PATH}"
export PYTHONUNBUFFERED=1
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=cluster/scripts/09_oriented_anchor_metagene.py
LOG=${LOG:-cluster/oriented_metagene.txt}
{
  "$PYTHON" -u "$SCRIPT"
} 2>&1 | tee "$LOG"
```

### `cooler_balance.sb` — Phase 6 prerequisite (HPC SLURM)

ICE-balance mcools at a given resolution (writes the `weight` column required by `cooltools.expected_cis`). Runs ctrl + mut concurrently with 2 nproc each.

```bash
#!/bin/bash
#SBATCH --account=csd940 --partition=shared --cpus-per-task=4 --mem=64G --time=02:00:00
#SBATCH --output=cluster/logs/cooler_balance_%j.out
MCOOL_CTRL="$1"; MCOOL_MUT="$2"; RES="${3:-10000}"
source ~/.bashrc; conda activate cluster
cooler balance --nproc 2 "${MCOOL_CTRL}::/resolutions/${RES}" &
cooler balance --nproc 2 "${MCOOL_MUT}::/resolutions/${RES}"  &
wait
```

Usage:
```bash
sbatch cluster/scripts/cooler_balance.sb \
  /expanse/.../ctrl_merged.mcool \
  /expanse/.../mut_merged.mcool \
  10000
```

### `sync_from_hpc.sh` — rsync HPC data

```bash
#!/bin/bash
# Sync BigWigs and mcools from SDSC Expanse (ssh alias: expanse)
# Run from repo root: bash cluster/scripts/sync_from_hpc.sh

REMOTE="expanse"
LOCAL_BW="cluster/data/bigwigs"
LOCAL_MCOOL="cluster/data/mcools"
mkdir -p "$LOCAL_BW" "$LOCAL_MCOOL"

echo "=== Syncing BigWigs ==="
# H2AK119ub (not in peaks/bigwigs/)
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H2AK119ubCtrl.bw" "$LOCAL_BW/"
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H2AK119ubMut.bw" "$LOCAL_BW/"
# H3K27me1
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H3K27me1Ctrl.bw" "$LOCAL_BW/"
rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/bigwigs/H3K27me1Mut.bw" "$LOCAL_BW/"

echo "=== Syncing mcools (Phase 6 only — large files, ~2-5GB each) ==="
echo "Uncomment below when ready for cooltools pileup:"
# rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/data/cool/250402/ctrl_merged.mcool" "$LOCAL_MCOOL/"
# rsync -avP "${REMOTE}:/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn/data/cool/250402/mut_merged.mcool" "$LOCAL_MCOOL/"

echo "Done. Files synced to:"
ls -lh "$LOCAL_BW"/ 2>/dev/null || echo "  (no bigwigs yet)"
ls -lh "$LOCAL_MCOOL"/ 2>/dev/null || echo "  (no mcools yet)"
```

**Note on BigWig paths in `sync_from_hpc.sh`:** The remote paths above assume `/expanse/lustre/projects/csd940/zalibhai/bigwigs/` matches what's on the user's Mac at `/Users/zakiralibhai/sdsc/bigwigs/`. If the HPC paths differ, update the script before running. The H3K27ac/H3K27me3 BigWigs are already on disk at `peaks/bigwigs/macs2.narrow.aug18.dedup/` and don't need syncing.

---

## Files to Create

| File | Phase | Description |
|------|-------|-------------|
| **Analysis scripts** | | |
| `cluster/scripts/01_build_loop_count_file.py` | 1.1 | Merged loops → Popay format count file |
| `cluster/scripts/02_build_mm10_gene_annotation.R` | 1.2 | mm10 promoter BED from TxDb |
| `cluster/scripts/03_chromhmm_segmentation.sh` | 2 | BinarizeBed + LearnModel pipeline |
| `cluster/scripts/04_clustering.py` | 3 | Elbow plot + Cluster 3.0 k-means + visualization (incl. `--min-ratio` / `--max-ratio` for outlier removal) |
| `cluster/scripts/utils/multi_format_output.py` | 3+ | Context-manager savefig wrapper + figure_subfolder helper (PNG + PDF + SVG + JPG, subfolder per figure); shared by Phase 4-5 |
| `cluster/scripts/05_grouped_analyses.py` | 4 | All Phase 4 downstream analyses (ChromHMM, ChIP, DiffBind, gene annotation) |
| `cluster/scripts/06_deeptools_metagene.py` | 5 | Per-cluster anchor BED prep + bed_pileup (8 BigWigs x 6 clusters, single combined heatmap) |
| `cluster/scripts/07_cooltools_pileup.py` | 6 | Loads combined-clusters.txt → cooler-convention bedpe_dict; calls mcool_pileup with viewframe built inline from cooler chromsizes (no UCSC fetch) |
| `cluster/scripts/07_cooltools_pileup.sb` | 6 | SLURM wrapper for Expanse: csd940 / shared / 16 cpus / 64G / 4h, two-root paths |
| `cluster/scripts/cooler_balance.sb` | 6 (prereq) | NEW: SLURM wrapper to ICE-balance mcools at a given resolution (concurrent ctrl+mut). Required when source mcool only has KR weights (no `weight` column from cooler balance). |
| `cluster/scripts/08_summary_figures.py` | 7 | Lab-meeting summary figures: cluster dashboard (6-panel), anchor-vs-span mechanism (clust5/clust6), feature summary heatmap. Reuses Phase 4 artifacts; live K119ub query. |
| `cluster/scripts/09_oriented_anchor_metagene.py` | 8 | Strand-aware anchor metagene: anchor1=+/anchor2=- → deepTools native flip → exterior/interior asymmetry-ready matrix (8 BigWigs × 6 clusters) |
| `cluster/scripts/quantify_orientation_asymmetry.py` | 8 | Reads oriented_anchors_values, splits ±5kb window into ext/int halves, runs Wilcoxon signed-rank → asymmetry_quantification.tsv (48 rows) |
| `cluster/scripts/visualize_orientation_asymmetry.py` | 8 | K27me3-focused dual-panel figure (per-cluster metagene profiles + ext/int bar chart with significance) |
| **Runner scripts** | | |
| `cluster/scripts/run_phase1.sh` | 1 | Data prep + directory setup |
| `cluster/scripts/run_phase2.sh` | 2 | ChromHMM segmentation |
| `cluster/scripts/run_phase3.sh` | 3 | Clustering — accepts `K FILTER_PCT MIN_RATIO MAX_RATIO` positional + `LOG` env var |
| `cluster/scripts/run_phase4.sh` | 4 | Downstream analyses — accepts `ANALYSES` positional (e.g. `4.4` or `all`) + `LOG` env var |
| `cluster/scripts/run_phase5.sh` | 5 | deepTools metagene — prepends cluster env bin to PATH, sets PYTHONUNBUFFERED=1, invokes `python3 -u` |
| `cluster/scripts/run_summary.sh` | 7 | Summary figures — minimal driver, tees stdout to `cluster/phase8_summary.txt` (legacy filename, predates phase renumbering) |
| `cluster/scripts/run_oriented_metagene.sh` | 8 | Oriented metagene step 1 only — same PATH/PYTHONUNBUFFERED setup as Phase 5; tees to `cluster/oriented_metagene.txt`. Steps 2–3 run manually. |

## Files to Modify

| File | Changes |
|------|---------|
| `cluster/plotting.py` | Add `import os`; fix custom_params path (L24-25); remove gene ylim overrides (L132-139); fix joint() logX bug (L502: `ycol`→`xcol`) |
| `cluster/cooltools_called.py` | Add `import os`; replace `plt.style.use('seaborn-poster')` with `'seaborn-v0_8-poster'` (L16); fix custom_params path (L18-19); change default genome to mm10 (L27) |
| `cluster/chromHMM_heatmap.py` | Fix custom_params path (L9-10) |
| `cluster/deeptools_plotting.py` | Add `import os`; fix custom_params path (L8-9); remove HA-NIPBL ylim override (L108) |
| `cluster/cluster_tools.py` | Fix elbow() savefig/show order (L18-19) |
| `cluster/deepTools_pipeline.py` | Add genome param to bam_coverage() with mm10 size (L15) |
| `cluster/bedpe_analysis.py` | Fix default gene path (L67); add FPKM_df None guard (L142-159) |

## Verification Plan

1. **Phase 0:** ✅ All Popay modules import without error (`plotting`, `cluster_tools`, `chromHMM_heatmap`, `deeptools_plotting`, `bedpe_analysis`, `cooltools_called`, `deepTools_pipeline`)
2. **Phase 1:** ✅ `late_merged_loop_counts.txt` has 39,344 rows, 8 columns; `mm10_knownGene_pp.bed` has 24,515 entries
3. **Phase 2:** ✅ Segmentation BED covers all 21 standard chroms (chr1-19 + chrX + chrY), 335,569 segments; 12-state emissions matrix biologically interpretable (Active_Promoter through Quiescent)
4. **Phase 3:** ✅ Elbow plot shows bend at k=4-5; v2 combined-clusters.txt has 38,948 rows, all 6 clusters > 600 loops, captures 99%+ of edgeR differential calls; multi-format outputs (PNG + PDF + SVG + JPG) per figure subfolder
5. **Phase 4.4 (KEY result):** ✅ ChromHMM anchor vs span heatmaps show differential Polycomb enrichment across clusters — clust5 anchor 6.59x + span 3.03x (extrusion impediment); clust6 anchor 2.09x + span 0.94x (sensitivity model)
6. **Phase 4.8:** ✅ Chi-squared p ≈ 0 for cluster x differential status (clust5=97% up, clust6=78% down, clust1=100% unchanged)
7. **Phase 5:** ✅ deepTools metagene shows clust5 anchors enriched for H3K27me3 + H2AK119ub + H3K27me1 with low H3K27ac (Polycomb signature); clust1-3 K27ac-dominant; clust4-5 paired vmax shows visible mut increases in repressive marks. Heatmap PDF/PNG/SVG/JPG all present at `cluster/bap1_late/figures/deeptools/histone_anchors/`
8. **Phase 6:** ✅ Cooltools pileup completed on Expanse (job `48577675`, 18:43 PDT 2026-04-27). Required ICE-balancing the mcools first (new `cooler_balance.sb` script, job `48574993`, ~52 min) and replacing the §6.3 make_arms_robust patch with an inline cooler-chromsizes viewframe. Output: `cluster/bap1_late/cooltools/obs_exp_contacts/obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}`. Visual confirmation: clust5 mut > ctrl central pixel (gain), clust6 mut < ctrl (loss), clust1 ctrl ≈ mut (unchanged).
9. **Phase 7:** ✅ Three composite summary figures produced — cluster dashboard (6-panel), anchor-vs-span mechanism (clust5/clust6), feature summary heatmap. Output: `cluster/bap1_late/figures/summary_figures/{dashboard,mechanism,heatmap}/`. `feature_values.tsv` records the underlying 6×12 feature matrix for reuse.
10. **Phase 8:** ✅ Oriented anchor metagene + asymmetry quantification completed. KEY result: K27me3 INTERIOR-enriched at clust5 (Ext/Int=0.91, p=0.004 ctrl/mut) and clust4 (0.92–0.93, p<5e-5) — Polycomb-domain compaction at gained loops. clust6 (lost loops): NO K27me3 asymmetry (Ext/Int=1.02, p=0.52) — symmetric Polycomb at anchor-disrupted loops, NOT a side-specific extrusion barrier. Output: `cluster/bap1_late/figures/deeptools/oriented_anchors/{oriented_anchors_values,asymmetry_quantification.tsv,orientation_asymmetry_figure}.*`.