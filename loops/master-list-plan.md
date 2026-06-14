# Master Annotated Loop List + CUT&RUN Signal Juggernaut

## Context

We're assembling the paper. We need two flexible, self-contained loop tables (TSVs) for the
**late/adult (250402)** timepoint that fold together "previous information" plus the recent
analyses, so that downstream figures/filters can be driven off one master file instead of
re-joining scattered outputs.

**Two deliverables (both keyed on `loop_id`, joinable to each other):**

1. **Sheet 1 — master annotated loop list:** prior loop annotation + span length + k-means
   cluster (incl. 6-short/6-long) + class (CRE/structural/…) + `SE_anchor` / `shared_anchor` /
   `DEG_status` flags.
2. **Sheet 2 — CUT&RUN/ATAC signal juggernaut:** per-anchor and per-span signal for K119ub,
   K27ac, K27me3, K4me3, ATAC, CTCF, plus KO−WT deltas. Repetition for shared/hub anchors is
   acceptable.

**Spine = the broader significant set, not the 2,910 logFC-filtered set.** Per your guidance,
we filter `merged_all_results.tsv` to **FDR < 0.05** → **7,981 loops** (4,253 up_in_mutant,
3,728 down_in_mutant), and *keep `logFC`* so the |logFC|>0.3 cut can be applied later when a
given analysis wants it. The 2,910 set is just this set further filtered by logFC.

> Note: every `.md`/`README` in the repo is AI-generated; this plan is grounded only in the
> actual scripts and data files I read and the row counts I verified by hand.

---

## Verified data sources (all exist, paths confirmed)

| Purpose | Path | Key facts |
|---|---|---|
| **Spine** (all loops + stats) | `loops/outputs/250402-late_outputs/merged_loops/merged_all_results.tsv` | 39,344 loops; `significant`==TRUE ⇔ FDR<0.05 → **7,981**; has `logFC`, per-res FDR/logFC, `direction`, BEDPE coords |
| k-means clusters | `cluster/outputs/bap1_late/cluster3/k-6/data/combined-clusters.txt` | 38,948 loops; cols `GROUP chr1 x1 x2 chr2 y1 y2 …`; join on `(chr1,start1,end1,chr2,start2,end2)`; **7,935/7,981 match**, 46 → NA |
| clust6 short/long split | derived: `loop_footprint = end2 − start1`, threshold **800,000 bp** (`<` short, `≥` long) | matches `cluster/scripts/10_clust6_subgroup_asymmetry.py` |
| Chromatin-state annotation logic | `loops/scripts/annotate_loops_extended.R` (functions `classify_anchor_type_extended`, `classify_loop_type_extended`) | 8 anchor categories; pure functions reused (see Path note below) |
| Class beds (CRE/structural) | `peaks/CTCF.bed` (32,487), `peaks/beds/H3K27acCerebellumLate2.bed` (15,105), `cluster/data/mm10_knownGene_pp.bed` (TSS±750bp) | exact inputs from `cluster/scripts/05_grouped_analyses.py` §4.2 |
| Super-enhancers | `peaks/Superenhancers_P60.bed` (1,046, adult/P60) | 4-col BED; `Superenhancers_encode.bed` (52) available as alt |
| Shared/"switched" anchors | re-derived on the 7,981 set (see Decisions); reference: `loops/output/shared_anchor_analysis/late/tables/shared_anchors.tsv` (212, from 2,910) | maxgap 10 kb; lost=down, gained=up |
| DEG list | `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` | cols `ensembl_gene_id`(=symbol), `log2FoldChange`, `padj` (per `load_rnaseq_degs()`) |
| DEG gene bodies | `biomodal/downstream/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz` | cols `Chromosome Start End Annotation Name`; **chrom has NO `chr` prefix** → prepend; `Name`=symbol → direct match |
| Signal bigwigs (5 marks ctrl+mut) | `/Users/zakiralibhai/sdsc/bigwigs/` → `H2AK119ub{Ctrl,Mut}.bw`, `H3K27ac{Ctrl,Mut}.bw`, `H3K27me3{Ctrl,Mut}.bw`, `H3K4me3{Ctrl,Mut}.bw`, `ATAC{ctrl,mut}.bw` | all non-zero; pre-normalized tracks |
| CTCF bigwig | `/Users/zakiralibhai/sdsc/bigwigs/CTCF.bw` (433 MB) | **ctrl-only, ENCODE (different Bl/6 adult mice)** → no mut/delta; paper must state this |
| Signal method (reuse) | `loops/scripts/preprocess_k119ub_anchor_signal.R` (`compute_mean_signal`: `import.bw(which=)` + `coverage(weight="score")` + `viewMeans`) | mean signal per region; cross-validate against its output `peaks/k119ub_anchor_signal.tsv` |

**Path note (important):** `loops/peaks/` does not exist, so `annotate_loops_extended.R`'s
mixed relative paths are CWD-fragile. New scripts run **from the repo root** with explicit
paths; the two *pure* classification functions are reproduced inline (with attribution), not
sourced.

---

## Build — Sheet 1 (`build_master_loop_table.R`, run from repo root)

1. **Spine:** read `merged_all_results.tsv`; keep rows where `significant==TRUE` (FDR<0.05) →
   7,981. Carry forward stats: `loop_id, chr1,start1,end1,chr2,start2,end2, coord_string,
   logFC, logCPM, F, PValue, FDR, significant, exploratory, category, direction, resolution_kb,
   n_resolutions_detected, resolutions_list, is_multi_resolution, FDR_5kb/10kb/25kb,
   logFC_5kb/10kb/25kb, kept_from_resolution, n_overlaps`.
2. **Previous information (re-derive on all 7,981):** build anchor1/anchor2 GRanges; load the 7
   late peak beds (`peaks/beds/H3K27ac/H3K27me3/H3K4me1/H3K4me3 …Late*.bed`,
   `peaks/beds/Bivalent_Cerebellum_Late.bed`, `peaks/CTCF.bed`, `peaks/ctcf_motifs_mm10.bed`);
   `countOverlaps>0` → 7 overlap flags/anchor; TSS distance via
   `TxDb.Mmusculus.UCSC.mm10.knownGene`; apply the reproduced `classify_anchor_type_extended`
   (8 categories, motif-driven CTCF_Site, TSS≤2kb) and `classify_loop_type_extended`. Yields
   `anchor{1,2}_{H3K27ac,H3K27me3,H3K4me1,H3K4me3,Bivalent_Promoter,CTCF,CTCF_motif}_overlap`,
   `anchor{1,2}_distance_to_tss`, `anchor{1,2}_type`, `loop_type`. Also add nearest gene name
   per anchor from gencode vM25 (cheap, since it's loaded for DEG).
3. **Span length:** `span_length = start2 − end1` (interior; all 7,981 are intra-chromosomal
   with span_length ≥ 30 kb, so no NA/edge cases). Also add `loop_footprint = end2 − start1`
   (used only for the clust6 split).
4. **Cluster:** left-join `combined-clusters.txt` on the 6 coord columns →
   `cluster` ∈ {clust2…clust6, NA(46)}. `cluster_label` = `cluster`, except clust6 →
   `clust6_short` (`loop_footprint < 800000`) / `clust6_long` (`≥ 800000`).
5. **Class (CRE/structural/mixed/unclassified):** reproduce `05_grouped_analyses.py` §4.2 exactly
   — per anchor `CTCF` = overlap `peaks/CTCF.bed`; `EorP` = overlap H3K27ac-late **or**
   `mm10_knownGene_pp.bed`. With `CTCF_sum`/`EorP_sum` over the two anchors (0–2):
   structural = `CTCF_sum==2 & EorP_sum<2`; CRE = `CTCF_sum<2 & EorP_sum==2`;
   mixed = both `==2`; else unclassified. Output `loop_class` (+ component flags for audit).
6. **`SE_anchor`:** TRUE if anchor1 or anchor2 overlaps `Superenhancers_P60.bed`.
7. **`shared_anchor`:** re-derive on the 7,981 set — cluster all significant-loop anchors at
   maxgap 10 kb; a locus is "shared" if it contains ≥1 down-loop anchor **and** ≥1 up-loop
   anchor (the loop-switching definition from `shared_anchor_analysis.R`). Flag a loop TRUE if
   either anchor sits in a shared locus.
8. **`DEG_status`:** parse the adult XLSX (`readxl`), keep `padj < 0.05` (no logFC cut, per your
   answer) → DEG symbols; map to gene bodies via gencode vM25 (prepend `chr`); TRUE if anchor1
   or anchor2 overlaps any DEG gene body. Also emit `DEG_genes` (matched symbols) and report the
   symbol→coordinate match rate (fail-loud, no silent drops).
9. Write `master_annotated_loops.tsv`.

---

## Build — Sheet 2 (`build_anchor_span_signal.R`, run from repo root)

1. Read the 7,981 spine (loop_id + coords). Define three region sets per loop: **anchor1**
   (`chr1,start1,end1`), **anchor2** (`chr2,start2,end2`), **span** (`chr1, end1, start2` =
   interior).
2. For each bigwig, compute **mean signal per region** via the reused `compute_mean_signal`
   pattern (`import.bw(which=)` + `viewMeans`) over the unique region set, then map back to each
   loop's anchor1/anchor2/span (hub anchors repeat — accepted).
   - 5 marks × {ctrl, mut}: K119ub, K27ac, K27me3, K4me3, ATAC.
   - CTCF: `CTCF.bw` ctrl only.
3. Per region × mark: `…_ctrl`, `…_mut`, `…_delta = mut − ctrl`,
   `…_log2fc = log2((mut+pc)/(ctrl+pc))` (adaptive pseudocount = 5th-pctile of non-zero, as in
   the existing script). CTCF: `…_CTCF_ctrl` only (no mut/delta).
4. Write `anchor_span_signal.tsv`.

**Column block (× regions anchor1/anchor2/span):** `{mark}_ctrl, {mark}_mut, {mark}_delta,
{mark}_log2fc` for the 5 marks + `CTCF_ctrl` → 3 × (5×4 + 1) = **63 signal columns** + ids.
Runs locally (bigwigs are on this Mac); ~unique 20k regions × 11 bigwig reads.

---

## Decisions & assumptions (please confirm at approval; all are parameters)

- **DEG = padj<0.05 only**, no |log2FC| cut (existing pipeline used |log2FC|>0.3). Switchable.
- **`shared_anchor` re-derived on the 7,981 set** (maxgap 10 kb), not the 2,910-based
  `shared_anchors.tsv` — consistent with the broader-spine philosophy. The original 212-anchor
  file remains available if you'd rather just overlap that.
- **Signal = raw per-condition mean** (bigwigs already normalized); `delta=mut−ctrl` + `log2fc`.
  No cross-condition median-ratio scaling by default (the K119ub script does it; can be toggled).
- **SE source = `Superenhancers_P60.bed`** (adult). `_encode` (52) available as alternative.
- **46 loops will have `cluster=NA`** (dropped by the clustering's count/ratio filter) — class,
  SE, DEG, span, and signal are still computed for them; nothing is fabricated.
- **CTCF caveat:** ctrl-only ENCODE track from different (Bl/6 adult cerebellum) mice — surfaced
  as a single `CTCF_ctrl` column per region, no delta; to be stated explicitly in the paper.
- Two distinct "sizes" are kept and named to avoid confusion: `span_length` (interior, the
  requested span) vs `loop_footprint` (end2−start1, used only for the clust6 split).

---

## New files

- `loops/scripts/build_master_loop_table.R` → `loops/output/master_loop_table/late/master_annotated_loops.tsv`
- `loops/scripts/build_anchor_span_signal.R` → `loops/output/master_loop_table/late/anchor_span_signal.tsv`
- (intermediate) `…/late/significant_loops.tsv`

R env: the Bioconductor stack used by `loops/` (`rtracklayer`, `GenomicRanges`, `InteractionSet`,
`TxDb.Mmusculus.UCSC.mm10.knownGene`, `readxl`, `tidyverse`) — your `mariner_env` R. I'll give
you the exact `Rscript` commands; I won't run them.

## Verification

1. **Counts:** master = 7,981 rows; signal = 7,981 rows; `cluster` non-NA = 7,935, NA = 46.
2. **Spot-check** `loop_32` (chr8:11,635,000–11,640,000 … 11,695,000–11,700,000): expect
   `cluster=clust5` (matches `combined-clusters.txt` row 1), `direction=up_in_mutant`.
3. **Cross-validate signal:** compare Sheet 2's K119ub per-anchor ctrl/mut means against the
   existing `peaks/k119ub_anchor_signal.tsv` for shared anchors (should track closely).
4. **Sanity distributions:** print `loop_class`, `cluster_label`, and TRUE-rates for
   `SE_anchor`/`shared_anchor`/`DEG_status`; per-mark non-zero signal fractions; DEG symbol→coord
   match rate. Eyeball before committing thresholds.
5. Confirm both TSVs are <100 MB and join cleanly on `loop_id`.

---

## Status (session handoff)

Purpose of this work: produce the two TSVs above (`master_annotated_loops.tsv` and
`anchor_span_signal.tsv`) for the late/250402 timepoint, in `loops/output/master_loop_table/late/`.

Both scripts are written:
- `loops/scripts/build_master_loop_table.R` (Sheet 1) — **has been run.**
- `loops/scripts/build_anchor_span_signal.R` (Sheet 2) — **not yet completed.**

Sheet 2 is being moved to HPC (Expanse). State of that:
- `build_anchor_span_signal.R` `CONFIG$bigwig_dir` is now set to the HPC path
  `/expanse/lustre/projects/csd940/zalibhai/bigwigs` (override locally with `--bigwig-dir`).
- That HPC bigwig dir contains the 5 marks (K119ub, K27ac, K27me3, K4me3, ATAC) ctrl+mut, but
  **not `CTCF.bw`** — CTCF.bw lives only on the Mac at `/Users/zakiralibhai/sdsc/bigwigs/CTCF.bw`
  and must be copied into the HPC bigwig dir before the run.
- Run from the HPC repo mirror `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/`, in the
  `mariner_env` R (has `rtracklayer`):
  `Rscript loops/scripts/build_anchor_span_signal.R 2>&1 | tee loops/output/master_loop_table/late/build_signal.log`
- Dependency to confirm present on the HPC mirror: the spine file
  `loops/outputs/250402-late_outputs/merged_loops/merged_all_results.tsv`.

Open thread: a local Sheet 2 run errored during `rtracklayer::import.bw()` on the Mac copy of
`H2AK119ubMut.bw` (`zlib data error`); moving the run to HPC against the Expanse bigwig copies is
the current approach.
