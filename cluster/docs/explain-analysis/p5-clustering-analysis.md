Please read the @cluster/CONTEXT-CLUSTER.md and @cluster/PLAN-CLUSTER.md - ultrathink about the given context, and proceed to implement this analysis, beginning with planning out p5. We have completed @cluster/phase1.txt @cluster/phase2.txt and @cluster/phase3_v2.txt along with @cluster/phase4.txt. So that it is easier for your context, I also recopied the plan into: @cluster/plan-p1.md (prior work done) and @cluster/plan-p2.md - feel free to update the latter when you are done in the same format as p1 with the completion state (never deleting things but properly noting). We have the results of phase4 and structure in the p2 file.
This serves as important biologically significant work; make sure to pay close attention for maximum scientific accuracy and read the important reference files to adhere to (like in @cluster/). Lastly make sure to follow existing patterns in @cluster/scripts/ - ultrathink and let me know if you have any questions with the planning tool! 
Note, if you launch subagents, make sure you are using sonnet rather than haiku (which fails in this env).

---
## Phase 5: deepTools Metagene Analysis

H3K27me3/H3K27ac/H2AK119ub signal across loop anchors +/-5kb per cluster.

### 5.1 — Prepare per-cluster anchor BED files

For each cluster, extract both anchors (deduplicated) into 3-column BED:
`cluster/bap1_late/figures/deeptools_input/{clust1..clustN}_anchors.bed`

### 5.2 — Run deepTools bed_pileup

Uses `deepTools_pipeline.bed_pileup()` which shells out to `computeMatrix reference-point` + `deeptools_plotting.heatmap_plot()`.

**BigWig dict — need user to confirm/sync these files:**
- H3K27me3: ctrl + mut (from `peaks/bigwigs/macs2.narrow.aug18.dedup/`)
- H3K27ac: ctrl + mut (same directory)
- H2AK119ub: ctrl + mut (user has at `/Users/zakiralibhai/sdsc/bigwigs/H2AK119ub{Ctrl,Mut}.bw` — available locally)

**Blacklist:** `tads/mm10-blacklist.v2.bed` (verified on disk)

**Output:** `cluster/bap1_late/figures/deeptools/histone_anchors.{png,pdf}`

---

All three figures rendered beautifully. Let me review the key biological signals:

**Dashboard (Fig 1):** The biological gradient is immediately clear — logFC goes from negative (red, left) to positive (blue, right). Panel C (KEY) visually nails the two mechanisms: clust5 has both anchor (6.6x) AND span (3.0x) bars above the 1.0 background line, while clust6 has only anchor (2.1x) above with span at background (0.9x). Loop size, classification, and K119ub all track the gradient.

**Mechanism (Fig 2):** The side-by-side comparison is striking. In clust5 (gained), Polycomb and Bivalent_Enhancer have massive solid anchor bars AND visible hatched span bars. In clust6 (lost), those same states have modest anchor-only signal while the span bars are tiny — but all active states (Active_Promoter, Active_Promoter_Flank, Strong_Enhancer) are highly enriched at anchors. The gene annotations show developmental TFs (Ebf3, Foxa2) for gained loops and histone genes for lost loops.

**Heatmap (Fig 3):** The z-score gradient flows cleanly from top (loss = red in loss columns, blue in gain columns) to bottom (gain = opposite). The 87% Polycomb proportion for clust5 and the 6.59x/3.03x Polycomb anchor/span for clust5 jump out in deep red. The K119ub mut column shows clust6 (strong loss) has the highest mutant K119ub signal (1.58), tying the upstream mechanism.

All 3 figures x 4 formats (PNG/PDF/SVG/JPG) = 12 files produced, plus the raw values TSV. The script ran in under 2 minutes including the fresh K119ub BigWig queries.

---

**Summary:** Created `cluster/scripts/08_summary_figures.py` and `cluster/scripts/run_summary.sh`. Three lab-meeting-ready composite figures in `cluster/bap1_late/figures/summary_figures/{dashboard,mechanism,heatmap}/`, each in 4 formats. The dashboard integrates all 6 dimensions (logFC, direction, Polycomb anchor vs span, loop size, classification, K119ub) into a single figure with biological gradient ordering and descriptive cluster labels. The mechanism figure directly contrasts the two models (domain compaction vs anchor disruption) with gene examples. The heatmap provides an at-a-glance z-scored summary of all 12 features.


----


## Image 1: Cluster × Differential Status (Analysis 4.8)

**How it was generated:** `05_grouped_analyses.py` sub-analysis 4.8. It joins `combined-clusters.txt` (the k-means cluster assignments from Phase 3) with `late_merged_loop_metadata.tsv` (the edgeR differential testing results from the main pipeline) on loop coordinates. Then builds a contingency table — for each cluster, what fraction of its loops were called `down_in_mutant`, `unchanged`, or `up_in_mutant` by edgeR (FDR<0.05, |logFC|>0.3). Plotted via `plotting.stacked()`. Chi-squared test: chi2 = 38,987, p ≈ 0.

**What it means:** The k-means clustering (which only sees raw Hi-C contact counts, NOT the edgeR p-values) independently recovers the same directional signal that edgeR found with full statistical testing. This validates that the clusters are biologically real, not noise:
- **clust1** = 100% unchanged (the largest, most stable loops)
- **clust5** = 97% gained in mutant (BAP1-KO creates new contacts here)
- **clust6** = 78% lost in mutant (BAP1-KO destroys contacts here)
- **clust3/4** = moderate loss/gain populations

## Image 2: Loop Classification (Analysis 4.2)

**How it was generated:** `05_grouped_analyses.py` sub-analysis 4.2. For each loop, it checks whether anchors overlap CTCF peaks (`peaks/CTCF.bed`), enhancer peaks (`H3K27acCerebellumLate2.bed`), or promoters (`mm10_knownGene_pp.bed`) using `bioframe.overlap`. Classification rules:
- **structural** = CTCF at both anchors (no enhancer/promoter)
- **CRE** = enhancer or promoter at both anchors (no CTCF)
- **mixed** = one CTCF anchor + one CRE anchor
- **unclassified** = everything else

Plotted as stacked bars per cluster. Chi-squared: chi2 = 1,773, p ≈ 0.

**What it means:** It shows what *type* of loop each cluster contains:
- **clust5** (gained) = 55% structural (CTCF-CTCF), only 3% CRE — these are new structural contacts forming within expanding Polycomb domains, not enhancer-promoter contacts
- **clust6** (lost) = 31% CRE, 29% structural — the lost loops are enriched for regulatory (enhancer-promoter) contacts, meaning BAP1 loss is breaking functional gene-regulatory loops
- **clust4** (moderate gain) = 52% structural — same structural bias as clust5

This connects to the abstract's claim: "long-range loops associated with repression were preferentially lost in exchange for shorter-range contacts." The lost CRE loops are the functional regulatory contacts being disrupted; the gained structural loops are the Polycomb-domain contacts replacing them.


---

