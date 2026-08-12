# Phase 1B: Figure evaluation decisions

This document captures what was decided in the Phase 1B figure evaluation chat. It builds on `phase1-decisions.md` (the three-beat structure, what's cut, the one-sentence story). A future Claude instance should read both documents together.

---

## Context

Zara uploaded 13 candidate figures, mostly from the Popay-style clustering analysis (Beat 2). The goal was to decide which figures carry the poster story, what format changes they need, and what gets cut.

## Beat 2 figures (locked)

Beat 2 ("clustering reveals two mechanisms") gets **three** poster figures. All will be remade in R for poster format.

### Figure A: Cluster feature heatmap (from Image 11)

**Original:** 6-row × 14-column z-score heatmap with raw value annotations across all 6 clusters.

**Poster version:** 2 rows only (clust5 "strong gain" and clust6 "strong loss"). The other 4 clusters are background; clust5 and clust6 are the ones with actual differential biology. 7 columns:

| Column | Why it stays |
|--------|-------------|
| Median logFC | Direction and magnitude |
| n (loop count) | Grounds comparison in sample size |
| Median size (kb) | Supports Beat 1 (distance shift) |
| % CRE | What kind of loops these are |
| Polycomb anchor | Half of the key mechanistic finding |
| Polycomb span | The other half; anchor-vs-span IS the punchline |
| K119ub mut | The driver mark; accumulation at lost-loop anchors |

**Dropped columns:** % Gained, % Lost (redundant with logFC), % Structural (inverse of % CRE), Bivalent anchor (supporting detail), K119ub ctrl (change matters more than baseline). Top gene labels optional.

**Role in the poster story:** This is the overview figure for Beat 2. It says: "We clustered 39,000 loops by their response to BAP1 loss and identified two extreme populations with distinct properties."

### Figure B: Anchor vs span enrichment (from Image 7)

**Original:** Side-by-side horizontal bar charts showing anchor (solid) vs span (hatched) enrichment for clust5 and clust6, using all 18 chromatin states from the 9-mark ChromHMM model.

**Poster version:** Collapse 18 states to ~4-5 grouped categories:
- **Polycomb/repressive:** K119ub_Only, Polycomb, Polycomb_K119ub, Repressed_Enhancer_K119ub, Heterochromatin
- **Active regulatory:** Active_Promoter, K9ac_Promoter, Active_Enhancer, Active_Enhancer_K9ac, Active_Enhancer_K119ub, Strong_Enhancer
- **Poised:** Poised_Enhancer, K119ub_Poised_Enhancer, Weak_Enhancer, ATAC_Enhancer
- **Structural:** CTCF_Open, Insulator
- (Quiescent can be dropped or shown as a thin bar)

Keep the anchor (solid) vs span (hatched) visual distinction. Keep the bottom labels ("Polycomb domain compaction" vs "Anchor disruption"). These are doing exactly the right work.

**Clust6 long question (decided):** Use pooled clust6, not clust6 long. Reasoning: pooled clust6 shows the "active vs repressed" contrast with clust5, which is the two-mechanism story. Clust6 long makes both panels Polycomb-heavy, turning the comparison into a span-only distinction that's harder for non-specialists to parse. The distance split can be mentioned verbally for judges ("when we split lost loops by size, the long ones show even stronger Polycomb anchor enrichment").

**Role in the poster story:** This is the mechanistic payoff for Beat 2. It says: "Gained loops have Polycomb across the entire domain (anchors AND span). Lost loops have active enhancers/promoters at anchors with no Polycomb enrichment in the span. Two different mechanisms."

### Figure C: Chromatin state composition of differential loops (from Image 9, reformatted)

**Original:** Four pie charts showing 8-category anchor type classification and CRE-vs-structural classification for gained vs lost loops. These are NOT clustering-derived; they compare all edgeR-differential loops directly.

**Poster version:** Produce both a grouped bar chart and a paired donut/pie chart with 4 collapsed categories, then compare which reads better at poster scale:
- **Polycomb:** Polycomb + Repressed_Promoter + Bivalent_Promoter
- **Active:** Active_Promoter + Active_Enhancer
- **Poised Enhancer:** keep as-is
- **CTCF/Structural:** CTCF_Site
- **Other:** thin slice if needed

The grouped bar version gives easier direct comparison; the donut version preserves parts-of-whole. Make both, decide later.

**Role in the poster story:** This bridges Beat 1 and Beat 2. After showing that lost loops are longer (Beat 1), this figure says: "Gained and lost loops sit at fundamentally different types of chromatin." Then the clustering (Figures A and B) explains the mechanism behind that difference. It's the "what" before the "how."

---

## What's cut from this batch

| Figure | Verdict | Where it goes |
|--------|---------|---------------|
| Images 3, 4 (orientation asymmetry profiles) | Paper only | Too many panels, too much context to explain ext/int |
| Images 5, 6 (18-state heatmaps) | Paper only | Raw data underlying Image 7; too many states for poster |
| Image 8 (clust6 short vs long) | Verbal Q&A | Third-order complexity; mention for judges if asked |
| Image 2 (asymmetry bar summary) | Verbal Q&A | Not central enough to justify word budget explaining it |
| Images 12, 13 (direction + loop type stacked bars) | Already in heatmap | Were uploaded as context; info captured in Figure A |

---

## Other figure decisions

**Syt1 locus (Image 1):** Crop to gene model + K119ub track + loop arcs only. Strip all other tracks (5mC, 5hmC, RNA-seq). Drop panels B and C entirely. This becomes a candidate Beat 1 supporting panel showing a concrete locus example of the distance shift. Zara will handle cropping in R.

**Beat 1 and Beat 3 figures:** Not evaluated in this chat. Beat 3 may be as small as a single sentence depending on decisions in a separate chat. Beat 1 figures (loop size distributions, shared anchor hubs) will be evaluated in a future session.

---

## Figure budget so far

The poster (36" × 48") can hold roughly 5-6 figures. Current allocation:

- Beat 1: 1-2 figures (TBD, future chat)
- Beat 1→2 transition: Figure C (chromatin composition)
- Beat 2: Figures A and B (heatmap + anchor-vs-span)
- Beat 3: probably 0-1 figures (may be text only)
- Schematic/model diagram: TBD (the BAP1 → K119ub → two mechanisms causal chain)
- Syt1 locus: small supporting panel if it fits after cropping

That's 4-6 total, which is right for this poster size.

---
---

# Claude Code prompts

These are copy-paste prompts for Claude Code to generate the poster-ready figures in R. Run them from the project root directory.

## Prompt 1: Two-row cluster feature heatmap

```
I need to remake a cluster feature summary heatmap for a poster figure. The original is a 6-row × 14-column z-score normalized heatmap with raw value annotations (generated by 08_summary_figures.py or similar). I need a simplified version.

Specifications:
- 2 rows only: clust5 ("Gained loops, n=667") and clust6 ("Lost loops, n=2,359")
- 7 columns: Median logFC, n, Median size (kb), % CRE, Polycomb anchor enrichment, Polycomb span enrichment, K119ub signal (mutant)
- Z-score normalize across the two rows so color shows relative contrast
- Annotate each cell with its raw value
- Use a diverging red-blue colormap (red = high z-score, blue = low)
- Clean, poster-ready styling: large axis labels (~18pt equivalent), no gridlines, clear row/column labels

The data sources are:
- Cluster assignments: look in the clustering output directory for the combined-clusters file (has logFC, cluster assignment, loop coordinates)
- ChromHMM enrichment: the anchor.txt and span.txt tables in the chromHMM output directory (Polycomb/E11 row for the 5-mark model)
- Loop classification: the loop classification output from Phase 4.2
- K119ub anchor signal: from Phase 4.3 anchor ChIP signal analysis
- Loop sizes: computed from loop coordinates in the cluster assignment file

Search the outputs/bap1_late/ directory tree to find these files. The clustering pipeline README (CLUSTER.md) documents the directory structure.
USER NOTE: @cluster/outputs/bap1_late/figures/summary_figures/heatmap
```

## Prompt 2: Simplified anchor vs span enrichment (collapsed states)

```
I need to remake the ChromHMM anchor-vs-span enrichment figure for clust5 vs clust6, but using the 9-mark/18-state model with states collapsed into ~5 groups.

Specifications:
- Two side-by-side panels: "Gained loops (clust5, n=667)" and "Lost loops (clust6, n=2,359)"
- Horizontal bar chart format (like the existing 18-state version)
- Solid bars = anchor enrichment, hatched/patterned bars = span enrichment
- Collapse the 18 states into these groups:
  * "Polycomb/repressive": K119ub_Only, Polycomb, Polycomb_K119ub, Repressed_Enhancer_K119ub, Heterochromatin
  * "Active regulatory": Active_Promoter, K9ac_Promoter, Active_Enhancer, Active_Enhancer_K9ac, Active_Enhancer_K119ub, Strong_Enhancer
  * "Poised": Poised_Enhancer, K119ub_Poised_Enhancer, Weak_Enhancer, ATAC_Enhancer
  * "Structural": CTCF_Open, Insulator
  * "Quiescent": Quiescent (can be dropped if enrichment is ~1.0 for all)
- For each group, compute a weighted-average fold enrichment across the member states (weighted by genome coverage of each state)
- Add dashed line at fold enrichment = 1.0 (genome baseline)
- Add text labels at bottom of each panel: "Polycomb domain compaction" (clust5) and "Anchor disruption" (clust6)
- Poster-ready styling: large labels, clean axes, no unnecessary gridlines

Data source: the 9-mark ChromHMM anchor and span enrichment tables. Look in the chromHMM_9mark or similar directory under outputs/bap1_late/. The tables have fold enrichment values for each state × cluster combination. The existing figure was generated by 13_mechanism_9mark.py.
USER_NOTE: @cluster/outputs/bap1_late/figures/summary_figures/mechanism_9mark
```

## Prompt 3: Chromatin composition of differential loops (both formats)

```
I need to make two versions of a figure comparing the chromatin state composition of gained vs lost differential loops. This is NOT clustering-based; it uses edgeR differential loop calls directly.

The original figure is a set of pie charts showing 8-category anchor type classification for "Strengthened (Gained in BAP1-KO)" vs "Weakened (Lost in BAP1-KO)" loops. I need both a grouped bar chart version and a paired donut chart version, each with categories collapsed to 4-5 groups.

Collapse the 8 anchor types to:
- "Polycomb": Polycomb + Repressed_Promoter + Bivalent_Promoter
- "Active": Active_Promoter + Active_Enhancer
- "Poised Enhancer": keep as-is
- "CTCF/Structural": CTCF_Site
- "Other": everything else

VERSION A (grouped bar chart):
- Two groups side by side: "Gained in BAP1-KO" and "Lost in BAP1-KO"
- Within each group, bars for each of the 4-5 categories
- Y-axis: percentage of loops
- Clear category colors consistent across groups
- Annotate percentages on or above each bar

VERSION B (paired donut chart):
- Two donut charts side by side
- Same 4-5 categories, same colors as Version A
- Annotate percentages on each slice
- Center of each donut can show n (loop count)

Both versions:
- Poster-ready styling, large labels
- Use a color palette that will be distinguishable at poster scale and for colorblind viewers

Data source: look for the loop anchor classification data. This likely comes from the loop classification analysis (Phase 4.2 of the clustering pipeline) but applied to the edgeR differential loop sets rather than the k-means clusters. The original pie charts were generated from anchor-level chromatin state assignments. If the pre-computed classification exists, use it; if not, the anchor type can be determined by overlapping loop anchors with the ChromHMM segmentation and taking the dominant non-quiescent state.
USER NOTE: @loops/output/loop_annotation_extended/late/plots/anchor_type_distribution_combined
```
