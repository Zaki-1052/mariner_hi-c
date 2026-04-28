Here's the full mapping. The 12 ChromHMM emission states were learned from your 5 ChIP-seq marks (H3K27ac, H3K4me3, H3K4me1, CTCF, H3K27me3) on the late/adult cerebellum, then manually named based on their emission profiles:

| State ID | Biological Name       | Marks ON (prob ~1.0)        | Marks OFF                          | Genome % |
| -------- | --------------------- | --------------------------- | ---------------------------------- | -------- |
| **E8**   | Active_Promoter       | H3K27ac + H3K4me3           | K4me1, CTCF, K27me3                | 0.35%    |
| **E5**   | Active_Promoter_Flank | H3K27ac + H3K4me3 + H3K4me1 | CTCF, K27me3                       | 0.20%    |
| **E6**   | Poised_Promoter       | H3K4me3 + H3K4me1           | K27ac, CTCF, K27me3                | 0.08%    |
| **E7**   | Weak_Promoter         | H3K4me3 only                | (slight K27me3 ~0.06)              | 0.11%    |
| **E9**   | Strong_Enhancer       | H3K27ac alone               | K4me3, K4me1, CTCF, K27me3 all off | 0.36%    |
| **E4**   | Active_Enhancer       | H3K27ac + H3K4me1           | K4me3, K27me3                      | 0.77%    |
| **E3**   | Poised_Enhancer       | H3K4me1 alone               | everything else off                | 3.67%    |
| **E12**  | Bivalent_Enhancer     | H3K4me1 + H3K27me3          | K27ac, K4me3                       | 0.31%    |
| **E11**  | Polycomb              | H3K27me3 alone              | everything else off                | 1.81%    |
| **E2**   | CTCF_Boundary         | H3K4me1 + CTCF              | K27ac, K4me3, K27me3               | 0.09%    |
| **E10**  | Insulator             | CTCF alone                  | everything else off                | 0.30%    |
| **E1**   | Quiescent             | none                        | all 5 marks off                    | 91.97%   |

The naming logic follows your project's existing 7-category anchor taxonomy but expands to 12 states because the emission matrix resolved finer distinctions:

- **Promoter states split 4 ways** based on whether K27ac is present (active vs poised), whether K4me1 flanks the TSS, and whether K27me3 co-occurs (weak/bivalent).
- **Enhancer states split 4 ways**: active (K27ac+K4me1), strong (K27ac-only, distal), poised (K4me1-only), bivalent (K4me1+K27me3).
- **Structural states split 2 ways**: CTCF+K4me1 = boundary element (often at active domain borders); CTCF alone = insulator.
- **E1 Quiescent** (92% of genome) = no marks, as expected for a 5-mark model.

These 12 states then feed directly into the Phase 4.4 anchor-vs-span enrichment analysis (the Popay Fig 2f equivalent) where the key finding was: **gained loops (clust5) show Polycomb enrichment at both anchors AND span** (extrusion impediment), while **lost loops (clust6) show Polycomb at anchors but NOT span** (sensitivity model).


---

All figures are under `cluster/bap1_late/`, each in its own subfolder with 4 formats (png, pdf, svg, jpg):

**ChromHMM (Phase 2 + 4.4 — the key anchor vs span enrichment)**
- `chromHMM/anchor.{png,pdf,svg,jpg}` — fold enrichment heatmap at loop anchors
- `chromHMM/span.{png,pdf,svg,jpg}` — fold enrichment heatmap across loop spans
- `chromHMM/learned_model/emissions_12.{png,svg}` — emission matrix
- `chromHMM/learned_model/transitions_12.{png,svg}` — transition matrix
- `chromHMM/learned_model/cerebellum_late_12_overlap.{png,svg}` — genome overlap
- `chromHMM/learned_model/cerebellum_late_12_RefSeqTSS_neighborhood.{png,svg}` — TSS enrichment
- `chromHMM/learned_model/cerebellum_late_12_RefSeqTES_neighborhood.{png,svg}` — TES enrichment

**Clustering (Phase 3, canonical v2)**
- `cluster3/elbow_plot/elbow_plot.*`
- `cluster3/k-6/figures/heatmap/heatmap.*`
- `cluster3/k-6/figures/lineplot/lineplot.*`
- `cluster3/k-6/figures/boxplot/boxplot.*`
- `cluster3/k-6/figures/stripplot/stripplot.*`

**Downstream analyses (Phase 4)**
- `figures/loop_size/loop_size.*` + `loop_size_strip.*`
- `figures/loop_classification/loop_classification.*`
- `figures/chromHMM_anchor/chromHMM_anchor.*` — stacked bar proportions
- `figures/cluster_differential_status/cluster_differential_status.*`
- `figures/ChIP_intersect/anchor_ChIP_box/anchor_ChIP_box.*` — raw ChIP signal
- `figures/ChIP_intersect/differential_binding_{K27ac,K27me3,K119ub}_fc{0p0,0p3}/` — DiffBind overlap (2 FC thresholds each)
- `figures/ChIP_intersect/ChIP_FC_{K27ac,K27me3,K119ub}_fc{0p0,0p3}/` — fold-change boxplots

**v1 diagnostic (preserved, not canonical)**
- `cluster3/elbow_plot_v1_no-ratio-bound/elbow_plot.*`
- `cluster3/k-6_v1_no-ratio-bound/figures/{heatmap,lineplot,boxplot,stripplot}/`


---


No, the clustering does NOT use ChromHMM states or edgeR logFC. It's simpler than that.

**What's actually being clustered:** The raw Hi-C contact count ratio — `mut_merge / ctrl_merge` — for each loop. That's it. Two numbers per loop, normalized so ctrl = 1.0, and k-means groups the 38,948 loops into 6 clusters based on how much their contact strength changed.

The flow is:

1. **Phase 3 (clustering):** Each loop's `[ctrl_merge, mut_merge]` counts are normalized to ctrl → `[1.0, ratio]`. Cluster 3.0 k-means on that 2D input produces 6 groups ranging from strong loss (clust6, median ratio 0.76) to strong gain (clust5, median ratio 1.34).

2. **Phase 4.4 (ChromHMM overlay — separate, downstream):** For each of those 6 clusters, extract the anchor and span coordinates, then run ChromHMM `OverlapEnrichment` to ask: "what chromatin states are enriched at the anchors vs the body of loops in this cluster?" This is where the emission states enter — as a characterization layer on top of clusters that were defined purely by contact change.

3. **edgeR logFC/FDR** is only used post-hoc in Phase 4.8 (the cluster x differential status crosstab) to validate that the contact-ratio clusters align with the statistical calls. They do — clust5 is 97% `up_in_mutant`, clust6 is 78% `down_in_mutant` — but those labels didn't inform the clustering.

So the clusters answer "which loops changed similarly in contact strength?", and then the ChromHMM overlay answers "what kind of chromatin is at those loops?" — two independent layers.