# Step 9b: Paired-Anchor Loop-ABC Analysis

## Position in the Pipeline

Step 9b sits between the independent cross-referencing (Step 9) and the enhancer subset analysis (Step 11). It asks a more stringent geometric question than Step 9:

- **Step 9 (independent overlap):** Do differential loop anchors and ABC-predicted enhancer-gene pairs share the same genomic region? (Answer: 51.4% concordance, essentially chance.)
- **Step 9b (paired-anchor):** Does a single Hi-C loop physically bridge an ABC-predicted enhancer at one anchor and its target gene's TSS at the other? (Answer: 89.7% concordance, p = 1.67e-48.)

The fundamental insight is that requiring **geometric co-localization** — enhancer and target gene TSS occupying opposite anchors of the same loop — transforms the loop-ABC relationship from noise to near-perfect agreement.

### Input Data

| Source | Rows | Description |
|--------|------|-------------|
| All Hi-C loops | 39,344 | Full loop set (visual background) |
| Differential loops | 2,910 | edgeR-called differential loops with logFC, FDR, direction |
| ABC E-G pairs | 180,423 | Distal enhancer-gene pairs with delta-ABC and delta-unnorm |
| mm10 TSS | 28,230 | Gene TSS coordinates |
| RNA-seq DE | 13,588 | Gene-level log2FC and padj |
| K119ub signal | ~15,000 | Per-enhancer H2AK119ub signal (ctrl/mut means, log2fc) |

---

## Matching Algorithm

The script tests two cases for every (ABC pair, loop) combination:
- **Case A:** Enhancer overlaps anchor 1, TSS overlaps anchor 2
- **Case B:** Enhancer overlaps anchor 2, TSS overlaps anchor 1

Overlaps are computed via `GenomicRanges::findOverlaps()` on all 39,344 loops, then annotated with differential loop metadata. Only the 2,910 differential loops are used for statistical analysis; the remaining matches provide visual context.

---

## Results

### Match Statistics

| Metric | Value |
|--------|-------|
| Total unique (ABC pair, loop) matches | 15,112 (from 8,339 unique loops) |
| **Differential loop matches** | **490** (from 304 unique differential loops) |
| Differential loops with >= 1 match | 304 / 2,910 (10.4%) |
| Unique E-G pairs confirmed | 465 / 180,423 (0.26%) |
| Unique genes in matched set | 353 |
| Unique enhancers in matched set | 348 |

**The coverage is narrow but deliberate.** Only 10.4% of differential loops have an ABC-predicted E-G pair where the enhancer and gene TSS sit at opposite anchors. This reflects the stringency of the geometric constraint: most loops connect regions that either (a) don't contain ABC-modeled enhancers, (b) contain enhancers whose ABC-predicted targets are not at the opposite anchor, or (c) involve promoter-promoter interactions outside the ABC framework. The 490 triplets that pass represent the subset where loop biology and ABC modeling are measuring the same physical E-G connection.

### Directional Concordance — The Central Result

| Method | Concordant | Total | Rate | p-value |
|--------|-----------|-------|------|---------|
| Step 9 independent overlap | — | — | 51.4% | (reference) |
| **Paired-anchor dABC** | **269** | **300** | **89.7%** | **1.67e-48** |
| **Paired-anchor d(AxC)** | **426** | **490** | **86.9%** | **1.16e-66** |

When a loop strengthens in KO (up_in_mutant), the ABC score at paired enhancer-gene connections almost always increases (gained). When a loop weakens (down_in_mutant), the ABC score decreases (lost). The jump from 51.4% (independent overlap) to 89.7% (paired-anchor) demonstrates that the geometric requirement is what was missing from the Step 9 cross-reference.

**Partial circularity caveat:** The ABC model uses Hi-C contact frequency as one of its two inputs (Activity x Contact). Loop logFC also measures Hi-C signal changes. So both metrics share a common data source (the .hic files), meaning some concordance is expected by construction. However, the ABC normalization (dividing by the sum of all Activity x Contact for a gene) should partially decouple them — if all contacts around a gene change proportionally, the normalized ABC score stays flat even though individual contacts changed. The 89.7% concordance suggests the local contact change at the specific E-G pair dominates the normalization denominator, which is biologically expected for strong loop-mediated interactions.

The unnormalized d(AxC) concordance (86.9%) is slightly lower than normalized dABC (89.7%), which is unexpected given that d(AxC) performed better in the Step 7 RNA-seq integration (65.3% vs 58.8%). The likely explanation: for these specific loop-anchored pairs, the ABC normalization correctly captures relative changes, while unnormalized scores include more "unchanged" pairs (d(AxC) uses a threshold of 0 rather than 0.01, so more pairs are directional, diluting the concordance rate).

### RNA-seq 3-Way Concordance — The Strongest Validation

Among 68 triplets where the ABC pair is directional AND the target gene is differentially expressed (padj < 0.05):

| Comparison | Concordant | Total | Rate | p-value |
|-----------|-----------|-------|------|---------|
| Loop-ABC 2-way | 63 | 68 | 92.6% | — |
| **Loop-ABC-DE 3-way** | **60** | **68** | **88.2%** | **1.69e-45** |
| Gene-level 3-way | 39 | 44 | 88.6% | — |

**This is the most important result in the entire analysis.** RNA-seq is fully independent of Hi-C data. When a loop weakens, the ABC score at the paired enhancer decreases, AND the target gene's expression decreases — all three agree 88% of the time. The chance expectation for 3-way agreement (all three same direction) is 12.5% (one of eight direction combinations). Observed 88.2% vs expected 12.5% gives p = 1.69e-45.

The sample size is small (68 triplets, 44 unique genes), which limits generalizability. But the effect is so extreme (7x above chance) and statistically unambiguous that it establishes the principle: **physically loop-connected enhancer-gene pairs with differential loop strength show concordant ABC and expression changes.**

### K119ub at Paired Enhancers

| ABC direction | N | Median K119ub log2fc |
|--------------|---|---------------------|
| Lost (dABC < -0.01) | 164 | **+0.321** |
| Gained (dABC > 0.01) | 136 | **-0.074** |

Spearman correlation: **rho = -0.401, p = 5.48e-13** (n = 490 quantifiable)

Enhancers that lose their E-G connections (negative dABC) tend to gain H2AK119ub (positive K119ub log2fc), while enhancers that gain connections show slight K119ub loss. This is consistent with the model: BAP1-KO causes K119ub accumulation, which disrupts enhancer-gene linkages.

**Context:** This rho = -0.401 is dramatically stronger than the rho = -0.035 seen in Step 11's tertile analysis. The difference: Step 11 examined all K119ub_Only enhancers regardless of whether they connect to a gene via a loop. Step 9b restricts to enhancers where the geometric constraint confirms a physical loop-mediated connection. The signal becomes much clearer when noise from non-loop-connected enhancers is removed.

### K119ub by Loop Type

K119ub signal varies significantly across loop types (Kruskal-Wallis p = 6.33e-18):

| Loop Type | Median K119ub log2fc | N |
|-----------|---------------------|---|
| Repressed_Promoter-Polycomb | **-0.335** | 39 |
| Repressed_Promoter-Poised_Enhancer | -0.178 | 32 |
| Other-Other | -0.116 | 15 |
| Active_Promoter-Polycomb | -0.084 | 15 |
| Poised_Enhancer-Other | -0.036 | 16 |
| Active_Promoter-Poised_Enhancer | +0.093 | 16 |
| Repressed_Promoter-Other | +0.112 | 23 |
| Active_Promoter-Repressed_Promoter | +0.138 | 19 |
| Active_Promoter-Active_Enhancer | +0.222 | 63 |
| **Active_Promoter-Active_Promoter** | **+0.231** | **133** |
| Active_Enhancer-Poised_Enhancer | +0.238 | 19 |
| Active_Promoter-Other | +0.316 | 47 |
| **Active_Enhancer-Active_Enhancer** | **+0.388** | 16 |

The gradient is biologically coherent:
- **Repressed/Polycomb loops lose K119ub** (negative log2fc): In BAP1-KO, these already-repressed loci may undergo Polycomb complex redistribution, reducing K119ub at some sites while it accumulates elsewhere.
- **Active loops gain K119ub** (positive log2fc): Active enhancer and promoter regions accumulate K119ub when BAP1 deubiquitinase is lost, consistent with the core BAP1 mechanism.

The Active_Promoter-Active_Promoter class (n=133) is the largest, with median K119ub gain of +0.231. These promoter-promoter loops connecting active genes may be particularly sensitive to K119ub accumulation.

### Distance-Dependent Concordance

| Distance | Concordant | Total | Rate | Binom p |
|----------|-----------|-------|------|---------|
| < 100 kb | 34 | 36 | 94.4% | — |
| 100-250 kb | 77 | 86 | 89.5% | — |
| 250-500 kb | 52 | 61 | 85.2% | — |
| 500 kb-1 Mb | 46 | 49 | 93.9% | — |
| > 1 Mb | 60 | 68 | 88.2% | — |

Concordance is uniformly high at all distances (85-94%), with no decay at longer range. This is important: it rules out a simple distance confound where short-range pairs might trivially agree due to proximity effects. Even at > 1 Mb, paired-anchor concordance is 88.2%.

Spearman distance vs |dABC|: rho = -0.207, p = 3.81e-6. This expected negative correlation confirms that more distant E-G pairs have smaller ABC score changes, consistent with the genomic distance decay of both contact frequency and enhancer influence.

### Effect Size

| Metric | Matched median | All-pairs median | Wilcoxon p |
|--------|---------------|-----------------|-----------|
| |dABC| | 0.0141 | 0.0102 | 3.0e-8 |
| |d(AxC)| | 0.1031 | 0.0611 | 2.3e-21 |

E-G pairs confirmed by paired-anchor loop overlap have ~40-70% larger absolute ABC changes than unfiltered pairs. This makes sense: loop-connected pairs have stronger physical contacts, so changes in loop strength translate to larger ABC perturbations.

### Gene Ontology and Pathway Enrichment

**GO Biological Process — Up_in_Mutant** (genes in strengthening loops):
- Biomineral/skeletal tissue development (q = 0.014): Sox9, Gli2, Tgfbr1, Dhrs3
- Embryonic organ morphogenesis (q = 0.014): Sox9, Foxf2, Myc, Gli2, Mafb
- Neuron fate commitment (q = 0.014): Sox9, Gli2, Lmo4, Tgfbr1, Nrg1
- Cell-substrate adhesion (q = 0.014): Dicer1, Lamc1, Cdk6, Hoxa7

These developmental and morphogenic programs are consistent with BAP1-KO de-repressing genes normally held in check by Polycomb-mediated repression in the cerebellum.

**GO Biological Process — Down_in_Mutant** (genes in weakening loops):
- Heparin proteoglycan biosynthesis (q = 2.2e-6): Angpt1, Ext1, Ext2, Ndst3, Slc10a7
- Negative regulation of DNA recombination (q = 1.2e-5): H1f4, H1f6, H1f5, H1f2, Bcl6, Polq, Cgas
- Heterochromatin formation (q = 1.6e-3): H2ac13, H2ac15, H2ac4, H2ac6 (+ 4 more)
- Negative regulation of gene expression, epigenetic (q = 1.6e-3): H2ac family, Bcl6, N6amt1, Ythdc1

The Down_in_Mutant enrichment is dominated by histone gene families (H2ac, H1f). These histone genes are clustered in the genome, and their downregulation in weakened loops suggests that BAP1-KO disrupts the 3D organization of histone gene loci. The epigenetic regulation terms are biologically coherent with BAP1's chromatin role.

**KEGG Pathways — All Down_in_Mutant** (8 terms):
- Systemic lupus erythematosus (q = 9.4e-18, 20 genes)
- Alcoholism (q = 2.5e-15, 20 genes)
- Neutrophil extracellular trap formation (q = 2.5e-15, 20 genes)

**Important caveat:** These three KEGG pathways share the same 20 genes — all from the histone H2ac/H1f families. SLE, alcoholism, and NET formation pathways in KEGG contain histone gene entries because nucleosome/chromatin biology is part of these disease mechanisms. This is a well-known KEGG annotation artifact where histone gene changes produce these specific pathway hits. The underlying biology (histone gene downregulation in BAP1-KO) is real, but the disease pathway labels should not be interpreted literally.

The remaining KEGG hits (viral carcinogenesis, chromatin remodeling, necroptosis, glycosaminoglycan biosynthesis) are more directly interpretable and consistent with chromatin regulation disruption.

---

## Relationship to Step 11 (Enhancer Subset Analysis)

Steps 9b and 11 examine overlapping biology from opposite angles, producing strikingly different conclusions:

| Aspect | Step 9b (Paired-Anchor) | Step 11 (Enhancer Subsets) |
|--------|------------------------|--------------------------|
| **Unit of analysis** | Loop-enhancer-gene triplets | Enhancers classified by epigenomic marks |
| **Geometric constraint** | Strict: enhancer at one anchor, TSS at the other | Loose: enhancer overlaps any loop anchor |
| **Loop-ABC concordance** | 89.7% (p = 1.67e-48) | Not directly tested at this level |
| **RNA-seq concordance** | 88.2% 3-way (p = 1.69e-45) | K119ub_Only: 49.4% (chance level) |
| **K119ub correlation** | rho = -0.401 (p = 5.5e-13) | rho = -0.035 (p = 0.01) |
| **Sample size** | 490 triplets (304 loops, 353 genes) | 2,479 K119ub_Only enhancers |
| **Conclusion** | Loop-connected E-G pairs are functionally concordant | K119ub_Only contact changes don't reach gene expression |

The apparent contradiction resolves as follows: Step 9b selects for the small fraction of enhancer-gene connections that are physically mediated by a specific, detectable Hi-C loop. These are the highest-confidence E-G connections. Step 11 examines all enhancers — the vast majority of which either (a) are not connected to genes by detectable loops, or (b) have ABC-predicted connections that rely on background contact or powerlaw estimates rather than specific loops. The K119ub effect on contacts is real but small (Step 11: median logFC = -0.054), and it only translates to gene expression changes at the subset of loci where the enhancer is directly loop-connected to its target (Step 9b).

---

## Limitations

1. **Small sample size.** 490 triplets (68 for 3-way concordance) from 304 differential loops means findings are based on ~10% of differential loops and 0.26% of ABC pairs. While statistically unambiguous, the generalizability to the full enhancer-gene landscape is limited.

2. **Shared Hi-C data source.** Both loop logFC and ABC contact scores derive from the same .hic files. Some concordance is expected by construction. The 3-way concordance with independent RNA-seq (88.2%) addresses this, but only for 68 triplets.

3. **Ascertainment bias.** The paired-anchor requirement selects for E-G pairs with strong, detectable loop connections. These are likely the most robust regulatory interactions in the genome — not representative of the typical enhancer-gene relationship, which may operate through weaker or more distributed contacts.

4. **Histone gene cluster effects.** The GO and KEGG enrichment for Down_in_Mutant is heavily influenced by clustered histone genes (H2ac, H1f families). While the biology is real, the enrichment statistics may be inflated by the non-independence of genes in physical clusters.

5. **No FDR stratification reported for exploratory loops.** The summary shows only the "significant" stratum (FDR < 0.05), suggesting all 300 directional matched pairs had FDR < 0.05. This means the analysis cannot assess whether sub-threshold loops show weaker concordance — they simply weren't represented in the matched set.

6. **ATAC-only ABC mode.** As with all ABC analyses in this project, the model used ATAC signal alone (no H3K27ac BAMs), which may underestimate enhancer activity and bias which E-G pairs pass the ABC threshold.

---

## Recommended Figures for Presentation

### Lead with Panel 11 (Multi-Panel Layout) — "The Complete Story in One Slide"

The `paired_anchor_panel/` figure is a 6-panel composite that puts the entire analysis on a single slide. The top row shows the three concordance bar charts (independent vs paired-anchor, FDR stratification, RNA-seq 3-way), and the bottom row shows the three scatter plots (loop logFC vs dABC, loop logFC vs d(AxC), K119ub vs dABC). This is the most efficient way to present Step 9b: the top-left bar (51.4% -> 89.7%) is the hook, and the bottom-right scatter (rho = -0.401) is the mechanistic punchline. If you have one slide for Step 9b, this is it.

### If You Have Two Slides

**Slide 1 — Panel 1 (Concordance Bar) + Panel 3 (logFC vs dABC scatter):**

Panel 1 is the single most impactful figure in the entire analysis. The jump from 51.4% (independent overlap, chance) to 89.7% (paired-anchor) is visually immediate — the grey bar barely clears the 50% dashed line, while the blue bar nearly reaches 90%. This one comparison transforms the narrative from "loops and ABC don't agree" to "loops and ABC agree when the geometry is correct." Pair it with Panel 3 (logFC vs dABC scatter, rho = 0.599) to show the continuous relationship underneath the categorical concordance. The blue/red coloring of differential matches against the grey background of non-differential matches is visually clear, and the Spearman annotation gives the quantitative backbone.

**Slide 2 — Panel 5 (RNA-seq 3-way concordance) + Panel 6 (K119ub scatter):**

Panel 5 shows that the concordance extends to the fully independent RNA-seq modality (88.2% 3-way agreement, where 12.5% is expected by chance). This is the strongest validation because RNA-seq shares no data source with Hi-C. Panel 6 (K119ub at paired enhancers, rho = -0.401) ties it back to the BAP1 mechanism: enhancers that lose E-G connections gain K119ub. The color separation (gained = red, positive dABC, lower K119ub vs lost = blue, negative dABC, higher K119ub) makes the anti-correlation visible without needing to parse the statistics.

### Why Not the Others?

- **Panel 2** (FDR-stratified concordance) only has one bar (all 300 pairs are FDR < 0.05) — visually redundant with Panel 1
- **Panel 4** (logFC vs d(AxC)) shows the same pattern as Panel 3 with slightly stronger rho (0.636 vs 0.599), but the unnormalized metric is harder to explain in a talk. Keep it for supplementary
- **Panel 7** (K119ub by loop type) is rich and biologically interesting (Repressed/Polycomb loops lose K119ub, Active loops gain K119ub), but has 13 categories and is too dense for a talk. Best suited for a paper figure or a targeted question about chromatin state
- **Panel 8** (distance concordance) is a robustness control — concordance is flat at 85-94% across all distance bins. Important for the paper, not for the talk
- **Panel 9** (GO dotplot) has many terms and works better in a paper supplement. The key biological insight (developmental/morphogenic terms for Up_in_Mutant, histone/epigenetic terms for Down_in_Mutant) can be communicated verbally
- **Panel 10** (KEGG dotplot) is dominated by the histone gene cluster artifact (SLE, alcoholism, NET formation all driven by the same 20 histone genes). Showing this in a talk risks derailing into the KEGG annotation caveat rather than the biology

---

## Output File Descriptions

### Plots (11 panels, each in PDF + SVG + JPG subfolders)

| # | File | Content |
|---|------|---------|
| 1 | `paired_anchor_concordance/` | Bar chart: independent overlap (51.4%) vs paired-anchor dABC (89.7%) vs d(AxC) (86.9%) |
| 2 | `fdr_stratified_concordance/` | Concordance by FDR stratum (significant only) |
| 3 | `logFC_vs_deltaABC/` | Scatter: loop logFC vs dABC with grey background (non-diff) and colored differential matches; Spearman annotation |
| 4 | `logFC_vs_delta_unnorm/` | Same as #3 but with unnormalized d(AxC) |
| 5 | `rnaseq_concordance/` | Bar chart: 2-way (92.6%) vs 3-way (88.2%) concordance among DE genes |
| 6 | `k119ub_at_paired_enhancers/` | Scatter: K119ub log2fc vs dABC, colored by ABC direction (rho = -0.401) |
| 7 | `k119ub_by_loop_type/` | Boxplots: K119ub log2fc stratified by loop type and ABC direction |
| 8 | `distance_concordance/` | Bar chart: concordance by enhancer-gene distance bin |
| 9 | `go_bp_dotplot/` | GO Biological Process enrichment, clustered by direction |
| 10 | `kegg_dotplot/` | KEGG pathway enrichment, clustered by direction |
| 11 | `paired_anchor_panel/` | Combined multi-panel layout (plots 1-6) |

Note: Root-level PDFs (e.g., `paired_anchor_concordance.pdf`) are legacy flat-file copies from a prior run. The canonical outputs are in the subfolders.

### TSV Files

| File | Rows | Description |
|------|------|-------------|
| `paired_anchor_all_matches.tsv` | 15,112 | All loop matches (differential + background), with loop coordinates, ABC scores, delta metrics, loop type, and differential status |
| `paired_anchor_matches.tsv` | 490 | Differential-loop matches only, with loop_id, loop_type, logFC, FDR, direction, enhancer/TSS coordinates, ABC scores, delta metrics |
| `paired_anchor_go_enrichment.tsv` | 150 | GO BP enrichment results (up + down clusters), with gene ratios, p-values, and gene lists |
| `paired_anchor_kegg_enrichment.tsv` | 8 | KEGG pathway enrichment (down cluster only), all histone-related |
| `paired_anchor_summary.txt` | — | Complete text summary of all statistical tests and results |

---

## Bottom Line

**What the data shows:** When the geometric constraint is satisfied — enhancer at one loop anchor, target gene TSS at the other — loop direction predicts ABC direction 89.7% of the time, and all three modalities (loop, ABC, RNA-seq) agree 88.2% of the time. K119ub gain at enhancers correlates strongly with E-G linkage loss (rho = -0.401). This subset of 490 triplets provides the clearest mechanistic picture: BAP1-KO causes K119ub accumulation at active enhancers, which weakens their physical loop connections to target genes, reducing enhancer-gene linkage scores and downregulating target gene expression.

**What the data does not show:** This applies to a small, highly curated subset (304/2,910 loops, 465/180,423 E-G pairs). The vast majority of enhancer-gene relationships in the genome do not have a detectable loop bridging them, and for those pairs, the concordance remains at chance (Step 9: 51.4%, Step 11 K119ub_Only: 49.4%). The paired-anchor analysis validates the mechanistic model at its strongest loci but does not demonstrate genome-wide functional impact.

**Interpretation:** The paired-anchor analysis is the "best case" demonstration: when a specific Hi-C loop connects an enhancer to a gene, BAP1-KO-driven changes are highly concordant across all three measurements. The majority of enhancer-gene relationships — which operate through weaker, distributed, or undetectable 3D contacts — do not show this concordance, suggesting that the functional impact of BAP1-KO on E-G linkage is concentrated at a minority of high-confidence loop-mediated connections rather than distributed broadly.
