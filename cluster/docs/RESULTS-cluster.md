# RESULTS.md

Results document for the Popay-style Hi-C loop clustering analysis adapted for BAP1-KO cerebellum. Covers the full Phase 0--8 pipeline: statistical outputs from each phase (Section 1), biological interpretation of the key findings (Section 2), and outstanding questions (Section 3).

**Data source of truth:** `docs/phase{1..5}.txt`, `docs/phase3_v2.txt`, `docs/phase4_test_4.4.txt`, `docs/phase8_summary.txt`, `docs/oriented_metagene.txt`, `docs/orientation.txt`, and the output tables under `bap1_late/`.

**Last updated:** 2026-05-02

---

## 1. Statistical Outputs by Phase

### 1.1 Phase 1 -- Data Preparation

**Loop count file** (`data/late_merged_loop_counts.txt`): 39,344 nonredundant loops x 8 columns (`chr1 x1 x2 chr2 y1 y2 ctrl_merge mut_merge`). Drawn from three resolutions: 25kb (16,890), 10kb (14,553), 5kb (7,901).

| Statistic | ctrl_merge | mut_merge |
|-----------|-----------|-----------|
| Mean      | 755       | 720       |
| Median    | 451       | 427       |
| Min       | 5.50      | 6.82      |
| 1%-tile   | 49        | --        |
| 99%-tile  | 4,527     | --        |
| Max       | 10,422    | 38,763    |

**Metadata sidecar** (`data/late_merged_loop_metadata.tsv`): 39,344 x 16. Direction split: unchanged = 31,363; up_in_mutant = 4,253; down_in_mutant = 3,728 (total significant at FDR<0.05 and |logFC|>0.3 = 7,981).

**Gene annotation** (`data/mm10_knownGene_pp.bed`): 24,515 promoter-proximal regions (TSS +/- 750 bp), standard chromosomes only (chr1-19, X, Y). 12,318 on + strand, 12,197 on - strand.

---

### 1.2 Phase 2 -- ChromHMM 12-State Segmentation

**Model:** 12 states learned from 5 histone/structural marks (H3K27ac, H3K4me3, H3K4me1, CTCF, H3K27me3) on 21 standard chromosomes using wildtype/control peak calls.

**Convergence:** 32 iterations, log-likelihood = -1,617,268.978, runtime = 93 sec (4 threads).

**Segmentation:** 335,569 merged segments covering 13,627,597 x 200 bp bins = 2,725,519,400 bp.

**12-state emission profile and genome distribution:**

| State | ID | Marks ON | Genome % | TSS fold | CpG fold |
|-------|----|----------|----------|----------|----------|
| Active_Promoter | E8 | K27ac + K4me3 | 0.35% | 45.07 | 83.76 |
| Active_Promoter_Flank | E5 | K27ac + K4me3 + K4me1 | 0.20% | 19.28 | 21.15 |
| Poised_Promoter | E6 | K4me3 + K4me1 | 0.08% | 19.47 | 20.66 |
| Weak_Promoter | E7 | K4me3 | 0.11% | 51.97 | 88.89 |
| Strong_Enhancer | E9 | K27ac | 0.36% | 7.18 | 12.70 |
| Active_Enhancer | E4 | K27ac + K4me1 | 0.77% | 3.00 | 2.33 |
| Poised_Enhancer | E3 | K4me1 | 3.67% | 2.32 | 1.60 |
| Bivalent_Enhancer | E12 | K4me1 + K27me3 | 0.31% | 23.33 | 49.75 |
| Polycomb | E11 | K27me3 | 1.81% | 2.93 | 3.79 |
| CTCF_Boundary | E2 | K4me1 + CTCF | 0.09% | 9.44 | 8.98 |
| Insulator | E10 | CTCF | 0.30% | 4.66 | 3.94 |
| Quiescent | E1 | none | 91.97% | 0.49 | 0.21 |

**Smoke-test OverlapEnrichment** (all 39,344 anchor1 regions): Active_Promoter = 5.31x, Active_Enhancer = 3.45x, Polycomb = 1.56x, Quiescent = 0.87x (depleted). Biologically coherent.

**No model degeneracy:** Every pair of states differs by at least one mark fully on vs fully off. No state >50% of segments.

---

### 1.3 Phase 3 -- K-Means Clustering

Two runs were performed and preserved:

**v1 (no ratio bounds):** 38,950 loops clustered with k=6. Produced a degenerate 2-loop cluster (clust1) driven by two outlier loops (`chr8:69440000` at mut/ctrl = 4.6 and `chr7:104225000` at 10.08), both flagged unchanged by edgeR (FDR > 0.25). Preserved at `bap1_late/cluster3/k-6_v1_no-ratio-bound/`.

**v2 (canonical, with ratio bounds 0.333-3.0):** 38,948 loops clustered with k=6. Used downstream by all subsequent phases.

| Cluster | n | % of total | Median mut/ctrl | Median logFC | edgeR direction (row %) |
|---------|---|-----------|-----------------|--------------|------------------------|
| clust1 | 12,298 | 31.6% | 0.93 | -0.03 | 100% unchanged |
| clust2 | 10,970 | 28.2% | 1.01 | +0.09 | 92% unchanged, 8% up |
| clust3 | 8,738 | 22.4% | 0.86 | -0.14 | 79% unchanged, 21% down |
| clust4 | 3,916 | 10.1% | 1.12 | +0.25 | 30% unchanged, 70% up |
| clust5 | 667 | 1.7% | 1.34 | +0.50 | 3% unchanged, **97% up** |
| clust6 | 2,359 | 6.1% | 0.76 | -0.31 | 21% unchanged, **78% down** |

**Differential capture:** 3,703 of 3,728 down_in_mutant (99.3%) and 4,232 of 4,253 up_in_mutant (99.5%) retained in v2.

**Resolution balance:** All 6 clusters draw from all 3 resolutions proportionally; no resolution dominates any cluster.

**Elbow plot:** SSE drops from 550 (k=1) to 100 (k=4) to 50 (k=6) to 30 (k=8); inflection at k=4-5. k=6 selected to match Popay precedent and sits at the inflection point.

---

### 1.4 Phase 4 -- Downstream Analyses (8 Sub-Analyses)

Total runtime: ~6 minutes on Mac. All 8 sub-analyses produced multi-format outputs (PNG + PDF + SVG + JPG).

#### 4.1 Loop Size

Median loop sizes per cluster (kb): clust1 = 200, clust2 = 190, clust3 = 300, clust4 = 290, clust5 = 350, **clust6 = 575**. Kruskal-Wallis p ~ 0; all 15 pairwise Wilcoxon comparisons significant after Bonferroni. clust6 (strong loss) has the longest median, consistent with the paper's finding that lost loops are systematically longer (median 625 kb overall).

#### 4.2 Loop Classification

chi2 = 1,773; p ~ 0; dof = 15.

| Type | clust1 | clust2 | clust3 | clust4 | clust5 | clust6 |
|------|--------|--------|--------|--------|--------|--------|
| Structural | 32.2% | 40.4% | 31.9% | 52.3% | **54.7%** | 29.3% |
| CRE | 24.7% | 16.6% | 28.3% | 7.2% | **2.8%** | **30.9%** |
| Mixed | 12.2% | 11.0% | 10.6% | 4.3% | 2.1% | 7.3% |
| Unclassified | 30.9% | 32.0% | 29.2% | 36.2% | 40.3% | 32.5% |

Gained loops (clust5) are predominantly structural (CTCF-CTCF); lost loops (clust6) are enriched for CRE (enhancer-promoter).

#### 4.3 Anchor ChIP Signal

All 8 mark x condition combinations (H3K27ac, H3K27me3, H2AK119ub, H3K27me1 -- ctrl and mut) show Kruskal-Wallis p << 0.001 across clusters. H2AK119ub_mut omnibus statistic (926) > H2AK119ub_ctrl (361), indicating cluster-level K119ub variance increases in mutant.

#### 4.4 ChromHMM Anchor vs Span Enrichment (KEY RESULT)

**Polycomb (E11) fold-enrichment by cluster:**

| Cluster | Anchor | Span | Anchor/Span ratio | Direction |
|---------|-------:|-----:|-----------------:|-----------|
| clust1 | 0.96 | 1.00 | 0.96 | 100% unchanged |
| clust2 | 1.80 | 1.23 | 1.46 | 92% unchanged |
| clust3 | 0.84 | 0.89 | 0.94 | 79% unchanged |
| clust4 | 3.91 | 1.69 | 2.31 | 70% up |
| **clust5** | **6.59** | **3.03** | **2.18** | **97% up** |
| **clust6** | **2.09** | **0.94** | **2.22** | **78% down** |

**Bivalent_Enhancer (E12) enrichment at clust5:** anchor = 7.91x, span = 2.21x.

Both gain and loss clusters show elevated anchor-Polycomb, but **only the gain cluster (clust5) shows elevated span-Polycomb**. Lost loops (clust6) have Polycomb at anchors with the span at genome baseline (~1.0x).

Active states (Active_Promoter, Active_Enhancer, Strong_Enhancer) show the inverse pattern: enriched at clust6 anchors (7.22x, 4.42x, 4.61x) but depleted at clust5 anchors (1.84x, 0.54x, 0.55x).

#### 4.5 ChromHMM Proportions (Stacked Bar)

Most-prominent state per anchor (Quiescent excluded):

| State | clust1 | clust2 | clust3 | clust4 | clust5 | clust6 |
|-------|--------|--------|--------|--------|--------|--------|
| Active_Promoter | 37.4% | 33.4% | 43.1% | 17.9% | 1.4% | 35.5% |
| Active_Enhancer | 29.8% | 27.7% | 29.4% | 5.6% | 0.0% | 25.5% |
| Polycomb | 8.8% | 17.5% | 6.7% | 60.6% | **87.0%** | 25.5% |
| Bivalent_Enhancer | 0.6% | 1.1% | 1.0% | 5.6% | 7.3% | 0.9% |

clust5 (strong gain) is overwhelmingly Polycomb-dominated (87% + 7% Bivalent = 94% repressive). clust6 (strong loss) is heterogeneous: 36% Active_Promoter, 26% Active_Enhancer, 26% Polycomb.

#### 4.6 Gene Annotation

Per-cluster gene lists: clust1 = 7,298 unique genes (16,217 loop-gene pairs); clust2 = 5,949; clust3 = 4,939; clust4 = 1,759; clust5 = 237; clust6 = 1,241. clust5's small gene count is consistent with its 87%-Polycomb signature (most anchors do not overlap promoters).

#### 4.7 DiffBind Relationship

All 3 marks show significant differential binding enrichment across clusters (chi2 >> 1,000, all p ~ 0):

| Mark | chi2 (|Fold|>0.3) | Significant peaks | Total peaks |
|------|--------------------|-------------------|-------------|
| H3K27ac | 3,488 | 11,647 | 25,669 |
| H3K27me3 | 3,009 | 7,103 | 18,324 |
| H2AK119ub | 1,693 | 20,512 | 41,392 |

Fold threshold 0.3 vs 0.0 produced structurally identical results (peak count differences <= 7%).

#### 4.8 Cluster x Differential Status Crosstab

chi2 = 38,987; p ~ 0; dof = 10. The k-means clustering (which sees only raw contact counts) independently recovers the edgeR differential calls: clust5 = 97.3% up_in_mutant, clust6 = 78.5% down_in_mutant, clust1 = 100% unchanged.

---

### 1.5 Phase 5 -- deepTools Metagene

Per-cluster anchor metagene profiles at +/- 5kb for 4 marks x ctrl/mut = 8 BigWigs. Runtime: 96 min on Mac.

**Anchor counts** (deduplicated per cluster): clust1 = 19,809; clust2 = 18,163; clust3 = 13,759; clust4 = 6,783; clust5 = 1,187; clust6 = 3,817. Total = 63,518 unique anchors (81.5% of 77,896 raw anchor entries). computeMatrix dropped 1,211 (1.9%) regions due to blacklist or BigWig gaps.

**Visual summary:** clust5 anchors show elevated H3K27me3 + H2AK119ub + H3K27me1 with low H3K27ac (Polycomb signature). clust1-3 show K27ac dominance. clust4-5 ctrl-vs-mut paired panels show visible mut increases in repressive marks.

---

### 1.6 Phase 6 -- Cooltools Pileup (HPC)

Per-cluster log2(obs/exp) Hi-C pileup at 10kb resolution, +/- 500kb flank, for ctrl_merged and mut_merged mcools. Final successful job: SLURM `48577675`, 35 min runtime on Expanse.

**Required a one-time prerequisite:** ICE balancing both mcools at 10kb (SLURM `48574993`, 52 min, 132 iterations per file) to add the `weight` column expected by `cooltools.expected_cis`.

**Visual result:** Central pixel of clust5 mut > ctrl (gain visible at aggregate level); clust6 mut < ctrl (loss visible); clust1 ctrl ~ mut (unchanged baseline). Background quadrants show log2(obs/exp) ~ 0 in both conditions.

---

### 1.7 Phase 7 -- Summary Figures

Three composite figures integrating Phase 4 outputs. Runtime: 12 sec.

**K119ub live query results** (49,556 unique anchors queried):

| Cluster | Direction | K119ub ctrl median | K119ub mut median |
|---------|-----------|--------------------|-------------------|
| clust6 | strong loss | 1.160 | **1.577** |
| clust3 | mod loss | 1.057 | 1.166 |
| clust1 | unchanged | 1.057 | 0.956 |
| clust2 | ~unchanged | 1.171 | 0.938 |
| clust4 | mod gain | 1.297 | 1.058 |
| clust5 | strong gain | 1.368 | 1.254 |

K119ub_mut is elevated at clust6 (lost-loop anchors) and clust3 (moderate-loss), consistent with BAP1 KO failing to remove K119ub at disrupted anchors. Gain clusters (clust4/5) show high K119ub in both ctrl and mut, indicating these are already Polycomb-marked loci.

---

### 1.8 Phase 8 -- Oriented Anchor Metagene + Asymmetry Quantification

Three-step pipeline testing whether histone marks are asymmetrically enriched on the exterior (away from loop body) vs interior (toward loop body) flanks of loop anchors.

**Oriented anchor counts** (deduplicated on chrom/start/end/strand): clust1 = 20,739 (930 hub-dual entries); clust2 = 18,752; clust3 = 14,424; clust4 = 6,887; clust5 = 1,198; clust6 = 3,964. computeMatrix retention: ~98%.

**H3K27me3 Ext/Int ratios and Wilcoxon signed-rank p-values:**

| Cluster | Direction | n | ctrl Ext/Int | ctrl p | mut Ext/Int | mut p |
|---------|-----------|---|--------------|--------|-------------|-------|
| clust6 | strong loss | 3,860 | 1.020 | 0.522 ns | 1.023 | 0.571 ns |
| clust3 | mod loss | 14,137 | 0.968 | 0.433 ns | 0.983 | 0.554 ns |
| clust1 | unchanged | 20,373 | 1.015 | 0.541 ns | 1.002 | 0.575 ns |
| clust2 | ~unchanged | 18,423 | 0.952 | 5.9e-6 | 0.956 | 0.016 |
| clust4 | mod gain | 6,749 | 0.928 | 5.1e-5 | 0.922 | 8.6e-6 |
| **clust5** | **strong gain** | **1,171** | **0.914** | **0.004** | **0.909** | **0.004** |

**H3K27me1 Ext/Int (complementary pattern):** clust4 ctrl = 1.076 (p = 0.001), clust5 ctrl = 1.066 (p = 0.040) -- exterior-enriched, opposite to K27me3. K27me1 marks the flanking chromatin outside the Polycomb core.

**H2AK119ub Ext/Int:** Weak interior enrichment at clust4 ctrl (0.953, p = 0.019); not significant at clust5 (p = 0.26).

---

### 1.9 Phase 11 -- Comprehensive Asymmetry (H2AK119ub, H3K27ac, PC1, Insulation)

Extended asymmetry analysis at two spatial scales: histone marks (H2AK119ub, H3K27ac) at +-5kb anchor-local, and compartment features (PC1 eigenvector, insulation score) at +-50kb domain-scale. PC1 computed via `cooltools eigs-cis` at 25kb (coarsened from 5kb mcool, ICE-balanced). Insulation via `cooltools insulation` at 10kb with 200kb diamond window. Focused on clust5, clust6, and clust6 short/long subgroups. Runtime: 91 min on Expanse (job 48713858).

**PC1 ctrl-mut correlation:** r = 0.9619 (n = 101,854 bins). High consistency confirms phasing-track sign orientation worked correctly.

**H2AK119ub and H3K27ac: No asymmetry** (all p > 0.05 for all 4 groups). K119ub is deposited uniformly at anchor regions without directional spreading. K27ac similarly symmetric.

**PC1 (compartment eigenvector):**

| Cluster | Direction | n | ctrl ext | ctrl int | ctrl p | mut ext | mut int | mut p |
|---------|-----------|---|---------|---------|--------|---------|---------|-------|
| clust5 | strong gain | 1,167 | -0.238 | -0.333 | 6.3e-44 *** | -0.113 | -0.164 | 6.7e-13 *** |
| clust6 | strong loss | 3,858 | 0.167 | 0.178 | 1.3e-4 *** | 0.103 | 0.109 | 0.032 * |
| clust6_short | loss <800kb | 2,502 | 0.269 | 0.302 | 1.0e-13 *** | 0.193 | 0.216 | 4.9e-7 *** |
| clust6_long | loss >=800kb | 1,474 | -0.003 | -0.029 | 0.003 ** | -0.050 | -0.075 | 0.005 ** |

Clust5 gained-loop anchors are in B compartment (negative PC1), with **interior more B** (deeper in heterochromatin). Clust6 lost-loop anchors are in A compartment (positive PC1), with interior slightly more A. Clust6_long straddles the A/B boundary.

**Insulation score (log2, 200kb diamond):**

| Cluster | Direction | n | ctrl ext | ctrl int | ctrl p | mut ext | mut int | mut p |
|---------|-----------|---|---------|---------|--------|---------|---------|-------|
| clust5 | strong gain | 1,171 | -0.275 | -0.175 | 3.5e-54 *** | -0.277 | -0.109 | 1.9e-89 *** |
| clust6 | strong loss | 3,859 | -0.059 | +0.098 | 2.9e-179 *** | -0.128 | +0.014 | 7.7e-151 *** |
| clust6_short | loss <800kb | 2,503 | +0.079 | +0.261 | 5.4e-141 *** | -0.001 | +0.155 | 1.2e-110 *** |
| clust6_long | loss >=800kb | 1,474 | -0.301 | -0.181 | 7.0e-52 *** | -0.352 | -0.229 | 1.2e-50 *** |

Clust5: exterior more insulated (boundary-like) than interior -- gained-loop anchors sit AT strong TAD boundaries, with the Polycomb domain interior being less insulated (self-associating). Clust6: dramatic sign flip across anchor (exterior insulated, interior open), marking a sharp boundary-to-open transition. Clust6_long: both sides boundary-like, exterior more so.

Note: Ext/Int ratio and asymmetry_index columns are unreliable for PC1 and insulation since those signals cross zero. Wilcoxon p-values and raw means are the correct metrics.

---

## 2. Biological Interpretation

This section synthesizes the statistical results above into biological conclusions. Interpretations are grounded in the measured data; speculative extensions are flagged as such.

### 2.1 Two Mechanisms in One Dataset

The central finding is that gained and lost chromatin loops in BAP1-KO cerebellum arise through distinct mechanisms, both identifiable from the ChromHMM anchor-vs-span enrichment analysis (Phase 4.4):

**Gained loops (clust5, 667 loops, 97% up_in_mutant):**
- Polycomb enriched at both anchors (6.59x) AND span (3.03x)
- 87% of anchors have Polycomb as the dominant ChromHMM state
- 55% structural (CTCF-CTCF), only 3% CRE
- Median size: 350 kb
- K27me3 is asymmetrically interior-enriched (Ext/Int = 0.91, p = 0.004)
- PC1: interior more B-compartment than exterior (ext = -0.24, int = -0.33, p = 6.3e-44)
- Insulation: exterior more boundary-like (ext = -0.28, int = -0.17, p = 3.5e-54)
- H2AK119ub and H3K27ac: no asymmetry (uniform at anchors)

These are new contacts forming **within expanding Polycomb domains**. Both anchors sit in heterochromatin, the chromatin between them is heterochromatic. The interior K27me3 asymmetry indicates the Polycomb domain extends preferentially inward from the anchors into the loop body, consistent with PRC-mediated domain compaction. The PC1 and insulation data confirm that clust5 anchors sit at **TAD boundaries at the edge of B-compartment Polycomb domains** -- the loop interior is deeper in heterochromatin (more B, less insulated / more self-associating) while the exterior faces the A/B boundary.

**Lost loops (clust6, 2,359 loops, 78% down_in_mutant):**
- Polycomb enriched at anchors (2.09x) but NOT at span (0.94x, genome baseline)
- 31% CRE, 29% structural -- the highest CRE fraction of any cluster
- Active_Promoter enriched at anchors (7.22x), Active_Enhancer (4.42x)
- Median size: 575 kb -- longest of any cluster
- K27me3 shows NO asymmetry (Ext/Int = 1.02, p = 0.52)
- K119ub_mut elevated (ctrl 1.16 -> mut 1.58)
- PC1: both sides A-compartment, interior slightly more A (ext = 0.17, int = 0.18, p = 1.3e-4)
- Insulation: dramatic sign flip across anchor (ext = -0.06, int = +0.10, p = 2.9e-179)
- H2AK119ub and H3K27ac: no asymmetry

These are existing active/CRE loops whose CTCF anchor sites become invaded by Polycomb upon BAP1 loss. The span remains euchromatic. The symmetric K27me3 at lost-loop anchors indicates Polycomb gain is not side-specific -- it does not create a euchromatin-to-heterochromatin boundary. Instead, the anchor region itself becomes heterochromatinized, disrupting CTCF binding. The insulation data reveals a sharp boundary-to-open transition at these anchors: the exterior is insulated while the interior (loop body) shows open contact structure, consistent with the anchor marking the edge of a structural domain that is being disrupted.

### 2.2 Mapping to the Dixon/Popay Framework

Jesse Dixon suggested testing whether Polycomb enrichment is at anchors vs spans of differential loops. Popay et al. (Nat Genet 2026) found a "mixed-dependency" cluster with K27me3 flanking but not at CTCF anchors -- an extrusion impediment signature where Polycomb blocks cohesin mid-journey.

Our BAP1-KO data shows a more complex picture than either the pure "sensitivity" (anchor-only) or "extrusion impediment" (span-only) model. The key distinction:

- **clust5 (gained)** has Polycomb at anchor AND span, with interior-biased asymmetry. This is consistent with **Polycomb-domain compaction** -- not extrusion blockade per se, but self-association of Polycomb-marked chromatin creating new short-range contacts.
- **clust6 (lost)** has Polycomb at anchor only, with symmetric distribution. This is consistent with the **sensitivity/anchor-disruption model** -- repressive chromatin directly occluding CTCF binding at anchor sites, analogous to the Bernstein IDH-mutant glioma mechanism.

Both mechanisms operate simultaneously in the same tissue, affecting different loop populations.

### 2.3 H2AK119ub as the Upstream Signal

K119ub data connects the two mechanisms to BAP1's enzymatic activity:

- At **lost-loop anchors** (clust6), K119ub_mut is elevated (ctrl 1.16 -> mut 1.58). BAP1 loss leads to K119ub accumulation at these specific anchor sites, which recruits PRC2 -> H3K27me3 -> CTCF displacement.
- At **gained-loop anchors** (clust5), K119ub is already high in both conditions (ctrl 1.37, mut 1.25). These are pre-existing Polycomb-marked loci where the chromatin domain is already repressive; BAP1 loss exacerbates domain expansion.
- The overall pattern is not a uniform mut-up shift but rather increased cluster-level K119ub variance in mutant (Kruskal-Wallis statistic: 926 mut vs 361 ctrl), consistent with locus-specific rather than genome-wide K119ub changes.

### 2.4 Loop Length and Classification Patterns

The loop length asymmetry -- lost loops are longer (median 575 kb for clust6), gained loops are shorter (350 kb for clust5) -- aligns with both the paper's abstract (lost median 625 kb vs gained 320 kb genome-wide) and the mechanistic interpretation:

- Long-range CRE loops are more vulnerable to anchor disruption because their anchors have greater cumulative exposure to Polycomb spreading
- New Polycomb-domain contacts form at short range because they reflect local chromatin compaction

The loop classification data reinforces this: gained loops are 55% CTCF-CTCF (structural contacts within Polycomb domains), while lost loops are 31% CRE (functional enhancer-promoter contacts being broken). This matches the abstract's claim: "long-range loops associated with repression were preferentially lost in exchange for shorter-range contacts."

### 2.5 Progressive Developmental Effect

The early timepoint (P12/250831) has too few differential loops (~200) for meaningful clustering. However, the directional progression -- from ~200 differential loops early to ~3,000 in adult -- is consistent with progressive Polycomb spreading over neurodevelopment. Stripenn data confirms: more lost stripes at P12, more gained in adult; bivalent promoter anchors 3x more prevalent at P12 (3.6%) than adult (1.2%).

### 2.6 Cross-Validation with Stripe and Compartment Data

The clustering results are internally consistent with other analyses performed on the same dataset:

- **Stripe analysis (Stripenn):** Gained stripes are enriched at Polycomb anchors (2.28x vs 0.84x for lost); stripe body genes in gained stripes are enriched for developmental TFs (Hox, Shh, Dlx5/6). Lost stripes concentrate at Active_Enhancer anchors. The stripe and loop populations share anchors non-randomly (permutation z-scores up to 8.2).
- **Compartment shifts:** Gained stripes have lower anchor PC1 (median 0.30, more B-like) than lost stripes (median 1.05). Gained stripes are 2.28x more likely to sit on compartment-switched bins.
- **Cooltools pileup:** Aggregate Hi-C contact strength visibly differentiates clust5 (gain) and clust6 (loss) at the central pixel, confirming the clusters separate at the raw signal level.

### 2.7 CpG Island Methylation at CTCF Anchors (Separate Analysis)

A complementary analysis overlaying biomodal 5mC/5hmC data onto CTCF loop anchors (Section 47 in the biomodal downstream pipeline) found:

- Dynamic CpG shores/shelves at lost CTCF anchors are hypermethylated (OR = 3.28, p < 2.2e-16) with concurrent hmC loss (OR = 2.08, p = 3.9e-6)
- Lost anchors are 2.18x more likely than gained to show coordinated mC-up/hmC-down -- consistent with TET blockade by Polycomb compaction
- The Cochran-Mantel-Haenszel test controlling for loop distance: OR = 2.87 [2.28-3.60], p < 2.2e-16. Strongest at short loops <200 kb (OR = 11.21), attenuates at >1 Mb (OR = 0.87, ns)
- Methylation effect correlates with loop loss: Spearman rho = -0.244 (p = 4.3e-9) between mC change and loop logFC

This suggests DNA hypermethylation at CpG-rich CTCF anchor sites may contribute to CTCF displacement at lost loops, paralleling the IDH-mutant glioma mechanism identified by Bernstein and colleagues.

### 2.8 CpG Island Ubiquitination (Separate Analysis)

A separate analysis (Section 48) testing whether K119ub accumulation drives CpG island methylation yielded a counterintuitive result:

- Hypermethylated CpG islands are **depleted** for K119ub overlap (OR = 0.48-0.49, p < 0.003)
- Hypomethylated CpG islands are **enriched** (OR = 1.91-2.50, p < 1e-8)
- 85% of hypermethylated CpG islands already had mC >= 0.20 in controls (mean baseline = 0.57) -- this is amplification of pre-existing methylation, not de novo methylation
- Hyper islands are depleted for H3K4me3 (OR = 0.11) and H3K27ac (OR = 0.43) -- they are at chromatin-inactive CpG islands

The methylation gains at CpG islands appear to occur independently of direct K119ub accumulation, suggesting an indirect mechanism.

---

## 3. Unanswered Questions and Further Directions

### 3.1 Questions Directly Addressable with Existing Data

1. **CTCF motif centering.** The CTEA April meeting requested re-centering anchor annotation around the actual CTCF motif within the 10 kb anchor window. CTCF motif BED is available (`peaks/ctcf_motifs_mm10.bed`). This would test whether Polycomb gain at clust6 anchors is precisely at the CTCF binding site vs spread across the wider anchor region.

2. **Superenhancer proximity to DEGs.** The CTEA meeting requested testing how frequently DEGs contact superenhancers within a 2 Mb window, sub-classified by K27ac change. SE BEDs exist (`peaks/Superenhancers_P60.bed`, 1,046 SEs).

3. **DEG-centric enhancer contact analysis.** Pull DEGs and expression-invariant genes as controls; compare contact strength to enhancers. This addresses a gap in the existing ABC analysis.

4. **Per-cluster GO/KEGG enrichment.** Phase 4.6 produced per-cluster gene lists but did not run functional enrichment. Testing whether clust5 genes are enriched for developmental/Polycomb targets and clust6 genes for synaptic/neural-function terms would tie the loop clustering to the transcriptional phenotype.

5. **Clust2 interior K27me3 enrichment.** clust2 (nominally "unchanged") shows significant K27me3 interior enrichment (Ext/Int = 0.952, p = 5.9e-6 in ctrl), though it is 92% unchanged by edgeR. This may indicate a pre-existing Polycomb signature at loops that have not yet reached the logFC threshold for differential calling but are biologically primed for remodeling.

### 3.2 Questions Requiring Additional Data or Computation

6. **Early timepoint clustering.** With only ~200 differential loops at P12, k-means is underpowered. However, a supervised approach (projecting P12 loops into the adult cluster space) could test whether the same two mechanisms are detectable in early development at lower magnitude.

7. **RAD21/cohesin ChIP-seq.** Popay normalizes ChIP signal to RAD21 at anchor peaks. Generating RAD21 ChIP-seq for BAP1-KO cerebellum would enable direct comparison with Popay's NIPBL-depletion results and test whether cohesin occupancy is reduced at clust6 anchors.

8. **Insulation score analysis.** ~~Dixon suggested Cooltools insulation scores to complement the loop-level analysis.~~ **DONE (Phase 11).** Insulation scores computed at 10kb/200kb diamond and tested for ext/int asymmetry. Clust5 anchors sit at strong TAD boundaries (p = 3.5e-54); clust6 anchors mark a sharp boundary-to-open transition (p = 2.9e-179).

9. **SNIPER subcompartment analysis.** Dixon suggested SNIPER for fine-grained subcompartment calling. This would test whether gained loops (clust5) specifically map to Polycomb-associated subcompartments (B2/B3).

10. **Distance-stratified anchor-vs-span analysis.** The current ChromHMM enrichment is computed across all loops per cluster regardless of size. Stratifying by loop length (e.g., <200 kb, 200-500 kb, 500 kb-1 Mb, >1 Mb) would test whether the anchor-vs-span pattern is consistent across length scales or driven by the longest/shortest loops.

### 3.3 Conceptual Questions

11. **Are clust5 loops truly cohesin-dependent?** The gained CTCF-CTCF loops within Polycomb domains could form through either cohesin-mediated extrusion (arrested within the Polycomb domain) or Polycomb-Polycomb self-association (phase separation). Acute cohesin degradation (as in Popay's NIPBL system) at clust5-specific loci would distinguish these, but is not available for this tissue.

12. **What drives the clust6 anchor-specific Polycomb gain?** The data shows K119ub elevation and Polycomb enrichment specifically at lost-loop anchors, but the mechanism by which Polycomb preferentially targets these CTCF sites (rather than spreading uniformly) is unclear. The CpG methylation data (Section 2.7) suggests DNA methylation may contribute, but the CpG island K119ub depletion (Section 2.8) argues against a simple K119ub -> methylation -> CTCF loss pathway.

13. **Is the progressive developmental expansion driven by Polycomb domain growth or cell-type changes?** The early-to-adult expansion from ~200 to ~3,000 differential loops could reflect progressive Polycomb spreading within the same cell type, or could reflect increasing representation of mature cell types (e.g., granule neurons) that have a different Polycomb landscape. Single-cell Hi-C or cell-type-specific BAP1 deletion would distinguish these.

---

## 4. File Inventory

### Analysis Scripts (`scripts/`)

| Script | Phase | Purpose |
|--------|-------|---------|
| `01_build_loop_count_file.py` | 1 | Merged loops -> Popay format + metadata sidecar |
| `02_build_mm10_gene_annotation.R` | 1 | mm10 TSS +/- 750 bp gene BED |
| `03_chromhmm_segmentation.sh` | 2 | BinarizeBed + LearnModel(k=12) |
| `04_clustering.py` | 3 | Elbow + Cluster 3.0 k-means + visualizations |
| `05_grouped_analyses.py` | 4 | 8 downstream sub-analyses orchestrator |
| `06_deeptools_metagene.py` | 5 | Anchor metagene (8 BigWigs x 6 clusters) |
| `07_cooltools_pileup.py` | 6 | Off-diagonal Hi-C pileup (HPC) |
| `08_summary_figures.py` | 7 | 3 composite lab-meeting figures |
| `09_oriented_anchor_metagene.py` | 8 | Strand-aware anchor metagene |
| `quantify_orientation_asymmetry.py` | 8 | Ext/Int Wilcoxon -> TSV |
| `visualize_orientation_asymmetry.py` | 8 | K27me3 dual-panel figure |
| `10_clust6_subgroup_asymmetry.py` | 9 | Clust6 short/long split + asymmetry |
| `11_histone_anchors_metagene.py` | 10 | Clean profile figure from Phase 5 matrix |
| `12_comprehensive_asymmetry.py` | 11 | H2AK119ub, H3K27ac, PC1, insulation asymmetry (HPC) |

### Key Output Files (`bap1_late/`)

| Path | What |
|------|------|
| `cluster3/k-6/data/combined-clusters.txt` | Canonical clustering (38,948 x 9) |
| `chromHMM/{anchor,span}.txt` | OverlapEnrichment fold-enrichment (12 states x 6 clusters) |
| `chromHMM/{anchor,span}.{png,pdf,svg,jpg}` | KEY result heatmaps |
| `figures/loop_size/` | Per-cluster loop size distributions |
| `figures/loop_classification/` | Structural/CRE/mixed/unclassified per cluster |
| `figures/chromHMM_anchor/` | ChromHMM proportions stacked bar |
| `figures/cluster_differential_status/` | Cluster x edgeR direction crosstab |
| `figures/ChIP_intersect/` | DiffBind + anchor ChIP signal |
| `figures/annotation/` | Per-cluster gene lists |
| `figures/deeptools/histone_anchors/` | Phase 5 combined heatmap |
| `figures/deeptools/oriented_anchors/` | Phase 8 strand-aware heatmap + asymmetry TSV + K27me3 figure |
| `figures/deeptools/clust6_subgroups/` | Phase 9 clust6 short/long asymmetry |
| `figures/deeptools/comprehensive_asymmetry/` | Phase 11 H2AK119ub, H3K27ac, PC1, insulation asymmetry |
| `figures/deeptools/comprehensive_asymmetry/bigwigs/` | PC1 + insulation BigWigs (computed from mcools) |
| `figures/summary_figures/{dashboard,mechanism,heatmap}/` | Phase 7 composites |
| `cooltools/obs_exp_contacts/` | Phase 6 pileup (ctrl + mut) |

### Run Logs (`docs/`)

| Log | Phase | Content |
|-----|-------|---------|
| `phase1.txt` | 1 | Data prep: count file + gene annotation |
| `phase2.txt` | 2 | ChromHMM BinarizeBed + LearnModel |
| `phase3.txt` | 3 | v1 diagnostic (degenerate cluster) |
| `phase3_v2.txt` | 3 | v2 canonical (6 well-populated clusters) |
| `phase4_test_4.4.txt` | 4 | KEY result smoke test |
| `phase4.txt` | 4 | Full 8-sub-analysis run |
| `phase5.txt` | 5 | deepTools metagene (~96 min) |
| `phase8_summary.txt` | 7 | Summary figures + K119ub query |
| `oriented_metagene.txt` | 8 | Oriented anchor metagene (~1.7 h) |

### HPC Logs (`logs/`)

| Log | Phase | Content |
|-----|-------|---------|
| `comprehensive_asymmetry_48713858.out` | 11 | H2AK119ub, H3K27ac, PC1, insulation asymmetry (~91 min) |
