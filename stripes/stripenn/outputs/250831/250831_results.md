# Early Timepoint (250831/P12) Differential Stripe Analysis — Stripenn

**Dataset**: BAP1-KO vs Wildtype Mouse Cerebellum (mm10)
**Timepoint**: Early (250831 / P12)
**Pipeline**: Stripenn v1.1.65 (Canny edge detection + pixel saturation)
**Comparison**: Mutant (BAP1-KO) vs Control (Wildtype), n=3 per condition
**Resolutions**: 5kb (primary) and 10kb (validation)

---

## 1. Summary Statistics

### Overall Counts

| Metric | 5kb | 10kb |
|--------|-----|------|
| Total union stripes | 4008 | 1879 |
| Significant (FDR<0.05) | 96 (2.4%) | 53 (2.8%) |
| Significant (FDR<0.10) | 224 (5.6%) | 162 (8.6%) |
| Lost (control-only) | 949 | 483 |
| Gained (mutant-only) | 776 | 401 |
| Unchanged (shared) | 2283 | 995 |

### Confidence Tiers (lost/gained stripes, 5kb)

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| High | 12 | 10 | 22 |
| Medium | 0 | 0 | 0 |
| Low | 937 | 766 | 1703 |

### Confidence Tiers (lost/gained stripes, 10kb)

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| High | 14 | 10 | 24 |
| Medium | 0 | 0 | 0 |
| Low | 469 | 391 | 860 |

### Cross-Resolution Concordance

| Metric | Value |
|--------|-------|
| Matched across resolutions | 1404 |
| 5kb-only stripes | 2604 |
| 10kb-only stripes | 475 |
| Both significant (FDR<0.05) | 22 / 1404 |
| Direction concordant (both sig) | 15 / 22 |
| logFC correlation (Pearson r) | 0.808 |
| Both concordant (direction-matched, any significance) | 759 |
| Both discordant | 645 |

### Effect Size Distribution (5kb, significant FDR<0.05 stripes only)

| Category | Count |
|----------|-------|
| Strong (\|logFC\| > 1.0) | 0 |
| Moderate (\|logFC\| 0.5–1.0) | 0 |
| Weak (\|logFC\| 0.3–0.5) | 0 |
| Minimal (\|logFC\| < 0.3) | 96 |

### edgeR Model Summary (5kb)

| Parameter | Value |
|-----------|-------|
| Common BCV | 0.020 |
| Median tagwise BCV | 0.018 |
| TMM normalization factor range | 0.9971 – 1.0028 |
| Median logFC (all stripes) | -0.003 |
| Median logFC (FDR<0.05 significant) | 0.104 |
| Median \|logFC\| (FDR<0.05 significant) | 0.117 |
| Directional consistency (lost/gained) | 40.9% |

### edgeR Model Summary (10kb)

| Parameter | Value |
|-----------|-------|
| Common BCV | 0.021 |
| Median tagwise BCV | 0.017 |
| TMM normalization factor range | 0.9979 – 1.0034 |
| Median logFC (all stripes) | -0.001 |
| Median logFC (FDR<0.05 significant) | 0.089 |
| Directional consistency (lost/gained) | 39.9% |

### Stripe Geometry (5kb, all union stripes)

| Metric | Lost | Gained | Unchanged |
|--------|------|--------|-----------|
| Median length (kb) | 385 | 377.5 | 450 |
| Mean length (kb) | 489.3 | 450.2 | — |
| Median width (kb) | 25 | 25 | 25 |
| N (3prime) | 477 | 377 | 1043 |
| N (5prime) | 472 | 399 | 1240 |

### ChIP-seq Anchor Annotation (5kb, all union stripes)

Anchor type distribution across all 4008 stripes, with breakdown by direction:

| Anchor Type | Lost (n=949) | Lost% | Gained (n=776) | Gained% | Unchanged (n=2283) | Unchanged% | All (n=4008) | All% |
|---|---|---|---|---|---|---|---|---|
| Poised_Enhancer | 319 | 33.6% | 272 | 35.1% | 821 | 36.0% | 1412 | 35.2% |
| Other | 323 | 34.0% | 252 | 32.5% | 643 | 28.2% | 1218 | 30.4% |
| Active_Promoter | 118 | 12.4% | 93 | 12.0% | 296 | 13.0% | 507 | 12.6% |
| Active_Enhancer | 102 | 10.7% | 90 | 11.6% | 287 | 12.6% | 479 | 12.0% |
| Repressed_Promoter | 47 | 5.0% | 42 | 5.4% | 124 | 5.4% | 213 | 5.3% |
| Bivalent_Promoter | 27 | 2.8% | 22 | 2.8% | 95 | 4.2% | 144 | 3.6% |
| Polycomb | 13 | 1.4% | 5 | 0.6% | 17 | 0.7% | 35 | 0.9% |

**ChIP mark overlap rates by direction (5kb, all union stripes):**

| ChIP Mark | Lost (n=949) | Lost% | Gained (n=776) | Gained% | Unchanged (n=2283) | Unchanged% |
|---|---|---|---|---|---|---|
| H3K27ac | 277 | 29.2% | 226 | 29.1% | 715 | 31.3% |
| H3K27me3 | 159 | 16.8% | 118 | 15.2% | 404 | 17.7% |
| H3K4me1 | 708 | 74.6% | 601 | 77.4% | 1849 | 81.0% |
| H3K4me3 | 185 | 19.5% | 159 | 20.5% | 483 | 21.2% |
| Bivalent (H3K4me3+H3K27me3) | 30 | 3.2% | 26 | 3.4% | 103 | 4.5% |

---

## 2. High-Priority Stripes for Validation

Top 10 stripes filtered to: direction IN (lost, gained) AND direction\_confidence == "high", sorted by FDR ascending. Gene name and anchor type now populated from visualization annotation.

| stripe\_id | chr | Anchor (pos1–pos2) | Extent (pos3–pos4) | Direction | dir\_type | logFC | FDR | Confidence | Resolution support | in\_10kb | Nearest gene | Anchor type |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| stripe\_2752 | chr4 | 125,490,001–125,510,000 | 125,490,001–125,970,000 | gained | 3prime | +0.1396 | 0.0155 | high | 5kb\_only | No | Grik3 | Bivalent\_Promoter |
| stripe\_3157 | chr6 | 51,845,001–51,865,000 | 51,465,001–51,865,000 | lost | 5prime | −0.1250 | 0.0173 | high | 5kb\_only | No | Skap2 | Active\_Enhancer |
| stripe\_2579 | chr4 | 3,595,001–3,620,000 | 3,595,001–3,935,000 | lost | 3prime | −0.1056 | 0.0290 | high | 5kb\_only | No | Tgs1 | Poised\_Enhancer |
| stripe\_3636 | chr8 | 64,010,001–64,040,000 | 64,010,001–64,740,000 | lost | 3prime | −0.1041 | 0.0307 | high | both\_concordant | Yes | Sgo2b | Poised\_Enhancer |
| stripe\_1504 | chr6 | 43,360,001–43,410,000 | 43,360,001–45,650,000 | lost | 3prime | −0.0876 | 0.0348 | high | 10kb\_only | Yes | Mir365-1 / Mrtfb | Bivalent\_Promoter |
| stripe\_1554 | chr6 | 111,040,001–111,110,000 | 111,040,001–112,190,000 | lost | 3prime | −0.0833 | 0.0366 | high | 10kb\_only | Yes | Arhgap31 | Other |
| stripe\_3805 | chr9 | 41,920,001–41,955,000 | 41,390,001–41,955,000 | gained | 5prime | +0.1046 | 0.0381 | high | 5kb\_only | No | Sorl1 | Poised\_Enhancer |
| stripe\_2321 | chr2 | 151,250,001–151,270,000 | 151,250,001–151,465,000 | gained | 3prime | +0.2112 | 0.0471 | high | 5kb\_only | No | 4930442J19Rik | Other |
| stripe\_2170 | chr2 | 70,545,001–70,570,000 | 70,120,001–70,570,000 | gained | 5prime | +0.0904 | 0.0481 | high | 5kb\_only | No | Gad1os / Gad1 | Bivalent\_Promoter |
| stripe\_3109 | chr6 | 21,215,001–21,235,000 | 21,215,001–21,705,000 | lost | 3prime | −0.1151 | 0.0494 | high | 5kb\_only | No | Kcnd2 | Active\_Promoter |

**Biological relevance notes:**
- **stripe\_2752 (Grik3, Bivalent\_Promoter):** Glutamate ionotropic receptor kainic acid type subunit 3. Bivalent promoter indicates H3K4me3+H3K27me3 co-occupancy — a Polycomb-regulated locus poised for activation. Gained in mutant.
- **stripe\_2170 (Gad1, Bivalent\_Promoter):** Glutamate decarboxylase 1 (GABAergic marker). Bivalent anchor; the nearest coding gene is Gad1 directly at the anchor. Gained in mutant.
- **stripe\_1504 (Mrtfb, Bivalent\_Promoter):** Myocardin-related transcription factor B. Bivalent promoter; the anchor overlaps Mir365-1 with Mrtfb 178kb away. Lost in mutant.
- **stripe\_3636 (Sgo2b, Poised\_Enhancer):** The only high-confidence stripe supported at both resolutions with concordant direction. Most reliable call in this timepoint. Lost in mutant.
- **stripe\_3109 (Kcnd2, Active\_Promoter):** Potassium voltage-gated channel subfamily D member 2. Active promoter anchor (H3K27ac+H3K4me1). Lost in mutant.

**Note on stripe\_1504 (Stripiness ctrl = −0.098):** Negative Stripiness score indicates Stripenn found no net stripe signal relative to background at this locus; the edgeR significance is driven by lower quantification scores in the mutant, not a classically stripe-shaped region in either condition. Treat with extra caution.

**Note on stripe\_2321 (Stripiness mut = −0.067):** Same caveat — absent or weak stripe signal in the mutant despite the "gained" call. The mutant has marginally higher contact counts by edgeR, but the Stripiness score does not confirm a stripe structure.

---

## 3. JuiceBox Coordinate Blocks

Each coordinate block uses 50 kb padding: `chr:min(pos1,pos3)−50000 – max(pos2,pos4)+50000`.

```
# LOST stripes — open CONTROL .hic first; confirm stripe is present, then switch to MUTANT to confirm loss
# Sorted by FDR ascending

stripe_3157  chr6:51415001-51915000    lost  5prime  logFC=-0.1250  FDR=0.0173  5kb_only         Skap2 / Active_Enhancer
stripe_2579  chr4:3545001-3985000      lost  3prime  logFC=-0.1056  FDR=0.0290  5kb_only         Tgs1 / Poised_Enhancer
stripe_3636  chr8:63960001-64790000    lost  3prime  logFC=-0.1041  FDR=0.0307  both_concordant  Sgo2b / Poised_Enhancer  [PRIORITY]
stripe_1504  chr6:43310001-45700000    lost  3prime  logFC=-0.0876  FDR=0.0348  10kb_only        Mrtfb / Bivalent_Promoter
stripe_1554  chr6:110990001-112240000  lost  3prime  logFC=-0.0833  FDR=0.0366  10kb_only        Arhgap31 / Other
stripe_3109  chr6:21165001-21755000    lost  3prime  logFC=-0.1151  FDR=0.0494  5kb_only         Kcnd2 / Active_Promoter

# GAINED stripes — open MUTANT .hic first; confirm stripe is present, then switch to CONTROL to confirm absence
# Sorted by FDR ascending

stripe_2752  chr4:125440001-126020000   gained  3prime  logFC=+0.1396  FDR=0.0155  5kb_only  Grik3 / Bivalent_Promoter  [PRIORITY]
stripe_3805  chr9:41340001-42005000     gained  5prime  logFC=+0.1046  FDR=0.0381  5kb_only  Sorl1 / Poised_Enhancer
stripe_2321  chr2:151200001-151515000   gained  3prime  logFC=+0.2112  FDR=0.0471  5kb_only  4930442J19Rik / Other    [highest |logFC|]
stripe_2170  chr2:70070001-70620000     gained  5prime  logFC=+0.0904  FDR=0.0481  5kb_only  Gad1 / Bivalent_Promoter
```

**Priority for manual inspection:** stripe\_3636 (chr8:63960001-64790000) is the only high-confidence stripe supported at both resolutions with concordant direction. This is the most reliable call in this timepoint. Among gained stripes, stripe\_2752 (Grik3, chr4:125440001-126020000) has the lowest FDR and a biologically relevant Bivalent\_Promoter anchor.

---

## 4. Biological Context

This is a factual summary of what the data shows at this timepoint.

- **Weak overall differential signal:** Only 96 stripes (2.4%) are significant at FDR<0.05 at 5kb; 53 (2.8%) at 10kb. These are the lowest significance rates compared to the late timepoint.
- **Direction bias toward loss:** 949 lost vs 776 gained at 5kb (ratio 1.22:1); 483 vs 401 at 10kb (ratio 1.20:1). More stripes are called control-only than mutant-only at this timepoint.
- **All effect sizes are minimal:** No significant stripe reaches the weak threshold (|logFC| > 0.3). The median logFC among significant stripes is 0.104 at 5kb. The maximum logFC among the top 10 high-confidence stripes is +0.211 (stripe\_2321).
- **Very low BCV:** Common BCV = 0.020 at 5kb and 0.021 at 10kb. This is substantially below typical bulk RNA-seq values (0.2–0.5), indicating very low between-replicate variability relative to mean count levels.
- **Low directional consistency:** 40.9% of lost/gained stripes at 5kb (39.9% at 10kb) show direction-consistent edgeR calls. This means that for the majority of detection-based lost/gained calls, the edgeR logFC sign does not reinforce the direction assigned by source (control-only / mutant-only).
- **Sparse cross-resolution support:** Only 22 stripes are significant at both resolutions simultaneously; 15 of those are direction-concordant. Only 1 stripe in the top-10 high-confidence list (stripe\_3636) is both\_concordant across resolutions.
- **No medium-confidence calls:** The medium confidence tier is empty at both resolutions — all lost/gained stripes are either high (22 at 5kb, 24 at 10kb) or low confidence.
- **Six lost stripes on chr6:** Of the 10 top high-confidence stripes, 4 are on chr6 (stripe\_3157, stripe\_1504, stripe\_1554, stripe\_3109). This concentration may reflect a regional chromatin feature or could be noise; the enrichment data (Section 6) does not clarify this.

---

## 5. Anchor Annotation Analysis

### Anchor Type Enrichment (Lost vs Gained vs Background)

Anchor type distributions across all 4008 stripes at 5kb. Background is the "unchanged" set (n=2283).

| Anchor Type | Lost% | Gained% | Unchanged% | Lost vs Background | Gained vs Background |
|---|---|---|---|---|---|
| Poised_Enhancer | 33.6% | 35.1% | 36.0% | Under-represented | Slightly under |
| Other | 34.0% | 32.5% | 28.2% | Over-represented | Over-represented |
| Active_Promoter | 12.4% | 12.0% | 13.0% | Similar | Similar |
| Active_Enhancer | 10.7% | 11.6% | 12.6% | Under-represented | Under-represented |
| Repressed_Promoter | 5.0% | 5.4% | 5.4% | Similar | Similar |
| Bivalent_Promoter | 2.8% | 2.8% | 4.2% | Under-represented | Under-represented |
| Polycomb | 1.4% | 0.6% | 0.7% | Over-represented (2×) | Similar |

Key observations:
- **Polycomb anchors are modestly enriched in lost stripes** (1.4% lost vs 0.7% unchanged). Absolute counts are very small (13 lost Polycomb vs 17 unchanged), so this cannot be interpreted as statistically robust.
- **Bivalent_Promoter anchors are depleted in both lost and gained** relative to background (2.8% vs 4.2%). This may reflect that bivalent loci are more stably organized at P12, or the small total count (n=144) makes this ratio unstable.
- **No strong directional enrichment of any anchor type** is apparent in this timepoint. The lost and gained distributions are highly similar to each other and to background.
- **Enrichment in high-confidence significant stripes:** Of the 12 high-confidence lost stripes, anchor types are Bivalent\_Promoter (2), Poised\_Enhancer (2), Active\_Promoter (1), Active\_Enhancer (1), Other (1); of the 10 high-confidence gained, Bivalent\_Promoter (2), Poised\_Enhancer (2), Other (1). Given these are only 22 stripes total, no statistical inference on enrichment is appropriate.

### ChIP-seq Mark Overlap Rates by Direction

| ChIP Mark | Lost (n=949) | Lost% | Gained (n=776) | Gained% | Unchanged (n=2283) | Unchanged% | Lost vs Unchanged | Gained vs Unchanged |
|---|---|---|---|---|---|---|---|---|
| H3K27ac | 277 | 29.2% | 226 | 29.1% | 715 | 31.3% | −2.1 pp | −2.2 pp |
| H3K27me3 | 159 | 16.8% | 118 | 15.2% | 404 | 17.7% | −0.9 pp | −2.5 pp |
| H3K4me1 | 708 | 74.6% | 601 | 77.4% | 1849 | 81.0% | −6.4 pp | −3.6 pp |
| H3K4me3 | 185 | 19.5% | 159 | 20.5% | 483 | 21.2% | −1.7 pp | −0.7 pp |
| Bivalent (H3K4me3+H3K27me3) | 30 | 3.2% | 26 | 3.4% | 103 | 4.5% | −1.3 pp | −1.1 pp |

(pp = percentage points relative to unchanged background)

Observations:
- **H3K4me1 shows the largest differential** relative to background: lost and gained stripes both have lower H3K4me1 overlap rates than unchanged stripes (−6.4 pp for lost, −3.6 pp for gained). H3K4me1 marks enhancers and poised regulatory elements, and its lower rate in directional stripes may reflect that differential stripes are somewhat less enriched at active enhancer regions.
- **No mark shows a directional difference between lost and gained** that is large in magnitude. The largest difference between lost and gained is for H3K27me3 (16.8% vs 15.2%, 1.6 pp) and H3K4me1 (74.6% vs 77.4%, 2.8 pp). These are likely not biologically meaningful given the overall very weak differential signal.
- **H3K27me3 overlap rates are similar across all groups** (~15–18%), suggesting that Polycomb-marked regions are not selectively depleted or gained in this timepoint.

### Stripiness Score Comparison (Shared Stripes)

From the visualization log: 2280 shared stripes with both ctrl and mutant Stripiness scores are available. Full correlation data is in `outputs/250831/visualizations/stripiness_250831.*` (visual output only; TSV not separately saved). The edgeR analysis is performed on `O_Sum_added` counts, not Stripiness scores directly. Stripes with large difference in Stripiness between conditions are not reported here as a separate table.

---

## 6. Pathway Enrichment Summary

Enrichment was performed using genes nearest to the stripe anchors in the high-confidence significant set. From the visualization log: early timepoint used **26 genes (lost direction)** and **19 genes (gained direction)** as input, with a background of 4344 genes (all stripe anchor-associated genes).

The gained direction produced **no significant GO or KEGG terms** at any threshold. All enrichment results below are for the lost direction only.

### GO Biological Process (Lost direction, early/P12)

49 significant terms for lost stripes (adj. p < 0.05). Top 5 by adjusted p-value:

| Rank | GO Term | Description | Gene Ratio | Fold Enrichment | Adj. p-value | Genes |
|---|---|---|---|---|---|---|
| 1 | GO:0007369 | gastrulation | 5/23 | 14.5× | 0.0161 | Inhba, Wnt5a, Prickle1, Klf4, Tenm4 |
| 2 | GO:0001702 | gastrulation with mouth forming second | 3/23 | 44.2× | 0.0161 | Wnt5a, Prickle1, Tenm4 |
| 2 | GO:2000648 | positive regulation of stem cell proliferation | 3/23 | 44.2× | 0.0161 | Ptprc, Wnt5a, Fgfr2 |
| 4 | GO:0098773 | skin epidermis development | 4/23 | 19.1× | 0.0163 | Inhba, Wnt5a, Klf4, Fgfr2 |
| 5 | GO:0097190 | apoptotic signaling pathway | 6/23 | 7.97× | 0.0201 | Ptprc, Inhba, Wnt5a, Zfp622, Klf4, Fgfr2 |

Additional significant terms (selected):
- GO:0043408 regulation of MAPK cascade (6/23, 5.85×, adj.p=0.0366) — Ptprc, Inhba, Wnt5a, Zfp622, Klf4, Fgfr2
- GO:0060759 regulation of response to cytokine stimulus (3/23, 17.1×, adj.p=0.0366) — Ptprc, Wnt5a, Klf4
- GO:0021879 forebrain neuron differentiation (3/23, 22.1×, adj.p=0.0354) — Inhba, Wnt5a, Fgfr2
- GO:0060070 canonical Wnt signaling pathway (4/23, 8.21×, adj.p=0.0427) — Wnt5a, Prickle1, Klf4, Fgfr2
- GO:0050877 nervous system process (8/23, 4.46×, adj.p=0.0296) — Gria1, Kcnd2, Tafa4, Tenm4, and others

**Key driver genes** (appearing in multiple terms): Wnt5a (non-canonical Wnt signaling), Fgfr2 (fibroblast growth factor receptor), Klf4 (Krüppel-like factor 4, stem cell transcription factor), Ptprc (CD45, immune cell marker), Inhba (activin subunit).

Note: The gene set (23 unique genes after background filtering from 26 input) is very small, so fold enrichments are high but confidence intervals are wide. The most striking terms — gastrulation, stem cell proliferation, skin epidermis development — likely reflect that these anchors are near developmentally regulated loci rather than indicating a specific biological process in P12 cerebellum.

### GO Cellular Component (Lost or Gained direction, early/P12)

**No significant GO CC terms** were found for either direction (confirmed in visualization log: "No significant GO CC terms").

### GO Molecular Function (Lost direction, early/P12)

1 significant term:

| GO Term | Description | Gene Ratio | Fold Enrichment | Adj. p-value | Genes |
|---|---|---|---|---|---|
| GO:0004984 | olfactory receptor activity | 3/24 | 33.4× | 0.0164 | Or51s1, Or51h5, Or51a8 |

All three olfactory receptor genes (Or51s1, Or51h5, Or51a8) are members of a cluster located near one another on chromosome 7. This enrichment reflects co-location of three olfactory receptor genes near a single stripe anchor region rather than a functional biology of olfaction in cerebellum.

### KEGG Pathways (Lost direction, early/P12)

2 significant pathways:

| Pathway | Description | Category | Gene Ratio | Fold Enrichment | Adj. p-value | Genes |
|---|---|---|---|---|---|---|
| mmu04550 | Signaling pathways regulating pluripotency of stem cells | Cellular Processes | 4/14 | 12.9× | 0.0094 | Inhba, Wnt5a, Klf4, Fgfr2 |
| mmu04740 | Olfactory transduction | Sensory system | 3/14 | 13.0× | 0.0349 | Or51s1, Or51h5, Or51a8 |

**No significant KEGG pathways** were found for the gained direction.

### Enrichment Summary and Caveats

- The gained direction is completely devoid of enrichment, consistent with the small gene set (19 genes) and the overall weak signal at this timepoint.
- The lost direction shows enrichment for **stem cell pluripotency signaling** (Wnt5a, Fgfr2, Klf4, Inhba) and **developmental patterning** (gastrulation, morphogenesis). Given the P12 cerebellar context, these likely reflect that the BAP1-KO disrupts chromatin organization at developmentally regulated loci that are still being resolved at this early postnatal stage.
- The olfactory receptor cluster terms (GO:0004984, mmu04740) are an artifact of three Or51-family genes near a single stripe anchor and should not be interpreted as cerebellum-specific biology.
- **No Polycomb-related GO terms** (e.g., PRC2, H2AK119ub1, chromatin silencing) appear in either direction. This is consistent with the overall weak signal: the 26 lost-stripe anchor genes do not specifically cluster near known Polycomb targets.
- No chromatin-remodeling or DNA repair GO terms appear, despite BAP1's known role in the PR-DUB complex.

---

## 7. Cross-Timepoint Comparison

### Significance Rates

| Metric | Early (250831/P12) | Late (250402/Adult) |
|--------|-------------------|---------------------|
| Union stripes (5kb) | 4008 | 7371 |
| Significant FDR<0.05 (5kb) | 96 (2.4%) | 2320 (31.5%) |
| High-confidence lost | 12 | 367 |
| High-confidence gained | 10 | 638 |
| Directional consistency | 40.9% | 60.7% |
| Common BCV | 0.020 | 0.012 |
| logFC correlation (5kb vs 10kb) | 0.808 | 0.850 |

The late timepoint shows 13× more significant stripes by count and 31.5% vs 2.4% by rate. This is a substantial difference. Both BCV values are very low, suggesting the difference is not driven by technical quality but by genuine differences in the magnitude and prevalence of differential stripe signal.

### Directional Bias

| Direction | Early (5kb) | Early% | Late (5kb) | Late% |
|-----------|-------------|--------|------------|-------|
| Lost (control-only) | 949 | 23.7% | 1528 | 20.7% |
| Gained (mutant-only) | 776 | 19.4% | 2052 | 27.8% |
| Unchanged (shared) | 2283 | 57.0% | 3784 | 51.3% |
| Lost:Gained ratio | — | 1.22:1 | — | 0.74:1 |

**The directional bias reverses between timepoints.** At P12 (early), more stripes are lost than gained (1.22:1). At adult (late), more stripes are gained than lost (0.74:1 lost:gained, i.e., gained outnumber lost). In absolute numbers: early has 173 more lost than gained; late has 524 more gained than lost.

### Anchor Type Comparison

| Anchor Type | Early All% | Late All% | Early Lost% | Late Lost% | Early Gained% | Late Gained% |
|---|---|---|---|---|---|---|
| Poised_Enhancer | 35.2% | 32.3% | 33.6% | 32.2% | 35.1% | 31.8% |
| Other | 30.4% | 26.1% | 34.0% | 25.3% | 32.5% | 26.8% |
| Active_Promoter | 12.6% | 17.4% | 12.4% | 17.3% | 12.0% | 17.3% |
| Active_Enhancer | 12.0% | 16.6% | 10.7% | 16.0% | 11.6% | 17.5% |
| Repressed_Promoter | 5.3% | 5.0% | 5.0% | 5.1% | 5.4% | 4.7% |
| Bivalent_Promoter | 3.6% | 1.2% | 2.8% | 1.6% | 2.8% | 1.6% |
| Polycomb | 0.9% | 1.5% | 1.4% | 1.9% | 0.6% | 1.4% |

Key cross-timepoint differences:
- **Active_Promoter and Active_Enhancer are more prevalent in the late timepoint** (+4.8 pp and +4.6 pp overall). This may reflect greater chromatin accessibility and gene activity at promoter and enhancer regions in adult cerebellum compared to P12.
- **Bivalent_Promoter is more prevalent in the early timepoint** (3.6% vs 1.2%), consistent with developmental biology: bivalent chromatin (H3K4me3+H3K27me3) is a hallmark of developmentally poised genes that resolve to active or repressed states during differentiation.
- **Polycomb anchors are more prevalent in the late timepoint** (1.5% vs 0.9%), which may reflect consolidation of repressive chromatin domains in adult tissue. In the late timepoint, Polycomb anchors are enriched in lost stripes (1.9%).

### Enrichment Concordance

Early (P12): Enrichment only in lost direction; driven by Wnt5a, Fgfr2, Klf4, Inhba (developmental signaling, stem cell pluripotency).

Late (Adult): Separate results document covers enrichment for the late timepoint. The late timepoint has 444 lost and 639 gained anchor genes, providing substantially more power for enrichment analysis.

No direct comparison of enrichment terms is possible without reading the late results, but the gene sets differ substantially in size (26 vs 444 for lost direction) and likely in biological content.

### Interpretation

The large difference in significance rates (2.4% vs 31.5%) is consistent with two scenarios that cannot be distinguished from this data alone: (1) BAP1-KO effects on chromatin stripe organization accumulate progressively and are minimal at P12 but substantial in adult; (2) the P12 Hi-C data has lower signal-to-noise for stripe detection (e.g., due to developmental heterogeneity in the P12 cerebellum). The reversal of directional bias (more lost at P12, more gained in adult) is a distinct observation that warrants specific investigation: it could indicate that the early loss of some stripes precedes the later formation of new ones, or that different biological processes dominate at each developmental stage.

---

## 8. Caveats and Confidence Assessment

- **All effect sizes are minimal:** No significant stripe reaches |logFC| > 0.3. The strongest call (stripe\_2321, logFC = +0.211) is below the weak threshold. Fold changes of this magnitude are at the limit of biological interpretability for chromatin contact data.
- **Low directional consistency (40.9%):** More than half of detection-based lost/gained calls have an edgeR logFC pointing in the opposite direction. This is consistent with the calls being dominated by detection noise rather than true differential signal.
- **Only 22 high-confidence significant stripes:** This is insufficient for robust enrichment analysis — any category-level statistics have very wide confidence intervals. The GO/KEGG enrichment (Section 6) uses only 26 lost-direction and 19 gained-direction genes; gained shows no enrichment at all.
- **Enrichment results driven by small gene clusters:** The two KEGG pathways (stem cell pluripotency, olfactory transduction) and the top GO MF term (olfactory receptor activity) are driven by 3–4 genes each. The olfactory receptor result (Or51s1/h5/a8) is almost certainly an artifact of gene cluster proximity to a stripe anchor.
- **Substantially weaker signal than late timepoint:** 31.5% vs 2.4% significant rate. The early timepoint may reflect a developmental stage where BAP1-KO effects on chromatin stripes are not yet established, or the differential effect is below the detection threshold.
- **Very low BCV (0.020):** While low BCV improves nominal power, it also makes the model sensitive to small systematic biases. The near-zero common dispersion warrants caution in interpreting any p-values.
- **Negative Stripiness scores:** stripe\_1504 (Stripiness ctrl = −0.098) and stripe\_2321 (Stripiness mut = −0.067) have negative Stripiness values in the condition they are called in, indicating Stripenn does not detect a classical stripe pattern at these loci. These calls should not be used as visual validation targets without first confirming the stripe structure visually in JuiceBox.
- **No medium-confidence tier:** The bimodal confidence distribution (high vs low, nothing in between) reflects the strict thresholds used for tier assignment. The 22 high-confidence calls are robustly classified; the 1703 low-confidence calls should be treated as exploratory only.
- **ChIP mark differences between lost/gained/unchanged are small** (largest: −6.4 pp for H3K4me1 in lost vs unchanged), and no directional mark differential exceeds 3 pp between lost and gained. These are not indicative of a specific chromatin state bias in differential stripes at this timepoint.
- **Anchor type enrichments do not survive careful interpretation:** Polycomb anchors are 2× enriched in lost stripes by percentage (1.4% vs 0.7%), but this represents 13 lost Polycomb stripes vs 17 unchanged — far too few for a meaningful enrichment test. No formal enrichment testing was performed.
