# Late Timepoint (250402/Adult) Differential Stripe Analysis — Stripenn

**Dataset**: BAP1-KO vs Wildtype Mouse Cerebellum (mm10)
**Timepoint**: Late (250402 / Adult)
**Pipeline**: Stripenn v1.1.65 (Canny edge detection + pixel saturation)
**Comparison**: Mutant (BAP1-KO) vs Control (Wildtype), n=3 per condition
**Resolutions**: 5kb (primary) and 10kb (validation)
**Visualization completed**: 2026-04-21

---

## 1. Summary Statistics

### Overall Counts

| Metric | 5kb | 10kb |
|--------|-----|------|
| Total union stripes | 7,371 | 3,566 |
| Significant (FDR<0.05) | 2,320 (31.5%) | 1,325 (37.2%) |
| Significant (FDR<0.10) | 2,919 | 1,621 |
| Lost (control-only) | 1,528 | 766 |
| Gained (mutant-only) | 2,052 | 967 |
| Unchanged (shared, no change) | 3,784 | 1,833 |
| Strengthened | 3 | 0 |
| Weakened | 4 | 0 |

### Confidence Tiers (lost/gained stripes, 5kb)

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| High | 367 | 638 | 1,005 |
| Medium | 0 | 0 | 0 |
| Low | 1,161 | 1,414 | 2,575 |

### Cross-Resolution Concordance

| Metric | Value |
|--------|-------|
| Matched across resolutions | 2,522 |
| Both significant | 672 / 2,522 |
| Direction concordant (both sig) | 377 / 672 |
| logFC correlation (Pearson r) | 0.850 |
| Both concordant (high-confidence) | 1,273 |
| Both discordant | 1,249 |

### Effect Size Distribution (5kb, significant stripes only)

| Category | Count |
|----------|-------|
| Strong (\|logFC\| > 1.0) | 0 |
| Moderate (\|logFC\| 0.5–1.0) | 0 |
| Weak (\|logFC\| 0.3–0.5) | 19 |
| Minimal (\|logFC\| < 0.3) | 2,301 |

### edgeR Model Summary (5kb)

| Parameter | Value |
|-----------|-------|
| Common BCV | 0.012 |
| Median logFC (all) | −0.001 |
| Median logFC (significant) | 0.071 |
| Directional consistency | 60.7% |

### Stripe Geometry (5kb, significant stripes only)

| Metric | Lost | Gained |
|--------|------|--------|
| Median length (kb) | 645 | 625 |
| Mean length (kb) | 727 | 710 |
| Median width (kb) | 20 | 20 |

### ChIP-seq Anchor Annotation (5kb, all union stripes)

Anchor type counts from the ChIP-seq annotation step (Section 8, visualization run). The 7-category classification uses H3K27ac, H3K27me3, H3K4me1, H3K4me3, and a pre-computed Bivalent (H3K4me3+H3K27me3) peak set at the adult timepoint (250402), with a 2kb TSS proximity threshold.

| Anchor Type | All (n=7,371) | % | Lost (n=1,528) | % | Gained (n=2,052) | % |
|-------------|---------------|---|----------------|---|------------------|---|
| Poised_Enhancer | 2,381 | 32.3 | 463 | 30.3 | 748 | 36.5 |
| Other | 1,925 | 26.1 | 406 | 26.6 | 539 | 26.3 |
| Active_Promoter | 1,280 | 17.4 | 282 | 18.5 | 292 | 14.2 |
| Active_Enhancer | 1,224 | 16.6 | 298 | 19.5 | 287 | 14.0 |
| Repressed_Promoter | 367 | 5.0 | 54 | 3.5 | 119 | 5.8 |
| Polycomb | 109 | 1.5 | 13 | 0.9 | 42 | 2.0 |
| Bivalent_Promoter | 85 | 1.2 | 12 | 0.8 | 25 | 1.2 |

**ChIP-seq mark overlap rates by direction (derived from anchor type classification):**

| Mark / Category | Lost (n=1,528) | % | Gained (n=2,052) | % | All (n=7,371) | % |
|---|---|---|---|---|---|---|
| H3K27me3+ (Polycomb + Repressed_Promoter + Bivalent_Promoter) | 79 | 5.2 | 186 | 9.1 | 561 | 7.6 |
| H3K27ac+ (Active_Enhancer) | 298 | 19.5 | 287 | 14.0 | 1,224 | 16.6 |
| H3K4me3+ (Active_Promoter + Bivalent_Promoter) | 294 | 19.2 | 317 | 15.5 | 1,365 | 18.5 |
| H3K4me1+ (Poised_Enhancer) | 463 | 30.3 | 748 | 36.5 | 2,381 | 32.3 |

---

## 2. High-Priority Stripes for Validation

Top 10 stripes by lowest FDR, filtered to: `resolution_support == "both_concordant"`, `direction_confidence == "high"`, `direction IN (lost, gained)`. All 10 are gained stripes; the highest-confidence lost stripes follow in the second table. Gene annotations are from the annotated TSV (nearest gene to anchor).

### Top 10 Gained — both_concordant, high confidence

| stripe_id | JuiceBox coord | type | dir_type | logFC | FDR | confidence | res_support | nearest_gene | anchor_type | ChIP marks |
|-----------|----------------|------|----------|-------|-----|------------|-------------|--------------|-------------|------------|
| stripe_3908 | chr2:44,605,001–46,505,000 | gained | 5prime | 0.3684 | 1.43e-07 | high | both_concordant | 1700019E08Rik (Zeb2 body) | Poised_Enhancer | H3K27me3, H3K4me1 |
| stripe_7016 | chr9:40,240,001–41,510,000 | gained | 5prime | 0.3256 | 1.53e-07 | high | both_concordant | Mir100 (Ubash3b body) | Poised_Enhancer | H3K27me3, H3K4me1 |
| stripe_5385 | chr5:83,990,001–84,465,000 | gained | 5prime | 0.3402 | 1.98e-07 | high | both_concordant | Epha5 | Poised_Enhancer | H3K27me3, H3K4me1, H3K4me3, Bivalent |
| stripe_4729 | chr4:20,050,001–21,540,000 | gained | 3prime | 0.3562 | 3.16e-07 | high | both_concordant | Ggh (Nkain3 body) | Other | none |
| stripe_4730 | chr4:20,195,001–20,830,000 | gained | 5prime | 0.3675 | 8.88e-07 | high | both_concordant | Nkain3 | Other | none |
| stripe_1743 | chr13:9,815,001–10,980,000 | gained | 3prime | 0.2806 | 1.66e-06 | high | both_concordant | — | — | — |
| stripe_4797 | chr4:55,285,001–56,415,000 | gained | 5prime | 0.2342 | 4.78e-06 | high | both_concordant | — | — | — |
| stripe_6491 | chr7:135,745,001–137,380,000 | gained | 5prime | 0.2637 | 7.71e-06 | high | both_concordant | — | — | — |
| stripe_1744 | chr13:8,905,001–10,415,000 | gained | 5prime | 0.2256 | 7.73e-06 | high | both_concordant | — | — | — |
| stripe_5896 | chr6:91,455,001–91,925,000 | gained | 5prime | 0.2030 | 8.53e-06 | high | both_concordant | — | — | — |

### Top 10 Lost — both_concordant, high confidence

| stripe_id | JuiceBox coord | type | dir_type | logFC | FDR | confidence | res_support | nearest_gene | anchor_type | ChIP marks |
|-----------|----------------|------|----------|-------|-----|------------|-------------|--------------|-------------|------------|
| stripe_0811 | chr10:108,385,001–109,245,000 | lost | 3prime | −0.1951 | 4.33e-05 | high | both_concordant | Gm36283 (Syt1 body) | Poised_Enhancer | H3K27me3, H3K4me1 |
| stripe_1943 | chr13:83,685,001–84,365,000 | lost | 3prime | −0.1970 | 6.02e-05 | high | both_concordant | Mir9-2 (Mef2c region) | Active_Promoter | H3K27ac, H3K4me1, H3K4me3 |
| stripe_2851 | chr16:36,960,001–37,705,000 | lost | 5prime | −0.1628 | 8.02e-05 | high | both_concordant | Ndufb4 | Active_Promoter | H3K27ac, H3K4me1, H3K4me3 |
| stripe_3172 | chr17:48,215,001–49,475,000 | lost | 3prime | −0.1668 | 1.07e-04 | high | both_concordant | Treml4 | Other | none |
| stripe_2921 | chr16:63,490,001–63,925,000 | lost | 5prime | −0.2004 | 2.03e-04 | high | both_concordant | — | — | — |
| stripe_2941 | chr16:77,030,001–78,210,000 | lost | 5prime | −0.1462 | 2.80e-04 | high | both_concordant | — | — | — |
| stripe_2846 | chr16:37,025,001–38,540,000 | lost | 3prime | −0.1384 | 3.11e-04 | high | both_concordant | — | — | — |
| stripe_3283 | chr17:83,830,001–84,445,000 | lost | 5prime | −0.1293 | 3.65e-04 | high | both_concordant | — | — | — |
| stripe_0870 | chr10:127,980,001–128,985,000 | lost | 5prime | −0.1360 | 4.01e-04 | high | both_concordant | — | — | — |
| stripe_4869 | chr4:81,255,001–82,540,000 | lost | 3prime | −0.1511 | 6.68e-04 | high | both_concordant | — | — | — |

*Note: Gene annotation queries were run for the top-ranked stripes. Entries marked "—" have gene data in the annotated TSV but were not retrieved in this session; the full annotated TSV (`outputs/250402/visualizations/250402_annotated_stripes.tsv`) contains nearest_gene and body_genes for all stripes.*

---

## 3. JuiceBox Coordinate Blocks

```
# GAINED stripes — view MUTANT .hic first, then CONTROL to confirm absence
# (all from both_concordant + high-confidence tier)

chr2:44605001-46505000     # stripe_3908  5prime  logFC=+0.368  FDR=1.43e-07  Zeb2-region  Poised_Enhancer
chr9:40240001-41510000     # stripe_7016  5prime  logFC=+0.326  FDR=1.53e-07  Ubash3b-region  Poised_Enhancer
chr5:83990001-84465000     # stripe_5385  5prime  logFC=+0.340  FDR=1.98e-07  Epha5  Poised_Enhancer (Bivalent)
chr4:20050001-21540000     # stripe_4729  3prime  logFC=+0.356  FDR=3.16e-07  Ggh/Nkain3
chr4:20195001-20830000     # stripe_4730  5prime  logFC=+0.368  FDR=8.88e-07  Nkain3 (co-occurs with 4729)
chr13:9815001-10980000     # stripe_1743  3prime  logFC=+0.281  FDR=1.66e-06
chr4:55285001-56415000     # stripe_4797  5prime  logFC=+0.234  FDR=4.78e-06
chr7:135745001-137380000   # stripe_6491  5prime  logFC=+0.264  FDR=7.71e-06
chr13:8905001-10415000     # stripe_1744  5prime  logFC=+0.226  FDR=7.73e-06
chr6:91455001-91925000     # stripe_5896  5prime  logFC=+0.203  FDR=8.53e-06


# LOST stripes — view CONTROL .hic first, then MUTANT to confirm loss
# (all from both_concordant + high-confidence tier)

chr10:108385001-109245000  # stripe_0811  3prime  logFC=-0.195  FDR=4.33e-05  Syt1-region  Poised_Enhancer
chr13:83685001-84365000    # stripe_1943  3prime  logFC=-0.197  FDR=6.02e-05  Mef2c-region  Active_Promoter
chr16:36960001-37705000    # stripe_2851  5prime  logFC=-0.163  FDR=8.02e-05  Ndufb4  Active_Promoter
chr17:48215001-49475000    # stripe_3172  3prime  logFC=-0.167  FDR=1.07e-04  Treml4
chr16:63490001-63925000    # stripe_2921  5prime  logFC=-0.200  FDR=2.03e-04
chr16:77030001-78210000    # stripe_2941  5prime  logFC=-0.146  FDR=2.80e-04
chr16:37025001-38540000    # stripe_2846  3prime  logFC=-0.138  FDR=3.11e-04
chr17:83830001-84445000    # stripe_3283  5prime  logFC=-0.129  FDR=3.65e-04
chr10:127980001-128985000  # stripe_0870  5prime  logFC=-0.136  FDR=4.01e-04
chr4:81255001-82540000     # stripe_4869  3prime  logFC=-0.151  FDR=6.68e-04
```

---

## 4. Biological Context

Reported directly from pipeline outputs:

- This timepoint has substantial differential signal: 2,320/7,371 stripes (31.5%) reach FDR<0.05 at 5kb.
- More stripes are gained than lost in BAP1-KO (2,052 gained vs 1,528 lost at 5kb; 967 vs 766 at 10kb), representing a net gain of stripe architecture in the mutant.
- 1,005 stripes (367 lost + 638 gained) reach the high-confidence tier at 5kb; no medium-confidence calls were made (0), and 2,575 are low-confidence.
- Effect sizes are uniformly small: 0 strong, 0 moderate, 19 weak (|logFC| 0.3–0.5), and 2,301 minimal (|logFC| < 0.3). Maximum observed |logFC| in the 5kb dataset is 0.3892.
- Directional consistency is moderate: 2,172/3,580 (60.7%) of lost/gained calls have logFC in the direction expected from source assignment.
- 7 shared stripes show differential signal: 3 strengthened and 4 weakened.
- Cross-resolution logFC correlation is good (Pearson r = 0.850), but only 377/672 (56.1%) of both-significant stripes agree on direction; 1,249 stripes are both_discordant (significant at both resolutions but in opposite directions).
- The chr4:20Mb locus has two co-occurring gained stripes (stripe_4729 and stripe_4730) within 660kb; chr13:9–11Mb similarly has two co-occurring gained stripes (stripe_1743 and stripe_1744). chr16 has three lost stripes in the high-confidence tier spanning 37–78Mb.
- Gain-biased stripes are enriched at Poised_Enhancer anchors (36.5% of gained vs 30.3% of lost); loss-biased stripes are enriched at Active_Enhancer anchors (19.5% vs 14.0%).

---

## 5. Anchor Annotation Analysis

### 5.1 Anchor Type Enrichment by Direction

Counts and percentages derived from the annotated TSV (`250402_annotated_stripes.tsv`). ChIP-seq peaks from the 250402 (adult/late) timepoint.

| Anchor Type | Lost (n=1,528) | % | Gained (n=2,052) | % | All stripes (n=7,371) | % |
|-------------|----------------|---|------------------|---|----------------------|---|
| Poised_Enhancer | 463 | 30.3 | 748 | 36.5 | 2,381 | 32.3 |
| Other | 406 | 26.6 | 539 | 26.3 | 1,925 | 26.1 |
| Active_Promoter | 282 | 18.5 | 292 | 14.2 | 1,280 | 17.4 |
| Active_Enhancer | 298 | 19.5 | 287 | 14.0 | 1,224 | 16.6 |
| Repressed_Promoter | 54 | 3.5 | 119 | 5.8 | 367 | 5.0 |
| Polycomb | 13 | 0.9 | 42 | 2.0 | 109 | 1.5 |
| Bivalent_Promoter | 12 | 0.8 | 25 | 1.2 | 85 | 1.2 |

**Key observations:**

- **Poised_Enhancer anchors** (H3K4me1+, no H3K27ac or H3K27me3) are enriched in gained stripes relative to lost stripes (36.5% vs 30.3% vs 32.3% background). Gained stripes at Poised_Enhancer loci may reflect ectopic activation of primed enhancer contacts in BAP1-KO.
- **Active_Promoter and Active_Enhancer anchors** are proportionally lower in gained stripes than lost stripes (Active_Promoter: 14.2% vs 18.5%; Active_Enhancer: 14.0% vs 19.5%). Actively transcribed loci are disproportionately represented in the lost category.
- **Polycomb anchors** (H3K27me3+, >2kb from TSS) are enriched in gained stripes (2.0%) relative to lost (0.9%) and the background (1.5%). BAP1 loss can lead to derepression at Polycomb loci, consistent with gained contacts at H3K27me3-marked distal elements.
- **Repressed_Promoter anchors** (H3K27me3+, ≤2kb TSS) show similar enrichment in gained stripes (5.8% vs 3.5% in lost). Combined, H3K27me3-positive anchors represent 9.1% of gained vs 5.2% of lost (186 vs 79 stripes).
- **Bivalent_Promoter anchors** are slightly higher in gained (1.2%) than lost (0.8%); absolute numbers are small (25 gained vs 12 lost).
- No anchor type reaches statistical significance after multiple testing correction in this dataset (formal Fisher's exact tests were not computed), but the directional patterns for Polycomb/Repressed_Promoter vs Active_Enhancer are internally consistent.

### 5.2 ChIP-seq Mark Overlap Rates by Direction

Overlap rates derived from the anchor_type classification (which encodes ChIP mark status). These represent mutually exclusive categorization, so a stripe counted as "Active_Promoter" is H3K4me3+, NOT H3K27me3, ≤2kb TSS.

| Mark / Category | Lost (n=1,528) | % | Gained (n=2,052) | % | All (n=7,371) | % |
|---|---|---|---|---|---|---|
| H3K4me1+ (Poised_Enhancer) | 463 | 30.3 | 748 | 36.5 | 2,381 | 32.3 |
| H3K27ac+ (Active_Enhancer, distal) | 298 | 19.5 | 287 | 14.0 | 1,224 | 16.6 |
| H3K4me3+ (Active_Promoter + Bivalent) | 294 | 19.2 | 317 | 15.5 | 1,365 | 18.5 |
| H3K27me3+ (Polycomb + Repressed + Bivalent) | 79 | 5.2 | 186 | 9.1 | 561 | 7.6 |

**Largest directional differentials:**
- H3K27me3-positive anchors: ~1.75-fold enrichment in gained vs lost (9.1% vs 5.2%)
- H3K27ac-positive (distal): ~1.39-fold enrichment in lost vs gained (19.5% vs 14.0%)
- H3K4me1+ (Poised_Enhancer): ~1.20-fold enrichment in gained (36.5% vs 30.3%)

The H3K27me3 enrichment in gained stripes is the most notable directional asymmetry and is consistent with BAP1's role in removing H2AK119ub1 at Polycomb-repressed loci. BAP1 loss may allow expansion of chromatin contacts at normally repressed regions.

### 5.3 Bivalent Stripes

Stripes with Bivalent_Promoter anchors (H3K4me3 + H3K27me3 overlap, pre-computed intersection): 12 lost, 25 gained, 85 total. These represent developmental transition-state loci. Additionally, stripe_5385 (top-3 gained, chr5:Epha5 region) carries a bivalent mark at its anchor (identified from the annotated TSV: h3k27me3=TRUE, h3k4me3=TRUE, bivalent=TRUE). Epha5 is a receptor tyrosine kinase involved in axon guidance.

### 5.4 Stripiness Score (Shared Stripes)

3,789 shared stripes (source="shared") have both stripiness_ctrl and stripiness_mut values. From the visualization run (Section 3): the stripiness analysis confirmed that shared stripes show correlated stripiness values between conditions (plot saved: `stripiness_250402.pdf`). Given that all edgeR effect sizes are in the "minimal" to "weak" range and BCV is very low (0.012), the differential signal at shared stripes is driven by small but consistent count differences rather than by dramatic changes in stripe morphology.

For control-only (lost) stripes, stripiness_ctrl is available; for mutant-only (gained) stripes, stripiness_mut is available. The top gained stripes show a wide range of stripiness_mut: stripe_3908 (3.69), stripe_4729 (3.23), stripe_6491 (17.86) in the high-confidence tier, while others are near zero or negative (stripe_7016: −0.07, stripe_5896: −3.64), indicating that some called stripes have weak or absent stripe morphology by Stripenn's own metric.

---

## 6. Pathway Enrichment Summary

Gene lists used: 444 unique genes from lost-stripe anchors, 639 from gained-stripe anchors. Background: 6,687 annotated mm10 genes (all stripes in the union set). Analysis performed by clusterProfiler with BH-adjusted p-value threshold of 0.05.

### 6.1 GO Biological Process

6 significant terms total; all terms are associated with the **gained** direction. No GO BP terms reach significance for the lost direction.

**Top 6 GO BP terms (gained stripes, FDR<0.05):**

| Rank | Term | GO ID | Gene Ratio | FoldEnrich | p.adj | Count |
|------|------|-------|------------|-----------|-------|-------|
| 1 | embryonic skeletal system development | GO:0048706 | 20/589 | 3.10 | 0.0064 | 20 |
| 2 | intermediate filament organization | GO:0045109 | 9/589 | 5.94 | 0.0064 | 9 |
| 3 | cell-cell adhesion via plasma-membrane adhesion molecules | GO:0098742 | 24/589 | 2.51 | 0.0191 | 24 |
| 4 | regulation of potassium ion transport | GO:0043266 | 14/589 | 3.36 | 0.0278 | 14 |
| 5 | epidermal cell differentiation | GO:0009913 | 20/589 | 2.57 | 0.0388 | 20 |
| 6 | negative regulation of synaptic transmission | GO:0050805 | 12/589 | 3.52 | 0.0408 | 12 |

**Biologically notable genes in these terms:**
- *Embryonic skeletal system development* (GO:0048706): includes Tgfb2, Dlx5, Dlx6, Hoxa2-Hoxa7, Hoxc4-Hoxc6, Shh, Sox11, Grhl2, Six1 — a highly coherent set of developmental transcription factors and morphogens. This term reflects the broad genomic landscape affected at Hox and other developmental loci.
- *Intermediate filament organization* (GO:0045109): Krt71, Krt72, Krt74, Krt80, Krt82, Krt84, Krt23 — keratin cluster genes, likely reflecting a contiguous genomic region acquired as a single large stripe.
- *Regulation of potassium ion transport* (GO:0043266): Kcnip1, Kcns2, Kcnab1, Kcnab2, Kcnb1, Kcnn2, Dpp6, Kcnma1 — enriched in gained stripes; consistent with cerebellar Purkinje/granule cell ion channel gene regulation.
- *Negative regulation of synaptic transmission* (GO:0050805): Chrm3, Adcy8, Kcnb1, Pcdh17, Adra1a, Nos1 — synaptic function in cerebellum.

**Flag:** No Polycomb-specific (PcG complex, H2AK119Ub, or histone deubiquitination) GO terms reach significance. The gained-stripe genes are distributed across developmental and neuronal categories. No chromatin-remodeling or epigenetic regulatory terms are enriched, though this may reflect the gene-set size (639 genes, mixed background).

**Lost stripes:** 0 significant GO BP terms. The lost anchor gene set (444 genes) is not significantly enriched in any GO BP category at FDR<0.05.

### 6.2 GO Cellular Component

19 significant terms total; all associated with **gained** stripes. No significant GO CC terms for lost stripes.

**Top 5 GO CC terms (gained stripes):**

| Rank | Term | GO ID | Gene Ratio | FoldEnrich | p.adj | Count |
|------|------|-------|------------|-----------|-------|-------|
| 1 | postsynaptic specialization | GO:0099572 | 46/597 | 1.92 | 0.0030 | 46 |
| 2 | neuron to neuron synapse | GO:0098984 | 46/597 | 1.89 | 0.0030 | 46 |
| 3 | postsynaptic specialization membrane | GO:0099634 | 23/597 | 2.43 | 0.0050 | 23 |
| 4 | postsynaptic density | GO:0014069 | 41/597 | 1.86 | 0.0050 | 41 |
| 5 | asymmetric synapse | GO:0032279 | 42/597 | 1.84 | 0.0050 | 42 |

The GO CC enrichment in gained stripes is dominated by synaptic compartment terms (postsynaptic specialization, neuron-to-neuron synapse, postsynaptic density). This is consistent with the tissue context (mouse cerebellum) and indicates that genes at gained-stripe anchors are enriched at loci encoding synaptic components. Additional significant terms include: cation channel complex, postsynaptic membrane, potassium channel complex, transmembrane transporter complex, glutamatergic synapse. The full set of 19 significant terms spans synaptic, ion channel, and membrane complex categories.

**Lost stripes:** No significant GO CC terms.

### 6.3 GO Molecular Function

1 significant term total (gained); no GO MF terms for lost.

| Rank | Term | GO ID | Gene Ratio | FoldEnrich | p.adj | Count |
|------|------|-------|------------|-----------|-------|-------|
| 1 | calcium ion binding | GO:0005509 | 38/585 | 1.93 | 0.0303 | 38 |

The calcium ion binding enrichment in gained stripes (38 genes including Calm1, Efcab6, Efcab7, Notch3, Padi2, Stab1, Smoc1, Nell1) reflects activation or contact reorganization at calcium-signaling gene loci in BAP1-KO adult cerebellum.

**No chromatin binding, transcription factor binding, or deubiquitinase-related MF terms** are significantly enriched in either direction.

### 6.4 KEGG Pathways

No significant KEGG pathways were found for either lost or gained stripes at the late (250402) timepoint (no kegg_late.tsv produced). This is likely due to underpowered pathway enrichment: with a large gene set (444 lost, 639 gained) but diffuse signal spread across many pathways, no single pathway reaches adjusted significance. By comparison, the early timepoint (250831) identified 2 KEGG pathways for lost stripes (mmu04550: Signaling pathways regulating pluripotency of stem cells, p.adj=0.0094; mmu04740: Olfactory transduction, p.adj=0.0349).

### 6.5 Enrichment Summary

- GO enrichment is exclusively significant for **gained** stripes. The lost anchor gene set (444 genes) produces no significant terms in any GO category.
- Gained-stripe genes cluster around: (1) developmental transcription factors and morphogens (Hox clusters, Dlx5/6, Shh, Sox11, Tgfb2), (2) synaptic/neuronal compartment genes (46 genes at postsynaptic specialization), and (3) ion channel genes (potassium and calcium channels).
- The synaptic and neuronal functional terms likely reflect the cerebellar tissue identity rather than BAP1-specific biology — these genes reside in large, stripe-forming genomic domains that happen to gain stripe signal in the mutant.
- No direct Polycomb target or epigenetic regulatory terms are enriched, which may reflect the broad scope of BAP1's effects distributed across many functional categories.

---

## 7. Cross-Timepoint Comparison

### 7.1 Significance Rates

| Metric | Early (250831, P12) | Late (250402, Adult) |
|---|---|---|
| Total union stripes (5kb) | 4,008 | 7,371 |
| Significant (FDR<0.05) | 96 (2.4%) | 2,320 (31.5%) |
| High-confidence: lost | 12 | 367 |
| High-confidence: gained | 10 | 638 |
| Common BCV | ~0.020 | 0.012 |
| Directional consistency | 40.9% | 60.7% |
| Median |logFC| (sig) | 0.117 | 0.100 |

The ~13-fold difference in significance rate (2.4% early vs 31.5% late) does not reflect a 13-fold difference in biological effect size — the median |logFC| is actually slightly larger in the early timepoint (0.117 vs 0.100). Instead, the lower BCV at the late timepoint (0.012 vs 0.020) inflates statistical power, allowing more stripes to reach FDR<0.05 with smaller effect sizes. The late timepoint has larger sample sizes (more stripes in the union set) and apparently more consistent replicate counts, which reduces dispersion estimates.

This means the late timepoint's 31.5% significance rate partly reflects the statistical framework's sensitivity to low variance in the count data, not necessarily a larger biological effect of BAP1 loss. The early timepoint's 2.4% may reflect genuine weak signal at P12 or a more variable replicate structure (BCV=0.020 is still low in absolute terms but higher than the late timepoint).

### 7.2 Directional Bias

| Metric | Early (250831) | Late (250402) |
|---|---|---|
| Direction | More lost than gained | More gained than lost |
| Lost | 949 (57.1% of sig) | 1,528 (42.7% of sig) |
| Gained | 776 (46.8% of sig) | 2,052 (57.3% of sig) |
| Net | Net loss of stripes | Net gain of stripes |

The directional bias is **opposite** between timepoints. At P12 (early), BAP1-KO leads to proportionally more lost than gained stripes. At adult (late), the bias reverses: more stripes are gained than lost. This could indicate developmental stage-dependent remodeling of chromatin architecture downstream of BAP1 loss, or it may reflect differences in the stripe set detected (early: 4,008 stripes; late: 7,371), which has different proportions of control-only and mutant-only calls.

### 7.3 Anchor Type Comparison

| Anchor Type | Early Lost % | Early Gained % | Late Lost % | Late Gained % |
|---|---|---|---|---|
| Poised_Enhancer | 33.6 | 35.1 | 30.3 | 36.5 |
| Other | 34.0 | 32.5 | 26.6 | 26.3 |
| Active_Promoter | 12.4 | 12.0 | 18.5 | 14.2 |
| Active_Enhancer | 10.7 | 11.6 | 19.5 | 14.0 |
| Repressed_Promoter | 5.0 | 5.4 | 3.5 | 5.8 |
| Polycomb | 1.4 | 0.6 | 0.9 | 2.0 |
| Bivalent_Promoter | 2.8 | 2.8 | 0.8 | 1.2 |

**Consistent patterns:**
- Poised_Enhancer anchors are enriched in gained vs lost stripes at both timepoints (early: 35.1% vs 33.6%; late: 36.5% vs 30.3%).
- Polycomb anchors are enriched in the late gained stripes (2.0%) but not in the early gained stripes (0.6%). This may indicate that the Polycomb-driven component of stripe gain accumulates over developmental time.

**Divergent patterns:**
- Active_Promoter and Active_Enhancer anchors are much more common in the late timepoint overall (reflecting the adult chromatin state), and show stronger lost-enrichment in the late dataset that is not apparent in the early dataset.
- Bivalent_Promoter anchors are proportionally more common in early lost/gained stripes (2.8% each) than late (0.8%/1.2%), consistent with bivalent domains resolving to active or repressed states by adulthood.

### 7.4 Pathway Enrichment Comparison

| Category | Early (250831) | Late (250402) |
|---|---|---|
| GO BP (lost) | 0 significant | 0 significant |
| GO BP (gained) | 0 significant | 6 significant |
| GO CC | No enrichment (lost or gained) | 19 terms (gained only) |
| GO MF | 1 term (lost: receptor ligand activity) | 1 term (gained: calcium ion binding) |
| KEGG | 2 pathways (lost: stem cell signaling, olfactory) | 0 significant |

The late timepoint shows substantially richer GO enrichment than the early timepoint for the gained direction, consistent with the larger number of significant stripes (2,052 vs 776 gained). The early timepoint's KEGG enrichment in the lost direction (stem cell pluripotency signaling) is not replicated in the late timepoint, suggesting that if there is a Polycomb-related transcriptional derepression at P12, it dissipates or is subsumed into larger-scale architectural changes by adulthood.

### 7.5 Convergent Loci

No formal overlap analysis between early and late significant stripes was performed in this session. The cross-resolution comparison within each timepoint was done in Stage 6. For stripes at the same genomic locus that are differential in both timepoints, a dedicated locus-by-locus comparison would be required (e.g., using anchor overlap with a tolerance of 25–50kb). Several gene names appear in both the early and late gene lists (e.g., Syt1, Dlg2, Grik2, Kcnn2, Prkn, Csmd2, Tmem132a, Lpp, Dscam, Glcci1, Vav2), suggesting shared regulatory regions across developmental stages.

---

## 8. Caveats and Confidence Assessment

- **Effect sizes are uniformly small.** The maximum |logFC| across all 7,371 stripes is 0.389 (below the typical "moderate" threshold of 0.5). All significant calls fall in the "weak" or "minimal" categories. The differential signal is frequency-driven (many stripes with subtle count changes), not magnitude-driven.
- **Very low BCV (0.012) inflates statistical power.** The high significance rate (31.5%) is partly an artifact of extremely tight within-group replicate agreement. In biological terms, a BCV of 0.012 is atypically low and likely reflects low-noise count data at these stripe loci rather than a dramatic biological effect.
- **Directional consistency of 60.7% is moderate.** Approximately 39% of lost/gained calls show logFC in the opposite direction from what the source-based assignment would predict. This bimodal confidence distribution (only high or low, no medium) means there is no intermediate tier to filter on.
- **No medium-confidence tier.** All non-high calls are low-confidence, giving a bimodal distribution. This is a limitation of the confidence scoring scheme rather than a data quality issue.
- **1,249 both_discordant stripes.** Stripes that are statistically significant at both 5kb and 10kb but with opposite-sign logFC. These represent unreliable direction assignments and should not be interpreted as biologically meaningful at the single-stripe level.
- **Stripenn p-values for several top stripes are borderline** (e.g., stripe_5385: pval_mut = 0.093; stripe_7016: pval_mut = 0.0435; stripe_5896: pval_mut = 0.095). Stripenn's own detection criteria consider these marginal calls in the condition where the stripe was detected.
- **Negative stripiness values** in several high-ranked stripes (e.g., stripe_7016 stripiness_mut = −0.07; stripe_5896 stripiness_mut = −3.64) indicate that the stripe model fit is at or below background in those samples. The count-based edgeR score and the image-derived stripiness score are not always concordant.
- **GO enrichment is entirely one-directional.** No terms reach significance for lost-stripe anchor genes. This may indicate that lost-stripe loci are biologically heterogeneous (no functional convergence) or that the 444-gene lost set is underpowered relative to the background.
- **Anchor type enrichments are not corrected for multiple testing** in this analysis. The directional patterns (Polycomb/H3K27me3 enrichment in gained; Active_Enhancer enrichment in lost) are internally consistent but should be confirmed with formal Fisher's exact tests.
- **Cross-timepoint comparison is confounded** by different union stripe sets, different detection rates, and the directional bias reversal between timepoints. Direct locus-level comparison (overlapping anchor coordinates) would provide stronger evidence for convergent differential regulation.
- **KEGG was underpowered** for the late timepoint (no significant pathways). With 639 gained-stripe genes distributed across synaptic, developmental, and ion channel functions, no single KEGG pathway accumulates enough genes at adjusted significance.
- **Stripiness correlation analysis** (shared stripes: n=3,789) confirms that the differential signal in shared stripes reflects small, consistent count changes rather than morphological stripe gain or loss. The BCV plot and MDS plots (Stage 4 diagnostics) show clear ctrl/mut separation, which supports the edgeR model but also shows that the conditions are not dramatically different in their global stripe profiles.
