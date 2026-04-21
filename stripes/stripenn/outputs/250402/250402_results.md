# Late Timepoint (250402/Adult) Differential Stripe Analysis — Stripenn

**Dataset**: BAP1-KO vs Wildtype Mouse Cerebellum (mm10)
**Timepoint**: Late (250402 / Adult)
**Pipeline**: Stripenn v1.1.65 (Canny edge detection + pixel saturation)
**Comparison**: Mutant (BAP1-KO) vs Control (Wildtype), n=3 per condition
**Resolutions**: 5kb (primary) and 10kb (validation)

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

### Stripe Geometry (5kb, all union stripes)

| Metric | Lost | Gained | Unchanged |
|--------|------|--------|-----------|
| Median length (kb) | 645 | 625 | 665 |
| Median width (kb) | 20 | 20 | 20 |
| N (3prime) | 779 | 1,008 | 1,918 |
| N (5prime) | 749 | 1,044 | 1,866 |

### ChIP-seq Anchor Annotation
*(Pending visualization script — requires HPC run for full ChIP-seq peak overlap analysis)*

---

## 2. High-Priority Stripes for Validation

Top 10 stripes by lowest FDR, filtered to: `resolution_support == "both_concordant"`, `direction_confidence == "high"`, `direction IN (lost, gained)`. All 10 are gained stripes; the highest-confidence lost stripes follow in Section 3.

### Top 10 Gained — both_concordant, high confidence

| stripe_id | JuiceBox coord | type | dir_type | logFC | FDR | confidence | res_support | strip_ctrl | strip_mut | pval_ctrl | pval_mut | len_kb | wid_kb | Rationale |
|-----------|----------------|------|----------|-------|-----|------------|-------------|------------|-----------|-----------|----------|--------|--------|-----------|
| stripe_3908 | chr2:44,605,001–46,505,000 | gained | 5prime | 0.3684 | 1.43e-07 | high | both_concordant | NA | 3.69 | NA | 0.028 | 1,800 | 20 | Highest FDR significance; long 1.8Mb stripe with strong stripiness in mutant |
| stripe_7016 | chr9:40,240,001–41,510,000 | gained | 5prime | 0.3256 | 1.53e-07 | high | both_concordant | NA | −0.07 | NA | 0.044 | 1,170 | 30 | Second-lowest FDR; wide anchor (30kb); stripiness near zero suggests weak but emergent signal |
| stripe_5385 | chr5:83,990,001–84,465,000 | gained | 5prime | 0.3402 | 1.98e-07 | high | both_concordant | NA | 0.24 | NA | 0.093 | 375 | 25 | Compact stripe (375kb); borderline Stripenn p-value but FDR highly significant |
| stripe_4729 | chr4:20,050,001–21,540,000 | gained | 3prime | 0.3562 | 3.16e-07 | high | both_concordant | NA | 3.23 | NA | 0.022 | 1,390 | 25 | Long 3prime stripe with high stripiness; neighboring stripe_4730 on same chromosome |
| stripe_4730 | chr4:20,195,001–20,830,000 | gained | 5prime | 0.3675 | 8.88e-07 | high | both_concordant | NA | 1.08 | NA | 0.058 | 535 | 20 | Adjacent to stripe_4729 at chr4:20Mb locus — co-occurring gained stripes in BAP1-KO |
| stripe_1743 | chr13:9,815,001–10,980,000 | gained | 3prime | 0.2806 | 1.66e-06 | high | both_concordant | NA | −0.50 | NA | 0.083 | 1,065 | 35 | Wide anchor (35kb); neighboring stripe_1744 on same chromosome |
| stripe_4797 | chr4:55,285,001–56,415,000 | gained | 5prime | 0.2342 | 4.78e-06 | high | both_concordant | NA | −0.70 | NA | 0.044 | 1,030 | 20 | Long stripe; three gained stripes on chr4 in top 10 |
| stripe_6491 | chr7:135,745,001–137,380,000 | gained | 5prime | 0.2637 | 7.71e-06 | high | both_concordant | NA | 17.86 | NA | 0.012 | 1,535 | 20 | Highest stripiness_mut in top 10 (17.86); long 1.5Mb stripe |
| stripe_1744 | chr13:8,905,001–10,415,000 | gained | 5prime | 0.2256 | 7.73e-06 | high | both_concordant | NA | 5.13 | NA | 0.034 | 1,410 | 20 | Co-occurs with stripe_1743 at chr13:9–11Mb; paired stripes flanking same locus |
| stripe_5896 | chr6:91,455,001–91,925,000 | gained | 5prime | 0.2030 | 8.53e-06 | high | both_concordant | NA | −3.64 | NA | 0.095 | 370 | 30 | Compact stripe; borderline Stripenn p-value; negative stripiness indicates edge case |

### Top 10 Lost — both_concordant, high confidence

| stripe_id | JuiceBox coord | type | dir_type | logFC | FDR | confidence | res_support | strip_ctrl | strip_mut | pval_ctrl | pval_mut | len_kb | wid_kb | Rationale |
|-----------|----------------|------|----------|-------|-----|------------|-------------|------------|-----------|-----------|----------|--------|--------|-----------|
| stripe_0811 | chr10:108,385,001–109,245,000 | lost | 3prime | −0.1951 | 4.33e-05 | high | both_concordant | 1.56 | NA | 0.054 | NA | 760 | 20 | Highest-confidence lost stripe; clear control-only with positive stripiness |
| stripe_1943 | chr13:83,685,001–84,365,000 | lost | 3prime | −0.1970 | 6.02e-05 | high | both_concordant | −8.96 | NA | 0.048 | NA | 580 | 20 | Strong edgeR significance; negative stripiness is an edge detection artifact |
| stripe_2851 | chr16:36,960,001–37,705,000 | lost | 5prime | −0.1628 | 8.02e-05 | high | both_concordant | −31.62 | NA | 0.036 | NA | 645 | 15 | Narrow anchor (15kb); two lost stripes at chr16:37Mb locus (with stripe_2846) |
| stripe_3172 | chr17:48,215,001–49,475,000 | lost | 3prime | −0.1668 | 1.07e-04 | high | both_concordant | 2.71 | NA | 0.080 | NA | 1,160 | 25 | Long 1.16Mb stripe; largest lost stripe in top tier |
| stripe_2921 | chr16:63,490,001–63,925,000 | lost | 5prime | −0.2004 | 2.03e-04 | high | both_concordant | 1.35 | NA | 0.046 | NA | 335 | 20 | Compact; highest |logFC| among lost (0.200); two lost stripes on chr16:63–78Mb |
| stripe_2941 | chr16:77,030,001–78,210,000 | lost | 5prime | −0.1462 | 2.80e-04 | high | both_concordant | 1.63 | NA | 0.072 | NA | 1,080 | 20 | Second chr16 lost stripe in same region; 77–78Mb |
| stripe_2846 | chr16:37,025,001–38,540,000 | lost | 3prime | −0.1384 | 3.11e-04 | high | both_concordant | −7.33 | NA | 0.068 | NA | 1,415 | 20 | Longest lost stripe in top tier (1.4Mb); flanks stripe_2851 at chr16:37Mb |
| stripe_3283 | chr17:83,830,001–84,445,000 | lost | 5prime | −0.1293 | 3.65e-04 | high | both_concordant | −1.59 | NA | 0.088 | NA | 515 | 30 | Wide anchor (30kb); second chr17 lost stripe in top tier |
| stripe_0870 | chr10:127,980,001–128,985,000 | lost | 5prime | −0.1360 | 4.01e-04 | high | both_concordant | −0.37 | NA | 0.058 | NA | 905 | 25 | Second chr10 lost stripe; 128Mb region distinct from stripe_0811 |
| stripe_4869 | chr4:81,255,001–82,540,000 | lost | 3prime | −0.1511 | 6.68e-04 | high | both_concordant | 1.14 | NA | 0.087 | NA | 1,185 | 20 | Long stripe; chr4 has both gained (top 3) and lost stripes in different regions |

---

## 3. JuiceBox Coordinate Blocks

```
# GAINED stripes — view MUTANT .hic first, then CONTROL to confirm absence
# (all from both_concordant + high-confidence tier)

chr2:44605001-46505000     # stripe_3908  5prime  logFC=+0.368  FDR=1.43e-07
chr9:40240001-41510000     # stripe_7016  5prime  logFC=+0.326  FDR=1.53e-07
chr5:83990001-84465000     # stripe_5385  5prime  logFC=+0.340  FDR=1.98e-07
chr4:20050001-21540000     # stripe_4729  3prime  logFC=+0.356  FDR=3.16e-07
chr4:20195001-20830000     # stripe_4730  5prime  logFC=+0.368  FDR=8.88e-07
chr13:9815001-10980000     # stripe_1743  3prime  logFC=+0.281  FDR=1.66e-06
chr4:55285001-56415000     # stripe_4797  5prime  logFC=+0.234  FDR=4.78e-06
chr7:135745001-137380000   # stripe_6491  5prime  logFC=+0.264  FDR=7.71e-06
chr13:8905001-10415000     # stripe_1744  5prime  logFC=+0.226  FDR=7.73e-06
chr6:91455001-91925000     # stripe_5896  5prime  logFC=+0.203  FDR=8.53e-06


# LOST stripes — view CONTROL .hic first, then MUTANT to confirm loss
# (all from both_concordant + high-confidence tier)

chr10:108385001-109245000  # stripe_0811  3prime  logFC=-0.195  FDR=4.33e-05
chr13:83685001-84365000    # stripe_1943  3prime  logFC=-0.197  FDR=6.02e-05
chr16:36960001-37705000    # stripe_2851  5prime  logFC=-0.163  FDR=8.02e-05
chr17:48215001-49475000    # stripe_3172  3prime  logFC=-0.167  FDR=1.07e-04
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

---

## 5. Caveats

- Effect sizes are uniformly small — the maximum |logFC| across all 7,371 stripes is 0.389 (below the typical "moderate" threshold of 0.5). All significant calls fall in the "weak" or "minimal" categories.
- Directional consistency of 60.7% is moderate — approximately 39% of lost/gained calls show logFC in the opposite direction from what the source-based assignment (control_only → lost, mutant_only → gained) would predict.
- The common BCV of 0.012 is very low, which inflates statistical power and likely drives the high significance rate (31.5%). Low BCV means very small count differences between replicates are declared significant.
- Cross-resolution: 1,249 stripes are "both_discordant" — significant at both 5kb and 10kb but with opposite-sign logFC. This is a concern for interpreting the lost/gained calls and suggests some direction assignments are unreliable.
- Despite the high significance rate, the small fold changes mean these are subtle quantitative differences in stripe intensity, not dramatic structural gains or losses.
- No medium-confidence tier stripes exist (0 medium); all non-high calls are low-confidence, giving a bimodal confidence distribution.
- Stripenn p-values for several top stripes are borderline (e.g., stripe_5385: pval_mut = 0.093; stripe_5896: pval_mut = 0.095), meaning the stripe was marginal by Stripenn's own detection criteria even in the condition where it was called.
- ChIP-seq anchor annotation is pending and not included in this document.
