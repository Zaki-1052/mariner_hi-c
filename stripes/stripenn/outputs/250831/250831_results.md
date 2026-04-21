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
| Median width (kb) | 25 | 25 | 25 |
| N (3prime) | 477 | 377 | 1043 |
| N (5prime) | 472 | 399 | 1240 |

### ChIP-seq Anchor Annotation

*(Pending visualization script — requires HPC run for full ChIP-seq peak overlap analysis)*

---

## 2. High-Priority Stripes for Validation

Top 10 stripes filtered to: direction IN (lost, gained) AND direction\_confidence == "high", sorted by FDR ascending.

| stripe\_id | chr | Anchor (pos1–pos2) | Extent (pos3–pos4) | Direction | dir\_type | logFC | FDR | Confidence | Resolution support | in\_10kb | pval\_ctrl | pval\_mut | Stripiness ctrl | Stripiness mut | Length (kb) | Width (kb) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| stripe\_2752 | chr4 | 125,490,001–125,510,000 | 125,490,001–125,970,000 | gained | 3prime | +0.1396 | 0.0155 | high | 5kb\_only | No | NA | 0.093 | NA | 1.112 | 480 | 20 |
| stripe\_3157 | chr6 | 51,845,001–51,865,000 | 51,465,001–51,865,000 | lost | 5prime | −0.1250 | 0.0173 | high | 5kb\_only | No | 0.067 | NA | 3.478 | NA | 400 | 20 |
| stripe\_2579 | chr4 | 3,595,001–3,620,000 | 3,595,001–3,935,000 | lost | 3prime | −0.1056 | 0.0290 | high | 5kb\_only | No | 0.078 | NA | 1.963 | NA | 340 | 25 |
| stripe\_3636 | chr8 | 64,010,001–64,040,000 | 64,010,001–64,740,000 | lost | 3prime | −0.1041 | 0.0307 | high | both\_concordant | Yes | 0.038 | NA | 3.784 | NA | 730 | 30 |
| stripe\_1504 | chr6 | 43,360,001–43,410,000 | 43,360,001–45,650,000 | lost | 3prime | −0.0876 | 0.0348 | high | 10kb\_only | Yes | 0.088 | NA | −0.098 | NA | 2290 | 50 |
| stripe\_1554 | chr6 | 111,040,001–111,110,000 | 111,040,001–112,190,000 | lost | 3prime | −0.0833 | 0.0366 | high | 10kb\_only | Yes | 0.071 | NA | 2.034 | NA | 1150 | 70 |
| stripe\_3805 | chr9 | 41,920,001–41,955,000 | 41,390,001–41,955,000 | gained | 5prime | +0.1046 | 0.0381 | high | 5kb\_only | No | NA | 0.046 | NA | 2.674 | 565 | 35 |
| stripe\_2321 | chr2 | 151,250,001–151,270,000 | 151,250,001–151,465,000 | gained | 3prime | +0.2112 | 0.0471 | high | 5kb\_only | No | NA | 0.074 | NA | −0.067 | 215 | 20 |
| stripe\_2170 | chr2 | 70,545,001–70,570,000 | 70,120,001–70,570,000 | gained | 5prime | +0.0904 | 0.0481 | high | 5kb\_only | No | NA | 0.050 | NA | 3.183 | 450 | 25 |
| stripe\_3109 | chr6 | 21,215,001–21,235,000 | 21,215,001–21,705,000 | lost | 3prime | −0.1151 | 0.0494 | high | 5kb\_only | No | 0.077 | NA | 0.417 | NA | 490 | 20 |

**Note on stripe\_1504 (Stripiness ctrl = −0.098):** Negative Stripiness scores indicate Stripenn found no net stripe signal relative to background at this locus; the edgeR significance is driven by lower quantification scores in the mutant, not a classically stripe-shaped region in either condition. Treat with extra caution.

**Note on stripe\_2321 (Stripiness mut = −0.067):** Same caveat — absent or weak stripe signal in the mutant despite the "gained" call. The mutant has marginally higher contact counts by edgeR, but the Stripiness score does not confirm a stripe structure.

---

## 3. JuiceBox Coordinate Blocks

Each coordinate block uses 50 kb padding: `chr:min(pos1,pos3)−50000 – max(pos2,pos4)+50000`.

```
# LOST stripes — open CONTROL .hic first; confirm stripe is present, then switch to MUTANT to confirm loss
# Sorted by FDR ascending

stripe_3157  chr6:51415001-51915000    lost  5prime  logFC=-0.1250  FDR=0.0173  5kb_only
stripe_2579  chr4:3545001-3985000      lost  3prime  logFC=-0.1056  FDR=0.0290  5kb_only
stripe_3636  chr8:63960001-64790000    lost  3prime  logFC=-0.1041  FDR=0.0307  both_concordant  [PRIORITY]
stripe_1504  chr6:43310001-45700000    lost  3prime  logFC=-0.0876  FDR=0.0348  10kb_only
stripe_1554  chr6:110990001-112240000  lost  3prime  logFC=-0.0833  FDR=0.0366  10kb_only
stripe_3109  chr6:21165001-21755000    lost  3prime  logFC=-0.1151  FDR=0.0494  5kb_only

# GAINED stripes — open MUTANT .hic first; confirm stripe is present, then switch to CONTROL to confirm absence
# Sorted by FDR ascending

stripe_2752  chr4:125440001-126020000   gained  3prime  logFC=+0.1396  FDR=0.0155  5kb_only
stripe_3805  chr9:41340001-42005000     gained  5prime  logFC=+0.1046  FDR=0.0381  5kb_only
stripe_2321  chr2:151200001-151515000   gained  3prime  logFC=+0.2112  FDR=0.0471  5kb_only  [highest |logFC|]
stripe_2170  chr2:70070001-70620000     gained  5prime  logFC=+0.0904  FDR=0.0481  5kb_only
```

**Priority for manual inspection:** stripe\_3636 (chr8:63960001-64790000) is the only high-confidence stripe supported at both resolutions with concordant direction. This is the most reliable call in this timepoint.

---

## 4. Biological Context

This is a factual summary of what the data shows at this timepoint.

- **Weak overall differential signal:** Only 96 stripes (2.4%) are significant at FDR<0.05 at 5kb; 53 (2.8%) at 10kb. These are the lowest significance rates compared to the late timepoint.
- **Direction bias toward loss:** 949 lost vs 776 gained at 5kb (ratio 1.22:1); 483 vs 401 at 10kb (ratio 1.20:1). More stripes are called control-only than mutant-only at this timepoint.
- **All effect sizes are minimal:** No significant stripe reaches the weak threshold (|logFC| > 0.3). The median logFC among significant stripes is 0.104 at 5kb. The maximum logFC among the top 10 high-confidence stripes is +0.211 (stripe\_2321).
- **Very low BCV:** Common BCV = 0.020 at 5kb and 0.021 at 10kb. This is substantially below typical bulk RNA-seq values (0.2–0.5), indicating very low between-replicate variability relative to mean count levels. This may reflect high replicate uniformity but also limits detection of small fold changes.
- **Low directional consistency:** 40.9% of lost/gained stripes at 5kb (39.9% at 10kb) show direction-consistent edgeR calls. This means that for the majority of detection-based lost/gained calls, the edgeR logFC sign does not reinforce the direction assigned by source (control-only / mutant-only).
- **Sparse cross-resolution support:** Only 22 stripes are significant at both resolutions simultaneously; 15 of those are direction-concordant. Only 1 stripe in the top-10 high-confidence list (stripe\_3636) is both\_concordant across resolutions.
- **No medium-confidence calls:** The medium confidence tier is empty at both resolutions — all lost/gained stripes are either high (22 at 5kb, 24 at 10kb) or low confidence.
- **Six independent lost stripes on chr6:** Of the 10 top high-confidence stripes, 4 are on chr6 (stripe\_3157, stripe\_1504, stripe\_1554, stripe\_3109). This concentration may reflect a regional chromatin feature or could be noise; without enrichment analysis it cannot be interpreted further.

---

## 5. Caveats

- **All effect sizes are minimal:** No significant stripe reaches |logFC| > 0.3. The strongest call (stripe\_2321, logFC = +0.211) is below the weak threshold. Fold changes of this magnitude are at the limit of biological interpretability for chromatin contact data.
- **Low directional consistency (40.9%):** More than half of detection-based lost/gained calls have an edgeR logFC pointing in the opposite direction. This is consistent with the calls being dominated by detection noise rather than true differential signal.
- **Only 22 high-confidence significant stripes:** This is insufficient for robust enrichment analysis (GO/KEGG, ChIP-seq overlap) — any category-level statistics would have very wide confidence intervals.
- **Substantially weaker signal than late timepoint:** The late timepoint (250402/P56) shows 31.5% significant stripes at 5kb vs 2.4% here. The early timepoint may reflect a developmental stage where BAP1-KO effects on chromatin stripes are not yet established, or the differential effect is too small to detect reliably.
- **Very low BCV (0.020):** While low BCV improves nominal power, it also makes the model sensitive to small systematic biases (e.g., library-size differences, technical variation). The TMM normalization factors are very close to 1 (range 0.997–1.003), suggesting the samples are already well-balanced, but the near-zero common dispersion warrants caution.
- **Negative Stripiness scores:** stripe\_1504 (Stripiness ctrl = −0.098) and stripe\_2321 (Stripiness mut = −0.067) have negative Stripiness values in the condition they are called in, indicating Stripenn does not detect a classical stripe pattern at these loci. These calls should not be used as visual validation targets without first confirming the stripe structure visually in JuiceBox.
- **No medium-confidence tier:** The bimodal confidence distribution (high vs low, nothing in between) reflects the strict thresholds used for tier assignment. The 22 high-confidence calls are robustly classified; the 1703 low-confidence calls should be treated as exploratory only.
