# Early Timepoint Differential Stripe Analysis

**Dataset**: BAP1-KO vs Wildtype Mouse Cerebellum (mm10)
**Timepoint**: Early (250831)
**Input file**: `early_stripes.bedpe`
**Analysis date**: 2024-12

---

## 1. Summary Statistics

### Overall Counts

| Metric | Count |
|--------|-------|
| **Total differential stripes** | 19 |
| **Lost stripes** (control-only) | 9 (47%) |
| **Gained stripes** (mutant-only) | 10 (53%) |

### Distribution by Direction Confidence

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| High | 0 | 0 | 0 |
| Medium | 9 | 10 | 19 |
| Low | 0 | 0 | 0 |

**Note**: All stripes have medium direction confidence. This is expected given:
- No stripes achieved FDR < 0.10 (required for "high" confidence)
- All stripes show logFC in the expected direction relative to source

### Distribution by Detection Confidence

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| High | 0 | 0 | 0 |
| Medium | 9 | 8 | 17 |
| Low | 0 | 2 | 2 |

The 2 low-confidence gained stripes (stripe_0247, stripe_0257) notably have the **highest logFC values** (1.0005 and 0.8767), suggesting these may represent real but weakly-detected signals.

### Distribution by Anchor Type

| Anchor Type | Lost | Gained | Total |
|-------------|------|--------|-------|
| **Promoter** | 5 | 6 | 11 (58%) |
| **Active_Enhancer** | 1 | 0 | 1 (5%) |
| **Other** | 3 | 4 | 7 (37%) |

**Observation**: The majority of differential stripes anchor at promoters, consistent with EP-stripe biology. The single Active_Enhancer-anchored stripe (stripe_0128) is lost in mutant.

### ChIP-seq Annotation Breakdown

| ChIP Mark | Lost | Gained | Total |
|-----------|------|--------|-------|
| **H3K27ac** | 6 (67%) | 7 (70%) | 13 (68%) |
| **H3K27me3** | 3 (33%) | 3 (30%) | 6 (32%) |
| **H3K4me1** | 0 | 0 | 0 |

**Note**: H3K4me1 data is not available for early timepoint (only H3K27ac and H3K27me3).

**Polycomb-marked stripes (H3K27me3+)**:
- Lost: Unc13d, Lingo2, Lrrc27 (3 stripes)
- Gained: Clca4b, Garin5a, 6820431F20Rik (3 stripes)

---

## 2. High-Priority Stripes for Validation

Prioritized by: (1) Quagga p-value, (2) in_10kb validation, (3) anchor annotation, (4) gene relevance.

### Priority 1: stripe_0056 (LOST)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr11:116,025,000-116,520,000` |
| **Gene** | *Unc13d* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K27me3+ |
| **Quagga p-val** | 6.55e-32 (exceptional) |
| **in_10kb** | TRUE |
| **logFC** | -0.27 |

**Rationale**: By far the strongest Quagga detection (p=10^-32). *Unc13d* (Munc13-4) is a synaptic vesicle priming protein with roles in exocytosis. The dual H3K27ac/H3K27me3 marks suggest a bivalent or Polycomb-poised state. Loss in BAP1-KO may indicate disrupted Polycomb boundary regulation.

---

### Priority 2: stripe_0265 (LOST)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr8:64,485,000-65,015,000` |
| **Gene** | *Msmo1* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+ |
| **Quagga p-val** | 1.23e-10 |
| **in_10kb** | TRUE |
| **logFC** | -0.30 |

**Rationale**: Strong Quagga detection at promoter of *Msmo1* (methylsterol monooxygenase 1), a cholesterol biosynthesis enzyme. Loss of stripe at metabolic gene promoter may reflect BAP1-dependent chromatin organization changes.

---

### Priority 3: stripe_0091 (LOST)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr14:122,085,000-122,545,000` |
| **Gene** | *Clybl* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+ |
| **Quagga p-val** | 7.31e-09 |
| **in_10kb** | TRUE |
| **logFC** | -0.25 |
| **FDR** | 0.437 (lowest among lost stripes) |

**Rationale**: Second-lowest FDR among lost stripes. *Clybl* (citrate lyase beta-like) is a metabolic enzyme. Good Quagga p-value with 10kb validation.

---

### Priority 4: stripe_0255 (LOST)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr7:138,975,000-139,550,000` |
| **Gene** | *Lrrc27* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K27me3+ |
| **Quagga p-val** | 1.30e-08 |
| **in_10kb** | TRUE |
| **logFC** | -0.23 |

**Rationale**: Bivalent H3K27ac/H3K27me3 marks suggest Polycomb-poised promoter. Strong detection with 10kb validation. *Lrrc27* (leucine-rich repeat containing 27) function less characterized but pattern fits BAP1-Polycomb biology.

---

### Priority 5: stripe_0247 (GAINED) - Highest logFC

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr7:104,450,000-105,350,000` |
| **Gene** | *Trim30e-ps1* (0 bp to TSS) |
| **Anchor type** | Other |
| **ChIP marks** | None |
| **Quagga p-val** | 9.68e-12 |
| **in_10kb** | FALSE |
| **logFC** | **1.00** (2-fold increase) |
| **FDR** | 0.437 |

**Rationale**: **Highest logFC in dataset** representing 2-fold signal increase in mutant. Despite low detection_confidence (no 10kb validation), the exceptional logFC and strong Quagga p-value warrant visual inspection. Pseudogene anchor but may represent ectopic activation.

---

### Priority 6: stripe_0257 (GAINED) - Second-highest logFC

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr8:20,240,000-21,170,000` |
| **Gene** | *6820431F20Rik* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K27me3+ |
| **Quagga p-val** | 2.28e-13 |
| **in_10kb** | FALSE |
| **logFC** | **0.88** (~1.8-fold increase) |

**Rationale**: Second-highest logFC with very strong Quagga detection. Bivalent anchor (H3K27ac+H3K27me3) at an annotated promoter suggests Polycomb-regulated locus gaining activity. Pattern consistent with BAP1-KO derepression.

---

### Priority 7: stripe_0230 (GAINED)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr7:44,420,000-44,915,000` |
| **Gene** | *Garin5a* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K27me3+ |
| **Quagga p-val** | 4.58e-05 |
| **in_10kb** | TRUE |
| **logFC** | 0.20 |

**Rationale**: Bivalent marks at promoter, validated at 10kb resolution. *Garin5a* (GTPase activating RPAP interacting 5A) has roles in development.

---

### Priority 8: stripe_0074 (GAINED)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr13:65,940,000-66,410,000` |
| **Gene** | *Zfp997* (0 bp to TSS) |
| **Anchor type** | Other |
| **ChIP marks** | None |
| **Quagga p-val** | 2.69e-10 |
| **in_10kb** | TRUE |
| **logFC** | 0.35 |

**Rationale**: Strong Quagga p-value with 10kb validation. *Zfp997* is a zinc finger protein (potential transcription factor). New stripe at TF promoter may indicate regulatory activation.

---

### Priority 9: stripe_0136 (GAINED) - Mbd2 proximity

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr18:70,965,000-72,930,000` |
| **Gene** | *Mbd2* (446 kb to TSS) |
| **Anchor type** | Other |
| **ChIP marks** | None |
| **Quagga p-val** | 0.96 (weak) |
| **in_10kb** | FALSE |
| **logFC** | **0.43** |
| **FDR** | 0.437 |

**Rationale**: While Quagga detection is weak, this stripe has the third-highest logFC (0.43) and is near *Mbd2* (methyl-CpG binding domain protein 2). Mbd2 is an epigenetic reader relevant to Polycomb biology. Worth checking despite low detection confidence.

---

### Priority 10: stripe_0009 (GAINED)

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr1:119,485,000-120,000,000` |
| **Gene** | *3830432H09Rik* (0 bp to TSS) |
| **Anchor type** | Other |
| **ChIP marks** | None |
| **Quagga p-val** | 1.14e-09 |
| **in_10kb** | TRUE |
| **logFC** | 0.25 |

**Rationale**: Strong Quagga detection validated at 10kb. Worth inspecting as additional gained stripe example.

---

## 3. JuiceBox Coordinate Blocks

### LOST Stripes (check for stripe disappearance in mutant)

```
# Top LOST stripes - load CONTROL .hic files to see stripe, MUTANT to confirm loss

# Priority 1: Unc13d - exceptional Quagga p-val, bivalent marks
chr11:116,025,000-116,520,000

# Priority 2: Msmo1 - metabolic gene, strong detection
chr8:64,485,000-65,015,000

# Priority 3: Clybl - lowest FDR among lost
chr14:122,085,000-122,545,000

# Priority 4: Lrrc27 - bivalent Polycomb marks
chr7:138,975,000-139,550,000

# Additional lost stripes
chr10:54,445,000-55,095,000    # Man1a region
chr17:89,540,000-90,545,000    # Gm4719 Active_Enhancer
chr6:59,555,000-60,170,000     # Gm5570
chr4:35,480,000-38,795,000     # Lingo2 (very large anchor)
chr8:21,435,000-21,850,000     # Defa20
```

### GAINED Stripes (check for new stripe appearance in mutant)

```
# Top GAINED stripes - load MUTANT .hic files to see stripe, CONTROL to confirm absence

# Priority 5: Trim30e-ps1 - HIGHEST logFC (2-fold)
chr7:104,450,000-105,350,000

# Priority 6: 6820431F20Rik - second-highest logFC, bivalent
chr8:20,240,000-21,170,000

# Priority 7: Garin5a - bivalent, 10kb validated
chr7:44,420,000-44,915,000

# Priority 8: Zfp997 - strong detection, 10kb validated
chr13:65,940,000-66,410,000

# Priority 9: Mbd2 region - epigenetic relevance
chr18:70,965,000-72,930,000

# Priority 10: 3830432H09Rik - good detection
chr1:119,485,000-120,000,000

# Additional gained stripes
chr13:84,145,000-85,405,000    # Tmem161b
chr17:19,865,000-20,300,000    # Fpr-rs6
chr3:144,625,000-145,530,000   # Clca4b
chr8:51,570,000-53,460,000     # 1190028D05Rik
```

---

## 4. Biological Interpretation

### 4.1 Gained Stripes and Activation Signatures

**H3K27ac enrichment at gained stripes**: 7/10 gained stripes (70%) have H3K27ac at anchor.

This is consistent with BAP1-KO biology: loss of BAP1 deubiquitinase activity leads to increased H2AK119ub1, but the downstream effects on gene activation are complex. The presence of active marks at gained stripes suggests these may represent **ectopic enhancer-promoter contacts at newly activated loci**.

Notable gained stripes with activation signatures:
- stripe_0257 (*6820431F20Rik*): H3K27ac+H3K27me3 (bivalent) with logFC=0.88
- stripe_0230 (*Garin5a*): H3K27ac+H3K27me3 (bivalent)
- stripe_0179 (*Clca4b*): H3K27ac+H3K27me3 (bivalent)

**Pattern**: 3/10 gained stripes have bivalent (H3K27ac+H3K27me3) signatures, suggesting Polycomb-poised promoters gaining stripe contacts.

### 4.2 Lost Stripes and Polycomb Associations

**H3K27me3 at lost stripes**: 3/9 lost stripes (33%) have H3K27me3 marks.

Lost stripes with Polycomb signatures:
- stripe_0056 (*Unc13d*): H3K27me3+ at promoter - best-detected stripe overall
- stripe_0183 (*Lingo2*): H3K27me3+ at promoter
- stripe_0255 (*Lrrc27*): H3K27me3+ at promoter

**Interpretation**: Loss of stripes at Polycomb-marked loci could indicate:
1. Disruption of Polycomb domain boundaries in BAP1-KO
2. Changes in cohesin-mediated loop extrusion at Polycomb-regulated regions
3. Altered chromatin compaction affecting stripe formation

### 4.3 Anchor Type Patterns

| Direction | Promoter | Active_Enhancer | Other |
|-----------|----------|-----------------|-------|
| Lost | 56% (5/9) | 11% (1/9) | 33% (3/9) |
| Gained | 60% (6/10) | 0% | 40% (4/10) |

**Observation**: The single Active_Enhancer-anchored stripe is **lost** in mutant (stripe_0128 near *Gm4719*). This could suggest enhancer-linked stripes are more sensitive to BAP1-KO perturbation.

Both lost and gained stripes are predominantly promoter-anchored (~58%), consistent with EP-stripe biology where stripes emanate from regulatory elements.

### 4.4 Interesting Gene Associations

**Metabolic genes with lost stripes:**
- *Msmo1* - cholesterol biosynthesis
- *Clybl* - citrate metabolism
- *Man1a* - glycoprotein processing

**Signaling/developmental genes:**
- *Unc13d* (lost) - synaptic vesicle exocytosis
- *Lingo2* (lost) - leucine-rich repeat protein
- *Garin5a* (gained) - GTPase regulation

**Epigenetic readers near differential stripes:**
- *Mbd2* region (gained) - methyl-CpG binding protein, relevant to chromatin regulation

---

## 5. Caveats and Confidence Assessment

### 5.1 Global FDR Concern

**All stripes have FDR > 0.43**. The lowest FDR is 0.437 (stripes 0091, 0136, 0247). This means:
- None pass conventional significance thresholds (FDR < 0.05 or even 0.10)
- Quantitative support for differential signal is weak
- Detection (source) is the primary evidence for condition-specificity

**Recommendation**: Treat this as an **exploratory analysis**. Visual validation in JuiceBox is essential before biological conclusions.

### 5.2 Direction Consistency

All stripes show consistent direction between source and logFC:
- Lost stripes: all have logFC < 0 (as expected)
- Gained stripes: all have logFC > 0 (as expected)

This is reassuring despite high FDR - the direction of effect is consistent.

### 5.3 Low Detection Confidence Stripes

Two gained stripes (stripe_0247, stripe_0257) have **low detection confidence** but the **highest logFC values** (1.00 and 0.88):

| Stripe | logFC | Detection Conf. | Quagga p-val | in_10kb |
|--------|-------|-----------------|--------------|---------|
| stripe_0247 | 1.00 | low | 9.7e-12 | FALSE |
| stripe_0257 | 0.88 | low | 2.3e-13 | FALSE |

These have excellent Quagga p-values but no 10kb validation. The strong quantitative signal (logFC) suggests they may be **real but resolution-sensitive features**. Prioritize for JuiceBox inspection.

### 5.4 Stripes Requiring Additional Scrutiny

1. **stripe_0183 (Lingo2)** - Exceptionally large anchor (3150 kb), which may indicate a merged or complex feature rather than a single stripe.

2. **stripe_0128 (Gm4719)** - Only Active_Enhancer-anchored stripe. Weak Quagga p-value (0.95) despite biological interest.

3. **stripe_0136 (Mbd2 region)** - Weak Quagga p-value (0.96) but epigenetically relevant gene proximity. Visual inspection will determine if this is signal or noise.

### 5.5 Limitations of This Dataset

1. **Sequencing depth** (~300M reads) is below optimal for stripe detection per Quagga paper benchmarks
2. **No CTCF data** available to classify CTCF vs EP stripes
3. **No H3K4me1** at early timepoint limits enhancer annotation
4. **Small total count** (19 stripes) limits statistical power for edgeR analysis

---

## 6. Recommended Validation Workflow

### Step 1: Visual Inspection in JuiceBox

Load pairs of .hic files (control merged vs mutant merged) at each priority coordinate.

For **LOST stripes**: Expect to see clear stripe in control, absent/weakened in mutant.
For **GAINED stripes**: Expect to see clear stripe in mutant, absent/weakened in control.

### Step 2: Cross-Reference with 10kb Calls

For stripes without 10kb validation (in_10kb=FALSE), check if any signal is visible at 10kb resolution.

### Step 3: Integration with ChIP-seq Tracks

Overlay H3K27ac and H3K27me3 BigWig tracks to confirm anchor annotations align with visual stripe anchors.

### Step 4: Prioritize for Downstream Analysis

After validation, create final high-confidence set for:
- Gene ontology analysis of anchor-proximal genes
- Integration with expression data (if available)
- Comparison with late timepoint results

---

## 7. Quick Reference: Top 5 Stripes per Direction

### Top 5 LOST (for control→mutant comparison)

| Rank | Stripe | Gene | Quagga p-val | Key Features |
|------|--------|------|--------------|--------------|
| 1 | stripe_0056 | *Unc13d* | 6.5e-32 | Bivalent marks, 10kb validated |
| 2 | stripe_0265 | *Msmo1* | 1.2e-10 | Metabolic gene, 10kb validated |
| 3 | stripe_0091 | *Clybl* | 7.3e-09 | Lowest FDR, 10kb validated |
| 4 | stripe_0255 | *Lrrc27* | 1.3e-08 | Bivalent, 10kb validated |
| 5 | stripe_0128 | *Gm4719* | 0.95 | Only Active_Enhancer |

### Top 5 GAINED (for mutant→control comparison)

| Rank | Stripe | Gene | logFC | Key Features |
|------|--------|------|-------|--------------|
| 1 | stripe_0247 | *Trim30e-ps1* | 1.00 | Highest logFC (2-fold) |
| 2 | stripe_0257 | *6820431F20Rik* | 0.88 | Bivalent, strong Quagga |
| 3 | stripe_0136 | *Mbd2* region | 0.43 | Epigenetic relevance |
| 4 | stripe_0074 | *Zfp997* | 0.35 | TF, 10kb validated |
| 5 | stripe_0230 | *Garin5a* | 0.20 | Bivalent, 10kb validated |

---

*Analysis generated for BAP1-KO differential stripe study*
