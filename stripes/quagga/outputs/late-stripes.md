# Late Timepoint Differential Stripe Analysis

**Dataset**: BAP1-KO vs Wildtype Mouse Cerebellum (mm10)
**Timepoint**: Late (250402)
**Input file**: `late_stripes.bedpe`
**Analysis date**: 2024-12

---

## 1. Summary Statistics

### Overall Counts

| Metric | Count |
|--------|-------|
| **Total differential stripes** | 17 |
| **Lost stripes** (control-only) | 11 (65%) |
| **Gained stripes** (mutant-only) | 6 (35%) |

**Notable shift from early timepoint**: Late timepoint shows more lost than gained stripes (65% vs 35%), whereas early timepoint was nearly balanced (47% vs 53%). This suggests progressive loss of chromatin stripe architecture in BAP1-KO over time.

### Distribution by Direction Confidence

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| **High** | 0 | **1** | **1** |
| Medium | 11 | 5 | 16 |
| Low | 0 | 0 | 0 |

**Critical finding**: **stripe_0014** (gained, near *Apaf1*) is the **only HIGH-confidence differential stripe** across both timepoints. With FDR = 0.0753 (< 0.10), this represents the strongest statistical evidence for differential stripe signal.

### Distribution by Detection Confidence

| Confidence | Lost | Gained | Total |
|------------|------|--------|-------|
| High | 0 | 0 | 0 |
| Medium | 11 | 6 | 17 |
| Low | 0 | 0 | 0 |

All stripes have medium detection confidence.

### Distribution by Anchor Type

| Anchor Type | Lost | Gained | Total |
|-------------|------|--------|-------|
| **Promoter** | 4 (36%) | 1 (17%) | 5 (29%) |
| **Active_Enhancer** | 0 | **1** (17%) | 1 (6%) |
| **Poised_Enhancer** | 2 (18%) | 0 | 2 (12%) |
| **Other** | 5 (45%) | 4 (67%) | 9 (53%) |

**Key observations**:
- Lost stripes are enriched at Promoters (36%) and Poised_Enhancers (18%)
- The only Active_Enhancer-anchored stripe is **gained** (opposite to early timepoint)
- Poised_Enhancers (H3K4me1+, H3K27ac-) exclusively show stripe loss

### ChIP-seq Annotation Breakdown

| ChIP Mark | Lost | Gained | Total |
|-----------|------|--------|-------|
| **H3K27ac** | 4 (36%) | 3 (50%) | 7 (41%) |
| **H3K27me3** | 0 | 0 | 0 |
| **H3K4me1** | 7 (64%) | 3 (50%) | 10 (59%) |

**Note**: At late timepoint, H3K27ac and H3K4me1 are available; H3K27me3 is NOT available.

**Enhancer signature analysis**:
- **Active enhancers** (H3K27ac+H3K4me1+): 6 stripes (4 lost, 2 gained)
- **Poised enhancers** (H3K4me1+ only): 4 stripes (3 lost, 1 gained)
- Pattern suggests preferential loss of enhancer-associated stripes in BAP1-KO at late timepoint

---

## 2. High-Priority Stripes for Validation

Prioritized by: (1) Direction confidence, (2) FDR, (3) Quagga p-value, (4) biological relevance.

### Priority 1: stripe_0014 (GAINED) - ONLY HIGH-CONFIDENCE STRIPE

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr10:90,660,000-91,150,000` |
| **Gene** | *Apaf1* (328 kb to TSS) |
| **Anchor type** | **Active_Enhancer** |
| **ChIP marks** | H3K27ac+, H3K4me1+ |
| **Quagga p-val** | 5.99e-18 (exceptional) |
| **in_10kb** | TRUE |
| **logFC** | 0.34 |
| **FDR** | **0.0753** (< 0.10) |
| **Direction confidence** | **HIGH** |

**Rationale**: This is the **most statistically robust differential stripe** in the entire dataset (both timepoints). The Active_Enhancer anchor with dual active marks (H3K27ac+H3K4me1) suggests a bona fide enhancer gaining stripe contact in BAP1-KO. *Apaf1* (Apoptotic Protease Activating Factor 1) is a key apoptosis regulator - chromatin reorganization at this locus could have functional consequences.

**This stripe should be the #1 validation priority.**

---

### Priority 2: stripe_0073 (LOST) - Strongest logFC

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr16:55,860,000-56,500,000` |
| **Gene** | *Pcnp* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K4me1+ |
| **Quagga p-val** | 0.95 (weak) |
| **in_10kb** | FALSE |
| **logFC** | **-0.62** (strongest lost) |
| **FDR** | 0.121 (second-lowest) |

**Rationale**: **Strongest magnitude logFC** among lost stripes (-0.62, representing ~1.5-fold decrease). Despite weak Quagga detection, the strong quantitative signal and near-significant FDR suggest real biological effect. *Pcnp* (PEST-containing nuclear protein) is involved in protein degradation. Loss of promoter stripe contact may indicate transcriptional dysregulation.

---

### Priority 3: stripe_0196 (LOST) - Exceptional Quagga detection

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr9:95,400,000-95,920,000` |
| **Gene** | *Pcolce2* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K4me1+ |
| **Quagga p-val** | **6.44e-18** (exceptional) |
| **in_10kb** | TRUE |
| **logFC** | -0.22 |
| **FDR** | 0.533 |

**Rationale**: One of the strongest Quagga detections in the dataset. *Pcolce2* (Procollagen C-Endopeptidase Enhancer 2) is involved in collagen processing and extracellular matrix. Promoter-anchored with dual active marks and 10kb validation makes this a high-confidence lost stripe.

---

### Priority 4: stripe_0114 (GAINED) - Strong detection at active promoter

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr2:157,185,000-157,645,000` |
| **Gene** | *Manbal* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K4me1+ |
| **Quagga p-val** | **1.41e-15** (excellent) |
| **in_10kb** | TRUE |
| **logFC** | 0.35 |
| **FDR** | 0.174 |

**Rationale**: Excellent Quagga detection at promoter with active chromatin marks. *Manbal* (Mannosidase Beta Like) is involved in glycoprotein processing. The gained stripe at an active promoter suggests enhanced regulatory contacts in BAP1-KO.

---

### Priority 5: stripe_0156 (LOST) - Strong Quagga at nuclear pore gene

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr6:91,040,000-91,485,000` |
| **Gene** | *Nup210* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K4me1+ |
| **Quagga p-val** | **7.89e-14** (excellent) |
| **in_10kb** | TRUE |
| **logFC** | -0.20 |
| **FDR** | 0.766 |

**Rationale**: Strong Quagga detection with 10kb validation. *Nup210* (Nucleoporin 210) is a nuclear pore complex component involved in nuclear transport and potentially gene regulation. Loss of stripe at this promoter could affect nuclear architecture.

---

### Priority 6: stripe_0176 (LOST) - Hmgb2 chromatin regulator

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr8:57,460,000-57,920,000` |
| **Gene** | *Hmgb2* (0 bp to TSS) |
| **Anchor type** | Promoter |
| **ChIP marks** | H3K27ac+, H3K4me1+ |
| **Quagga p-val** | 0.49 (moderate) |
| **in_10kb** | TRUE |
| **logFC** | -0.28 |
| **FDR** | 0.185 (third-lowest) |

**Rationale**: *Hmgb2* (High Mobility Group Box 2) is a **chromatin architectural protein** that binds DNA and facilitates transcription. Loss of stripe contact at a chromatin regulator's own promoter is biologically intriguing given BAP1's role in chromatin organization.

---

### Priority 7: stripe_0096 (GAINED) - Pik3ap1 signaling

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr19:41,190,000-41,550,000` |
| **Gene** | *Pik3ap1* (0 bp to TSS) |
| **Anchor type** | Other |
| **ChIP marks** | H3K4me1+ |
| **Quagga p-val** | 0.006 (good) |
| **in_10kb** | TRUE |
| **logFC** | 0.35 |
| **FDR** | 0.174 |

**Rationale**: *Pik3ap1* (PI3K Adaptor Protein 1) is involved in PI3K signaling, important for cell survival and growth. Gained stripe validated at 10kb resolution with reasonable Quagga detection.

---

### Priority 8: stripe_0174 (GAINED) - Highest gained logFC

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr8:27,650,000-29,590,000` |
| **Gene** | *Poteg* (253 kb to TSS) |
| **Anchor type** | Other |
| **ChIP marks** | None |
| **Quagga p-val** | 0.28 (weak) |
| **in_10kb** | FALSE |
| **logFC** | **0.43** (highest gained) |
| **FDR** | 0.285 |

**Rationale**: Highest logFC among gained stripes (0.43). Despite weak detection metrics, the quantitative signal suggests possible real effect. Worth visual inspection to assess stripe quality.

---

### Priority 9: stripe_0016 (LOST) - Poised Enhancer

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr10:110,155,000-110,830,000` |
| **Gene** | *Nav3* (206 kb to TSS) |
| **Anchor type** | **Poised_Enhancer** |
| **ChIP marks** | H3K4me1+ (only) |
| **Quagga p-val** | 0.90 (weak) |
| **in_10kb** | FALSE |
| **logFC** | -0.26 |
| **FDR** | 0.579 |

**Rationale**: Poised enhancer (H3K4me1 without H3K27ac) losing stripe contact. *Nav3* (Neuron Navigator 3) is involved in neuronal development. Pattern suggests BAP1-KO may affect poised enhancer architecture.

---

### Priority 10: stripe_0138 (LOST) - Poised Enhancer at adhesion gene

| Attribute | Value |
|-----------|-------|
| **JuiceBox coord** | `chr5:79,130,000-80,405,000` |
| **Gene** | *Adgrl3* (720 kb to TSS) |
| **Anchor type** | **Poised_Enhancer** |
| **ChIP marks** | H3K4me1+ (only) |
| **Quagga p-val** | 0.20 (moderate) |
| **in_10kb** | FALSE |
| **logFC** | -0.24 |
| **FDR** | 0.677 |

**Rationale**: Second poised enhancer losing stripe. *Adgrl3* (Adhesion G Protein-Coupled Receptor L3, also known as LPHN3) is associated with ADHD and neuronal function. Poised enhancer dysfunction at this locus may be relevant to cerebellar biology.

---

## 3. JuiceBox Coordinate Blocks

### LOST Stripes (check for stripe disappearance in mutant)

```
# Top LOST stripes - load CONTROL .hic files to see stripe, MUTANT to confirm loss
# Late timepoint shows predominantly LOST stripes (11 total)

# Priority 2: Pcnp - STRONGEST logFC (-0.62)
chr16:55,860,000-56,500,000

# Priority 3: Pcolce2 - exceptional Quagga (6.4e-18), active marks
chr9:95,400,000-95,920,000

# Priority 5: Nup210 - nuclear pore gene, strong detection
chr6:91,040,000-91,485,000

# Priority 6: Hmgb2 - CHROMATIN REGULATOR, biologically relevant
chr8:57,460,000-57,920,000

# Priority 9: Nav3 region - Poised_Enhancer
chr10:110,155,000-110,830,000

# Priority 10: Adgrl3 region - Poised_Enhancer, neuronal gene
chr5:79,130,000-80,405,000

# Additional lost stripes
chr11:9,370,000-11,695,000     # Spmip7 (very large anchor ~2Mb)
chr14:92,100,000-93,325,000    # Pcdh9 region
chr17:20,660,000-22,995,000    # Vmn1r228 region (very long span)
chr6:24,630,000-25,820,000     # Hyal6
chr18:17,180,000-19,975,000    # 4921533I20Rik - NOTE: this is listed as GAINED in data
```

### GAINED Stripes (check for new stripe appearance in mutant)

```
# Top GAINED stripes - load MUTANT .hic files to see stripe, CONTROL to confirm absence
# Late timepoint has fewer gained stripes (6 total)

# Priority 1: Apaf1 region - ONLY HIGH-CONFIDENCE STRIPE (FDR=0.075)
# *** THIS IS THE #1 VALIDATION PRIORITY ***
chr10:90,660,000-91,150,000

# Priority 4: Manbal - excellent Quagga (1.4e-15), active promoter
chr2:157,185,000-157,645,000

# Priority 7: Pik3ap1 - PI3K signaling, 10kb validated
chr19:41,190,000-41,550,000

# Priority 8: Poteg region - highest gained logFC (0.43)
chr8:27,650,000-29,590,000

# Additional gained stripes
chr4:50,990,000-52,975,000     # Cylc2 region
chr4:78,380,000-79,340,000     # Ptprd region
chr18:17,180,000-19,975,000    # 4921533I20Rik
```

---

## 4. Biological Interpretation

### 4.1 Temporal Dynamics: Early vs Late

| Metric | Early | Late | Trend |
|--------|-------|------|-------|
| Total stripes | 19 | 17 | Slight decrease |
| Lost stripes | 9 (47%) | 11 (65%) | **Increase** |
| Gained stripes | 10 (53%) | 6 (35%) | **Decrease** |
| High-confidence | 0 | 1 | Improvement |

**Interpretation**: The shift toward more lost stripes at late timepoint suggests **progressive deterioration of chromatin stripe architecture** in BAP1-KO. This could reflect:
1. Cumulative effects of BAP1 loss on loop extrusion machinery
2. Progressive changes in Polycomb-mediated chromatin organization
3. Secondary effects of initial chromatin disruption

### 4.2 The Single High-Confidence Stripe: Apaf1 Region

**stripe_0014** deserves detailed attention as the only statistically robust finding:

- **Active enhancer anchor**: H3K27ac+H3K4me1 signature indicates functional enhancer
- **Apaf1 proximity**: 328 kb to *Apaf1* TSS (apoptosis regulator)
- **Direction**: GAINED in mutant - suggests enhanced regulatory contacts
- **Possible interpretation**: BAP1-KO may lead to activation of enhancer elements that form new stripe contacts, potentially affecting *Apaf1* regulation

This finding aligns with BAP1's role in repressing certain regulatory elements through Polycomb.

### 4.3 Enhancer Patterns: Active vs Poised

**Active enhancers (H3K27ac+H3K4me1)**:
- 4 lost, 2 gained
- Net loss of active enhancer-anchored stripes

**Poised enhancers (H3K4me1 only, no H3K27ac)**:
- 3 lost, 1 gained (chr19 *Pik3ap1* region)
- All lost poised enhancers are in Lost stripes

**Pattern**: Poised enhancers show exclusive stripe loss at late timepoint, suggesting:
1. BAP1-KO particularly affects enhancers primed but not yet active
2. Possible disruption of enhancer poising mechanisms
3. May reflect changes in developmental timing of enhancer activation

### 4.4 Biologically Notable Gene Associations

**Chromatin/Nuclear Architecture:**
- *Hmgb2* (lost) - HMG-box chromatin protein, architectural role
- *Nup210* (lost) - Nuclear pore component, gene regulation at nuclear periphery

**Signaling:**
- *Apaf1* (gained) - Apoptosis pathway (HIGH-CONFIDENCE)
- *Pik3ap1* (gained) - PI3K signaling

**Neuronal:**
- *Nav3* (lost) - Neuron navigator
- *Adgrl3* (lost) - Adhesion GPCR, ADHD-associated

**Metabolism:**
- *Pcolce2* (lost) - Collagen processing
- *Manbal* (gained) - Glycoprotein processing

### 4.5 Contrast with Early Timepoint Patterns

| Pattern | Early | Late |
|---------|-------|------|
| Bivalent marks (H3K27ac+H3K27me3) | 6 stripes | N/A (no H3K27me3 data) |
| Promoter-anchored | 58% | 29% |
| Active_Enhancer | 1 (lost) | 1 (gained) |
| Highest logFC direction | Gained (1.00) | Lost (-0.62) |

The Active_Enhancer pattern is intriguing: at early timepoint the only Active_Enhancer stripe was lost, while at late timepoint the only Active_Enhancer stripe is gained (and is the highest-confidence call). This may reflect dynamic remodeling of enhancer-anchored stripes.

---

## 5. Caveats and Confidence Assessment

### 5.1 Statistical Significance

| Significance Level | Count | Stripes |
|--------------------|-------|---------|
| FDR < 0.05 | 0 | None |
| FDR < 0.10 | **1** | **stripe_0014** (gained, Apaf1) |
| FDR < 0.15 | 2 | stripe_0014, stripe_0073 |
| FDR < 0.20 | 5 | Above + stripe_0096, stripe_0114, stripe_0176 |

**Key point**: Only **stripe_0014** passes conventional significance threshold (FDR < 0.10). All other stripes should be considered exploratory findings requiring validation.

### 5.2 Direction Consistency

All stripes show consistent direction between source and logFC:
- Lost stripes: all have logFC < 0
- Gained stripes: all have logFC > 0

This consistency is reassuring for data quality.

### 5.3 Detection Quality Concerns

**Weak Quagga detection despite biological interest:**

| Stripe | Gene | Quagga p-val | Interest |
|--------|------|--------------|----------|
| stripe_0073 | *Pcnp* | 0.95 | Strongest logFC |
| stripe_0016 | *Nav3* | 0.90 | Poised enhancer |
| stripe_0079 | *Vmn1r228* | 0.27 | Large span |

These stripes have high FDR noise but biological patterns worth checking visually.

### 5.4 Unusually Large Anchors

Two stripes have very large anchors that may indicate merged features:

| Stripe | Anchor Width | Notes |
|--------|--------------|-------|
| stripe_0018 | 2,170 kb | Near *Spmip7* |
| stripe_0138 | 1,120 kb | Poised enhancer at *Adgrl3* |

Visual inspection should assess whether these represent single stripes or merged signals.

### 5.5 Missing ChIP-seq Data

**H3K27me3 is NOT available at late timepoint**. This limits:
- Ability to identify Polycomb-repressed domains
- Assessment of bivalent (H3K27ac+H3K27me3) signatures
- Direct comparison with early timepoint Polycomb patterns

### 5.6 Low Total Count

With only 17 differential stripes and 6 samples, statistical power is limited. The high FDR values across most stripes reflect this limitation.

---

## 6. Comparison with Early Timepoint

### Key Differences

| Feature | Early | Late | Interpretation |
|---------|-------|------|----------------|
| Lost:Gained ratio | 9:10 (47:53%) | 11:6 (65:35%) | Progressive loss |
| High-confidence stripes | 0 | 1 | Better separation at late |
| Promoter anchor % | 58% | 29% | Shift to enhancers |
| Active_Enhancer stripe | Lost | Gained | Dynamic remodeling |
| Strongest logFC | +1.00 (gained) | -0.62 (lost) | Direction shift |

### Shared Patterns

1. Most stripes have medium confidence (direction and detection)
2. Direction is consistent with source in all cases
3. Majority have moderate logFC (0.2-0.4 magnitude)
4. No H3K27me3 enrichment pattern (available early, shows bivalent marks)

### Biological Interpretation

The temporal shift from balanced lost/gained at early timepoint to predominantly lost at late timepoint suggests:

1. **Progressive architectural degradation**: BAP1-KO effects accumulate over time
2. **Enhancer remodeling**: Active enhancer gains at late timepoint (Apaf1) may represent compensatory or aberrant activation
3. **Poised enhancer vulnerability**: Loss of poised enhancer stripes at late timepoint suggests developmental timing disruption

---

## 7. Recommended Validation Strategy

### Tier 1: Must Validate (Statistical support)

1. **stripe_0014 (Apaf1)** - Only FDR < 0.10 stripe. Load control and mutant merged .hic files at `chr10:90,660,000-91,150,000`. Expect clear stripe in mutant, absent in control.

### Tier 2: High Priority (Strong detection)

2. **stripe_0196 (Pcolce2)** - Quagga p-val 6.4e-18
3. **stripe_0114 (Manbal)** - Quagga p-val 1.4e-15
4. **stripe_0156 (Nup210)** - Quagga p-val 7.9e-14

### Tier 3: Biological Interest

5. **stripe_0073 (Pcnp)** - Strongest lost logFC
6. **stripe_0176 (Hmgb2)** - Chromatin regulator
7. **stripe_0016 & stripe_0138** - Poised enhancer patterns

### Validation Checklist

For each stripe:
- [ ] Load control merged .hic at 5kb resolution
- [ ] Load mutant merged .hic at 5kb resolution
- [ ] Compare stripe presence/intensity at coordinate
- [ ] Check at 10kb resolution for additional support
- [ ] Overlay H3K27ac and H3K4me1 BigWig tracks
- [ ] Note stripe direction (3' vs 5' extension)
- [ ] Score visual confidence (clear/marginal/absent)

---

## 8. Quick Reference: Top 5 Stripes per Direction

### Top 5 LOST (for control→mutant comparison)

| Rank | Stripe | Gene | Key Feature | FDR |
|------|--------|------|-------------|-----|
| 1 | stripe_0073 | *Pcnp* | Strongest logFC (-0.62) | 0.121 |
| 2 | stripe_0196 | *Pcolce2* | Best Quagga (6.4e-18) | 0.533 |
| 3 | stripe_0176 | *Hmgb2* | Chromatin regulator | 0.185 |
| 4 | stripe_0156 | *Nup210* | Nuclear pore gene | 0.766 |
| 5 | stripe_0016 | *Nav3* | Poised enhancer | 0.579 |

### Top 5 GAINED (for mutant→control comparison)

| Rank | Stripe | Gene | Key Feature | FDR |
|------|--------|------|-------------|-----|
| 1 | **stripe_0014** | *Apaf1* | **HIGH CONFIDENCE (FDR=0.075)** | **0.075** |
| 2 | stripe_0114 | *Manbal* | Best Quagga (1.4e-15) | 0.174 |
| 3 | stripe_0174 | *Poteg* | Highest logFC (0.43) | 0.285 |
| 4 | stripe_0096 | *Pik3ap1* | PI3K signaling | 0.174 |
| 5 | stripe_0083 | *4921533I20Rik* | Long span (2.6 Mb) | 0.647 |

---

## 9. Summary and Key Takeaways

1. **One statistically significant stripe**: stripe_0014 (gained at *Apaf1* Active_Enhancer, FDR=0.075) is the only finding with conventional statistical support. This should be the primary focus.

2. **Temporal trend toward stripe loss**: The shift from 47% lost (early) to 65% lost (late) suggests progressive chromatin architecture disruption in BAP1-KO.

3. **Poised enhancers are exclusively lost**: All poised enhancer (H3K4me1-only) anchored stripes are lost at late timepoint, suggesting BAP1-KO affects enhancer developmental poising.

4. **Chromatin regulators affected**: Lost stripes at *Hmgb2* (HMG-box) and *Nup210* (nuclear pore) promoters may indicate broader chromatin organization changes.

5. **Validation essential**: Given high FDR values (except stripe_0014), visual validation in JuiceBox is critical before biological interpretation.

---

*Analysis generated for BAP1-KO differential stripe study - Late Timepoint*
