# Differential Chromatin Stripe Analysis: Context Document

## Document Purpose

This document summarizes all experimental context, data availability, analytical decisions, and pipeline architecture for differential chromatin stripe analysis in the BAP1-KO mouse cerebellum Hi-C project. It is intended to provide complete context for implementation of the analysis pipeline.

**This document contains only information and decisions that have been explicitly discussed and confirmed.** No new analytical decisions or recommendations are included.

---

## 1. Experimental Context

### 1.1 Biological System

| Parameter | Value |
|-----------|-------|
| **Organism** | Mouse |
| **Reference Genome** | mm10 |
| **Tissue** | Cerebellum |
| **Perturbation** | BAP1 knockout (KO) |
| **Control** | Wildtype littermates |

### 1.2 BAP1 Biology

BAP1 (BRCA1-Associated Protein 1) is a Polycomb regulator that functions as a deubiquitinase for H2AK119ub1. Key biological considerations:

- Loss of BAP1 may affect Polycomb-repressed domains
- Potential effects on enhancer-promoter contacts at derepressed loci
- Expected changes in chromatin architecture at Polycomb domain boundaries
- Possible emergence of ectopic EP-stripes at newly activated enhancers

### 1.3 Experimental Design

**Samples per condition per timepoint:**
- Control (wildtype): 3 biological replicates (ctrl_M1, ctrl_M2, ctrl_M3)
- Mutant (BAP1-KO): 3 biological replicates (mut_M1, mut_M2, mut_M3)
- Merged files: Pooled reads from all 3 replicates per condition

**Timepoints:**
- Early timepoint (250831)
- Late timepoint (250402)

**Sequencing depth:** ~300 million reads per individual replicate (slightly below the ~500M ideal identified in the Quagga paper for optimal stripe detection)

### 1.4 Available Orthogonal Data

| Timepoint | Available ChIP-seq |
|-----------|-------------------|
| Early | H3K27ac, H3K27me3 |
| Late | H3K27ac, H3K4me1 |

CTCF Cut&Run is in progress by another lab member but not yet available.

---

## 2. Stripe Calling Tool: Quagga

### 2.1 Tool Background

Quagga (Feng et al., Genome Research 2025) is a chromatin stripe detection tool that:
- Uses image processing (Gabor filtering, Gaussian blur) for candidate stripe identification
- Applies Poisson statistics for statistical verification
- Determines stripe length via p-value maximization along the stripe body
- Outputs bedpe format compatible with JuiceBox visualization

### 2.2 Stripe Biology (from Quagga paper)

**CTCF-stripes:**
- Anchor overlaps occupied CTCF motif
- Result of CTCF/cohesin-mediated loop extrusion
- Median length ~300kb
- Enriched for CTCF and RAD21 at anchor

**EP-stripes (Enhancer-Promoter):**
- CTCF-deficient at anchor
- Overlaps enhancer-promoter annotations
- Median length 41kb (E-P) to 71kb (P-P)
- Enriched for H3K4me2/3 and H3K27ac at anchor

**Indeterminate stripes:**
- Neither CTCF nor EP classification
- May represent false positives or low-mappability regions

### 2.3 Bedpe Output Format

```
#chr1    x1         x2         chr2    y1         y2         P_val
chr1     50015000   50060000   chr1    50060000   50645000   0.0138
```

**Column interpretation:**
- Columns 1-3 (chr1, x1, x2): Stripe anchor (narrow width, where CTCF/cohesin or EP anchor sits)
- Columns 4-6 (chr2, y1, y2): Stripe span (long length, the region the stripe extends into)
- Column 7: P-value from Poisson statistics (detection confidence, NOT differential)

**Stripe direction inference:**
- When y1 ≥ x2: Span extends downstream (3' stripe)
- When y2 ≤ x1: Span extends upstream (5' stripe)
- Direction inferred by comparing anchor width (x2-x1) vs span length (y2-y1)

### 2.4 Key Finding from Paper on Sequencing Depth

From Quagga paper: "Hi-C data stripe detection effective at 250 million filtered reads or more" with optimal detection at higher depths. The paper's benchmarking was done on datasets with billions of reads; detection at ~300M reads (our depth) will be less sensitive.

---

## 3. Quagga Parameters: Decision to Rerun

### 3.1 Original Parameters Used

```python
HIC_PARAMS_ORIGINAL = {
    "reference_genome": "mm10",
    "norm": "balanced",
    "threshold": 0.15,
    "resolution": 5000,
    "max_range": 7500000,
    "min_length": 500000,      # 500kb - PROBLEMATIC
    "min_distance": 500000,    # 500kb
    "max_width": 50000,
    "window_size": 7,
    "nstrata_blank": 10,
    "sigma": 1,
    "rel_height": 0.3,
}
```

### 3.2 Problem Identified

The `min_length=500kb` parameter was taken from the README "Recommended Hi-C settings" but is inconsistent with what the paper actually showed:
- Paper Figure 2A shows Quagga detecting stripes as short as 50kb
- README example code uses `min_length=200kb`
- Paper's EP-stripe analysis found median lengths of 41-71kb
- 500kb minimum excludes ALL EP-stripes and many CTCF-stripes

### 3.3 Revised Parameters (Confirmed)

```python
HIC_PARAMS_REVISED = {
    "reference_genome": "mm10",
    "norm": "balanced",
    "threshold": 0.15,
    "resolution": 5000,
    "max_range": 7500000,
    "min_length": 200000,      # CHANGED: 200kb
    "min_distance": 200000,    # CHANGED: 200kb
    "max_width": 50000,
    "window_size": 7,
    "nstrata_blank": 10,
    "sigma": 1,
    "rel_height": 0.3,
}
```

**Rationale for 200kb (not lower):**
- Captures CTCF-stripes (median ~300kb, distribution extends to 200kb)
- Captures longer EP-stripes (P-P median 71kb, distribution extends higher)
- Matches the README example code
- Conservative given sequencing depth (~300M reads)
- True short EP-stripes (<100kb) are difficult to detect reliably at this depth anyway
- Downstream filtering will handle false positives

### 3.4 Expected Impact of Rerun

| Resolution | Original (500kb min) | Expected (200kb min) |
|------------|---------------------|----------------------|
| 5kb merged | 45-52 stripes | ~120-180 stripes |
| Individual 5kb | 11-20 stripes | ~30-60 stripes |

### 3.5 Resolutions to Run

**Primary analysis:** 5kb resolution
- Better anchor localization
- Sharper CTCF/histone mark peaks at anchors (per paper)

**Validation:** 10kb resolution
- Used only for confidence annotation (not for adding stripes)
- Stripes detected at both 5kb and 10kb are higher confidence

---

## 4. Initial Stripe Counts (Pre-Rerun, for Reference)

### 4.1 At 5kb Resolution (min_length=500kb)

**Early timepoint:**
```
Individual replicates:
  ctrl_M1: 11 stripes
  ctrl_M2: 12 stripes
  ctrl_M3: 14 stripes
  mut_M1: 20 stripes
  mut_M2: 15 stripes
  mut_M3: 16 stripes

Merged files:
  ctrl_merged: 52 stripes
  mut_merged: 45 stripes
```

**Late timepoint:**
```
Individual replicates:
  ctrl_M1: 5 stripes
  ctrl_M2: 11 stripes
  ctrl_M3: 12 stripes
  mut_M1: 3 stripes
  mut_M2: 7 stripes
  mut_M3: 1 stripe

Merged files:
  ctrl_merged: 37 stripes
  mut_merged: 39 stripes
```

### 4.2 At 10kb Resolution (min_length=500kb)

**Early timepoint:**
```
Merged files:
  ctrl_merged: 365 stripes
  mut_merged: 293 stripes
```

**Late timepoint:**
```
Merged files:
  ctrl_merged: 182 stripes
  mut_merged: 183 stripes
```

### 4.3 Key Observations

- Merged files have substantially more calls due to pooled reads (expected per paper)
- Individual replicates have low counts due to depth limitations
- 10kb calls are ~5-7x more numerous than 5kb calls
- 10kb likely includes more false positives (TAD boundary artifacts)

---

## 5. Differential Analysis Strategy

### 5.1 Core Challenge

**Quagga's p-value ≠ Differential p-value**

Quagga's p-value measures enrichment of the stripe region over its local flanking background within a single sample. It is a detection confidence score, not a differential metric.

For differential analysis, we need to quantify stripe signal at the same genomic regions across all samples and perform proper statistical testing.

### 5.2 Chosen Approach: Anchor-Based Quantification

**Rationale:** Similar to the established loop differential analysis workflow, adapted for stripe geometry.

**Key adaptation:** Stripes are converted from line geometry (anchor + extended span) to point-like interactions suitable for mariner-based extraction by sampling at a fixed distance into the span.

---

## 6. Pipeline Architecture

### 6.1 Phase 1: Detection & Union Set Creation

**Input:**
- ctrl_merged_stripes.bedpe (5kb, 200kb min_length - from rerun)
- mut_merged_stripes.bedpe (5kb, 200kb min_length - from rerun)
- Individual replicate bedpe files (for confidence scoring)
- 10kb merged files (for validation annotation)

**Process:**
1. Load merged 5kb stripes for each condition
2. Intersect anchors with 50kb tolerance
3. Classify each stripe by presence/absence
4. Create UNION set of all unique stripe anchors
5. Annotate with confidence information

**Stripe Classification:**
- `control_only`: Anchor in ctrl merged but not mut merged
- `mutant_only`: Anchor in mut merged but not ctrl merged  
- `shared`: Anchor in both merged files (within 50kb tolerance)

**Special case - Direction switches:**
If same anchor position but different span direction between conditions, treat as TWO SEPARATE STRIPES (not shared). These represent functionally different/dysregulated features.

**Annotations to add:**
- `source`: "control_only" | "mutant_only" | "shared"
- `pval_ctrl`: Quagga p-value if called in control (NA if not)
- `pval_mut`: Quagga p-value if called in mutant (NA if not)
- `n_ctrl_reps`: Count of individual ctrl replicates with overlapping stripe (0-3)
- `n_mut_reps`: Count of individual mut replicates with overlapping stripe (0-3)
- `in_10kb`: Boolean, whether stripe validated by 10kb calls
- `confidence`: "high" | "medium" | "low"

**Confidence scoring logic:**
```
n_reps = n_ctrl_reps (for control_only) or n_mut_reps (for mutant_only) 
         or max(n_ctrl_reps, n_mut_reps) (for shared)
pval = pval_ctrl or pval_mut (whichever is present/lower)

confidence = case_when(
  (n_reps >= 2) & in_10kb ~ "high",
  (n_reps >= 1) | in_10kb | (pval < 1e-10) ~ "medium",
  TRUE ~ "low"
)
```

**Output:** `unified_stripes.rds` - GRanges or data.frame with all stripe coordinates and annotations

### 6.2 Phase 2: Quantification via Mariner

**Input:**
- unified_stripes from Phase 1
- 6 .hic files per timepoint (ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3)

**Stripe-to-GInteractions Conversion:**

Stripes must be converted to point-like interactions for mariner extraction:
- **anchor1**: Stripe anchor region (narrow side from bedpe: chr1, x1, x2)
- **anchor2**: Sampling point at 100kb into span direction from anchor center

**Sampling distance:** 100kb into the span
- With 200kb minimum stripe length, 100kb is always within the stripe body
- Consistent with the 200kb windows used in paper's histone modification pileups
- Accounts for direction (3' vs 5' stripes)

**Mariner workflow:**
```r
# 1. Convert to GInteractions
gi <- as_ginteractions(stripe_coords, ...)

# 2. Bin to resolution
binned <- assignToBins(gi, binSize = 5000, pos1 = "center", pos2 = "center")

# 3. Create buffer regions (5×5 matrix)
buffered <- pixelsToMatrices(binned, buffer = 2)

# 4. Extract matrices from .hic files
matrices <- pullHicMatrices(
  x = buffered,
  files = hic_files,
  binSize = 5000,
  norm = "KR",              # KR normalization available and preferred
  matrix = "observed",
  onDisk = TRUE,
  h5File = "stripe_matrices.h5"
)

# 5. Aggregate to counts (sum 5×5 = 25 pixels per stripe per sample)
counts[i, j] <- sum(matrices[, , i, j])
```

**Output:** `stripe_counts.tsv` - Stripes × samples count matrix

### 6.3 Phase 3: Differential Analysis (edgeR)

**Input:**
- stripe_counts.tsv
- Sample metadata (condition assignment)

**Process:**
1. Create DGEList with count matrix
2. **Skip low-count filtering** (retain all stripes due to small total feature count)
3. TMM normalization
4. Estimate dispersion
5. Fit GLM with ~condition design
6. Extract results with topTags()

**Output columns:**
- stripe_id
- logFC (positive = higher in mutant/BAP1-KO)
- logCPM
- LR (likelihood ratio)
- PValue
- FDR

### 6.4 Phase 4: Integration & Final Output

**Merge:** Phase 1 annotations + Phase 3 statistics

**Final output columns:**
- chr, anchor_x1, anchor_x2, span_y1, span_y2 (original bedpe coordinates)
- source: control_only | mutant_only | shared
- pval_ctrl, pval_mut (Quagga detection p-values)
- confidence: high | medium | low (Phase 1 detection confidence)
- n_ctrl_reps, n_mut_reps
- in_10kb
- logFC, FDR (edgeR differential statistics)
- direction: lost | gained | unchanged | strengthened | weakened
- direction_confidence: high | medium | low (tiered confidence for direction)
- direction_consistent: TRUE | FALSE | NA (whether logFC matches source)

**Direction classification logic (tiered confidence):**

Detection (source) is the PRIMARY evidence for condition-specific stripes.
FDR and logFC direction provide confidence tiers, not classification gates.

```
# Direction assignment (all condition-specific stripes get classified)
direction = case_when(
  source == "control_only" ~ "lost",       # By detection
  source == "mutant_only" ~ "gained",      # By detection
  source == "shared" & FDR < 0.05 & logFC > 0.3 ~ "strengthened",
  source == "shared" & FDR < 0.05 & logFC < -0.3 ~ "weakened",
  source == "shared" ~ "unchanged"
)

# Confidence tier (separate column)
direction_confidence = case_when(
  # Lost confidence tiers
  source == "control_only" & FDR < 0.10 & logFC < 0 ~ "high",
  source == "control_only" & logFC < -0.2 ~ "medium",
  source == "control_only" ~ "low",
  # Gained confidence tiers
  source == "mutant_only" & FDR < 0.10 & logFC > 0 ~ "high",
  source == "mutant_only" & logFC > 0.2 ~ "medium",
  source == "mutant_only" ~ "low",
  # Shared
  source == "shared" & FDR < 0.05 ~ "high",
  source == "shared" & FDR < 0.10 ~ "medium",
  TRUE ~ "low"
)

# Directional consistency flag
direction_consistent = case_when(
  source == "control_only" & logFC < 0 ~ TRUE,
  source == "mutant_only" & logFC > 0 ~ TRUE,
  source == "shared" ~ NA,
  TRUE ~ FALSE
)
```

**Confidence tier interpretation:**
- **High**: Quantitative support (FDR < 0.10 + correct direction)
- **Medium**: Directional support (logFC > 0.2 in expected direction)
- **Low**: Detection only (may lack quantitative support)

**Note:** Low directional consistency (~50%) may indicate noisy detection
or weak biological signal. This is tracked and reported in summaries.

---

## 7. Key Design Decisions (Confirmed)

### 7.1 Replicate Usage Strategy

**Decision:** Tiered confidence with merged files as primary detection source

- Merged files provide sensitivity for stripe detection
- Individual replicates provide confidence scoring
- Single output list with confidence annotation (not separate lists)

### 7.2 Resolution Strategy

**Decision:** 5kb as primary, 10kb for validation only

- 5kb used for all detection and quantification
- 10kb stripes used only to annotate `in_10kb` confidence flag
- No union of 5kb and 10kb stripe sets

### 7.3 Anchor Overlap Definition

**Decision:** Anchor overlap only (Option A)

- Stripes with anchors overlapping within 50kb tolerance are considered "the same stripe"
- Span differences ignored for matching purposes
- Anchor is the functional unit

### 7.4 Direction Switch Handling

**Decision:** Treat as separate stripes (Option B)

- If same anchor position but different span direction between conditions, these are TWO SEPARATE STRIPES
- Not classified as "shared"
- Functionally different/dysregulated features

### 7.5 Sampling Distance for Quantification

**Decision:** 100kb into span direction

- Fixed distance (not proportional)
- Direction-aware (respects 3' vs 5' stripe orientation)
- Within stripe body for all stripes ≥200kb

### 7.6 Normalization

**Decision:** KR (Knight-Ruiz) normalization

- KR normalization is available in .hic files (computed with juicer_tools)
- Preferred over VC when available

### 7.7 Count Filtering

**Decision:** Skip low-count filtering

- Retain all stripes for edgeR analysis
- Small total feature count (~100-200 stripes) means we should not filter

### 7.8 Analysis Scope

**Decision:** Control vs mutant at each timepoint (primary)

- Primary analysis: Differential stripes within each timepoint
- Secondary (lower priority): Early→late changes within genotype
- Early↔late comparison has limited statistical power without intensive APA

---

## 8. Data Availability Summary

### 8.1 Hi-C Data

| Sample | .hic File | Available |
|--------|-----------|-----------|
| ctrl_M1 | ✓ | Both timepoints |
| ctrl_M2 | ✓ | Both timepoints |
| ctrl_M3 | ✓ | Both timepoints |
| mut_M1 | ✓ | Both timepoints |
| mut_M2 | ✓ | Both timepoints |
| mut_M3 | ✓ | Both timepoints |

**Normalization available:** KR (Knight-Ruiz), VC (Vanilla Coverage)

**Resolutions available:** 5kb, 10kb (others likely available but not used)

### 8.2 Stripe Calls (Post-Rerun)

To be generated with revised parameters (min_length=200kb):
- 6 individual replicate bedpe files per timepoint per resolution
- 2 merged bedpe files per timepoint per resolution
- Total: 16 bedpe files per timepoint

### 8.3 ChIP-seq Data

| Timepoint | H3K27ac | H3K27me3 | H3K4me1 | CTCF |
|-----------|---------|----------|---------|------|
| Early | ✓ | ✓ | ✗ | In progress |
| Late | ✓ | ✗ | ✓ | In progress |

---

## 9. Analysis Goals

### 9.1 Primary Goal

Generate a high-confidence list of differential/dysregulated stripes between BAP1-KO and wildtype control at each timepoint.

### 9.2 Downstream Applications

1. Comparison with Cut&Run CTCF data (when available)
2. Integration with ChIP-seq (H3K27ac, H3K27me3, H3K4me1) for stripe classification
3. Hypothesis generation about BAP1-KO effects on chromatin architecture
4. Identification of regulatory E-P changes

### 9.3 Validation Plan

For top differential stripes:
1. 10kb resolution consistency check
2. JuiceBox manual visual inspection
3. CTCF/ChIP-seq overlap analysis
4. Length distribution analysis for CTCF vs EP classification
5. Gene proximity analysis

---

## 10. Biological Interpretation Framework

### 10.1 Expected Patterns Given BAP1 Biology

**Loss of stripes (control-only) may indicate:**
- Disrupted cohesin recruitment at Polycomb-regulated boundaries
- Altered chromatin compaction affecting loop extrusion processivity
- Loss of Polycomb domain boundaries

**Gain of stripes (mutant-only) may indicate:**
- Ectopic EP-stripes at derepressed enhancers
- Increased H2AK119ub1 → gene activation → new enhancer-promoter contacts
- Redistributed cohesin from lost boundaries

### 10.2 Stripe Classification Strategy

Using available ChIP-seq data:
- **CTCF-stripes**: Anchor overlaps CTCF peaks (when available)
- **Active stripes**: Anchor overlaps H3K27ac peaks
- **Polycomb-associated**: Anchor overlaps H3K27me3 peaks
- **Indeterminate**: No clear overlap

---

## 11. Technical Considerations

### 11.1 Mariner Adaptation for Stripes

Mariner was designed for point-like loop interactions. Adaptation for stripes requires:

1. Parsing bedpe format to extract anchor vs span
2. Inferring stripe direction from coordinates
3. Computing sampling point (anchor center + 100kb in span direction)
4. Converting to GInteractions format

### 11.2 Memory Management

Following established loop pipeline approach:
- HDF5-backed storage for extracted matrices
- On-disk processing to avoid RAM limitations
- Chunked extraction (blockSize parameter)

### 11.3 Coordinate Systems

- Bedpe files use 0-based coordinates
- Mariner handles conversion with `starts.in.df.are.0based = TRUE`
- Ensure consistency throughout pipeline

---

## 12. Files and Scripts (Existing Infrastructure)

### 12.1 Stripe Calling Script

`scripts/call_stripes.py` - Quagga wrapper (to be updated with revised parameters)

### 12.2 Loop Analysis Pipeline (Reference)

The existing loop differential analysis pipeline provides templates for:
- Mariner-based quantification (`scripts/prep_loops.R`, `scripts/extract_counts.R`)
- edgeR differential analysis (`scripts/edgeR.R`)
- APA visualization (`scripts/apa_analysis.R`)

### 12.3 Shared Resources

- Blacklist file: `/expanse/lustre/projects/csd940/ctea/HiC/dchic/250123blacklist.bed`
- Reference genome: mm10

---

## 13. Summary of What Will Be Implemented

1. **Quagga rerun** with min_length=200kb for all samples at 5kb and 10kb resolution

2. **Phase 1 script**: Detection and union set creation
   - Load bedpe files
   - Intersect anchors
   - Classify stripes
   - Annotate confidence

3. **Phase 2 script**: Quantification via mariner
   - Convert stripes to GInteractions
   - Extract 5×5 matrices at sampling points
   - Aggregate to count matrix

4. **Phase 3 script**: edgeR differential analysis
   - TMM normalization
   - GLM fitting
   - Statistical testing

5. **Phase 4 script**: Integration
   - Merge annotations and statistics
   - Classify direction (lost/gained/etc.)
   - Output final differential stripe list

---

## 14. Not Included in This Analysis

Per explicit decisions:
- APA-based aggregate quantification (not needed, using per-stripe quantification)
- Three separate confidence lists (single list with annotation column instead)
- Extensive early↔late comparison (secondary priority)
- Union of 5kb and 10kb stripe sets (10kb for validation only)
- min_length below 200kb (too permissive given sequencing depth)

---

*Document generated: December 2024*
*Project: BAP1-KO Differential Chromatin Stripe Analysis*
*Organism: Mouse (mm10)*