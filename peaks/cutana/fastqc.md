Here's the comprehensive FastQC report for your CUT&RUN data:

---

## FastQC Report — CUT&RUN (PUM1, old cells)

### Experiment Overview

**FastQC version**: 0.12.1 | **Read length**: 101 bp (PE) | **Encoding**: Sanger/Illumina 1.9

You have **6 samples** (12 FASTQ files, R1+R2 each):

| Sample | Condition | Mark | Total Read Pairs | %GC |
|--------|-----------|------|----------------:|----:|
| IgG | control | IgG | 24,028,048 | 48% |
| idx25_ctrl_1 | ctrl rep1 | H3K4me3 | 9,734,722 | 53% |
| idx26_mut_1 | mut rep1 | H3K4me3 | 6,672,993 | 51% |
| idx27_ctrl_2 | ctrl rep2 | H3K4me3 | 7,866,822 | 52% |
| idx28_mut_2 | mut rep2 | H3K4me3 | 7,683,131 | 55% |
| idx29_ctrl_1 | ctrl rep1 | H3K27ac | 13,959,976 | 47–48% |

**Missing**: `index_30_mut_1_H3K27ac` has only an HTML file — no `.txt` report. The R2 file also appears absent. This sample may have failed QC or processing was incomplete.

---

### Module Summary (pass / warn / fail)

| Module | IgG | K4me3 ctrl1 | K4me3 mut1 | K4me3 ctrl2 | K4me3 mut2 | K27ac ctrl1 |
|--------|:---:|:-----------:|:----------:|:-----------:|:----------:|:-----------:|
| Base quality | PASS | PASS | PASS | PASS | PASS | PASS |
| Tile quality | warn(R2) | PASS | PASS | PASS | PASS | PASS |
| Seq quality scores | PASS | PASS | PASS | PASS | PASS | PASS |
| Base content | **FAIL** | **FAIL** | **FAIL** | **FAIL** | **FAIL** | **FAIL** |
| GC content | **FAIL** | PASS | warn | warn | warn | PASS |
| N content | PASS | PASS | PASS | PASS | PASS | PASS |
| Length dist. | PASS | PASS | PASS | PASS | PASS | PASS |
| Duplication | **FAIL** | PASS | PASS | PASS | PASS | PASS |
| Overrepresented | PASS | PASS | PASS | warn | warn | warn |
| Adapter content | PASS | warn | warn | warn | warn | **FAIL** |

---

### Key Findings

#### 1. Sequence Quality — Excellent
All samples have mean Phred scores of **35–36+** across all positions (both R1 and R2). Medians are 37 across the board. No quality trimming is needed.

#### 2. Per Base Sequence Content — FAIL (all samples, expected)
Every sample fails this module. This is **normal and expected for CUT&RUN/CUT&Tag** data — the Tn5 transposase has a known insertion sequence bias in the first ~10 bp of reads. Not a concern.

#### 3. IgG Control — High Duplication (major concern)
- **Only 37–41% of reads are unique** (60%+ duplicated)
- Duplication breakdown: ~18% of sequences at >100x, ~22% at >500x
- GC content also fails, suggesting either:
  - **Over-sequencing of a very low-complexity library** (most likely — IgG controls in CUT&RUN often yield low library complexity)
  - Possible contamination
- At 24M reads, you have far more depth than needed for an IgG control, but the effective unique reads are only ~9M

#### 4. GC Content
- **IgG**: FAIL — bimodal or skewed GC distribution, consistent with low-complexity library
- **H3K4me3 samples**: warn for 3 of 4 — slight GC bias, likely reflecting enrichment at GC-rich promoter regions (expected for H3K4me3)
- **H3K27ac**: PASS — looks clean, GC ~47–48% (close to genome average)

#### 5. Adapter Contamination
- **H3K27ac ctrl_1 (idx29)**: **FAIL** — Illumina Universal Adapter detected at increasing levels along the read. Also significant PolyG signal in R2 (up to ~23%), indicating **short insert sizes** where reads run through the insert into adapters and then into poly-G (NextSeq/NovaSeq 2-color chemistry artifact)
- **All H3K4me3 samples**: warn — low-level adapter presence (TruSeq adapters)
- **Recommendation**: Run adapter trimming (e.g., `trim_galore` or `fastp`) before alignment, especially for the H3K27ac sample

#### 6. Overrepresented Sequences
- idx27, idx28, idx29 R1 reads flag TruSeq adapter sequences at 0.1–0.18% — minor, will be resolved by adapter trimming

#### 7. Read Depth
| Sample | Read pairs | Assessment |
|--------|-----------|------------|
| IgG | 24.0M (but ~9M unique) | Adequate for control |
| H3K4me3 samples | 6.7M–9.7M | On the lower side. CUT&RUN typically needs 5–8M for narrow marks, so you're borderline OK. May be tight for peak calling in the lower-depth mut_1 sample (6.7M) |
| H3K27ac ctrl_1 | 14.0M | Good depth for a broad mark |
| H3K27ac mut_1 | **Missing/incomplete** | Cannot assess |

---

### Recommendations

1. **Adapter trimming is essential** — Use `fastp` or `trim_galore` on all samples before alignment. This is critical for the H3K27ac sample (adapter fail + PolyG contamination in R2).

2. **Investigate the missing idx30 (mut_1 H3K27ac)** — You only have the HTML report for R1. Without this sample, you have no mutant replicate for H3K27ac.

3. **Mark duplicates during alignment** but consider whether to remove them — for CUT&RUN at these depths, PCR duplicates in the IgG should be removed, but for the enrichment samples, be cautious about removing "duplicates" that may be genuine enrichment signal.

4. **IgG control quality is poor** — The high duplication + GC bias means this IgG control will have limited utility for background subtraction. Consider whether you can use it, or if you'd be better off using a no-antibody control or IgG from a different batch.

5. **H3K4me3 mut_1 (idx26) depth** — At 6.7M reads, this is your shallowest sample. Monitor peak calling sensitivity for this sample vs. the others.