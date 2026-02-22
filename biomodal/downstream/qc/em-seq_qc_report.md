# EM-seq QC Report

**Run date**: 2026-02-19
**Reference genome**: mm10 + controls (`mm10+controls.fa`)
**Pipeline**: EM-seq (Nextflow), fastp → bwameth → MarkDuplicates → MethylDackel
**Lane**: L002

## Samples

| Sample | ID | Description |
|---|---|---|
| ctrl1 | 250504_A1_Bap1_P60_ctrl1_5mC_S44_L002 | Bap1 P60 control 1, 5mC |
| ctrl2 | 250504_B1_Bap1_P60_ctrl2_5mC_S45_L002 | Bap1 P60 control 2, 5mC |

---

## 1. Sequencing & Read Quality (fastp)

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Sequencing | PE 151 + 151 cycles | PE 151 + 151 cycles |
| Total reads (pre-filter) | 31,498,130 | 25,663,168 |
| Total read pairs | 15,749,065 | 12,831,584 |
| Reads passing filter | 31,498,130 (100%) | 25,663,168 (100%) |
| Low quality reads filtered | 0 | 0 |
| Q20 rate (post-filter) | 98.34% | 98.41% |
| Q30 rate (post-filter) | 94.50% | 94.67% |
| Mean read length (post-trim) | 147 bp | 147 bp |
| GC content (pre-filter) | 22.72% | 22.61% |
| GC content (post-filter) | 22.02% | 22.04% |
| Adapter-trimmed reads | 4,217,792 (13.4%) | 3,139,472 (12.2%) |
| Adapter-trimmed bases | 124,042,802 | 85,185,230 |
| fastp estimated duplication rate | 9.05% | 8.34% |
| Insert size peak (fastp) | 192 bp | 192 bp |

### Per-base Quality (FastQC)

Both samples: **pass**. Mean quality scores remain above Q37 across the full read length. Median is Q40 at all positions. The 10th percentile drops to Q24 at positions 145-149 and Q22 at 150-151.

---

## 2. FastQC Module Summary

| Module | ctrl1 | ctrl2 |
|---|---|---|
| Basic Statistics | pass | pass |
| Per base sequence quality | pass | pass |
| Per tile sequence quality | **fail** | **fail** |
| Per sequence quality scores | pass | pass |
| Per base sequence content | pass | pass |
| Per sequence GC content | **warn** | **warn** |
| Per base N content | pass | pass |
| Sequence Length Distribution | **warn** | **warn** |
| Sequence Duplication Levels | pass | pass |
| Overrepresented sequences | pass | pass |
| Adapter Content | pass | pass |

---

## 3. Per-tile Sequence Quality

Both samples fail this module. The same tiles are affected in both samples.

### Worst tiles (deviation from per-tile mean quality)

| Tile | Base Position | ctrl1 deviation | ctrl2 deviation |
|---|---|---|---|
| 2350 | 9 | -12.17 | -12.60 |
| 2245 | 9 | -9.82 | -9.95 |
| 2441 | 8 | -5.82 | -6.24 |
| 2395 | 95-99 | -4.40 | -4.54 |
| 2295 | 95-99 | -4.29 | -4.41 |
| 2287 | 50-54 | -3.90 | -3.76 |
| 2385 | 55-59 | -3.75 | -3.79 |
| 2282 | 150-151 | -3.74 | -3.74 |
| 2351 | 80-84 | -3.68 | -3.69 |
| 2242 | 120-124 | -3.67 | -3.84 |

### Tile quality summary

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Total tiles | 782 | 782 |
| Tiles with deviation < -2 | 39 (5.0%) | 38 (4.9%) |
| Tiles with deviation < -5 | 3 | 3 |
| Tiles with deviation < -10 | 1 | 1 |

All affected tiles are on surface 2 (tile IDs 2xxx). Most frequently affected base positions: 150-151 (8 tiles), 15-19 (6 tiles), 65-69 (4 tiles), 10-14 (4 tiles).

---

## 4. Per-sequence GC Content

Both samples: **warn**. The distribution peaks at 20% GC, consistent with C→T conversion of unmethylated cytosines in the mouse genome.

### Secondary peak at high GC (ctrl1)

| GC% | Read count |
|---|---|
| 90 | 1,538 |
| 95 | 5,648 |
| 96 | 8,940 |
| 97 | 7,517 |
| 98 | 11,703 |
| 99 | 28,950 |
| 100 | 45,094 |

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Reads with GC < 10% | 705,498 (2.19%) | 517,093 (2.00%) |
| Reads with GC 10-35% | 29,953,643 (92.78%) | 24,115,505 (93.17%) |
| Reads with GC > 80% | 124,345 (0.39%) | 74,061 (0.29%) |
| Reads with GC >= 90% | 115,758 (0.36%) | 49,087 (0.19%) |
| Reads with GC 95-100% | 107,850 (0.33%) | 44,321 (0.17%) |

---

## 5. Sequence Length Distribution

Both samples: **warn** (non-uniform length distribution due to adapter trimming).

### ctrl1

| Length bin | Read count | % of total |
|---|---|---|
| 30-34 bp | 2,602 | 0.01% |
| 35-49 bp | 116,492 | 0.37% |
| 50-79 bp | 451,959 | 1.42% |
| 80-99 bp | 468,531 | 1.47% |
| 100-129 bp | 1,263,432 | 3.97% |
| 130-149 bp | 2,156,922 | 6.77% |
| 150-152 bp | 27,382,332 | 85.98% |

Full-length reads (150-152 bp): 86.0%. Trimmed reads (<150 bp): 14.0%. Very short reads (<80 bp): 571,153 (1.79%).

### ctrl2

| Length bin | Read count | % of total |
|---|---|---|
| 30-34 bp | 2,168 | 0.01% |
| 35-49 bp | 58,601 | 0.23% |
| 50-79 bp | 198,924 | 0.77% |
| 80-99 bp | 273,182 | 1.06% |
| 100-129 bp | 980,843 | 3.80% |
| 130-149 bp | 1,874,172 | 7.25% |
| 150-152 bp | 22,609,066 | 87.56% |

Full-length reads (150-152 bp): 87.6%. Trimmed reads (<150 bp): 12.4%. Very short reads (<80 bp): 259,693 (1.01%).

---

## 6. Per-base Sequence Content

Both samples: **pass**. Base composition across the read body (ctrl1, representative):

| Position | A | T | G | C |
|---|---|---|---|---|
| 1 | 41.9% | 37.7% | 10.4% | 10.1% |
| 10-14 | 39.4% | 39.1% | 10.8% | 10.7% |
| 75-79 | 38.9% | 38.9% | 11.0% | 11.2% |
| 150-151 | 38.7% | 38.7% | 11.0% | 11.5% |

The genome-wide composition is shifted from the expected ~25% per base to ~39% A+T and ~11% G+C, reflecting C→T conversion of unmethylated cytosines. Composition is stable across read positions with no slope or endpoint distortion.

---

## 7. Alignment (samtools flagstat)

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Total reads (QC-passed) | 30,936,587 | 25,400,353 |
| Total reads (QC-failed) | 903,939 | 382,979 |
| Primary reads | 30,843,338 | 25,370,444 |
| Secondary alignments | 93,249 | 29,909 |
| Supplementary alignments | 0 | 0 |
| Mapped (QC-passed) | 30,933,112 (99.99%) | 25,398,188 (99.99%) |
| Primary mapped | 30,839,863 (99.99%) | 25,368,279 (99.99%) |
| Properly paired | 30,439,756 (98.69%) | 25,224,440 (99.42%) |
| Singletons | 3,475 (0.01%) | 2,165 (0.01%) |
| Mate on different chr | 76,216 | 40,454 |
| Mate on different chr (mapQ>=5) | 29,789 | 19,729 |

---

## 8. Picard Alignment Summary Metrics

Collected with `--IS_BISULFITE_SEQUENCED true`.

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| PCT_PF_READS | 97.92% | 98.86% |
| PCT_PF_READS_ALIGNED | 99.99% | 99.99% |
| PF_MISMATCH_RATE | 20.03% | 19.97% |
| PF_HQ_ERROR_RATE | 20.04% | 19.98% |
| PF_INDEL_RATE | 0.025% | 0.025% |
| PCT_READS_ALIGNED_IN_PAIRS | 99.99% | 99.99% |
| PCT_PF_READS_IMPROPER_PAIRS | 1.30% | 0.57% |
| STRAND_BALANCE | 0.500 | 0.500 |
| PCT_CHIMERAS | 1.01% | 0.29% |
| PCT_ADAPTER | 0.0002% | 0.0002% |
| PCT_SOFTCLIP (R1) | 0.12% | 0.10% |
| PCT_SOFTCLIP (R2) | 0.57% | 0.35% |

---

## 9. Duplication (Picard MarkDuplicates)

Run with `--BARCODE_TAG RX --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500`.

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Read pairs examined | 15,744,356 | 12,828,676 |
| Unpaired reads examined | 4,191 | 2,553 |
| Read pair duplicates | 2,400,696 | 1,791,436 |
| Read pair optical duplicates | 1,747,093 | 1,306,675 |
| Unpaired read duplicates | 1,351 | 737 |
| **PERCENT_DUPLICATION** | **15.25%** | **13.97%** |
| Optical duplicate fraction | 72.78% | 72.94% |
| Non-optical (PCR) duplicate pairs | 653,603 | 484,761 |
| PCR duplication rate (estimated) | ~4.2% | ~3.8% |
| Estimated library size | 145,176,493 | 133,061,635 |

---

## 10. Insert Size

Collected with Picard `CollectInsertSizeMetrics` on MAPQ-filtered reads.

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Median insert size | 236 bp | 236 bp |
| Mode insert size | 192 bp | 192 bp |
| Mean insert size | 260.4 bp | 259.3 bp |
| Standard deviation | 113.1 bp | 107.9 bp |
| Median absolute deviation | 66 bp | 63 bp |
| Min insert size | 2 bp | 2 bp |
| FR read pairs | 13,776,654 | 11,390,487 |
| RF read pairs | 687 | 498 |
| Tandem read pairs | 6 | 1 |

---

## 11. Spike-in Controls

### Read counts (samtools idxstats)

| Control | Ref length | ctrl1 mapped | ctrl2 mapped |
|---|---|---|---|
| phage_lambda (unmethylated) | 48,502 bp | 177,961 | 136,606 |
| plasmid_puc19c (methylated) | 2,686 bp | 22,960 | 18,990 |
| phage_T4 (5hmC) | 168,920 bp | 0 | 0 |
| phage_Xp12 (5mC) | 64,272 bp | 3 | 0 |

Spike-in reads as fraction of total mapped: lambda = 0.58% (ctrl1), 0.54% (ctrl2); pUC19 = 0.07% (ctrl1), 0.07% (ctrl2).

### Nonconverted read counts

Reads flagged as nonconverted (failed per-read conversion filter):

| Control | ctrl1 nonconverted / total | ctrl2 nonconverted / total |
|---|---|---|
| phage_lambda | 1,903 / 177,961 (1.07%) | 3,105 / 136,606 (2.27%) |
| plasmid_puc19c | 857 / 22,960 (3.73%) | 1,236 / 18,990 (6.51%) |

### Autosomal nonconverted read counts (per-chromosome)

| Chromosome | ctrl1 nonconverted | ctrl2 nonconverted |
|---|---|---|
| chr1 | 21,803 | 27,425 |
| chr2 | 21,858 | 28,119 |
| chr3 | 16,753 | 21,099 |
| chr4 | 17,617 | 22,801 |
| chr5 | 17,445 | 22,593 |
| chr6 | 17,808 | 22,049 |
| chr7 | 16,031 | 20,543 |
| chr8 | 15,262 | 19,169 |
| chr9 | 15,824 | 19,540 |
| chr10 | 13,180 | 17,312 |
| chr11 | 15,029 | 19,299 |
| chr12 | 14,253 | 18,205 |
| chr13 | 12,470 | 16,573 |
| chr14 | 13,830 | 18,180 |
| chr15 | 11,885 | 15,256 |
| chr16 | 11,192 | 14,076 |
| chr17 | 10,256 | 13,336 |
| chr18 | 10,233 | 13,141 |
| chr19 | 6,550 | 8,896 |
| chrX | 5,339 | 7,247 |
| chrY | 2,187 | 3,298 |
| chrM | 162 | 175 |

---

## 12. M-bias (MethylDackel)

### CHH context, chr1, read 1 (representative positions)

| Position | ctrl1 nMeth | ctrl1 nUnmeth | ctrl1 % meth | ctrl2 nMeth | ctrl2 nUnmeth | ctrl2 % meth |
|---|---|---|---|---|---|---|
| 1 | 858 | 66,961 | 1.27% | 914 | 55,510 | 1.62% |
| 10 | 787 | 67,762 | 1.15% | 835 | 56,644 | 1.45% |
| 50 | ~596 | ~67,185 | ~0.88% | ~727 | ~55,866 | ~1.29% |
| 100 | 596 | 67,185 | 0.88% | 727 | 55,866 | 1.28% |
| 120 | 613 | 64,359 | 0.94% | 709 | 54,710 | 1.28% |

CHH methylation across read positions is stable, with no pronounced endpoint bias. Background CHH methylation (false-positive rate): ctrl1 ~0.9-1.3%, ctrl2 ~1.3-1.6%.

---

## 13. Coverage Estimate

| Metric | ctrl1 | ctrl2 |
|---|---|---|
| Primary mapped read pairs | 15,419,931 | 12,684,139 |
| Estimated unique pairs (post-dedup) | ~13,019,235 | ~10,893,203 |
| Median insert size | 236 bp | 236 bp |
| Estimated unique bases | ~3.07 Gb | ~2.57 Gb |
| Mouse genome size (mm10) | 2.73 Gb | 2.73 Gb |
| Estimated unique coverage | ~1.1x | ~0.94x |

---

## 14. Tasmanian (Strand Bias)

Per-position nucleotide content by read and strand. Data collected at all positions for reads 1 and 2. The first few positions of read 1 show elevated off-diagonal signal (e.g., position 1: Nc_t = 85,121 in ctrl1 vs Nc_c = 3,475), consistent with C→T conversion on the expected strand. No anomalous cross-contamination patterns between I/C/N (incomplete/converted/non-converted) categories were observed in the first 30 positions of either sample.

---

*Generated 2026-02-21 from pipeline output in `em-seq_output/qc/`.*
