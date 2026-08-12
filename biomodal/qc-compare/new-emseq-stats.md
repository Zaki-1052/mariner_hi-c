```
rsync -avhPL --exclude='*.bam' --exclude='*.bai' --exclude='*.fastq.gz' --exclude='*.methylKit.gz' expanse:/expanse/lustre/projects/csd940/zalibhai/emseq/em-seq_output/ ~/sdsc/emseq/
```

This gets all the new QC outputs (stats/fastqc, flagstats, gc_bias, idxstats, insert_size, picard_alignment_metrics, tasmanian), the updated pipeline_info (execution traces from the resumed run), the M-bias plots from methylDackelExtracts/mbias/, and the nonconverted_counts.tsv files from bwameth_align.

Excludes: BAMs (9.4G), trimmed FASTQs, and the `.methylKit.gz` files (~550MB each, 9.9G total). Everything that comes through should be well under your 50-100MB threshold.


---


Now I have everything. Here's the complete report.

---

## EM-seq Pipeline QC Report — All 6 Samples

### Overview

|                 | Batch 1 (250504) | Batch 2 (260516)                         |
| --------------- | ---------------- | ---------------------------------------- |
| **Samples**     | ctrl1, ctrl2     | ctrl_3139, mut_3391, ctrl_3200, mut_3201 |
| **Read length** | 151 bp           | 101 bp                                   |
| **Read pairs**  | 15.4M, 12.7M     | 14.7M, 11.2M, 14.0M, 13.0M               |

---

### 1. Alignment (flagstats + Picard)

| Sample    | Primary Reads | Mapped     | Properly Paired | Singletons | Chimeras |
| --------- | ------------- | ---------- | --------------- | ---------- | -------- |
| ctrl1     | 30.8M         | **99.99%** | 98.69%          | 0.01%      | 1.01%    |
| ctrl2     | 25.4M         | **99.99%** | 99.42%          | 0.01%      | 0.29%    |
| ctrl_3139 | 29.4M         | **99.98%** | 99.31%          | 0.02%      | 0.42%    |
| mut_3391  | 22.3M         | **99.98%** | 99.40%          | 0.02%      | 0.40%    |
| ctrl_3200 | 28.0M         | **99.98%** | 99.68%          | 0.02%      | 0.13%    |
| mut_3201  | 26.0M         | **99.99%** | 99.76%          | 0.01%      | 0.10%    |

All samples have outstanding alignment. For context, the WGBS experiment on the same mice achieved only 70-78% alignment — EM-seq's enzymatic conversion preserves DNA integrity far better than bisulfite.

---

### 2. Duplication & Library Complexity

| Sample    | % Duplication | Optical Duplicates | Est. Library Size |
| --------- | ------------- | ------------------ | ----------------- |
| ctrl1     | 15.25%        | 1,747,093          | 145M              |
| ctrl2     | 13.97%        | 1,306,675          | 133M              |
| ctrl_3139 | 16.46%        | 1,365,944          | 80M               |
| mut_3391  | 13.68%        | 802,541            | 71M               |
| ctrl_3200 | 15.45%        | 1,202,560          | 81M               |
| mut_3201  | 13.50%        | 833,538            | 76M               |

Duplication rates 13.5-16.5% are normal. Batch 1 has larger estimated library sizes (133-145M vs 71-81M), likely because of longer reads. Batch 2 has fewer optical duplicates, consistent with the different sequencing platform/flow cell.

---

### 3. CpG Methylation (MethylDackel M-bias — the key metric)

| Sample    | CpG %     | CHH % (conversion proxy) |
| --------- | --------- | ------------------------ |
| ctrl1     | **76.1%** | 0.99%                    |
| ctrl2     | **76.4%** | 1.40%                    |
| ctrl_3139 | **79.7%** | 1.44%                    |
| mut_3391  | **79.0%** | 1.46%                    |
| ctrl_3200 | **79.4%** | 1.38%                    |
| mut_3201  | **79.7%** | 1.43%                    |

Batch 2 CpG methylation (79.0-79.7%) matches the WGBS measurements on the same mice (79.3-81.6%). Batch 1 is ~3% lower (76%) — could be biological (different animals, different litter) or technical (151bp reads = more end effects). CHH methylation at 1.0-1.5% across all samples is normal for mammalian somatic tissue and confirms enzymatic conversion is working.

---

### 4. Spike-in Conversion Controls

| Control               | Purpose               | Batch 1          | Batch 2         | Expected | Verdict                   |
| --------------------- | --------------------- | ---------------- | --------------- | -------- | ------------------------- |
| **Lambda** (unmeth)   | Conversion efficiency | CpG = 0.5-0.8%   | CpG = 0.3-0.9%  | ~0%      | **99.1-99.7% conversion** |
| **puc19c** (CpG-meth) | Protection efficiency | CpG = 99.4-99.5% | CpG = 97.1-100% | ~100%    | **97-100% protection**    |
| **phage_T4** (ghmC)   | Should not align      | 0 reads          | 0-2 reads       | 0        | Pass                      |
| **phage_Xp12** (5mC)  | Should not align      | 0-3 reads        | 0-4 reads       | 0        | Pass                      |

Lambda conversion at 99.1-99.7% is excellent. puc19c protection at 97-100% confirms the TET2+βGT protection step is working — methylated cytosines are correctly preserved.

One batch difference: batch 1 had 50-100x more lambda spike-in reads (137-178K vs 2.3-2.8K), so the nonconverted_counts flagged earlier (1,903/3,105 nonconverted reads in batch 1 vs 17-42 in batch 2) are proportionally a very small fraction (<1% of lambda reads in batch 1). The per-CpG-site conversion rates are comparable between batches.

---

### 5. Insert Size

| Sample    | Median | Mode | MAD |
| --------- | ------ | ---- | --- |
| ctrl1     | 236    | 192  | 66  |
| ctrl2     | 236    | 192  | 63  |
| ctrl_3139 | 253    | 185  | 71  |
| mut_3391  | 227    | 166  | 65  |
| ctrl_3200 | 235    | 170  | 67  |
| mut_3201  | 220    | 178  | 54  |

All in the 220-253bp median range with tight distributions (MAD 54-71). This is exactly where EM-seq fragments should be. Batch 1 is slightly more uniform (both at 236bp median).

---

### 6. Mismatch Rate (Picard, bisulfite-aware)

All samples show ~19.3-20.0% mismatch rate (PF_MISMATCH_RATE), which is expected for enzymatically-converted data — C→T conversion registers as "mismatches" against the reference. The fact that all samples are consistent confirms uniform conversion.

---

### 7. Head-to-Head: EM-seq vs WGBS (same 4 mice)

| Metric                   | EM-seq (batch 2) | WGBS       | Winner           |
| ------------------------ | ---------------- | ---------- | ---------------- |
| Alignment rate           | **99.98%**       | 70-78%     | EM-seq (+25-30%) |
| Duplication rate         | 13.5-16.5%       | 13.5-18.1% | Comparable       |
| CpG methylation          | 79.0-79.7%       | 79.3-81.6% | Concordant       |
| CHH methylation          | 1.38-1.46%       | 1.1-1.4%   | Comparable       |
| Properly paired          | 99.3-99.8%       | N/A        | —                |
| Genome coverage per read | ~30% more usable | —          | EM-seq           |

EM-seq produces concordant methylation measurements (within 2%) with ~30% more usable data per sequencing dollar due to the dramatically better alignment rate. At ~24M reads per sample, EM-seq gives you the equivalent genomic coverage of ~34M WGBS reads.

---

### Bottom Line

**The EM-seq experiment worked.** All quality metrics are excellent:
- 99.98-99.99% alignment
- 99.1-99.7% conversion efficiency (lambda)
- 97-100% protection efficiency (puc19c)
- CpG methylation concordant with WGBS
- Clean insert sizes, low chimeras, reasonable duplication

The one batch difference worth noting: batch 1 samples show 76% CpG methylation vs 79-80% for batch 2 and WGBS. This is a ~3% offset — unclear whether biological or technical, but it's consistent within batch 1 and the samples would need to be analyzed with batch as a covariate if combined.