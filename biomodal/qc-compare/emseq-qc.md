Here's a comprehensive summary of the EM-seq pilot run based on the available QC data.

---

## EM-seq Pilot Run Summary

### Samples

6 samples across 2 batches, all adult cerebellum:

| Sample              | Batch        | Condition | Mouse ID  |
| ------------------- | ------------ | --------- | --------- |
| 250504_A1_ctrl1     | 1 (May 2025) | ctrl      | P60 ctrl1 |
| 250504_B1_ctrl2     | 1 (May 2025) | ctrl      | P60 ctrl2 |
| 260516_A1_ctrl_3139 | 2 (May 2026) | ctrl      | 3139      |
| 260516_B1_mut_3391  | 2 (May 2026) | **mut**   | 3391      |
| 260516_C1_ctrl_3200 | 2 (May 2026) | ctrl      | 3200      |
| 260516_D1_mut_3201  | 2 (May 2026) | **mut**   | 3201      |

So you have **4 ctrl + 2 mut** — batch 1 is ctrl-only. The 260516 samples use the Math1-Cre naming convention.

### Pipeline Status

Ran June 25-26 on Expanse. **Core pipeline completed successfully** for all 6 samples:
- fastp trimming: 6/6 COMPLETED
- bwameth alignment: 44/44 chunks COMPLETED (split into 6-8 chunks per sample)
- merge + MarkDuplicates: 6/6 COMPLETED
- combine_nonconverted_counts: 6/6 COMPLETED
- **gc_bias: ABORTED** — only started for 1 sample, then pipeline stopped

The `gc_bias` abort means you won't have GC bias metrics, but all the aligned, deduplicated BAMs exist on Expanse. The `markduped_bams/` directory here is empty because the BAMs weren't rsync'd down.

### Sequencing Depth

| Sample              | Read Pairs (post-align) | ~Total Reads | Target Met? |
| ------------------- | ----------------------- | ------------ | ----------- |
| 250504_A1_ctrl1     | 15.7M                   | ~31.5M       | Yes         |
| 250504_B1_ctrl2     | 12.8M                   | ~25.7M       | Yes         |
| 260516_A1_ctrl_3139 | 14.8M                   | ~29.6M       | Yes         |
| 260516_B1_mut_3391  | 11.2M                   | ~22.5M       | Yes         |
| 260516_C1_ctrl_3200 | 14.0M                   | ~28.1M       | Yes         |
| 260516_D1_mut_3201  | 13.1M                   | ~26.1M       | Yes         |

All samples in the 22-32M read range — consistent with your ~24M shallow-seq target.

### Duplication Rates

| Sample              | % Dup     | Optical Dups | Non-optical Dups | Optical % of Total Dups |
| ------------------- | --------- | ------------ | ---------------- | ----------------------- |
| 250504_A1_ctrl1     | **15.3%** | 1.75M        | 0.65M            | 72.8%                   |
| 250504_B1_ctrl2     | **14.0%** | 1.31M        | 0.48M            | 72.9%                   |
| 260516_A1_ctrl_3139 | **16.5%** | 1.37M        | 1.07M            | 56.1%                   |
| 260516_B1_mut_3391  | **13.7%** | 0.80M        | 0.74M            | 52.2%                   |
| 260516_C1_ctrl_3200 | **15.5%** | 1.20M        | 0.97M            | 55.4%                   |
| 260516_D1_mut_3201  | **13.5%** | 0.83M        | 0.93M            | 47.3%                   |

Duplication rates of **13.5-16.5%** are quite good for EM-seq at this sequencing depth. For reference, WGBS/EM-seq pilot data under 30M reads typically runs 10-20% dups, so you're right in the sweet spot.

### Library Complexity

| Sample              | Estimated Library Size | Implication                                        |
| ------------------- | ---------------------- | -------------------------------------------------- |
| 250504_A1_ctrl1     | **145M**               | Excellent — deep-seq would yield high unique reads |
| 250504_B1_ctrl2     | **133M**               | Excellent                                          |
| 260516_A1_ctrl_3139 | **80M**                | Good                                               |
| 260516_B1_mut_3391  | **71M**                | Good                                               |
| 260516_C1_ctrl_3200 | **81M**                | Good                                               |
| 260516_D1_mut_3201  | **76M**                | Good                                               |

**Batch 1 has ~1.8x the library complexity of batch 2.** Both batches have plenty of headroom for deeper sequencing — at 25M reads you're sampling a fraction of either library — but batch 2 will start saturating (diminishing returns on unique reads) sooner. At typical deep-seq depths (~400M reads for your DUET), the 70-80M library sizes in batch 2 would hit high duplication (~80%+). At a more moderate 100-150M reads, they'd be fine.

### Notable Observations

**Batch differences in supplementary/secondary reads:**
- Batch 1 (250504): 120K-343K secondary/supplementary reads
- Batch 2 (260516): 12K-33K secondary/supplementary reads

Batch 1 has 10-25x more supplementary alignments. This could reflect different adapter contamination, chimeric reads, or different trimming stringency. Worth checking the fastp reports on Expanse to see if adapter content differs.

**Optical duplicate fraction differs by batch:**
- Batch 1: ~73% of dups are optical
- Batch 2: ~47-56% of dups are optical

This likely reflects different flowcell positions or sequencing runs (different clustering density, different patterned flowcell behavior).

### Bottom Line

**The EM-seq worked.** All 6 samples aligned successfully with reasonable duplication rates (13.5-16.5%) and good library complexity at this shallow depth. The pipeline completed through alignment and deduplication — only the gc_bias QC step failed (likely a job timeout or resource issue, not a data quality problem).

Two things to be aware of going forward:
1. **Batch 1 is ctrl-only** — if you're planning to use this for ctrl-vs-mut DMR calling, you only have n=2 mutants (both from batch 2), which confounds batch with genotype for the mutant arm. You'd need more mut reps or batch 1 mut samples to avoid that.
2. **The fastp QC reports and nonconverted_counts files are only on Expanse** (the local symlinks don't resolve). If you want to check conversion efficiency (the key EM-seq QC metric — tells you if the enzymatic conversion actually worked), you'll need to pull those nonconverted_counts TSVs from Expanse. That's the single most important number for confirming EM-seq actually worked vs just having good alignment stats.

Want me to help you pull those conversion efficiency files from Expanse, or investigate the gc_bias abort?


/expanse/lustre/projects/csd940/zalibhai/emseq/logs

```
rsync -avhP expanse:/expanse/lustre/projects/csd940/zalibhai/emseq/logs ~/sdsc/emseq/
```

---


Got it. Here's what happened and what to pull.

## The Abort: `tasmanian`, not `gc_bias`

The pipeline actually died on the **`tasmanian` (tasmanian-mismatch)** step, not `gc_bias`. At line 1869-1873 of the SLURM log:

```
ERROR ~ Error executing process > 'tasmanian (6)'
Caused by:
  Failed to create Conda environment
    command: conda create --yes --quiet --prefix ...
      bioconda::samtools=1.22 bioconda::tasmanian-mismatch=1.0.9
    status : 120
```

Conda exited with status 120 after successfully solving the environment and listing the package plan (80.3 MB of packages). The env creation appeared to complete the verify step but then failed — likely a **disk quota or filesystem timeout** on Expanse's lustre, since status 120 isn't a standard conda error. When `tasmanian` died, Nextflow killed the in-flight `gc_bias` task and never launched the remaining QC steps. That means **none** of the downstream QC ran:

| Step                 | Status             |
| -------------------- | ------------------ |
| methylDackel_extract | 0 of 6 (never ran) |
| methylDackel_mbias   | 0 of 6 (never ran) |
| gc_bias              | killed mid-run     |
| idx_stats            | never started      |
| flagstats            | never started      |
| fastqc               | 0 of 6 (never ran) |
| picard_metrics       | 0 of 6 (never ran) |
| insert_size_metrics  | 0 of 6 (never ran) |
| multiqc              | never started      |

The fix is straightforward: pre-create the tasmanian conda env before rerunning, or resume with `-resume`. You could also try:
```bash
conda create -n tasmanian bioconda::samtools=1.22 bioconda::tasmanian-mismatch=1.0.9
```
...in a separate job first to see if it's a quota/space issue.

## Pulling the files

The conversion efficiency data (nonconverted_counts) and fastp QC reports are on Expanse but not locally. Here's the rsync command to pull everything you need:

```bash
# From your Mac — pull the key QC files that the pipeline DID produce
rsync -avP zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/work_main/conda/ /dev/null  # skip this, just checking

# Pull the fastp JSON reports (6 files)
rsync -avP --include='*/' --include='*.fastp.json' --exclude='*' \
  zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/work_main/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/emseq/fastp_reports/

# Pull the nonconverted_counts TSVs (the combined ones from combine_nonconverted_counts)
# These are the final merged per-sample files, not the per-chunk ones
# Check where they landed first:
ssh zalibhai@login.expanse.sdsc.edu "find /expanse/lustre/projects/csd940/zalibhai/emseq/ -name '*.nonconverted_counts.tsv' -not -path '*/work_main/*' | head -20"
```

Actually, the cleanest approach is to find where the pipeline wrote its final outputs. The Nextflow `publishDir` would have copied results somewhere. Run this on Expanse to locate them:

```bash
# On Expanse — find all pipeline outputs
find /expanse/lustre/projects/csd940/zalibhai/emseq/ -maxdepth 3 -type d | grep -v work_main | sort

# Then rsync the whole output directory down
rsync -avP zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/output/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/emseq/output/
```

If the pipeline used a `publishDir` that isn't obvious, the fastp JSONs are in the work dirs. You can grab them directly:

```bash
# Grab fastp reports from the work dirs (6 files, small)
for hash in 80/504c448d d8/21e344f2 d4/17437fab d2/1bd0bd79 0d/5720b486 19/2629034d; do
  rsync -avP "zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/work_main/${hash}*/*.fastp.json" \
    /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/emseq/fastp_reports/
done
```

The conversion efficiency numbers should be in those nonconverted_counts TSVs — they tell you what fraction of reads had unconverted cytosines, which is the critical EM-seq QC metric.

## Source of the Batch Difference Stats

All from the **markdup log files** at `emseq/stats/markdups/*.markdups_log`. Specifically, the `## METRICS CLASS picard.sam.DuplicationMetrics` table in each file. The key columns I compared:

**`SECONDARY_OR_SUPPLEMENTARY_RDS`** (line with the metrics data in each file):
- `250504_A1_ctrl1.markdups_log`: **343,432**
- `250504_B1_ctrl2.markdups_log`: **120,874**
- `260516_A1_ctrl_3139.markdups_log`: **33,182**
- `260516_B1_mut_3391.markdups_log`: **27,412**
- `260516_C1_ctrl_3200.markdups_log`: **16,083**
- `260516_D1_mut_3201.markdups_log`: **12,626**

**`ESTIMATED_LIBRARY_SIZE`** (same metrics line):
- Batch 1: **145M** and **133M**
- Batch 2: **80M**, **71M**, **81M**, **76M**

**`READ_PAIR_OPTICAL_DUPLICATES`** / total duplicates (same line, derived):
- Batch 1: 72.8-72.9% optical
- Batch 2: 47.3-56.1% optical

These are all in the single metrics row after `## METRICS CLASS picard.sam.DuplicationMetrics` in each `.markdups_log` file. The supplementary reads difference (10-25x more in batch 1) could reflect different adapter contamination, library prep method differences, or different trimming behavior between the two batches.

---

Sure — here are the batch differences from the markdup logs side by side:

### Read Depth & Mapping

| Sample              | Batch | Condition | Read Pairs | Suppl/Secondary Reads | Unmapped |
| ------------------- | ----- | --------- | ---------- | --------------------- | -------- |
| 250504_A1_ctrl1     | 1     | ctrl      | 15,744,356 | 343,432               | 4,191    |
| 250504_B1_ctrl2     | 1     | ctrl      | 12,828,676 | 120,874               | 2,553    |
| 260516_A1_ctrl_3139 | 2     | ctrl      | 14,798,647 | 33,182                | 5,613    |
| 260516_B1_mut_3391  | 2     | mut       | 11,237,657 | 27,412                | 3,852    |
| 260516_C1_ctrl_3200 | 2     | ctrl      | 14,042,341 | 16,083                | 4,902    |
| 260516_D1_mut_3201  | 2     | mut       | 13,063,574 | 12,626                | 3,771    |

### Duplication & Library Complexity

| Sample              | Batch | % Dup  | Total Pair Dups | Optical Dups | % Optical | Est. Library Size |
| ------------------- | ----- | ------ | --------------- | ------------ | --------- | ----------------- |
| 250504_A1_ctrl1     | 1     | 15.25% | 2,400,696       | 1,747,093    | 72.8%     | **145,176,493**   |
| 250504_B1_ctrl2     | 1     | 13.97% | 1,791,436       | 1,306,675    | 72.9%     | **133,061,635**   |
| 260516_A1_ctrl_3139 | 2     | 16.46% | 2,435,540       | 1,365,944    | 56.1%     | 79,808,785        |
| 260516_B1_mut_3391  | 2     | 13.68% | 1,537,503       | 802,541      | 52.2%     | 70,558,966        |
| 260516_C1_ctrl_3200 | 2     | 15.45% | 2,169,624       | 1,202,560    | 55.4%     | 80,901,440        |
| 260516_D1_mut_3201  | 2     | 13.50% | 1,763,696       | 833,538      | 47.3%     | 76,271,788        |

### Batch Summary

| Metric                | Batch 1 (250504) | Batch 2 (260516) | Fold Difference             |
| --------------------- | ---------------- | ---------------- | --------------------------- |
| Suppl/secondary reads | 120K–343K        | 13K–33K          | **~10-25x more in batch 1** |
| Est. library size     | 133–145M         | 71–81M           | **~1.8x larger in batch 1** |
| Optical dup fraction  | 72.8–72.9%       | 47.3–56.1%       | **~1.3x higher in batch 1** |
| Overall dup rate      | 14.0–15.3%       | 13.5–16.5%       | Comparable                  |


---

You're right, I apologize for the scattered response. Let me start clean.

First, run this on Expanse to see where the pipeline wrote its outputs:

```bash
find /expanse/lustre/projects/csd940/zalibhai/emseq/ -maxdepth 3 -type d | grep -v work_main | sort
```

Send me the output and I'll give you the exact rsync command.


---

The output directory is `em-seq_output/` and it has a `stats/nonconverted_counts/` subdirectory plus `fastp/`. Run this to pull the QC files you need:

```bash
rsync -avhP zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/em-seq_output/stats/nonconverted_counts/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/emseq/stats/nonconverted_counts/
```

That gets the conversion efficiency data. If you also want the fastp JSONs (read quality, adapter content, trimming stats):

```bash
rsync -avP --include='*.fastp.json' --exclude='*' \
  zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/em-seq_output/fastp/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/emseq/fastp_reports/
```

Send me the nonconverted_counts files once they're down and I'll parse the conversion efficiency for you.


```
rsync -avhPL --include='*.fastp.json' --exclude='*' zalibhai@login.expanse.sdsc.edu:/expanse/lustre/projects/csd940/zalibhai/emseq/em-seq_output/fastp/ \
  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/emseq/fastp_reports/
```

---

Here's what the data shows.

## Fastp Read Quality

| Sample              | Batch | Raw Reads | Read Length | Q20   | Q30   | GC%   | Adapter-Trimmed |
| ------------------- | ----- | --------- | ----------- | ----- | ----- | ----- | --------------- |
| 250504_A1_ctrl1     | 1     | 31.5M     | **151 bp**  | 98.2% | 94.1% | 22.0% | 4.2M (13.4%)    |
| 250504_B1_ctrl2     | 1     | 25.7M     | **151 bp**  | 98.3% | 94.4% | 22.0% | 3.1M (12.2%)    |
| 260516_A1_ctrl_3139 | 2     | 29.6M     | **101 bp**  | 98.9% | 96.3% | 22.0% | 0.7M (2.4%)     |
| 260516_B1_mut_3391  | 2     | 22.5M     | **101 bp**  | 98.9% | 96.4% | 22.2% | 0.6M (2.5%)     |
| 260516_C1_ctrl_3200 | 2     | 28.1M     | **101 bp**  | 99.0% | 96.5% | 22.1% | 0.7M (2.5%)     |
| 260516_D1_mut_3201  | 2     | 26.1M     | **101 bp**  | 99.1% | 96.8% | 22.2% | 0.7M (2.5%)     |

The read lengths differ: **batch 1 = 151bp, batch 2 = 101bp**. That explains the higher adapter trimming and supplementary reads in batch 1 — longer reads read into adapter more often. GC content is 22% across all samples, which is the expected EM-seq signature (C→T conversion tanks the GC). Quality is excellent everywhere, slightly better in batch 2.

## Conversion Efficiency (Spike-in Controls)

This is the critical metric. The nonconverted_counts files report reads with too many unconverted cytosines per reference sequence. The spike-in controls tell you if the enzymatic conversion worked:

| Control            | Purpose                                | 250504_A1 | 250504_B1 | 260516_A1 | 260516_B1 | 260516_C1 | 260516_D1 |
| ------------------ | -------------------------------------- | --------- | --------- | --------- | --------- | --------- | --------- |
| **phage_lambda**   | Unmethylated (should convert)          | **1,903** | **3,105** | 42        | 17        | 18        | 27        |
| **plasmid_puc19c** | Unmethylated CpG-rich (should convert) | **857**   | **1,236** | 5         | 7         | 6         | 3         |
| phage_T4           | ghmC (protected)                       | 0         | 0         | 0         | 0         | 0         | 0         |
| phage_Xp12         | 5mC (should convert in EM-seq)         | 0         | 0         | 0         | 0         | 0         | 0         |
|                    |                                        |           |           |           |           |           |           |

**This is a major batch difference.** Batch 1 (250504) has **~50-100x more nonconverted reads** on the unmethylated controls compared to batch 2 (260516). Lambda nonconverted: 1,903-3,105 in batch 1 vs 17-42 in batch 2. Same pattern for puc19c.

This means **batch 1 had substantially worse conversion efficiency** — more reads retained cytosines where they should have been fully converted. Batch 2's conversion looks excellent (near-zero lambda failures).

To calculate the actual conversion rate you'd need the total reads mapping to lambda (not just the nonconverted ones) — that would come from the flagstats/idxstats output, which didn't run because of the tasmanian crash. But the relative difference is unambiguous.

## Autosomal Nonconverted Counts

For genomic chromosomes, the nonconverted counts are similar across batches (scaling with read depth), which is expected — those reflect real biological methylation being detected, not conversion failures:

| Chromosome | 250504_A1 | 250504_B1 | 260516_A1 | 260516_B1 | 260516_C1 | 260516_D1 |
| ---------- | --------- | --------- | --------- | --------- | --------- | --------- |
| chr1       | 21,803    | 27,425    | 20,481    | 18,219    | 18,667    | 21,135    |
| chr2       | 21,858    | 28,119    | 20,830    | 17,927    | 18,823    | 20,443    |
| chrM       | 162       | 175       | 136       | 127       | 62        | 89        |

## Summary

The EM-seq **worked** — conversion looks great in batch 2. Batch 1 has noticeably worse conversion on the spike-in controls, likely due to different kit lots, protocol versions, or library prep conditions (also different read lengths and higher adapter contamination). If you're combining these batches for downstream analysis, the conversion efficiency difference is something to account for — batch 1 will have slightly inflated methylation calls from unconverted reads.

---

Actually, **stop creating that named conda env** — Nextflow won't use it. Nextflow manages its own conda environments as prefix-based installs inside `work_main/conda/env-{hash}/`, not named envs. A `conda create -n tasmanian ...` goes to a completely different location.

Here's what actually happened and what to do:

## The Root Cause

The `tasmanian` process told Nextflow to `conda create --prefix /expanse/.../work_main/conda/env-4247555d91b9a90f55c383f15a5231f0`. The conda output shows it solved the env, prepared the transaction, and verified it — but there's no "Executing transaction: done" line. That means package extraction failed mid-way, likely from a **disk quota or lustre filesystem issue**. Exit code 120 isn't a standard conda error, which supports this.

The env directory at that prefix path now probably exists in a **broken/incomplete state**. If you just `-resume`, Nextflow might see the directory and assume the env is ready, then fail when trying to use it. Or `conda create` will refuse because the prefix already exists.

## What To Do (3 steps)

**Step 1:** Check your disk quota on Expanse to rule out space issues:
```bash
expanse-client user -p
```

**Step 2:** Remove the broken conda env prefix:
```bash
rm -rf /expanse/lustre/projects/csd940/zalibhai/emseq/work_main/conda/env-4247555d91b9a90f55c383f15a5231f0
```

**Step 3:** Resume the pipeline. You need to add `-resume` to the Nextflow command for the main pipeline (Step 2 in your wrapper). Nextflow will skip all 56 completed tasks (fastp, alignReads, mergeAndMarkDuplicates, combine_nonconverted_counts) from cache and only run the QC steps that never completed (methylDackel, gc_bias, flagstats, fastqc, picard_metrics, insert_size_metrics, tasmanian, multiqc).

I don't have your SLURM submission script — can you show me the file you ran with `sbatch`? I need to see the exact Nextflow command so I can tell you where to add `-resume` without breaking anything. It should be somewhere in `/expanse/lustre/projects/csd940/zalibhai/emseq/scripts/` or wherever you launched from.

---

Great question. The 4 batch-2 samples are the **same mice** in both experiments, so this is a pure chemistry-vs-chemistry comparison at similar shallow depths.

## Head-to-Head: Same 4 Samples

### Raw Reads → Usable Aligned Data

| Sample    | WGBS Raw Pairs | WGBS Align Rate | WGBS Dedup Reads | EM-seq Raw Pairs | EM-seq Align Rate | EM-seq Dedup Pairs |
| --------- | -------------- | --------------- | ---------------- | ---------------- | ----------------- | ------------------ |
| ctrl_3139 | 16.8M          | 72.4%           | 9.1M             | 14.8M            | ~99.9%            | ~12.4M             |
| ctrl_3200 | 13.8M          | 70.2%           | 7.4M             | 14.0M            | ~99.9%            | ~11.8M             |
| mut_3201  | **4.3M**       | 70.2%           | **2.3M**         | 13.1M            | ~99.9%            | ~11.3M             |
| mut_3391  | 13.6M          | 77.5%           | 8.3M             | 11.2M            | ~99.9%            | ~9.7M              |

The mut_3201 WGBS failure is striking — same mouse, same tissue, but 4.3M read pairs in WGBS vs 13.1M in EM-seq. That's a library prep or loading problem in the WGBS run, not biology.

### Efficiency: What fraction of your sequencing becomes usable data?

| Metric                                  | WGBS (Bismark) | EM-seq (bwameth) |
| --------------------------------------- | -------------- | ---------------- |
| Alignment rate                          | 70–78%         | ~99.9%           |
| Dup rate                                | 13.5–18.1%     | 13.5–16.5%       |
| **Usable fraction** (aligned × non-dup) | **~59–67%**    | **~83–86%**      |

**EM-seq converts ~25-30% more of your raw reads into unique aligned data.** At shallow depths where every read counts, that's the difference between borderline and usable.

### Conversion Efficiency

| Metric                          | WGBS     | EM-seq (batch 2)              |
| ------------------------------- | -------- | ----------------------------- |
| CHG methylation (should be ~0%) | 0.9–1.2% | N/A (MethylDackel didn't run) |
| CHH methylation (should be ~0%) | 1.1–1.4% | N/A                           |
| Lambda nonconverted reads       | N/A      | 17–42 (excellent)             |
| puc19c nonconverted reads       | N/A      | 3–7 (excellent)               |

WGBS conversion is good (CHG/CHH ~1% = ~99% conversion). EM-seq batch 2 spike-in controls look excellent (near-zero lambda failures), though without MethylDackel extraction we can't get the equivalent per-base metric. Both are adequate; EM-seq likely edges it slightly based on the spike-in counts.

### Estimated Coverage

| Sample    | WGBS Coverage | EM-seq Coverage (estimated) |
| --------- | ------------- | --------------------------- |
| ctrl_3139 | 0.65x         | ~0.93x                      |
| ctrl_3200 | 0.53x         | ~0.88x                      |
| mut_3201  | **0.16x**     | ~0.85x                      |
| mut_3391  | 0.59x         | ~0.73x                      |

EM-seq gives roughly **1.4–5x better coverage** per sample (5x for mut_3201 because WGBS essentially lost that sample).

### CpG Methylation (biology check)

| Metric          | WGBS       | EM-seq                                       |
| --------------- | ---------- | -------------------------------------------- |
| CpG methylation | 79.3–81.6% | Unknown (pipeline crashed before extraction) |

WGBS confirmed biologically correct CpG levels (~80%, expected for mammalian somatic tissue). EM-seq can't confirm yet — that's the one thing you're missing. Once you resume the pipeline and MethylDackel runs, you'll get this number.

## The Verdict

**EM-seq is clearly better at low read counts.** The reasons are fundamental to the chemistry:

1. **Bisulfite treatment destroys DNA.** High-temperature denaturation with sodium bisulfite causes depurination and strand breaks. This produces damaged fragments that won't align, which is why WGBS throws away ~25-30% of reads at the alignment step. EM-seq uses TET2 + APOBEC3A at mild conditions, preserving DNA integrity.

2. **More reads → more coverage → more power.** At ~24M reads, EM-seq gives ~0.8-0.9x coverage vs WGBS's ~0.5-0.65x. That's ~50-70% more coverage per sequencing dollar.

3. **No sample dropout.** All 4 EM-seq samples performed consistently. WGBS lost mut_3201 to severe undersequencing, killing your n=2 mutant arm.

4. **Dup rates are equivalent** (~14-17% for both), so library complexity isn't driving the difference — it's purely alignment efficiency from gentler chemistry.

One caveat: part of the alignment rate gap is the aligner (bwameth/BWA-MEM is more sensitive than Bismark/Bowtie2). But published benchmarks consistently show EM-seq gains 10-15% mapping rate over WGBS even with the same aligner. The chemistry difference is real.

**Bottom line:** if you're choosing one method to scale up for production sequencing, EM-seq will give you substantially more usable data per dollar, and the batch 2 conversion quality looks excellent. The one thing you need to confirm is that CpG methylation levels match the WGBS results once MethylDackel runs — if they do, it's an open-and-shut case for EM-seq.


---


Good — plenty of SU allocation, so it's not a quota issue. The conda env creation failure was likely a transient lustre filesystem hiccup.

Here's what to do. Modify the script to add `-resume` to both Nextflow commands, and skip Step 1 (already done). Run this:

```bash
cat > /expanse/lustre/projects/csd940/zalibhai/emseq/scripts/resume_emseq.sb << 'SBEOF'
#!/bin/bash
#SBATCH --job-name=emseq_resume
#SBATCH --output=logs/emseq_resume_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --account=csd940

echo "======================================"
echo "EM-seq Pipeline Resume"
echo "======================================"
echo "Started: $(date)"

source ~/.bashrc
conda activate emseq_env

WORK_DIR="/expanse/lustre/projects/csd940/zalibhai/emseq"
PIPELINE_DIR="${WORK_DIR}/EM-seq"

cd "${PIPELINE_DIR}"

# Step 1 already completed — skip entirely

# Step 2: Resume main pipeline (skips 56 completed tasks from cache)
echo "=== Resuming main EM-seq pipeline ==="
echo "Started: $(date)"

nextflow run main.nf \
  -resume \
  -profile conda \
  --genome 'mm10' \
  --ubam_dir "${PIPELINE_DIR}/ubam" \
  --email zalibhai@ucsd.edu \
  --flowcell 'L001' \
  --outputDir "${WORK_DIR}/em-seq_output" \
  -work-dir "${WORK_DIR}/work_main"

echo ""
echo "======================================"
echo "EM-seq Pipeline Complete"
echo "======================================"
echo "Finished: $(date)"
SBEOF
```

Then before submitting, clean up the broken conda env:

```bash
rm -rf /expanse/lustre/projects/csd940/zalibhai/emseq/work_main/conda/env-4247555d91b9a90f55c383f15a5231f0
```

Then submit:

```bash
sbatch /expanse/lustre/projects/csd940/zalibhai/emseq/scripts/resume_emseq.sb
```

Nextflow's `-resume` will hash-match all 56 completed tasks (fastp, alignReads, mergeAndMarkDuplicates, combine_nonconverted_counts) and skip straight to the QC steps that never ran. Should take well under an hour since the heavy alignment work is already cached.