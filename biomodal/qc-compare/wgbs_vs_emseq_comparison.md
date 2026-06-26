# WGBS vs EM-seq Pilot QC Comparison

Shallow-sequencing pilot of 4 adult BAP1 mouse cerebellum samples processed through both WGBS and EM-seq.

## Samples

| Mouse ID | Condition | Genotype |
|----------|-----------|----------|
| 3139 | ctrl | Bap1 f/f |
| 3200 | ctrl | Bap1 f/f |
| 3201 | mut | Bap1 f/f; Math1-Cre |
| 3391 | mut | Bap1 f/f; Math1-Cre |

Read length: 101 bp paired-end (both methods).

---

## Main Comparison

| Sample | Cond | Method | Total Reads | Q30 % | % Mapped | % Dup | Mean Cov | CpG % | Non-CG % |
|--------|------|--------|-------------|-------|----------|-------|----------|-------|----------|
| 3139 | ctrl | WGBS | 33,510,766 | 95.1 | 72.4 | 18.1 | 0.65x | 81.0 | 1.34 |
| 3139 | ctrl | EM-seq | 29,613,538 | 96.3 | 99.98 | 16.5 | ~0.90x | 78.4 | 1.39 |
| 3200 | ctrl | WGBS | 27,508,574 | 95.9 | 70.2 | 16.3 | 0.53x | 81.3 | 1.06 |
| 3200 | ctrl | EM-seq | 28,099,972 | 96.5 | 99.98 | 15.5 | ~0.87x | 78.2 | 1.31 |
| 3201 | mut | WGBS | 8,676,962 | 96.5 | 70.2 | 15.9 | 0.16x | 81.6 | 1.25 |
| 3201 | mut | EM-seq | 26,139,688 | 96.8 | 99.99 | 13.5 | ~0.82x | 78.4 | 1.34 |
| 3391 | mut | WGBS | 27,126,754 | 95.3 | 77.5 | 13.5 | 0.59x | 79.3 | 1.14 |
| 3391 | mut | EM-seq | 22,488,164 | 96.4 | 99.98 | 13.7 | ~0.71x | 77.8 | 1.38 |

See [Footnotes](#footnotes) for how each metric was computed and methodology differences between pipelines.

---

## Methylation by Context

| Sample | Cond | Method | CpG % | CHG % | CHH % | Non-CG % |
|--------|------|--------|-------|-------|-------|----------|
| 3139 | ctrl | WGBS | 81.0 | 1.2 | 1.4 | 1.34 |
| 3139 | ctrl | EM-seq | 78.4 | 1.10 | 1.48 | 1.39 |
| 3200 | ctrl | WGBS | 81.3 | 0.9 | 1.1 | 1.06 |
| 3200 | ctrl | EM-seq | 78.2 | 1.08 | 1.38 | 1.31 |
| 3201 | mut | WGBS | 81.6 | 1.1 | 1.3 | 1.25 |
| 3201 | mut | EM-seq | 78.4 | 1.11 | 1.41 | 1.34 |
| 3391 | mut | WGBS | 79.3 | 1.0 | 1.2 | 1.14 |
| 3391 | mut | EM-seq | 77.8 | 1.13 | 1.46 | 1.38 |

Non-CG % = (CHG_methylated + CHH_methylated) / (CHG_total + CHH_total), computed from raw methylated/unmethylated counts.

---

## Unique Aligned Read Pairs

Unique pairs = mapped pairs minus duplicates. % of raw = unique pairs / raw read pairs.

| Sample | Cond | WGBS Unique Pairs | EM-seq Unique Pairs | WGBS % of Raw | EM-seq % of Raw |
|--------|------|-------------------|---------------------|---------------|-----------------|
| 3139 | ctrl | 9,063,573 | 12,363,107 | 54.1% | 83.5% |
| 3200 | ctrl | 7,444,658 | 11,872,717 | 54.1% | 84.5% |
| 3201 | mut | 2,303,072 | 11,299,878 | 53.1% | 86.5% |
| 3391 | mut | 8,312,135 | 9,700,154 | 61.3% | 86.3% |

WGBS unique pairs from Bismark deduplication. EM-seq unique pairs = Picard READ_PAIRS_EXAMINED − READ_PAIR_DUPLICATES.

---

## Insert Size

| Sample | WGBS Median (bp) | EM-seq Median (bp) |
|--------|-------------------|---------------------|
| 3139 | 207 | 253 |
| 3200 | 165 | 235 |
| 3201 | 157 | 220 |
| 3391 | 201 | 227 |

WGBS from Qualimap. EM-seq from Picard InsertSizeMetrics.

---

## EM-seq Conversion Controls

Nonconverted read counts from spike-in references:

| Control | Purpose | 3139 | 3200 | 3201 | 3391 |
|---------|---------|------|------|------|------|
| phage_lambda | Unmethylated DNA | 42 | 18 | 27 | 17 |
| plasmid_puc19c | CpG-methylated DNA | 5 | 6 | 3 | 7 |
| phage_T4 | ghmC DNA | 0 | 0 | 0 | 0 |
| phage_Xp12 | 5mC DNA | 0 | 0 | 0 | 0 |

These are counts of reads flagged as having too many unconverted cytosines per spike-in reference. WGBS pipeline does not include spike-in controls; non-CG methylation rate serves as the conversion proxy.

Source: `emseq/stats/nonconverted_counts/260516_*.nonconverted_counts.for_agg.tsv`

---

## Footnotes

**Total Reads**: Both mates combined (R1 + R2), before alignment. WGBS from Trim Galore logs (R1 count × 2). EM-seq from fastp JSON (`total_reads`).

**Q30 %**: WGBS = % of reads whose per-read mean quality is ≥ Q30, computed from FastQC per-sequence quality score distributions (R1 + R2 weighted). EM-seq = % of individual bases with quality ≥ Q30, from fastp. These are different metrics from different tools. WGBS Q20 is 99.7% for all samples; EM-seq Q20 is 98.9–99.1%.

**% Mapped**: WGBS = Bismark unique pair alignment rate (uniquely-mapped pairs / trimmed input pairs). Uses Bowtie2 in bisulfite mode. EM-seq = primary reads mapped / total primary reads, from samtools flagstat. Uses bwameth (BWA-MEM backend). Different aligners with different stringency; the gap reflects both chemistry and aligner differences.

**% Dup**: WGBS = Bismark deduplicate_bismark (% of aligned pairs removed). EM-seq = Picard MarkDuplicates PERCENT_DUPLICATION (OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500, BARCODE_TAG=RX).

**Mean Cov**: WGBS = Qualimap `mean_coverage` on deduplicated BAMs (mapped_bases / genome_size). EM-seq = estimated as (primary_mapped − dup_reads) × 100 bp / 2.73 Gbp; prefixed with ~ because no coverage tool was run on the EM-seq BAMs.

**CpG %**: WGBS = Bismark methylation_extractor (post-dedup, M-bias corrected: `ignore_r2=2`, `ignore_3prime_r2=2`). EM-seq = summed from MethylDackel `combined_mbias.tsv` across all autosomal + sex chromosomes (excluding chrM, phage, plasmid controls), all read positions, no M-bias trimming. The ~2–3% gap between methods is partly from this M-bias correction difference (EM-seq end-positions tend to depress CpG values by ~1–2%), plus aligner/caller differences.

**Non-CG %**: WGBS = (CHG_meth + CHH_meth) / (CHG_total + CHH_total) from Bismark methylation_extractor counts. EM-seq = same formula from MethylDackel mbias counts. In mouse cerebellum, non-CG includes genuine neuronal mCA plus residual non-conversion.

---

## Source Files

All paths relative to `qc-compare/`.

### WGBS

| Metric | File |
|--------|------|
| Total reads | `wgbs/trimgalore/logs/*.trimming_report.txt` |
| Q30 | `wgbs/fastqc/zips/*_fastqc.zip` (per-sequence quality scores section) |
| % Mapped | `wgbs/bismark/alignments/logs/*_bismark_bt2_PE_report.txt` |
| % Dup | `wgbs/bismark/deduplicated/logs/*.deduplication_report.txt` |
| Mean coverage | `wgbs/multiqc/bismark/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt` |
| CpG/CHG/CHH % | `wgbs/multiqc/bismark/multiqc_data/multiqc_bismark_methextract.txt` |
| Non-CG | Derived from methextract raw counts (same file) |
| Insert size | Qualimap (same multiqc_qualimap file, `median_insert_size` column) |
| Unique pairs | Bismark deduplication reports (same dedup file) |

### EM-seq

| Metric | File |
|--------|------|
| Total reads | `emseq/fastp/260516_*.fastp.json` |
| Q30 | Same fastp JSON (`q30_rate` field, post-filtering) |
| % Mapped | `emseq/stats/flagstats/260516_*.flagstat` |
| % Dup | `emseq/stats/markdups/260516_*.markdups_log` (Picard metrics table) |
| Mean coverage | Estimated from flagstat + markdups |
| CpG/CHG/CHH % | `emseq/methylDackelExtracts/mbias/260516_*.combined_mbias.tsv` |
| Non-CG | Derived from MethylDackel mbias counts (same file) |
| Insert size | `emseq/stats/insert_size/260516_*.insertsize_metrics` (Picard) |
| Unique pairs | Derived from Picard markdups (READ_PAIRS_EXAMINED − READ_PAIR_DUPLICATES) |
| Conversion controls | `emseq/stats/nonconverted_counts/260516_*.nonconverted_counts.for_agg.tsv` |

---

## Pipeline Details

### WGBS
- **Pipeline**: nf-core/methylseq v4.2.0
- **Aligner**: Bismark v0.25.1 (Bowtie2 backend, directional, `--no-mixed --no-discordant --dovetail --maxins 500`)
- **Trimmer**: Trim Galore 0.6.10 / Cutadapt 4.9 (Q20 cutoff, min length 20 bp)
- **Dedup**: Bismark deduplicate_bismark
- **QC**: Qualimap 2.3, FastQC, MultiQC
- **Reference**: mm10 (UCSC)
- **Run date**: 2026-06-25

### EM-seq
- **Pipeline**: Custom Nextflow (nf-core/methylseq fork for EM-seq)
- **Aligner**: bwameth (BWA-MEM backend)
- **Trimmer**: fastp (adapter trimming only; no quality or length filtering applied — 0 reads removed)
- **Dedup**: Picard MarkDuplicates (BARCODE_TAG=RX, OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500)
- **Methylation caller**: MethylDackel (extract + mbias)
- **QC**: samtools flagstat, Picard CollectAlignmentSummaryMetrics, Picard CollectInsertSizeMetrics, FastQC (post-dedup), Tasmanian
- **Reference**: mm10 + spike-in sequences (phage_lambda, phage_T4, phage_Xp12, plasmid_puc19c)
- **Run date**: 2026-06-25 / 2026-06-26
