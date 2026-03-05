## Peak Calling Assessment

### Peak Counts & FRiP

| Sample | Peaks | FRiP | Assessment |
|---|---|---|---|
| IgG | **22** | 0.003 | Excellent negative control |
| K4me3_ctrl1 | 22,236 | 0.780 | Excellent |
| K4me3_mut1 | 20,764 | 0.707 | Excellent |
| K4me3_ctrl2 | 19,738 | 0.714 | Excellent |
| K4me3_mut2 | 22,148 | **0.868** | Excellent (notably high) |

**FRiP scores are outstanding.** ENCODE's minimum for ChIP-seq is 0.01; CUT&RUN typically exceeds that, but 0.70-0.87 is at the high end even for CUT&RUN. This confirms strong H3K4me3 enrichment across all samples.

**IgG called only 22 peaks with FRiP=0.003** — exactly what you want. The IgG is doing its job as a clean negative control.

**Peak counts (19.7K-22.2K)** are consistent across samples and in the expected range for H3K4me3 in mouse. No outlier behavior.

### Genomic Annotation — looks textbook for H3K4me3

| Annotation | Avg % of Peaks | Log2(obs/exp) | Expected? |
|---|---|---|---|
| **Promoter** | 16.3% | +3.8 to +4.0 | Yes — primary H3K4me3 target |
| **5'UTR** | 3.9% | +5.4 to +5.5 | Yes — strongest enrichment, right at TSS |
| **Exon** | 7.8% | +2.5 to +2.7 | Yes — H3K4me3 spreads into first exons |
| **Intron** | 54.4% | +0.6 to +0.7 | Neutral — introns are huge, so many peaks land there by genomic chance |
| **Intergenic** | 14.5% | **-2.0** | Yes — strongly depleted, as expected |
| **ncRNA** | 0.6% | +2.2 to +2.3 | Yes — many ncRNAs have active promoters |

The enrichment profile is exactly what H3K4me3 should look like: strong at promoters/5'UTR, depleted at intergenic regions. Highly consistent across all 4 samples — no sample is behaving as an outlier.

### Flags

1. **chrUn_JH584304 in top peaks** — All 22 IgG peaks and several top H3K4me3 peaks map to this unplaced scaffold. This is likely spurious mapping from residual E. coli DNA to a region with sequence homology on this small scaffold. Not concerning for your analysis, but when you generate final peak sets you may want to filter to canonical chromosomes (chr1-19, chrX, chrY).

2. **mut2 FRiP of 0.87 is notably higher than the other three (0.71-0.78).** Combined with its higher %GC (55%) from FastQC, this sample has the strongest enrichment signal. Peak count is consistent though (22K), so this isn't over-calling — it's likely genuinely stronger enrichment or slightly lower background. Not a problem, just worth noting for reproducibility.

3. **IgG "Uniquely Aligned Read Pairs" shows 13.7M in peak_caller_metrics** — this is ~2x the 6.9M unique pairs from alignment metrics, so this column is counting individual reads not pairs. Just a unit note for your records.

### Bottom line

Peak calling went very well. The data is clean, enrichment is strong and on-target, samples are concordant, and the IgG is a proper negative control. Ready for downstream differential analysis or consensus peak generation.

---



Good point, let me be straight about this.

**The MACS2 concern was partially valid but the results came out fine anyway.** Here's the honest breakdown:

### Where MACS2 did well
- FRiP 0.70-0.87 and 20-22K peaks per sample is genuinely excellent. No complaints about the core peak calling.
- Annotation enrichment profile is textbook H3K4me3. MACS2 clearly found real signal.
- IgG worked as a proper control despite its low complexity.

### Where SEACR might have done better (unknowable without running it)
- MACS2 assumes a Poisson/negative binomial background model calibrated on ChIP-seq noise. CUT&RUN has fundamentally lower and sparser background than ChIP-seq, which can cause MACS2 to either over-call in low-signal regions or be overly conservative depending on parameters.
- With your IgG being 54% E. coli, the effective IgG background after filtering is thin (~6.9M unique mouse reads). SEACR handles low-depth controls more gracefully since it uses the empirical signal distribution rather than modeling a noise distribution from the control.
- The chrUn_JH584304 artifact appearing in top peaks — SEACR's threshold approach may have been less susceptible to calling peaks on spurious scaffold mappings.

### The honest assessment for your PI

**For H3K4me3 specifically, MACS2 is fine.** H3K4me3 is the easiest mark to call peaks on — strong, sharp, abundant signal. You'd need a bad tool to mess it up. The results here validate that.

**The real test would be weaker or broader marks.** If you plan to run H3K27me3, H3K4me1, or especially H3K27ac CUT&RUN on this platform, the lack of SEACR (or SICER2 for broad marks) becomes a more meaningful limitation. MACS2 `--broad` mode is a poor substitute for purpose-built broad-mark callers.

**What I'd tell your PI:**
- The platform produced clean, high-quality results for this H3K4me3 test case
- MACS2-only is a real limitation if you'll be doing broad histone marks or weaker enrichments
- The E. coli carry-over in IgG and chrUn_JH584304 artifacts suggest the platform may not be filtering non-canonical chromosomes or spike-in contaminants — check if that's configurable
- The alignment metrics and QC are solid, so the upstream processing (trimming, alignment, dedup) is trustworthy
- Ask whether SEACR or SICER2 can be added as peak caller options — that would resolve the main gap

**Net: Good for sharp marks, jury's still out for your full experimental scope.**