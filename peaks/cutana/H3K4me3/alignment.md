## Alignment Metrics Summary

| Sample | Total Pairs | Unique Aligned | Unique Rate | Dup Rate | chrM | E. coli Rate |
|---|---|---|---|---|---|---|
| IgG | 23.5M | 6.9M | **29.1%** | 21.4% | 0.12% | **54.6%** |
| K4me3_ctrl1 | 9.5M | 7.6M | 80.2% | 12.5% | 0.0% | 0.01% |
| K4me3_mut1 | 6.5M | 5.3M | 81.8% | 9.2% | 0.0% | 0.01% |
| K4me3_ctrl2 | 7.7M | 6.1M | 78.8% | 11.8% | 0.01% | 0.01% |
| K4me3_mut2 | 7.5M | 6.3M | 84.0% | 10.7% | 0.0% | 0.01% |

### H3K4me3 samples — excellent

- **Unique alignment 79-84%** — very good for CUT&RUN
- **Duplication 9-12%** — notably lower than FastQC's sequence-level estimates (20-28%). Alignment-based dedup is more accurate; the libraries are higher-complexity than FastQC suggested
- **chrM ~0%** — essentially no mitochondrial contamination. This is ideal and indicates clean nuclear chromatin targeting
- **E. coli 0.01%** — negligible spike-in carry-over, meaning antibody enrichment is working very well
- **5.3-7.6M unique pairs per sample** — solid depth for H3K4me3 CUT&RUN peak calling

### IgG — the E. coli issue

The big finding: **54.6% of IgG reads align to E. coli**. This is the dominant reason for the low mouse alignment rate (29.1%) and explains the high duplication we saw in FastQC.

This is a known CUT&RUN artifact. The pA-MNase enzyme is produced in E. coli, and residual bacterial DNA carries over into the library. In H3K4me3 samples, the antibody pulls down abundant target chromatin that vastly outnumbers the E. coli contamination (0.01%). But in the IgG control, there's no specific target — the library is dominated by E. coli DNA and non-specific Tn5 cutting of accessible chromatin.

**After filtering out E. coli, you're left with ~6.9M unique mouse read pairs for IgG.** That's usable as peak-calling background, but it's worth knowing the effective depth is much lower than the 23.5M total reads suggested.

### Flags

1. **IgG usability** — 6.9M unique mouse pairs is adequate for background, but on the thin side. For SEACR (which is commonly used with CUT&RUN), this should be fine since it uses the IgG signal distribution rather than read depth per se. For MACS2, you may want to confirm coverage is sufficient.

2. **No E. coli spike-in normalization possible for H3K4me3** — at 0.01% E. coli, there are only ~500-1,000 E. coli read pairs per sample. That's far too few for reliable spike-in normalization. You'll need to use library-size or other normalization approaches instead.

3. **All samples look clean** — zero chrM, zero contamination, high alignment rates. No batch effects visible in the metrics.

---

**MACS2, not SICER2.** H3K4me3 is a narrow/sharp mark (punctate peaks at active promoters), and SICER2 is designed for broad diffuse marks (H3K27me3, H3K36me3, H3K9me3). Using SICER2 on H3K4me3 would merge nearby peaks and lose resolution.

That said, for CUT&RUN specifically, **SEACR** is worth considering over both:

| Tool | Best for | CUT&RUN suitability |
|---|---|---|
| **MACS2** | Narrow marks (H3K4me3, H3K27ac) | Good with parameter tweaks (`--nomodel`, adjusted `--extsize`) |
| **SICER2** | Broad marks (H3K27me3, H3K36me3) | Wrong mark type here |
| **SEACR** | Both narrow/broad CUT&RUN | Purpose-built for CUT&RUN's low background |

**SEACR advantages for your data:**
- Designed for CUT&RUN's sparse background (no assumptions about Poisson/negative binomial noise)
- Works directly with your IgG control or a numeric threshold
- Less sensitive to the IgG's low complexity (uses signal distribution, not per-base depth)

**MACS2 is also fine** — it's the most widely used and cited. Just needs `--nomodel --extsize 150` (or your measured fragment size) and `--keep-dup all` since you've already deduplicated.

If your PI's cloud platform offers both, you could run MACS2 narrow + SEACR and intersect for a high-confidence set. But if picking one: **SEACR** for CUT&RUN, **MACS2** if you want maximum comparability with published ChIP-seq datasets.