Now I have the full picture. Here's the landscape of downstream tools for differential non-CG analysis from EM-seq:

---

## What you have now (upstream → downstream gap)

**Upstream (done):** NEB EM-seq pipeline → MethylDackel → `.methylKit.gz` files per sample per context (CHH, CHG, CpG). These are per-site count files: chr, pos, strand, coverage, freqC, freqT. You also have the `non_cg_methylation.py` script that already does gene-body aggregation and the `ca_filter.py` for CA-context extraction.

**The gap:** No differential testing framework. Your existing code computes rates but doesn't test ctrl-vs-mut.

## Tool landscape for differential non-CG

### Tier 1: Purpose-built for this problem

**methylKit** (R/Bioconductor) — your most natural fit. It directly reads the `.methylKit.gz` files MethylDackel already produces. Key capabilities:
- `tileMethylCounts()` or `regionCounts()` to aggregate sites over custom regions (gene bodies, enhancers, TAD bins)
- `calculateDiffMeth()` uses logistic regression or Fisher's exact test per region
- Supports covariates (sex, batch) via the overdispersion-corrected logistic regression (`overdispersion="MN"`)
- Can filter by context (restrict to CA dinucleotides) if you pre-filter the methylKit input
- **Limitation:** designed for high-coverage CpG — at 1% mCA with 25x coverage, per-site tests are underpowered. The region aggregation is where it works.

**DSS** (R/Bioconductor) — beta-binomial hierarchical model, borrows information across sites. This is the strongest tool for low-signal contexts because the shrinkage estimator stabilizes noisy per-site rates:
- `DMLtest()` for site-level, `callDMR()` for regions
- Handles biological replicates natively
- No covariates in the standard version — would need to regress out sex/batch separately
- **Key advantage:** the Bayesian smoothing across nearby sites is exactly what mCA needs — individual sites are too sparse, but neighboring CA sites carry correlated signal

**dmrseq** (R/Bioconductor, built on bsseq) — two-stage approach specifically for DMR finding:
- Stage 1: identifies candidate regions by scanning for consistent methylation changes across samples
- Stage 2: permutation-based statistical test per candidate region
- Handles low coverage well because it evaluates regions, not sites
- **The best choice for genome-wide unbiased DMR discovery** without pre-defining regions
- Naturally works with non-CG contexts if you feed it CA-only input

### Tier 2: General frameworks you'd adapt

**edgeR/DESeq2 on aggregated counts** — the approach Luo et al. effectively used. For each gene body (or any region), sum methylated reads and total reads across all CA sites, then use a GLM on the ratio:
- Supports covariates (sex, batch) directly
- Well-characterized multiple testing correction
- You already have the aggregation code in `non_cg_methylation.py` — just need to output per-sample counts instead of pooled rates
- **This is what I'd recommend for the gene-body differential analysis** — it maps directly onto your DUET modality pipeline's approach and lets you compare CG and non-CG results in the same framework

**bsseq** (R/Bioconductor) — the foundational data structure. Store your per-site methylation as a BSseq object, then pipe into DSS, dmrseq, or custom analyses. Everything downstream works from this.

### Tier 3: Visualization/integration tools

**deepTools** `computeMatrix` + `plotHeatmap` — you already use this for the Hi-C/ChIP work. For mCA:
- Convert per-sample mCA bedGraphs to BigWig (you've already built this pipeline in `non_cg_analysis/bigwigs/`)
- `computeMatrix reference-point` over MeCP2 peaks, K119ub peaks, H3K27ac peaks, TAD boundaries
- Directly addresses questions 1, 3, 4 from your list

**GenomicRanges** + your existing `_shared_config.R` infrastructure — for the intersection analyses (mCA at K119ub-gained sites, at H3K27ac enhancers, etc.), you already have all the BED files loaded and the overlap framework in `_shared_config.R`.

## Mapping to your 5 questions

| Question                                | Tool                                         | Input                                                    | Approach                                                                                                                                                                                                       |
| --------------------------------------- | -------------------------------------------- | -------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1. K119ub-gained → non-CG gained?**   | GenomicRanges + edgeR                        | mCA gene-body counts + K119ub DiffBind                   | Subset genes to K119ub-gained (DiffBind FDR<0.05), test mCA differential. Fisher's overlap test. deepTools heatmap of mCA signal at K119ub-gained peaks                                                        |
| **2. Long gene-body aggregation**       | edgeR (counts) or methylKit (`regionCounts`) | Per-sample CA site counts aggregated per gene            | Weighted methylation level = Σ(meth_CA) / Σ(total_CA) per gene per sample → GLM with condition + sex + batch                                                                                                   |
| **3. TAD organization**                 | deepTools `computeMatrix` + custom R         | mCA BigWigs + TAD boundary BEDs from TADCompare          | `computeMatrix reference-point -R tad_boundaries.bed -S mCA_ctrl.bw mCA_mut.bw` → profile plot. Also: per-TAD mCA aggregation, test boundary-vs-interior                                                       |
| **4. Enhancer mCA near H3K27ac**        | GenomicRanges + edgeR                        | mCA counts in H3K27ac peak regions                       | `regionCounts()` over H3K27ac peaks (from `CHIP_PEAK_FILES`). Separately test active enhancers (H3K27ac+, >2kb from TSS)                                                                                       |
| **5. Increased hmC → increased non-CG** | Correlation analysis                         | DUET 5hmC gene-body levels + EM-seq mCA gene-body levels | Per-gene scatter: DUET hmC fold-change vs EM-seq mCA level. Test whether hmC-increased genes have higher mCA. This directly tests whether TET activity (high hmC = active TET) correlates with mCA maintenance |
| **CG comparison**                       | Merge tables                                 | DUET CG DMR results + EM-seq mCA results                 | Join on gene name. Test: do CG-hypermethylated genes also show mCA changes? Quadrant plot like your Section 05                                                                                                 |

## Recommended pipeline architecture

```
MethylDackel .methylKit.gz (per sample, CHH + CHG)
    ↓ [ca_filter.py — already exists, run on HPC]
CA-only .methylKit.gz (per sample)
    ↓ [bsseq::read.bismark() or methylKit::methRead()]
BSseq / methylKit object
    ↓ fork:
    ├── [dmrseq] → genome-wide DMR discovery (unbiased)
    ├── [regionCounts → edgeR GLM] → gene-body differential (aggregated)
    ├── [regionCounts → edgeR GLM] → enhancer/TAD/peak-region differential
    └── [BigWig → deepTools] → visualization at features
```

The CA-context filtering (`ca_filter.py`) is the critical first step — everything downstream should work on CA-only data, not pooled CHH/CHG. You already have this.

## What to build vs what exists

**Already built:** CA-context filtering, gene-body aggregation, BigWig conversion pipeline, genome-wide rate comparison, atlas comparison, all ChIP/ATAC/K119ub BED files loaded in `_shared_config.R`

**Need to build:**
1. Per-sample gene-body CA count extraction (modify `compute_emseq_per_gene()` to output per-sample counts, not pooled — ~30 lines)
2. R script wrapping edgeR/methylKit for differential testing on aggregated counts (~150 lines)
3. Region-based analysis script using GenomicRanges to intersect mCA with K119ub/H3K27ac/TAD boundaries (~200 lines)
4. DUET CG × EM-seq mCA comparison script (~100 lines)

The dmrseq genome-wide scan is the most useful "just run it and see what comes out" tool and requires minimal custom code — it's designed to find DMRs without pre-specifying regions.