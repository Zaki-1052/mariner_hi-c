

**Language > Specialist**: R > Computational Biologist / Epigenomics Integration Specialist

Both scripts ran successfully with no errors. Here's my analysis of the results:

## Section 20: Coordinated RNA-seq Integration

**Strong, clean results.** The key findings are biologically coherent:

- **5,708 coordinated mC up/hmC dn genes**, 5,203 (91.2%) matched to RNA-seq
- **Clear repression bias** in coordinated genes:
  - Tier 1: 742 Down vs 143 Up (5.2:1 ratio, Fisher's OR = 0.11, p = 2.1e-70)
  - Tier 2: 628 Down vs 28 Up (22.4:1 ratio, Fisher's OR = 0.02, p = 3.3e-106)
  - 83-87% show no significant expression change (expected for gene body methylation)
- **Spearman rho = -0.231** (p = 1e-63) -- combined methylation effect negatively correlates with expression, meaning larger mC up/hmC dn changes predict stronger repression
- **2x2 heatmap is striking**: mC Up + Expr Down = 847, mC Down + Expr Up = 278, Fisher's OR = 0.03 (p = 5.9e-118) -- extremely strong anti-correlation

**Interpretation**: mC up/hmC dn genes are overwhelmingly biased toward repression, consistent with TET-mediated demethylation block causing transcriptional silencing.

## Section 21: Discordant Gene Characterization

**Every single dimension is significant** -- discordant genes are genuinely distinct:

| Dimension | Key Finding | p-value |
|-----------|-------------|---------|
| \|mC diff\| | Disc smaller (0.017 vs 0.032) | <2.2e-16 |
| \|hmC diff\| | Disc smaller (0.012 vs 0.028) | <2.2e-16 |
| RNA-seq log2FC | Disc median +0.19 vs Coord -0.01 | <2.2e-16 |
| Expression Up/Down | **OR = 103.14** -- disc strongly upregulated | <2.2e-16 |
| Net ATAC | Disc median +1 vs Coord 0 -- more accessible | <2.2e-16 |
| ATAC Up/Down | OR = 14.55 -- disc enriched for ATAC Up | <2.2e-16 |
| K119ub gene body | Disc -0.031 vs Coord +0.033 -- opposite direction | 2.0e-06 |
| H3K27ac | Disc enriched for gained (101G/17L vs 1177G/823L) | 1.4e-09 |
| Chromatin state | Different distribution | 1.0e-04 |
| MeCP2 | Disc more MeCP2 Down relative to Up | 1.3e-13 |
| Loop involvement | Disc 15.5% vs Coord 5.1% -- 3x more at loops | <2.2e-16 |

**Biological picture of discordant genes**: These are genes where mC goes DOWN and hmC goes UP -- the opposite of the dominant BAP1-KO pattern. They have:
- **Smaller** methylation effect sizes (weaker changes)
- **Upregulated** expression (median log2FC = +0.19)
- **Gained** chromatin accessibility (ATAC Up enriched 14.5x)
- **Lost** K119ub signal (negative gb_log2fc) -- opposite of coordinated genes
- **Gained** H3K27ac (6:1 gained:lost ratio)
- **3x more** likely to be at Hi-C loop anchors

This profile suggests discordant genes are at loci where BAP1-KO *removes* repression rather than adding it -- possibly via indirect regulatory rewiring or compensatory mechanisms at loop-anchored, active regulatory regions.

## Minor warnings (non-issues):
1. **Many-to-many join warning** -- some duplicate gene names in DMR files (expected for overlapping annotations). Could silence with `relationship = "many-to-many"` if desired but doesn't affect results.
2. **No GO enrichment** -- likely because discordant genes (n=727 Entrez) tested against coordinated background (n=5,495) are too functionally similar for enrichment. This is a legitimate null result.
3. **Chromatin state**: As expected, gene body classification yields mostly Active_Promoter (3,950), Repressed_Promoter (618), and Other (2,114) -- the distal categories (Polycomb, Enhancer) are 0 since gene bodies contain their own TSS.

---

For a PI meeting, I'd recommend **3 plots** (or 2 if you want to be concise):

## Must-show

**20d (mC × expression heatmap)** -- This is your money plot. The 2×2 with OR = 0.03 (p = 5.9e-118) is immediately legible: 847 genes in mC Up / Expr Down vs 155 in mC Up / Expr Up. No ambiguity. It answers the question "does methylation gain cause repression?" with a single glance. Publication-ready.

**21a (discordant composite)** -- Normally 9 panels is too much, but here the point IS that every single dimension is significant. Your PI sees one figure and immediately understands: discordant genes aren't noise, they're a biologically distinct class with opposite K119ub (-0.03 vs +0.03), opposite ATAC (net +1 vs 0), opposite expression (+0.19 vs -0.01), and 3x loop enrichment. The consistency across all 9 panels is the argument.

## Strong supporting plot

**20a (stacked bar)** -- The Tier 2 panel is dramatic: 12.1% Down vs 0.5% Up among coordinated genes (22:1 ratio), while "Other mC DMR genes" and "All Genes" show roughly balanced Up/Down. This directly demonstrates the coordinated mC up/hmC dn pattern specifically predicts repression, not just methylation changes in general.

## If your PI wants only 1 plot per analysis

- Section 20: **20d** (heatmap) -- cleanest single-panel summary
- Section 21: **21a** (composite) -- the comprehensive evidence is the point

The scatter plots (20b, 21b, 21c) are good for the paper supplement but less impactful for a meeting since they require more explanation. The violin (20c) is clean but redundant with 20d.