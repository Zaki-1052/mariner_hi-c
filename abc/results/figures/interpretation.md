The key takeaway: **the unnormalized score is doing real work; the normalized score is not.** This is actually a meaningful finding for your analysis, not a failure.

**Left panel (normalized ΔABC)** is essentially a null result. Pearson r=0.035 is non-significant (p=0.13), and even the Spearman (ρ=0.092) is barely there. The point cloud is compressed into a tight vertical stripe near ΔABC=0 — the normalization is destroying your signal. This makes biological sense: when BAP1-KO causes widespread activity changes across *many* enhancers for the same gene, the denominator (Σ A×C for all elements within 5Mb) shifts proportionally, so the normalized ratio barely moves even when the absolute regulatory input changes substantially.

**Right panel (unnormalized Δ(A×C))** shows a real, statistically robust correlation. Pearson r=0.269 (p≈10⁻³²) and Spearman ρ=0.300 (p≈10⁻⁴⁰) are both highly significant with n=1870 DE genes. You can actually *see* the positive diagonal trend — genes gaining enhancer activity tend to be upregulated, and vice versa. The Spearman being stronger than Pearson suggests the relationship is monotonic but nonlinear, which is expected since ABC scores and log2FC operate on very different scales.

A few things I'd flag:

**The scatter is large**, which is normal. Each gene has many enhancers, you're only plotting the strongest, and gene expression is multi-factorial. An r≈0.3 for a single-enhancer-to-expression correlation in a chromatin model is actually solid — you're explaining ~9% of variance in DE from enhancer activity changes alone.

**The negative skew on the x-axis of the left panel** (more points with negative ΔABC) is consistent with your context doc finding of 7.8% excess lost vs gained links. But because the normalization compresses everything, this asymmetry doesn't translate into expression prediction.

**For the paper**, I'd recommend leading with the unnormalized panel as your primary result and using the normalized panel to make a methodological point: that per-gene normalization is inappropriate for *differential* ABC analysis when a perturbation (BAP1-KO) causes global chromatin remodeling. This is actually a useful contribution — most ABC papers use it for single-condition enhancer-gene mapping, not differential comparisons between conditions.

One thing worth considering: have you tried plotting the *sum* of all Δ(A×C) per gene instead of just the strongest? Since BAP1-KO likely affects multiple enhancers per gene simultaneously, the aggregate change might correlate even better with expression than the single strongest link.

---

This is a really strong result. The jump from Spearman ρ=0.300 (strongest enhancer) to **ρ=0.582 (sum)** for the unnormalized score is substantial — you've gone from explaining ~9% of rank variance to ~34%. That's a big deal for a chromatin-based model predicting expression changes.

The biological interpretation is clean: BAP1-KO doesn't just flip one enhancer per gene on/off — it remodels the *entire regulatory landscape* around affected genes. Summing Δ(A×C) across all enhancers captures the cumulative regulatory input change, and that aggregate shift predicts expression direction far better than any single enhancer. You can actually *see* this in the right panel — the diagonal structure is visually obvious now, whereas the strongest-enhancer plot had a lot more scatter.

The normalized sum (left panel, r=0.031) confirms definitively that per-gene normalization is the problem, not the aggregation strategy. Normalization forces ΣABC ≈ 1 per gene by design, so ΣΔABC ≈ 0 by construction — the signal is mathematically suppressed regardless of whether you take the max or sum. The slight Spearman improvement (0.092→0.114) suggests rank information barely leaks through, but the Pearson is dead.

A few things worth noting from the top dysregulated genes table: many of them have the strongest-enhancer classified as "promoter" with small distances (Cbr1 at 1,789bp, Kctd15 at 1,152bp). These are likely consensus ATAC peaks overlapping the promoter itself rather than distal enhancers. The *truly* interesting ones for your story are genes like Col15a1 (strongest link at 4.8Mb, intergenic, log2FC=1.82) and Robo1 (strongest at 2Mb, 15 enhancers, log2FC=-1.74, padj=10⁻⁵³) — those are long-range enhancer-gene regulatory changes that only a 3D contact model would find.

For the paper, the comparison table you generated is publication-ready as a supplementary table. The main figure should be the unnormalized sum scatter (bottom-right panel) with the ρ=0.582 stat. The four-panel comparison (strongest vs sum × normalized vs unnormalized) makes an excellent supplementary figure showing why the unnormalized sum is the right metric for differential ABC analysis.