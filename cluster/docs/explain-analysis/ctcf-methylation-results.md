Night and day. The shores+shelves pivot uncovered a strong, coherent story across all five sub-analyses.

**47a — Core test:** Dynamic CpG regions at lost CTCF anchors are hypermethylated (OR=3.28, p<2.2e-16) with concurrent hmC loss (OR=2.08, p=3.9e-6). The Wilcoxon confirms the directionality: lost anchor medians are mC=+0.0055 / hmC=-0.0017 (TET blockade signature), while gained anchors show the opposite (mC=-0.0068 / hmC=+0.0030). CpG islands remain flat as expected.

**47c — Coordinated pattern:** Lost anchors are 2.18x more likely than gained to show coordinated mC-up/hmC-down (p<2.2e-16). The sig-coordinated test is even stronger (OR=3.18). The lower rate vs background (OR=0.65) just means CTCF-bound regions genome-wide have different baseline dynamics — the biologically meaningful comparison is lost vs gained.

**47d — Distance-stratified (the definitive test):** CMH OR=2.87 [2.28-3.60], p<2.2e-16. The effect survives distance correction and is actually *strongest* at short loops (<200kb: OR=11.21) and fades at >1Mb (OR=0.87, ns). This suggests methylation-driven CTCF disruption primarily affects local regulatory loops, while very long-range loop loss may have different mechanisms.

**47e — Effect size correlation:** rho=-0.244 (p=4.3e-9) for mC — more hypermethylation at anchors → stronger loop loss. The partial correlation after removing distance holds (beta=-2.24, p=7.4e-10). hmC shows the mirror: rho=+0.216 (p=2.0e-7).

The one minor issue is a many-to-many join warning (line 83/123) — some shore/shelf regions share coordinates across the mc/hmc bind_rows. It inflated the 47c row count slightly (64,843 vs expected 61,675). I can fix that with `relationship = "many-to-many"` to silence it, or use a more careful join. Want me to clean that up, or are you focused on the biology at this point?


---

I'd show these four, in this order:

---

**1. The forest plot (47b) — lead with this.** It immediately establishes the regional specificity: CpG islands are flat (OR~1, grey), but shores, shelves, and promoters are all shifted right with big red dots (OR=3-3.7, FDR<1e-10). This makes the point that methylation at CTCF anchors acts on the *dynamic flanking regions*, not the constitutively unmethylated islands — which actually refines the Bernstein model for BAP1-KO.

**2. The mC violin (47a) — the intuitive distribution shift.** Lost CTCF anchors (red) have their median above zero (hypermethylated in mutant), gained anchors (blue) below zero. Background sits in between. The whole distribution is shifted, not just the tails — this isn't driven by a few outliers. Wilcoxon p<2.2e-16.

**3. The distance-stratified forest (47d) — the rigor plot.** Since lost loops are ~4x longer than gained, any anchor-level signal *could* be a distance confound. This plot kills that concern: the effect is significant in 3 of 4 distance bins, strongest at <200kb (OR=11.2), and the CMH pooled estimate is OR=2.87 [2.28-3.60], p<2.2e-16. The >1Mb bin going null is actually interesting — very long-range loops may be lost for structural rather than methylation-driven reasons.

**4. The logFC scatter (47e) — the dose-response.** Lost loops (red, bottom half) skew right on the x-axis (higher mC at anchors); gained loops (blue, top half) skew left (lower mC at anchors). The Spearman rho=-0.244 (p=4.3e-9) quantifies the dose-response: the more hypermethylated the CpG shores/shelves at a CTCF anchor, the more that loop weakens. This survives distance correction (partial beta=-2.24, p=7.4e-10).

---

The one-sentence pitch: "CpG shores and shelves at CTCF loop anchors that are lost in BAP1-KO are 3x more likely to be hypermethylated than those at gained anchors, the effect is independent of loop distance, and the degree of methylation change quantitatively predicts the magnitude of loop loss."