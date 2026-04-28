The run completed successfully — all 6 oriented BED6 files built correctly, `computeMatrix` finished, and the figure rendered with the right x-axis labels (`-5kb (exterior)` / `anchor` / `+5kb (interior)`). The "Skipping ._r..." messages are the plotting parser dropping blacklisted-region rows, same as Phase 5 — harmless.

## What the oriented figure shows

**H3K27me3 (the key columns):**
- **clust5** (cyan, 97% gained): K27me3 is high and **symmetric on both sides** of the anchor — the signal is broad and comparable on exterior and interior flanks. This cleanly matches your mentor's "Polycomb loop" prediction: the anchor is embedded within a Polycomb domain, not at a boundary.
- **clust4** (olive, 70% gained): Similar bilateral K27me3 pattern but lower magnitude. Same model, weaker effect.
- **clust1/2/3** (unchanged loops): Low, flat K27me3 baseline on both sides — no Polycomb involvement, as expected.
- **clust6** (brown, 78% lost): Moderate K27me3 signal, but the line plots are overlaid across all 6 clusters, making it hard to visually resolve whether clust6 specifically shows exterior > interior asymmetry at this scale.

**H3K27ac (context):**
- clust1/3/6 have the strongest K27ac peaks. Notably, the K27ac signal for clust6 appears to show a slight **interior bias** (right side of the plot) — consistent with the loop body being euchromatic while the exterior flank has repressive chromatin. This is a good sign.

**H2AK119ub / H3K27me1:**
- clust5 elevated in both, symmetric — Polycomb domain signature
- clust6 shows a visible ctrl→mut increase in K119ub, confirming these anchors gain ubiquitination upon BAP1 loss

## The honest limitation

The all-clusters-overlaid line plot makes it difficult to visually quantify clust6's asymmetry. The signal difference between exterior and interior for clust6 is subtle compared to the dramatic cluster-to-cluster differences. To rigorously answer your mentor's question, we should extract the mean signal in the exterior half (-5kb to 0) vs interior half (0 to +5kb) per cluster and compute an asymmetry ratio — that would give a single number rather than relying on visual inspection of overlapping lines.

Want me to write a quick quantification script that reads the `oriented_anchors_values` matrix and computes exterior vs interior mean signal per cluster per mark?

---

These results are striking — and they tell a different story than the simplest version of the mentor's hypothesis. Let me lay it out clearly.

## Results

### clust6 (lost loops, 78% down) — NOT asymmetric for K27me3

| Mark        | Exterior | Interior | Ext/Int | p-value |
| ----------- | -------- | -------- | ------- | ------- |
| K27me3_ctrl | 1.6966   | 1.6629   | 1.02    | 0.52 ns |
| K27me3_mut  | 2.4952   | 2.4397   | 1.02    | 0.57 ns |
| K27ac_ctrl  | 7.5889   | 7.4114   | 1.02    | 0.48 ns |

K27me3 at clust6 anchors is **symmetric** — no significant difference between exterior and interior (p > 0.5). The signal is equally distributed on both sides. This means clust6 anchors are NOT behaving as clean euchromatin-to-heterochromatin boundaries. The Polycomb enrichment seen in ChromHMM (2.09× at anchors, 0.94× at span) comes from the anchor region itself — not from an asymmetric flanking domain.

### clust5 (gained loops, 97% up) — INTERIOR-biased K27me3

| Mark        | Exterior | Interior | Ext/Int | p-value  |
| ----------- | -------- | -------- | ------- | -------- |
| K27me3_ctrl | 4.4652   | 4.8845   | 0.91    | 0.004 ** |
| K27me3_mut  | 3.8822   | 4.2705   | 0.91    | 0.004 ** |

This is the surprise: clust5 K27me3 is significantly **interior-biased** (Ext/Int = 0.91, p < 0.005). The Polycomb signal is ~9% stronger on the loop-body side than the exterior. This isn't symmetric embedding in a Polycomb domain — the domain extends preferentially *inward* into the loop. H2AK119ub shows the same interior bias (0.95–0.97 Ext/Int).

### clust4 (moderate gain, 70% up) — same pattern as clust5

K27me3 Ext/Int = 0.93 (p < 0.001). Interior-enriched. Consistent with Polycomb domain expansion preferentially into the loop body.

## What this means for the mentor's framework

1. **The extrusion impediment model (one side euchromatic, other heterochromatic) is NOT supported for clust6.** The K27me3 enrichment at lost-loop anchors is symmetric and modest — these anchors aren't sitting at sharp Polycomb domain boundaries. The loss of these loops may be more about direct anchor disruption (the K27ac drops from 7.59 to 4.77 ctrl→mut) than Polycomb-mediated extrusion blocking.

2. **The Polycomb domain compaction model for clust5 IS supported, with a refinement**: Polycomb is not uniformly distributed around the anchor — it's biased toward the loop interior. This is consistent with Polycomb-mediated compaction pulling the loop body together: the K27me3 domain extends inward from both anchors into the loop span, forming a self-reinforcing compacted domain.

3. **The two mechanisms are real but operate differently than the simple boundary model.** Lost loops (clust6) lose active chromatin marks symmetrically at anchors. Gained loops (clust5) gain Polycomb asymmetrically *into* the loop body.

The full quantification table is saved at `cluster/bap1_late/figures/deeptools/oriented_anchors/asymmetry_quantification.tsv`.


---

The figure came out clean. The profiles tell the story visually and the bar chart quantifies it.

**What the figure shows:**

The top row reveals a **monotonic gradient** from left to right. Lost/unchanged clusters (clust6, clust3, clust1) have flat, symmetric K27me3 profiles — the signal is the same on both sides of the anchor. As you move toward gained clusters (clust2 → clust4 → clust5), K27me3 progressively shifts toward the **interior** (loop-body side), with clust5 showing the strongest interior bias (E/I = 0.91, p < 0.005).

The bottom row makes the pattern unambiguous: bars at or above the dashed line = symmetric/exterior, bars below = interior-enriched. The significance stars land exclusively on the gained-loop clusters.

**How to represent this scientifically:**

This is actually a more nuanced finding than a binary "boundary vs. domain" result. The gradient suggests that K27me3 doesn't just passively coat loop anchors — it **extends preferentially into the loop body** at gained loops, consistent with Polycomb-mediated compaction drawing the interior together. At lost loops, there's no such directional bias.

For the mentor's framework specifically: clust6 anchors are **not** euchromatin-to-heterochromatin boundaries (no asymmetry). The Polycomb signal at those anchors is modest and symmetric — the loop loss is more likely driven by loss of active marks (K27ac drops from 7.59 to 4.77 ctrl→mut) than by Polycomb acting as an extrusion block.

Outputs saved at `cluster/bap1_late/figures/deeptools/oriented_anchors/orientation_asymmetry_figure.{png,pdf,svg,jpg}` alongside the quantification table (`asymmetry_quantification.tsv`).