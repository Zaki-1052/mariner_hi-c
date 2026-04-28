Now I have all the data. Here's how your mentor's framework maps onto the clusters:

---

Your mentor is laying out a clean diagnostic: **where is the H3K27me3 relative to the loop boundary?** That distinction separates two mechanisms. Let me walk through the data.

## The two clusters that matter most

### Clust5 (667 loops, 97% gained) — **Polycomb loop model**

|            | Polycomb enrichment | Bivalent enrichment |
| ---------- | ------------------- | ------------------- |
| **Anchor** | **6.59×**           | **7.91×**           |
| **Span**   | **3.03×**           | **2.21×**           |

Polycomb is enriched at **both** anchors AND span. The span is NOT euchromatic — it's 3× enriched for Polycomb. The anchor proportion data confirms this: **87% of clust5 anchors have Polycomb as their dominant ChromHMM state**. Active Enhancer is at 0.54× (depleted below genome average), Active Promoter at 1.84×.

By your mentor's logic: span is Polycomb-enriched → **these are loops within Polycomb domains**. The anchor is acting as an anchor for a Polycomb loop, not as an extrusion stopping site. These gained loops are consistent with **Polycomb-domain compaction** — BAP1 loss elevates H2AK119ub, Polycomb spreads, and new short-range contacts form within the expanding Polycomb territory.

This is further supported by:
- These loops are short (median 330 kb)
- Only 2.85% CRE loops (depleted for enhancer-promoter contacts)
- 54.7% structural (CTCF-CTCF) — Polycomb spreading into existing CTCF boundaries
- H2AK119ub is elevated in *both* ctrl and mut at these anchors (1.37 ctrl, 1.25 mut) — these are already Polycomb-marked loci

### Clust6 (2,359 loops, 78% lost) — **Anchor disruption / extrusion impediment**

|            | Polycomb enrichment | Bivalent enrichment |
| ---------- | ------------------- | ------------------- |
| **Anchor** | **2.09×**           | **2.48×**           |
| **Span**   | **0.94×**           | **0.95×**           |

Polycomb is enriched at anchors but the **span is at genome baseline** (~1.0×). The span is euchromatic — Active Enhancer at 1.12×, Active Promoter at 1.11×, Strong Enhancer at 1.12× in the span. Meanwhile, 25.5% of clust6 anchors have Polycomb as the dominant state, and 30.9% are CRE loops (the highest CRE fraction of any cluster).

By your mentor's logic: span is euchromatic, anchor has Polycomb → **this is where loop extrusion stops at a Polycomb site**. The Polycomb mark at the anchor acts as a boundary. These are active loops (enhancer-promoter, CRE) that traverse euchromatic territory but whose anchors get invaded by repressive chromatin upon BAP1 loss, disrupting the loop.

Supporting evidence:
- These loops are long (median 550 kb — longest of any cluster, p < 10⁻284 vs clust1/2)
- Highest Active Promoter enrichment at anchors (7.22×) and Active Enhancer (4.42×)
- H2AK119ub signal increases in mutant (ctrl 1.16 → mut 1.58), consistent with the anchor getting ub'd

## The asymmetry test your mentor suggests

Your mentor's "further evidence" — **one side euchromatic, the other heterochromatic** at the boundary — isn't directly testable from the current ChromHMM OverlapEnrichment output, which averages over all anchors per cluster. To test this, you would need to:

1. Take clust6 anchors individually
2. Extend each anchor ±N kb (e.g., ±10–20 kb)
3. Ask whether the ChromHMM state on the loop-interior side vs. loop-exterior side is asymmetric (euchromatic inside, Polycomb outside or vice versa)

The **Phase 5 deepTools metagene** data partially addresses this — the H3K27me3 ±5kb profile centered on clust6 anchors would show whether the signal is symmetric or skewed to one flank. This is the same analysis that Popay did (Fig 2g) where the mixed-dependency cluster showed K27me3 *flanking* the anchor, not at the anchor itself.

## Summary

|                         | Clust5 (gained)                                  | Clust6 (lost)                                         |
| ----------------------- | ------------------------------------------------ | ----------------------------------------------------- |
| Anchor chromatin        | Polycomb (6.59×)                                 | Polycomb (2.09×) + Active (7.22× promoter)            |
| Span chromatin          | Polycomb (3.03×)                                 | Euchromatic (~1.0×)                                   |
| Mentor's interpretation | **Polycomb loop** — anchor for a Polycomb domain | **Extrusion stop** — Polycomb boundary in euchromatin |
| Mechanism               | Polycomb domain compaction                       | Anchor disruption / extrusion impediment              |
| Loop size               | Short (330 kb)                                   | Long (550 kb)                                         |

**Both mechanisms are operating simultaneously**, just in different loop classes. The data cleanly separates them. The next step your mentor is pointing toward — testing chromatin asymmetry flanking the boundary — would nail whether clust6 anchors are truly euchromatin-to-heterochromatin transition points rather than just anchors that happen to overlap a Polycomb peak.