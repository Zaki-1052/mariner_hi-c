## What the Anchor vs. Span Results Would Tell Us

The ChromHMM analysis (Phase 4.4 — the Popay Fig 2f equivalent) asks: **where is Polycomb (H3K27me3) enriched relative to differential loops?** The two possible outcomes map to two distinct mechanistic models for how BAP1 loss disrupts 3D genome organization:

### Result 1: Polycomb enriched at **anchors** → Sensitivity Model

If H3K27me3 is concentrated specifically at the loop anchor points (where CTCF binds), this means BAP1 loss → H2AK119ub elevation → H3K27me3 spreading → **direct disruption of CTCF binding sites**. The repressive chromatin physically occludes CTCF, preventing it from anchoring loops.

This connects to Jesse Dixon's reference to the **Brad Bernstein IDH mutant glioma work** (dixon-meeting-summary L39): in IDH mutants, repressive chromatin gains knock out CTCF binding sites, collapsing loops. Same mechanism, different upstream cause.

Under this model:
- Lost loops break because their CTCF anchors become heterochromatinized
- Gained short CTCF-CTCF loops form at sites that escape repression
- The length dependence (lost loops longer, gained shorter) would be secondary — longer loops just have more anchors exposed to spreading Polycomb
- The CTEA meeting's question about centering on CTCF motifs (ctea-april-meeting-notes L4-5, L13-14) becomes critical — you'd want to show the chromatin state change is happening right at the motif

### Result 2: Polycomb enriched across the **span/body** → Extrusion Impediment Model

If H3K27me3 is spread across the loop body (the chromatin between anchors) but *not* concentrated at the anchors themselves, this means Polycomb domains act as **physical roadblocks to cohesin loop extrusion**. Cohesin starts extruding from a loading site, travels along chromatin, but gets stuck when it hits a Polycomb/H3K27me3-dense region.

This is exactly what Popay found for the "mixed-dependency" cluster in her NIPBL depletion data (CONTEXT-CLUSTER L128-129): anchors devoid of K27me3, but **strong K27me3 signal flanking the anchor**. The CTCF sites themselves are fine — it's the territory between them that's changed.

Under this model:
- Lost loops are long because cohesin can't extrude far enough to reach the distal anchor — it gets blocked by newly-established Polycomb domains mid-span
- Gained short loops form because cohesin stalls early, creating new contacts at closer range
- This directly explains the CTEA observation (ctea-april-meeting-notes L41-49): "cohesin trying to extrude/grow longer but more polycomb regions blocking it... extending along euchromatin stretch/span, but then hits polycomb K27me3 and stops"
- The length-dependent enrichment pattern (longer lost loops = more K27me3-associated) makes mechanistic sense: longer spans have more opportunity to accumulate Polycomb roadblocks
- K27ac-anchored loops showing NO length difference (p=0.41, CONTEXT-CLUSTER L33) fits perfectly — CRE loops are short and enhancer-driven, not extrusion-dependent

### Why This Distinction Matters for the Paper

The extrusion impediment model is the stronger story because it:
1. **Unifies the length and directionality findings** — long loops lost, short gained, with a physical mechanism
2. **Connects H2AK119ub mechanistically** — BAP1 loss → ub elevation → Polycomb spreading along chromatin → cohesin stalling. Your K119ub data (10x enrichment at disrupted anchors) provides the upstream trigger that Popay's NIPBL depletion study couldn't address
3. **Explains the progressive developmental phenotype** — Polycomb domains spread gradually over development, progressively impeding more extrusion events (200 differential loops early → 3,000 in adult)
4. **Matches the stripe data** — gained stripes enriched for H3K27me3 anchors, developmental TFs (Hox, Shh) — these are exactly the loci where Polycomb spreading would create new extrusion barriers

The sensitivity model is simpler but less interesting — it's essentially "repression kills CTCF binding" without explaining the length dependence or progressive nature.

### Most Likely Outcome

Based on the existing data, **both signals will probably be present** but the span enrichment should dominate in the clusters most sensitive to BAP1 loss. Jesse noted that Popay's mixed-dependency cluster showed exactly this hybrid pattern: anchors clear of K27me3, flanking regions loaded with it. Your data has the added advantage of H2AK119ub ChIP, which should show even clearer span enrichment since ub is the upstream signal that recruits PRC2 → H3K27me3.

---

## Stripes and the Extrusion Impediment Model

Stripes are the direct visual signature of **active loop extrusion in progress**. A stripe appears in Hi-C when cohesin loads at one site (usually a CTCF site, promoter, or super-enhancer) and extrudes outward, creating a line of contacts emanating from that anchor. The stripe extends until cohesin falls off or **hits a barrier**.

### The Connection

Under the extrusion impediment model, Polycomb spreading creates new barriers that **shorten or terminate stripes prematurely**:

1. **In wildtype**: Cohesin loads at a CTCF site, extrudes long distances, forming long loops and extended stripes
2. **In BAP1-KO**: H2AK119ub elevation → Polycomb (H3K27me3) spreads into euchromatic regions → cohesin stalls when it hits these new Polycomb domains → stripes get **truncated**, and new short-range contacts (gained loops) form at the stall points

This is exactly what the CTEA meeting was getting at (ctea-april-meeting-notes L41-49): "cohesin trying to extrude/grow longer but more polycomb regions blocking it... extending along euchromatin stretch/span, but then hits polycomb K27me3 and stops." The question "why loops vs stripes?" (L42) is answered: they're the same phenomenon at different scales — a stalled stripe creates a loop at the stall point.

### Your Stripenn Data Confirms This

From CONTEXT-CLUSTER Section 10, your completed Stripenn analysis shows:

- **More gained than lost stripes in adult** (2,052 gained vs 1,528 lost) — new short-range extrusion events replacing long-range ones
- **H3K27me3+ anchors 1.75x enriched in gained vs lost stripes** (9.1% vs 5.2%) — gained stripes emanate from Polycomb-associated sites
- **Active_Enhancer enriched in lost stripes** (19.5% vs 14.0%) — long-range enhancer-driven extrusion is preferentially disrupted
- **Gained stripes hit developmental TFs** (Hox, Shh, Dlx5/6) — exactly the loci where Polycomb domains are known to be large and expanding in BAP1-KO
- **Directional bias reverses between timepoints** (more lost at P12, more gained in adult) — progressive Polycomb spreading creates more barriers over time

### Stripes as Super-Enhancer Hubs

The CTEA meeting also noted that stripes can represent **super-enhancer hubs** — one SE contacting multiple promoters via extrusion (ctea-april-meeting-notes L29-36). If Polycomb spreading truncates these stripes:

- The SE loses contact with distal promoter targets
- Genes at the far end of the stripe lose enhancer input → downregulation
- Genes proximal to the SE retain contact → potentially upregulated relative to distal targets
- This would explain why lost stripes overlap with developmental/synaptic DEGs

### How Clustering Connects Stripes to the Anchor-vs-Span Question

The Popay clustering analysis will group loops by their sensitivity to BAP1 loss. If the extrusion impediment model is correct, you'd expect:

- **Most-sensitive clusters** (biggest loss in mutant): long loops with Polycomb enrichment **across the span** — these are the loops whose extrusion path got blocked
- **Stripes at those same loci**: truncated in mutant, with new gained stripes/loops appearing at shorter range
- **Least-sensitive clusters**: short CRE-type loops where extrusion distance is minimal, so Polycomb barriers don't matter

The stripe data provides independent validation — if the clustering shows span-enriched Polycomb in the most-sensitive loops AND the Stripenn data shows truncated stripes at those same anchors, that's two orthogonal lines of evidence for the same mechanism.

In short: **loops are the endpoints of extrusion, stripes are extrusion caught in action, and Polycomb barriers explain why both change in BAP1-KO**. The anchor-vs-span analysis tells you whether the barrier is at the endpoint (sensitivity) or along the path (impediment), and your stripe data already points toward impediment.

---

Great question — I think I was unclear. Let me reframe.

The key is **which loops you're looking at**. The ChromHMM analysis measures enrichment at the anchors and spans of the **original wildtype loops** (all 39,344), clustered by how much they change in BAP1-KO.

For a **lost** long loop (the main signal):
- The two anchors are CTCF sites — they're still there, still open chromatin
- Cohesin used to extrude from anchor A all the way to anchor B
- Now, Polycomb has spread into the chromatin **between** A and B (the span)
- Cohesin starts extruding, hits the new Polycomb domain **mid-journey**, and stalls — it never reaches anchor B
- Loop is lost because the extrusion path is blocked

So Polycomb is in the **span** of the lost loop — not at its anchors. The anchors are still CTCF. The roadblock is in the middle.

Your intuition — "cohesin hits Polycomb and stops, so Polycomb should be at the anchor" — is actually correct for **gained** loops. The new short-range contact forms exactly where cohesin stalled, at the **edge** of the Polycomb domain. That stall point becomes the anchor of the new gained loop. This is what Popay saw in her mixed-dependency cluster: "anchors devoid of K27me3 but strong K27me3 signal **flanking** the anchor" — the Polycomb is right next to where cohesin parked.

So both things are true simultaneously:
- **Lost loop spans** = Polycomb-invaded (the roadblock that killed the loop)
- **Gained loop anchors** = Polycomb-flanked (where cohesin stalled and made a new contact)

Same mechanism, different perspective depending on which set of loops you're measuring.