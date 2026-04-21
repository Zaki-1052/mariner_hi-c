Good — the annotated TSV has all the columns (h3k27ac, h3k27me3, h3k4me1, h3k4me3, bivalent, anchor_type) that the docs reference. Both documents are verified.

---

Here's your briefing:

## Key Findings

### Late timepoint (250402/Adult) — the main story

**The headline: BAP1-KO shows a net *gain* of chromatin stripes in adult cerebellum, with H3K27me3-marked (Polycomb) anchors preferentially gaining stripes.**

- 2,320/7,371 stripes (31.5%) are significant at FDR<0.05 — driven by very tight replicates (BCV=0.012), not by large fold-changes
- **More gained than lost** (2,052 vs 1,528), opposite to what you see at P12
- H3K27me3+ anchors are **1.75x enriched in gained vs lost** stripes (9.1% vs 5.2%) — this is the clearest directional asymmetry and is consistent with BAP1's Polycomb biology
- Active_Enhancer anchors go the other way: enriched in lost stripes (19.5% vs 14.0%)
- GO enrichment is **entirely one-directional** — only gained stripes show enrichment: developmental TFs (Hox clusters, Shh, Sox11, Dlx5/6, Tgfb2), synaptic genes (46 genes at postsynaptic specialization), and ion channels
- No GO enrichment for lost stripes at all; no KEGG pathways significant for either direction

### Early timepoint (250831/P12) — weak signal

- Only 2.4% significant — essentially underpowered
- Directional bias is **reversed**: more lost than gained (949 vs 776)
- Only 22 high-confidence stripes total; enrichment uses just 26 genes (lost) and 19 genes (gained)
- Lost-stripe enrichment hits Wnt5a/Fgfr2/Klf4 stem cell pluripotency pathway (KEGG p.adj=0.009) — interesting but based on 4 genes
- Bivalent_Promoter anchors are 3x more prevalent at P12 than adult (3.6% vs 1.2%), consistent with developmental bivalency resolving over time

## Critical Caveats to Flag

1. **All effect sizes are tiny** — max |logFC| = 0.39 across 7,371 stripes. Zero stripes reach "moderate" (0.5) fold-change. This is frequency-driven signal, not magnitude-driven.
2. **BCV of 0.012 inflates significance** — the 31.5% rate partly reflects the edgeR model being very sensitive when replicates are this tight. The "real" biological effect is subtle.
3. **1,249 both_discordant stripes** at the late timepoint — significant at both resolutions but with *opposite* direction. This means ~half of cross-resolution validated stripes disagree on gained vs lost.
4. **60.7% directional consistency** — ~39% of source-based direction calls (control-only → lost) have edgeR logFC pointing the wrong way.

## Plots to Look At

| Plot                         | Location                                                          | What it shows                                           |
| ---------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------- |
| **Volcano (combined)**       | `outputs/combined/volcano_combined/`                              | Late has the cloud of significant points; early is flat |
| **Anchor annotation (late)** | `outputs/250402/visualizations/anchor_annotation_250402/`         | The Polycomb/H3K27me3 enrichment in gained stripes      |
| **GO BP dotplot (late)**     | `outputs/combined/enrichment/go_bp_dotplot_late/`                 | Hox/developmental terms all in gained direction         |
| **GO CC dotplot (late)**     | `outputs/combined/enrichment/go_cc_dotplot_late/`                 | Synaptic compartment enrichment in gained               |
| **Replicate correlation**    | `outputs/250402/visualizations/replicate_correlation_250402_5kb/` | Shows how tight the replicates are (drives low BCV)     |
| **Cross-res comparison**     | `outputs/combined/cross_res_comparison/`                          | logFC correlation between 5kb and 10kb                  |
| **Direction breakdown**      | `outputs/combined/source_direction_summary/`                      | Lost vs gained across both timepoints                   |

## What to Tell Your PI

**Framing:** "Stripenn finds widespread but subtle differential stripe signal in adult BAP1-KO cerebellum. The data supports a model where BAP1 loss leads to a net gain of chromatin stripes, with gained stripes enriched at Polycomb-marked (H3K27me3+) and poised enhancer anchors — consistent with derepression at Polycomb target loci. GO enrichment in gained stripes hits developmental TFs (Hox, Dlx, Shh) and cerebellar synaptic genes."

**Honest caveats to include:** "All fold-changes are below 0.4 — this is a frequency signal (many stripes with tiny changes), not individual stripes dramatically appearing or disappearing. The very low BCV (0.012) means edgeR is declaring statistical significance at effect sizes that are at the edge of biological interpretability. JuiceBox validation of the top-ranked loci is the critical next step before drawing biological conclusions."

**The developmental angle:** "The directional bias reverses between P12 (more lost) and adult (more gained), and bivalent promoter anchors are 3x more prevalent at P12 — both consistent with progressive Polycomb-mediated chromatin remodeling downstream of BAP1 loss."