# Poster planning: Phase 1 decisions

This document captures what was decided in the Phase 1 planning chat. It supplements the main project prompt (`poster-PROMPT.md`) — don't repeat what's already there.

---

## The one-sentence story

BAP1 loss accumulates H2AK119ub, which collapses long-range Polycomb loops and disrupts active enhancer contacts through two distinct mechanisms — revealed by clustering the genome's 39,000 loops by their response.

---

## Three result beats (confirmed)

### Beat 1 — Distance shift
Lost loops are systematically longer than gained loops. At 212 shared anchor hubs, the same locus loses a distant partner and gains a nearer one. This is the most accessible finding and sets up the clustering mechanistic result.

### Beat 2 — Clustering reveals two mechanisms
The Popay-style clustering (k=6, 39,344 loops) is new work not in the submitted abstract and should be central to the poster. The key result is anchor-vs-span ChromHMM enrichment:
- **Clust5 (gained, n=667, 97% up):** Polycomb at anchors (6.59×) AND span (3.03×) → Polycomb domain compaction
- **Clust6 (lost, n=2,359, 78% down):** Polycomb at anchors (2.09×), span at baseline (0.94×) → anchor disruption of active/CRE loops

This answers Jesse Dixon's central question (anchor vs span?) and shows both mechanisms operate simultaneously in different loop populations.

### Beat 3 — K119ub as driver + temporal progression + gene expression
- Logistic regression OR = 10.7 for K119ub predicting loop loss
- Progressive: 165 differential loops at P12 → 2,910 in adult (18× amplification)
- ABC model concordance (88%) connects structural changes to gene expression
- Keep ABC presentation minimal — just enough to say "this matters for genes," not the Activity_Lost vs K119ub_Only mechanistic detail

---

## What's cut (and where it lives instead)

| Analysis | Why cut | Where it goes |
|----------|---------|---------------|
| TAD boundaries (merge/split/strength) | Requires too much context setup (what TADs are, what TADCompare does, what merge means) | Verbal Q&A |
| Compartments (44% genome shift) | PC1 visualizations not interpretable for non-specialists at poster scale | Mention as a number in conclusions text |
| Stripes | Real results (2,320 significant in adult, Polycomb-enriched gained stripes), but small effect sizes and the stripe-loop relationship needs too much comparative context | Conclusions text ("consistent with Polycomb-axis effects rather than cohesin/extrusion") or verbal Q&A |
| 9-mark / 18-state ChromHMM | Paper-level resolution; 5-mark model tells the same story more simply | Paper only |
| ABC mechanistic detail | Activity_Lost vs K119ub_Only threshold model is too nuanced for poster | Verbal Q&A for judges |
| DiffBind details, permutation tests | Supporting evidence for paper | Paper only |
| Individual locus examples | Contact maps require too much annotation at poster scale | *Exception:* Zara has a custom bigwig visualization tool for Syt1 — evaluate in next session whether it's clean enough to work as a small supporting panel |

---

## Open items for next session

1. **Evaluate clustering plots** — Zara will upload candidate figures. Assess which ones communicate the two-mechanism result at 36×48 poster scale without requiring ChromHMM expertise from the viewer.
2. **Evaluate the Syt1 bigwig visualization** — determine if it works as a small supporting panel or if it's too context-heavy.
3. **Begin Phase 2: layout and section planning** — spatial arrangement, word budget per section, figure placement, title drafting.

## Important (User) Note

NONE of the uploaded figures will be the "final" plot to show on the poster! Every single one can and will be redone in R to fit with the poster format and be simplified, made more focused, etc. Evaluate based on its conceptual role in the story and what aspect of the story it attempts to show; NOT based on whether it is "poster-ready". It can and will be made to fit, we are choosing WHAT to represent, and I can always fix HOW to represent them.
