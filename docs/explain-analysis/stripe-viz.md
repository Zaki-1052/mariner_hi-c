Please follow the @docs/stripenn-analysis-prompt.md to update both the late (250402) and early (250831) timepoint results documents - launch Sonnet 4.6 subagents to read the prompt in full, explore the results in outputs/ (make separate ones for each timepoint), and concurrently read said outputs from the final viz script to write/update the document describing them. Please make sure your prompt is specific and **requires** the full stripenn-analysis-prompt to be read and followed, and that their findings are appropriately supported by the data, not overly speculative, and based in the actual results. Lastly, briefly check their work after you are done, and then let me know in chat the key important findings and takeaways I need to know immediately, (overall, but mainly focused on the late timepoint (doc) since that has more stripes) which plots to look at that show this, and what to tell/show my PI upon update. If the statistical power or "realness" of the stripes was too low (ie through p-val or logfc) let me know/flag that and warn me about any other potential caveats you see. That should be sufficient to get me up to speed, thanks!

---

Good questions. Let me assess both based on what the script produces vs what already exists.

## Most Insightful Outputs from Stage 7

**Tier 1 — Biologically novel (didn't exist before):**

1. **Anchor type breakdown by direction** (Section 8, `anchor_annotation_*.{pdf,svg}`) — This is the single most important new result. BAP1 deubiquitinates H2AK119ub1, so the key biological question is whether lost stripes concentrate at Polycomb/Bivalent anchors and gained stripes at active/derepressed loci. The `{tp}_annotated_stripes.tsv` files will have this per-stripe.

2. **ChIP-seq overlap rates: lost vs gained** (Section 8, the `p_chip` panel) — H3K27me3 enrichment in lost stripes vs H3K27ac enrichment in gained stripes would directly support the Polycomb disruption model. This is the mechanistic evidence.

3. **GO/KEGG enrichment** (Section 10, `enrichment/go_bp_*.pdf`, `kegg_*.pdf`) — Whether anchor-proximal genes in lost vs gained stripes cluster in developmental/Polycomb-target pathways. The 250402 timepoint has enough genes (~hundreds of lost/gained) to power this; 250831 may fall below the 15-gene minimum and get skipped.

**Tier 2 — Validation/QC:**

4. **Stripiness scatter: ctrl vs mut** (Section 3) — Tests whether differential stripes have genuinely different stripiness scores between conditions, or if they're detection-frequency artifacts with similar signal strength. Given the uniformly small effect sizes (all |logFC| < 0.4), this is important for interpretation confidence.

5. **Replicate correlation heatmaps** (Section 7) — You already know BCV is very low (0.011-0.021), but the Spearman heatmaps will show whether that low BCV comes from genuinely tight replicates or from count inflation in the O_Sum_added metric.

6. **Early vs late combined panels** (Section 11, `combined/comparison_summary.*`) — The 250402 (31.5% sig) vs 250831 (2.4% sig) contrast is stark. The side-by-side will clarify whether early has consistent directionality at lower significance, or genuinely different biology.

**Tier 3 — Useful for JuiceBox follow-up:**

7. **Annotated 28-column BEDPEs** (Section 9) — The highconf/concordant tiers now include `nearest_gene`, `anchor_type`, and ChIP marks, making it much easier to prioritize which stripes to validate visually. These supersede the simple 15-column BEDPEs already generated.

## Analysis Prompt Assessment

The existing `stripenn-analysis-prompt.md` was designed for the **BEDPE files** and covers those columns well. But Stage 7 produces several data products the prompt doesn't address:

**Gaps to fill:**

| Missing from prompt                   | Where it lives                            | Why it matters                                                                                                                                                                 |
| ------------------------------------- | ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Stripiness scores (ctrl vs mut)       | `{tp}_annotated_stripes.tsv`              | Validates whether differential calls reflect signal differences or detection noise                                                                                             |
| GO/KEGG enrichment results            | `enrichment/go_bp_*.pdf` + TSVs           | Functional interpretation — the prompt asks for "biological interpretation" but has no section for pathway results                                                             |
| Replicate correlation structure       | `replicate_correlation_*.pdf`             | QC context for the edgeR results                                                                                                                                               |
| Early vs late comparison              | `combined/` directory                     | The prompt is single-timepoint; no section for cross-timepoint patterns                                                                                                        |
| `anchor_type` distribution statistics | `{tp}_annotated_stripes.tsv`              | The prompt asks for anchor_type breakdown but doesn't ask for enrichment *testing* (e.g., Fisher's exact: are lost stripes significantly enriched for Polycomb vs background?) |
| `h3k4me3` column                      | annotated TSV has it, simple BEDPEs don't | Prompt lists h3k27ac/h3k27me3/h3k4me1 but not h3k4me3 — relevant for Active_Promoter classification                                                                            |

**What I'd add to the prompt:**

A **Section 6: Anchor Annotation Analysis** that asks:
- Anchor type enrichment: is the distribution of anchor types in lost stripes significantly different from gained? From background (all stripes)?
- Which specific anchor types are over/under-represented in each direction?
- Bivalent promoter fate: are bivalent anchors preferentially lost or gained?

A **Section 7: Pathway Enrichment Summary** asking to report top GO BP and KEGG terms per direction if available.

A **Section 8: Cross-Timepoint Comparison** (when analyzing both) asking whether directional biases, anchor type enrichments, and pathway hits are consistent between early and late.

And update the column definitions table to include `h3k4me3` and `stripiness_ctrl`/`stripiness_mut`.

Want me to update the prompt file with these additions?