## Prompt for Claude: Stripenn Differential Stripe Analysis

You are analyzing differential chromatin stripe results from a BAP1-KO vs wildtype Hi-C experiment in mouse cerebellum (mm10). The attached BEDPE file contains stripes processed through a Stripenn detection + edgeR differential analysis pipeline with cross-resolution validation.

## Experimental Context

**Perturbation**: BAP1 knockout (Polycomb regulator, deubiquitinates H2AK119ub1)
**Expected biology**:
- Lost stripes may indicate disrupted Polycomb domain boundaries or cohesin loading sites
- Gained stripes may indicate ectopic enhancer-promoter contacts at derepressed loci
**Replicates**: n=3 per condition (ctrl_M1-M3, mut_M1-M3)
**Stripe caller**: Stripenn v1.1.65 (Canny edge detection + pixel saturation thresholds)
**Quantification**: Stripenn `score` command (O_Sum_added = total observed contacts per stripe)

## Column Definitions

| Column | Description |
|--------|-------------|
| chr1, x1, x2 | Stripe anchor (narrow, where CTCF/enhancer sits) |
| chr2, y1, y2 | Stripe span (extended region) |
| direction | `lost` (control-only) or `gained` (mutant-only) |
| direction_confidence | `high`, `medium`, or `low` (based on FDR + logFC agreement) |
| logFC | edgeR log2 fold-change (positive = higher in mutant) |
| FDR | Benjamini-Hochberg adjusted p-value |
| source | `control_only`, `mutant_only`, or `shared` |
| pval_ctrl / pval_mut | Stripenn detection p-value (lower = more confident call) |
| detection_confidence | Same as direction_confidence (Stripenn has no separate detection tier) |
| in_10kb | Whether stripe matched at 10kb resolution (from cross-resolution comparison) |
| nearest_gene | Closest gene to anchor |
| distance_to_tss | Distance to TSS in bp (0 = overlaps promoter) |
| anchor_type | `Active_Promoter`, `Repressed_Promoter`, `Bivalent_Promoter`, `Active_Enhancer`, `Poised_Enhancer`, `Polycomb`, or `Other` |
| h3k27ac, h3k27me3, h3k4me1, h3k4me3 | ChIP-seq peak overlap at anchor (TRUE/FALSE) |
| stripe_length_kb | Span length in kb |
| anchor_width_kb | Anchor width in kb |

**Additional columns in annotated TSV** (not in simple BEDPE, available in `{tp}_annotated_stripes.tsv`):

| Column | Description |
|--------|-------------|
| stripiness_ctrl | Stripiness score in control (higher = stronger stripe signal) |
| stripiness_mut | Stripiness score in mutant |
| bivalent | Overlap with pre-computed H3K4me3+H3K27me3 intersection peaks (TRUE/FALSE) |
| nearest_gene_id | Entrez gene ID for nearest gene |

## Stripenn-Specific Notes

- **Stripiness scores** (not in BEDPE but in source TSV): Measures the "stripe-ness" of the Hi-C signal. Negative values can occur after normalization and indicate absence of stripe signal.
- **direction_type** (3prime/5prime): Indicates the orientation of the stripe relative to the diagonal. Not directly strand but reflects geometric direction of signal extension.
- **resolution_support** (in cross_res_merged.tsv): `both_concordant` = highest confidence, `both_discordant` = significant at both but conflicting direction, `5kb_only` / `10kb_only` = single-resolution.
- **Effect sizes are uniformly small** in this dataset (all |logFC| < 0.4). The differential signal is driven by frequency (many stripes with small changes) rather than magnitude.

## Your Task

### 1. Summary Statistics
- Total stripes by direction (lost vs gained)
- Distribution by confidence tier
- Distribution by anchor_type
- ChIP-seq annotation breakdown (H3K27ac, H3K27me3, H3K4me1, H3K4me3 overlap rates by direction)

### 2. High-Priority Stripes for Validation

Identify the **top 8-10 stripes** to validate in JuiceBox, prioritizing by:
1. Cross-resolution validation (in_10kb == TRUE or resolution_support == "both_concordant")
2. Direction confidence (high > medium > low)
3. FDR (lower is better)
4. Biological relevance (Active_Promoter, Polycomb, or Bivalent_Promoter anchors)
5. logFC magnitude (even though all are small, relative ranking is informative)

For each priority stripe, provide:
- **JuiceBox coordinate string**: `chr:start-end` format with ~50kb padding
- **Stripe type**: Lost or gained
- **Anchor info**: Gene name, anchor_type, ChIP marks
- **Rationale**: Why this stripe is interesting to validate

### 3. JuiceBox Coordinate Blocks

Provide copy-paste ready coordinate blocks:

```
# LOST stripes (check for stripe disappearance in mutant)
chr:start-end   # gene: name
...

# GAINED stripes (check for new stripe appearance in mutant)
chr:start-end   # gene: name
...
```

Use format: `chr:anchor_start-span_end` with ~50kb padding on each side.

### 4. Biological Interpretation

Keep speculation minimal. Report what the data shows:
- Direction bias (more lost or gained?)
- Anchor type enrichment differences between lost vs gained
- ChIP mark overlap patterns
- Any notable gene associations

### 5. Anchor Annotation Analysis

Using the annotated TSV (`{tp}_annotated_stripes.tsv`) or the annotated BEDPEs:

**Anchor type enrichment:**
- Compare anchor_type distribution in lost stripes vs gained stripes vs all stripes (background)
- Flag any anchor types significantly over- or under-represented in a direction (Fisher's exact or chi-squared)
- Specifically assess: are Polycomb/Bivalent_Promoter anchors preferentially lost? Are Active_Enhancer anchors preferentially gained?

**ChIP-seq mark analysis:**
- Report H3K27me3, H3K27ac, H3K4me1, H3K4me3 overlap rates per direction as a table
- Compare lost vs gained overlap rates — which marks show the largest differential?
- Note any stripes with bivalent marks (H3K4me3 + H3K27me3) — these are key BAP1/Polycomb targets

**Stripiness score comparison** (if annotated TSV provided):
- For shared stripes: is stripiness_ctrl vs stripiness_mut correlated? Do differential stripes show genuinely different stripiness, or are effect sizes driven by count differences at similar stripe morphology?

### 6. Pathway Enrichment Summary

If GO/KEGG enrichment results are available (from `combined/enrichment/`):

**GO Biological Process** (`go_bp_{tp}.tsv` / `go_bp_dotplot_{tp}.*`):
- Report top 5 terms per direction (lost vs gained), with gene counts and adjusted p-values
- Flag any Polycomb-related, developmental, or chromatin-remodeling terms

**GO Cellular Component** (`go_cc_{tp}.tsv` / `go_cc_dotplot_{tp}.*`):
- Report top 3 terms per direction — are differential stripes enriched at specific subcellular structures (nuclear envelope, heterochromatin, etc.)?

**GO Molecular Function** (`go_mf_{tp}.tsv` / `go_mf_dotplot_{tp}.*`):
- Report top 3 terms per direction — any enrichment for transcription factor binding, chromatin binding, or deubiquitinase targets?

**KEGG Pathways** (`kegg_{tp}.tsv` / `kegg_dotplot_{tp}.*`):
- Report top 3 pathways per direction if significant

**General:**
- Note if enrichment was underpowered (too few genes) for either direction or timepoint
- Compare enrichment results across timepoints if both are available

### 7. Cross-Timepoint Comparison

When analyzing both early (250831/P12) and late (250402/adult) timepoints:

- Is the directional bias consistent? (Both show more lost? More gained? Or opposite?)
- Are anchor type enrichments consistent between timepoints?
- Do the same pathways appear in both timepoints' enrichment results?
- Identify any stripes at the same locus that are differential in both timepoints (convergent evidence)
- Note the large difference in significance rates (early ~2.4% vs late ~31.5%) and what this implies for developmental timing of BAP1's effect

### 8. Caveats and Confidence Assessment

Note any concerns about:
- Small effect sizes (all |logFC| < 0.4) — the differential signal is frequency-driven, not magnitude-driven
- Directional consistency rate
- Stripes with conflicting signals (e.g., both_discordant in cross-resolution)
- Low-confidence calls that appear biologically interesting
- Whether anchor type enrichments survive multiple testing correction
- Replicate correlation structure (tight within-group → trustworthy; inflated → cautious interpretation)

## Output Format

Use tables and structured lists. For JuiceBox coordinates, provide them in a code block that can be directly copy-pasted.
