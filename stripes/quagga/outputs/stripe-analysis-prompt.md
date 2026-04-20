## Prompt for Claude: Differential Stripe Analysis

You are analyzing differential chromatin stripe results from a BAP1-KO vs wildtype Hi-C experiment in mouse cerebellum (mm10). The attached bedpe file contains stripes that have been processed through a detection + quantification + edgeR differential analysis pipeline.

## Experimental Context

**Perturbation**: BAP1 knockout (Polycomb regulator, deubiquitinates H2AK119ub1)
**Expected biology**: 
- Lost stripes may indicate disrupted Polycomb domain boundaries
- Gained stripes may indicate ectopic enhancer-promoter contacts at derepressed loci
**Replicates**: n=3 per condition (ctrl_M1-M3, mut_M1-M3)

## Column Definitions

| Column | Description |
|--------|-------------|
| chr1, x1, x2 | Stripe anchor (narrow, where CTCF/enhancer sits) |
| chr2, y1, y2 | Stripe span (extended region) |
| direction | `lost` (control-only) or `gained` (mutant-only) |
| direction_confidence | `high`, `medium`, or `low` |
| logFC | edgeR log2 fold-change (positive = higher in mutant) |
| FDR | Benjamini-Hochberg adjusted p-value |
| source | `control_only`, `mutant_only`, or `shared` |
| pval_ctrl / pval_mut | Quagga detection p-value (lower = more confident call) |
| detection_confidence | Based on replicate support + 10kb validation |
| in_10kb | Whether stripe validated at 10kb resolution |
| nearest_gene | Closest gene to anchor |
| distance_to_tss | Distance to TSS (0 = overlaps promoter) |
| anchor_type | `Promoter`, `Active_Enhancer`, or `Other` |
| h3k27ac, h3k27me3, h3k4me1 | ChIP-seq peak overlap at anchor (TRUE/FALSE) |
| stripe_length_kb | Span length in kb |
| anchor_width_kb | Anchor width in kb |

## Your Task

### 1. Summary Statistics
- Total stripes by direction (lost vs gained)
- Distribution by confidence tier
- Distribution by anchor_type
- ChIP-seq annotation breakdown

### 2. High-Priority Stripes for Validation

Identify the **top 8-10 stripes** to validate in JuiceBox, prioritizing by:
1. Detection confidence (high > medium > low)
2. Quagga p-value (lower = better detection) OR logFC
3. Biological relevance (promoter/enhancer anchors, relevant ChIP marks)
4. Interesting genes (if recognizable)

For each priority stripe, provide:
- **JuiceBox coordinate string**: `chr:start-end` format for pasting
- **Stripe type**: Lost or gained
- **Anchor info**: Gene name, anchor_type, ChIP marks
- **Rationale**: Why this stripe is interesting to validate

### 3. JuiceBox Coordinate Blocks

Provide copy-paste ready coordinate blocks:

```

# LOST stripes (check for stripe disappearance in mutant)

chr1:start-end chr2:start-end ...

# GAINED stripes (check for new stripe appearance in mutant)

chr3:start-end chr4:start-end ...

```

Use format: `chr:anchor_start-span_end` with ~50kb padding on each side for context.

### 4. Biological Interpretation

Given BAP1's role as a Polycomb deubiquitinase:
- Do gained stripes show H3K27ac enrichment (activation)?
- Do lost stripes show H3K27me3 (Polycomb marks)?
- Are there interesting gene associations?
- Any patterns in anchor_type distribution between lost vs gained?

### 5. Caveats and Confidence Assessment

Note any concerns about:
- Stripes with conflicting signals (e.g., logFC direction vs source)
- Low-confidence calls that appear biologically interesting
- Stripes that may need additional validation

## Output Format

Use tables and structured lists. For JuiceBox coordinates, provide them in a code block that can be directly copy-pasted.
