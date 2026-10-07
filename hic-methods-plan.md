# Hi-C Methods Section — Writing Plan

## Context

Challana asked Zara to write the methods section for the BAP1 Hi-C paper. The section covers loop analysis, stripe detection, differential analysis, TAD/compartment analysis, GO/KEGG enrichment, and all software versions. This is a new methods section written from scratch.

PAPER.md is a separate complete paper (CUT&RUN/epigenetics). Its Methods section sets the register and style but is not the Hi-C paper.

The verified source of truth for every parameter, version, and count is `loops-stats.md` (repo root). Every number in that file was verified against scripts, configs, conda environments, or pipeline output logs.

**Process.** I discuss content and suggest text. Zara writes the methods. We work one section at a time in document order, following the process from the grant-writing session (grant-writing.txt).

## Resolved decisions

| Decision | Answer |
|----------|--------|
| Em dashes | Banned. Use commas, periods, colons, parentheses. |
| Upstream paragraph | Two versions produced (full + one-sentence brief). Challana chooses. |
| Pipeline identity | nf-core/hic (v2.1.0). Confirmed by nextflow config. |
| TADs and compartments | Included. Two systems: HOMER + TADCompare. |
| GO/KEGG enrichment | Included in Hi-C methods. |
| Register | Past tense, third person, passive voice. Match PAPER.md methods. |
| STE | Max 25 words per descriptive sentence. No semicolons. Simple verbs. |
| Anti-AI guide | Full banned list. No negative parallelisms, fix narratives, AI vocabulary, hedging. |

## Section structure

Each subsection follows PAPER.md convention: short title (2-5 words), one paragraph, ordered by analytical workflow. The source for every item below is `loops-stats.md` with the relevant line range noted.

---

### 1. Hi-C data processing (lines 7-15)

Two versions. Challana chooses.

**Full version (~8 sentences):** fastp (v1.0.1) trimming. nf-core/hic (v2.1.0), Nextflow (v25.04.6), Bowtie2 (v2.5.4) to mm10, Arima digestion. Contact maps at 10 bin sizes (1 Mb to 1 kb). TADs called with insulation and hicExplorer at 40/20/10 kb. Compartments at 1000/250/100 kb. SAMtools (v1.22). Valid pairs to Juicer v2.0 .hic via juicerpre (VC, VC_SQRT, GW_SCALE, SCALE, INTER_SCALE). Singularity on SDSC Expanse.

**Brief version (1-2 sentences):** Names nf-core/hic (v2.1.0), Bowtie2 (v2.5.4), mm10, Juicer v2.0 .hic files. Refers reader to supplemental for pipeline details.

---

### 2. Chromatin loop detection (lines 17-21)

Short. 3-4 sentences.

Content: HiCCUPS from Juicer Tools (v2.0.1). Singularity container, GPU nodes. Default parameters, no explicit flags. Three resolutions (5/10/25 kb). Six replicate .hic files processed separately.

---

### 3. Differential loop analysis (lines 23-97)

Two paragraphs. The most detailed section.

**Paragraph 1 — mariner processing (lines 23-57):**
mariner (v1.8.1). mergePairs: 10 kb radius, Manhattan distance, selectMax = TRUE. assignToBins at resolution, center positioning. pixelsToMatrices: buffer = 2 (5x5). pullHicMatrices: VC normalization, observed matrix, HDF5. Aggregation: sum of 25 pixels. chrX excluded.

**Paragraph 2 — edgeR QL-GLM (lines 59-87):**
edgeR (v4.4.0). filterByExpr: min.count = 5, min.total.count = 15, min.prop = 0.5. TMM via normLibSizes(). Design: ~group, no covariates. estimateDisp(robust = TRUE). glmQLFit(robust = TRUE). glmQLFTest(coef = 2). BH FDR. Primary: FDR < 0.05. Final stringent: |logFC| > 0.3 AND FDR < 0.03.

**Stays out:** Effect-size categories (results). QC thresholds (internal). Center-weighted aggregation (not used). "KR not available" (fix narrative).

---

### 4. Multi-resolution loop comparison (lines 99-103)

Short. 2-3 sentences.

Content: Matched across 5/10/25 kb. Same chromosome, both anchors within 10 kb. Closest by summed anchor distance.

---

### 5. Loop anchor annotation (lines 105-150)

One paragraph. ~6 sentences.

Content: 8-category priority system. Five histone marks (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent) + CTCF motifs. Priority order with rules (table at lines 121-132). TSS threshold: 2 kb (TxDb v3.10.0). Gene annotation: GENCODE vM25, nearest gene. Loop type: pairwise label. Loop class: structural/CRE/mixed/unclassified.

**Decision to make during writing:** Present the 8 categories as prose (inline) or reference a supplemental table. PAPER.md uses prose for similar classifications.

---

### 6. Aggregate peak analysis (lines 152-160)

Short. 2-3 sentences.

Content: GENOVA (v1.0.1). Merged per-condition contact maps. Not replicate-aware: single ctrl-vs-mutant comparison.

**Do NOT describe** the mariner-based apa_analysis.R. The paper uses GENOVA results.

---

### 7. Chromatin stripe detection (lines 289-331)

One paragraph. ~6 sentences.

Content: hictk (v2.1.4) conversion .hic to .mcool at 5 kb and 10 kb. Stripenn (v1.1.65.22) on merged .mcool files. KR normalization. Parameters: min length 10 bins, max width 8 bins, maxpixel percentiles 0.95-0.99, p-value 0.1, seed 42. chrX/Y/M removed post-hoc. Union set: 50 kb anchor tolerance, direction-aware matching. Replicate scoring: stripenn score on 6 replicates, O_Sum_added as count.

**Stays out:** The 250402 KR balancing fix (pipeline maintenance). Default seed number comparison (irrelevant). Canny edge detection algorithm description (Stripenn's internals, not our parameter).

---

### 8. Differential stripe analysis (lines 333-399)

One paragraph. ~8 sentences.

Content: edgeR QL-GLM, same design as loops. Two stated differences: (1) filterByExpr not applied, (2) calcNormFactors used. Same dispersion, fit, test. Primary: FDR < 0.05. Direction classification: lost/gained/strengthened/weakened/unchanged (state the rules). Cross-resolution: 50 kb tolerance, concordance = significant at both with same direction, Pearson logFC correlation.

Stripe annotation (lines 381-399): 7-category system (no CTCF_Site). Same marks, same 2 kb TSS threshold. Two differences from loops: (1) Bivalent_Promoter priority 1 vs 3, (2) Active_Promoter requires H3K27ac+.

**Stays out:** Confidence tiering details (QC metric). Directional consistency flag (QC). BCV and max |logFC| (results).

---

### 9. TAD analysis (lines 168-239)

Two paragraphs. Two independent systems.

**Paragraph 1 — HOMER TADs (lines 172-189):**
HOMER (v5.1). Tag directories from BAM files. Four steps: findtad per sample, mergetad consensus, scoretad inclusion ratios, getDiffExpression.pl differential. Volcano thresholds: relaxed (FDR < 0.15, |Diff| > 0.15), standard (FDR < 0.05, |Diff| > 0.30).

**Paragraph 2 — TADCompare boundaries (lines 191-239):**
TADCompare (v1.16.0). Resolution: 10 kb. Autosomes only (chr1-chr19). Input: sparse matrices from straw, NONE normalization. Six boundary classifications (Non-Differential, Strength Change, Shifted, Split, Merge, Complex). Blacklist: lab-specific (254 regions) + ENCODE mm10 (3,435 regions), union. Shift analysis: bedtools closest, minimum 10 kb. Boundary-loop cross-reference: Fisher's exact, permutation (1,000 iterations), distance thresholds 10/25/50/100 kb. Boundary ChIP classification: same 8-category system as loop annotation.

---

### 10. Compartment analysis (lines 181-187)

Short. 3-4 sentences. Can fold into section 9 or stand alone.

Content: HOMER PC1 at 25x50kb resolution per sample. Differential via getDiffExpression.pl. Same dual-threshold volcano. PC1 polarity auto-correction (genic bins validated for directionality). Genome percentage computed against mm10 assembled genome (2,730,871,774 bp).

---

### 11. Gene ontology and pathway enrichment (lines 162-166)

Short. 3-4 sentences.

Content: clusterProfiler (v4.16.0). Entrez IDs within 10 kb of up/down loop anchors. enrichGO (BP, CC, MF) and enrichKEGG. org.Mm.eg.db (v3.21.0). BH correction, p < 0.05 and q < 0.05.

---

## Writing process

1. Work through sections 1-11 in order.
2. Per section: I state content + stays-out + key parameters. We discuss. I suggest text. You edit. Move on.
3. One section per turn. Short adjacent sections (for example 4 + 5) can combine.
4. After all sections: consistency review (terminology, tense, style).

## Style rules (applied to every sentence)

- Past tense, third person, passive voice
- Max 25 words per sentence
- No em dashes. No semicolons.
- No AI vocabulary (delve, leverage, robust, crucial, seamless, etc.)
- No negative parallelisms
- No fix narratives (VC is the method, period)
- No hedging (may, could, might)
- No over-justification of standard methods
- No numbers substituting for plain facts
- "Pipeline" or tool name, never "framework"
- Parameters stated, not dramatized
- Result counts belong in Results, not Methods

## Verification (after all sections drafted)

1. Every parameter cross-checked against loops-stats.md
2. Every version cross-checked against software versions table (lines 480-539)
3. Sentence length audit (25-word max)
4. Banned word scan (anti-ai-writing-style-guide.md)
5. Tense consistency (past throughout)
6. Term consistency (each tool named identically everywhere)
7. Loop annotation = 8 categories. Stripe annotation = 7 categories. Both correct.
8. APA = GENOVA, not mariner apa_analysis.R
9. edgeR design = ~group, no covariates

## Unresolved items

1. **Sex composition of Hi-C samples.** Not in configs. Challana has lab records. Does not affect methods text (design has no sex covariate).
2. **Annotation table format.** Prose vs supplemental table for the 8-category system. Decide at section 5.
3. **Upstream version: full vs brief.** Challana chooses.
