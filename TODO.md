# Hi-C Remaining Analyses TODO

Reference: `CT-meeting-3.md` (meeting notes)

---

## Status Legend

- `[x]` = Done (script exists and has been run)
- `[~]` = Partially done (infrastructure exists, analysis incomplete)
- `[ ]` = Not started

---

## 1. Shared Anchor / Loop Switching Analysis

**Meeting notes:** "loop switching mechanism", "long range to short range polycomb loops", "if happening at same sites", "have to prove same exact sites that lost are same that gained", "shared anchor - might change", "long range anchors - with short range gained"

**Core question:** Do anchors that participate in lost long-range loops also participate in gained short-range loops? Prove loop "switching" at the same genomic sites.

### Tasks

- [x] **1a. Identify shared anchors between lost and gained loops.** For each anchor coordinate (with some tolerance window, e.g. 10kb), find anchors that appear in both a `down_in_mutant` loop AND an `up_in_mutant` loop. These are candidate "switching" anchors. -> `scripts/shared_anchor_analysis.R` (212 shared anchors identified)
- [x] **1b. Characterize shared anchors vs non-shared anchors.** Compare chromatin state (Polycomb, Active_Promoter, etc.), distance properties, and ChIP-seq marks between shared and non-shared anchors. -> `output/shared_anchor_analysis/{early,late}/` (Chi-square p=3.28e-31)
- [x] **1c. Show that lost loops at shared anchors are longer than gained loops at the same anchors.** Paired analysis: for each shared anchor, compare the distance of the lost loop vs the gained loop. -> Paired Wilcoxon p=1.17e-20, median lost=1.15Mb, median gained=340kb
- [x] **1d. Generate aggregate heatmaps (APA) for shared-anchor loop subsets.** Separate APA for: (a) lost long-range loops at shared anchors, (b) gained short-range loops at shared anchors. -> `scripts/apa_shared_anchors.R`, `scripts/apa_shared_anchors.sb`
- [x] **1e. Violin plot: expression of genes near shared anchors.** For genes proximal to shared anchors, show log2FC distribution (from RNA-seq) split by whether the anchor lost a long-range loop vs gained a short-range loop. -> `output/shared_anchor_analysis/{early,late}/plots/`

### Existing resources

- **Anchor coordinates:** `25042-late_outputs/merged_loops/characterized_loops.tsv` (columns: `anchor1_chr`, `anchor1_start`, `anchor1_end`, `anchor2_*`, `direction`, `loop_distance`)
- **Early equivalent:** `250831-early_outputs/merged_loops/characterized_loops.tsv`
- **Extended annotation (7 chromatin states):** `scripts/annotate_loops_extended.R` -> outputs in `outputs/loop_annotation_extended/{early,late}/extended_characterized_loops.tsv`
- **APA pipeline:** `scripts/apa_analysis.R` (extract Hi-C matrices, aggregate, compute P2LL scores)
- **RNA-seq DEG data:** `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` (late), `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx` (early)
- **Gene-boundary association pattern:** `tads/scripts/deg_tad_violin.R` (GREAT-style regulatory domains, can be adapted for loop anchors)

### Notes

- The `deg_tad_violin.R` script in `tads/scripts/` already implements GREAT-style gene association and violin plotting for TAD boundaries. The same pattern can be reused for loop anchors.
- Shared anchor identification needs a tolerance window since loops come from different resolutions (5kb, 10kb, 25kb).

---

## 2. RNA-seq Integration with Differential Loops

**Meeting notes:** "violin plot - genes around lost long range genes gained short", "are expression of those changing"

**Core question:** Are genes near differentially looped regions also differentially expressed?

### Tasks

- [ ] **2a. Associate genes with loop anchors using GREAT-style regulatory domains.** Map each loop anchor to its proximal genes (5kb upstream, 1kb downstream of TSS, 100kb max extension).
- [ ] **2b. Violin plot: log2FC of DEGs near lost vs gained loop anchors.** Analogous to `deg_tad_violin.R` but for loop anchors instead of TAD boundaries.
- [ ] **2c. Stratify by distance category.** Separate analysis for genes near long-range lost loops (>500kb) vs short-range gained loops (<500kb).
- [ ] **2d. Stratify by anchor chromatin state.** Show expression changes specifically for genes near Polycomb-classified anchors.

### Existing resources

- **RNA-seq data:** `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` (late/adult timepoint, contains columns: `ensembl_gene_id` [actually gene symbols], `log2FoldChange`, `padj`, `baseMean`)
- **RNA-seq data:** `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx` (early/young timepoint)
- **Reference implementation:** `tads/scripts/deg_tad_violin.R` (complete working example of GREAT-style gene association + violin plot + Wilcoxon test)
- **Gene annotations:** `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`
- **Loop data with gene annotations:** `characterized_loops.tsv` already has `anchor1_nearest_gene`, `anchor2_nearest_gene` columns

---

## 3. Polycomb Loop Story

**Meeting notes:** "k27me3 forming polycomb loops - with closer contact in mutant", "not interacting as much with closer contact", "polycomb story", "as lose k27me3 - gaining contacts - hopefully", "characterize anchors lost vs gained that are polycomb associated", "thru the genes, heatmap"

**Core question:** Are H3K27me3-marked (Polycomb) anchors specifically involved in the long-to-short loop switching? Do they gain closer contacts while losing distal ones?

### Tasks

- [x] **3a. CDF and density plots filtered by H3K27me3 overlap.** Three subsets: K27me3-anchored (at least one anchor), K27me3-both (both anchors), bivalent. -> `scripts/loop_distance_k27me3_filtered.R`
- [x] **3b. ChIP-seq mark trends across loop distance.** H3K27ac decreasing, H3K27me3 increasing with distance. -> `scripts/chip_distance_analysis.R`
- [x] **3c. Polycomb-specific shared anchor analysis.** Among shared anchors (task 1a), filter to those classified as "Polycomb" or "Repressed_Promoter" in extended annotation. Show that these are the anchors driving the switching. -> `scripts/polycomb_shared_anchor_analysis.R`, outputs in `output/shared_anchor_analysis/{early,late}/polycomb_specific/`
- [ ] **3d. APA heatmaps for Polycomb-anchored loops.** Aggregate contact signal at Polycomb-classified loop anchors, split by lost vs gained.
- [ ] **3e. Gene body heatmap for Polycomb-associated loops.** For genes at Polycomb anchors that switch from long to short loops, show a heatmap of ChIP-seq signal (H3K27me3) across the gene body.
- [x] **3f. Assess differential H3K27me3/H2AK119ub at Polycomb shared anchors.** -> `scripts/diff_chip_polycomb_enrichment.R`, `output/diff_chip_polycomb_enrichment/` Determine whether K27me3 and/or H2AK119ub signal is changing within the Polycomb-anchored shared anchor regions identified in 3c.
  - **Initial check (easy):** Overlap diffbind differential peak regions with: (a) long-range lost loops, (b) gained short-range loops, (c) unchanged loops. Caveat: summit=400bp may miss broader domains → expect false negatives, but worth a shot.
  - **Input:** Differential peaks from diffbind in `peaks/new/`. Individual peak files also on instance at `/data2/rs_256/Func_annotation_v2/subtracted_bedfiles`. BigWigs in `heatmaps/` for aggregate heatmaps/visualization.
- [ ] **3g. Aggregate ChIP-seq signal heatmaps at loop anchors.** Follow-up to 3f using bigWigs (preferred for visualizing broader H3K27me3/H2AK119ub domains beyond 400bp summits).

### Existing resources

- **K27me3 filtered analysis:** `scripts/loop_distance_k27me3_filtered.R` -> `output/loops_k27me3_filtered/{early,late}/`
- **K27me3 global analysis:** `scripts/loop_distance_k27me3_global.R` -> `output/loops_k27me3_global/{early,late}/`
- **ChIP distance analysis:** `scripts/chip_distance_analysis.R` -> `output/chip_distance/{early,late}/`
- **H3K27me3 peak files:** `peaks/beds/H3K27me3CerebellumLate1.bed`, `peaks/beds/H3K27me3CerebellumEarly1.bed`
- **Bivalent peak files:** `peaks/beds/Bivalent_Cerebellum_Late.bed`, `peaks/beds/Bivalent_Cerebellum_Early.bed`
- **Extended annotation with Polycomb category:** `scripts/annotate_loops_extended.R` (classifies anchors as Polycomb when H3K27me3+ and >2kb from TSS)
- **APA pipeline:** `scripts/apa_analysis.R`
- **Polycomb shared anchor analysis:** `scripts/polycomb_shared_anchor_analysis.R` -> `output/shared_anchor_analysis/{early,late}/polycomb_specific/`
- **Differential H3K27me3 peaks (diffbind, summit=400bp):** `peaks/new/adult_K27me3_down.bed`, `peaks/new/adult_K27me3_up.bed` (late), `peaks/new/P12_H3K27me3_down.bed`, `peaks/new/P12_H3K27me3_up.bed` (early)
- **Differential H2AK119ub peaks (diffbind, summit=400bp):** `peaks/new/H2AK119ub_down.bed`, `peaks/new/H2AK119ub_up.bed`

---

## 4. Anchor-Specific ChIP-seq Subsetting

**Meeting notes:** "do cdf density - for each histone modification", "k27ac on one site vs both sites (super enhancer)", "+ k4me3 - enhancer+promoter loop"

**Core question:** How does the distance distribution differ when filtering loops by specific ChIP-seq mark combinations at anchors?

### Tasks

- [x] **4a. CDF/density by H3K27me3.** One anchor vs both anchors. -> `scripts/loop_distance_k27me3_filtered.R`
- [x] **4b. CDF/density by H3K27ac.** One anchor (enhancer-only loop) vs both anchors (super-enhancer loop). -> `scripts/loop_distance_mark_filtered.R --marks H3K27ac` (92% of super-enhancer loops are lost in BAP1-KO)
- [x] **4c. CDF/density for H3K27ac + H3K4me3 (enhancer-promoter loops).** Filter to loops where one anchor is H3K27ac+ (enhancer) and the other is H3K4me3+ (promoter). Compare lost vs gained distance distributions. -> `scripts/loop_distance_ep_filtered.R`
- [x] **4d. CDF/density for each individual histone modification separately.** All 5 marks (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent) with one-anchor and both-anchor filters. -> `scripts/loop_distance_mark_filtered.R --marks all --timepoint both` outputs in `output/loops_mark_filtered/{early,late}/{mark}/`

### Existing resources

- **Template script:** `scripts/loop_distance_k27me3_filtered.R` (complete CDF + density + stats pipeline, can be generalized)
- **ChIP-seq peaks:** `peaks/beds/H3K27acCerebellumLate2.bed`, `H3K4me1CerebellumLate1.bed`, `H3K4me3CerebellumLate2.bed` (and Early equivalents)
- **Extended annotation:** `outputs/loop_annotation_extended/{early,late}/extended_characterized_loops.tsv` has all overlap columns (`anchor1_H3K27ac_overlap`, `anchor1_H3K27me3_overlap`, `anchor1_H3K4me1_overlap`, `anchor1_H3K4me3_overlap`, `anchor1_Bivalent_Promoter_overlap`, and anchor2 equivalents)
- **Multi-format output utility:** `scripts/utils/multi_format_output.R`

---

## 5. CTCF / Cohesin / Loop Extrusion Analysis

**Meeting notes:** "ctcf: not there anymore - losing checkpoint", "extrusion - rad21 faulty", "rad21 likely defunct", "more stripe level defects", "jesse dickson at salk - ctcf", "rad21 - ask for assay"

**Core question:** Is CTCF boundary function compromised in BAP1-KO? Is cohesin (RAD21) extrusion faulty, leading to stripe-level defects?

### Tasks

- [x] **5a. CTCF overlap analysis at lost vs gained loop anchors.** Are lost loops more likely to have CTCF at their anchors? Are gained loops CTCF-depleted? -> `scripts/loop_distance_mark_filtered.R --marks ctcf` outputs in `output/loops_mark_filtered/late/ctcf/`
- [x] **5b. CDF/density for CTCF-anchored loops.** Filter to loops with CTCF at anchors, compare lost vs gained distance distributions. -> `output/loops_mark_filtered/late/ctcf/` (one-anchor and both-anchor filters)
- [x] **5c. Cross-reference with stripe analysis.** Do regions with lost CTCF-anchored loops show stripe defects? -> `scripts/ctcf_stripe_crossref.R`, `output/ctcf_stripe_crossref/{early,late}/` (Late: 17/890 CTCF loops at lost stripes, Fisher's p=1.0, OR=1.14; no significant enrichment - loop loss and stripe defects appear independent)
- [ ] **5d. Obtain RAD21 ChIP-seq data.** No RAD21 data currently in repo. Need to request assay or find public data.
- [ ] **5e. (If RAD21 data obtained) Overlap RAD21 with loop anchors and TAD boundaries.**

### Existing resources

- **CTCF ChIP-seq peaks:** `peaks/CTCF.bed` (32,487 peaks, Late timepoint)
- **CTCF motif predictions:** `peaks/ctcf_motifs_mm10.bed` (genome-wide motif scan)
- **Extended annotation already includes CTCF:** `scripts/annotate_loops_extended.R` has CTCF_Site category (8-category version in `peaks/loop_annotation_extended/`)
- **CTCF loop distance analysis:** `output/loops_mark_filtered/late/ctcf/` (CDF/density plots for one-anchor and both-anchor CTCF loops)
- **Stripe pipeline:** `stripes/scripts/phase1_detection.R` through `phase4_integration.R`, with outputs in `stripes/outputs/{early,late}/`
- **RAD21 data:** NOT AVAILABLE - need to request (meeting note: "rad21 - ask for assay", "rao - ask about tet enzyme assay")

### Action items (external)

- Contact lab about RAD21 ChIP-seq assay
- Ask Rao about tet enzyme assay
- Look into Jesse Dickson (Salk) CTCF work for methodology reference

---

## 6. TAD Boundary Integration

**Meeting notes:** "tad boundary changing?"

### Tasks

- [x] **6a. TADCompare differential boundary analysis.** -> `tads/scripts/02_run_tadcompare.R` through `05_filter_blacklist.R`
- [x] **6b. DEG violin plot at TAD boundaries.** -> `tads/scripts/deg_tad_violin.R`
- [~] **6c. Cross-reference differential loops with differential TAD boundaries.** Are differential loops preferentially located near differential TAD boundaries? Overlap analysis between `characterized_loops.tsv` anchors and `tads/results/{tp}/final/tadcompare_final_filtered.tsv` boundaries. -> `tads/scripts/boundary_loop_crossref.R` (Late: 69.6% concordance, p<0.001; Merge boundaries enriched in lost loops OR=0.32, Strength Change enriched in gained loops OR=2.11). **Permutation test needs redo with regioneR/regioneReloaded.**
- [x] **6d. Proper permutation analysis for boundary-loop enrichment.** *(CTea handling)* Redo permutation testing using regioneR/regioneReloaded (Bioconductor). Current implementation is basic; proper approach shuffles one BED file across genome ≥1000 times (10,000 for final), measures overlap, builds null distribution. Consider restricting background set if asking about co-enrichment of specific chromatin features (e.g., euchromatin-associated markers → restrict to euchromatin).

### Existing resources

- **TAD results:** `tads/results/{early,late}/final/tadcompare_final_filtered.tsv` (differential boundaries with Enriched_In, Type, Gap_Score)
- **TAD visualization:** `tads/scripts/tad_visualizations.R` (40+ plots across 9 subdirectories)
- **Loop data:** `25042-late_outputs/merged_loops/characterized_loops.tsv`, `250831-early_outputs/merged_loops/characterized_loops.tsv`
- **Permutation analysis reference:** regioneReloaded vignette (https://www.bioconductor.org/packages/release/bioc/vignettes/regioneReloaded/inst/doc/regioneReloaded.html) - provides nice visuals for group BED file comparisons

---

## 7. ABC Model / Enhancer-Gene Linkage

**Meeting notes:** "tie enhancers to respective genes and link to expression changes", "activity by contact model (abc model)", "use abc model - formula", "enhancer gene linkages dysregulated tied to diff expressed genes", "activity times contact", "for every gene - difference in enhancer promoter contacts", "for promoter diff contact value", "differential gene expression", "super enhancers to gene due k27ac", "find a tool for this"

**Reference paper:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6886585/#F3

**Core question:** Can we use the ABC model (Activity x Contact) to link enhancers to genes, then show that enhancer-gene linkages are dysregulated in BAP1-KO and tied to differentially expressed genes?

### Tasks

- [ ] **7a. Implement ABC model scoring.** ABC score = (enhancer Activity) x (Contact frequency). Activity = H3K27ac signal at enhancer. Contact = Hi-C contact frequency between enhancer and promoter (from loop logFC or raw counts).
- [ ] **7b. For every gene, compute delta enhancer-promoter contact.** Using the differential loop data, calculate change in contact between each promoter and its linked enhancers.
- [ ] **7c. Correlate delta contacts with differential gene expression.** Plot delta E-P contact vs log2FC from RNA-seq.
- [ ] **7d. Identify dysregulated enhancer-gene linkages.** Find E-P pairs where both contact and expression change significantly.
- [ ] **7e. Super-enhancer to gene linkages.** Identify loops where both anchors are H3K27ac+ (super-enhancer signature), link to target genes, assess expression changes.

### Existing resources

- **H3K27ac peaks (enhancer activity proxy):** `peaks/beds/H3K27acCerebellumLate2.bed`, `H3K27acCerebellumEarly2.bed`
- **Loop contact data with logFC:** `characterized_loops.tsv` (has `logFC`, `anchor1_H3K27ac_overlap`, `anchor2_H3K27ac_overlap`, `anchor1_type`, `loop_type`)
- **RNA-seq:** `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`, `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx`
- **Gene annotations:** Already in characterized_loops.tsv (`anchor1_nearest_gene`, `anchor2_nearest_gene`, `anchor1_distance_to_tss`)
- **Reference paper:** Fulco et al. 2019 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6886585/) - ABC model methodology

### Notes

- Meeting notes: "h3k27ac hi-chip - not working" so Hi-C + H3K27ac peak overlap is the alternative approach
- Meeting notes distinguish two goals:
  - **1st goal:** Link change in E-P contacts to differentially expressed genes
  - **2nd goal:** Tie change in delta contacts to ubiquitinated histone (H2AK119ub)
- **H2AK119ub ChIP-seq data now available:** Differential peaks from diffbind in `peaks/new/` (up/down in mutant)
- Meeting hypothesis: "ub is buffer to stop k27ac contact" - once ubiquitination threshold is reached, contacts form

---

## 8. Ubiquitination (H2AK119ub) Integration

**Meeting notes:** "contact k27ac and ub - mediated by which?", "hypothesis: ub is buffer to stop k27ac contact", "once reach threshold - now forms contact", "correlate - enhancer gene linkage level", "how much ubiquitinated histone", "2nd: tie change in delta contacts to ubiquitinated histone"

### Tasks

- [x] **8a. Obtain H2AK119ub ChIP-seq data.** Differential peaks now available in `peaks/new/`.
- [ ] **8b. Overlap H2AK119ub with loop anchors.** Classify anchors by ubiquitination status. Use differential peaks from diffbind (up/down in mutant).
- [ ] **8c. Correlate delta E-P contacts with ubiquitination levels.** Test hypothesis that ubiquitination buffers H3K27ac-mediated contact formation.

### Existing resources

- **H2AK119ub differential peaks (diffbind):** `peaks/new/H2AK119ub_down.bed`, `peaks/new/H2AK119ub_up.bed` (summit=400bp, may underestimate broader domains)
- **Framework for ChIP overlap:** `scripts/annotate_loops_extended.R` (can add H2AK119ub as additional mark)

---

## Summary: Priority Order

Based on meeting notes emphasis and data availability:

| Priority | Analysis | Data Available? | Script Exists? |
|----------|----------|----------------|----------------|
| 1 | Shared anchor / loop switching (Section 1) | Yes | **COMPLETE** (`scripts/shared_anchor_analysis.R`, `scripts/apa_shared_anchors.R`) |
| 2 | RNA-seq integration with loops (Section 2) | Yes | No (TAD version exists as template) |
| 3 | Polycomb loop story completion (Section 3) | Yes | Partially (3c done: `scripts/polycomb_shared_anchor_analysis.R`, 3d-3f pending) |
| 4 | Per-mark CDF/density subsetting (Section 4) | Yes | **COMPLETE** (`scripts/loop_distance_mark_filtered.R` all 5 marks + CTCF) |
| 5 | CTCF analysis (Section 5) | Yes (CTCF peaks), No (RAD21) | Partially (5a-5b done, 5c-5e pending) |
| 6 | TAD-loop cross-reference (Section 6) | Yes | Partially (permutation needs redo with regioneR/regioneReloaded) |
| 7 | ABC model / E-P linkage (Section 7) | Yes (H2AK119ub now available) | No |
| 8 | H2AK119ub integration (Section 8) | Yes (differential peaks in `peaks/new/`) | No |

---

## Data File Reference

### Loop data (primary inputs for new analyses)

| Timepoint | File | Key columns |
|-----------|------|-------------|
| Late | `25042-late_outputs/merged_loops/characterized_loops.tsv` | 57 cols: coordinates, logFC, FDR, direction, distance, 5-mark ChIP overlaps (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent), 7-category anchor_type, gene annotations, loop_type |
| Early | `250831-early_outputs/merged_loops/characterized_loops.tsv` | Same schema (57 cols) |

### ChIP-seq peaks

| Mark | Late | Early |
|------|------|-------|
| H3K27ac | `peaks/beds/H3K27acCerebellumLate2.bed` | `peaks/beds/H3K27acCerebellumEarly2.bed` |
| H3K27me3 | `peaks/beds/H3K27me3CerebellumLate1.bed` | `peaks/beds/H3K27me3CerebellumEarly1.bed` |
| H3K4me1 | `peaks/beds/H3K4me1CerebellumLate1.bed` | `peaks/beds/H3K4me1CerebellumEarly1.bed` |
| H3K4me3 | `peaks/beds/H3K4me3CerebellumLate2.bed` | `peaks/beds/H3K4me3CerebellumEarly2.bed` |
| Bivalent | `peaks/beds/Bivalent_Cerebellum_Late.bed` | `peaks/beds/Bivalent_Cerebellum_Early.bed` |
| CTCF | `peaks/CTCF.bed` (Late only) | `peaks/ctcf_motifs_mm10.bed` (motif-based) |

### Differential ChIP-seq peaks (diffbind, summit=400bp)

| Mark | Down in mutant | Up in mutant |
|------|----------------|--------------|
| H3K27me3 (Late/Adult) | `peaks/new/adult_K27me3_down.bed` | `peaks/new/adult_K27me3_up.bed` |
| H3K27me3 (Early/P12) | `peaks/new/P12_H3K27me3_down.bed` | `peaks/new/P12_H3K27me3_up.bed` |
| H2AK119ub | `peaks/new/H2AK119ub_down.bed` | `peaks/new/H2AK119ub_up.bed` |

**Note:** Summit width of 400bp may underestimate broader H3K27me3/H2AK119ub domains. Consider aggregate signal heatmaps for these marks.

### RNA-seq

| Timepoint | File | Key columns |
|-----------|------|-------------|
| Late (Adult) | `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` | `ensembl_gene_id` (gene symbols), `log2FoldChange`, `padj`, `baseMean` |
| Early (Young) | `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx` | Same schema |

### Template scripts (reusable patterns)

| Pattern | Script | Reusable for |
|---------|--------|-------------|
| GREAT-style gene association + violin plot | `tads/scripts/deg_tad_violin.R` | Sections 1e, 2b |
| CDF + density + statistical tests | `scripts/loop_distance_k27me3_filtered.R` | Sections 4b-4d |
| Generalized mark-filtered analysis | `scripts/loop_distance_mark_filtered.R` | Any ChIP mark filtering |
| Extended anchor classification | `scripts/annotate_loops_extended.R` | Sections 1b, 3c |
| Polycomb-specific shared anchor | `scripts/polycomb_shared_anchor_analysis.R` | Section 3c, 3f |
| APA heatmaps | `scripts/apa_analysis.R` | Sections 1d, 3d |
| Multi-format output | `scripts/utils/multi_format_output.R` | All new scripts |
| Permutation analysis (proper) | regioneR/regioneReloaded (Bioconductor) | Section 6d, any enrichment testing |

### Permutation Analysis Best Practices (regioneR/regioneReloaded)

- Shuffle one BED file across genome N times (≥1000 for exploratory, 10,000 for final publication)
- Measure overlap between shuffled and target BED file each iteration
- Build null distribution, compare actual overlap
- **Background restriction:** If asking about co-enrichment of euchromatin-associated markers, restrict background to euchromatin regions (same logic for heterochromatin, etc.)
- regioneReloaded vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/regioneReloaded/inst/doc/regioneReloaded.html

---

## 9. Biomodal Section 13 Follow-Up: ATAC Chromatin State + Loop Anchor Accessibility

**Context:** Section 13 (`biomodal/downstream/scripts/viz_sections/section_13_atac_chromatin_and_loops.R`) established that ATAC Up peaks are enriched at Polycomb/Repressed_Promoter regions while ATAC Down peaks concentrate at Active_Enhancer/H3K27ac regions. Loop-ATAC concordance is 39.5% overall, with Active_Enhancer-Active_Enhancer loops showing 80% concordance.

### Tasks

- [ ] **9a. Anchor-level (unpooled) ATAC overlap rates.** Show ATAC overlap separately for anchor1 vs anchor2 instead of pooling. Confirms whether the concordance signal is symmetric or driven by one anchor position.
- [ ] **9b. Test Polycomb-enriched ATAC Up peaks at DMR hypermethylated regions.** Section 13a showed ATAC Up peaks are unexpectedly enriched at Polycomb regions (12.7% vs 0.5%). Section 12 showed hypomethylated DMRs overlap ATAC Up (50.5%). Overlap the Polycomb-classified ATAC Up peaks with hypermethylated DMRs to test whether the Polycomb accessibility gain is linked to methylation changes.
- [ ] **9c. Permutation test for loop anchor ATAC enrichment vs genomic background.** Current ATAC overlap rates at loop anchors (32.6% ATAC Up, 14.0% ATAC Down) lack a null comparison. Shuffle ATAC peaks across the genome (regioneR) and compare observed vs expected overlap at loop anchors to confirm the concordance is above chance.

### Existing resources

- **Section 13 script:** `biomodal/downstream/scripts/viz_sections/section_13_atac_chromatin_and_loops.R`
- **Section 13 outputs:** `biomodal/downstream/plots/visualizations/13{a-f}_*/`, tables in `tables/atac_chromatin_*.tsv`, `tables/loop_anchor_*.tsv`, `tables/loop_atac_*.tsv`
- **ATAC peaks:** `peaks/atac_seq/ATAC_up.bed` (7,620), `peaks/atac_seq/ATAC_down.bed` (3,744)
- **Loop annotations:** `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv` (2,910 loops, pre-computed 8-category anchor types)
- **DMR data:** loaded via `_shared_config.R` (`mc_dmr`, `hmc_dmr`)
- **Permutation framework:** regioneR/regioneReloaded (see Section 6d best practices)
