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
- [~] **3c. Polycomb-specific shared anchor analysis.** Among shared anchors (task 1a), filter to those classified as "Polycomb" or "Repressed_Promoter" in extended annotation. Show that these are the anchors driving the switching.
- [ ] **3d. APA heatmaps for Polycomb-anchored loops.** Aggregate contact signal at Polycomb-classified loop anchors, split by lost vs gained.
- [ ] **3e. Gene body heatmap for Polycomb-associated loops.** For genes at Polycomb anchors that switch from long to short loops, show a heatmap of ChIP-seq signal (H3K27me3) across the gene body.

### Existing resources

- **K27me3 filtered analysis:** `scripts/loop_distance_k27me3_filtered.R` -> `output/loops_k27me3_filtered/{early,late}/`
- **K27me3 global analysis:** `scripts/loop_distance_k27me3_global.R` -> `output/loops_k27me3_global/{early,late}/`
- **ChIP distance analysis:** `scripts/chip_distance_analysis.R` -> `output/chip_distance/{early,late}/`
- **H3K27me3 peak files:** `peaks/beds/H3K27me3CerebellumLate1.bed`, `peaks/beds/H3K27me3CerebellumEarly1.bed`
- **Bivalent peak files:** `peaks/beds/Bivalent_Cerebellum_Late.bed`, `peaks/beds/Bivalent_Cerebellum_Early.bed`
- **Extended annotation with Polycomb category:** `scripts/annotate_loops_extended.R` (classifies anchors as Polycomb when H3K27me3+ and >2kb from TSS)
- **APA pipeline:** `scripts/apa_analysis.R`

---

## 4. Anchor-Specific ChIP-seq Subsetting

**Meeting notes:** "do cdf density - for each histone modification", "k27ac on one site vs both sites (super enhancer)", "+ k4me3 - enhancer+promoter loop"

**Core question:** How does the distance distribution differ when filtering loops by specific ChIP-seq mark combinations at anchors?

### Tasks

- [x] **4a. CDF/density by H3K27me3.** One anchor vs both anchors. -> `scripts/loop_distance_k27me3_filtered.R`
- [~] **4b. CDF/density by H3K27ac.** One anchor (enhancer-only loop) vs both anchors (super-enhancer loop). Existing `chip_distance_analysis.R` does trends but not separate CDF/density plots per mark subset.
- [ ] **4c. CDF/density for H3K27ac + H3K4me3 (enhancer-promoter loops).** Filter to loops where one anchor is H3K27ac+ (enhancer) and the other is H3K4me3+ (promoter). Compare lost vs gained distance distributions.
- [ ] **4d. CDF/density for each individual histone modification separately.** Repeat the K27me3-filtered analysis pattern for H3K27ac, H3K4me1, H3K4me3, Bivalent marks.

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

- [ ] **5a. CTCF overlap analysis at lost vs gained loop anchors.** Are lost loops more likely to have CTCF at their anchors? Are gained loops CTCF-depleted?
- [ ] **5b. CDF/density for CTCF-anchored loops.** Filter to loops with CTCF at anchors, compare lost vs gained distance distributions.
- [ ] **5c. Cross-reference with stripe analysis.** Do regions with lost CTCF-anchored loops show stripe defects? Stripe pipeline already exists.
- [ ] **5d. Obtain RAD21 ChIP-seq data.** No RAD21 data currently in repo. Need to request assay or find public data.
- [ ] **5e. (If RAD21 data obtained) Overlap RAD21 with loop anchors and TAD boundaries.**

### Existing resources

- **CTCF ChIP-seq peaks:** `peaks/CTCF.bed` (32,487 peaks, Late timepoint)
- **CTCF motif predictions:** `peaks/ctcf_motifs_mm10.bed` (genome-wide motif scan)
- **Extended annotation already includes CTCF:** `scripts/annotate_loops_extended.R` has CTCF_Site category (8-category version in `peaks/loop_annotation_extended/`)
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
- [ ] **6c. Cross-reference differential loops with differential TAD boundaries.** Are differential loops preferentially located near differential TAD boundaries? Overlap analysis between `characterized_loops.tsv` anchors and `tads/results/{tp}/final/tadcompare_final_filtered.tsv` boundaries.

### Existing resources

- **TAD results:** `tads/results/{early,late}/final/tadcompare_final_filtered.tsv` (differential boundaries with Enriched_In, Type, Gap_Score)
- **TAD visualization:** `tads/scripts/tad_visualizations.R` (40+ plots across 9 subdirectories)
- **Loop data:** `25042-late_outputs/merged_loops/characterized_loops.tsv`, `250831-early_outputs/merged_loops/characterized_loops.tsv`

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
- H2AK119ub ChIP-seq data is NOT currently available in the repo
- Meeting hypothesis: "ub is buffer to stop k27ac contact" - once ubiquitination threshold is reached, contacts form

---

## 8. Ubiquitination (H2AK119ub) Integration

**Meeting notes:** "contact k27ac and ub - mediated by which?", "hypothesis: ub is buffer to stop k27ac contact", "once reach threshold - now forms contact", "correlate - enhancer gene linkage level", "how much ubiquitinated histone", "2nd: tie change in delta contacts to ubiquitinated histone"

### Tasks

- [ ] **8a. Obtain H2AK119ub ChIP-seq data.** Not currently in repo.
- [ ] **8b. (Once data obtained) Overlap H2AK119ub with loop anchors.** Classify anchors by ubiquitination status.
- [ ] **8c. Correlate delta E-P contacts with ubiquitination levels.** Test hypothesis that ubiquitination buffers H3K27ac-mediated contact formation.

### Existing resources

- **H2AK119ub data:** NOT AVAILABLE - need to generate or obtain
- **Framework for ChIP overlap:** `scripts/annotate_loops_extended.R` (can add H2AK119ub as additional mark)

---

## Summary: Priority Order

Based on meeting notes emphasis and data availability:

| Priority | Analysis | Data Available? | Script Exists? |
|----------|----------|----------------|----------------|
| 1 | Shared anchor / loop switching (Section 1) | Yes | **COMPLETE** (`scripts/shared_anchor_analysis.R`, `scripts/apa_shared_anchors.R`) |
| 2 | RNA-seq integration with loops (Section 2) | Yes | No (TAD version exists as template) |
| 3 | Polycomb loop story completion (Section 3) | Yes | Partially |
| 4 | Per-mark CDF/density subsetting (Section 4) | Yes | Partially (K27me3 done, others not) |
| 5 | ABC model / E-P linkage (Section 7) | Yes (except H2AK119ub) | No |
| 6 | CTCF analysis (Section 5) | Yes (CTCF peaks), No (RAD21) | No |
| 7 | TAD-loop cross-reference (Section 6) | Yes | No |
| 8 | H2AK119ub integration (Section 8) | No | No |

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
| Extended anchor classification | `scripts/annotate_loops_extended.R` | Sections 1b, 3c |
| APA heatmaps | `scripts/apa_analysis.R` | Sections 1d, 3d |
| Multi-format output | `scripts/utils/multi_format_output.R` | All new scripts |
