# TODO: New Analyses Action Items

Post-April 2026 CTEA meeting planning document. Derived from `ctea-april-meeting-notes.md`, `docs/dixon-meeting-summary.md`, reference papers (Popay et al. Nat Genet 2026; Yu et al. Commun Biol 2026), and the Tessa Popay email (`cluster/popay.txt`).

**Last updated:** 2026-04-21

---

## 1. Superenhancer Hub Analysis

**Source:** CTEA meeting notes, Dixon meeting

**Goal:** Determine whether differential stripes overlap super-enhancer (SE) domains, and whether genes within SE hubs show correlated changes in enhancer contact frequency and H3K27ac.

### 1a. SE-DEG Proximity Analysis

- Pull DEGs and a matched set of invariant genes as comparison
- For each DEG, identify SEs within a 2 Mb window
- Quantify contact frequency (APA) between DEGs and SEs vs. DEGs and regular enhancers
- Sub-classify contact frequency changes by H3K27ac change at the enhancer

**Input files:**
- `peaks/Superenhancers_P60.bed` (1,046 SEs, P60 cerebellum)
- `peaks/Superenhancers_encode.bed` (52 SEs, ENCODE/Bing Ren)
- `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` (late DEGs)
- `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx` (early DEGs)
- `peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt` (differential K27ac)

**New work:** No SE-DEG proximity script exists yet.

### 1b. Differential Stripes x Superenhancers

- Overlap differential stripes (from stripenn pipeline) with SE beds
- Check: are strengthened stripes overlapping gained SEs? Weakened stripes overlapping lost SEs?
- Overlap both P60 and ENCODE SE sets and compare

**Input files:**
- Stripenn outputs: `stripes/stripenn/outputs/{tp}/{tp}_stripes_*.bedpe`
- SE BEDs above

### 1c. SE Hub Gene Classification

- Classify genes within SE hubs (stripes as SE hubs: enhancers modifying multiple promoters)
- Different genes proximal to same enhancer -- are they coordinately regulated?

**Meeting note:** "stripes can be super-enhancer hubs" / "enhancers modifying multiple promoters"

---

## 2. Loop Anchor CTCF Motif Analysis

**Source:** CTEA meeting notes, Dixon meeting slide 13 discussion

**Goal:** Refine loop anchor annotation by centering on CTCF/Rad21 motifs within differential anchors. Current 10 kb anchors are large; the actual motif-level signal may be more informative.

### 2a. Re-center Annotation on CTCF Motif

- For each 10 kb loop anchor, identify the CTCF motif(s) within it
- Center the chromatin state annotation around the CTCF region specifically
- Most anchors should contain a CTCF motif; flag those that don't

**Key question from meeting:** "Are changes happening at the actual motif?"

### 2b. CTCF Motif Filter for Loops

- Set a filter: at least one anchor must have a CTCF motif (or be paired with one)
- Re-run the loop type distribution (`annotate_loops_extended.R`) with this filter
- Compare filtered vs. unfiltered distributions
- Jesse noted: "to have a loop, only need 1 anchor with CTCF; still a loop, not as stable; how you get a stripe"

**Input files:**
- `peaks/ctcf_motifs_mm10.bed` (genome-wide CTCF DNA motifs)
- `peaks/CTCF.bed` (CTCF ChIP-seq peaks)
- `output/loop_annotation_extended/{early,late}/extended_characterized_loops.tsv`

**Existing script to modify:** `scripts/annotate_loops_extended.R`

### 2c. Anchor Type at the Other End

- For loops where one anchor has a CTCF motif, characterize what's at the other end
- Enhancer? Promoter? Another CTCF? This informs loop function (structural vs. regulatory)

---

## 3. DEG-Centric Enhancer Contact Analysis

**Source:** CTEA meeting notes

**Goal:** Focus specifically on differentially expressed genes and quantify their contact strength with enhancers, rather than starting from loops.

### 3a. DEG vs. Invariant Contact Strength

- Pull DEGs and invariant genes (matched for expression level/chromosome)
- Measure contact strength between just those genes and nearby enhancers
- Compare: do DEGs show more differential enhancer contact than invariants?

**Meeting note:** "maybe will clear up ABC plot" -- this is intended to produce a cleaner signal than the genome-wide ABC scatter.

### 3b. K27ac Sub-classification

- Among DEG-enhancer contacts, sub-classify by K27ac change at the enhancer
- How frequently does a DEG contact a SE vs. regular enhancer?
- Does contact frequency change correlate with K27ac change magnitude?

**Existing partial work:**
- ABC pipeline (`abc/`) already computed delta-ABC scores for 180K E-G pairs
- `abc/results/gene_level_summary.tsv` has per-gene strongest enhancer stats
- `output/deg_loop_violin/late/` has DEG-loop intersection violin plots
- This task extends those analyses with enhancer-specific (not loop-specific) framing

---

## 4. Tessa Popay Anchor vs. Span Analysis (ChromHMM)

**Source:** Dixon meeting (high priority), Popay email (received 2026-04-21)

**Goal:** Determine whether Polycomb enrichment (K27me3) is at the anchors specifically or across the loop body/span. Body enrichment supports an extrusion impediment model; anchor-specific enrichment supports a sensitivity model.

### 4a. Adapt Popay Pipeline

- Popay's code: `https://github.com/tpopay/HiC-clustering/tree/main`
- Local copy in `cluster/` directory (notebooks + tools)
- Key notebook: `cluster/grouped_loops_figures.ipynb` > "ChromHMM relationship" section
- Uses `ChromHMM OverlapEnrichment` on anchors vs. spans separately
- Need to install ChromHMM: `https://compbio.mit.edu/ChromHMM/`
- Need to install Cluster 3.0: `http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm`
- Conda env specs: `cluster/environment_mac.yml` or `cluster/environment_linux.yml`

### 4b. Generate ChromHMM Segmentation for mm10 Cerebellum

- Popay's example uses hTERT RPE-1 12-state ChromHMM
- We need equivalent for mm10 cerebellum using our histone marks
- OR: use our existing 7-category annotation from `annotate_loops_extended.R` as a proxy
- Available marks: H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent, K119ub (per timepoint)

### 4c. Run Anchor vs. Span for Differential Loops

- Input: differential loops separated by direction (up/down in mutant)
- Format as BEDPE with first 6 columns (see Popay's Google Drive example)
- Run ChromHMM OverlapEnrichment for anchors and spans separately
- Generate heatmaps using `cluster/chromHMM_heatmap.py`

### 4d. ChIP-seq Metagene Plots at Loop Anchors

- Dixon's specific suggestion: metagene-style ChIP-seq plots from loop anchors
- Plot K27me3, K27ac, K119ub signal across anchor AND flanking regions
- Look for K27me3 enrichment flanking CTCF sites in gained loops
- This is Fig 2f/2g of Popay et al. (Nat Genet 2026) adapted to our data
- Use `cluster/deepTools_pipeline.py` `bed_pileup()` function with our bigWigs

**Key question:** Is K27me3 flanking the gained CTCF-CTCF loops (extrusion impediment), or at the anchors themselves (sensitivity to repressive state)?

### 4e. K-means Clustering of Loops (Popay Style)

- Cluster differential loops by their contact strength across conditions (Cluster 3.0)
- `cluster/HiC_cluster3.ipynb` provides the full workflow: elbow plot, k-means, cluster match
- Input: loop count file (6 columns coords + balanced counts per condition)
- Use our existing count matrices from `outputs/res_*kb/` or extract with cooltools
- Classify clusters by chromatin state composition, loop size, expression of anchor genes

---

## 5. Stripe Analysis Follow-ups

**Source:** CTEA meeting notes, Dixon meeting

### 5a. Stripenn Stage 7 Visualization (run on HPC)

- Script written, pending execution
- `cd /expanse/.../stripes/stripenn && sbatch scripts/stripenn_visualizations.sb`
- Produces annotated 28-col BEDPEs, volcano plots, ChIP-seq annotation, GO/KEGG enrichment

### 5b. Stripes as Boundary Explanation

- Dixon hypothesis: "differential boundaries" may actually be stripes/flares (directional cohesin loading)
- Check: are differential boundary-proximal DEGs enriched at stripe-like features?
- If stripes explain the boundary-expression association, that's more interpretable

### 5c. Polycomb Flanking Stripe Anchors

- From meeting notes: "polycomb regions flanking anchor sites?" / "span distance has no polycomb?" / "extending along euchromatin stretch/span, but then hits polycomb K27me3 and stops"
- Test: at differential stripe anchors, is K27me3 enriched at the termini (blocking extrusion)?
- This connects to the anchor-vs-span analysis (Section 4)

### 5d. Compare Stripenn vs. Quagga Results

- Stripes called by both methods = highest confidence set
- Both pipelines complete for 250402 and 250831
- Quagga: `stripes/quagga/`, Stripenn: `stripes/stripenn/`
- Cross-reference by anchor overlap with tolerance

---

## 6. Compartment Analysis Improvements

**Source:** Dixon meeting (medium priority)

### 6a. PC1 Scatter Plot (not Volcano)

- Jesse: plot PC1_mutant vs. PC1_WT as a scatter, not a volcano
- A change of +/-1 in PC1 doesn't necessarily mean a flip
- Partial work exists: `output/compartment_analysis/compartment_pc1_scatter/`
- Verify this is what Dixon described; may need refinement

### 6b. SNIPER Subcompartment Calling

- Tool: https://github.com/ma-compbio/SNIPER (deep learning, needs ~1B reads -- confirmed sufficient)
- Jesse recommended over DCIC (which we tried with iffy results)
- Goal: look for B1<->B2 or A1<->A2 switches rather than full A<->B flips
- Validate subcompartment calls against our control histone modification data
- H3K36me3, H4K20me, K79me have different subcompartment associations

**Status:** Not started. No SNIPER files in repo.

### 6c. Correlate Differential TADs with K119ub

- Jesse's direct suggestion: which TADs gaining interaction density are at K119ub-enriched regions?
- Increased TAD ID in B compartment (heterochromatin gaining density)
- Decreased TAD ID in A compartment (euchromatin losing density)
- Check with K119ub status to test model
- `output/tad_k119ub_analysis` exists

---

## 7. Insulation Score Investigation

**Source:** Dixon meeting (high priority)

### 7a. Cooltools Pileup at Differential Boundaries

- Jesse flagged unusual insulation profiles: low flanking -> rise -> dip instead of normal high -> dip
- Current analysis used Genova; cooltools is more standard
- Popay's `cluster/cooltools_called.py` has infrastructure but points to her local paths -- adapt for our data
- Cooltools conda env spec at: `cf-pipeline-scripts/conda_envs/conda_cooltools.yml`

### 7b. Browser Validation of Boundary Calls

- Browse specific loci in genome browser to visually check what TADcompare calls "boundaries"
- Determine if sub-TAD structure is being called rather than real boundaries
- Jesse noted the unusual profiles could be stripes/flares rather than true boundary changes

---

## 8. ABC Analysis Cleanup

**Source:** Dixon meeting

### 8a. Filter Per-Enhancer Plot

- Current per-enhancer scatter has 90% of data as a blob near zero
- Filter to only sites with significantly differential K119ub
- Use `peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt` for filtering

### 8b. Separate Activity vs. Contact

- Plot activity changes and contact changes separately (not multiplied as ABC score)
- Determine if the effect is mostly activity-driven (as suggested by current data)
- Relates to slide 32 discussion: contact changes are minimal, activity changes are clear

**Existing files:**
- `abc/results/delta_abc_all_pairs.tsv` (180K E-G pairs with activity + contact columns)
- `abc/results/k119ub_enhancer_signal.tsv` (K119ub at enhancers)

---

## 9. Loop Length / Extrusion Mechanism

**Source:** Dixon meeting, Yu et al. (Commun Biol 2026)

### 9a. Polycomb at Anchor vs. Body (Extrusion Impediment Test)

- Lost loops are longer than gained (median 625 kb vs 320 kb, KS p = 5.49e-48)
- Key question: is Polycomb repression at the anchors or across the loop body?
- Body -> impediment to extrusion; anchors -> CTCF sites sensitive to repressive state
- This is the same question as Section 4d but from the loop (not stripe) perspective

### 9b. CpG Island / GC Content at Affected CTCF Anchors

- Jesse asked: are CTCF sites that lose long-range loops in CpG island regions?
- Connects to Brad Bernstein IDH mutant glioma work (repressive chromatin knocking out CTCF sites)
- Test: GC content / CpG island density at CTCF anchors of lost vs. gained loops

### 9c. Polycomb-Cohesin Extrusion Model (from Reference Papers)

Popay et al. (Nat Genet 2026) findings to contextualize with our data:
- "Mixed-dependency" chromatin loops (K-means cluster) are associated with repressive chromatin
- Loop anchors devoid of K27me3 but strong K27me3 signal flanking the peak
- STAG1 stabilizes persistent loops; STAG2 at less-dependent loops
- NIPBL-dependent genes have stripe-like features emanating from TSS

Yu et al. (Commun Biol 2026) findings:
- NIPBL KD: loss of E/P loops, gain of Polycomb domain contacts
- Larger loops more sensitive to NIPBL loss (>400 kb showing most reduction)
- Lost Ctrl-specific loops are longer; gained si-NIPBL loops are shorter
- Loss of cohesin-mediated loops and gain of PRC compartmental contacts

**Relevance:** Both papers show the same pattern we see in BAP1-KO: long active loops lost, short repressive contacts gained. Our K119ub data provides a mechanistic link (H2AK119ub elevation -> Polycomb spreading -> extrusion impediment).

---

## 10. Developmental Comparison

**Source:** Dixon meeting (medium priority)

- Compare adult mutant Hi-C to P12 wild-type Hi-C
- Does the mutant adult look more like immature neurons?
- Could support a "blocked developmental remodeling" narrative
- Jesse's remaining puzzle: BAP1 is more expressed early but phenotype is stronger late

---

## 11. Reference Papers to Read / Cite

| Paper | Relevance | Status |
|-------|-----------|--------|
| Popay et al. Nat Genet 2026 (s41588-026-02516-y) | NIPBL/cohesin loop dynamics, anchor vs span ChromHMM, K-means clustering | In repo, read |
| Yu et al. Commun Biol 2026 (s42003-026-09838-x) | NIPBL KD, Polycomb-driven loop collapse, extrusion model | In repo, read |
| Luthi et al. Nat Commun 2025 | Referenced in meeting notes, context unknown | NOT in repo, need to retrieve |
| Brad Bernstein IDH mutant glioma | Repressive chromatin knocks out CTCF binding | NOT in repo, need to retrieve |

---

## Priority Ordering

### Immediate (for paper)

1. **Popay anchor-vs-span analysis** (Section 4) -- Tessa just shared code today; highest scientific value for distinguishing extrusion impediment vs. sensitivity model
2. **Stripenn Stage 7** (Section 5a) -- just needs `sbatch`, script is written (already run)
3. **CTCF motif re-centering** (Section 2) -- refines existing annotation, meeting consensus
4. **SE-DEG proximity** (Section 1a) -- direct meeting request, straightforward analysis
5. **DEG-centric enhancer contacts** (Section 3) -- fills gap in current ABC analysis

### Soon (strengthens story)

6. **ABC cleanup** (Section 8) -- quick filter + replot
7. **Loop body Polycomb analysis** (Section 9a) -- mechanistic insight
8. **Stripe-boundary connection** (Section 5b) -- interpretive
9. **Insulation score with cooltools** (Section 7) -- validates/corrects TADcompare results

### Medium-term

10. **SNIPER subcompartments** (Section 6b) -- HPC computation needed
11. **Developmental comparison** (Section 10) -- interpretive framing
12. **CpG/GC content at CTCF anchors** (Section 9b) -- mechanistic detail
13. **Stripenn vs. Quagga comparison** (Section 5d) -- confidence set

---

## Infrastructure Notes

### Popay Pipeline Setup

```bash
# Create conda env from Popay's specs
conda env create -f cluster/environment_mac.yml   # or environment_linux.yml for HPC

# Install Cluster 3.0
# http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm

# Install ChromHMM
# https://compbio.mit.edu/ChromHMM/

# Adapt cooltools_called.py paths from Popay's local to our HPC data
```

### Data Already Available

| Data | Location | Status |
|------|----------|--------|
| Differential loops (early + late) | `outputs/*/merged_loops/non_redundant_loops.tsv` | Ready |
| Loop annotations (8-category) | `output/loop_annotation_extended/` | Ready |
| Stripenn differential stripes | `stripes/stripenn/outputs/` | Ready (simple BEDPEs) |
| DEG lists (early + late) | `tads/*.xlsx` | Ready |
| ChIP-seq peaks (5 marks x 2 tp) | `peaks/beds/` | Ready |
| DiffBind results (4 marks) | `peaks/diffbind/` | Ready |
| Superenhancer BEDs | `peaks/Superenhancers_*.bed` | Ready |
| CTCF motifs (mm10) | `peaks/ctcf_motifs_mm10.bed` | Ready |
| ABC E-G pairs | `abc/results/delta_abc_all_pairs.tsv` | Ready |
| K119ub enhancer signal | `abc/results/k119ub_enhancer_signal.tsv` | Ready |
| Hi-C .hic files | HPC: `/expanse/.../juicerpre/*.hic` | On HPC |
| mcool files (stripenn) | HPC: `/expanse/.../stripes/stripenn/data/cool/` | On HPC |
| Compartment PC1 values | `tads/tad-pc-analysis/` | Ready |
