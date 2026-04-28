# Popay-Style Hi-C Loop Clustering for BAP1-KO Cerebellum

Adaptation of the [Popay et al. (Nat Genet 2026)](https://www.nature.com/articles/s41588-026-02516-y) Hi-C loop clustering pipeline for BAP1-KO mouse cerebellum. Answers the central mechanistic question: **is Polycomb (H3K27me3) enrichment at the anchors of differential loops, or across the loop body/span?**

The anchor-vs-span distinction separates two models of how BAP1 loss disrupts 3D genome organization:
- **Anchor enrichment** = sensitivity model (CTCF sites directly disrupted by repressive chromatin)
- **Span enrichment** = extrusion impediment model (Polycomb spreading blocks cohesin loop extrusion)

## Key Result

Both mechanisms operate simultaneously in different loop populations:

| Population | Polycomb Anchor | Polycomb Span | Mechanism |
|------------|----------------:|-------------:|-----------|
| **Gained loops** (clust5, n=667, 97% up) | **6.59x** | **3.03x** | Polycomb-domain compaction |
| **Lost loops** (clust6, n=2,359, 78% down) | **2.09x** | 0.94x | Anchor disruption |

Gained loops sit within expanding Polycomb domains (both anchors and span are heterochromatic). Lost loops have Polycomb gain specifically at anchor sites while the span remains euchromatic. The oriented metagene analysis (Phase 8) further refines this: K27me3 is asymmetrically interior-enriched at gained-loop anchors (Ext/Int = 0.91, p = 0.004) but symmetric at lost-loop anchors (Ext/Int = 1.02, p = 0.52).

---

## Quick Start

```bash
# All phases run from cluster/ (runners cd themselves)
cd cluster

# Phase 1: Data prep (~1 min)
bash scripts/run_phase1.sh

# Phase 2: ChromHMM segmentation (~2 min)
bash scripts/run_phase2.sh
# MANUAL: edit outputs/bap1_late/chromHMM/12state_rename_cerebellum.txt

# Phase 3: K-means clustering (~2 min)
bash scripts/run_phase3.sh 6 0.01 0.333 3.0

# Phase 4: All downstream analyses (~6 min)
bash scripts/run_phase4.sh all

# Phase 5: deepTools anchor metagene (~90 min)
bash scripts/run_phase5.sh

# Phase 6: Cooltools pileup (HPC only — see docs)
# Phase 7: Summary figures (~1 min)
bash scripts/run_summary.sh

# Phase 8: Oriented anchor metagene (~1.7 h)
bash scripts/run_oriented_metagene.sh
```

### Prerequisites

- **Conda env:** `cluster` (Python 3.8). Create from `modules/cluster.yml`.
- **Cluster 3.0:** k-means engine, must be on PATH as `cluster`. Install from [bonsai.hgc.jp](http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm).
- **ChromHMM v1.27:** Installed at `ChromHMM/ChromHMM.jar`. Wrapper script `ChromHMM/chromhmm` must be on PATH.
- **System R 4.5.2** with Bioconductor packages (Phase 1 only): `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`.
- **BigWig files** at `/Users/zakiralibhai/sdsc/bigwigs/` (H3K27ac, H3K27me3, H2AK119ub, H3K27me1 -- ctrl + mut).

---

## Pipeline Architecture

```
Phase 0: Bug fixes in inherited Popay modules (DONE)
Phase 1: 01_build_loop_count_file.py  ->  data/late_merged_loop_counts.txt (39,344 loops)
         02_build_mm10_gene_annotation.R  ->  data/mm10_knownGene_pp.bed
              |---------------------|
Phase 2: ChromHMM               Phase 3: K-means
  03_chromhmm_segmentation.sh     04_clustering.py
  (BinarizeBed -> LearnModel)     (elbow -> Cluster 3.0 -> sort)
              |----------|--------|
Phase 4: 05_grouped_analyses.py (8 sub-analyses)
  4.1 loop size  |  4.2 loop classification  |  4.3 anchor ChIP signal
  4.4 ChromHMM anchor vs span (KEY RESULT)   |  4.5 ChromHMM proportions
  4.6 gene annotation  |  4.7 DiffBind       |  4.8 cluster x differential
Phase 5: 06_deeptools_metagene.py (anchor +/-5kb per cluster x 8 BigWigs)
Phase 6: 07_cooltools_pileup.py (off-diagonal Hi-C pileup -- HPC only)
Phase 7: 08_summary_figures.py (3 composite lab-meeting figures)
Phase 8: 09_oriented_anchor_metagene.py (strand-aware anchor metagene)
       + quantify_orientation_asymmetry.py (Wilcoxon ext/int -> TSV)
       + visualize_orientation_asymmetry.py (K27me3 dual-panel figure)
```

Phases 2 and 3 are independent (both depend on Phase 1). Phase 4 requires both 2 and 3. Phases 5-8 depend on Phase 3's clustering output.

---

## Scope

**Timepoint:** Late/adult (250402) only. The early timepoint (P12/250831) has too few differential loops (~200) for meaningful clustering.

**Input:** 39,344 nonredundant merged loops across 3 resolutions (5kb: 7,901; 10kb: 14,553; 25kb: 16,890). Replicates averaged to 2 columns (ctrl_merge, mut_merge).

**Organism:** Mouse (mm10). **Conditions:** BAP1-KO mutant vs wildtype control, n=3 per condition.

---

## Canonical Clustering (k=6, v2)

Loops clustered by normalized contact-count ratio (mut/ctrl). Biological ordering for figures: loss -> unchanged -> gain.

| Cluster | n | Direction | Median logFC | Key biology |
|---------|---|-----------|-------------|-------------|
| clust6 | 2,359 | 78% down | -0.31 | Strong loss -- anchor disruption (Polycomb anchor 2.09x / span 0.94x) |
| clust3 | 8,738 | 21% down | -0.14 | Moderate loss |
| clust1 | 12,298 | 100% unchanged | -0.03 | High-signal baseline |
| clust2 | 10,970 | 92% unchanged | +0.09 | Moderate baseline |
| clust4 | 3,916 | 70% up | +0.25 | Moderate gain |
| clust5 | 667 | 97% up | +0.50 | Strong gain -- Polycomb domain compaction (anchor 6.59x / span 3.03x) |

The clustering captures 99.3% of edgeR's down_in_mutant calls and 99.5% of up_in_mutant calls. chi2 = 38,987 for cluster x differential status (p ~ 0).

---

## ChromHMM 12-State Model

Learned from 5 wildtype ChIP-seq marks on 21 standard chromosomes. 335,569 merged segments, converged in 32 iterations.

| ID | State | Marks ON | Genome % |
|----|-------|----------|----------|
| E8 | Active_Promoter | K27ac + K4me3 | 0.35% |
| E5 | Active_Promoter_Flank | K27ac + K4me3 + K4me1 | 0.20% |
| E6 | Poised_Promoter | K4me3 + K4me1 | 0.08% |
| E7 | Weak_Promoter | K4me3 | 0.11% |
| E9 | Strong_Enhancer | K27ac | 0.36% |
| E4 | Active_Enhancer | K27ac + K4me1 | 0.77% |
| E3 | Poised_Enhancer | K4me1 | 3.67% |
| E12 | Bivalent_Enhancer | K4me1 + K27me3 | 0.31% |
| E11 | Polycomb | K27me3 | 1.81% |
| E2 | CTCF_Boundary | K4me1 + CTCF | 0.09% |
| E10 | Insulator | CTCF | 0.30% |
| E1 | Quiescent | none | 91.97% |

---

## Key Adaptations from Popay

| Aspect | Popay (hTERT RPE-1) | This project (BAP1-KO cerebellum) |
|--------|---------------------|-----------------------------------|
| Genome | hg38 | mm10 |
| Conditions | DMSO, dTAG 4h, dTAG 24h | ctrl_merge, mut_merge |
| Perturbation | Acute NIPBL depletion | Genetic BAP1 KO (developmental) |
| Loops | RAD21 ChIP-based, 16,860 | HiCCUPS + mariner, 39,344 merged |
| ChromHMM | 12-state from RPE-1 marks | 12-state from 5 cerebellum marks |
| ChIP normalization | Normalized to RAD21 | Raw RPKM (no RAD21 data) |
| Loop classification | CTCF + RAD21 = structural | CTCF alone = structural |
| DiffBind marks | RAD21 only | K27ac, K27me3, K119ub |

---

## Output Structure

```
outputs/bap1_late/
  cluster3/
    elbow_plot/                         # Elbow plot for k selection
    k-6/data/combined-clusters.txt      # Canonical cluster assignments (38,948 x 9)
    k-6/figures/{heatmap,lineplot,boxplot,stripplot}/
  chromHMM/
    learned_model/                      # ChromHMM segmentation + emissions
    12state_rename_cerebellum.txt        # E1-E12 -> biological name mapping
    {anchor,span}.txt                   # KEY: fold-enrichment tables (12 states x 6 clusters)
    {anchor,span}.{png,pdf,svg,jpg}     # KEY: anchor-vs-span heatmaps
    {anchor,span}_input/                # Per-cluster BEDs for OverlapEnrichment
  figures/
    loop_size/                          # Per-cluster size distributions
    loop_classification/                # Structural/CRE/mixed/unclassified
    chromHMM_anchor/                    # ChromHMM proportions stacked bar
    cluster_differential_status/        # Cluster x edgeR direction crosstab
    ChIP_intersect/                     # DiffBind overlap + anchor ChIP signal
    annotation/                         # Per-cluster gene lists
    deeptools/
      histone_anchors/                  # Phase 5: 8-BigWig x 6-cluster metagene
      oriented_anchors/                 # Phase 8: strand-aware metagene + asymmetry
    summary_figures/
      dashboard/                        # 6-panel cluster summary
      mechanism/                        # clust5 vs clust6 mechanism comparison
      heatmap/                          # z-scored feature heatmap + feature_values.tsv
  cooltools/
    obs_exp_contacts/                   # Phase 6: Hi-C pileup (ctrl + mut)
```

All figures are saved in 4 formats (PNG + PDF + SVG + JPG) via the `multi_format_savefig()` context manager.

---

## Scripts

| Script | Phase | Purpose | Runtime |
|--------|-------|---------|---------|
| `01_build_loop_count_file.py` | 1 | Merged loops -> Popay-format count file + metadata sidecar | ~7 sec |
| `02_build_mm10_gene_annotation.R` | 1 | mm10 TSS +/- 750 bp gene BED (24,515 genes) | ~5 sec |
| `03_chromhmm_segmentation.sh` | 2 | Stage BEDs + BinarizeBed + LearnModel(k=12) | ~93 sec |
| `04_clustering.py` | 3 | Elbow + Cluster 3.0 k-means + 4 visualizations | ~30 sec |
| `05_grouped_analyses.py` | 4 | 8 downstream sub-analyses orchestrator | ~6 min |
| `06_deeptools_metagene.py` | 5 | Per-cluster anchor BEDs + deepTools bed_pileup | ~96 min |
| `07_cooltools_pileup.py` | 6 | Off-diagonal Hi-C pileup (HPC only) | ~35 min |
| `08_summary_figures.py` | 7 | 3 composite lab-meeting figures | ~12 sec |
| `09_oriented_anchor_metagene.py` | 8 | Strand-aware anchor metagene | ~1.7 h |
| `quantify_orientation_asymmetry.py` | 8 | Ext/Int Wilcoxon signed-rank -> TSV | <1 min |
| `visualize_orientation_asymmetry.py` | 8 | K27me3 dual-panel figure | <1 min |

Runner scripts in `scripts/run_*.sh` handle PATH setup, environment activation, and log capture.

---

## Data Dependencies

### In-Repo (`data/`)

| File | Rows x Cols | Source |
|------|-------------|--------|
| `late_merged_loop_counts.txt` | 39,344 x 8 | Phase 1 from merged BEDPE + per-resolution count matrices |
| `late_merged_loop_metadata.tsv` | 39,344 x 16 | Phase 1 edgeR stats sidecar (logFC, FDR, direction, resolution) |
| `mm10_knownGene_pp.bed` | 24,515 x 7 | Phase 1 from TxDb.Mmusculus.UCSC.mm10.knownGene |

### External (Not in Repo)

| Data | Location | Notes |
|------|----------|-------|
| BigWigs (4 marks x ctrl/mut) | `/Users/zakiralibhai/sdsc/bigwigs/` | H3K27ac, H3K27me3, H2AK119ub, H3K27me1. Do NOT use `peaks/bigwigs/macs2.narrow.aug18.dedup/` (0-byte mut files). |
| ChIP peak BEDs | `peaks/CTCF.bed`, `peaks/beds/*.bed` | H3K27ac (15,105), H3K27me3 (15,809), H3K4me1 (113,781), H3K4me3 (6,581), CTCF (32,487) |
| DiffBind results | `peaks/diffbind/` | K27ac, K27me3, K119ub differential binding |
| mcool files (HPC) | `/expanse/.../cool/250402/` | ctrl_merged.mcool, mut_merged.mcool -- Phase 6 only |
| Blacklist | `tads/mm10-blacklist.v2.bed` | mm10 ENCODE blacklist v2 |

---

## Provenance

The original Popay pipeline was cloned from [github.com/tpopay/HiC-clustering](https://github.com/tpopay/HiC-clustering) and shared directly by Tessa Popay (postdoc, Dixon Lab, Salk Institute). Eight Python modules in `modules/` are the adapted Popay library with bug fixes applied in Phase 0 (hardcoded paths, deprecated matplotlib styles, hg38 genome defaults, FPKM crash guard).

### Collaboration Context

Jesse Dixon (Salk Institute) suggested this analysis during the 2026-04-10 meeting. His postdoc Tessa Popay performed the equivalent analysis for NIPBL depletion in hTERT RPE-1 cells (Popay et al., Nat Genet 2026, Fig 2f). Tessa shared her code and environment specs via email (2026-04-21) and offered ongoing help. She did her PhD on HCF-1, which associates with BAP1 in the PR-DUB complex.

---

## Documentation

| File | Contents |
|------|----------|
| `docs/CONTEXT-cluster.md` | Full biological context, meeting notes, Popay paper summary, data inventory |
| `docs/PLAN-p1.md` | Phases 0-3 plan with corrections and verification |
| `docs/PLAN-p2.md` | Phases 4-6 plan with corrections and verification |
| `docs/PLAN-p3.md` | Phases 7-8 plan with corrections and verification |
| `docs/RESULTS-cluster.md` | Statistical outputs, biological interpretation, open questions |
| `docs/phase{1..5}.txt` | Phase run logs (stdout/stderr capture) |
| `docs/phase3_v2.txt` | Canonical v2 clustering run log |
| `docs/phase8_summary.txt` | Summary figures run log |
| `docs/oriented_metagene.txt` | Phase 8 step 1 run log |
| `CLAUDE.md` | AI assistant context (environment, gotchas, module API) |

---

## Environment Details

| Component | Version / Path |
|-----------|---------------|
| Python | 3.8.18 (`/opt/homebrew/anaconda3/envs/cluster/bin/python3`) |
| pandas | 1.5.3 |
| matplotlib | 3.7.5 |
| bioframe | 0.6.1 |
| cooltools | 0.6.1 |
| deeptools | 3.5.5 |
| Cluster 3.0 | `/usr/local/bin/cluster` |
| ChromHMM | v1.27 (`ChromHMM/ChromHMM.jar`) |
| R | 4.5.2 (system, Phase 1 only) |
| bedtools | 2.31.1 (`/opt/homebrew/bin/bedtools`) |
| Java | 25 (for ChromHMM) |

**Important:** `conda run -n cluster python3` does NOT reliably activate the environment. Always use the absolute interpreter path.
