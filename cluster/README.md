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

Gained loops sit within expanding Polycomb domains (both anchors and span are heterochromatic). Lost loops have Polycomb gain specifically at anchor sites while the span remains euchromatic.

### 9-Mark ChromHMM Expansion (Phase 2b)

Adding H2AK119ub, ATAC, H3K9ac, and H3K9me3 to the model (18 states) resolves K119ub-specific biology:

| State | clust5 anchor | clust5 span | clust6 anchor | clust6 span |
|-------|:---:|:---:|:---:|:---:|
| **Polycomb_K119ub** (K27me3 + K119ub) | **8.20x** | 2.57x | **4.14x** | 0.97x |
| **Repressed_Enhancer_K119ub** (K4me1 + K119ub + K27me3) | **9.02x** | 2.16x | 2.82x | 0.96x |
| **Active_Enhancer_K119ub** (active marks + K119ub) | 2.10x | 0.56x | **5.15x** | 1.06x |

Repressed_Enhancer_K119ub (9.02x) is the highest enrichment in the dataset -- fully PRC1+PRC2-silenced enhancers at gained-loop anchors. Active_Enhancer_K119ub (5.15x at clust6) represents the transitional state where BAP1 has failed to remove K119ub from active enhancers at lost-loop anchors.

The oriented metagene analysis (Phase 8) further refines this: K27me3 is asymmetrically interior-enriched at gained-loop anchors (Ext/Int = 0.91, p = 0.004) but symmetric at lost-loop anchors (Ext/Int = 1.02, p = 0.52). Comprehensive asymmetry (Phase 11) confirms this with orthogonal signals: gained-loop anchors sit at TAD boundaries (insulation p = 3.5e-54) at the edge of B-compartment domains (PC1 p = 6.3e-44), while lost-loop anchors reside in A-compartment euchromatin with no directional Polycomb spreading.

---

## Quick Start

All runner scripts accept an optional `CLUSTER_CONF` env var pointing to a config file. For the BAP1-KO late timepoint:

```bash
cd cluster
export CLUSTER_CONF=scripts/config/late.conf   # optional — sets all paths/params

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

# Phase 9: Clust6 subgroup asymmetry (~2-4 min)
bash scripts/run_clust6_subgroups.sh

# Phase 10: Histone anchors metagene profiles (~1-2 min)
bash scripts/run_phase5.sh   # or invoke Python directly

# Phase 11: Comprehensive asymmetry (HPC only — needs mcools)
sbatch scripts/12_comprehensive_asymmetry.sb

# Phase 2b: 9-mark ChromHMM expansion (HPC for segmentation, Mac for downstream)
sbatch scripts/03b_chromhmm_9mark.sb intersect
sbatch scripts/03b_chromhmm_9mark.sb union
# MANUAL: write {15,18}state_rename_cerebellum.txt after inspecting emissions
bash scripts/run_phase4_9mark.sh intersect 18
```

### Prerequisites

- **Conda env:** `cluster` (Python 3.8). Create from `modules/cluster.yml`.
- **Cluster 3.0:** k-means engine, must be on PATH as `cluster`. Install from [bonsai.hgc.jp](http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm).
- **ChromHMM v1.27:** Installed at `ChromHMM/ChromHMM.jar`. Wrapper script `ChromHMM/chromhmm` must be on PATH.
- **System R 4.5.2** with Bioconductor packages (Phase 1 only): `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`.
- **BigWig files** (4 marks x ctrl/mut). Default: `/Users/zakiralibhai/sdsc/bigwigs/`. Set `CLUSTER_BIGWIG_DIR` in your `.conf` to override.

---

## Configuration

The pipeline is parameterized via shell `.conf` files sourced by runner scripts. Config files live in `scripts/config/`:

| File | Purpose |
|------|---------|
| `late.conf` | BAP1-KO late/adult (250402) -- all values populated |
| `early.conf` | BAP1-KO early/P12 (250831) -- partially populated |
| `template.conf` | Documented template for new projects |

Config files define env vars read by Python scripts via `os.environ.get(VAR, default)`. All variables have BAP1-KO defaults, so existing runs produce identical output without any config.

Shared Python utilities live in `scripts/utils/pipeline_config.py` -- BigWig dict construction, cluster ordering, biological labels, deepTools header parsing, condition column names, and other constants previously duplicated across scripts.

### Using for a New Project

```bash
# 1. Copy the template
cp scripts/config/template.conf scripts/config/myproject.conf

# 2. Fill in required fields:
#    - CLUSTER_TIMEPOINT_LABEL, CLUSTER_TIMEPOINT_ID, CLUSTER_OUT_DIR, CLUSTER_CELL_NAME
#    - CLUSTER_CTRL_REPS, CLUSTER_MUT_REPS, CLUSTER_COND1_COL, CLUSTER_COND2_COL
#    - CLUSTER_MERGED_BEDPE, CLUSTER_COUNTS_TEMPLATE, CLUSTER_EDGER_TEMPLATE
#    - CLUSTER_PEAK_* (ChIP-seq peak BEDs for ChromHMM)
#    - CLUSTER_BIGWIG_DIR (and optionally CLUSTER_BW_* per-file overrides)
#    - CLUSTER_PYTHON, CLUSTER_BIN (machine-specific binary paths)

# 3. Run the pipeline
export CLUSTER_CONF=scripts/config/myproject.conf
bash scripts/run_phase1.sh
bash scripts/run_phase2.sh
# ... inspect emissions, create rename file ...
bash scripts/run_phase3.sh 6 0.01 0.333 3.0
bash scripts/run_phase4.sh all

# 4. After inspecting clustering results, add biological annotation:
#    CLUSTER_BIO_ORDER="clust4,clust2,clust1,clust3,clust5,clust6"
#    CLUSTER_BIO_NAME_clust5="Strong gain"
#    CLUSTER_BIO_NAME_clust6="Strong loss"
#    Then re-run Phase 7 for publication-quality summary figures.
```

The pipeline supports different organisms, condition counts, mark sets, replicate numbers, and machines via the `.conf` system. The only architectural constraint is exactly 2 conditions (for the mut/ctrl ratio used in k-means). Popay's 3-condition comparison (DMSO/4h/24h) would require a separate adaptation.

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
Phase 9: 10_clust6_subgroup_asymmetry.py (clust6 short/long split + asymmetry)
Phase 10: 11_histone_anchors_metagene.py (clean profile figure from Phase 5 matrix)
Phase 11: 12_comprehensive_asymmetry.py (H2AK119ub, H3K27ac, PC1, insulation --
           HPC only, computes PC1/insulation BigWigs from mcools)
Phase 2b: 03b_chromhmm_9mark_segmentation.sh (9-mark ChromHMM, 15+18 states --
           HPC; adds H2AK119ub, ATAC, H3K9ac, H3K9me3)
         + run_phase4_9mark.sh (Phase 4.4/4.5 rerun with 9-mark env vars)
```

Phases 2 and 3 are independent (both depend on Phase 1). Phase 4 requires both 2 and 3. Phases 5-11 depend on Phase 3's clustering output. Phase 11 additionally requires mcools on HPC.

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
      clust6_subgroups/                 # Phase 9: clust6 short/long asymmetry
      comprehensive_asymmetry/          # Phase 11: H2AK119ub, H3K27ac, PC1, insulation
        bigwigs/                        #   PC1 + insulation BigWigs (computed from mcools)
    summary_figures/
      dashboard/                        # 6-panel cluster summary
      mechanism/                        # clust5 vs clust6 mechanism comparison
      heatmap/                          # z-scored feature heatmap + feature_values.tsv
  chromHMM_9mark_intersect/              # Phase 2b: 9-mark model (18-state selected)
    learned_model_18/                   #   Segmentation + emissions (18 states x 9 marks)
    18state_rename_cerebellum.txt        #   E1-E18 -> biological name mapping
    {anchor,span}_18.txt                #   Fold-enrichment tables (18 states x 6 clusters)
    {anchor,span}_18.{png,pdf,svg,jpg}  #   Anchor-vs-span heatmaps
  chromHMM_9mark_union/                 # Phase 2b: union consensus (for comparison)
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
| `10_clust6_subgroup_asymmetry.py` | 9 | Clust6 short/long split + oriented asymmetry | ~2-4 min |
| `11_histone_anchors_metagene.py` | 10 | Clean profile figure from Phase 5 matrix | ~1-2 min |
| `12_comprehensive_asymmetry.py` | 11 | H2AK119ub, H3K27ac, PC1, insulation asymmetry (HPC) | ~90 min |
| `03b_chromhmm_9mark_segmentation.sh` | 2b | 9-mark ChromHMM (BinarizeBed + LearnModel k=15,18) | ~20 min |
| `03b_chromhmm_9mark.sb` | 2b | SLURM wrapper (creates consensus BEDs + runs segmentation) | ~20 min |
| `run_phase4_9mark.sh` | 2b | Phase 4.4+4.5 with 9-mark env var overrides | ~15 sec |

Runner scripts in `scripts/run_*.sh` handle PATH setup, environment activation, log capture, and optional `.conf` sourcing via `CLUSTER_CONF`.

### Shared Utilities (`scripts/utils/`)

| Utility | Purpose |
|---------|---------|
| `multi_format_output.py` | `multi_format_savefig()` context manager -- auto-emits SVG/PDF/JPG alongside PNG |
| `pipeline_config.py` | Shared configuration functions -- BigWig dict, cluster ordering, biological labels, deepTools header parser, condition columns, constants |

### Config Files (`scripts/config/`)

| Config | Purpose |
|--------|---------|
| `late.conf` | BAP1-KO late timepoint -- all variables populated |
| `early.conf` | BAP1-KO early timepoint -- partially populated |
| `template.conf` | Documented template for new projects |

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

**Important:** `conda run -n cluster python3` does NOT reliably activate the environment. Use the absolute interpreter path or set `CLUSTER_PYTHON` in your `.conf` file. All binary paths (`CLUSTER_PYTHON`, `CLUSTER_RSCRIPT`, `CLUSTER_BIN`) are configurable via `.conf` files and fall back to `$(command -v ...)` if unset.
