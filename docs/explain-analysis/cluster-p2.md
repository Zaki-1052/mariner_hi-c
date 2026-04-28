**Phase 2 complete** in 93 seconds (LearnModel converged in 32 iterations on 21 standard chroms / 13.6M bins).

**Deliverables:**
- `cluster/scripts/03_chromhmm_segmentation.sh` + `run_phase2.sh` (executable)
- `cluster/bap1_late/chromHMM/learned_model/cerebellum_late_12_segments.bed` — 335,569 segments
- `cluster/bap1_late/chromHMM/learned_model/emissions_12.txt`
- `cluster/bap1_late/chromHMM/12state_rename_cerebellum.txt` ← key Phase 4.4 input
- `cluster/phase2.txt` (full run log)

**12-state vocabulary** (ordered active→repressive→structural→quiescent):

```
E8  Active_Promoter        K27ac+K4me3        canonical TssA, 45× TSS-enriched
E5  Active_Promoter_Flank  K27ac+K4me3+K4me1  TssAFlnk, 19× TSS
E6  Poised_Promoter        K4me3+K4me1
E7  Weak_Promoter          K4me3 + slight K27me3   most TSS-specific (52×)
E9  Strong_Enhancer        K27ac alone
E4  Active_Enhancer        K27ac+K4me1
E3  Poised_Enhancer        K4me1 alone
E12 Bivalent_Enhancer      K4me1+K27me3       ← mechanistic key for BAP1 paper
E11 Polycomb               K27me3 alone
E2  CTCF_Boundary          K4me1+CTCF
E10 Insulator              CTCF alone
E1  Quiescent              92% of genome, depleted at loop anchors (0.87×)
```

**Smoke test (OverlapEnrichment of all 39,344 anchor1 BEDs against the rename file)**: Active_Promoter most enriched at loop anchors (**5.31×**), Quiescent depleted (0.87×), CTCF/promoter/enhancer all show coherent enrichment. Phase 4.4 plumbing is validated.

**Mechanistic note for the paper**: E12 Bivalent_Enhancer (K4me1+K27me3) is exactly the state expected from BAP1-driven K27me3 spreading onto active enhancers. Watching its anchor-vs-span enrichment across loop sensitivity clusters in Phase 4.4 will directly test the extrusion-impediment model.

**Next available phases**:
- **Phase 3** (k-means clustering) — independent of Phase 2 outputs, only reads `cluster/data/late_merged_loop_counts.txt`. Phase 4.4 needs both Phase 2 *and* Phase 3.
- **Phase 4** can start once Phase 3 produces `combined-clusters.txt`.

Want me to proceed with Phase 3 next?


---

## Phase 3 — DONE ✅

**Result summary (v2 canonical):**

| Cluster | n      | %     | mut/ctrl med | logFC med | dominant direction      |
| ------- | ------ | ----- | ------------ | --------- | ----------------------- |
| clust1  | 12,298 | 31.6% | 0.93         | -0.03     | 100% unchanged          |
| clust2  | 10,970 | 28.2% | 1.01         | +0.09     | 92% unchanged           |
| clust3  | 8,738  | 22.4% | 0.86         | -0.14     | 21% **lost**            |
| clust4  | 3,916  | 10.1% | 1.12         | +0.25     | 70% **gained**          |
| clust5  | 667    | 1.7%  | 1.34         | +0.50     | 97% **gained** (strong) |
| clust6  | 2,359  | 6.1%  | 0.76         | -0.31     | 79% **lost** (strong)   |

Captures **99.3%** of edgeR's down_in_mutant calls and **99.5%** of up_in_mutant. Smallest cluster (clust5) has 667 loops — viable for ChromHMM enrichment in Phase 4.4.

**Files preserved on disk (nothing deleted):**
- v1 (diagnostic): `cluster/bap1_late/cluster3/{elbow_plot,k-6}_v1_no-ratio-bound/` + `cluster/phase3.txt`
- v2 (canonical, used by Phase 4+): `cluster/bap1_late/cluster3/{elbow_plot,k-6}/` + `cluster/phase3_v2.txt`

**PLAN-CLUSTER.md updated** with `Phase 3 — DONE (2026-04-27)` section in same format as Phases 1-2: status, corrections applied (filter percentile, binary path, ratio bounds rationale, sort_clusters semantics, multi-format output, Y_range, Python 3.8 typing), 3.1/3.2/3.3 method descriptions, verification block with v1↔v2 comparison table, files-created table.

Phase 3 is complete and Phase 4 can now use `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt` as input.