All Phase 4 deliverables verified. Final summary:

**Implementation:**
- `cluster/scripts/05_grouped_analyses.py` — 8 sub-analyses orchestrated via `--analyses` CLI
- `cluster/scripts/run_phase4.sh` — driver mirroring `run_phase3.sh`
- `cluster/scripts/utils/multi_format_output.py` — patched to intercept `Figure.savefig` (covers `plotting.stacked()`)

**Outputs (`cluster/bap1_late/`):** 21 figures × 4 formats (PNG+PDF+SVG+JPG) + 12 stats files + 6 gene annotation TSVs + 2 ChromHMM enrichment matrices. Total runtime ~6 min.

**KEY biological finding (4.4 ChromHMM anchor vs span):**

| Cluster       | Polycomb Anchor | Polycomb Span | Direction          |
| ------------- | --------------: | ------------: | -----------------: |
| clust5 (gain) | **6.59×**       | **3.03×**     | 97% up_in_mutant   |
| clust6 (loss) | **2.09×**       | 0.94×         | 78% down_in_mutant |

Two mechanisms in the same paper:
- **Gained loops:** Polycomb at anchor AND span → Polycomb-Polycomb domain expansion / extrusion impediment
- **Lost loops:** Polycomb at anchor only → anchor-disruption / sensitivity model

**Supporting results:**
- 4.1: clust6 median loop = 575 kb (longest, matches abstract's "lost loops 625 kb")
- 4.2: clust5 = 55% structural (CTCF-CTCF), clust6 = 31% CRE → matches abstract
- 4.5: clust5 = 87% Polycomb-dominated anchors
- 4.8: χ² = 38,986 (p≈0) — clustering tracks edgeR direction cleanly
- 4.3: All 4 ChIP marks differential across clusters (Kruskal-Wallis p ≪ 0.001)
- 4.7: All 3 DiffBind marks chi² ≫ 1000; 0.3 vs 0.0 thresholds give very similar plots

**Documentation:**
- `cluster/plan-p2.md` — full Phase 4 — DONE record (corrections, verification, file inventory)
- `cluster/PLAN-CLUSTER.md` — DONE banner with summary + pointer to plan-p2.md
- `cluster/phase4.txt` — full run log (16 KB)

Phase 5 (deepTools metagene) is unblocked — note that PLAN-CLUSTER.md §4.3 / §5 still references `peaks/bigwigs/macs2.narrow.aug18.dedup/` paths that include the 0-byte mut files; Phase 5 should source from `/Users/zakiralibhai/sdsc/bigwigs/` exclusively.

Want me to /schedule a follow-up agent to plan Phase 5 deepTools metagene + Phase 6 cooltools pileup once you've reviewed the 4.4 heatmaps?

---


Good catch — the actual results are more nuanced than either clean model I described. Let me walk through what your data is actually showing.

## What I predicted vs. what happened

I said the extrusion impediment model would show Polycomb enriched in the **spans** of lost loops (roadblock mid-journey). But your clust6 (lost loops) shows the **opposite pattern**:

|                     | Anchor | Span               |
| ------------------- | ------ | ------------------ |
| **clust6 (lost)**   | 2.09×  | 0.94× (background) |
| **clust5 (gained)** | 6.59×  | 3.03×              |

Lost loops have Polycomb **at anchors only** — that's the sensitivity model, not extrusion impediment. The CTCF sites themselves are getting heterochromatinized, not the territory between them.

## What this means for each population

**Lost loops (clust6) — anchor disruption:** BAP1 loss → H2AK119ub elevation → H3K27me3 deposited specifically at CTCF anchor sites → CTCF can't bind → loop collapses. The span is fine — it's the endpoints that break. This is the Brad Bernstein IDH-mutant mechanism Jesse referenced: repressive chromatin directly knocking out CTCF binding.

**Gained loops (clust5) — Polycomb domain compaction:** Polycomb everywhere (anchors AND span) means these are new contacts forming **within expanding Polycomb domains**. Both anchors sit in heterochromatin, the chromatin between them is heterochromatin. These aren't extrusion-dependent loops in the normal sense — they're Polycomb-Polycomb contacts, likely driven by phase separation or PRC-mediated self-association of repressed chromatin. The fact that clust5 is 55% structural (CTCF-CTCF) means CTCF sites exist within these Polycomb domains, but the driving force is the Polycomb compaction, not cohesin extrusion.

## Why this is actually a stronger story than either model alone

You have **two mechanisms in the same dataset**, affecting different loop populations:

1. **Active loops die by anchor disruption** — long-range loops anchored at CTCF sites that become repressed. Polycomb at the anchor, span is clean. These are the CRE-enriched clust6 loops (31% CRE, median 575kb).

2. **Repressive loops form by domain compaction** — short-range Polycomb-Polycomb contacts created as heterochromatic domains expand and self-associate. Polycomb saturates both anchors and span. These are the structural clust5 loops (55% CTCF-CTCF, 97% up_in_mutant).

This connects directly to the CTEA meeting's observation: "extending along euchromatin stretch/span, but then hits polycomb K27me3 and stops." The Polycomb isn't blocking extrusion mid-span (impediment) — it's **replacing** the active chromatin landscape at specific anchor sites (sensitivity) while simultaneously **creating new compacted domains** that generate their own short-range contacts.

The H2AK119ub data ties it together: ub elevation is the upstream signal that both (a) seeds H3K27me3 at specific CTCF anchors, killing long loops, and (b) drives Polycomb domain expansion, creating new short contacts. Same cause, two consequences.