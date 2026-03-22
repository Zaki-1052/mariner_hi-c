# Hi-C Paper — Figure Panel Layouts (Version 2)

---

## Figure 1: Progressive Intensification of 3D Genome Phenotype
*TADs, boundaries, representative compartments*

```
┌─────────────────────────────────┬─────────────────────────────────┐
│                                 │                                 │
│        P13 Contact Map          │       Adult Contact Map         │
│                                 │                                 │
├────────────────┬────────────────┼────────────────┬────────────────┤
│ Insulation trk │                │ Insulation trk │                │
├────────────────┴────────────────┼─────────────────────────────────┤
│  Diff TAD ID     ┌───────────┐ │  P13 GENOVA                     │
│  P13 Volcano     │ Diff TAD  │ │  TAD Compare Boundaries         │
│                  │ ID Adult  │ ├─────────────────────────────────┤
│                  │ Volcano   │ │  Adult                          │
│                  └───────────┘ │  Boundaries                     │
├──────────────┬──────┬──────────┼─────────────────────────────────┤
│              │      │ Diff TAD │  Perm Test                      │
│  Diff Cmpt   │  ◇   │ P13      │  diff TAD / diff boundaries     │
│  Region      │ pie  ├──────────┤  in K27me3, K27ac,              │
│  Contact Map │chart │ Diff TAD │  K119ub, A&B cmpt               │
│              │      │ Adult    │                                 │
├──────────────┴──────┴──────────┴─────────────────────────────────┤
│                                                                   │
│  ┌─────┐  ── weakened PC1 ──►                                     │
│  │ P13 │  ── strengthened ──►  ┌───────┐                          │
│  └─────┘  ── flipped ───────►  │ Adult │                          │
│           ── no change ─────►  └───────┘                          │
│  (% of genome with differential PC1)                              │
└───────────────────────────────────────────────────────────────────┘
```

---

## Figure 2: Chromatin Loop Rewiring
*Differential loop analysis, APA, annotation w/ ctrl peaks*

```
┌──────────────────────┬───────────────────────────────────────────┐
│  P13 Volcano         │                                           │
│  (loops)             │            APA Analysis                   │
├──────────────────────┤            P13 & Adult                    │
│  Adult Volcano       │            (GENOVA)                       │
│  (loops)             │                                           │
├────────────┬─────────┴──┬────────────────────────────────────────┤
│            │            │                                        │
│  Lost Loop │ Gained Loop│   Anchor Annotation                    │
│  Contact   │ Contact    │   P13 & Adult                          │
│  Map       │ Map        │   (genomic anchor distribution)        │
├──────────┬─┴────────────┴────────────────────────────────────────┤
│ K27ac trk│  K27me3 track                                        │
├──────────┴───────────────────────────────────────────────────────┤
│                                                                   │
│              # Loops by Category (bar plot)                       │
│                                                                   │
├───────────────────────────────────────────────────────────────────┤
│                                                                   │
│           Loop logFC by Category (loop strength changes)          │
│                                                                   │
├───────────────┬───────────────┬─────────────────┬─────────────────┤
│  P13 eCDF     │  Adult eCDF   │ Adult K27me3    │ Adult K27ac     │
│  all loops    │  all loops    │ -anchored eCDF  │ -anchored eCDF  │
│  (lost v gain)│  (lost v gain)│                 │                 │
└───────────────┴───────────────┴─────────────────┴─────────────────┘
```

---

## Figure 3: Integration with Epigenetic Changes

```
┌───────────────────────────────────────────────────────────────────┐
│                                                                   │
│          ARA Plots for ATAC, K27me3, K27ac, K119ub               │
│                                                                   │
├───────────────────────────────┬───────────────────────────────────┤
│                               │                                   │
│  Histone Mark O/E at          │   Perm Test                       │
│  Diff TAD Boundaries          │   DPA peaks in structures         │
│                               │   (adult diff peaks → adult str.) │
├───────────────────────────────┼───────────────────────────────────┤
│                               │                                   │
│   △ ctrl     △ mut            │   Perm Test                       │
│  ─────────────────            │   P12 marks predict               │
│  (contact maps at             │   adult structures                │
│   diff boundary)              │                                   │
├───────────────────────────────┴───────────────────────────────────┤
│                                                                   │
│             Integrated Locus (probably Syt1)                      │
│                                                                   │
├───────────────────────────────────────────────────────────────────┤
│  + Possibly H2Az variant dynamics (if published by summer)        │
└───────────────────────────────────────────────────────────────────┘
```

---

## Figure 4: Enhancer ABC Analysis

```
┌───────────────────────────────┬───────────────────────────────────┐
│                               │                                   │
│   ABC–RNA Correlation         │   Concordance Analysis            │
│   (unnormalized score)        │   ○ ○  (pie charts)              │
│                               │   957 DEGs                        │
├───────────────────────────────┴───────────────────────────────────┤
│           Combined ctrl ATAC peaks (DB)                           │
├──────────────────┬──────────────────┬─────────────────────────────┤
│                  │                  │                              │
│   K27ac Volcano  │  ATAC Volcano    │   K119ub Volcano            │
│                  │                  │                              │
├──────────────────┴─┬────────────────┴─────────────────────────────┤
│                    │                                              │
│  Model for         │   Loop logFC          ~~~~~~~~~~~~~~~~       │
│  Subsetting        │   Δ Activity          ~~~~~~~~~~~~~~~~       │
│  Enhancers         │   Δ Contacts          ~~~~~~~~~~~~~~~~       │
│                    │   Gene logFC          ~~~~~~~~~~~~~~~~       │
├────────────────────┤                                              │
│  ABC Category      │                                              │
│  vs K119ub         │                                              │
│  Significance      │                                              │
└────────────────────┴──────────────────────────────────────────────┘
```

---

## Figure 5: Model and Functional Implications

```
┌───────────────────────────────────────────────────────────────────┐
│                                                                   │
│       Violin Plots: DEG logFC for genes near                     │
│       diff boundaries, at diff loops, with lost AxC connections   │
│                                                                   │
├───────────────────────────────┬───────────────────────────────────┤
│                               │                                   │
│                               │                                   │
│   Network Analysis            │   Heatmap (top 50 genes)          │
│   nodes = genes w/ combined   │   cols: logFC genes, logFC AxC    │
│     structural changes        │                                   │
│   node size = AxC Δ score     │                                   │
│   node color = gene logFC     │                                   │
│                               │                                   │
│                               │                                   │
├───────────────────────────────┴───────────────────────────────────┤
│  Conceptual Model                                                 │
│                        P13                Adult                    │
│              ┌─────────────────┬─────────────────┐                │
│              │ Epigenetic Layer│ Epigenetic Layer │                │
│              │        ↓       │        ↓         │                │
│              │ Structural Layer│ Structural Layer │                │
│              │        ↓       │        ↓         │                │
│              │ Functional (DEG)│ Functional (DEG)│                │
│              └─────────────────┴─────────────────┘                │
└───────────────────────────────────────────────────────────────────┘
```
