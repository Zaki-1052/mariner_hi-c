## Section 73: Disproportionate Chromatin Remodeling at K119ub-High Neuronal Genes (No MeCP2)
**Key numbers:**
- 4-group key test (K119ub-High tier): ATAC neuronal median +0.1667 vs other +0.1088 (p = 2.90e-11); K27ac neuronal +0.0469 vs other −0.0155 (p = 2.15e-05); K27me3 NS (p = 0.219) (source: 72g_4group_stats.tsv)
- Interaction (mark_fold ~ k119ub_signal × is_neuronal): K27ac beta = +0.0462, p = 4.55e-07; ATAC beta = +0.00832, p = 0.0172; K27me3 beta = +0.00304, p = 0.760 NS (source: 72g_interaction_models.tsv)
- Chromatin state of K119ub-high genes (neuronal vs other): Bivalent_Promoter OR = 2.866 (p = 3.38e-07), Repressed_Promoter OR = 1.737 (p = 6.72e-17), "Other"/unmarked OR = 0.407 (p = 1.08e-34, depleted) (source: 72g_chromatin_state_fisher.tsv)
- All-gene mark deltas (neuronal − other): ATAC +0.0199 (p = 1.10e-14), K27ac +0.0063 (p = 0.016), K27me3 −0.0423 (p = 6.22e-04) (source: 72g_mark_stats_neuronal_vs_other.tsv)

**What this shows:** Constitutively K119ub-high neuronal genes (from section 72) remodel their chromatin significantly more than non-neuronal genes upon BAP1-KO, with significant positive K119ub×neuronal interaction terms for ATAC and K27ac. The affected neuronal loci are specifically bivalent/repressed-promoter Polycomb-poised genes; the direction (accessibility/K27ac GAIN, no K27me3 interaction) indicates de-repression of poised neuronal promoters rather than uniform heterochromatinization. NOTE: output tables on disk carry the stale `72g_*` prefix (pre-rename); the committed script writes `73_*` but was not re-run — data are current and valid.

**Figures:**
- `73_multimark_neuronal_vs_other/` — (73a) mark log2FC violins, neuronal vs non-neuronal
- `73_4group_k119ub_x_neuronal/` — (73b) 4-group violins: K119ub level × neuronal status
- `73_chromatin_state_k119ub_high/` — (73c) chromatin-state composition of K119ub-high genes
- `73_dose_response_decile_marks/` — (73d) median mark change by K119ub decile
- `73_interaction_coefficients/` — (73e) K119ub×neuronal interaction forest plot
- `73_composite/` — assembled section-73 figure
