# Grant Figure SVG Index

Best existing SVGs for 4 grant figures, organized by theme.
All paths relative to `downstream/plots/visualizations/`.

Run the composite script to generate new patchwork versions:
```bash
cd downstream/
Rscript scripts/viz_sections/section_93_grant_figures.R 2>&1 | tee logs/section_93.txt
```

---

## Grant Figure 1: Bidirectional mC↑/hmC↓ Phenotype

**Story:** BAP1-KO causes coordinated 5mC increase + 5hmC decrease at gene bodies (78.7% of co-significant genes). Total modified C is flat = TET-block signature.

| Panel | Section | Key stat | SVG path |
|-------|---------|----------|----------|
| Dual volcano | 04 | 10,775 mC sig (70% hyper) / 11,484 hmC sig (83% down) | `04_volcano_plots/04_volcano_plots.svg` |
| Quadrant scatter | 05a | 6,589 coordinated (78.7%) | `05a_mc_hmc_scatter/05a_mc_hmc_scatter.svg` |
| Effect-size density | 07 | median mC +2.0% / hmC -1.8% | `07_effect_size_distributions/07_effect_size_distributions.svg` |
| Syt1 browser track | 46 | +17.3% mC / -15.0% hmC | `46_genome_browser/Syt1_locus_poster_v2_aligned/Syt1_locus_poster_v2_aligned.svg` |
| **Prior composite** | old-figs | 7-panel Figure 1 | `old-figs/figure1_methylation_phenotype/figure1_methylation_phenotype.svg` |

**Also useful:**
- `05b_top_coordinated_genes/05b_top_coordinated_genes.svg` — labeled top genes
- `05c_syt1_detail/05c_syt1_detail.svg` — Syt1 zoomed scatter
- `46_genome_browser/Syt1_locus_compact.svg` — smaller browser version
- `46_genome_browser/Zbtb20_locus.svg` — second-best example gene

---

## Grant Figure 2: Euchromatin Silencing

**Story:** H2AK119ub gain is the dominant predictor of hypermethylation (OR=4.71), and that hypermethylation is restricted to active, A-compartment chromatin.

### Chromatin State Classification System (7 categories)

Used throughout all figures referencing chromatin context. Shared with the Hi-C loop pipeline.

| State | Definition | Color |
|-------|-----------|-------|
| Active_Promoter | H3K4me3+ AND NOT H3K27me3, within 2kb of TSS | Red |
| Repressed_Promoter | H3K27me3+ AND NOT H3K27ac, within 2kb of TSS | Purple |
| Bivalent_Promoter | H3K4me3 + H3K27me3 overlap (pre-computed), within 2kb of TSS | Magenta |
| Polycomb | H3K27me3+ AND >2kb from TSS | Green |
| Active_Enhancer | H3K27ac+ AND >2kb from TSS | Blue |
| Poised_Enhancer | H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3, >2kb from TSS | Orange |
| Unmarked | No histone mark overlap | Gray |

Priority order: Active_Promoter > Repressed_Promoter > Bivalent_Promoter > Polycomb > Active_Enhancer > Poised_Enhancer > Unmarked (first match wins when a gene overlaps multiple marks).

| Panel | Section | Key stat | SVG path |
|-------|---------|----------|----------|
| Chromatin-state direction bar | 10b/10f | Active_Promoter 93% hyper, Repressed 94% hypo | `10f_chromatin_stacked_presentation/10f_chromatin_stacked_presentation.svg` |
| K119ub logistic forest plot | 33e | K119ub OR=4.71 (dominant) | `33_multi_mark_diffbind/33e_logistic_regression_forest/33e_logistic_regression_forest.svg` |
| A-compartment enrichment | 29e | mC-hyper in A: OR=13.6 | `29_ab_compartment_methylation/29e_dmr_direction_stacked_bar/29e_dmr_direction_stacked_bar.svg` |
| Compartment composite | 29g | A/B summary | `29_ab_compartment_methylation/29g_composite_compartment_summary/29g_composite_compartment_summary.svg` |
| Cross-mark correlation | 33b | 4-mark heatmap | `33_multi_mark_diffbind/33b_cross_mark_correlation_heatmap/33b_cross_mark_correlation_heatmap.svg` |
| **Prior composite** | old-figs | 5-panel Figure 2 | `old-figs/figure2_k119ub_chromatin_geography/figure2_k119ub_chromatin_geography.svg` |

**Also useful:**
- `10b_chromatin_by_methylation_direction/10b_chromatin_by_methylation_direction.svg`
- `10a_chromatin_state_distribution/10a_chromatin_state_distribution.svg`
- `33_multi_mark_diffbind/33d_methylation_vs_mark_scatters/33d_methylation_vs_mark_scatters.svg`

---

## Grant Figure 3: Heterochromatin Protection

**Story:** Polycomb/heterochromatin targets are EXCLUDED from hypermethylation (OR=0.063) and enriched in hypomethylation. Heterochromatin is protected — counter to naive predictions.

| Panel | Section | Key stat | SVG path |
|-------|---------|----------|----------|
| Per-state hyper rate | 30d | Active_Promoter 72%, Polycomb 1.8% | `30_polycomb_enrichment/30d_per_state_hypermethylation_rate/30d_per_state_hypermethylation_rate.svg` |
| Polycomb exclusion forest | 30b | Polycomb x hyper OR=0.063 | `30_polycomb_enrichment/30b_fisher_forest_plot/30b_fisher_forest_plot.svg` |
| Polycomb stacked bar | 30a | Polycomb vs non-Polycomb direction split | `30_polycomb_enrichment/30a_polycomb_vs_non_polycomb_stacked_bar/30a_polycomb_vs_non_polycomb_stacked_bar.svg` |
| Polycomb composite | 30f | 6-panel summary | `30_polycomb_enrichment/30f_composite_polycomb_summary/30f_composite_polycomb_summary.svg` |
| B-compartment mC violin | 29a | B-compartment methylation distribution | `29_ab_compartment_methylation/29a_mc_violin_by_compartment/29a_mc_violin_by_compartment.svg` |
| B-compartment hmC violin | 29b | B-compartment 5hmC distribution | `29_ab_compartment_methylation/29b_hmc_violin_by_compartment/29b_hmc_violin_by_compartment.svg` |

**Also useful (Q2 discordant = heterochromatin de-repression):**
- `21a_discordant_composite/21a_discordant_composite.svg` — Q2 discordant characterization composite
- `21a_panel_chromatin/21a_panel_chromatin.svg` — chromatin state of Q2 genes
- `21a_panel_k119ub/21a_panel_k119ub.svg` — K119ub at Q2 (losing K119ub)
- `21a_panel_k27ac/21a_panel_k27ac.svg` — K27ac at Q2 (gaining)

---

## Grant Figure 4: MeCP2 Reads Chromatin State, Not Methylation

**Story:** MeCP2 redistribution tracks H2AK119ub (chromatin R²=0.246), not CG methylation (R²=0.017). At 359 genes where MeCP2 rises with NO methylation change, K119ub is gained (OR=3.15).

| Panel | Section | Key stat | SVG path |
|-------|---------|----------|----------|
| R² comparison bars | 62a | CG 0.017 vs Chromatin 0.246 | `62_multifeature_chromatin_regression/62a_model_r2_comparison/62a_model_r2_comparison.svg` |
| Variance partition | 62b | CG-unique 1.5%, Chromatin-unique 24.3% | `62_multifeature_chromatin_regression/62b_variance_partition/62b_variance_partition.svg` |
| K119ub at MeCP2-no-meth | 67a | 72.8% vs 45.9% gain K119ub, OR=3.15 | `67a_k119ub_at_mecp2_no_meth/67a_k119ub_at_mecp2_no_meth.svg` |
| K119ub gained fraction | 67b | Grouped bar of gained % | `67b_k119ub_gained_fraction/67b_k119ub_gained_fraction.svg` |
| Mirror-image profiles | 60e | MeCP2-Up/Down per mark | `60_mecp2_status_epigenetic_profiles/60e_composite/60e_composite.svg` |
| Summary heatmap | 60b | Mark changes by MeCP2 status | `60_mecp2_status_epigenetic_profiles/60b_summary_heatmap/60b_summary_heatmap.svg` |
| Regression composite | 62g | Full regression summary | `62_multifeature_chromatin_regression/62g_composite/62g_composite.svg` |
| **Prior composite** | old-figs | 5-panel Figure 5 | `old-figs/figure5_mecp2_chromatin_reader/figure5_mecp2_chromatin_reader.svg` |

**Also useful:**
- `59a_k119ub_vs_mecp2/59a_k119ub_vs_mecp2.svg` — K119ub vs MeCP2 scatter
- `67c_mecp2_vs_k119ub_scatter/67c_mecp2_vs_k119ub_scatter.svg` — MeCP2 fold vs K119ub at 359 genes
- `62_multifeature_chromatin_regression/62d_standardized_coefficients/62d_standardized_coefficients.svg`

---

## Copy Command

Copies the single best SVG per panel into `originals/` for easy browsing:

```bash
cd downstream/plots/visualizations
mkdir -p 93_grant_figures/originals
cp 04_volcano_plots/04_volcano_plots.svg 93_grant_figures/originals/fig1_dual_volcano.svg
cp 05a_mc_hmc_scatter/05a_mc_hmc_scatter.svg 93_grant_figures/originals/fig1_quadrant_scatter.svg
cp 07_effect_size_distributions/07_effect_size_distributions.svg 93_grant_figures/originals/fig1_effect_size_density.svg
cp 46_genome_browser/Syt1_locus_poster_v2_aligned/Syt1_locus_poster_v2_aligned.svg 93_grant_figures/originals/fig1_syt1_browser.svg
cp 10f_chromatin_stacked_presentation/10f_chromatin_stacked_presentation.svg 93_grant_figures/originals/fig2_chromatin_direction.svg
cp 33_multi_mark_diffbind/33e_logistic_regression_forest/33e_logistic_regression_forest.svg 93_grant_figures/originals/fig2_k119ub_forest.svg
cp 29_ab_compartment_methylation/29e_dmr_direction_stacked_bar/29e_dmr_direction_stacked_bar.svg 93_grant_figures/originals/fig2_a_compartment.svg
cp 29_ab_compartment_methylation/29g_composite_compartment_summary/29g_composite_compartment_summary.svg 93_grant_figures/originals/fig2_compartment_composite.svg
cp 30_polycomb_enrichment/30d_per_state_hypermethylation_rate/30d_per_state_hypermethylation_rate.svg 93_grant_figures/originals/fig3_per_state_hyper.svg
cp 30_polycomb_enrichment/30b_fisher_forest_plot/30b_fisher_forest_plot.svg 93_grant_figures/originals/fig3_polycomb_exclusion.svg
cp 30_polycomb_enrichment/30a_polycomb_vs_non_polycomb_stacked_bar/30a_polycomb_vs_non_polycomb_stacked_bar.svg 93_grant_figures/originals/fig3_polycomb_bar.svg
cp 30_polycomb_enrichment/30f_composite_polycomb_summary/30f_composite_polycomb_summary.svg 93_grant_figures/originals/fig3_polycomb_composite.svg
cp 21a_discordant_composite/21a_discordant_composite.svg 93_grant_figures/originals/fig3_q2_discordant.svg
cp 62_multifeature_chromatin_regression/62a_model_r2_comparison/62a_model_r2_comparison.svg 93_grant_figures/originals/fig4_r2_comparison.svg
cp 62_multifeature_chromatin_regression/62b_variance_partition/62b_variance_partition.svg 93_grant_figures/originals/fig4_variance_partition.svg
cp 67a_k119ub_at_mecp2_no_meth/67a_k119ub_at_mecp2_no_meth.svg 93_grant_figures/originals/fig4_k119ub_linchpin.svg
cp 67b_k119ub_gained_fraction/67b_k119ub_gained_fraction.svg 93_grant_figures/originals/fig4_k119ub_gained.svg
cp 60_mecp2_status_epigenetic_profiles/60e_composite/60e_composite.svg 93_grant_figures/originals/fig4_mirror_profiles.svg
cp 60_mecp2_status_epigenetic_profiles/60b_summary_heatmap/60b_summary_heatmap.svg 93_grant_figures/originals/fig4_mecp2_heatmap.svg
cp old-figs/figure1_methylation_phenotype/figure1_methylation_phenotype.svg 93_grant_figures/originals/prior_fig1_composite.svg
cp old-figs/figure2_k119ub_chromatin_geography/figure2_k119ub_chromatin_geography.svg 93_grant_figures/originals/prior_fig2_composite.svg
cp old-figs/figure5_mecp2_chromatin_reader/figure5_mecp2_chromatin_reader.svg 93_grant_figures/originals/prior_fig4_composite.svg
```

---

## New Composites (from section_93_grant_figures.R)

After running the script, these appear in `93_grant_figures/`:

| File | Content |
|------|---------|
| `grant_fig1_bidirectional/` | 3-panel composite: dual volcano + quadrant + density |
| `grant_fig2_euchromatin/` | 3-panel composite: chromatin bar + K119ub forest + A-compartment |
| `grant_fig3_heterochromatin/` | 3-panel composite: per-state lollipop + Polycomb forest + B-compartment |
| `grant_fig4_mecp2_chromatin/` | 3-panel composite: R² bars + K119ub linchpin + mirror profiles |
| `93a_dual_volcano/` | Individual panel |
| `93b_quadrant_scatter/` | Individual panel |
| `93c_effect_size_density/` | Individual panel |
| `93d_chromatin_direction_bar/` | Individual panel |
| `93e_k119ub_forest/` | Individual panel |
| `93f_a_compartment_enrichment/` | Individual panel |
| `93g_per_state_hyper_rate/` | Individual panel |
| `93h_polycomb_exclusion_forest/` | Individual panel |
| `93i_b_compartment_enrichment/` | Individual panel |
| `93j_mecp2_r2_comparison/` | Individual panel |
| `93k_k119ub_mecp2_no_meth/` | Individual panel |
| `93l_mirror_profiles/` | Individual panel |
