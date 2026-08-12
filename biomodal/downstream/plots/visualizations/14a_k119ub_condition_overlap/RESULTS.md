## Section 14: H2AK119ub Correlation at DMRs
**Key numbers:**
- Condition-overlap (ctrl vs mut) Fisher OR = 0.700, p = 4.55e-16 (source: k119ub_condition_dmr_overlap_summary.tsv)
- Hyper DMRs: 35.4% overlap mutant K119ub vs 27.8% control (gain in mutant) (source: k119ub_condition_dmr_overlap_summary.tsv)
- Differential (gained/lost) Fisher OR = 4.46, p = 5.80e-97 (from 6,172 gained, 6,103 lost peaks) (source: k119ub_differential_dmr_overlap_summary.tsv)
- Hyper DMRs 19.6% overlap K119ub-Gained vs 10.4% Lost; hypo DMRs 12.2% Gained vs 29.0% Lost (source: k119ub_differential_dmr_overlap_summary.tsv)

**What this shows:** At the DMR level the predicted convergence holds strongly — hypermethylated sites gain K119ub and hypomethylated sites lose it (differential OR=4.46), the central evidence that K119ub accumulation is upstream of the methylation switch. (Section 17 re-examines whether this survives gene-level accounting.)

**Figures:**
- `14a_k119ub_condition_overlap` — ctrl-vs-mut overlap % by DMR direction
- `14b_k119ub_differential_overlap` — gained/lost overlap % by direction
- `14c_k119ub_coordinated_genes` — K119ub-mut box + top-20 coordinated bar
- `14d_mc_vs_k119ub_scatter` — gene-level mC% vs net-K119ub scatter
- `14e_k119ub_integration_heatmap` — 2×2 mC × K119ub O/E heatmap
