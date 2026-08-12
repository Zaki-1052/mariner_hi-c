# BAP1-KO Cerebellum: Subcompartment Analysis (ML/cmpts)

CALDER2 + SNIPER subcompartment calling pipeline for BAP1-knockout mouse cerebellum Hi-C data (mm10). Decomposes coarse A/B compartment calls into four subclasses (A.1, A.2, B.1, B.2) to distinguish true compartment flips from quantitative within-compartment weakening.

## Background

BAP1 is a Polycomb antagonist (H2AK119ub deubiquitinase). Its knockout in mouse cerebellum produces genome-wide A-compartment expansion at the subcompartment level: constitutive heterochromatin (B.2) shrinks, while the most active subcompartment (A.1) grows. This submodule quantifies that reorganization using two independent methods and integrates the results with loop, stripe, and epigenomic data from the broader project. For full biological context and detailed result interpretation, see the results guides in `docs/`.

## Three-Track Pipeline

All tracks completed as of 2026-06-02. Track A is the primary analysis (CALDER2 reference-guided subcompartment calling). Track B is secondary validation (SNIPER neural network, retrained from scratch on mm10). Track C connects subcompartments to the rest of the project.

```
Track A: CALDER2 (Primary)
  A0: Environment setup (interactive, one-time)
  A1: Run CALDER2 ── 4 jobs ── {250402,250831} x {ctrl,mut}_merged
       |
  A2: Differential subcompartments ── 2 jobs (1 per timepoint)
       |
       +--- A3: Epigenomic validation (1 job)  ──+
       |                                         |-- parallel after A2
       +--- A4: HOMER integration (1 job)     ──+

Track B: SNIPER (Secondary Validation)
  B0: Environment setup (Python 3.6 + TF-GPU 1.12)
  B1: Adapt for mm10 / generate crop maps (2 jobs)
  B2: Generate training labels from CALDER2 (2 jobs) ── depends on A2
  B3: Train SNIPER (2 GPU jobs, 1 per timepoint)
  B4: Apply SNIPER to all samples (4 GPU jobs)
  B5: Concordance analysis (1 job)

Track C: Integration (local, after A + B)
  C1: Combined comparison (developmental vs genotype effects)
  C2: Publication enrichment heatmap (ComplexHeatmap)
  C3: Loop/stripe x subcompartment integration
  C4: HOMER A/B decomposition (iceberg framing)
  C5: SNIPER concordant transitions (high-confidence subset)
```

## Key Results

| Finding | Result | Key Figure |
|---------|--------|------------|
| BAP1-KO opens the genome | A-compartment grows +3-4 pp; B.2 shrinks most. B-to-A transitions outnumber A-to-B by 3.7-4.9:1 | `outputs/calder2/late/250402_transition_heatmap/` |
| Most HOMER "switches" are not true flips | Only 2.8% of genome undergoes genuine A<->B flip; ~71% of HOMER-significant bins are within-compartment weakening or stable | `outputs/integration/homer_decomposition/plots/c4_iceberg_stacked/` |
| Polycomb loops form inside B compartment | Gained-loop cluster (clust5) is 70% B compartment; B.1 enrichment OR=3.22, p=2.96e-63 | `outputs/integration/loops_stripes/plots/c3_loop_subcompartment_enrichment/` |
| Developmental reorganization is ~2x larger than genotype effect | 32.7% bins change P12->adult (ctrl) vs 15.3% ctrl->BAP1-KO (late) | `outputs/calder2/combined/combined_stability_barplot/` |
| SNIPER independently validates CALDER2 | Cohen's kappa 0.64-0.69 across all conditions; only 3% accuracy drop in unseen mutant | `outputs/sniper/plots/250402_per_class_concordance/` |
| 707 high-confidence concordant transitions | Both tools agree on exact transition; all nearest-neighbor; B.2->B.1 dominant (29% confirmation) | `outputs/integration/sniper_concordant/plots/c5_concordant_heatmap_late/` |

## Running the Pipeline

### Prerequisites

All HPC jobs ran on SDSC Expanse (SLURM, account `csd940`). Track C scripts run locally. One-time environment setup required:

| Task | Command |
|------|---------|
| Track A env | `bash scripts/A0_setup_calder2_env.sh` (interactive, Expanse login node) |
| Track B env | See `scripts/B0_addnorm_kr.sb` and PLAN.md B0 section (interactive + gpu-debug test) |
| CALDER2 R deps | `Rscript scripts/utils/install_calder2_deps.R` (inside `calder2_env`) |

**Two-root path convention:**

| Root | Path | Contents |
|------|------|----------|
| CODE_DIR | `/expanse/.../mariner_hi-c/ML/cmpts` | Scripts, config, outputs <100MB (GitHub-synced) |
| DATA_DIR | `/expanse/.../sniper` | Large intermediates: .Rds, .h5, juicer dumps (NOT synced) |

### Track A (CALDER2, HPC)

```bash
cd /expanse/.../mariner_hi-c/ML/cmpts
bash scripts/run_full_calder2.sh                        # both timepoints, A1-A4
bash scripts/run_full_calder2.sh --tp early             # early (250831) only
bash scripts/run_full_calder2.sh --tp early --start A2  # resume from A2
```

Stages chain via `--dependency=afterok`. A3 and A4 run in parallel after A2. Estimated runtimes: A1 ~2-3 min/sample, A2 ~10s, A3 ~2.8 min, A4 ~13s.

### Track B (SNIPER, HPC)

```bash
# Crop maps (once per timepoint):
sbatch scripts/B1_generate_cropmap.sb 250402
sbatch scripts/B1_generate_cropmap.sb 250831

# Labels from CALDER2 (depends on A2):
sbatch scripts/B2_run.sb 250402
sbatch scripts/B2_run.sb 250831

# Train (GPU, depends on B1+B2):
bash scripts/B3_submit_train.sh

# Apply to all samples (GPU, depends on B3):
bash scripts/B4_submit_apply.sh

# Concordance (depends on B4):
sbatch scripts/B5_run.sb
```

### Track C (Integration, local)

```bash
cd ML/cmpts
Rscript scripts/C1_combined_comparison.R . .
Rscript scripts/C2_epigenomic_enrichment.R . .
Rscript scripts/C3_loops_stripes_integration.R . .
Rscript scripts/C4_homer_decomposition.R . .
Rscript scripts/C5_sniper_concordant_transitions.R . .
```

C3 requires loop clustering outputs (`cluster/outputs/`). C4 requires A4 outputs. C5 requires B5 and A3 outputs.

## Output Structure

```
outputs/
├── calder2/
│   ├── late/                              # 250402 results
│   │   ├── 250402/                        # raw CALDER2 per sample
│   │   │   ├── ctrl_merged/sub_compartments/
│   │   │   └── mut_merged/sub_compartments/
│   │   ├── 250402_subcompartment_labels_100kb.tsv   # primary output (A2)
│   │   ├── 250402_transition_matrix.tsv             # 4x4 counts
│   │   ├── 250402_bin_signals.tsv                   # epigenomic signals (A3)
│   │   ├── 250402_transition_heatmap/               # key figure
│   │   └── 250402_enrichment_combined/              # key figure
│   ├── early/                             # 250831 (same structure)
│   ├── combined/                          # C1: developmental comparison
│   └── enrichment/                        # C2: publication heatmaps
├── sniper/
│   ├── mm10_labels_{tp}.mat               # B2: training labels
│   ├── models_{tp}/                       # B3: 6 trained .h5 models per TP
│   ├── predictions/{tp}/{sample}/         # B4: predictions.bed + .mat
│   ├── tsvs/                              # B5: concordance tables
│   └── plots/                             # B5: confusion heatmaps, concordance
└── integration/
    ├── late/ + early/                     # A4: HOMER x CALDER2 joined tables + figures
    ├── loops_stripes/                     # C3: loop/stripe x subcompartment
    ├── homer_decomposition/               # C4: iceberg fractions, effect sizes
    └── sniper_concordant/                 # C5: high-confidence transitions
```

## Conda Environments

| Environment | Used by | Key packages |
|-------------|---------|-------------|
| `calder2_env` | A1-A4, B5, C1-C5 | R + CALDER2, data.table, rtracklayer, ComplexHeatmap, ggplot2 |
| `sniper_env` | B1-B4 | Python 3.6, TensorFlow-GPU 1.12.0, scipy, numpy, pandas |

`sniper_env` requires Python 3.6 exactly (TF 1.12 does not support 3.7+) and CUDA 9.0 + cuDNN 7.x (installed via conda, no `module load cuda` needed).

## Configuration

`config/sniper_config.yaml` centralizes HPC paths, timepoint mappings, BigWig filenames, and SLURM resource settings. Scripts currently hard-code their constants sections; the YAML is a reference and future migration target.

Key constants: 100kb resolution, autosomes chr1-19 only, 4 subcompartment classes (A.1/A.2/B.1/B.2, no B.3 in mouse), timepoints 250402=late/adult and 250831=early/P12, samples ctrl_merged and mut_merged.

## Documentation

| Document | Purpose |
|----------|---------|
| `docs/subcompartment_results_guide.md` | Track A (CALDER2) results: biology background, 3 key findings, figure-by-figure guide |
| `docs/sniper_integration_results_guide.md` | Track B (SNIPER) + Track C (integration) results: 4 additional findings, figure guide |
| `PLAN.md` | Complete implementation log: every script, decision, bug fix, and numerical result |
| `CLAUDE.md` | AI coding assistant instructions |

## Notes on mm10 Adaptation

SNIPER was originally designed for human hg19 (22 autosomes, hardcoded dimensions). Five changes were required for mm10: chromosome count reduction (22->19), asymmetric odd/even split (10 odd + 9 even vs 11+11), matrix dimension recalculation (12,808 x 11,831), removal of B.3 class (4 instead of 5 output neurons), and DAE epoch correction (upstream shipped 10 epochs, paper specifies 25). The mm10 module files are in `repos/SNIPER/utilities/` and are copies (not patches) of the originals to preserve the upstream repo.
