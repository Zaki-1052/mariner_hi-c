# LoopBin: Unsupervised Chromatin Loop Clustering

Deep learning-based classification of chromatin loops in BAP1-KO vs wildtype mouse cerebellum (mm10) using the [LoopBin](https://github.com/YajieZhu018/LoopBin) framework. LoopBin applies a Variational Autoencoder with Deep Embedding (VaDE) to learn joint representations of Hi-C contact matrices and epigenetic binding profiles, then clusters loops into functionally distinct subtypes via Gaussian Mixture Models.

## Biological Context

Chromatin loops identified by micro-C connect distal regulatory elements (enhancers) to gene promoters, but not all loops are functionally equivalent. By integrating Hi-C interaction strength with H3K27ac (active) and H3K27me3 (repressive) histone modification profiles at loop anchors, LoopBin separates loops into clusters that reflect distinct chromatin states. In the BAP1-KO system, this reveals which loop classes are most affected by loss of H2AK119ub deubiquitinase activity.

## Model Architecture

```
                     Input Features (~384 per loop)
                ┌────────────┴────────────────┐
         Micro-C Matrix              Epigenetic Signals
        (16x16 = 256 features)      (~128 features)
                │                          │
                └──────────┬───────────────┘
                           │
              ┌────────────▼────────────┐
              │  Autoencoder Pretraining │
              │  500 → 500 → 2000 → 10  │  (200 epochs, BCE loss)
              └────────────┬────────────┘
                           │
              ┌────────────▼────────────┐
              │       VaDE Training      │
              │  VAE + GMM (10 latent)   │  (2000 epochs)
              │  Recon + KL + Clustering │
              └────────────┬────────────┘
                           │
              ┌────────────▼────────────┐
              │     Cluster Labels       │
              │  + Membership Probs      │
              └─────────────────────────┘
```

**Encoder:** Input -> 500 -> 500 -> 2000 -> 10 (latent dim)
**Decoder:** 10 -> 2000 -> 500 -> 500 -> Input (sigmoid)
**Initialization:** K-means with elbow method (KneeLocator) on pretrained latent space, followed by GMM fitting

## Requirements

```bash
conda env create --file LoopBin/loopbin.yml
conda activate loopbin
```

Key dependencies: TensorFlow 2.5.0, scikit-learn 1.0.2, cooler 0.9.3, NumPy 1.19.5, pandas 1.3.5, kneed, bigWigAverageOverBed (UCSC tools)

## Input Data

| File | Format | Description |
|---|---|---|
| Loop coordinates | `.bedpe` | Chromatin loop anchor pairs (from SIP/HiCCUPS/etc.) |
| Hi-C matrix | `.mcool` | Micro-C contact matrix (used at 8 kb resolution) |
| H3K27ac | `.bw` | Active enhancer/promoter histone mark |
| H3K27me3 | `.bw` | Polycomb repressive mark |

For this project, loops come from the merged non-redundant set (2,413 differential loops), the mcool is `ctrl_merged.mcool`, and epigenetic marks are condition-merged bigWigs from CUT&Tag.

## Running the Pipeline

All steps are orchestrated through `main.py -f FUNCTION`. The SLURM wrapper `run_loopbin_mm10.sb` runs the full pipeline end-to-end on Expanse.

### Quick start (HPC)

```bash
cd /expanse/lustre/projects/csd940/zalibhai/loop-class
sbatch ML/LoopBin/run_loopbin_mm10.sb
```

### Step-by-step execution

#### Step 1: Preprocess bigWigs to bedGraphs

Converts bigWig signal to per-chromosome bedGraph files at 8 kb resolution using `bigWigAverageOverBed`. Parallelized across chromosomes (chr1-19, chrX).

```bash
# Via main.py (hg38 default)
python main.py -f preprocess -b BIGWIG_FILE -g OUTPUT_FOLDER -n PROTEIN_NAME

# Via mm10-adapted shell script (used in this project)
bash src/preprocess_local_mm10.sh H3K27ac.bw H3K27ac output_dir LoopBin/src/fn
bash src/preprocess_local_mm10.sh H3K27me3.bw H3K27me3 output_dir LoopBin/src/fn
```

Output: `{chr}_{mark}8K.bedgraph` per chromosome per mark (40 files for 2 marks x 20 chromosomes).

#### Step 2: Process input data

Extracts 16x16 Hi-C submatrices centered on each loop and samples epigenetic signal at both anchors. Produces a feature matrix with ~384 features per loop.

```bash
# Extract features
python main.py -f process \
    -l loops.bedpe \
    -c data.mcool \
    -g bedgraph_dir/ \
    -r 8 \
    -u output_dir/ctrl \
    -p H3K27ac,H3K27me3

# Normalize (log-max-min) and merge conditions
python main.py -f normalize \
    -e ctrl \
    -u output_dir/
```

Output: `merged_log_data.npy` — the training input matrix.

#### Step 3: Pretrain autoencoder

Trains a symmetric autoencoder (500-500-2000-10-2000-500-500) with binary cross-entropy loss. The 10-dimensional latent space is used for K-means initialization and GMM fitting.

```bash
python main.py -f pretrain \
    -d processed_dir/merged_log_data.npy \
    -u pretrain_dir/
```

Output: saved AE model, `elbow.pdf` (K-means inertia curve), `sil.pdf` (silhouette scores), `model_cluster.pkl` (GMM initialization).

#### Step 4: Train VaDE model

Joint variational autoencoder + GMM training. The loss combines reconstruction error, KL divergence from the GMM prior, and clustering loss. Checkpoints every 200 epochs.

```bash
python main.py -f train \
    -num 10 \
    -d processed_dir/merged_log_data.npy \
    -if_pre True \
    -pre pretrain_dir/ \
    -ep 2000 \
    -u train_dir/ \
    -p H3K27ac,H3K27me3
```

Output: `labels.npy`, `labels_loops.bedpe`, `prob.npy`, per-cluster heatmaps, t-SNE plot, loss curve.

#### Step 5: Predict clusters (re-apply trained model)

```bash
python main.py -f cluster \
    -d merged_log_data.npy \
    -m train_dir/ \
    -u output_dir/
```

#### Step 6: Merge small clusters (optional)

Merge clusters containing < 2% of loops into their nearest neighbor:

```bash
python main.py -f merge \
    -k 1,5 \
    -d merged_log_data.npy \
    -u train_dir/ \
    -p H3K27ac,H3K27me3
```

## Results

### Cluster Distribution

Training on 2,413 control-condition loops with 10 initial clusters produced:

| Cluster | Count | % of Total |
|---|---|---|
| 0 | 617 | 25.5% |
| 9 | 505 | 20.9% |
| 6 | 452 | 18.7% |
| 3 | 252 | 10.4% |
| 8 | 247 | 10.2% |
| 7 | 139 | 5.8% |
| 4 | 109 | 4.5% |
| 2 | 81 | 3.4% |
| 1 | 6 | 0.2% |
| 5 | 5 | 0.2% |

Clusters 1 and 5 (< 2% each) are candidates for merging.

### Key Output Files

```
loopbin_output/
├── bedgraphs/                    # Per-chrom bedGraphs (40 files)
│   ├── chr1_H3K27ac8K.bedgraph
│   └── chr1_H3K27me38K.bedgraph
├── processed/
│   ├── ctrl/                     # Raw extracted features
│   ├── loop_file_analysis.bedpe  # Loop coordinates used
│   └── merged_log_data.npy       # Normalized feature matrix (training input)
├── pretrain/
│   ├── saved_model.pb            # Pretrained AE weights
│   ├── model_cluster.pkl         # GMM initialization params
│   ├── elbow.pdf                 # K-means elbow plot
│   └── sil.pdf                   # Silhouette score plot
└── train/
    ├── saved_model.pb            # Trained VaDE weights
    ├── labels.npy                # Cluster assignments (2413,)
    ├── prob.npy                  # Membership probabilities (2413, 10)
    ├── labels_loops.bedpe        # Loop coordinates + cluster label (col 7)
    ├── recon.npy                 # Reconstructed features
    ├── all_clusters.pdf          # Per-cluster heatmap overview
    ├── cluster_{0-9}.pdf         # Individual cluster heatmaps
    ├── tsne_perplex100.pdf       # t-SNE visualization (perplexity=100)
    ├── pie.pdf                   # Cluster size distribution
    ├── loss.pdf                  # Training loss curve
    └── cluster_merged/           # Post-merge results
        ├── labels.npy
        └── all_clusters.pdf
```

The primary output for downstream analysis is `labels_loops.bedpe`: a 7-column file where columns 1-6 are BEDPE loop coordinates and column 7 is the cluster assignment.

## Post-hoc Analysis Scripts

Located in `LoopBin/scripts/`, these scripts characterize the biological meaning of each cluster:

### Cluster number optimization (`interfere_cluster_number/`)

| Script | Purpose |
|---|---|
| `01_tSNE_plot_pretrained_latent_space.py` | Visualize AE latent space via t-SNE |
| `02_evaluate_GMM_bic_aic.py` | BIC/AIC model selection for GMM |
| `05_calculate_sil_scores.py` | Silhouette scoring across k values |

### DEG integration (`degs/`)

| Script | Purpose |
|---|---|
| `00_extract_TSS.sh` | Extract TSS from gene annotations |
| `01_annotate_genes_loops.sh` | Assign nearest genes to loop anchors |
| `02_geneID_to_geneName.sh` | Convert Ensembl IDs to gene symbols |
| `03_add_fkpm_to_loops.sh` | Annotate loops with FPKM expression |
| `04_violinplot_pkfm.py` | Per-cluster expression violin plots |
| `05_barplot_perc_loops_with_promoters.py` | Promoter-loop fraction by cluster |

### Condition-specific loops (`unique/`)

| Script | Purpose |
|---|---|
| `01_unique.sh` | Identify condition-unique loops |
| `02_barplot_perc.py` | Cluster distribution bar plots |
| `05_violinplot_log2FC.py` | Per-cluster log2FC violin plots |

### Cross-condition analysis (`intersect/`)

| Script | Purpose |
|---|---|
| `01_intersect.sh` | Intersect loops between conditions |
| `02_heatmap.py` | Condition overlap heatmaps |
| `04_plot_length.py` | Loop length distributions by cluster |

### Cluster shifting (`cluster_shift/`)

| Script | Purpose |
|---|---|
| `01_pair_cluster.sh` | Match cluster identities across conditions |
| `02_heatmap.py` | Cluster correspondence heatmaps |
| `03_plot_nmi.py` | Normalized mutual information across runs |

### Consensus clustering (`merge_clusters/`)

`01_find_consensu_clusters.py` — Aggregates cluster labels from 5 independent runs using co-occurrence matrices and hierarchical clustering to identify stable clusters robust to initialization.

### Other utilities

- `probability/` — Violin plots of cluster membership probability; filter loops by confidence
- `metaplot/` — Format loop coordinates for metaplot visualization (APA-style)

## Hyperparameters

| Parameter | Value | Notes |
|---|---|---|
| Hi-C resolution | 8 kb | Submatrix = 16x16 bins = 128 kb window |
| Latent dimension | 10 | Shared by AE and VaDE |
| AE pretraining epochs | 200 | BCE loss |
| VaDE training epochs | 2000 | Recon + KL + clustering loss |
| Batch size | 256 | |
| Learning rate | 0.002 | Decays 0.9^(epoch//10), min 0.0002 |
| Initial clusters | 10 | VaDE auto-prunes redundant clusters |
| Small cluster threshold | 2% | Candidates for post-hoc merging |
| Random seed | 73 | Fixed for reproducibility |
| t-SNE perplexity | 100 | For visualization |

## mm10 Adaptation

The original LoopBin was developed for hg38. Adaptations for mm10:

- `src/fn/preprocessing_mm10.sh` — mm10 chromosome sizes (chr1-19, chrX; no chrY)
- `src/preprocess_local_mm10.sh` — Batched parallel preprocessing (3 batches to limit concurrent jobs)
- `run_loopbin_mm10.sb` — End-to-end SLURM wrapper with mm10-specific paths and only 2 epigenetic marks (H3K27ac, H3K27me3) instead of the original 4 (CTCF, H3K27ac, H3K27me3, SMC1A)

## Project Structure

```
ML/
├── LoopBin/
│   ├── main.py                         # CLI entry point (-f flag dispatch)
│   ├── loopbin.yml                     # Conda environment spec
│   ├── run_loopbin_mm10.sb             # SLURM end-to-end pipeline (mm10)
│   ├── src/
│   │   ├── model/
│   │   │   ├── ae.py                   # Autoencoder (pretraining)
│   │   │   └── vade_model.py           # VaDE: VAE + GMM joint model
│   │   ├── fn/
│   │   │   ├── processing.py           # Feature extraction (Hi-C + epigenetics)
│   │   │   ├── function.py             # K-means/GMM init, file validation
│   │   │   ├── init_data.py            # Data loading utilities
│   │   │   ├── preprocessing_mm10.sh   # Per-chrom bigWig → bedGraph (mm10)
│   │   │   └── preprocessing.sh        # Per-chrom bigWig → bedGraph (hg38)
│   │   ├── plot/
│   │   │   └── plotting.py             # t-SNE, heatmaps, cluster viz
│   │   ├── preprocess_local_mm10.sh    # Parallel preprocessing driver (mm10)
│   │   ├── preprocess_local.sh         # Parallel preprocessing driver (hg38)
│   │   └── jobs/                       # Reference SLURM job templates
│   └── scripts/                        # Post-hoc analysis (see above)
│       ├── interfere_cluster_number/   # Cluster count optimization
│       ├── degs/                       # DEG annotation at loops
│       ├── unique/                     # Condition-specific loops
│       ├── intersect/                  # Cross-condition overlap
│       ├── cluster_shift/              # Cluster identity matching
│       ├── merge_clusters/             # Consensus clustering
│       ├── probability/                # Membership confidence
│       └── metaplot/                   # APA-style visualization prep
├── loopbin_output/                     # Pipeline outputs (see Results)
├── clusters.md                         # Cluster count summary table
└── logs/                               # SLURM job logs
```

## References

Zhu, Y. & Bel, A. LoopBin: an unsupervised neural network for chromatin loop clustering. [GitHub](https://github.com/YajieZhu018/LoopBin)
