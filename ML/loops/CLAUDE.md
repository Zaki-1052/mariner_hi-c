# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**LoopBin** is an unsupervised deep learning pipeline for clustering chromatin loops based on genome interaction (micro-C) and epigenetic protein binding profiles (Cut&Tag data). It uses a Variational Autoencoder with Deep Embedding (VADE) model to learn genomic features in a latent space and separates them into clusters via Gaussian Mixture Models.

**Authors:** Yajie Zhu, Alexis Bel

## Installation & Environment

```bash
# Create and activate environment
conda env create --file LoopBin/loopbin.yml
conda activate loopbin
```

**Key dependencies:** TensorFlow 2.5.0, scikit-learn 1.0.2, cooler 0.9.3, NumPy 1.19.5, pandas 1.3.5

## Pipeline Commands

The pipeline is controlled via `main.py` with function flags:

```bash
python LoopBin/main.py -f FUNCTION [options]
```

### Step 1: Preprocess BigWig Files
Convert bigWig epigenetic files to bedGraph format (parallelized by chromosome):
```bash
python main.py -f preprocess -b BIGWIG_FILE -g OUTPUT_FOLDER -n PROTEIN_NAME
# -n: CTCF, H3K27ac, H3K27me3, or SMC1A
```

### Step 2: Process Input Data
Extract Hi-C matrices and epigenetic signals from loops:
```bash
python main.py -f process -l LOOPS.bedpe -c DATA.mcool -g BEDGRAPH_FOLDER -r NUM_CPUS -u OUTPUT_FOLDER
```

Merge and normalize across conditions:
```bash
python main.py -f normalize -e cond1,cond2,cond3 -u OUTPUT_FOLDER
```

### Step 3: Pretrain Autoencoder
```bash
python main.py -f pretrain -d processed_data.npy -u OUTPUT_FOLDER
```

### Step 4: Train VADE Model
```bash
python main.py -f train -num NUM_CLUSTERS -d data.npy -if_pre True -pre AE_MODEL -ep NUM_EPOCHS -u OUTPUT_FOLDER -p CTCF,H3K27ac,H3K27me3,SMC1A
```

### Step 5: Predict Clusters
```bash
python main.py -f cluster -d data.npy -m TRAINED_MODEL -u OUTPUT_FOLDER
```

### Step 6: Merge Small Clusters (Optional)
Merge clusters with <2% of loops:
```bash
python main.py -f merge -k 2,3 -d data.npy -u OUTPUT_FOLDER -p CTCF,H3K27ac,H3K27me3,SMC1A
```

### Additional Functions
```bash
# Cross-validation generalizability analysis (4-10 clusters)
python main.py -f calculateg -d data.npy -if_pre True -pre AE_MODEL -ep EPOCHS -u OUTPUT

# Calculate NMI between 5 independent runs
python main.py -f nmi -d MODEL_PATH_PREFIX -u OUTPUT_FOLDER
```

## Architecture

### Directory Structure
```
ML/
├── LoopBin/
│   ├── main.py                    # Entry point (argument parsing + orchestration)
│   ├── src/
│   │   ├── model/
│   │   │   ├── vade_model.py      # VADE neural network (VAE + GMM)
│   │   │   └── ae.py              # Autoencoder for pretraining
│   │   ├── fn/
│   │   │   ├── processing.py      # Hi-C matrix extraction + epigenetic signals
│   │   │   ├── function.py        # Utility functions + K-means setup
│   │   │   └── init_data.py       # Data loading/formatting
│   │   └── plot/
│   │       └── plotting.py        # Visualization (t-SNE, heatmaps, pie charts)
│   └── scripts/                   # Post-analysis scripts
│       ├── interfere_cluster_number/  # Cluster number optimization
│       ├── unique/                    # Condition-specific loop analysis
│       ├── intersect/                 # Intersection-based analysis
│       ├── cluster_shift/             # Cluster shifting between conditions
│       └── merge_clusters/            # Consensus clustering
└── loopbin_output/                # Pipeline outputs
    ├── bedgraphs/                 # Preprocessed epigenetic data
    ├── processed/                 # Feature matrices
    ├── pretrain/                  # AE model + GMM initialization
    └── train/                     # VADE model + cluster results
```

### Neural Network Architecture

**Autoencoder (pretraining):**
- Encoder: Input → 500 → 500 → 2000 → 10 (latent)
- Decoder: 10 → 2000 → 500 → 500 → Output
- Loss: Binary Cross-Entropy

**VADE Model:**
- Same encoder/decoder structure with variational sampling
- GMM layer learns cluster means, covariances, mixing weights
- Loss: Reconstruction + KL divergence + Clustering loss
- Output: Cluster labels + membership probabilities

### Data Flow

1. **Input:** Loop BEDPE + .mcool + BigWig files (CTCF, H3K27ac, H3K27me3, SMC1A)
2. **Preprocessing:** BigWig → BedGraph (8kb resolution, per chromosome)
3. **Processing:** Extract 16×16 submatrices (256 micro-C features) + epigenetic signals
4. **Normalization:** Log-max-min normalize, merge conditions
5. **Training:** AE pretraining → GMM initialization → VADE training
6. **Output:** Cluster labels, probabilities, reconstructions, visualizations

### Feature Matrix Structure
- Micro-C: 256 features (16×16 matrix flattened)
- Epigenetic: ~128 features (depends on proteins used)
- Total: ~384 features per loop

## Key Implementation Details

**Random seed:** Fixed at 73 for reproducibility in VADE training

**Training hyperparameters:**
- Batch size: 256
- Learning rate: 0.002 with decay (0.9^(epoch//10), min 0.0002)
- AE pretraining: 200 epochs
- VADE training: 1000 epochs (configurable)
- Checkpoints saved every 200 epochs

**K-means initialization:** Uses KneeLocator (elbow method) to determine initial cluster count from pretrained encoder

**Small cluster threshold:** 2% of total loops

## Output Files

After training:
```
output_folder/
├── labels.npy           # Cluster assignments (n_loops,)
├── prob.npy             # Membership probabilities (n_loops, n_clusters)
├── labels_loops.bedpe   # Original loops + cluster column
├── recon.npy            # Reconstructed features
├── checkpoints/         # Model snapshots every 200 epochs
├── cluster_{i}/         # Per-cluster heatmaps
├── cluster_all.pdf      # All clusters overview
├── pie_cluster.pdf      # Cluster distribution
├── tsne_latent.pdf      # t-SNE visualization
└── loss.pdf             # Training loss curves
```

## Analysis Scripts

Located in `LoopBin/scripts/`:

**Cluster number optimization:**
- `01_tSNE_plot_pretrained_latent_space.py` - Visualize learned representations
- `02_evaluate_GMM_bic_aic.py` - BIC/AIC scores for GMM selection
- `05_calculate_sil_scores.py` - Silhouette coefficient scoring

**Differential analysis:**
- `unique/02_barplot_perc.py` - Cluster percentage distributions
- `unique/05_violinplot_log2FC.py` - Visualize log2FC distributions
- `intersect/02_heatmap.py` - Condition overlap heatmaps

**Cluster shifting:**
- `cluster_shift/01_pair_cluster.sh` - Match clusters across conditions
- `cluster_shift/03_plot_nmi.py` - NMI heatmaps

**Consensus clustering:**
- `merge_clusters/01_find_consensu_clusters.py` - Aggregate 5 runs using co-occurrence matrices + hierarchical clustering

## Data Requirements

1. **Loop file** (.bedpe): Chromatin loop coordinates
2. **Hi-C matrix** (.mcool): Micro-C data at 8kb resolution
3. **Epigenetic BigWig files:** CTCF, H3K27ac, H3K27me3, SMC1A (or custom)

Example data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178593
