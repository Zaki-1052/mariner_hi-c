# HiC-clustering (Original Popay README)

> **Note:** This is the original README from the [Popay et al. HiC-clustering repo](https://github.com/tpopay/HiC-clustering). For the BAP1-KO cerebellum adaptation, see [`cluster/README.md`](../README.md).

Code used to cluster the Hi-C data and subsequent downstream analyses from Popay et al (https://www.nature.com/articles/s41588-026-02516-y)

Use the environment_mac.yml or environment_linux.yml files to create a conda environment with the required installations.
To perform clustering, Cluster 3.0 will also need to be installed independently: http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm
To utilize ChromHMM, this will also need to be install independently: https://compbio.mit.edu/ChromHMM/

Primary analysis can be performed using HiC_cluster3.ipynb (for clustering and cluster match) and grouped_loops_figures.ipynb

Example data is available here, with the exception of mcool files: https://drive.google.com/drive/folders/19EHKotAyxoccdB1J6QquAyIixMOqNJGu?usp=drive_link
