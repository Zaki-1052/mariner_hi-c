# R Package Dependencies for Mariner Pipeline

## Summary

This document lists all R packages required by the mariner differential chromatin loop analysis pipeline.

**Total packages: 36**
- CRAN: 12 packages
- Bioconductor: 24 packages

---

## Installation

Run the automated installation script:

```bash
Rscript install_R_packages.R
```

Or install manually:

```r
# Install BiocManager
install.packages("BiocManager")

# Install all packages
BiocManager::install(c(
  # CRAN packages
  "yaml", "ggplot2", "dplyr", "tibble", "tidyr", "pheatmap",
  "patchwork", "RColorBrewer", "viridis", "VennDiagram", "Matrix", "scales",

  # Bioconductor packages
  "mariner", "InteractionSet", "GenomicRanges", "strawr", "DelayedArray",
  "HDF5Array", "edgeR", "rtracklayer", "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db", "GenomicFeatures", "ChIPseeker", "clusterProfiler",
  "enrichplot", "DOSE", "SummarizedExperiment", "EnhancedVolcano"
))
```

---

## Package List by Category

### 1. Configuration & Data Manipulation (CRAN)

| Package | Purpose | Used In |
|---------|---------|---------|
| yaml | Configuration file parsing | edgeR.R |
| dplyr | Data manipulation | edgeR.R, visualizations.R, apa_analysis.R |
| tibble | Modern data frames | edgeR.R |
| tidyr | Data tidying | compare_resolutions.R, apa_analysis.R |
| Matrix | Sparse matrix operations | aggregate.R |

### 2. Core Visualization (CRAN)

| Package | Purpose | Used In |
|---------|---------|---------|
| ggplot2 | Core plotting framework | edgeR.R, compare_resolutions.R, downstream_analysis.R, visualizations.R, apa_analysis.R |
| patchwork | Plot composition | compare_resolutions.R, downstream_analysis.R, visualizations.R |
| pheatmap | Heatmap plotting | compare_resolutions.R, visualizations.R, apa_analysis.R |
| RColorBrewer | Color palettes | visualizations.R, apa_analysis.R |
| viridis | Perceptually uniform colors | apa_analysis.R |
| VennDiagram | Venn diagram generation | compare_resolutions.R |
| scales | Scale functions for ggplot2 | visualizations.R (via `scales::`) |

### 3. Bioconductor Core Infrastructure

| Package | Purpose | Used In |
|---------|---------|---------|
| GenomicRanges | Genomic interval operations | All scripts |
| InteractionSet | Hi-C interaction data structures | All scripts |
| SummarizedExperiment | Experiment data containers | apa_analysis.R |
| DelayedArray | Out-of-memory array operations | extract_counts.R, aggregate.R, apa_analysis.R |
| HDF5Array | HDF5-backed arrays | extract_counts.R, aggregate.R, apa_analysis.R |

### 4. Hi-C Analysis (Bioconductor)

| Package | Purpose | Used In |
|---------|---------|---------|
| mariner | Hi-C loop analysis framework | prep_loops.R, extract_counts.R, aggregate.R, downstream_analysis.R, visualizations.R, apa_analysis.R |
| strawr | .hic file I/O (Straw API) | extract_counts.R, visualizations.R, apa_analysis.R |

### 5. Differential Analysis (Bioconductor)

| Package | Purpose | Used In |
|---------|---------|---------|
| edgeR | Differential loop analysis (GLM) | aggregate.R, edgeR.R |
| limma | Linear models (edgeR dependency) | (Indirect) |

### 6. Genome Annotation (Bioconductor)

| Package | Purpose | Used In |
|---------|---------|---------|
| rtracklayer | Import/export genomic data | downstream_analysis.R |
| GenomicFeatures | Gene/transcript models | downstream_analysis.R |
| TxDb.Mmusculus.UCSC.mm10.knownGene | Mouse mm10 gene annotations | downstream_analysis.R, visualizations.R |
| org.Mm.eg.db | Mouse gene ID mappings | downstream_analysis.R, visualizations.R |
| AnnotationDbi | Annotation database interface | (Dependency) |

### 7. Functional Enrichment (Bioconductor)

| Package | Purpose | Used In |
|---------|---------|---------|
| ChIPseeker | ChIP-seq peak annotation | visualizations.R |
| clusterProfiler | GO/KEGG enrichment analysis | visualizations.R |
| enrichplot | Enrichment visualization | visualizations.R |
| DOSE | Disease ontology enrichment | visualizations.R |
| GO.db | Gene Ontology database | (Dependency) |
| KEGG.db | KEGG pathway database | (Dependency) |

### 8. Enhanced Visualization (Bioconductor)

| Package | Purpose | Used In |
|---------|---------|---------|
| EnhancedVolcano | Publication-quality volcano plots | visualizations.R |

---

## Script-Specific Dependencies

### prep_loops.R
- mariner
- InteractionSet
- GenomicRanges

### extract_counts.R
- mariner
- InteractionSet
- strawr
- DelayedArray
- HDF5Array

### aggregate.R
- mariner
- InteractionSet
- HDF5Array
- DelayedArray
- Matrix
- edgeR

### edgeR.R
- edgeR
- yaml
- GenomicRanges
- InteractionSet
- ggplot2
- dplyr
- tibble

### compare_resolutions.R
- tidyverse (ggplot2, dplyr, tidyr, etc.)
- GenomicRanges
- InteractionSet
- VennDiagram
- pheatmap
- patchwork

### downstream_analysis.R
- mariner
- InteractionSet
- GenomicRanges
- rtracklayer
- tidyverse (ggplot2, dplyr, tidyr)
- patchwork
- TxDb.Mmusculus.UCSC.mm10.knownGene
- org.Mm.eg.db
- GenomicFeatures

### visualizations.R
- ggplot2
- patchwork
- pheatmap
- RColorBrewer
- EnhancedVolcano
- GenomicRanges
- InteractionSet
- ChIPseeker
- TxDb.Mmusculus.UCSC.mm10.knownGene
- org.Mm.eg.db
- clusterProfiler
- enrichplot
- DOSE
- mariner (optional, for APA)
- strawr (optional, for APA)
- tidyverse

### apa_analysis.R
- mariner
- InteractionSet
- strawr
- HDF5Array
- SummarizedExperiment
- ggplot2
- pheatmap
- RColorBrewer
- viridis
- dplyr
- tidyr

---

## Version Requirements

- **R version**: ≥ 4.3.0 (for Bioconductor 3.18+)
- **Bioconductor**: 3.18 (R 4.3.x) or 3.19 (R 4.4.x)

The installation script automatically selects the appropriate Bioconductor version based on your R version.

---

## Troubleshooting

### Common Issues

**1. HDF5 library not found**
```bash
# macOS (Homebrew)
brew install hdf5

# Ubuntu/Debian
sudo apt-get install libhdf5-dev

# Then reinstall HDF5Array
BiocManager::install("HDF5Array", force = TRUE)
```

**2. Cairo graphics library missing (for PDF output)**
```bash
# macOS
brew install cairo

# Ubuntu/Debian
sudo apt-get install libcairo2-dev
```

**3. XML library missing (for enrichment analysis)**
```bash
# macOS
brew install libxml2

# Ubuntu/Debian
sudo apt-get install libxml2-dev
```

**4. strawr compilation issues**
```bash
# Ensure you have a C++ compiler
# macOS: Install Xcode Command Line Tools
xcode-select --install

# Ubuntu/Debian
sudo apt-get install build-essential
```

**5. Memory issues with large datasets**
- Increase R memory limit: `memory.limit(size = 100000)` (Windows)
- Use HDF5-backed arrays (already implemented in pipeline)
- Run on HPC with more RAM

---

## Testing Installation

After installation, test that all packages load correctly:

```r
# Test script
test_packages <- c(
  "mariner", "InteractionSet", "GenomicRanges", "strawr",
  "HDF5Array", "edgeR", "ggplot2", "ChIPseeker",
  "clusterProfiler", "EnhancedVolcano"
)

for (pkg in test_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("✓ %s: OK\n", pkg))
  } else {
    cat(sprintf("✗ %s: FAILED\n", pkg))
  }
}
```

---

## Additional Resources

- **Bioconductor**: https://bioconductor.org/
- **mariner documentation**: https://bioconductor.org/packages/mariner
- **edgeR User's Guide**: https://bioconductor.org/packages/edgeR
- **strawr**: https://github.com/aidenlab/straw/tree/master/R
- **EnhancedVolcano**: https://bioconductor.org/packages/EnhancedVolcano

---

## Citation

If you use this pipeline, please cite the relevant packages:

**mariner:**
- Kramer et al. (2022) mariner: Explore the Hi-Cs. Bioconductor.

**edgeR:**
- Robinson MD, McCarthy DJ, Smyth GK (2010). "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics, 26(1), 139-140.

**clusterProfiler:**
- Yu G, Wang L, Han Y, He Q (2012). "clusterProfiler: an R package for comparing biological themes among gene clusters." OMICS: A Journal of Integrative Biology, 16(5), 284-287.

**strawr:**
- Straw paper (Durand et al. 2016, Cell Systems)

---

Last updated: 2025-11-10
