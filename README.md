
# Z6Vx8

<!-- badges: start -->
[![R-CMD-check](https://github.com/annaxluo/Z6Vx8/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/annaxluo/Z6Vx8/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

**Z6Vx8** provides utilities for analyzing single-nucleus RNA sequencing data, including:

- Multi-fold clustering
- Consensus NMF decomposition
- Differential expression tests (pseudobulk and GLMM-based methods)
- Gene set processing for pathway analysis
- Visualization functions

## Installation

### 1. Install Bioconductor Dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.18")  # For R 4.3.x

BiocManager::install(c(
    "AnnotationDbi",
    "BiocSingular",
    "biomaRt",
    "bluster",
    "clusterProfiler",
    "DESeq2",
    "edgeR",
    "limma",
    "org.Hs.eg.db",
    "ReactomePA",
    "scater",
    "scran",
    "SingleCellExperiment"
))
```

### 2. Install RcppML from GitHub

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("zdebruine/RcppML")
```

### 3. Install Z6Vx8

```r
remotes::install_github("annaxluo/Z6Vx8")
```

## Requirements

- R (>= 4.3.0)
- Bioconductor 3.18

## License

GPL (>= 3)

## Author

Anna Luo
