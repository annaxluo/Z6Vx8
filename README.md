
# Z6Vx8

<!-- badges: start -->
[![R-CMD-check](https://github.com/annaxluo/Z6Vx8/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/annaxluo/Z6Vx8/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**Z6Vx8** provides utilities for analyzing single-nucleus RNA sequencing data, including:

- Multi-fold clustering
- Consensus non-negative matrix (NMF) decomposition
- Differential expression tests (pseudobulk and GLMM-based methods)
- Gene set processing for pathway analysis
- Plotting functions

## Installation

Z6Vx8 requires R 4.4.x.

The package uses several Bioconductor packages. Some dependencies are required for the core package functionality. Some are optional and only needed for specific functions.

### 1. Install Bioconductor

```r
if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install(version = "3.20") # For R 4.4.x
```

### 2. Install core dependencies

These packages are required for installing and using the main Z6Vx8 functionality.

```r
BiocManager::install(c(
  "AnnotationDbi",
  "BiocSingular",
  "biomaRt",
  "bluster",
  "DESeq2",
  "edgeR",
  "limma",
  "scater",
  "scran",
  "SingleCellExperiment",
  "SummarizedExperiment"
))

install.packages(c(
  "dplyr",
  "ggplot2",
  "ggrastr",
  "igraph",
  "magrittr",
  "Matrix",
  "pheatmap",
  "RANN",
  "RColorBrewer",
  "stringr",
  "tibble",
  "tidyr"
))
```

### 3. Install developer version of RcppML

Z6Vx8 uses `RcppML` for NMF decomposition. The developer versions '0.5.6' and '0.5.8' have been tested. 

```r
if(!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

devtools::install_github("zdebruine/RcppML")
```

### 4. Install optional dependencies for specific functions

Some Z6Vx8 functions require additional packages. These packages are listed in `Suggests` and are not required for installing the package, but they must be installed before running the corresponding functions.

#### Pathway analysis functions

The function `process_pathway` currently requires `clusterProfiler` to retrieve Gene Ontology (GO) and KEGG pathways, and `ReactomePA` to retrieve the Reactome pathways. Install these packages before using. In the future this will be replaced by internal functions to retrieve pathways to avoiding installing additional packages. When processing SynGO pathways, `process_pathway` requires the `org.Hs.eg.db` package to convert human gene symbol or Entrez IDs to Ensembl IDs. Install this package before using.  

- `clusterProfiler`
- `ReactomePA`
- `org.Hs.eg.db`

Install them with:

```r
BiocManager::install(c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db"
))
```

#### Seurat-related functions

The function `compute_hvg` that computes the highly-variable genes via Variance Stabilizing Transformation requires `FindVariableFeatures` from `Seurat`. In the future this will be replaced by an internal function to avoid installing additional packages. 

```r
install.packages("Seurat")
```

#### GLMM-based differential expression functions

GLMM-based differential expression function `DE_glmm` require:

- `nebula`

```r
devtools::install_github("lhe17/nebula")
```

Install this before running GLMM-based differential expression tests.

#### DESeq2 glmGamPoi fitType

When using the function `DE_pseudobulk`, dispersion estimates for DESeq2 can use the Gamma-Poisson distribution via the package `glmGamPoi`. Install this package before using `DESeq2.fitType = "glmGamPoi"`:

```r
BiocManager::install("glmGamPoi")
```

#### Plotting stacked violin plots with `theme_cowplot`

The function `plot_stackedViolin` requires the package `cowplot` to use `theme_cowplot`: 

```r
install.packages("cowplot")
```

### 5. Install Z6Vx8

```r
devtools::install_github("annaxluo/Z6Vx8")
```

## Requirements

- R >= 4.4.0 and < 4.5.0
- Bioconductor 3.20 for R 4.4.x

## License

GPL (>= 3)

## Author

Anna Luo




