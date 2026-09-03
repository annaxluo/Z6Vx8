
# Z6Vx8

<!-- badges: start -->
[![R-CMD-check](https://github.com/annaxluo/Z6Vx8/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/annaxluo/Z6Vx8/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**Z6Vx8** includes utilities for analyzing single-nucleus RNA sequencing data in my project, including:

- Multi-fold clustering
- Consensus non-negative matrix (NMF) decomposition
- Differential expression tests (pseudobulk and GLMM-based methods)
- Gene set processing for pathway analysis
- Plotting functions

## Table of contents

- [Installation](#installation)
  - [Install Bioconductor](#1-install-bioconductor)
  - [Install core dependencies](#2-install-core-dependencies)
  - [Install developer version of RcppML](#3-install-developer-version-of-rcppml)
  - [Install optional dependencies for specific functions](#4-install-optional-dependencies-for-specific-functions)
  - [Install Z6Vx8](#install-z6vx8)
- [Using the consensus non-negative matrix factorization cNMF functions](#using-the-consensus-non-negative-matrix-factorization-cnmf-functions)
- [Requirements](#requirements)

## Installation

Z6Vx8 requires R 4.4.x.

The package uses several Bioconductor packages. Some dependencies are required for the core functions. Some are optional and only needed for specific functions.

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

## Using the consensus non-negative matrix factorization (cNMF) functions

The package contains utilities for consensus NMF (cNMF) decomposition, adapted from Kotliar et al. (2019) (https://github.com/dylkot/cNMF/blob/main/src/cnmf/cnmf.py) with the following modifications:

- Removing outlier components/clusters (based on silhouette score and cluster size).
- Computing consensus NMF as one closest to the cluster centroid to preserve sparseness of components.
- Removing cNMFs that capture mainly random effects (samples, batches, etc.) based on Wilcoxon tests.

The cNMF functions require the `RcppML` package. 

### To compute cNMF for a dense or sparse log2-normalized expression matrix: 

#### 1. Run NMF for multiple iterations
Use `run_NMF_iter` with a list of random seeds to run NMF for multiple interations: 
```R
data("logcounts_mat", package = "Z6Vx8") 

k_used <- 10
nmf_model_dir <- "/path/to/NMF/output/dir"
seed_list <- c(123, 456, 789)

nmf_w_list <- run_NMF_iter(
    expr_mat = logcounts_mat,
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list)
```

#### 2. (Optional) visualize prefiltering clustergram via `plot_clustergram`
Plot a heatmap showing clustering of the NMF components before filtering: 
```R
hm_pre_filt <- plot_clustergram(
    nmf_w_list, 
    used_k, 
    plot_title="pre-filter")
```

#### 3. Find mean distance to `L` nearest neighbors
Find average distance to `L` nearest neighbors, where `L` is controlled by 
the parameter `local_neighborhood_size`. This function will return a matrix of 
mean distance to `L` nearest neighbors in the first field, and a histogram of 
the mean distance values in the second field with a vertical cutoff line specified 
by `distance_histogram_threshold`: 
```R
filter_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.4)
    
# plot 
filter_result$pl_ 
```

#### 4. Filter components based on mean distance to nearest neighbors
Based on the histogram in Step 3, select and apply an appropriate threshold on 
mean distance to `L` nearest neighbors, and filter components.
```R
nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 0.5)
```

#### 5. (Optional) visualize filtered clustergram via `plot_clustergram`
```R
hm_post_filt <- plot_clustergram(
    nmf_w_list_filtered, 
    used_k, 
    plot_title="post-filter")
```

#### 6. Optimzie clustering parameters.
Next, optimize the parameters for `leiden` clustering on the filtered components. Parameters to optimize include the number of nearest neighbors (`snn_k_list`) and the weight scheme (`snn_type_list`) for the SNN graph, and the resolution for Leiden clustering (`leiden_resolution_list`). This function returns a list with four fields, including the optimization outputs (`clust_leiden_opt`), the silhouette width (`sil_df`), mean silhouette width (`sil_df_mean`), and a plot of the mean silhouette width (`mean_silhouette_plot`). Choose appropriate parameters according to the mean silhouette width. 
```R
opt_results <- cluster_opt(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k_list = c(3, 5, 7),
    snn_type_list = c("rank"),
    leiden_resolution_list = c(0.4, 0.7))
```

#### 7. Cluster NMF components using optimized parameters
Cluster the filtered NMF components using optimized parameters: 
```R
cluster_result <- cluster_components(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k = 3,
    snn_type = "rank",
    leiden_res = 0.7)
```

#### 8. Apply post-clustering filtering criteria to remove outlier components (based on silhouette and cluster size)
Filter clusters based on silhouette width (`sil_approx_thresh`) and cluster size (`cluster_size_thresh`) to remove low-quality clusters: 
```R
sil_filtered <- post_clustering_filter(
    cluster_result$sil.approx,
    sil_approx_thresh = 0,
    cluster_size_thresh = 2)
```

#### 9. Compute consensus NMF components: median, or prototypic component
If `cnmf_type="median"`, consensus NMFs (cNMFs) will be computed as the median of the components in each cluster. If `cnmf_type="prototype"`, the component with the smallest Euclidean distance to the centroid of a cluster will be used as the cNMF for that cluster.  
```R
cnmf_w <- consensus_nmf(
    nmf_w_list_filtered = nmf_w_list_filtered,
    sil.approx.filtered = sil_filtered,
    cnmf_type = "median")
```

#### 10. Compute usage of cNMFs.
Next, project the data onto the cNMF components. The returned model will have the cNMF weight matrix **W** (features × k), the coefficient matrix **H** (k × samples), and the diagonal **d**, such that: 
**A ≈ W · diag(d) · H**, where **A** is the input data matrix. 
```R
nmf_model <- consensus_nmf_usage(
    expr_mat = logcounts_mat,
    consensus_nmf = cnmf_w)
```

#### 11. Filter cNMFs that capture mainly random effects based on Wilcoxon tests.
Since some cNMF components will capture sample- or batch-specific effects, I apply a Wilcoxon test on the coefficients of each cNMF component to detect those samples: 
```R
# Create mock grouping (e.g., batch labels)
n_cells <- ncol(logcounts_mat)
grouping <- factor(rep(c("A", "B"), length.out = n_cells))

result <- filter_cnmf(
    nmf_model = nmf_model,
    grouping = grouping,
    fdr_cutoff = 0.10)
```
Where `result` is data.frame containing the Wilcoxon statistics. These cNMFs could be removed from downstream analysis as a correction for batch- and sample-specific effects. 


## Requirements

- R >= 4.4.0 and < 4.5.0
- Bioconductor 3.20 for R 4.4.x

## License

GPL (>= 3)

## Author

Anna Luo




