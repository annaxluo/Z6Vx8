# utilities for clustering
# 1. cluster_one_fold
# 2. a function to compute DEscore
# 3. a function to automatically assignment cell type based on Z test
# 4. a function to compute cluster similarity based on the Adjusted Rand Index

#' Cluster cells at multiple leiden resolutions, and returns an evaluation of the clustering
#' by computing a `deScore`.
#'
#' @param expr_mat gene-by-cell expression matrix, log2 normalized counts ("logcounts").
#' @param leiden_resolution_list numeric vector, a list of resolutions for leiden clustering.
#' @param init_pca_k A numeric, initial PCA dimension. Default is 50.
#' @param pca_var_cutoff A numeric, maximum variance explained in retained PCs. Default is 0.80.
#' @param snn_type A string, the type of weighting scheme for shared neighbors in kNN graphs.
#' @param snn_k An integer, the number of nearest neighbors to consider during graph construction.
#' @param used_cells A vector, list of indices for cells to use. If NULL, use all cells.
#' @param verbose A logical, whether to print progress messages. Default is FALSE.
#' @param ... Additional arguments passed to `compute_hvg` (e.g., `fdr_top_cutoff`, `n_top_genes`)
#' and `compute_deScores` (e.g., `min_cluster_size`, `min_max_de_score`, `min_mean_de_score`).
#'
#' @returns A list of clustering outputs.
#'
#' @importFrom magrittr %>%
#' @export
cluster_one_fold <- function(expr_mat,
                             leiden_resolution_list,
                             init_pca_k = 50,
                             pca_var_cutoff = 0.80,
                             snn_type = c("rank", "number", "jaccard"),
                             snn_k = 10,
                             used_cells = NULL,
                             verbose = FALSE,
                             ...){

  extra_args <- list(...)

  if (verbose) print("Computing HVGs.")

  if(is.null(used_cells)){
    used_cells = seq(ncol(expr_mat))
  }

  hvg_args <- list(
    expr_mat = expr_mat,
    used_cells = used_cells
  )
  hvg_param_names <- c("fdr_top_cutoff", "n_top_genes")
  hvg_args <- c(hvg_args, extra_args[names(extra_args) %in% hvg_param_names])

  used_genes <- do.call(compute_hvg, hvg_args)

  expr_mat_sub <- expr_mat[used_genes, used_cells]

  # PCA
  if (verbose) print("computing PCA..")
  # initial PCA
  pca_comp_0 <- scater::calculatePCA(x = expr_mat_sub,
                                     ncomponents = init_pca_k,
                                     ntop = nrow(expr_mat_sub),
                                     subset_row = NULL,
                                     scale = TRUE,
                                     transposed = FALSE,
                                     BSPARAM = BiocSingular::IrlbaParam())
  pca_var_df <- data.frame(
    x=seq(init_pca_k),
    y=attributes(pca_comp_0)[["varExplained"]]/sum(attributes(pca_comp_0)[["varExplained"]]))
  pca_var_df$cumsum <- cumsum(pca_var_df$y)
  pca_k <- max(pca_var_df$x[pca_var_df$cumsum <= pca_var_cutoff])

  pca_comp <- scater::calculatePCA(x = expr_mat_sub,
                                   ncomponents = pca_k,
                                   ntop = nrow(expr_mat_sub),
                                   subset_row = NULL,
                                   scale = TRUE,
                                   transposed = FALSE,
                                   BSPARAM = BiocSingular::IrlbaParam())

  # make tmp sce
  sce <- SingleCellExperiment::SingleCellExperiment(
    assay = list("logcounts"=expr_mat_sub)
  )
  SingleCellExperiment::reducedDim(sce, "PCA") <- pca_comp

  # clustering
  if (verbose) print("Clustering...")
  # construct SNN graph
  graphmat <- scran::buildSNNGraph(sce,
                                   k = snn_k,
                                   use.dimred = "PCA",
                                   type = snn_type)

  clust_leiden_opt <- lapply(leiden_resolution_list, function(res_){
    igraph::cluster_leiden(graphmat,
                           objective_function = "CPM",
                           resolution_parameter = res_,
                           n_iterations = 5)
  })
  names(clust_leiden_opt) <- paste0(snn_type, ".snn_k", snn_k, ".res", seq(length(clust_leiden_opt)))

  # compute DEscores
  if (verbose) print("Computing DE scores....")

  de_param_names <- c("min_cluster_size", "min_max_de_score", "min_mean_de_score")
  de_extra_args <- extra_args[names(extra_args) %in% de_param_names]

  if (verbose) print("Computing DE scores....")
  DEscores <- lapply(seq(length(clust_leiden_opt)), function(ii){
    de_args <- c(
      list(membership = clust_leiden_opt[[ii]]$membership, expr_mat = expr_mat_sub),
      de_extra_args
    )
    do.call(compute_deScores, de_args)
  })
  DEscores <- do.call("rbind.data.frame", DEscores)

  if (verbose) print("Return results......")
  clust_opt_df <- data.frame(
    clust_name = names(clust_leiden_opt),
    clust_resolution = leiden_resolution_list,
    nb_clusters = sapply(clust_leiden_opt, function(ret_){ ret_$nb_clusters })
  )
  clust_opt_df <- cbind(clust_opt_df, DEscores)
  clust_opt_df <- clust_opt_df %>%
    dplyr::mutate(n.valid.clust.rate = n.valid.clust / nb_clusters)

  # outputs
  outs <- list("used_cells"=used_cells,
               "used_genes"=used_genes,
               "pca_k"=pca_k,
               "pca_comp"=pca_comp,
               "graphmat"=graphmat,
               "clust_leiden_opt"=clust_leiden_opt,
               "DEscores"=DEscores,
               "clust_opt_df"=clust_opt_df)
  return(outs)
}


#' Compute the `deScore` between clusters for a list of cluster memberships using pairwise
#' Wilcoxon tests.
#'
#' @param membership A vector specifying cluster memberships.
#' @param expr_mat A matrix of gene expression, log2 normalized counts ("logcounts").
#' @param min_cluster_size An integer, minimum cluster size to be considered valid.
#' @param min_max_de_score A numeric, minimum maximum pairwise deScore for a cluster to be considered valid.
#' @param min_mean_de_score A numeric, minimum mean pairwise deScore for a cluster to be considered valid.
#'
#' @returns A data.frame summarizing mean deScores for the cluster memberships.
#'
#' @importFrom magrittr %>%
#' @export
compute_deScores <- function(membership,
                             expr_mat,
                             min_cluster_size = 50,
                             min_max_de_score = 5,
                             min_mean_de_score = 1){
  grp_ <- factor(membership)
  nb_clusters <- length(levels(grp_))
  # exclude clusters with few cells
  exclude_idx <- levels(grp_)[table(grp_) < min_cluster_size]

  if(nb_clusters > 1 & length(setdiff(levels(grp_), exclude_idx)) > 1){
    DE.res <- scran::pairwiseWilcox(expr_mat,
                                    groups = grp_,
                                    exclude = exclude_idx,
                                    log.p = FALSE)

    de.score <- sapply(DE.res$statistics, function(res_){
      fdr_ <- res_$FDR
      fdr_[fdr_ == 0] <- .Machine$double.xmin
      mean(-log10(fdr_))
    })

    res_df <- DE.res$pairs %>% as.data.frame() %>%
      dplyr::mutate(de.score = de.score)

    # max de.score for each cluster
    de.score.summ <- res_df %>%
      dplyr::group_by(first) %>%
      dplyr::summarize(max.de.score = max(de.score),
                       min.de.score = min(de.score),
                       mean.de.score = mean(de.score))

    # valid cluster
    valid_clust <- de.score.summ$first[de.score.summ$max.de.score >= min_max_de_score &
                                         de.score.summ$mean.de.score >= min_mean_de_score]

    # DEscore for the valid clusters
    res_df_valid <- res_df %>% dplyr::filter(first %in% valid_clust & second %in% valid_clust)
    de.score.summ.valid <- res_df_valid %>%
      dplyr::group_by(first) %>%
      dplyr::summarize(max.de.score = max(de.score),
                       min.de.score = min(de.score),
                       mean.de.score = mean(de.score))

    ret_df <- data.frame(
      mean.max.de.score.all = mean(de.score.summ$max.de.score),
      mean.min.de.score.all = mean(de.score.summ$min.de.score),
      mean.mean.de.score.all = mean(de.score.summ$mean.de.score),
      mean.max.de.score.valid = mean(de.score.summ.valid$max.de.score),
      mean.min.de.score.valid = mean(de.score.summ.valid$min.de.score),
      mean.mean.de.score.valid = mean(de.score.summ.valid$mean.de.score),
      n.valid.clust = length(valid_clust),
      valid_clust = paste(valid_clust, collapse=",")
    )
    return(ret_df)
  }

  ret_df <- data.frame(
    mean.max.de.score.all = NA,
    mean.min.de.score.all = NA,
    mean.mean.de.score.all = NA,
    mean.max.de.score.valid = NA,
    mean.min.de.score.valid = NA,
    mean.mean.de.score.valid = NA,
    n.valid.clust = NA,
    valid_clust = NA
  )
  return(ret_df)
}


#' Automatically assign cell type based on a Z test.
#'
#' @param cL_ A matrix. A matrix of cell-type scores with clusters in the row and
#' cell types in the column. Matrix should be row-standardized (z scores). Output
#' of `sctype_score` function in Ianevski et al. (2021). Nat Commun.
#' @param p.adj.cutoff A numeric. Threshold for detecting a single cell type for
#' a cluster.
#' @param fold.diff A numeric. The factor by which the second smallest adjusted
#' p value must be greater than the smallest adjusted p value, for the cell type
#' with the smallest adjusted p value to be used to annotate the cluster.
#' @param p.adj.range A numeric. Use cell types with adjusted p value within this
#' range beyond the smallest adjusted p to annotate the cluster.
#'
#' @returns a vector of characters specifying the assigned cell type for each cluster
#'
#' @export
assign_cluster_celltypes <- function(cL_,
                                     p.adj.cutoff=0.10,
                                     fold.diff=3,
                                     p.adj.range=0.10){
  # perform z test (one-sided)
  cell_type_assignment <- apply(cL_, 1, function(z){
    pval <- 1 - pnorm(z)
    padj <- p.adjust(pval, method = "fdr")

    if(sum(padj < p.adj.cutoff) == 1){
      return(names(padj)[which.min(padj)])
    }else{ # padj==0
      sorted_type <- names(z)[sort.int(z,index.return=T,decreasing=T)$ix]
      #return(paste0(first_type, ";", second_type))

      # if padj of the second type is at least `fold.diff` times that of the first type, return first type
      if(padj[[sorted_type[[2]]]] >= fold.diff*padj[[sorted_type[[1]]]]){
        return(sorted_type[[1]])
      }else{
        # return cell-types with padj no more than `p.adj.range` greater than the minimum
        used_types <- names(padj)[padj <= p.adj.range + min(padj)]
        return(paste(used_types, collapse = ","))
      }
    }
  })
  return(cell_type_assignment)
}

#' Compute adjusted Rand Index between two clusters to evalute cluster similarity.
#'
#' @param c1 A vector containing cluster labels for a clustering
#' @param c2 A vector containing cluster labels for another clustering
#' @returns a numeric between 0 and 1.
#'
#' @export
adj.rand.index <- function(c1, c2){
  n <- length(c1)
  if( length(c2) != n ) stop("Clusterings must be the same length.")
  t1  <- table(c1)
  t2  <- table(c2)
  t12 <- table(c1,c2)
  expected <- sum(choose(t1, 2)) * sum(choose(t2, 2)) / choose(n, 2)
  numerator <- sum(choose(t12, 2)) - expected
  denominator <- 0.5*(sum(choose(t1, 2)) + sum(choose(t2, 2))) - expected
  numerator / denominator
}
