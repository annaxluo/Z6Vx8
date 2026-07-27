# utilities for consensus NMF decomposition
#
# Adapted from Kotliar et al. (2019) (https://github.com/dylkot/cNMF/blob/main/src/cnmf/cnmf.py)
# with the following modifications:
# a. Removing outlier components/clusters (based on silhouette score, cluster size).
# b. Computing consensus NMF as one closest to the cluster centroid to preserve
#    sparseness of components.
# c. Removing cNMFs that capture mainly random effects (samples, batches, etc.)
#    based on Z test.
#
# Using "RcppML".
#
# Steps:
# 1. Run NMF for multiple iterations
# 2. (Optional) visualize prefiltering clustergram via `plot_clustergram`
# 3. Find mean distance to `L` nearest neighbors
# 4. Filter components based on mean distance to nearest neighbors
# 5. (Optional) visualize filtered clustergram via `plot_clustergram`
# 6. Optimzie clustering parameters.
# 7. Cluster NMF components using optimized parameters
# 8. Apply post-clustering filtering criteria to remove outlier components (based
#    on silhouette and cluster size)
# 9. Compute consensus NMF components: median, or prototypic component
# 10. Compute usage of consensus NMFs.
# 11. Filter cNMFs that capture mainly random effects based on Z tests.

# helper functions --------------------------------
#' Access fields from RcppML NMF models across RcppML versions
#'
#' RcppML models may be returned either as an S4 object with slots `@w`, `@h`,
#' and `@d`, or as a list-like object with fields `$w`, `$h`, and `$d`.
#'
#' @param model An RcppML NMF model object.
#' @param field Character scalar. One of `"w"`, `"h"`, or `"d"`.
#'
#' @returns The requested model field.
#'
#' @keywords internal
.get_RcppML_model_field <- function(model, field = c("w", "h", "d")){

  field <- match.arg(field)

  # newer RcppML development versions: list-like model
  if(is.list(model) && field %in% names(model)){
    return(model[[field]])
  }

  # older RcppML versions: "nmf" object with slots
  if(methods::is(model, "nmf") && field %in% methods::slotNames(model)){
    return(methods::slot(model, field))
  }

  stop(
    sprintf(
      "Could not access field '%s' from RcppML model. Expected either model$%s or model@%s.",
      field, field, field
    ),
    call. = FALSE
  )
}

#' Call RcppML::project across RcppML versions.
#'
#' Older RcppML versions used argument `data`; newer versions use argument `A`.
#'
#' @param w NMF weight matrix.
#' @param A Input matrix.
#'
#' @returns Projected usage/loadings matrix.
#'
#' @keywords internal
.rcppml_project <- function(w, A){

  project_formals <- names(formals(RcppML::project))

  if ("A" %in% project_formals) {
    return(RcppML::project(w = w, A = A))
  }

  if ("data" %in% project_formals) {
    return(RcppML::project(w = w, data = A))
  }

  stop(
    "Unsupported RcppML::project() API. Expected argument `A` or `data`.",
    call. = FALSE
  )
}

#' Call RcppML::mse across RcppML versions
#'
#' Older RcppML versions used argument `data`; newer versions use argument `A`.
#'
#' @param w NMF weight matrix.
#' @param d NMF diagonal/scaling vector.
#' @param h NMF usage/loadings matrix.
#' @param A Input matrix.
#'
#' @returns Mean squared error.
#'
#' @keywords internal
.rcppml_mse <- function(w, d, h, A){

  mse_formals <- names(formals(RcppML::mse))

  if("A" %in% mse_formals){
    return(RcppML::mse(w = w, d = d, h = h, A = A))
  }

  if("data" %in% mse_formals){
    return(RcppML::mse(w = w, d = d, h = h, data = A))
  }

  stop(
    "Unsupported RcppML::mse() API. Expected argument `A` or `data`.",
    call. = FALSE
  )
}


#' Run NMF decomposition for multiple iterations.
#'
#' @param expr_mat A dense or sparse matrix of features in rows and samples in
#' columns. Prefer `matrix` or `Matrix::dgCMatrix`. Log2 normalized counts.
#' @param k_used An integer. Rank used.
#' @param nmf_model_dir A string. Directory to store nmf models.
#' @param seed_list A vector of random seeds numeric.
#'
#' @returns Logical. TRUE if no error.
#'
#' @export
run_NMF_iter <- function(expr_mat,
                         k_used,
                         nmf_model_dir,
                         seed_list) {

  if(!requireNamespace("RcppML", quietly = TRUE)){
    stop("Package 'RcppML' is required.", call. = FALSE)
  }

  if(!requireNamespace("Matrix", quietly = TRUE)){
    stop("Package 'Matrix' is required.", call. = FALSE)
  }

  # for CMD check
  if (!"package:Matrix" %in% search()) {
    suppressPackageStartupMessages(
      library(Matrix)
    )
  }

  # input validation
  if(!is.matrix(expr_mat) && !inherits(expr_mat, "dgCMatrix")){
    stop("expr_mat must be a matrix or dgCMatrix.", call. = FALSE)
  }

  if(!is.numeric(k_used) || length(k_used) != 1 || k_used < 1){
    stop("k_used must be a positive integer.", call. = FALSE)
  }

  if(length(seed_list) < 1){
    stop("seed_list must contain at least one seed.", call. = FALSE)
  }

  dir.create(nmf_model_dir, recursive = TRUE, showWarnings = FALSE)

  for(seed_ in seed_list){

    out_fn <- file.path(nmf_model_dir,paste0("nmf-model-iter_k", k_used, "_seed", seed_, ".rds"))

    nmf_model <- tryCatch(
      {
        RcppML::nmf(
          expr_mat,
          k = k_used,
          tol = 1e-05,
          maxit = 100,
          verbose = FALSE,
          seed = seed_
        )
      },
      error = function(e) {
        stop(
          sprintf("NMF failed for seed %s: %s", seed_, conditionMessage(e)),
          call. = FALSE
        )
      }
    )

    saveRDS(nmf_model, out_fn)

    if(!file.exists(out_fn)){
      stop("NMF model was not written to file: ", out_fn, call. = FALSE)
    }

    message("model saved to: ", basename(out_fn))
  }

  TRUE
}

#' Load model weights and perform L2 normalization
#'
#' @param k_used An integer. Rank used.
#' @param nmf_model_dir A string. Directory to store nmf models.
#' @param seed_list A vector of random seeds (numeric)
#'
#' @returns A list of normalized NMF weights.
#'
#' @export
load_model_weights <- function(k_used,
                               nmf_model_dir,
                               seed_list){

  nmf_w_list <- do.call("cbind", lapply(seed_list, function(seed_){
    nmf_model <- readRDS(file.path(nmf_model_dir,
                                   paste0("nmf-model-iter_k", k_used, "_seed", seed_, ".rds")))

    nmf_w <- .get_RcppML_model_field(nmf_model, "w")

    # L2 normalization
    nmf_w_norm <- nmf_w %*% diag(1 / sqrt(colSums(nmf_w^2)))
    colnames(nmf_w_norm) <- paste0("s", seed_, "_nmf", seq(k_used))
    nmf_w_norm
  }))

  return(nmf_w_list)
}

#' Plot clustergram before or after filtering.
#'
#' @param nmf_w_list A list of normalized NMF weights.
#' @param used_k An integer. Rank used.
#' @param plot_title A string. Plot title
#'
#' @returns A pheatmap object.
#' @export
plot_clustergram <- function(nmf_w_list,
                             used_k,
                             plot_title){

  # plot clustergram before filtering.
  dist_mat <- dist(t(nmf_w_list), method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'ward.D2')

  # force k clusters
  cut_avg <- cutree(hclust_avg, k=used_k)
  nn_list <- lapply(unique(cut_avg), function(k_){
    names(cut_avg)[cut_avg == k_]
  })

  # aggregate distance matrix
  dist_mat2 <- as.matrix(dist_mat)

  cluster_dist_mat <- do.call("rbind", lapply(nn_list, function(cl1){
    rowMeans(do.call("rbind", lapply(nn_list, function(cl2){
      ret_ <- sapply(cl1, function(c_){
        idx_i <- which(colnames(dist_mat2) == c_)
        idx_j <- which(colnames(dist_mat2) %in% cl2)
        mean(dist_mat2[idx_i, idx_j])
      })
    })))
  }))

  colnames(cluster_dist_mat) <- seq(used_k)
  rownames(cluster_dist_mat) <- seq(used_k)
  pl <- pheatmap::pheatmap(cluster_dist_mat,
                           cluster_rows = T,
                           cluster_cols = T,
                           show_rownames = T,
                           show_colnames = T,
                           main = plot_title)

  return(pl)
}

#' Find average distance to `L` nearest neighbors, and filter components.
#'
#' @param nmf_w_list A list of normalized NMF weights.
#' @param used_k An integer. Rank used.
#' @param local_neighborhood_size A numeric, specifies the local neighborhood size,
#' and controls the number of nearest neighbors to consider.
#' @param distance_histogram_threshold A numeric. A vertical cutoff to plot in the
#' histogram of mean distance to `L` nearest neighbor.
#'
#' @returns A list of length 2, including `graph_mat_mean` and `histogram`.
#'
#' @export
mean_nn_distance <- function(nmf_w_list,
                             used_k,
                             local_neighborhood_size=0.3,
                             distance_histogram_threshold=0.3){

  L <- as.integer(local_neighborhood_size * (ncol(nmf_w_list) / used_k))

  # find NNs
  knn <- RANN::nn2(t(nmf_w_list), k = L + 1) # plus 1
  target <- as.integer(t(knn$nn.idx))
  source <- rep(seq_len(nrow(knn$nn.idx)), each = ncol(knn$nn.idx))

  knn_dist_mat <- as.numeric(t(knn$nn.dists))
  graph_mat <- cbind(source, target, knn_dist_mat)
  rm_ <- graph_mat[,1] == graph_mat[,2]
  graph_mat <- graph_mat[!rm_,] %>% as.data.frame()
  graph_mat_mean <- graph_mat %>%
    dplyr::group_by(source) %>%
    dplyr::summarize(knn_dist_mean = mean(knn_dist_mat))

  # histogram of distance distribution
  pl_ <- graph_mat_mean %>%
    ggplot2::ggplot(ggplot2::aes(x = knn_dist_mean)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            colour = "black",
                            fill = "white", bins = 100) +
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::geom_vline(xintercept = distance_histogram_threshold) +
    ggplot2::labs(title = paste0("Histogram: mean distance to ", L, " nearest neighbors"),
                  x = "mean distance",
                  y = "density") +
    ggplot2::theme_bw()

  outs <- list("graph_mat_mean"=graph_mat_mean,
               "histogram"=pl_)
  return(outs)
}


#' Apply threshold on mean distance to `L` nearest neighbors, and filter components.
#'
#' @param nmf_w_list A list of normalized NMF weights.
#' @param graph_mat_mean A matrix of mean distance to `L` nearest neighbors.
#' @param density_threshold A numeric. Threshold on mean distance to `L` NNs.
#'
#' @returns A list of filtered NMF weights.
#'
#' @export
filter_components <- function(nmf_w_list,
                              graph_mat_mean,
                              density_threshold=0.30){
  # filter components and recompute aggregate distance matrix
  graph_mat_mean <- graph_mat_mean %>%
    dplyr::mutate(is_valid = knn_dist_mean < density_threshold)

  used_comp_idx <- graph_mat_mean$source[graph_mat_mean$is_valid]

  nmf_w_list_filtered <- nmf_w_list[,used_comp_idx]

  return(nmf_w_list_filtered)
}

#' Optimizing clustering NMF components via `leiden` clustering
#' @param nmf_w_list_filtered A list of fitlered NMF components
#' @param snn_k_list A list of integers, specifying the number of nearest
#' neighbors to consider during SNN graph construction.
#' @param snn_type_list A list of string, specifying the type of weighting
#' scheme to use for shared neighbors.
#' @param leiden_resolution_list A list of numeric, specifying resolutions of
#' leiden clustering.
#'
#' @returns A list of length 4, including the optimization outputs (`clust_leiden_opt`),
#' the silhouette width (`sil_df`), mean silhouette width (`sil_df_mean`), and a
#' plot of the mean silhouette width (`mean_silhouette_plot`).
#'
#' @export
cluster_opt <- function(nmf_w_list_filtered,
                        snn_k_list=c(5, 10, 15, 20, 25),
                        snn_type_list=c("rank", "number", "jaccard"),
                        leiden_resolution_list=c(0.4, 0.7, 1.0)){

  # opt
  clust_leiden_opt <- do.call("c", lapply(snn_type_list, function(type_){
    do.call("c", lapply(snn_k_list, function(k_){
      graphmat <- scran::buildSNNGraph(nmf_w_list_filtered,
                                       k = k_,
                                       d = 10,
                                       type = type_) # optimize `k`, `type`
      ret_ <- lapply(leiden_resolution_list, function(res_){
        igraph::cluster_leiden(graphmat,
                               resolution_parameter = res_,
                               n_iterations = 5)
      })
      names(ret_) <- paste0(type_, ".", "res", leiden_resolution_list*10, ".snn_k", k_)
      ret_
    }))
  }))

  # evaluate: silhouette width
  sil_df <- lapply(clust_leiden_opt, function(clust_leiden){
    if(clust_leiden$nb_clusters == 1){
      n_samples <- ncol(nmf_w_list_filtered)
      return(S4Vectors::DataFrame(
        cluster = factor(clust_leiden$membership),
        other = factor(rep(NA, n_samples)),
        width = rep(NA_real_, n_samples)
      ))
    }

    bluster::approxSilhouette(x = t(nmf_w_list_filtered),
                              clusters = factor(clust_leiden$membership))
  })

  names(sil_df) <- names(clust_leiden_opt)

  # mean silhouette width
  sil_df_mean <- do.call("rbind.data.frame", mapply(function(df_, opt_name){
    data.frame(opt = opt_name,
               n_cluster = max(as.numeric(df_$cluster)),
               silhouette_mean = mean(df_$width), # keep NA for single cluster
               type = stringr::str_split_i(opt_name, "\\.", 1),
               snn_k = as.numeric(stringr::str_extract_all(opt_name, "\\d+")[[1]][[2]]),
               resolution = as.numeric(stringr::str_extract_all(opt_name, "\\d+")[[1]][[1]])/10)
  }, sil_df, names(sil_df), SIMPLIFY = F))

  # mean silhouette plot
  pl_ <- ggplot2::ggplot(sil_df_mean,
                         ggplot2::aes(x = snn_k,
                                      y = silhouette_mean,
                                      color = factor(resolution))) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(alpha = .5, linewidth = 1) +
    ggplot2::scale_x_continuous(name = "SNN graph: k") +
    ggplot2::scale_y_continuous(name = "mean silhouette width") +
    ggplot2::scale_color_discrete(name = "leiden resolution") +
    ggplot2::facet_wrap(~type) +
    ggplot2::ggtitle("leiden clustering opt: filtered NMF components") +
    ggplot2::theme_classic()

  outs <- list("clust_leiden_opt"=clust_leiden_opt,
               "clust_leiden_opt_silhouette"=sil_df,
               "clust_leiden_opt_mean_silhouette"=sil_df_mean,
               "mean_silhouette_plot"=pl_)

  return(outs)
}


#' Cluster NMF components using optimized parameters.
#'
#' @param nmf_w_list_filtered A list of filtered NMF weights.
#' @param snn_k An integer, specifying the number of nearest neighbors to
#' consider during graph construction.
#' @param snn_type A string, specifying the type of weighting scheme to use
#' for shared neighbors.
#' @param leiden_res A numeric, specifying the leiden resolution parameter to use.
#'
#' @returns A list of length 4, including the leiden clustering outputs (`clust_leiden`),
#' the observed and expected edge weights (`edge_weight_leiden`), the modularity
#' (`mod_leiden`), and silhouette width (`sil.approx`).
#'
#' @export
cluster_components <- function(nmf_w_list_filtered,
                               snn_k,
                               snn_type,
                               leiden_res){

  graphmat <- scran::buildSNNGraph(nmf_w_list_filtered,
                                   k = snn_k,
                                   d = 10,
                                   type = snn_type) # optimize `k`, `type`

  # leiden
  clust_leiden <- igraph::cluster_leiden(graphmat,
                                         resolution_parameter = leiden_res,
                                         n_iterations = 5)

  # evaluation metrics -------------------
  # observed and expected edge weights
  edge_weight_leiden <- bluster::pairwiseModularity(
    graphmat, factor(clust_leiden$membership),
    get.weights = TRUE
  )

  # modularity
  mod_leiden <- bluster::pairwiseModularity(
    graphmat, factor(clust_leiden$membership),
    get.weights = FALSE
  )

  # silhouette width
  sil.approx <- bluster::approxSilhouette(x = t(nmf_w_list_filtered),
                                          clusters = factor(clust_leiden$membership))

  outs <- list("clust_leiden"=clust_leiden,
               "edge_weight_leiden"=edge_weight_leiden,
               "mod_leiden"=mod_leiden,
               "sil.approx"=sil.approx)

  return(outs)
}


#' Apply post-clustering filtering on NMF components.
#'
#' @param sil.approx A data.frame of silhouette width.
#' @param sil_approx_thresh A numeric, cutoff for filtering clusters with silhouette
#' width lower than this value.
#' @param cluster_size_thresh An integer, minimum cluster size to keep.
#'
#' @returns A data.frame of stats (mean silhouette) of filtered clusters.
#'
#' @importFrom magrittr %>%
#' @export
post_clustering_filter <- function(sil.approx,
                                   sil_approx_thresh=0,
                                   cluster_size_thresh=50){

  sil.approx.filtered <- sil.approx %>% as.data.frame() %>%
    dplyr::mutate(silhouette_valid = width > sil_approx_thresh)

  valid_clust_idx <- which(table(sil.approx.filtered$cluster[sil.approx.filtered$silhouette_valid]) > cluster_size_thresh)

  sil.approx.filtered <- sil.approx.filtered %>%
    dplyr::mutate(cluster_valid = cluster %in% valid_clust_idx)

  sil.approx.filtered <- sil.approx.filtered %>% dplyr::filter(silhouette_valid & cluster_valid)
  sil.approx.filtered$cluster <- droplevels(sil.approx.filtered$cluster)
  sil.approx.filtered$other <- droplevels(sil.approx.filtered$other)

  return(sil.approx.filtered)
}

#' Compute consensus NMFs
#'
#' @param nmf_w_list_filtered A list of filtered NMF weights.
#' @param sil.approx.filtered A data.frame of mean silhouette scores for valid
#' clusters. Output of `post_clustering_filter`
#' @param cnmf_type A string specifying the type of cNMF to compute. One of `median`
#' (using median of components in a cluster), or `prototype` (using component
#' closest to the cluster centroid).
#'
#' @returns A matrix of consensus NMF weights, or NULL.
#'
#' @export
consensus_nmf <- function(nmf_w_list_filtered,
                          sil.approx.filtered,
                          cnmf_type=c("median", "prototype")){

  if(!cnmf_type %in% c("median", "prototype")){
    stop("invalid consensus NMF method.")
  }

  # median
  if(cnmf_type == "median"){
    consensus_nmf <- do.call("cbind", lapply(levels(sil.approx.filtered$cluster), function(cl_){
      used_comp_idx <- rownames(sil.approx.filtered)[sil.approx.filtered$cluster == cl_]
      apply(nmf_w_list_filtered[,used_comp_idx], 1, median)
    }))
    # rename cNMFs
    colnames(consensus_nmf) <- paste0("cnmf", seq(length(levels(sil.approx.filtered$cluster))))
    rownames(consensus_nmf) <- rownames(nmf_w_list_filtered)
  }

  if(cnmf_type == "prototype"){
    consensus_nmf <- do.call("cbind", lapply(levels(sil.approx.filtered$cluster), function(cl_){
      used_comp_idx <- rownames(sil.approx.filtered)[sil.approx.filtered$cluster == cl_]
      centroid_ <- apply(nmf_w_list_filtered[,used_comp_idx], 1, mean) # mean
      # Euclidean distance
      distance_ <- sqrt(colSums((sweep(nmf_w_list_filtered[,used_comp_idx], 1, centroid_, "-"))^2))
      prototype_idx <- used_comp_idx[which.min(distance_)]
      nmf_w_list_filtered[,prototype_idx]
    }))
    # rename cNMFs
    colnames(consensus_nmf) <- paste0("cnmf", seq(length(levels(sil.approx.filtered$cluster))))
    rownames(consensus_nmf) <- rownames(nmf_w_list_filtered)
  }

  return(consensus_nmf)
}


#' Compute consensus NMF usage
#'
#' @param expr_mat An expression matrix with genes in the row and cells in the column.
#' @param consensus_nmf A matrix of consensus NMF weights.
#'
#' @returns A list of length 3, including the weights (`cnmf_w`), usage (`cnmf_h`), and
#' the diagonal (`cnmf_d`) of the consensus NMF.
#'
#' @export
consensus_nmf_usage <- function(expr_mat,
                                consensus_nmf){

  if(!requireNamespace("RcppML", quietly = TRUE)){
    stop("Package 'RcppML' is required.", call. = FALSE)
  }

  if(!requireNamespace("Matrix", quietly = TRUE)){
    stop("Package 'Matrix' is required.", call. = FALSE)
  }

  # for CMD check
  if(!"package:Matrix" %in% search()){
    suppressPackageStartupMessages(
      library(Matrix)
    )
  }

  # input validation
  if(!is.matrix(expr_mat) && !inherits(expr_mat, "dgCMatrix")){
    stop("expr_mat must be a matrix or dgCMatrix.", call. = FALSE)
  }

  if(!is.matrix(consensus_nmf) && !inherits(consensus_nmf, "dgCMatrix")){
    stop("consensus_nmf must be a matrix or dgCMatrix.", call. = FALSE)
  }

  if(nrow(expr_mat) != nrow(consensus_nmf)){
    stop("expr_mat and consensus_nmf must have the same number of rows.", call. = FALSE)
  }

  ret_h <- .rcppml_project(
    w = consensus_nmf,
    A = expr_mat
  )

  cnmf_d <- rowSums(ret_h)
  names(cnmf_d) <- colnames(consensus_nmf)
  cnmf_h <- t(ret_h) %*% diag(1 / cnmf_d)
  cnmf_h <- t(cnmf_h)
  rownames(cnmf_h) <- colnames(consensus_nmf)

  error <- .rcppml_mse(
    w = consensus_nmf,
    d = cnmf_d,
    h = cnmf_h,
    A = expr_mat
  )

  outs <- list("cnmf_w"=consensus_nmf,
               "cnmf_d"=cnmf_d,
               "cnmf_h"=cnmf_h,
               "error"=error)

  return(outs)
}


#' Filter cNMFs that capture mainly random effects. Test random effect by Wilcoxon
#' tests.
#'
#' @param nmf_model A list, the NMF model object, include `cnmf_w`, `cnmf_d` and `cnmf_h`.
#' @param grouping A vector. A vector of string or factor, specifying the random factor
#' for the examples (column of the `h` matrix).
#' @param fdr_cutoff A numeric. The cutoff applied to FDR to detect outlier.
#'
#' @returns A list of data.frame. Each data.frame contains the Wilcoxon test statistics
#' for all cNMFs for a level in `grouping`.
#'
#' @export
filter_cnmf <- function(nmf_model,
                        grouping,
                        fdr_cutoff=.10){

  h <- nmf_model$cnmf_h
  grp_levels <- unique(sort(as.character(grouping)))

  ret_ <- lapply(grp_levels, function(g_){
    grp_ <- factor(grouping==g_, levels=c("TRUE", "FALSE"))
    DE.res <- scran::pairwiseWilcox(nmf_model$cnmf_h,
                                    groups = grp_,
                                    log.p = FALSE)
    df_ <- DE.res$statistics[[1]]
    df_$is.outlier <- df_$FDR <= fdr_cutoff
    df_
  })

  names(ret_) <- grp_levels
  return(ret_)
}







