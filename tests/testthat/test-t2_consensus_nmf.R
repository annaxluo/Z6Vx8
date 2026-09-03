# test cNMF functions
# set up test data
data("logcounts_mat", package = "Z6Vx8")

# output dir
nmf_model_dir <- file.path(tempdir(), "nmf_models")
dir.create(nmf_model_dir, recursive = TRUE, showWarnings = FALSE)

# params
k_used <- 10
seed_list <- c(123, 456, 789)

# tests ------------------------------------------------------------
# 1. test run_NMF_iter
test_that("run_NMF_iter runs NMF decomposition for multiple iterations", {

  # Run NMF iterations
  result <- run_NMF_iter(
    expr_mat = logcounts_mat,
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  expect_true(result)

  # output files
  expected_files <- sapply(seed_list, function(s){
    file.path(nmf_model_dir, paste0("nmf-model-iter_k", k_used, "_seed", s, ".rds"))
  })

  for(fn in expected_files){
    expect_true(file.exists(fn), info = paste("file should exist:", fn))
  }

  # load one model
  model <- readRDS(expected_files[1])
  model_w <- .get_RcppML_model_field(model, "w")

  expect_true(is.matrix(model_w) || inherits(model_w, "Matrix"))
  expect_equal(ncol(model_w), k_used)
  expect_equal(nrow(model_w), nrow(logcounts_mat))
})


# 2. Test load_model_weights
test_that("load_model_weights loads and L2-normalizes NMF weights", {

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  expect_equal(nrow(nmf_w_list), nrow(logcounts_mat))
  expect_equal(ncol(nmf_w_list), k_used * length(seed_list))

  # L2 normalization: unit norm for each column
  col_norms <- sqrt(colSums(nmf_w_list^2))
  expect_true(all(abs(col_norms - 1) < 1e-6),
              info = "All columns should be L2-normalized to unit length")

  # check column names
  expect_true(all(grepl("^s\\d+_nmf\\d+$", colnames(nmf_w_list))),
              info = "Column names should follow pattern 's{seed}_nmf{k}'")
})

# 3. Test plot_clustergram
test_that("plot_clustergram creates a valid heatmap visualization", {

  # Load weights for testing
  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  pl <- plot_clustergram(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    plot_title = "Test Clustergram"
  )

  expect_s3_class(pl, "pheatmap")
  expect_true(!is.null(pl$tree_row), info = "Should have row dendrogram")
  expect_true(!is.null(pl$tree_col), info = "Should have column dendrogram")
})

# 4. Test mean_nn_distance
test_that("mean_nn_distance applies density threshold to filter components", {

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  # compute mean distances
  filter_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.4
  )

  expect_s3_class(filter_result[[1]], "tbl_df")
  expect_s3_class(filter_result[[2]], "ggplot")

})

# 5. Test filter_components
test_that("filter_components filters based on density threshold", {

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  # apply filter
  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 0.5
  )

  # Check output
  expect_true(is.matrix(nmf_w_list_filtered))
  expect_true(ncol(nmf_w_list_filtered) <= ncol(nmf_w_list),
              info = "Filtered matrix should have fewer or equal columns")
  expect_equal(nrow(nmf_w_list_filtered), nrow(nmf_w_list),
               info = "Row count should remain unchanged")
})

# 6. Test cluster_opt
test_that("cluster_opt optimizes clustering parameters", {

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 1.0 # keep all
  )

  # Run with reduced parameter space for faster testing
  result <- cluster_opt(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k_list = c(3, 5, 7),
    snn_type_list = c("rank"),
    leiden_resolution_list = c(0.4, 0.7)
  )

  expect_type(result, "list")
  expect_named(result, c("clust_leiden_opt", "clust_leiden_opt_silhouette",
                         "clust_leiden_opt_mean_silhouette", "mean_silhouette_plot"))

  expect_type(result$clust_leiden_opt, "list")
  expect_true(length(result$clust_leiden_opt) > 0)

  expect_type(result$clust_leiden_opt_silhouette, "list")
  expect_s3_class(result$clust_leiden_opt_mean_silhouette, "data.frame")

  expect_s3_class(result$mean_silhouette_plot, "ggplot")
})

# 7. Test cluster_components
test_that("cluster_components clusters with specified parameters", {
  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 1.0 # keep all
  )

  result <- cluster_components(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k = 3,
    snn_type = "rank",
    leiden_res = 0.7
  )

  expect_type(result, "list")
  expect_named(result, c("clust_leiden", "edge_weight_leiden", "mod_leiden", "sil.approx"))

  expect_true(is.matrix(result$mod_leiden))
  expect_s4_class(result$sil.approx, "DataFrame")
  expect_true("width" %in% names(result$sil.approx))
})

# 8. Test post_clustering_filter
test_that("post_clustering_filter returns correct structure", {

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 1.0 # keep all
  )

  cluster_result <- cluster_components(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k = 3,
    snn_type = "rank",
    leiden_res = 0.7
  )

  result <- post_clustering_filter(cluster_result$sil.approx,
                                   sil_approx_thresh = 0,
                                   cluster_size_thresh = 3) # filter some

  expect_s3_class(result, "data.frame")

  expect_true("silhouette_valid" %in% colnames(result))
  expect_true("cluster_valid" %in% colnames(result))
  expect_true("width" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))

  expect_true(all(result$width > 0))
  expect_true(all(result$silhouette_valid))
  expect_true(all(result$cluster_valid))

  expect_true(nrow(result) < nrow(cluster_result$sil.approx))
})

# 9. Test consensus_nmf
test_that("consensus_nmf with median method returns correct structure", {

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 1.0 # keep all
  )

  cluster_result <- cluster_components(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k = 3,
    snn_type = "rank",
    leiden_res = 0.7
  )

  sil_filtered <- post_clustering_filter(
    cluster_result$sil.approx,
    sil_approx_thresh = 0,
    cluster_size_thresh = 2  # keep all
  )

  result <- consensus_nmf(
    nmf_w_list_filtered = nmf_w_list_filtered,
    sil.approx.filtered = sil_filtered,
    cnmf_type = "median"
  )

  expect_true(is.matrix(result))
  n_clusters <- length(levels(sil_filtered$cluster))
  expect_equal(ncol(result), n_clusters)
  expect_equal(nrow(result), nrow(nmf_w_list_filtered))

  expect_true(all(grepl("^cnmf\\d+$", colnames(result))))
  expect_equal(colnames(result), paste0("cnmf", seq(n_clusters)))

  expect_equal(rownames(result), rownames(nmf_w_list_filtered))
})

# 10. Test consensus_nmf_usage
test_that("consensus_nmf_usage returns correct structure", {

  # load `logcounts_mat`

  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 1.0
  )

  cluster_result <- cluster_components(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k = 3,
    snn_type = "rank",
    leiden_res = 0.7
  )

  sil_filtered <- post_clustering_filter(
    cluster_result$sil.approx,
    sil_approx_thresh = 0,
    cluster_size_thresh = 2
  )

  cnmf_w <- consensus_nmf(
    nmf_w_list_filtered = nmf_w_list_filtered,
    sil.approx.filtered = sil_filtered,
    cnmf_type = "median"
  )

  result <- consensus_nmf_usage(
    expr_mat = logcounts_mat,
    consensus_nmf = cnmf_w
  )

  expect_type(result, "list")
  expect_named(result, c("cnmf_w", "cnmf_d", "cnmf_h", "error"))

  expect_equal(nrow(result$cnmf_w), nrow(logcounts_mat))
  expect_equal(ncol(result$cnmf_w), ncol(cnmf_w))

  expect_true(is.numeric(result$cnmf_d))
  expect_equal(length(result$cnmf_d), ncol(cnmf_w))
  expect_equal(names(result$cnmf_d), colnames(cnmf_w))

  expect_equal(nrow(result$cnmf_h), ncol(cnmf_w))
  expect_equal(ncol(result$cnmf_h), ncol(logcounts_mat))
  expect_equal(rownames(result$cnmf_h), colnames(cnmf_w))
})

# 11. Test filter_cnmf
test_that("filter_cnmf returns correct structure", {

  # load `logcounts_mat`
  nmf_w_list <- load_model_weights(
    k_used = k_used,
    nmf_model_dir = nmf_model_dir,
    seed_list = seed_list
  )

  distance_result <- mean_nn_distance(
    nmf_w_list = nmf_w_list,
    used_k = k_used,
    local_neighborhood_size = 2,
    distance_histogram_threshold = 0.5
  )

  nmf_w_list_filtered <- filter_components(
    nmf_w_list = nmf_w_list,
    graph_mat_mean = distance_result$graph_mat_mean,
    density_threshold = 1.0
  )

  cluster_result <- cluster_components(
    nmf_w_list_filtered = nmf_w_list_filtered,
    snn_k = 3,
    snn_type = "rank",
    leiden_res = 0.7
  )

  sil_filtered <- post_clustering_filter(
    cluster_result$sil.approx,
    sil_approx_thresh = 0,
    cluster_size_thresh = 2
  )

  cnmf_w <- consensus_nmf(
    nmf_w_list_filtered = nmf_w_list_filtered,
    sil.approx.filtered = sil_filtered,
    cnmf_type = "median"
  )

  nmf_model <- consensus_nmf_usage(
    expr_mat = logcounts_mat,
    consensus_nmf = cnmf_w
  )

  # Create mock grouping (e.g., batch labels)
  n_cells <- ncol(logcounts_mat)
  grouping <- factor(rep(c("A", "B"), length.out = n_cells))

  result <- filter_cnmf(
    nmf_model = nmf_model,
    grouping = grouping,
    fdr_cutoff = 0.10
  )

  expect_type(result, "list")

  expect_equal(length(result), 2)
  expect_true(all(c("A", "B") %in% names(result)))

  for(nm in names(result)){
    expect_s4_class(result[[nm]], "DFrame")
    expect_true("is.outlier" %in% colnames(result[[nm]]))
    expect_true("FDR" %in% colnames(result[[nm]]))
    expect_true("AUC" %in% colnames(result[[nm]]))

    expect_equal(nrow(result[[nm]]), ncol(cnmf_w))
  }
})

# clean up
test_that("cleanup temporary files", {

  files_to_remove <- list.files(nmf_model_dir, pattern = "nmf-model-iter", full.names = TRUE)
  if(length(files_to_remove) > 0){
    file.remove(files_to_remove)
  }
  expect_true(TRUE)
})



