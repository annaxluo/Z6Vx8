# test functions in `t1_clustering.R`

# test loading test data
test_that("test data can be loaded", {

  data("logcounts_mat", package = "Z6Vx8")
  expect_true(is.matrix(logcounts_mat) || inherits(logcounts_mat, "Matrix"))
  expect_true(nrow(logcounts_mat)==3973)
  expect_true(ncol(logcounts_mat)==5000)
})

# test compute_deScores function
test_that("compute_deScores returns correct structure", {
  data("logcounts_mat", package = "Z6Vx8")

  # create mock memberships
  n_cells <- ncol(logcounts_mat)
  membership <- sample(1:3, n_cells, replace = TRUE)

  result <- compute_deScores(
    membership = membership,
    expr_mat = logcounts_mat,
    min_cluster_size = 10,
    min_max_de_score = 5,
    min_mean_de_score = 1
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == 1)
  expected_cols <- c("mean.max.de.score.all", "mean.min.de.score.all",
                     "mean.mean.de.score.all", "mean.max.de.score.valid",
                     "mean.min.de.score.valid", "mean.mean.de.score.valid",
                     "n.valid.clust", "valid_clust")
  expect_true(all(expected_cols %in% colnames(result)))
})

# rest compute_deScores with single cluster
test_that("compute_deScores handles single cluster correctly", {
  data("logcounts_mat", package = "Z6Vx8")

  # one cluster
  n_cells <- ncol(logcounts_mat)
  membership <- rep(1, n_cells)

  result <- compute_deScores(
    membership = membership,
    expr_mat = logcounts_mat,
    min_cluster_size = 10,
    min_max_de_score = 5,
    min_mean_de_score = 1
  )

  expect_s3_class(result, "data.frame")
  expect_true(is.na(result$mean.max.de.score.all)) # shoudl return NA's
})

# test cluster_one_fold function
test_that("cluster_one_fold returns correct output structure", {

  testthat::skip_if_not_installed("Seurat")

  data("logcounts_mat", package = "Z6Vx8")

  result <- cluster_one_fold(
    expr_mat = logcounts_mat,
    leiden_resolution_list = c(0.01, 0.05, 0.10),
    init_pca_k = 50,
    pca_var_cutoff = 0.80,
    snn_type = "rank",
    snn_k = 10,
    used_cells = NULL,
    verbose = FALSE
  )

  expect_type(result, "list")
  expected_elements <- c("used_cells", "used_genes", "pca_k", "pca_comp",
                         "graphmat", "clust_leiden_opt", "DEscores", "clust_opt_df")
  expect_true(all(expected_elements %in% names(result)))

  # specific outputs
  expect_true(is.numeric(result$pca_k))
  expect_true(result$pca_k > 0 && result$pca_k <= 50)
  expect_s3_class(result$clust_opt_df, "data.frame")
  expect_equal(nrow(result$clust_opt_df), 3)
  expect_true("n.valid.clust.rate" %in% colnames(result$clust_opt_df))
})


# test assign_cluster_celltypes function
test_that("assign_cluster_celltypes runs and returns correct structure", {
  # mock cell-type score matrix (z-scores): 3 x 4
  set.seed(123)
  cL_ <- matrix(rnorm(12), nrow = 3, ncol = 4)
  rownames(cL_) <- paste0("cluster_", 1:3)
  colnames(cL_) <- c("A", "B", "C", "D")

  # make one cluster have a clear cell type signal
  cL_[1, "A"] <- 4  # high z-score for A in cluster 1

  result <- assign_cluster_celltypes(
    cL_ = cL_,
    p.adj.cutoff = 0.10,
    fold.diff = 3,
    p.adj.range = 0.10
  )

  expect_type(result, "character")
  expect_equal(length(result), 3)  # one assignment per cluster
})

