# tests for `t0_basics.R`

test_that("size factor computed without throwing error", {
  set.seed(123)
  test_mat <- matrix(runif(100), nrow = 10, ncol = 10)
  expect_no_error(estimate_sf_sparse(test_mat))
})

test_that("test compute highly variable genes", {

  set.seed(123)
  data("logcounts_mat", package = "Z6Vx8")
  used_cells <- rep(TRUE, ncol(logcounts_mat))

  expect_no_error(compute_hvg(logcounts_mat,
                              used_cells,
                              fdr_top_cutoff=.5))
})
