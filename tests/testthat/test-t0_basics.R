# tests for `t0_basics.R`

test_that("size factor computed without throwing error", {
  set.seed(123)
  test_mat <- matrix(runif(100), nrow = 10, ncol = 10)
  expect_no_error(estimate_sf_sparse(test_mat))
})

test_that("compute highly variable genes without throwing error", {

  testthat::skip_if_not_installed("Seurat")

  set.seed(123)
  expr_mat <- Matrix::rsparsematrix(nrow = 3000, ncol = 100, density = 0.30)
  used_cells <- rep(TRUE, times = 100)

  expect_no_error(compute_hvg(expr_mat,
                              used_cells,
                              fdr_top_cutoff=.5))
})
