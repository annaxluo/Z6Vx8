# test functions for GLMM-based differential expression tests

# skip if "nebula" is not installed.

testthat::skip_if_not_installed("nebula")

# test invalid inputs to DE_glmm -------------------------------------------
# 1. group_colname does not exist in colData
test_that("DE_glmm errors when group_colname not in colData", {

  data("sce", package = "Z6Vx8")

  expect_error(
    DE_glmm(
      sce = sce,
      group_colname = "nonexistent_column",
      random_var_colname = "dataset_id",
      normalization_method = "RLE",
      test_method = "NBGMM"
    ),
    "'group_colname'.*not found in colData"
  )
})

# 2. random_var_colname does not exist in colData
test_that("DE_glmm errors when random_var_colname not in colData", {

  data("sce", package = "Z6Vx8")

  expect_error(
    DE_glmm(
      sce = sce,
      group_colname = "condition",
      random_var_colname = "nonexistent_random",
      normalization_method = "RLE",
      test_method = "NBGMM"
    ),
    "'random_var_colname'.*not found in colData"
  )
})


# test DE_glmm function -------------------------------------------------------
# helper function to prepare test data
get_valid_sce <- function(n_genes = 100){
  data("sce", package = "Z6Vx8")

  # filter cell type
  dat_sub <- sce[,sce[["celltype"]]=="cell_type1"]
  if(n_genes < nrow(dat_sub)){
    dat_sub <- dat_sub[1:n_genes, ]
  }

  # sort by dataset_id to ensure cells are grouped by subject
  cell_order <- order(dat_sub[["dataset_id"]])
  dat_sub <- dat_sub[, cell_order]
  return(dat_sub)
}

#' Helper function to test output structures.
#' @param result The output from DE_glmm function
#' @param test_method The test method used
check_DE_glmm_output <- function(result, test_method = NULL){

  expect_type(result, "list")

  # check required elements in output
  expected_elements <- c("params", "design", "res", "res_tbl")
  expect_true(
    all(expected_elements %in% names(result)),
    info = paste0("Missing elements in output",
                  ifelse(!is.null(test_method), paste0(" for ", test_method), ""),
                  ": ", paste(setdiff(expected_elements, names(result)), collapse = ", "))
  )

  # check design is a matrix
  expect_true(is.matrix(result$design))

  # check res_tbl structure
  expect_s3_class(result$res_tbl, "data.frame")
  expect_true(
    nrow(result$res_tbl) > 0,
    info = "res_tbl should have at least one row"
  )

  # check required columns in res_tbl
  expected_cols <- c("gene_id", "gene_symbol", "log2FoldChange", "pvalue", "fdr")
  expect_true(
    all(expected_cols %in% colnames(result$res_tbl)),
    info = paste0("Missing columns in res_tbl: ",
                  paste(setdiff(expected_cols, colnames(result$res_tbl)), collapse = ", "))
  )
}

# 1. test NBGMM method with RLE normalization
test_that("DE_glmm works with NBGMM method and RLE normalization", {

  sce_valid <- get_valid_sce()

  result <- DE_glmm(
    sce = sce_valid,
    group_colname = "condition",
    random_var_colname = "dataset_id",
    normalization_method = "RLE",
    test_method = "NBGMM",
    ncore = 1
  )

  # common structure checks
  check_DE_glmm_output(result, "NBGMM")

  # method-specific checks
  expect_equal(result$params$test_method, "NBGMM")
  expect_equal(result$params$normalization_method, "RLE")
})

# 2. test NBGMM method with TMM normalization
test_that("DE_glmm works with NBGMM method and TMM normalization", {

  sce_valid <- get_valid_sce()

  result <- DE_glmm(
    sce = sce_valid,
    group_colname = "condition",
    random_var_colname = "dataset_id",
    normalization_method = "TMM",
    test_method = "NBGMM",
    ncore = 1
  )

  # common structure checks
  check_DE_glmm_output(result, "NBGMM")

  # method-specific checks
  expect_equal(result$params$test_method, "NBGMM")
  expect_equal(result$params$normalization_method, "TMM")
})

# 3. test NBGMM method with colSums normalization
test_that("DE_glmm works with NBGMM method and colSums normalization", {

  sce_valid <- get_valid_sce()

  result <- DE_glmm(
    sce = sce_valid,
    group_colname = "condition",
    random_var_colname = "dataset_id",
    normalization_method = "colSums",
    test_method = "NBGMM",
    ncore = 1
  )

  # common structure checks
  check_DE_glmm_output(result, "NBGMM")

  # method-specific checks
  expect_equal(result$params$test_method, "NBGMM")
  expect_equal(result$params$normalization_method, "colSums")
})

# 4. test NBLMM-reml method
test_that("DE_glmm works with NBLMM-reml method", {

  sce_valid <- get_valid_sce()

  result <- DE_glmm(
    sce = sce_valid,
    group_colname = "condition",
    random_var_colname = "dataset_id",
    normalization_method = "RLE",
    test_method = "NBLMM-reml",
    ncore = 1
  )

  # common structure checks
  check_DE_glmm_output(result, "NBLMM-reml")

  # method-specific checks
  expect_equal(result$params$test_method, "NBLMM-reml")
})

# 5. test PMM method
test_that("DE_glmm works with PMM method", {

  sce_valid <- get_valid_sce()

  result <- DE_glmm(
    sce = sce_valid,
    group_colname = "condition",
    random_var_colname = "dataset_id",
    normalization_method = "RLE",
    test_method = "PMM",
    ncore = 1
  )

  # common structure checks
  check_DE_glmm_output(result, "PMM")

  # method-specific checks
  expect_equal(result$params$test_method, "PMM")
})

# 6. test custom parameters (cpc, mincp)
test_that("DE_glmm respects custom cpc and mincp parameters", {

  sce_valid <- get_valid_sce()

  result <- DE_glmm(
    sce = sce_valid,
    group_colname = "condition",
    random_var_colname = "dataset_id",
    normalization_method = "RLE",
    test_method = "NBGMM",
    cpc = 0.01,
    mincp = 10,
    ncore = 1
  )

  # common structure checks
  check_DE_glmm_output(result, "NBGMM")

  expect_equal(result$params$cpc, 0.01)
  expect_equal(result$params$mincp, 10)
})

