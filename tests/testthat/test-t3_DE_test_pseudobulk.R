# test functions for pseudobulk differential expression tests

# test loading aggregated data
test_that("aggregated data can be loaded", {

  data("summed_sce", package = "Z6Vx8")
  expect_true(is(summed_sce, "SingleCellExperiment"))
  expect_true(nrow(summed_sce)==400)
  expect_true(ncol(summed_sce)==12)
})

# test filter_samples_pseudobulk function
test_that("test filtering samples for pseudobulk methods", {

  data("summed_sce", package = "Z6Vx8")

  current <- summed_sce[,summed_sce$celltype=="cell_type1"]

  # filter sampels
  sce_valid <- filter_samples_pseudobulk(current,
                                         group_by_var="condition",
                                         n_cells_dataset=39,
                                         n_datasets_celltype=2)

  # expect that one sample is filtered in each condition
  expect_true(sum(sce_valid$condition=="condition1")==3)
  expect_true(sum(sce_valid$condition=="condition2")==2)

  # expect all samples are filtered
  sce_valid2 <- filter_samples_pseudobulk(current,
                                          group_by_var="condition",
                                          n_cells_dataset=40,
                                          n_datasets_celltype=3)
  expect_null(sce_valid2)
})

# test DE_pseudobulk function ----------------------------------------------
# helper function to prepare test data
get_valid_sce <- function(){
  data("summed_sce", package = "Z6Vx8")
  # filter samples
  current <- summed_sce[,summed_sce$celltype=="cell_type1"]
  current_valid <- filter_samples_pseudobulk(current,
                                             group_by_var="condition",
                                             n_cells_dataset=10,
                                             n_datasets_celltype=2)

  return(current_valid)
}

#' helper function to check DE_pseudobulk output structure
#' @param result The output from DE_pseudobulk function
#' @param test_method The test method used (for informative error messages)
check_DE_pseudobulk_output <- function(result, test_method = NULL){
  # check return type
  expect_type(result, "list")

  # check required elements in output
  expected_elements <- c("pseudobulk_obj", "params", "design", "fit", "res", "res_tbl")
  expect_true(
    all(expected_elements %in% names(result)),
    info = paste0("Missing elements in output",
                  ifelse(!is.null(test_method), paste0(" for ", test_method), ""),
                  ": ", paste(setdiff(expected_elements, names(result)), collapse = ", "))
  )

  # check res_tbl structure
  expect_s3_class(result$res_tbl, "data.frame")
  expect_true(
    nrow(result$res_tbl) > 0,
    info = "res_tbl should have at least one row"
  )

  # check required columns in res_tbl
  expected_cols <- c("gene_id", "log2FoldChange", "pvalue", "fdr")
  expect_true(
    all(expected_cols %in% colnames(result$res_tbl)),
    info = paste0("Missing columns in res_tbl: ",
                  paste(setdiff(expected_cols, colnames(result$res_tbl)), collapse = ", "))
  )
}

# 1. test "DESeq2_pseudobulk" method with RLE (default) normalization
test_that("DE_pseudobulk works with DESeq2_pseudobulk method", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "RLE",
    test_method = "DESeq2_pseudobulk"
  )

  # common structure checks
  check_DE_pseudobulk_output(result, "DESeq2_pseudobulk")

  # method-specific checks
  expect_equal(result$params$test_method, "DESeq2_pseudobulk")
  expect_equal(result$params[["DESeq2.test"]], "LRT")
  expect_equal(result$params$normalization_method, "RLE")
})

# 2. test "DESeq2_pseudobulk" method with TMM normalization
test_that("DE_pseudobulk DESeq2 works with TMM normalization", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "TMM",
    test_method = "DESeq2_pseudobulk"
  )

  # common structure checks
  check_DE_pseudobulk_output(result, "DESeq2_pseudobulk")

  # method-specific checks
  expect_equal(result$params$test_method, "DESeq2_pseudobulk")
  expect_equal(result$params[["DESeq2.test"]], "LRT")
  expect_equal(result$params$normalization_method, "TMM")
})

# 3. test "DESeq2_pseudobulk" method with edgeR
test_that("DE_pseudobulk works with edgeR_pseudobulk method", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "RLE",
    test_method = "edgeR_pseudobulk"
  )

  # common structure checks
  check_DE_pseudobulk_output(result, "edgeR_pseudobulk")

  # method-specific checks
  expect_equal(result$params$test_method, "edgeR_pseudobulk")
  expect_equal(result$params$normalization_method, "RLE")
  expect_equal(result$params[["edgeR.test"]], "edgeRQLF") # default

  # check design matrix exists
  expect_true(is.matrix(result$design))
})

# 4. test "DESeq2_pseudobulk" method with edgeR using teh edgeRLRT test
test_that("DE_pseudobulk edgeR works with edgeRLRT test", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "TMM",
    test_method = "edgeR_pseudobulk",
    edgeR.test = "edgeRLRT"
  )

  # common structure checks
  check_DE_pseudobulk_output(result, "edgeR_pseudobulk")

  expect_equal(result$params$test_method, "edgeR_pseudobulk")
  expect_equal(result$params$normalization_method, "TMM")
  expect_equal(result$params[["edgeR.test"]], "edgeRLRT")

})

# 5. test "DESeq2_pseudobulk" method with limma
test_that("DE_pseudobulk works with limma_pseudobulk method", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "TMM",
    test_method = "limma_pseudobulk"
  )

  # common structure checks
  check_DE_pseudobulk_output(result, "limma_pseudobulk")

  # check params
  expect_equal(result$params$test_method, "limma_pseudobulk")
  expect_equal(result$params$normalization_method, "TMM")

  # check design matrix
  expect_true(is.matrix(result$design))

  # check fit object
  expect_true(!is.null(result$fit))
})


# 6. test gene filtering parameters
test_that("DE_pseudobulk accepts custom min.count and min.total.count", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "RLE",
    test_method = "edgeR_pseudobulk",
    min.count = 5,
    min.total.count = 10
  )

  expect_type(result, "list")
  expect_equal(result$params[["min.count"]], 5)
  expect_equal(result$params[["min.total.count"]], 10)
})

# 7. test invalid inputs: group_colname
test_that("DE_pseudobulk throws error for invalid group_colname", {
  current_valid <- get_valid_sce()

  expect_error(
    DE_pseudobulk(
      sce = current_valid,
      group_colname = "nonexistent_column",
      test_method = "edgeR_pseudobulk"
    ),
    "not found"
  )
})

# 8. test invalid inputs: test_method
test_that("DE_pseudobulk throws error for invalid test_method", {
  current_valid <- get_valid_sce()

  expect_error(
    DE_pseudobulk(
      sce = current_valid,
      group_colname = "condition",
      test_method = "invalid_method"
    )
  )
})

# 9. test output consistency across methods
test_that("all test methods return results with similar structure", {
  current_valid <- get_valid_sce()

  methods <- c("DESeq2_pseudobulk", "edgeR_pseudobulk", "limma_pseudobulk")

  results <- lapply(methods, function(method) {
    DE_pseudobulk(
      sce = current_valid,
      group_colname = "condition",
      normalization_method = "TMM",
      test_method = method
    )
  })
  names(results) <- methods

  # all should return lists
  expect_true(all(sapply(results, is.list)))

  # all res_tbl should have common columns
  common_cols <- c("gene_id", "gene_symbol", "log2FoldChange", "pvalue", "fdr")
  for (method in methods) {
    expect_true(
      all(common_cols %in% colnames(results[[method]]$res_tbl)),
      info = paste("Method:", method)
    )
  }

  # all res_tbl should have same number of genes after filtering
  n_genes <- sapply(results, function(r) nrow(r$res_tbl))
  expect_true(length(unique(n_genes)) == 1)
})

# test diagnostic plot ------------------------------------------
# 1. test diagnostic_plots with DESeq2_pseudobulk
test_that("diagnostic_plots works with DESeq2_pseudobulk output", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "RLE",
    test_method = "DESeq2_pseudobulk"
  )

  # test outputput
  tmp_dir <- tempdir()
  expect_no_error(diagnostic_plots(result, tmp_dir))
  expect_true(diagnostic_plots(result, tmp_dir))
  expect_true(file.exists(file.path(tmp_dir, "plot_DESeq2_DispEsts.pdf")))

  # cleanup
  unlink(file.path(tmp_dir, "plot_DESeq2_DispEsts.pdf"))
})


# 2. test diagnostic_plots with edgeR_pseudobulk
test_that("diagnostic_plots works with edgeR_pseudobulk output", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "RLE",
    test_method = "edgeR_pseudobulk"
  )

  # test output
  tmp_dir <- tempdir()
  expect_no_error(diagnostic_plots(result, tmp_dir))
  expect_true(diagnostic_plots(result, tmp_dir))
  expect_true(file.exists(file.path(tmp_dir, "plot_MD_samples.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "plot_MDS.pdf")))
  expect_true(file.exists(file.path(tmp_dir, "plot_diagnostics.pdf")))

  # cleanup
  unlink(file.path(tmp_dir, "plot_MD_samples.pdf"))
  unlink(file.path(tmp_dir, "plot_MDS.pdf"))
  unlink(file.path(tmp_dir, "plot_diagnostics.pdf"))
})

# 3. test diagnostic_plots with limma_pseudobulk
test_that("diagnostic_plots works with limma_pseudobulk output", {
  current_valid <- get_valid_sce()

  result <- DE_pseudobulk(
    sce = current_valid,
    group_colname = "condition",
    normalization_method = "RLE",
    test_method = "limma_pseudobulk"
  )

  tmp_dir <- tempdir()
  expect_no_error(diagnostic_plots(result, tmp_dir))
  expect_true(diagnostic_plots(result, tmp_dir))
  expect_true(file.exists(file.path(tmp_dir, "plot_diagnostics.pdf")))

  # cleanup
  unlink(file.path(tmp_dir, "plot_diagnostics.pdf"))
})

