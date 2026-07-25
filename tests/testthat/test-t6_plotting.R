# test plotting functions.

data("sce_reducedDim", package = "Z6Vx8")

marker_dict <- list(
  "A" = c("gene_102", "gene_198", "gene_299"),
  "B" = c("gene_400", "gene_500")
)

# test plot_reduced_dimension -------------------------------------
# 1. Discrete colData variable
test_that("plot_reduced_dimension returns ggplot object with discrete colData variable", {

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = "celltype",
    color_by_label = "Cell Type",
    plot_title = "UMAP by Cell Type",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2"
  )

  expect_s3_class(pl, "ggplot")
  expect_equal(pl$labels$title, "UMAP by Cell Type")
  expect_equal(pl$labels$colour, "Cell Type")
  expect_equal(pl$labels$x, "UMAP1")
  expect_equal(pl$labels$y, "UMAP2")
})

# 2. Continuous colData variable
test_that("plot_reduced_dimension works with continuous colData variable", {

  sce_test <- sce_reducedDim
  SummarizedExperiment::colData(sce_test)$continuous_var <- runif(ncol(sce_test))

  pl <- plot_reduced_dimension(
    sce = sce_test,
    color_by = "continuous_var",
    color_by_label = "Continuous Variable",
    plot_title = "UMAP with continuous color",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2"
  )

  expect_s3_class(pl, "ggplot")
})


# 3. Test with gene expression
test_that("plot_reduced_dimension works with gene expression color", {

  gene_name <- rownames(sce_reducedDim)[1]

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = gene_name,
    color_by_label = gene_name,
    plot_title = paste("Expression of", gene_name),
    assay.type = "logcounts",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2"
  )

  expect_s3_class(pl, "ggplot")
})

# 4. Test with numeric vector as color_by
test_that("plot_reduced_dimension works with numeric vector as color_by", {

  numeric_vals <- runif(ncol(sce_reducedDim))

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = numeric_vals,
    color_by_label = "Random Values",
    plot_title = "UMAP with numeric vector",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2"
  )

  expect_s3_class(pl, "ggplot")
})


# 5. Test PCA
test_that("plot_reduced_dimension works with PCA reduced dimension", {

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = "celltype",
    color_by_label = "Cell Type",
    plot_title = "PCA by Cell Type",
    reduced_dim_name = "PCA",
    reduced_dim_1 = "PC1",
    reduced_dim_2 = "PC2"
  )

  expect_s3_class(pl, "ggplot")
})

# 6. Test show_label = TRUE with discrete variable
test_that("plot_reduced_dimension adds text labels when show_label = TRUE", {

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = "celltype",
    color_by_label = "Cell Type",
    plot_title = "UMAP with labels",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2",
    show_label = TRUE
  )

  expect_s3_class(pl, "ggplot")
})

# 7. Test show_label is ignored for continuous variable
test_that("plot_reduced_dimension ignores show_label for continuous variable", {

  sce_test <- sce_reducedDim
  SummarizedExperiment::colData(sce_test)$continuous_var <- runif(ncol(sce_test))

  pl <- plot_reduced_dimension(
    sce = sce_test,
    color_by = "continuous_var",
    color_by_label = "Continuous",
    plot_title = "Test",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2",
    show_label = TRUE  # Should be ignored for continuous
  )

  expect_s3_class(pl, "ggplot")
})


# 8. Test split_by
test_that("plot_reduced_dimension splits plot with split_by parameter", {

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = "celltype",
    color_by_label = "Cell Type",
    plot_title = "UMAP split by dataset",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2",
    split_by = "dataset_id"
  )

  expect_s3_class(pl, "ggplot")
})

# 9. Test rasterize = FALSE
test_that("plot_reduced_dimension uses geom_point when rasterize = FALSE", {

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = "celltype",
    color_by_label = "Cell Type",
    plot_title = "Test",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2",
    rasterize = FALSE
  )

  expect_s3_class(pl, "ggplot")
})

# 10. Test rasterize = TRUE
test_that("plot_reduced_dimension uses rasterized points when rasterize = TRUE", {

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = "celltype",
    color_by_label = "Cell Type",
    plot_title = "Test",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2",
    rasterize = TRUE,
    rastr_dpi = 300
  )

  expect_s3_class(pl, "ggplot")
  layer_classes <- sapply(pl$layers, function(l) class(l$geom)[1])
  expect_true(any(grepl("Rast|rast", layer_classes, ignore.case = TRUE)) ||
                "GeomPoint" %in% layer_classes,
              info = "Should have rasterized or point layer when rasterize=TRUE")
})


# 11. Test error on invalid color_by variable (not in colData or rownames)
test_that("plot_reduced_dimension errors on invalid color_by variable", {

  expect_error(
    plot_reduced_dimension(
      sce = sce_reducedDim,
      color_by = "nonexistent_variable",
      color_by_label = "Test",
      plot_title = "Test",
      reduced_dim_name = "UMAP",
      reduced_dim_1 = "UMAP1",
      reduced_dim_2 = "UMAP2"
    ),
    regexp = "invalid"
  )
})

# 12. Test error on invalid split_by variable
test_that("plot_reduced_dimension errors on invalid split_by variable", {

  expect_error(
    plot_reduced_dimension(
      sce = sce_reducedDim,
      color_by = "celltype",
      color_by_label = "Cell Type",
      plot_title = "Test",
      split_by = "nonexistent_column",
      reduced_dim_name = "UMAP",
      reduced_dim_1 = "UMAP1",
      reduced_dim_2 = "UMAP2"
    ),
    regexp = "invalid"
  )
})

# 13. Test error when gene specified but assay.type not in data
test_that("plot_reduced_dimension errors when assay.type not found", {

  gene_name <- rownames(sce_reducedDim)[1]

  expect_error(
    plot_reduced_dimension(
      sce = sce_reducedDim,
      color_by = gene_name,
      color_by_label = gene_name,
      plot_title = "Test",
      assay.type = "nonexistent_assay",
      reduced_dim_name = "UMAP",
      reduced_dim_1 = "UMAP1",
      reduced_dim_2 = "UMAP2"
    ),
    regexp = "assay.*not in data"
  )
})

# 14. Test multiple conditions: split_by with gene expression
test_that("plot_reduced_dimension handles split_by with gene expression coloring", {

  gene_name <- rownames(sce_reducedDim)[1]

  pl <- plot_reduced_dimension(
    sce = sce_reducedDim,
    color_by = gene_name,
    color_by_label = gene_name,
    plot_title = "Gene Expression by Dataset",
    split_by = "condition",
    reduced_dim_name = "UMAP",
    reduced_dim_1 = "UMAP1",
    reduced_dim_2 = "UMAP2"
  )

  expect_s3_class(pl, "ggplot")
  expect_s3_class(pl$facet, "FacetWrap")
  expect_equal(pl$labels$colour, gene_name)
})

# test plot_stackedViolin --------------------------------------------
# 1. test when row_ensembl = TRUE
test_that("plot_stackedViolin returns ggplot object with row_ensembl = TRUE", {

  pl <- plot_stackedViolin(
    sce = sce_reducedDim,
    marker_dict = marker_dict,
    cluster_lab = "celltype",
    row_ensembl = TRUE,
    plot_title = "Test Stacked Violin"
  )

  expect_s3_class(pl, "ggplot")
  expect_equal(pl$labels$title, "Test Stacked Violin")
})

# 2. Test with row_ensembl = FALSE (gene symbols as rownames)
test_that("plot_stackedViolin works with row_ensembl = FALSE", {

  sce_test <- sce_reducedDim
  rownames(sce_test) <- SummarizedExperiment::rowData(sce_test)$gene_symbol

  pl <- plot_stackedViolin(
    sce = sce_test,
    marker_dict = marker_dict,
    cluster_lab = "celltype",
    row_ensembl = FALSE,
    plot_title = "Test with gene symbols as rownames"
  )

  expect_s3_class(pl, "ggplot")
})

# 3. Test filtering of invalid genes
test_that("plot_stackedViolin filters out genes not in data", {

  marker_dict2 <- list(
    "A" = c("gene_299", "gene_2", "nonexistent_gene"),
    "B" = c("gene_3", "also_nonexistent")
  )

  pl <- plot_stackedViolin(
    sce = sce_reducedDim,
    marker_dict = marker_dict2,
    cluster_lab = "celltype",
    row_ensembl = TRUE
  )

  expect_s3_class(pl, "ggplot")
  n_features <- length(unique(pl$data$feature))
  expect_equal(n_features, 3)
})

# 4. Test error when logcounts assay is missing
test_that("plot_stackedViolin errors when logcounts assay is missing", {

  sce_no_logcounts <- sce_reducedDim
  SummarizedExperiment::assays(sce_no_logcounts) <- list(
    counts = SummarizedExperiment::assay(sce_reducedDim, "logcounts")
  )

  expect_error(
    plot_stackedViolin(
      sce = sce_no_logcounts,
      marker_dict = marker_dict,
      cluster_lab = "celltype"
    ),
    regexp = "logcounts not in sce"
  )
})

# 5. Test error when cluster_lab is invalid
test_that("plot_stackedViolin errors on invalid cluster_lab", {

  expect_error(
    plot_stackedViolin(
      sce = sce_reducedDim,
      marker_dict = marker_dict,
      cluster_lab = "nonexistent_column"
    ),
    regexp = "invalid `cluster_lab` variable"
  )
})


# 6. Test with different cluster_lab variables
test_that("plot_stackedViolin works with different cluster_lab variables", {

  pl <- plot_stackedViolin(
    sce = sce_reducedDim,
    marker_dict = marker_dict,
    cluster_lab = "condition"
  )

  expect_s3_class(pl, "ggplot")

  pl2 <- plot_stackedViolin(
    sce = sce_reducedDim,
    marker_dict = marker_dict,
    cluster_lab = "dataset_id"
  )

  expect_s3_class(pl2, "ggplot")
})

# 7. Test with single marker gene
test_that("plot_stackedViolin works with single gene", {

  marker_dict1 <- list("A" = c("gene_1"))

  pl <- plot_stackedViolin(
    sce = sce_reducedDim,
    marker_dict = marker_dict1,
    cluster_lab = "celltype"
  )

  expect_s3_class(pl, "ggplot")
  n_features <- length(unique(pl$data$feature))
  expect_equal(n_features, 1)
})

# test plot_marker_heatmap ------------------------------------------------
# 1. test with row_ensembl=TRUE
test_that("plot_marker_heatmap creates pheatmap with row_ensembl=TRUE", {

  pl <- plot_marker_heatmap(
    sce = sce_reducedDim,
    marker_dict = marker_dict,
    cluster_lab = "celltype",
    row_ensembl = TRUE,
    plot_title = "Test Heatmap"
  )

  expect_s3_class(pl, "pheatmap")
  expect_true(!is.null(pl$gtable))
})

# 2. heatmap with row_ensembl=FALSE
test_that("plot_marker_heatmap creates pheatmap with row_ensembl=FALSE", {
  sce_symbols <- sce_reducedDim
  rownames(sce_symbols) <- SummarizedExperiment::rowData(sce_symbols)$gene_symbol

  pl <- plot_marker_heatmap(
    sce = sce_symbols,
    marker_dict = marker_dict,
    cluster_lab = "celltype",
    row_ensembl = FALSE,
    plot_title = "Test Heatmap (symbols)"
  )

  expect_s3_class(pl, "pheatmap")
  expect_true(!is.null(pl$gtable))
})

# 3. heatmap with scale_expr=FALSE
test_that("plot_marker_heatmap works without scaling when scale_expr=FALSE", {

  pl <- plot_marker_heatmap(
    sce = sce_reducedDim,
    marker_dict = marker_dict,
    cluster_lab = "celltype",
    scale_expr = FALSE
  )

  expect_s3_class(pl, "pheatmap")
})

