# test pathway analysis functions

# skip if not installed
testthat::skip_if_not_installed("clusterProfiler")
testthat::skip_if_not_installed("ReactomePA")
testthat::skip_if_not_installed("org.Hs.eg.db")

# Test GO pathway processing
# 1. GO_BP pathways
test_that("process_pathway processes GO_BP pathways and returns TRUE", {

  test_output_dir <- tempdir()

  result <- process_pathway(
    pathway_output_path = test_output_dir,
    pathway_db_ver = "test_v1",
    pathway_str = "GO_BP",
    ensembl_version = "110"
  )

  expect_true(result)

  pathway_list <- readRDS(file.path(test_output_dir, "pathways_GOBP_test_v1.rds"))
  expect_type(pathway_list, "list")
  expect_true(length(pathway_list) > 0)
  expect_false(is.null(attr(pathway_list, "pathway_name")))
  expect_type(pathway_list[[1]], "character")

  unlink(file.path(test_output_dir, "pathways_GOBP_test_v1.rds"))
}) # slow connection


# Test Reactome
test_that("process_pathway processes Reactome pathways and returns TRUE", {

    test_output_dir <- tempdir()

  result <- process_pathway(
    pathway_output_path = test_output_dir,
    pathway_db_ver = "test_v1",
    pathway_str = "Reactome",
    ensembl_version = "110"
  )

  expect_true(result)
  pathway_list <- readRDS(file.path(test_output_dir, "pathways_Reactome_test_v1.rds"))
  expect_type(pathway_list, "list")
  expect_true(length(pathway_list) > 0)
  expect_false(is.null(attr(pathway_list, "pathway_name")))
  expect_type(pathway_list[[1]], "character")

  unlink(file.path(test_output_dir, "pathways_Reactome_test_v1.rds"))
})

# test KEGG
test_that("process_pathway processes KEGG pathways and returns TRUE", {

  test_output_dir <- tempdir()

  result <- process_pathway(
    pathway_output_path = test_output_dir,
    pathway_db_ver = "test_v1",
    pathway_str = "KEGG",
    ensembl_version = "110"
  )

  expect_true(result)
  pathway_list <- readRDS(file.path(test_output_dir, "pathways_KEGG_test_v1.rds"))
  expect_type(pathway_list, "list")
  expect_true(length(pathway_list) > 0)
  expect_false(is.null(attr(pathway_list, "pathway_name")))
  expect_type(pathway_list[[1]], "character")

  unlink(file.path(test_output_dir, "pathways_KEGG_test_v1.rds"))
})

# test SynGO

# test mapping ------------------------------------------
# 1. HGNC inputs
test_that("map_human_gene_to_mouse_gene_symbol returns data.frame with correct columns", {

  human_genes <- c("TP53", "BRCA1")

  result <- map_human_gene_to_mouse_gene_symbol(
    gene_list = human_genes,
    input_type = "gene_symbol",
    ensembl_version = "110"
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("input", "human_ensembl_gene_id", "mouse_gene_symbol"))
  expect_type(result$input, "character")
  expect_type(result$human_ensembl_gene_id, "character")
  expect_type(result$mouse_gene_symbol, "character")
})

# 2. Ensemb id
test_that("map_human_gene_to_mouse_gene_symbol works with Ensembl IDs",{

  # ENSG00000141510 = TP53, ENSG00000012048 = BRCA1
  human_ensembl <- c("ENSG00000141510", "ENSG00000012048")

  result <- map_human_gene_to_mouse_gene_symbol(
    gene_list = human_ensembl,
    input_type = "ensembl_id",
    ensembl_version = "110"
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(any(result$input %in% human_ensembl))
})

# 3. Entrez id
test_that("map_human_gene_to_mouse_gene_symbol works with Entrez IDs", {

  # 7157 = TP53, 672 = BRCA1
  human_entrez <- c("7157", "672")

  result <- map_human_gene_to_mouse_gene_symbol(
    gene_list = human_entrez,
    input_type = "entrez_id",
    ensembl_version = "110"
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(any(result$input %in% human_entrez))
})

# test plotting function -------------------------------
# helper function
#' create mock GSEA results data frame for testing
#' @param n_rows Number of pathways
#' @param all_significant If TRUE, all pathways have fdr < 0.05
#' @param seed Random seed for reproducibility
create_mock_gsea_df <- function(n_rows = 20, all_significant = FALSE, seed = 123) {

  set.seed(seed)

  if (all_significant) {
    fdr_vals <- runif(n_rows, 0.001, 0.04)
    pval_vals <- fdr_vals * runif(n_rows, 0.5, 0.9)
  } else {
    # Mix of significant and non-significant
    fdr_vals <- c(runif(n_rows %/% 2, 0.001, 0.04),
                  runif(n_rows - n_rows %/% 2, 0.06, 0.5))
    pval_vals <- fdr_vals * runif(n_rows, 0.5, 0.9)
  }

  data.frame(
    pathway_id = paste0("GO:", sprintf("%07d", 1:n_rows)),
    pathway_name = paste0("pathway ", 1:n_rows),
    pathway_type = rep(c("GO_BP", "GO_CC", "GO_MF", "Reactome"), length.out = n_rows),
    NES = rnorm(n_rows, mean = 0, sd = 2),
    pval = pval_vals,
    padj = fdr_vals * 1.1,
    fdr = fdr_vals,
    gene_ratio = runif(n_rows, 0.01, 0.3),
    setSize = sample(20:500, n_rows, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# 1. test gsea_dotplot function
test_that("gsea_dotplot returns a ggplot object on success", {

  df <- create_mock_gsea_df(n_rows = 20, all_significant = TRUE)

  pl <- gsea_dotplot(df, title_str = "test")

  expect_true(ggplot2::is_ggplot(pl))
})

# 2. test filtering
test_that("gsea_dotplot uses n_min when few pathways pass threshold", {

  df <- create_mock_gsea_df(n_rows = 15)
  # 3 pathways significant, but n_min=5 should show 5
  df$fdr <- c(rep(0.01, 3), rep(0.5, 12))
  df$pval <- df$fdr * 0.8

  pl <- gsea_dotplot(df, title_str = "test", n_min = 5, p_value_threshold = 0.05)
  expect_true(ggplot2::is_ggplot(pl))

  plot_data <- ggplot2::ggplot_build(pl)$data[[1]]
  expect_equal(nrow(plot_data), 5)
})
