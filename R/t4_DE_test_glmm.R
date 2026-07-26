# utilities for Differential Expression tests: GLMM-based tests
# references:
# 1. Pseudobulk
# https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R
# https://github.com/csoneson/conquer_comparison/tree/master/scripts
# 2. GLMM
# https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20Tweedie%20GLMM.Rmd

#' Check validity of arguments for DE_glmm function.
#'
#' @param sce A SingleCellExperiment object
#' @param group_colname A string specifying the grouping column name
#' @param random_var_colname A string specifying the random variable column name
#' @param normalization_method A string specifying the normalization method
#' @param test_method A string specifying the test method
#' @param extra_params A list of additional parameters passed via `...`
#'
#' @returns Logical.
#'
#' @keywords internal
check_args_glmm <- function(sce,
                            group_colname,
                            random_var_colname,
                            normalization_method,
                            test_method,
                            extra_params = list()){

  # 1. check that sce is a SingleCellExperiment
  if(!inherits(sce, "SingleCellExperiment")){
    stop("'sce' must be a SingleCellExperiment object.")
  }

  # 2. check that group_colname is in colData of sce
  sce_coldata_names <- colnames(SummarizedExperiment::colData(sce))
  if(!group_colname %in% sce_coldata_names){
    stop(paste0("'group_colname' (\"", group_colname, "\") not found in colData of 'sce'."))
  }

  # 3. check that random_var_colname is in colData of sce
  if(!random_var_colname %in% sce_coldata_names){
    stop(paste0("'random_var_colname' (\"", random_var_colname, "\") not found in colData of 'sce'."))
  }

  # 4. check normalization_method is valid
  valid_norm_methods <- c("RLE", "TMM", "colSums")
  normalization_method <- match.arg(normalization_method, valid_norm_methods)

  # 5. check test_method is valid
  valid_test_methods <- c("NBGMM", "NBLMM-reml", "PMM")
  test_method <- match.arg(test_method, valid_test_methods)

  # 6. check that extra_params arguments are valid
  valid_extra_args <- c(
    "cpc",
    "mincp",
    "ncore"
  )

  if(length(extra_params) > 0){
    extra_arg_names <- names(extra_params)
    invalid_args <- extra_arg_names[!extra_arg_names %in% valid_extra_args]
    if(length(invalid_args) > 0){
      stop(paste0("invalid argument(s) in '...': ", paste(invalid_args, collapse = ", "), "."))
    }
  }

  return(TRUE)
}


#' Perform DE test using GLMM methods
#'
#' @param sce A SingleCellExperiment object of the single-cell data
#' @param group_colname A string. The grouping variable
#' @param random_var_colname A string. Name of the column that specifies the random variable.
#' @param normalization_method A string. Normalization method to use. One of "RLE", "TMM", and "colSums".
#' @param test_method A string. Distribution model to use. One of "NBGMM", "NBLMM-reml", and "PMM".
#' @param ... Additional parameters, including `cpc`, `mincp`, and `ncore`.
#'
#' @returns A list.
#'
#' @importFrom magrittr %>%
#'
#' @export
DE_glmm <- function(sce,
                    group_colname,
                    random_var_colname,
                    normalization_method = c("RLE", "TMM", "colSums"),
                    test_method=c("NBGMM", "NBLMM-reml", "PMM"),
                    ...){

  params <- list(...)

  # check argument validity
  args_check <- check_args_glmm(
    sce = sce,
    group_colname = group_colname,
    random_var_colname = random_var_colname,
    normalization_method = normalization_method,
    test_method = test_method,
    extra_params = params
  )

  if(!args_check){
    stop("invalid inputs.")
  }

  if(!requireNamespace("nebula", quietly = TRUE)){
    stop(
      "Package 'nebula' is required for DE_glmm(). ",
      "Install it with remotes::install_github('lhe17/nebula').",
      call. = FALSE
    )
  }

  params[["normalization_method"]] <- normalization_method
  params[["test_method"]] <- test_method

  # change grouping variable name
  sce[[".grouping"]] <- sce[[group_colname]]

  # make design matrix
  design <- model.matrix(
    ~ .grouping,
    data=droplevels(SummarizedExperiment::colData(sce))
  )

  # compute normalization factor
  norm.factor <- Matrix::colSums(SingleCellExperiment::counts(sce))
  if(params[["normalization_method"]] %in% c("RLE", "TMM")){
    norm.factor <- norm.factor * edgeR::calcNormFactors(
      SingleCellExperiment::counts(sce),
      method=params[["normalization_method"]]
    )
  }# else: use colSums

  # fit model
  if(!"cpc" %in% names(params))
    params[["cpc"]] <- .005
  if(!"mincp" %in% names(params))
    params[["mincp"]] <- 20
  if(!"ncore" %in% names(params))
    params[["ncore"]] <- 1

  if(params[["test_method"]]=="NBGMM"){
    # fit model
    nebula_result <- nebula::nebula(
      count = SingleCellExperiment::counts(sce),
      id = SummarizedExperiment::colData(sce)[[random_var_colname]],
      offset = norm.factor,
      pred = design,
      cpc = params[["cpc"]],
      mincp = params[["mincp"]],
      ncore = params[["ncore"]]
    )
  }else if(params[["test_method"]]=="NBLMM-reml"){
    nebula_result <- nebula::nebula(
      count = SingleCellExperiment::counts(sce),
      id = SummarizedExperiment::colData(sce)[[random_var_colname]],
      model = 'NBLMM',
      reml = 1,
      offset = norm.factor,
      pred = design,
      cpc = params[["cpc"]],
      mincp = params[["mincp"]],
      ncore = params[["ncore"]]
    )
  }else{ # PMM
    nebula_result <- nebula::nebula(
      count = SingleCellExperiment::counts(sce),
      id = SummarizedExperiment::colData(sce)[[random_var_colname]],
      model = 'PMM',
      offset = norm.factor,
      pred = design,
      cpc = params[["cpc"]],
      mincp = params[["mincp"]],
      ncore = params[["ncore"]]
    )
  }

  # results table
  res_tbl <- nebula_result$summary
  rownames(res_tbl) <- res_tbl$gene

  contrast_name <- paste0("p_", colnames(design)[[2]])
  colnames(res_tbl)[colnames(res_tbl) == "gene_id"] <- "gene_id_test"

  res_tbl <- res_tbl %>%
    tibble::add_column(gene_id = rownames(.), .before = 1) %>%
    tibble::add_column(gene_symbol = SummarizedExperiment::rowData(sce)[rownames(.), "gene_symbol"], .after = 1) %>%
    dplyr::mutate(fdr = p.adjust(.data[[contrast_name]], "fdr"))

  res_tbl <- res_tbl[order(res_tbl[[contrast_name]]),]

  # updata column names
  names(res_tbl)[which(names(res_tbl) == paste0("logFC_", colnames(design)[[2]]))] <- "log2FoldChange"
  names(res_tbl)[which(names(res_tbl) == paste0("p_", colnames(design)[[2]]))] <- "pvalue"

  # return outputs
  outs <- list("params" = params,
               "design" = design,
               "res" = nebula_result,
               "res_tbl" = res_tbl)

  return(outs)
}
