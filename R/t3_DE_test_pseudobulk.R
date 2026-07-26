# utilities for Differential Expression tests: pseudobulk tests
# references:
# 1. Pseudobulk
# https://github.com/neurorestore/DE-analysis/blob/master/R/functions/run_DE.R
# https://github.com/csoneson/conquer_comparison/tree/master/scripts
# 2. GLMM
# https://github.com/kdzimm/PseudoreplicationPaper/blob/master/Type_1_Error/Type%201%20-%20Tweedie%20GLMM.Rmd

#' Filter samples for differential expression test.
#'
#' @param sce SingleCellExperiment object of the aggregated data, output of `scater::aggregateAcrossCells`.
#' @param group_by_var A string. The grouping variable.
#' @param n_cells_dataset A numeric. The minimum number of cells per dataset for
#' the celltype to be considered valid
#' @param n_datasets_celltype A numeric. The mininum number of replicates (samples)
#' per condition in the `group_by_var` variable
#'
#' @returns A SingleCellExperiment object representing the data with valid samples
#' if the data is valid. Or NULL if the data is invalid.
#'
#' @importFrom magrittr %>%
#' @export
filter_samples_pseudobulk <- function(sce,
                                      group_by_var,
                                      n_cells_dataset,
                                      n_datasets_celltype){

  # filter invalid samples
  grouping_cols <- c(group_by_var)

  validity_df <- SummarizedExperiment::colData(sce) %>% as.data.frame() %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
    dplyr::summarize(n_datasets = dplyr::n(),
                     n_valid_datasets = sum(ncells >= n_cells_dataset))

  is_valid <- all(validity_df$n_valid_datasets >= n_datasets_celltype)
  if(!is_valid){
    return(NULL)
  }

  sce_valid <- sce[,sce[["ncells"]] >= n_cells_dataset]
  SummarizedExperiment::colData(sce_valid) <- droplevels(SummarizedExperiment::colData(sce_valid))

  return(sce_valid)
}

#' Check validity of arguments for DE_pseudobulk function
#'
#' @param sce A SingleCellExperiment object
#' @param group_colname A string specifying the grouping column name
#' @param normalization_method A string specifying the normalization method
#' @param test_method A string specifying the test method
#' @param extra_params A list of additional parameters passed via `...`
#'
#' @returns Logical.
#'
#' @keywords internal
check_args_pseudobulk <- function(sce,
                                  group_colname,
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
    stop(paste0("'group_colname' (\"", group_colname, "\") not found in colData of 'sce'. "))
  }

  # 3. check normalization_method is valid
  valid_norm_methods <- c("RLE", "TMM")
  normalization_method <- match.arg(normalization_method, valid_norm_methods)

  # 4. check test_method is valid
  valid_test_methods <- c("DESeq2_pseudobulk", "edgeR_pseudobulk", "limma_pseudobulk")
  test_method <- match.arg(test_method, valid_test_methods)

  # 5. check that extra_params arguments are valid
  valid_extra_args <- c(
    "min.count",
    "min.total.count",
    "DESeq2.test",
    "DESeq2.fitType",
    "DESeq2.sfType",
    "DESeq2.lfcShrink_type",
    "edgeR.test"
  )

  if(length(extra_params) > 0){
    extra_arg_names <- names(extra_params)
    invalid_args <- extra_arg_names[!extra_arg_names %in% valid_extra_args]
    if(length(invalid_args) > 0){
      stop(paste0("invalid argument(s) in '...': ", paste(invalid_args, collapse = ", "), ". "))
    }
  }

  return(TRUE)
}

#' Perform differential expression test using pseudobulking methods.
#'
#' @param sce An SingleCellExperiment object of the data. `sce` should be prefiltered to contain only valid samples
#' @param group_colname A string specifying the column (a factor) in the colData of `sce` to be used as grouping vairable
#' @param normalization_method A string specifying the normalization to use. One of `RLE`, and `TMM`
#' @param test_methood A string specifying the pseudobulk method to use.
#' @param ... Additional arguments. Valid options include: `min.count`, `min.total.count`,
#' `DESeq2.test`, `DESeq2.fitType`, `DESeq2.sfType`, `DESeq2.lfcShrink_type`, `edgeR.test`.
#'
#' @returns A list of outputs.
#'
#' @importFrom magrittr %>%
#'
#' @export
DE_pseudobulk <- function(sce,
                          group_colname,
                          normalization_method = c("RLE", "TMM"),
                          test_method = c("DESeq2_pseudobulk",
                                          "edgeR_pseudobulk",
                                          "limma_pseudobulk"),
                          ...){

  params <- list(...)

  # check variable validity
  args_check <- check_args_pseudobulk(
    sce = sce,
    group_colname = group_colname,
    normalization_method = normalization_method,
    test_method = test_method,
    extra_params = params
  )

  if(!args_check){
    stop("invalid inputs.")
  }

  # prepare argument list
  params[["normalization_method"]] <- normalization_method
  params[["test_method"]] <- test_method

  # change grouping variable name
  sce[[".grouping"]] <- sce[[group_colname]]

  # create required data object
  if(test_method == "DESeq2_pseudobulk"){

    pb <- DESeq2::DESeqDataSetFromMatrix(
      countData = SingleCellExperiment::counts(sce),
      colData = SummarizedExperiment::colData(sce),
      design = ~ .grouping
    )
  }else if(test_method %in% c("edgeR_pseudobulk", "limma_pseudobulk")){
    pb <- edgeR::DGEList(
      counts = SingleCellExperiment::counts(sce),
      sample = SummarizedExperiment::colData(sce)
    )
  }

  # filter genes
  if(!"min.count" %in% names(params))
    params[["min.count"]] <- 10
  if(!"min.total.count" %in% names(params))
    params[["min.total.count"]] <- 20

  keep <- edgeR::filterByExpr(pb,
                              group=sce[[group_colname]],
                              min.count = params[["min.count"]],
                              min.total.count = params[["min.total.count"]])
  pb <- pb[keep,]

  # compute normalization factor
  if(test_method == "DESeq2_pseudobulk"){
    if(normalization_method == "TMM"){
      tmp <- edgeR::DGEList(
        counts = DESeq2::counts(pb),
        sample = SummarizedExperiment::colData(pb)
      )
      tmp <- edgeR::calcNormFactors(tmp, method = "TMM")
      tmm_factors <- tmp$samples$norm.factors * tmp$samples$lib.size
      tmm_factors <- tmm_factors / exp(mean(log(tmm_factors)))
      DESeq2::sizeFactors(pb) <- tmm_factors
    }
    # normalization_method == "RLE": DESeq2 default

  }else if(test_method %in% c("edgeR_pseudobulk", "limma_pseudobulk")){
    pb <- edgeR::calcNormFactors(pb, method = normalization_method)
  }

  # fit model
  if(test_method == "DESeq2_pseudobulk"){
    # default parameters
    if(!"DESeq2.test" %in% names(params))
      params[["DESeq2.test"]] <- "LRT"

    if(!"DESeq2.fitType" %in% names(params)){
      params[["DESeq2.fitType"]] <- if (requireNamespace("glmGamPoi", quietly = TRUE)) {
        "glmGamPoi"
      } else {
        "parametric"
      }
    }

    if(!"DESeq2.sfType" %in% names(params))
      params[["DESeq2.sfType"]] <- "ratio"

    # if glmGamPoi is requested
    if(identical(params[["DESeq2.fitType"]], "glmGamPoi") &&
       !requireNamespace("glmGamPoi", quietly = TRUE)){
      stop(
        "DESeq2.fitType = 'glmGamPoi' requires the Bioconductor package 'glmGamPoi'. ",
        "Install it with BiocManager::install('glmGamPoi'), or use ",
        "DESeq2.fitType = 'parametric'.",
        call. = FALSE
      )
    }

    if(params[["DESeq2.test"]] == "LRT"){
      pb <- DESeq2::DESeq(pb,
                          test = params[["DESeq2.test"]],
                          fitType = params[["DESeq2.fitType"]],
                          sfType = params[["DESeq2.sfType"]],
                          reduced = ~1,
                          minReplicatesForReplace = Inf)
    }else{ # Wald test
      pb <- DESeq2::DESeq(pb,
                          test = params[["DESeq2.test"]],
                          fitType = params[["DESeq2.fitType"]],
                          sfType = params[["DESeq2.sfType"]])
    }
    design <- DESeq2::design(pb) # retrieve design for outputs

    # extract results
    contrast_name <- DESeq2::resultsNames(pb)[[2]]
    res <- DESeq2::results(pb,
                           name = contrast_name,
                           alpha = 0.05) # use 0.1 (default)?

    # apply lfcShrink if fitType is not glmGamPoi
    if(!params[["DESeq2.fitType"]] == "glmGamPoi"){
      if(!"DESeq2.lfcShrink_type" %in% names(params))
        params[["DESeq2.lfcShrink_type"]] <- "normal"

      res <- DESeq2::lfcShrink(pb, coef = 2, type = params[["DESeq2.lfcShrink_type"]], res = res)
    }

    # results table
    res_tbl <- res %>% data.frame() %>%
      tibble::add_column(gene_id = rownames(.), .before = 1) %>%
      tibble::add_column(gene_symbol = SummarizedExperiment::rowData(sce)[rownames(.),"gene_symbol"], .after = 1) %>%
      dplyr::mutate(fdr = p.adjust(pvalue, method = "fdr")) %>%
      dplyr::arrange(pvalue)

    # `fit` slot
    fit <- NULL

  }else if(test_method == "edgeR_pseudobulk"){
    # design
    design <- model.matrix(~ .grouping, pb$samples)
    # dispersion
    pb <- edgeR::estimateDisp(pb, design)

    # fit model
    if(!"edgeR.test" %in% names(params))
      params[["edgeR.test"]] <- "edgeRQLF" # options: edgeRQLF, edgeRLRT

    fit_func <- ifelse(params[["edgeR.test"]] == "edgeRQLF", edgeR::glmQLFit, edgeR::glmFit)
    test_func <- ifelse(params[["edgeR.test"]] == "edgeRQLF", edgeR::glmQLFTest, edgeR::glmLRT)

    fit <- fit_func(pb, design, robust=TRUE)
    res <- test_func(fit, coef=ncol(design))

    # results table
    res_tbl <- res$table %>%
      tibble::add_column(gene_id = rownames(.), .before = 1) %>%
      tibble::add_column(gene_symbol = SummarizedExperiment::rowData(sce)[rownames(.), "gene_symbol"], .after = 1) %>%
      dplyr::mutate(fdr = p.adjust(PValue, method = "fdr")) %>%
      dplyr::arrange(PValue)

    # update column names
    names(res_tbl)[which(names(res_tbl) == "logFC")] <- "log2FoldChange"
    names(res_tbl)[which(names(res_tbl) == "PValue")] <- "pvalue"

  }else{ # limma-voom
    # design
    design <- model.matrix(~ .grouping, pb$samples)

    v.beta <- limma::voomWithQualityWeights(pb, design)
    fit <- limma::lmFit(v.beta)
    fit <- limma::eBayes(fit, robust=TRUE)

    res <- limma::topTable(fit, sort.by="none", n=Inf, coef=ncol(design))

    res_tbl <- res %>%
      tibble::add_column(gene_id = rownames(.), .before = 1) %>%
      tibble::add_column(gene_symbol = SummarizedExperiment::rowData(sce)[rownames(.), "gene_symbol"], .after = 1) %>%
      dplyr::mutate(fdr = p.adjust(P.Value, method = "fdr")) %>%
      dplyr::arrange(P.Value)

    # update column names
    names(res_tbl)[which(names(res_tbl) == "logFC")] <- "log2FoldChange"
    names(res_tbl)[which(names(res_tbl) == "P.Value")] <- "pvalue"
    names(res_tbl)[which(names(res_tbl) == "adj.P.Val")] <- "adj.pvalue"
  }

  # return outputs
  outs <- list("pseudobulk_obj" = pb,
               "params" = params,
               "design" = design,
               "fit" = fit,
               "res" = res,
               "res_tbl" = res_tbl)
  return(outs)
}


#' Generate diagnostic plots for given outputs from the `DE_pseudobulk` function
#'
#' @param outs outputs from the `DE_pseudobulk` function
#' @param output_dir A string. Output directory to save plots
#'
#' @returns logical.
#'
#' @export
diagnostic_plots <- function(outs, output_dir){

  test_method <- outs$params$test_method
  pseudobulk_obj <- outs$pseudobulk_obj
  params <- outs$params

  if(test_method=="DESeq2_pseudobulk"){ # DESeq2 methods
    out_fn <- file.path(output_dir, "plot_DESeq2_DispEsts.pdf")
    pdf(out_fn, width = 5, height = 5)
    DESeq2::plotDispEsts(pseudobulk_obj)
    dev.off()

  }else if(test_method=="edgeR_pseudobulk"){ # edgeR methods
    # mean-difference
    out_fn1 <- file.path(output_dir, "plot_MD_samples.pdf")
    pdf(out_fn1, width = 10, height = 6)
    par(mfrow=c(2,ceiling(ncol(pseudobulk_obj)/2)))
    for (i in seq_len(ncol(pseudobulk_obj))) {
      limma::plotMD(pseudobulk_obj, column=i)
    }
    dev.off()

    # MDS plot
    out_fn2 <- file.path(output_dir, "plot_MDS.pdf")
    pdf(out_fn2, width = 7, height = 6)
    limma::plotMDS(edgeR::cpm(pseudobulk_obj, log=TRUE),
                   col=ifelse(pseudobulk_obj$samples$.grouping==levels(pseudobulk_obj$samples$.grouping)[[1]],
                              "#000004", "#3B719F"))
    dev.off()

    # plot dispersion
    out_fn3 <- file.path(output_dir, "plot_diagnostics.pdf")
    pdf(out_fn3, width = 10, height = 3.5)
    par(mfrow=c(1, 3))
    edgeR::plotBCV(pseudobulk_obj)
    edgeR::plotSmear(outs$res)
    if(params[["edgeR.test"]]=="edgeRQLF")
      edgeR::plotQLDisp(outs$fit)
    dev.off()

  }else if(test_method=="limma_pseudobulk"){
    # plot dispersion
    out_fn1 <- file.path(output_dir, "plot_diagnostics.pdf")
    pdf(out_fn1, width = 8, height = 3.5)
    ret_ <- limma::voomWithQualityWeights(pseudobulk_obj, outs$design, plot = T)
    dev.off()

  }else{
    return(FALSE)
  }

  return(TRUE)
}

