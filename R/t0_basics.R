# utilities

#' Compute size factor for a counts matrix
#' @param counts A Matrix or inheriting from Matrix.
#' @param round_exprs A logical.
#' @returns A numeric vector
#' @export
estimate_sf_sparse <- function(counts,
                               round_exprs=TRUE){
  if(round_exprs)
    counts <- round(counts)

  cell_total <- Matrix::colSums(counts)
  sfs <- cell_total / exp(mean(log(cell_total)))

  sfs[is.na(sfs)] <- 1
  sfs
}

#' Compute highly-variable genes for input matrix with given criteria
#' @param expr_mat dgCMatrix, expression matrix, normalized, log2-transformed
#' @param used_cells logical, used cells
#' @param fdr_top_cutoff numeric (0-1), FDR cutoff for top variable gene selection based
#' @param n_top_genes numeric, number of genes with top variability. Ignored if `fdr_top_cutoff` is not null
#' @returns A logical vector
#' @export
compute_hvg <- function(expr_mat,
                        used_cells,
                        fdr_top_cutoff=.5,
                        n_top_genes=3000){

  dat <- expr_mat[,used_cells]

  # Fit loess regression
  hvf.info <- Seurat::FindVariableFeatures(
    object = dat,
    selection.method = "vst",
    loess.span = 0.3,
    verbose = FALSE
  )

  if(is.null(fdr_top_cutoff)){
    hvg_ <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE),,drop = FALSE]
    hvf.info$is.hvg <- FALSE
    hvf.info[rownames(hvg_)[1:n_top_genes],]$is.hvg <- TRUE
  }else{
    # compute z scores
    z_scores <- compute_z_score(vec = hvf.info$vst.variance.standardized)
    pval <- 1 - pnorm(z_scores)
    padj <- p.adjust(pval, method = "fdr")
    hvf.info$is.hvg <- FALSE
    hvf.info[padj < fdr_top_cutoff,]$is.hvg <- TRUE
  }

  # return hvgs (logical vector)
  hvf.info <- hvf.info[rownames(expr_mat),]

  return(hvf.info$is.hvg)
}


#' Compute z score for a vector
#' Adapted from: https://github.com/AllenInstitute/scrattch.hicat/blob/master/R/vg.R#L122
#' This version of z scoring will differ from base R's scale()
#'
#' @param vec a numeric vector
#' @returns a numeric vector of z-scores
compute_z_score <- function(vec){

  ###fit normal with 25% to 75%
  IQR <- quantile(vec,
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))

  (vec  - m) / delta
}


#' Make mock aggregated data
#'
#' @param ncol numeric, number of columns
#' @param nrow numeric, number of rows
#' @returns A SingleCellExperiment object
#' @export
make_summed_sce <- function(ncol_=500, nrow_=400){

  data("sce", package = "Z6Vx8")

  summed <- scater::aggregateAcrossCells(
    sce,
    statistics="sum",
    ids=SummarizedExperiment::colData(sce)[,c("celltype", "dataset_id")]
  )

  col_data <- SummarizedExperiment::colData(summed)
  dup_cols <- duplicated(colnames(col_data))
  SummarizedExperiment::colData(summed) <- col_data[, !dup_cols]

  return(summed)
}

