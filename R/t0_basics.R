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
  hvf.info <- .FindVariableFeatures_vst(
    object = dat,
    loess.span = 0.3,
    clip.max = "auto"
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

  # fit normal with 25% to 75%
  IQR <- quantile(vec,
                  probs = c(0.25, 0.75),
                  na.rm = TRUE)
  m <- mean(IQR)
  delta <- (IQR[2] - IQR[1]) / (qnorm(0.75) - qnorm(0.25))

  (vec - m) / delta
}

#' Make mock aggregated data
#'
#' @param ncol numeric, number of columns
#' @param nrow numeric, number of rows
#' @returns A SingleCellExperiment object
#'
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

#' Replementation of the Seurat's `FindVariableFeatures` to apply Variance
#' Stabilizing Transformation on a matrix or dgCMatrix object.
#' @param counts A Matrix or inheriting from Matrix.
#' @param loess.span A numeric. Loess span parameter for VST
#' @param clip.max A numeric or "auto". The maximum value allowed after standardization.
#'
#' @returns A data.frame, specifying the mean, variance, expected variance
#' based on the mean-variance trend, and the standardized expected variance for
#' each gene.
#'
#' @keywords internal
.FindVariableFeatures_vst <- function(object,
                                      loess.span = 0.3,
                                      clip.max = "auto"){

  if(!requireNamespace("Matrix", quietly = TRUE)){
    stop("Package 'Matrix' is required.")
  }

  n_cells <- ncol(object)
  n_features <- nrow(object)

  if(clip.max == "auto"){
    clip.max <- sqrt(n_cells)
  }

  # sample variance per row
  row_means <- Matrix::rowMeans(object)
  row_sumsq <- Matrix::rowSums(object^2)
  row_vars <- (row_sumsq - n_cells * row_means^2) / (n_cells - 1)
  row_vars[row_vars < 0] <- 0

  hvf.info <- data.frame(
    mean = row_means,
    variance = row_vars,
    variance.expected = 0,
    variance.standardized = 0,
    row.names = rownames(object)
  )

  not.const <- hvf.info$variance > 0 & hvf.info$mean > 0

  fit <- stats::loess(
    formula = log10(variance) ~ log10(mean),
    data = hvf.info[not.const, , drop = FALSE],
    span = loess.span
  )

  hvf.info$variance.expected[not.const] <- 10^fit$fitted

  # variance of standardized values after clipping standardized values at clip.max.
  expected_sd <- sqrt(hvf.info$variance.expected)

  variance.standardized <- numeric(n_features)

  for(i in seq_len(n_features)){
    if(expected_sd[i] == 0 || is.na(expected_sd[i])){
      variance.standardized[i] <- 0
      next
    }

    row_start <- object@p[1:n_cells] + 1L
    row_end <- object@p[2:(n_cells + 1L)]

    ## Easier and clear standalone implementation:
    x <- as.numeric(object[i, ])

    z <- (x - row_means[i]) / expected_sd[i]
    z[z > clip.max] <- clip.max

    variance.standardized[i] <- stats::var(z)
  }

  hvf.info$variance.standardized <- variance.standardized

  colnames(hvf.info) <- paste0("vst.", colnames(hvf.info))

  return(hvf.info)
}

