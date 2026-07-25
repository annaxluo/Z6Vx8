# utilities for plotting

#' Plot 2D scatter plot for a reduced dimension
#'
#' @param sce A SingleCellExperiment object.
#' @param color_by A string specifying color of the data point. One of the
#' column variables of colData of sce (discrete or continuous), gene name,
#' or a numeric vector.
#' @param color_by_label A string. Labels for color.
#' @param plot_title A string. Title of plot.
#' @param show_label Logical. Show cluster labels over data point or not.
#' @param reduced_dim_name A string. Name of the reduced dimension to use.
#' @param reduced_dim_1 A string. Column name of first reduced dimension.
#' @param reduced_dim_2 A string. Column name of second reduced dimension.
#' @param assay.type A string. The type of `assay` to use when plotting expression
#' of a gene.
#' @param split_by A string. Name of a column in `colData` used to split the scatter
#' plot into subplots.
#' @param alpha A numeric. Alpha value for data points.
#' @param size A numeric. Size value for data points.
#' @param rasterize Logical. To rasterize the scatter plot or not.
#' @param rastr_dpi An integer. dpi value for rasterizing.
#'
#' @returns A `ggplot` object.
#'
#' @importFrom magrittr %>%
#' @export
plot_reduced_dimension <- function(sce,
                                   color_by,
                                   color_by_label,
                                   plot_title,
                                   show_label = FALSE,
                                   reduced_dim_name = "UMAP",
                                   reduced_dim_1 = "V1",
                                   reduced_dim_2 = "V2",
                                   assay.type = "logcounts",
                                   split_by = NULL,
                                   alpha = .5,
                                   size = .1,
                                   rasterize = FALSE,
                                   rastr_dpi = 600){

  if(!is.null(split_by) && !(split_by %in% colnames(SummarizedExperiment::colData(sce)))){
    stop("invalid `split_by` variable")
  }

  if(is.character(color_by)){
    if(!(color_by %in% colnames(SummarizedExperiment::colData(sce)) | color_by %in% rownames(sce)))
      stop("invalid `color_by` variable")
    else if(color_by %in% rownames(sce) & !(assay.type %in% SummarizedExperiment::assayNames(sce)))
      stop(paste0("assay ", assay.type, " not in data"))
  }

  if(!is.character(color_by)){
    color_by_val <- color_by
    show_label = FALSE
  }else if(color_by %in% colnames(SummarizedExperiment::colData(sce))){
    color_by_val <- SummarizedExperiment::colData(sce)[[color_by]]
    cluster_lab <- NULL
    if(is.factor(color_by_val) | is.character(color_by_val)){
      cluster_lab <- SingleCellExperiment::reducedDim(sce, reduced_dim_name)[,c(1,2)] %>%
        cbind(SummarizedExperiment::colData(sce)) %>%
        as.data.frame() %>%
        dplyr::group_by_at(color_by) %>%
        dplyr::summarise_at(dplyr::vars(dplyr::one_of(reduced_dim_1, reduced_dim_2)), list(mean))

      colnames(cluster_lab) <- c("group", "x", "y")
    }else{
      show_label = FALSE # for continuous variable in colData
    }
  }else{ # plot gene expression
    color_by_val <- SummarizedExperiment::assay(sce, assay.type)[color_by,]
    show_label = FALSE
  }

  plot_data <- SingleCellExperiment::reducedDim(sce, reduced_dim_name)[,c(1,2)] %>%
    as.data.frame() %>%
    dplyr::mutate(color_by_val = color_by_val)

  pl_ <- ggplot2::ggplot(plot_data)

  if(rasterize){
    pl_ <- pl_ +
      ggrastr::geom_point_rast(
        ggplot2::aes_string(x=reduced_dim_1,
                            y=reduced_dim_2,
                            color="color_by_val"),
        alpha=alpha,
        size=size,
        raster.dpi=rastr_dpi
      )
  }else{
    pl_ <- pl_ +
      ggplot2::geom_point(
        ggplot2::aes_string(x=reduced_dim_1,
                            y=reduced_dim_2,
                            color="color_by_val"),
        alpha=alpha,
        size=size
      )
  }

  pl_ <- pl_ +
    ggplot2::labs(color=color_by_label,
                  x=paste0(reduced_dim_name, "1"),
                  y=paste0(reduced_dim_name, "2"),
                  title=plot_title) +
    ggplot2::theme_classic()

  if(show_label & is.null(split_by)){
    pl_ <- pl_ + ggplot2::geom_text(data=cluster_lab,
                                    ggplot2::aes(label=group, x, y),
                                    size=2)
  }

  if(!is.null(split_by)){
    if(split_by %in% colnames(SummarizedExperiment::colData(sce))){
      pl_ <- pl_ +
        ggplot2::facet_wrap(~ SummarizedExperiment::colData(sce)[,split_by], scales="free") +
        ggplot2::scale_x_continuous(limits = c(ggplot2::layer_scales(pl_)$x$range$range[[1]],
                                               ggplot2::layer_scales(pl_)$x$range$range[[2]])) +
        ggplot2::scale_y_continuous(limits = c(ggplot2::layer_scales(pl_)$y$range$range[[1]],
                                               ggplot2::layer_scales(pl_)$y$range$range[[2]])) +
        ggplot2::theme(axis.line = ggplot2::element_line(),
                       strip.background = ggplot2::element_blank())
    }
  }
  return(pl_)
}


#' Internal function to compute marker data for plotting stacked violin or heatmap.
#'
#' @param sce A SingleCellExperiment object.
#' @param marker_dict A list of markers in a format of list("A"=c('gene1','gene2',...)).
#' as the grouping variable.
#' @param row_ensembl Logical. Column names of `sce` is Ensembl ID or not.
#'
#' @returns A `dgCMatrix` matrix of expression, with each row being a sample (cell),
#' and each column a marker gene.
plot_data_ <- function(sce,
                       marker_dict,
                       row_ensembl=TRUE){

  marker_list <- stack(marker_dict) %>% as.data.frame()
  colnames(marker_list) <- c("gene_symbol", "cell_type")
  marker_list <- marker_list %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  if(row_ensembl){
    marker_list$ensembl_id <- rownames(sce)[match(marker_list$gene_symbol,
                                                  SummarizedExperiment::rowData(sce)$gene_symbol)]
    idx_valid <- sapply(marker_list$ensembl_id, function(s){ s %in% rownames(sce) })
    marker_list <- marker_list[idx_valid,]
    plot_data <- Matrix::t(SummarizedExperiment::assay(sce, "logcounts")[marker_list$ensembl_id,,drop=FALSE])
  }else{
    idx_valid <- sapply(marker_list$gene_symbol, function(s){ s %in% rownames(sce) })
    marker_list <- marker_list[idx_valid,]
    plot_data <- Matrix::t(SummarizedExperiment::assay(sce, "logcounts")[marker_list$gene_symbol,,drop=FALSE])
  }

  return(list("marker_list"=marker_list, "plot_data"=plot_data))
}

#' Plot stacked violin plots for genes
#'
#' @param sce A SingleCellExperiment object.
#' @param marker_dict A list of markers in a format of list("A"=c('gene1','gene2',...)).
#' @param cluster_lab A string. Specifying a column in colData(sce) to be used
#' as the grouping variable.
#' @param row_ensembl Logical. Column names of `sce` is Ensembl ID or not.
#' @param plot_title A string. Plot title.
#'
#' @returns A `ggplot` object.
#'
#' @importFrom magrittr %>%
#'
#' @export
plot_stackedViolin <- function(sce,
                               marker_dict,
                               cluster_lab,
                               row_ensembl=TRUE,
                               plot_title="stacked violin plot",
                               x_lab="cluster",
                               y_lab="expression level (log)"){

  # check "logcounts" in assay
  if(!"logcounts" %in% SummarizedExperiment::assayNames(sce)){
    stop("logcounts not in sce.")
  }

  if(!cluster_lab %in% colnames(SummarizedExperiment::colData(sce))){
    stop("invalid `cluster_lab` variable")
  }

  outs <- plot_data_(sce, marker_dict, row_ensembl)

  plot_df <- outs$plot_data %>%
    as.matrix() %>%
    as.data.frame()

  marker_list <- outs$marker_list

  colnames(plot_df) <- marker_list$gene_symbol
  plot_df$cell_id <- rownames(plot_df)
  plot_df$cluster <- SummarizedExperiment::colData(sce)[,cluster_lab]
  if(!is.factor(plot_df$cluster)){
    plot_df$cluster <- factor(plot_df$cluster)
  }

  plot_df <- plot_df %>%
    tidyr::pivot_longer(cols = -c(cell_id, cluster), names_to = "feature",
                        values_to = "logcounts")
  plot_df$feature <- factor(plot_df$feature, levels = marker_list$gene_symbol)

  # x label order
  xlab_order <- data.frame(
    label = levels(plot_df$cluster),
    order = seq(length(levels(plot_df$cluster)))
  )
  plot_df$cluster_order <- xlab_order$order[match(plot_df$cluster, xlab_order$label)]
  plot_df$cluster_order <- factor(plot_df$cluster_order)

  # plot
  pl_ <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = cluster_order,
                                      y = logcounts,
                                      fill = feature)) +
    ggplot2::geom_violin(scale = "width", adjust = 1, trim = TRUE) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                position="right",
                                labels = function(x)
                                  c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    ggplot2::scale_x_discrete(labels = levels(plot_df$cluster)) +
    ggplot2::facet_grid(rows = ggplot2::vars(feature), scales = "free", switch = "y") +
    cowplot::theme_cowplot(font_size = 12) +
    ggplot2::theme(legend.position = "none",
                   panel.spacing = ggplot2::unit(0, "lines"),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.background = ggplot2::element_rect(fill = NA, color = "black"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"),
                   strip.text.y.left = ggplot2::element_text(angle = 0),
                   axis.text.x = ggplot2::element_text(size = 8, angle = 45)) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab)

  return(pl_)
}


#' Plot heatmap of markers
#' @param sce A SingleCellExperiment object.
#' @param marker_dict A list of markers in a format of list("A"=c('gene1','gene2',...)).
#' @param cluster_lab A string. Specifying a column in colData(sce) to be used
#' as the grouping variable.
#' @param row_ensembl Logical. Column names of `sce` is Ensembl ID or not.
#' @param scale_expr Logical. Scale expression for each gene or not.
#' @param cluster_rows Logical. Cluster heatmap rows (genes) or not.
#' @param cluster_cols Logical. Cluster heatmap columns (groups) or not.
#' @param plot_title A string. Plot title.
#'
#' @returns A `pheatmap` object.
#'
#' @importFrom magrittr %>%
#'
#' @export
plot_marker_heatmap <- function(sce,
                                marker_dict,
                                cluster_lab,
                                row_ensembl=TRUE,
                                scale_expr=TRUE,
                                cluster_rows=TRUE,
                                cluster_cols=TRUE,
                                plot_title="heatmap"){

  # check "logcounts" in assay
  if(!"logcounts" %in% SummarizedExperiment::assayNames(sce)){
    stop("logcounts not in sce.")
  }

  if(!cluster_lab %in% colnames(SummarizedExperiment::colData(sce))){
    stop("invalid `cluster_lab` variable")
  }

  outs <- plot_data_(sce, marker_dict, row_ensembl)

  plot_mat <- outs$plot_data %>%
    as.matrix()

  # mean expr
  expr_mat <- do.call("rbind", lapply(levels(sce[[cluster_lab]]), function(c_){
    rr_ <- plot_mat[sce[[cluster_lab]]==c_,]
    if(nrow(rr_) > 1){
      colMeans(rr_)
    }else{
      rr_
    }
  }))

  rownames(expr_mat) <- levels(sce[[cluster_lab]])
  colnames(expr_mat) <- SummarizedExperiment::rowData(sce)[colnames(expr_mat),"gene_symbol"]

  if(scale_expr){
    expr_mat <- scale(expr_mat)
  }

  # plot heatmap
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(101)
  limits_ <- ceiling(max(abs(expr_mat), na.rm=T))
  breaks <- seq(-limits_, limits_, length.out = 101)

  pl_ <- pheatmap::pheatmap(expr_mat,
                            color = colors,
                            breaks = breaks,
                            cluster_rows = cluster_rows,
                            cluster_cols = cluster_cols,
                            main = plot_title)

  return(pl_)
}

