# utilities for running pathway analysis
# 1. functions to retrieve pathways (GO, KEGG, Reactome, SynGO)
# 2. plot function: dot plot

#' Process a pathway (Gene Ontology, Reactome, KEGG, SynGO)
#'
#' @param pathway_str A String, specifies type of pathways (one of GO_BP, GO_CC,
#'                    GO_MF, Reactome, KEGG, SynGO)
#' @param pathway_output_path A String, output pathway
#' @param pathway_db_str A String, version of processed pathways
#'
#' @returns TRUE if no error.
#'
#' @export
process_pathway <- function(pathway_output_path,
                            pathway_db_ver,
                            pathway_str=c("GO_BP","GO_CC","GO_MF","Reactome","KEGG","SynGO"),
                            ensembl_version="110",
                            syngo_path=NULL){

  # SynGO
  if(pathway_str == "SynGO" & !is.null(syngo_path)){
    pathway_list <- extract_SynGO_pathways(syngo_path, ensembl_version=ensembl_version)
    path_fn <- file.path(pathway_output_path,
                         paste0("pathways_", pathway_str, "_", pathway_db_ver, ".rds"))
    saveRDS(pathway_list, path_fn)
    return(TRUE)
  }

  # GO terms
  if(pathway_str %in% c("GO_BP", "GO_CC", "GO_MF")){
    env_ <- clusterProfiler:::get_GO_data('org.Mm.eg.db',
                                          ont=stringr::str_split_i(pathway_str, "_", 2),
                                          keytype="ENSEMBL")
    type_ <- "GO"
  }

  # Reactome
  if(pathway_str == "Reactome"){
    env_ <- ReactomePA:::get_Reactome_DATA("mouse")
    type_ <- "Reactome"
  }

  # KEGG
  if(pathway_str == "KEGG"){
    species <- clusterProfiler:::organismMapper("mouse") # "mmu"
    env_ <- clusterProfiler:::prepare_KEGG(species)
    type_ <- "KEGG"
  }

  # extract pathways from environment
  pathway_list <- extract_list_from_env(env_, type=type_, ensembl_version=ensembl_version)
  path_fn <- file.path(pathway_output_path,
                       paste0("pathways_", sub("_", "", pathway_str), "_", pathway_db_ver, ".rds"))
  saveRDS(pathway_list, path_fn)

  return(TRUE)
}


#' Extracts pathways from an environment containing the pathway list.
#'
#' @param used_env An environment containing the pathway list.
#' @param type A string, specifying the type of pathway. One of "GO", "Reactome",
#' and "KEGG"
#' @param ensembl_version A string, specifying the Ensembl version to use.
#'
#' @returns A list. List of pathways (names and Ensembl IDs).
extract_list_from_env <- function(used_env,
                                  type=c("GO", "Reactome", "KEGG"),
                                  ensembl_version="110"){
  #require(org.Mm.eg.db)
  # list of all unique genes
  pathway_list_0 <- get("PATHID2EXTID", envir=used_env)
  all_genes_input <- unique(unlist(pathway_list_0))

  version_df <- biomaRt::listEnsemblArchives()
  host_ <- version_df$url[match("110", version_df$version)]

  # convert to MGI gene symbols
  ensembl <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "mmusculus_gene_ensembl",
    host = host_
  )

  if(type == "GO"){
    filter_type = "ensembl_gene_id"
  }else{
    filter_type = "entrezgene_id"
  }

  mapping_ <- biomaRt::getBM(
    attributes = c(filter_type, "mgi_symbol"),
    filters = filter_type,
    values = all_genes_input,
    mart = ensembl
  )

  # map genes
  pathway_list_1 <- lapply(pathway_list_0, function(p_){
    unique(mapping_[["mgi_symbol"]][mapping_[[filter_type]] %in% p_])
  })
  names(pathway_list_1) <- names(pathway_list_0)

  # remove empty lists
  pathway_list_1 <- pathway_list_1[!sapply(pathway_list_1, length)==0]

  list_1 <- get("PATHID2NAME", envir=used_env)
  list_names <- list_1[names(pathway_list_1)]
  if(type == "KEGG"){
    # remove  "- Mus musculus (house mouse)"
    list_names <- stringr::str_replace(list_names, " - Mus musculus \\(house mouse\\)", "")
    names(list_names) <- names(pathway_list_1)
  }

  # add original pathway id
  attr(pathway_list_1, "pathway_name") <- list_names
  return(pathway_list_1)

}

#' Extract SynGO pathways
#'
#' @param syngo_path path to SynGO downloads
#' @returns A list. List of SynGO pathways in gene symbol
extract_SynGO_pathways <- function(syngo_path, ensembl_version="110"){

  syngo_annotation <- openxlsx::read.xlsx(syngo_path)

  # extract pathways
  pathway_list_0 <- lapply(syngo_annotation$ensembl_id, function(gl_){
    trimws(stringr::str_split(gl_, ",")[[1]])
  })

  # map human gene symbols to mouse gene symbols
  human_gene_symbols <- unique(do.call("c", pathway_list_0))
  mapping_ <- map_human_gene_to_mouse_gene_symbol(human_gene_symbols,
                                                  input_type="ensembl_id",
                                                  ensembl_version=ensembl_version)

  pathway_list_1 <- lapply(pathway_list_0, function(p_){
    unique(mapping_[["mouse_gene_symbol"]][mapping_[["input"]] %in% p_])
  })
  names(pathway_list_1) <- paste0("SynGO_", syngo_annotation$id)

  # remove empty lists
  pathway_list_1 <- pathway_list_1[!sapply(pathway_list_1, length)==0]

  # add pathway attributes
  attr(pathway_list_1, "domain") <- syngo_annotation$domain
  p_names <- syngo_annotation$name
  names(p_names) <- names(pathway_list_1)
  attr(pathway_list_1, "pathway_name") <- p_names
  attr(pathway_list_1, "pathway_shortname") <- syngo_annotation$shortname

  return(pathway_list_1)
}

#' Map human genes (Ensembl ID or gene symbol) to mouse gene symbol
#'
#' @param gene_list A list of human genes
#' @param input_type A string. Type of the input human genes. One of `ensembl_id`,
#' `entrez_id`, and `gene_symbol`
#' @param ensembl_version A string. Ensembl version used.
#'
#' @returns A data.frame. The first columns is the input human gene, and the
#' second the associated mouse Ensembl ID
map_human_gene_to_mouse_gene_symbol <- function(gene_list,
                                                input_type=c("ensembl_id","entrez_id","gene_symbol"),
                                                ensembl_version="110"){

  # convert input to ensembl IDs
  if(input_type == "ensembl_id"){
    keytype_ <- "ENSEMBL"
  }else if(input_type == "entrez_id"){
    keytype_ <- "ENTREZID"
  }else{
    keytype_ <- "SYMBOL"
  }

  human_ensembl_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                             keys=gene_list,
                                             keytype=keytype_,
                                             column="ENSEMBL",
                                             multiVals="first")

  # remove NA's
  human_ensembl_ids <- human_ensembl_ids[!is.na(human_ensembl_ids)]

  # map to mouse gene symbols --------------------
  version_df <- biomaRt::listEnsemblArchives()
  host_ <- version_df$url[match("110", version_df$version)]

  # convert to MGI gene symbols
  human <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl",
    host = host_
  )

  mouse_map <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "mmusculus_homolog_associated_gene_name"
    ),
    filters = "ensembl_gene_id",
    values = unique(unname(human_ensembl_ids)),
    mart = human
  )

  mouse_map <- mouse_map[mouse_map$mmusculus_homolog_associated_gene_name != "",]

  input_map <- data.frame(
    input = names(human_ensembl_ids),
    ensembl_gene_id = unname(human_ensembl_ids)
  )

  genesV2 <- merge(input_map, mouse_map, by = "ensembl_gene_id", all.x = FALSE)

  colnames(genesV2) <- c("human_ensembl_gene_id", "input", "mouse_gene_symbol")

  genesV2 <- genesV2[, c("input","human_ensembl_gene_id","mouse_gene_symbol")]

  return(genesV2)
}


#' Plot enrichment analysis results
#'
#' @param df A data.frame. GSEA outputs
#' @param title_str  A string. Plot title.
#' @param p_value_str A string. Column name in the data.frame to be used as the probability value.
#' @param p_value_threshold A numeric. Probability value threshold.
#' @param n_limit A numeric. Maxiumum number of pathways to plot. If `NULL`, plot
#' all pathways.
#' @param n_min A numeric. Minimum number of pathways to plot.
#' @param y_order_by A string. Variable name for x axis and ordering.
#' @param color_by A string. Variable name for color.
#' @param size_by A string. Variable name for dot size.
#' @param pathway_label A string. Variable name for pathway labels.
#' @param pathway_type_label A string. Variable name for pathway label annotations.
#' @param x_axis_limits A nuermic vector of size 2, specifying the x axis limits.
#' @param colorbar_limits A nuermic vector of size 2, specifying the color legend limits.
#'
#' @returns A ggplot2 object or NULL
#'
#' @importFrom magrittr %>%
#'
#' @export
gsea_dotplot <- function(df,
                         title_str,
                         p_value_str="fdr",
                         p_value_threshold=0.05,
                         n_limit=NULL,
                         n_min=5,
                         y_order_by="NES",
                         color_by="fdr",
                         size_by="gene_ratio",
                         pathway_label="pathway_name",
                         pathway_type_label="pathway_type",
                         x_axis_limits=NULL,
                         colorbar_limits=NULL){

  # using padj2 or padj
  if(!p_value_str %in% colnames(df) | !y_order_by %in% colnames(df) |
     !size_by %in% colnames(df) | !pathway_label %in% colnames(df)){
    stop("Invalud slots.")
  }

  if(nrow(df) < n_min){
    print("Too few pathways.")
    return(NULL)
  }

  df <- df %>% as.data.frame()

  # gene ratio
  if(is.null(n_limit)){
    n_limit <- min(100, nrow(df))
    if(sum(df[,p_value_str] < p_value_threshold) < n_min){
      df2 <- df %>%
        dplyr::arrange(pval) %>%
        dplyr::slice(1:min(nrow(df), n_min))
    }else{
      df2 <- df %>% dplyr::filter(df[,p_value_str] < p_value_threshold)
      if(nrow(df2) > n_limit){
        df2 <- df2 %>%
          dplyr::arrange(pval) %>%
          dplyr::slice(1:n_limit)
      }
    }
  }else{
    df2 <- df %>% dplyr::filter(df[,p_value_str] < p_value_threshold)
    n_limit <- min(n_limit, nrow(df2))
    df2 <- df2 %>%
      dplyr::arrange(pval) %>%
      dplyr::slice(1:n_limit)
  }

  # handle case where no pathways pass threshold
  if (nrow(df2) == 0) {
    message("No pathways pass the threshold (", p_value_str, " < ",
            p_value_threshold, ").")
    return(NULL)
  }

  # deal with labels
  # decreasing order: p values. Increasing order: other
  is_decreasing <- ifelse(y_order_by %in% c("pval", "padj", "fdr"), FALSE, TRUE)
  idx <- order(df2[[y_order_by]], decreasing=is_decreasing)

  # wrap text
  df2$pathway_label <- paste0(substr(df2[,pathway_label], 1, 60), " (",
                              df2[,pathway_type_label], ")")

  df2$description <- factor(df2$pathway_label,
                            levels=rev(unique(df2$pathway_label[idx])))

  y_order_by_aes <- ifelse(y_order_by %in% c("pval", "padj", "fdr"),
                           paste0("-log10(", y_order_by, ")"),
                           y_order_by)
  x_lab <- stringr::str_replace(y_order_by_aes, "_", " ")

  color_by_aes <- ifelse(color_by %in% c("pval", "padj", "fdr"),
                         paste0("-log10(", color_by, ")"),
                         color_by)

  p1 <- ggplot2::ggplot(df2,
                        ggplot2::aes_string(x=y_order_by_aes,
                                            y="description",
                                            size=size_by,
                                            color=color_by_aes)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(limits=x_axis_limits) +
    ggplot2::ylab(NULL) +
    ggplot2::xlab(x_lab) +
    ggplot2::ggtitle(title_str) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9, face="bold"),
                   axis.title=ggplot2::element_text(size=9),
                   legend.title=ggplot2::element_text(size=9)) +
    ggplot2::guides(size=ggplot2::guide_legend(order=1))

  if(is.factor(df2[[color_by]])){
    p1 <- p1 + ggplot2::scale_color_discrete(name=color_by_aes)
  }else{
    p1 <- p1 +
      ggplot2::scale_color_continuous(low="blue", high="red",
                                      name=color_by_aes,
                                      limits=colorbar_limits) +
      ggplot2::guides(color=ggplot2::guide_colorbar(order=2, reverse=FALSE))
  }

  return(p1)
}




