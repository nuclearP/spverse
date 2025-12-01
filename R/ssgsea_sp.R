#' @title ssgsea_sp
#' @description
#' Estimates ssgsea enrichment scores for the given sp object
#' @param object an sp object \linkS4class{sp}.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Perform ssgsea enrichment analysis using the slot clean_data, and store the returned result in ssgsea_list[["ssgsea_list_clean"]].
#'   \item RID : Perform ssgsea enrichment analysis using the slot RID_bootstrap, and store the returned result in ssgsea_list[["ssgsea_list_RID"]].}
#' @param pathway a file path to a gmt format file, or a list storing vectors of different gene or protein IDs.
#' @param fromType The ID type of the row names in the selected expression matrix is default to NULL,
#' and it is used when converting the row names of the selected expression matrix to other ID type.
#' @param toType Corresponds to the fromType parameter, specifying the converted ID type.
#' @param species The name of the species-specific annotation package, such as "org.Hs.eg.db".
#' @details The ID type of the selected expression matrix should be consistent with those contained in the input pathway.
#' @export
#' @rdname ssgsea_sp
ssgsea_sp <- function(object,name="clean",pathway,fromType=NULL,toType=NULL,species=NULL){

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
    expr_df <- object@clean_data
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
    expr_df <- object@RID_bootstrap
  }

  if(is.list(pathway)){
    pathway_list <- list()
    pathway_name <- names(pathway)
    for (i in pathway_name) {
      gene_ids <- pathway[[i]]
      my_geneset <- GeneSet(
        geneIds = gene_ids,
        setName = i
      )
      pathway_list[[i]] <- my_geneset
    }
    gsc <- GeneSetCollection(pathway_list)
  }else{
    gsc <- getGmt(pathway)
  }

  if((!is.null(fromType)) & (!is.null(toType))){

    ids <- row.names(expr_df)

    if(!is.null(species)){
      id_conversion <- bitr(
        geneID = ids,
        fromType = fromType,
        toType = toType,
        OrgDb = species)
    }else{
      stop("You have not specified the name of the species-specific annotation package.")
    }

    expr_df[[fromType]] <- ids
    expr_df <- full_join(expr_df,id_conversion)
    expr_df[[fromType]] <- NULL
    expr_df <- expr_df %>% na.omit() %>% group_by(.data[[toType]]) %>%
      summarise(across(where(is.numeric),mean))
    row_name <- expr_df[[toType]]
    expr_df[[toType]] <- NULL
    row.names(expr_df) <- row_name
  }

  expr_matrix <- as.matrix(expr_df)

  gsvapar <- GSVA::ssgseaParam(exprData = expr_matrix, geneSets = gsc, normalize = T)
  scores <- GSVA::gsva(gsvapar)
  object@ssgsea_list[[paste0("ssgsea_",name)]] <- scores
  return(object)
}
