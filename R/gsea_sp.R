#' @title gsea_sp
#' @description
#' Estimates gsea enrichment scores for the given sp object
#' @param object an sp object \linkS4class{sp}.
#' @param protein_rank Named vector of gene-level stats.
#' @param dim_reduction_name The dimensionality reduction objects stored in the slot dim_reductions of an sp-type object
#' can be either PCA-type data or PLS-DA-type data. such as "pca_clean".
#' @param dim_value A positive integer value, or a vector composed of positive integers.
#' When using the dimensionality reduction data contained within the sp object as the basis for estimating gsea enrichment scores,
#' the ranking must be based on the contribution of genes across different dimensions of the dimensionality reduction data. Therefore,
#' it is necessary to select dimension information: if PCA data is chosen, the value in the parameter must not exceed the number of samples;
#' if PLS-DA data is chosen, the value in the parameter must be either not exceed 2.
#' @param pathway a file path to a gmt format file, or a list storing vectors of different gene or protein IDs.
#' @param fromType The ID type of the row names in the selected expression matrix is default to NULL,
#' and it is used when converting the row names of the selected expression matrix to other ID type.
#' @param toType Corresponds to the fromType parameter, specifying the converted ID type.
#' @param species The name of the species-specific annotation package, such as "org.Hs.eg.db".
#' @details If the protein_rank parameter is not NULL, the dimensionality reduction data contained
#' in the sp object will not be used, and the parameters dim_reduction_name and dim_value
#' will be ignored.
#' @export
#' @rdname gsea_sp
gsea_sp <- function(object,protein_rank=NULL,dim_reduction_name,dim_value=c(1),pathway,fromType=NULL,toType=NULL,species=NULL){

  if(!is.list(pathway)){
    pathway <- getGmt(pathway)
    pathway <- geneIds(pathway)
  }

  if(is.null(protein_rank)){##dim_reduction
    if(grepl("plsda",dim_reduction_name)){##plsda
      if(any(dim_value != 1 & dim_value != 2)){
        stop("PLS-DA only allows the calculation of GSEA scores using the gene rankings from the first and second dimensions.")
      }

      df_plsda <- object@dim_reductions[[dim_reduction_name]][["loadings"]][["X"]]
      ids <- row.names(df_plsda)

      plsda_list <- list()
      dim_value <- sort(dim_value)
      for (i in dim_value) {

        if((!is.null(fromType)) & (!is.null(toType))){##ID conversion

          if(!is.null(species)){
            protein_rank <- df_plsda[,i]
            protein_rank <- as.data.frame(protein_rank)
            protein_rank[[fromType]] <- ids
            colnames(protein_rank) <- c("value",fromType)

            id_conversion <- bitr(
              geneID = ids,
              fromType = fromType,
              toType = toType,
              OrgDb = species)

            mix_df <- full_join(protein_rank,id_conversion)
            mix_df[[fromType]] <- NULL
            mix_df <- mix_df %>% na.omit() %>% group_by(.data[[toType]]) %>%
              summarise(mean_value = mean(value))

            ranked_genes <- mix_df[["mean_value"]]
            names(ranked_genes) <- mix_df[[toType]]
            ranked_genes <- sort(ranked_genes,decreasing = F)
          }else{
            stop("You have not specified the name of the species-specific annotation package.")
          }

        }else{
          protein_rank <- df_plsda[,i]
          ranked_genes <- sort(protein_rank,decreasing = F)
        }

        fgsea_result <- fgsea(
          pathways = pathway,
          stats = ranked_genes,
          minSize = 1,
          maxSize = 500,
          nperm = 10000)

        fgsea_result[["Dims"]] <- paste0("comp",i)

        plsda_list[[i]] <- fgsea_result
      }
      plsda <- do.call(rbind, plsda_list)
      object@GSEA_list[[dim_reduction_name]] <- plsda
    }else{##pca

      pca_fit <- object@dim_reductions[[dim_reduction_name]]
      var <- get_pca_var(pca_fit)
      cors <- var[["cor"]]
      k <- dim(cors)[2]
      if(any(dim_value>k)){
        stop("The maximum dimension must not exceed"," ",k,".")
      }
      ids <- row.names(cors)

      pca_list <- list()
      dim_value <- sort(dim_value)

      for (i in dim_value) {

        if((!is.null(fromType)) & (!is.null(toType))){##ID conversion

          if(!is.null(species)){
            protein_rank <- cors[,i]
            protein_rank <- as.data.frame(protein_rank)
            protein_rank[[fromType]] <- ids
            colnames(protein_rank) <- c("value",fromType)

            id_conversion <- bitr(
              geneID = ids,
              fromType = fromType,
              toType = toType,
              OrgDb = species)

            mix_df <- full_join(protein_rank,id_conversion)
            mix_df[[fromType]] <- NULL
            mix_df <- mix_df %>% na.omit() %>% group_by(.data[[toType]]) %>%
              summarise(mean_value = mean(value))

            ranked_genes <- mix_df[["mean_value"]]
            names(ranked_genes) <- mix_df[[toType]]
            ranked_genes <- sort(ranked_genes,decreasing = F)
          }else{
            stop("You have not specified the name of the species-specific annotation package.")
          }

        }else{
          protein_rank <- cors[,i]
          ranked_genes <- sort(protein_rank,decreasing = F)
        }

        fgsea_result <- fgsea(
          pathways = pathway,
          stats = ranked_genes,
          minSize = 1,
          maxSize = 500,
          nperm = 10000)

        fgsea_result[["Dims"]] <- paste0("comp",i)

        pca_list[[i]] <- fgsea_result
      }
      pca <- do.call(rbind, pca_list)
      object@GSEA_list[[dim_reduction_name]] <- pca
    }

  }else{##protein_rank

    if((!is.null(fromType)) & (!is.null(toType))){##ID conversion

      if(!is.null(species)){
        ids <- names(protein_rank)

        protein_rank <- as.data.frame(protein_rank)
        protein_rank[[fromType]] <- ids
        colnames(protein_rank) <- c("value",fromType)

        id_conversion <- bitr(
          geneID = ids,
          fromType = fromType,
          toType = toType,
          OrgDb = species)

        mix_df <- full_join(protein_rank,id_conversion)
        mix_df[[fromType]] <- NULL
        mix_df <- mix_df %>% na.omit() %>% group_by(.data[[toType]]) %>%
          summarise(mean_value = mean(value))

        ranked_genes <- mix_df[["mean_value"]]
        names(ranked_genes) <- mix_df[[toType]]
        ranked_genes <- sort(ranked_genes,decreasing = F)

      }else{
        stop("You have not specified the name of the species-specific annotation package.")
      }

    }else{
      ranked_genes <- sort(protein_rank,decreasing = F)
    }

    fgsea_result <- fgsea(
      pathways = pathway,
      stats = ranked_genes,
      minSize = 1,
      maxSize = 500,
      nperm = 10000
    )
    object@GSEA_list[["non_dim_reduction"]] <- fgsea_result
  }
  return(object)
}
