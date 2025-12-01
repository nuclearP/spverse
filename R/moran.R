#' @title moran_sp
#' @description
#' Calculate the Moran's I of proteins or biological pathways
#' @param object an sp object \linkS4class{sp}.
#' @param ssgsea There are three options in total:
#' \itemize{ \item no :
#'   Default value, Calculate the Moran's I of features in the slot clean_data.
#'   \item clean : Calculate the Moran's I of the matrix named ssgsea_clean in ssgsea_list.
#'   \item RID : Calculate the Moran's I of the matrix named ssgsea_RID in ssgsea_list.}
#' @param scaled a logical indicating whether the Moran's I should be scaled so that it varies between -1 and +1  (default to TRUE).
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param alternative a character string specifying the alternative hypothesis that is tested
#' against the null hypothesis of no phylogenetic correlation; must be of one "two.sided",
#' "less", or "greater", or any unambiguous abbrevation of these.
#' @returns A data frame contains various elements. See also \code{\link[ape]{Moran.I}}
#' @export
#' @rdname moran_sp
moran_sp <- function(object,ssgsea="no",scaled = T,na.rm=T,alternative ="two.sided"){

  ssgsea <- match.arg(ssgsea, c("no","clean","RID"))

  if(ssgsea == "no"){
    df <- object@clean_data
  } else if(ssgsea == "clean"){
    df <- object@ssgsea_list[["ssgsea_clean"]] %>% as.data.frame()
  } else{
    df <- object@ssgsea_list[["sgsea_RID"]] %>% as.data.frame()
  }

  if(is.null(object@adj_matrix_list)){
    stop("Spatial adjacency matrix is missing.")
  }

  individual <- names(object@adj_matrix_list)

  moran_list <- list()

  for(i in individual){

    adj_matrix <- object@adj_matrix_list[[i]]

    df_new <- df[,colnames(adj_matrix)]

    df_new <- dplyr::select(df,colnames(adj_matrix))

    df_new <- t(df_new) %>% as.data.frame()

    moran_result <- map_df(df_new,~Moran.I(.x,adj_matrix,scaled = scaled,na.rm=na.rm,alternative = alternative))
    moran_result <- as.data.frame(moran_result)
    moran_result[["p"]] <- colnames(df_new)
    moran_result[["samples"]] <- i
    moran_list[[i]] <- moran_result
  }
  moran_list <- do.call(rbind, moran_list)

  return(moran_list)
}
