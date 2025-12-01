#' @title run_harmony_sp
#' @description
#' Performs the Harmony algorithm directly on the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param dimn Select the expression matrix of the first n dimensions from the
#' principal component analysis result for harmony.
#' @param group Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @param theta parameters for \code{\link[harmony]{RunHarmony.default}}
#' @param sigma parameters for \code{\link[harmony]{RunHarmony.default}}
#' @return The sp object with the slot corrective_PCA.
#' @export
#' @rdname run_harmony_sp
run_harmony_sp <- function(object,dimn,group,theta = 2.5,sigma = 1.2) {

  if(is.null(object@dim_reductions[["pca_clean"]])){
    stop("The pca_clean is missing in the sp object.")
  }

  df <- object@dim_reductions[["pca_clean"]][["x"]] # 提取PC score

  df <- df[,1:dimn]

  sam = object@sample_features[c("samples",group)]

  sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["pca_clean"]][["x"]]))

  meta <- sam[[group]]

  names(meta) <- row.names(df)

  set.seed(1)

  fi <- RunHarmony(df,meta=meta,theta = theta,sigma = sigma) %>% t() %>% as.data.frame()

  object@corrective_PCA <- fi

  return(object)
}
