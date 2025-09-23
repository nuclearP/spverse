#' @title run_harmony_sp
#' @description
#' Performs the Harmony algorithm directly on the given sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @export
#' @rdname run_harmony_sp
setGeneric("run_harmony_sp",function(object,...) standardGeneric("run_harmony_sp"))

#' @importFrom dplyr filter
#' @importFrom harmony RunHarmony
#' @param dimn Select the expression matrix of the first n dimensions from the
#' principal component analysis result for harmony.
#' @param group Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @param theta parameters for \code{\link[harmony]{RunHarmony.default}}
#' @param sigma parameters for \code{\link[harmony]{RunHarmony.default}}
#' @rdname run_harmony_sp
#' @exportMethod run_harmony_sp
setMethod(
  f ="run_harmony_sp",
  signature = signature(object= "sp"),
  function(object,dimn,group,theta = 2.5,sigma = 1.2) {

    if(is.null(object@dim_reductions[["pca_RAW"]])){
      stop("The pca_RAW is missing in the sp object.")
    }

    df1 <- object@dim_reductions[["pca_RAW"]][["x"]] # 提取PC score

    df1 <- df1[,1:dimn]

    sam = object@sample_features[c("samples",group)]

    sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["pca_RAW"]][["x"]]))

    meta <- sam[[group]]

    names(meta) <- row.names(df1)

    set.seed(1)

    fi <- RunHarmony(df1,meta=meta,theta = theta,sigma = sigma) %>% t() %>% as.data.frame()

    object@corrective_PCA <- fi

    return(object)
  }
)
