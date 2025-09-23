#' @title pca_bar_sp
#' @description
#' Visualize the variances of pca dimensions for the given sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @export
#' @rdname pca_bar_sp
setGeneric("pca_bar_sp",function(object,...) standardGeneric("pca_bar_sp"))

#' @importFrom factoextra fviz_eig
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 labs
#' @param addlabels logical value. If TRUE, labels are added at the top of bars or
#' points showing the information retained by each dimension.
#' @param ncp a numeric value specifying the number of dimensions to be shown.
#' The attribute of this column should be a character vector or a factor.
#' @param ggtheme function, ggplot2 theme name.
#' @param name There are two options in total:
#' \itemize{ \item RAW :
#'   Visualize the PCA results of the original data.
#'   \item RID : Visualize the PCA results of the data with individual differences removed.}
#' @rdname pca_bar_sp
#' @exportMethod pca_bar_sp
setMethod(
  f ="pca_bar_sp",
  signature = signature(object= "sp"),
  function(object,addlabels = TRUE,ncp = 15,ggtheme = theme_classic(),name = "RAW") {

    name <- match.arg(name, c("RAW","RID"))

    if(name == "RAW"){

      if(is.null(object@dim_reductions[["pca_RAW"]])){
        stop("The pca_RAW is missing in the sp object.")
      }

      p <- fviz_eig(object@dim_reductions[["pca_RAW"]], addlabels =addlabels,ncp = ncp,ggtheme = ggtheme)

      p <- p+labs(title = NULL)

    }else{

      if(is.null(object@dim_reductions[["pca_RID"]])){
        stop("The pca_RID slot is missing in the sp object.")
      }

      p <- fviz_eig(object@dim_reductions[["pca_RID"]], addlabels = addlabels,ncp = ncp,ggtheme = ggtheme)

      p <- p+labs(title = NULL)
    }
    return(p)
  }
)
