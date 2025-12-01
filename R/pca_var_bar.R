#' @title pca_bar_sp
#' @description
#' Visualize the variances of pca dimensions for the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param addlabels logical value. If TRUE, labels are added at the top of bars or
#' points showing the information retained by each dimension.
#' @param ncp a numeric value specifying the number of dimensions to be shown.
#' The attribute of this column should be a character vector or a factor.
#' @param ggtheme ggplot2 theme name.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Visualize the PCA results of the slot clean_data.
#'   \item RID : Visualize the PCA results of the slot RID_bootstrap.}
#' @return ggplot object.
#' @export
#' @rdname pca_bar_sp
pca_bar_sp <- function(object,addlabels = TRUE,ncp = 15,ggtheme = theme_classic(),name = "clean") {

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(is.null(object@dim_reductions[["pca_clean"]])){
      stop("The pca_clean is missing in the sp object.")
    }
    p <- fviz_eig(object@dim_reductions[["pca_clean"]], addlabels =addlabels,ncp = ncp,ggtheme = ggtheme)
    p <- p+labs(title = NULL)
  }else{
    if(is.null(object@dim_reductions[["pca_RID"]])){
      stop("The pca_RID is missing in the sp object.")
    }
    p <- fviz_eig(object@dim_reductions[["pca_RID"]], addlabels = addlabels,ncp = ncp,ggtheme = ggtheme)
    p <- p+labs(title = NULL)
  }
  return(p)
}
