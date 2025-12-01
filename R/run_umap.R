#' @title run_umap
#' @description
#' Perform the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique.
#' @param x an prcomp object of class PCA.
#' @param dimn Select the expression matrix of the first n dimensions from the
#' pca result for UMAP dimensional reduction.
#' @return return object of class umap.
#' @export
#' @rdname run_umap
run_umap <- function(x,dimn){

  df <- x[["x"]]
  fi <- df[,1:dimn]

  set.seed(1)
  umap_fit <- fi %>%
    umap()

  return(umap_fit)
}

#' @title run_umap_sp
#' @description
#' Performs Uniform Manifold Approximation and Projection (UMAP) dimensional reduction on the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param dimn Select the expression matrix of the first n dimensions from the
#' pca result for UMAP dimensional reduction.
#' @param name There are three options in total:
#' \itemize{ \item clean :
#'   The PCA results of the original data (the slot clean_data).
#'   \item RID : The PCA results of the data with individual differences removed through the bootstrap method.
#'    \item harmony : The PCA results of the data with individual differences removed through the harmony method.}
#' @details
#' UMAP dimensionality reduction is based on PCA results.
#' @seealso \code{\link{run_umap}}
#' @export
#' @rdname run_umap_sp
run_umap_sp <- function(object,dimn,name="clean") {
  name <- match.arg(name, c("clean","RID","harmony"))
  if(name == "clean"){
    if(is.null(object@dim_reductions[["pca_clean"]])){
      stop("The pca_clean is missing in the sp object.")
    }
    umap_fit <- run_umap(object@dim_reductions[["pca_clean"]],dimn = dimn)
    object@dim_reductions[["umap_clean"]] <- umap_fit
  }else if (name == "RID"){
    if(is.null(object@dim_reductions[["pca_RID"]])){
      stop("The pca_RID slot is missing in the sp object.")
    }
    umap_fit <- run_umap(object@dim_reductions[["pca_RID"]],dimn = dimn)
    object@dim_reductions[["umap_RID"]] <- umap_fit
  }else{
    if(is.null(object@corrective_PCA)){
      stop("The corrective_PCA slot is missing in the sp object.")
    }
    set.seed(1234)
    umap_fit <- object@corrective_PCA %>% t() %>%
      umap()
    object@dim_reductions[["umap_harmony"]] <- umap_fit
  }
  return(object)
}
