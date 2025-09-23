#' @title run_umap
#' @description
#' Perform the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
#' @param x an prcomp object of class PCA.
#' @param ... Arguments passed to other methods.
#' @return return object of class umap.
#' @export
#' @rdname run_umap
run_umap <- function(x,...) {
  UseMethod("run_umap")
}

#' @importFrom umap umap
#' @param dimn Select the expression matrix of the first n dimensions from the
#' principal component analysis for UMAP dimensional reduction.
#' @method run_umap prcomp
#' @export
#' @rdname run_umap
run_umap.prcomp <- function(x,dimn,...){

  df1 <- x[["x"]]

  fi <- df1[,1:dimn]

  set.seed(1234)
  umap_fit <- fi %>%
    umap()

  return(umap_fit)
}

#' @title run_umap_sp
#' @description
#' Performs Uniform Manifold Approximation and Projection (UMAP) dimensional reduction on the given sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @export
#' @rdname run_umap_sp
setGeneric("run_umap_sp",function(object,...) standardGeneric("run_umap_sp"))

#' @param dimn Select the expression matrix of the first n dimensions from the
#' principal component analysis for UMAP dimensional reduction.
#' @param name There are two options in total:
#' \itemize{ \item RAW :
#'   Visualize the PCA results of the original data.
#'   \item RID : Visualize the PCA results of the data with individual differences removed.}
#' @seealso \code{\link{run_umap.prcomp}}
#' @rdname run_umap_sp
#' @exportMethod run_umap_sp
setMethod(
  f ="run_umap_sp",
  signature = signature(object= "sp"),
  function(object,dimn,name="RAW") {
    name <- match.arg(name, c("RAW","RID"))

    if(name == "RAW"){

      if(is.null(object@dim_reductions[["pca_RAW"]])){
        stop("The pca_RAW is missing in the sp object.")
      }

      umap_fit <- run_umap(object@dim_reductions[["pca_RAW"]],dimn = dimn)

      object@dim_reductions[["umap_RAW"]] <- umap_fit

    }else{

      if(is.null(object@dim_reductions[["pca_RID"]])){
        stop("The pca_RID slot is missing in the sp object.")
      }

      umap_fit <- run_umap(object@dim_reductions[["pca_RID"]],dimn = dimn)

      object@dim_reductions[["umap_RID"]] <- umap_fit

    }
    return(object)
  }
)
