#' @include utilities.R
NULL

#' @title run_pca
#' @description
#' Performs a principal components analysis on the given protein matrix and
#' returns the results as an object of class prcomp.
#' @param x A protein expression matrix, where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param ... Arguments passed to other methods.
#' @return returns a list with class "prcomp".
#' @export
#' @rdname run_pca
run_pca <- function(x,...) {
  UseMethod("run_pca")
}

#' @importFrom dplyr mutate
#' @importFrom dplyr slice_max
#' @importFrom dplyr select
#' @importFrom tidyr drop_na
#' @importFrom dplyr mutate
#' @importFrom stats prcomp
#' @param filter_by Select features with large fluctuations for principal components analysis,
#' The standard deviation (sd) or coefficient of variation (cv) can be used for feature selection.
#' @param rank Select features with the largest fluctuations based on percentage.
#' By default, select the top 50\% of features for principal components analysis.
#' @details
#' If your protein expression matrix is in logarithmic form (e.g., after log2 transformation),
#' we recommend using the standard deviation as the screening criterion.
#' This is because the difference between data represents a relative change in fold (multiples),
#' rather than the absolute value of abundance.
#' @method run_pca data.frame
#' @export
#' @rdname run_pca
run_pca.data.frame <- function(x,filter_by = "sd",rank=0.5,...){

  filter_by <- match.arg(filter_by, c("sd", "cv"))

  if(filter_by == "sd"){

    filters = apply(x,1,removena_sd)

    x_fi <- mutate(x,sds = filters) %>%
      slice_max(sds,n=round((nrow(x)*rank))) %>%
      dplyr::select(-sds)
  }
  if(filter_by == "cv"){

    filters = apply(x,1,removena_cv)

    x_fi <- mutate(x,cvs = filters) %>%
      slice_max(cvs,n=round((nrow(x)*rank))) %>%
      dplyr::select(-cvs)
  }

  nrow_x_fi <- nrow(x_fi)

  message(nrow_x_fi," variables (proteins) will be used for principal components analysis.")

  set.seed(1234)
  pca_fit <- x_fi %>% t() %>%
    scale() %>%
    t() %>% as.data.frame() %>% drop_na() %>% t() %>%
    prcomp()

  return(pca_fit)
}

#' @title run_pca_sp
#' @description
#' Performs a principal components analysis on the given sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @export
#' @rdname run_pca_sp
setGeneric("run_pca_sp",function(object,...) standardGeneric("run_pca_sp"))

#' @param filter_by Select features with large fluctuations for principal components analysis,
#' The standard deviation (sd) or coefficient of variation (cv) can be used for feature selection.
#' @param rank Select features with the largest fluctuations based on percentage.
#' By default, select the top 50\% of features for principal components analysis.
#' @param name There are two options in total:
#' \itemize{ \item RAW :
#'   Perform PCA dimensionality reduction using cleaned_data lot, and store the returned result in "dim_reductions[["pca_RAW"]]".
#'   \item RID : Perform PCA dimensionality reduction using RID_bootstrap lot, and store the returned result in "dim_reductions[["pca_RID"]]".}
#' @seealso \code{\link{run_pca.data.frame}}
#' @rdname run_pca_sp
#' @exportMethod run_pca_sp
setMethod(
  f ="run_pca_sp",
  signature = signature(object= "sp"),
  function(object,filter_by = "sd",rank=0.5,name="RAW") {

    name <- match.arg(name, c("RAW","RID"))

    if(name == "RAW"){

      if(nrow(object@cleaned_data) < 1){
        stop("The cleaned_data slot is missing in the sp object.")
      }

      object@dim_reductions[["pca_RAW"]] <- run_pca(object@cleaned_data,filter_by = filter_by,rank=rank)

      return(object)

    }else{

      if(nrow(object@RID_bootstrap)  < 1){
        stop("The RID_bootstrap slot is missing in the sp object.")
      }

      object@dim_reductions[["pca_RID"]] <- run_pca(object@RID_bootstrap,filter_by = filter_by,rank=rank)
    }
    return(object)
  }
)
