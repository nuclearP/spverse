#' @include utilities.R
NULL

#' @title run_pca
#' @description
#' Performs a principal components analysis on the given protein matrix and
#' returns the results as an object of class prcomp.
#' @param x A protein expression matrix (data.frame format), where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param filter_by Select features with large fluctuations for principal components analysis,
#' The standard deviation (sd) or coefficient of variation (cv) can be used for feature selection.
#' @param rank A numeric value between 0 and 1. Select features with the largest fluctuations based on percentage.
#' By default, select the top 50\% of features for principal components analysis.
#' @details
#' If your protein expression matrix is in logarithmic form (e.g., after log2 transformation),
#' we recommend using the standard deviation as the screening criterion.
#' This is because the difference between data represents a relative change in fold (multiples),
#' rather than the absolute value of abundance.
#' @return returns a list with class "prcomp".
#' @rdname run_pca
run_pca <- function(x,filter_by = "sd",rank=0.5){

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

  set.seed(1)
  pca_fit <- x_fi %>% t() %>%
    scale() %>%
    t() %>% as.data.frame() %>% na.omit() %>% t() %>%
    prcomp()

  return(pca_fit)
}

#' @title run_pca_sp
#' @description
#' Performs a principal components analysis on the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param filter_by Select features with large fluctuations for principal components analysis,
#' The standard deviation (sd) or coefficient of variation (cv) can be used for feature selection.
#' @param rank A numeric value between 0 and 1. Select features with the largest fluctuations based on percentage.
#' By default, select the top 50\% of features for principal components analysis.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Perform PCA dimensionality reduction using the slot clean_data, and store the returned result in "dim_reductions[["pca_clean"]]".
#'   \item RID : Perform PCA dimensionality reduction using the slot RID_bootstrap, and store the returned result in "dim_reductions[["pca_RID"]]".}
#' @seealso \code{\link{run_pca}}
#' @export
#' @rdname run_pca_sp
run_pca_sp <- function(object,filter_by = "sd",rank=0.5,name="clean") {

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
    object@dim_reductions[["pca_clean"]] <- run_pca(object@clean_data,filter_by = filter_by,rank=rank)
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
    object@dim_reductions[["pca_RID"]] <- run_pca(object@RID_bootstrap,filter_by = filter_by,rank=rank)
  }
  return(object)
}
