#' @title mvi
#' @description
#' Impute missing values of protein expression matrix.
#' @param x A protein expression matrix, where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param ... Arguments passed to other methods.
#' @export
#' @return A protein expression matrix with missing values imputed.
#' @rdname mvi
mvi <- function(x,...) {
  UseMethod("mvi")
}

#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom rrcovNA impSeqRob
#' @importFrom rrcovNA impSeq
#' @importFrom DreamAI impute.RegImpute
#' @importFrom SeqKnn SeqKNN
#' @import glmnet
#' @param method Methods for missing value imputation: There are four methods in total,
#' namely impseqrob (Robust sequential imputation), impseq (Sequential imputation),
#' grr (Glmnet Ridge Regression), and seqknn (Sequential knn).
#' where the row names are protein or gene IDs and the column names are sample names.
#' @param nmr Non-missing rate. Using the non-missing rate as the threshold,
#' filter out proteins below this value before performing missing value imputation.
#' @method mvi data.frame
#' @export
#' @rdname mvi
mvi.data.frame <- function(x,method="impseqrob",nmr=0,...) {

  method <- match.arg(method, c("impseqrob", "impseq", "grr", "seqknn"))

  if(nmr < 0 || nmr > 1){
    stop("ir should be greater than 0 and less than 1.")
  }
  if(nmr == 1){
    stop("The value range of parameter mvr should be less than 1.")
  }
  if(0.5 < nmr & nmr < 1){
    warning("You have filtered out too many proteins.")
  }

  nmrp <- rowMeans(!is.na(x))

  nmrp <- nmrp[nmrp > nmr]

  x[["p"]] <- row.names(x)

  x <- x %>% dplyr::filter(p %in% names(nmrp)) %>% dplyr::select(-p)

  if(method == "impseqrob"){
    bridge<-impSeqRob(x,alpha=0.9)
    x <- bridge[["x"]]
  }
  if(method == "impseq"){
    x<-impSeq(x)
  }
  if(method == "grr"){
    x<-impute.RegImpute(data=as.matrix(x), fillmethod = "row_mean", maxiter_RegImpute = 10,conv_nrmse = 1e-03)
  }
  if(method == "seqknn"){
    x<- SeqKNN(x,k = 10)
  }

  x <- as.data.frame(x)

  return(x)

}

#' @title mvi_sp
#' @description
#' Impute missing values of sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @rdname mvi_sp
#' @export
setGeneric("mvi_sp",function(object,...) standardGeneric("mvi_sp"))

#' @returns The sp object obtained after performing missing value imputation on cleaned_data.
#' @param method Methods for missing value imputation: There are four methods in total,
#' namely impseqrob (Robust sequential imputation), impseq (Sequential imputation),
#' grr (Glmnet Ridge Regression), and seqknn (Sequential knn).
#' where the row names are protein or gene IDs and the column names are sample names.
#' @param nmr Non-missing rate. Using the non-missing rate as the threshold,
#' filter out proteins below this value before performing missing value imputation.
#' @seealso \code{\link{mvi.data.frame}}
#' @rdname mvi_sp
#' @exportMethod mvi_sp
setMethod(
  f ="mvi_sp",
  signature = signature(object= "sp"),
  function(object,method="impseqrob",nmr=0) {
    if(nrow(object@cleaned_data) < 1){
      stop("The cleaned_data slot is missing in the sp object.")
    }
    object@cleaned_data <- mvi(object@cleaned_data,method = method,nmr=nmr)

    return(object)

  }
)
