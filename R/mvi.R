#' @title mvi
#' @description
#' Impute missing values of protein expression matrix.
#' @param x A protein expression matrix (data.frame format), where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param method Methods for missing value imputation: There are four methods in total,
#' namely impseqrob (Robust sequential imputation), impseq (Sequential imputation),
#' grr (Glmnet Ridge Regression), and seqknn (Sequential knn).
#' where the row names are protein or gene IDs and the column names are sample names.
#' @param nmr Non-missing rate. Using the non-missing rate as the threshold,
#' filter out proteins below this value before performing missing value imputation.
#' @export
#' @return A protein expression matrix with missing values imputed.
#' @rdname mvi
mvi <- function(x,method="impseqrob",nmr=0) {

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

  x <- x[x$`p` %in% names(nmrp),]
  x <- x[, colnames(x) != "p"]

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
#' Impute missing values for sp object slot clean_data.
#' @param object an sp object \linkS4class{sp}.
#' @param method Methods for missing value imputation: There are four methods in total,
#' namely impseqrob (Robust sequential imputation), impseq (Sequential imputation),
#' grr (Glmnet Ridge Regression), and seqknn (Sequential knn).
#' where the row names are protein or gene IDs and the column names are sample names.
#' @param nmr Non-missing rate. Using the non-missing rate as the threshold,
#' filter out proteins below this value before performing missing value imputation.
#' @returns The sp object obtained after performing missing value imputation on clean_data.
#' @seealso \code{\link{mvi}}
#' @rdname mvi_sp
#' @export
mvi_sp <- function(object,method="impseqrob",nmr=0){

  if(nrow(object@clean_data) < 1){
    stop("The clean_data slot is missing in the sp object.")
  }

  object@clean_data <- mvi(object@clean_data,method = method,nmr=nmr)

  return(object)
}
