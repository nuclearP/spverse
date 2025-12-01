#' @title removeRowsAllNa
#' @description Delete rows that are all NA.
#' @param object a data.frame or matrix.
#' @export
#' @examples
#' \dontrun{
#'  x <- matrix(c(rep(NA,3),rep(1,6)),nrow = 3)
#'  x <- removeRowsAllNa(x)
#' }
removeRowsAllNa  <- function(object){
  if (!is.data.frame(object) & !is.matrix(object)){
    stop("the input data must be data.frame or matrix")
  }
  object[apply(object, 1, function(object) any(!is.na(object))),]
}

removena_cv = function(x){
  y=na.omit(x)
  return(sd(y)/mean(y))
}

removena_sd = function(x){
  y=na.omit(x)
  return(sd(y))
}
