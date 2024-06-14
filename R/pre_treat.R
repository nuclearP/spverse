#' removeRowsAllNa
#'
#' Delete rows that are all NA
#'
#' @param object a data.frame or matrix
#' @export
#' @examples
#' \dontrun{
#'  x <- matrix(c(rep(NA,3),rep(1,6)),nrow = 3)
#'  x <- removeRowsAllNa(x)
#' }
removeRowsAllNa  <- function(object){
  if (!is.data.frame(object) || !is.matrix(object)){
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

#' remove_outliers
#'
#' remove missing values based on IQR principles
#'
#' @param x numeric vector
#' @param n the coefficient of IQR
#' @param na.rm logical; if true, any NA and NaN's are removed from x before the quantiles are computed
#' @param ... further arguments passed to or from other methods
#' @export
remove_outliers <- function(x,n,na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 3 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}

