#' @title ppd (Probability of protein detection).
#' @description
#' Calculate the detection rate of proteins in different groups.
#' @param x A protein expression matrix, where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param ... Arguments passed to other methods.
#' @export
#' @return A data frame containing the detection rate of different proteins in each group.
#' @rdname ppd
ppd <- function(x,...) {
  UseMethod("ppd")
}

#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr across
#' @importFrom tidyselect everything
#' @param group information of samples, which is a character vector or factor.
#' @method ppd data.frame
#' @export
#' @rdname ppd
#' @examples
#' \dontrun{
#'  group = rep(c("groupA", "groupA", "gruupB"), each = 15)
#'
#'  set.seed(123)
#'  data <- data.frame(
#'  pro1 = sample(c(1:10, NA), 45, replace = TRUE),
#'  pro2 = sample(c(20:30, NA), 45, replace = TRUE),
#'  pro3 = sample(c(50:60, NA), 45, replace = TRUE)
#'  )
#'
#'  data <- t(data) %>% as.data.frame()
#'
#'  ppd(data,group)
#'
#' }
ppd.data.frame <- function(x,group,...) {

  if(dim(x)[2] != length(group)){
    stop("The grouping information you provided is incorrect.")
  }

  x <- t(x) %>% as.data.frame() %>% mutate(group = group)

  result <- x %>%
    group_by(group) %>%
    summarise(across(everything(), ~ mean(!is.na(.))))

  return(result)
}

#' @title ppd_sp
#' @description
#' calculate the detection rate of proteins in different groups for sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @returns A data frame containing the detection rate of different proteins in each group.
#' @export
#' @rdname ppd_sp
setGeneric("ppd_sp",function(object,...) standardGeneric("ppd_sp"))

#' @param group Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @seealso \code{\link{ppd.data.frame}}
#' @rdname ppd_sp
#' @exportMethod ppd_sp
setMethod(
  f ="ppd_sp",
  signature = signature(object= "sp"),
  function(object,group) {
    result <- ppd(object@rawdata,object@sample_features[[group]])
    return(result)
  })
