#' @title ppd (Probability of protein detection).
#' @description
#' Calculate the detection rate of proteins in different groups.
#' @param x A protein expression matrix (data.frame format), where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param group Group information of samples, which is a character vector or factor.
#' @return A data frame containing the detection rate of different proteins in each group.
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
ppd <- function(x,group) {

  if(dim(x)[2] != length(group)){
    stop("The grouping information you provided is incorrect.")
  }

  x <- t(x) %>% as.data.frame()
  x[["group"]] = group

  result <- x %>%
    group_by(group) %>%
    summarise(across(everything(), ~ mean(!is.na(.))))

  return(result)
}

#' @title ppd_sp
#' @description
#' Calculate the detection rate of proteins in different groups for sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @returns A data frame containing the detection rate of different proteins in each group.
#' @export
#' @rdname ppd_sp
#' @seealso \code{\link{ppd}}
ppd_sp <- function(object,group) {
  result <- ppd(object@rawdata,object@sample_features[[group]])
  return(result)
}
