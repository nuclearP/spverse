#' @title filter_sample
#' @description
#' Filter out samples where the number of identified proteins
#' is less than (25th percentile minus n times IQR).
#' @param x A protein expression matrix (data.frame format), where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param times The coefficient of the IQR value, with a default value of 1.
#' @param group Group information of samples, which is a character vector or factor.
#' @returns A data.frame with Some Columns (Samples) Filtered Out.
#' @export
#' @rdname filter_sample
filter_sample <- function(x,times = 1,group=NULL) {

  if(is.null(group)){
    df <- colSums(!is.na(x)) %>%
      as.data.frame() %>%
      magrittr::set_names("values") %>%
      dplyr::mutate(sample = colnames(x))

    filtered_df <- df %>%
      dplyr::filter(values >= (quantile(values, 0.25, na.rm = TRUE) -
                                 times * IQR(values, na.rm = TRUE)))
  }else{
    df <- colSums(!is.na(x)) %>%
      as.data.frame() %>%
      magrittr::set_names("values") %>%
      dplyr::mutate(sample = colnames(x),group = group)

    filtered_df <- df %>%
      dplyr::group_by(group) %>%
      dplyr::filter(values >= (quantile(values, 0.25, na.rm = TRUE) -
                                 times * IQR(values, na.rm = TRUE)))

  }

  x <- dplyr::select(x,filtered_df$sample) %>% removeRowsAllNa()

  return(x)
}

#' @title filter_sample_sp
#' @description
#' Filter out samples of sp object where the number of identified proteins
#' is less than (25th percentile minus n times IQR).
#' @param object an sp object \linkS4class{sp}.
#' @param times The coefficient of the IQR value, with a default value of 1.
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @returns Assign a value to the slot clean_data of the sp object.
#' @seealso \code{\link{filter_sample}}
#' @export
#' @rdname filter_sample_sp
filter_sample_sp <- function(object,times = 1,group = NULL) {
  if(is.null(group)){
    object@clean_data <- filter_sample(object@rawdata,times = times)
  }else{
    object@clean_data <- filter_sample(object@rawdata,times = times,object@sample_features[[group]])
  }
  return(object)
}
