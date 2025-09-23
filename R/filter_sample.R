#' @title filter_sample
#' @description
#' Filter out samples where the number of identified proteins
#' is less than (25th percentile minus n times IQR).
#' @param x A protein expression matrix, where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param ... Arguments passed to other methods.
#' @returns A data.frame with Some Columns (Samples) Filtered Out.
#' @export
#' @rdname filter_sample
filter_sample <- function(x,...) {
  UseMethod("filter_sample")
}

#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom magrittr set_names
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @param times The coefficient of the IQR value, with a default value of 1.
#' @param groups Group information of samples, which is a character vector or factor.
#' @method filter_sample data.frame
#' @export
#' @rdname filter_sample
filter_sample.data.frame <- function(x,times = 1,groups=NULL,...) {
  if(is.null(groups)){
    df <- colSums(!is.na(x)) %>%
      as.data.frame() %>%
      set_names("values") %>%
      dplyr::mutate(sample = colnames(x))

    filtered_df <- df %>%
      dplyr::filter(values >= (quantile(values, 0.25, na.rm = TRUE) -
                                 times * IQR(values, na.rm = TRUE)))

  }else{
    df <- colSums(!is.na(x)) %>%
      as.data.frame() %>%
      magrittr::set_names("values") %>%
      dplyr::mutate(sample = colnames(x),groups = groups)

    filtered_df <- df %>%
      dplyr::group_by(groups) %>%
      dplyr::filter(values >= (quantile(values, 0.25, na.rm = TRUE) -
                                 times * IQR(values, na.rm = TRUE)))

  }

  x <- dplyr::select(x,filtered_df$sample) %>% removeRowsAllNa()

}

#' @title filter_sample_sp
#' @description
#' Filter out samples of sp object where the number of identified proteins
#' is less than (25th percentile minus n times IQR).
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @returns Assign a value to the cleaned_data of the sp object.
#' @rdname filter_sample_sp
#' @export
setGeneric("filter_sample_sp",function(object,...) standardGeneric("filter_sample_sp") )

#' @param times The coefficient of the IQR value, with a default value of 1.
#' @param groups Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @seealso \code{\link{filter_sample.data.frame}}
#' @rdname filter_sample_sp
#' @exportMethod filter_sample_sp
setMethod(
  f ="filter_sample_sp",
  signature = signature(object= "sp"),
  function(object,times = 1,groups = NULL) {

    if(is.null(groups)){
      object@cleaned_data <- filter_sample(object@rawdata,times = times)
      return(object)
    }else{
      object@cleaned_data <- filter_sample(object@rawdata,times = times,object@sample_features[[groups]])
      return(object)
    }
  })
