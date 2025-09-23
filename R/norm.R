#' @title normlization_sp
#' @description normlize the cleaned_data slot of the sp object by Quantiles.
#' @importFrom limma normalizeQuantiles
#' @param object An \linkS4class{sp}.
#' @return A dataframe contains data that has been quality controlled and normalized.
#' @name normlization_sp
#' @rdname normlization_sp
#' @export
setGeneric("normlization_sp",function(object) standardGeneric("normlization_sp") )

#' @rdname normlization_sp
#' @exportMethod normlization_sp
setMethod("normlization_sp","sp",function(object){

  if(nrow(object@cleaned_data) < 1){
    stop("The cleaned_data slot is missing in the sp object.")
  }
  object@cleaned_data <- object@cleaned_data %>% normalizeQuantiles()

  return(object)
}
)
