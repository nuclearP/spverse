#' @title normlization_sp
#' @description Normlize the slot clean_data of the sp object by Quantiles method.
#' @param object an sp object \linkS4class{sp}.
#' @rdname normlization_sp
#' @export
normlization_sp <- function(object){
  if(nrow(object@clean_data) < 1){
    stop("The clean_data slot is missing in the sp object.")
  }
  object@clean_data <- object@clean_data %>% normalizeQuantiles()
  return(object)
}
