#' the sp class
#'
#' class designed for one experiment
#'
#' @slot rawdata data.frame.
#' @slot features data.frame.
#' @slot scaledata data.frame.
#' @slot dim list.
#' @exportClass sp
sp <- setClass(
  "sp",
  slots = c(
    rawdata = "data.frame",
    features = "data.frame",
    scaledata = "data.frame",
    dim = "list"
  )
)

#' The spatials class
#'
#' An S4 object that stores multiple sp objects, but only one sp object is activated at the same time
#'
#' @slot experiment list.
#' @slot active character.
#' @exportClass spatilas
spatials <- setClass(
  "spatials",
  slots = c(
    experiment = "list",
    active = "character"
  )
)
