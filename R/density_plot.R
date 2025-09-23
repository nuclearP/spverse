#' @title density_plot
#' @description
#' Plot the probability density plot of the protein expression matrix.
#' @param x A protein expression matrix, where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname density_plot
density_plot <- function(x,...) {
  UseMethod("density_plot")
}

#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @importFrom tidyr drop_na
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @method density_plot data.frame
#' @export
#' @rdname density_plot
density_plot.data.frame <- function(x,...) {

  c <- tidyr::pivot_longer(x,everything())

  c <- tidyr::drop_na(c)

  p <- ggplot(c,mapping = aes(x=value,color=name))+
    geom_density()+
    theme_classic()+
    theme(legend.position = "none")+
    labs(y = "Density",x=NULL)

  return(p)

}

#' @title density_plot_sp
#' @description
#' Plot the probability density plot of sp object.
#' @param object An \linkS4class{sp}.
#' @return ggplot object.
#' @export
#' @rdname density_plot_sp
setGeneric("density_plot_sp",function(object) standardGeneric("density_plot_sp") )

#' @rdname density_plot_sp
#' @exportMethod density_plot_sp
setMethod(
  f ="density_plot_sp",
  signature = signature(object= "sp"),
  function(object) {
    if(nrow(object@cleaned_data) < 1){
      stop("The cleaned_data slot is missing in the sp object.")
    }
    density_plot(object@cleaned_data)
  }
)
