#' @title density_plot
#' @description
#' Plot the probability density plot of the protein expression matrix.
#' @param x A protein expression matrix (data.frame format), where the row names are protein
#' or gene IDs and the column names are sample names.
#' @return ggplot object.
#' @export
#' @rdname density_plot
density_plot <- function(x) {

  c <- stack(x)

  c <- na.omit(c)

  p <- ggplot(c,mapping = aes(x=values,color=ind))+
    geom_density()+
    theme_classic()+
    theme(legend.position = "none")+
    labs(y = "Density",x=NULL)

  return(p)
}

#' @title density_plot_sp
#' @description
#' Plot the probability density plot of sp object.
#' @param object an sp object \linkS4class{sp}.
#' @return ggplot object.
#' @export
#' @rdname density_plot_sp
density_plot_sp <- function(object) {
  if(nrow(object@clean_data) < 1){
    stop("The cleaned_data slot is missing in the sp object.")
  }
  density_plot(object@clean_data)
}
