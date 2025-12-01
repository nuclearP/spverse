#' @title umap_point
#' @description
#' Create a ggplot2-based scatter plot for Uniform Manifold Approximation and Projection (UMAP) dimensional reduction.
#' @param x object of class umap.
#' @param group Group information of samples, which is a character vector or factor.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and NULL.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @return ggplot object.
#' @export
#' @rdname umap_point
umap_point <- function(x,group,pal=NULL,shape=NULL,legend_row=1){

  if(is.null(group)){
    stop("Data grouping information is missing.")
  }
  if(!is.factor(group)){
    group <- factor(group,levels = unique(group))
  }
  if(!is.null(pal) & (length(unique(group)) != length(pal))){
    stop("The number of variable categories you entered does not match the number of colors.")
  }
  if(is.null(pal)){
    warning("The default color scheme will be used.")
  }

  umap_df <- x[["layout"]] %>%
    as.data.frame()%>%
    mutate(sample = group)

  if(is.null(shape)){
    p <- ggplot(umap_df,aes(x = V1,y = V2,color = sample))+
      geom_point(size = 1.5)+
      # theme_classic()+
      labs(x = "UMAP1",
           y = "UMAP2",
           color="Group")
  }else if(shape == "chull") {
    p <- ggplot(umap_df,aes(x = V1,y = V2,color = sample))+
      geom_point(size = 1.5)+
      geom_encircle(aes(group = sample,fill = sample),expand=0,spread=0.5,s_shape=0.9,color = "gray50",alpha = 0.25,show.legend = F)+
      # theme_classic()+
      labs(x = "UMAP1",
           y = "UMAP2",
           color="Group")
  }else{
    stop("You have entered an incorrect shape parameter.")
  }

  if(!is.null(pal)){
    p <- p +
      theme_dr(xlength = 0.2,
               ylength = 0.2,
               arrow = grid::arrow(length = unit(0.1, "inches"),
                                   ends = 'last', type = "closed"))+
      theme(panel.grid = element_blank())+
      theme(legend.position = "bottom")+
      scale_color_manual(values = pal,guide = guide_legend(nrow = legend_row))+
      scale_fill_manual(values = pal)
  }else{
    p <- p +
      theme_dr(xlength = 0.2,
               ylength = 0.2,
               arrow = grid::arrow(length = unit(0.1, "inches"),
                                   ends = 'last', type = "closed"))+
      theme(panel.grid = element_blank())+
      theme(legend.position = "bottom")+
      scale_color_discrete(guide = guide_legend(nrow = legend_row))
  }

  return(p)
}

#' @title umap_point_sp
#' @description
#' Create a ggplot2-based scatter plot of umap for sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and ellipse.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @param name There are three options in total:
#' \itemize{ \item clean :
#'   TVisualize the umap results of the original data. (the umap_clean).
#'   \item RID : Visualize the umap results of the data with individual differences removed through the bootstrap method. (umap_RID)
#'    \item harmony : Visualize the umap results of the data with individual differences removed through the harmony method. (umap_harmony)}
#' @return ggplot object.
#' @export
#' @rdname umap_point_sp
#' @seealso \code{\link{umap_point}}
umap_point_sp <- function(object,group,pal=NULL,shape=NULL,legend_row=1,name="clean") {

  name <- match.arg(name, c("clean","RID","harmony"))

  if(name == "clean"){
    if(is.null(object@dim_reductions[["umap_clean"]])){
      stop("The umap_clean is missing in the sp object.")
    }
    sam = object@sample_features[c("samples",group)]
    sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["umap_clean"]][["layout"]]))
    p <- umap_point(object@dim_reductions[["umap_clean"]],group = sam[[group]],pal=pal,shape=shape,legend_row=legend_row)
  }else if(name == "RID"){
    if(is.null(object@dim_reductions[["umap_RID"]])){
      stop("The ump_RID slot is missing in the sp object.")
    }
    sam = object@sample_features[c("samples",group)]
    sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["umap_RID"]][["layout"]]))
    p <- umap_point(object@dim_reductions[["umap_RID"]],group = sam[[group]],pal=pal,shape=shape,legend_row=legend_row)
  }else{
    if(is.null(object@dim_reductions[["umap_harmony"]])){
      stop("The umap_harmony slot is missing in the sp object.")
    }
    sam = object@sample_features[c("samples",group)]
    sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["umap_harmony"]][["layout"]]))
    p <- umap_point(object@dim_reductions[["umap_harmony"]],group = sam[[group]],pal=pal,shape=shape,legend_row=legend_row)
  }

  return(p)
}

