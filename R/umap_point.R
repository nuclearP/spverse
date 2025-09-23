#' @title umap_point
#' @description
#' Create a ggplot2-based scatter plot for Uniform Manifold Approximation and Projection (UMAP) dimensional reduction.
#' @param x object of class umap.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname umap_point
umap_point <- function(x,...){
  UseMethod("umap_point")
}

#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggalt geom_encircle
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_discrete
#' @importFrom ggplot2 guide_legend
#' @importFrom tidydr theme_dr
#' @importFrom ggplot2 scale_color_discrete
#' @importFrom ggplot2 element_blank
#' @param group Group information of samples, which is a character vector or factor.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and NULL.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @method umap_point umap
#' @export
#' @rdname umap_point
umap_point.umap <- function(x,group,pal=NULL,shape="chull",legend_row=1,...){

  if(is.null(group)){
    stop("Data grouping information is missing.")
  }
  if(!is.factor(group)){
    group <- factor(group,levels = unique(group))
  }
  if(!is.null(pal) & (length(unique(group)) != length(pal))){
    stop("The number of variable categories you entered does not match the number of colors.")
  }
  shape <- match.arg(shape, c("chull",NULL))
  if(is.null(pal)){
    warning("The default color scheme will be used.")
  }

  umap_df <- x[["layout"]] %>%
    as.data.frame()%>%
    mutate(sample = group)

  if(shape == "chull"){
    p <- ggplot(umap_df,aes(x = V1,y = V2,color = sample))+
      geom_point(size = 1.5)+
      geom_encircle(aes(group = sample,fill = sample),expand=0,spread=0.5,s_shape=0.9,color = "gray50",alpha = 0.25,show.legend = F)+
      # theme_classic()+
      labs(x = "UMAP1",
           y = "UMAP2",
           color="Group")
  }else{
    p <- ggplot(umap_df,aes(x = V1,y = V2,color = sample))+
      geom_point(size = 1.5)+
      # theme_classic()+
      labs(x = "UMAP1",
           y = "UMAP2",
           color="Group")
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
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname umap_point_sp
setGeneric("umap_point_sp",function(object,...) standardGeneric("umap_point_sp") )

#' @param group Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and ellipse.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @param name There are two options in total:
#' \itemize{ \item RAW :
#'   Visualize the umap results of the original data.
#'   \item RID : Visualize the umap results of the data with individual differences removed.}
#' @seealso \code{\link{umap_point.umap}}
#' @rdname umap_point_sp
#' @exportMethod umap_point_sp
setMethod(
  f ="umap_point_sp",
  signature = signature(object= "sp"),
  function(object,group,pal=NULL,shape="chull",legend_row=1,name="RAW") {

    name <- match.arg(name, c("RAW","RID"))

    if(name == "RAW"){

      if(is.null(object@dim_reductions[["umap_RAW"]])){
        stop("The umap_RAW is missing in the sp object.")
      }

      sam = object@sample_features[c("samples",group)]

      sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["umap_RAW"]][["layout"]]))

      p <- umap_point(object@dim_reductions[["umap_RAW"]],group = sam[[group]],pal=pal,shape=shape,legend_row=legend_row)
    }else{

      if(is.null(object@dim_reductions[["umap_RID"]])){
        stop("The pca_RID slot is missing in the sp object.")
      }

      sam = object@sample_features[c("samples",group)]

      sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["umap_RID"]][["layout"]]))

      p <- umap_point(object@dim_reductions[["umap_RID"]],group = sam[[group]],pc=pc,pal=pal,shape=shape,legend_row=legend_row)
    }
    return(p)
  }
)

