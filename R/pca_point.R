#' @title pca_point
#' @description
#' Create a ggplot2-based scatter plot of PCA outputs prcomp.
#' @param x an prcomp object of class PCA.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname pca_point
pca_point <- function(x,...){
  UseMethod("pca_point")
}

#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggalt geom_encircle
#' @importFrom ggplot2 stat_ellipse
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_color_discrete
#' @importFrom ggplot2 guide_legend
#' @param group Group information of samples, which is a character vector or factor.
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and ellipse.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @method pca_point prcomp
#' @export
#' @rdname pca_point
pca_point.prcomp <- function(x,group,pc=c(1,2),pal=NULL,shape="chull",legend_row=1,...){

  if(is.null(group)){
    stop("Data grouping information is missing.")
  }
  if(length(pc)!=2 & !is.numeric(pc)){
    stop("Please enter a vector containing two numerical values to select the principal components for plotting.")
  }
  if(!is.factor(group)){
    group <- factor(group,levels = unique(group))
  }
  if(!is.null(pal) & (length(unique(group)) != length(pal))){
    stop("The number of variable categories you entered does not match the number of colors.")
  }
  shape <- match.arg(shape, c("chull","ellipse"))
  if(is.null(pal)){
    warning("The default color scheme will be used.")
  }

  summ1 <- summary(x)
  xlab1 <- paste0("PC",pc[1],"(",round(summ1$importance[2,pc[1]]*100,2),"%)")
  ylab1 <- paste0("PC",pc[2],"(",round(summ1$importance[2,pc[2]]*100,2),"%)")

  df1 <- x[["x"]]
  df1 <- as.data.frame(df1) %>% dplyr::select(pc)
  df1 <- dplyr::mutate(df1,group = group)

  if(shape == "chull"){
    p.pca1 <- ggplot(data = df1,aes(x = df1[,1],y = df1[,2],color=group))+
      geom_encircle(aes(group = group,fill = group),expand=0,spread=0.5,s_shape=0.9,color = "gray50",alpha = 0.25,show.legend = F)+
      geom_point(size = 1.5)+
      labs(x = xlab1,y = ylab1,color = "Group")+
      theme_classic()+
      theme(legend.position = "bottom")
  }else{
    p.pca1 <- ggplot(data = df1,aes(x = df1[,1],y = df1[,2],color=group))+
      stat_ellipse(aes(fill = group),
                   type = "norm",geom = "polygon",alpha = 0.25,level= 0.1,show.legend = F)+ # 添加置信椭圆
      geom_point(size = 1.5)+
      labs(x = xlab1,y = ylab1,color = "Group",fill="Group")+
      theme_classic()+
      theme(legend.position = "bottom")
  }

  if(!is.null(pal)){
    p.pca1 <- p.pca1 +
      scale_color_manual(values = pal,guide = guide_legend(nrow = legend_row))+
      scale_fill_manual(values = pal)
  }else{
    p.pca1 <- p.pca1 +
      scale_color_discrete(guide = guide_legend(nrow = legend_row))
  }

  return(p.pca1)
}

#' @title pca_point_sp
#' @description
#' Create a ggplot2-based scatter plot of PCA for sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname pca_point_sp
setGeneric("pca_point_sp",function(object,...) standardGeneric("pca_point_sp") )

#' @param group Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and ellipse.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @param name There are two options in total:
#' \itemize{ \item RAW :
#'   Visualize the PCA results of the original data.
#'   \item RID : Visualize the PCA results of the data with individual differences removed.}
#' @seealso \code{\link{pca_point.prcomp}}
#' @rdname pca_point_sp
#' @exportMethod pca_point_sp
setMethod(
  f ="pca_point_sp",
  signature = signature(object= "sp"),
  function(object,group,pc=c(1,2),pal=NULL,shape="chull",legend_row=1,name="RAW") {

    name <- match.arg(name, c("RAW","RID"))

    if(name == "RAW"){

      if(is.null(object@dim_reductions[["pca_RAW"]])){
        stop("The pca_RAW is missing in the sp object.")
      }

      sam = object@sample_features[c("samples",group)]

      sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["pca_RAW"]][["x"]]))

      p <- pca_point(object@dim_reductions[["pca_RAW"]],group = sam[[group]],pc=pc,pal=pal,shape=shape,legend_row=legend_row)
    }else{

      if(is.null(object@dim_reductions[["pca_RID"]])){
        stop("The pca_RID slot is missing in the sp object.")
      }

      sam = object@sample_features[c("samples",group)]

      sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["pca_RID"]][["x"]]))

      p <- pca_point(object@dim_reductions[["pca_RID"]],group = sam[[group]],pc=pc,pal=pal,shape=shape,legend_row=legend_row)
    }
    return(p)
  }
)

