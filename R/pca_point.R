#' @title pca_point
#' @description
#' Create a ggplot2-based scatter plot of PCA outputs prcomp.
#' @param x an prcomp object of class PCA.
#' @param group Group information of samples, which is a character vector or a factor.
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and ellipse.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @return ggplot object.
#' @export
#' @rdname pca_point
pca_point <- function(x,group,pc=c(1,2),pal=NULL,shape="chull",legend_row=1){

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
  }else if(shape == "ellipse"){
    p.pca1 <- ggplot(data = df1,aes(x = df1[,1],y = df1[,2],color=group))+
      stat_ellipse(geom = "polygon",level = 0.95,
                   linetype = 2,size=0.5,
                   aes(fill=group),
                   alpha=0.2,
                   show.legend = F)+
      geom_point(size = 1.5)+
      labs(x = xlab1,y = ylab1,color = "Group",fill="Group")+
      theme_classic()+
      theme(legend.position = "bottom")
  }else{
    p.pca1 <- ggplot(data = df1,aes(x = df1[,1],y = df1[,2],color=group))+
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
#' @param object an sp object \linkS4class{sp}.
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param shape There are two options available: chull (convex hull) and ellipse.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Visualize the PCA results of the slot clean_data.
#'   \item RID : Visualize the PCA results of the slot RID_bootstrap.}
#' @return ggplot object.
#' @export
#' @seealso \code{\link{pca_point}}
#' @rdname pca_point_sp
pca_point_sp <- function(object,group,pc=c(1,2),pal=NULL,shape="chull",legend_row=1,name="clean") {

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(is.null(object@dim_reductions[["pca_clean"]])){
      stop("The pca_clean is missing in the sp object.")
    }
    sam = object@sample_features[c("samples",group)]
    sam <- dplyr::filter(sam,sam[["samples"]] %in% row.names(object@dim_reductions[["pca_clean"]][["x"]]))
    p <- pca_point(object@dim_reductions[["pca_clean"]],group = sam[[group]],pc=pc,pal=pal,shape=shape,legend_row=legend_row)
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

