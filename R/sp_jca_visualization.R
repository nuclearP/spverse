#' @title sp_visualization
#' @description Visualize the global spatial proteome.
#' @param object an sp object \linkS4class{sp}.
#' @param variale The name of the variable to be displayed can be the name of a column
#' in the sample_features data frame, the name of a row in the clean_data data frame,
#' or the name of a pathway in the ssgsea results.
#' @param nrow Number of rows on the dividing surface.
#' @param ssgsea There are three options in total:
#' \itemize{ \item no :
#'   Default value, indicating that the variable to be visualized is not part of the ssgsea results.
#'   \item clean : The variable to be visualized is derived from the matrix named ssgsea_clean in ssgsea_list.
#'   \item RID : The variable to be visualized is derived from the matrix named ssgsea_RID in ssgsea_list.}
#' @param continuous_variable A logical value, defaulting to TRUE, which determines whether the variable
#' is continuous or discrete.
#' @param fills If the continuous_variable parameter is TRUE, input a vector containing three colors;
#' if continuous_variable is FALSE, input a number of colors corresponding to the number of categories
#' of the variable to be visualized.
#' @returns ggplot object.
#' @export
#' @rdname sp_visualization
sp_visualization <- function(object,variale,nrow=1,ssgsea="no",continuous_variable=T,fills = c("#2556A6","#ffffff","#EE2A29")){

  ssgsea <- match.arg(ssgsea, c("no","clean","RID"))

  sample_features <- object@sample_features

  if(ssgsea == "no"){
    pros <- object@clean_data %>% t() %>% as.data.frame()
    pros[["samples"]] <- row.names(pros)
    final <- full_join(sample_features,pros)
    final <- final[!is.na(final[[variale]]),]
  } else if(ssgsea == "clean"){
    pathways <- object@ssgsea_list[["ssgsea_clean"]] %>% t() %>% as.data.frame()
    pathways[["samples"]] <- row.names(pathways)
    pathways <- pathways[,c("samples",variale)]
    final <- full_join(sample_features,pathways)
    final <- final[!is.na(final[[variale]]),]
  } else{
    pathways <- object@ssgsea_list[["sgsea_RID"]] %>% t() %>% as.data.frame()
    pathways[["samples"]] <- row.names(pathways)
    pathways <- pathways[,c("samples",variale)]
    final <- full_join(sample_features,pathways)
    final <- final[!is.na(final[[variale]]),]
  }

  if(isTRUE(continuous_variable)){
    p <- ggplot(data =final,aes(x=x,y=y))+
      geom_tile(data =final,aes(x=x,y=y,fill = .data[[variale]]),color="black",size=0.8)+
      facet_wrap(~individual,nrow = nrow)+
      labs(fill=variale)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      theme(axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())+
      scale_fill_gradientn(colours =colorRampPalette(fills)(50))+
      coord_fixed()
  }else{
    p <- ggplot(data =final,aes(x=x,y=y))+
      geom_tile(data =final,aes(x=x,y=y,fill = .data[[variale]]),color="black",size=0.8)+
      facet_wrap(~individual,nrow = nrow)+
      labs(fill=variale)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      theme(axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())+
      scale_fill_manual(values = fills)+
      coord_fixed()
  }
  return(p)
}

#' @title jca_visualization
#' @description
#' Visualize the output results of the jca_sp function.
#' @param x The result returned by the jca_sp function can be a matrix or a list.
#' @param colors A vector containing three colors.
#' @param nrow Number of rows on the dividing surface. If x is a matrix, then this parameter is useless.
#' @export
#' @rdname jca_visualization
jca_visualization <- function(x,colors = c("#2556A6","#ffffff","#EE2A29"),nrow = 1){

  if(is.list(x) == F){

    x <- as.data.frame(x)
    x <- mutate(x,pair = row.names(x))
    x <- x[-nrow(x),]
    x <- dplyr::select(x,4,5)
    x <- separate(x,col = "pair",sep = ":",into = c("a","b"))

    p <- ggplot(x,mapping = aes(x=a,y=b,fill=`z-value`))+
      geom_tile(color="black",size=0.8)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      labs(x=NULL,y=NULL)+
      scale_fill_gradientn(colours =colorRampPalette(colors)(50))+
      coord_fixed()
  }else{

    list_name <- names(x)
    new_list <- list()

    for (i in list_name) {
      y <- as.data.frame(x[[i]])
      y[["pair"]] <- row.names(y)
      y <- y[-nrow(y),]
      y <- y[,c(4,5)]
      y <- separate(y,col = "pair",sep = ":",into = c("a","b"))
      y[["group"]] <- i
      new_list[[i]] <- y
    }

    y <- rbindlist(new_list)

    p <- ggplot(y,mapping = aes(x=a,y=b,fill=`z-value`))+
      geom_tile(color="black",size=0.8)+
      theme_bw()+
      theme(panel.grid = element_blank())+
      labs(x=NULL,y=NULL)+
      scale_fill_gradientn(colours =colorRampPalette(colors)(50))+
      facet_wrap(~group,nrow = nrow)+
      coord_fixed()
  }
  return(p)
}
