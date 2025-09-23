#' @title pca_var
#' @description
#' Visualize the correlation between features and different principal components using ggplot2.
#' @param x an prcomp object of class PCA.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname pca_var
pca_var <- function(x,...){
  UseMethod("pca_var")
}

#' @importFrom factoextra get_pca_var
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_segment
#' @importFrom grid arrow
#' @importFrom grid unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_classic
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param n The number of features to display, including n positively correlated
#' features and n negatively correlated features.
#' @method pca_var prcomp
#' @export
#' @rdname pca_var
pca_var.prcomp <- function(x,pc = c(1,2),n=10,...){

  pca_var <- get_pca_var(x)

  correlation <- pca_var[["cor"]]

  correlation <- as.data.frame(correlation)

  any_false <- pc %in% c(1:ncol(correlation))

  any_false <- any(any_false == F)

  if(any_false){
    stop("The dimension you entered does not exist.")
  }

  df <- correlation[,pc]

  df["name"] <- row.names(df)

  if(nrow(df) < (n*2) ){
    message("The number of features is insufficient, so all features will be displayed.")
  }else{
    sorted_df <- df[order(df[,pc[1]]),]
    min_n <- head(sorted_df, n)
    max_n <- tail(sorted_df, n)
    df <- rbind(min_n,max_n)
  }

  p <- ggplot(df, aes(x = df[,1],y = df[,2])) +
    geom_segment(
      aes(x = 0, y = 0, xend = df[,1], yend = df[,2]),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      color = "gray50"
    ) +
    geom_text_repel(
      aes(label = name)
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    # 添加水平虚线 y=0
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    xlim(-1,1) +
    ylim(-1,1) +
    labs(
      x = paste0("PC",pc[1]),
      y = paste0("PC",pc[2])
    ) +
    theme_classic()

  return(p)
}

#' @title pca_var_sp
#' @description
#' Visualize the correlation between features and different principal components using ggplot2 for sp object.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname pca_var_sp
setGeneric("pca_var_sp",function(object,...) standardGeneric("pca_var_sp") )

#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param n The number of features to display, including n positively correlated
#' @param name There are two options in total:
#' \itemize{ \item RAW :
#'   Visualize the PCA results of the original data.
#'   \item RID : Visualize the PCA results of the data with individual differences removed.}
#' @seealso \code{\link{pca_var.prcomp}}
#' @rdname pca_var_sp
#' @exportMethod pca_var_sp
setMethod(
  f ="pca_var_sp",
  signature = signature(object= "sp"),
  function(object,pc=c(1,2),n=10,name="RAW"){

    name <- match.arg(name, c("RAW","RID"))

    if(name == "RAW"){
      if(is.null(object@dim_reductions[["pca_RAW"]])){
        stop("The pca_RAW is missing in the sp object.")
      }

      p <- pca_var(object@dim_reductions[["pca_RAW"]],pc=pc,n=n,name=name)

    }else{

      if(is.null(object@dim_reductions[["pca_RID"]])){
        stop("The pca_RID slot is missing in the sp object.")
      }

      p <- pca_var(object@dim_reductions[["pca_RID"]],pc=pc,n=n,name=name)

    }
    return(p)
  }

)

