#' @title pca_var
#' @description
#' Visualize the correlation between features and different principal components using ggplot2.
#' @param x An prcomp object of class PCA.
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param n The number of features to display, including n positively correlated
#' features and n negatively correlated features.
#' @return ggplot object.
#' @export
#' @rdname pca_var
pca_var <- function(x,pc = c(1,2),n=10){

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
#' Visualize the correlation between features and different principal components using ggplot2 for the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param pc A vector containing two numerical values, used to select two principal
#' components for display, with the first and second dimensions as the default.
#' @param n The number of features to display, including n positively correlated
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Visualize the PCA results of the slot clean_data.
#'   \item RID : Visualize the PCA results of the slot RID_bootstrap.}
#' @return ggplot object.
#' @export
#' @seealso \code{\link{pca_var}}
#' @rdname pca_var_sp
pca_var_sp <- function(object,pc=c(1,2),n=10,name="clean"){

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){

    if(is.null(object@dim_reductions[["pca_clean"]])){
      stop("The slot pca_clean is missing in the sp object.")
    }
    p <- pca_var(object@dim_reductions[["pca_clean"]],pc=pc,n=n)

  }else{

    if(is.null(object@dim_reductions[["pca_RID"]])){
      stop("The slot pca_RID is missing in the sp object.")
    }
    p <- pca_var(object@dim_reductions[["pca_RID"]],pc=pc,n=n)

  }
  return(p)
}

