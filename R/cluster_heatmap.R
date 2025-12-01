#' @title heatmap_sp
#' @description
#' Generate a heatmap to visualize the clustering outcomes of the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param name \itemize{ \item clean :
#'   Visualize the pca results from the clean_data.
#'   \item RID : Visualize the pca results of the data with individual differences removed through the bootstrap method (RID_bootstrap).
#'   \item harmony : Visualize the pca results of the data with individual differences removed through the harmony method (corrective_PCA).}
#' @param dimn Select the expression matrix of the first n dimensions from the
#' pca results for visualization.
#' @param cluster_name Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param col Colors corresponding to value. A list contains named vectors.
#' @param show_row_names Whether show row names.
#' @param show_column_names Whether show column names.
#' @return A Heatmap-class object.
#' @export
#' @rdname heatmap_sp
heatmap_sp <- function(object,name = "clean",dimn,cluster_name,col=NULL,show_row_names=F,show_column_names=F) {

  name <- match.arg(name, c("clean","RID","harmony"))

  if(name == "clean"){
    if(is.null(object@dim_reductions[["pca_clean"]])){
      stop("The pca_clean is missing in the sp object.")
    }
    fi <- object@dim_reductions[["pca_clean"]][["x"]]
    fi <- t(fi[,1:dimn])
  }else if(name == "RID"){
    if(is.null(object@dim_reductions[["pca_RID"]])){
      stop("The pca_RID is missing in the sp object.")
    }
    fi <- object@dim_reductions[["pca_RID"]][["x"]]
    fi <- t(fi[,1:dimn])
  }else{
    if(is.null(object@corrective_PCA)){
      stop("The corrective_PCA slot is missing in the sp object.")
    }
    fi <- object@corrective_PCA
  }

  sample_features <- object@sample_features
  sample_features <- sample_features[,c("individual",cluster_name)]
  sample_features <- na.omit(sample_features)

  if(is.null(col)){
    top_anno <- HeatmapAnnotation(
      samples = sample_features[["individual"]],
      cluster = sample_features[[cluster_name]]
    )
  }else{
    top_anno <- HeatmapAnnotation(
      samples = sample_features[["individual"]],
      cluster = sample_features[[cluster_name]],
      col = col
    )
  }

  set.seed(1)

  zhc <- Heatmap(fi,
                 name = "Intensity",
                 show_row_names = show_row_names,
                 show_column_names = show_column_names,
                 column_split = sample_features[[cluster_name]],
                 show_row_dend = F,
                 top_annotation = top_anno,
                 column_title =NULL,
                 column_dend_height=unit(15,"mm"),
                 heatmap_legend_param = list(
                   legend_height=unit(2,"cm"), legend_direction="vertical",title_gp = gpar(fontsize = 10),
                   labels_gp =  gpar(fontsize = 8),title_position = "topleft"
                 ))


  return(zhc)
}
