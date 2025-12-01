#' @title cmeans_sp
#' @description
#' Fuzzy C-Means clustering for proteins. Perform
#' clustering on different features based on the results of sample clustering.
#' @param object an sp object \linkS4class{sp}.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Perform C-Means clustering using the slot clean_data.
#'   \item RID : Perform C-Means clustering using the slot RID_bootstrap.}
#' @param variable Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param ncluster number of clusters.
#' @returns Return the C-Means clustering results of features,
#' including a data frame feature_cluster containing category
#' information for different features, and an expression matrix of
#' different features in different variables.
#' @export
#' @rdname cmeans_sp
cmeans_sp <- function(object,name = "clean",variable,ncluster=10){

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
  }

  dat <- cmeans_data(object,name = name,variable=variable)
  dat <- new('ExpressionSet',exprs = dat)
  dat <- filter.std(dat,min.std=0,visu = F)
  dat <- standardise(dat)

  set.seed(1234)
  cl <- mfuzz(dat,c=ncluster,m=mestimate(dat))

  feature_cluster <- cbind(cl$cluster, cl$membership) %>% as.data.frame()
  colnames(feature_cluster)[1] <- 'catagory'
  feature_cluster[["features"]] = row.names(feature_cluster)

  zfeature <- dat@assayData$exprs %>% as.data.frame()
  zfeature[["features"]] <- row.names(zfeature)

  object@cmeans_list[[paste0("cmeans","_",name)]] <- list(feature_cluster,zfeature)

  names(object@cmeans_list[[paste0("cmeans","_",name)]]) <- c("feature_cluster","zfeature")

  return(object)
}

#' @title cmeans_elbowplot_sp
#' @description
#' Calculate the minimum centroid distances for a series of cluster
#' numbers and visualize the results to estimate the optimal number of clusters.
#' @param object an sp object \linkS4class{sp}.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Perform C-Means clustering using the slot clean_data.
#'   \item RID : Perform C-Means clustering using the slot RID_bootstrap.}
#' @param variable Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @export
#' @rdname cmeans_elbowplot_sp
cmeans_elbowplot_sp <- function(object,name = "clean",variable){

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
  }

  dat <- cmeans_data(object,name = name,variable=variable)
  dat <- new('ExpressionSet',exprs = dat)
  dat <- filter.std(dat,min.std=0,visu = F)
  dat <- standardise(dat)
  Dmin(dat,m=mestimate(dat))
}

#' @title cmeans_heatmap_sp
#' @description
#' Generate a heatmap to visualize the clustering outcomes of the C-Means clustering results.
#' @param object an sp object \linkS4class{sp}.
#' @param name \itemize{ \item clean :
#'   Visualize the C-Means clustering results from the slot clean_data.
#'   \item RID : Visualize the C-Means clustering results from the slot RID_bootstrap.}
#' @param min_mem The membership values range between 0 and 1 and are used for feature selection.
#' @param cluster_name Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param col Colors corresponding to value. A list contains named vectors.
#' @param show_row_names Whether show row names.
#' @param show_column_names Whether show column names.
#' @return A Heatmap-class object.
#' @export
#' @rdname cmeans_heatmap_sp
cmeans_heatmap_sp <- function(object,name = "clean",min_mem = 0.5,cluster_name,col=NULL,show_row_names=F,show_column_names=F){

  name <- match.arg(name, c("clean","RID"))
  if(name == "clean"){
    feature_cluster <- object@cmeans_list[["cmeans_clean"]][["feature_cluster"]]
    zfeature <- object@cmeans_list[["cmeans_clean"]][["zfeature"]]
    df_matrix <- object@clean_data
  }else{
    feature_cluster <- object@cmeans_list[["cmeans_RID"]][["feature_cluster"]]
    zfeature <- object@cmeans_list[["cmeans_RID"]][["zfeature"]]
    df_matrix <- object@RID_bootstrap
  }

  target_cols <- feature_cluster[, -c(1, ncol(feature_cluster))]
  filter_row <- apply(target_cols, 1, function(x) any(x > min_mem, na.rm = TRUE))
  df_filtered <- feature_cluster[filter_row, ]
  df_filtered <- df_filtered[, c(1,ncol(df_filtered))]

  zfeature_filtered <- merge(df_filtered,zfeature,all.x = TRUE)

  df_matrix <- t(df_matrix) %>% as.data.frame() %>% dplyr::select(all_of(zfeature_filtered[["features"]])) %>% t()

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

  zhc <- Heatmap(df_matrix,
                 name = "Intensity",
                 show_row_names = show_row_names,
                 show_column_names = show_column_names,
                 column_split = sample_features[[cluster_name]],
                 row_split = zfeature_filtered[["catagory"]],
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

#' @title cmeans_lineplot_sp
#' @description
#' Generate a heatmap to visualize the clustering outcomes of the C-Means clustering results.
#' @param object an sp object \linkS4class{sp}.
#' @param name \itemize{ \item clean :
#'   Visualize the C-Means clustering results from the slot clean_data.
#'   \item RID : Visualize the C-Means clustering results from the slot RID_bootstrap.}
#' @param min_mem The membership values range between 0 and 1 and are used for feature selection.
#' @param nrow Number of rows on the dividing surface.
#' @return ggplot2 object.
#' @export
#' @rdname cmeans_lineplot_sp
cmeans_lineplot_sp <- function(object,name = "clean",min_mem = 0.5,nrow=1){

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    feature_cluster <- object@cmeans_list[["cmeans_clean"]][["feature_cluster"]]
    zfeature <- object@cmeans_list[["cmeans_clean"]][["zfeature"]]
  }else{
    feature_cluster <- object@cmeans_list[["cmeans_RID"]][["feature_cluster"]]
    zfeature <- object@cmeans_list[["cmeans_RID"]][["zfeature"]]
  }

  target_cols <- feature_cluster[, -c(1, ncol(feature_cluster))]
  filter_row <- apply(target_cols, 1, function(x) any(x > min_mem, na.rm = TRUE))
  df_filtered <- feature_cluster[filter_row, ]
  df_filtered <- df_filtered[, c(1,ncol(df_filtered))]

  zfeature_filtered <- merge(df_filtered,zfeature,all.x = TRUE)

  df1 <- zfeature_filtered
  df1[["features"]] <- NULL
  df1 <- aggregate(. ~ catagory, data = df1, FUN = mean)
  df1 <- pivot_longer(df1,-catagory)
  df1[["group"]] = "1"
  df1[["name"]] <- factor(df1[["name"]],levels = unique(df1[["name"]]))

  df2 <- pivot_longer(zfeature_filtered,-c(features,catagory))
  df2[["name"]] <- factor(df2[["name"]],levels = unique(df2[["name"]]))

  p <- ggplot(df2,mapping=aes(x=name,y=value,group=features))+
    geom_line(color="#D2D2D2")+
    geom_line(df1,mapping=aes(x=name,y=value,group=group),size=0.5,color="#C00000")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 10))+
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(size = 1),
      legend.position = "none"
    )+
    facet_wrap(~catagory,nrow = nrow)+
    labs(x=NULL)

  return(p)
}


cmeans_data <- function(object,name = "clean",variable){
  if(name == "clean"){
    df <- object@clean_data
  }else{
    df <- object@RID_bootstrap
  }
  variable_df <- object@sample_features[,c("samples",variable)] %>% na.omit()
  if(ncol(df) != nrow(variable_df)){
    stop("The number of samples contained in the selected expression matrix is inconsistent with the number of non-NA values in the selected column of sample_features.")
  }
  df <- t(df) %>% as.data.frame()
  df[["samples"]] = row.names(df)
  df <- full_join(df,variable_df)
  df[["samples"]] <- NULL
  df <- df %>% group_by(.data[[variable]]) %>% summarise(across(where(is.numeric),mean,na.rm=T))
  row_name <- df[[variable]]
  df[[variable]] <- NULL
  row.names(df) <- row_name
  df <- t(df)
}
