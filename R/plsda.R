#' @include utilities.R
NULL

#' @title run_plsda
#' @description
#' Performs plsda analysis on the given protein matrix and
#' returns the results as an object of class prcomp.
#' @param x A protein expression matrix (data.frame format), where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param filter_by Select features with large fluctuations for plsda analysis,
#' The standard deviation (sd) or coefficient of variation (cv) can be used for feature selection.
#' @param rank A numeric value between 0 and 1. Select features with the largest fluctuations based on percentage.
#' By default, select the top 50\% of features for principal components analysis.
#' @param group Group information of samples, which is a character vector or factor.
#' @details
#' If your protein expression matrix is in logarithmic form (e.g., after log2 transformation),
#' we recommend using the standard deviation as the screening criterion.
#' This is because the difference between data represents a relative change in fold (multiples),
#' rather than the absolute value of abundance.
#' @return returns a list with class "plsda".
#' @rdname run_plsda
run_plsda <- function(x,filter_by = "sd",rank=0.5,group){

  filter_by <- match.arg(filter_by, c("sd", "cv"))

  if(filter_by == "sd"){

    filters = apply(x,1,removena_sd)

    x_fi <- mutate(x,sds = filters) %>%
      slice_max(sds,n=round((nrow(x)*rank))) %>%
      dplyr::select(-sds)
  }
  if(filter_by == "cv"){

    filters = apply(x,1,removena_cv)

    x_fi <- mutate(x,cvs = filters) %>%
      slice_max(cvs,n=round((nrow(x)*rank))) %>%
      dplyr::select(-cvs)
  }

  nrow_x_fi <- nrow(x_fi)

  message(nrow_x_fi," variables (proteins) will be used for plsda analysis.")

  x_fi <- t(x_fi)

  set.seed(1)
  plsda_fit <- plsda(x_fi,group,ncomp=2)

  return(plsda_fit)
}

#' @title run_plsda_sp
#' @description
#' Performs plsda analysis on the sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param filter_by Select features with large fluctuations for principal components analysis,
#' The standard deviation (sd) or coefficient of variation (cv) can be used for feature selection.
#' @param rank Select features with the largest fluctuations based on percentage.
#' By default, select the top 50\% of features for principal components analysis.
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Perform plsda dimensionality reduction using the slot clean_data, and store the returned result in dim_reductions[["plsda_clean"]].
#'   \item RID : Perform plsda dimensionality reduction using the slot RID_bootstrap, and store the returned result in dim_reductions[["pca_RID"]].}
#' @seealso \code{\link{run_plsda}}
#' @export
#' @rdname run_plsda_sp
run_plsda_sp <- function(object,filter_by = "sd",rank=0.5,group,name="clean") {

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
    sam = object@sample_features[c("samples",group)]
    sam <- dplyr::filter(sam,sam[["samples"]] %in% colnames(object@clean_data))
    object@dim_reductions[["plsda_clean"]] <- run_plsda(object@clean_data,filter_by = filter_by,rank=rank,group=sam[[group]])
    return(object)
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
    sam = object@sample_features[c("samples",group)]
    sam <- dplyr::filter(sam,sam[["samples"]] %in% colnames(x@RID_bootstrap))
    object@dim_reductions[["plsda_RID"]] <- run_plsda(object@RID_bootstrap,filter_by = filter_by,rank=rank,group=sam[[group]])
  }
  return(object)
}

#' plsda_point_sp
#' Create a ggplot2-based scatter plot of Plsda for sp object.
#' @param object \linkS4class{sp}.
#' @param col The color palette to be used for coloring and filling by groups.
#' @param name There are two options in total:
#' \itemize{ \item clean :
#'   Visualize the plsda results of the original data.
#'   \item RID : Visualize the plsda results of the data with individual differences removed.}
#' @return ggplot object.
#' @export
#' @rdname plsda_point_sp
plsda_point_sp <- function(object,name,col=NULL){

  name <- match.arg(name, c("clean","RID"))

  if(name == "clean"){
    if(is.null(object@dim_reductions[["plsda_clean"]])){
      stop("The plsda_clean is missing in the sp object.")
    }
    df_plsda <- object@dim_reductions[["plsda_clean"]] %>% unclass()
    df = as.data.frame(df_plsda$variates$X)
    df[["group"]] = df_plsda[["Y"]]
    df[["samples"]] = rownames(df)

    explain = df_plsda$prop_expl_var$X
    x_lable <- round(explain[1],digits=3)
    y_lable <- round(explain[2],digits=3)
  }else{
    if(is.null(object@dim_reductions[["plsda_RID"]])){
      stop("The plsda_RID slot is missing in the sp object.")
    }
    df_plsda <- object@dim_reductions[["plsda_RID"]] %>% unclass()
    df = as.data.frame(df_plsda$variates$X)
    df[["group"]] = df_plsda[["Y"]]
    df[["samples"]] = rownames(df)

    explain = df_plsda$prop_expl_var$X
    x_lable <- round(explain[1],digits=3)
    y_lable <- round(explain[2],digits=3)
  }

  if(is.null(col)){
    p1<-ggplot(df,aes(x=comp1,y=comp2,
                      color=group,shape=group))+
      theme_bw()+
      geom_point()+
      theme(panel.grid = element_blank())+
      geom_vline(xintercept = 0,lty="dashed")+
      geom_hline(yintercept = 0,lty="dashed")+
      labs(x=paste0("P1 (",x_lable*100,"%)"),
           y=paste0("P2 (",y_lable*100,"%)"))+
      stat_ellipse(data=df,
                   geom = "polygon",level = 0.95,
                   linetype = 2,size=0.5,
                   aes(fill=group),
                   alpha=0.2,
                   show.legend = T)
  }else{
    p1<-ggplot(df,aes(x=comp1,y=comp2,
                      color=group,shape=group))+
      theme_bw()+
      geom_point()+
      theme(panel.grid = element_blank())+
      geom_vline(xintercept = 0,lty="dashed")+
      geom_hline(yintercept = 0,lty="dashed")+
      labs(x=paste0("P1 (",x_lable*100,"%)"),
           y=paste0("P2 (",y_lable*100,"%)"))+
      stat_ellipse(data=df,
                   geom = "polygon",level = 0.95,
                   linetype = 2,size=0.5,
                   aes(fill=group),
                   alpha=0.2,
                   show.legend = T)+
      scale_color_manual(values = col)+
      scale_fill_manual(values = col)
  }
  return(p1)
}
