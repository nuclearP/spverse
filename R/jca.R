#' @title jca_sp
#' @description
#' Given the global counts in each cluster,
#' expected counts and variances are calculated under non-free sampling, and a z-value reported.
#' @param object an sp object \linkS4class{sp}.
#' @param cluster cluster information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object.
#' @param union default FALSE, if TRUE all individuals' neighbours list will be consolidated into one large list.
#' @return if union=F, A matrix with class jcmulti with row and column names for observed and expected
#' counts, variance, and z-value. if union=T, Return a list containing the jcmulti matrix for each individual.
#'  expected counts, variance, and z-value.
#' @export
#' @rdname jca_sp
jca_sp <- function(object,cluster,union = F){

  if(is.null(object@adj_matrix_list)){
    stop("Spatial adjacency matrix is missing.")
  }

  individuals <- names(object@adj_matrix_list)
  df <- object@sample_features
  df <- df[,c("individual",cluster)] %>% na.omit()

  if(union == T){
    k <- factor(df[[cluster]])
    jca_df <- joincount.multi(k,object@combine_listw)
    return(jca_df)
  }else{
    jca_list <- list()
    for (i in individuals) {
      df_new <- df[df$`individual` == i, ]
      k <- factor(df_new[[cluster]])
      jca_df <- joincount.multi(k,object@list_neighbour[[i]])
      jca_list[[i]] <- jca_df
    }
    return(jca_list)
  }
}
