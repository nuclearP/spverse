#' @title jaccard_sp
#' @description
#' Calculate the Jaccard index between a spot and its neighbors based on spatial
#'  adjacency relationships, and assign the results to a column in the sample_features slot.
#' @param object an sp object \linkS4class{sp}.
#' @export
#' @rdname jaccard_sp
jaccard_sp <- function(object){

  if(is.null(object@list_neighbour)){
    stop("The slot list_neighbour is missing.")
  }

  individual <- names(object@list_neighbour)

  fi <- list()

  for (i in individual) {
    fi[[i]] <- object@list_neighbour[[i]][["neighbours"]]
  }

  combined_list <- do.call(c, fi)

  df <- object@rawdata
  df_new <- df %>% dplyr::select(colnames(object@clean_data))

  pros <- list()

  for (i in (1:length(combined_list))) {
    pro <- df_new %>% dplyr::select(all_of(i)) %>% na.omit()
    pros[[i]] <- row.names(pro)
  }

  jaccards <- numeric(length(combined_list))

  for (i in (1:length(combined_list))) {
    p_1 <- pros[[i]]
    p_others <- combined_list[[i]]
    if(any(p_others == 0)){
      jaccards[[i]] <- NA
    }else{
      jaccard_one_sample <- numeric(length(p_others))
      for (j in (1:length(p_others))) {
        jaccard_one_sample[[j]] <- jaccard(p_1,pros[[p_others[[j]]]])
      }
      jaccards[[i]] <- mean(jaccard_one_sample,na.rm=T)
    }

  }

  fi <- data.frame(samples = colnames(df_new),jaccard_index = jaccards)

  object@sample_features <- full_join(object@sample_features,fi)

  return(object)
}

#' @title correlation_sp
#' @description
#' Calculate the correlation index between a spot and its neighbors based on spatial
#'  adjacency relationships, and assign the results to a column in the sample_features slot.
#' @param object an sp object \linkS4class{sp}.
#' @param use \code{\link[stats]{cor}}
#' @param method \code{\link[stats]{cor}}
#' @export
#' @rdname correlation_sp
correlation_sp <- function(object,use="pairwise.complete.obs",method="pearson"){

  if(is.null(object@list_neighbour)){
    stop("The slot list_neighbour is missing.")
  }

  individual <- names(object@list_neighbour)

  fi <- list()

  for (i in individual) {
    fi[[i]] <- object@list_neighbour[[i]][["neighbours"]]
  }

  combined_list <- do.call(c, fi)

  df <- object@rawdata

  df_new <- df %>% dplyr::select(colnames(object@clean_data))

  cor_matrix <- cor(df_new,use=use,method=method)

  corrs <- numeric(length(combined_list))

  for (i in (1:length(combined_list))) {
    cor_others <- combined_list[[i]]
    if(any(cor_others == 0)){
      corrs[[i]] <- NA
    }else{
      cor_one_sample <- numeric(length(cor_others))
      for (j in (1:length(cor_others))) {
        cor_one_sample[[j]] <- cor_matrix[j,i]
      }
      corrs[[i]] <- mean(cor_one_sample,na.rm=T)
    }
  }

  fi <- data.frame(samples = colnames(df_new),correlation = corrs)
  object@sample_features <- full_join(object@sample_features,fi)

  return(object)
}


jaccard <- function(a, b) {
  return(length(intersect(a, b)) / length(union(a, b)))
}
