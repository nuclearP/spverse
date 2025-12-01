#' @title cluster_sp
#' @description
#' Execute the ConsensusClusterPlus function for the given sp object.
#' @param object an sp object \linkS4class{sp}.
#' @param dimn Select the expression matrix of the first n dimensions from the
#' pca results for ConsensusClusterPlus function.
#' @param name There are three options in total:
#' \itemize{ \item clean :
#'   The pca results from the clean_data
#'   \item RID : The pca results of the data with individual differences removed through the bootstrap method (RID_bootstrap).
#'   \item harmony : The pca results of the data with individual differences removed through the harmony method (corrective_PCA).}
#' @param maxK integer value. maximum cluster number to evaluate.
#' @param seed optional numerical value. sets random seed for reproducible results.
#' @param reps integer value. number of subsamples.
#' @param pItem numerical value. proportion of items to sample.
#' @param pFeature numerical value. proportion of features to sample.
#' @param clusterAlg character value. cluster algorithm. 'hc' hierarchical (hclust), 'pam' for paritioning
#'  around medoids, 'km' for k-means upon data matrix, or a function that returns a clustering. See example
#'  and vignette for more details.
#' @return A folder containing clustering results.
#' @export
#' @rdname cluster_sp
#' @seealso \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
cluster_sp <- function(object,name = "clean",dimn,maxK = 10,seed = 1234,reps = 5000,pItem = 0.8,pFeature = 1,clusterAlg = 'km') {

  name <- match.arg(name, c("clean","RID","harmony"))

  if(name == "clean"){
    if(is.null(object@dim_reductions[["pca_clean"]])){
      stop("The pca_clean is missing in the sp object.")
    }
    x <- object@dim_reductions[["pca_clean"]]
    df1 <- x[["x"]]
    fi <- df1[,1:dimn]
    fi <- t(fi)
  }else if(name == "RID"){
    if(is.null(object@dim_reductions[["pca_RID"]])){
      stop("The pca_RID slot is missing in the sp object.")
    }
    x <- object@dim_reductions[["pca_RID"]]
    df1 <- x[["x"]]
    fi <- df1[,1:dimn]
    fi <- t(fi)
  }else{
    if(is.null(object@corrective_PCA)){
      stop("The corrective_PCA slot is missing in the sp object.")
    }
    fi <- as.matrix(object@corrective_PCA)
  }

  cluster <- ConsensusClusterPlus(
    d = fi, #数据矩阵，列是样本，行是变量
    maxK = maxK, #要评估的最大聚类簇数量
    seed = seed, #指定随机数种子，用于子样本处理等随机的过程
    reps = reps, #抽取的子样本数量
    pItem = pItem, #抽样样本的比例
    pFeature = pFeature, #抽样变量的比例
    clusterAlg = clusterAlg, #选择聚类算法
    title = name,
    plot = 'pdf', #聚类簇评估结果的输出格式
    writeTable=T
  )
}

#' @title cluster_read_sp
#' @description
#' Select the results from the clustering files and add a new column to the sample feature data frame of the sp object.
#' @param object \linkS4class{sp}.
#' @param name keep consistency with the setting of the name parameter in the cluster_sp function.
#' @param n_cluster integer value. number of clusters.
#' @param col_name the name of the newly generated column.
#' @export
#' @rdname cluster_read_sp

cluster_read_sp <- function(object,name = "clean",n_cluster,col_name) {

  file <- paste0(getwd(),"/",name,"/",name,".k=",n_cluster,".consensusClass.csv")
  xx <- read.csv(file,header = F)
  xx[["V2"]] <- paste0("cluster",xx[["V2"]])
  k <- xx[["V2"]]

  colnames(xx) <- c("samples",col_name)
  object@sample_features <- full_join(object@sample_features,xx)
  return(object)
}
