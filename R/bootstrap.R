#' @include utilities.R
NULL

#' @title bootstrap
#' @description
#' Use the bootstrap method to calculate the confidence intervals of each protein
#' in different groups, and filter out proteins with larger thresholds to control
#' the impact of individual differences on data analysis.
#' @param x A protein expression matrix, where the row names are protein
#' or gene IDs and the column names are sample names.
#' @param ... Arguments passed to other methods.
#' @return A protein expression matrix after removing some proteins with large fluctuations..
#' @export
#' @rdname bootstrap
bootstrap <- function(x,...) {
  UseMethod("bootstrap")
}

#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @param group The grouping information of samples, where each group contains samples from different individuals.
#' @param dif The range of the confidence interval.
#' @param R The number of bootstrap replicates.
#' @details
#' Since proteomics data is generally log-transformed, the range of the confidence interval—i.e.,
#' the difference between the two ends of the confidence interval—can be interpreted as the fold
#' change between these two ends.
#' @method bootstrap data.frame
#' @export
#' @rdname bootstrap
bootstrap.data.frame <- function(x,group,dif = 1.5,R=500,...){

  fis <- list()

  df <- data.frame(samples = colnames(x),group = group)

  for (i in unique(group)) {
    sams <- dplyr::filter(df,group == i)
    sub_samples <- dplyr::select(x,sams$samples) %>% t() %>% as.data.frame()

    sds <- apply(sub_samples,2,removena_sd) %>% as.data.frame()
    colnames(sds) <- "sd"
    sd0 <- dplyr::filter(sds,sd == 0)
    sd_no0 <- dplyr::filter(sds,sd != 0)

    sub_samples <- dplyr::select(sub_samples,row.names(sd_no0))
    fi <- boot_df(sub_samples,R=R)
    fi <- dplyr::filter(fi,difference < log2(dif))

    fi_proteins <-  fi[["Protein"]]
    fis[[i]] <- c(fi_proteins,row.names(sd0))
  }

  common_proteins <- Reduce(f = intersect, x = fis)

  x_fi <- x %>%  mutate(p=row.names(x)) %>% dplyr::filter(p %in% common_proteins) %>% dplyr::select(-p)

  return(x_fi)
}

#' @title bootstrap_sp
#' @description
#' Use the bootstrap method to calculate the confidence intervals of each protein
#' in different groups, and filter out proteins with larger thresholds of sp object to control
#' the impact of individual differences on data analysis.
#' @param object An \linkS4class{sp}.
#' @param ... Arguments passed to other methods.
#' @return ggplot object.
#' @export
#' @rdname bootstrap_sp
setGeneric("bootstrap_sp",function(object,...) standardGeneric("bootstrap_sp") )

#' @param group Group information of samples, which is a string corresponding.
#' to one element in the column name "sample_features" of the sp object
#' The attribute of this column should be a character vector or a factor.
#' @param dif The range of the confidence interval.
#' @param R The number of bootstrap replicates.
#' @seealso \code{\link{bootstrap.data.frame}}
#' @rdname bootstrap_sp
#' @exportMethod bootstrap_sp
setMethod(
  f ="bootstrap_sp",
  signature = signature(object= "sp"),
  function(object,group,dif = 1.5,R=500) {

    if(nrow(object@cleaned_data) < 1){
      stop("The cleaned_data slot is missing in the sp object.")
    }

    sam = object@sample_features[c("samples",group)]

    sam <- dplyr::filter(sam,sam[["samples"]] %in% colnames(object@cleaned_data))

    object@RID_bootstrap <- bootstrap(object@cleaned_data,group = sam[[group]],dif = dif,R=R)

    return(object)

  }
)


mean_func <- function(data, indices) {
  sample_data <- data[indices]
  return(mean(sample_data, na.rm = TRUE)) # 处理可能的缺失值
}

# 使用apply函数处理所有列
bootstrap_ci <- function(column,R=500) {
  print(R)
  boot_result <- boot(column, statistic = mean_func, R = R)
  boot_ci <- boot.ci(boot_result, conf = 0.95, type = "perc")
  return(c(
    Lower_CI = boot_ci$percent[4],
    Upper_CI = boot_ci$percent[5],
    Mean = mean(column, na.rm = TRUE)
  ))
}

boot_df <- function(x,R=500){

  # 应用函数到所有列
  results_matrix <- apply(x, 2, bootstrap_ci,R=R)

  # 转换为数据框
  results_df <- as.data.frame(t(results_matrix))
  results_df$Protein <- rownames(results_df)
  rownames(results_df) <- NULL

  # 重新排列列
  results_df <- results_df[, c("Protein", "Lower_CI", "Upper_CI", "Mean")]
  results_df <- mutate(results_df,difference = Upper_CI - Lower_CI)
}
