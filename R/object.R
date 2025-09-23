#' the sp class
#'
#' class designed for one experiment.
#'
#' @slot rawdata data.frame.
#' @slot sample_features data.frame.
#' @slot protein_features data.frame.
#' @slot cleaned_data data.frame.
#' @slot dim_reductions list.
#' @slot corrective_PCA data.frame.
#' @slot RID_bootstrap data.frame.
#' @export

sp <- setClass(
  "sp",
  slots = list(
    rawdata = "data.frame",
    sample_features = "data.frame",
    protein_features = "data.frame",
    cleaned_data = "data.frame",
    dim_reductions = "list",
    corrective_PCA = "data.frame",
    RID_bootstrap  = "data.frame"
  )
)

#' @method [ sp
#' @export
`[.sp` <- function(x, i, j, ...) {
  y <- x@sample_feature[i, j, ...]
  return(y)
}

##' @method $ sp
##' @export
`$.sp` <-  function(x, name) {
  x@rawdata[, name]
}

#' creatsp
#'
#' creat sp class
#'
#' @param pathway The path where the original data is stored
#' @param filetype File extension
#'
#' @return a sp class
#' @export
creatsp <- function(
  pathway,
  filetype = "xls"
){
  rawdata <- paste0(pathway,"/","rawdata.",filetype)
  location <- paste0(pathway,"/","sample_feature.",filetype)
  if(!file.exists(rawdata)){
    stop("Original data not found")
  }
  #rawdata
  rawdata <- fread(rawdata) %>% as.data.frame()
  rownames(rawdata) <- rawdata$p
  rawdata <- subset(rawdata,select = -p)
  #sample_features
  sample_features <- colSums(!is.na(rawdata)) %>% as.data.frame()
  colnames(sample_features) <- "num"
  sample_features$samples <- row.names(sample_features)
  sample_features$means <- apply(rawdata,2,mean,na.rm = T)
  sample_features$cvs <- apply(rawdata,2,removena_cv)
  #protein_feature
  protein_features <- rowSums(!is.na(rawdata)) %>% as.data.frame()
  colnames(protein_features) <- "num"
  protein_features$samples <- row.names(protein_features)
  protein_features$means <- apply(rawdata,1,mean,na.rm = T)
  protein_features$cvs <- apply(rawdata,1,removena_cv)

  #creat_object
  if(!file.exists(location)){
    object <- new(Class = "sp",rawdata = rawdata,sample_features = sample_features,protein_features = protein_features)
  }else{
    location <- fread(location)
    sample_features <- merge(location,sample_features,by = "samples")
    object <- new(Class = "sp",rawdata = rawdata,sample_features = sample_features,protein_features = protein_features)
  }
  return(object)
}

##' @importFrom methods show
#' @exportMethod show
setMethod(
  f = 'show',
  signature = 'sp',
  definition = function(object) {
    cat(
      class(object)[1],
      'data with',
      ncol(object@rawdata),
      'samples for',
      nrow(object@rawdata), 'features\n'
    )
}
)

#' @exportMethod dim
setMethod(
  f = 'dim',
  signature = 'sp',
  definition = function(x){
    dim(x@rawdata)
  }
)


