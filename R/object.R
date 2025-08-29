#' @include pre_treat.R
#'
NULL

#' the sp class
#'
#' class designed for one experiment
#'
#' @slot rawdata data.frame.
#' @slot sample_features data.frame.
#' @slot protein_features data.frame.
#' @slot scaledata data.frame.
#' @slot dims list.
#' @slot corrective_PCA data.frame.
#' @slot trend list.
#' @slot group factor.
#'
#' @return
#' @export
#'
#' @examples
sp <- setClass(
  "sp",
  slots = c(
    rawdata = "data.frame",
    sample_features = "data.frame",
    protein_features = "data.frame",
    scaledata = "data.frame",
    dims = "list",
    corrective_PCA = "data.frame",
    trend = "list",
    group = "factor"
  )
)

#' @method [ sp
#' @export
`[.sp` <- function(x, i, j, ...) {
  y <- x@feature[i, j, ...]
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
  location_ <- paste0(pathway,"/","location.",filetype)
  if(!file.exists(rawdata)){
    stop("Original data not found")
  }
  #rawdata
  rawdata <- fread(rawdata) %>% as.data.frame()
  rownames(rawdata) <- rawdata$p
  rawdata <- subset(rawdata,select = -rawdata$p)
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
  protein_features$means <- apply(rawdata,2,mean,na.rm = T)
  protein_features$cvs <- apply(rawdata,2,removena_cv)
  #scaledata
  int.scaledata <- data.frame()
  #creat_object
  if(!file.exists(location)){
    object <- new(Class = "sp",rawdata = rawdata,sample_features = sample_features,protein_features = protein_features,sclaedata = int.scaledata,dim = list())
  }else{
    location <- fread(location)
    sample_features <- merge(location,sample_features,by = "samples")
    object <- new(Class = "sp",rawdata = rawdata,sample_features = sample_features,protein_features = protein_features,sclaedata = int.scaledata,dim = list())
  }
  return(object)
}

#' @export
setGeneric("mean_cv",function(object) standardGeneric("mean_cv") )

#' @title mean_cv
#' @description Calculate mean and cv
#' @param object  A parameter object of one of the following classe sp or spatials
#' @return a dataframe containing the mean and standard deviation of each sample
#' @name mean_cv
#' @rdname mean_cv
NULL

#' @rdname mean_cv
#' @exportMethod mean_cv
setMethod("mean_cv","sp",function(object){
  feature <- colSums(!is.na(object@rawdata)) %>% as.data.frame()
  colnames(feature) <- "num"
  feature$samples <- row.names(feature)
  feature$means <- apply(object@rawdata,2,mean,na.rm = T)
  feature$cvs <- apply(object@rawdata,2,removena_cv)
  if(is.null(object@features)){
    object@features <- feature
  }else{
    object@features <- merge(object@features,feature,by = "samples")
  }
  return(object)
}
)


#' @export
setGeneric("normlization",function(object,...) standardGeneric("normlization") )

#' @title normlization
#' @description normalize raw data
#' @importFrom limma normalizeQuantiles
#' @param object  A parameter object of one of the following classe sp or spatials
#' @param meancut  A parameter object of one of the following classe sp or spatials
#' @param cvcut  A parameter object of one of the following classe sp or spatials
#' @return A dataframe contains data that has been quality controlled and normalized
#' @name normlization
#' @rdname normlization
NULL

#' @rdname normlization
#' @exportMethod normlization
setMethod("normlization","sp",function(object,meancut,cvcut){
  if(!is.numeric(meancut)){
    stop("meancut should be a numerical value")
  }
  if(!is.numeric(cvcut)){
    stop("cvcut should be a numerical value")
  }
  selected <- subset(object@features,object@features$means < meancut & object@features$cvs < cvcut)
  scaledata <- subset(object@rawdata,select = selected$samples)
  if(nrow(scaledata) < 1){
    stop("All data is filtered")
  }
  scaledata <- removeRowsAllNa(scaledata) %>% normalizeQuantiles()
  object@scaledata <- scaledata
  return(object)
}
)

#' the spatials class
#'
#' An S4 object that stores multiple sp objects, but only one sp object is activated at the same time
#'
#' @slot experiment list.
#' @slot active character.
#' @exportClass spatials
spatials <- setClass(
  "spatials",
  slots = c(
    experiment = "list",
    active = "character"
  )
)

#' @name show
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

setMethod(
  f = 'dim',
  signature = 'sp',
  definition = function(object){
    dim(object@rawdata)
  }
)


