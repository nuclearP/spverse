#' @include pre_treat.R
#'
NULL

#' the sp class
#'
#' class designed for one experiment
#'
#' @slot rawdata data.frame.
#' @slot features data.frame.
#' @slot scaledata data.frame.
#' @slot dim list.
#' @exportClass sp
sp <- setClass(
  "sp",
  slots = c(
    rawdata = "data.frame",
    features = "data.frame",
    scaledata = "data.frame",
    dim = "list"
  ),
  prototype = list(features = NULL,scaledata = NULL,dim = NULL)
)

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
  location <- paste0(pathway,"/","location.",filetype)
  if(!file.exists(rawdata)){
    stop("Original data not found")
  }
  rawdata <- fread(rawdata)
  row.names(rawdata) <- rawdata$p
  rawdata <- rawdata[-p]

  if(!file.exists(location)){
    object <- new(Class = "sp",rawdata = rawdata)
  }else{
    location <- fread(location)
    object <- new(Class = "sp",rawdata = rawdata,features = location)
  }
  object
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
  feature$means <- apply(object@rawdata,1,mean,na.rm = T)
  feature$cvs <- apply(object@rawdata,1,removena_cv)
  if(is.null(object@features)){
    object@features <- feature
  }else{
    object@features <- merge(object@features,feature,by = "samples")
  }
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
  selected <- subset(object@features,means < meancut,cvs < cvcut)
  scaledata <- object@rawdata[selected$samples]
  if(nrow(scaledata) < 1){
    stop("All data is filtered")
  }
  scaledata <- removeRowsAllNa(scaledata) %>% normalizeQuantiles(scaledata)
  object@scaledata <- scaledata
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


