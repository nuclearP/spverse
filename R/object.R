#' the sp class
#'
#' class designed for one experiment.
#'
#' @slot rawdata data.frame.
#' @slot sample_features data.frame.
#' @slot protein_features data.frame.
#' @slot clean_data data.frame.
#' @slot dim_reductions list.
#' @slot corrective_PCA data.frame.
#' @slot RID_bootstrap data.frame.
#' @slot cmeans_list "list".
#' @slot ssgsea_list list.
#' @slot GSEA_list list.
#' @slot adj_matrix_list list.
#' @slot list_neighbour list.
#' @slot combine_listw ANY.
#' @slot deconvolution_result data.frame.
#' @export

sp <- setClass(
  "sp",
  slots = list(
    rawdata = "data.frame",
    sample_features = "data.frame",
    protein_features = "data.frame",
    clean_data = "data.frame",
    dim_reductions = "list",
    corrective_PCA = "data.frame",
    RID_bootstrap  = "data.frame",
    cmeans_list = "list",
    ssgsea_list = "list",
    GSEA_list = "list",
    adj_matrix_list = "list",
    list_neighbour = "list",
    combine_listw = "ANY",
    deconvolution_result = "data.frame"
  ),
  validity = function(object) {
    if (is.null(object@combine_listw)) {
      return(TRUE)
    }
    if (!inherits(object@combine_listw, "listw")) {
      return("combine_listw must be of the 'listw' type (S3) from the spdep package or NULL")
    }
    TRUE
  }
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
  feature <- paste0(pathway,"/","feature.",filetype)
  location <- paste0(pathway,"/","location.xlsx")

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

  feature <- fread(feature) %>% as.data.frame()
  sample_features <- merge(feature,sample_features,by = "samples")

  #creat_object
  if(!file.exists(location)){
    object <- new(Class = "sp",rawdata = rawdata,sample_features = sample_features,protein_features = protein_features)
  }else{
    location <- paste0(pathway,"/","location.xlsx")
    location <- creat_coordinate(location)
    sample_features <- merge(location,sample_features,by = "samples")
    order_idx <- match(colnames(rawdata),sample_features$samples)
    sample_features <- sample_features[order_idx, ]
    object <- new(Class = "sp",rawdata = rawdata,sample_features = sample_features,protein_features = protein_features)
  }
  return(object)
}

#' @importFrom readxl excel_sheets
#' @importFrom readxl read_excel
#' @importFrom stats na.omit
#' @importFrom data.table rbindlist
creat_coordinate <- function(pathway){

  sheet_names <- excel_sheets(pathway)
  all_sheets <- lapply(sheet_names, read_excel, path = pathway,col_names = F)
  names(all_sheets) <- sheet_names

  fi <- list()
  for (i in 1:length(sheet_names)) {
    matrix_sheet <- as.matrix(all_sheets[[i]])
    valid_cells <- !is.na(matrix_sheet)
    ####location
    valid_indices <- which(valid_cells, arr.ind = TRUE) %>% as.data.frame()
    valid_indices[["row"]] <- (max(valid_indices$row)+1) - valid_indices$row
    valid_indices["samples"] <- na.omit(as.vector(matrix_sheet))
    # valid_indices[["individual"]] <- sheet_names[i]
    names(valid_indices) <- c("y","x","samples")
    fi[[i]] <- valid_indices
  }
  fi <- rbindlist(fi)
}


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

