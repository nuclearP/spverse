#' @title deconvolution_sp
#' @description
#' Perform deconvolution for the sp object using the CIBERSORT method.
#' @param object an sp object \linkS4class{sp}.
#' @param slot There are two options in total:
#' \itemize{ \item clean_data :
#'   Perform cell deconvolution using the slot clean_data.
#'   \item RID_bootstrap : Perform cell deconvolution using the slot RID_bootstrap.}
#' @param sig_matrix A matrix of expression profile of cells.
#' @param fromType The ID type of the row names in the selected expression matrix is default to NULL,
#' and it is used when converting the row names of the selected expression matrix to other ID type.
#' @param toType Corresponds to the fromType parameter, specifying the converted ID type.
#' @param organ Filter features with high correlation between transcriptome and proteome
#' in the mix_matrix based on the organ types of the samples in the matrix, which will be
#' used for deconvolution. If your sig_matrix is derived from RNA,
#' feature filtering is required prior to deconvolution.
#' @details
#' \itemize{ \item
#' The object and slot parameters are used to select the mix_matrix,
#' while the fromType and toType parameters are employed to convert
#' the IDs of the matrix for compatibility with the IDs of sig_matrix.
#'   \item When the organ parameter is not NULL, fromType cannot be NULL.}
#' @export
#' @rdname deconvolution_sp
#' @return Assign a value to the slot deconvolution_result of the sp object.
deconvolution_sp <- function(object,slot="clean_data",sig_matrix,fromType=NULL,toType=NULL,organ = NULL){

  slot <- match.arg(slot, c("clean_data","RID_bootstrap"))

  if(slot == "clean_data"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
    mix_matrix <- object@clean_data
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
    mix_matrix <- object@RID_bootstrap
  }

  if(!is.null(organ)){
    if(!(organ %in% c("Bladder","Bone","Breast",
                      "Central Nervous System","Cervix","Endometrium",
                      "Esophagus","Haematopoietic and Lymphoid","Head and Neck",
                      "Kidney","Large Intestine","Liver",
                      "Lung","Ovary","Pancreas" ,
                      "Peripheral Nervous System","Prostate","Skin",
                      "Soft tissue","Stomach","Thyroid"))){
      stop("The organ type you entered is invalid. Please refer to the function usage instructions.")
    }
    if(is.null(fromType)){
      stop("When the organ parameter is not NULL, fromType cannot be NULL.")
    }
    rds_path <- system.file("extdata", "HighCor.xlsx", package = "spverse")
    higr_cor <- readxl::read_xlsx(rds_path)
    target_features <- higr_cor[[organ]] %>% na.omit()
    if(fromType == "SYMBOL"){
      mix_matrix[["features"]] <- row.names(mix_matrix)
      mix_matrix[mix_matrix[["features"]] %in% target_features, ]
      mix_matrix[["features"]] <- NULL
    }else{
      ids = row.names(mix_matrix)
      id_conversion <- bitr(
        geneID = ids,
        fromType = fromType,
        toType = "SYMBOL",
        OrgDb = species)
      id_conversion[id_conversion[["SYMBOL"]] %in% target_features, ]
      slected_features <- id_conversion[,1] %>% na.omit()
      mix_matrix[["features"]] <- row.names(mix_matrix)
      mix_matrix[mix_matrix[["features"]] %in% target_features, ]
      mix_matrix[["features"]] <- NULL
    }
  }

  if((!is.null(fromType)) & (!is.null(toType))){
    mix_matrix[[fromType]] <- row.names(mix_matrix)
    ids <- row.names(mix_matrix)
    id_conversion <- bitr(
      geneID = ids,
      fromType = fromType,
      toType = toType,
      OrgDb = "org.Hs.eg.db")
    mix_matrix <- full_join(mix_matrix,id_conversion)
    mix_matrix[[fromType]] <- NULL
    mix_matrix <- mix_matrix %>% na.omit() %>% group_by(.data[[toType]]) %>%
      summarise(mean_value = mean(value))
    to_id <- mix_matrix[[toType]]
    mix_matrix[[toType]] <- NULL
    mix_matrix <- as.data.frame(mix_matrix)
    row.names(mix_matrix) <- to_id
  }

  results <- cibersort(sig_matrix = sig_matrix, mixture_file = mix_matrix,perm = 100)
  object@deconvolution_result <- results
  return(object)
}

#' @title deconvolution_view
#' @description
#' Create a ggplot2-based pie chart for the deconvolution result.
#' @param object an sp object \linkS4class{sp}.
#' @param nrow Number of rows on the dividing surface.
#' @param pal The color palette to be used for coloring and filling by groups.
#' @param legend_row The number of rows in the legend. If there are too many types of groups,
#' you can appropriately increase the number of rows to improve the aesthetics of the plot.
#' @return ggplot object.
#' @export
#' @rdname deconvolution_view
deconvolution_view <- function(object,nrow=1,pal=NULL,legend_row=1){

  df_dec <- object@deconvolution_result
  df_dec <- df_dec[,1:(ncol(df_dec)-3)]
  df_dec[["samples"]] <- row.names(df_dec)
  if(!is.null(pal)){
    if(length(pal) != ncol(df_dec)){
      stop("The length of the color vector you entered is invalid. Please input a color vector of length ",ncol(df_dec),".")
    }
  }

  locs <- object@sample_features
  locs <- locs[,c("x","y","samples","individuals")] %>% ma.omit()

  df <- merge(locs,df_dec,by="samples")
  df_long <- pivot_longer(df,-c("x","y","samples","individuals"),names_to = c("type"),values_to = "per")

  if(is.null(pal)){
    p <- ggplot()+
      geom_tile(data =final,aes(x=x,y=y),fill = "white",color = "black",size=0.8)+
      geom_jjPointPie(data = df_long,aes(x = x,y = y,pievar = per,fill = type,group = samples),width = 0.8,color = 'white',line.size = .05)+
      facet_wrap(~individuals,nrow = nrow)+
      theme(panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank()
      )+
      guides(fill=guide_legend(nrow = legend_row))+
      coord_fixed()
  }else{
    p <- ggplot()+
      geom_tile(data =final,aes(x=x,y=y),fill = "white",color = "black",size=0.8)+
      geom_jjPointPie(data = df_long,aes(x = x,y = y,pievar = per,fill = type,group = samples),width = 0.8,color = 'white',line.size = .05)+
      facet_wrap(~individuals,nrow = nrow)+
      scale_fill_manual(values = pal)
      theme(panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank()
      )+
      guides(fill=guide_legend(nrow = legend_row))+
      coord_fixed()
  }

  return(p)
}
