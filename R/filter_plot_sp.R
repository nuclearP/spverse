#' @title filter_plot_sp
#' @description
#' Draw violin plots of protein identification numbers and filter dashed lines.
#' @param object an sp object \linkS4class{sp}.
#' @param times The coefficient of the IQR value, with a default value of 1.
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param nrow Number of rows on the dividing surface.
#' @return ggplot object.
#' @export
#' @rdname filter_plot_sp
filter_plot_sp <- function(object,times = 1,group=NULL,nrow = 1) {

  if(is.null(group)){
    df <- colSums(!is.na(object@rawdata)) %>%
      as.data.frame() %>%
      magrittr::set_names("values") %>%
      dplyr::mutate(sample = colnames(object@rawdata),
                    hline = (quantile(values, 0.25, na.rm = TRUE) -
                               times * IQR(values, na.rm = TRUE)),
                    name = "X")

    ggplot(df,aes(x= name, y=values,color= name))+
      geom_violin(aes(fill=name),color="white",alpha=0.6)+
      geom_boxplot(width=0.05,position = position_dodge(0.9),color="white",outliers = F)+
      geom_jitter(size = 0.5,alpha = 1,width = 0.2)+
      geom_hline(aes(yintercept = hline), linetype = "dashed", color = "gray",size = 0.8)+
      labs(x = NULL,y = "Num of Proteins")+
      theme_classic()+
      theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank())

  }else{

    df <- colSums(!is.na(object@rawdata)) %>%
      as.data.frame() %>%
      magrittr::set_names("values") %>%
      dplyr::mutate(sample = colnames(object@rawdata),group = object@sample_features[[group]]) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(hlines = (quantile(values, 0.25, na.rm = TRUE) -
                                times * IQR(values, na.rm = TRUE)))

    ggplot(df, aes(x=group, y=values,color=group))+
      geom_violin(aes(fill=group),color="white",alpha=0.6)+
      geom_boxplot(width=0.05,position = position_dodge(0.9),color="white",outliers = F)+
      geom_jitter(size = 0.5,alpha = 1,width = 0.2)+
      geom_hline(aes(yintercept = hlines), linetype = "dashed", color = "gray",size = 0.8)+
      facet_wrap(~group,scales = "free",nrow = nrow)+
      labs(x = NULL,y = "Num of Proteins")+
      theme_classic()+
      theme(legend.position = "none")
  }
}
