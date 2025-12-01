#' @title dep_sp
#' @description
#' Calculate differentially expressed proteins among different groups in spatial proteomics (sp) objects
#' @param object an sp object \linkS4class{sp}.
#' @param slot There are two options in total:
#' \itemize{ \item clean_data :
#'   Perform differential calculation using the slot clean_data.
#'   \item RID_bootstrap : Perform differential calculation using the slot RID_bootstrap.}
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param catagory1 Corresponding to one or more groups in the group parameter,
#' this argument is used to select samples for comparison.
#' @param catagory2 Corresponding to one or more groups in the group parameter,
#' this argument is used to select samples for comparison.
#' @param method There are three options in total:
#' \itemize{ \item wilcox :
#'   Differentially expressed proteins were identified using the wilcox test.
#'   \item t_test : Differentially expressed proteins were identified using the t_test test.
#'   \item moderated_t : Differentially expressed proteins were identified using the wilcox moderated t_test from the limma package..}
#' @return A data frame containing statistical information on differentially expressed proteins.
#' @export
#' @rdname dep_sp
dep_sp <- function(object,slot = "clean_data",group,catagory1,catagory2,method = "wilcox"){

  slot <- match.arg(slot, c("clean_data","RID_bootstrap"))

  method <- match.arg(method, c("wilcox","t_test","moderated_t"))

  sam = object@sample_features[c("samples",group)] %>% na.omit()
  sam1 <- sam[sam[[group]]%in%catagory1,]
  sam1 <- sam1[["samples"]]
  sam2 <- sam[sam[[group]]%in%catagory2,]
  sam2 <- sam2[["samples"]]

  if (length(sam1) < 3 | length(sam2) < 3) {
    stop("If the number of samples in either the treatment group or the control group is less than 3, the analysis cannot proceed!")
  }

  cat("Sample size validation passed!\n")
  cat("the number of samples in the treatment group : ", length(sam1), "\n")
  cat("the number of samples in the control group : ", length(sam2), "\n")

  if(slot == "clean_data"){
    if(nrow(object@clean_data) < 1){
      stop("The clean_data slot is missing in the sp object.")
    }
    df <- object@clean_data
  }else{
    if(nrow(object@RID_bootstrap)  < 1){
      stop("The RID_bootstrap slot is missing in the sp object.")
    }
    df <- object@RID_bootstrap
  }

  df_treatmeant <- df[,sam1]
  df_control <- df[,sam2]

  if(method == "wilcox"){

    catagory1 <- paste0(catagory1,collapse="_")
    catagory2 <- paste0(catagory2,collapse="_")

    df_treatmeant <- df_treatmeant %>% t() %>% as.data.frame()
    df_control <- df_control %>% t() %>% as.data.frame()

    out1 <- union(as.numeric(which(apply(df_treatmeant,2,var)==0)),as.numeric(which(apply(df_control,2,var)==0)))
    out2 <- union(as.numeric(which(colSums(!is.na(df_treatmeant))<3)),as.numeric(which(colSums(!is.na(df_control))<3)))
    out <- union(out1,out2)

    if(length(out) != 0){
      df_treatmeant <- df_treatmeant[,-out]
      df_control <- df_control[,-out]
    }

    fc <- map2_dbl(df_treatmeant,df_control,~fc_sp(.x,.y))
    dep_wilcox <- map2(df_treatmeant,df_control,~wilcox.test(.x,.y))
    dep <- data.frame(names(dep_wilcox),map_dbl(dep_wilcox,"p.value"),fc,rep(paste0(catagory1,"/",catagory2),length(dep_wilcox)))
    colnames(dep) <- c("features","pvalue","logFC","case_control")
    dep[["fdr"]] <-p.adjust(dep[["pvalue"]],method = "BH")
  }else if(method == "t_test"){

    catagory1 <- paste0(catagory1,collapse="_")
    catagory2 <- paste0(catagory2,collapse="_")

    df_treatmeant <- df_treatmeant %>% t() %>% as.data.frame()
    df_control <- df_control %>% t() %>% as.data.frame()

    out1 <- union(as.numeric(which(apply(df_treatmeant,2,var)==0)),as.numeric(which(apply(df_control,2,var)==0)))
    out2 <- union(as.numeric(which(colSums(!is.na(df_treatmeant))<3)),as.numeric(which(colSums(!is.na(df_control))<3)))
    out <- union(out1,out2)

    if(length(out) != 0){
      df_treatmeant <- df_treatmeant[,-out]
      df_control <- df_control[,-out]
    }

    fc <- map2_dbl(df_treatmeant,df_control,~fc_sp(.x,.y))
    dep_t <- map2(df_treatmeant,df_control,~t.test(.x,.y, var.equal = T))
    dep <- data.frame(names(dep_t),map_dbl(dep_t,"p.value"),fc,rep(paste0(catagory1,"/",catagory2),length(dep_t)))
    colnames(dep) <- c("features","pvalue","logFC","case_control")
    dep[["fdr"]] <-p.adjust(dep[["pvalue"]],method = "BH")
  }else{
    catagory1 <- paste0(catagory1,collapse="_")
    catagory2 <- paste0(catagory2,collapse="_")
    group_list = c(rep(catagory1,length(sam1)),rep(catagory2,length(sam2)))
    design <- model.matrix(~0+factor(group_list))
    colnames(design)=levels(factor(group_list))
    contrast.matrix<-makeContrasts(contrasts = paste0(catagory1,"-",catagory2),levels = design)
    fit <- lmFit(cbind(df_treatmeant,df_control),design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    dep <-topTable(fit2, coef=1, n=Inf) %>% na.omit()
    dep[["features"]] <- row.names(dep)
    dep[["case_control"]] <- paste0(catagory1,"/",catagory2)
    dep <- dep[,c(7,4,1,8,5)]
    colnames(dep) <- c("features","pvalue","logFC","case_control","fdr")
  }
  return(dep)
}

#' @title dep_all_sp
#' @description
#' Calculate differentially expressed proteins among all groups in a specific column
#' of sample_features within the spatial proteomics (sp) object.
#' @param object an sp object \linkS4class{sp}.
#' @param slot There are two options in total:
#' \itemize{ \item clean_data :
#'   Perform differential calculation using the slot clean_data.
#'   \item RID_bootstrap : Perform differential calculation using the slot RID_bootstrap.}
#' @param group Group information of samples, which is a string corresponding
#' to one element in the column name of the sp object slot "sample_features".
#' The attribute of this column should be a character vector or a factor.
#' @param method There are three options in total:
#' \itemize{ \item wilcox :
#'   Differentially expressed proteins were identified using the wilcox test.
#'   \item t_test : Differentially expressed proteins were identified using the t_test test.
#'   \item moderated_t : Differentially expressed proteins were identified using the wilcox moderated t_test from the limma package..}
#' @return A data frame containing statistical information on differentially expressed proteins.
#' @export
#' @rdname dep_sp
dep_all_sp <- function(object,slot = "clean_data",group,method = "wilcox"){

  slot <- match.arg(slot, c("clean_data","RID_bootstrap"))

  method <- match.arg(method, c("wilcox","t_test","moderated_t","anova"))

  sam = object@sample_features[c("samples",group)] %>% na.omit()
  category_count <- table(sam[[group]])
  valid_categories <- names(category_count)[category_count >= 3]
  sam_filtered <- sam[sam[[group]] %in% valid_categories, ]
  if(nrow(sam_filtered) < 1){
    stop("Since each group contains fewer than three samples, differential protein analysis cannot be performed.")
  }

  if(method != "anova"){
    all_catagory <- unique(sam_filtered[[group]])
    deps <- list()
    for (i in all_catagory) {
      c1 <- i
      c2 <- all_catagory[all_catagory != i]
      deps[[i]] <- dep_sp(object,slot = slot,group=group,catagory1 = c1,catagory2 = c2,method = method)
    }
    deps = do.call(rbind,deps)
  }else{

    if(slot == "clean_data"){
      if(nrow(object@clean_data) < 1){
        stop("The clean_data slot is missing in the sp object.")
      }
      df <- object@clean_data
    }else{
      if(nrow(object@RID_bootstrap)  < 1){
        stop("The RID_bootstrap slot is missing in the sp object.")
      }
      df <- object@RID_bootstrap
    }
    df[["features"]] <- row.names(df)
    df <- df %>% pivot_longer(-features, names_to = "samples", values_to = "expression")
    df <- merge(df,sam_filtered,by = "samples", all.x = TRUE)
    colnames(df) <- c("samples","features","expression","group")
    deps <- anova_sp(df)
  }
  return(deps)
}

####Calculate the fold change####
fc_sp <- function(x,y){
  fc <- mean(x,na.rm=T) - mean(y,na.rm=T)
  return(fc)
}

####anova_mix####
anova_sp <- function(df){

  anova_results <- df %>%
    group_by(features) %>%
    do(anova_per_feature(.)) %>%
    ungroup()

  sig_proteins <- anova_results %>%
    filter(anova_p_adj < 0.05) %>%
    pull(features)

  sig_expr <- df %>%
    filter(features %in% sig_proteins)

  tukey_results <- sig_expr %>%
    group_by(features) %>%
    do(tukey_per_protein(.)) %>%
    ungroup()

  final_results <- tukey_results %>%
    left_join(anova_results[, c("features", "anova_p_adj")])

  colnames(final_results) <- c("features","contrast","diff","lwr","upr","p_adj","sig","p_adj_all")

  return(final_results)

}

####as a whole####
anova_per_feature <- function(protein_data) {
  # Levene
  levene_test <- leveneTest(expression ~ group, data = protein_data)
  levene_p <- levene_test$`Pr(>F)`[1]

  # ANOVA：ANOVA,Welch ANOVA
  if (levene_p > 0.05) {
    anova_res <- aov(expression ~ group, data = protein_data)
    anova_summary <- summary(anova_res)
    anova_p <- anova_summary[[1]]$`Pr(>F)`[1] # ANOVA P
  } else {
    anova_res <- oneway.test(expression ~ group, data = protein_data, var.equal = FALSE)
    anova_p <- anova_res$p.value # Welch ANOVA P
  }

  return(data.frame(
    levene_p = levene_p,
    anova_p = anova_p,
    anova_p_adj = p.adjust(anova_p, method = "fdr") # FDR
  ))
}

####Each group####
tukey_per_protein <- function(protein_data) {
  # ANOVA
  anova_res <- aov(expression ~ group, data = protein_data)
  # Tukey HSD
  tukey_res <- TukeyHSD(anova_res)$group %>% as.data.frame()
  tukey_res[["contrast"]] <- row.names(tukey_res)
  tukey_res <- tukey_res %>% mutate(
      sig = ifelse(`p adj` < 0.05, "significant", "ns")
    )
  return(tukey_res)
}

