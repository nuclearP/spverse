rm(list = ls())

library(tidyverse)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(ggsci)
library(sva)
library(ComplexHeatmap)
library(umap)
library(circlize)
library(Seurat)
library(limma)
library(NbClust)
library(cluster)
library(factoextra)
library(jjAnno)
library(ggpubr)
library(CellChat)
library(RColorBrewer)
library(corrplot)

###############rank###############
load(file = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/dep/dep.Rdata")

comdata <- test@comdata
cluster <- test@cluster

comdata <- t(comdata) %>% as.data.frame() %>% mutate(clusters = paste0("cluster",cluster))

comdata <- comdata %>% group_by(clusters) %>% summarise(across(where(is.numeric),mean))
comdata <- column_to_rownames(comdata,"clusters") %>% t() %>% as.data.frame()

c <- cor(comdata)

col_fun=colorRamp2(c(0.65,0.7,0.75,0.8,0.85,0.9,0.95,1), colorRampPalette(rev(brewer.pal(11, 'Spectral')))(8) )

zhc <- Heatmap(c,  name = "Pearson correlation",
               show_row_names = T,
               show_column_names = T,
               col = col_fun,
               # column_dend_height=unit(30,"mm"),
               # row_dend_width = unit(30, "mm"),
               # show_row_dend = F,
               row_title_gp = gpar(fontsize=10),
               # row_names_max_width = max_text_width(colnames(gsvas)),
               column_title =NULL,column_title_gp =  gpar(fontsize=5),
               
               heatmap_legend_param = list(
                 legend_height=unit(4,"cm"), legend_direction="vertical",title_gp = gpar(fontsize = 12),
                 labels_gp =  gpar(fontsize = 12),title_position = "topleft"
                 # title_position = "topcenter"
               )
)

png(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cluster_cor/","cluster_cor.png"),width = 4000,height =3500,res = 300)
print(draw(zhc,merge_legend = TRUE))
dev.off()

pdf(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cluster_cor/","cluster_cor.pdf"),width = 9,height =8)
draw(zhc,merge_legend = TRUE)
dev.off()


corrplot(c,
         type = 'lower',
         method = 'color',
         # order = "AOE",
         # method = "color",            # 图案形状 "square"方框,"circle"圆, "ellipse"椭圆, "number"数字, "shade"阴影花纹, "color"颜色方框, "pie饼图"
         # type = "lower",               # "full", "lower", "upper"
         col = colorRampPalette(brewer.pal(9,"YlGnBu")[1:9])(50),
         addCoef.col = "black",
         bg = "white",                 # 背景颜色
         col.lim = c(min(c),max(c)),
         tl.srt = 45, 
         # tl.offset = 0.8, 
         cl.ratio = 0.08,
         cl.cex = 1.5,
         # tl.pos = "d",
         tl.col = "black",
         # number.cex = 2,
         # add = T,                      # 是否在原来的图层上添加图形
         # diag = T                # 是否显示主对角
         # addCoefasPercent = F,         # 是否把相关性数值改为百分数
         is.corr = F
         # order = "original",
         # tl.col = "black",
         # outline = "black"
         
)

###################cluster_np_cor###################

x <- test@rawdata
y <- test@description

x <- t(x) %>% as.data.frame() %>% mutate(pn = y$pn,cluster = y$cluster)

for (i in 1:9) {
  z <- dplyr::filter(x,cluster == i)
  k <- z %>% dplyr::select(-cluster,-pn) %>% t()
  corel <- cor(k,method = "pearson",use = "pairwise.complete.obs")
  
  col_fun=colorRamp2(c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1), colorRampPalette(rev(brewer.pal(11, 'Spectral')))(9) )
  
  top_anno <- HeatmapAnnotation(
    
    postive_negative = z$pn,
    col = list(postive_negative=c("P"="#B24745FF","N"="#6A6599FF")),
    show_legend = T,
    # how_annotation_name = F,
    annotation_legend_param = list(
      legend_height=unit(4,"cm"), legend_direction="vertical",title_gp = gpar(fontsize = 12),
      labels_gp =  gpar(fontsize = 12),title_position = "topleft"
    )
  )
  
  left_anno <- rowAnnotation(
    
    postive_negative = z$pn,
    col = list(postive_negative=c("P"="#B24745FF","N"="#6A6599FF")),
    show_legend = F,
    # how_annotation_name = F,
    annotation_legend_param = list(
      legend_height=unit(4,"cm"), legend_direction="vertical",title_gp = gpar(fontsize = 12),
      labels_gp =  gpar(fontsize = 12),title_position = "topleft"
    )
  )
  
  set.seed(1)
  
  
  zhc <- Heatmap(corel,  name = "Pearson correlation",
                 show_row_names = F,
                 show_column_names = F,
                 col = col_fun,
                 column_dend_height=unit(30,"mm"),
                 row_dend_width = unit(30, "mm"),
                 # show_row_dend = F,
                 top_annotation =top_anno,
                 left_annotation = left_anno,
                 row_title_gp = gpar(fontsize=10),
                 # row_names_max_width = max_text_width(colnames(gsvas)),
                 column_title =NULL,column_title_gp =  gpar(fontsize=5),
                 
                 heatmap_legend_param = list(
                   legend_height=unit(4,"cm"), legend_direction="vertical",title_gp = gpar(fontsize = 12),
                   labels_gp =  gpar(fontsize = 12),title_position = "topleft"
                   # title_position = "topcenter"
                 )
  )
  
  png(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cluster_cor/","cluster_",i,"formal.png"),width = 4000,height =3500,res = 300)
  print(draw(zhc,merge_legend = TRUE))
  dev.off()
  
  pdf(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cluster_cor/","cluster_",i,"formal.pdf"),width = 9,height =8)
  draw(zhc,merge_legend = TRUE)
  dev.off()
  
}

###################cluster_np_cor###################
rm(list = ls())
load(file = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/dep/dep.Rdata")
x <- test@rawdata
y <- test@description

x <- t(x) %>% as.data.frame() %>% mutate(cluster = y$cluster)

