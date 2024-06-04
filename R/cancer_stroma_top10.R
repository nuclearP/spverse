rm(list = ls())

library(tidyverse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggsci)

removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}

cs <- fread("D:/software/R/Rdata/fly2/lfq/ov_cs/cs_lfq.txt")
x <- fread("D:/software/R/Rdata/fly2/lfq/ov_cs/msstat_input.tsv")

###################cs_pca#################

cs <- separate(cs,"V1",c("p","b"),sep = ";") %>% 
  dplyr::select(-b)

cs <- column_to_rownames(cs,"p")


x <- dplyr::select(x,c("Condition","BioReplicate","Run"))

x <- distinct(x)

x <- mutate(x,sample = paste0(x$Condition,x$BioReplicate))

colnames(cs) <- x$sample

y <- colSums(!is.na(cs)) %>% as.data.frame() %>% set_names(c("p")) %>% dplyr::filter(p > 4000)

cs <- dplyr::select(cs,row.names(y))

cs <- removeRowsAllNa(cs)

cs[is.na(cs)] <- min(cs,na.rm = T)

set.seed(1234)
umap_fit <- cs %>% t() %>% 
  scale() %>% 
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  mutate(sample = colnames(cs))%>%
  inner_join(x, by="sample")

umaps <- ggplot(umap_df,aes(x = V1, 
                            y = V2, 
                            color = Condition))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP",
       color="Condition")+
  scale_color_aaas()+
  theme_bw()

umaps

ggsave("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/umaps_pn.png",umaps)
ggsave("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/umaps_pn.pdf",umaps)


######################pca######################
set.seed(1234)
pca_fit <- cs %>% t() %>% 
  scale() %>%
  prcomp()
# 查看成分重要性
summary(pca_fit)

# 可视化PC1和PC2
summ1 <- summary(pca_fit)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

df1 <- pca_fit$x # 提取PC score
df1 <- as.data.frame(df1) 
df1 <- mutate(df1,sample=row.names(df1)) %>% 
  inner_join(x,by="sample")


p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color=Condition))+
  # stat_ellipse(aes(fill = pca$x),
  #              type = "norm",geom = "polygon",alpha = 0.25,color=NA)+ # 添加置信椭圆
  geom_point()+
  
  labs(x = xlab1,y = ylab1,color = "patients",title = "PCA")+
  theme_bw()+
  scale_color_aaas()

p.pca1

ggsave(p.pca1,filename = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/PCA_pn.png")
ggsave(p.pca1,filename = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/PCA_pn.pdf")

################鉴定数################

cs <- fread("D:/software/R/Rdata/fly2/lfq/ov_cs/cs_lfq.txt")
x <- fread("D:/software/R/Rdata/fly2/lfq/ov_cs/msstat_input.tsv")

cs <- separate(cs,"V1",c("p","b"),sep = ";") %>% 
  dplyr::select(-b)

cs <- column_to_rownames(cs,"p")


x <- dplyr::select(x,c("Condition","BioReplicate","Run"))

x <- distinct(x)

x <- mutate(x,sample = paste0(x$Condition,x$BioReplicate))

colnames(cs) <- x$sample

y <- colSums(!is.na(cs)) %>% as.data.frame()

colnames(y) <- "proteins"

y <- mutate(y,sample = row.names(y)) %>% left_join(x)

y <- filter(y,y$proteins > 4000)

p1=ggplot(y, aes(x=Condition, y=proteins,color=Condition))+
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.3,size=1)+
  geom_boxplot(aes(),notch = F,size=0.7)+
  geom_jitter(size = 2.5,alpha = 0.6,width = 0.2)+ 
  theme_bw()+scale_color_npg()+scale_fill_npg()+

  theme(axis.text=element_text(colour='black',size=11),axis.title = element_text(colour='black',size=15))
p1
ggsave(filename = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/proteins.png")
ggsave(filename = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/proteins.pdf")


y <- arrange(y,y$Condition)
cs <- dplyr::select(cs,y$sample)
z <- y %>% group_by(Condition) %>% count()


c <- pivot_longer(cs,everything())
c <- drop_na(c)




p <- ggplot(c,mapping =aes(x=value))+
  geom_density(mapping = aes(color= name))+
  scale_colour_manual(values = c(rep("#208D85",z[1,2]),rep("#208A09",z[2,2]),rep("#272E6A",z[3,2]) ,  rep("#8A9FD1",z[4,2])))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 25,face = "bold"))+
  theme(axis.text = element_text(size = 20,face = "bold"))+theme(axis.ticks = element_line(linewidth = 1.5),axis.line = element_line(linewidth = 1.5),axis.ticks.length = unit(0.2, "cm") )+
  annotate("segment",x = 11.5 ,y= 0.4,xend = 12,yend = 0.4,colour = "#208D85",linewidth = 2)+
  annotate("text",x = 12.6 ,y= 0.4,colour = "black",label= "Neg_cancer",size = 4)+
  annotate("segment",x = 13.2 ,y= 0.4,xend = 13.7,yend = 0.4,colour = "#272E6A",linewidth = 2)+
  annotate("text",x = 14.3 ,y= 0.4,colour = "black",label= "Neg_stroma",size = 4)+
  annotate("segment",x = 11.5 ,y= 0.37,xend = 12,yend = 0.37,colour = "#208A42",linewidth = 2)+
  annotate("text",x = 12.6 ,y= 0.37,colour = "black",label= "Pos_cancer",size = 4)+
  annotate("segment",x = 13.2 ,y= 0.37,xend = 13.7,yend = 0.37,colour = "#8A9FD1",linewidth = 2)+
  annotate("text",x = 14.3 ,y= 0.37,colour = "black",label= "Pos_stroma",size = 4)+
  labs(y = "Density")

# xlim(c(0,15))

p

ggsave("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/distribution.png",height = 7,width = 9)

ggsave("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/distribution.pdf",height = 7,width = 9)
################top10_mapping######################

rm(list = ls())

load(file = "D:/software/R/Rdata/fly2/iq_test/lfq_kmean/dep/dep.Rdata")

pro_gene <- read.delim("D:/software/R/Rdata/ch/p_g.xls")

x <- test@comdata %>% mutate(p = row.names(.))

x <- left_join(x,pro_gene) %>% dplyr::select(-"p") %>% filter(!is.na(gene)) %>% group_by(gene) %>%
  summarise(across(where(is.numeric),mean)) %>% column_to_rownames("gene")

x <- x %>%  t() %>% as.data.frame() %>% mutate(samples = row.names(.))

y <- data.frame(sample=1:100,x=rep(1:10,10),y=rep(10:1,each=10))

zhc <- list()

for (i in c("N","P")) {
  for (j in 1:5) {
    z <- y %>% mutate(y,samples = paste0(i,j,"_",sample))
    zhc[[paste0(i,j)]] <- z 
  }
}

zhc <- rbindlist(zhc)

load("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/lfq_data/cancer_stroma_dep.Rdata")
depn_c <- dplyr::filter(depn_c,`adj.P.Val` < 0.0001) %>% slice_max(logFC,n=10) %>% pull(p) 
depn_s <- dplyr::filter(depn_s,`adj.P.Val` < 0.0001) %>% slice_min(logFC,n=10) %>% pull(p) 
depp_c <- dplyr::filter(depp_c,`adj.P.Val` < 0.0001) %>% slice_max(logFC,n=10) %>% pull(p) 
depp_s <- dplyr::filter(depp_s,`adj.P.Val` < 0.0001) %>% slice_min(logFC,n=10) %>% pull(p) 
marker <- c(depn_c,depn_s,depp_c,depp_s)

x <- dplyr::select(x,any_of(unique(marker)),samples)

fi <- full_join(zhc,x) %>% mutate(face = str_extract(samples,"[A-Z0-9]+"))

set.seed(1234)
pca_fit <- test@comdata %>% t() %>% 
  scale() %>%
  prcomp(center = T)
# 查看成分重要性
summary(pca_fit)

# 可视化PC1和PC2
summ1 <- summary(pca_fit)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

df1 <- pca_fit$x # 提取PC score
df1 <- as.data.frame(df1) 
df1 <- mutate(df1,samples=row.names(df1)) %>% 
  inner_join(x,by="samples")


for(i in 5:29){
  z1 <- ggplot(fi,mapping = aes(x,y))+
    # geom_tile(data =y,aes(x=x,y=y),fill = "grey50",color = "black",size=0.8)+
    geom_tile(data =fi,aes(x=x,y=y,fill=fi[[i]]),color="black",size=0.6)+
    facet_wrap(~face,nrow = 2)+
    # scale_fill_gradientn(colours =colorRampPalette(c("#2556A6","#ffffff","#EE2A29"))(50))+
    # scale_fill_manual(values = c("#E8C5E0","#FFD217","#002E9F"))+
    # scale_fill_gradientn(colours =colorRampPalette(c("#2556A6","#ffffff","#EE2A29"))(50),limit=c(0.5,0.95),values = scales::rescale(c(0.5, 0.75, 0.9)))+
    # scale_fill_gradientn(colours = brewer.pal(11,'YlGnBu'))+
    scale_fill_viridis(option = "H")+
    labs(fill = colnames(fi)[i])+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_x_continuous(name = c(1:10),breaks = 1:10)+
    scale_y_continuous(name = c(1:10),breaks = 1:10)+
    theme(axis.ticks = element_blank(),axis.title = element_blank(),axis.text = element_text(size = 9),axis.line = element_blank(),
          strip.background = element_rect(fill = c("#e7e7e7")), #修改分页背景
          strip.text.x = element_text(size = 8),legend.position = "bottom")+
    guides(fill = guide_colorbar(barwidth= 13))+
    coord_fixed()
  z1
  ggsave(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/top10/",colnames(fi)[i],"_map.png"),width = 9,height = 7)
  ggsave(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/top10/",colnames(fi)[i],"_map.pdf"),width = 9,height = 7)
  
  
  p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color=df1[[967+i]]))+
    # stat_ellipse(aes(fill = pca$x),
    #              type = "norm",geom = "polygon",alpha = 0.25,color=NA)+ # 添加置信椭圆
    geom_point()+
    
    labs(x = xlab1,y = ylab1,color =colnames(fi)[i] ,title = "PCA")+
    theme_bw()+
    scale_color_viridis(option = "H")
  
  p.pca1
  ggsave(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/top10/",colnames(fi)[i],"_pca.png"))
  ggsave(paste0("D:/software/R/Rdata/fly2/iq_test/lfq_kmean/cancer_scores/top10/",colnames(fi)[i],"_pca.pdf"))
}
