setwd("/data1/zhuzn/wangsn/reviewer/")
rm(list=ls())

# library -----------------------------------------------------------------
library(Seurat)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(scater)
library(cowplot)
library(Matrix)
library(reshape2)
library(tidyverse)
library(SingleCellExperiment)
#library(destiny)
library(umap)
library(ggthemes)
library(patchwork)
library(scales)
library(magrittr)
library(ggsci)
library(magick)
library(devtools)
library(stringr)
library(scales)

# figureS1A 50 SASP_tissue-合并 ----------------------------------------------------------------------
rm(list=ls())
library(ggridges)
library(ggpubr)
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined1 <- all_sample.combined
data.list_Tissue <- SplitObject(all_sample.combined1, split.by = "Tissue")
sample_list <- unique(all_sample.combined@meta.data$Tissue)
Sasp <- data.frame()
Sasp2 <- data.frame()
Sasp3 <- data.frame()
for (Tissue in sample_list){
  all_sample.combined <- data.list_Tissue[[Tissue]]
  Idents(all_sample.combined)
  DefaultAssay(all_sample.combined) <- "RNA"
  DefaultAssay(all_sample.combined)
  sample.integrated <- all_sample.combined
  SASP_gene_set <- read.delim('/data1/zhuzn/wangsn/SASP/SASP_V2.list',col.names = F)
  gene <- as.list(SASP_gene_set)
  sample.integrated <- AddModuleScore(object = sample.integrated,features = gene,ctrl = 100, name = 'SASP_Features',replace = T)
  
  names(sample.integrated@meta.data)
  sasp <- sample.integrated@meta.data[,c(5,19)]
  
  sasp2 <-aggregate(sasp$SASP_Features1,list(sasp$Age),mean)
  names(sasp) <- c('Group','Gene_set_score')
  sasp$Tissue <- rep(Tissue)
  names(sasp2) <- c('Group','mean')
  sasp2$Tissue <- rep(Tissue)
  # sasp$Group <- factor(sasp$Group,levels = c("young","old"))
  sasp$Group <- factor(sasp$Group,levels = c("old","young"))
  sasp3 <- compare_means(Gene_set_score~Group, data=sasp,method = "t.test")
  sasp3$Tissue <- rep(Tissue)
  
  Sasp <- rbind(Sasp,sasp)
  Sasp2 <- rbind(Sasp2,sasp2)
  Sasp3 <- rbind(Sasp3,sasp3)
}
library(dplyr)
library(forcats)
p2 <- Sasp %>%
  mutate(TissueFct = fct_rev(as.factor(Tissue))) %>%
  ggplot(aes(y = TissueFct)) +
  geom_density_ridges(aes(x = Gene_set_score, fill = paste(TissueFct, Group)),color = "grey", from = -0.2, to = 0.3) +
  labs(x = "Gene_set_score",y = "",title = "Tissue",subtitle = "",caption = "") +
  scale_y_discrete(expand = c(0.01, 0)) +scale_x_continuous(expand = c(0.01, 0)) +
  scale_fill_cyclical(
    breaks = c("Bone_marrow young", "Bone_marrow old","Mesenteric_lymph young", "Mesenteric_lymph old","PBMC young", "PBMC old","Spleen young", "Spleen old"),
    labels = c(`Bone_marrow young` = "young",`Bone_marrow old` = "old",`Mesenteric_lymph young` = "young",`Mesenteric_lymph old` = "old",`PBMC young` = "young", `PBMC old` = "old",`Spleen young` = "young", `Spleen old` = "old"),
    values = c('#35D168','#238b45','#DF96C3','#df65b0','#93CAF8','#2171b5','#CC8F6B','#cc4c02'),
    name = "Age", guide = "legend") +
  # geom_vline(data = Sasp2, aes(xintercept = mean, color=Tissue),linetype='dashed')+
  theme_ridges()+
  theme(axis.text = element_text(size = 20),legend.text = element_text(size = 18),legend.key.size = unit(20,"pt"),legend.key = element_rect(colour='grey32'))
p2
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figureS1/',"figure.",'S1A.pdf')
ggsave(plotfile, plot=p2, dpi = 1000, width = 10, height = 7)

# figureS1 FindMarkers-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$Tissue)
diff_gene <- FindAllMarkers(all_sample.combined,logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                       'NOT'))
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","tissue_others_tissue","_DE.txt"),quote = FALSE,sep = "\t",row.names = F)

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$Age)

data.list_tissue <- SplitObject(all_sample.combined, split.by = "Tissue")
sample_list <- unique(all_sample.combined@meta.data$Tissue)
for (tissue in sample_list){
  tissue_data <- data.list_tissue[[tissue]]
  
  diff_gene <- FindMarkers(tissue_data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
  diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
  diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                         ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                         'NOT'))
  write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/",tissue,"_old_vs_young.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
}

# figureS1BC DEGs数量_UpSetR-------------------------------------------------------------------------
rm(list=ls())
sample_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
Bone_marrow <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Bone_marrow","_old_vs_young.findmark.txt"))
Mesenteric_lymph <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Mesenteric_lymph","_old_vs_young.findmark.txt"))
PBMC <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","PBMC","_old_vs_young.findmark.txt"))
Spleen <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Spleen","_old_vs_young.findmark.txt"))

#install.packages('UpSetR')
library('UpSetR')
#取出准备取交集的数据集们
listinput_UP <- list("Bone_marrow" = Bone_marrow[Bone_marrow$change_V1 == "UP",6],
                     "Mesenteric_lymph" = Mesenteric_lymph[Mesenteric_lymph$change_V1 == "UP",6],
                     "PBMC" = PBMC[PBMC$change_V1 == "UP",6],
                     "Spleen" = Spleen[Spleen$change_V1 == "UP",6])
p1<-upset(fromList(listinput_UP),nsets = length(listinput_UP), order.by = "freq",
          sets = c("Spleen","PBMC","Mesenteric_lymph","Bone_marrow"),keep.order = TRUE,
          mainbar.y.label="Intersection of UP genes", matrix.color = "red",
          sets.bar.color = c('#238b45','#df65b0','#2171b5','#cc4c02'),
          line.size = 0.7,sets.x.label = "UP DEGs",text.scale = c(1.2, 1.1, 1.1, 1.1, 1.2, 1))
p1
pdf(file = paste0('/data1/zhuzn/wangsn/reviewer_4000/figureS1/',"figure.",'S1B.pdf'), p1, width = 7, height = 5)
print(p1)
dev.off()

#取出准备取交集的数据集们
listinput_DOWN <- list("Bone_marrow" = Bone_marrow[Bone_marrow$change_V1 == "DOWN",6],
                       "Mesenteric_lymph" = Mesenteric_lymph[Mesenteric_lymph$change_V1 == "DOWN",6],
                       "PBMC" = PBMC[PBMC$change_V1 == "DOWN",6],
                       "Spleen" = Spleen[Spleen$change_V1 == "DOWN",6])
p2<-upset(fromList(listinput_DOWN),nsets = length(listinput_DOWN), order.by = "freq",
          sets = c("Spleen","PBMC","Mesenteric_lymph","Bone_marrow"),keep.order = TRUE,
          mainbar.y.label="Intersection of DOWN genes", matrix.color = "red",
          sets.bar.color = c('#238b45','#df65b0','#cc4c02','#2171b5'),
          line.size = 0.7,sets.x.label = "DOWN DEGs",text.scale = c(1.2, 1.1, 1.1, 1.1, 1.2, 1))
p2
pdf(file = paste0('/data1/zhuzn/wangsn/reviewer_4000/figureS1/',"figure.",'S1C.pdf'), p2, width = 7, height = 5)
print(p2)
dev.off()

# figureS1D DE_union（上调-共同-特异）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library(reshape2)
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
Tissue <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figureS1/tissue_others_tissue_DE.txt", sep="") %>% .[.$change_V1 == "UP",]
Tissue$gene <- gsub("\\..*","",Tissue$gene)
Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Bone_marrow","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Mesenteric_lymph <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Mesenteric_lymph","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Mesenteric_lymph <- Mesenteric_lymph[!grepl("ENSMMUG", Mesenteric_lymph$gene),]
PBMC <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","PBMC","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
PBMC <- PBMC[!grepl("ENSMMUG", PBMC$gene),]
Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Spleen","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Spleen <- Spleen[!grepl("ENSMMUG", Spleen$gene),]

DE <- rbind(Bone_marrow,Mesenteric_lymph,PBMC,Spleen)
gene_DE <- data.frame(table(DE$gene));gene_DE <- gene_DE[gene_DE$Freq == "4",]
gene_union_Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Bone_marrow$tissue <- rep("Bone_marrow")
gene_union_Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Mesenteric_lymph$tissue <- rep("Mesenteric_lymph")
gene_union_PBMC <- PBMC[PBMC$gene %in% gene_DE$Var1,c(2,6)]
gene_union_PBMC$tissue <- rep("PBMC")
gene_union_Spleen <- Spleen[Spleen$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Spleen$tissue <- rep("Spleen")
gene_union <- rbind(gene_union_Bone_marrow,gene_union_Mesenteric_lymph,gene_union_PBMC,gene_union_Spleen)

overlap_FC <- gene_union
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$tissue <- factor(overlap_FC$tissue,levels = c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen"))
overlap_FC1 <- overlap_FC

Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]

Tissue_Bone_marrow <- Tissue[Tissue$cluster == "Bone_marrow",]
Bone_marrow1 <- Bone_marrow1[Bone_marrow1$gene %in% intersect(Bone_marrow1$gene,Tissue_Bone_marrow$gene),]
Tissue_Mesenteric_lymph <- Tissue[Tissue$cluster == "Mesenteric_lymph",]
Mesenteric_lymph1 <- Mesenteric_lymph1[Mesenteric_lymph1$gene %in% intersect(Mesenteric_lymph1$gene,Tissue_Mesenteric_lymph$gene),]
Tissue_PBMC <- Tissue[Tissue$cluster == "PBMC",]
PBMC1 <- PBMC1[PBMC1$gene %in% intersect(PBMC1$gene,Tissue_PBMC$gene),]
Tissue_Spleen <- Tissue[Tissue$cluster == "Spleen",]
Spleen1 <- Spleen1[Spleen1$gene %in% intersect(Spleen1$gene,Tissue_Spleen$gene),]

diff_gene <- c(Bone_marrow1$gene,Mesenteric_lymph1$gene,PBMC1$gene,Spleen1$gene)
diff_gene_merge <- unique(c(as.character(gene_DE$Var1),diff_gene))
diff_gene_overlap_FC <- data.frame()
DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
for (j in DE_list){
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/",j,"_old_vs_young.findmark.txt"), header=T)
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
  diff_gene1$tissue <- rep(j)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}
diff_gene_overlap_FC <- diff_gene_overlap_FC[!grepl("ENSMMUG", diff_gene_overlap_FC$gene),]
overlap_FC <- dcast(diff_gene_overlap_FC,tissue~diff_gene_overlap_FC$gene,value.var = "avg_log2FC",fill = "0")
overlap_FC <- melt(overlap_FC,id.vars = "tissue",variable.name='gene',value.name="avg_log2FC")
colnames(overlap_FC) <- c("tissue","gene","avg_log2FC")
str(overlap_FC)
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$tissue <- factor(overlap_FC$tissue,levels = DE_list)
overlap_FC2 <- overlap_FC
overlap_FC2 <- overlap_FC2[,c(3,2,1)]
overlap_FC <- rbind(overlap_FC1,overlap_FC2)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene_merge))
overlap_FC[overlap_FC > 1] <- 1
colors <- c("#154889","#4359A7","#39BBF9","#73CEF9","#C3E8F9","#EAF0F5","#F7ADAD","#F76C6C","#F73E3E","#F70000")
show_col(colors)
p <- ggplot(overlap_FC,aes(x=tissue,y=gene,fill=avg_log2FC))+
  geom_raster()+
  # scale_fill_gradientn(colors = c("#154889","#1A58A7","#4592E5","#eff3ff","#F76C6C","#F70000","#981216"))+
  # scale_fill_gradientn(colors = c("#154889","#39BBF9","#73CEF9","#C3E8F9","#EAF0F5","#F7ADAD","#F76C6C","#F76C6C","#F70000"))+
  # scale_fill_gradient2(low="#154889",mid="#EAF0F5",high ="#F70000",midpoint = 0.1)+
  scale_fill_gradient2(low="#4592E5",mid="#eff3ff",high ="#F70000",midpoint = 0)+
  theme_bw() + 
  labs(x = "",y = "",title = "Tisssue_old_vs_young_UP_setdiff_union_merge")+ 
  theme(legend.title = element_text(size = 15),legend.text = element_text(size = 20),
        panel.grid = element_blank(),plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 20),axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.border = element_blank())
# theme(panel.grid = element_blank(),plot.title = element_text(size = 10),axis.text.y = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","figureS1D.pdf"),plot = p,width = 10,height = 15)

# figureS1E DE_union（下调-共同-特异）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

Tissue <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figureS1/tissue_others_tissue_DE_old.txt", sep="") %>% .[.$change_V1 == "DOWN",]
Tissue$gene <- gsub("\\..*","",Tissue$gene)

Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Bone_marrow","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Mesenteric_lymph <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Mesenteric_lymph","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
Mesenteric_lymph <- Mesenteric_lymph[!grepl("ENSMMUG", Mesenteric_lymph$gene),]
PBMC <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","PBMC","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
PBMC <- PBMC[!grepl("ENSMMUG", PBMC$gene),]
Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Spleen","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
Spleen <- Spleen[!grepl("ENSMMUG", Spleen$gene),]

DE <- rbind(Bone_marrow,Mesenteric_lymph,PBMC,Spleen)
gene_DE <- data.frame(table(DE$gene));gene_DE <- gene_DE[gene_DE$Freq == "4",]
gene_union_Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Bone_marrow$tissue <- rep("Bone_marrow")
gene_union_Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Mesenteric_lymph$tissue <- rep("Mesenteric_lymph")
gene_union_PBMC <- PBMC[PBMC$gene %in% gene_DE$Var1,c(2,6)]
gene_union_PBMC$tissue <- rep("PBMC")
gene_union_Spleen <- Spleen[Spleen$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Spleen$tissue <- rep("Spleen")
gene_union <- rbind(gene_union_Bone_marrow,gene_union_Mesenteric_lymph,gene_union_PBMC,gene_union_Spleen)

overlap_FC <- gene_union
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$tissue <- factor(overlap_FC$tissue,levels = c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen"))
overlap_FC1 <- overlap_FC

Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]

Tissue_Bone_marrow <- Tissue[Tissue$cluster == "Bone_marrow",]
Bone_marrow1 <- Bone_marrow1[Bone_marrow1$gene %in% intersect(Bone_marrow1$gene,Tissue_Bone_marrow$gene),]
Tissue_Mesenteric_lymph <- Tissue[Tissue$cluster == "Mesenteric_lymph",]
Mesenteric_lymph1 <- Mesenteric_lymph1[Mesenteric_lymph1$gene %in% intersect(Mesenteric_lymph1$gene,Tissue_Mesenteric_lymph$gene),]
Tissue_PBMC <- Tissue[Tissue$cluster == "PBMC",]
PBMC1 <- PBMC1[PBMC1$gene %in% intersect(PBMC1$gene,Tissue_PBMC$gene),]
Tissue_Spleen <- Tissue[Tissue$cluster == "Spleen",]
Spleen1 <- Spleen1[Spleen1$gene %in% intersect(Spleen1$gene,Tissue_Spleen$gene),]

diff_gene <- c(Bone_marrow1$gene,Mesenteric_lymph1$gene,PBMC1$gene,Spleen1$gene)
diff_gene_merge <- c(as.character(gene_DE$Var1),diff_gene)
diff_gene_overlap_FC <- data.frame()
DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
for (j in DE_list){
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/",j,"_old_vs_young.findmark.txt"), header=T)
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
  diff_gene1$tissue <- rep(j)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}
diff_gene_overlap_FC <- diff_gene_overlap_FC[!grepl("ENSMMUG", diff_gene_overlap_FC$gene),]
overlap_FC <- dcast(diff_gene_overlap_FC,tissue~diff_gene_overlap_FC$gene,value.var = "avg_log2FC",fill = "0")
overlap_FC <- melt(overlap_FC,id.vars = "tissue",variable.name='gene',value.name="avg_log2FC")
colnames(overlap_FC) <- c("tissue","gene","avg_log2FC")
str(overlap_FC)
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$tissue <- factor(overlap_FC$tissue,levels = DE_list)
overlap_FC2 <- overlap_FC
overlap_FC2 <- overlap_FC2[,c(3,2,1)]
overlap_FC <- rbind(overlap_FC1,overlap_FC2)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene_merge))
overlap_FC[overlap_FC < -1] <- -1
colors <- c("#154889","#1A58A7","#4592E5","#eff3ff","#eff3ff","#F76C6C","#F70000")
show_col(colors)
p <- ggplot(overlap_FC,aes(x=tissue,y=gene,fill=avg_log2FC))+
  geom_raster()+
  # scale_fill_gradientn(colors = c("#154889","#1A58A7","#4592E5","#eff3ff","#F76C6C","#F70000"))+
  # scale_fill_gradientn(colors = c("#981216","#F70000","#F76C6C","#eff3ff","#4592E5","#1A58A7","#154889"))+
  scale_fill_gradient2(low="#4592E5",mid="#eff3ff",high ="#F70000",midpoint = 0)+
  # scale_fill_gradient2(low="#154889",mid="#EAF0F5",high ="#F70000",midpoint = -0.1)+
  theme_bw() + 
  labs(x = "",y = "",title = "Tisssue_old_vs_young_UP_setdiff_union_merge")+ 
  theme(legend.title = element_text(size = 15),legend.text = element_text(size = 20),
        panel.grid = element_blank(),plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 20),axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.border = element_blank())
# theme(panel.grid = element_blank(),plot.title = element_text(size = 10),axis.text.y = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","figureS1E.pdf"),plot = p,width = 10,height = 15)

# figureS1D DE_union（上调-GO富集）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Bone_marrow","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Mesenteric_lymph <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Mesenteric_lymph","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Mesenteric_lymph <- Mesenteric_lymph[!grepl("ENSMMUG", Mesenteric_lymph$gene),]
PBMC <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","PBMC","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
PBMC <- PBMC[!grepl("ENSMMUG", PBMC$gene),]
Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Spleen","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Spleen <- Spleen[!grepl("ENSMMUG", Spleen$gene),]
DE <- rbind(Bone_marrow,Mesenteric_lymph,PBMC,Spleen)
gene_DE <- data.frame(table(DE$gene));gene_DE <- gene_DE[gene_DE$Freq == "4",]
gene_union_Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Bone_marrow$tissue <- rep("Bone_marrow")
gene_union_Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Mesenteric_lymph$tissue <- rep("Mesenteric_lymph")
gene_union_PBMC <- PBMC[PBMC$gene %in% gene_DE$Var1,c(2,6)]
gene_union_PBMC$tissue <- rep("PBMC")
gene_union_Spleen <- Spleen[Spleen$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Spleen$tissue <- rep("Spleen")
gene_union <- rbind(gene_union_Bone_marrow,gene_union_Mesenteric_lymph,gene_union_PBMC,gene_union_Spleen)
colnames(gene_union)[2] <- "SYMBOL"
df_id <- bitr(gene_union$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
gene_union <- merge(gene_union,df_id,by = "SYMBOL",all=F)
gene <- unique(gene_union[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","union","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
Bone_marrow1 <- Bone_marrow1[order(Bone_marrow1$avg_log2FC,decreasing = T),]
colnames(Bone_marrow1)[6] <- "SYMBOL"
df_id <- bitr(Bone_marrow1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Bone_marrow1 <- merge(Bone_marrow1,df_id,by = "SYMBOL",all=F)
gene <- unique(Bone_marrow1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","Bone_marrow","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
Mesenteric_lymph1 <- Mesenteric_lymph1[order(Mesenteric_lymph1$avg_log2FC,decreasing = T),]
colnames(Mesenteric_lymph1)[6] <- "SYMBOL"
df_id <- bitr(Mesenteric_lymph1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Mesenteric_lymph1 <- merge(Mesenteric_lymph1,df_id,by = "SYMBOL",all=F)
gene <- unique(Mesenteric_lymph1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","Mesenteric_lymph","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
PBMC1 <- PBMC1[order(PBMC1$avg_log2FC,decreasing = T),]
colnames(PBMC1)[6] <- "SYMBOL"
df_id <- bitr(PBMC1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
PBMC1 <- merge(PBMC1,df_id,by = "SYMBOL",all=F)
gene <- unique(PBMC1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","PBMC","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]
Spleen1 <- Spleen1[order(Spleen1$avg_log2FC,decreasing = T),]
colnames(Spleen1)[6] <- "SYMBOL"
df_id <- bitr(Spleen1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Spleen1 <- merge(Spleen1,df_id,by = "SYMBOL",all=F)
gene <- unique(Spleen1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","Spleen","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

# figureS1 D DE_union（上调-GO富集）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
DE_list <- c("union","Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO <- data.frame()
for (celltype in DE_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/",celltype,"_setdiff_GO_UP.tab"),stringsAsFactors = F)
  # GO_change <- GO1[c(1:5),]
  GO_change <- GO1
  GO_change$celltype <- rep(celltype)
  GO <- rbind(GO,GO_change)
  GO <- GO[GO$pvalue < 0.05,]
  write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","figure.S1D_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
}

# figureS1D DE_union（下调-GO富集）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Bone_marrow","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Mesenteric_lymph <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Mesenteric_lymph","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
Mesenteric_lymph <- Mesenteric_lymph[!grepl("ENSMMUG", Mesenteric_lymph$gene),]
PBMC <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","PBMC","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
PBMC <- PBMC[!grepl("ENSMMUG", PBMC$gene),]
Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","Spleen","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "DOWN",]
Spleen <- Spleen[!grepl("ENSMMUG", Spleen$gene),]
DE <- rbind(Bone_marrow,Mesenteric_lymph,PBMC,Spleen)
gene_DE <- data.frame(table(DE$gene));gene_DE <- gene_DE[gene_DE$Freq == "4",]
gene_union_Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Bone_marrow$tissue <- rep("Bone_marrow")
gene_union_Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Mesenteric_lymph$tissue <- rep("Mesenteric_lymph")
gene_union_PBMC <- PBMC[PBMC$gene %in% gene_DE$Var1,c(2,6)]
gene_union_PBMC$tissue <- rep("PBMC")
gene_union_Spleen <- Spleen[Spleen$gene %in% gene_DE$Var1,c(2,6)]
gene_union_Spleen$tissue <- rep("Spleen")
gene_union <- rbind(gene_union_Bone_marrow,gene_union_Mesenteric_lymph,gene_union_PBMC,gene_union_Spleen)
colnames(gene_union)[2] <- "SYMBOL"
df_id <- bitr(gene_union$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
gene_union <- merge(gene_union,df_id,by = "SYMBOL",all=F)
gene <- unique(gene_union[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","union","_setdiff_GO_DOWN.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
Bone_marrow1 <- Bone_marrow1[order(Bone_marrow1$avg_log2FC,decreasing = T),]
colnames(Bone_marrow1)[6] <- "SYMBOL"
df_id <- bitr(Bone_marrow1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Bone_marrow1 <- merge(Bone_marrow1,df_id,by = "SYMBOL",all=F)
gene <- unique(Bone_marrow1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","Bone_marrow","_setdiff_GO_DOWN.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
Mesenteric_lymph1 <- Mesenteric_lymph1[order(Mesenteric_lymph1$avg_log2FC,decreasing = T),]
colnames(Mesenteric_lymph1)[6] <- "SYMBOL"
df_id <- bitr(Mesenteric_lymph1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Mesenteric_lymph1 <- merge(Mesenteric_lymph1,df_id,by = "SYMBOL",all=F)
gene <- unique(Mesenteric_lymph1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","Mesenteric_lymph","_setdiff_GO_DOWN.tab"), quote = FALSE,sep="\t",row.names = FALSE)

PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
PBMC1 <- PBMC1[order(PBMC1$avg_log2FC,decreasing = T),]
colnames(PBMC1)[6] <- "SYMBOL"
df_id <- bitr(PBMC1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
PBMC1 <- merge(PBMC1,df_id,by = "SYMBOL",all=F)
gene <- unique(PBMC1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","PBMC","_setdiff_GO_DOWN.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]
Spleen1 <- Spleen1[order(Spleen1$avg_log2FC,decreasing = T),]
colnames(Spleen1)[6] <- "SYMBOL"
df_id <- bitr(Spleen1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Spleen1 <- merge(Spleen1,df_id,by = "SYMBOL",all=F)
gene <- unique(Spleen1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/","Spleen","_setdiff_GO_DOWN.tab"), quote = FALSE,sep="\t",row.names = FALSE)

# figureS1 E DE_union（下调-GO富集）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
DE_list <- c("union","Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO <- data.frame()
for (celltype in DE_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/GO/",celltype,"_setdiff_GO_DOWN.tab"),stringsAsFactors = F)
  # GO_change <- GO1[c(1:5),]
  GO_change <- GO1
  GO_change$celltype <- rep(celltype)
  GO <- rbind(GO,GO_change)
  GO <- GO[GO$pvalue < 0.05,]
  write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS1/","figure.S1E_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
}
