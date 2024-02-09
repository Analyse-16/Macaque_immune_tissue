setwd("/data1/zhuzn/wangsn/reviewer_4000/")
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

#### Monocyte 单核亚型-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "integrated"
unique(all_sample.combined@meta.data$celltype)
Cell.combined<-subset(all_sample.combined, subset = celltype  %in% "Monocyte")
set.seed(1)
Cell.combined <- RunPCA(Cell.combined, npcs = 30, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 42, reduction = "pca", dims = 1:20)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:20)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.2,0.4,0.8,1))
p1 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE)
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE)
p3 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.8", label = TRUE)
p4 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.1", label = TRUE)
p <- p<-(p1|p2)/(p3|p4)
p
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_cluster.RData")

rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Monocyte/Cell.combined.celltype.RData")
Cell.combined@meta.data$Tissue <-  gsub("Lymph_node","Mesenteric_lymph",Cell.combined@meta.data$Tissue)
DefaultAssay(Cell.combined) <- "RNA"
colnames(Cell.combined@meta.data)
cell <- Cell.combined@meta.data[,c(4,5,18)]
cell %<>% {.$id<-rownames(.);.}
cell$id <-gsub("Lymph_node","Mesenteric_lymph",cell$id)
cell$id <-gsub("-",".",cell$id)

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
colnames(Cell.combined@meta.data)
cell1 <- Cell.combined@meta.data[,c(4,5)]
cell1 %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(cell1) <- c("Tissue","Age","id","Celltype_old")
cell1[is.na(cell1)] <- "NON"
cell2 <- cell1[cell1$Celltype_old == "NON",]
cell2 <- separate(cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined <- AddMetaData(Cell.combined,cell1$Celltype_old,col.name = "Celltype_old")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "Celltype_old", label = TRUE,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_cluster.RData")

rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
gene <- c("CD14","FCGR3")
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.0.4)

p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.0.4',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black")) + RotatedAxis()
p

rm(list=ls())
celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_celltype.txt")
unique(celltype$celltype)
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_cluster.RData")
colnames(Cell.combined@meta.data)
table(Cell.combined@meta.data$integrated_snn_res.0.4)
cluster <- Cell.combined@meta.data[,c(1,10)]
colnames(cluster) <- c("orig.ident","cluster")
library(plyr)
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "celltype")
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_celltype.RData")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype", label = TRUE,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p

#### Macrophage 巨噬亚型-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "integrated"
unique(all_sample.combined@meta.data$celltype)
Cell.combined<-subset(all_sample.combined, subset = celltype  %in% "Macrophage")
set.seed(1)
Cell.combined <- RunPCA(Cell.combined, npcs = 30, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 42, reduction = "pca", dims = 1:20)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:20)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.2,0.4,0.8,1))
p1 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE)
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE)
p3 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.8", label = TRUE)
p4 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.1", label = TRUE)
p <- p<-(p1|p2)/(p3|p4)
p
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_cluster.RData")

rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Macrophage/Cell.combined.celltype.RData")
Cell.combined@meta.data$Tissue <-  gsub("Lymph_node","Mesenteric_lymph",Cell.combined@meta.data$Tissue)
DefaultAssay(Cell.combined) <- "RNA"
colnames(Cell.combined@meta.data)
cell <- Cell.combined@meta.data[,c(4,5,18)]
cell %<>% {.$id<-rownames(.);.}
cell$id <-gsub("Lymph_node","Mesenteric_lymph",cell$id)
cell$id <-gsub("-",".",cell$id)

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
colnames(Cell.combined@meta.data)
cell1 <- Cell.combined@meta.data[,c(4,5)]
cell1 %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(cell1) <- c("Tissue","Age","id","Celltype_old")
cell1[is.na(cell1)] <- "NON"
cell2 <- cell1[cell1$Celltype_old == "NON",]
cell2 <- separate(cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined <- AddMetaData(Cell.combined,cell1$Celltype_old,col.name = "Celltype_old")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "Celltype_old", label = TRUE,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_cluster.RData")

rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
gene <- c("CXCL10","TLR2","TLR4","CD80","CD86","CD163","C1QB","C1QC","MRC1","H4C3","PCNA","MKI67")
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.0.4)

p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.0.4',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black")) + RotatedAxis()
p

rm(list=ls())
celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_celltype.txt")
unique(celltype$celltype)
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_cluster.RData")
colnames(Cell.combined@meta.data)
table(Cell.combined@meta.data$integrated_snn_res.0.4)
cluster <- Cell.combined@meta.data[,c(1,10)]
colnames(cluster) <- c("orig.ident","cluster")
library(plyr)
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "celltype")
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_celltype.RData")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype", label = TRUE,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p

# figure6A  髓系细胞(单核、巨噬、巨核、树突、粒细胞、HSC) ############# --------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V8/all_sample.combined.celltype.RData")
all_sample.combined2 <- subset(all_sample.combined, subset = celltype  %in% c("Megakaryocyte"))
all_sample.combined2@meta.data$id <- colnames(all_sample.combined2)
all_sample.combined2$id <-gsub("-",".",all_sample.combined2$id)
all_sample.combined2$id <-gsub("Lymph_node","Mesenteric_lymph",all_sample.combined2$id)
table(all_sample.combined2$celltype)
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined@meta.data$id <- colnames(all_sample.combined)
all_sample.combined$celltype <- as.character(all_sample.combined$celltype)
table(all_sample.combined$celltype)
all_sample.combined$Celltype <- ifelse(all_sample.combined$id %in% all_sample.combined2$id,all_sample.combined2$celltype,all_sample.combined$celltype)
table(all_sample.combined$Celltype)
all_sample.combined$celltype <- all_sample.combined$Celltype
table(all_sample.combined$celltype)

colnames(all_sample.combined@meta.data)
unique(all_sample.combined@meta.data$celltype)
Cell.combined1 <- subset(all_sample.combined, subset = celltype  %in% c("Monocyte","Macrophage","Megakaryocyte","Dendritic cell","Granulocyte","HSC"))
all_cell <- Cell.combined1@meta.data[,c(4,5,18)]

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_celltype.RData")
colnames(Cell.combined@meta.data)
cell_Mon <- Cell.combined@meta.data[,c(4,5,18)]
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_celltype.RData")
colnames(Cell.combined@meta.data)
cell_Mac <- Cell.combined@meta.data[,c(4,5,18)]
cell <- rbind(cell_Mon,cell_Mac)
cell  %<>% {.$id<-rownames(.);.} 

all_cell %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell) <- c("Tissue","Age","celltype","id","Celltype")
all_cell[is.na(all_cell)] <- "NON"
all_cell$celltype <- as.character(all_cell$celltype)
table(all_cell$Celltype)
all_cell$Celltype <- ifelse(all_cell$Celltype == "NON" ,all_cell$celltype,all_cell$Celltype)
unique(all_cell$Celltype)
Cell.combined1 <- AddMetaData(Cell.combined1,all_cell$Celltype,col.name = "sub_celltype")
table(Cell.combined1@meta.data$celltype)
table(Cell.combined1@meta.data$sub_celltype)
save(Cell.combined1, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")

### 重新聚类
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
Cell.combined <- Cell.combined1
ncol(Cell.combined)
DefaultAssay(Cell.combined) <- "integrated"
Cell.combined <- RunPCA(Cell.combined, npcs = 20, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 1, reduction = "pca", dims = 1:10)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:10)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.2,0.4,0.8,1))
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
Cell.combined$sub_celltype <- factor(Cell.combined$sub_celltype,levels = c("CD14_Mon","In-termed_Mon","CD16_Mon","M1",
                                                                             "M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC"))
colors <- c("#2092E2","#008A8A","#7756A1","#FF6666","#74B404","#B3DE69","#DBBAFA","#FFD95A","#FFAF18")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "sub_celltype",cols = colors,label = F)
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6A.pdf",plot = p,width = 10,height = 8)

Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
Cell.combined@meta.data$Tissue_Age <- factor(Cell.combined@meta.data$Tissue_Age,levels = c("Bone_marrow_young","Mesenteric_lymph_young","PBMC_young","Spleen_young",
                                                                                               "Bone_marrow_old","Mesenteric_lymph_old","PBMC_old","Spleen_old"))

p <- DimPlot(Cell.combined, reduction = "umap", group.by = "sub_celltype",
             split.by = "Tissue_Age",ncol = 4,cols = colors,label = F)+
  theme_bw()+
  labs(title = "")+ 
  theme(strip.text.x = element_blank()) +
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),panel.border = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/","figure.S6B",".pdf"),p,width = 12,height = 6)

# figure6B FeaturePlot-------------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
DefaultAssay(Cell.combined) <- "RNA"
gene <- c("CD14","FCGR3","CD86","CD163","PPBP","CD83","ITGAM","MCL1")
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$sub_celltype)
unique(Cell.combined$sub_celltype)
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'sub_celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)

data$Cluster <- factor(data$Cluster,levels = rev(c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC")))
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  # scale_color_gradient2(low="#e0e0e0",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  scale_color_gradientn(colors = c("grey","#eff3ff","#F76C6C","#981216"))+
  theme_classic()+labs(x = '', y = '')+
  theme(legend.position = "top")+
  theme(axis.text=element_text(size=13, color="black")) + RotatedAxis()
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6B.pdf",p,width = 5,height = 4)

# figure6C FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
colnames(Cell.combined@meta.data)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,21)]
celltype_tissue$sub_celltype <- as.character(celltype_tissue$sub_celltype)
unique(celltype_tissue$sub_celltype)
cell_list <-  c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$sub_celltype == j),]))
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               FC = as.character(cell[which(cell$Age == "old"),"Freq"]/cell[which(cell$Age == "young"),"Freq"]))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = "0",
                               FC = as.character(cell[which(cell$Age == "old"),"Freq"]))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = "0",
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               FC = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }
  }
  # freq_data_tissue$old_per <- freq_data_tissue$old_cell/sum(as.numeric(freq_data_tissue$old_cell))
  # freq_data_tissue$young_per <- freq_data_tissue$young_cell/sum(as.numeric(freq_data_tissue$young_cell))
  freq_data <- rbind(freq_data,freq_data_tissue)
}
freq_data$log2FC <- log2(as.numeric(freq_data$FC))

###fisher检验
fisher_data <- data.frame()
for (tissue in unique(celltype_tissue$Tissue)){
  test_data <- freq_data[freq_data$Tissue == tissue,]
  for (sub_celltype in cell_list){
    test_data1 <- test_data[test_data$sub_celltype == sub_celltype,c(3,4)]
    test_data2 <- test_data[test_data$sub_celltype != sub_celltype,]
    test_data3 <- data.frame(old_cell = sum(as.numeric(test_data2$old_cell)),young_cell = sum(as.numeric(test_data2$young_cell)))
    test_data4 <- rbind(test_data1,test_data3)
    test_data4$old_cell <- as.numeric(test_data4$old_cell);test_data4$young_cell <- as.numeric(test_data4$young_cell)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,sub_celltype = sub_celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))
library("ggplot2")
unique(freq_data$sub_celltype)
freq_data$sub_celltype <- factor(freq_data$sub_celltype,levels = c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#2092E2","#008A8A","#7756A1","#FF6666","#74B404","#B3DE69","#DBBAFA","#FFD95A","#FFAF18")
show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
freq_data$log2FC[!is.finite(freq_data$log2FC)] <- 0
p <- ggplot(data=freq_data, aes(x=sub_celltype, y=log2FC, fill = sub_celltype, width=0.8))+
  ylim(-2, 6)+
  facet_wrap(.~Tissue,scales= "free",nrow=2)+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+ 
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=sub_celltype, y=log2FC,label=fisher_data$p.value1),size=6,position= position_dodge(0.6))
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6C.pdf",plot = p1,width = 7,height = 7)
dev.off()

# figure6C  占比饼图 ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
colnames(Cell.combined@meta.data)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,21)]
celltype_tissue$sub_celltype <- as.character(celltype_tissue$sub_celltype)
unique(celltype_tissue$sub_celltype)
cell_list <-  c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  tissue1 <- tissue[tissue$Age == "young",]
  tissue2 <- tissue[tissue$Age == "old",]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$sub_celltype == j),]))
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(nrow(cell) == 0){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell_per = "0",
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell_per = "0",
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }
  }
  freq_data <- rbind(freq_data,freq_data_tissue)
}
str(freq_data)
freq_data$old_cell_per <- as.numeric(freq_data$old_cell_per);freq_data$young_cell_per <- as.numeric(freq_data$young_cell_per)
freq_data1 <- melt(freq_data,id.vars = c("Tissue","sub_celltype"))
colnames(freq_data1) <- c("Tissue","sub_celltype","Age","cell_per")
freq_data1$Age <- gsub("young_cell_per","Young",freq_data1$Age)
freq_data1$Age <- gsub("old_cell_per","Aged",freq_data1$Age)
freq_data2 <- spread(freq_data1,sub_celltype,cell_per)
freq_data2$x=c(4,2,8,6,12,10,16,14)
freq_data2$y=rep(2)
freq_data2$region=as.character(1:8)
colnames(freq_data2)
freq_data3 <- freq_data2[,c(1,2,12,13,3,8,4,9,10,11,5,6,7,14)]
colors <- c("#2092E2","#008A8A","#7756A1","#FF6666","#74B404","#B3DE69","#DBBAFA","#FFD95A","#FFAF18")
p <- ggplot()+
  geom_scatterpie(data=freq_data2,
                  aes(x,y,group=region,r=0.9),
                  cols = cell_list)+
  coord_equal()+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values = colors)
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6C1.pdf",plot = p,width = 15,height = 7)
dev.off()

# figure6C percent_FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
colnames(Cell.combined@meta.data)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,19)]
celltype_tissue$sub_celltype <- as.character(celltype_tissue$sub_celltype)
unique(celltype_tissue$sub_celltype)
cell_list <-  c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  tissue1 <- tissue[tissue$Age == "young",]
  tissue2 <- tissue[tissue$Age == "old",]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$sub_celltype == j),]))
    
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(nrow(cell) == 0){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = "0",
                               young_cell = "0",
                               old_cell_per = "0",
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = "0",
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = "0",
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               old_cell_per = "0",
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }
  }
  freq_data <- rbind(freq_data,freq_data_tissue)
}
str(freq_data)
freq_data$old_cell_per <- as.numeric(freq_data$old_cell_per);freq_data$young_cell_per <- as.numeric(freq_data$young_cell_per)
freq_data$log2FC <- freq_data$old_cell_per-freq_data$young_cell_per

###fisher检验
fisher_data <- data.frame()
for (tissue in unique(celltype_tissue$Tissue)){
  test_data <- freq_data[freq_data$Tissue == tissue,]
  for (sub_celltype in cell_list){
    test_data1 <- test_data[test_data$sub_celltype == sub_celltype,c(3,4)]
    test_data2 <- test_data[test_data$sub_celltype != sub_celltype,]
    test_data3 <- data.frame(old_cell = sum(as.numeric(test_data2$old_cell)),young_cell = sum(as.numeric(test_data2$young_cell)))
    test_data4 <- rbind(test_data1,test_data3)
    test_data4$old_cell <- as.numeric(test_data4$old_cell);test_data4$young_cell <- as.numeric(test_data4$young_cell)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,sub_celltype = sub_celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))
library("ggplot2")
unique(freq_data$sub_celltype)
freq_data$sub_celltype <- factor(freq_data$sub_celltype,levels = c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#2092E2","#008A8A","#7756A1","#FF6666","#74B404","#B3DE69","#DBBAFA","#FFD95A","#FFAF18")
show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
freq_data$log2FC[!is.finite(freq_data$log2FC)] <- 0
p <- ggplot(data=freq_data, aes(x=sub_celltype, y=log2FC, fill = sub_celltype, width=0.8))+
  facet_wrap(.~Tissue,scales= "free",nrow=2)+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  scale_y_continuous(limits=c(-0.1,0.1))+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+ 
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,colour = "black"),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=sub_celltype, y=log2FC,label=fisher_data$p.value1),size=6,position= position_dodge(0.6))
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6C.pdf",plot = p1,width = 7,height = 7)
dev.off()

# figure6D 髓系细胞DEG (组织-细胞类型-old-young) -----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Age)
data.list_tissue <- SplitObject(Cell.combined, split.by = "Tissue")
tissue_list <- unique(Cell.combined$Tissue)
for (tissue in tissue_list){
  tissue = "Mesenteric_lymph"
  Cell.combined2 <- data.list_tissue[[tissue]]
  table(Cell.combined2$sub_celltype)
  data.list_celltype <- SplitObject(Cell.combined2, split.by = "sub_celltype")
  for (celltype in unique(Cell.combined2$sub_celltype)){
    celltype = "HSC"
    Cell.combined3 <- data.list_celltype[[celltype]]
    d <- data.frame(table(Cell.combined3$Age))
    if(nrow(d) == 1){
      print(paste0(tissue,"_",celltype," no young or old cells"))
    }
    else{
      a = d[1,2];b = d[2,2]
      if(a<3 | b <3){
        print(paste0(tissue,"_",celltype,"Cells fewer than 3 cells")) 
      }else{
        diff_gene <- FindMarkers(Cell.combined3,ident.1 = "old", ident.2 = "young",min.cells.group = 0,logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
        diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
        diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                               ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                      ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                               'NOT'))
        diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
        celltype <- gsub(" ","_",celltype)
        write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE/",tissue,"_",celltype,"_DE.txt"),quote = FALSE,sep = "\t",row.names = F) 
      } 
    }
  }
}

# figure6D inflammatory.gene FC （组织细胞类型）----------------------------------------------------
rm(list=ls())
DE_list <- list.files("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE/")
DE_list <- gsub("_DE.txt","",DE_list)
DE_list1 <- c("Bone_marrow_CD14_Mon","Bone_marrow_In-termed_Mon","Bone_marrow_CD16_Mon","Bone_marrow_M1","Bone_marrow_M2",
             "Bone_marrow_Megakaryocyte","Bone_marrow_Dendritic_cell","Bone_marrow_Granulocyte","Bone_marrow_HSC",
             "Mesenteric_lymph_CD14_Mon","Mesenteric_lymph_In-termed_Mon","Mesenteric_lymph_CD16_Mon",
             "Mesenteric_lymph_M1","Mesenteric_lymph_M2",
             "Mesenteric_lymph_Megakaryocyte","Mesenteric_lymph_Dendritic_cell","Mesenteric_lymph_Granulocyte","Mesenteric_lymph_HSC",
             "PBMC_CD14_Mon","PBMC_In-termed_Mon","PBMC_CD16_Mon","PBMC_M1","PBMC_M2","PBMC_Megakaryocyte","PBMC_Dendritic_cell","PBMC_Granulocyte","PBMC_HSC",
             "Spleen_CD14_Mon","Spleen_In-termed_Mon","Spleen_CD16_Mon","Spleen_M1","Spleen_M2","Spleen_Megakaryocyte",
             "Spleen_Dendritic_cell","Spleen_Granulocyte","Spleen_HSC")

inflammatory_gene <- read.csv("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/inflammatory_gene1.txt", sep="")
# inflammatory_gene <- read.csv("/data1/zhuzn/wangsn/SASP/myeloid_SASP/chemokine_gene.txt", sep="")
# inflammatory_gene <- read.csv("/data1/zhuzn/wangsn/SASP/myeloid_SASP/cytokines_+_gene.txt", sep="")
# inflammatory_gene <- read.csv("/data1/zhuzn/wangsn/SASP/myeloid_SASP/cytokines_-_gene.txt", sep="")
# inflammatory_gene <- read.csv("/data1/zhuzn/wangsn/SASP/myeloid_SASP/antigen presentation_gene.txt", sep="")

genes <- unique(inflammatory_gene$gene)
# genes <- c("CCL1","CCL11","CCL14","CCL2","CCL20","CCL24","CCL3","CCL4","CCL5",
#           "CXCL1","CXCL10","CXCL12","CXCL16","CXCL2","CXCL3","CXCL4","CXCL8","CXCL9","XCL1","XCL2")
diff_gene_overlap_FC <- data.frame()
for (DE in DE_list){
  DE_data <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE/",DE,"_DE.txt"))
  DE_data1 <- DE_data[DE_data$gene %in% genes,]
  if(nrow(DE_data1) < 1){
    print(paste0("no genes"))
  }
  else{
    DE_data1 <- DE_data1[,c(6,2,1)];DE_data1$celltype <- rep(DE)
    
    diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,DE_data1)
  }
}
diff_gene_overlap_FC <- diff_gene_overlap_FC[!grepl("ENSMMUG", diff_gene_overlap_FC$gene),]
# write.table(diff_gene_overlap_FC,file = paste0("./inflammatory_gene_FC.txt"),quote = FALSE,sep = "\t",row.names = F)

overlap_FC <- diff_gene_overlap_FC
str(overlap_FC)
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$p_val <- as.numeric(overlap_FC$p_val)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(genes))
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = DE_list)
overlap_FC$group <- rep("inflammatory")
overlap_FC$group <- ifelse(overlap_FC$gene %in% c("CCL3","CCL5","CCL2") ,"chemokine",overlap_FC$group)
overlap_FC$group <- ifelse(overlap_FC$gene %in% c("TNF","IFNG") ,"cytokines",overlap_FC$group)
overlap_FC$group <- ifelse(overlap_FC$gene %in% c("PSMB8","TAPBP") ,"antigen presentation",overlap_FC$group)
overlap_FC$group <- factor(overlap_FC$group,levels = c("inflammatory","chemokine","cytokines","antigen presentation"))
colors <- c("#1A58A7","#4592E5","#6baed6","#c6dbef","#eff3ff","#fcbba1","#F76C6C","#F70000","#981216")
colors <- c("#1A58A7","#4592E5","#6baed6","#89B9D6","#A4C4D6","#c6dbef","#eff3ff","#F2F5FF","#FCC8B3","#fcbba1","#F76C6C","#F70000","#981216")
show_col(colors)
overlap_FC <- overlap_FC %>%
  mutate(text = case_when(
    p_val > 0.05 ~ "", 
    p_val <= 0.05 & p_val >= 0.01 ~ "*",
    p_val <= 0.01 ~ "**"))
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = DE_list1)
overlap_FC$avg_log2FC[overlap_FC$avg_log2FC > 2] = 2
overlap_FC$avg_log2FC[overlap_FC$avg_log2FC < -2] = -2

p <- ggplot(overlap_FC,aes(x=celltype,y=gene,fill=avg_log2FC))+
  geom_raster()+
  facet_grid(group~.,scales = "free",space = "free")+
  scale_fill_gradientn(colors = colors)+
  # geom_text(aes(label=text),col ="black",size = 5) +
  theme_bw() + 
  labs(title = "Myeloid_Cell_FC")+ 
  theme(strip.text.y = element_text(size = 15, angle = 0)) +
  theme(legend.title = element_text(size = 20),legend.text = element_text(size = 20),
        panel.grid = element_blank(),plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle = 90,vjust = 0.5,hjust = 1),axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),panel.border = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6D.pdf",plot = p,width = 10,height = 7.5)

# figure6E 脾脏巨核细胞diff（小提琴图）-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V8/all_sample.combined.celltype.RData")
all_sample.combined2 <- subset(all_sample.combined, subset = celltype  %in% c("Megakaryocyte"))
all_sample.combined2@meta.data$id <- colnames(all_sample.combined2)
all_sample.combined2$id <-gsub("-",".",all_sample.combined2$id)
all_sample.combined2$id <-gsub("Lymph_node","Mesenteric_lymph",all_sample.combined2$id)
table(all_sample.combined2$celltype)
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined@meta.data$id <- colnames(all_sample.combined)
all_sample.combined$celltype <- as.character(all_sample.combined$celltype)
all_sample.combined$Celltype <- ifelse(all_sample.combined$id %in% all_sample.combined2$id,all_sample.combined2$celltype,all_sample.combined$celltype)
table(all_sample.combined$Celltype)
all_sample.combined$celltype <- all_sample.combined$Celltype

data.list_tissue <- SplitObject(all_sample.combined, split.by = "Tissue")
Cell.combined <- data.list_tissue[["Spleen"]]
data.list_celltype <- SplitObject(Cell.combined, split.by = "celltype")
Cell.combined <- data.list_celltype[["Megakaryocyte"]]
colnames(Cell.combined@meta.data)
Cell.combined@meta.data[which(Cell.combined@meta.data$Age == "young"),5] <- "Young"
Cell.combined@meta.data[which(Cell.combined@meta.data$Age == "old"),5] <- "Aged"
Cell.combined$Age <- factor(Cell.combined$Age,levels = c("Young","Aged"))
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Age)
diff_gene <- FindMarkers(Cell.combined,ident.1 = "Aged", ident.2 = "Young",logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE/Spleen_Megakaryocyte_DE.txt"),quote = FALSE,sep = "\t",row.names = F) 

DE <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE/Spleen_Megakaryocyte_DE.txt") %>% .[.$change_V1 == "UP",]
gene <- c("COX6C","CCL5","KLRD1","MIF","SAMHD1","TNFRSF1B")
str(Cell.combined@meta.data)
p1 <- VlnPlot(Cell.combined, cols = c('#39A8E3',"#D75156"),features = gene, pt.size = 0,ncol = 3) + NoLegend()
p1
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6E.pdf"),plot = p1,width = 5,height = 6)
dev.off()

# figure6F 脾脏巨核细胞diff（富集）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
DE <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE/Spleen_Megakaryocyte_DE.txt") %>% .[.$change_V1 == "UP",]
colnames(DE)[6] = "SYMBOL"
DE <- merge(x=DE,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_result<-as.data.frame(ego_BP@result)
GO <- c("purine deoxyribonucleotide catabolic process","deoxyribonucleoside triphosphate catabolic process",
        "positive regulation by host of viral process","glyoxylate metabolic process",
        "ATP metabolic process")
GO  <- ego_BP_result [ego_BP_result$Description %in% GO,]
library(DOSE)
GO$GeneRatio1 <- parse_ratio(GO$GeneRatio)
p <- ggplot(GO,aes(reorder(Description,GeneRatio1),GeneRatio1,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_colour_gradient(low="#810000",high="#F56F6F")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
  geom_point(size = 2.0,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6F.pdf"),plot = p,width = 5,height = 3.5)
dev.off()

# figure6 髓系CellChat-------------------------------------------------------------------------
rm(list=ls())
library(CellChat)
library(patchwork)
library(tidyverse)
library(ggalluvial)
load(file = "/data1/zhuzn/wangsn/V8/all_sample.combined.celltype.RData")
all_sample.combined2 <- subset(all_sample.combined, subset = celltype  %in% c("Fibroblast","Megakaryocyte"))
all_sample.combined2@meta.data$id <- colnames(all_sample.combined2)
all_sample.combined2$id <-gsub("-",".",all_sample.combined2$id)
all_sample.combined2$id <-gsub("Lymph_node","Mesenteric_lymph",all_sample.combined2$id)
table(all_sample.combined2$celltype)
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined@meta.data$id <- colnames(all_sample.combined)
all_sample.combined$celltype <- as.character(all_sample.combined$celltype)
table(all_sample.combined$celltype)
all_sample.combined$Celltype <- ifelse(all_sample.combined$id %in% all_sample.combined2$id,all_sample.combined2$celltype,all_sample.combined$celltype)
table(all_sample.combined$Celltype)
all_sample.combined$celltype <- all_sample.combined$Celltype
table(all_sample.combined$celltype)

all_sample.combined@meta.data$celltype <- factor(all_sample.combined@meta.data$celltype,
                                                  levels = c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                             "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"))
unique(all_sample.combined$celltype)
data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
data_age <- data.list_age[["old"]]
# data_age <- data.list_age[["young"]]

data.list_tissue <- SplitObject(data_age, split.by = "Tissue")
Cell.combined <- data.list_tissue[["Spleen"]]
unique(Cell.combined$celltype)
scRNA <- Cell.combined
##提取表达矩阵和细胞分类信息
# CellChat要求输入标准化后的表达数据
data.input = scRNA@assays[["RNA"]]@data
meta = subset(scRNA@meta.data, select = "celltype")
##创建cellchat对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "celltype") 
groupSize <- as.numeric(table(cellchat@idents)) # 后面有用
##设置参考数据库
# 选择合适的物种，可选CellChatDB.human, CellChatDB.mouse
CellChatDB <- CellChatDB.human  
# 使用"Secreted Signaling"用于细胞通讯分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# 将数据库传递给cellchat对象
cellchat@DB <- CellChatDB.use 
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.mouse)
##推测细胞通讯网络
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cells <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
mat <- cellchat@net$weight
i = 9
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
mat2 <- mat2[cells,cells]
colnames(mat)
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
show_col(colors)
p <- netVisual_circle(mat2, color.use = colors,vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat2)[i])
print(p)
dev.off()

cellchat@netP$pathways
pathways.show <- c("MIF") 
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# netVisual_bubble(cellchat, signaling = c("CCL"), remove.isolate = FALSE)
p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.use = colors,color.heatmap = "Reds")
print(p)
dev.off()
saveRDS(cellchat, file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/cellchat_Myeloid_old.rds")
# saveRDS(cellchat, file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/cellchat_Myeloid_young.rds")
dev.off()

# figure6GH CellChat结果-------------------------------------------------------------------------
rm(list=ls())
library(CellChat)
library(patchwork)
library(tidyverse)
library(ggalluvial)
cellchat1 <- readRDS("/data1/zhuzn/wangsn/reviewer_4000/figure6/cellchat_Myeloid_old.rds")
cellchat2 <- readRDS("/data1/zhuzn/wangsn/reviewer_4000/figure6/cellchat_Myeloid_young.rds")

groupSize <- as.numeric(table(cellchat1@idents))
cells <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
cellchat1@net$weight1 <- cellchat1@net$weight-cellchat2@net$weight
mat <- cellchat1@net$weight1
i = 9
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
colnames(mat)
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
show_col(colors)
pdf(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6G.pdf"), width = 5, height = 5)
p <- netVisual_circle(mat2, color.use = colors,vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat2)[i])
print(p)
dev.off()

cellchat1@netP$pathways
cellchat2@netP$pathways
cellchat1@netP$prob <- cellchat1@netP$prob[,,c(1,2)]-cellchat2@netP$prob[,,c(1,2)]
cellchat1@netP$pathways <- c("MIF","RESISTIN")
pathways.show <- c("MIF") 
# Heatmap
par(mfrow=c(1,1))
cellchat1@netP$prob
netVisual_heatmap(cellchat1, signaling = pathways.show)
pdf(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure6H.pdf"), width = 4, height = 3)
cells <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
p <- netVisual_heatmap(cellchat1, signaling = pathways.show, color.use = colors)
print(p)
dev.off()

# 伪时间-------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
# monocle3 
library(ggplot2)
library(Seurat)
library(scater)
library(Matrix)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(monocle3)
library(ggplot2)
library(dplyr)
data.list_age <- SplitObject(Cell.combined, split.by = "Age")
young_data <- data.list_age[["young"]]
cell_data <- Cell.combined

unique(cell_data@meta.data$sub_celltype)
# DefaultAssay(cell_data) <- "integrated"
# data <- as.matrix(cell_data@assays[["integrated"]]@scale.data)

DefaultAssay(cell_data) <- "RNA"
data <- as(as.matrix(cell_data@assays$RNA@counts), 'sparseMatrix') 
# metadata
cell_metadata <- as.data.frame(cell_data@meta.data)
cell_metadata <- cell_metadata[colnames(data),]
# gene_annotation
gene_annotation <- data.frame(gene_short_name=rownames(data))
rownames(gene_annotation) <- gene_annotation$gene_short_name
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 20,norm_method = c("none"), verbose = T)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method ="UMAP",max_components = 2,cores = 1, verbose = T,umap.fast_sgd = FALSE)
cds <- cluster_cells(cds, reduction_method ="UMAP",
                     #resolution=0.01, 
                     verbose = T)
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
cds <- learn_graph(cds,verbose = T)

# cds = order_cells(cds)
get_earliest_principal_node <- function(cds, NSC_type="HSC"){
  cell_ids <- which(colData(cds)[, "sub_celltype"] == NSC_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
# save(cds,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/","myeloid_cds_inter",".RData"))
save(cds,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/","myeloid_cds_RNA",".RData"))

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/myeloid_cds_RNA.RData")
cds$Tissue_Age <- paste0(cds$Tissue,"_",cds$Age)
cds@colData@listData[["Tissue_Age"]] <- factor(cds@colData@listData[["Tissue_Age"]],levels = c("Bone_marrow_young","Bone_marrow_old",
                                                                                               "Mesenteric_lymph_young","Mesenteric_lymph_old",
                                                                                               "PBMC_young","PBMC_old",
                                                                                               "Spleen_young","Spleen_old"))
# cds@colData@listData[["Age"]] <- factor(cds@colData@listData[["Age"]],levels = c("young","old"))
unique(cds@colData@listData[["sub_celltype"]])
cds@colData@listData[["sub_celltype"]] <- factor(cds@colData@listData[["sub_celltype"]],levels = c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2",
                                                                                                   "Megakaryocyte","Dendritic cell","Granulocyte","HSC"))
colors <- c("#1f78b4","#33a02c","#ff7f00","#8366A8","black","#e7298a","#a65628","#6DB1B7","#ef3b2c")
show_col(colors)
p1<-plot_cells(cds, reduction_method ="UMAP",color_cells_by="pseudotime", group_label_size = 6,
               cell_size = 0.3,graph_label_size = 3,label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE)+ 
  labs(title="")+
  theme(plot.title=element_text(size=16,color="black",hjust = 0.5))
p1
# ggsave(paste0("/data1/zhuzn/wangsn/V7/result/","figure_4D_1",".pdf"),p1,width = 6,height = 5)

p2<-plot_cells(cds, reduction_method ="UMAP",color_cells_by="sub_celltype", group_label_size = 6,
               cell_size = 0.5,graph_label_size = 3,
               label_branch_points = F,label_roots = F, label_leaves = F,
               label_cell_groups=F,label_groups_by_cluster=FALSE)+
  scale_color_manual(values=colors)+
  labs(title="")+
  theme(plot.title=element_text(size=16,color="black",hjust = 0.5),legend.position="right")
p2
p4<-plot_cells(cds, reduction_method ="UMAP",color_cells_by="sub_celltype", group_label_size = 6,
               cell_size = 0.5,graph_label_size = 3,
               label_branch_points = F,label_roots = F, label_leaves = F,
               label_cell_groups=F,label_groups_by_cluster=FALSE)+
  scale_color_manual(values=colors)+
  labs(title="")+
  facet_wrap(~Tissue_Age,ncol = 4)+
  theme(plot.title=element_text(size=16,color="black",hjust = 0.5),legend.position="right")
p4
p <- (p1|p2)/p4
p


# 髓系HSC两支差异基因 ---------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")
table(Cell.combined$sub_celltype)
colnames(Cell.combined@meta.data)
Cell.combined$sub_celltype <- as.character(Cell.combined$sub_celltype)
Cell.combined@meta.data[which(Cell.combined$sub_celltype %in% c("Dendritic cell","In-termed_Mon","CD16_Mon","M1",
                                                                 "M2","CD14_Mon","Megakaryocyte")),19] <- "Others"
table(Cell.combined$sub_celltype)
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$sub_celltype)
data.list_tissue <- SplitObject(Cell.combined, split.by = "Tissue")
sample_list <- unique(Cell.combined@meta.data$Tissue)
for (tissue in sample_list){
  tissue_data <- data.list_tissue[[tissue]]
  d <- data.frame(table(tissue_data$Age))
  if(nrow(d) == 1){
    print(paste0(celltype," no young or old cells"))
  }
  else{
    a = d[1,2];b = d[2,2]
    if(a<3 | b <3){
      print(paste0(celltype,"Cells fewer than 3 cells")) 
    }else{
      diff_gene <- FindAllMarkers(tissue_data,logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
      diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
      diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                             ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                    ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                             'NOT'))
      diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
      write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/",tissue,"DE.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
    } 
  }
}

# 伪时间映射细胞类型 ---------------------------------------------------------------
rm(list=ls())
library(monocle3)
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/myeloid_cds_RNA.RData")
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Myeloid_Cell.RData")

Cell.combined <- AddMetaData(Cell.combined,pseudotime,col.name = "Pseudotime")
colnames(Cell.combined@meta.data)
mat1=as.data.frame(Cell.combined@meta.data[,c(4,5,21,26)])
mat2=Embeddings(Cell.combined,"umap")
mat3=merge(mat2,mat1,by="row.names")

p <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=Pseudotime),size = 0.01)+
  # scale_color_gradientn(colors = c("#140789","#38049A","#7701A8","#9B169F","#AB2494","#BC3587","#CD4A76",
  #                                  "#EB7655","#F79143","#F9983E","#FDAE32","#FBD324","#F0F921"))+
  scale_color_gradientn(colors = c("#140789","#38049A","#9B169F","#AB2494",
                                   "#F9983E","#FDAE32","#FBD324","#F0F921"))+
  theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),panel.border = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/figure.S6A.pdf",plot = p,dpi=1000,width = 5,height = 5)

# figure4F 伪时间基因表达-------------------------------------------
rm(list=ls())
# t <- c("PBMC","Spleen","Mesenteric_lymph","Bone_marrow")
# DE.findmark <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/Mesenteric_lymphDE.findmark.txt")

# Track_genes_sig <- c("LTB","MAMU-DRA")#Bone_marrow
# Track_genes_sig <- c("DEFA1B","CST3")#Spleen
# Track_genes_sig <- c("HOPX","CTSB")#PBMC
Track_genes_sig <- c("CMA1","VMO1")#Mesenteric_lymph

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/myeloid_cds_RNA.RData")
cells <- data.frame(cds@colData@listData)
cells1 <- cells[cells$Tissue == "Mesenteric_lymph",]
cds1 <- cds[, rownames(cells1)]
cds1@colData$sub_celltype <- factor(cds1@colData$sub_celltype,levels = c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Dendritic cell","Granulocyte","HSC"))
colors <- c("#2092E2","#008A8A","#7756A1","#FF6666","#74B404","#B3DE69","#DBBAFA","#FFD95A","#FFAF18")
show_col(colors)
#基因表达趋势图
p <-plot_genes_in_pseudotime(cds1[Track_genes_sig], color_cells_by="sub_celltype",min_expr=0.5, ncol = 2)+
  scale_color_manual(values=colors)+
  theme_bw()+
  theme(strip.text.x = element_blank()) +
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),panel.border = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/","figure.S6C","Mesenteric_lymph",".pdf"),p,width = 10,height = 4)

# HSC组织之间的差异基因、每个组织年老比年轻的差异基因 ------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
unique(all_sample.combined@meta.data$celltype)
cell_data<-subset(all_sample.combined, subset = celltype  %in% "HSC")
DefaultAssay(cell_data) <- "RNA"
cell_data <- SetIdent(object = cell_data, value = cell_data@meta.data$Tissue)
diff_gene <- FindAllMarkers(cell_data, logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                       'NOT'))
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","HSC","_tissue.txt"),quote = FALSE,sep = "\t",row.names = F)

cell_data <- SetIdent(object = cell_data, value = cell_data@meta.data$Age)
data.list_tissue <- SplitObject(cell_data, split.by = "Tissue")
sample_list <- unique(cell_data@meta.data$Tissue)
unique(cell_data$celltype)
for (tissue in sample_list){
  tissue_data <- data.list_tissue[[tissue]]
  d <- data.frame(table(tissue_data$Age))
  if(nrow(d) == 1){
    print(paste0(" no young or old cells"))
  }
  else{
    a = d[1,2];b = d[2,2]
    if(a<3 | b <3){
      print(paste0("Cells fewer than 3 cells")) 
    }else{
      diff_gene <- FindMarkers(tissue_data, ident.1 = "old", ident.2 = "young",min.cells.group = 0,logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
      diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
      diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                             ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                    ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                             'NOT'))
      write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/",tissue,"_old_vs_young.txt"),quote = FALSE,sep = "\t",row.names = F)
    } 
  }
}

# HSC骨髓、脾脏富集 --------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_old_vs_young.txt"), sep="") %>% .[.$change_V1 %in% c("UP","DOWN"),]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
colnames(Bone_marrow)[6] <- "SYMBOL"
df_id <- bitr(Bone_marrow$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Bone_marrow <- merge(Bone_marrow,df_id,by = "SYMBOL",all=F)
gene <- unique(Bone_marrow[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_GO.tab"), quote = FALSE,sep="\t",row.names = FALSE)


Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_old_vs_young.txt"), sep="") %>% .[.$change_V1 %in% c("UP","DOWN"),]
Spleen <- Spleen[!grepl("ENSMMUG", Spleen$gene),]
colnames(Spleen)[6] <- "SYMBOL"
df_id <- bitr(Spleen$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Spleen <- merge(Spleen,df_id,by = "SYMBOL",all=F)
gene <- unique(Spleen[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_GO.tab"), quote = FALSE,sep="\t",row.names = FALSE)

# HSC（上调-特异）热图 -------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
cell <- c("HSC")

tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/",cell,"_tissue.txt")) %>% .[.$change_V1 == "UP",]
DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
Bone_marrow <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
DE1 <- tissue_DE[tissue_DE$cluster == "Bone_marrow",]
Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% intersect(DE1$gene,Bone_marrow$gene),]
Mesenteric_lymph <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Mesenteric_lymph","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
DE1 <- tissue_DE[tissue_DE$cluster == "Mesenteric_lymph",]
Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% intersect(DE1$gene,Mesenteric_lymph$gene),]
PBMC <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","PBMC","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
DE1 <- tissue_DE[tissue_DE$cluster == "PBMC",]
PBMC <- PBMC[PBMC$gene %in% intersect(DE1$gene,PBMC$gene),]
Spleen <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
DE1 <- tissue_DE[tissue_DE$cluster == "Spleen",]
Spleen <- Spleen[Spleen$gene %in% intersect(DE1$gene,Spleen$gene),]

Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
Bone_marrow1 <- Bone_marrow1[order(Bone_marrow1$avg_log2FC,decreasing = T),]

Bone_marrow2 <- Bone_marrow1
colnames(Bone_marrow2)[6] <- "SYMBOL"
df_id <- bitr(Bone_marrow2$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Bone_marrow2 <- merge(Bone_marrow2,df_id,by = "SYMBOL",all=F)
gene <- unique(Bone_marrow2[,'ENTREZID'])
ego_BP1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP1, OrgDb = org.Mmu.eg.db);ego_BP_result1<-as.data.frame(ego_BP1@result)
ego_BP_result1<-ego_BP_result1[ego_BP_result1$pvalue<0.05,]
ego_BP_result1$Tissue <-rep("Bone_marrow")

Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
Mesenteric_lymph1 <- Mesenteric_lymph1[order(Mesenteric_lymph1$avg_log2FC,decreasing = T),]

Mesenteric_lymph2 <- Mesenteric_lymph1
colnames(Mesenteric_lymph2)[6] <- "SYMBOL"
df_id <- bitr(Mesenteric_lymph2$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Mesenteric_lymph2 <- merge(Mesenteric_lymph2,df_id,by = "SYMBOL",all=F)
gene <- unique(Mesenteric_lymph2[,'ENTREZID'])
ego_BP2 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP2 <- setReadable(ego_BP2, OrgDb = org.Mmu.eg.db);ego_BP_result2<-as.data.frame(ego_BP2@result)
ego_BP_result2<-ego_BP_result2[ego_BP_result2$pvalue<0.05,]
ego_BP_result2$Tissue <-rep("Mesenteric_lymph")

PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
PBMC1 <- PBMC1[order(PBMC1$avg_log2FC,decreasing = T),]

PBMC2 <- PBMC1
colnames(PBMC2)[6] <- "SYMBOL"
df_id <- bitr(PBMC2$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
PBMC2 <- merge(PBMC2,df_id,by = "SYMBOL",all=F)
gene <- unique(PBMC2[,'ENTREZID'])
ego_BP3 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP3 <- setReadable(ego_BP3, OrgDb = org.Mmu.eg.db);ego_BP_result3<-as.data.frame(ego_BP3@result)
ego_BP_result3<-ego_BP_result3[ego_BP_result3$pvalue<0.05,]
ego_BP_result3$Tissue <-rep("PBMC")

Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]
Spleen1 <- Spleen1[order(Spleen1$avg_log2FC,decreasing = T),]

Spleen2 <- Spleen1
colnames(Spleen2)[6] <- "SYMBOL"
df_id <- bitr(Spleen2$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Spleen2 <- merge(Spleen2,df_id,by = "SYMBOL",all=F)
gene <- unique(Spleen2[,'ENTREZID'])
ego_BP4 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP4 <- setReadable(ego_BP4, OrgDb = org.Mmu.eg.db);ego_BP_result4<-as.data.frame(ego_BP4@result)
ego_BP_result4<-ego_BP_result4[ego_BP_result4$pvalue<0.05,]
ego_BP_result4$Tissue <-rep("Spleen")
ego_BP_result <- rbind(ego_BP_result1,ego_BP_result2,ego_BP_result3,ego_BP_result4)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/",cell,"_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)


diff_gene <- c(Bone_marrow1$gene,Mesenteric_lymph1$gene,PBMC1$gene,Spleen1$gene)
diff_gene <- diff_gene[!grepl("ENSMMUG", diff_gene)]
diff_gene_overlap_FC <- data.frame()
for (i in DE_list){
  old_vs_young <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/",i,"_old_vs_young.txt"))
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
  diff_gene2 <- tissue_DE[tissue_DE$cluster == i,c(2,7)]
  diff_gene1 <- merge(diff_gene1,diff_gene2,all.x=TRUE,by = "gene")
  colnames(diff_gene1) <- c("gene","log2FC.old_young","log2FC.tissue")
  diff_gene1$tissue <- rep(i)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}
diff_gene_overlap_FC[is.na(diff_gene_overlap_FC)] <- "0"
overlap_FC <- diff_gene_overlap_FC
str(overlap_FC)
overlap_FC$log2FC.tissue <- as.numeric(overlap_FC$log2FC.tissue)
overlap_FC$tissue <- factor(overlap_FC$tissue,levels = DE_list)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
p<-ggplot(overlap_FC,aes(x = tissue,y=gene,size=log2FC.old_young,colour=log2FC.tissue,ylab=''))+
  geom_point()+
  # coord_flip() +
  scale_size_continuous(range=c(0,6))+
  scale_color_gradientn(colors = c("blue","#3636FF","#2EB5E9","#B3B3FF","#B05AD4","#C297D4","grey","#eff3ff","#F76C6C","#981216"))+
  # scale_color_gradient2(low = "blue",mid = "grey", high = "#981216",midpoint = 1)+
  theme_classic()+
  labs(x = '', y = '',title = cell)+
  theme(axis.text.x=element_text(size=15, color="black",angle = 45,vjust = 1,hjust = 1),legend.position = "right",
        axis.text.y=element_text(size=15, color="black"))# + RotatedAxis():45度
p

p <- ggplot(overlap_FC,aes(x=tissue,y=gene,fill=log2FC.old_young))+
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

# HSC DE_union（上调-共同-特异）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library(reshape2)
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
Tissue <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/HSC_tissue.txt", sep="") %>% .[.$change_V1 == "UP",]
Tissue$gene <- gsub("\\..*","",Tissue$gene)
Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_old_vs_young.txt"), sep="") %>% .[.$change_V1 == "UP",]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Mesenteric_lymph <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Mesenteric_lymph","_old_vs_young.txt"), sep="") %>% .[.$change_V1 == "UP",]
Mesenteric_lymph <- Mesenteric_lymph[!grepl("ENSMMUG", Mesenteric_lymph$gene),]
PBMC <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","PBMC","_old_vs_young.txt"), sep="") %>% .[.$change_V1 == "UP",]
PBMC <- PBMC[!grepl("ENSMMUG", PBMC$gene),]
Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_old_vs_young.txt"), sep="") %>% .[.$change_V1 == "UP",]
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
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/",j,"_old_vs_young.txt"), header=T)
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
# overlap_FC[overlap_FC > 1] <- 1
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
# ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","figureS1D.pdf"),plot = p,width = 10,height = 15)

# HSC 组织之间的差异基因富集-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library(reshape2)
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
Tissue <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/HSC_tissue.txt", sep="") %>% .[.$change_V1 == "UP",]
Tissue$gene <- gsub("\\..*","",Tissue$gene)
colnames(Tissue)[7] <- "SYMBOL"
df_id <- bitr(Tissue$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Tissue <- merge(Tissue,df_id,by = "SYMBOL",all=F)
gene <- unique(Tissue[Tissue$cluster == "Bone_marrow",'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_GO_UP.txt"), quote = FALSE,sep="\t",row.names = FALSE)
gene <- unique(Tissue[Tissue$cluster == "Spleen",'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_GO_UP.txt"), quote = FALSE,sep="\t",row.names = FALSE)

Tissue <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/HSC_tissue.txt", sep="") %>% .[.$change_V1 == "DOWN",]
Tissue$gene <- gsub("\\..*","",Tissue$gene)
colnames(Tissue)[7] <- "SYMBOL"
df_id <- bitr(Tissue$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Tissue <- merge(Tissue,df_id,by = "SYMBOL",all=F)
gene <- unique(Tissue[Tissue$cluster == "Bone_marrow",'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_GO_DOWN.txt"), quote = FALSE,sep="\t",row.names = FALSE)
gene <- unique(Tissue[Tissue$cluster == "Spleen",'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_GO_DOWN.txt"), quote = FALSE,sep="\t",row.names = FALSE)

Tissue <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/HSC_tissue.txt", sep="")
Tissue$gene <- gsub("\\..*","",Tissue$gene)
colnames(Tissue)[7] <- "SYMBOL"
df_id <- bitr(Tissue$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Tissue <- merge(Tissue,df_id,by = "SYMBOL",all=F)
gene <- unique(Tissue[Tissue$cluster == "Bone_marrow",'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
gene <- unique(Tissue[Tissue$cluster == "Spleen",'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# HSC 组织表达基因富集-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library(reshape2)
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
unique(all_sample.combined@meta.data$celltype)
cell_data<-subset(all_sample.combined, subset = celltype  %in% "HSC")
data.list <- SplitObject(cell_data, split.by = "Tissue")

data <- data.list[["Bone_marrow"]]
data1 <- data.frame(data@assays[["RNA"]]@counts)
data1$sum <- rowSums(data1[,])
data1 <- data1[data1$sum != 0,]
Tissue <- data.frame(gene=rownames(data1),exp=data1[,"sum"])
colnames(Tissue)[1] <- "SYMBOL"
df_id <- bitr(Tissue$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Tissue <- merge(Tissue,df_id,by = "SYMBOL",all=F)
gene <- unique(Tissue[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result1<-as.data.frame(ego_BP1@result)
ego_BP_result1$Group <- rep("BP")
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "MF",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result2<-as.data.frame(ego_BP1@result)
ego_BP_result2$Group <- rep("MF")
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "CC",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result3<-as.data.frame(ego_BP1@result)
ego_BP_result3$Group <- rep("CC")
ego_BP_result <- rbind(ego_BP_result1,ego_BP_result2,ego_BP_result3)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

data <- data.list[["Spleen"]]
data1 <- data.frame(data@assays[["RNA"]]@counts)
data1$sum <- rowSums(data1[,])
data1 <- data1[data1$sum != 0,]
Tissue <- data.frame(gene=rownames(data1),exp=data1[,"sum"])
colnames(Tissue)[1] <- "SYMBOL"
df_id <- bitr(Tissue$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Tissue <- merge(Tissue,df_id,by = "SYMBOL",all=F)
gene <- unique(Tissue[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result1<-as.data.frame(ego_BP1@result)
ego_BP_result1$Group <- rep("BP")
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "MF",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result2<-as.data.frame(ego_BP1@result)
ego_BP_result3$Group <- rep("MF")
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "CC",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result3<-as.data.frame(ego_BP1@result)
ego_BP_result3$Group <- rep("CC")
ego_BP_result <- rbind(ego_BP_result1,ego_BP_result2,ego_BP_result3)
ego_BP_result <- ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Sleen","_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# 富集结果图 -------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
GO <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/Bone_marrow_GO.txt")
GO_lists <- c("GO:0002763","GO:0022409","GO:0045639","GO:0097529")
GO  <- GO[GO$ID %in% GO_lists,]
library(DOSE)
GO$GeneRatio1 <- parse_ratio(GO$GeneRatio)
p <- ggplot(GO,aes(reorder(Description,GeneRatio1),GeneRatio1,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_colour_gradient(low="#810000",high="#F56F6F")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
  geom_point(size = 2.0,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figureS6D1.pdf"),plot = p,width = 5,height = 5)
dev.off()

rm(list=ls())
GO <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/Sleen_GO.txt")
GO_lists <- c("GO:0002764","GO:0002252","GO:0002262","GO:0016887")
GO  <- GO[GO$ID %in% GO_lists,]
library(DOSE)
GO$GeneRatio1 <- parse_ratio(GO$GeneRatio)
p <- ggplot(GO,aes(reorder(Description,GeneRatio1),GeneRatio1,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_colour_gradient(low="#810000",high="#F56F6F")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
  geom_point(size = 2.0,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figureS6E1.pdf"),plot = p,width = 5,height = 5)
dev.off()

# 通路基因表达 ------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
unique(all_sample.combined@meta.data$celltype)
cell_data<-subset(all_sample.combined, subset = celltype  %in% "HSC")
data.list_tissue <- SplitObject(cell_data, split.by = "Tissue")

Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Bone_marrow","_old_vs_young.txt"), sep="") %>% .[.$change_V1 %in% c("UP","DOWN"),]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Cell.combined <- data.list_tissue[["Bone_marrow"]]
Cell.combined@meta.data[which(Cell.combined@meta.data$Age == "young"),5] <- "Young"
Cell.combined@meta.data[which(Cell.combined@meta.data$Age == "old"),5] <- "Aged"
Cell.combined$Age <- factor(Cell.combined$Age,levels = c("Young","Aged"))
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Age)
# gene <- c("ACIN1","ATP6AP1","CCR1","CD74","CREB1","EVI2B","FES","FOS","HAX1","HCLS1","IFNG","IL23A","IL34","KLF10","LEF1","NOTCH2",
#           "RB1","RUNX1","TESC","TGFB1","TMEM64","TNF","TNFSF11","TRAF6","TRIB1","AIF1","ALOX5","ANXA1","AP3D1","BAD","BCL10","BMP7","CARD11",
# "CBFB","CD28","CD40LG","CD44","CD46","CD47","CD74","CD80","CD83","CYRIB","DENND6A","DHPS","DPP4","DUSP10","EFNB1","EPHB6","F11R",
# "FLOT1","FOXP3","GATA3","HAS2","HES1","HLX","ICOS","IFNG","IL15","IL1B","IL23A","IL7R","ITPKB","JAK1","KIFAP3","LCK","LEF1","MALT1",
# "NCK1","NCK2","NCKAP1L","NKAP","NLRP3","PDCD1LG2","PODXL","PYCARD","RELA","RHOH","RIPK2","SART1","SHB","SLC7A1","SOCS1","SOCS5","SYK",
# "TFRC","TGFB1","TGFBR2","TNF","TNFSF11","TNFSF13B","TRAF6","VCAM1","VSIR","XBP1","ZAP70","ACIN1","ACVR1B","ANKRD54","ATP6AP1","CCR1",
# "CD74","CREB1","EVI2B","FES","FOS","GATA2","HAX1","HCLS1","HIF1A","HMGB2","IFNG","IL23A","IL34","ISG15","JAG1","KAT7","KLF10","LEF1",
# "MAPK14","NCKAP1L","NOTCH2","PITHD1","PRKDC","RB1","RUNX1","SCIN","TESC","TGFB1","TMEM64","TNF","TNFSF11","TRAF6","TRIB1","ADGRE2",
# "AIF1","AKIRIN1","ANXA1","C1QBP","C5AR2","CCL17","CCL2","CCL24","CCL3","CCL4L1","CCL5","CCR1","CCR2","CD300H","CD74","CD9","CD99L2",
# "CREB3","CSF3R","CXCL1","CXCL12","CXCL17","CXCL6","CXCL8","CXCR1","DNM1L","DPP4","DUSP1","EMILIN1","FCER1G","HRH1","IL1R1","IL23A",
# "IRAK4","JAGN1","JAML","LGALS3","LGMN","MAPK1","MCU","MOSPD2","MPP1","MSMP","MST1","NCKAP1L","NUP85","P2RY12","PDE4B","PLA2G7","PLCB1",
# "PPBP","PPIA","PTGER4","RAC2","RHOG","RHOH","S100A8","S100A9","SLC37A4","SPI1","SWAP70","SYK","TNFSF11","VEGFB")
# gene1 <- intersect(gene,Bone_marrow$gene)
gene1 <- c("FES","FOS","FLOT1","PLA2G7")
p1 <- VlnPlot(Cell.combined, cols = c('#39A8E3',"#D75156"),features = gene1, pt.size = 0,ncol = 2) + NoLegend()
p1
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figureS6D2.pdf"),plot = p1,width = 4,height = 6)
dev.off()

Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/DE_HSC/","Spleen","_old_vs_young.txt"), sep="") %>% .[.$change_V1 %in% c("UP","DOWN"),]
Spleen <- Spleen[!grepl("ENSMMUG", Spleen$gene),]
Cell.combined <- data.list_tissue[["Spleen"]]
Cell.combined@meta.data[which(Cell.combined@meta.data$Age == "young"),5] <- "Young"
Cell.combined@meta.data[which(Cell.combined@meta.data$Age == "old"),5] <- "Aged"
Cell.combined$Age <- factor(Cell.combined$Age,levels = c("Young","Aged"))
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Age)
# gene1 <- intersect(gene,Spleen$gene)
gene1 <- c("TYROBP","KHDRBS1","LAT2","NFAM1")
p1 <- VlnPlot(Cell.combined, cols = c('#39A8E3',"#D75156"),features = gene1, pt.size = 0,ncol = 2) + NoLegend()
p1
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/figureS6E2.pdf"),plot = p1,width = 4,height = 6)
dev.off()

# 粒细胞亚型 -------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "integrated"
celltype_lists <- unique(all_sample.combined@meta.data$celltype)
Cell.combined<-subset(all_sample.combined, subset = celltype  %in% "Granulocyte")
set.seed(1)
Cell.combined <- RunPCA(Cell.combined, npcs = 30, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 42, reduction = "pca", dims = 1:10)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:10)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.1,0.2,0.3,0.4))
p1 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE, label.size = 3)
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE, label.size = 3)
p <- (p1/p2)
p
DefaultAssay(Cell.combined) <- "RNA"
celltype_marker <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/celltype_marker.txt")
gene <- celltype_marker$Gene
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.0.4',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
p1
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  # scale_color_manual(values=color16)+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black")) + RotatedAxis()
p

celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/celltype.txt")
unique(celltype$celltype)
colnames(Cell.combined@meta.data)
table(Cell.combined@meta.data$integrated_snn_res.0.4)
cluster <- Cell.combined@meta.data[,c(1,22)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "sub_celltype")
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/Granulocyte.RData")

Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/Granulocyte.RData")
colors <- c("#140789","#38049A","#9B169F","#AB2494",
  "#F9983E","#FDAE32","#FBD324","#F0F921")
colors <- c("#2092E2","#008A8A","#7756A1","#FF6666","#74B404","#B3DE69","#DBBAFA","#FFD95A","#FFAF18")
colors <- c("#140789","#9B169F","#F9983E")

# show_col(colors)
p <- DimPlot(Cell.combined, reduction = "umap", pt.size = 0.1,group.by = "sub_celltype",cols = colors,
             label.size = 4,label = F,repel = T)+
  # theme(legend.position = 'none') +
  xlab("") + ylab("")
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/figure.S61.pdf",plot = p,dpi=1000,width = 5,height = 4)

# figureS4B 粒细胞 features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/Granulocyte.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$sub_celltype)
celltype_marker <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/celltype_marker.txt")
gene <- celltype_marker$Gene
gene <- c("ITGAM","ITGB2","CD55","CEACAM8","PTPRC","CD19","CD22","CD63","ENPP3","FCER1A")
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'sub_celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
data$Cluster <- factor(data$Cluster,levels = rev(c("Neutrophils","Eosinophils","Basophils")))
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  # scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  scale_color_gradientn(colors = c("grey","grey","#eff3ff","#F76C6C","#981216"))+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black"),
        # legend.position = "none",
        axis.text.x=element_text(angle = 90,hjust = 1,vjust = 1))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/figure.S62.pdf",plot = p,dpi=1000,width = 4,height = 3)

# figure粒细胞 FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/Granulocyte.RData")
colnames(Cell.combined@meta.data)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,23)]
celltype_tissue$sub_celltype <- as.character(celltype_tissue$sub_celltype)
unique(celltype_tissue$sub_celltype)
cell_list <-  c("Neutrophils","Eosinophils","Basophils")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$sub_celltype == j),]))
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               FC = as.character(cell[which(cell$Age == "old"),"Freq"]/cell[which(cell$Age == "young"),"Freq"]))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(nrow(cell) == 0){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = "0",
                               young_cell = "0",
                               FC = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)  
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = "0",
                               FC = as.character(cell[which(cell$Age == "old"),"Freq"]))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,sub_celltype = j,
                               old_cell = "0",
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               FC = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }
  }
  # freq_data_tissue$old_per <- freq_data_tissue$old_cell/sum(as.numeric(freq_data_tissue$old_cell))
  # freq_data_tissue$young_per <- freq_data_tissue$young_cell/sum(as.numeric(freq_data_tissue$young_cell))
  freq_data <- rbind(freq_data,freq_data_tissue)
}
freq_data$log2FC <- log2(as.numeric(freq_data$FC))


# ###fisher检验
fisher_data <- data.frame()
tissue_list <- unique(celltype_tissue$Tissue)
tissue_list <- c("Spleen","Bone_marrow")
cell_list <-  c("Neutrophils","Eosinophils","Basophils")
for (tissue in tissue_list){
  test_data <- freq_data[freq_data$Tissue == tissue,]
  for (celltype in cell_list){
    test_data1 <- test_data[test_data$sub_celltype == celltype,c(3,4)]
    test_data2 <- test_data[test_data$sub_celltype != celltype,]
    test_data3 <- data.frame(old_cell = sum(as.numeric(test_data2$old_cell)),young_cell = sum(as.numeric(test_data2$young_cell)))
    test_data4 <- rbind(test_data1,test_data3)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,celltype = celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))
library("ggplot2")
unique(freq_data$sub_celltype)
freq_data$sub_celltype <- factor(freq_data$sub_celltype,levels = c("Neutrophils","Eosinophils","Basophils"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#140789","#9B169F","#F9983E")
show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
# freq_data$Tissue <- gsub("Lymph_node","Mesenteric_lymph",freq_data$Tissue)
freq_data$log2FC[!is.finite(freq_data$log2FC)] <- 0
p <- ggplot(data=freq_data, aes(x=sub_celltype, y=log2FC, fill = sub_celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+ 
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"))+
  xlab("") + ylab("")
p
ggsave(plot = p, paste0('/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/','FC.pdf'), width = 6, height = 3,dpi = 600)
dev.off()

# figure5C  占比饼图 ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/Granulocyte.RData")
colnames(Cell.combined@meta.data)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,23)]
celltype_tissue$sub_celltype <- as.character(celltype_tissue$sub_celltype)
unique(celltype_tissue$sub_celltype)
cell_list <-  c("Neutrophils","Eosinophils","Basophils")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  tissue1 <- tissue[tissue$Age == "young",]
  tissue2 <- tissue[tissue$Age == "old",]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$sub_celltype == j),]))
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(nrow(cell) == 0){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell_per = "0",
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell_per = "0",
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }
  }
  freq_data <- rbind(freq_data,freq_data_tissue)
}
str(freq_data)
freq_data$old_cell_per <- as.numeric(freq_data$old_cell_per);freq_data$young_cell_per <- as.numeric(freq_data$young_cell_per)
freq_data1 <- melt(freq_data,id.vars = c("Tissue","celltype"))
colnames(freq_data1) <- c("Tissue","celltype","Age","cell_per")
freq_data1$Age <- gsub("young_cell_per","Young",freq_data1$Age)
freq_data1$Age <- gsub("old_cell_per","Aged",freq_data1$Age)
freq_data2 <- spread(freq_data1,celltype,cell_per)
freq_data2$x=c(4,2,8,6,12,10,16,14)
freq_data2$y=rep(2)
freq_data2$region=as.character(1:8)
colnames(freq_data2)
freq_data3 <- freq_data2[,c(1,2,6,7,3,4,5,8)]
colors <- c("#140789","#9B169F","#F9983E")
p <- ggplot()+
  geom_scatterpie(data=freq_data2,
                  aes(x,y,group=region,r=0.9),
                  cols = cell_list)+
  coord_equal()+
  theme_void()+
  # theme(legend.position = "none")+
  scale_fill_manual(values = colors)
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/pie.pdf",plot = p,width = 15,height = 7)
dev.off()


# GZMB 表达--------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/Granulocyte.RData")
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$sub_celltype)
data.list_age <- SplitObject(Cell.combined, split.by = "Age")
age_list <- unique(Cell.combined$Age)
data_mean <- data.frame()
for (age in age_list){
  data_age <- data.list_age[[age]]
  data.list_tissue <- SplitObject(data_age, split.by = "Tissue")
  tissue_list <- unique(data_age$Tissue)
  for (tissue in tissue_list){
    data_tissue <- data.list_tissue[[tissue]]
    data.list_celltype <- SplitObject(data_tissue, split.by = "sub_celltype")
    celltype_list <- unique(data_tissue$sub_celltype)
    for (celltype in celltype_list){
      data_celltype <- data.list_celltype[[celltype]]
      DefaultAssay(data_celltype) <- "RNA"
      data <- data.frame(GetAssayData(data_celltype, slot = "data"));data$gene <- rownames(data)
      data <- data[data$gene == "GZMB",]
      data1 <- data[,-ncol(data)]
      sum(as.numeric(data1))
      if (sum(as.numeric(data1)) == 0){
        data2 <- data.frame(avg.exp=0,Cluster=celltype,Age=age,Tissue=tissue)
        data_mean <- rbind(data_mean,data2)
      }else if(sum(as.numeric(data1)) != 0){
        data2 <- data.frame(avg.exp=sum(as.numeric(data1))/(ncol(data)-1),Cluster=celltype,Age=age,Tissue=tissue)
        data_mean <- rbind(data_mean,data2) 
      }
    }
  }
}
write.table(data_mean,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure6/","GZMB.txt"),quote = FALSE,sep = "\t",row.names = F)

rm(list=ls())
data_mean <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure6/GZMB.txt")
young <- data_mean[data_mean$Age == "young",]
young <- young[order(young$Tissue, young$Cluster),]
old <- data_mean[data_mean$Age == "old",]
old <- old[order(old$Tissue, old$Cluster),]
table(old$Tissue)
data <- cbind(young,old)
colnames(data) <- c("avg.exp_young","Cluster","Age","Tissue","avg.exp_old","Cluster","Age","Tissue")
data <- data[,c(4,2,1,5)]
data$fc <- data$avg.exp_old - data$avg.exp_young

tg.sub <- data[,c(1,2,5)]
dt <- reshape2::recast(tg.sub,Tissue ~ Cluster)
dt[is.na(dt)] <- 0
dt.pct <- dt
for (j in 2:ncol(dt)) {
  if (sum(dt[,j] != 0)) {
    dt.pct[,j] <- dt[,j]/sum(dt[,j])
  }
}
tissue_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
dt.pct$Tissue <- factor(dt.pct$Tissue,levels = tissue_list)
library(ggradar)
p1 <- ggradar(dt,axis.label.size=6.5,
              # values.radar = c(0, ceiling(max(dt[,-1]))/2, ceiling(max(dt[,-1]))),
              values.radar = c(0, 2, 4),
              grid.min = 0,
              grid.mid = 2,
              grid.max = 4,
              group.colours = c('#23B363','#DF63CC','#339FBA','#CA7D44')) +
  ggtitle("GZMB mean expression (Granulocyte)")+
  theme(legend.position = "none")
p1
ggsave(plot = p1, paste0('/data1/zhuzn/wangsn/reviewer_4000/figure6/Granulocyte/','GZMB.pdf'), width = 14, height = 16,dpi = 600)
