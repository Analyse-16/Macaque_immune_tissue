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

# ############# T细胞分类 #############---------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "integrated"
unique(all_sample.combined@meta.data$celltype)
Cell.combined<-subset(all_sample.combined, subset = celltype  %in% "T Cell")
Cell.combined <- RunPCA(Cell.combined, npcs = 30, verbose = T,seed.use = 42)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 2, reduction = "pca", dims = 1:10)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:10)
Cell.combined <- FindClusters(Cell.combined, resolution = c(1.2,1.4,1.6,1.8,2))
p1 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.1.4", label = TRUE,raster=FALSE)
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.1.6", label = TRUE,raster=FALSE)
p3 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.1.8", label = TRUE,raster=FALSE)
p4 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.2", label = TRUE,raster=FALSE)
p <- p<-(p1|p2)/(p3|p4)
p
save(Cell.combined, file="/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T_cluster.RData")

#添加细胞类型
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
all_cell <- Cell.combined@meta.data[,c(4,5,18)]
all_cell %<>% {.$id<-rownames(.);.}
all_cell$id <-gsub("-",".",all_cell$id)
all_cell$id <-gsub("Lymph_node","Mesenteric_lymph",all_cell$id)

load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T_cluster.RData")
colnames(Cell.combined@meta.data)
all_cell1 <- Cell.combined@meta.data[,c(4,5)]
all_cell1 %<>% {.$id<-rownames(.);.} %>% join(.,all_cell[,c("id","celltype")],by=c("id"="id"))
d <- intersect(all_cell$id,all_cell1$id)
colnames(all_cell1) <- c("Tissue","Age","id","Celltype_old")
str(all_cell1)
all_cell1$Celltype_old <- as.character(all_cell1$Celltype_old)
all_cell1[is.na(all_cell1)] <- "NON"
all_cell2 <- all_cell1[all_cell1$Celltype_old == "NON",]
all_cell2 <- separate(all_cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined <- AddMetaData(Cell.combined,all_cell1$Celltype_old,col.name = "Celltype_old")
unique(Cell.combined@meta.data$Celltype_old)
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "Celltype_old", label = TRUE,raster=FALSE)
p
save(Cell.combined, file="/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T_cluster.RData")

### CD4-CD8
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
gene <- c("CD4","CD8A","CD8B")
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.0.6)
###FeaturePlot
library(patchwork)
p <- FeaturePlot(Cell.combined, features = gene,min.cutoff = 0, max.cutoff = 4,reduction = "umap",cols = c("grey", "red"),
                 pt.size = 0.01,raster=FALSE,label = F)+
  plot_annotation(theme = theme(plot.title = element_text(size = 20,hjust = 0.5,face = "bold")))
p

celltype <- read.delim("/data1/zhuzn/wangsn/data_QC_test/celltype_4000_T.txt")
colnames(Cell.combined@meta.data)
table(Cell.combined@meta.data$integrated_snn_res.1.2)
cluster <- Cell.combined@meta.data[,c(1,10)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "celltype")
p1 <- DimPlot(Cell.combined, reduction = "umap", group.by = "Celltype_old", label = TRUE,raster=FALSE)
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype", label = TRUE,raster=FALSE)
# Cell.combined$celltype1 <- ifelse(Cell.combined$Celltype_old == "NON" ,Cell.combined$celltype,Cell.combined$Celltype_old)
# p3 <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype1", label = TRUE,raster=FALSE)
p <- p1|p2
p
save(Cell.combined, file="/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T.RData")

Cell.combined_CD4 <- subset(Cell.combined, subset = celltype == "CD4")
save(Cell.combined_CD4, file="/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_CD4.RData")
Cell.combined_CD8 <- subset(Cell.combined, subset = celltype == "CD8")
save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_CD8.RData")

# CD8重新聚类分簇-------------------------------------------------------------------------
#添加细胞类型
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
all_cell <- Cell.combined@meta.data[,c(4,5,18)]
all_cell %<>% {.$id<-rownames(.);.}
all_cell$id <-gsub("-",".",all_cell$id)
all_cell$id <-gsub("Lymph_node","Mesenteric_lymph",all_cell$id)

load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
all_cell1 <- Cell.combined_CD8@meta.data[,c(4,5)]
all_cell1 %<>% {.$id<-rownames(.);.} %>% join(.,all_cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell1) <- c("Tissue","Age","id","Celltype_old")
all_cell1[is.na(all_cell1)] <- "NON"
all_cell2 <- all_cell1[all_cell1$Celltype_old == "NON",]
all_cell2 <- separate(all_cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined_CD8 <- AddMetaData(Cell.combined_CD8,all_cell1$Celltype_old,col.name = "Celltype_old")
unique(Cell.combined_CD8@meta.data$Celltype_old)
save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8_cluster.RData")

# 重新聚类
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8_cluster.RData")
DefaultAssay(Cell.combined_CD8) <- "integrated"
unique(Cell.combined_CD8@meta.data$celltype)
unique(Cell.combined_CD8$Celltype_old)
Cell.combined_CD8 <- RunPCA(Cell.combined_CD8, npcs = 30, verbose = T)
ElbowPlot(Cell.combined_CD8)
Cell.combined_CD8 <- RunUMAP(Cell.combined_CD8, seed.use = 2, reduction = "pca", dims = 1:17)
Cell.combined_CD8 <- FindNeighbors(Cell.combined_CD8, reduction = "pca", dims = 1:17)
Cell.combined_CD8 <- FindClusters(Cell.combined_CD8, resolution = c(1,1.2,1.4,1.6,1.8,2))
p1 <- DimPlot(Cell.combined_CD8, reduction = "umap",group.by = "integrated_snn_res.2",pt.size = 0.01,label = T)
p1
Cell.combined_CD8@meta.data$Celltype_old <- factor(Cell.combined_CD8@meta.data$Celltype_old,levels = c("Naive","Tem","Tcm","Tex","CTL/TRM",
                                                                                                       "TRM","CTL","NON"))
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$Celltype_old)
Cell.combined_CD8$Tissue_Age <- paste0(Cell.combined_CD8$Tissue,"_",Cell.combined_CD8$Age)
colors <- c("#46ADE3","#92DF4D","#339933","#E3393B","#FB7473","#FF6600","#B064D6","black")
# Cell.combined_CD8 <- RunUMAP(Cell.combined_CD8, seed.use = 222, reduction = "pca", dims = 1:15)
p1 <- DimPlot(Cell.combined_CD8, reduction = "umap",cols = colors,group.by = "Celltype_old",
              # split.by = "Age",ncol = 2,
              pt.size = 0.01,label = F)
p1
save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8_cluster.RData")

# test
for (i in c(10:30)){
  Cell.combined_CD8 <- RunUMAP(Cell.combined_CD8, seed.use = 2,reduction = "pca", dims = 1:i)
  p <- DimPlot(Cell.combined_CD8, reduction = "umap", group.by = "Celltype_old",pt.size = 0.01,
               label = T, repel = TRUE, label.size = 5,raster=F)+
    theme(legend.position = "right")
  p
  ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS4/round_4000/",i,".png"),p,width = 10,height = 8)
}

# figureS4 CD8 features_Bubble_cluster-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8_cluster.RData")
Cell.combined <- Cell.combined_CD8
DefaultAssay(Cell.combined) <- "RNA"

Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.1.2)
cell_markers <- read.delim("/data1/zhuzn/wangsn/reviewer/CD8/celltype_marker.txt")
gene <- cell_markers$Gene
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.1.2',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
library(plyr)
data1 <- join(data,cell_markers)
getAnywhere(join)
unique(data1$celltype)
p <- ggplot(data=data1,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`))+ geom_point() +
  facet_wrap(.~celltype,scales= "free",ncol= 3)+
  scale_size_continuous(range=c(1,5))+
  scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ffcc29",midpoint = 0)+
  theme_bw() +
  labs(title = "",x = '', y = '')+theme(panel.grid = element_blank()) +
  theme(axis.text=element_text(size=11, color="black"), plot.title = element_text(angle = 90)) + RotatedAxis()
p
dev.off()

# figureS4A CD8 添加细胞类型信息 --------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8_cluster.RData")
celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS4/celltype.txt")
colnames(Cell.combined_CD8@meta.data)
cluster <- Cell.combined_CD8@meta.data[,c(1,10)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined_CD8 <- AddMetaData(Cell.combined_CD8,cluster_celltype$celltype,col.name = "celltype")
Cell.combined_CD8@meta.data$Age <- factor(Cell.combined_CD8@meta.data$Age,levels = c("young","old"))
Cell.combined_CD8@meta.data$celltype <- factor(Cell.combined_CD8@meta.data$celltype,levels = c("Naive","Tem","Tcm","Tex","CTL/TRM",
                                                                                               "TRM","CTL"))
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype)
Cell.combined_CD8$Tissue_Age <- paste0(Cell.combined_CD8$Tissue,"_",Cell.combined_CD8$Age)
colors <- c("#46ADE3","#92DF4D","#339933","#E3393B","#FB7473","#FF6600","#B064D6")
p1 <- DimPlot(Cell.combined_CD8, reduction = "umap",cols = colors,group.by = "celltype",
              split.by = "Age",ncol = 2,
              pt.size = 0.01,label = F)
p1
save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")

# 挑选一部分cluster14细胞(分为CD4) --------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
table(Cell.combined_CD8@meta.data$integrated_snn_res.1.2)
colnames(Cell.combined_CD8@meta.data)
# 获取该簇中的细胞索引
cells <- Cell.combined_CD8@meta.data[,c(1,10)]
cells$id <- rownames(cells)
cells_in_cluster <- cells[which(cells$integrated_snn_res.1.2 == 14),3]
# 例如，选择1000个细胞
random_cells <- sample(cells_in_cluster, size = 1000)
data.filt <- subset(Cell.combined_CD8,cells=setdiff(WhichCells(Cell.combined_CD8),random_cells))
colors <- c("#46ADE3","#92DF4D","#339933","#E3393B","#FB7473","#FF6600","#B064D6")
p1 <- DimPlot(data.filt, reduction = "umap",cols = colors,group.by = "celltype",
              split.by = "Age",ncol = 2,
              pt.size = 0.01,label = F)
p1
Cell.combined_CD41 <- subset(Cell.combined_CD8,cells=intersect(WhichCells(Cell.combined_CD8),random_cells))
save(Cell.combined_CD41, file="/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD41.RData")

Cell.combined_CD8 <- data.filt
save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")

# figureS4A,3A,4A CD8 celltype umap --------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
Cell.combined_CD8@meta.data$celltype <- factor(Cell.combined_CD8@meta.data$celltype,levels = c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL"))
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype)
colors <- c("#46ADE3","#92DF4D","#339933","#E3393B","#FB7473","#FF6600","#B064D6")
colors <- c("#FF6600","#FB7473","#339933","#E3393B","#92DF4D","#46ADE3","#B064D6")
show_col(colors)
p <- DimPlot(Cell.combined_CD8, reduction = "umap", cols = colors, raster=FALSE,repel = F, label = T)+xlab("") + ylab("") + 
  theme(axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none',
        axis.text = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4A.pdf",plot = p,dpi=1000,width = 8,height = 8)
# save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")

data.list_age <- SplitObject(Cell.combined_CD8, split.by = "Age")
young_data <- data.list_age[["young"]]
old_data <- data.list_age[["old"]]
show_col(colors)
p1 <- DimPlot(young_data, reduction = "umap",cols = colors, raster=FALSE,repel = F, label = F)+xlab("") + ylab("") + 
  theme(axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none',
        axis.text = element_blank())
p2 <- DimPlot(old_data, reduction = "umap", cols = colors, raster=FALSE,repel = F, label = F)+xlab("") + ylab("") + 
  theme(axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none',
        axis.text = element_blank())
p <- p1/p2
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure3/figure.3-4A.pdf",plot = p,dpi=1000,width = 8,height = 16)

# figureS4B CD8 features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
DefaultAssay(Cell.combined_CD8) <- "RNA"
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype)

gene <- c("CCR7","CCR6","GZMK","CD69","HAVCR2","CTLA4","ITGAE","EOMES","RUNX3","GZMB")
# gene <- c("CCR7","CCR6","GZMK","CD69","CTLA4","ITGAE","EOMES")
p1<-DotPlot(object = Cell.combined_CD8, features = gene, group.by  = 'celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
data$Cluster <- factor(data$Cluster,levels = rev(c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL")))
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  # scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  scale_color_gradientn(colors = c("grey","grey","#eff3ff","#F76C6C","#981216"))+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black"),axis.text.x=element_text(angle = 90,hjust = 1,vjust = 1),legend.position = "none")
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4B.pdf",plot = p,dpi=1000,width = 6.3,height = 5.3)

# figureS4C FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
celltype_tissue <- Cell.combined_CD8@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
unique(celltype_tissue$celltype)
cell_list <-  c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$celltype == j),]))
    freq_data1 <- data.frame(Tissue = i,celltype = j,
                             old_cell = cell[which(cell$Age == "old"),"Freq"],
                             young_cell = cell[which(cell$Age == "young"),"Freq"],
                             FC = as.character(cell[which(cell$Age == "old"),"Freq"]/cell[which(cell$Age == "young"),"Freq"]))
    freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
  }
  freq_data_tissue$old_per <- freq_data_tissue$old_cell/sum(freq_data_tissue$old_cell)
  freq_data_tissue$young_per <- freq_data_tissue$young_cell/sum(freq_data_tissue$young_cell)
  freq_data <- rbind(freq_data,freq_data_tissue)
}
freq_data$log2FC <- log2(as.numeric(freq_data$FC))

###fisher检验
fisher_data <- data.frame()
for (tissue in unique(celltype_tissue$Tissue)){
  test_data <- freq_data[freq_data$Tissue == tissue,]
  for (celltype in cell_list){
    test_data1 <- test_data[test_data$celltype == celltype,c(3,4)]
    test_data2 <- test_data[test_data$celltype != celltype,]
    test_data3 <- data.frame(old_cell = sum(test_data2$old_cell),young_cell = sum(test_data2$young_cell))
    test_data4 <- rbind(test_data1,test_data3)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,celltype = celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))
library("ggplot2")
unique(freq_data$celltype)
freq_data$celltype <- factor(freq_data$celltype,levels = c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#FF6600","#FB7473","#339933","#E3393B","#92DF4D","#46ADE3","#B064D6")
show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
freq_data$Tissue <- gsub("Lymph_node","Mesenteric_lymph",freq_data$Tissue)
p <- ggplot(data=freq_data, aes(x=celltype, y=log2FC, fill = celltype, width=0.8))+
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
p1<-p+geom_text(aes(x=celltype, y=log2FC,label=p.value1),size=6,position= position_dodge(0.6))
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4C.pdf",plot = p1,dpi=1000,width = 10,height = 4)
dev.off()

# figure5C  占比饼图 ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
celltype_tissue <- Cell.combined_CD8@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
unique(celltype_tissue$celltype)
cell_list <-  c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  tissue1 <- tissue[tissue$Age == "young",]
  tissue2 <- tissue[tissue$Age == "old",]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$celltype == j),]))
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
freq_data3 <- freq_data2[,c(1,2,10,11,5,7,6,8,4,9,3,12)]
colors <- c("#FF6600","#FB7473","#339933","#E3393B","#92DF4D","#46ADE3","#B064D6")
p <- ggplot()+
  geom_scatterpie(data=freq_data2,
                  aes(x,y,group=region,r=0.9),
                  cols = cell_list)+
  coord_equal()+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values = colors)
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figureS4C1.pdf",plot = p,width = 15,height = 7)
dev.off()

# figureS4C percent_FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
celltype_tissue <- Cell.combined_CD8@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
unique(celltype_tissue$celltype)
cell_list <-  c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  tissue1 <- tissue[tissue$Age == "young",]
  tissue2 <- tissue[tissue$Age == "old",]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$celltype == j),]))
    
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = (cell[which(cell$Age == "young"),"Freq"])/nrow(tissue1))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(nrow(cell) == 0){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell = "0",
                               young_cell = "0",
                               old_cell_per = "0",
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = "0",
                               old_cell_per = (cell[which(cell$Age == "old"),"Freq"])/nrow(tissue2),
                               young_cell_per = "0")
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
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
  for (celltype in cell_list){
    test_data1 <- test_data[test_data$celltype == celltype,c(3,4)]
    test_data2 <- test_data[test_data$celltype != celltype,]
    test_data3 <- data.frame(old_cell = sum(test_data2$old_cell),young_cell = sum(test_data2$young_cell))
    test_data4 <- rbind(test_data1,test_data3)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,celltype = celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))
library("ggplot2")
unique(freq_data$celltype)
freq_data$celltype <- factor(freq_data$celltype,levels = c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#46ADE3","#92DF4D","#339933","#E3393B","#FB7473","#FF6600","#B064D6")
show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
freq_data$Tissue <- gsub("Lymph_node","Mesenteric_lymph",freq_data$Tissue)
p <- ggplot(data=freq_data, aes(x=celltype, y=log2FC, fill = celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+ 
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,colour = "black"),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=celltype, y=log2FC,label=p.value1),size=6,position= position_dodge(0.6))
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4C.pdf",plot = p1,dpi=1000,width = 10,height = 4)
dev.off()

#  figureS4D CD8 （ggplot）-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
Cell.combined_CD8<-AddMetaData(Cell.combined_CD8,Cell.combined_CD8@reductions$umap@cell.embeddings,col.name = colnames(Cell.combined_CD8@reductions$umap@cell.embeddings))
colnames(Cell.combined_CD8@meta.data)

Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype)
DefaultAssay(Cell.combined_CD8) <- "RNA"
data.list_age <- SplitObject(Cell.combined_CD8, split.by = "Age")
# genes <- c("S100A4","PDCD1","CCL3")
genes <- c("GADD45A","BHLHE40","IL18")
young_data <- data.list_age[["young"]]
mat1=as.data.frame(young_data[["RNA"]]@data["GADD45A",])

colnames(mat1)="exp"
mat2=Embeddings(young_data,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat3$exp[mat3$exp > 4] = 4
colors=rev(c("#A10000","#B30000","#CE1212","#F03B3B","#F05A5A","#F07979","#FFAFAF20"))
p1 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p1
# ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4D.pdf",plot = p1,width = 6,height = 6)#图例文件

for (gene in genes[2:3]){
  mat1=as.data.frame(young_data[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(young_data,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat3$exp[mat3$exp > 4] = 4
  p12 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p1 <- p1|p12
}
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4D_young2.pdf",plot = p1,width = 17,height = 6)

old_data <- data.list_age[["old"]]
mat1=as.data.frame(old_data[["RNA"]]@data["GADD45A",])

colnames(mat1)="exp"
mat2=Embeddings(old_data,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat3$exp[mat3$exp > 4] = 4
p2 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
for (gene in genes[2:3]){
  mat1=as.data.frame(old_data[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(old_data,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat3$exp[mat3$exp > 4] = 4
  p22 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p2 <- p2|p22
}
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS4/figure.S4D_old2.pdf",plot = p2,width = 17,height = 6)

#  figureS4D CD8 （ggplot-TEST）-------------------------------------------------------------------------
rm(list=ls())
diff_young <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_young.txt")
diff_young <- diff_young[diff_young$change_V1 == "UP",]
diff_young <- diff_young[!grepl("ENSMMUG", diff_young$gene),]
diff_old <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_old.txt")
diff_old <- diff_old[diff_old$change_V1 == "UP",]
diff_old <- diff_old[!grepl("ENSMMUG", diff_old$gene),]
diff_old_young <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_old_vs_young.txt")
diff_old_young <- diff_old_young[diff_old_young$change_V1 == "UP",]
diff_old_young <- diff_old_young[!grepl("ENSMMUG", diff_old_young$gene),]
overlap <- Reduce(intersect, list(diff_young$gene,diff_old$gene,diff_old_young$gene))
genes <- c("GZMB","PDCD1","CX3CR1","ADGRG1","RGS9","ZEB2","S100A4","CCL5","CCL3","GADD45A","BHLHE40","IL18",
           "APP","CTNNB1","MAPK1","ARF1","JUNB","LEF1","CCR7","S100A11","MAF","IKZF2","NFKBIA","S100A9","S100A8","JUN","FOS")
g <- setdiff(overlap,genes)

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
p1 <- FeaturePlot(Cell.combined_CD8, features = genes,reduction = "umap",cols = c("grey", "red"),raster=FALSE,min.cutoff = 0, max.cutoff = 5,pt.size = 0.01,label = TRUE)+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p1

colnames(Cell.combined_CD8@meta.data)
Cell.combined_CD8<-AddMetaData(Cell.combined_CD8,Cell.combined_CD8@reductions$umap@cell.embeddings,col.name = colnames(Cell.combined_CD8@reductions$umap@cell.embeddings))
colnames(Cell.combined_CD8@meta.data)

Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype)
DefaultAssay(Cell.combined_CD8) <- "RNA"
data.list_age <- SplitObject(Cell.combined_CD8, split.by = "Age")
genes <- c("PDCD1","S100A4","CCL3","GADD45A","BHLHE40","IL18")
young_data <- data.list_age[["young"]]
mat1=as.data.frame(young_data[["RNA"]]@data["PDCD1",])

colnames(mat1)="exp"
mat2=Embeddings(young_data,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat3$exp[mat3$exp > 4] = 4
colors=rev(c("#B30000","#CE1212","#F03B3B","#F05A5A","#F07979","#FFAFAF","#FFDBDB","#FFE6E6"))
show_col(colors)
p1 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p1

for (gene in genes[2:6]){
  mat1=as.data.frame(young_data[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(young_data,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat3$exp[mat3$exp > 4] = 4
  p12 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p1 <- p1|p12
}

old_data <- data.list_age[["old"]]
mat1=as.data.frame(old_data[["RNA"]]@data["PDCD1",])

colnames(mat1)="exp"
mat2=Embeddings(old_data,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat3$exp[mat3$exp > 4] = 4
p2 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
for (gene in genes[2:6]){
  mat1=as.data.frame(old_data[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(old_data,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat3$exp[mat3$exp > 4] = 4
  p22 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p2 <- p2|p22
}


