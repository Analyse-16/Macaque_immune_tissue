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

# figure5 A 4000阈值数据 B细胞 ------------------------------------------------------------
#添加细胞类型
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/B Cell/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
all_cell <- Cell.combined@meta.data[,c(4,5,23)]
all_cell %<>% {.$id<-rownames(.);.}
all_cell$id <-gsub("-",".",all_cell$id)
all_cell$id <-gsub("Lymph_node","Mesenteric_lymph",all_cell$id)

load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "integrated"
unique(all_sample.combined@meta.data$celltype)
Cell.combined_B<-subset(all_sample.combined, subset = celltype  %in% "B Cell")

all_cell1 <- Cell.combined_B@meta.data[,c(4,5)]
all_cell1 %<>% {.$id<-rownames(.);.} %>% join(.,all_cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell1) <- c("Tissue","Age","id","Celltype_old")
all_cell1[is.na(all_cell1)] <- "NON"
all_cell2 <- all_cell1[all_cell1$Celltype_old == "NON",]
all_cell2 <- separate(all_cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined_B <- AddMetaData(Cell.combined_B,all_cell1$Celltype_old,col.name = "Celltype_old")
unique(Cell.combined_B@meta.data$Celltype_old)

Cell.combined_B@meta.data$Celltype_old <- factor(Cell.combined_B@meta.data$Celltype_old,levels = c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell","NON"))
Cell.combined_B <- SetIdent(object = Cell.combined_B, value = Cell.combined_B@meta.data$Celltype_old)
colors <- c("#B3DE69","#FDB462","#1084D3","#FCCDE5","#8DD3C7","#80B1D3","#FB8072","black")
p <- DimPlot(Cell.combined_B, reduction = "umap", cols = colors, raster=FALSE,repel = F, label = T)+xlab("") + ylab("") + 
  theme(axis.ticks = element_blank(),axis.line = element_blank(),
        # legend.position = 'none',
        axis.text = element_blank())
p
save(Cell.combined_B, file="/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B.RData")

Cell.combined_B <- RunPCA(Cell.combined_B, npcs = 30, verbose = T)
Cell.combined_B <- RunUMAP(Cell.combined_B, seed.use = 2, reduction = "pca", dims = 1:20)
Cell.combined_B <- FindNeighbors(Cell.combined_B, reduction = "pca", dims = 1:20)
Cell.combined_B <- FindClusters(Cell.combined_B, resolution = c(0.4))
p <- DimPlot(Cell.combined_B, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE,raster=FALSE)
p
Cell.combined_B$Tissue_Age <- paste0(Cell.combined_B$Tissue,"_",Cell.combined_B$Age)
p1 <- DimPlot(Cell.combined_B, reduction = "umap",cols = colors,group.by = "Celltype_old",label = TRUE)
p1
save(Cell.combined_B, file="/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B.RData")

# figure5 A B细胞亚型分类 -----------------------------------------------------------------
rm(list = ls())
celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure5/celltype_4000_B.txt")
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B.RData")
p <- DimPlot(Cell.combined_B, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE,raster=FALSE)
p
colnames(Cell.combined_B@meta.data)
table(Cell.combined_B@meta.data$integrated_snn_res.0.4)
cluster <- Cell.combined_B@meta.data[,c(1,10)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined_B <- AddMetaData(Cell.combined_B,cluster_celltype$celltype,col.name = "celltype")
save(Cell.combined_B, file="/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")

Cell.combined_B@meta.data$celltype <- factor(Cell.combined_B@meta.data$celltype,levels = c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell"))
colors <- c("#B3DE69","#FDB462","#1084D3","#FCCDE5","#8DD3C7","#80B1D3","#FB8072")
show_col(colors)
p <- DimPlot(Cell.combined_B, reduction = "umap", group.by = "celltype", cols = colors,label = TRUE,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p

# figure5 A tissue.young_old_umap-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
unique(Cell.combined_B$Tissue_Age)
Cell.combined_B@meta.data$Tissue_Age <- factor(Cell.combined_B@meta.data$Tissue_Age,levels = c("Bone_marrow_young","Mesenteric_lymph_young","PBMC_young","Spleen_young",
                                                                                               "Bone_marrow_old","Mesenteric_lymph_old","PBMC_old","Spleen_old"))

Cell.combined_B@meta.data$celltype <- factor(Cell.combined_B@meta.data$celltype,levels = c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell"))
Cell.combined_B <- SetIdent(object = Cell.combined_B, value = Cell.combined_B@meta.data$celltype)
colors <- c("#B3DE69","#FDB462","#1084D3","#FCCDE5","#8DD3C7","#80B1D3","#FB8072")
show_col(colors)
p <- DimPlot(Cell.combined_B, reduction = "umap",cols = colors,group.by = "celltype",
              split.by = "Tissue_Age",ncol = 4,pt.size = 0.01,label = F)
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure5/","figure5A.pdf"),plot = p,width = 13,height = 7.5)

# figure5 A B细胞统计 FC(all_tissue.young_old)-------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
colnames(Cell.combined_B@meta.data)
celltype_tissue <- Cell.combined_B@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
cell_list <- c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell")
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
freq_data$celltype <- factor(freq_data$celltype,levels = c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#B3DE69","#FDB462","#1084D3","#FCCDE5","#8DD3C7","#80B1D3","#FB8072")
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]

p <- ggplot(data=freq_data, aes(x=celltype, y=log2FC, fill = celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+ 
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,colour = "black"),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"),
        panel.grid = element_blank(),strip.background = NULL,
        strip.text = element_text(size = 15))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=celltype, y=log2FC,label=p.value1),size=6,position= position_dodge(0.6))
p1
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure5/',"figure5A",'_freq.pdf')
ggsave(plotfile, plot=p1, dpi = 1000, width = 10, height = 4)
dev.off()

# figure5C  占比饼图 ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
colnames(Cell.combined_B@meta.data)
celltype_tissue <- Cell.combined_B@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
unique(celltype_tissue$celltype)
cell_list <- c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell")
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
freq_data3 <- freq_data2[,c(1,2,10,11,9,8,5,4,3,7,6,12)]
colors <- c("#B3DE69","#FDB462","#1084D3","#FCCDE5","#8DD3C7","#80B1D3","#FB8072")
p <- ggplot()+
  geom_scatterpie(data=freq_data2,
                  aes(x,y,group=region,r=0.9),
                  cols = cell_list)+
  coord_equal()+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values = colors)
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure5/figure5_freq1.pdf",plot = p,width = 15,height = 7)
dev.off()

# figure5 A B细胞统计 percent_FC(all_tissue.young_old)-------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
colnames(Cell.combined_B@meta.data)
celltype_tissue <- Cell.combined_B@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
cell_list <- c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell")
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
freq_data$celltype <- factor(freq_data$celltype,levels = c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#B3DE69","#FDB462","#1084D3","#FCCDE5","#8DD3C7","#80B1D3","#FB8072")
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]

p <- ggplot(data=freq_data, aes(x=celltype, y=log2FC, fill = celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+ 
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,colour = "black"),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"),
        panel.grid = element_blank(),strip.background = NULL,
        strip.text = element_text(size = 15))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=celltype, y=log2FC,label=p.value1),size=6,position= position_dodge(0.6))
p1
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure5/',"figure5A",'_freq.pdf')
ggsave(plotfile, plot=p1, dpi = 1000, width = 10, height = 4)
dev.off()

# figure5 B features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
DefaultAssay(Cell.combined_B) <- "RNA"
Cell.combined_B <- SetIdent(object = Cell.combined_B, value = Cell.combined_B@meta.data$celltype)
gene <- c("FCN2","RUNX1","CD34","CIITA","CD40","IGHM","TCL1A","HES4","KLK1","S100A10","S100A4","ITGAX","HOPX","MZB1","JCHAIN","IRF4")
p1<-DotPlot(object = Cell.combined_B, features = gene, group.by  = 'celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
data$Cluster <- factor(data$Cluster,levels = rev(c("Pro_B","Pre_B","Naive BC","Memory BC","ABC","plasmablast","Plasma Cell")))
colors = c("#1A58A7","#4592E5","#FFFFFF","#eff3ff","#F76C6C","#F70000","#981216")
show_col(colors)
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  # scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  scale_color_gradientn(colors = c("grey","#eff3ff","#F76C6C","#981216"))+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black"),axis.text.x=element_text(angle = 90,hjust = 1,vjust = 1))
p
ggsave(paste0('/data1/zhuzn/wangsn/reviewer_4000/figure5/',"figure5B.pdf"),plot = p,width = 6,height = 3)

# figure5 C B细胞 naive分簇 --------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
# p <- DimPlot(Cell.combined_B, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE,raster=FALSE)
# p
# Cell.combined_B@meta.data$group <- rep("Other BC")
# colnames(Cell.combined_B@meta.data)
# Cell.combined_B@meta.data[which(Cell.combined_B@meta.data$integrated_snn_res.0.4 == "7"),20] <- "PDCD4low Naive"
# Cell.combined_B@meta.data[which(Cell.combined_B@meta.data$integrated_snn_res.0.4 %in% c("0","12")),20] <- "PDCD4high Naive"
# Cell.combined_B@meta.data$group <- factor(Cell.combined_B@meta.data$group,levels = c("PDCD4high Naive","PDCD4low Naive","Other BC"))
# save(Cell.combined_B, file="/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")

colors <- c("#ABFF56","#C972E4","#D9D9D9")
p1 <- DimPlot(Cell.combined_B, reduction = "umap", group.by = "group",cols = colors,label = F)+
  theme(legend.position = 'none')+ggtitle("")+
  theme(legend.position = 'none',axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p1
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure5/","figure5C1.tiff"),plot = p1,width = 6,height = 6)

# figure5 C cluster7 vs cluster0,12 DEG (骨髓） ------------------------------------------------
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
DefaultAssay(Cell.combined_B) <- "RNA"
data.list_tissue <- SplitObject(Cell.combined_B, split.by = "Tissue")
tissue_data <- data.list_tissue[["Bone_marrow"]]
tissue_data <- SetIdent(object = tissue_data, value = tissue_data@meta.data$integrated_snn_res.0.4)
diff_gene <- FindMarkers(tissue_data, ident.1 = "7", ident.2 = c("0","12"), logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_BM_cluster7v0_findmark.txt",quote = FALSE,sep = "\t",row.names = F)

# figure5 C cluster7 old vs young DEG (骨髓）------------------------------------------------
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_B_celltype.RData")
DefaultAssay(Cell.combined_B) <- "RNA"
Cell.combined1 <- subset(Cell.combined_B, subset = (integrated_snn_res.0.4  %in% "7"))
Cell.combined1 <- SetIdent(object = Cell.combined1, value = Cell.combined1@meta.data$Age)
table(Cell.combined1$Age)
diff_gene <- FindMarkers(Cell.combined1, ident.1 = "old", ident.2 = "young", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt",quote = FALSE,sep = "\t",row.names = F)

# figure5 C cluster7  vs 其他细胞------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "RNA"
colnames(all_sample.combined@meta.data)
all_cell <- all_sample.combined@meta.data[,c(4,5,18)]
### 添加亚型信息
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
colnames(Cell.combined_B@meta.data)
cell <- Cell.combined_B@meta.data[,c(4,5,20)]
colnames(cell) <- c("Tissue","Age","celltype")
cell  %<>% {.$id<-rownames(.);.} 

all_cell %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell) <- c("Tissue","Age","celltype","id","Celltype")
all_cell$celltype <- as.character(all_cell$celltype)
all_cell$Celltype <- as.character(all_cell$Celltype)
all_cell[is.na(all_cell)] <- "NON"
all_cell$Celltype <- ifelse(all_cell$Celltype == "NON" ,all_cell$celltype,all_cell$Celltype)
unique(all_cell$Celltype)
all_sample.combined <- AddMetaData(all_sample.combined,all_cell$Celltype,col.name = "sub_celltype")
table(all_sample.combined@meta.data$celltype)
table(all_sample.combined@meta.data$sub_celltype)
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$sub_celltype)
save(all_sample.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figure5/DE.RData")

diff_gene <- FindMarkers(all_sample.combined, ident.1 = "PDCD4low Naive",logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_vs_others.txt",quote = FALSE,sep = "\t",row.names = F)

# figure5 C B细胞 cluster7 vs 其他细胞 cluster7 old vs young DEG --------------------------------------------------------------
rm(list=ls())
DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_vs_others.txt", header=T)
DEG1 <- DEG1[DEG1$change_V1 == "UP",]
# DEG1 <- DEG1[DEG1$change_V1 == "DOWN",]
DEG1 <- DEG1[,c(1,2,6,7)];colnames(DEG1) <- c("p_val_7vother","avg_log2FC_7vother","gene","change_7vother")
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt", header=T)
DEG2 <- DEG2[DEG2$change_V1 == "UP",]
# DEG2 <- DEG2[DEG2$change_V1 == "DOWN",]
DEG2 <- DEG2[,c(1,2,6,7)];colnames(DEG2) <- c("p_val_7old_young","avg_log2FC_7old_young","gene","change_7old_young")
intersect(DEG1$gene,DEG2$gene)
DEG <- merge(DEG1,DEG2,all=F)
DEG <- DEG[!grepl("ENSMMUG", DEG$gene),]
write.table(DEG,file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/PDCD4_other_overlap_UP.txt",quote = FALSE,sep = "\t",row.names = F)

# figure5 C B细胞 cluster7 vs cluster0 cluster7 old vs young overlap --------------------------------------------------------------
rm(list=ls())
DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_BM_cluster7v0_findmark.txt", header=T)
# DEG1 <- DEG1[DEG1$change_V1 != "NEUTRAL",]
DEG1 <- DEG1[DEG1$change_V1 == "DOWN",]
DEG1 <- DEG1[,c(1,2,6,7)];colnames(DEG1) <- c("p_val_7v0","avg_log2FC_7v0","gene","change_7v0")
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt", header=T)
# DEG2 <- DEG2[DEG2$change_V1 != "NEUTRAL",]
DEG2 <- DEG2[DEG2$change_V1 == "DOWN",]
DEG2 <- DEG2[,c(1,2,6,7)];colnames(DEG2) <- c("p_val_7old_young","avg_log2FC_7old_young","gene","change_7old_young")
intersect(DEG1$gene,DEG2$gene)
DEG <- merge(DEG1,DEG2,all=F)
write.table(DEG,file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_overlap_DOWN.txt",quote = FALSE,sep = "\t",row.names = F)

rm(list=ls())
library(VennDiagram)
DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_BM_cluster7v0_findmark.txt", header=T) %>% .[.$change_V1 %in% c("UP","DOWN"),]
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt", header=T) %>% .[.$change_V1 %in% c("UP","DOWN"),]
intersect(DEG1$gene,DEG2$gene)
venn<- venn.diagram(list(DEG1=DEG1$gene,
                         DEG2=DEG2$gene),
                    filename=NULL,fill = c("cornflowerblue", "green"),
                    #col = "black",
                    col = "transparent", 
                    alpha = 0.4, cat.cex = 1.5,rotation.degree = 0)
pdf("/data1/zhuzn/wangsn/reviewer_4000/figure5/figure5C2.pdf",height = 5,width = 5)
grid.draw(venn)
dev.off()

# figure5 D overlap_小提琴图------------------------------------------------
rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/Cell.combined_B_celltype.RData")
DefaultAssay(Cell.combined_B) <- "RNA"
Cell.combined_B <- subset(Cell.combined_B, subset = (Tissue == "Bone_marrow" & celltype == "Naive BC"))
unique(Cell.combined_B$group)
DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_BM_cluster7v0_findmark.txt", header=T) %>% .[.$change_V1 == "DOWN",]
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt", header=T) %>% .[.$change_V1 == "DOWN",]

genes <- intersect(DEG1$gene,DEG2$gene)
genes <- genes[!grepl("ENSMMUG", genes)]
genes <- c("PDCD4","LDHA","CAST","EAF2")
# genes <- c("PDCD4","EAF2")
tsv <- data.frame(GetAssayData(Cell.combined_B,slot = 'data'),check.names=FALSE)
data1 <- data.frame(t(tsv[rownames(tsv) %in% genes,]))
data1$cell <- rownames(data1)
colnames(Cell.combined_B@meta.data)
metadata <- Cell.combined_B@meta.data[,c(5,20)]
metadata$cell <- rownames(metadata)
data <- merge(data1,metadata,by="cell")
data$Age <- factor(data$Age,levels = c("young","old"))
library(devtools)
library(ggunchained)
plot_list = list()
for (i in genes) {
  i <- gsub("-",".",i)
  data_1 <- data[,c(i,"Age","group")]
  colnames(data_1) <- c("expression","group","celltype")
  comparisons <- list(c("PDCD4low Naive","PDCD4high Naive"))
  p = data_1 %>% ggplot(aes(x=celltype,y=expression,fill=group)) +#修改x轴坐标   
    geom_split_violin(trim=F,color="#C972E4") + #绘制分半的小提琴图
    # geom_signif(comparisons = comparisons,map_signif_level = T, textsize = 6, test = wilcox.test, step_increase = 0.2) +
    # geom_boxplot(width=0.2,position = position_dodge(0.3))+#添加箱图
    scale_fill_manual(values = c('#39A8E3',"#D75156"))+ #设置填充的颜色
    scale_y_continuous(limits = c(-1,4))+
    theme_bw()+ #背景变为白色
    theme(axis.text.x = element_text(size = 15, angle = 90,hjust = 1,vjust = 0.5),
          panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
          panel.grid.major = element_blank(),   #不显示网格线（去掉横网格线）
          panel.grid.minor = element_blank())+  #不显示网格线（去掉竖网格线）
    labs(title=i)+
    # ylim(-1,4)+
    ylab("expression")+xlab("")#设置x轴和y轴的标题
  plot_list[[i]] = p
}
library(cowplot)
plot_grid(plotlist = plot_list,align = "h",ncol = 4)
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure5/","figure5D.pdf"),width = 18,height = 7)

# figure5 E cluster0 vs cluster8 cluster8 old vs young overlap_富集 --------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
# DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_BM_cluster7v0_findmark.txt", header=T) %>% .[.$change_V1 == "DOWN",]
# DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt", header=T) %>% .[.$change_V1 == "DOWN",]

DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_BM_cluster7v0_findmark.txt", header=T) %>% .[.$change_V1 == "UP",]
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/reviewer_4000/figure5/B_Cell_cluster7_old_young.txt", header=T) %>% .[.$change_V1 == "UP",]

DEG <- data.frame(intersect(DEG1$gene,DEG2$gene))
colnames(DEG) = "SYMBOL"
DEG <- merge(x=DEG,y=geneinfo,by="SYMBOL")
gene <- unique(DEG[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_result<-as.data.frame(ego_BP@result)#8107
ego_BP_result<-ego_BP_result[ego_BP_result$pvalue < 0.05,]
ego_BP_result$group <- rep("DOWN")
# write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure5/","B_DOWN_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure5/","B_UP_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

library("ggplot2")
GO <- c("GO:0048006","GO:0035442","GO:0048007","GO:0033632","GO:0019883")
GO <- c("GO:0002250","GO:0002504","GO:0002520","GO:0019724","GO:0030890","GO:0050871")
GO  <- ego_BP_result[ego_BP_result$ID %in% GO,]

library(DOSE)
GO$GeneRatio1 <- parse_ratio(GO$GeneRatio)
p <- ggplot(GO,aes(reorder(Description,GeneRatio1),GeneRatio1,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_colour_gradient(low="#3E1AFD",high="#C1B5FC")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
  geom_point(size = 2.0,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure5/","figure5E.pdf"),plot = p,width = 5.5,height = 5)

# figure3D overlap TF targets FC 富集结果图（按上下调分开） --------------------------------------------------
rm(list=ls())
BHLHE40_GO_BP_UP <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure5/B_UP_GO.txt")
BHLHE40_GO_BP_UP$Group <- rep("UP")
BHLHE40_GO_BP_DOWN <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure5/B_DOWN_GO.txt")
BHLHE40_GO_BP_DOWN$Group <- rep("DOWN")
GO <- rbind(BHLHE40_GO_BP_UP,BHLHE40_GO_BP_DOWN)
GO_lists <- c("GO:2000405","GO:0007162","GO:0002683","GO:0032682","GO:0006954",
            "GO:0002250","GO:0002504","GO:0002520","GO:0019724","GO:0030890","GO:0050871")
dat <- GO[GO$ID %in% GO_lists,]

dat1 <- dat[dat$Group == "UP",] %>% arrange(-Count)
dat1$Description <- factor(dat1$Description, levels = dat1$Description)
p1 <- ggplot(dat1,aes(x=Count, y=Description)) +
  geom_bar(fill="#C41013", alpha = 0.6,stat = "identity") + 
  theme_classic() + labs(x="",y="",title = "Upregulated Targets GO") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = 'black',size = 12)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_text(aes(x=0.1,y=Description,label=Description), hjust = 0)
p1

dat2 <- dat[dat$Group == "DOWN",] %>% arrange(Count)
dat2$Description <- factor(dat2$Description, levels = dat2$Description)
p2 <- ggplot(dat2,aes(x=Count, y=Description)) +
  geom_bar(fill="#3076AD", alpha = 0.8,stat = "identity") + 
  theme_classic() + labs(x="",y="",title = "Downregulated Targets GO") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = 'black',size = 12)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_text(aes(x=0.1,y=Description,label=Description), hjust = 0)
p2
p <- p1/p2
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure5figure5F.pdf",plot = p,width = 6,height = 6)

# figure3D overlap TF targets FC 富集结果图（按上下调分开） --------------------------------------------------
rm(list=ls())
ego_BP_result_UP <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure5/B_UP_GO.txt")
ego_BP_result_UP$Group <- rep("UP")
c <- c("GO:2000405","GO:0007162","GO:0002683","GO:0032682","GO:0006954")
top_go_up <- ego_BP_result_UP[ego_BP_result_UP$ID %in% c,]
top_go_up <- top_go_up[order(top_go_up$Count,decreasing = T),]
top_go_up$Number <- c(5:1)

ego_BP_result_DOWN <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure5/B_DOWN_GO.txt")
ego_BP_result_DOWN$Group <- rep("DOWN")
c <- c("GO:0002250","GO:0002504","GO:0002520","GO:0019724","GO:0030890","GO:0050871")
top_go_down <- ego_BP_result_DOWN[ego_BP_result_DOWN$ID %in% c,]
top_go_down <- top_go_down[order(top_go_down$Count,decreasing = F),]
top_go_down$Number <- c(-1:-6)
top_go <- rbind(top_go_up,top_go_down)
top_go$number = factor(rev(1:nrow(top_go)))
top_go$Count = as.numeric(ifelse(top_go$Group == 'DOWN' ,paste0("-",top_go$Count),top_go$Count))

p <- ggplot(data=top_go, aes(x=number, y=Count,fill = Number)) + geom_bar(stat="identity", width=0.8) + 
  scale_x_discrete(labels=rev(top_go$Description)) +
  # scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
  scale_fill_gradient2(low="#52B088",mid="#eff3ff",high ="#8152E7",midpoint = 0)+
  coord_flip() +theme_bw() + xlab("") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme(axis.text=element_text(size=20,face = "plain", color="black"),panel.grid = element_blank(),
        axis.title=element_text(size=10),legend.text=element_text(size=10),legend.title = element_text(size=20))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure5/figure5F.pdf",plot = p,width = 14,height = 6)

# 实验统计图 -------------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
data <- data.frame(exp = c(75,73,70,48,51,50),Aged = c(rep("Young",3),rep("Aged",3)))
data$Aged <- factor(data$Aged,levels = c("Young","Aged"))
p <- ggplot(data=data, aes(x=Aged,y=exp,color=Aged,fill=Aged))+
  # geom_boxplot()+
  stat_summary(fun = mean,geom = "bar",aes(color=Aged,fill=Aged,width=0.5))+
  stat_summary(geom = "errorbar",fun.min=min,fun.max=max,width=0.2)+
  theme_bw()+
  scale_color_manual(values=c("blue","red"))+
  scale_fill_manual(values=c("white","white"))+
  theme(legend.title=element_blank())+labs(x="",y="PDCD4+/TCL1+ B Cell",caption="")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure5/figure.5H.pdf",plot = p,width = 5,height = 6)



