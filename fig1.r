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

# figure1  主图重新定义细胞类型 --------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
ncol(all_sample.combined)
table(all_sample.combined$Tissue)
unique(all_sample.combined@meta.data$celltype)
all_sample.combined@meta.data$celltype <- factor(all_sample.combined@meta.data$celltype,
                                                 levels = c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                            "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"))

colnames(all_sample.combined@meta.data)
table(all_sample.combined@meta.data$Tissue)
d <- data.frame(table(all_sample.combined@meta.data$celltype))
d$freq <- d$Freq/sum(d$Freq)*100
rownames(d) <- d$Var1
d <- d[c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"),]
d$number <- c(1:11)
colnames(d) <- c("celltype","Freq","freq","number")
celltype_cells <- all_sample.combined@meta.data[,c(4,5,18)]
celltype_cells <- join(celltype_cells,d) 
celltype_cells$celltype1 <- paste0(celltype_cells$number,":",celltype_cells$celltype,"(",round(celltype_cells$freq,2),"%)")
all_sample.combined <- AddMetaData(all_sample.combined,celltype_cells$celltype1,col.name = "celltype1")
all_sample.combined <- AddMetaData(all_sample.combined,celltype_cells$number,col.name = "number")
# save(all_sample.combined,file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")

all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$celltype1)
unique(all_sample.combined@meta.data$celltype1)
all_sample.combined@meta.data$celltype1 <- factor(all_sample.combined@meta.data$celltype1,
                                                  levels = c("1:T Cell(52.16%)","2:NK(2.98%)","3:B Cell(29.01%)",
                                                             "4:Monocyte(2.73%)","5:Macrophage(1.32%)","6:Dendritic cell(1.11%)",
                                                             "7:Fibroblast(0.25%)",
                                                             "8:Erythrocyte(0.99%)","9:Megakaryocyte(0.21%)",
                                                             "10:Granulocyte(5.46%)","11:HSC(3.77%)"))
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")

data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
young_data <- data.list_age[["young"]]
p1 <- DimPlot(young_data, reduction = "umap", group.by = "number",pt.size = 0.01,
              cols = colors, label = F, repel = TRUE,raster=F)+xlab("") + ylab("") + 
  ggtitle("")+
  theme(legend.position = 'none',axis.ticks = element_blank(),axis.line = element_blank(),
        axis.text = element_blank())
p1
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure1/","figure.1A_young.tiff"),plot = p1,dpi = 1000,width = 16,height = 16)

old_data <- data.list_age[["old"]]
p2 <- DimPlot(old_data, reduction = "umap", group.by = "number",pt.size = 0.01,
              cols = colors, label = F, repel = TRUE,raster=F)+xlab("") + ylab("") + 
  ggtitle("")+
  theme(legend.position = 'none',axis.ticks = element_blank(),axis.line = element_blank(),
        axis.text = element_blank())
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure1/","figure.1A_old.tiff"),plot = p2,dpi = 1000,width = 16,height = 16)

p <- DimPlot(all_sample.combined, reduction = "umap", group.by = "number",pt.size = 0.01, 
             cols = colors, label = T, repel = TRUE, label.size = 10,raster=F)+
  ggtitle("")+xlab("") + ylab("") + 
  theme(legend.position = 'none',axis.ticks = element_blank(),axis.line = element_blank(),
        axis.text = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure1/figure.1B.pdf",dpi = 1000,plot = p,width = 12,height = 12)

table(all_sample.combined@meta.data$Age)
all_sample.combined@meta.data$Age <- factor(all_sample.combined@meta.data$Age,levels = c("young","old"))
colors <- c('#39A8E3',"#D75156")
show_col(colors)
p <- DimPlot(all_sample.combined, reduction = "umap", group.by = "Age",pt.size = 0.05,
             cols = colors,
             raster=FALSE)+
  ggtitle("")+
  theme(legend.position = 'none',axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure1/figure.1B1.tiff",dpi = 1000,plot = p,width = 8,height = 8)

# figure1C  features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$celltype)
cell_markers <- read.delim("/data1/zhuzn/wangsn/V4/celltype/all_sample_celltype_marker2.txt")
unique(cell_markers$Cell.type)
gene <- as.character(unique(cell_markers$Gene))
p1<-DotPlot(object = all_sample.combined, features = gene, group.by  = 'celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
unique(cell_markers$Cell.type)
data$Cluster <- factor(data$Cluster,levels = rev(c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                   "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")))

p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradientn(colors = c("grey","grey","#eff3ff","#F76C6C","#981216"))+
  theme_classic()+
  labs(x = '', y = '')+
  theme(axis.text.x=element_text(size=20, color="black",angle = 90,vjust = 0.5,hjust = 1),legend.position = "top",
        axis.text.y=element_text(size=22, color="black"))# + RotatedAxis():45度
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure1/","figure.","1C.pdf"),p,dpi = 1000,width = 11,height = 6)

# figure1D  FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined$Tissue_Age <- paste0(all_sample.combined$Tissue,"_",all_sample.combined$Age)
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
p1 <- DimPlot(all_sample.combined, reduction = "umap", group.by = "celltype",pt.size = 0.01,
              split.by = "Tissue_Age",ncol = 4,
              cols = colors, label = F, repel = TRUE,raster=F)+xlab("") + ylab("") + 
  ggtitle("")+
  theme(legend.position = 'none',axis.ticks = element_blank(),axis.line = element_blank(),
        axis.text = element_blank())
p1

colnames(all_sample.combined@meta.data)
celltype_tissue <- all_sample.combined@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
cell_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
               "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
freq_data <- data.frame()
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  freq_data_tissue <- data.frame()
  for (j in cell_list){
    cell <- data.frame(table(tissue[which(tissue$celltype == j),]))
    
    if (nrow(cell) == 2){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = cell[which(cell$Age == "young"),"Freq"],
                               FC = as.character(cell[which(cell$Age == "old"),"Freq"]/cell[which(cell$Age == "young"),"Freq"]))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(nrow(cell) == 0){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                             old_cell = "0",
                             young_cell = "0",
                             FC = "0")
    freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "old"){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
                               old_cell = cell[which(cell$Age == "old"),"Freq"],
                               young_cell = "0",
                               FC = as.character(cell[which(cell$Age == "old"),"Freq"]))
      freq_data_tissue <- rbind(freq_data_tissue,freq_data1)
    }else if(cell[,2] == "young"){
      freq_data1 <- data.frame(Tissue = i,celltype = j,
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
  for (celltype in cell_list){
    test_data1 <- test_data[test_data$celltype == celltype,c(3,4)]
    test_data2 <- test_data[test_data$celltype != celltype,]
    test_data3 <- data.frame(old_cell = sum(as.numeric(test_data2$old_cell)),young_cell = sum(as.numeric(test_data2$young_cell)))
    test_data4 <- rbind(test_data1,test_data3)
    test_data4$old_cell <- as.numeric(test_data4$old_cell);test_data4$young_cell <- as.numeric(test_data4$young_cell)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,celltype = celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))

library("ggplot2")
unique(freq_data$celltype)
freq_data$celltype <- factor(freq_data$celltype,levels = c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                           "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")

show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
freq_data$log2FC[!is.finite(freq_data$log2FC)] <- 0
p <- ggplot(data=freq_data, aes(x=celltype, y=log2FC, fill = celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+
  scale_y_continuous(limits=c(-4,4))+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,colour = "black"),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"),
        panel.grid = element_blank(),strip.background = NULL,
        strip.text = element_text(size = 15))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=celltype, y=log2FC,label=fisher_data$p.value1),size=6,position= position_dodge(0.6))
p1
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1D_FC.pdf')
ggsave(plotfile, plot=p1, dpi = 1000, width = 8, height = 4)
dev.off()

# figure1D  占比饼图 ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
colnames(all_sample.combined@meta.data)
celltype_tissue <- all_sample.combined@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
cell_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
               "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
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
freq_data3 <- freq_data2[,c(1,2,14,15,13,12,3,11,9,4,6,5,10,7,8,16)]
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
p <- ggplot()+
  geom_scatterpie(data=freq_data2,
                  aes(x,y,group=region,r=0.9),
                  cols = cell_list)+
  coord_equal()+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values = colors)
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure1/figure1D_FC1.pdf",plot = p,width = 15,height = 7)
dev.off()

# figure1D  percent_FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
# 设置scipen参数，禁止科学计数法显示小数
options(scipen = 999)
colnames(all_sample.combined@meta.data)
celltype_tissue <- all_sample.combined@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
cell_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
               "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
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
# freq_data$log2FC <- log2(as.numeric(freq_data$FC))

###fisher检验
fisher_data <- data.frame()
for (tissue in unique(celltype_tissue$Tissue)){
  test_data <- freq_data[freq_data$Tissue == tissue,]
  for (celltype in cell_list){
    test_data1 <- test_data[test_data$celltype == celltype,c(3,4)]
    test_data2 <- test_data[test_data$celltype != celltype,]
    test_data3 <- data.frame(old_cell = sum(as.numeric(test_data2$old_cell)),young_cell = sum(as.numeric(test_data2$young_cell)))
    test_data4 <- rbind(test_data1,test_data3)
    test_data4$old_cell <- as.numeric(test_data4$old_cell);test_data4$young_cell <- as.numeric(test_data4$young_cell)
    rownames(test_data4) <- c("cell","others")
    fisher_data1 <- data.frame(Tissue = tissue,celltype = celltype,p.value = fisher.test(test_data4)$p.value)
    fisher_data <- rbind(fisher_data,fisher_data1)
  }
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))

library("ggplot2")
unique(freq_data$celltype)
freq_data$celltype <- factor(freq_data$celltype,levels = c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                           "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")

show_col(colors)
freq_data <- cbind(freq_data,fisher_data)
freq_data <- freq_data[,-c(9,10)]
freq_data$log2FC[!is.finite(freq_data$log2FC)] <- 0
p <- ggplot(data=freq_data, aes(x=celltype, y=log2FC, fill = celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.7)) +
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "",x = '', y = 'log2FC')+
  scale_y_continuous(breaks = seq(-0.26, 0.22, by = 0.02))+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.text.y = element_text(size = 18,colour = "black"),legend.position = "top",
        legend.text = element_text(size = 10),legend.key.size = unit(10,"pt"),
        panel.grid = element_blank(),strip.background = NULL,
        strip.text = element_text(size = 15))+
  xlab("") + ylab("")
p
p1<-p+geom_text(aes(x=celltype, y=log2FC,label=fisher_data$p.value1),size=4,position= position_dodge(0.6))
p1
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1D_FC.pdf')
ggsave(plotfile, plot=p1, dpi = 1000, width = 8, height = 8)
dev.off()

library(gg.gap)
p2 = gg.gap(plot = p,
            segments = c(0.08,0.2),
            tick_width = 0.02,
            rel_heights = c(0.25,0,0.1),# 设置分隔为的三个部分的宽度
            ylim = c(-0.26, 0.22))
p2
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1D_+.pdf')
ggsave(plotfile, plot=p2, dpi = 1000, width = 8, height = 8)
dev.off()
p3 = gg.gap(plot = p,
            segments = c(-0.2,-0.1),
            tick_width = 0.02,
            rel_heights = c(0.25,0,0.1),# 设置分隔为的三个部分的宽度
            ylim = c(-0.26, 0.22))
p3
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1D_-.pdf')
ggsave(plotfile, plot=p3, dpi = 1000, width = 8, height = 8)
dev.off()

# figure1D  细胞数量(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
colnames(all_sample.combined@meta.data)
celltype_tissue <- all_sample.combined@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
unique(celltype_tissue$celltype)
cell_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast",
               "Erythrocyte","Megakaryocyte","Granulocyte","HSC")
cell_data <- data.frame(Tissue = "Tissue",Age = "Age",celltype = "celltype",Freq = "Freq")
for (i in unique(celltype_tissue$Tissue)){
  tissue <- celltype_tissue[which(celltype_tissue$Tissue == i),]
  for (j in cell_list){
    cell_data1 <- data.frame(table(tissue[which(tissue$celltype == j),]))
    cell_data <- rbind(cell_data,cell_data1)
  }
}
cell_data <- cell_data[-1,]

library("ggplot2")
unique(cell_data$celltype)
show_col(colors)
cell_data$Freq <- as.numeric(cell_data$Freq)
unique(cell_data$celltype)

cell_data$Celltype <- paste0(cell_data$celltype,"_",cell_data$Age)
cell_data$Celltype <- factor(cell_data$Celltype,levels = c("T Cell_young","T Cell_old","NK_young","NK_old","B Cell_young","B Cell_old",
                                                           "Monocyte_young","Monocyte_old","Macrophage_young","Macrophage_old",
                                                           "Dendritic cell_young","Dendritic cell_old",
                                                           "Fibroblast_young","Fibroblast_old","Erythrocyte_young","Erythrocyte_old",
                                                           "Megakaryocyte_young","Megakaryocyte_old","Granulocyte_young","Granulocyte_old","HSC_young","HSC_old"))
# colors <- c("#6699CC","#6699CC","#1f78b4","#1f78b4","#33a02c","#33a02c","#ff7f00","#ff7f00","#b2df8a","#b2df8a","#beaed4","#beaed4",
#             "#fb9a99","#fb9a99","#e78ac3","#e78ac3","#e6ab02","#e6ab02","#ef3b2c","#ef3b2c","#21888F","#21888F","#83639D","#83639D")
colors <- c("#99B3CC","#6699CC","#80B5DA","#1f78b4","#90D98A","#33a02c","#FFB266","#ff7f00","#C9DFB5","#b2df8a","#C5B9D4","#beaed4",
            "#FBC1C0","#fb9a99","#E7C2D9","#e78ac3","#E6C35D","#e6ab02","#EF837A","#ef3b2c","#88D9DD","#21888F","#90809D","#83639D")
show_col(colors)
#ggplot2 作图
library(gg.gap)
p <- ggplot(data=cell_data, aes(x=Celltype, y=Freq, fill = Celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.8)) +
  scale_fill_manual(values = colors, guide = "none")+ 
  theme_bw() + 
  # scale_y_log10()+
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0, angle = 270))+
  #geom_text(aes(label=Freq),position = position_stack(vjust = 1.0),size = 1.5)+
  xlab("") + ylab("cell numbers")
p
p2 = gg.gap(plot = p,
            segments = c(7500,10500),
            tick_width = 500,
            rel_heights = c(0.25,0,0.1),# 设置分隔为的三个部分的宽度
            ylim = c(0,23000))
p2
# ggsave(paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1D_FC.pdf'), plot=p1, dpi = 1000, width = 8, height = 4)
dev.off()

# figure1D  细胞数量百分比(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
colnames(all_sample.combined@meta.data)
celltype_tissue <- all_sample.combined@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
unique(celltype_tissue$celltype)
cell_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast",
               "Erythrocyte","Megakaryocyte","Granulocyte","HSC")
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
freq_data <- freq_data[,c(1,2,5,6)]
colnames(freq_data) <- c("Tissue","celltype","old","young")
cell_data<-melt(freq_data,
                id.vars = c("Tissue","celltype"))
colnames(cell_data) <- c("Tissue","celltype","Age","Freq")
library("ggplot2")
unique(cell_data$celltype)
show_col(colors)
cell_data$Freq <- as.numeric(cell_data$Freq)
unique(cell_data$celltype)

cell_data$Celltype <- paste0(cell_data$celltype,"_",cell_data$Age)
cell_data$Celltype <- factor(cell_data$Celltype,levels = c("T Cell_young","T Cell_old","NK_young","NK_old","B Cell_young","B Cell_old",
                                                           "Monocyte_young","Monocyte_old","Macrophage_young","Macrophage_old",
                                                           "Dendritic cell_young","Dendritic cell_old",
                                                           "Fibroblast_young","Fibroblast_old","Erythrocyte_young","Erythrocyte_old",
                                                           "Megakaryocyte_young","Megakaryocyte_old","Granulocyte_young","Granulocyte_old","HSC_young","HSC_old"))
# colors <- c("#6699CC","#6699CC","#1f78b4","#1f78b4","#33a02c","#33a02c","#ff7f00","#ff7f00","#b2df8a","#b2df8a","#beaed4","#beaed4",
#             "#fb9a99","#fb9a99","#e78ac3","#e78ac3","#e6ab02","#e6ab02","#ef3b2c","#ef3b2c","#21888F","#21888F","#83639D","#83639D")
colors <- c("#99B3CC","#6699CC","#80B5DA","#1f78b4","#90D98A","#33a02c","#FFB266","#ff7f00","#C9DFB5","#b2df8a","#C5B9D4","#beaed4",
            "#FBC1C0","#fb9a99","#E7C2D9","#e78ac3","#E6C35D","#e6ab02","#EF837A","#ef3b2c","#88D9DD","#21888F","#90809D","#83639D")
show_col(colors)
#ggplot2 作图
library(gg.gap)
p <- ggplot(data=cell_data, aes(x=Celltype, y=Freq, fill = Celltype, width=0.8))+
  facet_grid(.~Tissue,scales= "free")+
  geom_bar(stat="identity",position=position_dodge(0.8)) +
  scale_fill_manual(values = colors, guide = "none")+ 
  theme_bw() + 
  # scale_y_log10()+
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0, angle = 270))+
  #geom_text(aes(label=Freq),position = position_stack(vjust = 1.0),size = 1.5)+
  xlab("") + ylab("cell percent")
p
p2 = gg.gap(plot = p,
            segments = c(7500,10500),
            tick_width = 500,
            rel_heights = c(0.25,0,0.1),# 设置分隔为的三个部分的宽度
            ylim = c(0,23000))
p2
# ggsave(paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1D_FC.pdf'), plot=p1, dpi = 1000, width = 8, height = 4)
dev.off()

# figure1E  SASP ----------------------------------------------------------------------
rm(list=ls())
library(ggridges)
library(ggpubr)
setwd('/data1/zhuzn/wangsn/SASP/')
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined1 <- all_sample.combined
unique(all_sample.combined1@meta.data$celltype)
data.list_celltype <- SplitObject(all_sample.combined1, split.by = "celltype")
sample_list <- unique(all_sample.combined@meta.data$celltype)
Sasp <- data.frame()
Sasp2 <- data.frame()
Sasp3 <- data.frame()
for (celltype in sample_list){
  all_sample.combined <- data.list_celltype[[celltype]]
  Idents(all_sample.combined)
  DefaultAssay(all_sample.combined) <- "RNA"
  DefaultAssay(all_sample.combined)
  sample.integrated <- all_sample.combined
  SASP_gene_set <- read.table("/data1/zhuzn/wangsn/SASP/SASP_V2.list", quote="\"", comment.char="")
  gene <- as.list(SASP_gene_set)
  sample.integrated <- AddModuleScore(object = sample.integrated,features = gene,ctrl = 100, name = 'SASP_Features',replace = T)
  
  names(sample.integrated@meta.data)
  sasp <- sample.integrated@meta.data[,c(5,19)]
  
  sasp2 <-aggregate(sasp$SASP_Features1,list(sasp$Age),mean)
  names(sasp) <- c('Group','Gene_set_score')
  sasp$celltype <- rep(celltype)
  names(sasp2) <- c('Group','mean')
  sasp2$celltype <- rep(celltype)
  # sasp$Group <- factor(sasp$Group,levels = c("young","old"))
  sasp$Group <- factor(sasp$Group,levels = c("old","young"))
  sasp3 <- compare_means(Gene_set_score~Group, data=sasp,method = "t.test")
  sasp3$celltype <- rep(celltype)
  
  Sasp <- rbind(Sasp,sasp)
  Sasp2 <- rbind(Sasp2,sasp2)
  Sasp3 <- rbind(Sasp3,sasp3)
}
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
colors <- c("#7756A1","black","#BEBADA","black","#1084D3","black","#73C508","black","#CFE6A7","black","#144F2B","black",
            "#FF348B","black","#B7D5F3","black","#A91C51","black","#E9F033","black","#FFAF18","black")

show_col(colors)
library(dplyr)
library(forcats)
Sasp1 <- Sasp
Sasp1$Celltype <- paste0(Sasp1$celltype," ",Sasp1$Group)
Sasp1$celltype <- factor(Sasp1$celltype,levels = rev(c("T Cell","NK","B Cell","Monocyte","Macrophage",
                                                       "Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")))
Sasp1$Celltype <- factor(Sasp1$Celltype,levels = c("T Cell old","T Cell young","NK old","NK young",
                                                   "B Cell old","B Cell young",
                                                   "Monocyte old","Monocyte young","Macrophage old","Macrophage young",
                                                   "Dendritic cell old","Dendritic cell young",
                                                   "Fibroblast old","Fibroblast young",
                                                   "Erythrocyte old","Erythrocyte young",
                                                   "Megakaryocyte old","Megakaryocyte young",
                                                   "Granulocyte old","Granulocyte young","HSC old","HSC young"))
p <- Sasp1 %>%
  ggplot(aes(y = celltype)) +
  geom_density_ridges(aes(x = Gene_set_score, color = Celltype),fill = "white") +
  labs(x = "Gene_set_score",y = "",title = "celltype",subtitle = "",caption = "") +
  theme_ridges()+
  scale_color_cyclical(values = colors,
                       name = "Age", guide = "legend") +
  theme(axis.text = element_text(size = 20),
        legend.position = "none",
        legend.text = element_text(size = 18),legend.key.size = unit(20,"pt"),legend.key = element_rect(colour='grey32'))
p
plotfile = paste0('/data1/zhuzn/wangsn/reviewer_4000/figure1/',"figure.",'1E_SASP.pdf')
ggsave(plotfile, plot=p, dpi = 1000, width = 7.5, height = 7)
