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

# figure1  主图重新定义细胞类型 --------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
unique(all_sample.combined@meta.data$celltype1)
all_sample.combined@meta.data$celltype <- factor(all_sample.combined@meta.data$celltype,
                                                 levels = c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                            "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"))
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$celltype1)
all_sample.combined@meta.data$celltype1 <- factor(all_sample.combined@meta.data$celltype1,
                                                  levels = c("1:T Cell(53.25%)","2:NK(2.93%)","3:B Cell(27.94%)",
                                                             "4:Monocyte(2.51%)","5:Macrophage(1.15%)","6:Dendritic cell(0.89%)","7:Fibroblast(0.49%)",
                                                             "8:Erythrocyte(0.95%)","9:Megakaryocyte(0.55%)",
                                                             "10:Granulocyte(6.76%)","11:HSC(2.58%)"))
library(RColorBrewer)
display.brewer.all()
colors <- brewer.pal(n = 10, name = "Set3")
colors
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
p <- DimPlot(all_sample.combined, reduction = "umap", group.by = "number",pt.size = 0.01, split.by = "Age",
             cols = colors, label = T, repel = TRUE, label.size = 10,raster=F)+
  theme(legend.position = "right")
p
data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
young_data <- data.list_age[["young"]]
p1 <- DimPlot(young_data, reduction = "umap", group.by = "number",pt.size = 0.01,
              cols = colors, label = F, repel = TRUE,raster=F)+xlab("") + ylab("") + 
  ggtitle("")+
  theme(legend.position = 'none',axis.ticks = element_blank(),axis.line = element_blank(),
        axis.text = element_blank())
p1
ggsave(paste0("/data1/zhuzn/wangsn/reviewer/figure1/","figure.1A_young.jpeg"),plot = p1,dpi = 1000,width = 16,height = 16)

old_data <- data.list_age[["old"]]
p2 <- DimPlot(old_data, reduction = "umap", group.by = "number",pt.size = 0.01,
              cols = colors, label = F, repel = TRUE,raster=F)+xlab("") + ylab("") + 
  ggtitle("")+
  theme(legend.position = 'none',axis.ticks = element_blank(),axis.line = element_blank(),
        axis.text = element_blank())
ggsave(paste0("/data1/zhuzn/wangsn/reviewer/figure1/","figure.1A_old.jpeg"),plot = p2,dpi = 1000,width = 16,height = 16)

p <- DimPlot(all_sample.combined, reduction = "umap", group.by = "number",pt.size = 0.01, 
             cols = colors, label = T, repel = TRUE, label.size = 10,raster=F)+
  theme(legend.position = "right")
p
ggsave("/data1/zhuzn/wangsn/reviewer/figure1/figure.1B.pdf",dpi = 1000,plot = p,width = 12,height = 12)

table(all_sample.combined@meta.data$Age)
all_sample.combined@meta.data$Age <- factor(all_sample.combined@meta.data$Age,levels = c("young","old"))
colors <- c('#39A8E3',"#D75156")
show_col(colors)
p <- DimPlot(all_sample.combined, reduction = "umap", group.by = "Age",pt.size = 0.05,
             cols = colors,
             raster=FALSE)+
  ggtitle("")+
  theme(axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer/figure1/figure.1B1.pdf",dpi = 1000,plot = p,width = 8,height = 7.5)

# figure1  features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
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
ggsave(paste0("/data1/zhuzn/wangsn/reviewer/figure1/","figure.","1C.pdf"),p,dpi = 1000,width = 11,height = 6)

# figure1  FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
colnames(all_sample.combined@meta.data)
celltype_tissue <- all_sample.combined@meta.data[,c(4,5,23)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
cell_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
               "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
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
freq_data$celltype <- factor(freq_data$celltype,levels = c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                           "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC"))
freq_data$log2FC <- as.numeric(as.character(freq_data$log2FC))
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")

show_col(colors)
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
plotfile = paste0('/data1/zhuzn/wangsn/reviewer/figure1/',"figure.",'1D_FC.pdf')
ggsave(plotfile, plot=p1, dpi = 1000, width = 10, height = 4)
dev.off()

# figure1  SASP ----------------------------------------------------------------------
rm(list=ls())
library(ggridges)
library(ggpubr)
setwd('/data1/zhuzn/wangsn/SASP/')
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
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
  sasp <- sample.integrated@meta.data[,c(5,26)]
  
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
        # legend.position = "none",
        legend.text = element_text(size = 18),legend.key.size = unit(20,"pt"),legend.key = element_rect(colour='grey32'))
p
plotfile = paste0('/data1/zhuzn/wangsn/reviewer/figure1/',"figure.",'1_SASP.pdf')
ggsave(plotfile, plot=p, dpi = 1000, width = 7.5, height = 7)

# figure2  DEGs_tissue_celltype ----------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
ncol(all_sample.combined)
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$Age)
data.list_tissue <- SplitObject(all_sample.combined, split.by = "Tissue")

for (tissue in unique(all_sample.combined@meta.data$Tissue)){
  tissue_data <- data.list_tissue[[tissue]]
  data.list_celltype <- SplitObject(tissue_data, split.by = "celltype")
  for (celltype in unique(all_sample.combined@meta.data$celltype)){
    celltype_data <- data.list_celltype[[celltype]]
    d <- data.frame(table(celltype_data$Age))
    if(nrow(d) == 1){
      print(paste0(tissue,"_",celltype," no young or old cells"))
    }
    else{
      a = d[1,2];b = d[2,2]
      if(a<3 | b <3){
        print(paste0(tissue,"_",celltype,"Cells fewer than 3 cells")) 
      }else{
        diff_gene <- FindMarkers(celltype_data,ident.1 = "old", ident.2 = "young",logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
        diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
        diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                               ifelse(diff_gene$avg_log2FC > 0.5 ,'UP',
                                                      ifelse(diff_gene$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                               'NOT'))
        celltype <- gsub(" ","_",celltype)
        write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
      } 
    }
  }
}

# figure2 A  DEGs数量 --------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
tissue_list <- c("Lymph_node","PBMC","Spleen","Bone_marrow")
celltype_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                   "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
DEG_data <- data.frame(p_val = "p_val",gene = "gene",change_V1 = "change_V1",celltype = "celltype")
for (tissue in tissue_list){
  for (celltype in celltype_list){
    celltype <- gsub(" ","_",celltype)
    if (file.exists(paste0("./DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"))){
      dif_data <- read.delim(paste0("./DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
      dif_data <- dif_data[which(dif_data$change_V1 != "NOT"),]
      DEG_data1 <- dif_data[,c(1,6,7)]
      DEG_data1$celltype <- rep(paste0(tissue,"_",celltype))
      DEG_data <- rbind(DEG_data,DEG_data1)
    }
    else{
      next
    }
  }
}
DEG_data <- DEG_data[-1,]
DEG_data$sig_data <- "sig"
class(DEG_data$celltype)

for (i in unique(DEG_data$gene)){
  DEG_data[which(DEG_data$gene == i),5] <- nrow(DEG_data[DEG_data$gene == i,])
}
range(DEG_data$sig_data)
D <- data.frame(table(DEG_data$sig_data))
D$Var1 <- factor(D$Var1,levels = c(1:30))

p <- ggplot(data = D, mapping = aes(x = Var1, y = Freq)) +
  geom_bar(stat = 'identity', color = "black",fill = 'lightblue',position="dodge",width =0.6)+theme_bw()+
  coord_flip()+
  labs(x = 'Unique Cell Identities', y = 'DiffEx Genes')+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(size = 1,colour = "black"),
        axis.text=element_text(size=20,color="black"),axis.title=element_text(size=20,color="black"))
p

library(gg.gap)
p <- ggplot(D, aes(Var1, Freq,group = 1))+
  geom_line(size=1,color='#3398CC')+geom_point(size=10,color='#D1484F')+
  labs(x = 'Unique Cell Identities', y = 'DiffEx Genes')+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(size = 1,colour = "black"),
        axis.text=element_text(size=20,color="black"),axis.title=element_text(size=20,color="black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer/figure2/figure.2A.pdf",plot = p,dpi=1000,width = 15,height = 10)

# figure2 A 50% 差异基因 ----------------------------------------------------------------
rm(list=ls())
library("ggplot2")
tissue_list <- c("Lymph_node","PBMC","Spleen","Bone_marrow")
celltype_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                   "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
DEG_data <- data.frame(gene = "gene",change_V1 = "change_V1",celltype = "celltype")

for (tissue in tissue_list){
  for (celltype in celltype_list){
    celltype <- gsub(" ","_",celltype)
    if (file.exists(paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"))){
      dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
      dif_data$change_V1 = as.factor(ifelse(dif_data$p_val < 0.05,
                                             ifelse(dif_data$avg_log2FC > 0.25 ,'UP',
                                                    ifelse(dif_data$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                             'NOT'))
      dif_data <- dif_data[which(dif_data$change_V1 != "NOT"),]
      DEG_data1 <- dif_data[,c(6,7)]
      DEG_data1$celltype <- rep(paste0(tissue,"_",celltype))
      DEG_data <- rbind(DEG_data,DEG_data1)
    }
    else{
      next
    }
  }
}
DEG_data <- DEG_data[-1,]
DEG_data <- DEG_data[DEG_data$change_V1 != "NOT",]
DEG_data$sig_data <- "sig"
class(DEG_data$celltype)
for (i in unique(DEG_data$gene)){
  DEG_data$sig_data[which(DEG_data$gene == i)] <- nrow(DEG_data[DEG_data$gene == i,])
}
change_data <- DEG_data[,c(1,2,4)]
change_data <- unique(change_data)
change_data$sig_data <- as.numeric(change_data$sig_data)
change_data$change_percent <- (change_data$sig_data/30)*100
gene_list <- unique(change_data$gene)

tissue = "Lymph_node"
celltype = "T_Cell"
dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
gene_FC <- data.frame(p_val = "p_val",avg_log2FC = "avg_log2FC",gene = "gene",change_V1 = "change_V1",celltype = "celltype")

for (tissue in tissue_list){
  for (celltype in celltype_list){
    celltype <- gsub(" ","_",celltype)
    if (file.exists(paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"))){
      dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
      for (gene in gene_list){
        if(gene %in%  dif_data$gene) {
          gene_dif1 <- dif_data[dif_data$gene == gene,c(1,2,6,7)]
          gene_dif1$celltype <- paste0(tissue,"_",celltype)
          gene_FC <- rbind(gene_FC,gene_dif1)
        }else{
          gene_dif1 <- data.frame(p_val = "0",avg_log2FC = "0",gene = gene,change_V1 = "NOT")
          gene_dif1$celltype <- paste0(tissue,"_",celltype)
          gene_FC <- rbind(gene_FC,gene_dif1)
        }
      }
    }
    else{
      next
    }
  }
}
FC <- gene_FC
FC <- FC[-1,]
write.table(FC, file = paste0("/data1/zhuzn/wangsn/V4/DE/","sig","FC.tab"), quote = FALSE,sep="\t",row.names = FALSE)

FC <- read.delim("/data1/zhuzn/wangsn/V4/DE/sigFC.tab")
FC$avg_log2FC <- as.numeric(FC$avg_log2FC)
FC$p_val <- as.numeric(FC$p_val)
# FC3 <- FC[FC$p_val < 0.05 & FC$change_V1 != "NOT",]

# FC_mean <- aggregate(avg_log2FC~gene,FC,mean)
FC_mean <- aggregate(FC[,1:2],by=list(gene=FC$gene),FUN=mean)
range(FC_mean$avg_log2FC)
# -1.214696  1.444049
range(FC_mean$p_val)
# 2.567196e-234  1.665017e-01
# FC_mean$avg_log2FC <- 6*(FC_mean$avg_log2FC-min)/(max-mix)-3
nrDEG_mean <- merge(change_data, FC_mean, by = "gene", all = TRUE)
nrDEG_mean$sig_gene <- nrDEG_mean$change_V1
nrDEG_mean$sig_gene <- ifelse(nrDEG_mean$change_percent > 50,ifelse(nrDEG_mean$avg_log2FC > 0, 'UP', 'DOWN' ),'FALSE')

gene_sig <- nrDEG_mean[nrDEG_mean$change_percent > 50 & nrDEG_mean$p_val < 0.05,]
gene_sig <- gene_sig[,-2]
gene_sig <- unique(gene_sig)
write.table(gene_sig, file = paste0("/data1/zhuzn/wangsn/V4/DE/","sig","_50.tab"), quote = FALSE,sep="\t",row.names = FALSE)

# figure2 B heatmap(交集TF) -----------------------------------------------------------------
rm(list=ls())
sig_50 <- read.delim("/data1/zhuzn/wangsn/V4/DE/sig_50.tab")
sig_50 <- sig_50[!grepl("ENSMMUG", sig_50$gene),]
TF <- read.table("/data1/zhuzn/wangsn/V6/pyscenic/data/hs_hgnc_tfs.txt", quote="\"", comment.char="")
sig_TF <- intersect(sig_50$gene,TF$V1)
colnames(sig_50)

sig_gene <- sig_50
sig_gene <- sig_gene[order(sig_gene[,6],decreasing=T),]
# sig_50$gene <- as.factor(sig_50$gene)
gene_list <- unique(sig_gene$gene)
tissue_list <- c("Lymph_node","PBMC","Spleen","Bone_marrow")
celltype_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                   "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
gene_dif <- data.frame(p_val = "p_val",avg_log2FC = "avg_log2FC",gene = "gene",change_V1 = "change_V1",
                       celltype = "celltype")
for (tissue in tissue_list){
  for (celltype in celltype_list){
    dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/tissue_celltype/",tissue,"_",celltype,"_old_vs_young.findmark.txt"),stringsAsFactors = F)
    for (gene in gene_list){
      if(gene %in%  dif_data$gene) {
        gene_dif1 <- dif_data[dif_data$gene == gene,c(1,2,6,8)]
        gene_dif1$celltype <- paste0(tissue,"_",celltype)
        gene_dif <- rbind(gene_dif,gene_dif1)
      }else{
        gene_dif1 <- data.frame(p_val = "0",avg_log2FC = "0",gene = gene,change_V1 = "NOT",celltype = celltype)
        gene_dif1$celltype <- paste0(tissue,"_",celltype)
        gene_dif <- rbind(gene_dif,gene_dif1)
      }
    }
  }
}
gene_dif <- gene_dif[-1,]
gene_dif$gene <- factor(gene_dif$gene,levels = gene_list)
gene_dif$avg_log2FC <- as.numeric(gene_dif$avg_log2FC)
gene_dif$celltype <-  gsub("Lymph_node","Mesenteric_lymph",gene_dif$celltype)
library(ggplot2)
p <- ggplot(gene_dif,aes(gene,celltype)) + geom_tile(aes(fill=avg_log2FC),color = "white") +
  guides(fill=guide_colorbar("avg_Expression_FC")) +
  scale_fill_gradient2(limits=c(-4,6),breaks = seq(-4, 6, by = 1),low="#52B088",mid="#eff3ff",high ="#8152E7",midpoint = 0)+
  # scale_fill_gradient2(limits=c(-4,6),breaks = seq(-4, 6, by = 1),low="#154889",mid="#eff3ff",high ="#AB2524",midpoint = 0)+
  labs(x = "",y = "",title = "") +
  theme(axis.text.x = element_text(size = 15, hjust = 1, vjust = 0.5, angle = 90,colour = "black"),
        axis.text.y = element_text(size = 15,colour = "black"))
p
ggsave(paste0('./DE/','sig','_heatmap.jpeg'), plot=p, dpi = 1000,width = 15, height = 10)
ggsave("/data1/zhuzn/wangsn/reviewer/figure2/figure2B.pdf",plot = p,dpi=1000,width = 15,height = 10)
dev.off()

# figure2 C 富集 -----------------------------------------------------------------
rm(list=ls())
sig_50 <- read.delim("/data1/zhuzn/wangsn/V4/DE/sig_50.tab")
# sig_50 <- read.delim("/data1/zhuzn/wangsn/V4/DE/sig_0.tab")

sig_50 <- sig_50[!grepl("ENSMMUG", sig_50$gene),]
SASP_gene_set <- read.delim('/data1/zhuzn/wangsn/SASP/SASP_V2.list',col.names = F)
in_gene <- intersect(sig_50$gene,SASP_gene_set$FALSE.)

#GO
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

sig_gene <- sig_50[sig_50$Var1 == "UP",]
colnames(sig_gene)[2] = "SYMBOL"
DE <- merge(x=sig_gene,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ekegg <- enrichKEGG(gene,organism="mcc",pvalueCutoff=0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Mmu.eg.db)
ekegg_result<-as.data.frame(ekegg@result)

ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_result_UP<-as.data.frame(ego_BP1@result)#8107
ego_BP_result_UP$group <- rep("UP")
c <- c("activation of innate immune response","negative regulation of natural killer cell mediated immunity",
       "negative regulation of T cell mediated cytotoxicity","negative regulation of leukocyte mediated cytotoxicity")
top_go_up <- ego_BP_result_UP[ego_BP_result_UP$Description %in% c,]
top_go_up <- top_go_up[order(top_go_up$Count,decreasing = T),]
top_go_up$Number <- c(4:1)

sig_gene <- sig_50[sig_50$Var1 == "DOWN",]
colnames(sig_gene)[2] = "SYMBOL"
DE <- merge(x=sig_gene,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_result_DOWN<-as.data.frame(ego_BP1@result)#8107
ego_BP_result_DOWN$group <- rep("DOWN")
c <- c("ATP metabolic process","telomerase holoenzyme complex assembly","electron transport chain","protein stabilization","protein folding")
top_go_down <- ego_BP_result_DOWN[ego_BP_result_DOWN$Description %in% c,]
top_go_down <- top_go_down[order(top_go_down$Count,decreasing = F),]
top_go_down$Number <- c(-1:-5)
top_go <- rbind(top_go_up,top_go_down)
top_go$number = factor(rev(1:nrow(top_go)))
top_go$Count = as.numeric(ifelse(top_go$group == 'DOWN' ,paste0("-",top_go$Count),top_go$Count))

GO_BP <- rbind(ego_BP_result_UP,ego_BP_result_DOWN)
write.table(GO_BP, file = paste0("./","figure.2C.txt"), quote = FALSE,sep="\t",row.names = FALSE)

p <- ggplot(data=top_go, aes(x=number, y=Count,fill = Number)) + geom_bar(stat="identity", width=0.8) + 
  scale_x_discrete(labels=rev(top_go$Description)) +
  scale_fill_gradient2(low="#52B088",mid="#eff3ff",high ="#8152E7",midpoint = 0)+
  coord_flip() +theme_bw() + xlab("") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme(axis.text=element_text(size=20,face = "plain", color="black"),panel.grid = element_blank(),
        axis.title=element_text(size=10),legend.text=element_text(size=10),legend.title = element_text(size=20))
p
ggsave("/data1/zhuzn/wangsn/reviewer/figure2/figure.2C.pdf",plot = p,dpi=1000,width = 12,height = 5)

# figure2 D GZMB 雷达图 ----------------------------------------------------------------
# devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
del_ery <- WhichCells(all_sample.combined,expression = celltype  %in% "Erythrocyte")
ncol(del_ery)
all_sample.combined <- subset(all_sample.combined,cells=setdiff(WhichCells(all_sample.combined),del_ery))

data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
age_list <- unique(all_sample.combined$Age)
data_mean <- data.frame()
for (age in age_list){
  data_age <- data.list_age[[age]]
  data.list_tissue <- SplitObject(data_age, split.by = "Tissue")
  tissue_list <- unique(data_age$Tissue)
  for (tissue in tissue_list){
    data_tissue <- data.list_tissue[[tissue]]
    
    data.list_celltype <- SplitObject(data_tissue, split.by = "celltype")
    celltype_list <- unique(data_tissue$celltype)
    for (celltype in celltype_list){
      data_celltype <- data.list_celltype[[celltype]]
      DefaultAssay(data_celltype) <- "RNA"
      data <- data.frame(GetAssayData(data_celltype, slot = "data"));data <- data[rownames(data) == "GZMB",]
      data1 <- data.frame(avg.exp=rowMeans(data),Cluster=celltype,Age=age,Tissue=tissue)
      data_mean <- rbind(data_mean,data1)
    }
  }
}
young <- data_mean[data_mean$Age == "young",]
young <- young[order(young$Tissue, young$Cluster),]
old <- data_mean[data_mean$Age == "old",]
old <- old[order(old$Tissue, old$Cluster),]
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
dt.pct$Tissue <- factor(dt.pct$Tissue,levels = tissue_list)
library(ggradar)
p1 <- ggradar(dt,axis.label.size=6.5,
              # values.radar = c(0, ceiling(max(dt[,-1]))/2, ceiling(max(dt[,-1]))),
              values.radar = c(0, 1, 2),
              group.colours = c('#23B363','#DF63CC','#339FBA','#CA7D44'),
              grid.min = 0,
              grid.mid = 1,
              grid.max = 2) + ggtitle(paste0("old-young",' mean exp'))
p1
ggsave(plot = p1, paste0('/data1/zhuzn/wangsn/reviewer/figure2/','figure2D.pdf'), width = 14, height = 16,dpi = 600)

# figure3 D CD8 36个TF的分析 ----------------------------------------------------------------
### AUC热图
rm(list=ls())
diff_young <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/CD8_cluster7_DE_6-11_young.txt")
diff_young$change_V1 = as.factor(ifelse(diff_young$p_val < 0.05,
                                       ifelse(diff_young$avg_log2FC > 0.5 ,'UP',
                                              ifelse(diff_young$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                       'NOT'))
diff_young <- diff_young[diff_young$change_V1 == "UP",]
# diff_young <- diff_young[diff_young$change_V1 != "NOT",]

diff_old <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/CD8_cluster7_DE_6-11_old.txt")
diff_old$change_V1 = as.factor(ifelse(diff_old$p_val < 0.05,
                                       ifelse(diff_old$avg_log2FC > 0.5 ,'UP',
                                              ifelse(diff_old$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                       'NOT'))
diff_old <- diff_old[diff_old$change_V1 == "UP",]
# diff_old <- diff_old[diff_old$change_V1 != "NOT",]

diff_old_young <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/CD8_cluster_DE_6-11_old_vs_young.txt")
diff_old_young$change_V1 = as.factor(ifelse(diff_old_young$p_val < 0.05,
                                       ifelse(diff_old_young$avg_log2FC > 0.5 ,'UP',
                                              ifelse(diff_old_young$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                       'NOT'))
diff_old_young <- diff_old_young[diff_old_young$change_V1 == "UP",]
# diff_old_young <- diff_old_young[diff_old_young$change_V1 != "NOT",]

regulonAUC <- read.csv("/data1/zhuzn/wangsn/V6/pyscenic/CD8/cluster_BH/AmnionStep3Result.csv", row.names=1)
colnames(regulonAUC) <- gsub("\\.","",colnames(regulonAUC))
TF_list <- unique(colnames(regulonAUC))

overlap <- Reduce(intersect, list(diff_young$gene,diff_old$gene,diff_old_young$gene,TF_list))
#韦恩图
# install.packages("VennDiagram",repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
library(VennDiagram)
venn<- venn.diagram(list(diff_young=diff_young$gene,
                         diff_old=diff_old$gene,
                         diff_old_young=diff_old_young$gene,
                         TF_list=TF_list),
                    filename=NULL,fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
                    #col = "black",
                    col = "transparent", 
                    alpha = 0.4, cat.cex = 1.5,rotation.degree = 0)
pdf("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/figure3C.pdf",height = 4.5,width = 5)
grid.draw(venn)
dev.off()

# figure3 D overlap TF targets FC 富集（按上下调分开） --------------------------------------------------
rm(list=ls())
diff_old_young <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/CD8_cluster_DE_6-11_old_vs_young.txt")
diff_old_young <- diff_old_young[!grepl("ENSMMUG", diff_old_young$gene),]
diff_old_young <- diff_old_young[abs(diff_old_young$avg_log2FC) > 0.5,]
diff_old_young_UP <- diff_old_young[diff_old_young$change_V1 == "UP",c(2,6)]
diff_old_young_DOWN <- diff_old_young[diff_old_young$change_V1 == "DOWN",c(2,6)]

targets_all <- read.delim("/data1/zhuzn/wangsn/V6/pyscenic/CD8/cluster0_11/AmnionStep1Result.tsv")
targets_all <- targets_all[!grepl("ENSMMUG", targets_all$target),]
# TF_lists <- c("BHLHE40","KLF4")
TF = "BHLHE40"
TF = TF_lists
targets <- targets_all[targets_all$TF %in% TF,]

colnames(targets)[2] <- "gene"
targets %<>% join(.,diff_old_young[,c(2,8,6)],by=c("gene"="gene"))
targets <- na.omit(targets)
targets <- targets[targets$change_V1 != "NOT",]
# targets <- targets[targets$importance > 1,]
targets_UP <- targets[targets$change_V1 == "UP",]
targets_DOWN <- targets[targets$change_V1 == "DOWN",]
table(targets$change_V1)
write.table(targets,file = paste0("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/","TF_targets_DE_FC.txt"),quote = FALSE,sep = "\t",row.names = F)
write.table(targets_UP$gene,file = paste0("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/","TF_targets_DE_UP.txt"),quote = FALSE,sep = "\t",row.names = F)
write.table(targets_DOWN$gene,file = paste0("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/","TF_targets_DE_DOWN.txt"),quote = FALSE,sep = "\t",row.names = F)

###富集
rm(list = ls())
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
targets <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/","TF_targets_DE_FC.txt"))
targets1 <- targets[targets$change_V1 == "UP",]
# targets1 <- targets[targets$change_V1 == "DOWN",]
colnames(targets1)[2] = "SYMBOL"
targets1 <- merge(x=targets1,y=geneinfo,by="SYMBOL")
gene <- unique(targets1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_result<-as.data.frame(ego_BP1@result)#8107
ego_BP_result<-ego_BP_result[ego_BP_result$pvalue < 0.05,]
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/","BHLHE40_GO_BP_UP.txt"), quote = FALSE,sep="\t",row.names = FALSE)
# write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer/figure3_BHLHE40/","BHLHE40_GO_BP_DOWN.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# figure4 E 富集(Tcm单独聚类）------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
data.list_age <- SplitObject(Cell.combined, split.by = "Age")
age_data <- data.list_age[["old"]]
data.list_celltype <- SplitObject(age_data, split.by = "celltype")
celltype_data <- data.list_celltype[["Tcm"]]
celltype_data <- SetIdent(object = celltype_data, value = celltype_data@meta.data$celltype1)
diff_gene <- FindMarkers(celltype_data, ident.1 = "Tcm_2", ident.2 = "Tcm_1", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.5 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/DE/","Tcm","_old_1_2.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)

rm(list=ls())
DE <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/DE/Tcm_old_1_2.findmark.txt")
DE <- DE[DE$change_V1 == "UP",]
# DE <- DE[DE$change_V1 == "DOWN",]
colnames(DE)[6] <- "SYMBOL"
df_id <- bitr(DE$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
DE <- merge(DE,df_id,by = "SYMBOL",all=F)
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/DE/","Tcm","_GO_UP.txt"),quote = FALSE,sep = "\t",row.names = F)
# write.table(ego_BP_result,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/DE/","Tcm","_GO_DOWN.txt"),quote = FALSE,sep = "\t",row.names = F)

# GO富集结果合并图 
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
DE_list <- c("UP","DOWN")
GO <- data.frame()
for (DE in DE_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/DE/Tcm_GO_",DE,".txt"),stringsAsFactors = F)
  GO1$DE <- rep(DE)
  GO <- rbind(GO,GO1)
}
GO <- GO[GO$pvalue < 0.05,]
write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer/","Tcm_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

ID <- c("GO:0006954","GO:0031347","GO:0034341","GO:0002222","GO:0140194","GO:0007253","GO:0043922","GO:0002460","GO:0034599",
        "GO:0002250","GO:0098609","GO:0032113","GO:0060550","GO:0050870","GO:0007155","GO:0002250","GO:0002823","GO:0002819")
GO <- GO[GO$ID %in% ID,]

GO$Description <- factor(GO$Description,levels = unique(GO$Description))
GO$Description <- gsub("GO:","",GO$Description)
GO$DE <- factor(GO$DE,levels = unique(GO$DE))
p <- ggplot(GO,aes(Description,DE,color = -log10(pvalue), size=Count)) + geom_point()+
  coord_flip() +
  # scale_colour_gradient(low="grey",high="#CB1217")+
  scale_colour_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#F0BCBC")))+
  geom_point()+scale_size(range=c(5,10))+
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels=function(x) str_wrap(x, width=60))+
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20, color = "black"),plot.title = element_text(size = 30),
        axis.text.y = element_text(size = 25, color = "black"),
        legend.position = "right")
p
unique(GO$Description)
ggsave("/data1/zhuzn/wangsn/reviewer/figure.4E.pdf",plot = p,dpi=1000,width = 12.5,height = 10)
dev.off()

# figure5 A B细胞统计 FC(all_tissue.young_old)-------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/B Cell/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,23)]
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

show_col(colors)
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
plotfile = paste0('/data1/zhuzn/wangsn/reviewer/',"figure.",'B_freq.pdf')
ggsave(plotfile, plot=p1, dpi = 1000, width = 10, height = 4)
dev.off()

# figure5 C B细胞 cluster0 vs cluster8 cluster8 old vs young overlap_富集 --------------------------------------------------------------
rm(list=ls())
DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/V4/sub_celltype/B_Cell_BM_cluster8v0_findmark.txt", header=T)
DEG1$change_V1 = as.factor(ifelse(DEG1$p_val < 0.05,ifelse( DEG1$avg_log2FC > 0.5 ,'UP',ifelse( DEG1$avg_log2FC < -0.25 ,'DOWN','NOT')),'NOT'))
DEG1 <- DEG1[DEG1$change_V1 == "DOWN",]
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/V4/sub_celltype/B_Cell_BM_cluster8_findmark.txt", header=T)
DEG2$change_V1 = as.factor(ifelse(DEG2$p_val < 0.05,ifelse( DEG2$avg_log2FC > 0.5 ,'UP',ifelse( DEG2$avg_log2FC < -0.25 ,'DOWN','NOT')),'NOT'))
DEG2 <- DEG2[DEG2$change_V1 == "DOWN",]
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
library(VennDiagram)
intersect(DEG1$gene,DEG2$gene)
venn<- venn.diagram(list(DEG1=DEG1$gene,
                         DEG2=DEG2$gene),
                    filename=NULL,fill = c("cornflowerblue", "green"),
                    #col = "black",
                    col = "transparent", 
                    alpha = 0.4, cat.cex = 1.5,rotation.degree = 0)
pdf("/data1/zhuzn/wangsn/reviewer/figure5-B_DOWN_overlap.pdf",height = 5,width = 5)
grid.draw(venn)
dev.off()

# figure5 E cluster0 vs cluster8 cluster8 old vs young overlap_富集 --------------------------------------------------------------
rm(list=ls())
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
DEG1 <- read.delim(file = "/data1/zhuzn/wangsn/V4/sub_celltype/B_Cell_BM_cluster8v0_findmark.txt", header=T)
DEG1$change_V1 = as.factor(ifelse(DEG1$p_val < 0.05,ifelse( DEG1$avg_log2FC > 0.5 ,'UP',ifelse( DEG1$avg_log2FC < -0.25 ,'DOWN','NOT')),'NOT'))
DEG1 <- DEG1[DEG1$change_V1 == "DOWN",]
DEG2 <- read.delim(file = "/data1/zhuzn/wangsn/V4/sub_celltype/B_Cell_BM_cluster8_findmark.txt", header=T)
DEG2$change_V1 = as.factor(ifelse(DEG2$p_val < 0.05,ifelse( DEG2$avg_log2FC > 0.5 ,'UP',ifelse( DEG2$avg_log2FC < -0.25 ,'DOWN','NOT')),'NOT'))
DEG2 <- DEG2[DEG2$change_V1 == "DOWN",]

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
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer/","B_DOWN_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

library("ggplot2")
GO <- c("GO:0048006","GO:0035442","GO:0048007","GO:0033632","GO:0019883")
GO  <- ego_BP_result[ego_BP_result$ID %in% GO,]

library(DOSE)
GO$GeneRatio1 <- parse_ratio(GO$GeneRatio)
p <- ggplot(GO,aes(reorder(Description,GeneRatio1),GeneRatio1,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_colour_gradient(low="#3E1AFD",high="#C1B5FC")+
  scale_x_discrete(labels=function(x) str_wrap(x, width=42))+
  geom_point(size = 2.0,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer/","B Cell_GO_DOWN.pdf"),plot = p,width = 5.5,height = 3)

# figure6 髓系CellChat-------------------------------------------------------------------------
rm(list=ls())
library(CellChat)
library(patchwork)
library(tidyverse)
library(ggalluvial)
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
unique(all_sample.combined$celltype)
data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
# data_age <- data.list_age[["old"]]
data_age <- data.list_age[["young"]]

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

cells <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte")
mat <- cellchat@net$weight
i = 9
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[7, ] <- mat[7, ]
mat2 <- mat2[cells,cells]
colnames(mat)
colors <- c('#1f78b4','#a6cee3',"#09963B","#fb9a99","#ff7f00","#008A8A","#1D1D1D","#D31C8D","#e31a1c",'#6a3d9a')
show_col(colors)
# pdf(paste0("/data1/zhuzn/wangsn/V7/result/figure_4F3.pdf"), width = 5, height = 5)
p <- netVisual_circle(mat2, color.use = colors,vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat2)[i])
print(p)
dev.off()

cellchat@netP$pathways
pathways.show <- c("MIF") 
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# netVisual_bubble(cellchat, signaling = c("CCL"), remove.isolate = FALSE)

# pdf(paste0("/data1/zhuzn/wangsn/V7/result/figure_4F5.pdf"), width = 4, height = 3)
colors <- rev(c("#8DD3C7","#FB8072","#80B1D3","#36D936","#FDB462","#B4CBEB","#FCCDE5","#BC80BD","#B3DE69","#BEBADA"))
p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.use = colors,color.heatmap = "Reds")
print(p)
dev.off()
# saveRDS(cellchat, file = "/data1/zhuzn/wangsn/reviewer/cellchat_Myeloid_old.rds")
saveRDS(cellchat, file = "/data1/zhuzn/wangsn/reviewer/cellchat_Myeloid_young.rds")
dev.off()

# figure6 CellChat结果-------------------------------------------------------------------------
rm(list=ls())
library(CellChat)
library(patchwork)
library(tidyverse)
library(ggalluvial)
cellchat1 <- readRDS("/data1/zhuzn/wangsn/reviewer/cellchat_Myeloid_old.rds")
cellchat2 <- readRDS("/data1/zhuzn/wangsn/reviewer/cellchat_Myeloid_young.rds")

groupSize <- as.numeric(table(cellchat1@idents))
cells <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte")
cellchat1@net$weight1 <- cellchat1@net$weight-cellchat2@net$weight
mat <- cellchat1@net$weight1
i = 9
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[7, ] <- mat[7, ]
mat2 <- mat2[cells,cells]
colnames(mat)
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#DBBAFA","#FF348B","#B7D5F3","#B3DE69","#FFD95A")
show_col(colors)
pdf(paste0("/data1/zhuzn/wangsn/reviewer/figure_6G.pdf"), width = 5, height = 5)
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
pdf(paste0("/data1/zhuzn/wangsn/reviewer/figure_6H.pdf"), width = 4, height = 3)
cells <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell","Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte")
colors <- c("#1084D3","#DBBAFA","#B7D5F3","#FF348B","#FFD95A","#CFE6A7","#B3DE69","#73C508","#BEBADA","#7756A1")
p <- netVisual_heatmap(cellchat1, signaling = pathways.show, color.use = colors)
print(p)
dev.off()

# figureS1 FindMarkers-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$Age)

data.list_tissue <- SplitObject(all_sample.combined, split.by = "Tissue")
sample_list <- unique(all_sample.combined@meta.data$Tissue)
for (tissue in sample_list){
  tissue_data <- data.list_tissue[[tissue]]
  
  diff_gene <- FindMarkers(tissue_data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
  diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
  diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                         ifelse(diff_gene$avg_log2FC > 0.5 ,'UP',
                                                ifelse(diff_gene$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                         'NOT'))
  write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/V4/DE/",tissue,"_old_vs_young.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
}

# figureS1 DEGs数量_UpSetR-------------------------------------------------------------------------
rm(list=ls())
sample_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
Bone_marrow <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/","Bone_marrow","_old_vs_young.findmark.txt"))
Mesenteric_lymph <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/","Mesenteric_lymph","_old_vs_young.findmark.txt"))
PBMC <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/","PBMC","_old_vs_young.findmark.txt"))
Spleen <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/","Spleen","_old_vs_young.findmark.txt"))

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
pdf(file = paste0('./result/',"figure.",'2A.pdf'), p1, width = 7, height = 5)
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
pdf(file = paste0('./result/',"figure.",'2B.pdf'), p2, width = 7, height = 5)
print(p2)
dev.off()

# figureS1DE_union（上调-GO富集-共同）-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
Bone_marrow <- read.csv(paste0("/data1/zhuzn/wangsn/V4/DE/","Bone_marrow","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Bone_marrow <- Bone_marrow[!grepl("ENSMMUG", Bone_marrow$gene),]
Mesenteric_lymph <- read.csv(paste0("/data1/zhuzn/wangsn/V4/DE/","Mesenteric_lymph","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
Mesenteric_lymph <- Mesenteric_lymph[!grepl("ENSMMUG", Mesenteric_lymph$gene),]
PBMC <- read.csv(paste0("/data1/zhuzn/wangsn/V4/DE/","PBMC","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
PBMC <- PBMC[!grepl("ENSMMUG", PBMC$gene),]
Spleen <- read.csv(paste0("/data1/zhuzn/wangsn/V4/DE/","Spleen","_old_vs_young.findmark.txt"), sep="") %>% .[.$change_V1 == "UP",]
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
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/V4/DE/GO/","union","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
Bone_marrow1 <- Bone_marrow1[order(Bone_marrow1$avg_log2FC,decreasing = T),]
write.table(Bone_marrow1,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/setdiff/","Bone_marrow","_tissue_setdiff_UP.txt"),quote = FALSE,sep = "\t",row.names = T,col.names = T)
colnames(Bone_marrow1)[6] <- "SYMBOL"
df_id <- bitr(Bone_marrow1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Bone_marrow1 <- merge(Bone_marrow1,df_id,by = "SYMBOL",all=F)
gene <- unique(Bone_marrow1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/V4/DE/GO/","Bone_marrow","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
Mesenteric_lymph1 <- Mesenteric_lymph1[order(Mesenteric_lymph1$avg_log2FC,decreasing = T),]
write.table(Mesenteric_lymph1,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/setdiff/","Mesenteric_lymph","_tissue_setdiff_UP.txt"),quote = FALSE,sep = "\t",row.names = T,col.names = T)
colnames(Mesenteric_lymph1)[6] <- "SYMBOL"
df_id <- bitr(Mesenteric_lymph1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Mesenteric_lymph1 <- merge(Mesenteric_lymph1,df_id,by = "SYMBOL",all=F)
gene <- unique(Mesenteric_lymph1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/V4/DE/GO/","Mesenteric_lymph","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
PBMC1 <- PBMC1[order(PBMC1$avg_log2FC,decreasing = T),]
write.table(PBMC1,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/setdiff/","PBMC","_tissue_setdiff_UP.txt"),quote = FALSE,sep = "\t",row.names = T,col.names = T)
colnames(PBMC1)[6] <- "SYMBOL"
df_id <- bitr(PBMC1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
PBMC1 <- merge(PBMC1,df_id,by = "SYMBOL",all=F)
gene <- unique(PBMC1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/V4/DE/GO/","PBMC","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]
Spleen1 <- Spleen1[order(Spleen1$avg_log2FC,decreasing = T),]
write.table(Spleen1,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/setdiff/","Spleen","_tissue_setdiff_UP.txt"),quote = FALSE,sep = "\t",row.names = T,col.names = T)
colnames(Spleen1)[6] <- "SYMBOL"
df_id <- bitr(Spleen1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
Spleen1 <- merge(Spleen1,df_id,by = "SYMBOL",all=F)
gene <- unique(Spleen1[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/V4/DE/GO/","Spleen","_setdiff_GO_UP.tab"), quote = FALSE,sep="\t",row.names = FALSE)

# figureS1 D DE_union（上调-共同-特异）GO富集-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
DE_list <- c("union","Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO <- data.frame()
for (celltype in DE_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/GO/",celltype,"_setdiff_GO_UP.tab"),stringsAsFactors = F)
  # GO_change <- GO1[c(1:5),]
  GO_change <- GO1
  GO_change$celltype <- rep(celltype)
  GO <- rbind(GO,GO_change)
  write.table(GO, file = paste0("./","figure.S1D.txt"), quote = FALSE,sep="\t",row.names = FALSE)
}

# figureS1 E DE_union（下调-共同-特异）GO富集-------------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
DE_list <- c("union","Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO <- data.frame()
for (celltype in DE_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/V4/DE/GO/",celltype,"_setdiff_GO_DOWN.tab"),stringsAsFactors = F)
  # GO_change <- GO1[c(1:5),]
  GO_change <- GO1
  GO_change$celltype <- rep(celltype)
  GO <- rbind(GO,GO_change)
  write.table(GO, file = paste0("./","figure.S1E.txt"), quote = FALSE,sep="\t",row.names = FALSE)
}

# figureS2 NK DEGs(按celltype做差异） ------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/NK/Cell.combined.celltype.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Age)
data.list_celltype <- SplitObject(Cell.combined, split.by = "celltype")
for (celltype in unique(Cell.combined$celltype)){
  Cell.combined3 <- data.list_celltype[[celltype]]
  d <- data.frame(table(Cell.combined3$Age))
  if(nrow(d) == 1){
    print(paste0(celltype," no young or old cells"))
  }
  else{
    a = d[1,2];b = d[2,2]
    if(a<3 | b <3){
      print(paste0(celltype,"Cells fewer than 3 cells")) 
    }else{
      diff_gene <- FindMarkers(Cell.combined3,ident.1 = "old", ident.2 = "young",logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
      diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
      diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                             ifelse(diff_gene$avg_log2FC > 0.5 ,'UP',
                                                    ifelse(diff_gene$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                             'NOT'))
      diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
      celltype <- gsub(" ","_",celltype)
      write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/V7/DE_NK/",celltype,"_DE.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
    } 
  }
}

# figureS2 NK不同亚型DEG_setdiff ------------------------------------------------------------
rm(list=ls())
NK1 <- read.delim("/data1/zhuzn/wangsn/V7/DE_NK/NK1_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK2 <- read.delim("/data1/zhuzn/wangsn/V7/DE_NK/NK2_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK3 <- read.delim("/data1/zhuzn/wangsn/V7/DE_NK/NK3_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
library(reshape2)
NK11 <- NK1[NK1$gene %in% Reduce(setdiff,list(NK1$gene,NK2$gene,NK3$gene)),]
NK11$celltype <- rep("NK11")
NK21 <- NK2[NK2$gene %in% Reduce(setdiff,list(NK2$gene,NK1$gene,NK3$gene)),]
NK21$celltype <- rep("NK21")
NK31 <- NK3[NK3$gene %in% Reduce(setdiff,list(NK3$gene,NK1$gene,NK2$gene)),]
NK31$celltype <- rep("NK31")

diff_gene <- c(NK11$gene,NK21$gene,NK31$gene)

celltype.lists <- c("NK1","NK2","NK3")
diff_gene_overlap_FC <- data.frame()
for (j in celltype.lists){
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/V7/DE_NK/",j,"_DE.findmark.txt"), header=T)
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
  diff_gene1$celltype <- rep(j)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}
diff_gene_overlap_FC <- diff_gene_overlap_FC[!grepl("ENSMMUG", diff_gene_overlap_FC$gene),]
overlap_FC <- dcast(diff_gene_overlap_FC,celltype~diff_gene_overlap_FC$gene,value.var = "avg_log2FC",fill = "0")
overlap_FC1 <- data.frame(t(overlap_FC));colnames(overlap_FC1) <- overlap_FC1[1,];overlap_FC1 <- overlap_FC1[-1,]
overlap_FC1 <- overlap_FC1[diff_gene,]
overlap_FC <- melt(overlap_FC,id.vars = "celltype",variable.name='gene',value.name="avg_log2FC")
colnames(overlap_FC) <- c("celltype","gene","avg_log2FC")
str(overlap_FC)
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = celltype.lists)
# overlap_FC[overlap_FC > 1] <- 1
# overlap_FC[overlap_FC < -1] <- -1

p <- ggplot(overlap_FC,aes(x=celltype,y=gene,fill=avg_log2FC))+
  geom_raster()+
  # scale_fill_gradientn(colors = c("grey90","grey90","#a50f15"))+
  scale_fill_gradient2(low="#4592E5",mid="#eff3ff",high ="#F70000",midpoint = 0)+
  theme_bw() + 
  labs(title = "")+ 
  theme(legend.title = element_text(size = 20),legend.text = element_text(size = 20),
        panel.grid = element_blank(),plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle = 90,vjust = 1,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.border = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer/figureS2","_NK_DEG_UP.pdf"),plot = p,width = 5,height = 6)
dev.off()

# figureS2 NK 不同细胞类型_DE_setdiff-富集-------------------------------------------------------------------------
rm(list=ls())
NK1 <- read.delim("/data1/zhuzn/wangsn/V7/DE_NK/NK1_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK2 <- read.delim("/data1/zhuzn/wangsn/V7/DE_NK/NK2_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK3 <- read.delim("/data1/zhuzn/wangsn/V7/DE_NK/NK3_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

NK11 <- NK1[NK1$gene %in% Reduce(setdiff,list(NK1$gene,NK2$gene,NK3$gene)),]
colnames(NK11)[6] = "SYMBOL"
DE <- merge(x=NK11,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_NK11<-as.data.frame(ego_BP1@result)
ego_BP_NK11$celltype <- rep("NK1")
ego_BP_NK11 <- ego_BP_NK11[ego_BP_NK11$pvalue < 0.05,]
# ego_BP_NK11<-ego_BP_NK11[c(1:5),]
NK21 <- NK2[NK2$gene %in% Reduce(setdiff,list(NK2$gene,NK1$gene,NK3$gene)),]
colnames(NK21)[6] = "SYMBOL"
DE <- merge(x=NK21,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_NK21<-as.data.frame(ego_BP1@result)
ego_BP_NK21$celltype <- rep("NK2")
ego_BP_NK21 <- ego_BP_NK21[ego_BP_NK21$pvalue < 0.05,]
# ego_BP_NK21<-ego_BP_NK21[c(1:5),]
NK31 <- NK3[NK3$gene %in% Reduce(setdiff,list(NK3$gene,NK1$gene,NK2$gene)),]
colnames(NK31)[6] = "SYMBOL"
DE <- merge(x=NK31,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_NK31<-as.data.frame(ego_BP1@result)
ego_BP_NK31$celltype <- rep("NK3")
ego_BP_NK31 <- ego_BP_NK31[ego_BP_NK31$pvalue < 0.05,]
# ego_BP_NK31<-ego_BP_NK31[c(1:5),]
GO <- rbind(ego_BP_NK11,ego_BP_NK21,ego_BP_NK31)
write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer/","NK_UP_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# figureS3 CD4 不同亚型DEG_setdiff ------------------------------------------------------------
rm(list=ls())
celltype.lists <- c("Naive","Tem","Tcm","TRM","CTL")
for (j in celltype.lists){
  diff_gene <- read.delim(paste0("/data1/zhuzn/wangsn/V7/DE_CD4/",j,"_DE.findmark.txt"))
  diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                         ifelse(diff_gene$avg_log2FC > 0.5 ,'UP',
                                                ifelse(diff_gene$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                         'NOT'))
  write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/V7/DE_CD4/",j,"_DE.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
}

Naive <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/Naive_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tem <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/Tem_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tcm <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/Tcm_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
TRM <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/TRM_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
CTL <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/CTL_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
library(reshape2)
Naive1 <- Naive[Naive$gene %in% Reduce(setdiff,list(Naive$gene,Tem$gene,Tcm$gene,TRM$gene,CTL$gene)),]
Naive1$celltype <- rep("Naive1")
Tem1 <- Tem[Tem$gene %in% Reduce(setdiff,list(Tem$gene,Naive$gene,Tcm$gene,TRM$gene,CTL$gene)),]
Tem1$celltype <- rep("Tem1")
Tcm1 <- Tcm[Tcm$gene %in% Reduce(setdiff,list(Tcm$gene,Naive$gene,Tem$gene,TRM$gene,CTL$gene)),]
Tcm1$celltype <- rep("Tcm1")
TRM1 <- TRM[TRM$gene %in% Reduce(setdiff,list(TRM$gene,Naive$gene,Tem$gene,Tcm$gene,CTL$gene)),]
TRM1$celltype <- rep("TRM1")
CTL1 <- CTL[CTL$gene %in% Reduce(setdiff,list(CTL$gene,Naive$gene,Tem$gene,Tcm$gene,TRM$gene)),]
CTL1$celltype <- rep("CTL1")

diff_gene <- c(Naive1$gene,Tcm1$gene,Tem1$gene,TRM1$gene,CTL1$gene)

celltype.lists <- c("Naive","Tcm","Tem","TRM","CTL")
diff_gene_overlap_FC <- data.frame()
for (j in celltype.lists){
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/V7/DE_CD4/",j,"_DE.findmark.txt"), header=T)
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
  diff_gene1$celltype <- rep(j)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}
diff_gene_overlap_FC <- diff_gene_overlap_FC[!grepl("ENSMMUG", diff_gene_overlap_FC$gene),]
overlap_FC <- dcast(diff_gene_overlap_FC,celltype~diff_gene_overlap_FC$gene,value.var = "avg_log2FC",fill = "0")
overlap_FC1 <- data.frame(t(overlap_FC));colnames(overlap_FC1) <- overlap_FC1[1,];overlap_FC1 <- overlap_FC1[-1,]
overlap_FC1 <- overlap_FC1[diff_gene,]
overlap_FC <- melt(overlap_FC,id.vars = "celltype",variable.name='gene',value.name="avg_log2FC")
colnames(overlap_FC) <- c("celltype","gene","avg_log2FC")
str(overlap_FC)
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = celltype.lists)
# overlap_FC[overlap_FC > 0.5] <- 0.5
# overlap_FC[overlap_FC < -0.5] <- -0.5

overlap_FC$celltype <- factor(overlap_FC$celltype,levels = c("Naive","Tcm","Tem","TRM","CTL"))
p <- ggplot(overlap_FC,aes(x=celltype,y=gene,fill=avg_log2FC))+
  geom_raster()+
  # scale_fill_gradientn(colors = c("grey90","grey90","#a50f15"))+
  scale_fill_gradient2(low="#4592E5",mid="#eff3ff",high ="#F70000",midpoint = 0)+
  theme_bw() + 
  labs(title = "")+ 
  theme(legend.title = element_text(size = 20),legend.text = element_text(size = 20),
        panel.grid = element_blank(),plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle = 90,vjust = 1,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.border = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer/","figureS3_CD4_DEG_UP.pdf"),plot = p,width = 4,height = 4)
dev.off()

# figureS3 CD4 不同细胞类型_DE_setdiff-富集-------------------------------------------------------------------------
rm(list=ls())
Naive <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/Naive_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tem <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/Tem_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tcm <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/Tcm_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
TRM <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/TRM_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
CTL <- read.delim("/data1/zhuzn/wangsn/V7/DE_CD4/CTL_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

Naive1 <- Naive[Naive$gene %in% Reduce(setdiff,list(Naive$gene,Tem$gene,Tcm$gene,TRM$gene,CTL$gene)),]
# colnames(Naive1)[6] = "SYMBOL"
# DE <- merge(x=Naive1,y=geneinfo,by="SYMBOL")
# gene <- unique(DE[,'ENTREZID'])
# ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
#                    pAdjustMethod = "BH",minGSSize = 1,
#                    pvalueCutoff = 0.05)
# ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
# ego_BP_Naive1<-as.data.frame(ego_BP1@result)
# ego_BP_Naive1$celltype <- rep("Naive")
Tem1 <- Tem[Tem$gene %in% Reduce(setdiff,list(Tem$gene,Naive$gene,Tcm$gene,TRM$gene,CTL$gene)),]
colnames(Tem1)[6] = "SYMBOL"
DE <- merge(x=Tem1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_Tem1<-as.data.frame(ego_BP1@result)
ego_BP_Tem1$celltype <- rep("Tem")
Tcm1 <- Tcm[Tcm$gene %in% Reduce(setdiff,list(Tcm$gene,Naive$gene,Tem$gene,TRM$gene,CTL$gene)),]
colnames(Tcm1)[6] = "SYMBOL"
DE <- merge(x=Tcm1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_Tcm1<-as.data.frame(ego_BP1@result)
ego_BP_Tcm1$celltype <- rep("Tcm")
TRM1 <- TRM[TRM$gene %in% Reduce(setdiff,list(TRM$gene,Naive$gene,Tem$gene,Tcm$gene,CTL$gene)),]
colnames(TRM1)[6] = "SYMBOL"
DE <- merge(x=TRM1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_TRM1<-as.data.frame(ego_BP1@result)
ego_BP_TRM1$celltype <- rep("TRM")
CTL1 <- CTL[CTL$gene %in% Reduce(setdiff,list(CTL$gene,Naive$gene,Tem$gene,Tcm$gene,TRM$gene)),]
colnames(CTL1)[6] = "SYMBOL"
DE <- merge(x=CTL1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_CTL1<-as.data.frame(ego_BP1@result)
ego_BP_CTL1$celltype <- rep("CTL")

# GO <- rbind(ego_BP_Naive1,ego_BP_Tcm1,ego_BP_Tem1,ego_BP_TRM1,ego_BP_CTL1)
GO <- rbind(ego_BP_Tcm1,ego_BP_TRM1,ego_BP_CTL1)
GO <- GO[GO$pvalue < 0.05,]
write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer/","CD4_UP_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# figureS4 CD8 features_Bubble-----------------------------------------------------------------------
rm(list = ls())
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined_CD8.after_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.0.4)
cell_markers <- read.delim("/data1/zhuzn/wangsn/reviewer/CD8/celltype_marker.txt")
gene <- cell_markers$Gene
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.0.4',col.min=-2,col.max = 2,
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

# figureS4 CD8 features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$celltype)

gene <- c("CCR7","CCR6","GZMK","CD69","HAVCR2","CTLA4","ITGAE","EOMES","RUNX3")
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'celltype',col.min=-2,col.max = 2,
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
ggsave("/data1/zhuzn/wangsn/reviewer/CD8/figure.S4B.pdf",plot = p,dpi=1000,width = 6.3,height = 5.3)

# SASP基因（问题6 TNF,P16,P21） -------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
colnames(all_sample.combined@meta.data)
all_sample.combined <- AddMetaData(all_sample.combined,all_sample.combined@reductions$umap@cell.embeddings,col.name = colnames(all_sample.combined@reductions$umap@cell.embeddings))
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$celltype)
DefaultAssay(all_sample.combined) <- "RNA"
data.list <- SplitObject(all_sample.combined, split.by = "Tissue")
unique(all_sample.combined$Tissue)
genes <- c("IL6","TNF","CDKN1A")
mat <- data.frame()
data1 <- data.list[["Lymph_node"]]
mat1=as.data.frame(data1[["RNA"]]@data["IL6",])
colnames(mat1)="exp"
mat2=Embeddings(data1,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat0 <- data.frame(tissue = "Lymph_node",gene = "IL6",exp = sum(mat3$exp))
mat <- rbind(mat,mat0)
mat3$exp[mat3$exp > 3] = 3
colors=rev(c("#B30000","#CE1212","#F03B3B","#F05A5A","#F07979","#FFAFAF","#FFDBDB","#FFE6E6"))
show_col(colors)
p1 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p1
for (gene in genes[2:3]){
  mat1=as.data.frame(data1[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(data1,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat0 <- data.frame(tissue = "Lymph_node",gene = gene,exp = sum(mat3$exp))
  mat <- rbind(mat,mat0)
  mat3$exp[mat3$exp > 3] = 3
  p12 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p1 <- p1|p12
}

data2 <- data.list[["PBMC"]]
mat1=as.data.frame(data2[["RNA"]]@data["IL6",])
colnames(mat1)="exp"
mat2=Embeddings(data2,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat0 <- data.frame(tissue = "PBMC",gene = "IL6",exp = sum(mat3$exp))
mat <- rbind(mat,mat0)
mat3$exp[mat3$exp > 3] = 3
p2 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
for (gene in genes[2:3]){
  mat1=as.data.frame(data2[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(data2,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat0 <- data.frame(tissue = "PBMC",gene = gene,exp = sum(mat3$exp))
  mat <- rbind(mat,mat0)
  mat3$exp[mat3$exp > 3] = 3
  p22 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p2 <- p2|p22
}

data3 <- data.list[["Spleen"]]
mat1=as.data.frame(data3[["RNA"]]@data["IL6",])
colnames(mat1)="exp"
mat2=Embeddings(data3,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat0 <- data.frame(tissue = "Spleen",gene = "IL6",exp = sum(mat3$exp))
mat <- rbind(mat,mat0)
mat3$exp[mat3$exp > 3] = 3
p3 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
for (gene in genes[2:3]){
  mat1=as.data.frame(data3[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(data3,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat0 <- data.frame(tissue = "Spleen",gene = gene,exp = sum(mat3$exp))
  mat <- rbind(mat,mat0)
  mat3$exp[mat3$exp > 3] = 3
  p33 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p3 <- p3|p33
}

data4 <- data.list[["Bone_marrow"]]
mat1=as.data.frame(data4[["RNA"]]@data["IL6",])
colnames(mat1)="exp"
mat2=Embeddings(data4,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat0 <- data.frame(tissue = "Bone_marrow",gene = "IL6",exp = sum(mat3$exp))
mat <- rbind(mat,mat0)
mat3$exp[mat3$exp > 3] = 3
p4 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
for (gene in genes[2:3]){
  mat1=as.data.frame(data4[["RNA"]]@data[gene,])
  colnames(mat1)="exp"
  mat2=Embeddings(data4,"umap")
  mat3=merge(mat2,mat1,by="row.names")
  mat0 <- data.frame(tissue = "Bone_marrow",gene = gene,exp = sum(mat3$exp))
  mat <- rbind(mat,mat0)
  mat3$exp[mat3$exp > 3] = 3
  p44 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
    scale_color_gradientn(limits = c(0,3),breaks = c(0,1.5,3),colors=colors)+theme_bw()+
    theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",
          panel.grid = element_blank())
  p4 <- p4|p44
}
p <- p1/p2/p3/p4
p
ggsave("/data1/zhuzn/wangsn/reviewer/SASP.png",plot = p,width = 24,height = 6)

# SASP基因和PBMC差异基因交集 -------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
ncol(all_sample.combined)
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$Age)
data.list_tissue <- SplitObject(all_sample.combined, split.by = "Tissue")

tissue_data <- data.list_tissue[["PBMC"]]
diff_gene <- FindMarkers(tissue_data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.5 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.5 ,'DOWN','NOT')),
                                       'NOT'))
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer/","PBMC","_.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)

rm(list=ls())
diff_gene <- read.delim("/data1/zhuzn/wangsn/reviewer/PBMC_.findmark.txt")
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,ifelse( diff_gene$avg_log2FC > 0 ,'UP','DOWN'),'NOT'))
diff_gene <- diff_gene[diff_gene$change_V1 != "NOT",]
SASP_gene_set <- read.table("/data1/zhuzn/wangsn/SASP/SASP_V2.list", quote="\"", comment.char="")
diff_gene <- diff_gene[diff_gene$gene %in% intersect(diff_gene$gene,SASP_gene_set$V1),]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer/","PBMC","_SASP_intersect.txt"),quote = FALSE,sep = "\t",row.names = F)

# GZMB表达 ------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/NK/Cell.combined.celltype.RData")
Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
p <- FeaturePlot(Cell.combined, features = "GZMB",reduction = "umap",cols = c("grey", "red"),pt.size = 0.01,label = F)+
  ggtitle("NK")+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p
p <- FeaturePlot(Cell.combined, features = "GZMB",reduction = "umap",cols = c("grey", "red"),pt.size = 0.01,label = F,
                 split.by = "Tissue_Age")+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p

load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
p <- FeaturePlot(Cell.combined, features = "GZMB",reduction = "umap",cols = c("grey", "red"),pt.size = 0.01,label = F)+
  ggtitle("CD8")+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p
p <- FeaturePlot(Cell.combined, features = "GZMB",reduction = "umap",cols = c("grey", "red"),pt.size = 0.01,label = F,
                 split.by = "Tissue_Age")+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p

load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/Granulocyte.RData")
Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
p <- FeaturePlot(Cell.combined, features = "GZMB",reduction = "umap",cols = c("grey", "red"),pt.size = 0.01,label = F)+
  ggtitle("Granulocyte")+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p
p <- FeaturePlot(Cell.combined, features = "GZMB",reduction = "umap",cols = c("grey", "red"),pt.size = 0.01,label = F,
                 split.by = "Tissue_Age")+
  plot_annotation(theme = theme(plot.title = element_text(size = 10,hjust = 0.5,face = "bold")))
p

# 质控标准 --------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/merge_all_sample_seurat.RData")
p1 <- VlnPlot(alldata,features = "nFeature_RNA", pt.size = 0) + NoLegend()+
  theme(axis.text.x = element_text(size = 15, angle = 65))+
  plot_annotation(title = paste0("Before_filter_",ncol(alldata)),theme = theme(plot.title = element_text(size = 20)))
p2 <- VlnPlot(alldata,features = "percent.mito", pt.size = 0) + NoLegend()+
  theme(axis.text.x = element_text(size = 15, angle = 65))
p <- p1/p2
p
dev.off()

del_all <- WhichCells(alldata,expression = (nFeature_RNA>4000 | nFeature_RNA<200 | percent.mito>5))
data.filt <- subset(alldata,cells=setdiff(WhichCells(alldata),del_all))
ncol(data.filt)
p1 <- VlnPlot(data.filt,features = "nFeature_RNA", pt.size = 0) + NoLegend()+
  theme(axis.text.x = element_text(size = 15, angle = 65))+
  plot_annotation(title = paste0("After_filter_",ncol(data.filt)),theme = theme(plot.title = element_text(size = 20)))
p2 <- VlnPlot(data.filt,features = "percent.mito", pt.size = 0) + NoLegend()+
  theme(axis.text.x = element_text(size = 15, angle = 65))
p <- p1/p2
p
dev.off()

# ############# CD8 celltype_umap #############---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
Cell.combined@meta.data$celltype <- factor(Cell.combined@meta.data$celltype,levels = c("Naive","Tem","Tcm","Tex","CTL/TRM","TRM","CTL"))
Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$celltype)

colors <- c("#46ADE3","#92DF4D","#339933","#E3393B","#FB7473","#FF6600","#B064D6")
show_col(colors)
p <- DimPlot(Cell.combined, reduction = "umap",cols = colors, 
              split.by = "Tissue_Age",ncol = 4,
              raster=FALSE,repel = F, label = F)+xlab("") + ylab("") + 
  theme(axis.ticks = element_blank(),axis.line = element_blank(),
        legend.position = 'none',
        axis.text = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer/figure.CD8.png",plot = p,dpi=1000,width = 20,height = 8)

# HSC features_Bubble_cluster ------------------------------------------------------------
rm(list = ls())
load(paste0("/data1/zhuzn/wangsn/V4/data/","merge_all_sample",".combined.after_cluster.RData"))
DefaultAssay(all_sample.combined) <- "RNA"
all_sample.combined <- SetIdent(object = all_sample.combined, value = all_sample.combined@meta.data$integrated_snn_res.1)
table(all_sample.combined@meta.data$integrated_snn_res.1)

cell_markers <- read.delim("/data1/zhuzn/wangsn/V4/celltype/all_sample_celltype_marker2.txt")
unique(cell_markers$Cell.type)
gene <- as.character(unique(cell_markers$Gene))
p1<-DotPlot(object = all_sample.combined, features = gene, group.by  = 'integrated_snn_res.1',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
library(magrittr)
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
data <- na.omit(data)
data$Cluster <- factor(data$Cluster,levels = c(0:35))
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ffcc29",midpoint = 0)+
  # scale_color_manual(values=color16)+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black"),axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1))
p

# 粒细胞亚型 -------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
DefaultAssay(all_sample.combined) <- "integrated"
celltype_lists <- unique(all_sample.combined@meta.data$celltype)
Cell.combined<-subset(all_sample.combined, subset = celltype  %in% "Granulocyte")
set.seed(1)
Cell.combined <- RunPCA(Cell.combined, npcs = 30, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 42, reduction = "pca", dims = 1:20)
Cell.combined <- RunTSNE(Cell.combined, seed.use = 2, reduction = "pca", dims = 1:20)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:20)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.1,0.2,0.3,0.4))
save(Cell.combined, file="/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/Granulocyte.RData")

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/Granulocyte.RData")
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

celltype <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/celltype.txt")
unique(celltype$celltype)
colnames(Cell.combined@meta.data)
table(Cell.combined@meta.data$integrated_snn_res.0.6)
cluster <- Cell.combined@meta.data[,c(1,13)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "sub_celltype")
save(Cell.combined, file="/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/Granulocyte/Granulocyte.RData")

Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
p <- DimPlot(Cell.combined, reduction = "umap", pt.size = 0.1,group.by = "sub_celltype", 
             # split.by = "Tissue_Age",ncol = 4,
              label.size = 4,label = F,repel = T)+
  # theme(legend.position = 'none') +
  xlab("") + ylab("")
p

# 淋巴结不同细胞类型DEG_setdiff ------------------------------------------------------------
rm(list=ls())
celltype.lists <- c("CD14_Mon","CD16_Mon","In-termed_Mon","M1","M2","Megakaryocyte","Granulocyte","Dendritic_cell")
for (j in celltype.lists){
  diff_gene <- read.delim(paste0("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_",j,"_DE.txt"))
  diff_gene$change_V = diff_gene$change_V1
  diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,ifelse( diff_gene$avg_log2FC > 0.5 ,'UP','DOWN'),'NOT'))
  diff_gene$change_V2 = as.factor(ifelse(diff_gene$p_val < 0.05,ifelse( diff_gene$avg_log2FC > 1 ,'UP','DOWN'),'NOT'))
  diff_gene$change_V3 = as.factor(ifelse(diff_gene$p_val < 0.05,ifelse( diff_gene$avg_log2FC > 1.5 ,'UP','DOWN'),'NOT'))
  diff_gene$change_V4 = as.factor(ifelse(diff_gene$p_val < 0.05,ifelse( diff_gene$avg_log2FC > 2 ,'UP','DOWN'),'NOT'))
  write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_",j,"_DE.txt"),quote = FALSE,sep = "\t",row.names = F)
}

CD14_Mon <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_CD14_Mon_DE.txt") %>% .[.$change_V1 == "UP",]
CD16_Mon <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_CD16_Mon_DE.txt") %>% .[.$change_V1 == "UP",]
In_termed_Mon <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_In-termed_Mon_DE.txt") %>% .[.$change_V1 == "UP",]
M1 <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_M1_DE.txt") %>% .[.$change_V1 == "UP",]
M2 <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_M2_DE.txt") %>% .[.$change_V1 == "UP",]
Megakaryocyte <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_Megakaryocyte_DE.txt") %>% .[.$change_V1 == "UP",]
Granulocyte <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_Granulocyte_DE.txt") %>% .[.$change_V1 == "UP",]
Dendritic_cell <- read.delim("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_Dendritic_cell_DE.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
CD14_Mon1 <- CD14_Mon[CD14_Mon$gene %in% Reduce(setdiff,list(CD14_Mon$gene,CD16_Mon$gene,In_termed_Mon$gene,M1$gene,M2$gene,Megakaryocyte$gene,Granulocyte$gene,Dendritic_cell$gene)),]
CD16_Mon1 <- CD16_Mon[CD16_Mon$gene %in% Reduce(setdiff,list(CD16_Mon$gene,CD14_Mon$gene,In_termed_Mon$gene,M1$gene,M2$gene,Megakaryocyte$gene,Granulocyte$gene,Dendritic_cell$gene)),]
In_termed_Mon1 <- In_termed_Mon[In_termed_Mon$gene %in% Reduce(setdiff,list(In_termed_Mon$gene,CD14_Mon$gene,CD16_Mon$gene,M1$gene,M2$gene,Megakaryocyte$gene,Granulocyte$gene,Dendritic_cell$gene)),]
M11 <- M1[M1$gene %in% Reduce(setdiff,list(M1$gene,CD14_Mon$gene,CD16_Mon$gene,In_termed_Mon$gene,M2$gene,Megakaryocyte$gene,Granulocyte$gene,Dendritic_cell$gene)),]
M21 <- M2[M2$gene %in% Reduce(setdiff,list(M2$gene,CD14_Mon$gene,CD16_Mon$gene,In_termed_Mon$gene,M1$gene,Megakaryocyte$gene,Granulocyte$gene,Dendritic_cell$gene)),]
Megakaryocyte1 <- Megakaryocyte[Megakaryocyte$gene %in% Reduce(setdiff,list(Megakaryocyte$gene,CD14_Mon$gene,CD16_Mon$gene,In_termed_Mon$gene,M1$gene,M2$gene,
                                                                            Granulocyte$gene,Dendritic_cell$gene)),]
Granulocyte1 <- Granulocyte[Granulocyte$gene %in% Reduce(setdiff,list(Granulocyte$gene,CD14_Mon$gene,CD16_Mon$gene,In_termed_Mon$gene,M1$gene,M2$gene,
                                                                      Megakaryocyte$gene,Dendritic_cell$gene)),]
Dendritic_cell1 <- Dendritic_cell[Dendritic_cell$gene %in% Reduce(setdiff,list(Dendritic_cell$gene,CD14_Mon$gene,CD16_Mon$gene,In_termed_Mon$gene,M1$gene,M2$gene,
                                                                               Granulocyte$gene,Megakaryocyte$gene)),]

diff_gene <- c(CD14_Mon1$gene,In_termed_Mon1$gene,CD16_Mon1$gene,M11$gene,M21$gene,Megakaryocyte1$gene,Granulocyte1$gene,Dendritic_cell1$gene)

celltype.lists <- c("CD14_Mon","In-termed_Mon","CD16_Mon","M1","M2","Megakaryocyte","Granulocyte","Dendritic_cell")
diff_gene_overlap_FC <- data.frame()
for (j in celltype.lists){
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/V4/sub_celltype/Myeloid_Cell/DE/tissue/Mesenteric_lymph_",j,"_DE.txt"), header=T)
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
  diff_gene1$celltype <- rep(j)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}
diff_gene_overlap_FC <- diff_gene_overlap_FC[!grepl("ENSMMUG", diff_gene_overlap_FC$gene),]
overlap_FC <- dcast(diff_gene_overlap_FC,celltype~diff_gene_overlap_FC$gene,value.var = "avg_log2FC",fill = "0")
overlap_FC1 <- data.frame(t(overlap_FC));colnames(overlap_FC1) <- overlap_FC1[1,];overlap_FC1 <- overlap_FC1[-1,]
overlap_FC1 <- overlap_FC1[diff_gene,]
overlap_FC <- melt(overlap_FC,id.vars = "celltype",variable.name='gene',value.name="avg_log2FC")
colnames(overlap_FC) <- c("celltype","gene","avg_log2FC")
str(overlap_FC)
overlap_FC$avg_log2FC <- as.numeric(overlap_FC$avg_log2FC)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = celltype.lists)

overlap_FC[overlap_FC > 1] <- 1
overlap_FC[overlap_FC < -1] <- -1

p <- ggplot(overlap_FC,aes(x=celltype,y=gene,fill=avg_log2FC))+
  geom_raster()+
  scale_fill_gradientn(colors = c("#154889","#1A58A7","#4592E5","#eff3ff","#F76C6C","#F70000","#981216"))+
  # scale_fill_gradient2(low="#4C99AD",mid="#eff3ff",high ="#D86080",midpoint = 0)+
  theme_bw() + 
  labs(title = "Myeloid_Cell_old_vs_young_UP_setdiff_merge")+ 
  theme(legend.title = element_text(size = 20),legend.text = element_text(size = 20),
        panel.grid = element_blank(),plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 20,angle = 90,vjust = 1,hjust = 1),axis.text.y = element_blank(),
        axis.ticks = element_blank(),panel.border = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer/figureS5_Mesenteric_lymph_DEG.pdf",plot = p,width = 10,height = 12)
dev.off()

# Tcm细胞数量 -----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/data/all_sample.combined.celltype.RData")
table(all_sample.combined@meta.data$Tissue)

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
Cell.combined$Tissue_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$Age)
Cell.combined$celltype_Age <- paste0(Cell.combined$celltype1,"_",Cell.combined$Age)
Cell.combined$Tissue_celltype <- paste0(Cell.combined$Tissue,"_",Cell.combined$celltype1)
Cell.combined$Tissue_celltype_Age <- paste0(Cell.combined$Tissue,"_",Cell.combined$celltype1,"_",Cell.combined$Age)
ncol(Cell.combined)
table(Cell.combined@meta.data$Tissue)
table(Cell.combined@meta.data$Tissue_Age)
table(Cell.combined@meta.data$celltype1)
table(Cell.combined@meta.data$celltype_Age)
table(Cell.combined@meta.data$Tissue_celltype)
table(Cell.combined@meta.data$Tissue_celltype_Age)
