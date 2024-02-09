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

# 重新聚类分簇-----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/NK/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
all_cell <- Cell.combined@meta.data[,c(4,5,18)]
all_cell %<>% {.$id<-rownames(.);.}
all_cell$id <-gsub("-",".",all_cell$id)
all_cell$id <-gsub("Lymph_node","Mesenteric_lymph",all_cell$id)

load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "integrated"
unique(all_sample.combined@meta.data$celltype)
Cell.combined<-subset(all_sample.combined, subset = celltype  %in% "NK")

all_cell1 <- Cell.combined@meta.data[,c(4,5)]
all_cell1 %<>% {.$id<-rownames(.);.} %>% join(.,all_cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell1) <- c("Tissue","Age","id","Celltype_old")
all_cell1[is.na(all_cell1)] <- "NON"
all_cell2 <- all_cell1[all_cell1$Celltype_old == "NON",]
all_cell2 <- separate(all_cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined <- AddMetaData(Cell.combined,all_cell1$Celltype_old,col.name = "Celltype_old")
unique(Cell.combined@meta.data$Celltype_old)

Cell.combined <- RunPCA(Cell.combined, npcs = 30, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 42, reduction = "pca", dims = 1:20)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:20)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.2,0.4,0.8,1.2))
p1 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE)
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE)
p3 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.8", label = TRUE)
p4 <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.1.2", label = TRUE)
table(Cell.combined@meta.data$Age)
p <- p<-(p1|p2)/(p3|p4)
p
Cell.combined$Celltype_old <- factor(Cell.combined$Celltype_old,levels = c("NK1","NK2","NK3","NON"))
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Celltype_old)
colors <- c("#80B1D3","#8DD3C7","#FF348B","grey80")
p <- DimPlot(Cell.combined, reduction = "umap", cols = colors, raster=FALSE,repel = F, label = T)+xlab("") + ylab("") + 
  theme(axis.ticks = element_blank(),axis.line = element_blank(),
        # legend.position = 'none',
        axis.text = element_blank())
p
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.after_cluster.RData")

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.after_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.0.4)
###FeaturePlot
gene <- c("NCAM1","IL7R","FCGR3","IL2RB")
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.1.2',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
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

###umap
rm(list=ls())
celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/celltype.txt")
unique(celltype$celltype)
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.after_cluster.RData")
colnames(Cell.combined@meta.data)
cluster <- Cell.combined@meta.data[,c(1,10)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "celltype")
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.celltype.RData")

# NK UMAP -----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
Cell.combined$celltype <- factor(Cell.combined$celltype,levels = c("NK1","NK2","NK3"))
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$celltype)
unique(Cell.combined@meta.data$celltype)
a <- data.frame(table(Cell.combined@meta.data$celltype))
# colors <- c("#80B1D3","#8DD3C7","#FF348B")
colors <- c("#FF348B","#8DD3C7","#1084D3")
show_col(colors)

p <- DimPlot(Cell.combined, reduction = "umap", pt.size = 0.1,group.by = "celltype", 
             cols = colors,
             label.size = 4,label = F,repel = T)+
  xlab("") + ylab("")+ggtitle("")+
  theme(legend.position = 'none',axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS2/figureS2A.tiff",plot = p,width = 4,height = 4)

# NK featureplot ----------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.celltype.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$celltype)
###FeaturePlot
gene <- c("NCAM1","IL7R","FCGR3","IL2RB")
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
data$Cluster <- factor(data$Cluster,levels = rev(c("NK1","NK2","NK3")))
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,3))+
  scale_color_gradientn(colors = c("grey","grey","#eff3ff","#F76C6C","#981216"))+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=11, color="black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS2/figureS2B.pdf",plot = p,width = 4,height = 2)

# 百分比(tissue_celltype) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
unique(Cell.combined@meta.data$Tissue)
celltype_tissue <- Cell.combined@meta.data[,c(4,5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)
young_data <- celltype_tissue[which(celltype_tissue$Age == "young"),]
old_data <- celltype_tissue[which(celltype_tissue$Age == "old"),]

cell_list <- unique(celltype_tissue$celltype)
freq_data <- data.frame(Age = "Age",celltype = "celltype",Freq = "Freq",cell_freq = "cell_freq")
for (j in cell_list){
  young_freq_data <- data.frame(Age = "young",celltype = j,Freq = as.character(nrow(young_data[which(young_data$celltype == j),])),
                                cell_freq = as.character(nrow(young_data[which(young_data$celltype == j),])/nrow(young_data)))
  old_freq_data <- data.frame(Age = "old",celltype = j,Freq = as.character(nrow(old_data[which(old_data$celltype == j),])),
                              cell_freq = as.character(nrow(old_data[which(old_data$celltype == j),])/nrow(old_data)))
  freq_data1 <- rbind(young_freq_data,old_freq_data)
  freq_data <- rbind(freq_data,freq_data1)
}
freq_data <- freq_data[-1,]
###分组堆叠柱状图
freq_data$cell_freq <- as.numeric(freq_data$cell_freq)
colors <- c("#80B1D3","#8DD3C7","#FF348B")
freq_data
p <- ggplot(data=freq_data, aes(x=Age, y=100 * cell_freq, fill = celltype, width=0.8))+scale_x_discrete(limits=rev)+
  geom_bar(stat="identity", position = "fill", width = 0.4, size = 0.25) +
  labs(x = '', y = 'cell_proportion')+
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "")+ 
  theme(panel.grid = element_blank()) +
  xlab("") + ylab("")
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figureS2/figureS2C.pdf",plot = p,width = 3,height = 4)
dev.off()

# FC(all_tissue.young_old) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.celltype.RData")
colnames(Cell.combined@meta.data)
unique(Cell.combined@meta.data$Tissue)
celltype_tissue <- Cell.combined@meta.data[,c(5,18)]
celltype_tissue$celltype <- as.character(celltype_tissue$celltype)

cell_list <- unique(celltype_tissue$celltype)
freq_data <- data.frame()
for (j in cell_list){
  cell <- data.frame(table(celltype_tissue[which(celltype_tissue$celltype == j),]))
  freq_data1 <- data.frame(celltype = j,
                           old_cell = cell[which(cell$Age == "old"),"Freq"],
                           young_cell = cell[which(cell$Age == "young"),"Freq"],
                           FC = as.character(cell[which(cell$Age == "old"),"Freq"]/cell[which(cell$Age == "young"),"Freq"]))
  freq_data <- rbind(freq_data,freq_data1)
}
freq_data$old_per <- freq_data$old_cell/sum(freq_data$old_cell)
freq_data$young_per <- freq_data$young_cell/sum(freq_data$young_cell)
freq_data$log2FC <- log2(as.numeric(freq_data$FC))

###fisher检验
fisher_data <- data.frame()
for (celltype in cell_list){
  test_data1 <- freq_data[freq_data$celltype == celltype,c(2,3)]
  test_data2 <- freq_data[freq_data$celltype != celltype,]
  test_data3 <- data.frame(old_cell = sum(test_data2$old_cell),young_cell = sum(test_data2$young_cell))
  test_data4 <- rbind(test_data1,test_data3)
  rownames(test_data4) <- c("cell","others")
  fisher_data1 <- data.frame(celltype = celltype,p.value = fisher.test(test_data4)$p.value)
  fisher_data <- rbind(fisher_data,fisher_data1)
}
fisher_data$p.value1 = as.factor(ifelse(fisher_data$p.value < 0.01,"**",
                                        ifelse(0.01 < fisher_data$p.value & fisher_data$p.value < 0.05,"*", "ns")))



# figureS2 NK DEGs(按celltype做差异） ------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS2/Cell.combined.celltype.RData")
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
                                             ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                    ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                             'NOT'))
      diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
      celltype <- gsub(" ","_",celltype)
      write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/",celltype,"_DE.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
    } 
  }
}

# figureS2 NK不同亚型DEG_setdiff ------------------------------------------------------------
rm(list=ls())
NK1 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/NK1_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK2 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/NK2_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK3 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/NK3_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

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
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/",j,"_DE.findmark.txt"), header=T)
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
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS2/","figureS2D.pdf"),plot = p,width = 5,height = 6)
dev.off()

# figureS2 NK 不同细胞类型_DE_setdiff-富集-------------------------------------------------------------------------
rm(list=ls())
NK1 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/NK1_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK2 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/NK2_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
NK3 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS2/DE/NK3_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

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
write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS2/","NK_UP_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
