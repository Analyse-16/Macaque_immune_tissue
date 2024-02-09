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


# 合并CD4细胞 -----------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD41.RData")
load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_CD4.RData")
cells_CD4 <- c(colnames(Cell.combined_CD4),colnames(Cell.combined_CD41))
load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T.RData")
Cell.combined_CD4 <- subset(Cell.combined, cells = cells_CD4)
save(Cell.combined_CD4, file="/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")

# 添加细胞类型 ------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD4/Cell.combined.celltype.RData")
Cell.combined@meta.data$celltype <- factor(Cell.combined@meta.data$celltype,levels =
                                             c("Naive","Tcm","Tem","TRM","CTL"))
library(RColorBrewer)
colors <- c("#B73532","#8DD3C7","#FCB8DB","#80B1D3","#FB6655")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype", cols = colors,label = F,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p

colnames(Cell.combined@meta.data)
all_cell <- Cell.combined@meta.data[,c(4,5,18)]
all_cell %<>% {.$id<-rownames(.);.}
all_cell$id <-gsub("-",".",all_cell$id)
all_cell$id <-gsub("Lymph_node","Mesenteric_lymph",all_cell$id)

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")
colnames(Cell.combined_CD4@meta.data)
all_cell1 <- Cell.combined_CD4@meta.data[,c(4,5)]
all_cell1 %<>% {.$id<-rownames(.);.} %>% join(.,all_cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell1) <- c("Tissue","Age","id","Celltype_old")
str(all_cell1)
all_cell1$Celltype_old <- as.character(all_cell1$Celltype_old)
all_cell1[is.na(all_cell1)] <- "NON"
all_cell2 <- all_cell1[all_cell1$Celltype_old == "NON",]
all_cell2 <- separate(all_cell2,id,into = c("V1","V2"),sep = "\\.",remove = F)
table(all_cell2$V2)
Cell.combined_CD4 <- AddMetaData(Cell.combined_CD4,all_cell1$Celltype_old,col.name = "Celltype_old")
unique(Cell.combined_CD4@meta.data$Celltype_old)
save(Cell.combined_CD4, file="/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")

# 分簇-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")
Cell.combined <- Cell.combined_CD4
DefaultAssay(Cell.combined) <- "integrated"
unique(Cell.combined@meta.data$celltype)
Cell.combined <- RunPCA(Cell.combined, npcs = 20, verbose = T)
Cell.combined <- RunUMAP(Cell.combined, seed.use = 42, reduction = "pca", dims = 1:10)
Cell.combined <- FindNeighbors(Cell.combined, reduction = "pca", dims = 1:10)
Cell.combined <- FindClusters(Cell.combined, resolution = c(0.4,0.8,1.2,1.5))
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "integrated_snn_res.0.8", label = TRUE)
p
Cell.combined@meta.data$Celltype_old <- factor(Cell.combined@meta.data$Celltype_old,levels =
                                             c("Naive","Tcm","Tem","TRM","NON"))
colors <- c("#B73532","#8DD3C7","#FCB8DB","#80B1D3","#FB6655","black")
p <- DimPlot(Cell.combined, reduction = "umap", group.by = "Celltype_old", cols = colors,label = F,label.size = 6)+
  theme(legend.text=element_text(colour= 'black',size=20))
p
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")

# test
for (i in c(10:30)){
  Cell.combined <- RunUMAP(Cell.combined, seed.use = 42,reduction = "pca", dims = 1:i)
  p <- DimPlot(Cell.combined, reduction = "umap", group.by = "Celltype_old",pt.size = 0.01,
               label = T, repel = TRUE, label.size = 5,raster=F)+
    theme(legend.position = "right")
  p
  ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/round_4000/",i,".png"),p,width = 10,height = 8)
}

# FeaturePlot-------------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$integrated_snn_res.0.8)
gene <- c("CCR7","CD69","CCR6","ITGAE","EOMES","RUNX3")
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'integrated_snn_res.0.4',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
# data$Cluster <- factor(data$Cluster,levels = c(0:17))
p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradient2(low="#330066",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  # scale_color_manual(values=color16)+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text=element_text(size=11, color="black")) + RotatedAxis()
p
dev.off()

# figureS3A CD4 添加细胞类型信息 --------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4_cluster.RData")
celltype <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/celltype.txt")
colnames(Cell.combined@meta.data)
cluster <- Cell.combined@meta.data[,c(1,21)]
colnames(cluster) <- c("orig.ident","cluster")
cluster_celltype <- join(cluster,celltype)
Cell.combined <- AddMetaData(Cell.combined,cluster_celltype$celltype,col.name = "celltype")
Cell.combined@meta.data$Age <- factor(Cell.combined@meta.data$Age,levels = c("young","old"))
Cell.combined@meta.data$celltype <- factor(Cell.combined@meta.data$celltype,levels =
                                                 c("Naive","Tcm","Tem","CTL/TRM"))
# colors <- c("#B73532","#8DD3C7","#FCB8DB","#80B1D3","#FB6655")
colors <- c("#80B1D3","#BD8AD6","#A5DF72","#B73532")
p1 <- DimPlot(Cell.combined, reduction = "umap",cols = colors,group.by = "celltype",
              split.by = "Age",ncol = 2,
              pt.size = 0.01,label = F)
p1
save(Cell.combined, file="/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4.RData")

#  figureS3A celltype_umap-----------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4.RData")
Cell.combined@meta.data$celltype <- factor(Cell.combined@meta.data$celltype,levels =
                                             c("Naive","Tcm","Tem","CTL/TRM"))
colors <- c(colors <- c("#80B1D3","#BD8AD6","#A5DF72","#DB6664"))
p <- DimPlot(Cell.combined, reduction = "umap",cols = colors,group.by = "celltype",
              pt.size = 0.01,label = F)+
  ggtitle("")+
  theme(legend.position = 'none',axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/","figureS3A.tiff"),plot = p,width = 4,height = 4)

# features_Bubble_celltype-----------------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4.RData")
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$celltype)
gene <- c("CCR7","CD69","CCR6","ITGAE","EOMES","RUNX3")
p1<-DotPlot(object = Cell.combined, features = gene, group.by  = 'celltype',col.min=-2,col.max = 2,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
data <- p1$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
unique(data$Cluster)
data$Cluster <- factor(data$Cluster,levels = rev(c("Naive","Tcm","Tem","CTL/TRM")))

p<-ggplot(data,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradientn(colors = c("grey","grey","#eff3ff","#F76C6C","#981216"))+
  # scale_color_gradient2(low="#e0e0e0",mid="#e5f5f9",high ="#ef3b2c",midpoint = 0)+
  theme_classic()+labs(x = '', y = '')+
  theme(axis.text.x=element_text(size=15, angle = 90,color="black"),
        axis.text.y=element_text(size=15, color="black")) + RotatedAxis()
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/","figureS3B.pdf"),plot = p,width = 5.4,height = 3)

# 百分比(celltype-柱状图) ---------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4.RData")
colnames(Cell.combined@meta.data)
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
freq_data$celltype <- factor(freq_data$celltype,levels = c("Naive","Tcm","Tem","CTL/TRM"))
colors <- c("#80B1D3","#BD8AD6","#A5DF72","#DB6664")
show_col(colors)
p <- ggplot(data=freq_data, aes(x=Age, y=100 * cell_freq, fill = celltype, width=0.8))+scale_x_discrete(limits=rev)+
  geom_bar(stat="identity", position = "fill", width = 0.4, size = 0.25) +
  labs(x = '', y = 'cell_proportion')+
  scale_fill_manual(values = colors)+
  theme_bw() + 
  labs(title = "")+ 
  theme(panel.grid = element_blank()) +
  xlab("") + ylab("")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/","figureS3C.pdf"),plot = p,width = 3,height = 4)
dev.off()

# DEGs(按celltype做差异） ------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS3/Cell.combined_CD4.RData")
colnames(Cell.combined@meta.data)
DefaultAssay(Cell.combined) <- "RNA"
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$Age)
data.list_celltype <- SplitObject(Cell.combined, split.by = "celltype")
celltype.lists <- unique(Cell.combined@meta.data$celltype)
for (j in celltype.lists){
  data_celltype <- data.list_celltype[[j]]
  d <- data.frame(table(data_celltype$Age))
  if(nrow(d) == 1){
    print(paste0(tissue,"_",j," no young or old cells"))
  }
  else{
    a = d[1,2];b = d[2,2]
    if(a<3 | b <3){
      print(paste0(tissue,"_",j,"Cells fewer than 3 cells")) 
    }else{
      diff_gene <- FindMarkers(data_celltype,ident.1 = "old", ident.2 = "young",logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
      diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
      diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                             ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                    ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                             'NOT'))
      diff_gene <- diff_gene[diff_gene$change_V1 != "NOT",]
      diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
      j <- gsub("/","_",j)
      write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/",j,"_DE.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
    } 
  }
}

# figureS3 CD4 不同亚型DEG_setdiff ------------------------------------------------------------
rm(list=ls())
Naive <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/Naive_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tem <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/Tem_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tcm <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/Tcm_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
CTL <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/CTL_TRM_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
library(reshape2)
Naive1 <- Naive[Naive$gene %in% Reduce(setdiff,list(Naive$gene,Tem$gene,Tcm$gene,CTL$gene)),]
Naive1$celltype <- rep("Naive1")
Tem1 <- Tem[Tem$gene %in% Reduce(setdiff,list(Tem$gene,Naive$gene,Tcm$gene,CTL$gene)),]
Tem1$celltype <- rep("Tem1")
Tcm1 <- Tcm[Tcm$gene %in% Reduce(setdiff,list(Tcm$gene,Naive$gene,Tem$gene,CTL$gene)),]
Tcm1$celltype <- rep("Tcm1")
CTL1 <- CTL[CTL$gene %in% Reduce(setdiff,list(CTL$gene,Naive$gene,Tem$gene,Tcm$gene)),]
CTL1$celltype <- rep("CTL1")

diff_gene <- c(Naive1$gene,Tcm1$gene,Tem1$gene,CTL1$gene)

celltype.lists <- c("Naive","Tcm","Tem","CTL_TRM")
diff_gene_overlap_FC <- data.frame()
for (j in celltype.lists){
  old_vs_young <- read.delim(file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/",j,"_DE.findmark.txt"), header=T)
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

overlap_FC$celltype <- factor(overlap_FC$celltype,levels = c("Naive","Tcm","Tem","CTL_TRM"))
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
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/","figureS3D.pdf"),plot = p,width = 4,height = 4)
dev.off()

# figureS3 CD4 不同细胞类型_DE_setdiff-富集-------------------------------------------------------------------------
rm(list=ls())
Naive <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/Naive_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tem <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/Tem_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
Tcm <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/Tcm_DE.findmark.txt") %>% .[.$change_V1 == "UP",]
CTL <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figureS3/DE/CTL_TRM_DE.findmark.txt") %>% .[.$change_V1 == "UP",]

library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

Naive1 <- Naive[Naive$gene %in% Reduce(setdiff,list(Naive$gene,Tem$gene,Tcm$gene,CTL$gene)),]
colnames(Naive1)[6] = "SYMBOL"
DE <- merge(x=Naive1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_Naive1<-as.data.frame(ego_BP1@result)
ego_BP_Naive1$celltype <- rep("Naive")
Tem1 <- Tem[Tem$gene %in% Reduce(setdiff,list(Tem$gene,Naive$gene,Tcm$gene,CTL$gene)),]
colnames(Tem1)[6] = "SYMBOL"
DE <- merge(x=Tem1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_Tem1<-as.data.frame(ego_BP1@result)
ego_BP_Tem1$celltype <- rep("Tem")
Tcm1 <- Tcm[Tcm$gene %in% Reduce(setdiff,list(Tcm$gene,Naive$gene,Tem$gene,CTL$gene)),]
colnames(Tcm1)[6] = "SYMBOL"
DE <- merge(x=Tcm1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_Tcm1<-as.data.frame(ego_BP1@result)
ego_BP_Tcm1$celltype <- rep("Tcm")
CTL1 <- CTL[CTL$gene %in% Reduce(setdiff,list(CTL$gene,Naive$gene,Tem$gene,Tcm$gene)),]
colnames(CTL1)[6] = "SYMBOL"
DE <- merge(x=CTL1,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_CTL1<-as.data.frame(ego_BP1@result)
ego_BP_CTL1$celltype <- rep("CTL")

GO <- rbind(ego_BP_Naive1,ego_BP_Tcm1,ego_BP_Tem1,ego_BP_CTL1)
GO <- GO[GO$pvalue < 0.05,]
write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figureS3/","CD4_UP_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
