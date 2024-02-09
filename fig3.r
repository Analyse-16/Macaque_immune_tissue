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

#  figure3B CD8 （ggplot）(3A,4A代码在S4图）-------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
Cell.combined_CD8<-AddMetaData(Cell.combined_CD8,Cell.combined_CD8@reductions$umap@cell.embeddings,col.name = colnames(Cell.combined_CD8@reductions$umap@cell.embeddings))
colnames(Cell.combined_CD8@meta.data)

Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype)
DefaultAssay(Cell.combined_CD8) <- "RNA"
data.list_age <- SplitObject(Cell.combined_CD8, split.by = "Age")
genes <- c("GZMB","CX3CR1","ADGRG1")
# genes <- c("CCL5","RGS9","ZEB2")
young_data <- data.list_age[["young"]]
mat1=as.data.frame(young_data[["RNA"]]@data["GZMB",])

colnames(mat1)="exp"
mat2=Embeddings(young_data,"umap")
mat3=merge(mat2,mat1,by="row.names")
mat3$exp[mat3$exp > 4] = 4
colors=rev(c("#A10000","#B30000","#CE1212","#F03B3B","#F05A5A","#F07979","#FFAFAF20"))
show_col(colors)
p1 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=exp),size = 0.01)+
  scale_color_gradientn(limits = c(0,4),breaks = c(0,2,4),colors=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p1

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
  p1 <- p1/p12
}
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure3/figure.3B_young1.pdf",plot = p1,width = 6,height = 17)

old_data <- data.list_age[["old"]]
mat1=as.data.frame(old_data[["RNA"]]@data["GZMB",])

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
  p2 <- p2/p22
}
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure3/figure.3B_old1.pdf",plot = p2,width = 6,height = 17)

# figure3C CD8 36个TF的分析(差异分析)(tcm重新聚类在图4) ----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
DefaultAssay(Cell.combined_CD8) <- "RNA"
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype1)
p <- DimPlot(Cell.combined_CD8, reduction = "umap", group.by = "celltype1", label = TRUE,label.size = 8)+
  theme(legend.text=element_text(colour= 'black',size=20))
p
unique(Cell.combined_CD8@meta.data$celltype1)
data.list_age <- SplitObject(Cell.combined_CD8, split.by = "Age")
###young
data_young <- data.list_age[["young"]]
diff_gene_young <- FindMarkers(data_young, ident.1 = c("CTL/TRM","TRM","Tex","XCL1high Tcm"), logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene_young <- diff_gene_young  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene_young$change = as.factor(ifelse(abs(diff_gene_young$avg_log2FC) >= 1,'Diff','NOT' ))
diff_gene_young$change_V1 = as.factor( ifelse( diff_gene_young$p_val < 0.05,
                                               ifelse( diff_gene_young$avg_log2FC > -0.25, 
                                                       ifelse( diff_gene_young$avg_log2FC >= 0.25, 'UP', 'NOT' ),
                                                       'DOWN' ),
                                               'NOT' ))
diff_gene_young <- diff_gene_young[order(diff_gene_young$avg_log2FC,decreasing = T),]
write.table(diff_gene_young,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","CD8","_cluster_young.txt"),quote = FALSE,sep = "\t",row.names = F)
###old
data_old <- data.list_age[["old"]]
diff_gene_old <- FindMarkers(data_old, ident.1 = c("CTL/TRM","TRM","Tex","XCL1high Tcm"), logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene_old <- diff_gene_old  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene_old$change = as.factor(ifelse(abs(diff_gene_old$avg_log2FC) >= 1,'Diff','NOT' ))
diff_gene_old$change_V1 = as.factor( ifelse( diff_gene_old$p_val < 0.05,
                                             ifelse( diff_gene_old$avg_log2FC > -0.25, 
                                                     ifelse( diff_gene_old$avg_log2FC >= 0.25, 'UP', 'NOT' ),
                                                     'DOWN' ),
                                             'NOT' ))
diff_gene_old <- diff_gene_old[order(diff_gene_old$avg_log2FC,decreasing = T),]
write.table(diff_gene_old,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","CD8","_cluster_old.txt"),quote = FALSE,sep = "\t",row.names = F)
###old/young
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
Idents(Cell.combined_CD8)
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype1)
Cell.combined_CD8 <- subset(Cell.combined_CD8, subset = (celltype1 %in% c("CTL/TRM","TRM","Tex","XCL1high Tcm")))
DefaultAssay(Cell.combined_CD8) <- "RNA"
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$Age)
diff_gene <- FindMarkers(Cell.combined_CD8, ident.1 = "old", ident.2 = "young", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change = as.factor(ifelse(abs(diff_gene$avg_log2FC) >= 1,'Diff','NOT' ))
diff_gene$change_V1 = as.factor( ifelse( diff_gene$p_val < 0.05,
                                         ifelse( diff_gene$avg_log2FC > -0.25, 
                                                 ifelse( diff_gene$avg_log2FC >= 0.25, 'UP', 'NOT' ),
                                                 'DOWN' ),
                                         'NOT' ))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","CD8","_cluster_old_vs_young.txt"),quote = FALSE,sep = "\t",row.names = F)

#  figure3C CD8 36个TF的分析(pySCENIC) ############# -----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype1)
Cell.combined_CD8_sel <- subset(Cell.combined_CD8, subset = (celltype1 %in% c("CTL/TRM","TRM","Tex","XCL1high Tcm")))
countData <- t(as.matrix(Cell.combined_CD8_sel@assays[["RNA"]]@counts))
write.csv(countData,paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/pyscenic/","amnionCountMatrix.csv"))

# /data1/zhuzn/wangsn/reviewer_4000/figure3/pyscenic/pyscenic.sh

# figure3C CD8 36个TF的分析 ----------------------------------------------------------------
### AUC热图
rm(list=ls())
diff_young <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_young.txt")
diff_young <- diff_young[diff_young$change_V1 == "UP",]
# diff_young <- diff_young[diff_young$change_V1 != "NOT",]

diff_old <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_old.txt")
diff_old <- diff_old[diff_old$change_V1 == "UP",]
# diff_old <- diff_old[diff_old$change_V1 != "NOT",]

diff_old_young <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_old_vs_young.txt")
diff_old_young <- diff_old_young[diff_old_young$change_V1 == "UP",]
# diff_old_young <- diff_old_young[diff_old_young$change_V1 != "NOT",]

regulonAUC <- read.csv("/data1/zhuzn/wangsn/reviewer_4000/figure3/pyscenic/AmnionStep3Result.csv", row.names=1)
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
pdf("/data1/zhuzn/wangsn/reviewer_4000/figure3/figure3C.pdf",height = 4.5,width = 5)
grid.draw(venn)
dev.off()

# figure3D overlap TF targets FC 富集（按上下调分开） --------------------------------------------------
rm(list=ls())
diff_old_young <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/CD8_cluster_old_vs_young.txt")
diff_old_young <- diff_old_young[!grepl("ENSMMUG", diff_old_young$gene),]
diff_old_young <- diff_old_young[abs(diff_old_young$avg_log2FC) > 0.25,]
diff_old_young_UP <- diff_old_young[diff_old_young$change_V1 == "UP",c(2,6)]
diff_old_young_DOWN <- diff_old_young[diff_old_young$change_V1 == "DOWN",c(2,6)]

targets_all <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/pyscenic/AmnionStep1Result.tsv")
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
write.table(targets,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","TF_targets_DE_FC.txt"),quote = FALSE,sep = "\t",row.names = F)
write.table(targets_UP$gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","TF_targets_DE_UP.txt"),quote = FALSE,sep = "\t",row.names = F)
write.table(targets_DOWN$gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","TF_targets_DE_DOWN.txt"),quote = FALSE,sep = "\t",row.names = F)

###富集
rm(list = ls())
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))
targets <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","TF_targets_DE_FC.txt"))
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
write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","BHLHE40_GO_BP_UP.txt"), quote = FALSE,sep="\t",row.names = FALSE)
# write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/","BHLHE40_GO_BP_DOWN.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# figure3D overlap TF targets FC 富集结果图（按上下调分开） --------------------------------------------------
rm(list=ls())
BHLHE40_GO_BP_UP <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/BHLHE40_GO_BP_UP.txt")
BHLHE40_GO_BP_UP$Group <- rep("UP")
BHLHE40_GO_BP_DOWN <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure3/BHLHE40_GO_BP_DOWN.txt")
BHLHE40_GO_BP_DOWN$Group <- rep("DOWN")
GO <- rbind(BHLHE40_GO_BP_UP,BHLHE40_GO_BP_DOWN)
GO_lists <- c("negative regulation of cellular protein metabolic process","inflammatory response",
        "regulation of interleukin-6 production","cellular response to molecule of bacterial origin",
        "I-kappaB kinase/NF-kappaB signaling",
        "T cell differentiation","T-helper cell differentiation",
        "alpha-beta T cell activation involved in immune response","alpha-beta T cell differentiation",
        "positive regulation of neutrophil apoptotic process")
dat <- GO[GO$Description %in% GO_lists,]

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
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure3/figure.3D_GO.pdf",plot = p,width = 6,height = 6)

# 细胞互作测试 CellChat------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
Cell.combined <- Cell.combined_CD8
# load(file = "/data1/zhuzn/wangsn/V4/sub_celltype/T Cell/CD8/Cell.combined.celltype.RData")
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype1",label = TRUE)
p2
colnames(Cell.combined@meta.data)
cells_data <- Cell.combined@meta.data[,c(4,5,22)]
unique(cells_data$celltype1)
cells_data$celltype <- "CD8 T"
cells_data[which(cells_data$celltype1 %in% c("Tcm_2","CTL/TRM","Tex","TRM")),4] <- "CD8 Taa"

load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
# load(file = "/data1/zhuzn/wangsn/V8/all_sample.combined.celltype.RData")
colnames(all_sample.combined@meta.data)
unique(all_sample.combined@meta.data$celltype)
all_cell <- all_sample.combined@meta.data[,c(4,5,18)]

cell <- cells_data
cell  %<>% {.$id<-rownames(.);.} 

all_cell %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell) <- c("Tissue","Age","celltype","id","Celltype")
all_cell[is.na(all_cell)] <- "NON"
all_cell$celltype <- as.character(all_cell$celltype)
all_cell$Celltype <- ifelse(all_cell$Celltype == "NON" ,all_cell$celltype,all_cell$Celltype)
all_cell[which(all_cell$Celltype == "T Cell"),5] <- "CD4 T"

unique(all_cell$Celltype)
all_sample.combined <- AddMetaData(all_sample.combined,all_cell$Celltype,col.name = "celltype")
table(all_sample.combined@meta.data$celltype)

data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
Cell.combined <- data.list_age[["old"]]
data.list_tissue <- SplitObject(Cell.combined, split.by = "Tissue")

for (t in unique(Cell.combined$Tissue)){
  scRNA <- data.list_tissue[[t]]
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
  saveRDS(cellchat, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/",t,"_Taa_old.rds"))
}

# CellChat结果-------------------------------------------------------------------------
rm(list=ls())
library(CellChat)
library(patchwork)
library(tidyverse)
library(ggalluvial)
t <- c("PBMC","Spleen","Mesenteric_lymph","Bone_marrow")
cellchat <- readRDS(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Mesenteric_lymph","_Taa_old.rds"))

groupSize <- as.numeric(table(cellchat@idents))
mat <- cellchat@net$weight
i = 4
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[i, ] <- mat[i, ]
colnames(mat)
p <- netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat2)[i])
print(p)
dev.off()

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

pathways <- cellchat@netP$pathways
pathways.show <- "MIF"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

netVisual_bubble(cellchat,
                 # sources.use = c(4),
                 targets.use = c(4),
                 # signaling = pathways,
                 remove.isolate = F)

# 细胞互作测试cellphoneDB ------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
Cell.combined <- Cell.combined_CD8
p2 <- DimPlot(Cell.combined, reduction = "umap", group.by = "celltype1",label = TRUE)
p2
colnames(Cell.combined@meta.data)
cells_data <- Cell.combined@meta.data[,c(4,5,22)]
unique(cells_data$celltype1)
cells_data$celltype <- "CD8 T"
cells_data[which(cells_data$celltype1 %in% c("Tcm_2","CTL/TRM","Tex","TRM")),4] <- "CD8 Taa"

load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
colnames(all_sample.combined@meta.data)
unique(all_sample.combined@meta.data$celltype)
all_cell <- all_sample.combined@meta.data[,c(4,5,18)]

cell <- cells_data
cell  %<>% {.$id<-rownames(.);.} 

all_cell %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell) <- c("Tissue","Age","celltype","id","Celltype")
all_cell[is.na(all_cell)] <- "NON"
all_cell$celltype <- as.character(all_cell$celltype)
all_cell$Celltype <- ifelse(all_cell$Celltype == "NON" ,all_cell$celltype,all_cell$Celltype)
all_cell[which(all_cell$Celltype == "T Cell"),5] <- "CD4 T"

unique(all_cell$Celltype)
all_sample.combined <- AddMetaData(all_sample.combined,all_cell$Celltype,col.name = "celltype")
table(all_sample.combined@meta.data$celltype)
# cellphoneDB
data.list_age <- SplitObject(all_sample.combined, split.by = "Age")
Cell.combined <- data.list_age[["old"]]
data.list_tissue <- SplitObject(Cell.combined, split.by = "Tissue")

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
for (t in unique(Cell.combined$Tissue)){
  dir.create(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/cellphoneDB/",t))
  scRNA <- data.list_tissue[[t]]
  data <- as_matrix(scRNA@assays$RNA@data)
  write.table(data,paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/cellphoneDB/",t,"/cellphonedb_count.txt"), sep='\t', quote=F)
  meta_data <- cbind(rownames(scRNA@meta.data), scRNA@meta.data[,'celltype', drop=F])  
  meta_data <- as.matrix(meta_data)
  meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
  write.table(meta_data,paste0("/data1/zhuzn/wangsn/reviewer_4000/figure3/cellphoneDB/",t,"/cellphonedb_meta.txt"), sep='\t', quote=F, row.names=F)
}

# cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --counts-data=gene_name
# cellphonedb 自己的绘图
# cellphonedb plot dot_plot
# cellphonedb plot heatmap_plot cellphonedb_meta_young.txt

# 网络图-------------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(scales)
# install.packages("/data1/zhuzn/soft/igraph_1.1.2.tar.gz", repos=NULL, type = "source")
library(igraph)
pvalues=read.table("/data1/zhuzn/wangsn/V4/cellphoneDB/young/out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
# pvalues=read.table("/data1/zhuzn/wangsn/V4/cellphoneDB/aging/out/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)

pvalues=pvalues[,12:dim(pvalues)[2]]
statdf=as.data.frame(colSums(pvalues < 0.05))
colnames(statdf)=c("number")
statdf$index=str_replace(rownames(statdf),"ell\\.","ell ")
statdf$index=str_replace(statdf$index,"e\\.","e ")
statdf$index=str_replace(statdf$index,"Fibroblast.","Fibroblast ")
statdf$index=str_replace(statdf$index,"CMP.","CMP ")
statdf$index=str_replace(statdf$index,"NK.","NK ")
statdf <- separate(statdf, index,into = c("indexa","indexb"),sep = " ",remove = T)
rankname=sort(unique(statdf$indexa)) 
A=c()
B=c()
C=c()
remaining=rankname
for (i in rankname[-11]) {
  remaining=setdiff(remaining,i)
  for (j in remaining) {
    count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
      statdf[statdf$indexb == i & statdf$indexa == j,"number"]
    A=append(A,i)
    B=append(B,j)
    C=append(C,count)
  }
}
statdf2=data.frame(indexa=A,indexb=B,number=C)
statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
statdf2=statdf2[statdf2$number > 0,] #过滤掉值为0的观测
#设置节点和连线的颜色
color1=c("#8DD3C7", "#FDB462", "#B3DE69", "#FCCDE5", "#fb8072", "#BC80BD", "#377eb8", "#4daf4a", "#984ea3","#a65628","#ff7f00","#e7298a")
show_col(color1)
dev.off()
names(color1)=rankname
color2=colorRampPalette(brewer.pal(9, "Reds")[3:7])(40) #将颜色分成多少份，取决于互作关系数目的最大值
names(color2)=1:40 #每一份颜色用对应的数字命名
#做网络图
##下面的四行代码相对固定
net <- graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
group <-  cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(group)))

E(net)$width <- E(net)$number / 2 #将数值映射到连线的宽度，有时还需要微调，这里除以2就是这个目的
E(net)$color <- color2[as.character(ifelse(E(net)$number > 40,40,E(net)$number))] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
E(net)$label = E(net)$number #连线的标注
E(net)$label.color <- "black" #连线标注的颜色
V(net)$label.color <- "black" #节点标注的颜色
V(net)$color <- color1[names(V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

#调整节点位置的线条角度
##如果没有这两行代码，节点位置的圆圈是向右的
loop.angle<-ifelse(coords[igraph::V(net),1]>0,-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]))
igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
#pdf("interaction.num.3.pdf",width = 6,height = 6)
p <- plot(net,
          edge.arrow.size = 0, #连线不带箭头
          edge.curved = 0, #连线不弯曲
          vertex.frame.color = "black", #节点外框颜色
          layout = coords,
          vertex.label.cex = 1, #节点标注字体大小
          vertex.size = 20) #节点大小
p
ggsave(paste0('./cellphoneDB/young/','net.png'), plot=p, dpi = 600,width = 12, height = 10)
# ggsave(paste0('./cellphoneDB/aging/','net.png'), plot=p, dpi = 600,width = 12, height = 10)
dev.off()

# 网络图1-------------------------------------------------------------------------
rm(list=ls())
# out='/data1/zhuzn/wangsn/V4/cellphoneDB/young/out/' ##  outs 文件放在这里了。
out='/data1/zhuzn/wangsn/V4/cellphoneDB/aging/out/'
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
netf<- "count_network.txt"
mynet <- read.delim(paste0(out,"count_network.txt"), check.names = FALSE)
# write.table(mynet, './cellphoneDB/young/cellphonedb_mynet_young.txt', sep='\t', quote=F,row.names = F)
write.table(mynet, './cellphoneDB/aging/cellphonedb_mynet_old.txt', sep='\t', quote=F,row.names = F)

###合并年轻年老
rm(list=ls())
cellphonedb_young <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/young/cellphonedb_mynet_young.txt")
cellphonedb_young$group <- rep("young")
cellphonedb_aging <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/aging/cellphonedb_mynet_old.txt")
cellphonedb_aging$group <- rep("aging")
cellphonedb_edge <- rbind(cellphonedb_young,cellphonedb_aging)
range(cellphonedb_edge$number)
cellphonedb_edge$number_rm <- as.factor(ifelse(cellphonedb_edge$number > 24 ,"T","F"))
cellphonedb_edge$number_rm1 <- as.character(cellphonedb_edge$number_rm)
cellphonedb_edge$number_rm1[which(cellphonedb_edge$group == "young" & cellphonedb_edge$number_rm == "T")] <- "YOUNG"
cellphonedb_edge$number_rm1[which(cellphonedb_edge$group == "aging" & cellphonedb_edge$number_rm == "T")] <- "OLD"
cellphonedb_edge$number1 <- cellphonedb_edge$number
cellphonedb_edge$number1[which(cellphonedb_edge$number_rm1 == "F")] <- "0"
write.table(cellphonedb_edge, './cellphoneDB/cellphonedb_mynet_merge.txt', sep='\t', quote=F,row.names = F)

### 挑选差异明显的互作关系
rm(list=ls())
cellphonedb_young <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/young/cellphonedb_mynet_young.txt")
cellphonedb_young$group <- rep("young")
cellphonedb_aging <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/aging/cellphonedb_mynet_old.txt")
cellphonedb_aging$group <- rep("aging")
cellphonedb_edge <- cbind(cellphonedb_young,cellphonedb_aging)
cellphonedb_edge <- cellphonedb_edge[,c(1,2,3,8)]
colnames(cellphonedb_edge) <- c("SOURCE","TARGET","young","old")
cellphonedb_edge <- cellphonedb_edge[cellphonedb_edge$old != "0",]
cellphonedb_edge$FC<-ifelse(cellphonedb_edge$old >= cellphonedb_edge$young,(cellphonedb_edge$old/cellphonedb_edge$young),
                            -(cellphonedb_edge$young/cellphonedb_edge$old))
cellphonedb_edge <- cellphonedb_edge[cellphonedb_edge$FC > 1.3 | cellphonedb_edge$FC < -3,]
cellphonedb_edge <- cellphonedb_edge[,-5]
cellphonedb_edge1 <- melt(cellphonedb_edge,value.name="number")
cellphonedb_edge1$FC<-ifelse(cellphonedb_edge$old > cellphonedb_edge$young,"UP","DOWN")
cellphonedb_edge1$FC1<-paste(cellphonedb_edge1$variable,cellphonedb_edge1$FC,sep = "_")
write.table(cellphonedb_edge1, './cellphoneDB/cellphonedb_mynet_FC.txt', sep='\t', quote=F,row.names = F)

rm(list=ls())
# cellphonedb_edge <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/young/cellphonedb_edge.txt")
cellphonedb_edge <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/aging/cellphonedb_edge.txt")
range(cellphonedb_edge$number)
cellphonedb_edge$group <- as.factor(ifelse(cellphonedb_edge$number > 20 ,"T","F"))
colnames(cellphonedb_edge)
cellphonedb_edge <- cellphonedb_edge[,c(2,7)]
colnames(cellphonedb_edge)[1] <- "edge"
cellphonedb_edge1 <- cellphonedb_edge %>% separate(edge,into = c("source","target"),sep = " \\(interacts with\\) ",remove = FALSE)
colnames(cellphonedb_edge1)[2] <- "celltype"
# write.table(cellphonedb_edge1, './cellphoneDB/young/cellphonedb_edge1.txt', sep='\t', quote=F,row.names = F)
write.table(cellphonedb_edge1, './cellphoneDB/aging/cellphonedb_edge1.txt', sep='\t', quote=F,row.names = F)

rm(list=ls())
# mynet <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/young/cellphonedb_mynet_young.txt")
mynet <- read.delim("/data1/zhuzn/wangsn/V4/cellphoneDB/aging/cellphonedb_mynet_old.txt")
range(mynet$count)
#0-57
mynet <- mynet[order(mynet[,1]),]
table(mynet$count)
mynet %>% filter(count>0) -> mynet  # 有零会报错
head(mynet)
net<- graph_from_data_frame(mynet)
plot(net)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # 设置网络布局
E(net)$width  <- E(net)$count/10  # 边点权重（粗细）
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 
net2 <- net

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
}
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7)
dev.off()

# png(paste0('./cellphoneDB/young/','net.png'), width = 500, height = 500)
png(paste0('./cellphoneDB/aging/','net.png'), width = 500, height = 500)
# graphics.off()
par("mar")
par(mar=c(1,1,1,1))
plot(net, edge.arrow.size=.1, 
     edge.curved=0.2, # 只是调了这个参数
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=1)
dev.off()

length(unique(mynet$SOURCE))
par(mfrow=c(3,4), mar=c(.3,.3,.3,.3))
# png(paste0('./cellphoneDB/young/','net1.png'), width = 2500, height = 500)
# png(paste0('./cellphoneDB/aging/','net1.png'), width = 2500, height = 500)
for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  # 故技重施
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.4,
       edge.label = E(net1)$count, # 绘制边的权重
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1
  )
}
dev.off()

# 点图-------------------------------------------------------------------------
rm(list=ls())
t <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
out='/data1/zhuzn/wangsn/reviewer_4000/figure3/cellphoneDB/Spleen/out/'
mypvals <- read.delim(paste0(out,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(out,"means.txt"), check.names = FALSE)

mymeans %>% dplyr::select("interacting_pair",starts_with("CD8 Taa"),ends_with("CD8 Taa"))  %>%  
  reshape2::melt() -> meansdf

colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals %>% dplyr::select("interacting_pair",starts_with("CD8 Taa"),ends_with("CD8 Taa"))%>%  
  reshape2::melt()-> pvalsdf

colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

summary((filter(pldf,means >1))$means)

pldf%>% filter(means >1) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 2)+ theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 1)) 
