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

# figure4B celltype_umap(Tcm单独聚类）(4A,3A在S4）-----------------------------------------------------------------
###Tcm单独聚类
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
DefaultAssay(Cell.combined_CD8) <- "integrated"
unique(Cell.combined_CD8@meta.data$celltype)
Cell.combined_Tcm<-subset(Cell.combined_CD8, subset = celltype  %in% "Tcm")
set.seed(1)
Cell.combined_Tcm <- RunPCA(Cell.combined_Tcm, npcs = 30, verbose = T)
Cell.combined_Tcm <- RunUMAP(Cell.combined_Tcm, seed.use = 42, reduction = "pca", dims = 1:20)
Cell.combined_Tcm <- FindNeighbors(Cell.combined_Tcm, reduction = "pca", dims = 1:20)
Cell.combined_Tcm <- FindClusters(Cell.combined_Tcm, resolution = c(0.1,0.2))
p1 <- DimPlot(Cell.combined_Tcm, reduction = "umap", group.by = "integrated_snn_res.0.1",label = TRUE)
p1
p2<-DotPlot(object = Cell.combined_Tcm, features = c("XCL1","IFNG"), group.by  = 'integrated_snn_res.0.1',
            cols = c("#21aba5", "#e84a5f"))
p2
Cell.combined_Tcm$celltype1 <- rep("XCL1high Tcm")
colnames(Cell.combined_Tcm@meta.data)
Cell.combined_Tcm@meta.data[which(Cell.combined_Tcm@meta.data$integrated_snn_res.0.1 == "1"),24] <- "IFNGhigh Tcm"
p3 <- DimPlot(Cell.combined_Tcm, reduction = "umap", group.by = "celltype1",label = TRUE)
p3
save(Cell.combined_Tcm, file="/data1/zhuzn/wangsn/reviewer_4000/figure4/Cell.combined_Tcm.RData")

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/Cell.combined_Tcm.RData")
colnames(Cell.combined_Tcm@meta.data)
cell_Tcm <- Cell.combined_Tcm@meta.data[,c(4,5,24)]
colnames(cell_Tcm) <- c("Tissue","Age","celltype")
cell_Tcm %<>% {.$id<-rownames(.);.} 

load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
colnames(Cell.combined_CD8@meta.data)
all_cell <- Cell.combined_CD8@meta.data[,c(4,5,18)]
all_cell %<>% {.$id<-rownames(.);.} %>% join(.,cell_Tcm[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell) <- c("Tissue","Age","celltype","id","Celltype")
all_cell[is.na(all_cell)] <- "NON"
all_cell$celltype <- as.character(all_cell$celltype)
all_cell$Celltype <- ifelse(all_cell$Celltype == "NON" ,all_cell$celltype,all_cell$Celltype)
unique(all_cell$Celltype)
Cell.combined_CD8 <- AddMetaData(Cell.combined_CD8,all_cell$Celltype,col.name = "celltype1")
p4 <- DimPlot(Cell.combined_CD8, reduction = "umap", pt.size = 0.25, group.by = "celltype1")
p4
save(Cell.combined_CD8, file="/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/Cell.combined_Tcm.RData")
colors <- c('#39A8E3',"#D75156")
Cell.combined_Tcm$Age <- factor(Cell.combined_Tcm$Age,levels=rev(c("young","old")))
table(Cell.combined_Tcm$Age)
p4 <- DimPlot(Cell.combined_Tcm, reduction = "umap", pt.size = 1, group.by = "Age",cols = rev(colors))+
  ggtitle("")+
  theme(legend.position = 'none',axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p4
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4B1.tiff",plot = p4,dpi=1000,width = 5,height = 4.5)

colors <- c("#6F6AF5","#D74048")
colors <- c("#6F6AF560","#D7404860")
colors <- c("#A3A0F9","#E48287")
p1 <- DimPlot(Cell.combined_Tcm, reduction = "umap", pt.size = 1, group.by = "celltype1",cols = colors, label = F)+
  ggtitle("")+
  theme(legend.position = 'none',axis.text = element_blank(),axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4B2.tiff",plot = p1,dpi=1000,width = 5,height = 4.5)

mat1=as.data.frame(Cell.combined_Tcm$celltype1)
colnames(mat1)="celltype1"
mat2=Embeddings(Cell.combined_Tcm,"umap")
mat3=merge(mat2,mat1,by="row.names")
colors=c("#6F6AF540","#D7404840")
p1 <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=celltype1),size = 1)+
  scale_color_manual(values=colors)+theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        panel.grid = element_blank())
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4B1.pdf",plot = p1,width = 6,height = 6)#图例文件



# figure4C SASP-----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/Cell.combined_Tcm.RData")
data.list_age <- SplitObject(Cell.combined_Tcm, split.by = "Age")

DefaultAssay(Cell.combined_Tcm) <- "RNA"
DefaultAssay(Cell.combined_Tcm)
sample.integrated <- Cell.combined_Tcm
SASP_gene_set <- read.delim('./SASP_V2.list',col.names = F)
gene <- as.list(SASP_gene_set)
sample.integrated <- AddModuleScore(
  object = sample.integrated,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'SASP_Features',
  replace = T
)#计算结果保存在pbmc@meta.data[["CD_Features1"]]  得到的score是在每个细胞中算出来的gene的表达均值。

names(sample.integrated@meta.data)
sasp <- sample.integrated@meta.data[,c(24,25)]

sasp2 <-aggregate(sasp$SASP_Features1,list(sasp$celltype1),mean)
names(sasp) <- c('Group','Gene_set_score')
names(sasp2) <- c('Group','mean')

library('ggpubr')
library(ggridges)
library(ggplot2)
sasp3 <- compare_means(Gene_set_score~Group, data=sasp,method = "t.test")
sasp$Group <- factor(sasp$Group,levels = c("XCL1high Tcm","IFNGhigh Tcm"))
sasp2$Group <- factor(sasp2$Group,levels = c("XCL1high Tcm","IFNGhigh Tcm"))
p<-ggplot(sasp, aes(x = Gene_set_score))+
  geom_density(aes(fill = Group), alpha=0.4)+
  # geom_signif(comparisons = list(c("young", "old")))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_vline(data = sasp2, aes(xintercept = mean, color=Group),linetype='dashed')+
  labs(x = "Gene set score",y = "Density",title = paste0("SASP gene set"),
       subtitle = paste0(sasp3$method,":P=",sasp3$p))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4C.pdf",plot = p,dpi=1000,width = 5,height = 4)

# figure4D Tcm差异分析 -----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/Cell.combined_Tcm.RData")
DefaultAssay(Cell.combined_Tcm) <- "RNA"
Cell.combined_Tcm <- SetIdent(object = Cell.combined_Tcm, value = Cell.combined_Tcm@meta.data$celltype1)
data.list_age <- SplitObject(Cell.combined_Tcm, split.by = "Age")
age_data <- data.list_age[["old"]]
diff_gene <- FindMarkers(age_data, ident.1 = "XCL1high Tcm", ident.2 = "IFNGhigh Tcm", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Tcm","_old_1_2.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)

# figure4D Tcm差异基因Top10基因表达气泡图 -----------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/Cell.combined_Tcm.RData")
colnames(Cell.combined_Tcm@meta.data)
Cell.combined_Tcm <- SetIdent(object = Cell.combined_Tcm, value = Cell.combined_Tcm@meta.data$celltype1)
cell_markers <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Tcm","_old_1_2.findmark.txt"))
cell_markers <- cell_markers[!grepl("ENSMMUG", cell_markers$gene),]
cell_markers <- cell_markers[cell_markers$change_V1 != "NOT",]
gene <- as.character(unique(cell_markers$gene[c(1:20,(nrow(cell_markers)-19):nrow(cell_markers))]))
gene <- c("XCL1","GZMB","CCL4L1","CCL5","ID2","KLRD1","CCL3","NFKBIA","GZMA","RGS1",
          "IFNG","LTB","CTSW","GZMK","SAMHD1","KLRB1","GZMM","THY1","IL7R","IGFBP2")
d <- cell_markers[cell_markers$gene %in% gene,]
p1<-DotPlot(object = Cell.combined_Tcm, features = gene, group.by  = 'celltype1',cols = c("grey", "#CE1212")) + 
  theme(legend.position = "top",axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
p1
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4D.pdf",plot = p1,dpi=1000,width = 6,height = 3)

# figure4E 富集(Tcm单独聚类）------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
DE <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure4/Tcm_old_1_2.findmark.txt")
# DE <- DE[DE$change_V1 == "UP",]
DE <- DE[DE$change_V1 == "DOWN",]
colnames(DE)[6] <- "SYMBOL"
df_id <- bitr(DE$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
DE <- merge(DE,df_id,by = "SYMBOL",all=F)
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db);ego_BP_result<-as.data.frame(ego_BP1@result)
# write.table(ego_BP_result,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Tcm","_GO_UP.txt"),quote = FALSE,sep = "\t",row.names = F)
write.table(ego_BP_result,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Tcm","_GO_DOWN.txt"),quote = FALSE,sep = "\t",row.names = F)

# GO富集结果合并图 
rm(list=ls())
library("ggplot2")
library( "clusterProfiler")
DE_list <- c("UP","DOWN")
GO <- data.frame()
for (DE in DE_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/Tcm_GO_",DE,".txt"),stringsAsFactors = F)
  GO1$DE <- rep(DE)
  GO <- rbind(GO,GO1)
}
GO <- GO[GO$pvalue < 0.05,]
write.table(GO, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Tcm_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

# ID <- c("GO:2000379","GO:0032635","GO:0038066","GO:0097190","GO:0043068","GO:0032101",
#         "GO:0006412","GO:0043043","GO:0042254","GO:0006364","GO:0022613")
ID <- c("positive regulation of reactive oxygen species metabolic process",
        "p38MAPK cascade","interleukin-6 production","apoptotic signaling pathway",
        "positive regulation of programmed cell death","regulation of response to external stimulus",
        "ribosome biogenesis","rRNA processing","ribonucleoprotein complex biogenesis","translation","peptide biosynthetic process")
GO <- GO[GO$Description %in% ID,]

GO$Description <- factor(GO$Description,levels = rev(ID))
GO$DE <- factor(GO$DE,levels = unique(GO$DE))
p <- ggplot(GO,aes(Description,DE,color = -log10(pvalue), size=Count)) + geom_point()+
  coord_flip() +
  # scale_colour_gradient(low="grey",high="#CB1217")+
  scale_colour_gradientn(colors=rev(c("#810000","#CE1212","#F05454")))+
  geom_point()+scale_size(range=c(5,10))+
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels=function(x) str_wrap(x, width=60))+
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 20, color = "black"),plot.title = element_text(size = 30),
        axis.text.y = element_text(size = 25, color = "black"),
        legend.position = "right")
p
unique(GO$Description)
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4E.pdf",plot = p,dpi=1000,width = 12.5,height = 10)
dev.off()

# figure4F 伪时间------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
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
cell_data <- Cell.combined_CD8

unique(cell_data@meta.data$celltype1)
DefaultAssay(cell_data) <- "RNA"
##创建CDS对象并预处理数据
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
get_earliest_principal_node <- function(cds, NSC_type="Naive"){
  cell_ids <- which(colData(cds)[, "celltype1"] == NSC_type)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
save(cds,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","RNA_monocle",".RData"))

# 伪时间映射细胞类型 ---------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/RNA_monocle.RData")
pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
Cell.combined_CD8 <- AddMetaData(Cell.combined_CD8,pseudotime,col.name = "Pseudotime")
colnames(Cell.combined_CD8@meta.data)
mat1=as.data.frame(Cell.combined_CD8@meta.data[,c(5,23)])
mat2=Embeddings(Cell.combined_CD8,"umap")
mat3=merge(mat2,mat1,by="row.names")

p <- mat3%>%ggplot(aes(UMAP_1,UMAP_2))+geom_point(aes(color=Pseudotime),size = 0.01)+
  # scale_color_gradientn(colors = c("#140789","#38049A","#7701A8","#9B169F","#AB2494","#BC3587","#CD4A76",
  #                                  "#EB7655","#F79143","#F9983E","#FDAE32","#FBD324","#F0F921"))+
  scale_color_gradientn(colors = c("#140789","#38049A","#9B169F","#AB2494",
                                   "#F9983E","#FDAE32","#FBD324","#F0F921"))+
  theme_bw()+
  theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
        panel.grid = element_blank(),panel.border = element_blank(),legend.position = "none")
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4A3.pdf",plot = p,dpi=1000,width = 4,height = 4)

# figure4F 伪时间基因表达-------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
DefaultAssay(Cell.combined_CD8) <- "RNA"
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype1)
diff_gene <- FindMarkers(Cell.combined_CD8, ident.1 = "XCL1high Tcm", verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
diff_gene <- diff_gene[!grepl("ENSMMUG", diff_gene$gene),]
diff_gene <- diff_gene[diff_gene$change_V1 != "NOT",]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","XCL1high Tcm_vs_others.txt"),quote = FALSE,sep = "\t",row.names = F)
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figureS4/Cell.combined_CD8.RData")
DefaultAssay(Cell.combined_CD8) <- "RNA"
Cell.combined_CD8 <- SetIdent(object = Cell.combined_CD8, value = Cell.combined_CD8@meta.data$celltype1)
diff_gene <- FindMarkers(Cell.combined_CD8, ident.1 = "IFNGhigh Tcm", verbose = T,only.pos = F)
diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NOT')),
                                       'NOT'))
diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
diff_gene <- diff_gene[!grepl("ENSMMUG", diff_gene$gene),]
diff_gene <- diff_gene[diff_gene$change_V1 != "NOT",]
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","IFNGhigh Tcm_vs_others.txt"),quote = FALSE,sep = "\t",row.names = F)

rm(list=ls())
DE1 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure4/XCL1high Tcm_vs_others.txt")
DE2 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure4/IFNGhigh Tcm_vs_others.txt")
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure4/RNA_monocle.RData")
Track_genes_sig <- c("NFKBIA","CCL3")
Track_genes_sig <- c("CCL3","FOSB","CCL4L1","CD74","DNAJA1","DNAJB1","RGS1","IFNG","XCL1")
Track_genes_sig <- c("CD74","NFKBIA")
cds@colData$celltype1 <- factor(cds@colData$celltype1,levels = c("Naive","Tex","Tem","IFNGhigh Tcm","XCL1high Tcm","TRM","CTL","CTL/TRM"))
# # colors <- c("#1f78b4","#33a02c","#ff7f00","#e31a1c","#fb9a99","#e78ac3","#83639D")
colors <- c("#B8A99D50","#B8A99D50","#B8A99D50","#6F6AF5","#D74048","#B8A99D50","#B8A99D50","#B8A99D50")
show_col(colors)
#基因表达趋势图
p <-plot_genes_in_pseudotime(cds[Track_genes_sig], color_cells_by="celltype1",min_expr=0.5, ncol = 2)+
  scale_color_manual(values=colors)
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure4/","Tcm_2_gene",".pdf"),p,width = 9,height = 3)

# 实验统计图 -------------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
data <- data.frame(IFNG = c(31.3,32.5,31.9,15.2,14.5,15.8),
                   XCL1 = rep(100,6),
                   Aged = c(rep("Young",3),rep("Aged",3)))
data$XCL1 <- data$XCL1 - data$IFNG
data <- melt(data,id.vars = "Aged")
data$Group <- paste0(data$Aged,"_",data$variable)
unique(data$Group)
data$Group <- factor(data$Group,levels = c("Young_IFNG","Young_XCL1","Aged_IFNG","Aged_XCL1"))
my_comparisons <- list(c("Young_IFNG","Young_XCL1"),c("Aged_IFNG","Aged_XCL1"))
p <- ggplot(data=data, aes(x=Group,y=value,color=variable,fill=variable))+
  stat_summary(fun = mean,geom = "bar",aes(color=variable,fill=variable,width=0.5))+
  stat_summary(geom = "errorbar",fun.min=min,fun.max=max,width=0.2)+
  # stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))+
  geom_signif(comparisons = my_comparisons,map_signif_level=T,test = "t.test")+
  # scale_y_continuous(limits =c(0, 15) ,expand = c(0,0))+ 
  theme_bw()+
  scale_color_manual(values=c("blue","red"))+
  scale_fill_manual(values=c("white","white"))+
  theme(legend.position = "none",legend.title=element_blank())+labs(x="",y="Relative RNA expression",caption="")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"))
p

# 实验统计图 -------------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
data <- data.frame(value = c(31.3,32.5,31.9,15.2,14.5,15.8),
                   Aged = c(rep("IFNG+Tcm",3),rep("XCL1+Tcm",3)))
data$Aged <- factor(data$Aged,levels = c("IFNG+Tcm","XCL1+Tcm"))
p <- ggplot(data=data, aes(x=Aged,y=value,color=Aged,fill=Aged))+
  # geom_boxplot()+
  stat_summary(fun = mean,geom = "bar",aes(color=Aged,fill=Aged,width=0.5))+
  stat_summary(geom = "errorbar",fun.min=min,fun.max=max,width=0.2)+
  theme_bw()+
  scale_color_manual(values=c("blue","red"))+
  scale_fill_manual(values=c("white","white"))+
  theme(legend.title=element_blank())+labs(x="",y="GZMK/Tcm(%)",caption="")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure4/figure.4G.pdf",plot = p,width = 3.5,height = 3.5)
