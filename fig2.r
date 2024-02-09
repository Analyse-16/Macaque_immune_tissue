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

# figure2  DEGs_tissue_celltype ----------------------------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
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
                                               ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                      ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                               'NOT'))
        celltype <- gsub(" ","_",celltype)
        write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),quote = FALSE,sep = "\t",row.names = F)
      } 
    }
  }
}

# figure2A 50% 差异基因 ----------------------------------------------------------------
rm(list=ls())
tissue_list <- c("Mesenteric_lymph","PBMC","Spleen","Bone_marrow")
celltype_list <- c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                   "Fibroblast","Erythrocyte","Megakaryocyte","Granulocyte","HSC")
DEG_data <- data.frame(gene = "gene",change_V1 = "change_V1",celltype = "celltype")

for (tissue in tissue_list){
  for (celltype in celltype_list){
    celltype <- gsub(" ","_",celltype)
    if (file.exists(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"))){
      dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
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
change_data$change_percent <- (change_data$sig_data/max(change_data$sig_data))*100
gene_list <- unique(change_data$gene)

# tissue = "Mesenteric_lymph"
# celltype = "T_Cell"
# dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
# gene_FC <- data.frame(p_val = "p_val",avg_log2FC = "avg_log2FC",gene = "gene",change_V1 = "change_V1",celltype = "celltype")
# 
# for (tissue in tissue_list){
#   for (celltype in celltype_list){
#     celltype <- gsub(" ","_",celltype)
#     if (file.exists(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"))){
#       dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",tissue,"_",celltype,"_.findmark.txt"),stringsAsFactors = F)
#       for (gene in gene_list){
#         if(gene %in%  dif_data$gene) {
#           gene_dif1 <- dif_data[dif_data$gene == gene,c(1,2,6,7)]
#           gene_dif1$celltype <- paste0(tissue,"_",celltype)
#           gene_FC <- rbind(gene_FC,gene_dif1)
#         }else{
#           gene_dif1 <- data.frame(p_val = "0",avg_log2FC = "0",gene = gene,change_V1 = "NOT")
#           gene_dif1$celltype <- paste0(tissue,"_",celltype)
#           gene_FC <- rbind(gene_FC,gene_dif1)
#         }
#       }
#     }
#     else{
#       next
#     }
#   }
# }
# FC <- gene_FC
# FC <- FC[-1,]
# write.table(FC, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/","sig","FC.tab"), quote = FALSE,sep="\t",row.names = FALSE)

FC <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/sigFC.tab")
FC$avg_log2FC <- as.numeric(FC$avg_log2FC)
FC$p_val <- as.numeric(FC$p_val)
FC_mean <- aggregate(FC[,1:2],by=list(gene=FC$gene),FUN=mean)
range(FC_mean$avg_log2FC)
range(FC_mean$p_val)
nrDEG_mean <- merge(change_data, FC_mean, by = "gene", all = TRUE)
nrDEG_mean$sig_gene <- nrDEG_mean$change_V1
nrDEG_mean$sig_gene <- ifelse(nrDEG_mean$change_percent > 50,ifelse(nrDEG_mean$avg_log2FC > 0, 'UP', 'DOWN' ),'FALSE')

D <- data.frame(table(nrDEG_mean$sig_data))
D$Var1 <- as.numeric(D$Var1)
max(D$Var1)
D$Var1 <- factor(D$Var1,levels = c(1:max(D$Var1)))

library(gg.gap)
p <- ggplot(D, aes(Var1, Freq,group = 1))+
  geom_line(size=1,color='#3398CC')+geom_point(size=10,color='#D1484F')+
  labs(x = 'Unique Cell Identities', y = 'DiffEx Genes')+
  # scale_y_continuous(breaks = seq(0,3000,100))+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(size = 1,colour = "black"),
        axis.text=element_text(size=20,color="black"),axis.title=element_text(size=20,color="black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure2/figure.2A.pdf",plot = p,dpi=1000,width = 15,height = 10)

gene_sig <- nrDEG_mean[nrDEG_mean$change_percent > 50 & nrDEG_mean$p_val < 0.1,]
gene_sig <- gene_sig[!grepl("ENSMMUG", gene_sig$gene),]
gene_sig <- gene_sig[,-2]
gene_sig <- unique(gene_sig)
write.table(gene_sig, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/","sig","_50.tab"), quote = FALSE,sep="\t",row.names = FALSE)

# figure2 B heatmap(交集TF) -----------------------------------------------------------------
rm(list=ls())
sig_50 <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/sig_50.tab")
TF <- read.table("/data1/zhuzn/wangsn/V6/pyscenic/data/hs_hgnc_tfs.txt", quote="\"", comment.char="")
sig_TF <- intersect(sig_50$gene,TF$V1)
colnames(sig_50)

sig_gene <- sig_50
sig_gene <- sig_gene[order(sig_gene[,6],decreasing=T),]
sig_50$gene <- as.factor(sig_50$gene)
gene_list <- unique(sig_gene$gene)
# f.list <- list.files("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/")
# f.list <- f.list[grep("findmark.txt$",f.list)]
tissue_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
celltype_list <- c("T_Cell","NK","B_Cell","Monocyte","Macrophage","Dendritic_cell","Fibroblast","Megakaryocyte","Granulocyte","HSC")
f.list <- c()
for (tissue in tissue_list){
  for (celltype in celltype_list){
    f.list1 <- paste0(tissue,"_",celltype,"_.findmark.txt")
    f.list <- c(f.list,f.list1)
  }
}
# f.list <- as.factor(f.list)
gene_dif <- data.frame(p_val = "p_val",avg_log2FC = "avg_log2FC",gene = "gene",change_V1 = "change_V1",
                       celltype = "celltype")
for (i in f.list){
  dif_data <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/",i),stringsAsFactors = F)
  for (gene in gene_list){
    if(gene %in%  dif_data$gene) {
      gene_dif1 <- dif_data[dif_data$gene == gene,c(1,2,6,7)]
      gene_dif1$celltype <- i
      gene_dif <- rbind(gene_dif,gene_dif1)
    }else{
      gene_dif1 <- data.frame(p_val = "0",avg_log2FC = "0",gene = gene,change_V1 = "NOT",celltype = i)
      gene_dif <- rbind(gene_dif,gene_dif1)
    }
  }
}
gene_dif <- gene_dif[-1,]
gene_dif$celltype <-gsub("_.findmark.txt","",gene_dif$celltype)
gene_dif$celltype <- factor(gene_dif$celltype,levels = unique(gene_dif$celltype))
gene_dif$gene <- factor(gene_dif$gene,levels = rev(gene_list))
gene_dif$avg_log2FC <- as.numeric(gene_dif$avg_log2FC)
library(ggplot2)
p <- ggplot(gene_dif,aes(gene,celltype)) + geom_tile(aes(fill=avg_log2FC),color = "white") +
  guides(fill=guide_colorbar("avg_Expression_FC")) +
  coord_flip() +
  scale_fill_gradient2(limits=c(-4,6),breaks = seq(-4, 6, by = 1),low="#52B088",mid="#eff3ff",high ="#8152E7",midpoint = 0)+
  # scale_fill_gradient2(limits=c(-4,6),breaks = seq(-4, 6, by = 1),low="#154889",mid="#eff3ff",high ="#AB2524",midpoint = 0)+
  labs(x = "",y = "",title = "") +
  theme(axis.text.x = element_text(size = 15, hjust = 1, vjust = 0.5, angle = 90,colour = "black"),
        axis.text.y = element_text(size = 15,colour = "black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure2/figure2B.pdf",plot = p,dpi=1000,width = 9,height = 10)
dev.off()

# figure2 C 富集 -----------------------------------------------------------------
rm(list=ls())
sig_gene <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_tissue_celltype/sig_50.tab")
SASP_gene_set <- read.delim('/data1/zhuzn/wangsn/SASP/SASP_V2.list',col.names = F)
in_gene <- intersect(sig_gene$gene,SASP_gene_set$FALSE.)

#GO
library("ggplot2")
library( "clusterProfiler")
library("org.Mmu.eg.db")
columns(org.Mmu.eg.db)
geneinfo = select(org.Mmu.eg.db, keys=keys(org.Mmu.eg.db), columns = c('ENTREZID',"SYMBOL"))

sig_gene <- sig_gene[sig_gene$sig_gene == "UP",]
colnames(sig_gene)[1] = "SYMBOL"
DE <- merge(x=sig_gene,y=geneinfo,by="SYMBOL")
gene <- unique(DE[,'ENTREZID'])
ego_BP <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",
                   pAdjustMethod = "BH",minGSSize = 1,
                   pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Mmu.eg.db)
ego_BP_result_UP<-as.data.frame(ego_BP1@result)#8107
ego_BP_result_UP$group <- rep("UP")
c <- c("cellular response to interferon-gamma","glycolytic process","granulocyte chemotaxis",
       "negative regulation of translation","regulation of viral life cycle")
top_go_up <- ego_BP_result_UP[ego_BP_result_UP$Description %in% c,]
top_go_up <- top_go_up[order(top_go_up$Count,decreasing = T),]
top_go_up$Number <- c(5:1)

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
c <- c("positive regulation of T cell chemotaxis","positive regulation of lymphocyte migration",
       "cellular response to chemokine","translation")
top_go_down <- ego_BP_result_DOWN[ego_BP_result_DOWN$Description %in% c,]
top_go_down <- top_go_down[order(top_go_down$Count,decreasing = F),]
top_go_down$Number <- c(-1:-4)
top_go <- rbind(top_go_up,top_go_down)
top_go$number = factor(rev(1:nrow(top_go)))
top_go$Count = as.numeric(ifelse(top_go$group == 'DOWN' ,paste0("-",top_go$Count),top_go$Count))

GO_BP <- rbind(ego_BP_result_UP,ego_BP_result_DOWN)
GO_BP <- GO_BP[GO_BP$pvalue < 0.05,]
write.table(GO_BP, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","figure.2C.txt"), quote = FALSE,sep="\t",row.names = FALSE)

p <- ggplot(data=top_go, aes(x=number, y=Count,fill = Number)) + geom_bar(stat="identity", width=0.8) + 
  scale_x_discrete(labels=rev(top_go$Description)) +
  scale_fill_gradient2(low="#52B088",mid="#eff3ff",high ="#8152E7",midpoint = 0)+
  coord_flip() +theme_bw() + xlab("") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))+
  theme(axis.text=element_text(size=20,face = "plain", color="black"),panel.grid = element_blank(),
        axis.title=element_text(size=10),legend.text=element_text(size=10),legend.title = element_text(size=20))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure2/figure.2C.pdf",plot = p,dpi=1000,width = 10,height = 5)

# figure2 D GZMB 雷达图----------------------------------------------------------------
rm(list = ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
all_sample.combined1 <- subset(all_sample.combined, subset = celltype  %in% c("T Cell","NK","B Cell","Monocyte","Macrophage","Dendritic cell",
                                                                              "Fibroblast","Megakaryocyte","Granulocyte","HSC"))
all_sample.combined1@meta.data$id <- colnames(all_sample.combined1)
all_sample.combined1$celltype <- as.character(all_sample.combined1$celltype)
load(file = "/data1/zhuzn/wangsn/V8/all_sample.combined.celltype.RData")
all_sample.combined2 <- subset(all_sample.combined, subset = celltype  %in% c("Fibroblast","Megakaryocyte"))
all_sample.combined2@meta.data$id <- colnames(all_sample.combined2)
all_sample.combined2$id <-gsub("-",".",all_sample.combined2$id)
all_sample.combined2$id <-gsub("Lymph_node","Mesenteric_lymph",all_sample.combined2$id)
table(all_sample.combined2$celltype)
all_sample.combined1$Celltype <- ifelse(all_sample.combined1$id %in% all_sample.combined2$id,all_sample.combined2$celltype,all_sample.combined1$celltype)
colors <- c("#7756A1","#BEBADA","#1084D3","#73C508","#CFE6A7","#144F2B","#FF348B","#B7D5F3","#A91C51","#E9F033","#FFAF18")
p1 <- DimPlot(all_sample.combined1, reduction = "umap", group.by = "celltype",cols = colors,pt.size = 0.01,label = F, repel = TRUE,raster=F)
p2 <- DimPlot(all_sample.combined1, reduction = "umap", group.by = "Celltype",cols = colors,pt.size = 0.01,label = F, repel = TRUE,raster=F)
p <- p1/p2
p
table(all_sample.combined1$Celltype)
all_sample.combined1$celltype <- all_sample.combined1$Celltype
data.list_age <- SplitObject(all_sample.combined1, split.by = "Age")
age_list <- unique(all_sample.combined1$Age)
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
table(old$Tissue)
data <- cbind(young,old)
colnames(data) <- c("avg.exp_young","Cluster","Age","Tissue","avg.exp_old","Cluster","Age","Tissue")
data <- data[,c(4,2,1,5)]
# data$fc <- data$avg.exp_old - data$avg.exp_young
data$fc <- data$avg.exp_old/data$avg.exp_young

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
dt$Granulocyte[!is.finite(dt$Granulocyte)] <- 0;dt$HSC[!is.finite(dt$HSC)] <- 0;dt$Megakaryocyte[!is.finite(dt$Megakaryocyte)] <- 0
p1 <- ggradar(dt,axis.label.size=6.5,
              # values.radar = c(0, 1.5, 3),grid.min = 0,grid.mid = 1.5,grid.max = 3,
              values.radar = c(0, 14, 28),grid.min = 0,grid.mid = 14,grid.max = 28,
              group.colours = c('#23B363','#DF63CC','#339FBA','#CA7D44')) +
  ggtitle("GZMB mean expression (Aged vs Young)")+
  theme(legend.position = "none")
p1
ggsave(plot = p1, paste0('/data1/zhuzn/wangsn/reviewer_4000/figure2/','figure2D.pdf'), width = 14, height = 16,dpi = 600)

# q-PCR柱状图-GZMB----------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
GZMB <- read.delim("/data1/zhuzn/wangsn/V8/q-PCR_data/GZMB.txt")
GZMB$Group <- paste0(GZMB$group,"_",GZMB$gene)
unique(GZMB$Group)
GZMB$Group <- factor(GZMB$Group,levels = unique(GZMB$Group))
my_comparisons <- list(c("Young_GZMB","Aged_GZMB"),c("Young_S100A8","Aged_S100A8"),c("Young_S100A4","Aged_S100A4"))
p <- ggplot(data=GZMB, aes(x=Group,y=exp,color=group,fill=group))+
  stat_summary(fun = mean,geom = "bar",aes(color=group,fill=group,width=0.5))+
  stat_summary(geom = "errorbar",fun.min=min,fun.max=max,width=0.2)+
  # stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))+
  geom_signif(comparisons = my_comparisons,map_signif_level=T,test = "t.test")+
  # scale_y_continuous(limits =c(0, 15) ,expand = c(0,0))+ theme_bw()+
  scale_color_manual(values=c("red","blue"))+
  scale_fill_manual(values=c("white","white"))+
  theme(legend.title=element_blank())+labs(x="",y="Relative RNA expression",caption="")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure2/figure.2E.pdf",plot = p,width = 5,height = 3)

# q-PCR柱状图-GZMB----------------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(ggsignif)
library(ggpubr)
GZMB <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/GZMB.txt")
GZMB$Group <- paste0(GZMB$group,"_",GZMB$gene)
unique(GZMB$Group)
GZMB$Group <- factor(GZMB$Group,levels = c("Control_TIMP1","GZMB_TIMP1","Control_NFKBIA","GZMB_NFKBIA","Control_P21","GZMB_P21",
                                           "Control_JUN","GZMB_JUN","Control_GZMB","GZMB_GZMB","Control_S100A4","GZMB_S100A4",
                                           "Control_S100A8","GZMB_S100A8"))
p <- ggplot(data=GZMB, aes(x=Group,y=exp,color=group,fill=group))+
  stat_summary(fun = mean,geom = "bar",aes(color=group,fill=group,width=0.5))+
  stat_summary(geom = "errorbar",fun.min=min,fun.max=max,width=0.2)+
  # stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",width = 0.25,position = position_dodge( .9))+
  # geom_signif(comparisons = my_comparisons,map_signif_level=T,test = "t.test")+
  # scale_y_continuous(limits =c(0, 15) ,expand = c(0,0))+ 
  theme_bw()+
  scale_color_manual(values=c("blue","red"))+
  scale_fill_manual(values=c("white","white"))+
  theme(legend.position = "none",legend.title=element_blank())+labs(x="",y="Relative RNA expression",caption="")+
  theme(panel.grid = element_blank(),panel.border = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.line = element_line(size=0.5, colour = "black"))
p
ggsave("/data1/zhuzn/wangsn/reviewer_4000/figure2/figure.2E.pdf",plot = p,width = 5,height = 4)
