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

# 免疫细胞组织之间的差异基因、每个组织年老比年轻的差异基因 ------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")
DefaultAssay(all_sample.combined) <- "RNA"
colnames(all_sample.combined@meta.data)
all_cell <- all_sample.combined@meta.data[,c(4,5,18)]
### 添加亚型信息
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Monocyte_celltype.RData")
colnames(Cell.combined@meta.data)
cell_Mon <- Cell.combined@meta.data[,c(4,5,18)]
load(file = "/data1/zhuzn/wangsn/reviewer_4000/figure6/Macrophage_celltype.RData")
colnames(Cell.combined@meta.data)
cell_Mac <- Cell.combined@meta.data[,c(4,5,18)]
load(file = "/data1/zhuzn/wangsn/data_QC_test/f4000/Cell.combined_T.RData")
colnames(Cell.combined@meta.data)
cell_T <- Cell.combined@meta.data[,c(4,5,18)]

cell <- rbind(cell_Mon,cell_Mac,cell_T)
cell  %<>% {.$id<-rownames(.);.} 

all_cell %<>% {.$id<-rownames(.);.} %>% join(.,cell[,c("id","celltype")],by=c("id"="id"))
colnames(all_cell) <- c("Tissue","Age","celltype","id","Celltype")
all_cell[is.na(all_cell)] <- "NON"
all_cell$celltype <- as.character(all_cell$celltype)
all_cell$Celltype <- ifelse(all_cell$Celltype == "NON" ,all_cell$celltype,all_cell$Celltype)
unique(all_cell$Celltype)
all_sample.combined <- AddMetaData(all_sample.combined,all_cell$Celltype,col.name = "sub_celltype")
table(all_sample.combined@meta.data$celltype)
table(all_sample.combined@meta.data$sub_celltype)
save(all_sample.combined, file="/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000_sub.RData")

rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000_sub.RData")
data.list_celltype <- SplitObject(all_sample.combined, split.by = "sub_celltype")
cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
for (cell in cell_list){
  cell_data <- data.list_celltype[[cell]]
  DefaultAssay(cell_data) <- "RNA"
  cell_data <- SetIdent(object = cell_data, value = cell_data@meta.data$Tissue)
  diff_gene <- FindAllMarkers(cell_data, logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
  diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                         ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                         'NOT'))
  write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_tissue.txt"),quote = FALSE,sep = "\t",row.names = F)
}

for (cell in cell_list){
  cell_data <- data.list_celltype[[cell]]
  cell_data <- SetIdent(object = cell_data, value = cell_data@meta.data$Age)
  data.list_tissue <- SplitObject(cell_data, split.by = "Tissue")
  sample_list <- unique(cell_data@meta.data$Tissue)
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
        diff_gene <- FindMarkers(tissue_data, ident.1 = "old", ident.2 = "young", logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
        diff_gene <- diff_gene  %>% {.$avg_log2FC<-as.numeric(.$avg_log2FC);.} %>% {.$gene<-rownames(.);.}
        diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                               ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                                      ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                               'NOT'))
        write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_",tissue,"_old_vs_young.txt"),quote = FALSE,sep = "\t",row.names = F)
      } 
    }
  }
}

# 免疫细胞组织细胞类型与其他免疫细胞DE ------------------------------------------------
rm(list=ls())
load(file = "/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000_sub.RData")
unique(all_sample.combined$sub_celltype)
cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
Cell.combined<-subset(all_sample.combined, subset = sub_celltype  %in% cell_list)
unique(Cell.combined$sub_celltype)
Cell.combined$sub_celltype1 <- paste0(Cell.combined$Tissue,"_",Cell.combined$sub_celltype)
unique(Cell.combined$sub_celltype1)
Cell.combined <- SetIdent(object = Cell.combined, value = Cell.combined@meta.data$sub_celltype1)
diff_gene <- FindAllMarkers(cell_data, logfc.threshold = 0,min.pct = 0,verbose = T,only.pos = F)
diff_gene$change_V1 = as.factor(ifelse(diff_gene$p_val < 0.05,
                                       ifelse(diff_gene$avg_log2FC > 0.25 ,'UP',
                                              ifelse(diff_gene$avg_log2FC < -0.25 ,'DOWN','NEUTRAL')),
                                       'NOT'))
write.table(diff_gene,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/tissue/","tissue1.txt"),quote = FALSE,sep = "\t",row.names = F)

# 免疫细胞点图-组织细胞类型差异 ------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
# cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
# cell_list <- c("CD16_Mon","In-termed_Mon")
cell_list <- c("M1","M2")
for (cell in cell_list){
  # tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"tissue.txt")) %>% .[.$change_V1 == "UP",]
  tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/tissue/tissue1.txt")) %>% .[.$change_V1 == "UP",]
  DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
  Bone_marrow <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","Bone_marrow","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("Bone_marrow","_",cell),]
  Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% intersect(DE1$gene,Bone_marrow$gene),]
  Mesenteric_lymph <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","Mesenteric_lymph","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("Mesenteric_lymph","_",cell),]
  Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% intersect(DE1$gene,Mesenteric_lymph$gene),]
  PBMC <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","PBMC","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("PBMC","_",cell),]
  PBMC <- PBMC[PBMC$gene %in% intersect(DE1$gene,PBMC$gene),]
  Spleen <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","Spleen","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("Spleen","_",cell),]
  Spleen <- Spleen[Spleen$gene %in% intersect(DE1$gene,Spleen$gene),]
  
  Bone_marrow1 <- Bone_marrow[Bone_marrow$gene %in% Reduce(setdiff,list(Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene,Spleen$gene)),]
  Bone_marrow1 <- Bone_marrow1[order(Bone_marrow1$avg_log2FC,decreasing = T),]
  
  Mesenteric_lymph1 <- Mesenteric_lymph[Mesenteric_lymph$gene %in% Reduce(setdiff,list(Mesenteric_lymph$gene,Bone_marrow$gene,PBMC$gene,Spleen$gene)),]
  Mesenteric_lymph1 <- Mesenteric_lymph1[order(Mesenteric_lymph1$avg_log2FC,decreasing = T),]
  
  PBMC1 <- PBMC[PBMC$gene %in% Reduce(setdiff,list(PBMC$gene,Bone_marrow$gene,Mesenteric_lymph$gene,Spleen$gene)),]
  PBMC1 <- PBMC1[order(PBMC1$avg_log2FC,decreasing = T),]
  
  Spleen1 <- Spleen[Spleen$gene %in% Reduce(setdiff,list(Spleen$gene,Bone_marrow$gene,Mesenteric_lymph$gene,PBMC$gene)),]
  Spleen1 <- Spleen1[order(Spleen1$avg_log2FC,decreasing = T),]

  diff_gene <- c(Bone_marrow1$gene,Mesenteric_lymph1$gene,PBMC1$gene,Spleen1$gene)
  diff_gene <- diff_gene[!grepl("ENSMMUG", diff_gene)]
  # diff_gene <- c("HBB","HBA","LTF","JCHAIN","PLAC8","C1QB","GZMB","LTA4H","CD163","LY6D","LILRA3","LRIF1")
  diff_gene_overlap_FC <- data.frame()
  for (i in DE_list){
    old_vs_young <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_",i,"_old_vs_young.txt"))
    diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
    diff_gene2 <- tissue_DE[tissue_DE$cluster == paste0(i,"_",cell),c(2,7)]
    diff_gene1 <- merge(diff_gene1,diff_gene2,all.x=TRUE,by = "gene")
    colnames(diff_gene1) <- c("gene","log2FC.old_young","log2FC.tissue")
    diff_gene1$celltype <- rep(i)
    diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
  }
  diff_gene_overlap_FC[is.na(diff_gene_overlap_FC)] <- "0"
  overlap_FC <- diff_gene_overlap_FC
  str(overlap_FC)
  overlap_FC$log2FC.tissue <- as.numeric(overlap_FC$log2FC.tissue)
  overlap_FC$celltype <- factor(overlap_FC$celltype,levels = DE_list)
  overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
  p<-ggplot(overlap_FC,aes(x = celltype,y=gene,size=log2FC.old_young,colour=log2FC.tissue,ylab=''))+
    geom_point()+
    # coord_flip() +
    scale_size_continuous(range=c(0,6))+
    scale_color_gradientn(colors = c("#00008A","#0000FF","#3636FF","#4D4DFF","#9B9BFF","#F7C0C0","#F78888","#F76C6C","#981216"))+
    # scale_color_gradient2(low = "blue",mid = "grey", high = "#981216",midpoint = 1)+
    theme_classic()+
    labs(x = '', y = '',title = cell)+
    theme(axis.text.x=element_text(size=15, color="black",angle = 45,vjust = 1,hjust = 1),legend.position = "right",
          axis.text.y=element_text(size=15, color="black"))# + RotatedAxis():45度
  p
  # ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/","immune_",cell,".png"),p,width = 10,height = 40)
  ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/","immune_",cell,".png"),p,width = 10,height = 55,limitsize = FALSE)
}

# 免疫细胞点图-富集 ------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
cell_list <- c("CD8")
for (cell in cell_list){
  # tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_tissue.txt")) %>% .[.$change_V1 == "UP",]
  tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/tissue/tissue1.txt")) %>% .[.$change_V1 == "UP",]
  DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
  Bone_marrow <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","Bone_marrow","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("Bone_marrow","_",cell),]
  Bone_marrow <- Bone_marrow[Bone_marrow$gene %in% intersect(DE1$gene,Bone_marrow$gene),]
  Mesenteric_lymph <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","Mesenteric_lymph","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("Mesenteric_lymph","_",cell),]
  Mesenteric_lymph <- Mesenteric_lymph[Mesenteric_lymph$gene %in% intersect(DE1$gene,Mesenteric_lymph$gene),]
  PBMC <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","PBMC","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("PBMC","_",cell),]
  PBMC <- PBMC[PBMC$gene %in% intersect(DE1$gene,PBMC$gene),]
  Spleen <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_","Spleen","_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  DE1 <- tissue_DE[tissue_DE$cluster == paste0("Spleen","_",cell),]
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
  write.table(ego_BP_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/",cell,"_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
}

# 免疫细胞点图-合并 ----------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
diff_gene <- c("S100A9","LTF","SOD2","LGALS3","IFI16","VSIR",
               "CCL3","IFNG","NFKBIA","TNFAIP3","DUSP1","ISG15","CCL4L1","DUSP10","CD83","ZFP36","ZC3H12A","IRF4","CD74",
               "AFAP1L2","CD226","PRF1","NCR3","SLC11A1","PPBP","SYT11","VCL","ADGRG1","CIB1","FLNA","CDC42","NRROS",
               "KLRD1","LAG3","CCL5","UBASH3B","FGR","ISG20","LYN","SWAP70","BST2",
               "CHMP7","IGKC","APOL2",
               "CX3CR1","RGS9","FCRL6","CSF1",
               "FOS","FOSB","GZMB","COTL1",
               "RTD1A","RASSF4","C1QC","JCHAIN","VMO1","LTA4H","LY6D","LILRA3",
               "DEFA1B","AHSP","CXCL9","UBD","CDKN1C","CFD",
               "PRMT9","RASGEF1B","IFI27","MT1E","GPX4","ARL2BP",
               "THBS1","SLC2A3","SERPINB9","RGS10","FCGR3","NKG7","CYBB","NDRG2",
               "RPS20","LST1","PPA1","RNF213","CD163","DUSP6","WARS1.1")

# diff_gene <- c("S100A8","LTF","NFKBIA","DUSP1","PPBP","CHMP7","IGKC","APOL2",
#               "PRTN3","CAMP","CCL3","IFNG","CX3CR1","RGS9","FCRL6","CSF1",
#               "SOD2","S100A9","FOS","FOSB","GZMB","PRF1","COTL1","FGR",
#               "RTD1A","RASSF4","C1QC","JCHAIN","VMO1","LTA4H","LY6D","LILRA3",
#               "DEFA1B","AHSP","CXCL9","UBD","CDKN1C","CFD","KLRD1","CCL5",
#               "HBB","HBA","PRMT9","RASGEF1B","IFI27","MT1E","GPX4","ARL2BP",
#               "THBS1","SLC2A3","SERPINB9","RGS10","FCGR3","NKG7","CYBB","NDRG2",
#               "RPS20","LST1","ISG15","PPA1","RNF213","CD163","DUSP6","WARS1.1")
  
# diff_gene <- c("S100A8","CAMP","COA1","DUSP1","IL7R","TPT1","IGKC","GIMAP6",
#                "LTF","MMP9","RGS1","JUN","AHNAK","CD52","CSF1","FCRL6",
#                "RTD1A","S100A9","NFKBIA","GADD45B","GTF3C1","RHOBTB3","ITGAX","CCL5",
#                "S100A6","CTSB","IL1B","JCHAIN","AIF1","LYN","ACAD10","PSAP",
#                "CTSD","GLUL","MAMU-DRA","UBD","LTA4H","CFD","CTSC","IFI27",
#                "LTAF","RETN","MAMU-DRB1","IDO1","CASP1","MT1E","FES","ODF3B",
#                "PLBD1","VCAN","IRF8","CIITA","FN1","IGSF6","GRN","RNASE2",
#                "GCA","LST1","ZFP36","SAT1","IFI30","TAK","FGR","BRI3")


diff_gene_overlap_FC <- data.frame()
for (cell in cell_list){
  tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_tissue.txt")) %>% .[.$change_V1 %in% c("UP"),]
  # tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/tissue/tissue1.txt")) %>% .[.$change_V1 %in% c("UP"),]
  DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
  for (i in DE_list){
    old_vs_young <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_",i,"_old_vs_young.txt")) %>% .[.$change_V1 %in% c("UP"),]
    diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
    diff_gene2 <- tissue_DE[tissue_DE$cluster == i,c(2,7)]
    # diff_gene2 <- tissue_DE[tissue_DE$cluster == paste0(i,"_",cell),c(2,7)]
    diff_gene2 <- diff_gene2[diff_gene2$gene %in% diff_gene,]
    diff_gene1 <- merge(diff_gene1,diff_gene2,all.x=TRUE,by = "gene")
    colnames(diff_gene1) <- c("gene","log2FC.old_young","log2FC.tissue")
    diff_gene1$celltype <- rep(cell)
    diff_gene1$Tissue <- rep(i)
    diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
  }
}
diff_gene_overlap_FC <- na.omit(diff_gene_overlap_FC)
# diff_gene_overlap_FC$log2FC.tissue[diff_gene_overlap_FC$log2FC.tissue > 3] = 3
# diff_gene_overlap_FC$log2FC.tissue[diff_gene_overlap_FC$log2FC.tissue < -3] = -3
overlap_FC <- diff_gene_overlap_FC
str(overlap_FC)
overlap_FC$log2FC.tissue <- as.numeric(overlap_FC$log2FC.tissue)
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = cell_list)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
# colors <- c("#00008A","#0000FF","#3636FF","#4D4DFF","#9B9BFF","#F7CECE","#F7C0C0","#F78888","#F76C6C","#981216")
colors <- c("#F7C0C0","#F78888","#F76C6C","#981216")
show_col(colors)
p<-ggplot(overlap_FC,aes(x = Tissue,y=gene,size=log2FC.old_young,colour=log2FC.tissue,ylab=''))+
  geom_point()+
  facet_wrap(.~celltype,ncol= 8)+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradientn(colors = colors)+
  theme_classic()+
  labs(x = '', y = '',title = cell)+
  theme(axis.text.x=element_text(size=12, color="black",angle = 45,vjust = 1,hjust = 1),legend.position = "right",
        axis.text.y=element_text(size=12, color="black"))# + RotatedAxis():45度
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","immune_gene12",".pdf"),p,width = 12,height = 15)

# 构造数据用于添加线段
df3<-data.frame(
  x = seq(1.5,4.5,1),
  xend = seq(1.4,4.5,1),
  y = -Inf,
  yend = Inf
)
df3
df4<-data.frame(
  x = -Inf,
  xend = Inf,
  y = seq(1.5,64.5,1),
  yend = seq(1.5,64.5,1)
)
p <- ggplot(data=overlap_FC,aes(x=Tissue,y=gene))+
  geom_point(aes(size=log2FC.old_young,
                 color=log2FC.tissue),
             shape=15)+
  # scale_color_manual(values = c(rep("#fe0000",4),rep("#009ccc",4)))+
  scale_color_gradientn(colors = colors)+
  facet_wrap(.~celltype,ncol= 8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="grey"),
        axis.ticks = element_blank())+
  geom_segment(data=df3,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey")+
  geom_segment(data=df4,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey")+
  scale_size_continuous(range = c(2,10))+
  # scale_y_discrete(position = "right")+
  labs(x=NULL,y=NULL)
p

# 免疫细胞通路-合并 ----------------------------------------------------------------------
rm(list=ls())
library("ggplot2")

tissue_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO0 <- data.frame()
for (tissue in tissue_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/CD8_",tissue,"_Gene_GO.txt"),stringsAsFactors = F)
  GO1$Tissue <- rep(tissue)
  GO1$cell <- rep("CD8")
  GO0 <- rbind(GO0,GO1)
}
colnames(GO0)
GO0 <- GO0[,-10]
cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
GO <- data.frame()
for (cell in cell_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/",cell,"_GO.txt"),stringsAsFactors = F)
  GO1$cell <- rep(cell)
  GO <- rbind(GO,GO1)
}
colnames(GO)
GO <- rbind(GO,GO0)
ID <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/GO2.txt")
unique(ID$ID)
GO_sel <- GO[GO$ID %in% ID$ID,]
GO_sel$cell <- factor(GO_sel$cell,levels = cell_list)
GO_sel$Description <- factor(GO_sel$Description,levels = rev(unique(ID$Description)))
GO_sel$GeneRatio1 <- parse_ratio(GO_sel$GeneRatio)
unique(GO_sel$Description)

colors <- rev(c("#F78888","#F76C6C","#F75656","#981216"))
p <- ggplot(GO_sel,aes(Description,Tissue,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  # facet_grid(vars(cell.x), vars(cell.y),scales="free")+
  facet_wrap(.~cell,ncol= 8)+
  scale_color_gradientn(colors = colors)+
  # scale_colour_gradient(low="#F7A4A4",high="#981216")+
  # scale_x_discrete(labels=function(x) str_wrap(x, width=60))+
  geom_point(size = 2,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_classic()+ theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 8,angle = 90,vjust = 1,hjust = 1,color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "right")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","immune","_GO.pdf"),p,width = 14,height = 13)

# 单独做CD8-点图 ------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
cell_list <- c("CD4","CD8","NK","M1","M2","CD14_Mon","CD16_Mon","In-termed_Mon")
diff_gene <- c("S100A8","CAMP","COA1","DUSP1","IL7R","TPT1","IGKC","GIMAP6",
               "LTF","MMP9","RGS1","JUN","AHNAK","CD52","CSF1","FCRL6",
               "RTD1A","S100A9","NFKBIA","GADD45B","GTF3C1","RHOBTB3","ITGAX","CCL5",
               "S100A6","CTSB","IL1B","JCHAIN","AIF1","LYN","ACAD10","PSAP",
               "CTSD","GLUL","MAMU-DRA","UBD","LTA4H","CFD","CTSC","IFI27",
               "LTAF","RETN","MAMU-DRB1","IDO1","CASP1","MT1E","FES","ODF3B",
               "PLBD1","VCAN","IRF8","CIITA","FN1","IGSF6","GRN","RNASE2",
               "GCA","LST1","ZFP36","SAT1","IFI30","TAK","FGR","BRI3")
diff_gene_overlap_FC <- data.frame()
for (cell in cell_list){
  # tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_tissue.txt")) %>% .[.$change_V1 %in% c("UP","DOWN"),]
  tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/tissue/tissue1.txt")) %>% .[.$change_V1 %in% c("UP"),]
  DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
  for (i in DE_list){
    old_vs_young <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_",i,"_old_vs_young.txt")) %>% .[.$change_V1 %in% c("UP"),]
    diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene,c(2,6)]
    # diff_gene2 <- tissue_DE[tissue_DE$cluster == i,c(2,7)]
    diff_gene2 <- tissue_DE[tissue_DE$cluster == paste0(i,"_",cell),c(2,7)]
    diff_gene2 <- diff_gene2[diff_gene2$gene %in% diff_gene,]
    diff_gene1 <- merge(diff_gene1,diff_gene2,all.x=TRUE,by = "gene")
    colnames(diff_gene1) <- c("gene","log2FC.old_young","log2FC.tissue")
    diff_gene1$celltype <- rep(cell)
    diff_gene1$Tissue <- rep(i)
    diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
  }
}
diff_gene_overlap_FC <- na.omit(diff_gene_overlap_FC)
diff_gene_overlap_FC$log2FC.tissue[diff_gene_overlap_FC$log2FC.tissue > 3] = 3
diff_gene_overlap_FC$log2FC.tissue[diff_gene_overlap_FC$log2FC.tissue < -3] = -3
overlap_FC <- diff_gene_overlap_FC
str(overlap_FC)
overlap_FC$log2FC.tissue <- as.numeric(overlap_FC$log2FC.tissue)
overlap_FC$celltype <- factor(overlap_FC$celltype,levels = cell_list)
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(diff_gene))
colors <- c("#00008A","#0000FF","#3636FF","#4D4DFF","#9B9BFF","#F7CECE","#F7C0C0","#F78888","#F76C6C","#981216")
# show_col(colors)
p<-ggplot(overlap_FC,aes(x = Tissue,y=gene,size=log2FC.old_young,colour=log2FC.tissue,ylab=''))+
  geom_point()+
  facet_wrap(.~celltype,ncol= 8)+
  scale_size_continuous(range=c(0,6))+
  scale_color_gradientn(colors = colors)+
  theme_classic()+
  labs(x = '', y = '',title = cell)+
  theme(axis.text.x=element_text(size=12, color="black",angle = 45,vjust = 1,hjust = 1),legend.position = "right",
        axis.text.y=element_text(size=12, color="black"))# + RotatedAxis():45度
p
# ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","immune_gene12",".png"),p,width = 20,height = 20)

# 构造数据用于添加线段
df3<-data.frame(
  x = seq(1.5,4.5,1),
  xend = seq(1.4,4.5,1),
  y = -Inf,
  yend = Inf
)
df3
df4<-data.frame(
  x = -Inf,
  xend = Inf,
  y = seq(1.5,64.5,1),
  yend = seq(1.5,64.5,1)
)
p <- ggplot(data=overlap_FC,aes(x=Tissue,y=gene))+
  geom_point(aes(size=log2FC.old_young,
                 color=log2FC.tissue),
             shape=15)+
  # scale_color_manual(values = c(rep("#fe0000",4),rep("#009ccc",4)))+
  scale_color_gradientn(colors = colors)+
  facet_wrap(.~celltype,ncol= 8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="grey"),
        axis.ticks = element_blank())+
  geom_segment(data=df3,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey")+
  geom_segment(data=df4,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey")+
  scale_size_continuous(range = c(2,10))+
  # scale_y_discrete(position = "right")+
  labs(x=NULL,y=NULL)
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","immune_FC",".png"),p,width = 15,height = 15)

# 单独做CD8-富集 ------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
GO <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/","CD8","_GO.txt"),stringsAsFactors = F)
d <- data.frame(table(GO$ID))
GO_sel <- GO[GO$ID %in% d[d$Freq == 1,1],]
write.table(GO_sel, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/","CD8","_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
# ID <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/GO2.txt")
# GO_sel <- GO[GO$ID %in% ID$ID,]
# GO_sel$cell <- factor(GO_sel$cell,levels = cell_list)
GO_sel$Description <- factor(GO_sel$Description,levels = rev(unique(GO_sel$Description)))
# GO_sel$GeneRatio1 <- parse_ratio(GO_sel$GeneRatio)
# unique(GO_sel$Description)
p <- ggplot(GO_sel,aes(Description,Tissue,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_colour_gradient(low="#F7A4A4",high="#981216")+
  # scale_x_discrete(labels=function(x) str_wrap(x, width=60))+
  geom_point(size = 2,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 8,angle = 90,vjust = 1,hjust = 1,color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "right")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","immune","_GO.png"),p,width = 15,height = 100,limitsize = F)

# 单独做CD8-指定基因富集 ------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
gene <- read.delim2("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/CD8_GO_Gene.txt")
library(tidyr)
gene1 <- gene %>% separate_rows(gene, sep = ",")
colnames(gene1)[1] <- "SYMBOL"
df_id <- bitr(gene1$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
gene1 <- merge(gene1,df_id,by = "SYMBOL",all=F)
gene <- unique(gene1[,'ENTREZID'])
ego_BP1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_BP1 <- setReadable(ego_BP1, OrgDb = org.Mmu.eg.db);ego_BP_result1<-as.data.frame(ego_BP1@result)
ego_BP_result1<-ego_BP_result1[ego_BP_result1$pvalue<0.05,]
ego_BP_result1$Group <-rep("BP")

ego_MF1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "MF",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_MF1 <- setReadable(ego_MF1, OrgDb = org.Mmu.eg.db);ego_MF_result1<-as.data.frame(ego_MF1@result)
ego_MF_result1<-ego_MF_result1[ego_MF_result1$pvalue<0.05,]
ego_MF_result1$Group <-rep("MF")

ego_CC1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "CC",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
ego_CC1 <- setReadable(ego_CC1, OrgDb = org.Mmu.eg.db);ego_CC_result1<-as.data.frame(ego_CC1@result)
ego_CC_result1<-ego_CC_result1[ego_CC_result1$pvalue<0.05,]
ego_CC_result1$Group <-rep("CC")
ego_result <- rbind(ego_BP_result1,ego_MF_result1,ego_CC_result1)
write.table(ego_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/","CD8","_Gene_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)

GO_lists <- c("GO:0002285","GO:0097696","GO:0002292","GO:0001817","GO:0030155",
              "GO:0002252","GO:0002381","GO:0031341","GO:0042129","GO:1990869",
              "GO:1903659","GO:0071706","GO:0002439","GO:0051712","")

# 单独做CD8-四个组织年老年轻差异基因富集 ------------------------------------------------------------------
rm(list=ls())
tissue_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
for (tissue in tissue_list){
  de <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/","CD8","_",tissue,"_old_vs_young.txt")) %>% .[.$change_V1 == "UP",]
  colnames(de)[6] <- "SYMBOL"
  df_id <- bitr(de$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mmu.eg.db")
  gene1 <- merge(de,df_id,by = "SYMBOL",all=F)
  gene <- unique(gene1[,'ENTREZID'])
  ego_BP1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
  ego_BP1 <- setReadable(ego_BP1, OrgDb = org.Mmu.eg.db);ego_BP_result1<-as.data.frame(ego_BP1@result)
  ego_BP_result1<-ego_BP_result1[ego_BP_result1$pvalue<0.05,]
  ego_BP_result1$Group <-rep("BP")
  
  ego_MF1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "MF",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
  ego_MF1 <- setReadable(ego_MF1, OrgDb = org.Mmu.eg.db);ego_MF_result1<-as.data.frame(ego_MF1@result)
  ego_MF_result1<-ego_MF_result1[ego_MF_result1$pvalue<0.05,]
  ego_MF_result1$Group <-rep("MF")
  
  ego_CC1 <- enrichGO(gene = gene,OrgDb=org.Mmu.eg.db,ont = "CC",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = 0.05)
  ego_CC1 <- setReadable(ego_CC1, OrgDb = org.Mmu.eg.db);ego_CC_result1<-as.data.frame(ego_CC1@result)
  ego_CC_result1<-ego_CC_result1[ego_CC_result1$pvalue<0.05,]
  ego_CC_result1$Group <-rep("CC")
  ego_result <- rbind(ego_BP_result1,ego_MF_result1,ego_CC_result1)
  write.table(ego_result, file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/","CD8","_",tissue,"_Gene_GO.txt"), quote = FALSE,sep="\t",row.names = FALSE)
}

# 单独做CD8-富集通路图-合并 -------------------------------------------------------------------
rm(list=ls())
# Bone_marrow <- read.delim("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/CD8_Bone_marrow_Gene_GO.txt")
# Bone_marrow1 <- Bone_marrow[grepl("inflam", Bone_marrow$Description),]
# Bone_marrow2 <- Bone_marrow[grepl("damage", Bone_marrow$Description),]
# Bone_marrow3 <- Bone_marrow[grepl("age-", Bone_marrow$Description),]
# Bone_marrow4 <- Bone_marrow[grepl("inhibitory MHC", Bone_marrow$Description),]
# Bone_marrow5 <- Bone_marrow[grepl("^antigen processing and presentation of peptide antigen via MHC class I", Bone_marrow$Description),]
# Bone_marrow6 <- Bone_marrow[grepl("regulation of immune effector process", Bone_marrow$Description),]
# Bone_marrow7 <- Bone_marrow[grepl("T cell mediated cytotoxicity", Bone_marrow$Description),]
# Bone_marrow8 <- Bone_marrow[grepl("positive regulation of immune effector process", Bone_marrow$Description),]
# Bone_marrow9 <- Bone_marrow[Bone_marrow$ID %in% GO_lists,]
# d <- rbind(Bone_marrow1,Bone_marrow2,Bone_marrow3,Bone_marrow4,Bone_marrow5,Bone_marrow6,Bone_marrow7,Bone_marrow8)
# GO_lists1 <- unique(d$ID)
GO_lists1 <- c("GO:0002544","GO:0006954","GO:0042771","GO:1902165","GO:0008630","GO:0043518","GO:2001021","GO:1904290","GO:2001020",
              "GO:0001315","GO:0062080","GO:0062082","GO:0002474","GO:0002699","GO:0002697","GO:0001913","GO:0001914")
GO_lists2 <- c("GO:0002285","GO:0097696","GO:0002292","GO:0001817","GO:0030155",
              "GO:0002252","GO:0002381","GO:0031341","GO:0042129","GO:1990869",
              "GO:1903659","GO:0071706","GO:0002439","GO:0051712")
GO_lists <- c(GO_lists1,GO_lists2)

library("ggplot2")
tissue_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO <- data.frame()
for (tissue in tissue_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/CD8_",tissue,"_Gene_GO.txt"),stringsAsFactors = F)
  GO1$Tissue <- rep(tissue)
  GO <- rbind(GO,GO1)
}

GO_sel <- GO[GO$ID %in% GO_lists,]
GO_sel$Description <- factor(GO_sel$Description,
                             levels = rev(c("inflammatory response","HLA-E specific inhibitory MHC class Ib receptor activity",
                                            "inhibitory MHC class Ib receptor activity",
                                            "positive regulation of killing of cells of other organism","T cell mediated cytotoxicity",
                                            "cellular response to chemokine","tumor necrosis factor superfamily cytokine production",
                                            "regulation of cell killing","regulation of cell adhesion",
                                            "positive regulation of immune effector process","regulation of immune effector process",
                                            "immune effector process",
                                            "regulation of cytokine production","regulation of T cell mediated cytotoxicity",
                                            "negative regulation of response to DNA damage stimulus",
                                            "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
                                            "regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
                                            "intrinsic apoptotic signaling pathway in response to DNA damage",
                                            "negative regulation of DNA damage response, signal transduction by p53 class mediator",
                                            "regulation of T cell proliferation",
                                            "regulation of response to DNA damage stimulus",
                                            "negative regulation of mitotic DNA damage checkpoint","chronic inflammatory response",
                                            "age-dependent response to reactive oxygen species",
                                            "antigen processing and presentation of peptide antigen via MHC class I",
                                            "lymphocyte activation involved in immune response")))
GO_sel$GeneRatio1 <- parse_ratio(GO_sel$GeneRatio)
unique(GO_sel$Description)
colors <- c("#F78888","#F76C6C","#F75656","#981216")
p <- ggplot(GO_sel,aes(reorder(Description,Tissue),Tissue,color = pvalue, size=Count)) +
  geom_point()+
  coord_flip() +
  scale_color_gradientn(colors = rev(colors))+
  # scale_colour_gradient(low="#F7CECE",high="#981216")+
  # scale_x_discrete(labels=function(x) str_wrap(x, width=60))+
  geom_point(size = 2,shape = 16)+
  labs(x = "", y = "", title = "") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(size = 8,angle = 90,vjust = 1,hjust = 1,color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "right")
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","CD8","_GO.pdf"),p,width = 8,height = 10)

# 单独做CD8-通路基因-点图-合并 ----------------------------------------------------------------------
rm(list=ls())
library("ggplot2")
GO_lists <- c("GO:0002544","GO:0006954","GO:0042771","GO:1902165","GO:0008630","GO:0043518","GO:2001021","GO:1904290","GO:2001020",
               "GO:0001315","GO:0062080","GO:0062082","GO:0002474","GO:0002699","GO:0002697","GO:0001913","GO:0001914",
               "GO:0002285","GO:0097696","GO:0002292","GO:0001817","GO:0030155",
               "GO:0002252","GO:0002381","GO:0031341","GO:0042129","GO:1990869",
               "GO:1903659","GO:0071706","GO:0002439","GO:0051712")
tissue_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
GO <- data.frame()
for (tissue in tissue_list){
  GO1 <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/CD8_",tissue,"_Gene_GO.txt"),stringsAsFactors = F)
  GO1$Tissue <- rep(tissue)
  GO <- rbind(GO,GO1)
}

GO_sel <- GO[GO$ID %in% GO_lists,]
d <- unique(GO_sel[,c(1,2)])
write.table(d,file = paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/immune_cell/GO3.txt"),quote = FALSE,sep = "\t",row.names = F)

diff_gene <- GO_sel[,c(2,8)]
diff_gene <- diff_gene %>% separate_rows(geneID, sep = "/")

diff_gene_overlap_FC <- data.frame()
cell = "CD8"
tissue_DE <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_tissue.txt")) %>% .[.$change_V1 %in% c("UP"),]
DE_list <- c("Bone_marrow","Mesenteric_lymph","PBMC","Spleen")
for (i in DE_list){
  old_vs_young <- read.delim(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/DE_immune/",cell,"_",i,"_old_vs_young.txt")) %>% .[.$change_V1 %in% c(c("UP")),]
  diff_gene1 <- old_vs_young[old_vs_young$gene %in% diff_gene$geneID,c(2,6)]
  diff_gene2 <- tissue_DE[tissue_DE$cluster == i,c(2,7)]
  diff_gene2 <- diff_gene2[diff_gene2$gene %in% diff_gene$geneID,]
  diff_gene1 <- merge(diff_gene1,diff_gene2,all.x=TRUE,by = "gene")
  colnames(diff_gene1) <- c("gene","log2FC.old_young","log2FC.tissue")
  diff_gene1$celltype <- rep(cell)
  diff_gene1$Tissue <- rep(i)
  diff_gene_overlap_FC <- rbind(diff_gene_overlap_FC,diff_gene1) 
}

diff_gene_overlap_FC <- na.omit(diff_gene_overlap_FC)
# diff_gene_overlap_FC$log2FC.tissue[diff_gene_overlap_FC$log2FC.tissue > 3] = 3
overlap_FC <- diff_gene_overlap_FC
overlap_FC$log2FC.tissue <- as.numeric(overlap_FC$log2FC.tissue)
gene_lists <- c("S100A9","LTF","SOD2","LGALS3","IFI16","VSIR",
                "CCL3","IFNG","NFKBIA","TNFAIP3","DUSP1","ISG15","CCL4L1","DUSP10","CD83","ZFP36","ZC3H12A","IRF4","CD74",
                "AFAP1L2","CD226","PRF1","NCR3","SLC11A1","PPBP","SYT11","VCL","ADGRG1","CIB1","FLNA","CDC42","NRROS",
                "KLRD1","LAG3","CCL5","UBASH3B","FGR","ISG20","LYN","SWAP70","BST2")
overlap_FC$gene <- factor(overlap_FC$gene,levels = rev(gene_lists))
colors <- c("#F78888","#F76C6C","#F75656","#981216")
p<-ggplot(overlap_FC,aes(x = Tissue,y=gene,size=log2FC.old_young,colour=log2FC.tissue,ylab=''))+
  geom_point()+
  scale_size_continuous(range=c(0,4))+
  scale_color_gradientn(colors = colors)+
  theme_classic()+
  labs(x = '', y = '',title = "")+
  theme(axis.text.x=element_text(size=12, color="black",angle = 90,vjust = 1,hjust = 1),legend.position = "right",
        axis.text.y=element_text(size=12, color="black"))# + RotatedAxis():45度
p
ggsave(paste0("/data1/zhuzn/wangsn/reviewer_4000/figure2/","CD8","_GO_gene.pdf"),p,width = 4,height = 8)
