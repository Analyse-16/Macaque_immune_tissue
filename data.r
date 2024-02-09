# ############# 4. Filter(200-4000) ############# --------------------------------------------------------------
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
load(file = "/data1/zhuzn/wangsn/V4/data/merge_all_sample_seurat.RData")
ncol(alldata)
unique(alldata@meta.data$orig.ident)
del_all <- WhichCells(alldata,expression = (nFeature_RNA>4000 | nFeature_RNA<200 | percent.mito>5))
data.filt <- subset(alldata,cells=setdiff(WhichCells(alldata),del_all))
ncol(data.filt)
## Normalize & Find variable features 
data.filt.list <- SplitObject(data.filt, split.by = "orig.ident")
data.filt.list <- lapply(X = data.filt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})
## Perform integration & cluster 
all_sample.anchors <- FindIntegrationAnchors(object.list = data.filt.list , dims = 1:20)
all_sample.combined <- IntegrateData(anchorset = all_sample.anchors, dims = 1:20)
save(all_sample.combined, file="/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000.RData")

DefaultAssay(all_sample.combined) <- "integrated"
all_sample.combined <- ScaleData(all_sample.combined,verbose = T)
all_sample.combined <- RunPCA(all_sample.combined, npcs = 30, verbose = T)
# t-SNE and Clustering
all_sample.combined <- RunUMAP(all_sample.combined, seed.use = -1, reduction = "pca", dims = 1:20)
all_sample.combined <- RunTSNE(all_sample.combined, seed.use = 2, reduction = "pca", dims = 1:20)
all_sample.combined <- FindNeighbors(all_sample.combined, reduction = "pca", dims = 1:20)
all_sample.combined <- FindClusters(all_sample.combined, resolution = 1)
p <- DimPlot(all_sample.combined, reduction = "umap", group.by = "integrated_snn_res.1", label = TRUE,raster=FALSE)
p
ggsave(paste0("/data1/zhuzn/wangsn/data_QC_test/","test_4000.pdf"),p,dpi = 1000,width = 9,height = 8)
save(all_sample.combined, file="/data1/zhuzn/wangsn/data_QC_test/Integrate_test_4000_cluster.RData")
