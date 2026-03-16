library(cowplot)
library(Matrix)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)
library(plyr)
library(tidyverse)
library(tidyft)
library("ComplexHeatmap")
library("clustifyr")
library("leiden")
library("limma")
library(data.table)
library("RColorBrewer")
library(scico)
# Plotting
library(ggraph)
library("clustree")
library(viridisLite)
library("viridis")
library(circlize)




wd1 <- "/scratch/gpfs/LEVINE/llemaire/singleCellAnalysis/26.01.26_iGmT"
setwd(wd1)
object <- "26.01.26_iGmT2.RData"
load(object)
print("Environment Loaded")

object <- "26.01.27_iGmT3.RData" ###Should be after load object, wrote over by load over
inf <- sessionInfo()
write(capture.output(inf),file = "26.01.27_sessionInfo_script3.txt")

# split the dataset into a list of layer within the same object

iGmT[["RNA"]] <- split(iGmT[["RNA"]], f = iGmT$orig.ident)

# run standard anlaysis workflow
iGmT <- FindVariableFeatures(iGmT, selection.method = "vst", nfeatures = 4000)
iGmT <- ScaleData(iGmT, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(iGmT))

print("Data are scaled")


iGmT <- RunPCA(iGmT)
DimPlot(iGmT, reduction = "pca")
ggsave("PCA.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

save.image(object)

### Cluster and project without integration

iGmT <- FindNeighbors(iGmT, dims = 1:30, reduction = "pca")
iGmT <- FindClusters(iGmT, resolution = 0.6, cluster.name = "unintegrated_clusters")
iGmT <- RunUMAP(iGmT, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(iGmT, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))
ggsave("Unintegrated.UMAP.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")

print("Unintegrated data are projected")

#Perform integration

iGmT <- IntegrateLayers(object = iGmT, method = CCAIntegration, orig.reduction = "pca",
                           new.reduction = "integrated.cca",
                           verbose = TRUE)


iGmT
print("Integration Done")
# re-join layers after integration

iGmT[["RNA"]] <- JoinLayers(iGmT[["RNA"]])

iGmT
print("Layers joined")
save.image(object)