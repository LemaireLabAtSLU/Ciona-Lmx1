library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(patchwork)
#library(plyr)
library(dbplyr)
library(pheatmap)
library(tidyft)
library(stringr)
library(tidyverse)
library(splitstackshape)
library(data.table)
library(viridis)
library("clustree")
library("ComplexHeatmap")
library("zoo")
library(jcolors)
library(circlize)
library(RColorBrewer)

library(slingshot)
library(SingleCellExperiment)
library(scales)
library(UpSetR)
library(msigdbr)
library(fgsea)
library(knitr)
library(gridExtra)
library(tradeSeq)
library(clusterExperiment)
library("DelayedMatrixStats")

#################### Script performed with SLURM

print("library loaded")

wd1 <- "/YourDirectory"
setwd(wd1)

object <- "blineage_3.RData"
load(object)
object <- "blineage_4.RData"


#Clustering of genes according to their expression pattern
##Clustering using RSEC, clusterExperiment
# Tranfer reduce dimension to sce object with tradeseq
reducedDims(npbb.sce.2)$PCA <- npbb.sce@int_colData@listData[["reducedDims"]]@listData[["NPBB2.PCA"]]

head(npbb.sce@int_colData@listData[["reducedDims"]]@listData[["NPBB2.PCA"]])
head(reducedDims(npbb.sce.2)$PCA)
reducedDimNames(npbb.sce.2)
listBuiltInReducedDims()
class(npbb.sce.2)

nPointsClus <- 40
clusPat1 <- clusterExpressionPatterns(npbb.sce.2, nPoints = nPointsClus, reduceMethod = "PCA", verbose = TRUE,
                                      genes = deg.lineage1, minSizes = 3,  ncore = 4)

print(primaryCluster(clusPat1$rsec))
print(tableClusters(clusPat1$rsec))

clusPat2 <- clusterExpressionPatterns(npbb.sce.2, nPoints = nPointsClus, reduceMethod = "PCA", verbose = TRUE,
                                      genes = deg.lineage2, , minSizes = 3,  ncore = 4)

print(primaryCluster(clusPat2$rsec))
print(tableClusters(clusPat2$rsec))

clusPat3 <- clusterExpressionPatterns(npbb.sce.2, nPoints = nPointsClus, reduceMethod = "PCA", verbose = TRUE,
                                      genes = deg.lineage3, minSizes = 3,  ncore = 4)

print(primaryCluster(clusPat3$rsec))
print(tableClusters(clusPat3$rsec))

save.image(object)

print("End Script")
