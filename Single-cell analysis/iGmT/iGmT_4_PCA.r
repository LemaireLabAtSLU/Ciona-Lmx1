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
library("leidenbase")
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

wd1 <- "/YourDirectory"
setwd(wd1)

object <- "iGmT3.RData"

load(object)

iGmT

object <- "iGmT4.RData"

#Choosing number of PCAs for clustering

ElbowPlot(iGmT, ndims = 40, reduction = "pca") ###Cannot check integrated.cca
ggsave("Elbowplot.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

iGmT@meta.data[["orig.ident"]] <- factor(iGmT@meta.data[["orig.ident"]], levels = ids1)
Idents(iGmT) <- "orig.ident"

DimPlot(iGmT, reduction ="pca", split.by="orig.ident")+NoLegend()

n.dims <- 16 ###Number of dimension for downstream analysis

iGmT <- FindNeighbors(iGmT, reduction = "integrated.cca", dims = 1:n.dims)

resolutions <- seq(1.0, 2.0, 0.2)

iGmT <- FindClusters(iGmT,  algorithm = 4, resolution = resolutions) ###(Resolution should be at least 1for leiden)

save.image(object)



pdf(file="cluster-tree.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(iGmT)
dev.off()

DimPlot(iGmT, reduction ="pca", group.by="Phase")

DimPlot(iGmT, reduction ="integrated.cca", cells.highlight = WhichCells(iGmT, idents = "24"),
        cols.highlight = "darkblue")

DimPlot(iGmT, reduction ="integrated.cca", group.by="Phase")


#Project cells into UMAP and tSNE space

iGmT <- RunTSNE(iGmT, reduction = "integrated.cca", dims = 1:n.dims) 

iGmT <- RunUMAP(iGmT, reduction = "integrated.cca",
                   dims = 1:n.dims, reduction.name = "UMAP")

stage <- c("iG","iG",
           "mG", "mG",
           "eN", "eN",
           "lN", "lN",
           "iT", "iT",
           "eT", "eT",
           "mT", "mT")

stage <- data.frame(ids1, stage)

colnames(stage)[1] <- "orig.ident"

metadata <- iGmT@meta.data
metadata <- rownames_to_column(metadata)
metadata <- merge(metadata, stage, by = "orig.ident")
metadata <- column_to_rownames(metadata)
md1 <- subset(metadata, select = "stage")
head(md1)
iGmT <- AddMetaData(iGmT, md1)

stage.order <- c("iG", "mG", "eN", "lN", "iT", "eT", "mT")
iGmT$stage <- factor(iGmT$stage,
                        levels = stage.order)

col1 <- c("#35608D", "#66CC33", "#E31A1C")
col2 <- c("lightgrey", "#330066")
col3 <- c("blue", "red")
nb.cols <- 34 #Number of clusters
col4 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 14 #Number of samples
col5 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 7 #Number of stage
col6 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
barplot(1:34, col = col4)






DimPlot(iGmT, reduction = "tsne", 
        group.by = "orig.ident", cols = alpha(col5, 0.5))+
  coord_fixed()
ggsave("tSNE_orig.ident.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "orig.ident", cols = alpha(col5, 0.5), label.size = 2, pt.size = 1)+
  coord_fixed()
ggsave("UMAP_orig.ident.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(iGmT, reduction = "UMAP",
        group.by = "stage",
        cols = col6,
        pt.size = 0.5)+
  coord_fixed()
ggsave("UMAP_stage.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")

DimPlot(iGmT, reduction = "tsne",
        group.by = "stage",
        cols = col6,
        pt.size = 0.5)+
  coord_fixed()
ggsave("tsne_stage.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")

DimPlot(iGmT, reduction = "UMAP",
        group.by = "Phase", split.by = "stage")
ggsave("UMAP_phase.stage.pdf", device= "pdf", width = 40, 
       height = 16, units = "cm")


save.image(object)

#### Check "essential markers"

FeaturePlot(iGmT, features = c("KY21:KY21.Chr7.872",#Usx _ Epidermis
                               "KY21:KY21.Chr11.1104", #thr _ Endoderm
                               "KY21:KY21.Chr12.6", #brachyury _ Notochord
                               "KY21:KY21.Chr9.100", #twist _ Mesenchyme
                               "KY21:KY21.Chr3.978", #MesP _ Tvcs
                               "KY21:KY21.Chr14.750", #Mrf _ Muscle
                               "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl _ germ cells
                               "KY21:KY21.Chr6.58", #ETR = Celf3
                               "KY21:KY21.Chr5.707", #Dmrt1 _a-lineage
                               "KY21:KY21.Chr12.803", #Tyr _Pigment cells
                             "KY21:KY21.Chr6.345", #Flvcr2_ A-lineage
                             "KY21:KY21.Chr2.1031",#Msx_ b-lineage
                             "KY21:KY21.Chr9.822", #Wnttun5 _ Rubber cells
                             "KY21:KY21.Chr1.1701" #hand2 _ Tvcs, Rubber cells
                             ), 
            reduction = "UMAP", ncol = 5, cols = col2,
            min.cutoff = "q05", max.cutoff = "q95")
ggsave("iGmT.umap_tissuemarkers.pdf", device= "pdf", width = 45, 
       height = 27, units = "cm") 

FeaturePlot(iGmT, features = c("KY21:KY21.Chr7.872",#Usx _ Epidermis
                               "KY21:KY21.Chr11.1104", #thr _ Endoderm
                               "KY21:KY21.Chr12.6", #brachyury _ Notochord
                               "KY21:KY21.Chr9.100", #twist _ Mesenchyme
                               "KY21:KY21.Chr3.978", #MesP _ Tvcs
                               "KY21:KY21.Chr14.750", #Mrf _ Muscle
                               "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl _ germ cells
                               "KY21:KY21.Chr6.58", #ETR = Celf3
                               "KY21:KY21.Chr5.707", #Dmrt1 _a-lineage
                               "KY21:KY21.Chr12.803", #Tyr _Pigment cells
                               "KY21:KY21.Chr6.345", #Flvcr2_ A-lineage
                               "KY21:KY21.Chr2.1031",#Msx_ b-lineage
                               "KY21:KY21.Chr9.822", #Wnttun5 _ Rubber cells
                               "KY21:KY21.Chr1.1701" #hand2 _ Tvcs, Rubber cells
), 
reduction = "tsne", ncol = 5, cols = col2,
min.cutoff = "q05", max.cutoff = "q95")


DimPlot(iGmT, reduction = "UMAP", 
        cols = alpha(col4, 0.5), group.by = "RNA_snn_res.1",
        pt.size = 0.5, label = TRUE, raster = FALSE)+
  coord_fixed()

Idents(iGmT) <- "RNA_snn_res.1"

DimPlot(iGmT, cells.highlight = WhichCells(iGmT, idents = "6"),
        cols.highlight = "darkblue", reduction = "UMAP")+
  coord_fixed()


FeaturePlot(iGmT, features = c("KY21:KY21.Chr5.707"), #Dmrt1 _a-lineage
reduction = "UMAP", cols = col2,
min.cutoff = "q05", max.cutoff = "q95")+
  coord_fixed()

###Select 1 as resolution

iGmT$seurat_clusters <- iGmT$RNA_snn_res.1

Idents(iGmT) <- "seurat_clusters"


DimPlot(iGmT, reduction = "UMAP", 
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE, raster = FALSE)+
  coord_fixed()
ggsave("UMAP_clusters.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(iGmT, reduction = "tsne", 
        cols = alpha(col4, 0.5),
        pt.size = 0.2, label = TRUE, raster = FALSE)+
  coord_fixed()
ggsave("tsne_clusters.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

save.image(object)
