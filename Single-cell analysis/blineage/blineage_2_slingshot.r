sessionInfo()

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

wd1 <- "/scratch/gpfs/LEVINE/llemaire/singleCellAnalysis/26.01.26_iGmT/26.02.03_NS/26.02.08_blineage/26.02.10_subset npbb"
setwd(wd1)
write(capture.output(sessionInfo()),file = "26.02.10sessionInfo.txt")

object <- "26.02.10blineage_2.RData"

npbb <- readRDS( file = "blineage2.rds")



################################################3##Slingshot: trajectory inference
col1b <- c("#35608D", "#66CC33", "#E31A1C")
col3 <- c("lightgrey", "#330066")
nb.cols <- 14
col4 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 5
col5 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 7
col6 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
col7 <- c("#A6CEE3" ,"#1F78B4","#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A", "#FB9A99", "#E31A1C","#FDBF6F", "#FF7F00","#B15928","#E8D079")
# col9 <- colorRampPalette(c("#A6CEE3" ,"#1F78B4"))(6)

col1 = viridis(8, option = "viridis")
pal <- colorRampPalette(col1)
col2 <- pal(4)

show_col(col7)

#Plot are not supporting proper  color, it has to be implemented via this function
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}



npbb.sce <- as.SingleCellExperiment(npbb, assay = "RNA")

head(npbb.sce@int_colData$reducedDims$NPBB2.UMAP)
head(npbb.sce@int_colData$reducedDims$NPBB2.PCA)

dim(npbb.sce)
class(npbb.sce)
head(npbb.sce$cell.type.4)
head(npbb.sce$npbb2.RNA_snn_res.0.9)


cell_colors <- cell_pal(npbb.sce$cell.type.4, brewer_pal("qual", "Paired"))
cell_colors2 <- cell_pal(npbb.sce$npbb2.RNA_snn_res.0.9, brewer_pal("qual", "Paired"))


plot(reducedDims(npbb.sce)$NPBB2.PCA, col = cell_colors2, 
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )

plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2, 
     asp = 1, pch = 16, xlab = "umap-1", ylab = "umap-2" )

##Perform slinggshot using PCAs

#REsolution npbb2_snn_res.0.9 with start cluster 4
set.seed(1)
lin1 <- getLineages(reducedDims(npbb.sce)$NPBB2.PCA, 
                    clusterLabels = colData(npbb.sce)$npbb2.RNA_snn_res.0.9,
                    start.clus = 4)


plot(reducedDims(npbb.sce)$NPBB2.PCA,
     col = cell_colors2,
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )
lines(SlingshotDataSet(lin1), type = 'lineages', lwd=2, col='black', show.constraints = TRUE) 

cur1 <- getCurves(lin1, approx_points = FALSE,
                  thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

pdf("slingshot_pca_start_cluster4_res0.9.pdf",
    width = 6,
    height = 6)
plot(reducedDims(npbb.sce)$NPBB2.PCA, col = cell_colors2,
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )
lines(SlingshotDataSet(cur1), type = 'curves', lwd=2, col='black') 
dev.off()


cur1_UMAP <- embedCurves(cur1, reducedDims(npbb.sce)$NPBB2.UMAP)

pdf("slingshot_pca_start_cluster4_res0.9_UMAP.pdf",
    width = 6,
    height = 6)
plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2" )
lines(SlingshotDataSet(cur1_UMAP), type = 'curves', lwd=2, col='black') 
dev.off()


#### If no start
set.seed(1)
lin2 <- getLineages(reducedDims(npbb.sce)$NPBB2.PCA, 
                    clusterLabels = colData(npbb.sce)$npbb2.RNA_snn_res.0.9)


plot(reducedDims(npbb.sce)$NPBB2.PCA,
     col = cell_colors2,
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )
lines(SlingshotDataSet(lin2), type = 'lineages', lwd=2, col='black', show.constraints = TRUE) 

cur2 <- getCurves(lin1, approx_points = FALSE)

plot(reducedDims(npbb.sce)$NPBB2.PCA, col = cell_colors2,
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )
lines(SlingshotDataSet(cur2), type = 'curves', lwd=2, col='black') 

cur2_UMAP <- embedCurves(cur2, reducedDims(npbb.sce)$NPBB2.UMAP)

plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2" )
lines(SlingshotDataSet(cur1_UMAP), type = 'curves', lwd=2, col='black') 


####Put start or not very similar with three lineages

##Look at pseudotime using curve1
nc <- 3
pt1 <- slingPseudotime(cur1)
nms <- colnames(pt1)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)

pdf("slingshot_pca_cur1start_cluster4_res0.9_pseudotime.pdf",
    width = 24,
    height = 8)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt1[,i], breaks = 100)]
  plot(reducedDims(npbb.sce)$NPBB2.PCA, col = colors, pch = 16, cex = 1, main = i,
       xlab = "PC-1", ylab = "PC-2", asp = 1,)
  lines(SlingshotDataSet(cur1), lwd = 2, col = 'black', type = 'curves')
}
dev.off()

pdf("slingshot_pca_start_cluster4_res0.9_pseudotime_UMAP.pdf",
    width = 24,
    height = 8)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt1[,i], breaks = 100)]
  plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = colors, pch = 16, cex = 1, main = i,
       xlab = "PC-1", ylab = "PC-2", asp = 1)
  lines(SlingshotDataSet(cur1_UMAP), lwd = 2, col = 'black', type = 'curves')
}
dev.off()

### Try with low resolution res0.3
set.seed(1)
lin3 <- getLineages(reducedDims(npbb.sce)$NPBB2.PCA, 
                    clusterLabels = colData(npbb.sce)$npbb2.RNA_snn_res.0.3)


plot(reducedDims(npbb.sce)$NPBB2.PCA,
     col = cell_colors2,
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )
lines(SlingshotDataSet(lin3), type = 'lineages', lwd=2, col='black', show.constraints = TRUE) 

cur3 <- getCurves(lin3, approx_points = FALSE)

plot(reducedDims(npbb.sce)$NPBB2.PCA, col = cell_colors2,
     asp = 1, pch = 16, xlab = "PC-1", ylab = "PC-2" )
lines(SlingshotDataSet(cur3), type = 'curves', lwd=2, col='black') 

cur3_UMAP <- embedCurves(cur3, reducedDims(npbb.sce)$NPBB2.UMAP)

plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2" )
lines(SlingshotDataSet(cur3_UMAP), type = 'curves', lwd=2, col='black') 


#### Does not work since mTB cells used as start.


#Use UMAP instead of PCA (However not recommended)

#npbb_snn_res.0.9 to have a proper start cluster
set.seed(1)
lin4 <- getLineages(reducedDims(npbb.sce)$NPBB2.UMAP,
                    clusterLabels = colData(npbb.sce)$npbb2.RNA_snn_res.0.9,
                    start.clus = 4)


plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2,
     asp = 0.5, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2")
lines(SlingshotDataSet(lin4), type = 'lineages', lwd=2, col='black', show.constraints = TRUE) 

cur4 <- getCurves(lin4, approx_points = FALSE,
                  thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

pdf("slingshot_umap_cluster4_res0.9.pdf",
    width = 6,
    height = 6)
plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2,
     asp = 1, pch = 16, xlab = "umap-1", ylab = "umap-2" )
lines(SlingshotDataSet(cur4), type = 'curves', lwd=2, col='black') 
dev.off()

cell_colors3 <- cell_pal(npbb.sce$stage, brewer_pal("qual","YlGnBu", direction=+1))

plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors3,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2" )
lines(SlingshotDataSet(cur4), type = 'curves', lwd=2, col='black') 

ggsave("slingshot_umap_cluster4_res0.9_stage.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

nc <- 3
pt4 <- slingPseudotime(cur4)
nms <- colnames(pt4)
nr <- ceiling(length(nms)/nc)
pal <- plasma(100, end = 0.95)
par(mfrow = c(nr, nc))

pdf("slingshot_umap_cluster4_res0.9_pseudotime_lineage.pdf",
    width = 6,
    height = 6)
for (i in nms) {
  colors <- pal[cut(pt4[,i], breaks = 100)]
  plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = colors, pch = 19, cex = 0.5, main = i,
       xlab = "UMAP-1", ylab = "UMAP-2", asp = 1)
  lines(SlingshotDataSet(cur4), lwd = 2, col = 'black', type = 'curves')
}
dev.off()


plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors2,
     asp = 1, pch = 19, xlab = "UMAP-1", ylab = "UMAP-2", main = "stage", cex = .5 )
lines(SlingshotDataSet(cur4), type = 'curves', lwd=2, col='black') 


show_col(viridis_pal()(7))
cell_colors4 <- viridis(7, option = "viridis")
pal <- colorRampPalette(cell_colors4)
cell_colors5 <- pal(7)
show_col(cell_colors5)

pdf("slingshot_umap_cluster4_res0.9_stage2.pdf",
    width = 6,
    height = 6)
plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = cell_colors5[colData(npbb.sce)$stage],
     asp = 1, pch = 19, xlab = "UMAP-1", ylab = "UMAP-2", main = "stage", cex = .5 )
legend("topright", legend = levels(factor(colData(npbb.sce)$stage)), col = cell_colors5,
       pch = 16, bty = "n", cex = 1)
lines(SlingshotDataSet(cur4), type = 'curves', lwd=2, col='black') 
dev.off()


nc = 3
nms <- colnames(pt4)
nr <- ceiling(length(nms)/nc)
pal <- plasma(100, end = 0.95)
pal2 <- plasma(5, end = 0.95)

pdf("slingshot_umap_cluster4_res0.9_pseudotime_lineage.pdf",
    width = 12.5,
    height = 6)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt4[,i], breaks = 100)]
  names(pal2) <- levels(cut(pt4[,i], breaks = 5))
  plot(reducedDims(npbb.sce)$NPBB2.UMAP, col = colors, pch = 19, cex = 0.5, main = i,
       xlab = "UMAP-1", ylab = "UMAP-2", asp = 1)
  legend("topright", legend = names(pal2), col = pal2, pch = 16, bty = "n", cex = 2 / 3)
  lines(SlingshotDataSet(cur4), lwd = 2, col = 'black', type = 'curves')
}
dev.off()

save.image(object)


################## Stop here on 26.02.10


## Transfer lin1 pseudotime (PCA with start cluster4) to seurat object
pt1.df <- as.data.frame(pt1)
head(pt1.df)

colnames(pt1.df)<- c("Lineage1.pca", "Lineage2.pca", "Lineage3.pca")
npbb <- AddMetaData(npbb, pt1.df)


DimPlot(npbb, reduction = "npbb2.umap", group.by = "Lineage3.pca", label = FALSE, na.value = "grey50")+
  NoLegend()+ coord_fixed()  + guides(fill = guide_colourbar())+ 
  scale_color_viridis(discrete = TRUE, option="plasma", guide = "colourbar")
DimPlot(npbb, reduction = "npbb2.umap", group.by = "stage", label = FALSE) +
  coord_fixed()+
  scale_color_viridis(discrete = TRUE, option="viridis")


## Transfer lin4 pseudotime (UMAP-start cluster4) to seurat object
pt4.df <- as.data.frame(pt4)
head(pt4.df)
colnames(pt4.df)<- c("Lineage1.umap", "Lineage2.umap", "Lineage3.umap")
npbb <- AddMetaData(npbb, pt4.df)


DimPlot(npbb, reduction = "npbb2.umap", group.by = "Lineage1.umap", label = FALSE, na.value = "grey50")+
  NoLegend()+
  coord_fixed()  +
  guides(fill = guide_colourbar())+ 
  scale_color_viridis(discrete = TRUE, option="plasma", guide = "colourbar")

DimPlot(npbb, reduction = "npbb2.umap", group.by = "Lineage2.umap", label = FALSE, na.value = "grey50")+
  NoLegend()+ coord_fixed()  + guides(fill = guide_colourbar())+ 
  scale_color_viridis(discrete = TRUE, option="plasma", guide = "colourbar")

DimPlot(npbb, reduction = "npbb2.umap", group.by = "Lineage3.umap", label = FALSE, na.value = "grey50")+
  NoLegend()+ coord_fixed()  + guides(fill = guide_colourbar())+ 
  scale_color_viridis(discrete = TRUE, option="plasma", guide = "colourbar")

save.image(object)

##### Done script on 26.02.11