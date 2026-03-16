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

wd1 <- "/YourDirectory"
setwd(wd1)


object <- "blineage_4.RData" #(Cluster gene expression pattern with 40 PCAs)
load(object)
object <- "blineage_5.RData"


clusterLabels3 <- primaryCluster(clusPat3$rsec) #Correspond to mergeClusters, first column of clusterMatrix
cUniq3 <- unique(clusterLabels3)
cUniq3a <- cUniq3[!cUniq3 == -1] # remove unclustered genes

names(clusterLabels3 )<- clusPat3[["rsec"]]@colData@rownames

head(clusterLabels3)
clusterLabels3 <- as.data.frame(clusterLabels3)
clusterLabels3 <- rownames_to_column(clusterLabels3)
colnames(clusterLabels3) [1] <- "gene"
clusterLabels3 <- tidyft::left_join(clusterLabels3,human.homo,
                                    by = "gene")
write.table(clusterLabels3, file="blineage_clusterLabels_lineage3.txt",
            quote=F, sep="\t", col.names=NA)

viridis(8, option = "viridis")
show_col(viridis(8, option = "viridis"))

clusterLabels3b <- primaryCluster(clusPat3$rsec)

q <- list()
for (xx in cUniq3a) {
  cId <- which(clusterLabels3b == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat3$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat3$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
                                  y = clusPat3$yhatScaled[geneId, ],
                                  lineage = rep(0:2, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = "none") +
    scale_color_manual(values = c("#FDE725FF", "#277f8EFF", "#440154FF"),
                      breaks = c("0", "1", "2"))  
  print(p)
  q[[xx]] <- p
}



pdf("geneCluster_lineage3.pdf",
    width = 25,
    height = 25)
do.call(grid.arrange,q)
dev.off()

r <- list()
for (xx in cUniq3a) {
  cId <- which(clusterLabels3b == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat3$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat3$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 3),
                                  y = clusPat3$yhatScaled[geneId, ],
                                  lineage = rep(0:2, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = "none") +
    scale_color_manual(values = c("#FFFFFF00", "#FFFFFF00", "#440154FF"),
                       breaks = c("0", "1", "2"))  
  print(p)
  r[[xx]] <- p
}

pdf("geneCluster_lineage3_only.pdf",
    width = 25,
    height = 25)
do.call(grid.arrange,r)
dev.off()


metrics_pt_lineage3 <- metrics_pt [ which(metrics_pt$gene %in% deg.lineage3), ]

metrics_pt_lineage3 <- tidyft::left_join(clusterLabels3, metrics_pt_lineage3,
                                         by = "gene")
write.table(metrics_pt_lineage3, file="metrics_pt_lineage3.txt",
            quote=F, sep="\t", col.names=NA)

save.image(object)

##Gene expression cascade of TF

# Find tf

tf <- read.csv("/YourDirectory/cionaGeneModel/TF_KY21.csv", stringsAsFactors=FALSE, fileEncoding="latin1")

head(tf)

tf1 <- subset(tf, select = c("Gene.Name_KY", "Meaning_KY", "KY.gene.model"))


ZFtf <- read.csv("/YourDirectory/cionaGeneModel/ZF-TF_KY21.csv",  
                 stringsAsFactors=FALSE, fileEncoding="latin1",
                 header = TRUE)

head(ZFtf)

ZFtf1 <- subset(ZFtf, select = c("Gene.Name", "Meaning_KY", "KY.gene.model"))
colnames(ZFtf1) [1] <- "Gene.Name_KY"

ZFtf2 <- subset(ZFtf1, ZFtf1$KY.gene.model %in% tf1$KY.gene.model)

ZFtf3 <- dplyr::anti_join(ZFtf1, ZFtf2, by = "KY.gene.model")

B<-c(1, 2, 251, 253, 254, 255, 256, 257)

ZFtf4<-ZFtf3[B,]


ZFtf4$KY21 <- paste("KY21:", ZFtf4$KY.gene.model, sep = "")
tf1$KY21 <- paste("KY21:", tf1$KY.gene.model, sep = "")

tf1 <- rbind(tf1, ZFtf4)

head(tf1)
tf.name <- tf1$KY21
tf.name <- na.omit(tf.name)

rm(tf1, ZFtf1, ZFtf2, ZFtf3, ZFtf4)

tf.deg.lineage3 <- deg.lineage3[which(deg.lineage3 %in% tf.name)]

tf.deg.lineage3 <- as.data.frame(tf.deg.lineage3)
colnames(tf.deg.lineage3) [1] <- "gene"
tf.deg.lineage3 <- tidyft::left_join(tf.deg.lineage3,human.homo,
                                    by = "gene")

### instead of dynamic TF, take all TF above a threshold
#### npbb.sce.2 has the fitGAm
#### pt4 and pt4.df has the lineage information

lineage3.cells <- pt4.df[!is.na(pt4.df$Lineage3.umap),]
lineage3.cells <- rownames_to_column(lineage3.cells)
lineage3.cells <- lineage3.cells$rowname

npbb.sce.lineage3 <- npbb.sce[, lineage3.cells]

geneFilter.lineage3 <- apply(assays(npbb.sce.lineage3)$counts,1,function(x){
  sum(x >= 4) >= 12})

head(geneFilter.lineage3)

geneFilter.lineage3 <- as.data.frame(geneFilter.lineage3)
table(geneFilter.lineage3)
geneFilter.lineage3 <- rownames_to_column(geneFilter.lineage3)
colnames(geneFilter.lineage3)[1] <- "gene"

geneFilter.lineage3 <-  geneFilter.lineage3[ which(geneFilter.lineage3$geneFilter.lineage3 == TRUE), ]
geneFilter.lineage3 <- geneFilter.lineage3$gene

tf.express <- subset(tf.name, tf.name %in% geneFilter.lineage3)

tf.express <- as.data.frame(tf.express)

colnames(tf.express)[1] <- "gene"

tf.express <- tidyft::left_join(tf.express,human.homo,
                                by = "gene")

### Smad2/3a 	KY21.Chr6.524
### Smad1/5/9 Chr2.623
### Zicl Chr6.26 - Chr6.31
### Lmx Chr9.606
### Pax3/7 Chr10.288
### Msx Chr2.1031




### based on mean smoother

yhatSmooth <- predictSmooth(npbb.sce.2, gene = tf.express$gene, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 101:150]))),
                       fontsize = 6,
                       cluster_cols = FALSE,
                       show_rownames = TRUE, 
                       show_colnames = FALSE,
                       cluster_rows = TRUE) 
pdf("GRC_predictsmooth.pdf",
    width = 8,
    height = 24)
heatSmooth
dev.off()

#Rename TF with common name

tf.express.name  <- structure(as.character(tf.express$ciona.name),
                                   names = as.character(tf.express$gene))

col8 <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(100)
pdf("GRC_predictsmooth_2.pdf",
    width = 50,
    height = 25)
Heatmap(t(scale(t(yhatSmooth[, 101:150]))),
        cluster_rows=TRUE, cluster_columns = FALSE,
        show_column_names = FALSE,
        col = col8,
        row_names_gp = gpar(fontsize = 8),
        width = ncol(yhatSmooth)*unit(0.4, "mm"),
        height = nrow(yhatSmooth)*unit(3, "mm"),
        row_labels = tf.express.name[rownames(yhatSmooth)])
dev.off()


#Order gene based on first appear
window <- 2
step <- 1

yhatSmooth2 <- yhatSmooth[, 101:150]

yhatSmooth2  <- yhatSmooth2 [order(apply(t(rollapply(t(yhatSmooth2 ),
                                                     width=window,
                                                     by=step,
                                                     FUN=mean)), 1, which.max)), ]

order3 <- rownames(yhatSmooth2 )

tf.express.order <- tf.express[match(order3, tf.express$gene), ]

write.table(tf.express.order, file="tf_expressed-order_lineage3.txt",
            quote=F, sep="\t", col.names=NA)

pdf("GRC_predictsmooth_2_order_KY21.pdf",
    width = 10,
    height = 20)
Heatmap(t(scale(t(yhatSmooth[, 101:150]))),
        cluster_rows=FALSE, row_order = order3,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        col = col8,
        # row_labels = tf.express.name[rownames(yhatSmooth)],
        row_names_gp = gpar(fontsize = 8),
     
        width = ncol(yhatSmooth)*unit(0.4, "mm"),
        height = nrow(yhatSmooth)*unit(3, "mm"))
dev.off()

pdf("GRC_predictsmooth_2_order_name.pdf",
    width = 10,
    height = 20)
Heatmap(t(scale(t(yhatSmooth[, 101:150]))),
        cluster_rows=FALSE, row_order = order3,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        col = col8,
        row_labels = tf.express.name[rownames(yhatSmooth)],
        row_names_gp = gpar(fontsize = 6),
        width = ncol(yhatSmooth)*unit(0.4, "mm"),
        height = nrow(yhatSmooth)*unit(2, "mm"))
dev.off()


############## Decrease scale range
col8 <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(101)
show_col(col8)

col9 <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                          "RdYlBu")))(31)
show_col(col9)
col10 =  viridis(7, option = "viridis")
barplot(1:13, col=col10)

col_fun = colorRamp2(seq(-2, 2, 0.04), c(col8))

col_fun(seq(-3, 3))

pdf("GRC_predictsmooth_scaled_order_KY21.pdf",
    width = 10,
    height = 20)
Heatmap(t(scale(t(yhatSmooth[, 101:150]))),
        cluster_rows=FALSE, row_order = order3,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        col = col_fun,
        # row_labels = tf.express.name[rownames(yhatSmooth)],
        row_names_gp = gpar(fontsize = 6),
        width = ncol(yhatSmooth)*unit(0.3, "mm"),
        height = nrow(yhatSmooth)*unit(1.8, "mm"))
dev.off()

pdf("GRC_predictsmooth_scaled_order_name.pdf",
    width = 10,
    height = 20)
Heatmap(t(scale(t(yhatSmooth[, 101:150]))),
        cluster_rows=FALSE, row_order = order3,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        col = col_fun,
        row_labels = tf.express.name[rownames(yhatSmooth)],
        row_names_gp = gpar(fontsize = 6),
        width = ncol(yhatSmooth)*unit(0.3, "mm"),
        height = nrow(yhatSmooth)*unit(1.8, "mm"))
dev.off()

### If no Z score
Heatmap(yhatSmooth[, 101:150],
        cluster_rows=TRUE, 
        cluster_columns = FALSE,
        show_column_names = FALSE,
        col = col8,
        row_labels = tf.express.name[rownames(yhatSmooth)],
        row_names_gp = gpar(fontsize = 5),
        width = ncol(yhatSmooth)*unit(0.4, "mm"),
        height = nrow(yhatSmooth)*unit(1.5, "mm"))

### based on fitted values (plotting takes a while to run)
yhatCell <- predictCells(npbb.sce.2, gene = tf.express$gene)

# order according to pseudotime only
Lineage3 <- pt4.df[complete.cases(pt4.df[ , 3]),]
Lineage3 <- subset(Lineage3, select = "Lineage3.umap")
colnames(Lineage3)[1] <- "pseudotime"
Lineage3 <- rownames_to_column(Lineage3)
oLineage3<- Lineage3[order(Lineage3$pseudotime),]
oLineage3cell <- oLineage3$rowname

yhatCellLineage3 <- yhatCell[,oLineage3cell]

pdf("GRC_predictcell.pdf",
    width = 8,
    height = 8)
pheatmap(t(scale(t(yhatCellLineage3))), cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames=FALSE)
dev.off() ##### No smoothing therefore uneasy to read

save.image(object)


####Add pseudotime lineages to seurat object

metadata <- npbb@meta.data

for (x in 1:nrow(metadata))
{
  metadata[x, "Lineage3b"] <- if (is.na(metadata[x, "Lineage3.umap"]))
  {FALSE} 
  else 
  {TRUE}
}

for (x in 1:nrow(metadata))
{
  metadata[x, "Lineage2b"] <- if (is.na(metadata[x, "Lineage2.umap"]))
  {FALSE} 
  else 
  {TRUE}
}

for (x in 1:nrow(metadata))
{
  metadata[x, "Lineage1b"] <- if (is.na(metadata[x, "Lineage1.umap"]))
  {FALSE} 
  else 
  {TRUE}
}

for (x in 1:nrow(metadata))
{
  metadata[x, "Lineage.pseudotime"] <- if (metadata[x, "Lineage1b"] == TRUE 
                                           & metadata[x, "Lineage2b"] == TRUE
                                           & metadata[x, "Lineage3b"] == TRUE)
  {"common lineage"} 
  else
  {if (metadata[x, "Lineage1b"] == TRUE 
       & metadata[x, "Lineage2b"] == TRUE)
  {"lineage 1-2"}
    else
    {if (metadata[x, "Lineage2b"] == TRUE 
         & metadata[x, "Lineage3b"] == TRUE)
    {"lineage 2-3"}
      else
      {if (metadata[x, "Lineage1b"] == TRUE 
           & metadata[x, "Lineage3b"] == TRUE)
      {"lineage 1-3"}
        else
        {if (metadata[x, "Lineage1b"] == TRUE)
  {"lineage 1"}
    else
          {if (metadata[x, "Lineage2b"] == TRUE)
          {"lineage 2"}
          else  
    {"lineage 3"}}
} } } }}

md3 <- subset(metadata, select= c("Lineage3b",
                                  "Lineage2b",
                                  "Lineage1b",
                                  "Lineage.pseudotime"))
npbb <- AddMetaData(npbb, md3)

DimPlot(npbb, reduction = "npbb2.umap", group.by = "Lineage.pseudotime",
        cols = col10)+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1)+
  labs( title = "Lineage")
ggsave("Umap_lineage.slingshot.pdf", device= "pdf", width = 14, 
       height = 10, units = "cm")

save.image(object)



Idents(npbb) <- "Lineage.pseudotime"

class(npbb$Lineage.pseudotime)

FeaturePlot(npbb, reduction = "npbb2.umap", features = "Lineage3.umap")+
  scale_color_viridis(discrete = FALSE, option="plasma", guide = "colourbar", na.value="lightgrey")+ 
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1)+
  labs( title = "Pseudotime")
ggsave("npbb_Umap_lineage3_pseudotime.pdf", device= "pdf", width = 14, 
       height = 10, units = "cm")
dev.off()

Idents(npbb) <- "Lineage3b"
npbb.lin3 <- subset(npbb, idents = "TRUE") 

Lineage3Data.1 <- FeaturePlot(npbb.lin3, reduction = "npbb2.umap", features = "Lineage3.umap")+
  scale_color_viridis(discrete = FALSE, option="plasma", guide = "colourbar") 

Lineage3Data.2  <- ggplot_build(Lineage3Data.1)

head(Lineage3Data.2$data[[1]])

pdf(file = "Umap_lineage3_trajectory.pdf", width = 6, 
    height = 6)
par(pty = "s") #plot area will be square
plot(x = Lineage3Data.2[["data"]][[1]][["x"]], 
     y = Lineage3Data.2[["data"]][[1]][["y"]],
     col = Lineage3Data.2[["data"]][[1]][["colour"]],
     pch = 19, xlab = "Umap-1", ylab = "Umap-2",
     xlim = c(-5,5),
     ylim = c(-6,6))
lines(SlingshotDataSet(cur4), type = 'curves', lwd=2, col='black')
dev.off()

save.image(object)

### Extract trajectory
 
cur4b <- slingCurves((SlingshotDataSet(cur4))) #object containing the coordinates for the 3 lineage trajectories

cur4b_lineage1 <- data.frame(cur4b[[1]]$s[cur4b[[1]]$ord, ]) # ordered coordinates from the lineage1 trajectory

cur4b_lineage2 <-  data.frame(cur4b[[2]]$s[cur4b[[2]]$ord, ]) # ordered coordinates from the lineage2 trajectory

cur4b_lineage3 <- data.frame(cur4b[[3]]$s[cur4b[[3]]$ord, ]) # ordered coordinates from the lineage3 trajectory
  

ggplot(cur4b_lineage1, aes( x = npbb2umap_1, y = npbb2umap_2)) +
  geom_path()


FeaturePlot(npbb, reduction = "npbb2.umap", features = "Lineage3.umap")+coord_fixed()+
  scale_color_viridis(discrete = FALSE, option="plasma", guide = "colourbar", na.value="#D3D3D380")+ 
    geom_path(data = cur4b_lineage1, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage2, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage3, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
 theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1)+
  labs( title = "Lineage 3 Pseudotime")
ggsave("npbb_Umap_lineage3_pseudotime-trajectories.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


FeaturePlot(npbb, reduction = "npbb2.umap", features = "Lineage1.umap")+coord_fixed()+
  scale_color_viridis(discrete = FALSE, option="plasma", guide = "colourbar", na.value="#D3D3D380")+ 
  geom_path(data = cur4b_lineage1, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage2, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage3, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1)+
  labs( title = "Lineage 1 Pseudotime")
ggsave("npbb_Umap_lineage1_pseudotime-trajectories.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

FeaturePlot(npbb, reduction = "npbb2.umap", features = "Lineage2.umap")+coord_fixed()+
  scale_color_viridis(discrete = FALSE, option="plasma", guide = "colourbar", na.value="#D3D3D380")+ 
  geom_path(data = cur4b_lineage1, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage2, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage3, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1)+
  labs( title = "Lineage 2 Pseudotime")
ggsave("npbb_Umap_lineage2_pseudotime-trajectories.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

FeaturePlot(npbb.lin3, reduction = "npbb2.umap", features = "Lineage3.umap")+
  scale_color_viridis(discrete = FALSE, option="plasma", guide = "colourbar", na.value="lightgrey")+ 
  #geom_path(data = cur4b_lineage1, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  #geom_path(data = cur4b_lineage2, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  geom_path(data = cur4b_lineage3, mapping = aes( x = npbb2umap_1, y = npbb2umap_2), linewidth = 1)+
  xlim(c(-3, 5))+
  ylim(c(-7, 5))+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1)+
  labs( title = "Lineage 3 Pseudotime")
ggsave("npbb-lineage3_only_Umap_pseudotime-trajectory.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

save.image(object)



save.image(object)



