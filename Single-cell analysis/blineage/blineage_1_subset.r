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

wd1 <- "/scratch/gpfs/LEVINE/llemaire/singleCellAnalysis/26.01.26_iGmT/26.02.03_NS/26.02.08_blineage/26.02.10_subset_npbb"
setwd(wd1)
write(capture.output(sessionInfo()),file = "26.02.10sessionInfo.txt")

object <- "26.02.10blineage.RData"

npbb <- readRDS( file = "blineage.rds") 

human.homo <- read.csv("/scratch/gpfs/LEVINE/llemaire/cionaGeneModel/HT.KY21Gene.2.gff3/23.11.01CionaHomolog/gene-homolog2.csv", header = TRUE, row.names = 1)
human.homo <- as.data.table(human.homo)
human.homo <- human.homo %>% mutate(human.homolog=coalesce(human.homolog,gene))

tf <- read.csv("/scratch/gpfs/LEVINE/llemaire/cionaGeneModel/TF_KY21.csv", stringsAsFactors=FALSE, fileEncoding="latin1")

head(tf)



tf1 <- subset(tf, select = c("Gene.Name_KY", "Meaning_KY", "KY.gene.model"))


ZFtf <- read.csv("/scratch/gpfs/LEVINE/llemaire/cionaGeneModel/ZF-TF_KY21.csv",  
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

tf.name <- tf1$KY21
tf.name <- na.omit(tf.name)

rm(tf1, ZFtf1, ZFtf2, ZFtf3, ZFtf4)

DimPlot(npbb, reduction = "NS.umap", group.by = "seurat_clusters")


col1 <- c( "lightgrey","#330066")

### Nerve cord-b marker
FeaturePlot(npbb, features = c("KY21:KY21.Chr5.707", #Dmrt1 _a-lineage
                             "KY21:KY21.Chr6.345", #Flvcr2_ A-lineage
                             "KY21:KY21.Chr2.1031"), #Msx_ b-lineage
            reduction = "NS.umap", ncol = 3, cols = col1,
            min.cutoff = "q05", max.cutoff = "q95")+
  coord_fixed()


FeaturePlot(npbb, features = c("KY21:KY21.Chr14.584",#Cdx
                             "KY21:KY21.Chr4.737", # Ank1
                             "KY21:KY21.Chr12.134", # Lmna
                             "KY21:KY21.Chr8.1024", # Kcnip4
                             "KY21:KY21.Chr1.1257", # Lama2
                             "KY21:KY21.Chr8.1169", # Wnt7
                             "KY21:KY21.Chr9.606", # Lmx
                             "KY21:KY21.Chr10.167", # Cbl
                             "KY21:KY21.Chr7.611"), # SLit1 Marker for trunk nerve cord-b
                             ncol = 3,
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

#### Rubber cells
FeaturePlot(npbb, features = c("KY21:KY21.Chr9.822", #Wnttun5
                               "KY21:KY21.Chr4.1174", #Wnt5
                               "KY21:KY21.Chr8.449", #Rubber cells
                               "KY21:KY21.Chr7.170",#"KH2012:KH.L141.73", #muscle Rubber cells 
                               "KY21:KY21.Chr1.1701",#hand2 
                               "KY21:KY21.Chr7.709"), #Hox12
            reduction = "NS.umap", ncol=3,
            min.cutoff = "q10", max.cutoff = "q90")



save.image(object)

## Check cell number
samplename = npbb@meta.data$orig.ident
filteredCell<- as.data.frame(table(samplename))

save.image(object)
## Can keep sample separated because at least than 30 cells


##Clustering b.npb (start workflow again)

# split the dataset into a list of layer within the same object  

Idents(npbb) <- "orig.ident"

npbb[["RNA"]] <- split(npbb[["RNA"]], f = npbb$orig.ident)

# run standard anlaysis workflow

npbb <- FindVariableFeatures(npbb, selection.method = "vst", nfeatures = 4000)

save.image(object)
#This meanpbb that signals separating non-cycling cells and cycling cells will be maintained, 
#but differences in cell cycle phase among proliferating cells (which are often uninteresting), 
#will be regressed out of the data

######Not done, since most of the cells are still cycling
#npbb$CC.Difference <- npbb$S.Score - npbb$G2M.Score

npbb <- ScaleData(npbb, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(npbb))

npbb <- RunPCA(npbb, reduction.name = "npbb.pca",
             reduction.key = "npbbPC_")

DimPlot(npbb, reduction = "npbb.pca")
ggsave("npbb_PCA.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

#save.image(object)
#Check unintegrated data

npbb <- FindNeighbors(npbb, dims = 1:30, reduction = "npbb.pca")
npbb <- FindClusters(npbb, resolution = 0.6, cluster.name = "unintegrated.npbb_clusters")
npbb <- RunUMAP(npbb, dims = 1:30, reduction = "pca", reduction.name = "npbb.umap.unintegrated")
DimPlot(npbb, reduction = "npbb.umap.unintegrated", group.by = c("orig.ident", "unintegrated.npbb_clusters"))
ggsave("npbb.Unintegrated.UMAP.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")



#Perform integration ####### Decrease the number of k.wieght since only 36 cells in datasets

npbb <- IntegrateLayers(object = npbb, method = CCAIntegration, orig.reduction = "npbb.pca",
                        new.reduction = "npbb.integrated.cca", k.weight = 35, #dims = 1:30 ,k.anchor = NA,  , k.filter = NA,
                      verbose = TRUE)


# We can also specify parameters such as `k.anchor` to increase the strength of integration

npbb.2 <- IntegrateLayers(object = npbb, method = CCAIntegration, orig.reduction = "npbb.pca",
                        new.reduction = "npbb.integrated.cca.2", k.weight = 35, k.anchor = 10, #dims = 1:30 , , k.filter = NA,
                        verbose = TRUE)


###### dims = 1:30 (default)
## Error in FindIntegrationAnchors(object.list = object.list, anchor.features = features,  : 
##Max dimension too large: objects 10 contain fewer than 10 cells. 
##Please specify a maximum dimensions that is less than the number of cells in any object (10).
##### k.weight = 100
## k.weight (100) is set larger than the number of cells 
## in the smallest object (36). Please choose a smaller k.weight. 
#### 
## Issue with datasets smaller than 30 cells
## https://github.com/satijalab/seurat/issues/4803
## Suggestion: Maybe you can merge the small objects into the bigger one.

# re-join layers after integration

npbb[["RNA"]] <- JoinLayers(npbb[["RNA"]])

npbb.2[["RNA"]] <- JoinLayers(npbb.2[["RNA"]])
npbb.2

save.image(object)

#Choosing number of PCAs for clustering

ElbowPlot(npbb.2, ndims = 40, reduction = "npbb.pca") ###Cannot check integrated.cca
ggsave("npbb.Elbowplot.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


pdf("PCAHeatmap.pdf", width = 20, height = 20)
DimHeatmap(npbb.2, dims = 1:25, cells = 500, balanced = TRUE, reduction = "npbb.pca")
dev.off()

DimPlot(npbb.2, reduction ="npbb.pca", split.by="orig.ident", ncol = 3)

n.dims <- 12 ###Number of dimension for downstream analysis


npbb.default.integration <- npbb
npbb <- npbb.2

rm(npbb.2)



################################ K.anchor  = 10
npbb <- FindNeighbors(npbb, reduction = "npbb.integrated.cca.2", dims = 1:n.dims,
                    graph.name = c("npbb.RNA_nn","npbb.RNA_snn"))



resolutions <- seq(0.3, 0.9, 0.2)

npbb <- FindClusters(npbb,  algorithm = 4, resolution = resolutions, random.seed = 1,
                   graph.name = "npbb.RNA_snn") 

save.image(object)

pdf(file="cluster-tree.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(npbb, prefix = "npbb.RNA_snn_res.")
dev.off()

DimPlot(npbb, reduction ="pca", group.by="Phase", split.by = "stage")

#Project cells into UMAP and tSNE space

npbb <- RunTSNE(npbb, reduction = "npbb.integrated.cca.2", dims = 1:n.dims,
                reduction.name = "npbb.tsne", reduction.key = "npbbtsne_") 

npbb <- RunUMAP(npbb, reduction = "npbb.integrated.cca.2",
                dims = 1:n.dims, reduction.name = "npbb.umap",
                reduction.key = "npbbumap_")

DimPlot(npbb, reduction ="npbb.umap", group.by=c("npbb.RNA_snn_res.0.3",
                                                 "npbb.RNA_snn_res.0.5",
                                                 "npbb.RNA_snn_res.0.7",
                                                 "npbb.RNA_snn_res.0.9"))

ggsave("npbb.umap.res0.3-0.9.pdf", device= "pdf", width = 40, 
       height = 40, units = "cm")



DimPlot(npbb, reduction ="npbb.umap", group.by=c("stage",
                                                 "Phase"))
ggsave("npbb.umap.stage.phase.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")


DimPlot(npbb, reduction ="npbb.tsne", group.by=c("stage",
                                                 "Phase"))
ggsave("npbb.tsne.stage.phase.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")


Idents(npbb) <- "stage"


DimPlot(npbb, reduction = "npbb.umap", 
        cells.highlight = list(WhichCells(npbb, idents = "iG")),
        cols.highlight = c("darkblue"), sizes.highlight = 0.3,
        pt.size = 0.2) +coord_fixed()


DimPlot(npbb, reduction = "npbb.umap", 
        cells.highlight = list(WhichCells(npbb, idents = "iG"),
                               WhichCells(npbb, idents = "mG")),
        cols.highlight = c("darkred","darkblue"), sizes.highlight = 0.3,
        pt.size = 0.2) +coord_fixed() +
  scale_color_manual(labels = c("Other clusters",
                                "mG",
                                "iG"), values = c("grey","darkred","darkblue")) +
  #labs(color = "Analyzed Cell")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("npbb.umap.stage.iG-mG.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


DimPlot(npbb, reduction ="npbb.umap", group.by=c("cell.type.2",
                                                 "cell.type.3"))+ coord_fixed()
ggsave("npbb.umap.cell.type.2-3.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")



DimPlot(npbb, reduction ="npbb.tsne", group.by=c("npbb.RNA_snn_res.0.3",
                                                 "npbb.RNA_snn_res.0.5",
                                                 "npbb.RNA_snn_res.0.7",
                                                 "npbb.RNA_snn_res.0.9"))


########################### if default k.anchor

npbb <- RunUMAP(npbb, reduction = "npbb.integrated.cca",
                dims = 1:n.dims, reduction.name = "npbb.default.umap",
                reduction.key = "npbbDefaultumap_")

DimPlot(npbb, reduction ="npbb.default.umap", group.by=c("stage",
                                                 "Phase"))
ggsave("npbb.default.umap.stage.phase.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

###################################


#### Select res 0.5 as resolution

Idents(npbb) <- "npbb.RNA_snn_res.0.5" 

npbb$seurat_clusters <- npbb$npbb.RNA_snn_res.0.5


#### Check marker and find if should remove cluster 5

npbb.markers <- FindAllMarkers(npbb, only.pos = TRUE,
                             min.pct = 0.5,
                             logfc.threshold = 0.5,
                             test.use = "roc",
                             slot = "data")

save.image(object)

#annotate markers

npbb.markers <- tidyft::left_join(npbb.markers,human.homo,
                                by = "gene")
write.table(npbb.markers, file="npbb.DEG_combined_res0.5.txt",
            quote=F, sep="\t", col.names=NA)

save.image(file = object)

################## cluster5

DotPlot(npbb, features = c("KY21:KY21.Chr6.58",#Etr
                         "KY21:KY21.Chr11.637", #Tub
                         "KY21:KY21.Chr7.1172", # Eya
                         "KY21:KY21.Chr1.2012", #Tff1
                         "KY21:KY21.Chr1.275", #Tff2
                         "KY21:KY21.Chr14.584", # Cdx
                         "KY21:KY21.Chr7.872" #Uxs1
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/6)

#### Cell are likely epithelial, remove

#### Cluster1 - no marker

cluster1.markers <- FindMarkers(npbb, ident.1 = 1, only.pos = TRUE,
                               min.pct = 0.5,
                               logfc.threshold = 0.5,
                               test.use = "roc",
                               slot = "data")
cluster1.markers <- rownames_to_column(cluster1.markers)
colnames(cluster1.markers)[1] <- "gene"

cluster1.markers <- tidyft::left_join(cluster1.markers,human.homo,
                                  by = "gene")
write.table(cluster1.markers, file="cluster1.DEG_combined_res0.5.txt",
            quote=F, sep="\t", col.names=NA)

DotPlot(npbb, features = c("KY21:KY21.Chr6.58",#Etr
                           "KY21:KY21.Chr11.637", #Tub
                           "KY21:KY21.Chr2.1031", # Msx
                           "KY21:KY21.Chr4.947", #Tbx2/3
                           "KY21:KY21.Chr9.606", #Lmx1
                           "KY21:KY21.Chr14.584", # Cdx
                           "KY21:KY21.Chr7.872" #Uxs1
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/3)

FeaturePlot(npbb, features = c("KY21:KY21.Chr14.584",#Cdx
                               "KY21:KY21.Chr9.606", # Lmx
                               "KY21:KY21.Chr9.822" #Wnttun5
                               ), 
            ncol = 3,
            reduction = "npbb.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

ggsave("npbb.UMAP.Cdx_Lmx_Wnttun5.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")



###### assign cell types

clusters <- c(1,2,3,4,5)
cell.type.4 <- c("precusor cells",
                 "rubber cells",
                 "posterior b-neural tube cells",
                 "anterior b-neural tube cells",
                 "epithelial cells")

cell.type.4 <- data.frame(clusters, cell.type.4)

metadata<- npbb@meta.data

head(metadata)
head(cell.type)


colnames(cell.type.4)[1] <- "seurat_clusters"

cell.type.4$seurat_clusters <- factor(cell.type.4$seurat_clusters)
class(cell.type.4$seurat_clusters)

md1 <- subset(metadata,select = "seurat_clusters")
head(md1)
md1 <- rownames_to_column(md1)

md1 <- left_join(md1, cell.type.4, by = "seurat_clusters")
md1 <- column_to_rownames(md1)

npbb <- AddMetaData(npbb, md1)



############################
col2 <- c("#35608D", "#66CC33", "#E31A1C")
col1 <- c("lightgrey", "#330066")
col3 <- c("#35608D", "#66CC33")
nb.cols <- 14
col4 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 5
col5 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 7
col6 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

DimPlot(npbb, reduction = "npbb.tsne", 
        group.by = "orig.ident", cols = alpha(col4, 0.5))+
  coord_fixed()
ggsave("npbb.tSNE_orig.ident.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


DimPlot(npbb, reduction = "npbb.umap", 
        group.by = "orig.ident", cols = alpha(col4, 0.5), label.size = 2, pt.size = 0.5)+
  coord_fixed()
ggsave("npbb.UMAP_orig.ident.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(npbb, reduction = "npbb.umap", group.by = "npbb.RNA_snn_res.0.9",
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()
ggsave("npbb.UMAP_clusters_res0.9.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 


DimPlot(npbb, reduction = "npbb.umap", group.by = "cell.type.4",
        cols = alpha(col4, 0.5),
        pt.size = 0.5)+
  coord_fixed()
ggsave("npbb.UMAP_cell.type.4.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(npbb, reduction = "npbb.umap", group.by = "stage",
        cols = alpha(col6, 0.5),
        pt.size = 0.5)+
  coord_fixed()
ggsave("npbb.UMAP_stage.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 


#### Rubber cells
FeaturePlot(npbb, features = c("KY21:KY21.Chr9.822", #Wnttun5
                               "KY21:KY21.Chr4.1174", #Wnt5
                               "KY21:KY21.Chr8.449", #Rubber cells
                               "KY21:KY21.Chr7.170",#"KH2012:KH.L141.73", #muscle Rubber cells 
                               "KY21:KY21.Chr1.1701",#hand2 
                               "KY21:KY21.Chr7.709"), #Hox12
            reduction = "npbb.umap", ncol=3,
            min.cutoff = "q10", max.cutoff = "q90")

#### Nerve cord b marker
FeaturePlot(npbb, features = c("KY21:KY21.Chr14.584",#Cdx
                               "KY21:KY21.Chr4.737", # Ank1
                               "KY21:KY21.Chr12.134", # Lmna
                               "KY21:KY21.Chr8.1024", # Kcnip4
                               "KY21:KY21.Chr1.1257", # Lama2
                               "KY21:KY21.Chr8.1169", # Wnt7
                               "KY21:KY21.Chr9.606", # Lmx
                               "KY21:KY21.Chr10.167", # Cbl
                               "KY21:KY21.Chr7.611"), # SLit1 Marker for trunk nerve cord-b
            ncol = 3,
            reduction = "npbb.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

save.image(file = object)

###### Add tissue type to npbb

DimPlot(npbb, reduction = "npbb.umap", group.by = "tissue.type.1",
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()

DimPlot(npbb, reduction = "npbb.umap", group.by = "cell.type.1",
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()

DimPlot(npbb, reduction = "npbb.umap", group.by = "cell.type.2",
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()

DimPlot(npbb, reduction = "npbb.umap", group.by = "cell.type.3",
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()


DimPlot(npbb, reduction = "npbb.umap", group.by = "lineage",
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()

clusters <- c(1,2,3,4,5)
tissue.type.2 <- c("b lineage",
                 "b lineage",
                 "b lineage",
                 "b lineage",
                 "epithelial cells")

tissue.type.2<- data.frame(clusters, tissue.type.2)
colnames(tissue.type.2)[1] <- "seurat_clusters"

tissue.type.2$seurat_clusters<- as.factor(tissue.type.2$seurat_clusters)
metadata <- npbb@meta.data

md2 <- subset(metadata,select = "seurat_clusters")
head(md2)
md2 <- rownames_to_column(md2)

md2 <- left_join(md2, tissue.type.2, by = "seurat_clusters")
md2 <- column_to_rownames(md2)

npbb <- AddMetaData(npbb, md2)

DimPlot(npbb, reduction = "npbb.umap", group.by = "tissue.type.2",
        cols = alpha(c(col4[1], col4[4]), 0.5),
        pt.size = 0.5)+
  coord_fixed()
ggsave("npbb.UMAP_tissue.type.2.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

FeaturePlot(npbb, features = "KY21:KY21.Chr7.872", #Uxs1
cols = col1,
reduction = "npbb.umap", 
min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
ggsave("npbb.UMAP_uxs1.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

save.image("26.02.10blineage_npbb.RData")

####### Remove epithelial cells

saveRDS(npbb, file = "blineage.clean.rds")

npbb <- readRDS( file = "blineage.clean.rds") 

head(Idents(npbb))

DimPlot(npbb, reduction = "npbb.umap", 
        cols = alpha(col4, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()

npbb.2 <- subset(x = npbb, idents = '5', invert = TRUE)

DimPlot(npbb.2, reduction = "npbb.umap", group.by = "cell.type.4",
        cols = alpha(col4, 0.5),
        pt.size = 0.5)+
  coord_fixed()
ggsave("npbb.2.UMAP_cell.type.4.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 


FeaturePlot(npbb.2, features = c("KY21:KY21.Chr14.584",#Cdx
                               "KY21:KY21.Chr9.606", # Lmx
                               "KY21:KY21.Chr9.822" #Wnttun5
                                ), 
ncol = 3, cols = col1,
reduction = "npbb.umap", 
min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

ggsave("npbb.UMAP.Cdx_Lmx_Wnttun5.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

rm(npbb)
### Reproject the data


## Check cell number
samplename = npbb.2@meta.data$orig.ident
filteredCell<- as.data.frame(table(samplename))

save.image(object)
## Can keep sample separated because at least than 30 cells


##Clustering b.npb (start workflow again)

# split the dataset into a list of layer within the same object  

Idents(npbb.2) <- "orig.ident"

npbb.2[["RNA"]] <- split(npbb.2[["RNA"]], f = npbb.2$orig.ident)

# run standard anlaysis workflow

npbb.2 <- FindVariableFeatures(npbb.2, selection.method = "vst", nfeatures = 4000)

save.image(object)
#This meanpbb that signals separating non-cycling cells and cycling cells will be maintained, 
#but differences in cell cycle phase among proliferating cells (which are often uninteresting), 
#will be regressed out of the data

######Not done, since most of the cells are still cycling
#npbb$CC.Difference <- npbb$S.Score - npbb$G2M.Score

npbb.2 <- ScaleData(npbb.2, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(npbb.2))

npbb.2 <- RunPCA(npbb.2, reduction.name = "npbb2.pca",
               reduction.key = "npbb2PC_")

DimPlot(npbb.2, reduction = "npbb2.pca", #dims = c(3, 4),
        cols = col4)
ggsave("npbb2_PCA.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

save.image(object)
#Check unintegrated data

#Perform integration ####### Decrease the number of k.wieght since only 36 cells in datasets

npbb.2b <- IntegrateLayers(object = npbb.2, method = CCAIntegration, orig.reduction = "npbb2.pca",
                        new.reduction = "npbb2.integrated.cca", k.weight = 34, #dims = 1:30 ,,  , k.filter = NA,
                        k.anchor = 10, verbose = TRUE)

npbb.2c <- IntegrateLayers(object = npbb.2, method = CCAIntegration, orig.reduction = "npbb2.pca",
                           new.reduction = "npbb2.integrated.cca", k.weight = 34, #dims = 1:30 ,,  , k.filter = NA,
                           k.anchor = 10, verbose = TRUE)

npbb.2d <- IntegrateLayers(object = npbb.2, method = RPCAIntegration, orig.reduction = "npbb2.pca",
                           new.reduction = "npbb2.integrated.cca", k.weight = 30, #dims = 1:30 ,,  , k.filter = NA,
                           k.anchor = 3, verbose = TRUE)

# We can also specify parameters such as `k.anchor` to increase/decrease the strength of integration, default is k.anchor = 5
#Try 3 but still a mass of cells
#Try 10, blops
##############Try with RPCA and low k.anchor, it works separate the different cell types.


###### dims = 1:30 (default)
## Error in FindIntegrationAnchors(object.list = object.list, anchor.features = features,  : 
##Max dimension too large: objects 10 contain fewer than 10 cells. 
##Please specify a maximum dimensions that is less than the number of cells in any object (10).
##### k.weight = 100
## k.weight (100) is set larger than the number of cells 
## in the smallest object (36). Please choose a smaller k.weight. 
#### 
## Issue with datasets smaller than 30 cells
## https://github.com/satijalab/seurat/issues/4803
## Suggestion: Maybe you can merge the small objects into the bigger one.

# re-join layers after integration

npbb.2b[["RNA"]] <- JoinLayers(npbb.2b[["RNA"]])
npbb.2c[["RNA"]] <- JoinLayers(npbb.2c[["RNA"]])
npbb.2d[["RNA"]] <- JoinLayers(npbb.2d[["RNA"]])

#save.image(object)





################################
#Project cells into UMAP 

#Choosing number of PCAs for projection and clustering

ElbowPlot(npbb.2d, ndims = 40, reduction = "npbb2.pca") ###Cannot check integrated.cca
ggsave("npbb2.Elbowplot.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


#pdf("PCAHeatmap.pdf", width = 20, height = 20)
DimHeatmap(npbb.2d, dims = 1:12, cells = 200, balanced = TRUE, reduction = "npbb2.pca")
#dev.off()

DimPlot(npbb.2d, reduction ="npbb2.pca", split.by="orig.ident", ncol = 3)
n.dims <- 10 ###Number of dimension for downstream analysis

#npbb.2b <- RunTSNE(npbb.2b, reduction = "npbb2.integrated.cca", dims = 1:n.dims,
#                  reduction.name = "npbb2.tsne", reduction.key = "npbb2tsne_") 

npbb.2b <- RunUMAP(npbb.2b, reduction = "npbb2.integrated.cca",
                  dims = 1:n.dims, reduction.name = "npbb2.umap",
                  reduction.key = "npbb2umap_")

npbb.2c <- RunUMAP(npbb.2c, reduction = "npbb2.integrated.cca",
                   dims = 1:n.dims, reduction.name = "npbb2.umap",
                   reduction.key = "npbb2umap_")

npbb.2d <- RunUMAP(npbb.2d, reduction = "npbb2.integrated.cca",
                   dims = 1:n.dims, reduction.name = "npbb2.umap",
                   reduction.key = "npbb2umap_")

DimPlot(npbb.2d, reduction ="npbb2.umap",
        group.by="cell.type.4", cols = col4)+
  coord_fixed()
ggsave("npbb2_UMAP_cell.type.4.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

DimPlot(npbb.2d, reduction ="npbb2.umap",
        group.by="stage")+
  scale_color_viridis(discrete = TRUE, option= "viridis", begin = 0.2,
                      end = 1, direction = -1)+
  coord_fixed()
ggsave("npbb2_UMAP_stage.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

#npbb.2d looks great, change to npbb.2

npbb.2 <- npbb.2d


######################## Find clusters

npbb.2 <- FindNeighbors(npbb.2, reduction = "npbb2.integrated.cca", dims = 1:n.dims,
                      graph.name = c("npbb2.RNA_nn","npbb2.RNA_snn"))

resolutions <- seq(0.3, 1.2, 0.3)

npbb.2 <- FindClusters(npbb.2,  algorithm = 4, resolution = resolutions, random.seed = 1,
                     graph.name = "npbb2.RNA_snn") 

save.image(object)

pdf(file="cluster-tree_npbb2.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(npbb.2, prefix = "npbb2.RNA_snn_res.")
dev.off()

DimPlot(npbb.2, reduction ="npbb2.pca", group.by="Phase", split.by = "stage")

#Project clusters to UMAP

DimPlot(npbb.2, reduction ="npbb2.umap", group.by=c("npbb2.RNA_snn_res.0.3",
                                                 "npbb2.RNA_snn_res.0.6",
                                                 "npbb2.RNA_snn_res.0.9",
                                                 "npbb2.RNA_snn_res.1.2"),
        cols = col4, label = TRUE, ncol = 4)+coord_fixed()
ggsave("npbb2.umap.res0.3-1.2.pdf", device= "pdf", width = 80, 
       height = 20, units = "cm")


Idents(npbb.2) <- "stage"

DimPlot(npbb.2, reduction = "npbb2.umap", 
        cells.highlight = list(WhichCells(npbb.2, idents = "iG")),
        cols.highlight = c("darkblue"), sizes.highlight = 0.3,
        pt.size = 0.2) +coord_fixed()


DimPlot(npbb.2, reduction = "npbb2.umap", 
        cells.highlight = list(WhichCells(npbb.2, idents = "iG"),
                               WhichCells(npbb.2, idents = "mG")),
        cols.highlight = c("darkred","darkblue"), sizes.highlight = 0.3,
        pt.size = 0.2) +coord_fixed() +
  scale_color_manual(labels = c("Other clusters",
                                "mG",
                                "iG"), values = c("grey","darkred","darkblue")) +
  #labs(color = "Analyzed Cell")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("npbb2.umap.stage.iG-mG.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


FeaturePlot(npbb.2, features = c("KY21:KY21.Chr14.584",#Cdx
                               "KY21:KY21.Chr9.606", # Lmx
                               "KY21:KY21.Chr9.822" #Wnttun5
), 
ncol = 3,
reduction = "npbb2.umap", cols = col1,
min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

ggsave("npbb2.UMAP.Cdx_Lmx_Wnttun5.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

FeaturePlot(npbb.2, reduction = "npbb2.umap", features = c("KY21:KY21.Chr9.606", #Lmx
                                                         "KY21:KY21.Chr10.215", #Hep-r.a
                                                         "KY21:KY21.Chr2.1031", #Msx
                                                         "KY21:KY21.Chr14.584"), #Cdx
            cols = col1, ncol = 4, min.cutoff = "q5",
            max.cutoff = "q95")+ coord_fixed()
ggsave("npbb2_Lmx_Hep-r.a_Msx_Cdx.pdf", device= "pdf", width = 80, 
       height = 20, units = "cm")


########### Add cell types (Also add in visualization 6 script)

npbb.2$seurat_clusters <- npbb.2$npbb2.RNA_snn_res.0.9
DimPlot(npbb.2, reduction = "npbb2.umap", group.by = "seurat_clusters")


Idents(npbb.2) <- "seurat_clusters"

clusters <- c(1,2,3,4,5,6)
cell.type.5 <- c("anterior b-neural tube cells",
                 "posterior b-neural tube cells",
                  "rubber cells",
                 "precusor cells",
                  "rubber cells",
                 "anterior b-neural tube cells")

cell.type.5 <- data.frame(clusters, cell.type.5)

metadata<- npbb.2@meta.data


colnames(cell.type.5)[1] <- "seurat_clusters"

cell.type.5$seurat_clusters <- factor(cell.type.5$seurat_clusters)
class(cell.type.5$seurat_clusters)

md1 <- subset(metadata,select = "seurat_clusters")
head(md1)
md1 <- rownames_to_column(md1)

md1 <- left_join(md1, cell.type.5, by = "seurat_clusters")
md1 <- column_to_rownames(md1)

npbb.2 <- AddMetaData(npbb.2, md1)

DimPlot(npbb.2, reduction = "npbb2.umap", group.by = "cell.type.5",
        cols = alpha(col4, 0.5),
        pt.size = 0.5)+
  coord_fixed()
ggsave("npbb.2.UMAP_cell.type.5.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 




saveRDS(npbb.2, file = "blineage2.rds")




save.image(object)

