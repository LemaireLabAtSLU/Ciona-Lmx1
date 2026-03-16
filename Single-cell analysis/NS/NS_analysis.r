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

wd1 <- "/scratch/gpfs/LEVINE/llemaire/singleCellAnalysis/26.01.26_iGmT/26.02.03_NS"
setwd(wd1)

write(capture.output(sessionInfo()),file = "26.02.03_sessionInfo.txt")


object <- "26.02.03.iGmT.NS.RData"

NS <- readRDS( file = "NS.rds")       ### include midline which will give rise to epidermal neurons 

DimPlot(NS, reduction = "UMAP", group.by = "seurat_clusters")
DimPlot(NS, reduction = "UMAP", group.by = "cell.type.1") + coord_fixed()

DimPlot(NS, reduction = "UMAP", split.by = "stage", group.by = "cell.type.1") + coord_fixed()

DimPlot(NS, reduction = "UMAP", split.by = "stage", group.by = "Phase") + coord_fixed()

#### Scale the data with regression over cell cycle

# split the dataset into a list of layer within the same object  

Idents(NS) <- "orig.ident"

NS[["RNA"]] <- split(NS[["RNA"]], f = NS$orig.ident)

# run standard anlaysis workflow

NS <- FindVariableFeatures(NS, selection.method = "vst", nfeatures = 4000)

save.image(object)

#This means that signals separating non-cycling cells and cycling cells will be maintained, 
#but differences in cell cycle phase among proliferating cells (which are often uninteresting), 
#will be regressed out of the data
#NS$CC.Difference <- NS$S.Score - NS$G2M.Score ######Not done, since most of the cells are still cycling

NS <- ScaleData(NS, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(NS))

NS <- RunPCA(NS, reduction.name = "NS.pca",
             reduction.key = "NSPC_")

DimPlot(NS, reduction = "NS.pca")
ggsave("NS_PCA.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

save.image(object)

### Cluster and project without integration

NS <- FindNeighbors(NS, dims = 1:30, reduction = "NS.pca")
NS <- FindClusters(NS, resolution = 0.6, cluster.name = "unintegrated.NS_clusters")
NS <- RunUMAP(NS, dims = 1:30, reduction = "pca", reduction.name = "NS.umap.unintegrated")
DimPlot(NS, reduction = "NS.umap.unintegrated", group.by = c("orig.ident", "unintegrated.NS_clusters"))
ggsave("NS.Unintegrated.UMAP.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")


#Perform integration

NS <- IntegrateLayers(object = NS, method = CCAIntegration, orig.reduction = "NS.pca", new.reduction = "NS.integrated.cca",
                           verbose = TRUE)

# re-join layers after integration

NS[["RNA"]] <- JoinLayers(NS[["RNA"]])

save.image(object)

#Choosing number of PCAs for clustering

ElbowPlot(NS, ndims = 40, reduction = "NS.pca") ###Cannot check integrated.cca
ggsave("NS.Elbowplot.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


pdf("PCAHeatmap.pdf", width = 20, height = 20)
DimHeatmap(NS, dims = 1:25, cells = 500, balanced = TRUE, reduction = "NS.pca")
dev.off()

DimPlot(NS, reduction ="pca", split.by="orig.ident", ncol = 3)

n.dims <- 20 ###Number of dimension for downstream analysis

NS <- FindNeighbors(NS, reduction = "NS.integrated.cca", dims = 1:n.dims,
                    graph.name = c("NS.RNA_nn","NS.RNA_snn"))

resolutions <- seq(0.8, 4.0, 0.4)

NS <- FindClusters(NS,  algorithm = 4, resolution = resolutions,
                   graph.name = "NS.RNA_snn", random.seed = 1) 

save.image(object)



pdf(file="cluster-tree.pdf", 
    width = 16, # The width of the plot in inches
    height = 12) # The height of the plot in inches
clustree(NS, prefix = "NS.RNA_snn_res.")
dev.off()

DimPlot(NS, reduction ="pca", group.by="Phase", split.by = "stage")

#### Select res 2.4 as resolution

Idents(NS) <- "NS.RNA_snn_res.2.4" 

NS$seurat_clusters <- NS$NS.RNA_snn_res.2.4

#Project cells into UMAP and tSNE space

NS <- RunTSNE(NS, reduction = "NS.integrated.cca", dims = 1:n.dims,
              reduction.name = "NS.tsne", reduction.key = "nstsne_") 

NS <- RunUMAP(NS, reduction = "NS.integrated.cca",
                   dims = 1:n.dims, reduction.name = "NS.umap",
              reduction.key = "NSumap_")

col1 <- c("#35608D", "#66CC33", "#E31A1C")
col2 <- c("lightgrey", "#330066")
col3 <- c("#35608D", "#66CC33")
nb.cols <- 14
col4 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 35
col5 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
nb.cols <- 24
col6 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

DimPlot(NS, reduction = "NS.tsne", 
        group.by = "orig.ident", cols = alpha(col4, 0.5))+
  coord_fixed()
ggsave("NS.tSNE_orig.ident.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")


DimPlot(NS, reduction = "NS.umap", 
        group.by = "orig.ident", cols = alpha(col4, 0.5), label.size = 2, pt.size = 0.5)+
  coord_fixed()
ggsave("NS.UMAP_orig.ident.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(NS, reduction = "NS.umap", 
        cols = alpha(col5, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()
ggsave("NS.UMAP_clusters_res2.4.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(NS, reduction = "NS.umap", 
        cols = alpha(col5, 0.5), split.by = "stage",
        pt.size = 0.5, label = TRUE, ncol = 4)+
  coord_fixed()
ggsave("NS.UMAP_stage_clusters_res2.4.pdf", device= "pdf", width = 48, 
       height = 32, units = "cm") 


DimPlot(NS, reduction = "NS.tsne", 
        cols = alpha(col5, 0.5),
        pt.size = 0.5, label = TRUE)+
  coord_fixed()
ggsave("NS.tsne_clusters_res2.4.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

DimPlot(NS, reduction = "NS.umap", 
        cols = alpha(col6, 0.5),
        pt.size = 0.5, label = TRUE, group.by = "cell.type.1")+
  coord_fixed()
ggsave("NS.umap_cell.type.1_iGmT.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 


DimPlot(NS, reduction = "NS.umap", 
        #cols = alpha(col6, 0.5),
        pt.size = 0.5, label = TRUE, group.by = "KHTissue")+
  coord_fixed()
ggsave("NS.umap_KHTissue_iGmT.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm") 

############# Lineage

FeaturePlot(NS, features = c("KY21:KY21.Chr5.707", #Dmrt1 _a-lineage
                             "KY21:KY21.Chr6.345", #Flvcr2_ A-lineage
                             "KY21:KY21.Chr2.1031"), #Msx_ b-lineage
            reduction = "NS.umap", ncol = 3, cols = col2,
            min.cutoff = "q05", max.cutoff = "q95")+
  coord_fixed()
ggsave("NS.iGmT.umap_Dmrt-Foxb-Msx.pdf", device= "pdf", width = 48, 
       height = 16, units = "cm") 


DotPlot(NS, features = c("KY21:KY21.Chr5.707", #Dmrt1 _a-lineage
                         "KY21:KY21.Chr6.345", #Flvcr2_ A-lineage
                         "KY21:KY21.Chr2.1031"), #Msx_ b-lineage
               cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/6)
ggsave("NS.iGmT.dotplot_Dmrt-Foxb-Msx.pdf", device= "pdf", width = 48, 
       height = 16, units = "cm") 

FeaturePlot(NS, features = "KY21:KY21.Chr6.345", #Flvcr2_ A-lineage
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr2.1031", #Msx_ b-lineage
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")



DimPlot(NS, group.by = "Phase", ncol = 5,
        reduction = "NS.umap", split.by = "orig.ident")+ coord_fixed()
ggsave("NS.iGmT.umap_Phase_replicate.pdf", device= "pdf", width = 48, 
       height = 16, units = "cm") 

save.image(object)

# Identify cluster types (in v5 no more difference between integrated and RNA assay slot)

Idents(NS) <- "seurat_clusters"

NS.markers <- FindAllMarkers(NS, only.pos = TRUE,
                                  min.pct = 0.5,
                                  logfc.threshold = 0.5,
                                  test.use = "roc",
                                  slot = "data")

save.image(object)

#annotate markers
human.homo <- read.csv("/scratch/gpfs/LEVINE/llemaire/cionaGeneModel/HT.KY21Gene.2.gff3/23.11.01CionaHomolog/gene-homolog2.csv", header = TRUE, row.names = 1)

human.homo <- as.data.table(human.homo)

human.homo <- human.homo %>% mutate(human.homolog=coalesce(human.homolog,gene))

NS.markers <- tidyft::left_join(NS.markers,human.homo,
                                     by = "gene")
write.table(NS.markers, file="NS.DEG_combined_res2.4.txt",
            quote=F, sep="\t", col.names=NA)

save.image(file = object)

#### Check marker expression

###Nervous system 
FeaturePlot(NS, features = "KY21:KY21.Chr6.58", #ETR = Celf3
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr2.1203", #syt1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr10.367", #PC2
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr11.637", #tub
            reduction = "NS.umap", 
            min.cutoff = "q15", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr14.252", #Nut1
            reduction = "NS.umap", 
            min.cutoff = "q15", max.cutoff = "q95")+coord_fixed()

DotPlot(NS, features = c("KY21:KY21.Chr6.58",#Etr
                         "KY21:KY21.Chr11.637", #Tub
                         "KY21:KY21.Chr14.252", # Nut1
                         "KY21:KY21.Chr2.1203", #Syt1
                         "KY21:KY21.Chr10.367", #PC2
                         "KY21:KY21.Chr9.941", # Syt9
                         "KY21:KY21.Chr7.872" #Uxs1
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)




FeaturePlot(NS, features = "KY21:KY21.Chr1.145", #sspop
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr11.1129", #Foxa
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

##Coronet cells

FeaturePlot(NS, features = "KY21:KY21.Chr11.539", #Ptf1a
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

####Pigment cells
FeaturePlot(NS, features = "KY21:KY21.Chr12.803",
            reduction = "NS.umap", max.cutoff = "q90", cols = col2) 
FeaturePlot(NS, features = "KY21:KY21.Chr5.389",
            reduction = "NS.umap", max.cutoff = "q90", cols = col2)


############ Hh2+ cells

FeaturePlot(NS, features = "KY21:KY21.Chr5.641", #Hh2
            reduction = "NS.umap", max.cutoff = "q90", cols = col2)
FeaturePlot(NS, features = "KY21:KY21.Chr12.575", #Slc23a1
            reduction = "NS.umap", max.cutoff = "q90", cols = col2)

DotPlot(NS, features = c("KY21:KY21.Chr5.641", #Hh2
                                                   "KY21:KY21.Chr12.575" #Slc23a1
                                                   ))+coord_flip()+ 
     theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
                     axis.text.y = element_text(size = 8),
                     plot.title = element_text(hjust = 0.5),
                     aspect.ratio=4/15)





#### Nerve cord
FeaturePlot(NS, features = "KY21:KY21.Chr1.567", #Collagen
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q100")

FeaturePlot(NS, features = "KY21:KY21.Chr7.1151", #Kdm8
            reduction = "NS.umap", split.by = "stage", 
            min.cutoff = "q05", max.cutoff = "q100")

FeaturePlot(NS, features = "KY21:KY21.Chr10.167", #Cbl
            reduction = "NS.umap", #split.by = "stage", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr9.606", #Lmx
            reduction = "NS.umap", split.by = "stage", 
            min.cutoff = "q05", max.cutoff = "q100")

FeaturePlot(NS, features = "KY21:KY21.Chr9.606", #Lmx
            reduction = "NS.umap", #split.by = "orig.ident", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

col7 <- rep(col2, 10)
DotPlot(NS, features = "KY21:KY21.Chr10.167", #Cbl
            group.by = "stage",
        cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/10)

DotPlot(NS, features = "KY21:KY21.Chr10.167", #Cbl
        group.by = "cell.type.1", 
        cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/8)

NS$stage_cell.type.1 <- paste(NS$stage, NS$cell.type.1, sep = "_")

DotPlot(NS, features = "KY21:KY21.Chr9.606", #lmx
        group.by = "stage_cell.type.1",
        cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/25)

DotPlot(NS, features = "KY21:KY21.Chr9.606", #lmx
        group.by = "orig.ident")+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/25)

DotPlot(NS, features = "KY21:KY21.Chr10.167", #cbl
        group.by = "stage_cell.type.1",
        cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/25)





FeaturePlot(NS, features = "KY21:KY21.Chr1.976", #Arx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q100")


### Nerve cord-b


FeaturePlot(NS, features = "KY21:KY21.Chr14.584", #Cdx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr4.737", # Ank1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr12.134", # Lmna
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr8.1024", # Kcnip4
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr1.1257", # Lama2
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr8.1169", # Wnt7
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr7.611", # SLit1 Marker for trunk nerve cord-b
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

#### PSC related
FeaturePlot(NS, features = "KY21:KY21.Chr13.415", #Sp8
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q100")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr2.892", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q100")+coord_fixed()

#### PSC
FeaturePlot(NS, features = "KY21:KY21.Chr10.318", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q100")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr1.2268", #bgcrystallin
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

#### Epidermis
FeaturePlot(NS, features = "KY21:KY21.Chr7.872", #Uxs1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr9.1071", #Umod
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

#### Photoreceptor
FeaturePlot(NS, features = "KY21:KY21.Chr1.838", # Opn1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

#### MHB????
FeaturePlot(NS, features = "KY21:KY21.Chr1.254", #Soxb1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()


##########
FeaturePlot(NS, features = "KY21:KY21.Chr8.654", #Foxd
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()



FeaturePlot(NS, features = "KY21:KY21.Chr7.360", # Dlx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr10.288", # Pax3/7
            reduction = "NS.umap", #split.by = "stage",
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr7.480", # FoxP
            reduction = "NS.umap", #split.by = "stage",
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()



########## Eminens neuron


FeaturePlot(NS, features = "KY21:KY21.Chr11.1031", # Prop
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
DotPlot(NS, features = "KY21:KY21.Chr11.1031")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/20)

###### Cholinergic neurons

FeaturePlot(NS, features = "KY21:KY21.Chr1.783", # vAcht
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()



###### GABAergic neurons

FeaturePlot(NS, features = "KY21:KY21.Chr7.1221", # Gad
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr2.793", # Viaat/vGat
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


###### Glutamatergic neurons

FeaturePlot(NS, features = "KY21:KY21.Chr3.1172", # vGlut
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")+coord_fixed()

### Endoderm
FeaturePlot(NS, features = "KY21:KY21.Chr11.1104", #thr
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr14.805", #sod1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

### Tectotial plate marker?
FeaturePlot(NS, features = "KY21:KY21.Chr1.2025", #tecta
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr8.832", #tecta
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

#### Pitx ANB

FeaturePlot(NS, features = "KY21:KY21.Chr4.347", #BMP
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr6.180", #Pitx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

### Motor ganglion neurons

FeaturePlot(NS, features = "KY21:KY21.Chr2.1436", #Kcnb1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr14.851", #Kcnb1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


#### VP and VPR neurons
 
FeaturePlot(NS, features = "KY21:KY21.Chr8.539", #Acta
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr6.340", #Vp
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr11.257", #C11.631
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")



###### aATENs neurons - 

FeaturePlot(NS, features = "KY21:KY21.Chr2.375", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


###### pATENs neurons - cluster

FeaturePlot(NS, features = "KY21:KY21.Chr7.62", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

########## 
FeaturePlot(NS, features = "KY21:KY21.Chr5.501", #Fgf8/17/18
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr7.709", #Hox12
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

############

FeaturePlot(NS, features = "KY21:KY21.Chr1.422", #Ebf
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
DotPlot(NS, features = "KY21:KY21.Chr1.422")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/20)


FeaturePlot(NS, features = "KY21:KY21.Chr14.647", #Otp
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
DotPlot(NS, features = "KY21:KY21.Chr14.647")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/20)
save.image(file = object)


####### engrailed ependymal cells

FeaturePlot(NS, features = "KY21:KY21.Chr7.541", # En
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr1.2091", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")



########## Pax2/5/8a Neck
FeaturePlot(NS, features = "KY21:KY21.Chr6.690", # Pax2/5/8a
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()


########### Mid hindbrain boundary

FeaturePlot(NS, features = "KY21:KY21.Chr3.942", # Rgs3
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr9.112", # Rgs8
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr9.113", # Rgs8
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr1.254", # Sox1/2/3
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr1.1871", # Sox14/15/22
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr5.501", # Fgf7/17/18
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

#### CESNs
FeaturePlot(NS, features = "KY21:KY21.Chr10.728", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr3.809", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr1.1", # Fgf11/12/13/14
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")


#### RTEN
FeaturePlot(NS, features = "KY21:KY21.Chr4.166", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr4.955", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")

FeaturePlot(NS, features = "KY21:KY21.Chr9.123", # Gnrh1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr4.1164", # Islet
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")

####### RTEN and CESNs
FeaturePlot(NS, features = "KY21:KY21.Chr12.875", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")
FeaturePlot(NS, features = "KY21:KY21.Chr8.242", # Atoh
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")

########## pATENs
FeaturePlot(NS, features = "KY21:KY21.Chr7.62", # Marker
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q90")

####### Photoreceptor
FeaturePlot(NS, features = "KY21:KY21.Chr1.838", # Opn1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr12.306", # Rax/Rx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr8.555", # Bsx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

########## Trunk nerve chord A
FeaturePlot(NS, features = "KY21:KY21.Chr7.1151", # Kdm8
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


FeaturePlot(NS, features = "KY21:KY21.Chr14.584", # Cdx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

####### engrailed ependymal cells

FeaturePlot(NS, features = "KY21:KY21.Chr7.541", # En
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr1.2091", # 
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")



########## Pax2/5/8a Neck
FeaturePlot(NS, features = "KY21:KY21.Chr6.690", # Pax2/5/8a
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


############## Mib+ RNs also called Lhx1+ GABAergic neurons

FeaturePlot(NS, features = "KY21:KY21.Chr11.888", # Mib
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q80")

FeaturePlot(NS, features = "KY21:KY21.Chr5.104", # Unc.a
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q80")


########### Antenna cells

FeaturePlot(NS, features = "KY21:KY21.Chr2.734", #AmpaR
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr2.1188", #Prod2 -OACCs
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr2.396", #Capn15 -OACCs
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


FeaturePlot(NS, features = "KY21:KY21.Chr10.826", #Alox5/lox5 -OACCs
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

### Motor ganglion neurons

FeaturePlot(NS, features = "KY21:KY21.Chr11.296", #Amd+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr14.610", #Amd+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr6.281", #Amd+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr8.285", #Amd+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr6.520", #Glra1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr1.1190", #Glra1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr2.1436", #Kcnb1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr14.851", #Kcnb1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr10.362", #Kcnb1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr13.438", #Kcnb1+ MG
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

##### Six3/6 Pro-aSV

FeaturePlot(NS, features = "KY21:KY21.Chr9.751", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

######## Eminens neurons


FeaturePlot(NS, features = "KY21:KY21.Chr11.1031", #Prop
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

########### aATENs


FeaturePlot(NS, features = "KY21:KY21.Chr2.375", #Rhag
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

######## Hox expression

FeaturePlot(NS, features = "KY21:KY21.Chr1.2230", #Hox5
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr1.841", #Hox1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr1.2091", #Hox3
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr1.2189", #Hox10 (part of cluster with Foxd+ cells but only at LTBII)
            reduction = "NS.umap", # split.by = "stage",
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr8.654", #foxd
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr7.709", #Hox12
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

########## hh2+ cells

FeaturePlot(NS, features = "KY21:KY21.Chr5.641", #Hh2
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr12.575", #Slc23a1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr2.1031", #Msx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

###########33 Trunk nerve cord A

FeaturePlot(NS, features = "KY21:KY21.Chr2.199", #Slc22a16
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(NS, features = "KY21:KY21.Chr4.737", #Ank1 (Shili reporter)
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr11.1129", #Foxa2
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr2.376", #Slc22a3
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr8.986", # Shili reporter
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr6.568", # (Nerve chord b? Shili reporter)
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr2.369", # Shili reporter
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr7.1151", # Shili reporter
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


###########33 Tail nerve cord A

FeaturePlot(NS, features = "KY21:KY21.Chr1.1326", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr1.2230", # Hox5
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr13.132", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(NS, features = "KY21:KY21.Chr11.1298", #Fgf3/7/10/22
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr3.881", #Efna.d
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr3.875", #Efna.b
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()


################# Pigment cells

FeaturePlot(NS, features = "KY21:KY21.Chr10.578", #Mitf - otolith
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr5.389", #Tyr
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr12.803", #Tyrp
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr2.1075", #C2.1036 - ocellus
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

##########BTNs

FeaturePlot(NS, features = "KY21:KY21.Chr1.2051", #Asic
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

DotPlot(NS, features = "KY21:KY21.Chr1.2051")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/15)

FeaturePlot(NS, features = "KY21:KY21.Chr8.1318", #Hmx
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

DotPlot(NS, features = "KY21:KY21.Chr8.1318")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/15)

DotPlot(NS, features = "KY21:KY21.Chr1.1164")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/15)

FeaturePlot(NS, features = "KY21:KY21.Chr12.81", #Grid1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

DotPlot(NS, features = "KY21:KY21.Chr12.81")+coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/15)

DotPlot(NS, features = "KY21:KY21.Chr3.561")+coord_flip()+ #Mdga2
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/15)


FeaturePlot(NS, features = "KY21:KY21.Chr14.647", #Otp
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

FeaturePlot(NS, features = "KY21:KY21.Chr4.1089", #Slc39a2
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()

#######Switch neurons
FeaturePlot(NS, features = "KY21:KY21.Chr2.1092", #Neurod/Ascal1
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()
FeaturePlot(NS, features = "KY21:KY21.Chr9.240", #Lhx2/9
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")+coord_fixed()



DimPlot(NS, reduction = "NS.umap", 
        cells.highlight = list(WhichCells(NS, idents = 4),
                               WhichCells(NS, idents = 31)),
        cols.highlight = c("darkred","darkblue"), sizes.highlight = 0.3,
        pt.size = 0.2) +coord_fixed()+
  scale_color_manual(labels = c("Other clusters",
                               "31",
                              "4"), values = c("grey","darkred","darkblue")) +
  #labs(color = "Analyzed Cell")+
  theme(plot.title = element_text(hjust = 0.5))
#  labs( y = "T-SNE_2" , x= "T-SNE_1", title = "Analyzed Cell Types") 
ggsave("NS_iGmT_UMAP_.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")



### Cluster 5

DotPlot(NS, features = c("KY21:KY21.Chr4.768",#Foxb
                         "KY21:KY21.Chr2.996", #Gsx
                         "KY21:KY21.Chr1.188", #Eph.a
                         "KY21:KY21.Chr2.824", #Fgf9/16/20
                         "KY21:KY21.Chr1.1192", #Myt1
                         "KY21:KY21.Chr4.720", #Otx
                         "KY21:KY21.Chr9.1050" #Pax6
                         ))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)

#### Row IV-VI - Palp
DotPlot(NS, features = c("KY21:KY21.Chr10.262",#Six3/6
                         "KY21:KY21.Chr3.541", #Six1/2
                         "KY21:KY21.Chr7.1172", #Eya
                         "KY21:KY21.Chr12.158" #FoxC
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)


##### Six3/6 Pro-aSV

FeaturePlot(NS, features = "KY21:KY21.Chr9.751", #
            reduction = "NS.umap", 
            min.cutoff = "q05", max.cutoff = "q95")


#### ANB

DotPlot(NS, features = c("KY21:KY21.Chr6.180",#Pitx
                         "KY21:KY21.Chr7.360", #Dlx
                         "KY21:KY21.Chr7.1172", #Eya
                         "KY21:KY21.Chr12.158" #FoxC
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)

### Cluster 11

DotPlot(NS, features = c("KY21:KY21.Chr5.1052",#Sall1
                        "KY21:KY21.Chr5.118", #Max
                        "KY21:KY21.Chr10.346" #Ets1/2.b
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)

### Cluster 12 Midline epidermis


DotPlot(NS, features = c("KY21:KY21.Chr10.288",#Pax3/7
                         "KY21:KY21.Chr2.1031", #Msx
                         "KY21:KY21.Chr14.584", # Cdx
                         "KY21:KY21.Chr4.1174", #Wnt5
                         "KY21:KY21.Chr7.361", #Dlx
                         "KY21:KY21.Chr7.709", #Hox12
                         "KY21:KY21.Chr5.695" #Klf2/4
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)

### Cluster 15 and 21 b lineage

DotPlot(NS, features = c("KY21:KY21.Chr10.288",#Pax3/7
                         "KY21:KY21.Chr2.1031", #Msx
                         "KY21:KY21.Chr14.584", # Cdx
                         "KY21:KY21.Chr8.1169", #Wnt7
                         "KY21:KY21.Chr4.947", #Tbx2/3
                         "KY21:KY21.Chr7.361", #Dlx
                         "KY21:KY21.Chr7.872", #Uxs1
                         "KY21:KY21.Chr1.1701", #Hand-r
                         "KY21:KY21.Chr7.709", #Hox12
                         "KY21:KY21.Chr9.606" #Lmx
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)


#### Cluster 17 - motor ganglion

DotPlot(NS, features = c("KY21:KY21.Chr9.1050",#Pax6 - MN1
                         "KY21:KY21.Chr13.457", #Lhx3/4
                         "KY21:KY21.Chr1.422", # Ebf
                         "KY21:KY21.Chr6.231", #Onecut
                         "KY21:KY21.Chr6.434", #Neurog
                         "KY21:KY21.Chr7.541", # En
                         "KY21:KY21.Chr4.1164", #Islet -MN2
                         "KY21:KY21.Chr5.29" #Nkx6 - MN1
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)


### Cluster 31

DotPlot(NS, features = c("KY21:KY21.Chr2.996",#Gsx
                         "KY21:KY21.Chr11.321", #Piwi
                         "KY21:KY21.Chr9.687", # Foxh.a
                         "KY21:KY21.Chr6.231", #Onecut
                         "KY21:KY21.Chr4.720", #Otx
                         "KY21:KY21.UAContig4.11", # Ptprq
                         "KY21:KY21.Chr8.555", # Bsx
                         "KY21:KY21.Chr12.306", # Rax/Rx
                         "KY21:KY21.Chr1.838", # Opn1
                         "KY21:KY21.Chr9.606" #Lmx
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=4/15)

#### CLuster 34
DotPlot(NS, features = c("KY21:KY21.Chr10.362", #Mnx
                         "KY21:KY21.Chr7.709", #Hox12
                         "KY21:KY21.Chr14.584", #Cdx
                         "KY21:KY21.Chr9.1050" #Pax6
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)

#### CLuster 32
DotPlot(NS, features = c("KY21:KY21.Chr9.606", #Lmx
                         "KY21:KY21.Chr1.422", #Ebf
                         "KY21:KY21.Chr6.231", #Onecut
                         "KY21:KY21.Chr6.434" #Neurog
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/15)

### Cluster29
DotPlot(NS, features = c("KY21:KY21.Chr4.720", #Otx
                         "KY21:KY21.Chr10.288",#Pax3/7
                         "KY21:KY21.Chr2.1031", #Msx
                         "KY21:KY21.Chr9.606" #Lmx
))+coord_flip()+ 
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=4/15)


#### Assign cell type based on clustering 2.4

cell.type <- read.csv("26.02.08_NS.cell.types.csv",header=TRUE)
metadata<- NS@meta.data

head(metadata)
head(cell.type)

cell.type1 <- subset(cell.type, select = c("cluster", "cell.type.1", "cell.type.2", "lineage"))
head(cell.type1)
colnames(cell.type1)[1] <- "seurat_clusters"
colnames(cell.type1)[2] <- "cell.type.3"
cell.type1$seurat_clusters <- factor(cell.type1$seurat_clusters)
class(cell.type1$seurat_clusters)

md2 <- subset(metadata,select = "seurat_clusters")
head(md2)
md2 <- rownames_to_column(md2)

md2 <- left_join(md2, cell.type1, by = "seurat_clusters")
md2 <- column_to_rownames(md2)

NS <- AddMetaData(NS, md2)

nb.cols <- 30
col8 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)


DimPlot(NS, reduction = "NS.umap", group.by = "cell.type.3", 
        cols = col8)+coord_fixed()
ggsave("26.02.08_NS.UMAP_cell.types.pdf", device= "pdf", width = 40, 
       height = 30, units = "cm") 
DimPlot(NS, reduction = "NS.umap", 
        group.by = "cell.type.2", cols = col8)+coord_fixed()
ggsave("26.02.08_NS.UMAP_future.cell.types.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

DimPlot(NS, reduction = "NS.umap", 
        group.by = "lineage", cols = col6)+coord_fixed()
ggsave("26.02.08_NS.UMAP_lineage.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

save.image(object)

NS$stage_lineage <- paste(NS$stage, NS$lineage, sep = "_")

DotPlot(NS, features = "KY21:KY21.Chr9.606", #lmx
        group.by = "stage_lineage", cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/25)

DotPlot(NS, features = "KY21:KY21.Chr10.167", #cbl
        group.by = "stage_lineage",
        cols = col2)+ coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/25)

save.image(object)




##################### Subset b-lineage

Idents(NS) <- "lineage"
head(Idents(NS))

blineage <- subset(NS, idents = "b-lineage" )

DimPlot(blineage, group.by = "stage", reduction = "NS.umap")

saveRDS(blineage, file = "blineage.rds")

Idents(NS) <- "seurat_clusters"
save.image(object)
