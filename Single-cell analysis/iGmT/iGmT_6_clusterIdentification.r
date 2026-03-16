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
object <- "iGmT5.RData"
load(object)
object <- "iGmT6.RData"


## KH tissue annotation (Cao et al, Nature, 2019)

KHclusterID <- read_delim("/YourDirectory/cellID_KH2012/ciona10stage.cluster.upload.new.txt")
KHclusterID[c("orig.ID","cell.name")] <- str_split_fixed(KHclusterID$NAME, "_", 2)

KHsample<- data.frame(table(KHclusterID$orig.ID))

colnames(KHsample) <- c("orig.ID", "cellNumber")

class(KHsample$orig.ID)

for (x in 1:nrow(KHsample))
{
  KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "C110.1")
  {"IG_rep1"}
  else
  {KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "C110.2")
  {"IG_rep2"}
  else
    {KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "earlyN.1")
    {"EN_rep1"}
    else
    {
      KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "earlyN.2")
      {"EN_rep2"}
      else
      {
        KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "ETB.1")
        {"ETB_rep1"} 
        else
        {
          KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "ETB.2")
          {"ETB_rep2"}
          else
          {
            KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "ITB.1")
            {"ITB_rep1"}
            else
            {
              KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "ITB.3")
              {"ITB_rep2"}
              else
              {
                KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "lateN.1")
                {"LN_rep1"}
                else
                {
                  KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "lateN.2")
                  {"LN_rep2"}
                  else
                  {
                    KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "LTB1.1")
                    {"LTBI_rep2"}
                    else
                      KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "LTB2.2")
                      {"LTBII_rep2"}
                    else
                    {
                      KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "lv.1")
                      {"larva_rep3"}
                      else
                      {
                        KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "lv.3")
                        {"larva_rep1"}
                        else
                        {
                          KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "lv.4")
                          {"larva_rep2"}
                          else
                            KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "midG.1")
                            {"MG_rep1"}
                          else
                          {
                            KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "midG.2")
                            {"MG_rep2"}
                            else
                            {
                              KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "MTB.1")
                              {"MTB_rep1"}
                              else
                              {
                                KHsample[x,"orig.ident"] <- if(KHsample[x,"orig.ID"] == "MTB.2")
                                {"MTB_rep2"}
                                else
                                {NA}
                              }
                            }
                              
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    }
  }
}

KHsample <- head(KHsample, -1)
KHsample2 <- KHsample[,-2]

KHclusterID <- KHclusterID[-1,]

KHclusterID <- KHclusterID[,-c(2,3)]

KHclusterID <- left_join(KHclusterID, KHsample2, by ="orig.ID")

KHclusterID$cell.name2 <- paste(KHclusterID$cell.name, KHclusterID$orig.ident, sep = "-") 
head(KHclusterID)
KHclusterID2 <- subset(KHclusterID, select = c("cell.name2", "Tissue Type"))
colnames(KHclusterID2) <- c("rowname", "KHTissue")
metadata <- iGmT@meta.data
metadata <- rownames_to_column(metadata)
head(metadata)
metadata <- left_join(metadata, KHclusterID2, by = "rowname")
metadata <- column_to_rownames(metadata)
md2 <- subset(metadata, select = c("KHTissue"))
head(md2)
iGmT <- AddMetaData(iGmT, md2)



col7 <- brewer.pal(8, "Paired")

barplot(1:8, col = col7)

DimPlot(iGmT, group.by = "KHTissue", cols = col7,
        raster=FALSE, reduction = "tsne") 
ggsave("tsne_KHtissue.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm")

DimPlot(iGmT, group.by = "KHTissue", raster=FALSE,
        cols = col7, reduction = "UMAP")
ggsave("UMAP_KHtissue.pdf", device= "pdf", width = 20, 
       height = 16, units = "cm")


save.image(object)

###investigate gene expression/marker genes

FeaturePlot(iGmT, features = c("KY21:KY21.Chr4.92", #CCNb1
                                  "KY21:KY21.Chr2.319"), #Ccna1
            ncol=3, cols = col2,
            min.cutoff = "q30", max.cutoff = "q95")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr12.306"), #Rx
            cols = col2,
            min.cutoff = "q05", max.cutoff = "q20")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr1.145"), #SSPO
            cols = col2,
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr13.415"), #Sp8
            cols = col2, reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr3.1055"), #Ryr
            cols = col2,  reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr1.838"), #Opn1/Opn1lw
            cols = col2,  reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr3.989", #Foxf
            reduction = "tsne", raster=FALSE,
            min.cutoff = "q05", max.cutoff = "q95")

DotPlot(iGmT, features = "KY21:KY21.Chr3.989") + #Foxf
  coord_flip()+
theme(axis.text.x = element_text(size = 6, angle=45, vjust=1, hjust=1),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      aspect.ratio=3/35)

DotPlot(iGmT, features = c("KY21:KY21.Chr7.872",#Usx
                           "KY21:KY21.Chr8.693", #FoxG
                           "KY21:KY21.Chr4.720", #Otx
                           "KY21:KY21.Chr5.695", #Klf
                           "KY21:KY21.Chr13.415")  ) + #Sp8/9
       
  coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=8/35)

DotPlot(iGmT, features = "KY21:KY21.Chr7.872") + #Usx
  coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=3/35)

DotPlot(iGmT, features = c("KY21:KY21.Chr9.606", #Lmx1
                           "KY21:KY21.Chr8.1169",#Wnt7
                           "KY21:KY21.Chr5.695", #Klf
                           "KY21:KY21.Chr2.1031", #Msx
                           "KY21:KY21.Chr6.345",#FLVCR2
                           "KY21:KY21.Chr7.1151")) + #Kdm8
  coord_flip()+
  theme(axis.text.x = element_text(size = 6, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=8/35)

  
FeaturePlot(iGmT, features = "KY21:KY21.Chr2.862", 
            reduction = "UMAP", raster=FALSE,
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr4.186", #"KH2012:KH.C4.697", #ATM 
                                  "KY21:KY21.Chr8.285", # "KH2012:KH.C8.738", #ATM TVCs sensory vesicle
                                  "KY21:KY21.Chr8.449", # "KH2012:KH.C8.796", #Rubber cells
                                  "KY21:KY21.Chr7.170",#"KH2012:KH.L141.73", #muscle Rubber cells 
                                  "KY21:KY21.Chr1.8"), #"KH2012:KH.C1.344"), #Stronger in b muscle
            reduction = "tsne", ncol=3,
            min.cutoff = "q10", max.cutoff = "q90")

#Rubber cells
FeaturePlot(iGmT, features = c("KY21:KY21.Chr9.822", #Wnttun5 
                                  "KY21:KY21.Chr8.449", #Rudder cells
                                  "KY21:KY21.Chr7.170",#"KH2012:KH.L141.73", #muscle Rubber cells 
                                  "KY21:KY21.Chr1.1701",#hand2 
                                  "KY21:KY21.Chr7.709"), #Hox12
            reduction = "tsne", ncol=3,
            min.cutoff = "q10", max.cutoff = "q90")

# -> Cluster 20

#Germ cells

pdf(file="Germ cells-UMAP", 
    width = 16, # The width of the plot in inches
    height = 8) # The height of the plot in inches
FeaturePlot(iGmT, features = c( "KY21:KY21.Chr11.316", #    "KH2012:KH.C11.383",# Brsk
                                   "KY21:KY21.Chr1.618", #"KH2012:KH.C1.755", #Pem1
                                   "KY21:KY21.Chr8.476", #"KH2012:KH.C8.370"?,
                                   "KY21:KY21.Chr4.1200" , #"KH2012:KH.S852.2", #Pem3
                                   "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl
                                   "KY21:KY21.Chr2.1118" ,#"KH2012:KH.C2.261"),
                                   "KY21:KY21.Chr1.230", #Pem4
                                   "KY21:KY21.Chr3.749", #Pem5
                                   "KY21:KY21.Chr3.118", #Pem6
                                   "KY21:KY21.Chr2.302", #Pen1
                                   "KY21:KY21.Chr14.518", #Pem12
                                   "KY21:KY21.Chr2.1111", #Pem13
                                   "KY21:KY21.Chr2.1153", #Epb4115
                                   "KY21:KY21.Chr8.1130", #Midn
                                   "KY21:KY21.Chr12.615", #Rev3l
                                   "KY21:KY21.Chr5.966", #Syde2
                                   "KY21:KY21.Chr4.491", # PTP-like
                                   "KY21:KY21.Chr9.1133", #Tdrd-r.a
                                   "KY21:KY21.Chr13.286", #Naa40
                                   "KY21:KY21.Chr4.118"),  #Pem2
            reduction = "UMAP", ncol=5,
            min.cutoff = "q10", max.cutoff = "q90")
dev.off()

FeaturePlot(iGmT, features = c( "KY21:KY21.Chr11.316", #    "KH2012:KH.C11.383",# Brsk
                                   "KY21:KY21.Chr1.618", #"KH2012:KH.C1.755", #Pem1
                                   "KY21:KY21.Chr8.476", #"KH2012:KH.C8.370"?,
                                   "KY21:KY21.Chr4.1200" , #"KH2012:KH.S852.2", #Pem3
                                   "KY21:KY21.Chr11.697" , #"KH2012:KH.L154.37",#Haus1 Ci-ZF114 // Zfl
                                   "KY21:KY21.Chr2.1118") ,#"KH2012:KH.C2.261"),
            reduction = "UMAP", ncol=3, cols = col2,
            min.cutoff = "q05", max.cutoff = "q80")
ggsave("UMAP-germcells.pdf", device= "pdf", width = 48, 
       height = 24, units = "cm")

save.image(object)

##Germ cells are cluster 19


##endoderm
FeaturePlot(iGmT, features = "KY21:KY21.Chr11.1104", #thr
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr14.805", #sod1
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")


## -> cluster 8, 13

##Epidermis
FeaturePlot(iGmT, features = "KY21:KY21.Chr1.2012", #tff1
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = c("KY21:KY21.Chr7.872"), #Uxs1
            cols = col2, reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = c("KY21:KY21.Chr4.1174"), #Wnt5
            cols = col2, reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr7.933"), #Btln9
            cols = col2, reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr2.382", #Pinhead, trunk ventral epidermis 
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr2.381", #Admp, trunk dorsal epidermis 
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr12.740", #Aqp8 Trunk epidermis 
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr1.258", #Efnb Tail epidermis 
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr1.1871", #Sox14/15/21
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

## -> cluster 1, 2, 3, 6, 7, 16 (Cluster 6 is midline epidermis)

##Notochord
FeaturePlot(iGmT, features = "KY21:KY21.Chr12.6", #brachyury
            reduction = "tsne", 
            min.cutoff = "q05", max.cutoff = "q95")
## ->cluster 5

##mesenchyme 
FeaturePlot(iGmT, features = "KY21:KY21.Chr9.100", #twist
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr1.1501", #Kdm8
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr5.716", #Kdm8
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr11.984", #Hlx
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

## ->cluster 11, 14, 15, 17

###Nervous system 
FeaturePlot(iGmT, features = "KY21:KY21.Chr6.58", #ETR = Celf3
            reduction = "tsne", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr2.1203", #syt1
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr11.637", #tub
            reduction = "UMAP", 
            min.cutoff = "q15", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr1.145", #sspop
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q100")

FeaturePlot(iGmT, features = "KY21:KY21.Chr11.1129", #Foxa
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q100")

FeaturePlot(iGmT, features = "KY21:KY21.Chr12.803",
            reduction = "tsne", max.cutoff = "q90", cols = col2) ##Pigment cells



## -> clusters 4, 6, 9, 10, 18, 20

##Tvc 
FeaturePlot(iGmT, features = "KY21:KY21.Chr1.1771", #Rhod/f
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr4.861", #Casq (cluster 35, 36, 46 = rubber cells)
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")
FeaturePlot(iGmT, features = "KY21:KY21.Chr3.978", #MesP
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q80")

# -> cluster12


##### Muscle

FeaturePlot(iGmT, features = "KY21:KY21.Chr11.474", #Mlx
            reduction = "UMAP", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr14.750", #Mrf
            reduction = "tsne", 
            min.cutoff = "q05", max.cutoff = "q95")

FeaturePlot(iGmT, features = "KY21:KY21.Chr1.1601",
            reduction = "UMAP", max.cutoff = "q90", cols = col2) ##Muscle

## -> cluster 12


#### If necessary, Explore cell cycle
Idents(iGmT) <- "seurat_clusters"
FeaturePlot(iGmT, features = c("KY21:KY21.Chr2.319", #ccn1a
                              "KY21:KY21.Chr4.92", #ccn1b
                              "KY21:KY21.Chr2.1280",#ki67
                              "KY21:KY21.Chr12.862", #Pcna
                              "KY21:KY21.Chr5.730"),#cdc25
            ncol = 3, reduction = "UMAP",
            min.cutoff = "q05", max.cutoff = "q80") 
ggsave("UMAP-cell cycle.pdf", device= "pdf", width = 72, 
       height = 24, units = "cm")

FeaturePlot(iGmT, features = c("KY21:KY21.Chr11.637", #tubb4b - NS
                                  "KY21:KY21.Chr9.690", #tubb4b
                                  "KY21:KY21.Chr11.636"),#tubb4b
            ncol = 3, min.cutoff = "q05", max.cutoff = "q80") 




#assign tissue types

tissue.type1 <- read.csv("iGmT.tissue.types.csv",header=TRUE)
tissue.type1 <- subset(tissue.type1, select = c("cluster","tissue.type.1","cell.type.1"))

metadata <- iGmT@meta.data

md3 <- subset(metadata,select = seurat_clusters)
md3 <- rownames_to_column(md3)
head(md3)
class(md3$seurat_clusters)
class(tissue.type1$cluster)
tissue.type1[,"cluster"]<-factor(tissue.type1[,"cluster"])
class(tissue.type1)
class(tissue.type1$cluster)

md3 <- left_join(md3,tissue.type1, by = c("seurat_clusters" = "cluster"))
md3 <- column_to_rownames(md3)

head(md3)

iGmT <- AddMetaData(iGmT, md3)

save.image(object)

DimPlot(iGmT, reduction = "tsne", cols = col7,
        group.by = "tissue.type.1", label.size = 2, pt.size = 0.5,
        label = FALSE)+
  coord_fixed()


#Paired palette

col8epi <-"#62A3CB"       #Blue
col8NS <- c("#6FBE58")#Green
col8endo <- "#EB4445" #Red
col8mes <-"#FDB55F"#Orange
col8TVCmusc2 <- c("#CAB2D6", "#6A3D9A") #Violet
col8TVCmusc <- "#6A3D9A"
col8germ<- "#B15928" #Brown
col8noto <- "#E8D079" #Yellow

barplot(1:2, col = c(col8germ, col8noto))

col8 <- c(col8epi, col8NS, col8endo, col8mes, col8TVCmusc, col8noto, col8germ)

iGmT$tissue.type.1 <- factor(x= iGmT$tissue.type.1, levels = c("Epidermis",
                                                       "Nervous System",
                                                       "Endoderm",
                                                       "Mesenchyme",
                                                       "Muscle and Heart Lineage",
                                                       "Notochord",
                                                       "Germline"))

iGmT$cell.type.1 <- factor(x= iGmT$cell.type.1, levels = c("Endoderm.1",
                                                            "Endoderm.2",
                                                           "Anterior Dorsal Epidermis.1",
                                                                 "Anterior Dorsal Epidermis.2",
                                                                 "Anterior Ventral Epidermis",
                                                                 
                                                                 "Midline Epidermis",
                                                                 "Posterior Epidermis.1",
                                                                 "Posterior Epidermis.2",
                                                           "Mesenchyme.1",
                                                                 "Mesenchyme.2",
                                                                 "Mesenchyme.3",
                                                                 "Mesenchyme.4",
                                                                 "Muscle cells and TVCs",
                                                                 
                                                                "A lineage Nervous System",
                                                                 "a lineage Nervous System.1",
                                                                 "a lineage Nervous System.2",
                                                                 "b lineage Nervous System",
                                                                 "Rubber cells",
                                                                 "Notochord",
                                                                 "Germ cells"))

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "tissue.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col8)+
  coord_fixed()
ggsave("UMAP.tissue.iGmT.pdf", device= "pdf", width = 30, 
       height = 20, units = "cm") 

DimPlot(iGmT, reduction = "tsne", 
        group.by = "tissue.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col8)+
  coord_fixed()
ggsave("tsne.tissue.iGmT.pdf", device= "pdf", width = 30, 
       height = 20, units = "cm") 


col9epi <- colorRampPalette(c("#A6CEE3" ,"#1F78B4"))(6)
col9NS <- colorRampPalette(c("#B2DF8A", "#33A02C"))(5)
col9endo <- colorRampPalette(c("#FB9A99", "#E31A1C"))(2)
col9mes <- colorRampPalette(c("#FDBF6F", "#FF7F00"))(4)
col9TVCmusc <- "#6A3D9A" #colorRampPalette(c("#CAB2D6", "#6A3D9A"))(2)
col9germ<- "#B15928"
col9noto <- "#E8D079" #c("#FFFF99","#E8D079")

col9<- c(col9endo, col9epi, col9mes, col9TVCmusc, col9NS, col9noto, col9germ)

barplot(1:20, col = col9)

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "cell.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col9)+
  coord_fixed()
ggsave("UMAP.cell.type.iGmT.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm") 


DimPlot(iGmT,  
        group.by = "cell.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col9)+
  coord_fixed()


DimPlot(iGmT, reduction = "tsne", 
        group.by = "cell.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col9)+
  coord_fixed()
ggsave("tsne.cell.type.1.iGmT.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm") 

save.image(object)

#### create annotation lineage nervous system or midline epidermis

metadata <- iGmT@meta.data

for (x in 1:nrow(metadata))
{
  if(metadata[x,"tissue.type.1"] == "Endoderm")
  {
  metadata[x,"tissue.type.2"] <- "Endoderm"
  metadata[x,"tissue.type.3"] <- "Endoderm"
  }
  else
  {
    if(metadata[x,"tissue.type.1"] == "Mesenchyme")
    {
      metadata[x,"tissue.type.2"] <- "Mesenchyme"
      metadata[x,"tissue.type.3"] <- "Mesenchyme"
    }
    else
    {
      if(metadata[x,"cell.type.1"] == "Anterior Dorsal Epidermis.1")
      {
        metadata[x,"tissue.type.2"] <- "Anterior Epidermis"
        metadata[x,"tissue.type.3"] <- "Epidermis"
      }
      else
      {
        if(metadata[x,"cell.type.1"] == "Anterior Dorsal Epidermis.2")
        {
          metadata[x,"tissue.type.2"] <- "Anterior Epidermis"
          metadata[x,"tissue.type.3"] <- "Epidermis"
        }
        else
        {
          if(metadata[x,"cell.type.1"] == "Anterior Ventral Epidermis")
          {
            metadata[x,"tissue.type.2"] <- "Anterior Epidermis"
            metadata[x,"tissue.type.3"] <- "Epidermis"
          }
          else
          {
            if(metadata[x,"cell.type.1"] == "Midline Epidermis")
            {
              metadata[x,"tissue.type.2"] <- "Midline Epidermis"
              metadata[x,"tissue.type.3"] <- "Epidermis"
            }
            else
            {
              if(metadata[x,"cell.type.1"] == "Posterior Epidermis.1")
              {
                metadata[x,"tissue.type.2"] <- "Posterior Epidermis"
                metadata[x,"tissue.type.3"] <- "Epidermis"
              }
              else
              {
                if(metadata[x,"cell.type.1"] == "Posterior Epidermis.2")
                {
                  metadata[x,"tissue.type.2"] <- "Posterior Epidermis"
                  metadata[x,"tissue.type.3"] <- "Epidermis"
                }
                else
                  {
                    if(metadata[x,"cell.type.1"] == "Posterior Epidermis.2")
                  {
                    metadata[x,"tissue.type.2"] <- "Posterior Epidermis"
                    metadata[x,"tissue.type.3"] <- "Epidermis"
                    }
                    else
                    {
                      if(metadata[x,"cell.type.1"] == "A lineage Nervous System")
                      {
                        metadata[x,"tissue.type.2"] <- "Nervous System"
                        metadata[x,"tissue.type.3"] <- "A lineage Nervous System"
                      }
                      else
                      {
                        if(metadata[x,"cell.type.1"] == "a lineage Nervous System.1")
                        {
                          metadata[x,"tissue.type.2"] <- "Nervous System"
                          metadata[x,"tissue.type.3"] <- "a lineage Nervous System"
                        }
                        else
                        {
                          if(metadata[x,"cell.type.1"] == "a lineage Nervous System.2")
                          {
                            metadata[x,"tissue.type.2"] <- "Nervous System"
                            metadata[x,"tissue.type.3"] <- "a lineage Nervous System"
                          }
                          else
                          { 
                          if(metadata[x,"cell.type.1"] == "b lineage Nervous System")
                          {
                            metadata[x,"tissue.type.2"] <- "Nervous System"
                            metadata[x,"tissue.type.3"] <- "b lineage Nervous System"
                          }
                          else
                            {
                              if(metadata[x,"cell.type.1"] ==  "Rubber cells")
                            {
                              metadata[x,"tissue.type.2"] <- "Nervous System"
                              metadata[x,"tissue.type.3"] <- "b lineage Nervous System"
                              }
                              else
                                if(metadata[x,"tissue.type.1"] ==  "Muscle and Heart Lineage")
                                {
                                  metadata[x,"tissue.type.2"] <- "Muscle and Heart Lineage"
                                  metadata[x,"tissue.type.3"] <- "Muscle and Heart Lineage"
                                }
                              else
                              {
                                if(metadata[x,"tissue.type.1"] ==  "Notochord")
                                {
                                  metadata[x,"tissue.type.2"] <- "Notochord"
                                  metadata[x,"tissue.type.3"] <- "Notochord"
                                }
                                else
                                {
                                  metadata[x,"tissue.type.2"] <- "Germ cells"
                                  metadata[x,"tissue.type.3"] <- "Germ cells"
                                }
                              }
                            }
                        }
                      }
                    }
                  }
              }
                
            }
          }
        }
      }
    }
  }
}}

md4 <- subset(metadata,select = c(tissue.type.2, tissue.type.3))

iGmT <- AddMetaData(iGmT, md4)

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "tissue.type.2", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col9)+
  coord_fixed()

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "tissue.type.3", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col9)+
  coord_fixed()

iGmT$tissue.type.2 <- factor(x= iGmT$tissue.type.2, levels = c("Anterior Epidermis",
                                                               "Midline Epidermis",
                                                           "Posterior Epidermis",
                                                           "Nervous System",
                                                           "Endoderm",
                                                           "Mesenchyme",
                                                           "Muscle and Heart Lineage",
                                                           "Notochord",
                                                           "Germ cells"))


iGmT$tissue.type.3 <- factor(x= iGmT$tissue.type.3, levels = c("Epidermis",
                                                               "a lineage Nervous System",
                                                               "b lineage Nervous System",
                                                               "A lineage Nervous System",
                                                               "Endoderm",
                                                               "Mesenchyme",
                                                               "Muscle and Heart Lineage",
                                                               "Notochord",
                                                               "Germ cells"))


col8epi.2 <-colorRampPalette(c("#B2DF8A", "#33A02C"))(3)
col8NS.2 <- "#62A3CB"
col8epi.3 <- "#6FBE58"
col8NS.3 <- colorRampPalette(c("#A6CEE3" ,"#1F78B4"))(3)

col8.2 <- c(col8epi.2, col8NS.2, col8endo, col8mes, col8TVCmusc, col8noto, col8germ)

col8.3 <- c(col8epi.3, col8NS.3, col8endo, col8mes, col8TVCmusc, col8noto, col8germ)
barplot(1:12, col = col8.3)

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "tissue.type.2", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col8.2)+
  coord_fixed()
ggsave("UMAP.tissue.type.2.iGmT.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm") 


DimPlot(iGmT, reduction = "UMAP", 
        group.by = "tissue.type.3", label.size = 2, pt.size = 0.2,
        label = FALSE, raster=FALSE, cols = col8.3)+
  coord_fixed()
ggsave("UMAP.tissue.type.3.iGmT.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm") 

#### create annotation lineage nervous system combine to stage

iGmT$tissue.type.3_stage <- paste(iGmT$tissue.type.3, iGmT$stage, sep = "_")

DimPlot(iGmT, reduction = "UMAP", 
        group.by = "tissue.type.3_stage", pt.size = 0.2,
        label = FALSE, raster=FALSE)+
  coord_fixed()

save.image(object)

#Subset Nervous system (based on cell types to include midline epidermis) from iG to mT: cluster: 4, 6, 9, 10, 18, 20

head(iGmT@meta.data)

head(Idents(iGmT))

NS <- subset(iGmT, idents = c(4, 6, 9, 10, 18, 20))  



DimPlot(NS, reduction = "UMAP", 
        group.by = "cell.type.1", label.size = 2, pt.size = 0.2,
        label = FALSE)+
  coord_fixed()
 

saveRDS(NS, file = "NS.rds")

save.image(object)


