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
setwd("/scratch/gpfs/LEVINE/llemaire/singleCellAnalysis/26.01.26_iGmT")


load(object)
object  <- "26.01.26_iGmT2.RData"

#### Normalized data for cell cycle score

# split the dataset into a list of layer within the same object

iGmT[["RNA"]] <- split(iGmT[["RNA"]], f = iGmT$orig.ident)

iGmT <- NormalizeData(iGmT) 

iGmT[["RNA"]] <- JoinLayers(iGmT[["RNA"]])

### Explore cell cycle and calculate phase

cellcycle <- read.csv("Ciona_CellCycleGenes.csv", header = TRUE)
s.genes <- cellcycle[cellcycle$role %in% "s.genes",]
s.genes <- s.genes$gene
g2m.genes <- cellcycle[cellcycle$role %in% "g2m.genes",]
g2m.genes <- g2m.genes$gene
iGmT <- CellCycleScoring(iGmT, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(iGmT, features = c("KY21:KY21.Chr2.319", #ccn1a
                                "KY21:KY21.Chr4.92", #ccn1b
                                "KY21:KY21.Chr2.1280",#ki67
                                "KY21:KY21.Chr12.862", #Pcna
                                "KY21:KY21.Chr5.730"),#cdc25
           ncol = 3)
ggsave("CellCycle.pdf", device= "pdf", width = 24, 
       height = 16, units = "cm")

Idents(iGmT) <- "orig.ident"

save.image(object)



