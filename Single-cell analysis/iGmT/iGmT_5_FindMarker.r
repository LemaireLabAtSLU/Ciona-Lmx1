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

######## This script is run via SLURM

wd1 <- "/YourDirectory"
setwd(wd1)
object <- "iGmT4.RData"
load(object)

print("Environment Loaded")

object <-"iGmT5.RData"

Idents(iGmT) <- "seurat_clusters"

iGmT.markers <- FindAllMarkers(iGmT, only.pos = TRUE,
                                  min.pct = 0.5,
                                  logfc.threshold = 0.5,
                                  test.use = "roc",
                                  slot = "data")

print("Markers Found")
save.image(object)


#annotate markers
human.homo <- read.csv("/scratch/gpfs/LEVINE/llemaire/cionaGeneModel/HT.KY21Gene.2.gff3/23.11.01CionaHomolog/gene-homolog2.csv", header = TRUE, row.names = 1)

human.homo <- as.data.table(human.homo)

human.homo <- human.homo %>% mutate(human.homolog=coalesce(human.homolog,gene))

iGmT.markers <- tidyft::left_join(iGmT.markers,human.homo,
                                     by = "gene")
write.table(iGmT.markers, file="iGmT_DEG_combined_res1.txt",
            quote=F, sep="\t", col.names=NA)

save.image(object)
print("Done")
