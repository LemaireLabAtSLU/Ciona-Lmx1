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
library(ggraph)
library("clustree")
library(viridisLite)
library("viridis")
library(circlize)

wd1 <- "/YourDirectory"

setwd(wd1)

inf <- sessionInfo()

write(capture.output(inf),file = "sessionInfo")

wd2 <- "/scratch/gpfs/LEVINE/llemaire/CellRanger"
wd3 <- "outs/filtered_feature_bc_matrix/"    ### use filtered matrices (empty droplets are removed) 
object <- "iGmT1.RData"

ids1 <- c("IG_rep1","IG_rep2",
          "MG_rep1", "MG_rep2",
          "EN_rep1", "EN_rep2",
          "LN_rep1", "LN_rep2",
          "ITB_rep1", "ITB_rep2",
          "ETB_rep1", "ETB_rep2",
          "MTB_rep1", "MTB_rep2") ### Not adding any other data files than the one published by Cao et al., 2019

data <- sapply(ids1, function(i){
  d10x <- Read10X(file.path(wd2,i,wd3))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

data <- do.call("cbind", data)

iGmT <- CreateAssayObject(counts = data) 

save.image(object)

iGmT <- CreateSeuratObject(counts =  iGmT, 
                            project = "iGmT", 
                            min.cells = 3, 
                            min.features = 200,   
                            names.field = 2,
                            names.delim = "\\-")

iGmT@meta.data[["orig.ident"]] <- factor(iGmT@meta.data[["orig.ident"]], levels = ids1)

head(Idents(iGmT))

samplename = iGmT@meta.data$orig.ident
unfilteredCell<- as.data.frame(table(samplename))

iGmT@meta.data[["orig.ident"]] <- factor(iGmT@meta.data[["orig.ident"]], levels = ids1)

Idents(iGmT) <- "orig.ident"

###Low-quality / dying cells often exhibit extensive mitochondrial contamination
###We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features


iGmT[["percent.mt"]] <- PercentageFeatureSet(iGmT, pattern = "^mt.")


VlnPlot(iGmT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("gene-RNA-raw.content-1.pdf", device= "pdf", width = 100, 
       height = 20, units = "cm")


VlnPlot(iGmT, features = "percent.mt") +
  geom_hline(yintercept=c(0.08, 20), color = "red", linewidth=2)
ggsave("percent.mt-raw.pdf", device= "pdf", width = 40, 
       height = 20, units = "cm")

VlnPlot(iGmT, features = "nFeature_RNA") +
  geom_hline(yintercept=1000, color = "red", linewidth=2)
ggsave("gene-raw.pdf", device= "pdf", width = 40, 
         height = 20, units = "cm")
  
VlnPlot(iGmT, features = "nCount_RNA") +
  geom_hline(yintercept=mean(iGmT$nCount_RNA) + 5*sd(iGmT$nCount_RNA), color = "red", linewidth=2)
ggsave("gene-raw_mean+5SD.pdf", device= "pdf", width = 40, 
         height = 20, units = "cm")
  
rm(data)

save.image(object)

##Samples quality control: 
# minimum 1000 for genes and less than mean +5SD for nCount, percent mito between 0.08 and 20.

iGmT <- SplitObject(iGmT, split.by = "orig.ident")

qc_thresholds <- list(
  "IG_rep1" = list(min_features=1000,
                   max_features=mean(iGmT[["IG_rep1"]]$nCount_RNA) + 5*sd(iGmT[["IG_rep1"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "IG_rep2" = list(min_features=1000,
                   max_features=mean(iGmT[["IG_rep2"]]$nCount_RNA) + 5*sd(iGmT[["IG_rep2"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "MG_rep1" = list(min_features=1000,
                   max_features=mean(iGmT[["MG_rep1"]]$nCount_RNA) + 5*sd(iGmT[["MG_rep1"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "MG_rep2" = list(min_features=1000,
                   max_features=mean(iGmT[["MG_rep2"]]$nCount_RNA) + 5*sd(iGmT[["MG_rep2"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "EN_rep1" = list(min_features=1000,
                    max_features=mean(iGmT[["EN_rep1"]]$nCount_RNA) + 5*sd(iGmT[["EN_rep1"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08),
  "EN_rep2" = list(min_features=1000,
                   max_features=mean(iGmT[["EN_rep2"]]$nCount_RNA) + 5*sd(iGmT[["EN_rep2"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "LN_rep1" = list(min_features=1000,
                   max_features=mean(iGmT[["LN_rep1"]]$nCount_RNA) + 5*sd(iGmT[["LN_rep1"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "LN_rep2" = list(min_features=1000,
                   max_features=mean(iGmT[["LN_rep2"]]$nCount_RNA) + 5*sd(iGmT[["LN_rep2"]]$nCount_RNA),
                   max_mt=20, min_mt=0.08),
  "ITB_rep1" = list(min_features=1000,
                    max_features=mean(iGmT[["ITB_rep1"]]$nCount_RNA) + 5*sd(iGmT[["ITB_rep1"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08),
  "ITB_rep2" = list(min_features=1000,
                    max_features=mean(iGmT[["ITB_rep2"]]$nCount_RNA) + 5*sd(iGmT[["ITB_rep2"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08),
  "ETB_rep1" = list(min_features=1000,
                    max_features=mean(iGmT[["ETB_rep1"]]$nCount_RNA) + 5*sd(iGmT[["ETB_rep1"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08),
  "ETB_rep2" = list(min_features=1000,
                    max_features=mean(iGmT[["ETB_rep2"]]$nCount_RNA) + 5*sd(iGmT[["ETB_rep2"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08),
  "MTB_rep1" = list(min_features=1000,
                    max_features=mean(iGmT[["MTB_rep1"]]$nCount_RNA) + 5*sd(iGmT[["MTB_rep1"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08),
  "MTB_rep2" = list(min_features=1000,
                    max_features=mean(iGmT[["MTB_rep2"]]$nCount_RNA) + 5*sd(iGmT[["MTB_rep2"]]$nCount_RNA),
                    max_mt=20, min_mt=0.08))


for (s in names(iGmT)) 
  {
  th <- qc_thresholds[[s]]
  iGmT[[s]] <- subset(iGmT[[s]], subset = nFeature_RNA > th$min_features & 
                                                  nCount_RNA < th$max_features & 
                                                  percent.mt <= th$max_mt &
                                                  percent.mt >= th$min_mt )
  }

iGmT <- merge(iGmT[[1]], y = iGmT[-1])

VlnPlot(iGmT, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)
ggsave("gene-RNA-filtered.content.pdf", device= "pdf", width = 100, 
       height = 20, units = "cm")

head(Idents(iGmT))







samplename = iGmT@meta.data$orig.ident
filteredCell<- as.data.frame(table(samplename))
colnames(filteredCell)[2] <- "Filtered"
colnames(unfilteredCell)[2] <- "Unfiltered"

cellNumber <- left_join(unfilteredCell, filteredCell, by = "samplename")

save.image(object)

rm(filteredCell)
rm(unfilteredCell)
rm(qc_thresholds)
rm(th)
rm(samplename)

object  <- "iGmT2.RData"
save.image(object)






