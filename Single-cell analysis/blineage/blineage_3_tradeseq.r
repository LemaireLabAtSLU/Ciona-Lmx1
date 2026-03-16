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
write(capture.output(sessionInfo()),file = "26.02.11sessionInfo.txt")

object <- "26.02.10blineage_2.RData"
load(object)
object <- "26.02.10blineage_3.RData"

#ggplot(npbb.sce, aes(x = Day, fill = is_even)) +
#  geom_density(alpha = .5) +
#  theme_minimal() +
#  scale_fill_brewer(type = "qual")

##Finding genes differentially expressed along pseudotime (TradeSeq)
##### This script mostly followed https://bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html

plotGeneCount(curve = cur4, clusters = npbb$cell.type.4,
              title = "Colored by cell types")+ xlab("UMAP-1") + ylab ("UMAP-2")

npbb.sce.count <- as.matrix(assays(npbb.sce)$counts) 

#Should I use logcounts or counts?
## In this tutorial it is counts and not logcounts
## https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html


w4 <- slingCurveWeights(cur4)

save.image(object)

set.seed(5)
icMat <- evaluateK(counts = npbb.sce.count,
                   pseudotime = pt4                      ,
                   cellWeights = w4,
                   nGenes = 300,
                   k = 3:10)

#select nknots = 7 (The middle panels show that the drop in AIC levels off if the number of knots is increased beyond 6
#, and we will choose that number of knots to fit the tradeSeq models)s

save.image(object)

set.seed(7)
pt4b <- slingPseudotime(cur4, na = FALSE) 
# to prevent the error: The pseudotimes contain NA values, and these cannot be used for GAM fitting.

#b.npb2 <- fitGAM(counts = npbb.sce.count,
#              pseudotime = pt4b, 
#              cellWeights = w4,
#              nknots = 6, verbose = TRUE)
#
##### too many genes, take too much time

##Most of the DE genes are only expressed in a few cells/1 cell. filtered them

geneFilter <- apply(assays(npbb.sce)$counts,1,function(x){
  sum(x >= 3) >= 5
})

### sum(x >= 3) >= 10
### Returns TRUE if the gene is expressed (≥3 counts)
### in at least 10 cells

head(geneFilter)
geneFilter <- as.data.frame(geneFilter)
table(geneFilter)
geneFilter <- rownames_to_column(geneFilter)
colnames(geneFilter)[1] <- "gene"
geneFilter <- rownames_to_column(geneFilter)

geneFilter1 <-  geneFilter[ which(geneFilter$geneFilter == TRUE), ]
geneFilter1 <- geneFilter1$gene

head(geneFilter1)
save.image(object)


npbb.sce.2 <- fitGAM(counts = npbb.sce.count,
                 pseudotime = pt4b, 
                 cellWeights = w4,
                 nknots = 7, verbose = TRUE, gene = geneFilter1)


plotGeneCount(curve = cur4,npbb.sce.count, clusters = npbb$cell.type.4,
              title = "Colored by cell types", models =npbb.sce.2 )+
  xlab("UMAP-1") + 
  ylab ("UMAP-2")+
  scale_color_manual(values = col4)+
  coord_fixed()
ggsave("fitGAM_UMAP_npbb_cell.type.4.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

plotGeneCount(curve = cur4,npbb.sce.count, clusters = npbb$npbb2.RNA_snn_res.0.9,
              title = "Colored by clusters", models =npbb.sce.2 )+
  xlab("UMAP-1") + 
  ylab ("UMAP-2")+
  scale_color_manual(values = col4)+
  coord_fixed()
ggsave("fitGAM_UMAP_npbb_res0.9.pdf", device= "pdf", width = 20, 
       height = 20, units = "cm")

table(rowData(npbb.sce.2)$tradeSeq$converged)

save.image(object) 

assocRes <- associationTest(npbb.sce.2, lineages = TRUE, contrastType = "consecutive")

#	
#The contrast used to impose the log2 fold-change threshold. Defaults to "start".
#Three options are possible: 
#- If "start", the starting point of each lineage is used to compare against all other points, and the fold-change threshold is applied on these contrasts.
#- If "end", the procedure is similar to "start", except that the reference point will now be the endpoint of each lineage rather than the starting point. 
#- If "consecutive", then consecutive points along each lineage will be used as contrasts. 
#This is the original way the associationTest was implemented and is kept for backwards compatibility. 
#If a fold change threshold has been set, we recommend users to use either the "start" or "end" options.
#
#### NB: number for DEGs for lineage 1 (rudder cells) is much higher with contrastType = "start" than with the consecutive contrary to lineage 3

head(assocRes)

deg.lineage1 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 1e-5)
]
head(deg.lineage1)
length(deg.lineage1)

deg.lineage2 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_2, "fdr") <= 1e-5)
]
head(deg.lineage2)
length(deg.lineage2)

deg.lineage3 <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_3, "fdr") <= 1e-5)
]
head(deg.lineage3)
length(deg.lineage3)



pdf("tradeseq.UMAP.summary.deg.pdf",
    width = 24,
    height = 12)
UpSetR::upset(fromList(list(lineage1 = deg.lineage1, 
                            lineage2 = deg.lineage2,
                            lineage3 = deg.lineage3)))
dev.off()




metrics_pt <- rownames_to_column(assocRes)
colnames(metrics_pt)[1] <- "gene"


metrics_pt$p.adjust1 <- p.adjust(metrics_pt$pvalue_1, "fdr")
metrics_pt$p.adjust2 <- p.adjust(metrics_pt$pvalue_2, "fdr")
metrics_pt$p.adjust3 <- p.adjust(metrics_pt$pvalue_3, "fdr")

################ Add human gene name most similar to KY21
human.homo <- read.csv("/scratch/gpfs/LEVINE/llemaire/cionaGeneModel/HT.KY21Gene.2.gff3/23.11.01CionaHomolog/gene-homolog2.csv", header = TRUE, row.names = 1)
human.homo <- as.data.table(human.homo)
human.homo <- human.homo %>% mutate(human.homolog=coalesce(human.homolog,gene))


metrics_pt <- tidyft::left_join(metrics_pt,human.homo,
                                by = "gene" )

write.table(metrics_pt, file="Associatio.test_metrics_pt_npbb_consecutive.txt",
            quote=F, sep="\t", col.names=NA)

pdf("Lmx over pseudotime.pdf",
    width = 6,
    height = 6)
plotSmoothers(npbb.sce.2, npbb.sce.count, gene = "KY21:KY21.Chr9.606",
              alpha = 1, border = TRUE)
dev.off()

#### To "hide" lineage plot them as transparent "#FFFFFF00"

plotSmoothers(npbb.sce.2, npbb.sce.count, gene = "KY21:KY21.Chr9.606",
              alpha = 1, border = TRUE, curvesCols = c("#FFFFFF00","#FFFFFF00", "#440154FF")) + 
 ggplot2::scale_color_manual(values = c( "#FFFFFF00", "#FFFFFF00", "#440154FF"))+
  scale_x_continuous(limits = c(0.0, 8.0))+
  scale_y_continuous(limits = c(0.0, 2.2))+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1.5)+
  labs( title = "Lmx1")
ggsave("lmx_pseudotime.lineage3.pdf", device= "pdf", width = 18, 
       height = 10, units = "cm")

plotSmoothers(npbb.sce.2, npbb.sce.count, gene = "KY21:KY21.Chr9.410",
              alpha = 1, border = TRUE, curvesCols = c("#FFFFFF00","#FFFFFF00", "#440154FF")) + 
  ggplot2::scale_color_manual(values = c("#FFFFFF00","#FFFFFF00", "#440154FF"))+
  scale_x_continuous(limits = c(0.0, 8.0))+
  scale_y_continuous(limits = c(0.0, 3.2))+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1.5)+
  labs( title = "Rac3")
ggsave("Rac3_pseudotime.lineage3.pdf", device= "pdf", width = 18, 
       height = 10, units = "cm")


plotSmoothers(npbb.sce.2, npbb.sce.count, gene = "KY21:KY21.Chr1.1257",
              alpha = 1, border = TRUE, curvesCols = c("#FFFFFF00","#FFFFFF00","#440154FF")) + 
  ggplot2::scale_color_manual(values = c( "#FFFFFF00","#FFFFFF00", "#440154FF"))+
  scale_x_continuous(limits = c(0.0, 8.0))+
  scale_y_continuous(limits = c(0.0, 3.2))+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1.5)+
  labs( title = "Laminin")
ggsave("Laminin_pseudotime.lineage3.pdf", device= "pdf", width = 18, 
       height = 10, units = "cm")

plotSmoothers(npbb.sce.2, npbb.sce.count, gene = "KY21:KY21.Chr10.167",
              alpha = 1, border = TRUE, curvesCols = c("#FFFFFF00","#FFFFFF00","#440154FF")) + 
  ggplot2::scale_color_manual(values = c( "#FFFFFF00","#FFFFFF00", "#440154FF"))+
  scale_x_continuous(limits = c(0.0, 8.0))+
  scale_y_continuous(limits = c(0.0, 3.2))+
  theme(axis.text.x = element_text(size = 10,  vjust=1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio=1/1.5)+
  labs( title = "Cbl")
ggsave("Cbl_pseudotime.lineage3.pdf", device= "pdf", width = 18, 
       height = 10, units = "cm")

save.image(object)

FeaturePlot(npbb, reduction = "npbb2.umap", features = c("KY21:KY21.Chr9.606",
                                                     "KY21:KY21.Chr9.410",
                                                     "KY21:KY21.Chr1.1257",
                                                     "KY21:KY21.Chr10.167"), #Cbl
            cols = col3, ncol = 4)+ coord_fixed()
ggsave("npbb2_Lmx_Rac3_Laminin.pdf", device= "pdf", width = 40, 
       height = 10, units = "cm")

save.image(object)

#########Discovering progenitor marker genes, not performed
#startRes <- startVsEndTest(npbb.sce.2, lineages = TRUE)
#oStart <- order(startRes$waldStat, decreasing = TRUE)
#sigGeneStart_2 <- names(npbb.sce.2)[oStart[2]]
#plotSmoothers(npbb.sce.2, assays(npbb.sce.2)$counts, gene = sigGeneStart_2)
#plotGeneCount(cur5, assays(npbb.sce.2)$counts, gene = sigGeneStart_2)
#save.image(object)

##############Between-lineage comparisons: Discovering differentiated cell type markers,not performed
##endRes <- diffEndTest(npbb.sce.2)
#o <- order(endRes$waldStat, decreasing = TRUE)
#sigGene3 <- names(npbb.sce.2)[o[3]]
#plotSmoothers(npbb.sce.2, assays(npbb.sce.2)$counts, "KH2012:KH.C9.439")
#plotGeneCount(cur5, assays(npbb.sce.2)$counts, gene = "KH2012:KH.C9.439")

###############Discovering genes with different expression patterns, not performed
#patternRes <- patternTest(npbb.sce.2)
#oPat <- order(patternRes$waldStat, decreasing = TRUE)
#head(rownames(patternRes)[oPat])
#plotSmoothers(npbb.sce.2, assays(npbb.sce.2)$counts, gene = rownames(patternRes)[oPat][5])

#######################Combining patternTest with diffEndTest results, not performed
#patternRes$Gene <- rownames(patternRes)
#patternRes$pattern <- patternRes$waldStat
#patternRes <- patternRes[, c("Gene", "pattern")]
#endRes$Gene <- rownames(endRes)
#endRes$end <- endRes$waldStat
#endRes <- endRes[, c("Gene", "end")]
#compare <- merge(patternRes, endRes, by = "Gene", all = FALSE)
#compare$tranpbbientScore <- 
#  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2
#ggplot(compare, aes(x = log(pattern), y = log(end))) +
#  geom_point(aes(col = tranpbbientScore)) +
#  labs(x = "patternTest Wald Statistic (log scale)",
#       y = "diffEndTest Wald Statistic (log scale)") +
#  scale_color_continuous(low = "yellow", high = "red") +
#  theme_classic()
#topTranpbbient <- compare[which.max(compare$tranpbbientScore), "Gene"]
#plotSmoothers(npbb.sce.2, assays(npbb.sce.2)$counts, gene = "KH2012:KH.C3.297")
#plotGeneCount(cur5, assays(npbb.sce.2)$counts, gene = topTranpbbient)
#head(
#  compare[order(compare$tranpbbientScore, decreasing = TRUE), "Gene"],
#  n = 10)

################Early drivers of differentiation, not performed

#plotGeneCount(curve = cur5, counts = assays(npbb.sce.2)$counts,
#              clusters = apply(slingClusterLabels(cur5), 1, which.max),
#              models = npbb.sce.2)
#earlyDERes <- earlyDETest( npbb.sce.2, knots = c(1, 2))
#oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
#head(rownames(earlyDERes)[oEarly])
#plotSmoothers(npbb.sce.2, assays(npbb.sce.2)$counts, gene = rownames(earlyDERes)[oEarly][6])
#save.image(object)

################################################# Beginning of slurm

#Clustering of genes according to their expression pattern
##Clustering using RSEC, clusterExperiment
# Tranfer reduce dimension to sce object wit tradeseq
reducedDims(npbb.sce.2)$PCA <- npbb.sce@int_colData@listData[["reducedDims"]]@listData[["NPBB2.PCA"]]

head(npbb.sce@int_colData@listData[["reducedDims"]]@listData[["NPBB2.PCA"]])
head(reducedDims(npbb.sce.2)$PCA)
reducedDimNames(npbb.sce.2)
listBuiltInReducedDims()
class(npbb.sce.2)

nPointsClus <- 20
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

save.image()

print("End Script")
################################################# End of slurm





