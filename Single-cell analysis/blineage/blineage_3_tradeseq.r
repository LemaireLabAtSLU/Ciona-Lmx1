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

object <- "blineage_2.RData"
load(object)
object <- "blineage_3.RData"


##Finding genes differentially expressed along pseudotime (TradeSeq)
##### This script mostly followed https://bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/tradeSeq.html

plotGeneCount(curve = cur4, clusters = npbb$cell.type.4,
              title = "Colored by cell types")+ xlab("UMAP-1") + ylab ("UMAP-2")

npbb.sce.count <- as.matrix(assays(npbb.sce)$counts) 

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
human.homo <- read.csv("/YourDirectory/cionaGeneModel/gene-homolog2.csv", header = TRUE, row.names = 1)
human.homo <- as.data.table(human.homo)
human.homo <- human.homo %>% mutate(human.homolog=coalesce(human.homolog,gene))


metrics_pt <- tidyft::left_join(metrics_pt,human.homo,
                                by = "gene" )

write.table(metrics_pt, file="blineage_Association.test_metrics_pt_npbb_consecutive.txt",
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

save.image(object)
