library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(ggpubr)
library(pheatmap)

addArchRGenome("mm10")
addArchRThreads(threads = 6) 

projDev = loadArchRProject("TMOnly//")

##########################################
### Confusion matrix: RNA vs ATAC clusters
##########################################

data <- readRDS("12072023_P60_wt_single_nuc_atac_pom.RDS")
data = RenameCells(object = data, new.names = paste0("P60_WT#", colnames(data)))
common.cells = intersect(colnames(data), projDev$cellNames)
data = subset(data, cells = common.cells)


cM <- confusionMatrix(paste0(projDev$CellType), paste0(projDev$Clusters2) )


cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
ggsave(p, file = "Confusion_AllTM_reclustering.png", width = 6, height = 6)
ggsave(p, file = "Confusion_AllTM_reclustering.pdf", width = 6, height = 6)


###########################################################
# heatmap of promotor openness vs gene expression
###########################################################

projDev <- addPeak2GeneLinks(
  ArchRProj = projDev,
  reducedDims = "IterativeLSI"
)

pals = c("#A5CDE2", "#E3191C", "#2F9E2A")
names(pals) = c("TM1", "TM2", "TM3")

p <- plotPeak2GeneHeatmap(ArchRProj = projDev, groupBy = "CellType", k = 2, 
                          palGroup = pals)

plotPDF(p,
        name = "TMOnly_gene2peak_heatmap_July4_Update.pdf",
        ArchRProj = projDev,
        addDOC = FALSE, width = 8, height = 10)


