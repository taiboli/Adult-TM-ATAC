library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(ggpubr)

addArchRGenome("mm10")
addArchRThreads(threads = 6) 

projDev = loadArchRProject("TMOnly//")


marker.genes = c("Cygb", "Gpc3", "Crym", "Fn1", "Gpc3", "Acta2")

p <- plotBrowserTrack(
  ArchRProj = projDev, 
  groupBy = "CellType", 
  geneSymbol = marker.genes, 
  # plotSummary = "geneTrack",
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes_TMOnly.pdf", 
        ArchRProj = projDev, 
        addDOC = FALSE, width = 5, height = 5)
