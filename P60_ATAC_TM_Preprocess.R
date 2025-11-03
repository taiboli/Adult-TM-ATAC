library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRGenome("mm10")
addArchRThreads(threads = 6) 

ArrowFiles = c("ArchR/P60_WT.arrow")
names(ArrowFiles) <- c("P60")

projDev <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR/",
  copyArrows = TRUE 
)

projDev <- filterDoublets(projDev)


### Dimensionality reduction

projDev <- addIterativeLSI(
  ArchRProj = projDev,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  force = T
)

### Clustering 
projDev <- addClusters(
  input = projDev,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)


### UMAP and visualization
projDev <- addUMAP(
  ArchRProj = projDev, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)


# p1 <- plotEmbedding(ArchRProj = projDev, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = projDev, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# ggAlignPlots(p1, p2, type = "h")
plotPDF(p1, name = "Plot-UMAP_All_P60.pdf", ArchRProj = projDev, addDOC = FALSE, width = 5, height = 5)

#################################
##### snRNA-seq integration
#################################
data <- readRDS("12072023_P60_wt_single_nuc_atac_pom.RDS")
data = RenameCells(object = data, new.names = paste0("P60_WT#", colnames(data)))

data = subset(data, CellType %in% c("TM1", "TM2", "TM3"))


common.cells = intersect(colnames(data), projDev$cellNames)

data = subset(data, cells = common.cells)


projDev = subsetArchRProject(
  ArchRProj = projDev,
  cells = common.cells,
  outputDirectory = "TMOnly//",
  dropCells = TRUE,
  logFile = NULL,
  threads = 1,
  force = TRUE
)


#################################

projDev <- addGeneIntegrationMatrix(
  ArchRProj = projDev, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = data,
  addToArrow = TRUE,
  force= TRUE,
  # groupList = groupList,
  groupRNA = "CellType",
  nameCell = "predictedCell_Integrated",
  nameGroup = "predictedGroup_Integrated",
  nameScore = "predictedScore_Integrated"
)



#################################
### Peak calling
#################################
######### Add cluster information

celltype.data = data.frame(cell = colnames(data),
                           type = as.character(Idents(data)))

add.data = celltype.data$type
names(add.data) = celltype.data$cell
add.data = add.data[projDev$cellNames]



projDev = addCellColData(
  ArchRProj = projDev,
  data = celltype.data$type,
  name = "CellType",
  cells = celltype.data$cell,
  force = FALSE
)



### pseudobulk replicates

projDev <- addGroupCoverages(ArchRProj = projDev, groupBy = "CellType", force = T)

pathToMacs2 <- findMacs2()


projDev <- addReproduciblePeakSet(
  ArchRProj = projDev, 
  groupBy = "CellType", 
  pathToMacs2 = pathToMacs2
)


projDev <- addPeakMatrix(projDev)


### Marker peaks only in cell types with sufficient cells 
# count.data <- as.data.frame(table(projDev.TOM$predictedGroup_Integrated))
# count.data$Var1 <- as.character(count.data$Var1)
# to.choose <- count.data$Var1[which(count.data$Freq > 100)]


markersPeaks <- getMarkerFeatures(
  ArchRProj = projDev.TOM, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup_Integrated",
  # useGroups = to.choose,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5",
  labelRows = F,
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "TOM_SeuratIntegration_Peak-Marker-Heatmap.pdf", width = 8, height = 6, 
        ArchRProj = projDev, addDOC = FALSE)



############################################
## Motif enrichment
#############################################

projDev <- addMotifAnnotations(ArchRProj = projDev, motifSet = "cisbp", name = "Motif",
                                   force = T)


enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projDev,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEM, name = "TOM_SeuratIntegration_Motifs-Enriched-Marker-Heatmap", 
        width = 10, height = 6, ArchRProj = projDev, addDOC = FALSE)



df = data.table(assay(enrichMotifs))
df$TF = rownames(enrichMotifs)
fwrite(df, file = "TMOnly/TOM_KeyTF_ArchR.csv", row.names = F, quote = F)



saveArchRProject(ArchRProj = projDev.TOM, outputDirectory = "TMOnly//", load = FALSE)


