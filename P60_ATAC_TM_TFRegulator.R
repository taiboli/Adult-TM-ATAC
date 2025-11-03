library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(ggpubr)
library(dplyr)

addArchRGenome("mm10")
addArchRThreads(threads = 6) 

projDev = loadArchRProject("TMOnly//")

### REMOVE TFs with no expression 
gene.expression = fread("RNA_expression_all_tm.csv")


#########################
## positive and negative
## regulators
#########################


seGroupMotif <- getGroupSE(ArchRProj = projDev, useMatrix = "MotifMatrix", groupBy = "CellType")
## just deviation z-scores
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
## max between clusters 
## Now fixed May 2025
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x], useNames = FALSE)
}) %>% Reduce("cbind", .) %>% rowMaxs(useNames = FALSE)

### correlated TF gene ~ activity
corGSM_MM <- correlateMatrices(
  ArchRProj = projDev,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)


corGIM_MM <- correlateMatrices(
  ArchRProj = projDev,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)


corGEM_MM <- correlateMatrices(
  ArchRProj = projDev,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)


### max delta 
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGEM_MM$maxDelta <- rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]


## positive TF regulators 
this.cor = corGSM_MM
this.label = "GeneScore"


this.cor <- this.cor[order(abs(this.cor$cor), decreasing = TRUE), ]
this.cor <- this.cor[which(!duplicated(gsub("\\-.*","",this.cor[,"MotifMatrix_name"]))), ]
this.cor$TFRegulator <- "NO"
this.cor$TFRegulator[which(this.cor$cor > 0.5 & this.cor$padj < 0.01 & this.cor$maxDelta > quantile(this.cor$maxDelta, 0.75))] <- "Positive"
this.cor$TFRegulator[which(this.cor$cor < -0.5 & this.cor$padj < 0.01 & this.cor$maxDelta > quantile(this.cor$maxDelta, 0.75))] <- "Negative"

plot.data <- data.frame(this.cor)
plot.data$label = ifelse(plot.data$TFRegulator %in% c("Positive", "Negative"), plot.data[[paste0(this.label, "Matrix_name")]], NA)

plot.data$Expression = gene.expression$TM_expression[match(plot.data$GeneScoreMatrix_name, gene.expression$Gene_ID)]



p <- ggplot(plot.data %>% filter(Expression > 0), aes(cor, maxDelta, color = TFRegulator, label = label)) +
  geom_point() + 
  ggrepel::geom_label_repel() +
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "Positive"="#e41a1c", "Negative" = "#377eb8")) +
  # geom_text(subset(data.frame(corGSM_MM), TFRegulator=="YES"), aes(cor, maxDelta, label = GeneScoreMatrix_name)) +
  xlab(paste0("Correlation To ", this.label)) +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

ggsave(p, file = paste0("TMOnly///Plots/TMOnly_PositiveTF_AllTM_", this.label, "_Matrix_update_Oct27_RemoveNoExpression.png"), width = 8, height = 6)
ggsave(p, file = paste0("TMOnly///Plots/TMOnly_PositiveTF_AllTM_", this.label, "_Matrix_update_Oct27_RemoveNoExpression.pdf"), width = 8, height = 6)
ggsave(p, file = paste0("TMOnly///Plots/TMOnly_PositiveTF_AllTM_", this.label, "_Matrix_update_Oct27_RemoveNoExpression.svg"), width = 8, height = 6)
