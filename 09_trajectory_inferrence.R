setwd("/home/jcfan/human_brain/HC_scrnaseq/pseudotime/")

cM <- confusionMatrix(paste0(HC_TI$maintype), paste0(HC_TI$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
HC_TI <- myAddGeneIntegrationMatrix(
  ArchRProj = HC_TI, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = scRNA,
  addToArrow = TRUE,
  force= TRUE,
  groupRNA = "tmptype",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
HC_TI <- addGroupCoverages(ArchRProj = HC_TI, groupBy = "maintype")
###IHC
{
head(HC_TI$IHC[!is.na(HC_TI$IHC)])
p <- plotTectory(HC_TI, Tectory = "IHC", colorBy = "cellColData", name = "IHC")

HC_TI <- addImputeWeights(HC_TI)
TMM  <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTectoryHeatmap(TMM, pal = paletteContinuous(set = "solarExtra"))
TGSM <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra"))
TGIM <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTectoryHeatmap(TGIM,  pal = paletteContinuous(set = "blueYellow"))
plotPDF(p1, name = "peakHeatmap_Motif_IHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
plotPDF(p2, name = "peakHeatmap_GeneScore_IHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
plotPDF(p3, name = "peakHeatmap_Gene_IHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
TPM  <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTectoryHeatmap(TPM, pal = myPal)
plotPDF(p4, name = "peakHeatmap_peak_IHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
corGSM_MM <- correlateTectories(TGIM, TMM)

idxToRemove <- grep(pattern = "deviations", x = corGSM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove > 1)){
  corGSM_MM[["correlatedMappings"]] <- corGSM_MM[["correlatedMappings"]][-idxToRemove,]
}
corGSM_MM[["correlatedMappings"]]
TGSM2 <- TGSM[corGSM_MM[["correlatedMappings"]]$name1, ]
TMM2 <- TMM[corGSM_MM[["correlatedMappings"]]$name2, ]
TCombined <- TGSM2
assay(TCombined, withDimnames=FALSE) <- t(apply(assay(TGSM2), 1, scale)) + t(apply(assay(TMM2), 1, scale))
combinedMat <- plotTectoryHeatmap(TCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(TGSM2))
ht1 <- plotTectoryHeatmap(TGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTectoryHeatmap(TMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
pdf("heatmap_IHC_combined.pdf", width = 12, height = 8)
ComplexHeatmap::draw(ht1 + ht2)
dev.off()
}
#####OHC
head(HC_TI$OHC[!is.na(HC_TI$OHC)])
p <- plotTectory(HC_TI, Tectory = "OHC", colorBy = "cellColData", name = "OHC")
myPal <- colorRampPalette(c("#2C7BB6", "#ABD9E9", "#F7FBFF", "#CC00FF", "#FF7F0E", "#FFFF00"))(100)

TMM  <- getTectory(ArchRProj = HC_TI, name = "OHC", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTectoryHeatmap(TMM, pal = paletteContinuous(set = "solarExtra"))
TGSM <- getTectory(ArchRProj = HC_TI, name = "OHC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra"))
TGIM <- getTectory(ArchRProj = HC_TI, name = "OHC", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTectoryHeatmap(TGIM,  pal = paletteContinuous(set = "blueYellow"))
plotPDF(p1, name = "peakHeatmap_Motif_OHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
plotPDF(p2, name = "peakHeatmap_GeneScore_OHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
plotPDF(p3, name = "peakHeatmap_Gene_OHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
TPM  <- getTectory(ArchRProj = HC_TI, name = "OHC", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTectoryHeatmap(TPM, pal = myPal)
plotPDF(p4, name = "peakHeatmap_peak_OHC", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
corGSM_MM <- correlateTectories(TGIM, TMM)

idxToRemove <- grep(pattern = "deviations", x = corGSM_MM[["correlatedMappings"]]$name2)
if(length(idxToRemove > 1)){
  corGSM_MM[["correlatedMappings"]] <- corGSM_MM[["correlatedMappings"]][-idxToRemove,]
}
corGSM_MM[["correlatedMappings"]]
TGSM2 <- TGSM[corGSM_MM[["correlatedMappings"]]$name1, ]
TMM2 <- TMM[corGSM_MM[["correlatedMappings"]]$name2, ]
TCombined <- TGSM2
assay(TCombined, withDimnames=FALSE) <- t(apply(assay(TGSM2), 1, scale)) + t(apply(assay(TMM2), 1, scale))
combinedMat <- plotTectoryHeatmap(TCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(TGSM2))
ht1 <- plotTectoryHeatmap(TGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTectoryHeatmap(TMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
pdf("heatmap_OHC_combined.pdf", width = 12, height = 8)
ComplexHeatmap::draw(ht1 + ht2)
dev.off()
write.csv(p1@ht_list[["MotifMatrix"]]@matrix,file = "OHC_motif.csv")
write.csv(p2@ht_list[["GeneScoreMatrix"]]@matrix,file="OHC_geneactivity(peak).csv")
write.csv(p3@ht_list[["GeneIntegrationMatrix"]]@matrix,file="OHC_gene.csv")
write.csv(p4@ht_list[["PeakMatrix"]]@matrix,file = "OHC_peak.csv")
#####IHC
TGSM <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra") ,returnMat = TRUE)
rowOrder<-rownames(p)
p2 <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra"), rowOrder = rowOrder)
TGIM <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTectoryHeatmap(TGIM,  pal = paletteContinuous(set = "blueYellow"),returnMatrix = TRUE)
rowOrder<-intersect(rownames(p3),rowOrder)
TGIM2<-TGIM[rowOrder,]
TGSM2<-TGSM[rowOrder,]
rowOrder2 <- sub("^[^:]+:", "", rowOrder)
tf_symbols <- unique(dorothea_hs$tf) 
isTF <- toupper(rowOrder2) %in% tf_symbols
rowOrder_TF <- rowOrder[isTF]
TGIM2<-TGIM[rowOrder_TF,]
TGSM2<-TGSM[rowOrder_TF,]
p3 <- plotTectoryHeatmap(TGIM2,  pal = paletteContinuous(set = "blueYellow"), rowOrder = rowOrder_TF,  varCutOff = 0)
p2 <- plotTectoryHeatmap(TGSM2,  pal = paletteContinuous(set = "horizonExtra"), rowOrder = rowOrder_TF,  varCutOff = 0)
plotPDF(p2, name = "peakHeatmap_GeneScore_IHC_ordered", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
plotPDF(p3, name = "peakHeatmap_Gene_IHC_ordered", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
####OHC
TGSM <- getTectory(ArchRProj = HC_TI, name = "OHC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra") ,returnMat = TRUE)
rowOrder<-rownames(p)
p2 <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra"), rowOrder = rowOrder)
TGIM <- getTectory(ArchRProj = HC_TI, name = "OHC", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTectoryHeatmap(TGIM,  pal = paletteContinuous(set = "blueYellow"),returnMatrix = TRUE)
rowOrder<-intersect(rownames(p3),rowOrder)
TGIM2<-TGIM[rowOrder,]
TGSM2<-TGSM[rowOrder,]
rowOrder2 <- sub("^[^:]+:", "", rowOrder)
tf_symbols <- unique(dorothea_hs$tf) 
isTF <- toupper(rowOrder2) %in% tf_symbols
rowOrder_TF <- rowOrder[isTF]
TGIM2<-TGIM[rowOrder_TF,]
TGSM2<-TGSM[rowOrder_TF,]
p3 <- plotTectoryHeatmap(TGIM2,  pal = paletteContinuous(set = "blueYellow"), rowOrder = rowOrder_TF,  varCutOff = 0)
p2 <- plotTectoryHeatmap(TGSM2,  pal = paletteContinuous(set = "horizonExtra"), rowOrder = rowOrder_TF,  varCutOff = 0)
plotPDF(p2, name = "peakHeatmap_GeneScore_OHC_ordered", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
plotPDF(p3, name = "peakHeatmap_Gene_OHC_ordered", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)

plotPDF(p2, name = "peakHeatmap_GeneScore_OHC_TF", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
write.csv(p2@ht_list[["GeneScoreMatrix"]]@matrix,file="OHC_geneactivity(TF).csv")

#####IHC
TGSM <- getTectory(ArchRProj = HC_TI, name = "IHC", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p <- plotTectoryHeatmap(TGSM,  pal = paletteContinuous(set = "horizonExtra") ,returnMat = TRUE)
rowOrder<-rownames(p)
rowOrder2 <- sub("^[^:]+:", "", rowOrder)
tf_symbols <- unique(dorothea_hs$tf) 
isTF <- toupper(rowOrder2) %in% tf_symbols
rowOrder_TF <- rowOrder[isTF]
TGSM2<-TGSM[rowOrder_TF,]
p2 <- plotTectoryHeatmap(TGSM2,  pal = paletteContinuous(set = "horizonExtra"), rowOrder = rowOrder_TF,  varCutOff = 0)
plotPDF(p2, name = "peakHeatmap_GeneScore_IHC_TF", width = 8, height = 8, ArchRProj = HC_TI, addDOC = FALSE)
write.csv(p2@ht_list[["GeneScoreMatrix"]]@matrix,file="IHC_geneactivity(TF).csv")


#######################################################

#monocle3
library(monocle3)
library(ggplot2)
library(dplyr)
library(ArchR)

cds_path1 <- getMonocleTrajectories(
  ArchRProj = projHeme5, 
  name = "PC/DC_path",
  useGroups = c("PC/DC_progenitor" , "PC"  ,   "DC"),
  principalGroup = "PC/DC_progenitor" ,
  groupBy = "final_type",
  embedding = "UMAP",
  clusterParams = list(k = 50),
  seed = 1
)

cds_path2 <- getMonocleTrajectories(
  ArchRProj = projHeme5, 
  name = "IBC/IPhC_path",
  useGroups = c("IBC/IPhC_progenitor"  ,"IPhC"   ,    "IBC"),
  principalGroup = "IBC/IPhC_progenitor" ,
  groupBy = "final_type",
  embedding = "UMAP",
  clusterParams = list(k = 50),
  seed = 1
)

cds_path3 <- getMonocleTrajectories(
  ArchRProj = projHeme5, 
  name = "HC_path",
  useGroups = c( "HC_progenitor" , "Intermediate_HC"  ,"IHC"  ,   "OHC") ,
  principalGroup = "HC_progenitor",
  groupBy = "final_type",
  embedding = "UMAP",
  clusterParams = list(k = 50),
  seed = 1
)


projHeme5 <- addMonocleTrajectory(
  ArchRProj = projHeme5,
  name = "PC/DC_path",
  useGroups = c("PC/DC_progenitor" , "PC"  ,   "DC"),
  groupBy = "final_type",
  monocleCDS = cds_path1,
  force = TRUE
)

projHeme5 <- addMonocleTrajectory(
  ArchRProj = projHeme5,
  name = "IBC/IPhC_path",
  useGroups = c("IBC/IPhC_progenitor"  ,"IPhC"   ,    "IBC"),
  groupBy = "final_type",
  monocleCDS = cds_path2,
  force = TRUE
)

projHeme5 <- addMonocleTrajectory(
  ArchRProj = projHeme5,
  name = "HC_path",
  useGroups = c( "HC_progenitor" , "Intermediate_HC"  ,"IHC"  ,   "OHC") ,
  groupBy = "final_type",
  monocleCDS = cds_path3,
  force = TRUE
)

p_cds_path1 <- plotTrajectory(projHeme5, trajectory = "PC/DC_path", colorBy = "cellColData", name = "PC/DC path", addArrow = FALSE)




