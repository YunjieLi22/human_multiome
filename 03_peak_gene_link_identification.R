#ArchR

library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(motifmatchr)
library(chromVAR)
library(GenomicRanges)
library(rhdf5)
library(HDF5Array)

ArrowFiles <- createArrowFiles(
 inputFiles = atacfiles,
 sampleNames = names(atacfiles),
 minTSS = 2, #Dont set this too high because you can always increase later
 minFrags = 500, 
 addTileMat = TRUE,
 addGeneScoreMat = TRUE
 )
projMulti1 <- ArchRProject(ArrowFiles = ArrowFiles)
seRNA <- import10xFeatureMatrix(
  input = rnafiles,
  names = names(rnafiles),
  strictMatch = TRUE
)
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE)
}
library(EnsDb.Hsapiens.v86)
length(which(getCellNames(projMulti1) %ni% colnames(seRNA)))
length(which( colnames(seRNA)%ni% getCellNames(projMulti1)))
cells.use <- intersect(na.omit(Cells(all)), projMulti1$cellNames)
projMulti2 <- subsetArchRProject(ArchRProj = projMulti1, cells =cells.use, outputDirectory = "Save-ProjMulti", force = TRUE)
projMulti2 <- addGeneExpressionMatrix(input = projMulti2, seRNA = seRNA, strictMatch = TRUE, force = TRUE)
projMulti2 <- addDoubletScores(projMulti2, force = TRUE)
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)
projMulti2 <- addCombinedDims(projMulti2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
metadata<-metadata[projMulti2@cellColData@rownames,]
cells_archr <- getCellNames(projMulti2)
ix <- match(cells_archr, rownames(metadata))
{v <- metadata$batch
  if (inherits(v, "integer64")) v <- as.numeric(v)  # bit64 -> numeric
  if (is.list(v)) {
    v <- vapply(v, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
  }
  if (is.factor(v)) v <- as.character(v)
  projMulti2<- addCellColData(
    ArchRProj = projMulti2,
    data  = v[ix],
    cells = cells_archr,          
    name  = "batch",
    force = TRUE
  )
}
{v <- metadata$maintype
  if (inherits(v, "integer64")) v <- as.numeric(v)  # bit64 -> numeric
  if (is.list(v)) {
    v <- vapply(v, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
  }
  if (is.factor(v)) v <- as.character(v)
  projMulti2 <- addCellColData(
    ArchRProj = projMulti2,
    data  = v[ix],
    cells = cells_archr,      
    name  = "maintype",
    force = TRUE
  )
}
cM_atac_rna <- confusionMatrix(paste0(projMulti2$maintype), paste0(projMulti2$maintype))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)
library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p_atac_rna
projMulti2 <- addGroupCoverages(ArchRProj = projMulti2, groupBy = "maintype", verbose = FALSE)
projMulti2 <- addReproduciblePeakSet(ArchRProj = projMulti2, groupBy = "maintype")
projMulti2 <- addPeakMatrix(ArchRProj = projMulti2)
projMulti2 <- addPeak2GeneLinks(ArchRProj = projMulti2, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")
p2g <- getPeak2GeneLinks(ArchRProj = projMulti2,returnLoops = FALSE)
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Save-ProjMulti2", overwrite = TRUE, load = TRUE)
#####
md <- getCellColData(projMulti2, select = "maintype") |> as.data.frame()
target_per_type <- 500
md$cell<-rownames(md)
spl <- split(md$cell, md$maintype)
keep_cells <- lapply(spl, function(v){
  if (length(v) <= target_per_type) v else sample(v, target_per_type)
})
keep_cells <- unique(unlist(keep_cells, use.names = FALSE))
all_cells <- rownames(getCellColData(projMulti2))
keep_cells <- intersect(keep_cells, all_cells)
proj_mt <- subsetArchRProject(
  ArchRProj       = projMulti2,
  cells           = keep_cells,
  outputDirectory = "Proj-Downsample-maintype",force = TRUE
)
proj_mt <- addPeak2GeneLinks(
  ArchRProj = proj_mt,
  reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = proj_mt,
  corCutOff = 0.45,
  returnLoops = FALSE
)
palGroup <- c(
  Neuron        = "#17becf",
  Glial_cell    = "#D62728",
  Endothelial_cell = "#FFFF00",
  Mesenchymal_cell = "#C1CDC1",
  Immunocyte    = "#8C564B",
  Epithelial_cell = "#FFA500",
  Melanocyte    = "#E377C2",
  Hair_cell     = "#008B3D",
  Pericyte      = "#0F4C5C",
  Roof_cell     = "#AA40FC",
  DIAPH3_cell   = "#B5BD61"
)
palGroup<-palGroup[desired_order]
proj_mt@cellColData$maintype <- factor(
  as.character(proj_mt@cellColData$maintype),
  levels = desired_order
)
proj_mt@cellColData$maintype <- droplevels(proj_mt@cellColData$maintype)

p <- plotPeak2GeneHeatmap.distal (ArchRProj = proj_mt, groupBy = "maintype",nPlot = 60000,palGroup = palGroup)
plotPDF(p, name = "Plot-heatmap_p2g", width = 12, height = 10, ArchRProj = proj_mt, addDOC = FALSE)
p <- plotPeak2GeneHeatmap.distal (ArchRProj = proj_mt, groupBy = "maintype",returnMatrices = TRUE,nPlot = 60000)
pairs$kmeans<-as.character(list)
write.csv(pairs,file="P2Gpairs.csv")
proj_mt <- saveArchRProject(ArchRProj = proj_mt, outputDirectory = "Proj-Downsample-maintype/", overwrite = TRUE, load = TRUE)

#####roof, floor
metadata<-REH@meta.data
keep_cells<-intersect( rownames(metadata),projMulti2$cellNames)
proj_epi<- subsetArchRProject(
  ArchRProj       = projMulti2,
  cells           = keep_cells,
  outputDirectory = "Proj-Downsample-pei",force = TRUE
)
metadata<-metadata[proj_epi@cellColData@rownames,]
cells_archr <- getCellNames(proj_epi)
ix <- match(cells_archr, rownames(metadata))
{v <- metadata$subtype
  if (inherits(v, "integer64")) v <- as.numeric(v)  # bit64 -> numeric
  if (is.list(v)) {
    v <- vapply(v, function(y) if (length(y)) as.character(y[[1]]) else NA_character_, character(1))
  }
  if (is.factor(v)) v <- as.character(v)
  proj_epi<- addCellColData(
    ArchRProj = proj_epi,
    data  = v[ix],
    cells = cells_archr,          
    name  = "subtype",
    force = TRUE
  )
}
md <- getCellColData(proj_epi, select = "subtype") |> as.data.frame()
target_per_type <- 300
md$cell<-rownames(md)
spl <- split(md$cell, md$subtype)
keep_cells <- lapply(spl, function(v){
  if (length(v) <= target_per_type) v else sample(v, target_per_type)
})
keep_cells <- unique(unlist(keep_cells, use.names = FALSE))
all_cells <- rownames(getCellColData(proj_epi))
keep_cells <- intersect(keep_cells, all_cells)
proj_mt2 <- subsetArchRProject(
  ArchRProj       = proj_epi,
  cells           = keep_cells,
  outputDirectory = "Proj-Downsample-subtype",force = TRUE
)
proj_mt2 <- addPeak2GeneLinks(
  ArchRProj = proj_mt2,
  reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = proj_mt2,
  corCutOff = 0.45,
  returnLoops = FALSE
)
p <- plotPeak2GeneHeatmap.distal (ArchRProj = proj_mt2, groupBy = "subtype",palGroup = palGroup)
plotPDF(p, name = "Plot-heatmap_p2g_subtype", width = 12, height = 10, ArchRProj = proj_mt2, addDOC = FALSE)
palGroup <- c(
  "MEC"          = "#b6bd61",
  "OSC"          = "#ffbb78",
  "HC"           = "#289e67",
  "SDC"          = "#c4b0d5",
  "LSDC"         = "#e377c3",
  "KOC"          = "#ab40fc",
  "ISC"      = "#d72929",
  "Cycling_RFC"  = "#1f76b3",
  "FF"           = "#ff7f0f",
  "MSC2"         = "#afc7e8",
  "MSC1"         = "#15becf",
  "CyclingFF"    = "#28a3ff",
  "RF"           = "#99df8a",
  "SPC"         = "#c49c94",
  "RMC"         = "#ff9896"
)
palGroup <- c(
  "MEC"            = "#C2C878",
  "OSC"            = "#FFC88E",
  "HC"             = "#13AC7C",
  "SDC"            = "#D0BEDC",
  "LSDC"           = "#EE8ECC",
  "KOC"            = "#C05BF9",
  "ISC"  = "#E54238",
  "Cycling_RFC"    = "#1DB2FD",
  "FF"             = "#FF9422",
  "MSC2"           = "#BCD2EB",
  "MSC1"           = "#00C8D8",
  "CyclingFF"      = "#1A89BF",
  "RF"             = "#A4E39F",
  "SPC"     = "#D1ADA6",
  "RMC"     = "#FFABA5"
)
keep_cells <- rownames(subset(metadata, !(subtype %in% c("Junct_cell1","Junct_cell2"))))
keep_cells<-intersect(keep_cells,proj_mt2$cellNames)
proj_mt_noJunct <- subsetArchRProject(
  ArchRProj        = proj_epi,
  cells            = keep_cells,
  outputDirectory  = "Save-proj_epi_noJunct",
  force            = TRUE
)
proj_mt_noJunct <- addPeak2GeneLinks(
  ArchRProj = proj_mt_noJunct,
  reducedDims = "LSI_Combined",useMatrix = "GeneExpressionMatrix"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = proj_mt_noJunct,
  corCutOff = 0.45,
  returnLoops = FALSE
)
ps <- metadata(p2g)$peakSet   # GRanges of peaks
gs <- metadata(p2g)$geneSet   # GRanges of genes
peak_all <- paste0(as.character(seqnames(ps)), ":", start(ps), "-", end(ps))
gene_all <- mcols(gs)$name
pairs <- data.table(
  peakName = peak_all[p2g$idxATAC],
  geneName = gene_all[p2g$idxRNA])

p <- plotPeak2GeneHeatmap.distal (ArchRProj = proj_mt_noJunct, groupBy = "subtype",palGroup = palGroup)
plotPDF(p, name = "Plot-heatmap_p2g", width = 12, height = 10, ArchRProj = proj_mt_noJunct, addDOC = FALSE)
p <- plotPeak2GeneHeatmap.distal (ArchRProj = proj_mt_noJunct, groupBy = "maintype",returnMatrices = TRUE)
list<-p@listData[["RNA"]]@listData[["kmeansId"]]
pairs$kmeans<-as.character(list)
write.csv(pairs,file = "EPI_P2G_k.csv")



#Cicero
library(cicero)
library(monocle3)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

m.cds <- as.cell_data_set(x = snRNA)
cicero_cds <- make_cicero_cds(m.cds, reduced_coordinates = reducedDims(m.cds)$UMAP)

genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome_df <- data.frame(
  V1 = names(genome),
  V2 = as.numeric(genome),
  stringsAsFactors = FALSE)
sample_genome <- genome_df[1:24,]

conns <- run_cicero(cicero_cds, sample_genome) 
head(conns)

ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)

Links(snRNA) <- links

CoveragePlot(snRNA, region = "chrX-133980149-133990649")#example

