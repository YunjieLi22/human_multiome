library(Seurat)
library(cowplot)
#library(feather)
library(dplyr)
library(Matrix)
library(Hmisc)
library(tidyr)
library(stringr)

library(matrixStats)
library(ggplot2)
library(igraph)

seurat_Human_sc=readRDS('seurat_Macaca_sc_raw_WuZhongJian.rds.gz')
seurat_Human_SN=readRDS('seurat_Macaca_SN_raw_WuZhongJian.rds.gz')
seurat_Mouse=readRDS('seurat_mouse_raw_WuZhongJian.rds.gz')
Idents(seurat_Human_sc) <- seurat_Human_sc$orig_new_type1
Idents(seurat_Human_SN) <- seurat_Human_SN$orig_new_type1
Idents(seurat_Mouse) <- seurat_Mouse$orig_new_type1

seurat_Human_sc=subset(seurat_Human_sc,idents=c('GC','MC','IPhC/IBC','RMC','ISC/IdC','KOC','RtC','SGN','MgC','PC/DC','HeC/OSC','RMC2','IC','HC','M'))
seurat_Human_SN=subset(seurat_Human_SN,idents=c('GC','MC','IPhC/IBC','ISC/IdC','KOC','SGN','MgC','HeC/OSC','RMC1','RMC2','IC','HC','M'))
seurat_Mouse=subset(seurat_Mouse,idents=c('GC','M','IC','HeC/OSC','KOC','MC','RMC','IPhC/IBC','ISC','RtC','SGN','PC/DC','HC','MgC'))

seurat_Mouse$orig.ident_species <- "mouse"
seurat_Human_SN$orig.ident_species <- "human_SN"
seurat_Human_sc$orig.ident_species <- "human_sc"

seurat_Mouse$orig.ident_species=paste(seurat_Mouse$orig.ident_species,seurat_Mouse$orig_orig.ident,sep='__')
seurat_Human_SN$orig.ident_species=paste(seurat_Human_SN$orig.ident_species,seurat_Human_SN$orig_orig.ident,sep='__')
seurat_Human_sc$orig.ident_species=paste(seurat_Human_sc$orig.ident_species,seurat_Human_sc$orig_orig.ident,sep='__')

all.data <- merge(x = seurat_Human_sc, y = list(seurat_Human_SN,seurat_Mouse), add.cell.ids = c('human_sc', "human_SN","Mouse"))

combined.list <- SplitObject(all.data, split.by = "orig.ident_species")

for (i in 1:length(combined.list)) {
  print(paste(i,'SCT'))
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = TRUE, variable.features.n = 3500)
}

combined.features <- SelectIntegrationFeatures(object.list = combined.list, 
    nfeatures = length(Reduce(union, lapply(combined.list,VariableFeatures))))
combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = combined.features)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT",
                                           anchor.features = combined.features, verbose = TRUE)
sample.combined <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                 verbose = TRUE)

sample.combined <- RunPCA(sample.combined)
ElbowPlot(sample.combined, ndims = 40)
sample.combined <- FindNeighbors(sample.combined, reduction = "pca", 
      dims = 1:40, nn.eps = 0,k.param=25)
sample.combined <- FindClusters(sample.combined, resolution = 1) 
sample.combined$seurat_clusters.new <- as.integer(sample.combined$seurat_clusters)
sample.combined <- RunUMAP(sample.combined, reduction = "pca", dims = 1:40)

saveRDS(sample.combined,'seurat_mouse_human.rds.gz')

