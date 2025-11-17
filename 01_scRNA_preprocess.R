#scRNA-seq preprocessing

##PCW12----------------------
data <- Read10X(data.dir = "/home/jcfan/human_brain/arcresult/humanP12wc1/outs/filtered_feature_bc_matrix/")
E12_1<-CreateSeuratObject(counts = data,project = "E12",min.cells = 3,min.features = 200)

data <- Read10X(data.dir = "/home/jcfan/human_brain/arcresult/humanP12wc2/outs/filtered_feature_bc_matrix/")
E12_2<-CreateSeuratObject(counts = data,project = "E12_2",min.cells = 3,min.features = 200)
colnames(E12_1)<-paste0("E12_1",colnames(E12_1))
colnames(E12_2)<-paste0("E12_2",colnames(E12_2))

E12<-merge(E12_1,E12_2)

E12[["percent.mt"]] <- PercentageFeatureSet(E12, pattern = "^MT-")
VlnPlot(E12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

E12 <- subset(
  E12,
  subset = nFeature_RNA >800 & 
    percent.mt < 5
)

E12 <- NormalizeData(E12)
E12 <- FindVariableFeatures(E12, nfeatures = 4000)
all.genes <- rownames(E12)
E12 <- ScaleData(E12,features = all.genes)
E12 <- RunPCA(E12, npcs = 50)

E12 <- FindNeighbors(E12, dims = 1:30)
E12 <- FindClusters(E12, resolution = 0.8)
E12 <- RunUMAP(E12, dims = 1:30)
DimPlot(E12, group.by = c("orig.ident","seurat_clusters"),label = T)
DimPlot(E12, group.by = "orig.ident")
saveRDS(E12,file = "E12.rds")

#mapping from snRNA
E11_rename@active.assay<-"RNA"
E11_rename <- NormalizeData(E11_rename) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
E12<-JoinLayers(E12)
anchors <- FindTransferAnchors(
  reference = E11_rename,
  query = E12,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = E11_rename$maintype,
  dims = 1:30
)

#addMetaData
E12 <- AddMetaData(
  object = E12,
  metadata = predictions
)

#visulization
#UMAP
DimPlot(E12, group.by = c("predicted.id"), label = TRUE) 
FeaturePlot(E12,features = "nFeature_RNA")
pancreas.ref<-RunUMAP(E11_rename, dims = 1:20, n.neighbors = 50,return.model=TRUE)
pancreas.query <- MapQuery(anchorset = anchors, reference = pancreas.ref, query = E12,
                           refdata = list(celltype = "maintype"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "maintype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query,reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
E12<-pancreas.query
saveRDS(E12,file = "E12.rds")
FeaturePlot(E12,features = c("LGR5","SOX2"),label = T)
FeaturePlot(E12,features = c("EPCAM","OTOG","GATA3","USH1C"),label = T,reduction = "ref.umap")
FeaturePlot(E12,features = c("OC90","LMX1A","OTX2","PAX2"),label = T,reduction = "ref.umap")

DimPlot(E12,label = T,reduction = "ref.umap")
E12@active.ident<-factor(E12$predicted.celltype)
#####reh分类
NEWREH@active.ident<-factor(NEWREH$batch)
reh11<-subset(NEWREH,idents=c("E11"))
DimPlot(reh11,reduction = "umap.cca")
ROOF_EPI<-subset(E12,idents = c("Roof_cell","Epithelial_cell"))
DimPlot(ROOF_EPI)
anchors <- FindTransferAnchors(
  reference = reh11,
  query = ROOF_EPI,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = reh11$subtype,
  dims = 1:30
)
ROOF_EPI <- AddMetaData(
  object = ROOF_EPI,
  metadata = predictions
)
DimPlot(ROOF_EPI)
ROOF_EPI$maintype<-ROOF_EPI@active.ident
ROOF_EPI$subtype<-ROOF_EPI$predicted.id
DimPlot(ROOF_EPI,group.by = "subtype")
ROOF_EPI <- NormalizeData(ROOF_EPI)
ROOF_EPI <- FindVariableFeatures(ROOF_EPI, nfeatures = 2000)
all.genes <- rownames(ROOF_EPI)
ROOF_EPI <- ScaleData(ROOF_EPI,features = all.genes)
ROOF_EPI <- RunPCA(ROOF_EPI, npcs = 50)
ROOF_EPI<-RunUMAP(ROOF_EPI,dims = 1:30)
DimPlot(ROOF_EPI,group.by = "subtype",label = T)
saveRDS(ROOF_EPI,file = "../RNA12w/REH12.rds")
ROOF_EPI@active.ident<-factor(ROOF_EPI$subtype)
FeaturePlot(ROOF_EPI,features = c("LGR5","SOX2","FGFR3","KLHL4"),label = T)
FeaturePlot(ROOF_EPI,features = c("ATOH1","PCP4","POU4F3","CCER2"),label = T)
FeaturePlot(ROOF_EPI,features = c("CDH20","TECTB","EDA"),label = T)
##SDC
SDC11<-subset(ROOF_EPI,idents=c("SDC"))
LSDC11<-subset(ROOF_EPI,idents=c("LSDC"))
DimPlot(SDC11)
DimPlot(ROOF_EPI)


##PCW13----------------------

data <- Read10X(data.dir = "/home/jcfan/human_brain/arcresult/humanP13wc1/outs/filtered_feature_bc_matrix/")
E13_1a<-CreateSeuratObject(counts = data,project = "E13_1a",min.cells = 3,min.features = 200)

data <- Read10X(data.dir = "/home/jcfan/human_brain/arcresult/humanP13wc2/outs/filtered_feature_bc_matrix/")
E13_1b<-CreateSeuratObject(counts = data,project = "E13_1b",min.cells = 3,min.features = 200)
colnames(E13_1a)<-paste0("E13_1a",colnames(E13_1a))
colnames(E13_1b)<-paste0("E13_1b",colnames(E13_1b))

E13<-merge(E13_1a,E13_1b)

E13[["percent.mt"]] <- PercentageFeatureSet(E13, pattern = "^MT-")
VlnPlot(E13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

E13 <- subset(
  E13,
  subset = nFeature_RNA >800 & 
    percent.mt < 5
)

E13 <- NormalizeData(E13)
E13 <- FindVariableFeatures(E13, nfeatures = 4000)
all.genes <- rownames(E13)
E13 <- ScaleData(E13,features = all.genes)
E13 <- RunPCA(E13, npcs = 50)

E13 <- FindNeighbors(E13, dims = 1:30)
E13 <- FindClusters(E13, resolution = 0.8)
E13 <- RunUMAP(E13, dims = 1:30)
DimPlot(E13, group.by = c("orig.ident","seurat_clusters"),label = T)
DimPlot(E13, group.by = "orig.ident")
save(E13,file = "E13.Rdata")

rm(E11_rename)
E14_rename@active.assay<-"RNA"
E14_rename <- NormalizeData(E14_rename) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
E13<-JoinLayers(E13)
anchors <- FindTransferAnchors(
  reference = E14_rename,
  query = E13,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = E14_rename$maintype,
  dims = 1:30
)

E13 <- AddMetaData(
  object = E13,
  metadata = predictions
)

DimPlot(E13, group.by = "predicted.id", label = TRUE) 
FeaturePlot(E13,features = "nFeature_RNA")

pancreas.ref<-RunUMAP(E14_rename, dims = 1:20, n.neighbors = 50,return.model=TRUE)
pancreas.query <- MapQuery(anchorset = anchors, reference = pancreas.ref, query = E13,
                           refdata = list(celltype = "maintype"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "maintype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query,reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
E13<-pancreas.query
saveRDS(E13,file = "E13.rds")
DimPlot(E13, group.by = c("predicted.id","seurat_clusters"),label = T)
E13_filtered<-subset(E13,idents=c("0","2","16","15","8","33","1","14","6","4","10"),invert=T)
DimPlot(E13_filtered, group.by = c("predicted.id","seurat_clusters"),label = T,reduction = "ref.umap")
VlnPlot(E13_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save(E13_filtered,file = "E13_filtered.Rdata")

E13_filtered@active.ident <- factor(E13_filtered$predicted.celltype)
E13_filtered<-FindVariableFeatures(E13_filtered,nfeatures = 2000)
all.genes <- rownames(E13_filtered)
E13_filtered <- ScaleData(E13_filtered, features = all.genes)
E13_filtered <- RunPCA(E13_filtered, features = VariableFeatures(object = E13_filtered))
E13_filtered <- RunUMAP(E13_filtered, dims = 1:10)
DimPlot(E13_filtered,label = T,group.by = "CHOIR_parent_clusters")
E13_filtered <- FindNeighbors(E13_filtered, dims = 1:10)
E13_filtered <- FindClusters(E13_filtered, resolution = 0.5)
DimPlot(E13_filtered,group.by = c("RNA_snn_res.0.5","predicted.id"),label = T)
FeaturePlot(E13_filtered)
E13_filtered <- NormalizeData(E13_filtered)
E13_filtered <- FindVariableFeatures(E13_filtered, nfeatures = 2000)
all.genes <- rownames(E13_filtered)
E13_filtered <- ScaleData(E13_filtered,features = all.genes)
E13_filtered <- RunPCA(E13_filtered, npcs = 50)
E13_filtered <- RunUMAP(E13_filtered, dims = 1:10)

NEWREH@active.ident <- factor(NEWREH$batch)
reh14 <- subset(NEWREH, idents = c("E14"))
DimPlot(reh14, reduction = "umap.cca")

#subtype classification（以E13_filtered为例）
ROOF_EPI <- subset(E13_filtered, idents = c("Roof_cell", "Epithelial_cell","Hair_cell"))
DimPlot(ROOF_EPI)

#integration and prediction
anchors <- FindTransferAnchors(
  reference = reh14,
  query = ROOF_EPI,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = reh14$subtype,
  dims = 1:30
)
ROOF_EPI <- AddMetaData(
  object = ROOF_EPI,
  metadata = predictions
)

ROOF_EPI$maintype <- ROOF_EPI@active.ident
ROOF_EPI$subtype <- ROOF_EPI$predicted.id
DimPlot(ROOF_EPI, group.by = "subtype", label = TRUE)

ROOF_EPI <- NormalizeData(ROOF_EPI)
ROOF_EPI <- FindVariableFeatures(ROOF_EPI, nfeatures = 2000)
all.genes <- rownames(ROOF_EPI)
ROOF_EPI <- ScaleData(ROOF_EPI, features = all.genes)
ROOF_EPI <- RunPCA(ROOF_EPI, npcs = 50)
ROOF_EPI <- RunUMAP(ROOF_EPI, dims = 1:30)
ROOF_EPI <- CHOIR(ROOF_EPI, 
                  n_cores = 16)
ROOF_EPI <- runCHOIRumap(ROOF_EPI,
                         reduction = "P0_reduction")
plotCHOIR(ROOF_EPI,group_by = "subtype")
plotCHOIR(ROOF_EPI)
ROOF_EPI <- FindNeighbors(ROOF_EPI, dims = 1:10)
ROOF_EPI <- FindClusters(ROOF_EPI, resolution = 0.5)
DimPlot(ROOF_EPI, label = TRUE,group.by = "subtype")


DimPlot(ROOF_EPI, group.by = c("maintype","subtype"), label = TRUE)
saveRDS(ROOF_EPI, file = "../RNA13w/REH13_filtered.rds")  # 修改文件名避免覆盖


##PCW15---------------------
data <- Read10X(data.dir = "/home/jcfan/human_brain/arcresult/humanP15wc1/outs/filtered_feature_bc_matrix/")
E15_1a<-CreateSeuratObject(counts = data,project = "E15_1a",min.cells = 3,min.features = 200)

data <- Read10X(data.dir = "/home/jcfan/human_brain/arcresult/humanP15wc2/outs/filtered_feature_bc_matrix/")
E15_1b<-CreateSeuratObject(counts = data,project = "E15_1b",min.cells = 3,min.features = 200)
colnames(E15_1a)<-paste0("E15_1a",colnames(E15_1a))
colnames(E15_1b)<-paste0("E15_1b",colnames(E15_1b))

E15<-merge(E15_1a,E15_1b)

E15[["percent.mt"]] <- PercentageFeatureSet(E15, pattern = "^MT-")
VlnPlot(E15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

E15 <- subset(
  E15,
  subset = nFeature_RNA >800 & 
    percent.mt < 5
)

E15 <- NormalizeData(E15)
E15 <- FindVariableFeatures(E15, nfeatures = 4000)
all.genes <- rownames(E15)
E15 <- ScaleData(E15,features = all.genes)
E15 <- RunPCA(E15, npcs = 50)

E15 <- FindNeighbors(E15, dims = 1:30)
E15 <- FindClusters(E15, resolution = 0.8)
E15 <- RunUMAP(E15, dims = 1:30)
DimPlot(E15, group.by = c("orig.ident","seurat_clusters"),label = T)
DimPlot(E15, group.by = "orig.ident")
save(E15,file = "E15.Rdata")

rm(E14_rename)
E16_rename@active.assay<-"RNA"
E16_rename <- NormalizeData(E16_rename) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
E15<-JoinLayers(E15)
anchors <- FindTransferAnchors(
  reference = E16_rename,
  query = E15,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(
  anchorset = anchors,
  refdata = E16_rename$maintype,
  dims = 1:30
)

E15 <- AddMetaData(
  object = E15,
  metadata = predictions
)


DimPlot(E15, group.by = c("predicted.id"), label = TRUE,reduction = "ref.umap",repel = T) 
FeaturePlot(E15,features = "nFeature_RNA")

pancreas.ref<-RunUMAP(E16_rename, dims = 1:20, n.neighbors = 50,return.model=TRUE)
pancreas.query <- MapQuery(anchorset = anchors, reference = pancreas.ref, query = E15,
                           refdata = list(celltype = "maintype"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "maintype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query,reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2
E15<-pancreas.query
saveRDS(E15,file = "E15.rds")
DimPlot(E15, group.by = c("maintype","seurat_clusters"),label = T)
E15_filtered<-subset(cc,idents=c("0","2","16","15","8","33","1","14","6","4","10"),invert=T)
DimPlot(E15_filtered, group.by = c("predicted.id","seurat_clusters"),label = T,reduction = "ref.umap")
VlnPlot(E15_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save(E15_filtered,file = "E15_filtered.Rdata")

FeaturePlot(E15,c("ATOH1", "PCP4", "POU4F3","CCER2"))
FeaturePlot(E13_filtered,c("ATOH1", "PCP4", "POU4F3","CCER2"))
FeaturePlot(E12,c("ATOH1", "PCP4", "POU4F3","CCER2"))

