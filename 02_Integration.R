
#integration of scRNA-seq & snRNA-seq

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidydr)

rna_only_obj <- CreateSeuratObject(
  counts = GetAssayData(all_peakcalled, assay = "RNA", slot = "counts"),
  meta.data = all_peakcalled@meta.data  # 保留所有metadata
)

if ("SCT" %in% Assays(all_peakcalled)) {
  rna_only_obj[["SCT"]] <- all_peakcalled[["SCT"]]
}

rna_only_obj@reductions <- all_peakcalled@reductions  # 全部保留
snRNA_all = rna_only_obj
snRNA_all$tag = "snRNA"
saveRDS(snRNA_all, file = "E:/human multiome/0626画图数据/snRNA_only.rds")

scRNA = readRDS("E:/human multiome/0626画图数据/scRNA_new.rds")


scRNA <- JoinLayers(scRNA, assay = "RNA", layer = "counts")
snRNA <- JoinLayers(snRNA, assay = "RNA", layer = "counts")
snRNA <- JoinLayers(snRNA, assay = "RNA", layer = "data")
scRNA <- JoinLayers(scRNA, assay = "RNA", layer = "data")


preprocess <- function(obj) {
  obj <- NormalizeData(obj, layer = "counts") 
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- ScaleData(obj)  
  obj <- RunPCA(obj)
  return(obj)
}

snRNA <- preprocess(snRNA)
scRNA <- preprocess(scRNA)

features <- SelectIntegrationFeatures(object.list = list(snRNA, scRNA))

anchors <- FindIntegrationAnchors(
  object.list = list(scRNA, snRNA),
  anchor.features = features, 
  reduction = "rpca",     
  dims = 1:30
)

combined <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  new.assay.name = "integrated",
  features.to.integrate = rownames(combined)  # 明确指定所有基因
)

combined <- FindVariableFeatures(combined, assay = "integrated")
combined <- ScaleData(combined, assay = "integrated")
combined <- RunPCA(
  combined,
  assay = "integrated",
  features = VariableFeatures(combined),
  verbose = FALSE
)
combined <- RunUMAP(combined, dims = 1:30)


combined$age = combined$batch

combined$age[which(combined$age == "E15")] = "3"
combined$age[which(combined$age == "E16")] = "3"
combined$age[which(combined$age == "E14")] = "2"
combined$age[which(combined$age == "E13")] = "2"
combined$age[which(combined$age == "E11")] = "1"
combined$age[which(combined$age == "E12")] = "1"

#batch removal
library(harmony)
combined <- RunHarmony(
  combined,
  group.by.vars = "age",
  assay.use = "RNA",       
  max.iter.harmony = 30    
)
combined = RunUMAP(combined, reduction = "harmony", reduction.name = "umap_harmony", dims = 1:30)



########################################################################


#integration of subtypes(glial cells, endothelial cells, immunocyte)


#subset scRNA-seq
sc_all = readRDS("E:/human multiome/0626画图数据/scRNA_new.rds")
glial <- subset(sc_all, cells = colnames(sc_all)[which(sc_all$maintype == "Glial")])
endothelial <- subset(sc_all, cells = colnames(sc_all)[which(sc_all$maintype == "Endothelial")])
neuron <- subset(sc_all, cells = colnames(sc_all)[which(sc_all$maintype == "Neuron")])

glial$tag <- "scRNA"
endothelial$tag <- "scRNA"
neuron$tag <- "scRNA"


#snRNA-seq
all_peakcalled <- readRDS("E:/human multiome/0626画图数据/multiome_final.rds")
DefaultAssay(all_peakcalled) <- "RNA"
Idents(all_peakcalled) <- all_peakcalled$maintype
sn_glial <- subset(all_peakcalled, cells = colnames(all_peakcalled)[which(all_peakcalled$maintype == "Glial_cell")])
sn_endothelial <- subset(all_peakcalled, cells = colnames(all_peakcalled)[which(all_peakcalled$maintype == "Endothelial_cell")])
sn_neuron <- subset(all_peakcalled, cells = colnames(all_peakcalled)[which(all_peakcalled$maintype == "Neuron")])

sn_glial$tag <- "snRNA"
sn_endothelial$tag<- "snRNA"
sn_neuron$tag <- "snRNA"

save(glial, endothelial, neuron, sn_endothelial, sn_glial, sn_neuron, file = "human multiome/data/glial_endo_neuron_subset.RData")

#common genes intersection(appled for glial cells,endothelial cells and immunocyte )
scRNA_subset = glial
snRNA_subset = sn_glial
common_genes <- intersect(rownames(scRNA_subset), rownames(snRNA_subset))
scRNA_subset <- scRNA_subset[common_genes, ]
snRNA_subset <- snRNA_subset[common_genes, ]

#merge
combined <- merge(scRNA_subset, snRNA_subset)

#batch effect removal
##harmony
library(harmony)
combined <- NormalizeData(combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
combined <- RunHarmony(combined, group.by.vars = "tag")
combined <- RunHarmony(combined, group.by.vars = "age")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)

combined <- RunHarmony(
object = combined,
group.by.vars = "batch",       
reduction = "pca",
reduction.save = "harmony_age",
verbose = TRUE
)
combined <- RunHarmony(
object = combined,
group.by.vars = "tag",    
reduction = "harmony_age",    
reduction.save = "harmony_final",
verbose = TRUE
)

DimPlot(combined, group.by = "tag") 
DimPlot(combined, group.by = "batch") 

ElbowPlot(combined, ndims = 50)  
combined <- FindNeighbors(
  combined, 
  reduction = "harmony_final",          
  dims = 1:20                 
)
combined <- FindClusters(
  combined,
  resolution = 0.1,          
  algorithm = 1,             
  random.seed = 42
)
combined <- RunUMAP(
  combined, 
  reduction = "harmony_final",         
  dims = 1:20,               
  seed.use = 42
)

combined <- FindClusters(
  combined,
  resolution = 0.05,          
  algorithm = 1,             
  random.seed = 42
)

combined <- RunUMAP(
combine,
reduction = "harmony",          # 若用Harmony则改为 "harmony"
dims = 1:20,                # 与FindNeighbors一致
seed.use = 42
)
combined_endo = combined


# 
DimPlot(combined, group.by = "seurat_clusters", label = TRUE)
DimPlot(combined, group.by = "new_batch", shuffle = TRUE)  # 理想情况应混合均匀

combined <- JoinLayers(combined, assay = "RNA", layer = "data")



#find markers
markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,            
  min.pct = 0.25,             
  logfc.threshold = 0.5       
)

top20 <- markers %>%  group_by(cluster) %>% slice_head(n = 20)

t <- markers
t$p_val_adj[t$p_val_adj == 0] <- min(t$p_val_adj[t$p_val_adj > 0])
t$combined_score <- abs(t$avg_log2FC) * -log10(t$p_val_adj) * t$pct.1

t <- t %>%
  group_by(cluster) %>%                
  arrange(desc(combined_score), .by_group = TRUE)

top20 <- t %>% group_by(cluster) %>% slice_max(n = 20, order_by = combined_score)

glial_markers = t
endo_markers  = t
neuron_markers  = t




#age re-annotated
combined$age <- combined$batch
combined$age[combined$age %in% c("E12", "E12_2")] = "E10"
combined$age[combined$age %in% c("E13_1a", "E13_1b")] = "E13"
combined$age[combined$age %in% c("E15_1a", "E15_1b")] = "E15"

combined$new_age = combined$age
combined$new_age[combined$new_age %in% c("E11", "E12") ] = "1"
combined$new_age[combined$new_age %in% c("E13", "E14") ] = "2"
combined$new_age[combined$new_age %in% c("E15", "E16") ] = "3"
saveRDS(combined, file = "0626画图数据/combine_endo.rds")
saveRDS(combined, file = "0626画图数据/combine_glial.rds")


combine_glial = combined
combined_endo = combined
combine_neuron = combined

DimPlot(combined, group.by = "age", label = TRUE)

combined$age_highlight <- ifelse(
combined$age %in% c("E12", "E13", "E15"),  
as.character(combined$age),
"Other"
)

DimPlot(combined, group.by = "age_highlight") +
scale_color_manual(
values = c("E12" = "orange", "E13" = "skyblue", "E15" = "green3", "Other" = "gray90"),
breaks = c("E12", "E13", "E15") 
)

combined$age_highlight2 <- ifelse(
combined$age %in% c("E11", "E14", "E16"),  
as.character(combined$age),
"Other"
)

DimPlot(combined, group.by = "age_highlight2") +
scale_color_manual(
values = c("E11" = "#FE8A6E", "E14" = "#7AD2E1", "E16" = "#1C9E76", "Other" = "gray90"),
breaks = c("E11", "E14", "E16")  
)

#top20 markers heatmap
cols6 <- c("#2269AD", "#6FAED0", "#D1E5F0", "#F9CBB1", "#DD6F59", "#9D1126")
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("#2269AD",  "#D1E5F0","white", "#F9CBB1",  "#9D1126"))

p <- DotPlot(combined, features = top20$gene) +theme(axis.text.x = element_text(angle = 90))+theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),  axis.title.y = element_blank())+scale_color_gradientn(colours = cols6)+theme(panel.border = element_rect(color = "black"), panel.spacing = unit(1, "mm"))

df <- p$data
df_sub <- df %>%
dplyr::select(features.plot, id, avg.exp.scaled) 
exp_mat <- df_sub %>%
  pivot_wider(
    names_from = id,
    values_from = avg.exp.scaled,
    values_fill = NA 
  ) %>%
  as.data.frame()
row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()
exp_mat <- t(exp_mat) 

Heatmap(exp_mat,
heatmap_legend_param=list(title="avg_expression"),
col= col_fun ,
cluster_rows = F,
cluster_columns = F,
row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 6), border = "black")










