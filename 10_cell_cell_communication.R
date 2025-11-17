
#niche_based cell-cell communication

library(ggsignif)
library(ggpubr)
library(cowplot)
library(dplyr)
library(NMF)
library(ggalluvial)
library(CellChat)
library(pheatmap)
library(qs)
library(patchwork)
library(presto)
library(stringr)
library(future)
plan(strategy = "multisession")

snRNA = readRDS("E:/human multiome/0626画图数据/multiome_final.rds")

#Floor niche
target_groups = c( "Endothelial_0" ,  "Endothelial_1",   "Endothelial_2" ,"Endothelial_4" ,  "Immunocyte_1",    "Immunocyte_2" ,   "Immunocyte_3" , "Mesenchymal_0" ,"Mesenchymal_1" ,  "Mesenchymal_3" ,"PC/DC_domain", "IBC/IPhC_domain",  "HC_domain",  "OSC"  , "KOC", "ISC/IdC" )
PCW16_down = subset(snRNA, subset = batch == "E16" & cellchat_type %in% target_groups)
PCW11_down = subset(snRNA, subset = batch == "E11" & cellchat_type %in% target_groups)


#SV niche
target_groups1 = c("Endothelial_0" ,  "Endothelial_1",   "Endothelial_2" ,"Endothelial_4" ,   "Immunocyte_1",    "Immunocyte_2" ,   "Immunocyte_3" , "Mesenchymal_0" ,"Mesenchymal_1" ,  "Mesenchymal_3",  "Mesenchymal_7", "MSC1" ,"MSC2" ,   "Melanocyte", "Cycling_RFC", "MEC" )
PCW16_roof = subset(snRNA, subset = batch == "E16" & cellchat_type %in% target_groups1)
PCW11_roof = subset(snRNA, subset = batch == "E11" & cellchat_type %in% target_groups1)


#RMC niche
target_groups2 = c( "Immunocyte_1",    "Immunocyte_2" ,   "Immunocyte_3" ,  "Mesenchymal_5", "RMC","SPC", "RF",  "Cycling_RFC" )
PCW16_rmc = subset(snRNA, subset = batch == "E16" & cellchat_type %in% target_groups2)
PCW11_rmc = subset(snRNA, subset = batch == "E11" & cellchat_type %in% target_groups2)


#SGN niche
target_groups3 = c("Endothelial_0" ,  "Endothelial_1",   "Endothelial_2" ,"Endothelial_4" , "Immunocyte_1",    "Immunocyte_2" ,   "Immunocyte_3" , "Glial_0" ,"Glial_1" ,"Glial_2"  ,"Glial_3"  , "Glial_4", "Mesenchymal_0" ,"Mesenchymal_1" ,  "Mesenchymal_3"   , "Mesenchymal_4" , "Mesenchymal_8", "SGN"  )
PCW16_sgn = subset(snRNA, subset = batch == "E16" & cellchat_type %in% target_groups3)
PCW11_sgn = subset(snRNA, subset = batch == "E11" & cellchat_type %in% target_groups3)



#Differential Cell Communication Analysis
##Floor niche
custom_order = c("Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_4",   "Immunocyte_1",  "Immunocyte_2",  "Immunocyte_3", 
                   "Mesenchymal_0" ,"Mesenchymal_1", "Mesenchymal_3","FF","HC_domain", "PC/DC_domain", "IBC/IPhC_domain", "ISC/IdC", "KOC", "OSC") 
  
Idents(PCW11_down) = factor(PCW11_down$cellchat_type, levels = custom_order)
Idents(PCW16_down) = factor(PCW16_down$cellchat_type, levels = custom_order)

##SGN niche
custom_order = c("Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_4", "Glial_0"  ,     "Glial_1" ,     
                 "Glial_2" ,      "Glial_3" ,      "Glial_4"   ,    "Immunocyte_1",  "Immunocyte_2",  "Immunocyte_3",
                 "Mesenchymal_4" ,"Mesenchymal_8","Mesenchymal_0" ,"Mesenchymal_1", "Mesenchymal_3","SGN" )
Idents(PCW11_sgn) = factor(PCW11_sgn$cellchat_type, levels = custom_order)
Idents(PCW16_sgn) = factor(PCW16_sgn$cellchat_type, levels = custom_order)

##SV niche
custom_order <- c("Mesenchymal_0", "Mesenchymal_1", "Mesenchymal_3", "Mesenchymal_7", "Melanocyte", "Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_4", "Immunocyte_1", "Immunocyte_2", "Immunocyte_3", "Cycling_RFC",  "MSC1", "MSC2","MEC")
Idents(PCW11_roof) = factor(PCW11_roof$cellchat_type, levels = custom_order)
Idents(PCW16_roof) = factor(PCW16_roof$cellchat_type, levels = custom_order)


##RMC niche
custom_order = c( "Immunocyte_1",  "Immunocyte_2" , "Immunocyte_3" , "Mesenchymal_5" ,"RF"  ,   "Cycling_RFC",  "RMC" ,"SPC")
Idents(PCW11_ves) = factor(PCW11_ves$cellchat_type, levels = custom_order)
Idents(PCW16_ves) = factor(PCW16_ves$cellchat_type, levels = custom_order)


                             
PCW11_down$samples = PCW11_down$batch
PCW11_roof$samples = PCW11_roof$batch
PCW11_sgn$samples = PCW11_sgn$batch
PCW11_ves$samples = PCW11_ves$batch


identity1 <- subset(PCW11_down@meta.data, select = c("cellchat_type","samples") )
identity1 <- subset(PCW11_roof@meta.data, select = c("cellchat_type","samples") )
identity1 <- subset(PCW11_sgn@meta.data, select = c("cellchat_type","samples") )
identity1 <- subset(PCW11_ves@meta.data, select = c("cellchat_type","samples") )


cellchat1 <- createCellChat(PCW11_down, meta = identity1, group.by = "cellchat_type")
cellchat1 <- createCellChat(PCW11_roof, meta = identity1, group.by = "cellchat_type")
cellchat1 <- createCellChat(PCW11_sgn, meta = identity1, group.by = "cellchat_type")
cellchat1 <- createCellChat(PCW11_ves, meta = identity1, group.by = "cellchat_type")


PCW16_down$samples = PCW16_down$batch
PCW16_roof$samples = PCW16_roof$batch
PCW16_sgn$samples = PCW16_sgn$batch
PCW16_ves$samples = PCW16_ves$batch


identity2 <- subset(PCW16_down@meta.data, select = c("cellchat_type","samples") )
identity2 <- subset(PCW16_roof@meta.data, select = c("cellchat_type","samples") )
identity2 <- subset(PCW16_sgn@meta.data, select = c("cellchat_type","samples") )
identity2 <- subset(PCW16_ves@meta.data, select = c("cellchat_type","samples") )


cellchat2 <- createCellChat(PCW16_down, meta = identity2, group.by = "cellchat_type")
cellchat2 <- createCellChat(PCW16_roof, meta = identity2, group.by = "cellchat_type")
cellchat2 <- createCellChat(PCW16_sgn, meta = identity2, group.by = "cellchat_type")
cellchat2 <- createCellChat(PCW16_ves, meta = identity2, group.by = "cellchat_type")


cellchat = cellchat1
cellchat = cellchat2

CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

options(future.globals.maxSize = 1024 * 1024^2)  # 1GB
future::plan("sequential")
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)

cellchat <- computeCommunProb(cellchat, raw.use = FALSE)  
cellchat <- filterCommunication(cellchat, min.cells = 10)  

df.net <- subsetCommunication(cellchat) 

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

object.list <- list(PCW11 = cellchat1, PCW16 = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



cluster_info <- data.frame(
  Cluster_ID = seq_along(levels(cellchat1@idents)),
  Cluster_Name = levels(cellchat1@idents)
)
print(cluster_info)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions - PCW11")
netVisual_circle(cellchat1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength - PCW11")

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat2@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions - PCW16")
netVisual_circle(cellchat2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength - PCW16")


##communication pattern heatmap
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8 ,height = 22,color.heatmap  = "Reds")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 22, color.heatmap  = "Reds")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 22, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 22, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


##signaling role scatterplot
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, dot.size = c(5,15),label.size = 5, dot.alpha = 0.9)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)



#PCW11_16 HC_SGN communication
##SDC "#c4b0d5"
##HC "#289e67"
##SGN "#50bdd1"

##PCW11
cluster_info <- data.frame(
  Cluster_ID = seq_along(levels(cellchat1@idents)),
  Cluster_Name = levels(cellchat1@idents)
)
print(cluster_info)

communication_df <- subsetCommunication(cellchat1)
filtered_com1 <- communication_df[
communication_df$prob > 0.05 &  
communication_df$source == "Neuron" &  
communication_df$target == "SDC",  
]

netVisual_chord_gene(cellchat1, sources.use = 1, targets.use = 2, lab.cex = 0.5,legend.pos.y = 30, color.use = c("#50bdd1", "#c4b0d5"), net = filtered_com1)


filtered_com2 <- communication_df[
  communication_df$prob > 0.05 &  
    communication_df$source == "SDC" &  
    communication_df$target == "Neuron",   
]

netVisual_chord_gene(cellchat1, sources.use = 2, targets.use = 1, lab.cex = 0.5,legend.pos.y = 30, color.use = c("#c4b0d5", "#50bdd1"), net = filtered_com2)

communication_df1 = communication_df

write.csv(filtered_com1, file = "data/cellchat/PCW11_SGN_to_SDC_topPairs.csv")
write.csv(filtered_com2, file = "data/cellchat/PCW11_SDC_to_SGN_topPairs.csv")
write.csv(communication_df1, file = "data/cellchat/PCW11_SDC_SGN.csv")


##PCw16
cluster_info <- data.frame(
  Cluster_ID = seq_along(levels(cellchat2@idents)),
  Cluster_Name = levels(cellchat2@idents)
)
print(cluster_info)

communication_df <- subsetCommunication(cellchat2)
filtered_com1 <- communication_df[
  communication_df$prob > 0.05 &  
    communication_df$source == "HC" &  
    communication_df$target == "SGN",  
]

netVisual_chord_gene(cellchat2, sources.use = 1, targets.use = 2, lab.cex = 0.5,legend.pos.y = 30, color.use = c("#289e67","#50bdd1"), net = filtered_com1)


filtered_com2 <- communication_df[
  communication_df$prob > 0.05 &  
    communication_df$source == "SGN" & 
    communication_df$target == "HC",  
]

netVisual_chord_gene(cellchat2, sources.use = 2, targets.use = 1, lab.cex = 0.5,legend.pos.y = 30, color.use = c( "#50bdd1", "#289e67"), net = filtered_com2)

communication_df2 = communication_df

write.csv(communication_df2, file = "data/cellchat/PCW16_HC_SGN.csv")
write.csv(filtered_com1, file = "data/cellchat/PCW16_HC_to_SGN_topPairs.csv")
write.csv(filtered_com2, file = "data/cellchat/PCW16_SGN_to_HC_topPairs.csv")



#cell-cell communication clustering
library(dplyr)
library(purrr)

##floor: M0- HC_domain; PC/DC_domain; IBC/IPhC_domain; ISC/IdC ; KOC ;OSC
##RMC : M5- RMC; SPC
##SV: M7 - MSC1; MSC2; MEC
##SGN: M4/M8 - SGN

target_groups = c("Mesenchymal_0"  , "HC_domain" ,"PC/DC_domain" , "IBC/IPhC_domain" ,"ISC/IdC"  ,  "KOC" , "OSC",    
                  "Mesenchymal_5","RMC", "SPC", "Mesenchymal_7", "MSC1" ,"MSC2" , "MEC" , "Mesenchymal_4" , "Mesenchymal_8"  ,  "SGN"  )
receptor_group = c( "HC_domain" ,"PC/DC_domain" , "IBC/IPhC_domain" ,"ISC/IdC"  ,  "KOC" , "OSC",    
                   "RMC", "SPC",  "MSC1" ,"MSC2" , "MEC" , "SGN"  )
cellchat_obj = subset(snRNA, subset =  cellchat_type %in% target_groups)
cellchat_receptor = subset(cellchat_obj, subset = cellchat_type %in% receptor_group)

receptor_markers <- FindAllMarkers(
cellchat_receptor,
only.pos = TRUE,            
min.pct = 0.25,             
logfc.threshold = 0.5       
)


cellchat_obj$samples = cellchat_obj$batch
identity <- subset(cellchat_obj@meta.data, select = c("cellchat_type","samples") )
cellchat <- createCellChat(cellchat_obj, meta = identity, group.by = "cellchat_type")
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
options(future.globals.maxSize = 1024 * 1024^2)  # 1GB
future::plan("sequential")
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10) 

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


df.net <- subsetCommunication(cellchat) #汇总所有细胞通信
df.net_subset1 <- subsetCommunication(cellchat, sources.use = "Mesenchymal_0", targets.use = c("HC_domain", "PC/DC_domain", "IBC/IPhC_domain","ISC/IdC", "KOC", "OSC"))
df.net_subset2 <- subsetCommunication(cellchat, sources.use = "Mesenchymal_5", targets.use = c("RMC", "SPC"))
df.net_subset3 <- subsetCommunication(cellchat, sources.use = "Mesenchymal_7", targets.use = c("MSC1", "MSC2", "MEC"))
df.net_subset4 <- subsetCommunication(cellchat, sources.use = c("Mesenchymal_4","Mesenchymal_8" ), targets.use = c("SGN"))

target1 = c("HC_domain", "PC/DC_domain", "IBC/IPhC_domain","ISC/IdC", "KOC", "OSC")
target2 = c("RMC", "SPC")
target3 = c("MSC1", "MSC2", "MEC")
target4 = c("SGN")

receptor1 =  filtered_receptors[filtered_receptors$cluster %in% target1,  ]
receptor2 = filtered_receptors[filtered_receptors$cluster %in% target2,  ]
receptor3 = filtered_receptors[filtered_receptors$cluster %in% target3,  ]
receptor4 = filtered_receptors[filtered_receptors$cluster %in% target4,  ]

df_sorted <- df.net_subset1 %>%
group_by(target) %>%
arrange(desc(prob), .by_group = TRUE) %>%
ungroup()
df.net_subset1 = df_sorted
df.net_subset2 = df_sorted
df.net_subset3 = df_sorted
df.net_subset4 = df_sorted

df_sorted = df.net_subset1

##receptors-specific selection
find_partial_matches <- function(df, marker_genes) {
  matched_rows <- data.frame()
  
  for(i in 1:nrow(df)) {
    receptor_gene <- df$receptor[i]
    
    matches <- sapply(marker_genes, function(marker) {
      grepl(marker, receptor_gene) | grepl(receptor_gene, marker)
    })
    
    if(any(matches)) {
      matched_markers <- marker_genes[matches]
      df$matched_markers[i] <- paste(matched_markers, collapse = ", ")
      matched_rows <- rbind(matched_rows, df[i, ])
    }
  }
  
  return(matched_rows)
}

df_matched_flexible <- find_partial_matches(df_sorted, receptor1$gene)


matched_genes_from_markers <- receptor1 %>%
  filter(sapply(gene, function(marker_gene) {
    any(sapply(df_sorted$receptor, function(receptor_gene) {
      grepl(marker_gene, receptor_gene) | grepl(receptor_gene, marker_gene)
    }))
  }))


complete_double_match_enhanced <- function(df_sorted, matched_genes) {
  df_sorted <- df_sorted %>%
    mutate(across(where(is.factor), as.character))
  
  matched_genes <- matched_genes %>%
    mutate(across(where(is.factor), as.character))
  
  df_sorted %>%
    rowwise() %>%
    mutate(
      is_matched = {
        target_match <- target %in% matched_genes$cluster
        if(target_match) {
          target_genes <- matched_genes$gene[matched_genes$cluster == target]
          if(grepl("_", receptor)) {
            components <- strsplit(receptor, "_")[[1]]
            any(components %in% target_genes)
          } else {
            receptor %in% target_genes
          }
        } else {
          FALSE
        }
      },
      matched_gene = {
        if(target %in% matched_genes$cluster) {
          target_genes <- matched_genes$gene[matched_genes$cluster == target]
          
          if(grepl("_", receptor)) {
            components <- strsplit(receptor, "_")[[1]]
            matching_components <- components[components %in% target_genes]
            if(length(matching_components) > 0) {
              paste(matching_components, collapse = ", ")
            } else {
              NA
            }
          } else {
            if(receptor %in% target_genes) {
              receptor
            } else {
              NA
            }
          }
        } else {
          NA
        }
      }
    ) %>%
    filter(is_matched) %>%
    select(-is_matched) %>%
    ungroup()
}

df_final_matched <- complete_double_match_enhanced(df_sorted, matched_genes_from_markers)
print("增强版匹配结果（包含匹配gene信息）:")
print(df_final_matched)

df_final_1 = df_final_matched
df_final_2 = df_final_matched
df_final_3 = df_final_matched
df_final_4 = df_final_matched
df_final = rbind(df_final_1, df_final_2, df_final_3, df_final_4)

df_final$unit = paste0(df_final$source, "-", df_final$target)
df_final$interaction = paste0(df_final$ligand, "|", df_final$receptor)



##chord plot for selected pathways

# print(cluster_info)
#Cluster_ID    Cluster_Name
#1           1   Mesenchymal_0
#2           2   Mesenchymal_5
#3           3   Mesenchymal_7
#4           4   Mesenchymal_4
#5           5   Mesenchymal_8
#6           6       HC_domain
#7           7 IBC/IPhC_domain
#8           8         ISC/IdC
#9           9             KOC
#10         10             OSC
#11         11    PC/DC_domain
#12         12             RMC
#13         13             SPC
#14         14             MEC
#15         15            MSC1
#16         16            MSC2
#17         17             SGN

color.use = c("#b7d0de", "#d055e5", "#7f7f7f", "#8ef9c1", "#fefd60",
              "#4fa97f", "#e193c9", "#d34f41", "#b25ff1","#f6ca96","#cdbfda",
              "#f3afa7", "#cbaea7",
              "#c3c882", "#5ac5d5", "#c0d1e9",
              "#5464dc")

pathways.show = c("ADGRL")#PTN #EPHA
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[22,]

netVisual_individual(cellchat, signaling = pathways.show,
pairLR.use = LR.show,  
sources.use = c(1,2,3,4,5),  
targets.use = c(6,7,8,9,10,14, 11, 12, 13, 15, 17, 16), 
layout = "chord", color.use = color.use)


netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show,  
                     sources.use = c(1), 
                     targets.use = c(6,7,8,9,10,11), 
                     layout = "chord", color.use = color.use)
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show,  
                     sources.use = c(2),  
                     targets.use = c(12,13),  
                     layout = "chord", color.use = color.use)
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, 
                     sources.use = c(3),  
                     targets.use = c(14,15,16),  
                     layout = "chord", color.use = color.use)
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, 
                     sources.use = c(4,5),  
                     targets.use = c(17), 
                     layout = "chord", color.use = color.use)