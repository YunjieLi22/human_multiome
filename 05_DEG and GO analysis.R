
#DEG analysis and GO analysis
library(tidyverse)
library(ggh4x)
library(ggfun)
library(ggnewscale)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)

all_markers = FindAllMarkers(newreh_merge, only.pos = TRUE,           
min.pct = 0.25,             
logfc.threshold = 0.5       
)

top_markers <- all_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.58) %>%  
  group_by(cluster) %>%                          
  arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE)

desired_cluster_order = levels(newreh_merge)
top_markers$cluster <- factor(top_markers$cluster, levels = desired_cluster_order)
ordered_markers <- top_markers %>%
  arrange(cluster)
gene_list <- split(ordered_markers$gene, ordered_markers$cluster)
sapply(gene_list, length)

go_results <- lapply(names(gene_list), function(cluster) {
  genes <- gene_list[[cluster]]
  eg <- bitr(genes, fromType = "SYMBOL", 
             toType = "ENTREZID", 
             OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene          = eg$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",        
                  pAdjustMethod = "BH",         
                  pvalueCutoff = 0.1,         
                  qvalueCutoff = 0.5,         
                  readable      = TRUE)   
  if (!is.null(ego)) {
    ego@result$Cluster <- cluster  
  }
  return(ego)
})


  
names(go_results) <- names(gene_list)

top_go <- bind_rows(lapply(go_results, function(x) {
  if (!is.null(x)) {
    x@result %>% 
      filter(p.adjust < 0.05) %>%  
      arrange(p.adjust) %>%        
      slice_head(n = 10)         
  }
}))

write.csv(top_go, file = "data/Figure2/top10_newreh_BP_GO.csv")

celltype1 = levels(newreh_merge)
celltype

#reorder gene_list
gene_list_ordered <- gene_list[levels(newreh_merge)]
sapply(gene_list_ordered, length)

newreh_GO <- read_excel("data/Figure2/newreh_GO.xlsx")

celltype = c( "RMC",  "RMC_1"   ,"RMC_2" ,   "SPC"   ,  "SPC_1" ,"SPC_2" ,    "RF" , "RF_1"  , "RF_2",   "Cycling_RFC","Cycling_RFC_1", "Cycling_RFC_2","MSC1"    ,  "MSC1_1","MSC1_2" ,   "MSC2"    , "MSC2_1" ,"MSC2_2"  ,  "FF"    ,   "FF_1" ,  "FF_2"  ,     "CyclingFF"  ,"CyclingFF_1" ,"CyclingFF_2" , "InC"   ,  "InC_1"  ,"InC_2" ,     "KOC"   ,  "KOC_1" , "KOC_2",  "LSDC" ,   "LSDC_1", "LSDC_2",       "SDC"    ,     "SDC_1" ,  "SDC_2",   "HC"    ,   "HC_1", "HC_2"  ,    "OSC"    , "OSC_1","OSC_2"  ,    "MEC", "MEC_1" ,"MEC_2"      )
celltype1 = levels( newreh_merge)

cell_data <- data.frame( CellType = celltype1,
                         MarkerNumber = c( 316 , 280 ,178 ,157 , 240,150 ,230 ,134 ,275,182,157 , 217, 391 , 182, 129
                         ))
cell_data$CellType <- factor(cell_data$CellType,levels = rev(celltype1))

go_data <- data.frame( GOItem = newreh_GO$Description ,Pvalue = newreh_GO$p.adjust)

cols1 = c("#ff9896", "#c49c94", "#99df8a", "#1f76b3","#15becf",   "#afc7e8","#ff7f0f", "#28a3ff","#d72929", "#ab40fc", "#e377c3","#c4b0d5",  "#289e67", "#ffbb78", "#b6bd61" )
color = rep(rev(cols1), each = 3)


go_data <- go_data %>%
  mutate(Pvalue = as.numeric(as.character(Pvalue))) %>%
  mutate(logP = -log10(Pvalue))
go_data$GOItem <- factor(go_data$GOItem, levels = rev(go_data$GOItem))
go_data$original_order <- 1:nrow(go_data)

p1 <-ggplot(cell_data,aes(x = MarkerNumber, y = CellType)) +
  geom_segment(aes(y = CellType, yend = CellType, x =0, xend = MarkerNumber), linewidth =0.8, color ='grey80') +
  geom_point(aes(color = CellType), size =6) +scale_color_manual(values = color) +scale_x_continuous(limits =c(0,500),expand =c(0,0)) +
  labs(title ="Cell marker numbers", x =NULL, y =NULL) +theme_bw(base_size =18) +theme(plot.title =element_text(hjust =0.5, face ="bold", size =18),
                                                                                       axis.text.x =element_text(color ="black"), axis.text.y =element_text(hjust =1,color ="black"),
                                                                                       panel.grid =element_blank(),legend.position ="none" )

p2 <- ggplot(go_data, aes(y = GOItem)) +
  geom_bar(aes(x = logP, fill = GOItem),
           stat = "identity", width = 0.5, 
           color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.3, label = GOItem), 
            hjust = 0, size = 5.5) +
  labs(title = "GO enrichment item", 
       x = "-log10(Pvalue)", y = NULL) +
  scale_fill_manual(values = color) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(), 
    panel.grid = element_blank() 
  )

p1 + p2 +plot_layout(widths =c(1,1.5))



