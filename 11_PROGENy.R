#pathway activity score with PROGENy

library(progeny)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

model <- progeny::model_human_full
model_100 <- model %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 200)

# Plot
ggplot(data=model_100, aes(x=weight, color=pathway, fill=pathway)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  facet_wrap(~ pathway, scales='free') +
  xlab('scores') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")


m_merge <- JoinLayers(m_merge, assay = "RNA", layer = "data")
m_merge <- JoinLayers(m_merge, assay = "RNA", layer = "counts")

m_merge <- progeny(m_merge, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)

DefaultAssay(m_merge) <- "progeny"
m_merge <- Seurat::ScaleData(m_merge)

progeny_scores_df <- 
  as.data.frame(t(GetAssayData(m_merge, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 


levels = levels(as.factor(m_merge$age_clusters))
Idents(m_merge)= factor(m_merge$age_clusters, levels = levels)

CellsClusters <- data.frame(Cell = names(Idents(m_merge)), 
                            CellType = as.character(Idents(m_merge)),
                            stringsAsFactors = FALSE)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

col_fun1 <- colorRamp2(c(-1, -0.5, 0, 0.5, 1),  c( "#2980d3", "#a0d3f4", "#fdfefe", "#fea888", "#d72c23"))
col_fun3 <- colorRamp2(c(-1,   0,  1),  c( "#4886bf",  "#fefefb", "#c07140"))


Heatmap(t(summarized_progeny_scores_df),
        heatmap_legend_param=list(title="pathway activity score"),
        col= col_fun1 ,
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 8), border = "black")

#specific pathways
mat = summarized_progeny_scores_df[ , c(2,6,9,10,11,14)]
order = c("age_1_Mesenchymal_0", "age_2_Mesenchymal_0", "age_3_Mesenchymal_0", "age_1_Mesenchymal_1", "age_2_Mesenchymal_1","age_3_Mesenchymal_1",
          "age_1_Mesenchymal_2", "age_2_Mesenchymal_2","age_3_Mesenchymal_2", "age_1_Mesenchymal_3", "age_2_Mesenchymal_3","age_3_Mesenchymal_3",
          "age_1_Mesenchymal_4", "age_2_Mesenchymal_4","age_3_Mesenchymal_4",  "age_1_Mesenchymal_5", "age_2_Mesenchymal_5","age_3_Mesenchymal_5",
          "age_1_Mesenchymal_6", "age_2_Mesenchymal_6","age_3_Mesenchymal_6", "age_1_Mesenchymal_7", "age_2_Mesenchymal_7","age_3_Mesenchymal_7",
          "age_1_Mesenchymal_8", "age_2_Mesenchymal_8","age_3_Mesenchymal_8" )

mat_ordered = mat[order, ]

df <- mat_ordered %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  select(cluster, value = 2)  

df$group <- rep(1:9, each = 3, length.out = nrow(df))


cols9 <- c("#A4C7D9", "#A30059", "#063FAC",
           "#71C33A", "#63FFAC", "#E14AEC",
           "#8FB0FF", "#7F7F7F", "#FFFF00")


pdf("pictures/Figure3/TNFa_pathway_score_activity.pdf", width = 3, height = 9) #EGFR, MAPK, PI3K, TGFb, TNFa, WNT
ggplot(df,
       aes(x = factor(cluster, levels = rev(unique(df$cluster))),
           y = value,
           color = factor(group))) +
  geom_segment(aes(xend = cluster, yend = 0), size = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = cols9, guide = "none") +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "Expression") +
  theme(
    panel.background = element_rect(fill = "grey90"),
    panel.grid.major.x = element_line(colour = "white", size = 0.6),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_blank()
  )
dev.off()