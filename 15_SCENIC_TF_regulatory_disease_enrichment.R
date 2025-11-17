sample.combined=readRDS('seurat_Macaca_SCT_WuZhongJian.rds.gz')

#===mouse network enrichment 
regulonTargetsInfo_mouse=readRDS('2.5_regulonTargetsInfo_gene.Rds')

all_HL_genes=readRDS('all_HL_genes.rds')
names(all_HL_genes)[which(names(all_HL_genes)=='Auditory.Neuropathy')]='Auditory.Neuropathy.Nonsyndromic'#'Alport.Syndrome'
names(all_HL_genes)[which(names(all_HL_genes)=='Modifiers')]='Modifiers.Nonsyndromic'#'Alport.Syndrome'
names(all_HL_genes)

library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(qvalue)
 library(Seurat) 

background_genes=unique(c(unique(regulonTargetsInfo_mouse$gene_2),rownames(sample.combined)))
tf_targets=split(regulonTargetsInfo_mouse$gene_2,regulonTargetsInfo_mouse$TF_2)
func_sets=all_HL_genes[-c(8,9,10)]

enrich_result <- data.frame()

for (tf in names(tf_targets)) {
  tf_genes <- tf_targets[[tf]]
  n_tf_genes <- length(tf_genes)  
  
  for (func in names(func_sets)) {
    func_genes <- func_sets[[func]]
    n_func_genes <- length(func_genes)  
    n_background <- length(background_genes) 
    
    overlap_genes <- intersect(tf_genes, func_genes)
    n_overlap <- length(overlap_genes) 
    
    p_value <- 1 - phyper(
      q = n_overlap - 1,               
      m = n_func_genes,                
      n = n_background - n_func_genes,   
      k = n_tf_genes                   
    )
        result_row <- data.frame(
      TF = tf,
      Functional_Set = func,
      Background_Size = n_background,
      TF_Target_Size = n_tf_genes,
      Func_Set_Size = n_func_genes,
      Overlap_Size = n_overlap,
      Overlap_Genes = paste(overlap_genes, collapse = ", "),  
      P_Value = p_value
    )
    enrich_result <- rbind(enrich_result, result_row)
  }
}
#enrichment calculation
enrich_result <- enrich_result %>%
  group_by(TF) %>%
  mutate(FDR = p.adjust(P_Value, method = "fdr")) %>% 
  ungroup() %>%
  mutate(Significance = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01 ~ "**",
    FDR < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

#visualization
library(tidyr)
plot_data <- enrich_result %>%
  filter(P_Value < 0.05) %>%   
  mutate(Functional_Set_Simple = gsub("GO_", "", Functional_Set)) %>%
   mutate(Overlap_Ratio = Overlap_Size / TF_Target_Size)

 ggplot(plot_data, aes(x = TF, y = Functional_Set_Simple)) +
   geom_point(aes(size = Overlap_Size, color = -log10(P_Value)), alpha = 0.8) +
   scale_color_gradient(low = "orange", high = "red", name = "-log10(P_Value)") +
   scale_size(range = c(2, 8), name = "Overlap Genes") +
   theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  ) +
   labs(
    x = "Functional Set",
    y = "TF",
    title = "Significant Overlap Between TF Targets and disease (FDR<0.05)"
  )

#===human network enrichment
regulonTargetsInfo=read.delim2('eRegulon_filtered.tsv')

library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(qvalue) 
 library(Seurat) 

background_genes=unique(c(unique(regulonTargetsInfo$Gene),rownames(sample.combined)))
tf_targets=split(regulonTargetsInfo$Gene,regulonTargetsInfo$TF)
func_sets=all_HL_genes[-c(8,9,10)]

enrich_result <- data.frame()

for (tf in names(tf_targets)) {
  tf_genes <- tf_targets[[tf]]
  n_tf_genes <- length(tf_genes) 
  
  for (func in names(func_sets)) {
    func_genes <- func_sets[[func]]
    n_func_genes <- length(func_genes)  
    n_background <- length(background_genes) 
    
    overlap_genes <- intersect(tf_genes, func_genes)
    n_overlap <- length(overlap_genes) 
    
    p_value <- 1 - phyper(
      q = n_overlap - 1,               
      m = n_func_genes,                
      n = n_background - n_func_genes,  
      k = n_tf_genes                   
    )
    
    result_row <- data.frame(
      TF = tf,
      Functional_Set = func,
      Background_Size = n_background,
      TF_Target_Size = n_tf_genes,
      Func_Set_Size = n_func_genes,
      Overlap_Size = n_overlap,
      Overlap_Genes = paste(overlap_genes, collapse = ", "), 
      P_Value = p_value
    )
    enrich_result <- rbind(enrich_result, result_row)
  }
}

enrich_result <- enrich_result %>%
  group_by(TF) %>%
  mutate(FDR = p.adjust(P_Value, method = "fdr")) %>% 
  ungroup() %>%
  mutate(Significance = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01 ~ "**",
    FDR < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

#visualization
plot_data <- enrich_result %>%
  filter(P_Value < 0.05) %>%  
  mutate(Functional_Set_Simple = gsub("GO_", "", Functional_Set)) %>%
  mutate(Overlap_Ratio = Overlap_Size / TF_Target_Size)

ggplot(plot_data, aes(x = TF, y = Functional_Set_Simple)) +
  geom_point(aes(size = Overlap_Size, color = -log10(P_Value)), alpha = 0.8) +
  scale_color_gradient(low = "orange", high = "red", name = "-log10(P_Value)") +
  scale_size(range = c(2, 8), name = "Overlap Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "Functional Set",
    y = "TF",
    title = "Significant Overlap Between TF Targets and disease (FDR<0.05)"
  )


