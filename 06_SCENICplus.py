#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
os.chdir("/home/jcfan/human_brain/subtype_scenicplus/NEWREH/outs")


# In[3]:


import mudata
scplus_mdata = mudata.read("scplusmdata.h5mu")


# In[4]:


import scipy.sparse
import pandas as pd
if scipy.sparse.issparse(scplus_mdata["direct_gene_based_AUC"].X):
    matrix = scplus_mdata["direct_gene_based_AUC"].X.toarray()
else:
    matrix = scplus_mdata["direct_gene_based_AUC"].X

df_auc = pd.DataFrame(matrix, index=scplus_mdata["direct_gene_based_AUC"].obs_names, columns=scplus_mdata["direct_gene_based_AUC"].var_names)


# In[7]:


df_auc.to_csv("gene_based_AUC.csv")


# In[3]:


scplus_mdata.uns["direct_e_regulon_metadata"]


# In[4]:


scplus_mdata.uns["extended_e_regulon_metadata"]


# In[3]:


import scanpy as sc
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)


# In[6]:


eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]


# In[7]:


sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")


# In[8]:


sc.tl.umap(eRegulon_gene_AUC)


# In[9]:


sc.pl.umap(eRegulon_gene_AUC, color = "scATAC_counts:subtype")


# In[3]:


scplus_mdata


# In[9]:


from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
import numpy as np
import pandas as pd


# In[7]:


rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scATAC_counts:subtype",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)


# In[11]:


eregulons_to_plot = sorted(set(
    np.concatenate([rss.loc[g].nlargest(40).index.values for g in rss.index])
))




# In[19]:


meta["eRegulon_name"]


# In[13]:


plot_rss(
    data_matrix = rss,
    top_n = 3,
    num_columns = 5
)


# In[14]:


sc.pl.umap(eRegulon_gene_AUC, color = list(set([x for xs in [rss.loc[ct].sort_values()[0:2].index for ct in rss.index] for x in xs ])))


# In[5]:


from scenicplus.plotting.dotplot import heatmap_dotplot


# In[6]:


meta_key = "direct_e_regulon_metadata"
meta = scplus_mdata.uns[meta_key].copy()  


# In[11]:


scplus_mdata.uns


# In[7]:


TF_corr_col   = "TF_GEX_correlation"
Peak_corr_col = "Peak_ACC_correlation"
reg_type_col  = "regulation_type"           # activator / repressor
name_col      = "eRegulon_name"


# In[9]:


meta


# In[8]:


def sign_pair(row):
    tf_sign   = '+' if row[TF_corr_col]   > 0 else '-'
    peak_sign = '+' if row[Peak_corr_col] > 0 else '-'
    return f"{tf_sign}/{peak_sign}"

meta["sign_pair"] = meta.apply(sign_pair, axis=1)



heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scATAC_counts:subtype",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "TF",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    split_repressor_activator=False
    figsize = (25, 16),save="allTFdotplot.pdf"
)

