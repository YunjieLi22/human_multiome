# -*- coding: utf-8 -*-
import scdrs
import scanpy as sc
sc.set_figure_params(dpi=125)
from anndata import AnnData
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings
import sys

warnings.filterwarnings("ignore")

output_file = 'allsubtype_E16.pkl'
output_path2 =  'allsubtype_stat_E16.pkl'
    
H5AD_FILE = 'E16_2_merge_all.h5ad' #"/home/dongyu/project/xiongxian/temp/reuslt/0711singlecell_shuangzuxue/two_human_mouse_all_celltype_0726/integrate_library/spatial_data/1.SAW/different_bin/S16w_cochlea_bin50_scanpy_out.h5ad"
#COV_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.cov")
DATA_PATH = ''
#COV_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.cov")
GS_FILE = os.path.join(DATA_PATH, "all_HL_genes_scDRS.txt")

# Load .h5ad file, .cov file, and .gs file
#adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)

#load data
adata = sc.read_h5ad(H5AD_FILE)
#adata = adata.raw.to_adata()

adata.obs['predicted.celltype'] = adata.obs['subtype_means_cell2location']

#df_cov = pd.read_csv(COV_FILE, sep="\t", index_col=0)
df_gs = scdrs.util.load_gs(GS_FILE)

# normalization
# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)
# Log-transform the data
sc.pp.log1p(adata)

#==adjust cell ratio
scdrs.preprocess(adata, adj_prop='subtype_means_cell2location')

#Computing scDRS scores
dict_df_score = dict()
for trait in df_gs:
    gene_list, gene_weights = df_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=1000,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
import pickle
# save result
with open(output_file, "wb") as f:
    pickle.dump(dict_df_score, f)
    
import scanpy as sc
# Build the adjacency relationship between cells
sc.pp.pca(adata, n_comps=45) 
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

df_stats_all = dict()
#Downstream analyses for SCZ
for key in dict_df_score.keys():
    print(key)
    df_stats_all[key] = scdrs.method.downstream_group_analysis(
        adata=adata,
        df_full_score=dict_df_score[key],
        group_cols=["predicted.celltype"])#["orig_new_type1"]
    #scdrs.util.plot_group_stats(df_stats)
#display(df_stats["orig_new_type1"].style.set_caption("Group-level statistics for HP_0000407__EuropeanAndAfrican"))   

# save result
with open(output_path2, "wb") as f:
    pickle.dump(df_stats_all, f)
