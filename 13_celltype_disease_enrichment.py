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

warnings.filterwarnings("ignore")

H5AD_FILE = "seurat_macaque_two.h5ad"
DATA_PATH = ''
GS_FILE = os.path.join(DATA_PATH, "all_HL_genes_scDRS.txt")

adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
df_gs = scdrs.util.load_gs(GS_FILE)

adata.obs['orig_new_type1'] = adata.obs['orig_new_type1'].replace({
    'RMC1': 'RMC',
    'RMC2': 'RMC',
    'ISC/IdC': 'ISC'
})

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
df_cor = pd.DataFrame(data= adata.obs, columns=['nCount_RNA','nFeature_RNA'])
scdrs.preprocess(adata, cov=df_cor)
scdrs.preprocess(adata, adj_prop='orig_new_type1')

#calculate score
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

with open("human_two.pkl", "wb") as f:
    pickle.dump(dict_df_score, f)
    
import scanpy as sc
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)  

#calculate difference
df_stats_all = dict()
for key in dict_df_score.keys():
    print(key)
    df_stats_all[key] = scdrs.method.downstream_group_analysis(
        adata=adata,
        df_full_score=dict_df_score[key],
        group_cols=["orig_new_type1"])

with open("human_two_stats.pkl", "wb") as f:
    pickle.dump(df_stats_all, f)
