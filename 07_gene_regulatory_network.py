"""
Export per‑cell‑type TF‑Gene mapping tables (simple version) from a SCENIC+ object.

Workflow:
1. For each cell subtype, restrict to its eRegulons.
2. Use `create_nx_tables()` to obtain edge tables.
3. Extract the TF→Gene edge list (`TF2G`).
4. Keep only the `TF` and `Gene` columns, drop duplicates.
5. Write to `<celltype>_TF_Gene_table_simple.csv`.
"""

import pickle
from pathlib import Path

import pandas as pd
from scenicplus.networks import create_nx_tables
import scipy.sparse
import pandas as pd
import scanpy as sc
import anndata
import dill
import pickle
from scenicplus.networks import *
import networkx as nx
# In[74]:
import matplotlib.pyplot as plt
import matplotlib as mp
from matplotlib import colors, cm

# ---------------------------------------------------------------------
# 0. I/O paths – adjust if your directory structure changes
# ---------------------------------------------------------------------
BASE_DIR = Path("/home/jcfan/human_brain/subtype_scenicplus/NEWREH/outs")
SCPLUS_PKL = BASE_DIR / "scplus_obj.pkl"  # pickled SCENIC+ object

# ---------------------------------------------------------------------
# 1. Load SCENIC+ object
# ---------------------------------------------------------------------
print("[INFO] Loading SCENIC+ object…")
with SCPLUS_PKL.open("rb") as fh:
    scplus_obj = pickle.load(fh)
print("[INFO]   ✔ loaded")
import mudata
scplus_mdata = mudata.read("/home/jcfan/human_brain/subtype_scenicplus/NEWREH/outs/scplusmdata.h5mu")
scplus_obj.metadata_cell["subtype"]=scplus_mdata.obs["scATAC_counts:subtype"]

# ---------------------------------------------------------------------
# 2. Define the eRegulon lists for each cell subtype
# ---------------------------------------------------------------------
REGULONS_BY_SUBTYPE = {
    "HC": [
        "ATOH1_direct", "BACH2_direct", "GFI1_direct", "GTF2IRD1_direct",
        "IRX2_direct", "POU2F1_direct", "POU4F3_direct", "RFX7_direct",
        "VEZF1_direct", "YY1_direct",
    ],
    "SDC": [
        "GATA3_direct", "LIN28B_direct", "MXI1_direct", "SOX11_direct",
        "SOX2_direct", "TCF4_direct", "TPRS1_direct", "ZNF69_direct",
    ],
    "LSDC": [
        "EBF1_direct", "EBF3_direct", "ISL1_direct", "SOX5_direct",
        "SOX6_direct", "TFCP2_direct", "THRB_direct",
    ],
    "ISC": [
        "CHD1_direct", "EGR1_direct", "ELF1_direct", "FOS_direct",
        "NFIB_direct", "NFIX_direct", "RB1_direct", "ZNF385D_direct",
    ],
    "KOC": [
        "NCOA1_direct", "SIX1_direct", "SOX6_direct", "SREBF2_direct",
        "THRB_direct",
    ],
    "OSC": [
        "ETV1_direct", "ETV5_direct", "GATA2_direct", "GATA3_direct",
        "LIN28B_direct", "PBX1_direct", "SMAD2_direct", "SOX4_direct",
        "SREBF2_direct", "TFDP2_direct",
    ],
    "MSC1": [
        "JUN_direct", "MTF2_direct", "PAX2_direct", "TEAD1_direct",
        "ZFP64_direct", "ZNF407_direct", "ZNF518A_direct",
    ],
    "MSC2": [
        "KLF3_direct", "MITF_direct", "NCOA1_direct", "PDLIM5_direct",
        "ZNF407_direct",
    ],
    "RMC": [
        "IRX2_direct", "LMX1A_direct", "MEIS1_direct", "MEIS2_direct",
    ],
    "SPC": [
        "NFIA_direct", "OTX2_direct", "PDLIM5_direct", "GATA2_direct",
    ],
}

# ---------------------------------------------------------------------
# 3. Iterate over cell subtypes and write TF–Gene tables
# ---------------------------------------------------------------------
OUTPUT_DIR = BASE_DIR / "tf_gene_tables"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for subtype, regulons in REGULONS_BY_SUBTYPE.items():
    print(f"[INFO] Processing {subtype} (n_regulons = {len(regulons)})…")

    # 3.1 Build nx_tables restricted to this regulon set
    nx_tables = create_nx_tables(
        scplus_obj,
        eRegulon_metadata_key="eRegulon_metadata",
        subset_eRegulons=regulons,
        add_differential_gene_expression=False,
        add_differential_region_accessibility=False,
    )

    # 3.2 Extract simple TF‑Gene list, remove duplicates
    tf_gene = (
        nx_tables["Edge"]["TF2G"].copy()
        .drop_duplicates()
        .sort_values(["TF", "Gene"])
    )

    # 3.3 Save
    out_path = OUTPUT_DIR / f"{subtype}_TF_Gene_table_simple.csv"
    tf_gene.to_csv(out_path, index=False)
    print(
        f"[INFO]   → saved {len(tf_gene):,} rows to {out_path.relative_to(BASE_DIR)}"
    )

print("[DONE] All subtype tables exported.")
