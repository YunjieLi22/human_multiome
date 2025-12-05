# A Human Fetal Cochlear Cell Atlas Reveals a Regulatory Blueprint for Spatial Patterning

## Preamable
- This directory contains illustration of the overall workflow of creating human fetal cochlear cell atlas
- Related codes illustrate the workflow of processing single-cell multiome, scRNA-seq and stereo-seq data
- Samples of the sequencing data were all from human cochleae out of PCW11, PCW14 and PCW16.
- All the scripts were runned on Linux and R.4.5.1, unless specially mentioned.
- The original data can be assessed on the ArrayExpress (sn-Multiome：E-MTAB-16331; scRNA-seq: E-MTAB-16336; ).
- Code contributors: Jiacheng Fan, Yunjie Li and Yu Dong.



## Analysis Workflow and Script Descriptions
- This repository contains the key computational pipeline for building a human fetal cochlear cell atlas from multi-modal single-cell sequencing data. Each following script performs a specific analytical step.

| Script Name | Language | Description |
| :--- | :--- | :--- |
| **`01_scRNA_preprocess.R`** | R | Performs initial quality control, filtering, normalization, and feature selection on the raw single-cell RNA-seq data. |
| **`02_Integration.R`** | R | Integration of scRNA-seq and snRNA-seq datasets to correct for technical variance and enable joint analysis. |
| **`03_peak_gene_link_identification.R`** | R | Processes single-cell multiome (ATAC + RNA) data to identify statistically significant links between chromatin accessibility peaks and potential target genes. |
| **`04_peak_annotation_and_motif_enrichment.R`** | R | Annotates identified ATAC-seq peaks to genomic features (promoters, enhancers, etc.) and performs motif enrichment analysis to predict binding transcription factors. |
| **`05_DEG and GO analysis.R`** | R | Identifies differentially expressed genes (DEGs) between cell types or conditions and performs Gene Ontology (GO) enrichment analysis to interpret biological functions. |
| **`06_SCENiCplus.py`** | Python | Runs the SCENIC+ pipeline, which extends SCENIC to multiome data to infer enhancer-driven gene regulatory networks (GRNs) and cis-regulatory interactions. |
| **`07_gene_regulatory_network.py`** | Python | Constructs and analyzes cell-type-specific gene regulatory networks based on transcription factor activity and target gene expression as a downstream analysis of SCENIC+. |
| **`08_perturbation_analysis.py`** | Python | Simulates in silico perturbations (e.g., TF knockout) on the regulatory networks to predict key regulators and network stability. |
| **`09_trajectory_inference.R`** | R | Performs pseudotime analysis and trajectory inference to model cellular differentiation or developmental processes within the cochlea. |
| **`10_cell_cell_communication.R`** | R | Infers intercellular communication networks by analyzing ligand-receptor interactions between different cell types within each functional niche. |
| **`11_PROGENy.R`** | R | Applies the PROGENy method to infer activity of downstream signaling pathways from gene expression data in a cell-type-specific manner. |
| **`12_human_mouse_integration.R`** | R | Integrates and compares the human fetal cochlear atlas with relevant mouse datasets to assess conservation and identify species-specific features. |
| **`13_celltype_disease_enrichment.py`** | Python | Maps cell-type-specific gene signatures to known human disease (e.g., hearing loss) associations using enrichment analysis. |
| **`14_Stereo-seq_disease_enrichment.py`** | Python | Performs spatial enrichment analysis on Stereo-seq data to correlate disease-associated genes or signatures with specific anatomical regions of the cochlea. |
| **`15_SCENIC_TF_regulatory_disease_enrichment.R`** | R | Integrates SCENIC-inferred transcription factor regulons with disease gene sets to identify TFs potentially linked to cochlear pathologies. |

## Usage Notes
- Dependencies: Ensure all required R packages (e.g., Seurat, Signac, Monocle3) and Python libraries (e.g., scanpy, pySCENIC+) are installed before execution.
- Input Data: The pipeline is designed to start from processed count matrices (for scRNA-seq) and fragment files (for ATAC-seq). Raw FASTQ files require separate alignment and counting steps not included here.
- Customization: Parameters for filtering, integration, and analysis thresholds may need adjustment for different datasets. Key parameters are defined within the scripts.


