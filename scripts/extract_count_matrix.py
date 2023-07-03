#!/usr/bin/env python

# > June 2023                                                                                                              
# > Script : Extract_count_matrixes.R                                                                                                      
# > Function : extract csv of gene expression from Tabula sapiens data               
# @ COLAJANNI Antonin                                                          
################################################################################

# # Loading and filtering features from Single cell dataset 
import pandas as pd        
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display
from collections import Counter

path_data = "/shared/projects/microbiome_translocation/data/Tabula_sapiens_immune_all/"

adata = sc.read(path_data + "immune_all_cell.h5ad")

blood_adata = adata[adata.obs.tissue_in_publication.isin(["Blood"]),:]

df = blood_adata.to_df()

# Renaming genes and cell id
gene_name = np.asarray(blood_adata.var.feature_name).tolist()
gene_id = np.asarray(blood_adata.var.index).tolist()
ID = blood_adata.obs.index.tolist()
types = np.asarray(blood_adata.obs.cell_type).tolist()

df = df.rename(columns=dict(zip(gene_id,gene_name)),index=dict(zip(ID,types)))
df_no_zero = df.loc[:, (df.sum(axis=0) > 0)]
df_no_zero.to_csv(path_data+"50000cell_immune_expressed.csv", sep='\t')


file = open(path_data+'50k_cell_labels.txt','w')
for cells in df.index.to_list():
	file.write(cells+"\n")
file.close()

file = open(path_data+'50k_cells_geneset/50k_cell_genes.txt','w')
for cells in df_no_zero.columns.to_list():
	file.write(cells+"\n")
file.close()
