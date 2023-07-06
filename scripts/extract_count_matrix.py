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
from functions.utils import *
from functions.function_h5ad import *
from functions.function_SingleCell import *



path_data = "/shared/projects/microbiome_translocation/data/Tabula_sapiens_immune_all/"

#################################################################
# Tabula Sapiens dataset

adata = sc.read(path_data + "immune_all_cell.h5ad")
blood_adata = adata[adata.obs.tissue_in_publication.isin(["Blood"]),:]

# Renaming genes and cell id and transform into dataframe
df = adata_to_expr_df(blood_adata)

df_no_zero = df.loc[:, (df.sum(axis=0) > 0)]

# Saving dataframe, label, genes
df_no_zero.to_csv(path_data+"Tabula_Sapiens.csv", sep='\t')

file = open(path_data+'50k_cells_geneset/50k_cell_labels.txt','w')
for cells in df_no_zero.index.to_list():
	file.write(cells+"\n")
file.close()

file = open(path_data+'50k_cells_geneset/50k_cell_genes.txt','w')
for cells in df_no_zero.columns.to_list():
	file.write(cells+"\n")
file.close()


################################################################# 
# Liu et al 2021 Cell

path_data = "/shared/projects/microbiome_translocation/data/scRNAseq/Liu_2021_cell/"

adaptative = sc.read(path_data+"adaptative.h5ad")
adaptative = adaptative[adaptative.obs.disease.isin(["normal"]),:]
adaptative_df = adata_to_expr_df(adaptative)

innate = sc.read(path_data+"innate.h5ad")
innate = innate[innate.obs.disease.isin(["normal"]),:]
innate_df = adata_to_expr_df(innate)

# Merging cells from innate and adaptative immunity
merged_data = pd.concat([innate_df,adaptative_df])

indexes = merged_data.index.to_list()

label_list = ["CD4","CD8","Monocyte","memory_B_cell","naive_B_cell","NK_cell"]

# Removing cells that does not interest us
new_index = replace_label(indexes, 
              pattern_list=["CD4","CD8","monocyte","memory B cell","naive B cell","natural"],
              replace_list=label_list,
              replace_other='other' )

# Remove cell with the label 'other'
bool_list = []
for item in new_index:
    if item in label_list:
        bool_list.append(True)
    else:
        bool_list.append(False)
        
merged_data.index = new_index
merged_data = merged_data[bool_list]
merged_data_no_zero = merged_data.loc[:, (merged_data.sum(axis=0) > 0)]

# Saving dataframe, label, genes
merged_data_no_zero.to_csv(path_data+"Liu2021_Cell.csv", sep='\t')

file = open(path_data+'Liu2021_index.txt','w')
for cells in merged_data_no_zero.index.to_list():
	file.write(cells+"\n")
file.close()

file = open(path_data+'Liu2021_genes.txt','w')
for genes in merged_data_no_zero.columns.to_list():
	file.write(genes+"\n")
file.close()