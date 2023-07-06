################################################################################                                                                           
# > December 2022                                                                                                             
# > Script : function_h5ad.py                                                                                                           
# > Function : function to extract informations and count matrixes for Single-seq data (.h5ad)                               
# @ COLAJANNI Antonin                                                          
################################################################################



import scanpy as sc
import numpy as np
import pandas as pd


"""
For a given cell type and a tissue of an AnnData object, count the number of cell that have a non-null expression for each gene

Parameters
----------
    adata : AnnData object
        scRNAseq data type by scanpy

    tissue : character string

    cell_type : character string

    arg_filter1 : character string
        one element of adata.obs

    arg_filter2 : character string
        one element of adata.obs

    filter_cell : character string
        string to filter a specific type of cell

Returns
----------
    list of int
"""
def Count_expressing_cell(adata, tissue, cell_type, arg_filter1 = "tissue_in_publication", arg_filter2 = "cell_type"):
    # axis=0 : columnwise (in this format, we have cell X genes)
    N = np.asarray( (adata[adata.obs[arg_filter2].isin([cell_type]) & 
                           adata.obs[arg_filter1].isin([tissue]),:].X != 0).sum(axis = 0)).ravel()
    return N 


"""
For a given cell type and a tissue of an AnnData object, computes the expression sum for each gene

Parameters
----------
    adata : AnnData object
        scRNAseq data type by scanpy

    tissue : character string

    cell_type : character string

    arg_filter1 : character string
        one element of adata.obs

    arg_filter2 : character string
        one element of adata.obs

    filter_cell : character string
        string to filter a specific type of cell

Returns
----------
    list of float
"""
def Expression_sum(adata, tissue, cell_type, arg_filter1 = "tissue_in_publication", arg_filter2 = "cell_type"):

    Total = np.asarray(adata[adata.obs[arg_filter2].isin([cell_type]) & 
                           adata.obs[arg_filter1].isin([tissue]),:].X.sum(axis=0)).ravel()
    return Total

"""
Returns the element of a list matching at least partially a string

Parameters
----------
    string : character string

    c_list : list

Returns
----------
    list of character strings
"""
def str_detect(string, c_list) : 
    return [item for item in c_list if string in item]


"""
Convert np.inf into np.nan in vector

Parameters
----------
    vector : list

Returns
----------
    list
"""
def replace_NaN(vector):
    vector = vector.astype('float')
    vector[vector == np.inf] = np.nan
    return vector

"""
Convert AnnData object to a python Dictonnary of tissue specific, cellular specific mean expression 

Parameters
----------
    adata : AnnData object
        scRNAseq data type by scanpy
    arg_filter1 : character string
        one element of adata.obs

    arg_filter2 : character string
        one element of adata.obs

    filter_cell : character string
        string to filter a specific type of cell

    filter_tissue : character string
        string to filter a specific tissue

    zero_filtering : boolean
Returns
----------
    dict of dataframes
"""
def Anndata_to_mean_expr_df(adata, arg_filter1 = "tissue_in_publication", 
                        arg_filter2 = "cell_type", 
                        filter_cell=None, 
                        filter_tissue=None,
                        zero_filtering = False) :
    
    tissu_dict = {}
    ## Filtering on tissue
    if filter_tissue != None : 
        filter_tissue = str_detect(filter_tissue, list(adata.obs[arg_filter1].cat.categories) )
    else : 
        filter_tissue = adata.obs[arg_filter1].cat.categories

    ## Filtering on cell
    if filter_cell != None : 
        filter_cell = str_detect(filter_cell, list(adata.obs[arg_filter2].cat.categories) )
    else : 
        filter_cell = adata.obs[arg_filter2].cat.categories

    ## 1st for loop
    for tissue in filter_tissue:
        
        tissu_dict[tissue] = pd.DataFrame(columns=adata.var.feature_name,
                                          index=filter_cell)

        
        for cell in filter_cell :
            if zero_filtering :
                #Total of expression value
                Sum = Expression_sum(adata, tissue, cell)
                #Number of cell that have a non null expression
                N = Count_expressing_cell(adata, tissue, cell)

                tissu_dict[tissue].loc[cell] = replace_NaN(Sum / N)
                # Removing artifact of number divided by 0
                tissu_dict[tissue].replace(np.inf,np.nan,inplace=True)
        
            else : 
                expr_mat = adata[adata.obs[arg_filter2].isin([cell]) & 
                                           adata.obs[arg_filter1].isin([tissue]),:].X
        
                if expr_mat.shape[0] != 0 :
                    tissu_dict[tissue].loc[cell] = expr_mat.mean(0)
        
        tissu_dict[tissue] = tissu_dict[tissue].T

    return(tissu_dict)

"""
Convert AnnData object to a python Dictonnary of tissue specific, cellular specific, mean expression where 0 values where not accounted

Parameters
----------
    adata : AnnData object
        scRNAseq data type by scanpy
    arg_filter1 : character string
        one element of adata.obs

    arg_filter2 : character string
        one element of adata.obs

    filter_cell : character string
        string to filter a specific type of cell

Returns
----------
    dict of dataframes
"""
def get_expr_per_CellOrgan(adata, arg_filter1 = "tissue_in_publication", arg_filter2 = "cell_type", filter_cell=None) :
   ##### OLD 
    tissu_dict = {}
    
    for tissue in adata.obs[arg_filter1].cat.categories:
        # Index : genes
        tissu_dict[tissue] = pd.DataFrame(index=adata.var.feature_name)
        
        for cell in adata.obs[arg_filter2].cat.categories:  
            #Total of expression value
            Sum = Expression_sum(adata, tissue, cell)
            #Number of cell that have a non null expression
            N = Count_expressing_cell(adata, tissue, cell)

            tissu_dict[tissue][cell] = Sum / N
        # Removing artifact of number divided by 0
        tissu_dict[tissue].replace(np.inf,np.nan,inplace=True)
        
    return(tissu_dict)

"""
Convert AnnData object to a python dataframe with expression from all different cell types

Parameters
----------
    adata : AnnData object
        scRNAseq data type by scanpy
        
    rename_cell : boolean
        
    rename_gene : boolean string

Returns
----------
    pandas.dataframe
"""
def adata_to_expr_df(adata, rename_cell=True, rename_gene=True):

    df = adata.to_df()

    # Renaming genes and cell id
    gene_name = np.asarray(adata.var.feature_name).tolist()
    gene_id = np.asarray(adata.var.index).tolist()
    ID = adata.obs.index.tolist()
    types = np.asarray(adata.obs.cell_type).tolist()

    if rename_cell and rename_gene :
        df = df.rename(columns=dict(zip(gene_id,gene_name)),index=dict(zip(ID,types)))
    elif rename_cell and not rename_gene : 
        df = df.rename(index=dict(zip(ID,types)))
    elif not rename_cell and rename_gene : 
        df = df.rename(columns=dict(zip(gene_id,gene_name)))

    return(df)
