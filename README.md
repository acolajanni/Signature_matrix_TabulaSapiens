# Signature matrix built from Tabula Sapiens single cell dataset

# Work in progress

**Be Carefull, the current signature matrix is built on 16487 cells, not 50115 as it should be, results will differ**

**Contact:**
Antonin Colajanni: antonin.colajanni@u-bordeaux.fr

# Aim

This page describes the composition and the method of the reference matrix built from Tabula Sapiens single cell dataset. It contains: 
- T CD4 +
- T CD8 +
- Memory B cell
- Naive B cell
- Monocyte
- Natural killer cell
- Neutrophil
- Macrophage
- Platelet
- Plasma cell (plasmocyte)



## Heatmap of the signature matrix

*work in progress*


## File Description: 

The [.html document](/markdown/feature_selection_v2.html) describes the filtration steps before and after the permutaions steps. The code to produce the figures is the [*feature_selection.R*](/scripts/Feature_selection.R) script (à ajouter quand il sera fini), in the same [folder](/scripts/) as the other scripts.


The [results](/results) folder contains the .csv with the permutation importance of each variable for the 50 permutations. The obtained features after selection for each comparison, from the broader selection to the final one are available [here](/genesets).

At last, the evaluated classification performance are contained in the folder [performance](/performance).

## Data
THE TABULA SAPIENS CONSORTIUM. 2022. « The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans ». Science 376 (6594): eabl4896. https://doi.org/10.1126/science.abl4896.

**URL:** 
https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5 - *Tabula Sapiens - Immune* (filtered on blood tissues) 

Alternatively, the data *Tabula Sapiens - Blood* can be directly downloaded.


## Cross Validation Data
Liu, Can, Andrew J. Martins, William W. Lau, Nicholas Rachmaninoff, Jinguo Chen, Luisa Imberti, Darius Mostaghimi, et al. 2021. « Time-Resolved Systems Immunology Reveals a Late Juncture Linked to Fatal COVID-19 ». Cell 184 (7): 1836-1857.e22. https://doi.org/10.1016/j.cell.2021.02.018.

**URL:** https://cellxgene.cziscience.com/collections/ed9185e3-5b82-40c7-9824-b2141590c7f0




## Method description 

<p align="center">
<img src="/doc/Diapo_pathseq-TabulaSapiens.drawio.png" height="500">


After the obtention of the normalised count table, a first step of feature engineering is done. First, for each celltype the mean expression of it and without the said celltype is computed. This way, a fold change between the background and the given cell can be calculated.

Then, the number of cell expressing each gene for each celltype is also computed in order to filter on the percentage of cell expressing a gene. 

With these two parameters, a filtration can be done. First if we compare one celltype vs all other, the filtration is done on the absolute value of the log2 Fold Change being superior to 1, meaning, a gene must be, at least twice as expressed in the given celltype compared to be background and conversly and it must be expressed in 25% of the cells of this given celltype. Alternatively, the curve of the absolute Fold change is generated, and all the genes with a fold change higher than the inflexion point of the curve* are selected, independantly of the proportion of cell expressing it. 

For the alternative, a gene can be a signature for a given celltype if it is expressed everywhere except in the said celltype. And, for the 25% threshold, it is supposed to be relatively generous to avoid selecting false negative gene. Moreover, for some celltypes, I aggregated several sub celltypes, it is the case for the T.CD4 cells for example, that corresponds to the aggregation of 
*CD4-positive, alpha-beta T cell*, *CD4-positive, alpha-beta memory T cell* and *naive thymus-derived CD4-positive, alpha-beta T cell* that are not equally represented.

*using the kneedle algorithm contained is the eponym package in R 


For the Celltype vs Celltype comparison, it is done between certain celltypes, known to be relatively similar, to help the model finding disciminative genes that could have been obscured by the comparison with the Fold change computed with the background mean expression. The same threshold the Fold change is used, however, it has been found that a threshold at 25% of cell expressing a gene was too stringent, so the choice has been made to select genes expressed in at least 2 cells in the given celltypes. 


The chosen algorithm to realise permutation on, and compute variable importance is the [Random Forest implemented in the scikit-learn package](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html#sklearn.ensemble.RandomForestClassifier). The forests are set to contain 1000 trees, with the otehr parameters set to default.
This classifier is then embeded in the [permutation_importance](https://scikit-learn.org/stable/modules/generated/sklearn.inspection.permutation_importance.html#sklearn.inspection.permutation_importance) method which will permute the expression value of a given gene between groups to compare the effect on the decrease in mean squared error (MSE) relative to a situation without permutation. This way, each gene will have it's values permutated 50 times. The mean decrease in MSE for the 50 permutations is the permutation Variable Importance (VI) used to select the genesets for each celltype. 


The other possibility for feature selection was to use a logarithmic regression with a lasso penalty, where the low coefficients will shrink towards zero for the less informative variables. However, given the number of features, the model had to be used several times to converge towards a relatively low number of selected features which is not ideal. Moreover, in this kind of model, each variable is observed independantly as it will be far too complex to model the interactions between each gene for the selection. On the contrary with the Random Forest, for each nodes of a decision tree an interaction is modelized. For example if gene A > 2.5 & gene B > 4 & Gene C < ...   

The edges of the decision trees models the interaction between the expression of several genes.


To evaluate the quality of the prediction with the selected features, the single cell transcriptomics data from Liu et al, 2021, Cell.

Since we already used the Tabula Sapiens data more than once to select the features and compute the variable importances, this dataset has been chosen, as it contains a lot of cells: 8439 after filtration: Only the cells belonging to healthy patients are kept. We want to compare the predictive value of our set of gene with the state of the art reference matrix LM22, used by default by Cibersort*, a software that predicts cell proportion in RNAseq data based on the expression values. This way, only the celltype present in the data and our genesets and LM22 are kept : 
- T CD4 +
- T CD8 +
- Memory B cell
- Naive B cell
- Monocyte
- Natural killer cell

Neutrophils, macrophages and plasma cells being absent of the dataset, and platelet cells not in LM22 matrix.

Chen, Binbin, Michael S. Khodadoust, Chih Long Liu, Aaron M. Newman, et Ash A. Alizadeh. 2018. « Profiling tumor infiltrating immune cells with CIBERSORT ». Methods in molecular biology (Clifton, N.J.) 1711: 243‑59. https://doi.org/10.1007/978-1-4939-7493-1_12.


## Pipeline

All the mentionned script are available [here](/scripts/).


As shown in the previous diagram, for a given comparison, two series of permutations can be done if it remains more than 250 features after filtration with permutation variable importance.
It is also important to note that the filtration after permutation is not automated. It must be launched manually after each of those steps. 


<p align="center">
<img align="center" src="/doc/Tabula_sapiens_pipeline.drawio.png" height="750">
