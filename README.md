# Signature matrix built from Tabula Sapiens single cell dataset


**Be Carefull, the current signature matrix is built on 16487 cells, not 50115 as it should be, results will differ**


**Contact:**
Antonin Colajanni: antonin.colajanni@u-bordeaux.fr

## Data
THE TABULA SAPIENS CONSORTIUM. 2022. « The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans ». Science 376 (6594): eabl4896. https://doi.org/10.1126/science.abl4896.

**URL:** 
https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5 - *Tabula Sapiens - Immune* (filtered on blood tissues) 

Alternatively, the data *Tabula Sapiens - Blood* can be directly downloaded.


## Heatmap of the signature matrix

*to add*


## File Description: 

The [.html document](/markdown/feature_selection_v2.html) describes the filtration steps before and after the permutaions steps. The code to produce the figures is the [*feature_selection.R*](/scripts/Feature_selection.R) script (à ajouter quand il sera fini), in the same [folder](/scripts/) as the other scripts.


[Permutation results](/results) contains the .csv with the permutation importance of each variable for the 50 permutation. The obtained features after selection for each comparison, from the broader selection to the final one are available [here](/genesets).


## Method description 

<p align="center">
<img src="/doc/TabulaSapiens.drawio.png" height="500">


After the obtention of the normalised count table, a first step of feature engineering is done. First, for each celltype the mean expression of it and without the said celltype is computed. This way, a fold change between the background and the given cell can be calculated.

Then, the number of cell expressing each gene for each celltype is also computed in order to filter on the percentage of cell expressing a gene. 

With these two parameters, a filtration can be done. First if we compare one celltype vs all other, the filtration is done on the absolute value of the log2 Fold Change being superior to 1, meaning, a gene must be, at least twice as expressed in the given celltype compared to be background and conversly and it must be expressed in 25% of the cells of this given celltype. Alternatively, the curve of the absolute Fold change is generated, and all the genes with a fold change higher than the inflexion point of the curve* is selected, independantly of the proportion of cell expressing it. 

For the alternative, a gene can be a signature for a given celltype if it is expressed everywhere except in the said celltype. And, for the 25% threshold, it is supposed to be relatively generous to avoid selecting false negative gene. Moreover, for some celltypes, I aggregated several sub celltypes, it is the case for the T.CD4 cells for example, that corresponds to the aggregation of 
*CD4-positive, alpha-beta T cell*, *CD4-positive, alpha-beta memory T cell* and *naive thymus-derived CD4-positive, alpha-beta T cell* that are not equally represented.

*using the kneedle algorithm contained is the eponym package in R 


For the Celltype vs Celltype comparison, it is done between certain celltypes, known to be relatively similar, to help the model finding disciminative genes that could have been obscured by the comparison with the Fold change computed with the background mean expression. The same threshold the Fold change is used, however, it has been found that a threshold at 25% of cell expressing a gene was too stringent, so the choice has been made to select genes expressed in at least 2 cells in the given celltypes. 



## Pipeline

All the mentionned script are available [here](/scripts/).


As shown in the previous diagram, for a given comparison, two series of permutations can be done if it remains more than 250 features after filtration with permutation variable importance.
It is also important to note that the filtration after permutation is not automated. It must be launched manually after each of those steps. 

<p align="center">
<img align="center" src="/doc/Diapo_pathseq-Tabula_sapiens_pipeline.drawio.png" height="750">
