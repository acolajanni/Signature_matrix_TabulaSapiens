# Signature matrix built from Tabula Sapiens single cell dataset

**Contact:**
Antonin Colajanni: antonin.colajanni@u-bordeaux.fr

## Data
THE TABULA SAPIENS CONSORTIUM. 2022. « The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans ». Science 376 (6594): eabl4896. https://doi.org/10.1126/science.abl4896.

**URL:** 
https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5

Tabula Sapiens - Immune 


## File Description: 

[Signature matrix constitution](/markdown/): document .html describing the filtration step for the data



## Method description 

<img src="/doc/Diapo_pathseq-TabulaSapiens.drawio.png" height="750">


## Pipeline

All the mentionned script are available [here](/scripts/).


As shown in the previous diagram, for a given comparison, two series of permutations can be done if it remains more than 250 features after filtration with permutation variable importance.
It is also important to note that the filtration after permutation is not automated. It must be launched manually after each of those steps. 

<img src="/doc/Diapo_pathseq-Tabula_sapiens_pipeline.drawio.png" height="750">
