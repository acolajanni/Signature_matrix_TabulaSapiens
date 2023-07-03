# > June 2023                                                                                                              
# > Script : Process_matrix_50k.R                                                                                                      
# > Function : Computaing necessary steps for the feature filtration                  
# @ COLAJANNI Antonin                                                          
################################################################################

library(stringr)
library(UpSetR)
library(tidyverse)

################################################################################
# CLUSTER version
################################################################################
# P A T H

# distant
path = file.path("/shared/projects/microbiome_translocation/")
data_dir = file.path(path, "data/Tabula_sapiens_immune_all/")
result_dir = file.path(path, "results/Tabula_sapiens/")
setwd(path)

################################################################################
# F U N C T I O N S

create_dir = function(directory){
  if ( ! file.exists(directory)){
    dir.create(directory)  }
  return(directory) }

################################################################################
# I M P O R T
mat = read.csv(file.path(data_dir,"50000cell_immune_expressed.csv"), header=TRUE,sep='\t')#,nrows = 1000 )
cellTypes = readLines(file.path(data_dir,"50k_cell_labels.txt"))
mat$cell_id = cellTypes
mat$X = NULL
#row.names(mat) = cellTypes

original_labels = cellTypes

# Rename cells
Count_cell = as.data.frame(table(cellTypes))
Count_cell

table(cellTypes[str_detect(cellTypes, "NK")])
table(cellTypes[str_detect(cellTypes, "CD4")])
table(cellTypes[str_detect(cellTypes, "CD8")])
table(cellTypes[str_detect(cellTypes, "naive B")])
table(cellTypes[str_detect(cellTypes, "memory B")])
table(cellTypes[str_detect(cellTypes, "neutr")])
table(cellTypes[str_detect(cellTypes, "monocyte")])
table(cellTypes[str_detect(cellTypes, "plasma cell")])
table(cellTypes[str_detect(cellTypes, "erythrocyte")])
table(cellTypes[str_detect(cellTypes, "macrophage")])


reduced_CellTypes = ifelse(str_detect(cellTypes,"CD8" ), "CD8", cellTypes )
reduced_CellTypes = ifelse(str_detect(cellTypes,"CD4"), "CD4", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"NK"), "NK_cell", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"naive B"), "naive_B_cell", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"memory B"), "memory_B_cell", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"neutro"), "Neutrophil", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"monocyte"), "Monocyte", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"plasma cell"), "Plasmocyte", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"erythrocyte"), "Erythrocyte", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"macrophage"), "Macrophage", reduced_CellTypes)
reduced_CellTypes = ifelse(str_detect(cellTypes,"platelet"), "Platelet", reduced_CellTypes)

Count_cell = as.data.frame(table(reduced_CellTypes))

reduced_CellTypes = ifelse( ! reduced_CellTypes %in% c("NK_cell",'CD4','CD8',
                                                       "naive_B_cell" , "memory_B_cell" , 
                                                       "Erythrocyte","Neutrophil", "Monocyte",
                                                       "Plasmocyte","Macrophage","Platelet"  ), 
                            "other", reduced_CellTypes)

Count_cell = as.data.frame(table(reduced_CellTypes))

writeLines(reduced_CellTypes,file.path(data_dir,"Celltypes_of_interest_50k.txt"))

################################################################################
# V A R I A N C E   O N   M E A N   E X P R   G R O U P W I S E
mat = read.csv(file.path(data_dir,"50000cell_immune_expressed.csv"), header=TRUE,sep='\t')#,nrows = 1000 )
mat$X = NULL
cellTypes = readLines(file.path(data_dir,"Celltypes_of_interest_50k.txt"))
mat$cell_id = cellTypes


# Aggregate rows by group by using mean function
groupwise_mean = aggregate(mat[colnames(mat)[colnames(mat) != 'cell_id']], list(mat$cell_id), FUN=mean)


# Saving it in a txt file
groupwise_mean = as.data.frame(t(groupwise_mean))
colnames(groupwise_mean) = groupwise_mean[1,]
groupwise_mean = groupwise_mean[-1,]

write.table(groupwise_mean, file.path(data_dir,"50k_cells_geneset/GroupWise_mean_50k.txt"),quote=FALSE, sep="\t", row.names=TRUE)


################################################################################
# Mean vs all

reduced_CellTypes = y

means_list = list()
for (label in cells){
  print(label)
  mat$cell_id = ifelse(reduced_CellTypes != label, yes = paste("other_vs_",label), no=reduced_CellTypes)
  groupwise_mean = aggregate(mat[colnames(mat)[colnames(mat) != 'cell_id']], list(mat$cell_id), FUN=mean)
  groupwise_mean = as.data.frame(t(groupwise_mean))
  colnames(groupwise_mean) = groupwise_mean[1,]
  groupwise_mean = groupwise_mean[-1,]
  means_list[[label]] = groupwise_mean }

means_list = lapply(means_list, function(x){
  x$genes = row.names(x)
  return(x) })



groupwise_mean = reduce(means_list, full_join, by="genes")
col = colnames(groupwise_mean)[colnames(groupwise_mean)!='genes']
groupwise_mean = groupwise_mean[c('genes',col)]
groupwise_mean[col] = dplyr::mutate_all(groupwise_mean[col], 
                                        function(x) as.numeric(as.character(x) )) 


write.table(groupwise_mean, file.path(data_dir,"50k_cells_geneset/Groupwise_mean_vs_50k.txt"),quote=FALSE, sep="\t", 
            row.names=FALSE)

################################################################################
# count number of cells expressing each genes vs all + proportion 

# Counting 0 values for each cell type
y = readLines(paste0(data_dir,"Celltypes_of_interest_50k.txt"))
mat$cell_id = y

cells = unique(y)
cells = cells[! cells %in% c("other", "Erythrocyte")]

nonzero <- function(x) sum(x != 0)

count_nonzero = lapply(cells, function(x) {
  print(x)
  tmp_mat = mat[mat$cell_id == x , ]
  tmp_mat = as.data.frame(t(plyr::numcolwise(nonzero)(tmp_mat)))
  colnames(tmp_mat) = paste0(x, "_n_expr")
  return( tmp_mat )
})

names(count_nonzero) = cells
count_nonzero_df = as.data.frame(count_nonzero)

write.table(count_nonzero_df, file.path(data_dir,"50k_cells_geneset/Groupwise_nexpr_expressed.txt"),quote=FALSE, sep="\t", row.names=TRUE)


count_nonzero_df$genes = row.names(count_nonzero_df)
count_nonzero_df$genes = NULL

###  Computing non zero proportion
n_cell = as.data.frame(t(as.data.frame(table(y))))
colnames(n_cell) = n_cell[1,]
n_cell = n_cell[-c(1),]

count_nonzero_prop = lapply(colnames(count_nonzero_df), function(x) {
  x2 = str_remove(x, pattern = "_n_expr")
  print(x)
  return(count_nonzero_df[[x]] / as.numeric(n_cell[[x2]]) ) })

names(count_nonzero_prop) = paste0(str_remove(colnames(count_nonzero_df), pattern = "_n_expr"), '_prop_expr' )
#count_nonzero_prop$B_cell_prop_expr = NULL
count_nonzero_prop = as.data.frame(count_nonzero_prop)
row.names(count_nonzero_prop) = row.names(count_nonzero_df)

count_nonzero_prop = read.table(file.path(data_dir,"50k_cells_geneset/Groupwise_nprop_expressed.txt"),sep="\t", header = TRUE)


# comparing proportion of cells expressing: one cell vs all others

count_non_null_list = list()

non_null_count = as.data.frame(t(count_nonzero_df))
non_null_count$cell_id= row.names(non_null_count)


count_non_null_list = list()

for (label_cell in cells){
  print(label_cell)
  label = paste0(label_cell,"_n_expr")
  non_null_count$cell_id = row.names(non_null_count)
  non_null_count$cell_id = ifelse(non_null_count$cell_id != label, yes = paste("other_vs_",label), no=label)
  groupwise_sum = aggregate(non_null_count[colnames(non_null_count)[colnames(non_null_count) != 'cell_id']], list(non_null_count$cell_id), FUN=sum)
  groupwise_sum = as.data.frame(t(groupwise_sum))
  colnames(groupwise_sum) = groupwise_sum[1,]
  groupwise_sum = groupwise_sum[-1,]
  count_non_null_list[[label_cell]] = groupwise_sum
}

tmp = as.data.frame(count_non_null_list)
count_non_null = as.data.frame(apply(tmp, 2, function(x) as.numeric(as.character(x)) ))
row.names(count_non_null) = row.names(count_non_null_list$NK_cell)


tmp_names = strsplit(colnames(count_non_null), ".",T)

tmp_names = lapply(tmp_names, function(x) x[-1])
tmp_names = lapply(tmp_names, function(x) paste(x, collapse = '') )
tmp_names = unlist(tmp_names)
tmp_names = str_remove(tmp_names, "_n_expr")
colnames(count_non_null) = tmp_names


write.table(count_non_null, file.path(data_dir,"50k_cells_geneset/Groupwise_count_non_null_vs_50k.txt"),quote=FALSE, sep="\t", 
            row.names=FALSE)


n_cell=as.data.frame(table(y))
row.names(n_cell) = n_cell$y
n_cell$y = NULL
n_cell = t(n_cell)
n_cell = as.data.frame(n_cell)
n_cell$total = rowSums(n_cell)



other_vs_df = list()
other_vs_df = lapply(cells, function(label){
  new_col = paste0("other_vs_",label)
  other_vs_df[[new_col]] = n_cell$total - n_cell[[label]]
  return(other_vs_df)
})
other_vs_df = as.data.frame(t(as.data.frame(unlist(other_vs_df))))
n_cell = cbind(n_cell,other_vs_df)


count_non_null_plus1 = count_non_null
count_non_null_plus1[count_non_null_plus1 == 0] = 1


count_cell_prop = lapply(colnames(count_non_null_plus1), function(feature) {
  print(feature)
  tot = n_cell[[feature]]
  tmp[[feature]] = as.numeric(count_non_null_plus1[[feature]]) / tot
})
names(count_cell_prop) = colnames(count_non_null_plus1) 
count_cell_prop = as.data.frame(count_cell_prop)
row.names(count_cell_prop) = row.names(count_non_null)

write.table(count_cell_prop, file.path(data_dir,"50k_cells_geneset/Groupwise_count_non_null_vs_50k.txt"),quote=FALSE, sep="\t", 
            row.names=TRUE)
