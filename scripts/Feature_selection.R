# > June 2023                                                                                                              
# > Script : Feature_selection.R                                                                                                    
# > Function : Selecting and filtering features before and after the permutation step    
# > This script is a transcript of the Rmarkdown file in /markdown
# @ COLAJANNI Antonin                                                          
################################################################################
	
# I M P O R T	
library(stringr)	
library (plyr) 	
library(ggplot2)	
library(viridis)	
library(UpSetR)	
library(ggrepel)	
library(factoextra)	
library(FactoMineR)	
library(viridis)	
library(RColorBrewer)	
library(dplyr)	
library(pheatmap)	
library(DT) 	
library(cluster)	
library(grid)	
source("/home/acolajanni/Documents/work/GITHUB/Signature_matrix_TabulaSapiens/scripts/functions/functions_R.R")

# V A R I A B L E S	
path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation"	
local_path = "/home/acolajanni/Documents/work/"	
#setwd(main.dir)	
	
#path = "/home/acolajanni/Documents/work/"	
data_dir = file.path(path, "data/Tabula_sapiens_immune_all")	
	
result_dir = file.path(path, "results/Tabula_sapiens/RF_filter_50k/")	
#setwd(path)	
	
git_dir = file.path(local_path,"GITHUB/Signature_matrix_TabulaSapiens/")	
git_geneset = file.path(git_dir,"genesets/")	




#' 	
#' 	
#' # Cell Types	
#' 	
#' 	
	
cellTypes = read.table(file.path(data_dir,"50k_cell_labels.txt"), header=FALSE, sep='\t')	
label_count = cellTypes$V1	
label_count = as.data.frame(table(label_count) )	
	
	
ggplot( label_count, aes(x = reorder(label_count,Freq), y = Freq, fill=Freq ) )  +	
  geom_col(position = position_dodge(0), width = 0.75) +	
    geom_text(color="black",hjust=-0.2, size=3.8, angle=0,	
            aes(y = ,label = Freq))+	
  coord_flip(ylim = c(0,12000)  ) +	
  labs(fill = "Geneset") +	
  ylab('Frequency') + xlab('Cell type')+	
  ggtitle(paste0("Frequency of cell types in Tabula Sapiens immune (Blood) dataset"))+	
  scale_fill_viridis(discrete = FALSE, alpha=1, direction = -1, option='E', name="Absolute frequency") +	
  theme_linedraw()+  	
  theme( plot.title = element_text(size=15),	
         axis.text.y = element_text(size = 12),	
         axis.text.x = element_text(size = 12))	
	
#' 	
#' 	
#' 	
	
cellTypes = read.table(file.path(data_dir,"Celltypes_of_interest_50k.txt"), header=FALSE, sep='\t')	
label_count = cellTypes$V1	
label_count = as.data.frame(table(label_count) )	
	
ggplot( label_count, aes(x = reorder(label_count,Freq), y = Freq, fill=Freq ) )  +	
  geom_col(position = position_dodge(0), width = 0.75) +	
    geom_text(color="black",hjust=-0.2, size=3.8, angle=0,	
            aes(y = ,label = Freq))+	
  coord_flip(ylim = c(0,17000)  ) +	
  labs(fill = "Geneset") +	
  ylab('Frequency') + xlab('Cell type')+	
  ggtitle(paste0("Frequency of cell types in Tabula Sapiens immune (Blood) dataset"))+	
  scale_fill_viridis(discrete = FALSE, alpha=1, direction = -1, option='E', name="Absolute frequency") +	
  theme_linedraw()+  	
  theme( plot.title = element_text(size=15),	
         axis.text.y = element_text(size = 12),	
         axis.text.x = element_text(size = 12))	
	
	
	
	
n_cells = label_count	
row.names(n_cells) = label_count$label_count	
n_cells$label_count = NULL	
n_cells = as.data.frame(t(n_cells))	
n_cells$total = sum(n_cells[1,])	
#' 	
#' 	
#' 	
data_dir_tmp = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/data/scRNAseq/Liu_2021_cell/"	
	
cellTypes = readLines(file.path(data_dir_tmp,"Liu2021_index.txt"))	
label_count = as.data.frame(table(cellTypes) )	
	
	
ggplot( label_count, aes(x = reorder(cellTypes,Freq), y = Freq, fill=Freq ) )  +	
  geom_col(position = position_dodge(0), width = 0.75) +	
    geom_text(color="black",hjust=-0.2, size=3.8, angle=0,	
            aes(y = ,label = Freq))+	
  coord_flip(ylim = c(0,40000)  ) +	
  labs(fill = "Geneset") +	
  ylab('Frequency') + xlab('Cell type')+	
  ggtitle(paste0("Frequency of cell types in Liu et al. 2021, Cell. Healthy patients, Innate and Adaptative  "))+	
  scale_fill_viridis(discrete = FALSE, alpha=1, direction = -1, option='E', name="Absolute frequency") +	
  theme_linedraw()+  	
  theme( plot.title = element_text(size=15),	
         axis.text.y = element_text(size = 12),	
         axis.text.x = element_text(size = 12))	
	
#' 	
#' 	
#' 	
#' 	
data_dir_tmp = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/data/scRNAseq/Liu_2021_cell/"	
	
cellTypes_innate = read.table(file.path(data_dir_tmp,"Innate_healthy_labels.txt"), header=FALSE, sep='\t')	
label_count = cellTypes_innate$V1	
label_count = as.data.frame(table(label_count) )	
colnames(label_count) = c("cells","Freq")	
	
	
cellTypes_adaptative = read.table(file.path(data_dir_tmp,"Adaptative_healthy_labels.txt"), header=FALSE, sep='\t')	
label_count2 = cellTypes_adaptative$V1	
label_count2 = as.data.frame(table(label_count2) )	
colnames(label_count2) = c("cells","Freq")	
	
label_count = rbind(label_count,label_count2)	
	
ggplot( label_count, aes(x = reorder(cells,Freq), y = Freq, fill=Freq ) )  +	
  geom_col(position = position_dodge(0), width = 0.75) +	
    geom_text(color="black",hjust=-0.2, size=3.8, angle=0,	
            aes(y = ,label = Freq))+	
  coord_flip(ylim = c(0,40000)  ) +	
  labs(fill = "Geneset") +	
  ylab('Frequency') + xlab('Cell type')+	
  ggtitle(paste0("Frequency of cell types in Liu et al. 2021, Cell. Healthy patients, Innate and Adaptative  "))+	
  scale_fill_viridis(discrete = FALSE, alpha=1, direction = -1, option='E', name="Absolute frequency") +	
  theme_linedraw()+	
  theme( plot.title = element_text(size=15),	
         axis.text.y = element_text(size = 12),	
         axis.text.x = element_text(size = 12))	
	
#' 	
#' 	
#' 	
#' # Filtering genes - Cell vs Cell {.tabset}	
#' 	
#' 	
data_dir = file.path(data_dir,"50k_cells_geneset")	
	
groupwise_mean = read.table(file.path( data_dir,"Groupwise_count_non_null_vs_50k.txt"), sep='\t', header = TRUE)	
groupwise_prop_full = read.table(file.path( data_dir,"Groupwise_count_non_null_vs_50k.txt"), sep='\t', header = TRUE)	
	
	
gene_label = readLines(file.path(data_dir,"50k_cell_genes.txt"))	
groupwise_mean$genes = gene_label	
groupwise_prop_full$genes = gene_label	
row.names(groupwise_mean) = gene_label	
row.names(groupwise_prop_full) = gene_label	
	
	
	
# Compute median expr over all datatable	
col = colnames(groupwise_mean)[!str_detect(colnames(groupwise_mean), "vs|genes" )]	
expr_median = compute_median_expr(groupwise_mean,col )	
	
	
colnames(groupwise_prop_full) = paste0(colnames(groupwise_prop_full),"_prop_expr")	
groupwise_prop_full$genes = row.names(groupwise_prop_full)	
groupwise_mean$genes = row.names(groupwise_mean)	
  	
groupwise_mean = merge(groupwise_mean, groupwise_prop_full, by='genes')	
groupwise_mean$genes_prop_expr = NULL	
	
#' 	
#' 	
#' ## CD4 vs CD8	
#' 	
#' **61 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_cd48 = groupwise_mean[str_detect(colnames(groupwise_mean),"CD4|CD8")]	
FC_df_cd48$FC_cd48 = log2( (FC_df_cd48$CD4+expr_median ) / (FC_df_cd48$CD8+expr_median) )	
	
elbow=kneedle_threshold(FC_df_cd48, FC_df_cd48$FC_cd48)	
elbow$plot	
cd48 = Rprop_plot(FC_df_cd48, FC_x = "FC_cd48", y="CD4_prop_expr", 	
                title = "Log Fold Change of CD4 vs CD8 depending on the proportion of CD4 cells expressing the genes", 	
                #add_label = TRUE, 	
                FC_threshold = 1, 	
                prop_threshold = 2/n_cells$CD4, 	
                prop_threshold_to_show = paste0("2 / ",as.character(n_cells$CD4)),	
                highlight = "NCR1",	
                FC_thresh2 = elbow$threshold) 	
	
cd84 = Rprop_plot(FC_df_cd48, FC_x = "FC_cd48", y="CD8_prop_expr", 	
                title = "Log Fold Change of CD4 vs CD8 depending on the proportion of CD8 cells expressing the genes", 	
                #add_label = TRUE, 	
                FC_threshold = 1, 	
                prop_threshold = 2/n_cells$CD8,	
                prop_threshold_to_show = paste0("2 / ",as.character(n_cells$CD8)),	
                highlight = "NCR1",	
                FC_thresh2 = elbow$threshold) 	
	
	
cowplot::plot_grid(cd84$plot, cd48$plot, labels = 'AUTO', nrow = 1)	
	
	
	
cd84_new = row.names(cd84$df[cd84$df$selection != "Not selected",])	
cd48_new = row.names(cd48$df[cd48$df$selection != "Not selected",])	
	
	
cd48_genes = unique(c(cd84_new,cd48_new))	
	
writeLines(cd48_genes, file.path(data_dir,'CD4vsCD8_3435.txt') )	
writeLines(cd48_genes, file.path(git_geneset,'CD4vsCD8_3435.txt') )	
	
#' 	
#' 	
#' ### 3435 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "CD4-CD8/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_cd48 = Permutation_list$`RF_permutation_First_CD4 vs CD8`	
other_cd48$genes = row.names(other_cd48)	
	
cd48_computation = importance_barplot(other_cd48, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 3435 genes selected in CD4 vs CD8 classification")	
                                        	
cd48_computation$plot	
	
top61_cd48_genes = cd48_computation$ordered_permutation_dataframe[cd48_computation$ordered_permutation_dataframe$mean_importance >0,]$genes 	
	
writeLines(top61_cd48_genes, file.path(data_dir,'CD4vsCD8_61.txt') )	
writeLines(top61_cd48_genes, file.path(git_geneset,'CD4vsCD8_61.txt') )	
	
#' 	
#' 	
#' ## NK vs CD8	
#' 	
#' **79 Features selected**	
#' 	
#' 	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_NKcd8 = groupwise_mean[str_detect(colnames(groupwise_mean),"NK|CD8")]	
FC_df_NKcd8$FC_NKCD8 = log2( (FC_df_NKcd8$NK_cell+expr_median ) / (FC_df_NKcd8$CD8+expr_median) )	
FC_df_NKcd8$Avalue = 0.5* ( log2(FC_df_NKcd8$NK_cell+1) + log2(FC_df_NKcd8$CD8+1) )	
	
elbow=kneedle_threshold(FC_df_NKcd8, FC_df_NKcd8$FC_NKCD8)	
elbow$plot	
	
nkcd8 = Rprop_plot(FC_df_NKcd8, FC_x = "FC_NKCD8", y="NK_cell_prop_expr", 	
                title = "Log Fold Change of NK vs CD8 depending on the proportion of NK cells expressing the genes", 	
                #add_label = TRUE, 	
                FC_threshold = 1, 	
                prop_threshold = 2/n_cells$NK_cell,	
                prop_threshold_to_show = paste0("2 / ",as.character(n_cells$NK_cell)),	
                #highlight = test,	
                FC_thresh2 = elbow$threshold) 	
	
cd8nk = Rprop_plot(FC_df_NKcd8, FC_x = "FC_NKCD8", y="CD8_prop_expr", 	
                title = "Log Fold Change of NK vs CD8 depending on the proportion of CD8 cells expressing the genes", 	
                #add_label = TRUE, 	
                FC_threshold = 1, 	
                prop_threshold = 2/n_cells$CD8,	
                prop_threshold_to_show = paste0("2 / ",as.character(n_cells$CD8)),	
                #highlight = test,	
                FC_thresh2 = elbow$threshold) 	
	
	
	
cowplot::plot_grid(nkcd8$plot, cd8nk$plot, labels = 'AUTO', nrow = 1)	
	
nkcd8_new = row.names(nkcd8$df[nkcd8$df$selection != "Not selected",])	
cd8nk_new = row.names(cd8nk$df[cd8nk$df$selection != "Not selected",])	
	
cd8nk_genes = unique(c(nkcd8_new,cd8nk_new))	
	
writeLines(cd8nk_genes, file.path(data_dir,'CD8vsNK_1028.txt') )	
writeLines(cd8nk_genes, file.path(git_geneset,'CD8vsNK_1028.txt') )	
	
#' 	
#' 	
#' ### 1028 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "NK-CD8/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
	
nkcd8 = Permutation_list$`RF_permutation_First_NK vs CD8`	
nkcd8$genes = row.names(nkcd8)	
	
	
nkcd8_computation = importance_barplot(nkcd8, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 1028 genes selected in NK vs CD8 classification")	
                                        	
nkcd8_computation$plot+	
   theme(axis.text.y=element_blank(),	
        axis.ticks.y=element_blank()	
        )	
	
top247_nkcd8_genes = nkcd8_computation$ordered_permutation_dataframe[nkcd8_computation$ordered_permutation_dataframe$mean_importance >0,]$genes 	
writeLines(top247_nkcd8_genes, file.path(data_dir,'NKvsCD8_247.txt') )	
writeLines(top247_nkcd8_genes, file.path(git_geneset,'NKvsCD8_247.txt') )	
	
#' 	
#' 	
#' ### 247 features	
#' 	
#' 	
	
	
Permutation_list$`RF_permutation_First_NK vs CD8`= NULL	
nkcd8_2 = lapply(Permutation_list, function(x) {return(x[1:5]) } ) 	
nkcd8_2 = rlist::list.cbind(nkcd8_2)	
mean_imp = rowMeans(nkcd8_2)	
nkcd8_2$genes = row.names(nkcd8_2) 	
row_stdev <- apply(nkcd8_2, 1, sd, na.rm=TRUE)	
	
nkcd8_2$mean_importance = mean_imp	
nkcd8_2$std_importance = row_stdev	
	
	
nkcd8_computation = importance_barplot(nkcd8_2, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 247 genes selected in NK vs CD8 classification")	
                                        	
nkcd8_computation$plot	
	
top79_nkcd8_genes = nkcd8_computation$ordered_permutation_dataframe[nkcd8_computation$ordered_permutation_dataframe$mean_importance >0,]$genes 	
writeLines(top79_nkcd8_genes, file.path(data_dir,'NKvsCD8_79.txt') )	
writeLines(top79_nkcd8_genes, file.path(git_geneset,'NKvsCD8_79.txt') )	
	
#' 	
#' 	
#' 	
#' ## Naive vs Memory Bcells	
#' 	
#' **42 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_naiveMem = groupwise_mean[str_detect(colnames(groupwise_mean),"naive|memory")]	
	
FC_df_naiveMem$FC_nmB = log2( 	
  (FC_df_naiveMem$memory_B_cell+expr_median ) / (FC_df_naiveMem$naive_B_cell+expr_median) )	
	
elbow=kneedle_threshold(FC_df_naiveMem, FC_df_naiveMem$FC_nmB)	
elbow$plot	
	
naivemem = Rprop_plot(FC_df_naiveMem, FC_x = "FC_nmB", y="naive_B_cell_prop_expr", 	
                title = "Log Fold Change of Naive vs Memory Bcells depending on the proportion of naive Bcells expressing the genes", 	
                #add_label = TRUE, 	
                FC_threshold = 1, 	
                prop_threshold = 2/n_cells$naive_B_cell,	
                prop_threshold_to_show = paste0("2 / ",as.character(n_cells$naive_B_cell)),	
                FC_thresh2 = elbow$threshold) 	
	
memnaive = Rprop_plot(FC_df_naiveMem, FC_x = "FC_nmB", y="memory_B_cell_prop_expr", 	
                title = "Log Fold Change of Naive vs Memory Bcells depending on the proportion of memory Bcells expressing the genes", 	
                #add_label = TRUE, 	
                FC_threshold = 1, 	
                prop_threshold = 2/n_cells$memory_B_cell, 	
                prop_threshold_to_show = paste0("2 / ",as.character(n_cells$memory_B_cell)),	
                FC_thresh2 = elbow$threshold) 	
	
	
cowplot::plot_grid(naivemem$plot, memnaive$plot, labels = 'AUTO', nrow = 1)	
	
	
memnaive_new = row.names(memnaive$df[memnaive$df$selection != "Not selected",])	
naivemem_new = row.names(naivemem$df[naivemem$df$selection != "Not selected",])	
nmBcell_genes = unique(c(memnaive_new,naivemem_new))	
	
	
writeLines(nmBcell_genes, file.path(data_dir,'naiveVSmemory_B_3582.txt') )	
writeLines(nmBcell_genes, file.path(git_geneset,'naiveVSmemory_B_3582.txt') )	
	
#' 	
#' 	
#' ### 3582 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Memory-Naive_Bcell/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_nm = Permutation_list$`RF_permutation_First_Naive vs Memory`	
other_nm$genes = row.names(other_nm)	
	
nm_computation = importance_barplot(other_nm, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 3582 genes selected in Memory vs Naive Bcells classification")	
                                        	
nm_computation$plot	
	
top42_naive_mem = nm_computation$ordered_permutation_dataframe[nm_computation$ordered_permutation_dataframe$mean_importance >0,]$genes 	
	
writeLines(top42_naive_mem, file.path(data_dir,'naiveVSmemory_B_42.txt') )	
writeLines(top42_naive_mem, file.path(git_geneset,'naiveVSmemory_B_42.txt') )	
	
#' 	
#' # Filtering genes - Cell vs all {.tabset}	
#' 	
#' ## NK vs other	
#' 	
#' **178 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_NK = groupwise_mean[str_detect(colnames(groupwise_mean),"NK")]	
FC_df_NK = Get_value_RA(FC_df_NK, expr_median)	
FC_df_NK = FC_prop_expr(FC_df_NK)	
	
	
nk = RA_plot(FC_df_NK, "Avalue_NK_cell", "FC_NK_cell", title = "RA plot NK vs other", add_label = FALSE)	
	
	
	
elbow=kneedle_threshold(FC_df_NK, FC_df_NK$FC_NK_cell)	
elbow$plot	
	
nk = Rprop_plot(FC_df_NK,FC_x = "FC_NK_cell", y="NK_cell_prop_expr", 	
                title = "Log Fold Change of NK vs other depending on the proportion of NK cells expressing the genes", 	
                highlight = "KLRB1",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
	
nk$plot	
	
	
nk_new = row.names(nk$df[nk$df$selection != "Not selected",])	
	
	
writeLines(nk_new, file.path(data_dir,'NK_926.txt') )	
writeLines(nk_new, file.path(git_geneset,'NK_926.txt') )	
	
#' 	
#' 	
#' ### 926 features	
#' 	
#' 	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "NK/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_NK = Permutation_list$`RF_permutation_First_NK vs other`	
other_NK$genes = row.names(other_NK)	
	
NK_computation = importance_barplot(other_NK, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 926 genes selected in Other vs NK cell classification", 	
                                    conf_interval = TRUE) 	
                                        	
NK_computation$plot  +	
   theme(axis.text.y=element_blank(),	
        axis.ticks.y=element_blank()	
        )	
	
top178_nk = NK_computation$ordered_permutation_dataframe[NK_computation$ordered_permutation_dataframe$mean_importance >0,]$genes 	
writeLines(top178_nk, file.path(data_dir,'NK_178.txt') )	
writeLines(top178_nk, file.path(git_geneset,'NK_178.txt') )	
	

#' 	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/NK_cell/" ),	
  til10_path = file.path(path_heatmap, "TIL10/NK_cell/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/NK_cell/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","NK_cell"))	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/NK_cell/" ),	
  til10_path = file.path(path_heatmap, "TIL10/NK_cell/" ),	
  lm22_path = file.path(path_heatmap, "LM22/NK_cell/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","NK_cell"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
cf = rbind(cf_tbs,cf_Liu2021)	
	
comp = build_confusion_heatmap(cf,  	
                               c("Other_cells","NK_cell"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - NK cell Classification") 	

#' 	
#' ## CD8 vs other	
#' 	
#' **79 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_CD8 = groupwise_mean[str_detect(colnames(groupwise_mean),"CD8")]	
FC_df_CD8 = Get_value_RA(FC_df_CD8, expr_median)	
FC_df_CD8 = FC_prop_expr(FC_df_CD8)	
	
	
elbow=kneedle_threshold(FC_df_CD8, FC_df_CD8$FC_CD8) 	
elbow$plot	
	
cd8 = Rprop_plot(FC_df_CD8,FC_x = "FC_CD8", y="CD8_prop_expr", 	
                title = "Log Fold Change of CD8 vs other depending on the proportion of CD8 cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
cd8$plot	
	
	
cd8_new = row.names(cd8$df[cd8$df$selection != "Not selected",])	
	
	
writeLines(cd8_new, file.path(data_dir,'CD8_639.txt') )	
writeLines(cd8_new, file.path(git_geneset,'CD8_639.txt') )	
	
#' 	
#' 	
#' ### 639 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "CD8/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
	
other_CD8 = Permutation_list$`RF_permutation_First_CD8 vs other`	
other_CD8$genes = row.names(other_CD8)	
	
cd8_computation = importance_barplot(other_CD8, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 639 genes selected in Other vs CD8 classification")	
                                        	
cd8_computation$plot+ 	
  theme(axis.text.y = element_text(size=10))	
	
	
top79_cd8 = cd8_computation$ordered_permutation_dataframe[cd8_computation$ordered_permutation_dataframe$mean_importance >0,]$genes 	
writeLines(top79_cd8, file.path(data_dir,'CD8_79.txt') )	
writeLines(top79_cd8, file.path(git_geneset,'CD8_79.txt') )	
	

#' 	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/CD8/" ),	
  til10_path = file.path(path_heatmap, "TIL10/CD8/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/CD8/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","CD8"))	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/CD8/" ),	
  til10_path = file.path(path_heatmap, "TIL10/CD8/" ),	
  lm22_path = file.path(path_heatmap, "LM22/CD8/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","CD8"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
cf = rbind(cf_tbs,cf_Liu2021)	
	
comp = build_confusion_heatmap(cf,  	
                               c("Other_cells","CD8"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - CD8 Classification") 	
#' 	
#' 	
#' ## CD4 vs other	
#' 	
#' **114 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_CD4 = groupwise_mean[str_detect(colnames(groupwise_mean),"CD4")]	
FC_df_CD4 = Get_value_RA(FC_df_CD4, expr_median)	
FC_df_CD4 = FC_prop_expr(FC_df_CD4)	
	
	
elbow=kneedle_threshold(FC_df_CD4, FC_df_CD4$FC_CD4) 	
elbow$plot	
	
cd4 = Rprop_plot(FC_df_CD4,FC_x = "FC_CD4", y="CD4_prop_expr", 	
                title = "Log Fold Change of CD4 vs other depending on the proportion of CD4 cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
cd4$plot	
	
	
cd4_new = row.names(cd4$df[cd4$df$selection != "Not selected",])	
	
	
writeLines(cd4_new, file.path(data_dir,'CD4_1402.txt') )	
writeLines(cd4_new, file.path(git_geneset,'CD4_1402.txt') )	
	
#' 	
#' 	
#' ### 1402 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "CD4/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_CD4 = Permutation_list$`RF_permutation_First_CD4 vs other`	
other_CD4$genes = row.names(other_CD4)	
	
cd4_computation = importance_barplot(other_CD4, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 1402 genes selected in Other vs CD4 classification")	
                                        	
cd4_computation$plot+ 	
  theme(axis.text.y = element_text(size=0))	
	
	
CD4_top319 = cd4_computation$ordered_permutation_dataframe[cd4_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
	
writeLines(CD4_top319, file.path(data_dir,'CD4_319.txt') )	
writeLines(CD4_top319, file.path(git_geneset,'CD4_319.txt') )	
	
#' 	
#' 	
#' ### 319 features	
#' 	
#' 	
	
other_CD4 = Permutation_list$`RF_permutation_Second_CD4 vs other`	
other_CD4$genes = row.names(other_CD4)	
	
cd4_computation = importance_barplot(other_CD4, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 319 genes selected in Other vs CD4 classification")	
                                        	
cd4_computation$plot+ 	
  theme(axis.text.y = element_text(size=8))	
	
	
CD4_top114 = cd4_computation$ordered_permutation_dataframe[cd4_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
writeLines(CD4_top114, file.path(data_dir,'CD4_114.txt') )	
writeLines(CD4_top114, file.path(git_geneset,'CD4_114.txt') )	
	
#' 	
#' 	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/CD4/" ),	
  til10_path = file.path(path_heatmap, "TIL10/CD4/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/CD4/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","CD4"))	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/CD4/" ),	
  til10_path = file.path(path_heatmap, "TIL10/CD4/" ),	
  lm22_path = file.path(path_heatmap, "LM22/CD4/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","CD4"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
cf = rbind(cf_tbs,cf_Liu2021)	
	
comp = build_confusion_heatmap(cf,  	
                               c("Other_cells","CD4"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - CD4 Classification") 	
#' 	
#' 	
#' ## Naive Bcell vs other	
#' 	
#' **29 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_naiveB = groupwise_mean[str_detect(colnames(groupwise_mean),"naive")]	
FC_df_naiveB = Get_value_RA(FC_df_naiveB, expr_median)	
FC_df_naiveB = FC_prop_expr(FC_df_naiveB)	
	
	
elbow=kneedle_threshold(FC_df_naiveB, FC_df_naiveB$naive_B_cell, min = 1) 	
	
elbow$plot	
	
naiveB = Rprop_plot(FC_df_naiveB,FC_x = "FC_naive_B_cell", y="naive_B_cell_prop_expr", 	
                title = "Log Fold Change of Naive B cell vs other depending on the proportion of CD8 cells expressing the genes",	
                add_label = FALSE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
naiveB$plot	
	
	
naiveB_genes = row.names(naiveB$df[naiveB$df$selection != "Not selected",])	
	
writeLines(naiveB_genes, file.path(data_dir,'naiveB_6098.txt') )	
writeLines(naiveB_genes, file.path(git_geneset,'naiveB_6098.txt') )	
	
#' 	
#' 	
#' ### 6098 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Naive_Bcell/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
	
other_naiveB = Permutation_list$`RF_permutation_First_Naive_Bcell vs other`	
other_naiveB$genes = row.names(other_naiveB)	
	
naiveB_computation = importance_barplot(other_naiveB, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 6098 genes selected in Other vs naive Bcell classification",	
                                        conf_interval = FALSE)	
                                        	
naiveB_computation$plot + 	
  theme(axis.text.y = element_text(size=10))	
	
naiveB_top29 = naiveB_computation$ordered_permutation_dataframe[naiveB_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
writeLines(naiveB_top29, file.path(data_dir,'naiveB_29.txt') )	
writeLines(naiveB_top29, file.path(git_geneset,'naiveB_29.txt') )	
	
#' 	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/naive_B_cell/" ),	
  til10_path = file.path(path_heatmap, "TIL10/naive_B_cell/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/naive_B_cell/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","naive_B_cell"))	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/naive_B_cell/" ),	
  til10_path = file.path(path_heatmap, "TIL10/naive_B_cell/" ),	
  lm22_path = file.path(path_heatmap, "LM22/naive_B_cell/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","naive_B_cell"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
cf = rbind(cf_tbs,cf_Liu2021)	
	
comp = build_confusion_heatmap(cf,  	
                               c("Other_cells","naive_B_cell"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - naive B cell Classification") 	
#' 	
#' 	
#' ## Memory Bcell vs other	
#' 	
#' **53 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_memB = groupwise_mean[str_detect(colnames(groupwise_mean),"memory")]	
FC_df_memB = Get_value_RA(FC_df_memB, expr_median)	
FC_df_memB = FC_prop_expr(FC_df_memB)	
	
	
elbow=kneedle_threshold(FC_df_memB, FC_df_memB$memory_B_cell,min = 1) 	
elbow$plot	
memB = Rprop_plot(FC_df_memB,FC_x = "FC_memory_B_cell", y="memory_B_cell_prop_expr", 	
                title = "Log Fold Change of Memory B cell vs other depending on the proportion of CD8 cells expressing the genes",	
                highlight = "NCR1",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
memB$plot	
	
memB_genes = row.names(memB$df[memB$df$selection != "Not selected",])	
	
writeLines(memB_genes, file.path(data_dir,'memoryB_7208.txt') )	
writeLines(memB_genes, file.path(git_geneset,'memoryB_7208.txt') )	
	

#' 	
#' 	
#' ### 7208 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Memory_Bcell/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	

other_memB = Permutation_list$`RF_permutation_First_Memory_Bcell vs other`	
other_memB$Genes = row.names(other_memB) 	
  	
memB_computation = importance_barplot(other_memB, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 7208 genes selected in Other vs nnaive Bcell classification")	
                                        	
memB_computation$plot + 	
  theme(axis.text.y = element_text(size=9))	
	
top53_memB = memB_computation$ordered_permutation_dataframe[memB_computation$ordered_permutation_dataframe$mean_importance > 0 ,]$genes	
	
writeLines(top53_memB, file.path(data_dir,'memoryB_53.txt') )	
writeLines(top53_memB, file.path(git_geneset,'memoryB_53.txt') )	
	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/memory_B_cell/" ),	
  til10_path = file.path(path_heatmap, "TIL10/memory_B_cell/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/memory_B_cell/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","memory_B_cell"))	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/memory_B_cell/" ),	
  til10_path = file.path(path_heatmap, "TIL10/memory_B_cell/" ),	
  lm22_path = file.path(path_heatmap, "LM22/memory_B_cell/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","memory_B_cell"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
cf = rbind(cf_tbs,cf_Liu2021)	
	
comp = build_confusion_heatmap(cf,  	
                               c("Other_cells","memory_B_cell"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - memory B cell Classification") 	
#' 	
#' 	
#' 	
#' ## Neutrophils vs other	
#' 	
#' **84 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_neutro = groupwise_mean[str_detect(colnames(groupwise_mean),"Neutro")]	
FC_df_neutro = Get_value_RA(FC_df_neutro, expr_median)	
FC_df_neutro = FC_prop_expr(FC_df_neutro)	
	
	
elbow=kneedle_threshold(FC_df_neutro, FC_df_neutro$FC_Neutrophil) 	
elbow$plot	
	
neutro = Rprop_plot(FC_df_neutro,FC_x = "FC_Neutrophil", y="Neutrophil_prop_expr", 	
                title = "Log Fold Change of Neutrophil vs other depending on the proportion of Neutrophil cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
neutro$plot	
	
	
neutro_new = row.names(neutro$df[neutro$df$selection != "Not selected",])	
	
	
writeLines(neutro_new, file.path(data_dir,'Neutro_825.txt') )	
writeLines(neutro_new, file.path(git_geneset,'Neutro_825.txt') )	
	
#' 	
#' 	
#' ### 825 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Neutrophil/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_neutro = Permutation_list$`RF_permutation_First_Neutrophil vs other`	
other_neutro$genes = row.names(other_neutro)	
	
neutro_computation = importance_barplot(other_neutro, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 825 genes selected in Other vs neutrophil classification")	
                                        	
neutro_computation$plot + 	
  theme(axis.text.y = element_text(size=10))	
	
neutro_top84 = neutro_computation$ordered_permutation_dataframe[neutro_computation$ordered_permutation_dataframe$mean_importance>0,]$genes	
writeLines(neutro_top84, file.path(data_dir,'Neutro_84.txt') )	
writeLines(neutro_top84, file.path(git_geneset,'Neutro_84.txt') )	
	
#' 	
#' 	
#' 	
#' ### Performances - No neutrophils in Liu2021 data	
#' 	
#' 	
	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/Neutrophil/" ),	
  til10_path = file.path(path_heatmap, "TIL10/Neutrophil/" ),	
  lm22_path = file.path(path_heatmap, "LM22/Neutrophil/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","Neutrophil"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
	
comp = build_confusion_heatmap(cf_tbs,  	
                               c("Other_cells","Neutrophil"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Neutrophils Classification") 	
#' 	
#' 	
#' ## Plasmocyte vs other	
#' 	
#' **19 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_Plasmo = groupwise_mean[str_detect(colnames(groupwise_mean),"Plasmo")]	
FC_df_Plasmo = Get_value_RA(FC_df_Plasmo, expr_median)	
FC_df_Plasmo = FC_prop_expr(FC_df_Plasmo)	
	
	
elbow=kneedle_threshold(FC_df_Plasmo, FC_df_Plasmo$FC_Plasmocyte) 	
elbow$plot	
	
plasmo = Rprop_plot(FC_df_Plasmo,FC_x = "FC_Plasmocyte", y="Plasmocyte_prop_expr", 	
                title = "Log Fold Change of Plasmocyte vs other depending on the proportion of Plasmocyte cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
plasmo$plot	
	
	
plasmo_new = row.names(plasmo$df[plasmo$df$selection != "Not selected",])	
	
	
writeLines(plasmo_new, file.path(data_dir,'Plasmocyte_5172.txt') )	
writeLines(plasmo_new, file.path(git_geneset,'Plasmocyte_5172.txt') )	
	
#' 	
#' 	
#' ### 5172 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Plasmocyte/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
	
other_plasmo = Permutation_list$`RF_permutation_First_Plasmocyte vs other`	
	
plasmo_computation = importance_barplot(other_plasmo, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 5172 genes selected in Other vs Plasmocyte classification")	
                                        	
plasmo_computation$plot + 	
  theme(axis.text.y = element_text(size=10))	
	
plasmo_top19 = plasmo_computation$ordered_permutation_dataframe[plasmo_computation$ordered_permutation_dataframe$mean_importance>0,]$genes	
writeLines(plasmo_top19, file.path(data_dir,'Plasmocyte_19.txt') )	
writeLines(plasmo_top19, file.path(git_geneset,'Plasmocyte_19.txt') )	
	
#' 	
#' 	
#' ### Performances - No Plasmocytes in Liu2021 data / TIL10	
#' 	
#' 	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf1 = confusion_matrix_setup(	
    directory = file.path(path_heatmap, "Tabula_Sapiens_50k/Plasmocyte/" ),	
    pattern2= "/confusion_heatmap_",	
    cells = c("Other_cells","Plasmocyte"),	
    signature_matrix = "Tabula Sapiens")	
	
cf2 = confusion_matrix_setup(	
    directory = file.path(path_heatmap, "LM22/Plasmocyte/" ),	
    pattern2= "/confusion_heatmap_",	
    cells = c("Other_cells","Plasmocyte"),	
    signature_matrix = "LM22")	
	
cf_tbs = rbind(cf1,cf2)	
cf_tbs$dataset = "Tabula Sapiens data"	
	
comp = build_confusion_heatmap(cf_tbs,  	
                               c("Other_cells","Plasmocyte"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Plasmocytes Classification") 	
#' 	
#' 	
#' 	
#' ## Monocyte vs other	
#' 	
#' **174 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_mono = groupwise_mean[str_detect(colnames(groupwise_mean),"Mono")]	
FC_df_mono = Get_value_RA(FC_df_mono, expr_median)	
FC_df_mono = FC_prop_expr(FC_df_mono)	
	
	
elbow=kneedle_threshold(FC_df_mono, FC_df_mono$FC_Monocyte) 	
elbow$plot	
	
mono = Rprop_plot(FC_df_mono,FC_x = "FC_Monocyte", y="Monocyte_prop_expr", 	
                title = "Log Fold Change of Monocyte vs other depending on the proportion of Monocyte cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
mono$plot	
	
	
mono_new = row.names(mono$df[mono$df$selection != "Not selected",])	
	
	
writeLines(mono_new, file.path(data_dir,'Monocyte_4022.txt') )	
writeLines(mono_new, file.path(git_geneset,'Monocyte_4022.txt') )	
	
#' 	
#' 	
#' ### 4022 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Monocyte/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
	
other_mono = Permutation_list$`RF_permutation_First_Monocyte vs other`	
	
mono_computation = importance_barplot(other_mono, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 4022 genes selected in Other vs Monocyte classification")	
                                        	
mono_computation$plot + 	
  theme(axis.text.y = element_text(size=0))	
	
	
mono_top1870 = mono_computation$ordered_permutation_dataframe[mono_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
writeLines(mono_top1870, file.path(git_geneset,'Monocyte_1870.txt') )	
writeLines(mono_top1870, file.path(data_dir,'Monocyte_1870.txt') )	
	
#' 	
#' 	
#' 	
#' ### 1870 features	
#' 	
#' 	
	
other_mono = Permutation_list$`RF_permutation_second_Monocyte vs other`	
	
	
	
mono_computation = importance_barplot(other_mono, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 4022 genes selected in Other vs Monocyte classification")	
                                        	
mono_computation$plot + 	
  theme(axis.text.y = element_text(size=0))	
	
	
mono_top369 = mono_computation$ordered_permutation_dataframe[mono_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
writeLines(mono_top369, file.path(git_geneset,'Monocyte_369.txt') )	
writeLines(mono_top369, file.path(data_dir,'Monocyte_369.txt') )	
	
#' 	
#' 	
#' 	
#' ### 369 features	
#' 	
#' 	
	
other_mono = Permutation_list$`RF_permutation_third_Monocyte vs other`	
	
	
mono_computation = importance_barplot(other_mono, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 4022 genes selected in Other vs Monocyte classification")	
                                        	
mono_computation$plot + 	
  theme(axis.text.y = element_text(size=5))	
	
	
mono_top174 = mono_computation$ordered_permutation_dataframe[mono_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
writeLines(mono_top174, file.path(git_geneset,'Monocyte_174.txt') )	
writeLines(mono_top174, file.path(data_dir,'Monocyte_174.txt') )	
	
#' 	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/Monocyte/" ),	
  til10_path = file.path(path_heatmap, "TIL10/Monocyte/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/Monocyte/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","Monocyte"))	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/Monocyte/" ),	
  til10_path = file.path(path_heatmap, "TIL10/Monocyte/" ),	
  lm22_path = file.path(path_heatmap, "LM22/Monocyte/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","Monocyte"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
cf = rbind(cf_tbs,cf_Liu2021)	
	
comp = build_confusion_heatmap(cf,  	
                               c("Other_cells","Monocyte"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Monocyte Classification") 	
#' 	
#' 	
#' ## Macrophage vs other	
#' 	
#' **68 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_macro = groupwise_mean[str_detect(colnames(groupwise_mean),"Macrophage")]	
FC_df_macro = Get_value_RA(FC_df_macro, expr_median)	
FC_df_macro = FC_prop_expr(FC_df_macro)	
	
	
elbow=kneedle_threshold(FC_df_macro, FC_df_macro$FC_Macrophage, sens=1, min=1) 	
elbow$plot	
macro = Rprop_plot(FC_df_macro,FC_x = "FC_Macrophage", y="Macrophage_prop_expr", 	
                title = "Log Fold Change of Macrophage vs other depending on the proportion of Macrophage cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
macro$plot	
	
	
macro_new = row.names(macro$df[macro$df$selection != "Not selected",])	
	
	
writeLines(macro_new, file.path(data_dir,'Macrophage_369.txt') )	
writeLines(macro_new, file.path(git_geneset,'Macrophage_369.txt') )	
	
#' 	
#' 	
#' 	
#' 	
#' ### 369 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Macrophage/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_macro = Permutation_list$`RF_permutation_First_Macrophage vs other`	
other_macro$genes = row.names(other_macro)	
	
macro_computation = importance_barplot(other_macro, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 369 genes selected in Other vs Macrophage classification")	
                                        	
macro_computation$plot + 	
  theme(axis.text.y = element_text(size=7))	
	
	
macro_top68 = macro_computation$ordered_permutation_dataframe[macro_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
	
writeLines(macro_top68, file.path(data_dir,'Macrophage_68.txt') )	
writeLines(macro_top68, file.path(git_geneset,'Macrophage_68.txt') )	
	
#' 	
#' 	
#' 	
#' ### Performances - No Macrophages in Liu2021 data	
#' 	
#' 	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/Macrophage/" ),	
  til10_path = file.path(path_heatmap, "TIL10/Macrophage/" ),	
  lm22_path = file.path(path_heatmap, "LM22/Macrophage/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = c("Other_cells","Macrophage"))	
	
cf_tbs$dataset = "Tabula Sapiens data"	
	
comp = build_confusion_heatmap(cf_tbs,  	
                               c("Other_cells","Macrophage"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
	
        	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Macrophage Classification") 	

#' 	
#' ## Platelet vs other	
#' 	
#' **40 Features selected**	
#' 	
#' 	
	
row.names(groupwise_mean)= groupwise_mean$genes	
FC_df_plat = groupwise_mean[str_detect(colnames(groupwise_mean),"Platelet")]	
FC_df_plat = Get_value_RA(FC_df_plat, expr_median)	
FC_df_plat = FC_prop_expr(FC_df_plat)	
	

	
	
elbow=kneedle_threshold(FC_df_plat, FC_df_plat$FC_Platelet,sens = 1) 	
elbow$plot	
plat = Rprop_plot(FC_df_plat,FC_x = "FC_Platelet", y="Platelet_prop_expr", 	
                title = "Log Fold Change of Platelet vs other depending on the proportion of Platelet cells expressing the genes",	
                add_label = TRUE, FC_threshold = 1, 	
                prop_threshold = 0.25, 	
                FC_thresh2 = elbow$threshold) 	
	
plat$plot	
	
	
platelet_new = row.names(plat$df[plat$df$selection != "Not selected",])	
	
	
writeLines(platelet_new, file.path(data_dir,'Platelet_510.txt') )	
writeLines(platelet_new, file.path(git_geneset,'Platelet_510.txt') )	
	
#' 	
#' 	
#' ### 510 features	
#' 	
#' 	
	
Permutation_list = get_permutation_list(	
  directory = file.path(result_dir, "Platelet/" ),	
  pattern2= "/RF_permutation/Permutation_result_")	
	
other_plat = Permutation_list$`RF_permutation_First_Platelet vs other`	
other_plat$genes = row.names(other_plat)	
	
plat_computation = importance_barplot(other_plat, 50, 	
                                        filter = "mean", 	
                                        title = "Mean feature importance of 50 permutations with 510 genes selected in Other vs Macrophage classification")	
                                        	
plat_computation$plot + 	
  theme(axis.text.y = element_text(size=10))	
	
	
plat_top40 = plat_computation$ordered_permutation_dataframe[plat_computation$ordered_permutation_dataframe$mean_importance > 0,]$genes	
	
writeLines(plat_top40, file.path(git_geneset,'Platelet_40.txt') )	
writeLines(plat_top40, file.path(data_dir,'Platelet_40.txt') )	
	
#' 	
#' 	
#' ### Performances - Platelet unique to Tabula Sapiens	
#' 	
#' 	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = confusion_matrix_setup(	
    directory = file.path(path_heatmap, "Tabula_Sapiens_50k-10celltypes/Platelet/" ),	
    pattern2= "/confusion_heatmap_",	
    cells = c("Other_cells","Platelet"),	
    signature_matrix = "Tabula Sapiens")	
	
cf_tbs$dataset = "Tabula Sapiens data"	
	
comp = build_confusion_heatmap(cf_tbs,  	
                               c("Other_cells","Platelet"),	
                               filter = FALSE,	
                               direction_color_palette = 1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Platelet Classification") 	
#' 	
#' 	
#' # Combining Features	
#' 	
#' 	
	
	
LM22 = read.table(file.path(path,"results/Tabula_sapiens/LM22/LM22_merged.txt"), sep="\t")	
TIL10 = read.table(file.path(path,"results/Tabula_sapiens/LM22/TIL10_binary_standardized.txt"), sep="\t")	
	
	
LM22_list = lapply( names(LM22), function(x) {	
  geneset = row.names(LM22[LM22[[x]] ==1 ,])	
  return(geneset) })	
names(LM22_list) = colnames(LM22)	
	
TIL10_list = lapply( names(TIL10), function(x) {	
  geneset = row.names(TIL10[TIL10[[x]] ==1 ,])	
  return(geneset) })	
names(TIL10_list) = colnames(TIL10)	
	
	
TabulaSapiens_geneset = list(	
  "NK_cell"   = unique(c(top79_nkcd8_genes,top178_nk)),	
  "Neutrophil"= neutro_top84,	
  "Monocyte"  = mono_top174,	
  "Plasmocyte"= plasmo_top19,	
  	
  "CD4"       = unique(c(CD4_top114, top61_cd48_genes)),	
  "CD8"       = unique(c(top61_cd48_genes, top79_cd8, top247_nkcd8_genes)), 	
  	
  "Macrophage"= macro_top68,	
  "Platelet"  = plat_top40,	
  "memory_B_cell"   = unique(c(top42_naive_mem, top53_memB)), 	
  "naive_B_cell"    = unique(c(naiveB_top29,top42_naive_mem))	
  )	
	
#save(TabulaSapiens_geneset, file = file.path(git_geneset,"TabulaSapiens_geneset_final.Rdata"))	
  	
geneset_flat = unique(unlist(TabulaSapiens_geneset))	
	
	
upset(fromList(TabulaSapiens_geneset),     	
      sets.bar.color = "#56B4E9", order.by = "freq", nsets = 50,	
      empty.intersections = NULL, set_size.show = TRUE, mb.ratio = c(0.7, 0.3),	
      text.scale = 1.75, set_size.scale_max = 400)	
	
plot_list = list()	
for (cell in names(TabulaSapiens_geneset)) {	
  if (cell != "Platelet" ){	
  	
  cell_list = list(	
    TabulaSapiens_geneset[[cell]],	
    TIL10_list[[cell]], 	
    LM22_list[[cell]] )	
  	
  names(cell_list) = c(paste0("Tabula Sapiens - ",cell),	
                       paste0("TIL10 - ",cell), 	
                       paste0("LM22 - ",cell) )	
  	
  plot_list[[cell]]= upset(fromList(cell_list),     	
      sets.bar.color = "#56B4E9", order.by = "freq", nsets = 50,	
      empty.intersections = NULL, set_size.show = TRUE, mb.ratio = c(0.7, 0.3),	
      text.scale = 1.75, set_size.scale_max = max(unlist(lapply(cell_list, function(x) length(x)))) * 1.25 )	
}}	
	
plot_list	
	
	
#' 	
#' 	
#' ### Performances	
#' 	
#' 	
 path_heatmap = file.path(path,"results/Liu2021/Signature_matrix_v2/")	
cells = c("CD4","CD8","Monocyte","NK_cell","memory_B_cell","naive_B_cell")	
	
cf_Liu2021 = merge_confusion_matrix( 	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/Multiclass/" ),	
  til10_path = file.path(path_heatmap, "TIL10/Multiclass/RF/" ),	
  lm22_path = file.path(path_heatmap, "LM22/Multiclass/RF/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = cells)	
	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
cf_Liu2021$dataset = "Liu 2021 Cell data"	
	
comp = build_confusion_heatmap(cf_Liu2021,	
                               c("NK_cell","CD8",	
                                 "CD4","Monocyte",	
                                 "memory_B_cell","naive_B_cell"),	
                               filter = FALSE,	
                               direction_color_palette = -1)	
	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold"))+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Multiclass Classification") 	
#' 	
#' 	
#' 	
#' 	
#' 	
cells = c("CD4","CD8","Erythrocyte","Macrophage",	
          "Monocyte","NK_cell","Neutrophil", 	
          "memory_B_cell","naive_B_cell")	
path_heatmap = file.path(path,"results/Tabula_sapiens/Signature_matrix_v2/")	
	
cf_tbs = merge_confusion_matrix(	
  tbs_path = file.path(path_heatmap, "Tabula_Sapiens_50k/Multiclass/" ),	
  til10_path = file.path(path_heatmap, "TIL10/Multiclass/" ),	
  lm22_path = file.path(path_heatmap, "LM22/Multiclass/" ),	
  pattern2= "/confusion_heatmap_",	
  cells = cells)	
	
cf_tbs$dataset = "Tabula Sapiens data"	
	
cf_tbs = cf_tbs[ cf_tbs$True_label != "Erythrocyte" &  cf_tbs$Predicted_label != "Erythrocyte",]	
	
balanced_acc1 = round(mean(cf_tbs[cf_tbs$True_label == cf_tbs$Predicted_label & cf_tbs$signature_matrix == "Tabula Sapiens",]$Percentage)*100,2)	
	
balanced_accLM22 = round(mean(cf_tbs[cf_tbs$True_label == cf_tbs$Predicted_label & cf_tbs$signature_matrix == "LM22",]$Percentage)*100,2)	
	
balanced_accTIL10 = round(mean(cf_tbs[cf_tbs$True_label == cf_tbs$Predicted_label & cf_tbs$signature_matrix == "TIL10",]$Percentage)*100,2)	
	
	
cf_tbs[cf_tbs$signature_matrix == "Tabula Sapiens",]$subtitle = 	
  paste0("Tabula Sapiens \n Balanced accuracy = ",balanced_acc1,"%")	
	
cf_tbs[cf_tbs$signature_matrix == "TIL10",]$subtitle = 	
  paste0("TIL10 \n Balanced accuracy = ",balanced_accTIL10,"%")	
	
cf_tbs[cf_tbs$signature_matrix == "LM22",]$subtitle = 	
  paste0("LM22 \n Balanced accuracy = ",balanced_accLM22,"%")	
  	
  	
comp = build_confusion_heatmap(cf_tbs,	
                               c("NK_cell","CD8","CD4",	
                                 "Monocyte","Macrophage","Neutrophil",	
                                 "memory_B_cell","naive_B_cell"),	
                               filter = FALSE,	
                               direction_color_palette = -1)	
comp$plot + theme_dark() + 	
  facet_wrap(. ~ dataset+subtitle, ncol=3) + 	
  theme(strip.text = element_text(size=14,face="bold"),	
        panel.background  = element_rect(fill = "white"),	
        panel.grid.major  = element_line(color = "white"),	
        panel.border = element_rect(color = "black",fill=NA),	
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),	
        plot.title = element_text(hjust = 0.5, size = 16,face="bold") )+	
        	
  ggtitle("Signature matrix prediction performance with Random Forest - Multiclass Classification") 	
#' 	
#' 	
#' 	
TabulaSapiens_binary = geneset_to_binary_df(TabulaSapiens_geneset)	
print(table(rowSums(TabulaSapiens_binary)))	
	
writeLines(geneset_flat, file.path(path,"/results/Tabula_sapiens/LM22/signature_genes_TabulaSapiens_50k.txt"))	
	
write.table(TabulaSapiens_binary, file.path(path,"/results/Tabula_sapiens/LM22/TabulaSapiens_binary.txt"),	
             quote = FALSE, row.names = TRUE, sep='\t')	
            	
            	
	
#' 	
#' 	
#' # Mean Expression heatmap	
#' 	
#' 	
	
heatmap = groupwise_mean[groupwise_mean$genes %in% geneset_flat,]	
heatmap = heatmap[,! str_detect(colnames(heatmap), "other|Erythrocyte|prop_expr")]	
	
	
row.names(heatmap) = heatmap$genes	
heatmap$genes = NULL	
	
write.table(heatmap, file.path(data_dir,"Signature_TabulaSapiens_977.txt"),	
          quote = FALSE, row.names = TRUE, sep='\t')	
write.table(heatmap, file.path(git_dir,"/Signature_matrix/Signature_TabulaSapiens_977_10celltypes.txt"),	
          quote = FALSE, row.names = TRUE, sep='\t')	
	
#' 	
#' 	
#' 	
	
	
genes_annotation = heatmap	
genes_annotation$NK_cell = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$NK_cell, "yes" , "no")	
genes_annotation$CD4 = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$CD4, "yes" , "no")	
genes_annotation$CD8 = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$CD8,"yes" , "no")	
genes_annotation$Neutrophil = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$Neutrophil,"yes" , "no")	
genes_annotation$naive_B_cell = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$naive_B_cell,"yes" , "no")	
genes_annotation$memory_B_cell = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$memory_B_cell,"yes" , "no")	
genes_annotation$Plasmocyte = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$Plasmocyte,"yes" , "no")	
genes_annotation$Monocyte = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$Monocyte,"yes" , "no")	
genes_annotation$Macrophage = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$Macrophage,"yes" , "no")	
genes_annotation$Platelet = ifelse(row.names(genes_annotation) %in% TabulaSapiens_geneset$Platelet,"yes" , "no")	
	
	
color_yes = "darkgreen"	
color_no = "#CCCCCC"	
	
my_colour = list(	
  B_cell = c(yes = color_yes, no = color_no),	
  naive_B_cell = c(yes = color_yes, no = color_no),	
  memory_B_cell = c(yes = color_yes, no = color_no),	
  Neutrophil = c(yes = color_yes, no = color_no),	
  CD4 = c(yes = color_yes, no = color_no), 	
  NK_cell = c(yes = color_yes, no = color_no) ,	
  CD8 = c(yes = color_yes, no = color_no),	
  Plasmocyte = c(yes = color_yes, no = color_no),	
  Monocyte = c(yes = color_yes, no = color_no),	
  Macrophage = c(yes = color_yes, no = color_no),	
  Platelet = c(yes = color_yes, no = color_no) )	
	
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=0.95, name="vp", just=c("right","top"))), action="prepend")	
pheatmap(t(heatmap), scale = "column", 	
         cutree_cols = 10,	
         angle_col = "45", main = "Mean expression heatmap of the 977 features used to predict 10 cell types",	
         clustering_distance_rows = 'correlation', 	
         clustering_distance_cols = 'correlation',	
         clustering_method = "ward.D2",	
         annotation_col = genes_annotation, annotation_colors = my_colour,	
         color = colorRampPalette(c("navy", "#DDDDDD", "firebrick"))(25),	
         show_colnames = FALSE)	
setHook("grid.newpage", NULL, "replace")	
grid.text("Genes", y=-0.02, gp=gpar(fontsize=16))	
	

	
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=1, name="vp", just=c("right","top"))), action="prepend")	
pheatmap(heatmap, scale = "row", 	
         cutree_rows = 10,	
         angle_col = "45", main = "Mean expression heatmap of the 977 features used to predict 10 cell types",	
         clustering_distance_rows = 'correlation', 	
         clustering_distance_cols = 'correlation',	
         clustering_method = "ward.D2",	
         annotation_row = genes_annotation, annotation_colors = my_colour,	
         color = colorRampPalette(c("navy", "#DDDDDD", "firebrick"))(25),	
         show_rownames = FALSE)	
setHook("grid.newpage", NULL, "replace")	
grid.text("Genes", y=0.35, x=0.9, gp=gpar(fontsize=16))	
