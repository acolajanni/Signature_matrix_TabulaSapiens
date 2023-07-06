################################################################################                                                                           
# > january-february 2023                                                                                                             
# > Script : extract_features_FromReducedDataset.py                                                                                              
# > Function : extracting features / important genes used to recognized cells with RF permut followed by LassoCV 
# @ COLAJANNI Antonin                                                          
################################################################################

from functions.utils import *
from functions.function_SingleCell import *
from functions.dataviz import *
from collections import Counter
import argparse
import sys

    
#################
# Path distant : SFBI
path_data="/shared/projects/microbiome_translocation/data/Tabula_sapiens_immune_all/"
path_res=create_dir("/shared/projects/microbiome_translocation/results/Tabula_sapiens/")

# Path local 
#path_data="/home/acolajanni/Documents/work/data/Tabula_sapiens_immune_all/"
#path_res="/home/acolajanni/Documents/work/results/Tabula_sapiens/"


# path_res=create_dir(path_res+"LassoCV_python/")
path_res=create_dir(path_res+"RF_filter_50k/")


path_res_CD8        = create_dir(path_res+"CD8/")
path_res_CD4        = create_dir(path_res+"CD4/")
path_res_CD48       = create_dir(path_res+"CD4-CD8/")
path_res_NKCD8       = create_dir(path_res+"NK-CD8/")
path_res_Naivebcell = create_dir(path_res+"Naive_Bcell/")
path_res_Membcell   = create_dir(path_res+"Memory_Bcell/")
path_res_NaiveMemory= create_dir(path_res+"Memory-Naive_Bcell/")

path_res_Plasmocyte = create_dir(path_res+"Plasmocyte/")
path_res_Monocyte   = create_dir(path_res+"Monocyte/")
path_res_NK         = create_dir(path_res+"NK/")
path_res_Neutrophil = create_dir(path_res+"Neutrophil/")
path_res_macrophage = create_dir(path_res+"Macrophage/")
path_res_platelet   = create_dir(path_res+"Platelet/")


# Get gene subset for other cell types : Bcell, plasmocyte, monocytes, neutrophils, NK, etc

geneset_plasmocyte = pd.read_csv(path_data+"50k_cells_geneset/Plasmocyte_5172.txt",header=None)[0].to_list()
geneset_Monocyte = pd.read_csv(path_data+"50k_cells_geneset/Monocyte_4022.txt",header=None)[0].to_list()
geneset_NK = pd.read_csv(path_data+"50k_cells_geneset/NK_926.txt",header=None)[0].to_list()
geneset_CD8 = pd.read_csv(path_data+"50k_cells_geneset/CD8_639.txt",header=None)[0].to_list()
geneset_CD4 = pd.read_csv(path_data+"50k_cells_geneset/CD4_1402.txt",header=None)[0].to_list()
geneset_neutrophil = pd.read_csv(path_data+"50k_cells_geneset/Neutro_825.txt",header=None)[0].to_list()

geneset_CD48 = pd.read_csv(path_data+"50k_cells_geneset/CD4vsCD8_3435.txt",header=None)[0].to_list()
geneset_NKCD8 = pd.read_csv(path_data+"50k_cells_geneset/CD8vsNK_1028.txt",header=None)[0].to_list()
geneset_NaiveMemory=pd.read_csv(path_data+"50k_cells_geneset/naiveVSmemory_B_3582.txt",header=None)[0].to_list()

geneset_MemBcell = pd.read_csv(path_data+"50k_cells_geneset/memoryB_7208.txt",header=None)[0].to_list()
geneset_NaiveBcell= pd.read_csv(path_data+"50k_cells_geneset/naiveB_6098.txt",header=None)[0].to_list()
geneset_macrophage = pd.read_csv(path_data+"50k_cells_geneset/Macrophage_369.txt",header=None)[0].to_list()
geneset_platelet = pd.read_csv(path_data+"50k_cells_geneset/Platelet_510.txt",header=None)[0].to_list()


geneset_CD4 = pd.read_csv(path_data+"50k_cells_geneset/CD4_319.txt",header=None)[0].to_list()
geneset_NKCD8 = pd.read_csv(path_data+"50k_cells_geneset/NKvsCD8_247.txt",header=None)[0].to_list()


# Import
mat=pd.read_csv(path_data+'Tabula_Sapiens.csv', sep='\t', index_col=0)
print("import done")

y = pd.read_csv(path_data+"Celltypes_of_interest_50k.txt",header=None)[0].to_list()
mat.index = y


################################################################################
RF = RandomForestClassifier(n_jobs=-1,
                            n_estimators=1000,
                            random_state=0,
                            verbose=1)

################################################################################


def RF_permutation_multiple_cell(cell_list, filename, n_perm=50, save_rate=5, model = RF, save_dir_extension="RF_permutation/", cell_type_path=path_data+"Celltypes_of_interest_50k.txt"):    
    for cell in cell_list:
        y = pd.read_csv(cell_type_path,header=None)[0].to_list()
         
        if cell == "naive_Bcell" :
            X = mat[geneset_NaiveBcell].reset_index(drop=True)
            geneset = geneset_NaiveBcell
            result_directory = path_res_Naivebcell
            y_binary = replace_label(y, ['naive_B_cell'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "memory_Bcell" :
            X = mat[geneset_MemBcell].reset_index(drop=True)
            geneset = geneset_MemBcell
            result_directory = path_res_Membcell
            y_binary = replace_label(y, ['memory_B_cell'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)

        elif cell == "NK":
            X = mat[geneset_NK].reset_index(drop=True)
            geneset = geneset_NK
            result_directory = path_res_NK
            y_binary = replace_label(y, ['NK'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "Neutrophil":
            X = mat[geneset_neutrophil].reset_index(drop=True)
            geneset = geneset_neutrophil
            result_directory = path_res_Neutrophil
            y_binary = replace_label(y, ['Neutrophil'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "Monocyte":
            X = mat[geneset_Monocyte].reset_index(drop=True)
            geneset = geneset_Monocyte
            result_directory = path_res_Monocyte
            y_binary = replace_label(y, ['Monocyte'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "Plasmocyte":
            X = mat[geneset_plasmocyte].reset_index(drop=True)
            geneset = geneset_plasmocyte
            result_directory = path_res_Plasmocyte
            y_binary = replace_label(y, ['Plasmocyte'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
        
        elif cell == "CD4":
            X = mat[geneset_CD4].reset_index(drop=True)
            geneset = geneset_CD4
            result_directory = path_res_CD4
            y_binary = replace_label(y, ['CD4'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "CD8":
            X = mat[geneset_CD8].reset_index(drop=True)
            geneset = geneset_CD8
            result_directory = path_res_CD8
            y_binary = replace_label(y, ['CD8'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
        
        elif cell == "Macrophage":
            X = mat[geneset_macrophage].reset_index(drop=True)
            geneset = geneset_macrophage
            result_directory = path_res_macrophage
            y_binary = replace_label(y, ['Macrophage'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "Platelet":
            X = mat[geneset_platelet].reset_index(drop=True)
            geneset = geneset_platelet
            result_directory = path_res_platelet
            y_binary = replace_label(y, ['Platelet'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "CD48":
            X = mat[geneset_CD48].reset_index(drop=True)
            X.index = y 

            geneset = geneset_CD48
            result_directory = path_res_CD48
            y_binary = y
            
            y_binary = [match for match in y if match in ['CD4','CD8']]
            CD8_index = X.index == "CD8"
            CD4_index = X.index == "CD4"
            CD48_index = np.logical_or(CD8_index, CD4_index)
            X = X[CD48_index]
            y_binary = replace_label(y_binary, ['CD8'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "NKCD8":
            X = mat[geneset_NKCD8].reset_index(drop=True)
            X.index = y 

            #X = X_CD48
            geneset = geneset_NKCD8
            result_directory = path_res_NKCD8
            y_binary = y
            
            y_binary = [match for match in y if match in ['NK_cell','CD8']]
            CD8_index = X.index == "CD8"
            NK_index = X.index == "NK_cell"
            NKCD8_index = np.logical_or(CD8_index, NK_index)
            X = X[NKCD8_index]
            y_binary = replace_label(y_binary, ['CD8'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        elif cell == "NaiveMemory":
            X = mat[geneset_NaiveMemory].reset_index(drop=True)
            X.index = y 

            geneset = geneset_NaiveMemory
            result_directory = path_res_NaiveMemory
            y_binary = y
            
            y_binary = [match for match in y if match in ['memory_B_cell','naive_B_cell']]
            naive_index = X.index == "memory_B_cell"
            memory_index = X.index == "naive_B_cell"
            naiveMemory_index = np.logical_or(naive_index, memory_index)
            X = X[naiveMemory_index]
            y_binary = replace_label(y_binary, ['memory_B_cell'], [1], 0) 
            _, le, _ = numeric_encoding(y_binary)
            
        else: 
            print("Error, cell type not recognized:",cell )
            return None 
        
        condition = f'{cell} vs other'
        le_name_mapping = {"other":0 , cell:1}
        
        if cell == "CD48":
           condition = f'CD4 vs CD8'
           le_name_mapping = {"CD4":0 , 'CD8':1} 
           
        elif cell == "NKCD8":
           condition = f'NK vs CD8'
           le_name_mapping = {"NK_cell":0 , 'CD8':1} 
        
        elif cell == "NaiveMemory":
           condition = f'Naive vs Memory'
           le_name_mapping = {"memory_B_cell":0 , 'naive_B_cell':1} 
        
        
        print( cell,"- FC Filter - Permutations")
        save_dir=create_dir(result_directory+save_dir_extension)
        X_train,X_test, y_train,y_test = train_test_split(X,y_binary, test_size=0.3, random_state=1)

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        if cell == "CD48":
            acc, clf, cm = analyse_model( y_test, y_pred )
            cell = "CD4"
            ROC_values = None
                        
        if cell == "NKCD8":
            acc, clf, cm = analyse_model( y_test, y_pred )
            cell = "NK_cell"
            ROC_values = None
            
        if cell == "NaiveMemory":
            acc, clf, cm = analyse_model( y_test, y_pred )
            cell = "memory_B_cell"
            ROC_values = None
                       
        else :
            acc, clf, cm = analyse_model(le.inverse_transform(y_test), le.inverse_transform(y_pred) )
            ROC_values = ROC_analysis(model, X_train, y_train, le_name_mapping, positive_label=cell)
        
        metrics = [acc,clf,cm]
        importance_df = get_importance(model, geneset )  
        save_metrics(
            acc=metrics[0], 
            clf=metrics[1],
            cm=metrics[2],
            imp_df = importance_df ,
            save_dir=save_dir , 
            filename=filename, 
            ROC_values=ROC_values)
            
        Full_model_evaluation([acc,clf,cm], importance_df, 
            f"Confusion heatmap {condition}", 
            f"Classification report heatmap {condition}", 
            f"Gene importance in classification {condition}",
            save_dir, f'{filename}_{condition}', 
            n_genes=len(importance_df[importance_df.importance.cumsum()<=0.95]) )


        n_perm_before_save = int( n_perm/save_rate )
        result = []
        c = 0
        for i in range(n_perm_before_save):
            importance_dict = permutation_importance(model, 
                                       X_test, y_test, 
                                       n_repeats=save_rate, random_state= n_perm+i , 
                                       n_jobs=-1, scoring="neg_mean_squared_error")
            
            result = result_from_permutations( 
                    importance_dict,
                    geneset,
                    sort=False,
                    filename=filename+"_"+condition+"_batch"+str(i+1),
                    save_dir=save_dir)
            
        c=0
        for df in result:
            df.drop(["mean_importance","std_importance"],axis=1,inplace=True)
            df.columns = [i for i in range(c,save_rate+c)]
            c+=save_rate
            
        result = pd.concat(result, axis=1, ignore_index=False)    
        result['mean_importance'] = result.mean(axis=1)
        result['std_importance'] = result.std(axis=1)
        
        result.to_csv(f'{save_dir}/Permutation_result_{filename}_{condition}.csv')


    return result


full_list = [
    "NaiveMemory",
    "CD48",
    "NKCD8", 
    "NK", 
    "Neutrophil",
    "naive_Bcell",
    "memory_Bcell",
    "Monocyte" 
    "Plasmocyte" ,
    "CD4" , 
    "CD8",
    "Macrophage",
    "Platelet"
    ]


RF_permutation_multiple_cell(full_list, "RF_permutation_First", 50, 5, RF)
RF_permutation_multiple_cell(["CD4","NKCD8","Monocyte"] , "RF_permutation_second", 50, 5, RF)
