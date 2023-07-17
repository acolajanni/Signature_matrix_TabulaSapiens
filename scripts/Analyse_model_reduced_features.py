################################################################################                                                                           
# > february 2023                                                                                                             
# > Script : Analyze_model_reduced_features.py                                                                                       
# > Function : Analyze model performance on a reduced number of features
# @ COLAJANNI Antonin                                                          
################################################################################

from functions.function_SingleCell import *
from functions.dataviz import *
from functions.utils import *
from collections import Counter



def run_on_reduced_dataset_from_dictionnary(dictionnary, X, path_result, 
                                            filename,OnevsAll=True, to_remove = [], 
                                            ROC=True,
                                            dataset = "Liu2021",
                                            color='rocket'):

    if not OnevsAll: 
        # cell class available in X 
        cells_avail = list(dict.fromkeys(X.index.to_list().copy()))
    
        # Celltypes in dictionnary that are also present in X
        signature_cell=[ p for p in list(dictionnary.keys()) for cell in cells_avail if p in cell]
    
        # Create a dictionnary['Multiclass'] = all unique genes
        res = []
        for celltype in signature_cell :             
            res.append(dictionnary[celltype] )
        dictionnary = {}
        dictionnary["Multiclass"] = list(dict.fromkeys([item for sublist in res for item in sublist]))
        
    for geneset in dictionnary:
        y = X.index.to_list().copy()
        print(geneset, "\n",  dictionnary[geneset] )
        
        if dictionnary[geneset] == [] : 
            continue
        
        if dataset == "Liu2021" and (geneset == "Neutro" or geneset == "Plasm") : # No neutrophil/plasma cell in Liu2021 data
            continue


        if OnevsAll : #2 class : 0 vs 1
            geneset_in_data = [match for match in dictionnary[geneset] if match in X.columns.to_list()]
            y_binary = replace_label(y, [geneset], [1], 0) # type: ignore
            _, le, _ = numeric_encoding(y_binary)
            le_name_mapping = {0:0 , 1:1}
            X_tmp = X[geneset_in_data]
            path_geneset = create_dir(path_result+geneset+"/")
            
        else : # Multiclasses                         
            X_tmp, y, geneset_in_data = get_existing_index_features(X, dictionnary[geneset], y, to_remove= to_remove)  
            y_binary, le, le_name_mapping = numeric_encoding(y)    
            print(len(X), "shape raw")
            print(len(X_tmp), "shape updated")
            path_geneset = create_dir(path_result+"/Multiclass/")

            
            
        gene_number = len(geneset_in_data)
        
            
        #run_on_reduced_dataset(X_tmp[geneset_in_data], y_binary, geneset_in_data, 
        #               le, le_name_mapping, LR, create_dir(path_geneset+"/LR/"), 
        #               filename+"_"+geneset+"_LR_"+str(gene_number)+"_features",
        #               ROC = ROC, color=color)
        
        run_on_reduced_dataset(X_tmp[geneset_in_data], y_binary, geneset_in_data, 
                       le, le_name_mapping, RF, create_dir(path_geneset),#+"/RF/"), 
                       filename+"_"+geneset+"_RF_"+str(gene_number)+"_features", 
                       ROC=ROC, color=color)
                
        print("-------------------------------------------------------")       
        
    return None

################################################################################

LR = LogisticRegressionCV(n_jobs=-1, verbose=3,
                             Cs=10, cv=5, max_iter = 2500,
                             random_state=0, penalty='l1',
                             scoring='neg_mean_squared_error', 
                             solver="saga")

RF = RandomForestClassifier(n_jobs=-1,
                            n_estimators=1000,
                            random_state=0 )

################################################################################


# Signature matrix assessments
path = "/shared/projects/microbiome_translocation/results/Tabula_sapiens/LM22/"

LM22      = pd.read_csv(path+"LM22_merged.txt", sep="\t")
Vallania  = pd.read_csv(path+"Vallania_Binary_reduced.txt", sep="\t")
Sigmatrix = pd.read_csv(path+"Sigmatrix_Binary_reduced.txt", sep="\t")
TIL10     = pd.read_csv(path+"TIL10_binary_standardized.txt", sep="\t")
TBS       = pd.read_csv(path+"TabulaSapiens_binary.txt", sep="\t")

#TBS = pd.read_csv(path+".txt", sep="\t")

# Retrieve genesets from binary matrix and make dictionnary : {CD4 : [IL7R, NKG7, ......], CD8 : [...]}
TIL10_genesets     = get_geneset_from_binary(TIL10)
Sigmatrix_genesets = get_geneset_from_binary(Sigmatrix)
Vallania_genesets  = get_geneset_from_binary(Vallania)
LM22_genesets      = get_geneset_from_binary(LM22)
TBS_genesets       = get_geneset_from_binary(TBS)

## GENESETS
path_data     = "/shared/projects/microbiome_translocation/data/Tabula_sapiens_immune_all/"
path_data_cv  = "/shared/projects/microbiome_translocation/data/scRNAseq/Liu_2021_cell/"

# Comment the dataset to not remove from the two
#dataset = 'Liu2021'
dataset = 'Tabula_Sapiens'

 
# Genesets to evaluate performances
Reference_gs = [
    #Sigmatrix_genesets, Vallania_genesets, 
    LM22_genesets, TIL10_genesets, 
    TBS_genesets   
    ]

Reference_names = [
    # "Sigmatrix", "Vallania", "LM22", 
    "LM22", "TIL10",
    #"Tabula_Sapiens_50k-10celltypes"
    "Tabula_Sapiens_50k"
    ]

for i in range(len(Reference_gs)): 
    
    if dataset == "Liu2021": 
        path_res=create_dir("/shared/projects/microbiome_translocation/results/Liu2021/Signature_matrix_v2/")
        mat=pd.read_csv(path_data_cv+'Liu2021_Cell.csv', sep='\t', index_col=0)#, nrows=2500)
        y = pd.read_csv(path_data_cv+"Liu2021_index.txt",header=None)[0].to_list()#[0:2500]
        mat.index = y
        
        # These cell types are absent from this dataset
        Reference_gs[i]["Neutrophil"] = []
        Reference_gs[i]["Plasmocyte"] = [] # not in TIL10
        Reference_gs[i]["Macrophage"] = [] 
        Reference_gs[i]["Platelet"]   = [] # not in TIL10 and LM22
        
    elif dataset == "Tabula_Sapiens" :
        path_res=create_dir("/shared/projects/microbiome_translocation/results/Tabula_sapiens/Signature_matrix_v2/")
        mat=pd.read_csv(path_data+'50000cell_immune_expressed.csv', sep='\t', index_col=0)
        y = pd.read_csv(path_data+"Celltypes_of_interest_50k.txt",header=None)[0].to_list()
        mat.index = y
        
        #Reference_gs[i]["Plasmocyte"] = [] # not in TIL10
        #Reference_gs[i]["Platelet"]   = [] # not in TIL10 and LM22
    
    
    path_result = create_dir(path_res+Reference_names[i]+"/")
    print(Reference_names[i])
    
    to_remove = []
    

    # Uncomment to Compare NK vs CD4 vs CD8
    #to_remove = ["Monocyte", "Plasmocyte", "Neutrophil", 
    #             "naive_B_cell", "memory_B_cell" ] 
        
    # Uncomment to Compare common celltypes in Tabula sapiens dataset 
    to_remove = ["Plasmocyte", "Platelet" ]  
            
    run_on_reduced_dataset_from_dictionnary(Reference_gs[i],
                                            mat, 
                                            path_result, 
                                            Reference_names[i]+"_Multiclass",
                                            OnevsAll=False, 
                                            ROC= False,
                                            to_remove = to_remove ,
                                            color='rocket_r')

    run_on_reduced_dataset_from_dictionnary(Reference_gs[i], 
                                            mat,
                                            path_result, 
                                            Reference_names[i],
                                            OnevsAll=True, 
                                            ROC= True,
                                            to_remove = to_remove,
                                            color='rocket')