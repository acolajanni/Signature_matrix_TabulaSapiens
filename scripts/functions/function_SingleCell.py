################################################################################                                                                           
# > December 2022                                                                                                             
# > Script : function.py                                                                                                           
# > Function : Used to get a reservoir of function                                 
# @ COLAJANNI Antonin                                                          
################################################################################

import pandas as pd        
import numpy as np

from sklearn.metrics import confusion_matrix, classification_report, accuracy_score, roc_auc_score, roc_curve, balanced_accuracy_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LassoCV, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.multiclass import OneVsOneClassifier
from sklearn.pipeline import Pipeline
from sklearn.inspection import permutation_importance
from .dataviz import ROC_curve_3classes, Coordinate_descent_Multiplot_logit, Full_model_evaluation
from .utils import *

import joblib
import os, re, os.path
from collections import Counter




"""
Extract labels and features from dataframe

Parameters
----------
    df : pandas.Dataframe
        df of expression with genes in columns cells in rows
        
Returns
----------
    label : list of character string
    features : list of character string
"""
def get_label_features(df):
    label = df.index.tolist()
    features = df.columns.tolist()
    return(label,features)



"""
Get the importance value for each feature
Parameters
----------
    model : sklearn model 

    features : list

Returns
----------
    importance_df : pandas.Dataframe
"""
def get_importance(model, features):
    #importance = zip(model.feature_names_in_, model.feature_importances_)   
    importance_df = pd.DataFrame({"genes":features,
                              "importance":model.feature_importances_})
    importance_df.sort_values(by=['importance'], inplace=True, ascending=False)
    
    return importance_df

def get_coefficients(model, features):

    if type(model) == GridSearchCV:
        coef = model.best_estimator_['model'].coef_[0]
    else : 
        coef = model.coef_[0]

    coef_df = pd.DataFrame(data={ 'genes': features, 'Coefficients': coef })

    return coef_df

"""
Reduce the number of features and predict label

Parameters
----------
    model : sklearn model 

    X_train : pandas.Dataframe
        Train dataset
    
    X_test : pandas.Dataframe
        Test dataset

    y_train : list
        list of label in the train dataset

    threshold : str or float, default=None

    n_feat : int, The maximum number of features to select.

Returns
----------
    X_important_train, X_important_test : pandas.Dataframe with only most important features
    selected_feat : character list, most important features
"""
def ROC_analysis(classifiers, X_test, y_test, le_name_mapping, positive_label):

    result_table = pd.DataFrame(columns=['classifiers', 'fpr','tpr','auc'])

    lab_number = le_name_mapping[positive_label]

    for cls in classifiers:      
        yproba = cls.predict_proba(X_test)    
        auc = roc_auc_score(y_test, yproba[::,1]) 
        fpr, tpr, _ = roc_curve(y_test,  yproba[:,lab_number], pos_label=lab_number)
                
        result_table = result_table.append({'classifiers':f'{cls.__class__.__name__}',
                                        'fpr':fpr, 
                                        'tpr':tpr, 
                                        'auc':auc}, ignore_index=True)

    result_table.set_index('classifiers', inplace=True)


    return result_table



def single_ROC(cls, X_test, y_test, le_name_mapping, positive_label):

    result_table = pd.DataFrame(columns=['classifiers', 'fpr','tpr','auc'])

    lab_number = le_name_mapping[positive_label]

    yproba = cls.predict_proba(X_test)    
    auc = roc_auc_score(y_test, yproba[::,1]) 
    fpr, tpr, _ = roc_curve(y_test,  yproba[:,lab_number], pos_label=lab_number)
                
    result_table = result_table.append({'classifiers':f'{cls.__class__.__name__}',
                                        'fpr':fpr, 
                                        'tpr':tpr, 
                                        'auc':auc}, ignore_index=True)

    result_table.set_index('classifiers', inplace=True)


    return result_table

"""
Analyse of a trained model based on its predictions

Parameters
----------
    y_test : pandas.Dataframe
        Train dataset
    
    y_pred : pandas.Dataframe
        Test dataset

Returns
----------
    acc : int, accuracy score
    clf : classification repport
    cm : numpy.array confusion matrix 
"""
def analyse_model(y_test, y_pred ):

    acc = balanced_accuracy_score(y_test, y_pred)
    clf = classification_report(y_test,y_pred, output_dict=True )
    cm = confusion_matrix(y_test, y_pred)
    return acc, clf, cm


def save_model(model, save_dir , filename):
    joblib.dump(model, f'{save_dir}{filename}.joblib')
    

def load_model(file):
    loaded_model = joblib.load(file)
    return loaded_model


def save_cm(cm, file):
    np.savetxt(file, cm, fmt='%d')

def load_cm(file):
    return np.loadtxt(file, dtype=int)

def save_clf(clf,file):
    clf_df = pd.DataFrame(clf)
    clf_df.to_csv(file)

def load_clf(file):
    return pd.read_csv(file, index_col=0)

def save_acc(acc,file):
    acc = pd.DataFrame({"accuracy":[acc]})
    acc.to_csv(file)

def load_acc(file):
    acc = load_clf(file)
    return np.float64(acc["accuracy"])

def save_importance(imp_df, file):
    if type(imp_df) == pd.core.frame.DataFrame : 
        imp_df.to_csv(file)
    else : 
        return None

def load_importance(file):
    imp_df = load_clf(file)
    return imp_df

def save_CV_results(CV_results, file):
    titles = ['Other-vs-CD4','Other-vs-CD8','CD4-vs-CD8']
    for index in range(len(CV_results)) :
        CV_results[index].to_csv(f'{file}_{titles[index]}.csv')        
        
def load_CV_results(path_folder) : 
    files = find_ALL_geneset_file(path_folder, pattern='CV-results')
    pattern = ['Other-vs-CD4','Other-vs-CD8','CD4-vs-CD8']
    CV_results = []
    
    for pat in pattern :
        file = [f for f in files if pat in f][0]
        df = pd.read_csv(file, index_col=0)
        CV_results.append(df)
        
    return CV_results

def save_ROC_values(result_table,file):
    result_df = pd.DataFrame(columns=["fpr","tpr","clf","auc"])
    for i in result_table.index:
        fpr = result_table["fpr"].loc[i]
        tpr = result_table["tpr"].loc[i]
        auc = result_table["auc"].loc[i]
    
        data = pd.DataFrame({"fpr":fpr, "tpr":tpr, "clf":[i]*len(fpr), "auc":auc })
        result_df = pd.concat([result_df,data])

    result_df.to_csv(file)
    

def load_ROC_values(file):
    tmp = pd.read_csv(file, index_col=0)
    result_table = pd.DataFrame(columns=['clf', 'fpr','tpr','auc'])
    for clf in set(tmp["clf"]) : 
        fpr = tmp["fpr"][tmp["clf"] == clf ].to_list()
        tpr = tmp["tpr"][tmp["clf"] == clf ].to_list()
        auc = tmp["auc"][tmp["clf"] == clf ].to_list()[0]
            
        result_table = result_table.append({'clf':f'{clf}','fpr':fpr,'tpr':tpr,'auc':auc},
                                           ignore_index=True)
    
    result_table.set_index('clf', inplace=True)
    return result_table


def save_metrics(acc, clf, cm, imp_df , save_dir, filename, ROC_values=None):
    save_acc(acc, f'{save_dir}accuracy_{filename}.txt' )
    save_cm(cm, f'{save_dir}confusion_heatmap_{filename}.txt' )
    save_clf(clf, f'{save_dir}classification_report_{filename}.txt' )
    save_importance(imp_df, f'{save_dir}importance_df_{filename}.txt')

    if type(ROC_values) == pd.core.frame.DataFrame : 
        save_ROC_values(ROC_values, f'{save_dir}ROC_values_{filename}.txt') 


"""
From one directory, get the path of files corresponding to accuracy,clf,confusion heatmap, importance_df
"""
def get_path_metrics(save_dir):
    files=os.listdir(save_dir)
    acc = get_one_path(save_dir,files,"accuracy") 
    clf = get_one_path(save_dir,files,"classification_report")
    cm = get_one_path(save_dir,files,"confusion_heatmap")
    importance_df = get_one_path(save_dir,files,"importance_df")

    try : 
        ROC_values = get_one_path(save_dir,files,"ROC_values")
    except : 
        ROC_values = None

    return acc, clf, cm, importance_df, ROC_values

"""
load accuracy,clf,confusion heatmap, importance_df from one directory
"""
def load_metrics(save_dir):
    acc_path, clf_path, cm_path, importance_df_path, ROC_path = get_path_metrics(save_dir)
    acc = load_acc(acc_path)
    clf = load_clf(clf_path)
    cm = load_cm(cm_path)
    importance_df = load_importance(importance_df_path)

    try :
        ROC_values=load_ROC_values(ROC_path)
    except :
        ROC_values=None

    return acc, clf, cm, importance_df, ROC_values


"""
Numeric encoding for qualitative label

Parameters
----------
    y_train : list
        list of label in the dataset
 
Returns
----------
    y_encoded : list of numerical value
    le : sklearn.preprocessing.LabelEncoder() object
    le_name_mapping : python dictionary with y label and their numeric
"""
def numeric_encoding(y_labels):
    le = LabelEncoder()
    le.fit(y_labels)
    y_encoded = le.transform(y_labels)
    le_name_mapping = dict(zip(le.classes_, le.transform(le.classes_)))
    return y_encoded, le, le_name_mapping




"""
Produce a dataframe From the result of permutations 

Parameters
----------
    importances_dict : output class from sklearn.inspection.permutation_importance
    
    features : list
        list of features used for permutations
        
    save_dir : string
        path to folder to save in

    filename : character string 

Returns
----------
    result_df : pd.DataFrame
"""
def result_from_permutations(importances_dict, features, filename, save_dir, sort=True ): 
    result_df = pd.DataFrame(importances_dict.importances)
    result_df["mean_importance"] = importances_dict.importances_mean
    result_df["std_importance"] = importances_dict.importances_std
    result_df.index = features

    if sort :
        result_df.sort_values(by="mean_importance", ascending=True, inplace=True)
    #result_df = result_df[result_df["mean_importance"] > 0]
    
    result_df.to_csv(f'{save_dir}/Permutation_result_{filename}.csv')
    
    return result_df



"""
Replace element from list based on element matching

Parameters
----------    
    y : list

    to_replace, replace_by : list of str
Returns
----------
    y : list with replaced element
    
"""
def remap_index(y, to_replace=["2"], replace_by=["1"]):
    if type(y) == np.ndarray :
        y = y.tolist()
    
    y = list(map(str, y))
    y = replace_label(y, to_replace, replace_by, None)
    
    return( list(map(int,y)) )


"""
Remove rows from dataframe / associated element from list corresponding to one of the three classes : CD4, CD8, other

Parameters
----------
    X : pd.DataFrame
    
    y : list

    to_remove : character string in ["CD8", "CD4" , "other" ]

Returns
----------
    X_binary : pd.DataFrame
    y_binary : list
"""
def remove_index(X, y, to_remove = "CD8"):
    
    y = np.asarray(y)
    if to_remove == "CD4" :
        # other : 0 / CD8 : 1
        bool_list = y != 1
        y_binary = y[y != 1]        
        y_binary = remap_index(y_binary, ['2'], ['1'])
        
    elif to_remove == "other" :
        # CD4 : 0 / CD8 : 1
        bool_list = y != 0 
        y_binary = y[y != 0]
        y_binary = remap_index(y_binary, ['1'], ['0'])
        y_binary = remap_index(y_binary, ['2'], ['1'])
    
    elif to_remove == "CD8" : 
        #  other : 0 / CD4 : 1
        bool_list = y != 2
        y_binary = y[y != 2]
        
    X_binary = X[bool_list]

    return(X_binary, y_binary)




"""
Streamline the analyses with One vs One classifier from fitting the model to saving the result  

Parameters
----------
    X_train : pandas.Dataframe
        Train dataset

    X_test : pandas.Dataframe
        Test dataset

    y_train : list
        list of label in the train dataset
 
    y_test : list
        list of label in the test dataset

    save_dir : string
        path to folder to save in

    filename : character string
    
    max_depth, n_estimator : int, argument of RandomForestClassifier
    
    n_perm : int, number of permutation
    
    othervsCD4, othervsCD8, CD4vsCD8 : Boolean, True to compute permutation for one of these classifiers
    
    OvO : default : None, sklearn.multiclass.OneVsOneClassifier type
    
    clf : default : string 'None', sklearn classifier object compatible with sklearn.multiclass.OneVsOneClassifier
    
    scorer : str, or list of str. 
        must be one of sklearn.metrics.get_scorer_names()


Returns
----------
    sklearn.multiclass.OneVsOneClassifier
    dict containing result of permutation for each classifier (other vs CD4, etc.) 
    
"""
def RF_permutation(X_train, y_train, X_test, y_test,filename, save_dir, n_estimator=1000,max_depth=None ,n_perm = 5, CD4vsCD8 = False, 
                    othervsCD8 = True, othervsCD4 = True, OvO = None, clf = "None", scorer="roc_auc_ovo_weighted", Save_metric=True, permutation=True):
        
    if type(OvO) != OneVsOneClassifier : 
        
        if clf == "None" :
        
            OvO = OneVsOneClassifier( n_jobs=-1, 
                                estimator=RandomForestClassifier(n_jobs=-1, 
                                                                n_estimators=n_estimator, max_depth=max_depth,
                                                                random_state=0) )
        else :
            OvO = OneVsOneClassifier( n_jobs=-1, 
                                     estimator= clf)
            
        OvO.fit(X_train, y_train)
        
    save_model(OvO, save_dir, f'RF_1vs1_{filename}' )   

    features = X_train.columns.to_list()
    
    result = {}
    
    index_list = []
    if othervsCD4 :
        index_list.append(0)
    if othervsCD8:
        index_list.append(1)
    if CD4vsCD8 : 
        index_list.append(2)
        
    
    for index in index_list:
        if   index == 0 : 
            condition = "OthervsCD4"
            X_binary, y_binary = remove_index(X_test, y_test, to_remove='CD8')
            le_name_mapping = {'CD4': 1, 'other': 0}
            pos_lab = 'CD4'
  
        elif index == 1 :
            condition = "OthervsCD8"
            X_binary, y_binary = remove_index(X_test, y_test, to_remove='CD4')
            le_name_mapping = {'CD8': 1, 'other': 0}
            pos_lab = 'CD8'

        elif index == 2 :
            condition = "CD4vsCD8"
            X_binary, y_binary = remove_index(X_test, y_test, to_remove='other')
            le_name_mapping = {'CD4': 0, 'CD8': 1}
            pos_lab = 'CD8'
        
        if Save_metric : 
            _, le, _ = numeric_encoding(y_binary)
            model = OvO.estimators_[index]
            y_pred = model.predict(X_binary)

            acc, clf, cm = analyse_model(le.inverse_transform(y_binary), le.inverse_transform(y_pred) )
            metrics = [acc,clf,cm]
            importance_df = get_importance(model, features )  
            save_metrics(
                acc=metrics[0], 
                clf=metrics[1],
                cm=metrics[2],
                imp_df = importance_df ,
                save_dir=save_dir , 
                filename=f'RF_9th_quantile_permutation', 
                ROC_values=ROC_analysis(model, X_binary, y_binary, le_name_mapping, positive_label=pos_lab))
            
            Full_model_evaluation([acc,clf,cm], importance_df, 
                        f"Confusion heatmap {condition}", 
                        f"Classification report heatmap {condition}", 
                        f"Gene importance in classification {condition}",
                        save_dir, f'{filename}_{condition}', 
                        n_genes=len(importance_df[importance_df.importance.cumsum()<=0.95]) )

        if not permutation :   
            return OvO, importance_df
        
        importance_dict = permutation_importance(OvO.estimators_[index], 
                                   X_binary, y_binary, 
                                   n_repeats=n_perm, random_state=0, 
                                   n_jobs=-1, scoring=scorer)
        
        
        if len(scorer) > 1 : 
            score_dict = {}
            for score in scorer : 
                score_dict[score] = result_from_permutations(
                    importance_dict[score], 
                    features,
                    filename=filename+"_"+condition+"_"+score,
                    save_dir=save_dir )
                
            result[ condition ] = score_dict
        
        else : 
            result[ condition ] = result_from_permutations( 
                importance_dict,
                features,
                filename=filename+"_"+condition,
                save_dir=save_dir)
        
    
    return OvO, result



"""
Extract informations from OvO individual estimators (typically a CV estimator like LassoCV or LogisticRegressionCV)

Parameters
----------
    OvO_model : sklearn.multiclass.OneVsOneClassifier
 
Returns
----------
    dictionnary with mse path, best alpha parameters, and all alphas
"""
def get_informations_from_OVO_estimator(OvO_model):
    results_df_dict = {}
    counter = 0
    for model in OvO_model.estimators_ : 
        
        cv_df = pd.DataFrame(model.mse_path_) 
        cv_df['alpha']= model.alphas_
        cv_df['best_alpha']= model.alpha_
        
        results_df_dict[counter] = cv_df
        counter += 1    
        
    return results_df_dict


"""
Save OvO estimators parameters for each estimators : regularization parameters and other parameters 

Parameters
----------
    OvO_model : sklearn.multiclass.OneVsOneClassifier
 
Returns
----------
    dictionnary with mse path, best alpha parameters, and all alphas
"""
def OvO_parameters(OvO_model, filename, save_dir) : 
    titles = ['Other-vs-CD4','Other-vs-CD8','CD4-vs-CD8']
    counter = 0
    for estimator in OvO_model.estimators_ : 
        param = estimator.get_params(False)
        param['min_C'] = min(estimator.Cs_)
        param['max_C'] = max(estimator.Cs_)
        param['best_C'] = estimator.C_[0] 
        df = pd.DataFrame.from_dict(param, orient='index', columns=['value'])
        df.to_csv(f'{save_dir}/Model-parameters_{titles[counter]}_{filename}.csv')
        counter += 1

"""
Lasso Feature selection in a 1v1 model for 3 class classification

Parameters
----------
    X_train : pandas.Dataframe
        Train dataset

    y_train : list
        list of label in the train dataset. Labels must absolutely be :
        0 : other-cell
        1 : CD4
        2 : CD8

    save_dir : string
        path to folder to save in

    filename : character string 

Returns
----------
    OvO : sklearn.multiclass.OneVsOneClassifier
    coef_df : pd.DataFrame with coefficient for each classifier trained in 1v1
    results_dict : dictionary of pd.DaraFrame with CV MSE and alphas/best_alpha for each classification in 1v1
"""
def lassoCV_Feature_reduction(X_train,y_train, filename, save_dir, n_alphas=100,verbose=3,cv=5, save=True, mod_type='linear'):
    
    if mod_type == 'linear' : 
        Model = Pipeline([('scaler',StandardScaler()),
                    ('model',OneVsOneClassifier( n_jobs=-1,
                         estimator=LassoCV(n_jobs=-1,
                                           verbose=verbose, 
                                           n_alphas=n_alphas, 
                                           random_state=0,
                                           cv=cv))) ]).fit(X_train,y_train)
        OvO = Model['model']

        coef_df = pd.DataFrame(
            {"coef_other-CD4":OvO.estimators_[0].coef_,
            "coef_other-CD8": OvO.estimators_[1].coef_,
            "coef_CD4-CD8": OvO.estimators_[2].coef_,
            "genes":X_train.columns })
            
        if save :
            results_dict = get_informations_from_OVO_estimator(OvO)
            save_CV_results(results_dict, file= f'{save_dir}/CV-results_{filename}' )    

    elif mod_type == 'logistic' : 
        Model = Pipeline([('scaler',StandardScaler()),
                    ('model',OneVsOneClassifier(n_jobs=-1,
                        estimator = LogisticRegressionCV(n_jobs=-1, verbose=verbose,
                                                         Cs=n_alphas, cv=cv, max_iter = 250,
                                                         random_state=0, penalty='l1',
                                                         scoring='neg_mean_squared_error', 
                                                         solver="saga",))) ]).fit(X_train,y_train)
        OvO = Model['model']
        Coordinate_descent_Multiplot_logit(OvO, save_dir=save_dir, filename=filename)

        coef_df = pd.DataFrame(
            {"coef_other-CD4":OvO.estimators_[0].coef_[0],
            "coef_other-CD8": OvO.estimators_[1].coef_[0],
            "coef_CD4-CD8": OvO.estimators_[2].coef_[0],
            "genes":X_train.columns })
        
        results_dict = None
        if save :
            OvO_parameters(OvO, filename, save_dir )

    save_model(OvO, save_dir, f'LassoCV_1v1_{filename}' )   
    coef_df = coef_df.loc[coef_df[ ["coef_other-CD4","coef_other-CD8","coef_CD4-CD8"] ].sum(axis=1) != 0 , ]
    coef_df.to_csv(f'{save_dir}/coef-LassoCV_{filename}.csv')
    
    return OvO , coef_df, results_dict


def run_on_reduced_dataset(X, y, feature, le, le_name_mapping, model, save_dir, filename, positive_label=1, ROC = True):
    
    X_train,X_test, y_train,y_test = train_test_split(X,y,test_size=0.3,random_state=1)
    model.fit(X_train,y_train)
    y_pred = model.predict(X_test)
    
    acc, clf, cm = analyse_model(le.inverse_transform(y_test), le.inverse_transform(y_pred) )
    metrics = [acc,clf,cm]
    
    if type(model) == RandomForestClassifier:
        coef_df = get_importance(model, feature)
        n_gene = len(coef_df.index.to_list())
        
    elif type(model) == LogisticRegressionCV:
        coef_df = get_coefficients(model,feature)
        n_gene = len(coef_df[coef_df.Coefficients.cumsum()<=0.95])
    
    if ROC : 
        save_metrics(acc    = metrics[0], 
             clf        = metrics[1],
             cm         = metrics[2],
             imp_df     = coef_df ,
             save_dir   = save_dir, 
             filename   = filename, 
             ROC_values = single_ROC(model, X_test, y_test, le_name_mapping, positive_label=positive_label) )
    else :
        save_metrics(acc    = metrics[0], 
             clf        = metrics[1],
             cm         = metrics[2],
             imp_df     = coef_df ,
             save_dir   = save_dir, 
             filename   = filename, 
             ROC_values = None )
    
    Full_model_evaluation([acc,clf,cm], coef_df, 
                            f"Confusion heatmap {filename}", 
                            f"Classification report heatmap {filename}", 
                            f"Gene importance in classification {filename}",
                            save_dir, f'{filename}' ,
                            n_genes = n_gene)
    
    
    #save_ROC_values(ROC_values, f'{save_dir}{filename}.txt')
    
    save_model(model, save_dir, f'{filename}' )  
    
    return model




def get_single_geneset(Binary_matrix, pattern):
    name = [match for match in Binary_matrix.columns.to_list() if pattern in match]
    if name == [] : 
        return name
    geneset = Binary_matrix[Binary_matrix[ name[0] ] != 0]
    geneset = geneset.index.to_list()
    return(geneset)

def get_geneset_from_binary(Binary_matrix) :
    patterns = ["NK","B","Neutro","Mono","Plasm","CD4","CD8"]
    genesets = {}
    for p in patterns : 
        genesets[p] = get_single_geneset(Binary_matrix, p)
    return(genesets)

def get_existing_index_features(mat, geneset, label, to_remove=["Erythrocyte"] ):
    geneset = [gene for gene in geneset if gene in mat.columns.to_list() ]
    X         = mat[geneset]
    X.index   = label
    X         = X[X.index != "other" ]
    y_uptated = [celltype for celltype in label if celltype != "other" ]
    
    for remove in to_remove:
        X = X[X.index != remove ]
        y_uptated = [celltype for celltype in y_uptated if celltype != remove ]
    
    return(X, y_uptated, geneset )