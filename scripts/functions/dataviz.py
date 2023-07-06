################################################################################                                                                           
# > December 2022                                                                                                             
# > Script : dataviz.py                                                                                                            
# > Function : list of function used to visualize data/results                                 
# @ COLAJANNI Antonin                                                          
################################################################################


#import matplotlib

import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
import scikitplot as skplt
import seaborn as sns 
import pandas as pd
import numpy as np
from collections import Counter


def conf_mat_percentage(cm, classes):
    cm = pd.DataFrame(cm, columns = classes, index = classes)
    cm_percent = cm.div(cm.sum(axis=1), axis=0)
    np.asarray(cm_percent).flatten()    
    return cm_percent


"""
Plot a confusion heatmap and save the image
Parameters
----------
    cm : pandas.dataframe 

    title : str
        plot title

    save_dir : str
        save directory

    filename : str

    colors : str
        color palette
    
    acc : int (optional)
        accuracy score
        0 <= int <= 1

Returns
----------
    None
"""
def confusion_heatmap(cm, classes, title, save_dir, filename ,acc = None, colors = "rocket"):
    #if int(cm.sum(axis=1).sum(axis=0)) == len(cm):
    #    precision = ".2%"
    #else :
    #    precision = ".5g"

    # Print label both in %age and true number
    cm_percent = conf_mat_percentage(cm, classes)

    group_counts = ["{0:0.0f}".format(value) for value in cm.flatten()]
    group_percentages = ["{0:.2%}".format(value) for value in np.asarray(cm_percent).flatten()]
    labels = [f"{v1}\n{v2}" for v1, v2, in zip(group_percentages,group_counts)]
    labels = np.asarray(labels).reshape(len(classes),len(classes))
    
    
    plt.figure(figsize = (12,12))
    g = sns.heatmap(cm_percent, annot=labels, fmt="", cmap=colors,
    #g = sns.heatmap(cm, annot=True, cmap = colors ,fmt=precision, 
           linewidths=2,
           square=True, cbar=True,
           cbar_kws={"shrink": 0.5})
    # rotation des lables / changement de taille de police
    plt.setp(g.get_yticklabels(), rotation=0, fontsize = 14)
    plt.setp(g.get_xticklabels(), rotation=45, fontsize = 14)

    if acc != None : 
        title = f'{title} \n Balanced accuracy : {acc:.2%}'

    g.set_xlabel('Predicted label', fontsize=16)
    g.set_ylabel('True label', fontsize=16)
    plt.title(title, fontsize = 16)
    plt.tight_layout()
    plt.savefig(save_dir+filename+".png", dpi=500)
    plt.clf()


"""
Plot a classification heatmap and save the image
Parameters
----------
    clf : dictionnary
        output of sklearn.metrics.classification_report

    title : str
        plot title

    save_dir : str
        save directory

    filename : str

    colors : str
        color palette

Returns
----------
    None
"""
def plot_clf(clf, title, save_dir, filename,colors = "gist_heat" ):

    plt.figure(figsize = (12,12))
    # .iloc[:-1, -3:] enlver le support et les moyennes
    g = sns.clustermap(pd.DataFrame(clf).iloc[:-1, :-3].T, annot=True, cmap = colors, col_cluster = False )

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=0, fontsize = 14)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize = 14)
    plt.title(title, fontsize = 18)
    plt.tight_layout()
    plt.savefig(save_dir+filename+".png", dpi=500)
    plt.clf()

"""
Plot Features importance barplot and save the image
Parameters
----------
    importance_df : dictionnary
        output of sklearn.metrics.classification_report

    title : str
        plot title

    save_dir : str
        save directory

    filename : str

    n_genes : int
        number of genes to keep in the plot

Returns
----------
    None
"""
def Importance_barplot(importance_df,title, save_dir, filename, n_genes=50):
    extension = ".png"
    if type(importance_df) != pd.core.frame.DataFrame : 
        return None

    importance_df = importance_df.head(n_genes)
    tot = sum(importance_df.head(n_genes).importance)

    plt.figure(figsize=(10,12))
    # make barplot and sort bars
    sns.barplot(y='genes',
            x="importance", 
            color = "steelblue",
            data=importance_df, 
            order=importance_df.sort_values('importance',ascending=False).genes)


    # set labels
    plt.xlabel("Importance value", size=15)
    plt.ylabel("Genes", size=15)
    plt.title(f'{title} \n Importance sum : {tot:.2%}', size=18)
    plt.tight_layout()
    plt.savefig(save_dir+filename+extension, dpi=500)
    plt.clf()


"""
Plot a Coefficients barplot and save the image (regression)
Parameters
----------
    importance_df : dictionnary
        output of sklearn.metrics.classification_report

    title : str
        plot title

    save_dir : str
        save directory

    filename : str
    
Returns
----------
    None
"""
def Coef_barplot(coef_df, title, save_dir, filename):
    extension = ".png"
    coef_df = coef_df.sort_values(by='Coefficients', ascending=False)
    if len(coef_df) > 1000 : 
        coef_df = coef_df[coef_df["Coefficients"] != 0] 
        title = title+" 0 filtered "

    plt.figure(figsize = (10,25))
    plt.barh(y=coef_df['genes'], width=coef_df['Coefficients'], color='#087E8B')
    plt.xticks(rotation='vertical', size = 10)

    plt.title(title+f'\n {len(coef_df)} features ', size = 16)
    plt.savefig(save_dir+filename+extension, dpi=500)
    plt.clf()

"""
Plot confusion heatmap, classification report, and barplot of importance

Parameters
----------
    importance_df : dictionnary
        output of sklearn.metrics.classification_report

    title : str
        plot title

    save_dir : str
        save directory

    filename : str

    n_genes : int
        number of genes to keep in the plot

Returns
----------
    None
"""
def Full_model_evaluation(evaluation_metrics, importance_df, title_confusion, title_report, title_barplot, save_dir, filename, n_genes=50, colors = "rocket" ):

    classes = []
    not_wanted = ['macro avg', 'weighted avg','accuracy']
    for key in evaluation_metrics[1].keys() : 
        if key not in not_wanted : 
            classes.append(key)

    confusion_heatmap(evaluation_metrics[2], classes ,title_confusion, save_dir, f'{filename}_conf_heatmap',evaluation_metrics[0],colors=colors )
    plot_clf(evaluation_metrics[1], title_report, save_dir, f'{filename}_classfication_report')
    try :
        Importance_barplot(importance_df,title_barplot, save_dir, f'{filename}_Top_50_importance_genes', n_genes= n_genes)
    except:
        Coef_barplot(importance_df, title_barplot, save_dir, f'{filename}_Coefficients')



"""
Barplot counting different label in a list

Parameters
----------
    label_list : list

    title : str
        plot title

    save_dir : str
        save directory

    filename : str

Returns
----------
    None
"""
def label_barplot(label_list, title, save_dir, filename):
    count = dict(Counter(label_list))
    count= pd.DataFrame({"Cell_types":count.keys(), "count":count.values()})
    plt.figure(figsize = (10,15))
    splot = sns.barplot(x = "Cell_types", y = "count", data=count)
    for g in splot.patches:
        splot.annotate(format(g.get_height(), '.0f'),
                   (g.get_x() + g.get_width() / 2., g.get_height()),
                   ha = 'center', va = 'center',
                   xytext = (0, 15),
                   textcoords = 'offset points',
                   size=13)
    plt.xlabel("Cell Types", size = 16)
    plt.ylabel("Occurence", size = 16)
    plt.xticks(rotation=90, size = 12)
    plt.tight_layout()
    plt.title(f'{title} \n Cell number : {len(label_list)}', size=18)
    plt.savefig(save_dir+filename, dpi=500)
    plt.clf()


"""
ROC curve for 3 classes classification

Parameters
----------
    y_true : y true label / list

    y_proba : output of clf.predict_proba

    label_encoder : fitted sklearn.preprocessing.LabelEncoder()

    save_dir : str
        save directory

    filename : str

Returns
----------
    None
"""
def ROC_curve_3classes(y_true, y_proba, label_encoder, save_dir, filename):
    extension = ".png"
    plot = skplt.metrics.plot_roc(label_encoder.inverse_transform(y_true),
                              y_proba, 
                              plot_micro=False, plot_macro=False,
                              figsize=(10,10))

    plot.figure.savefig( fname = f'{save_dir}{filename}_ROC_curve{extension}', dpi=500)
    plt.clf()

"""
ROC curve for 2 classes classification

Parameters
----------
    result_table : pd.DataFrame
        4 column : 
            classifiers : name of the used classifier
            fpr : false positive rate output of clf.predict_proba()
            tpr : false positive rate output of clf.predict_proba()
            auc : roc_auc_score() value

    save_dir : str
        save directory

    filename : str

Returns
----------
    None
"""
def ROC_curve(result_table, save_dir, filename, title=""):
    extension = ".png"

    plt.figure(figsize=(8,6))
    for i in result_table.index:
        plt.plot(result_table.loc[i]['fpr'], 
                result_table.loc[i]['tpr'], 
                label="{}, AUC={:.3f}".format(i, result_table.loc[i]['auc']))
    
    plt.plot([0,1], [0,1], color='orange', linestyle='--')

    plt.xticks(np.arange(0.0, 1.1, step=0.1))
    plt.xlabel("False Positive Rate", fontsize=15)

    plt.yticks(np.arange(0.0, 1.1, step=0.1))
    plt.ylabel("True Positive Rate", fontsize=15)

    plt.title(f'ROC Curve Analysis \n {title}', fontweight='bold', fontsize=15)
    plt.legend(prop={'size':13}, loc='lower right')

    plt.savefig(f'{save_dir}{filename}_ROC_curve{extension}', dpi=500)
    plt.clf()
    


"""
Function used internally to produce figures of MSe paths accross alpha

Parameters
----------
    results_dict : dictionary of pd.DaraFrame with CV MSE and alphas/best_alpha for each classification in 1v1

Returns
----------
    results : dictionary with 3 elements mse_paths,best_alpha,alphas that contains list of int
"""
def get_results_from_OvO_dict(result_dict) : 
    mse_col = result_dict[0].axes[1][0: len(result_dict[0].axes[1]) - 2]
    
    results = {'mse_path' : [],
               'best_alpha':[],
               'alphas' : []}
    
    for index in range(len(result_dict)):
        results['best_alpha'].append(result_dict[index]['best_alpha'][0])
        results['alphas'].append(result_dict[index]['alpha'])
        results['mse_path'].append(np.asarray(result_dict[index][mse_col]))

    return(results)


"""
Mean square error path plot accross alpha for Lasso Cross Validation

Parameters
----------
    mse_path : np.array
    
    alphas : list
    
    best_alpha : alpha value for which R² score is maximal
    
    ax : matplotlib.axes._subplots.AxesSubplot
    
    title : string
    
    y_label, xlabel,legend : boolean, whether or not to show labels/plot

Returns
----------
    None
"""
def Coordinate_descent_plot(mse_path, alphas, best_alpha ,ax, title, ylabel = True, xlabel=True, legend=True):
    #ax = plt.figure()
    
    ax.semilogx(alphas, mse_path, linestyle=":",linewidth=1)

    ax.plot(
        alphas,
        mse_path.mean(axis=-1),
        color="black",
        label="Average across the folds",
        linewidth=1,
    )
    ax.axvline(best_alpha, linestyle="--", color="red", label="alpha: CV estimate")
    ax.errorbar(alphas, 
                mse_path.mean(axis=-1), 
                yerr=mse_path.std(axis=-1), 
                linestyle='None', marker='.', color='black')

    
    if xlabel : 
        ax.set_xlabel(r"$\lambda$")
    if ylabel :
        ax.set_ylabel("Mean square error")
    if legend :
        ax.legend()
        
    _ = ax.set_title(title) 

    return ax

"""
Mean square error path plot accross alpha for Lasso Cross Validation for 3 plots : 
Used for 1v1 classifier in 3 class classification

Parameters
----------
    results_df_dict : dictionary of dataframes. Object obtained with lassoCV_Feature_reduction
        Dict with 3 index for each of the 3 trained model.
        Each dataframe has n_col for each cross-validation step, +1 column for alphas, +1 for best alpha value.
        
        
    save_dir : string. path to save the plot

    filename: string. name of the file to create
    
    save : boolean. True to save the plot

Returns
----------
    None
"""
def Coordinate_descent_Multiplot(results_df_dict, save_dir, filename, save=True ):
    
    results = get_results_from_OvO_dict(results_df_dict)
                             
    fig, axs = plt.subplots(3,1, sharey=False, sharex=True, figsize=(8,10) )
    titles = ['Other-vs-CD4','Other-vs-CD8','CD4-vs-CD8']

    fig.suptitle("Mean square error on each fold: coordinate descent")

    Coordinate_descent_plot(results["mse_path"][0], results["alphas"][0], results['best_alpha'][0], axs[0], titles[0],ylabel = False, xlabel=False, legend=True )
    Coordinate_descent_plot(results["mse_path"][1], results["alphas"][1], results['best_alpha'][1], axs[1], titles[1],ylabel = True, xlabel=False, legend=False)
    Coordinate_descent_plot(results["mse_path"][2], results["alphas"][2], results['best_alpha'][2], axs[2], titles[2],ylabel = False, xlabel=True, legend=False)
    
    if save :
        fig.figure.savefig(f"{save_dir}/Coordinate_descent_{filename}.png",dpi=900)
        
        
def Coordinate_descent_plot_logit(scores, Cs, best_C ,ax, title, ylabel = True, xlabel=True, legend=True):
    for i in range(len(scores[1])):
        ax.plot(Cs, scores[1][i], linestyle=":", linewidth=1)
    
    ax.plot(Cs, scores[1].mean(axis=0),
        color="black",
        label="Average across the folds",
        linewidth=1, )
    
    ax.axvline(best_C, linestyle="--", color="red", label=r"1 / $\lambda$ : CV estimate")
    ax.errorbar(Cs, scores[1].mean(axis=0), 
                yerr=scores[1].std(axis=0), 
                linestyle='None', marker='.', color='black')
    
    if xlabel : 
        ax.set_xlabel(r"C = 1 / $\lambda$")
    if ylabel :
        ax.set_ylabel("Negative mean squared error")
    if legend :
        ax.legend()

    _ = ax.set_title(title) 
    return ax

def Coordinate_descent_Multiplot_logit(OvO_model,  save_dir, filename, save=True ):
    
    model_0 = OvO_model.estimators_[0] 
    model_1 = OvO_model.estimators_[1] 
    model_2 = OvO_model.estimators_[2] 
     
    fig, axs = plt.subplots(3,1, sharey=False, sharex=True, figsize=(8,10) )
    
    for i in range(len(axs)) :
        axs[i].set_xscale('log')
        axs[i].xaxis.set_major_locator(mticker.LogLocator(numticks=999))
        axs[i].xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))

    titles = ['Other-vs-CD4','Other-vs-CD8','CD4-vs-CD8']
    fig.suptitle("Negative mean squared error on each fold: coordinate descent")

    Coordinate_descent_plot_logit(model_0.scores_, model_0.Cs_, model_0.C_, axs[0], titles[0],ylabel = False, xlabel=False, legend=True )
    Coordinate_descent_plot_logit(model_1.scores_, model_1.Cs_, model_1.C_, axs[1], titles[1],ylabel = True, xlabel=False, legend=False)
    Coordinate_descent_plot_logit(model_2.scores_, model_2.Cs_, model_2.C_, axs[2], titles[2],ylabel = False, xlabel=True, legend=False)
        
    if save :
        fig.figure.savefig(f"{save_dir}/Coordinate_descent_{filename}.png",dpi=900)