################################################################################                                                                           
# > february 2023                                                                                                             
# > Script : utils.py                                                                                                           
# > Function : Used to get a reservoir of utilitary functions
# @ COLAJANNI Antonin                                                          
################################################################################

import pandas as pd        
import numpy as np
import os, re, os.path


def save_geneset(selected_feat, save_dir, filename):
        reduced_N_feat = len(selected_feat) 
        with open(f'{save_dir}gene_set_{reduced_N_feat}_{filename}.txt', 'w') as filehandle:
                for gene in selected_feat:
                        filehandle.write('%s\n' % gene)
                        

def load_geneset(file):
    geneset = []
    with open(file, 'r') as filehandle:
        for line in filehandle:
            curr_genes = line[:-1]
            geneset.append(curr_genes)
    return(geneset)


def find_ALL_geneset_file(path, pattern = "gene_set|geneset|common_genes"):
    #pattern = "gene_set|geneset|common_genes"
    file_list = []
    for root, dirs, files in os.walk(path):
        for file in filter(lambda x: re.match(pattern, x), files):
            file_list.append(root+"/"+file )
            
    return file_list

def find_geneset_file(path, pattern):
    files = find_ALL_geneset_file(path)
    for pat in pattern : 
        files = [f for f in files if pat in f]
        
    #path = f'{path}/{files[0]}'
    return files[0]

"""
Extract one full path from a directory matching a pattern

Parameters
----------
    save_dir : string
        path to folder to save in

    filename : character string 

    pattern : character string or list of character string

Returns
----------
    character string, path
"""
def get_one_path(save_dir,file_list,pattern):
    if type(pattern) == str : 
        pattern = [pattern]
    for pat in pattern : 
        file_list = [f for f in file_list if pat in f]
    path = f'{save_dir}/{file_list[0]}'
    return path


"""
Replace elements in list if substring is contained in an element of the list 

Parameters
----------
    lab_list : list of character string

    pattern_list : list of str to match
        should be the same length as replace_list 

    replace_list : list of str to replace in lab_list of it matches one of pattern_list
        should be the same length as replace_pattern

    replace_other : string
        string to replace element if it doesn't match anything
        replace_other should be None if you don't want to replace those values

Returns
----------
    list of character string
"""
def replace_label(lab_list, pattern_list, replace_list, replace_other = "other_cell"):

    for i, item in enumerate(lab_list):
        c = 0

        for index, pattern in enumerate(pattern_list) : 

            if pattern in item :
                lab_list[i] = replace_list[index]
                continue 
            c+=1
            if replace_other != None and c == len(pattern_list):
                lab_list[i] = replace_other

    return(lab_list)



def create_dir(directory):
    # checking if the directory demo_folder 
    # exist or not.
    if not os.path.exists(directory):
        # if the demo_folder directory is not present 
        # then create it.
        os.makedirs(directory)
    return(directory)