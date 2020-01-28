# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:03:57 2020

@author: yjk12
"""

import pandas as pd;
import numpy as np;

def CLR_across_cell(matrix):
    matrix=np.log(1+matrix)-(np.log(1+matrix).sum(axis=1,keepdims=True))/matrix.shape[1]
    return matrix

def CLR_across_gene(matrix):
    matrix=np.log(1+matrix)-(np.log(1+matrix).sum(axis=0,keepdims=True))/matrix.shape[0]
    return matrix

def pos_CLR_across_cell(matrix):
    t=np.exp((np.log1p(matrix).sum(axis=1,keepdims=True))/matrix.shape[1])
    matrix=np.log1p(matrix/t)
    return matrix

def pos_CLR_across_gene(matrix):
    t=np.exp((np.log1p(matrix).sum(axis=0,keepdims=True))/matrix.shape[0])
    matrix=np.log1p(matrix/t)
    return matrix

### read file and transform to numpy objects
file=pd.read_csv('con_PBMC1+2+3.csv',index_col='Unnamed: 0');
mat=file.to_numpy();

### CLR transformation
origin=CLR_across_gene(mat)
positive=pos_CLR_across_gene(mat)

### Transform numpy objects to Dataframe 
origin=pd.DataFrame(origin,index=file.index,columns=file.columns)
positive=pd.DataFrame(positive,index=file.index,columns=file.columns)


### write to csv
origin.to_csv('CLR_PBMC.csv')
positive.to_csv('CLR_positive_PBMC.csv')