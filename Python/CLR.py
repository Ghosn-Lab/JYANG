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
file=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC2_3500.csv',index_col='Unnamed: 0');
file=file.T
preserve=file.apply(lambda x: x.sum(),axis=1)>0
file=file.loc[preserve,]
#file=file.drop(index="CD38_TotalSeqC")
mat=file.to_numpy();

### CLR transformation
origin=CLR_across_gene(mat)
positive=pos_CLR_across_gene(mat)

### Transform numpy objects to Dataframe 
origin=pd.DataFrame(origin,index=file.index,columns=file.columns)
positive=pd.DataFrame(positive,index=file.index,columns=file.columns)


### write to csv
origin.to_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/CLR_PBMC2_3500.csv')
#positive.to_csv('C:/Users/yjk12/Box/Data/supperseq/CLR_positive_Cite_BM3.csv')
