# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:42:14 2019

@author: yjk12
"""

import pandas as pd
import csv
import numpy as np
import scipy.io
import os
from scipy.sparse import csr_matrix

d = os.getcwd().split('\\')
d = d[-1]

with open('genes.tsv','r+') as f1:
    gene= [row[0] for row in csv.reader(f1, delimiter="\t")]
    for i in range(0,len(gene)-1):
        gene[i]=gene[i].split('-')
        gene[i]=gene[i][0]

with open('barcodes.tsv','r') as f2:
    barcodes = [row[0] for row in csv.reader(f2, delimiter="\t")]
    for i in range(0,len(barcodes)):
        barcodes[i] = barcodes[i]+'_'+str(d[-1])
 
mat= scipy.io.mmread('matrix.mtx')
mat=mat.tocsr()
mat=mat.todense()
mat=pd.DataFrame(data=mat,columns=barcodes,index=gene)
mat=mat.drop(['unmapped'],axis=0)
#mat=mat.drop(['-','unmapped'],axis=0)
mat.to_csv(d+'.csv',sep=',')