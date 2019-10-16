# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:38:09 2019

@author: yjk12
"""
import pandas as pd
import csv
import numpy as np
import scipy.io
filename= input('mtx file name:')
with open('genes.tsv','r') as f1:
    gene= [row[0] for row in csv.reader(f1, delimiter="\t")]
with open('barcodes.tsv','r') as f2:
    barcodes = [row[0] for row in csv.reader(f2, delimiter="\t")]
with open('genes.tsv','r') as f3:
    gene_types = [row[2] for row in csv.reader(f3, delimiter="\t")]   
mat= scipy.io.mmread(filename)     
mat=mat.tocsr()

gene_index=[]
gene_name=[]
for i in range(len(gene_types)):
    if gene_types[i].startswith('Anti'):
        gene_index.append(i)
        gene_name.append(gene[i])
        
matrix=np.zeros((len(gene_index),len(barcodes)))
for i in range(len(gene_index)):
    for j in range(len(barcodes)):
        matrix[i,j]=mat[gene_index[i],j]
BEI=pd.DataFrame(index=gene_name,data=matrix)
BEI.columns=barcodes
BEI.to_csv('ADT.csv',sep='\t')
