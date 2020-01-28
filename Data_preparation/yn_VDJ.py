# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 12:08:39 2020

@author: yjk12
"""

import pandas as pd
import numpy as np
import copy as cp
def modify_name(name,index):
    for i in range(len(name)):
        name[i]=name[i][0:-2]+'_'+str(index)
def remove_missing(vdj,adt,back_vdj):
    for i in back_vdj:
        if i not in adt:
            vdj.remove(i)
        
                

vdj1=pd.read_csv('PBMC1_filtered_contig_annotations.csv')
vdj2=pd.read_csv('PBMC2_filtered_contig_annotations.csv')
vdj3=pd.read_csv('PBMC3_filtered_contig_annotations.csv')

adt1=pd.read_csv('F:/Ghosn_lab/supperseq/CLR_PBMC.csv',index_col='Unnamed: 0')




vdj1=list(set(vdj1['barcode']))
vdj2=list(set(vdj2['barcode']))
vdj3=list(set(vdj3['barcode']))

modify_name(vdj1,1)
modify_name(vdj2,2)
modify_name(vdj3,3)

vdj=vdj1+vdj2+vdj3
back_vdj=cp.deepcopy(vdj)
remove_missing(vdj,adt1.columns.tolist(),back_vdj)

row = ['is_B']
col = adt1.columns.tolist()
VDJ = pd.DataFrame(np.zeros((1,len(col))),index=row,columns=col);
VDJ[vdj]=1

ADT=pd.concat([adt1,VDJ])

ADT.to_csv('F:/Ghosn_lab/supperseq/CLR_adt+vdj_PBMC_filtered.csv')




