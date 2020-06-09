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
        
                
unfiltered_sample1=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/VDJ SuPERR-seq/BM1_all_contig_annotations.csv')
unfiltered_sample2=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/VDJ SuPERR-seq/BM2_all_contig_annotations.csv')
unfiltered_sample3=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/VDJ SuPERR-seq/BM3_all_contig_annotations.csv')

sample1=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/VDJ SuPERR-seq/BM1_filtered_contig_annotations.csv')
sample2=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/VDJ SuPERR-seq/BM2_filtered_contig_annotations.csv')
sample3=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/VDJ SuPERR-seq/BM3_filtered_contig_annotations.csv')

adt1=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_BM1.csv',index_col='Unnamed: 0')
adt2=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_BM2.csv',index_col='Unnamed: 0')
adt3=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_BM3.csv',index_col='Unnamed: 0')

ADT=[adt1,adt2,adt3];
SAMPLE=[sample1,sample2,sample3]
Unfiltered_sample=[unfiltered_sample1,unfiltered_sample2,unfiltered_sample3]
VDJ=[];
EXTRA_VDJ=[];
for i in range(len(SAMPLE)):
    SAMPLE[i]=SAMPLE[i][["barcode","umis"]].groupby("barcode").sum()
    SAMPLE[i].columns=['Total_Ig_trans']
    Unfiltered_sample[i]=Unfiltered_sample[i][["barcode","umis"]].groupby("barcode").sum()
    vdj=SAMPLE[i].index.tolist()    
    extra_vdj=list(set(Unfiltered_sample[i].index.tolist())-set(vdj)) 
    Unfiltered_sample[i]=Unfiltered_sample[i].loc[extra_vdj]
    modify_name(vdj,i+1)
    modify_name(extra_vdj,i+1)
    SAMPLE[i].index=vdj;
    Unfiltered_sample[i].index=extra_vdj
    back_vdj=cp.deepcopy(vdj)
    back_extra_vdj=cp.deepcopy(extra_vdj)
    remove_missing(vdj,ADT[i].columns.tolist(),back_vdj)
    remove_missing(extra_vdj,ADT[i].columns.tolist(),back_extra_vdj)   
    VDJ.append(vdj)
    EXTRA_VDJ.append(extra_vdj)
    SAMPLE[i]=SAMPLE[i].loc[vdj]
    Unfiltered_sample[i]=Unfiltered_sample[i].loc[extra_vdj]


row1 = 'Productive_VDJ'
row2 = "Non_productive"

for i in range(len(ADT)):
    col = ADT[i].columns.tolist()
    New = pd.DataFrame(np.zeros((2,len(col))),index=[row1,row2],columns=col);
    New.loc["Productive_VDJ",VDJ[i]]=1
    New.loc["Non_productive",EXTRA_VDJ[i]]=1
    SAMPLE[i]=pd.concat([SAMPLE[i],Unfiltered_sample[i]])
    SAMPLE[i]=SAMPLE[i].T
    ADT[i]=pd.concat([ADT[i],New,SAMPLE[i]])
    ADT[i]=ADT[i].fillna(0)
    ADT[i].to_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature/DSB_adt+vdj_BM'+str(i+1)+'.csv')




