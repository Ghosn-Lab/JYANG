# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 10:33:12 2020

@author: yjk12
"""

import pandas as pd
import copy as cp

def modify_name(name,index):
    for i in range(len(name)):
        name[i]=name[i][0:-2]+'_'+str(index)
def remove_missing(vdj,adt,back_vdj):
    for i in back_vdj:
        if i not in adt:
            vdj.remove(i)

unfiltered=pd.read_csv('C:/Users/yjk12/Box/Data/VDJ SuPERR-seq/PBMC2_all_contig_annotations.csv')
adt=pd.read_csv('C:/Users/yjk12/Box/Data/supperseq/CLR_adt+vdj_PBMC2_filtered.csv',index_col='Unnamed: 0')
adt=adt.T

### B cells in filtered VDJ
fB=adt[adt["is_B"]==1].index.tolist()
nfB=unfiltered[["barcode","umis"]].groupby("barcode").sum()
nfB=nfB.index.tolist()

modify_name(nfB,2)

back_vdj=cp.deepcopy(nfB)
remove_missing(nfB,adt.index.tolist(),back_vdj)

### B cells only belong to unfiltered VDJ
ufB=list(set(nfB).difference(set(fB)))

### cells not belong to unfiltered VDJ
nB=list(set(adt.index.tolist()).difference(set(nfB)))


fB=adt[adt["is_B"]==1]
ufB=adt.loc[ufB]
nB=adt.loc[nB]

fB.to_csv("C:/Users/yjk12/Box/Data/VDJ SuPERR-seq/B_cell/B_PBMC2.csv")
ufB.to_csv("C:/Users/yjk12/Box/Data/VDJ SuPERR-seq/B_unfiltered_only/unfiltered_B_PBMC2.csv")
nB.to_csv("C:/Users/yjk12/Box/Data/VDJ SuPERR-seq/Non_B_unfiltered/non_B_PBMC2.csv")
