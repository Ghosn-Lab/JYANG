# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 13:27:54 2020

@author: yjk12
"""
import os
import pandas as pd
import numpy as np
os.chdir('F:/Ghosn_lab/supperseq')

file1=pd.read_csv('F:/Ghosn_lab/supperseq/PBMC1.csv',index_col='Unnamed: 0')
file2=pd.read_csv('F:/Ghosn_lab/supperseq/PBMC2.csv',index_col='Unnamed: 0')
file3=pd.read_csv('F:/Ghosn_lab/supperseq/PBMC3.csv',index_col='Unnamed: 0')

file4=pd.read_csv('F:/Ghosn_lab/supperseq/BM1.csv',index_col='Unnamed: 0')
file5=pd.read_csv('F:/Ghosn_lab/supperseq/BM2.csv',index_col='Unnamed: 0')
file6=pd.read_csv('F:/Ghosn_lab/supperseq/BM3.csv',index_col='Unnamed: 0')
#def concatenate_matrix(frame1,frame2):

pbmc=pd.concat([file1,file2,file3],join='inner',axis=1)
BM=pd.concat([file2,file3],join='inner',axis=1)

pbmc.to_csv('con_PBMC.csv')
BM.to_csv('con_BM.csv')
