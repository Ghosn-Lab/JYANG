# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 13:27:54 2020

@author: yjk12
"""
import os
import pandas as pd
import numpy as np
os.chdir('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq')





file1=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_adt+vdj_PBMC1_filtered.csv',index_col='Unnamed: 0')
file2=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_adt+vdj_PBMC2_filtered.csv',index_col='Unnamed: 0')
file3=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_adt+vdj_PBMC3_filtered.csv',index_col='Unnamed: 0')

file4=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature/DSB_adt+vdj_PBMC1.csv',index_col='Unnamed: 0')
file5=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature/DSB_adt+vdj_BM2.csv',index_col='Unnamed: 0')
file6=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature/DSB_adt+vdj_BM3.csv',index_col='Unnamed: 0')

file1=file1.T
file2=file2.T
file3=file3.T

file4=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_adt+vdj_PBMC1.csv',index_col='Unnamed: 0')
file5=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_adt+vdj_PBMC2.csv',index_col='Unnamed: 0')
file6=pd.read_csv('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_adt+vdj_PBMC3.csv',index_col='Unnamed: 0')
file4=file4.T
file5=file5.T
file6=file6.T
#def concatenate_matrix(frame1,frame2):

pbmc=pd.concat([file1,file2,file3],join='inner',axis=1)
BM=pd.concat([file5,file6],join='inner',axis=1)
DSB_pbmc=pd.concat([file4,file5,file6],join='inner',axis=1)

pbmc.to_csv('DSB_adt+vdj_PBMC_con.csv')
pbmc.to_csv('con_PBMC_force.csv')
BM.to_csv('DSB_BM2+3_vdj.csv')
DSB_pbmc.to_csv("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature/DSB_adt+vdj_BM2+3.csv")
