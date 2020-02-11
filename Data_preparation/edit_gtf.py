# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:17:21 2020

@author: yjk12
"""

import os
import argparse

path='F:/Ghosn_lab/Devon/Adding sequence/Fluorescent reporters/'   ### Folder only contains new fastq files
files=os.listdir(path)

Dic={}
for i in files:
    file=open(path+i,'r')
    title=file.readline()
    title=title.split()[0]
    seq=file.read()
    file.close()
    seq=seq.replace('\n','')
    Dic[title]=len(seq)

source='Junkai'
feature='exon'
Start='1'
Score='.'
strand='+'
Frame='.'


#%% Write to gtf
file=open('add_gene.gtf','w')
for key in Dic.keys():
    Attributes= 'gene_id "'+key[1:]+'"; transcript_id "'+key[1:]+'";\n'
    out='\t'.join([key,source,feature,Start,str(len(seq)),Score,strand,Frame,Attributes])
    file.write(out)
file.close()
    
