#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 16:40:40 2019

@author: Laurentijn Tilleman
"""

# packages
import pandas as pd
import re

## panels
# IDT
idt = set(pd.read_csv('data/PanCancer V2 gene list.csv', sep=' *, *')['Gene List'])

# Roche
expanded = set(pd.read_csv('data/expanded.txt', sep=' *, *',header=None)[0])
surveillance = set(pd.read_csv('data/surveillance.txt', sep=' *, *',header=None)[0])
target = set(pd.read_csv('data/target.txt', sep=' *, *',header=None)[0])

# Foundation Medicine
oneL = set(pd.read_csv('data/FoundationOne.txt', sep=' *, *',header=None)[0])
oneCDx = set()
with open('data/FoundationOneCDx.txt','r') as fh:
    for line in fh.readlines():
        for gene in re.sub("\(.*?\)",'',line).replace('*','').split():
            oneCDx.add(gene)

# qiagen
qiagen = set(pd.read_csv("data/qiagen_compCancer_list.txt",header=None)[0])

# illumina
illumina = set(pd.read_csv('data/illuminaTSO500_list.txt',header=None)[0])

# Thermofisher
ampliseq = set(pd.read_csv('data/CCP.20131001.submitted.bed',sep='\t',
                skiprows=1,header=None)[5])

# panels
panels = [idt,expanded,surveillance,target,oneL,oneCDx,illumina,qiagen,ampliseq]
# panel names
panelNames = ['Total','IDT','Avenio Target','Avenio Expanded','Avenio Surveillance',
              'Foundation One Liquid','Foundation One CDx','Illumina','Qiagen',
              'AmpliSeq']

data = pd.DataFrame(index=panelNames[1:],columns = panelNames[1:])

for i,x in enumerate(panels):
    for j,y in enumerate(panels):
        data.iloc[i,j] = len(set(x) & set(y))

data.to_csv('results/panelCompare.csv')
