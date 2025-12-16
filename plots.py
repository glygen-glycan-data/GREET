#!/bin/env python3.12

import sys, os
import pandas as pd

df = pd.read_csv(sys.argv[1],sep='\t')

import matplotlib
from matplotlib.pyplot import *
import seaborn as sns

allsampletypes = set(df['SampleType'])
celltypes = set(filter(lambda s: s.startswith('Celltype:'),allsampletypes))
tissues = set(filter(lambda s: s.startswith('Tissue:'),allsampletypes))
celltypeintissue = allsampletypes.difference(celltypes).difference(tissues)
celltypesandtissues = celltypes.union(tissues)
endolymphatic = set(filter(lambda s: 'Endothelial cell (lymphatic)' in s,allsampletypes))
endovascular = set(filter(lambda s: 'Endothelial cell (vascular)' in s,allsampletypes))
mucouscell = set(filter(lambda s: 'Esophagus mucosa' in s,allsampletypes))
# Endothelial cell (lymphatic) in Esophagus muscularis
endomuscularis = set(filter(lambda s: 'Endothelial cell (lymphatic) in Esophagus muscularis' == s,allsampletypes))
endomuscularis = set(filter(lambda s: 'Endothelial cell (lymphatic) in Esophagus mucosa' == s,allsampletypes))

def fixcelltype(s):
    if ' in ' in s:
        return s.split(' in ')[0]
    return "All Cell-types"
 

def baseplot(df,gene,sampletypes,sortkey='SampleType',sortasc=True,hidexticks=False):
    df1 = df
    df1 = df1[df1['Genes'] == gene]
    df1 = df1[df1['SampleType'].isin(sampletypes)]
    df1 = df1.drop(columns=['GeneSet'])
    df1 = df1.drop_duplicates()
    df1 = df1.sort_values([sortkey],ascending=sortasc)
    # df1['SampleType'] = [ s.split(':',1)[-1] for s in df1['SampleType'] ]
    df1['SampleType'] = [ fixcelltype(s) for s in df1['SampleType'] ]
    df1['Recall'] = 100*df1['Recall']
    fig, ax1 = subplots()
    h1 = ax1.bar(df1['SampleType'], df1['Neg10LogPVal'],label='-10log10(p-Value)')
    h2 = ax1.plot(df1['SampleType'], df1['Recall'], 'o', color="red",label='% Recall')
    h3 = ax1.plot([-1,len(df1['SampleType'])],[50,50], 'g--', )
    axis('tight')
    ax1.set_ylim([0,100])
    # legend(loc='best')
    if not hidexticks:
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        subplots_adjust(bottom=0.5)
    else:
        xticks([])
    show()

def fixgenes(s,focus):
    gs = s.split(',')
    if len(gs) == 1:
        return gs[0] + " Alone"
    gs.remove(focus)
    if gs[0].startswith('ST6GALNAC'):
        return "XXXXX"
    return "w/ "+gs[0]

def otherplot(df,focus,sampletypes,sortkey,sortasc=True,hideticks=False):
    df1 = df
    genes = [ g for g in set(df1['Genes']) if focus in g.split(',') ]
    df1 = df1[df1['Genes'].isin(genes)]
    df1 = df1[df1['SampleType'].isin(sampletypes)]
    df1 = df1.drop(columns=['GeneSet'])
    df1 = df1.drop_duplicates()
    df1 = df1.sort_values([sortkey],ascending=sortasc)
    df1['Genes'] = [ fixgenes(s,focus) for s in df1['Genes'] ]
    df1 = df1[df1['Genes'] != 'XXXXX']
    df1['Recall'] = 100*df1['Recall']
    fig, ax1 = subplots()
    h1 = ax1.bar(df1['Genes'], df1['Neg10LogPVal'],label='-10log10(p-Value)')
    h2 = ax1.plot(df1['Genes'], df1['Recall'], 'o', color="red",label='% Recall')
    h3 = ax1.plot([-1,len(df1['SampleType'])],[50,50], 'g--', )
    ax1.set_ylim([0,100])
    if not hideticks:
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        subplots_adjust(bottom=0.5)
    else:
        xticks([])
    show()

gene = 'ST6GALNAC3'
otherplot(df,gene,endomuscularis,'Neg10LogPVal',False)
# baseplot(df,gene,celltypesandtissues)
# baseplot(df,gene,celltypeintissue,'Neg10LogPVal',False,True)
# baseplot(df,gene,mucouscell,'Neg10LogPVal',False)
