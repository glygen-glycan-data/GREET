#!/bin/env python3.12

import sys, os
import pandas as pd

df1 = pd.read_csv(sys.argv[1],sep='\t').pivot_table("Z Score","Tissue","Gene")
df2 = pd.read_csv(sys.argv[2],sep='\t').pivot_table("Z Score","Tissue","Gene")

def transform1(inds):
    newinds = []
    for i in inds:
        if ' in [' in i:
            i = " ".join(i.split(" in "))
        elif i == "Adipocyte in Breast":
            i = "Adipocyte [Breast]"
        elif i == "Immune (NK cell) in Heart":
            i = "Immune (NK cell) [Heart]"
        elif i == "Immune (mast cell) in Esophagus muscularis":
            i = "Immune (mast cell) [Esophagus muscularis]"
        newinds.append(i)
    return newinds

df1 = df1.reset_index()
df1['Tissue'] = transform1(df1['Tissue'])
df1 = df1.set_index("Tissue")

def transform2(inds):
    newinds = []
    for i in inds:
        if i.startswith("Tissue: "):
            i = i.split(None,1)[1]
        elif i.startswith("Cell-type: "):
            i = i.split(None,1)[1]
        newinds.append(i)
    return newinds

df2 = df2.reset_index()
df2['Tissue'] = transform2(df2['Tissue'])
df2 = df2.set_index("Tissue")

rows1 = list(df1.index)
rows2 = list(df2.index)

combined_rows = list(set(rows1).intersection(set(rows2)))

df1 = df1.loc[combined_rows]
df2 = df2.loc[combined_rows]

df = df1.merge(df2, left_index=True, right_index=True)

import matplotlib
from matplotlib.pyplot import *
import seaborn as sns

genes = ['B4GALNT3, CHST9','B4GALNT4, CHST9','B4GALNT3, CHST8','B4GALNT4, CHST8','CHST9','CHST8','B4GALNT3', 'B4GALNT4']

genes = list(filter(None,"""
B4GALNT3, ST6GALNAC3
B4GALNT4, ST6GALNAC3
B4GALNT3, ST6GALNAC5
B4GALNT4, ST6GALNAC5
ST6GALNAC3
ST6GALNAC5
B4GALNT3
B4GALNT4
""".splitlines()))



df1 = df[genes]
df1 = df1.reset_index()
df1 = df1.melt(id_vars=["Tissue"],value_vars=genes)
df1['Gene(s)'] = df1['Gene']
df1['Cell-type'] = df1['Tissue']
df1['Z-Score'] = df1['value']


ax = sns.displot(df1,x='Z-Score', hue="Gene(s)", kind="ecdf", complementary=True, stat='percent')
sns.move_legend(ax, "upper right")
plot((5,5),(0,100),'--',color='0.6')
ylim([0,20])
ylabel('Percent of Cell-types')
tight_layout()
savefig("ST6GALNAC3.png",dpi=150)
# show()


