
import pandas as pd
import numpy as np
import random

from collections import defaultdict

class DataIO(object):
    def __init__(self,config):
        self.config = config
        self.enzymes = config.get_enzymes()

    def glyco_enzymes(self):
        return self.enzymes.all_enzymes()

    def glyco_enzyme_genesets(self):
        seen = set()
        for name,geneset in self.enzymes.genesets():
            if len(set(geneset).intersection(self.glyco_genes)) == len(geneset):
                names = self.enzymes.geneset_names(geneset)
                key = ",".join(sorted(geneset))
                if key not in seen:
                    yield names,geneset
                    seen.add(key)

    def sample_non_glyco_genes(self,n):
        return random.sample(self.non_glyco_genes,n)

    def create_genebased_dataframe(self,genes,samples):
        df = self.gene_restricted_df(genes)
        df['Class'] = 0
        ind = self.samples_to_index(samples)
        df.loc[ind,'Class'] = 1
        return df

    def random_dfs(self,ngenes,randcnt):
        for i in range(randcnt):
            yield self.gene_restricted_df(self.sample_non_glyco_genes(ngenes))
 
    def create_random_dataframes(self,ngenes,samples,randcnt):
        for df in self.random_dfs(ngenes,randcnt):
            df['Class'] = 0
            ind = self.samples_to_index(samples)
            df.loc[ind,'Class'] = 1
            yield df

    def gene_restriction_plus_random_dfs(self,genes,randcnt):
        dfs = [self.gene_restricted_df(genes)]
        for i in range(randcnt):
            dfs.append(self.gene_restricted_df(self.sample_non_glyco_genes(len(genes))))
        return dfs

    def create_dataframes(self,genes,samples,randcnt):
        dfs = []
        for df in self.gene_restriction_plus_random_dfs(genes,randcnt):
            df['Class'] = 0
            ind = self.samples_to_index(samples)
            df.loc[ind,'Class'] = 1
            dfs.append(df)
        return dfs

import GTExData

class GTExTissueRNASeq(DataIO):

    def read(self,filename):
        self.data = GTExData.GTExData(filename)
        self.glyco_genes = sorted(set(self.data.genes()).intersection(set(self.glyco_enzymes())))
        self.non_glyco_genes = sorted(set(self.data.genes()).difference(self.glyco_genes))

        self.sampledf = self.data.df[['Tissue']]
        self.tissues = set()
        tissues = set()
        freq = defaultdict(set)
        for t in self.data.tissues:
            tissues.add(t)
            if ' - ' not in t:
                continue
            base,rest = t.split(' - ',1)
            freq[base].add(t)
        for base in freq:
            if len(freq[base]) > 1:
                tissues.add(base)
        for t in tissues:
            ind = self.tissue_toindex(t)
            if sum(ind) >= self.config.min_samples_per_sampletype():
                self.tissues.add(t)

    def tissue_toindex(self,t):
        if ' - ' in t:
            return (self.sampledf['Tissue']==t)
        ts = set()
        for ti in self.data.tissues:
            if t == ti.split(' - ',1)[0]:
                ts.add(ti)
        ind = pd.Series(False,index=self.sampledf.index)
        for ti in ts:
            ind = (ind | (self.sampledf['Tissue'] == ti))
        return ind
            
    def samples_to_index(self,samples):
        return self.tissue_toindex(samples)

    def topredict(self):
        for t in sorted(self.tissues):
            yield t,t

    def gene_restricted_df(self,genes):
        return self.data.df.loc[:,list(genes)]

import scanpy

class GTExCelltypeSCRNASeq(DataIO):
    
    def read(self,filename):
        self.data = scanpy.read_h5ad(filename)
        self.glyco_genes = sorted(set(self.data.var_names).intersection(set(self.glyco_enzymes())))
        self.non_glyco_genes = sorted(set(self.data.var_names).difference(self.glyco_genes))

        self.sampledf = self.data.obs[['Broad cell type','Tissue']]
        celltypes = set(self.sampledf.itertuples(index=False, name=None))
        freq = defaultdict(set)
        for ct,ti in celltypes:
            if ct == 'Unknown':
                continue
            freq[('celltype',ct)].add((ct,ti))
            freq[('tissue',ti)].add((ct,ti))
        for k,v in sorted(freq.items(),key=lambda t: -len(t[1])):
            if len(v) > 1:
                celltypes.add(k)
        self.celltypes = set()
        for ct in celltypes:
            if ct[0] == 'Unknown':
                continue
            ind = self.celltype_toindex(ct)
            if sum(ind) >= self.config.min_samples_per_sampletype():
                self.celltypes.add(ct)
        pass

    def cttoname(self,ct):
        if ct[0] not in ('celltype','tissue'):
            return "%s in %s"%ct
        return "%s:%s"%(ct[0].title(),ct[1])

    def topredict(self):
        for ct in sorted(self.celltypes):
            yield self.cttoname(ct),ct

    def samples_to_index(self,samples):
        return self.celltype_toindex(samples)

    def celltype_toindex(self,ctin):
        cts = set()
        for ct in set(self.sampledf.itertuples(index=False, name=None)):
            if ctin[0] not in ('celltype','tissue') and ctin == ct:
                cts.add(ct)
            elif ctin[0] == 'celltype' and ctin[1] == ct[0]:
                cts.add(ct)
            elif ctin[0] == 'tissue' and ctin[1] == ct[1]:
                cts.add(ct)
        ind = pd.Series(False,index=self.sampledf.index)
        for ct in cts:
            ind = (ind | ((self.sampledf['Broad cell type'] == ct[0]) & (self.sampledf['Tissue'] == ct[1])))
        return ind

    def gene_restricted_df(self,genes):
        colind = self.data.var_names.isin(genes)
        df = pd.DataFrame.sparse.from_spmatrix(self.data.X[:,colind], 
                                               index=self.data.obs_names,
                                               columns=self.data.var_names[colind])
        return df

class GTExCelltypeSCRNASeqThreshold(GTExCelltypeSCRNASeq):

    def __init__(self,*args,**kwargs):
        self.threshold = kwargs.get('threshold',0)
        if 'threshold' in kwargs:
            del kwargs['threshold']
        self.mode = kwargs.get('mode','GT')
        if 'mode' in kwargs:
            del kwargs['mode']
        assert self.mode in ('GE','GT')
        super().__init__(*args,**kwargs)

    def transform(self,df):
        df1 = df.copy()
        if mode == "GT":
            ind = (df1 > self.threshold)
        elif mode == "GE":
            ind = (df1 >= self.threshold)
        df1[ind] = 1
        df1[~ind] = 0
        return df1

    def gene_restricted_df(self,genes):
        return self.transform(super().gene_restricted_df(genes))
